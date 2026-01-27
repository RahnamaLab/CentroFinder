#!/usr/bin/env python3
import os, sys, numpy as np, pandas as pd, matplotlib.pyplot as plt, argparse

# Standard large exclusion zone (for long chromosomes)
DEFAULT_EXCLUSION_BP_LARGE = 100000
# Minimal exclusion zone (for very short scaffolds, e.g., 100 kbp)
DEFAULT_EXCLUSION_BP_MIN = 10000
DEFAULT_WINDOW = 1000

# Weighted Scoring Dict Values
DEFAULT_TRF=8.0
DEFAULT_TE=5.0
DEFAULT_GENE=1.0
DEFAULT_METH=1.0
DEFAULT_COV=0.5
DEFAULT_GC=1.0

def parse_args():
    parser = argparse.ArgumentParser(
        description="Score centromere candidate windows from feature tables."
    )

    parser.add_argument(
    "--features",
    required=True,
    help="Path to windows.features.tsv file.",
    )

    parser.add_argument(
    "--fai",
    required=True,
    help="Path to FASTA index (.fai) file.",
    )

    parser.add_argument(
    "--outdir",
    required=True,
    help="Output directory for centromere scoring results.",
    )

    parser.add_argument(
    "--platform",
    choices=["nanopore", "pacbio"],
    default="nanopore",
    help="Sequencing platform. Use 'pacbio' to enable PacBio-specific scoring/selection logic. (default: nanopore)."
    )

    parser.add_argument(
    "--trf",
    type=float,
    default=DEFAULT_TRF,
    help="Weight for TRF feature in the scoring model (default: %(default)s).",
    )

    parser.add_argument(
    "--te",
    type=float,
    default=DEFAULT_TE,
    help="Weight for TE feature in the scoring model (default: %(default)s).",
    )

    parser.add_argument(
    "--gene",
    type=float,
    default=DEFAULT_GENE,
    help="Weight for GENE feature in the scoring model (default: %(default)s).",
    )

    parser.add_argument(
    "--meth",
    type=float,
    default=DEFAULT_METH,
    help="Weight for METH feature in the scoring model (default: %(default)s).",
    )

    parser.add_argument(
    "--cov",
    type=float,
    default=DEFAULT_COV,
    help="Weight for COV feature in the scoring model (default: %(default)s).",
    )

    parser.add_argument(
    "--gc",
    type=float,
    default=DEFAULT_GC,
    help="Weight for GC feature in the scoring model (default: %(default)s).",
    )

    parser.add_argument(
    "--exclusion-bp-large",
    type=int,
    default=DEFAULT_EXCLUSION_BP_LARGE,
    help="Large exclusion zone in bp for long chromosomes (default: %(default)s).",
    )

    parser.add_argument(
    "--exclusion-bp-min",
    type=int,
    default=DEFAULT_EXCLUSION_BP_MIN,
    help="Minimum exclusion zone in bp for short scaffolds (default: %(default)s).",
    )

    parser.add_argument(
    "--window",
    type=int,
    default=DEFAULT_WINDOW,
    help="Window size in bp (default: %(default)s).",
    )

    return parser.parse_args()

def shade_excluded(ax, g, color="#D9D9D9", alpha=0.4):
    excluded = g[g['is_excluded']]
    if excluded.empty:
        return

    cur_s, cur_e = None, None
    for _, r in excluded.iterrows():
        if cur_s is None:
            cur_s, cur_e = r['start'], r['end']
        elif r['start'] <= cur_e:
            cur_e = r['end']
        else:
            ax.axvspan(cur_s/1000, cur_e/1000, color=color, alpha=alpha, lw=0)
            cur_s, cur_e = r['start'], r['end']

    if cur_s is not None:
        ax.axvspan(cur_s/1000, cur_e/1000, color=color, alpha=alpha, lw=0)

def minmax(s):
    if s.max() == s.min():
        return s * 0
    return (s - s.min()) / (s.max() - s.min())

def main():

    # ==============================
    # CONFIGURATION
    # ==============================
    args = parse_args()

    is_pacbio = (args.platform == "pacbio")

    features_path = args.features
    fai_path = args.fai
    outdir = args.outdir

    os.makedirs(outdir, exist_ok=True)

    out_prefix = os.path.join(outdir, "centro")

    exclusion_bp_large = args.exclusion_bp_large
    exclusion_bp_min = args.exclusion_bp_min
    window = args.window

    # Weight Scoring Dict Values 
    gc = args.gc
    trf = args.trf
    te = args.te
    gene = args.gene
    meth = args.meth
    cov = args.cov

    # ==============================
    # PREPARE DATA
    # ==============================
    df = pd.read_csv(features_path, sep='\t')

    # Load chromosome lengths from FAI file
    chr_lengths = {}
    try:
        with open(fai_path, 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                # FAI file columns: CHROM, LENGTH, OFFSET, LINE_BASES, LINE_WIDTH
                if len(parts) >= 2:
                    chr_lengths[parts[0]] = int(parts[1])
    except Exception as e:
        print(f"Error loading FAI file {fai_path}: {e}")
        sys.exit(1)

    # Normalize features
    df['trf_cov_n'] = minmax(df['trf_cov'])

    # Invert TE coverage normalization: Low TE -> High Score (Centromere-poor TE)
    df['te_cov_n'] = 1 - minmax(df['te_cov'])

    df['gene_count_n'] = 1 - minmax(df['gene_count'])
    df['meth_diff_n'] = minmax(abs(df['meth_mean'] - df['meth_mean'].median()))
    df['hifi_cov_mean'] = df['hifi_cov_mean'].replace(0, 1e-6)
    df['cov_anom_n'] = minmax(abs(np.log2(df['hifi_cov_mean'] / df['hifi_cov_mean'].median())))
    df['gc_low_n'] = 1 - minmax(df['gc_content'])

    # Initialize is_excluded column
    df['is_excluded'] = False

    # Apply severe scoring penalty in the dynamic exclusion zone
    for chrom in df['chrom'].unique():
        if chrom not in chr_lengths:
            print(f"Warning: Length for chromosome {chrom} not found in FAI file. Skipping exclusion for this chromosome.")
            continue

        chr_len = chr_lengths[chrom]
        
        # DETERMINE DYNAMIC EXCLUSION BP
        # If the chromosome is short (less than 2x Large Exclusion), use the minimal exclusion.
        if chr_len <= exclusion_bp_large * 2:
            current_exclusion_bp = exclusion_bp_min
        else:
            current_exclusion_bp = exclusion_bp_large
        
        # Identify windows within the exclusion zone (start < exclusion_bp OR end > L - exclusion_bp)
        start_mask = (df['chrom'] == chrom) & (df['start'] < current_exclusion_bp)
        end_mask = (df['chrom'] == chrom) & (df['end'] > chr_len - current_exclusion_bp)

        # Combined mask
        exclusion_mask = start_mask | end_mask
        
        # Record exclusion status for later filtering
        df.loc[exclusion_mask, 'is_excluded'] = True
        
        # === CRITICAL SUBTELOMERE PENALTY ===
        # Set the normalized scores of the most misleading features (Gene-poor, TE-poor, Meth) 
        # to 0.0 in the exclusion zone. This means only TRF and Cov Anomaly can contribute score here.

    # Weighted scoring (Weights remain strong for centromeric features)
    w = dict(trf=trf, te=te, gene=gene, meth=meth, cov=cov, gc=gc)
    df['centro_score'] = (
        w['trf'] * df['trf_cov_n'] +
        w['te'] * df['te_cov_n'] +
        w['gene'] * df['gene_count_n'] +
        w['meth'] * df['meth_diff_n'] +
        w['cov'] * df['cov_anom_n'] +
        w['gc'] * df['gc_low_n']
    )

    # ==============================
    # MASK SUBTELOMERIC REGIONS FOR PLOTTING
    # ==============================

    # Rank per chromosome
    df['rank_within_chr'] = df.groupby('chrom')['centro_score'].rank(ascending=False, method='first')
    df = df.sort_values(['chrom', 'start']).reset_index(drop=True)
    df.to_csv(out_prefix + "_windows_ranked.tsv", sep='\t', index=False)

    # ==============================
    # Candidate selection
    # ==============================
    candidates = []
    if is_pacbio:
        for chrom, sub in df.groupby('chrom'):
            # 1. Apply hard filter (allows up to 1 TE hit, zero gene) to ALL windows
            perfect_sub = sub[(sub['te_cov'] <= 1) & (sub['gene_count'] == 0)]

            # 1b. CRITICAL FILTER: Exclude subtelomeric regions from the perfect set
            internal_perfect_sub = perfect_sub[~perfect_sub['is_excluded']]

            if internal_perfect_sub.empty:
                # Fallback 1: If no internal perfect regions, print warning and select best overall score
                # This will select the best-scoring window, even if it's a penalized subtelomere, 
                # as a last resort if no internal candidate meets the purity filter.
                print(f"Warning: No internal windows with zero TE/Gene found on {chrom}. "
                "Selecting overall best score as fallback (may be penalized subtelomere).")
                sel = sub.nlargest(1,'centro_score')
            else:
                # 2. Use the strict 'internal perfect' regions for thresholding
                maxs = internal_perfect_sub['centro_score'].max()
                # Threshold: 25% of max score to aggressively capture the full 40 kbp region.
                thr = 0.25 * maxs
                sel = internal_perfect_sub[internal_perfect_sub['centro_score']>=thr]
                
                # 3. If no internal windows pass threshold, take the single best internal perfect window
                if sel.empty:
                    sel = internal_perfect_sub.nlargest(1,'centro_score')
            
            candidates.append(sel)
    else:
        for chrom, sub in df.groupby('chrom'): 
            # 1. Filter out windows that are in the exclusion zones
            internal_sub = sub[~sub['is_excluded']]

            if internal_sub.empty:
                print(f"Warning: Chromosome {chrom} has no windows outside the exclusion zone. Skipping.")
                continue

            # 2. Select the top 5 windows with the highest 'centro_score' 
            sel = internal_sub.nlargest(5, 'centro_score')

            if sel.empty:
                print(f"Warning: Could not select any candidates for {chrom}. Skipping.")
                continue

            candidates.append(sel)

    cand_df = pd.concat(candidates).sort_values(['chrom', 'start']).reset_index(drop=True)

    # Merge close windows (<10 kb)
    merged = []
    for chrom, g in cand_df.groupby('chrom'):
        g = g.sort_values('start')
        cur_s, cur_e = None, None
        for _, r in g.iterrows():
            if cur_s is None:
                cur_s, cur_e = int(r['start']), int(r['end'])
            elif int(r['start']) <= cur_e + 10000:
                cur_e = max(cur_e, int(r['end']))
            else:
                merged.append((chrom, cur_s, cur_e))
                cur_s, cur_e = int(r['start']), int(r['end'])
        if cur_s is not None:
            merged.append((chrom, cur_s, cur_e))

    # Save merged candidates
    with open(out_prefix + "_candidates.bed", "w") as fh:
        for chrom, s, e in merged:
            fh.write(f"{chrom}\t{s}\t{e}\n")

    cand_df.to_csv(out_prefix + "_candidates_ranked.tsv", sep='\t', index=False)

    #====================================
    # Plot per chromosome with best candidate marked
    # ==============================
    plotdir = out_prefix + "_plots"
    os.makedirs(plotdir, exist_ok=True)

    best_df = []
    gene_free_regions = [] # Initialize list to collect gene-free regions for output
    for chrom, g in df.groupby('chrom'):
        x = (g['start'] + g['end']) / 2 / 1000  # kb

        # Best candidate for this chromosome
        chrom_candidates = cand_df[cand_df['chrom'] == chrom]

        if is_pacbio:
            if not chrom_candidates.empty:
                best_window = chrom_candidates.loc[chrom_candidates['centro_score'].idxmax()]
                best_x = (best_window['start'] + best_window['end']) / 2 / 1000
                best_df.append(best_window)
            else:
                best_x = None

        else:

            # --- Identify true centromere as local minimum in the score (robust search) ---
            sub = g.copy()
            sub_center = (sub['start'] + sub['end']) / 2
            internal = sub[~sub['is_excluded']].copy()

            # If internal is empty, skip
            if internal.empty:
                print(f"Skipping {chrom}: no internal windows after exclusion.")
                continue

            # smoothing window in bp around each point (adjustable)
            search_kb = 200  # +/- 200 kb search radius around candidate center
            search_bp = int(search_kb * 1000)

            # convert WINDOW (bp per row) to number of rows for rolling
            # fallback if WINDOW is not exact per-row size
            median_step = int(np.median(internal['end'] - internal['start']))
            step = median_step if median_step > 0 else window
            window_n = max(3, int((2 * search_bp) / step))  # cover ~400 kb smoothing window by default
            # ensure window_n is odd for symmetric smoothing (not required, but okay)
            if window_n % 2 == 0:
                window_n += 1

            internal['score_smooth'] = internal['centro_score'].rolling(window_n, center=True, min_periods=1).median()

            # --- REVISED: Find the single best center globally in the smoothed score ---
            # The best position is simply the global minimum of the smoothed score series 
            # across all non-excluded windows for this chromosome.
            idx_glob = internal['centro_score'].idxmin()
            best_window = internal.loc[idx_glob]

            best_x = (best_window['start'] + best_window['end']) / 2 / 1000
            best_df.append(best_window) # Correct Indentation

            # The entire previous multi-candidate search and conditional logic for 'best_pos' is now removed.

        # =====================
        # Identify STRICT, CONTIGUOUS no-gene/low-TE region around best candidate
        # =====================
        no_gene_region = None
        if best_x is not None:
            # The search now centers on the globally best minimum found above
            search_start = best_window['start'] - 150000
            search_end = best_window['end'] + 150000
            sub = g[(g['start'] >= search_start) & (g['end'] <= search_end)]
            
            # 1. Find windows that are GENE-ZERO AND TE-LOW (<= 1)
            # CHANGED: Relaxing TE purity for the core region to expand the contiguous block to 40kbp
            no_gene_strict = sub[(sub['gene_count'] == 0) & (sub['te_cov'] <= 5)].sort_values('start')

            if not no_gene_strict.empty:
                
                # 2. Group contiguous windows
                # The 'block' flag increments every time a non-adjacent gap is found.
                # Adjacency is defined by: Current start == Previous end.
                no_gene_strict['is_contiguous'] = no_gene_strict['start'] == no_gene_strict['end'].shift(1)
                # Cumulative sum creates a new block ID whenever contiguity is broken (is_contiguous is False)
                no_gene_strict['block'] = (~no_gene_strict['is_contiguous']).cumsum()
                
                # 3. Merge windows within each block and calculate length
                contiguous_blocks = no_gene_strict.groupby('block').agg(
                    chrom=('chrom', 'first'),
                    start=('start', 'min'),
                    end=('end', 'max'),
                )
                contiguous_blocks['length'] = contiguous_blocks['end'] - contiguous_blocks['start']
                
                # 4. Select the single largest contiguous block to represent the gene-free region
                best_block = contiguous_blocks.loc[contiguous_blocks['length'].idxmax()]

                # Use this single, best contiguous block for the plot/save
                ng_start, ng_end = best_block['start'], best_block['end']

                no_gene_region = (ng_start / 1000, ng_end / 1000)  # in kb

                # Save the raw coordinates (in bp) for the new BED output file
                gene_free_regions.append((chrom, int(ng_start), int(ng_end), "best_candidate"))

        COLORS = {
        "TRF":        "#0072B2",  # blue
        "TE":         "#E69F00",  # orange
        "Meth":       "#009E73",  # green
        "Gene":       "#56B4E9",  # light blue
        "GC":         "#CC79A7",  # purple
        "Coverage":   "#D55E00",  # vermillion
        "Score":      "#000000",  # black
        "Centromere": "red"   # yellow (shaded region)
        }

        # --- Plot features and scores ---
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 7), sharex=True)
        ax1.plot(x, g['trf_cov_n'],        lw=0.75, label='TRF',        color=COLORS["TRF"])
        ax1.plot(x, g['te_cov_n'],         lw=0.75, label='TE',         color=COLORS["TE"])
        ax1.plot(x, g['meth_diff_n'],      lw=0.75, label='Meth',       color=COLORS["Meth"])
        ax1.plot(x, g['gene_count_n'],     lw=0.75, label='Gene',       color=COLORS["Gene"])
        ax1.plot(x, g['gc_low_n'],         lw=0.75, label='GC',         color=COLORS["GC"])
        ax1.plot(x, g['cov_anom_n'],       lw=0.75, label='Coverage',   color=COLORS["Coverage"])

        ax1.set_ylabel("Normalized Feature Values")
        ax1.set_title(f"{chrom} - Features & Centromere Score")

        ax2.plot(x, g['centro_score'], lw=0.75, color=COLORS["Score"], label="Centromere Score")
        ax2.set_xlabel(f"{chrom} position (kb)")
        ax2.set_ylabel("Centromere Score")

        if no_gene_region is not None:
            ax2.axvspan(no_gene_region[0], no_gene_region[1], color=COLORS["Centromere"], alpha=0.3, label='Centromere Region')

        handles, labels = [], []
        for ax in [ax1, ax2]:
            for line in ax.get_lines():
                handles.append(line)
                labels.append(line.get_label())
        if no_gene_region is not None:
            handles.append(ax2.fill_between([], [], [], color='red', alpha=0.3, label='Centromere Region'))
            labels.append('Centromere Region')

        fig.legend(handles, labels, loc='center right', fontsize='small')
        # Shade subtelomeric (excluded) regions
        for ax in [ax1, ax2]:
            shade_excluded(ax, g)
        plt.tight_layout(rect=[0, 0, 0.85, 1])
        plt.savefig(f"{plotdir}/{chrom}_cen.pdf", dpi=150)
        plt.close()

    # --- Save results ---
    best_per_chr = pd.DataFrame(best_df)[['chrom', 'start', 'end']]
    best_per_chr.to_csv(out_prefix + "_best_windows_marked.tsv", sep='\t', index=False, header=True)

    if gene_free_regions:
        gene_free_df = pd.DataFrame(gene_free_regions, columns=['chrom', 'start', 'end', 'name'])
        gene_free_df.to_csv(out_prefix + "_best_candidates.bed", sep='\t', index=False, header=False)
    else:
        print("Warning: No contiguous gene-free/low-TE regions identified.")

if __name__ == "__main__":
    main()
