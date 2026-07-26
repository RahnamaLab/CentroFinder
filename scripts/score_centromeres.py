#!/usr/bin/env python3
import os, sys, numpy as np, pandas as pd, matplotlib.pyplot as plt, argparse

# Standard large exclusion zone (for long chromosomes)
DEFAULT_EXCLUSION_BP_LARGE = 100000
# Minimal exclusion zone (for very mini chromosomes, e.g., 100 kbp)
DEFAULT_EXCLUSION_BP_MIN = 10000
DEFAULT_MINI_CHR_CUTOFF = 500000 # Added script fallback baseline
DEFAULT_WINDOW = 1000

# Weighted Scoring Dict Values
DEFAULT_TRF=4.0
DEFAULT_TE=8.0
DEFAULT_GENE=1.0
DEFAULT_METH=1.0
DEFAULT_COV=1.0
DEFAULT_GC=1.0

def parse_args():
    parser = argparse.ArgumentParser(
        description="Score centromere candidate windows strictly within gene-free islands."
    )
    parser.add_argument("--features", required=True, help="Path to windows.features.tsv file.")
    parser.add_argument("--fai", required=True, help="Path to FASTA index (.fai) file.")
    parser.add_argument("--outdir", required=True, help="Output directory for centromere scoring results.")
    parser.add_argument(
        "--platform",
        choices=["nanopore", "pacbio"],
        default="nanopore",
        help="Sequencing platform flag (kept for compatibility)."
    )
    parser.add_argument("--trf", type=float, default=DEFAULT_TRF, help="Weight for TRF feature (default: %(default)s).")
    parser.add_argument("--te", type=float, default=DEFAULT_TE, help="Weight for TE feature (default: %(default)s).")
    parser.add_argument("--gene", type=float, default=DEFAULT_GENE, help="Weight for GENE feature (default: %(default)s).")
    parser.add_argument("--meth", type=float, default=DEFAULT_METH, help="Weight for METH feature (default: %(default)s).")
    parser.add_argument("--cov", type=float, default=DEFAULT_COV, help="Weight for COV feature (default: %(default)s).")
    parser.add_argument("--gc", type=float, default=DEFAULT_GC, help="Weight for GC feature (default: %(default)s).")
    parser.add_argument("--exclusion-bp-large", type=int, default=DEFAULT_EXCLUSION_BP_LARGE, help="Large exclusion zone (default: %(default)s).")
    parser.add_argument("--exclusion-bp-min", type=int, default=DEFAULT_EXCLUSION_BP_MIN, help="Minimum exclusion zone (default: %(default)s).")
    parser.add_argument("--mini-chr-cutoff", type=int, default=DEFAULT_MINI_CHR_CUTOFF, help="Threshold length below which a mini chromosome is categorized as short.")
    parser.add_argument("--window", type=int, default=DEFAULT_WINDOW, help="Window size in bp (default: %(default)s).")
    parser.add_argument("--accum-windows", type=int, default=50, help="Number of windows to aggregate for regional accumulation (default: 50 = ~50kb).")
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
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    out_prefix = os.path.join(args.outdir, "centro")

    # ==============================
    # PREPARE DATA
    # ==============================
    df = pd.read_csv(args.features, sep='\t')

    try:
        chr_lengths = {}
        with open(args.fai, 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 2:
                    chr_lengths[parts[0]] = int(parts[1])
    except Exception as e:
        print(f"Error loading FAI file {args.fai}: {e}")
        sys.exit(1)

    # Normalize features 
    df['trf_cov_n'] = df.groupby('chrom')['trf_cov'].transform(minmax)
    df['te_cov_n'] = df.groupby('chrom')['te_cov'].transform(minmax)             
    df['gene_count_n'] = 1 - df.groupby('chrom')['gene_count'].transform(minmax) 
    
    # Absolute difference calculations
    df['meth_diff'] = df.groupby('chrom').apply(lambda x: abs(x['meth_mean'] - x['meth_mean'].median()), include_groups=False).reset_index(level=0, drop=True)
    df['meth_diff_n'] = df.groupby('chrom')['meth_diff'].transform(minmax)
    
    df['hifi_cov_mean'] = df['hifi_cov_mean'].replace(0, 1e-6)
    df['cov_anom'] = df.groupby('chrom').apply(lambda x: abs(np.log2(x['hifi_cov_mean'] / x['hifi_cov_mean'].median())), include_groups=False).reset_index(level=0, drop=True)
    df['cov_anom_n'] = df.groupby('chrom')['cov_anom'].transform(minmax)
    
    df['gc_low_n'] = 1 - df.groupby('chrom')['gc_content'].transform(minmax)
    df['is_excluded'] = False
    # ==========================================
    # Dynamic exclusion zone masking (Updated)
    # ==========================================
    for chrom in df['chrom'].unique():
        if chrom not in chr_lengths:
            print(f"Warning: Length for chromosome {chrom} not found in FAI file. Skipping.")
            continue

        chr_len = chr_lengths[chrom]
        
        # MODIFICATION 1: Evaluates scaffold length via configurable argument flag
        if chr_len <= args.mini_chr_cutoff:
            current_exclusion_bp = args.exclusion_bp_min
        else:
            current_exclusion_bp = args.exclusion_bp_large
        
        # MODIFICATION 2: Drop the chromosome if target exclusion zones clash or overflow 
        if chr_len <= (current_exclusion_bp * 2):
            print(f"Warning: Skipping {chrom} because its length ({chr_len} bp) is too short "
                  f"to accommodate symmetrical subtelomeric exclusion zones ({current_exclusion_bp} bp each).")
            df.loc[df['chrom'] == chrom, 'is_excluded'] = True
            continue

        # Map positions safely if length tests pass
        exclusion_mask = (df['chrom'] == chrom) & ((df['start'] < current_exclusion_bp) | (df['end'] > chr_len - current_exclusion_bp))
        df.loc[exclusion_mask, 'is_excluded'] = True
    # Calculate weighted final metrics
    w = dict(trf=args.trf, te=args.te, gene=args.gene, meth=args.meth, cov=args.cov, gc=args.gc)
    df['centro_score'] = (
        w['trf'] * df['trf_cov_n'] + w['te'] * df['te_cov_n'] +
        w['gene'] * df['gene_count_n'] + w['meth'] * df['meth_diff_n'] +
        w['cov'] * df['cov_anom_n'] + w['gc'] * df['gc_low_n']
    )
    
    # Zero out scores strictly inside the exclusion zones
    df.loc[df['is_excluded'], 'centro_score'] = 0.0

    df['rank_within_chr'] = df.groupby('chrom')['centro_score'].rank(ascending=False, method='first')
    df = df.sort_values(['chrom', 'start']).reset_index(drop=True)
    df.to_csv(out_prefix + "_windows_ranked.tsv", sep='\t', index=False)

    # ==============================
    # Candidate selection (Gene-Free Island Boundaries)
    # ==============================
    candidates = []
    
    for chrom, sub in df.groupby('chrom'): 
        internal_sub = sub[~sub['is_excluded']].copy().sort_values('start')
        if internal_sub.empty:
            print(f"Warning: Chromosome {chrom} has no windows outside exclusion zone. Skipping.")
            continue
        
        # 1. Calculate the rolling window profile across the chromosome
        internal_sub['accumulated_score'] = internal_sub['centro_score'].rolling(
            window=args.accum_windows, center=True, min_periods=1
        ).sum()
        
        max_accum_score = internal_sub['accumulated_score'].max()
        if max_accum_score == 0:
            continue
            
        best_accumulated_idx = internal_sub['accumulated_score'].idxmax()
        
        # 2. Map distinct "gene-free islands" across the chromosome
        internal_sub['is_gene_free'] = internal_sub['gene_count'] == 0
        internal_sub['island_id'] = (~internal_sub['is_gene_free']).cumsum()
        
        # Ensure our target anchor index sits safely inside a gene free window
        if not internal_sub.loc[best_accumulated_idx, 'is_gene_free']:
            gene_free_indices = internal_sub[internal_sub['is_gene_free']].index
            if gene_free_indices.empty:
                print(f"Warning: Chromosome {chrom} contains no gene-free windows at all. Skipping.")
                continue
            best_accumulated_idx = gene_free_indices[np.argmin(np.abs(gene_free_indices - best_accumulated_idx))]
            
        target_island_id = internal_sub.loc[best_accumulated_idx, 'island_id']
        
        # 3. Isolate the dataframe strictly down to this continuous gene-free domain
        island_windows = internal_sub[(internal_sub['island_id'] == target_island_id) & (internal_sub['is_gene_free'])].copy()
        
        # 4. Extract centromere boundaries ONLY within this island
        score_threshold = max_accum_score * 0.35 
        chrom_candidates = island_windows[island_windows['accumulated_score'] >= score_threshold]
        
        if not chrom_candidates.empty:
            candidates.append(chrom_candidates)

    if not candidates:
        print("Fatal Error: No candidates were found on any chromosome.")
        sys.exit(1)

    cand_df = pd.concat(candidates).sort_values(['chrom', 'start']).reset_index(drop=True)

    # Merge contiguous neighboring windows (tight tolerance to prevent expanding into flanking unmapped gaps)
    merged = []
    for chrom, g in cand_df.groupby('chrom'):
        g = g.sort_values('start')
        cur_s, cur_e = None, None
        for _, r in g.iterrows():
            if cur_s is None:
                cur_s, cur_e = int(r['start']), int(r['end'])
            elif int(r['start']) <= cur_e + 2000: # Tightened merging gap to ensure boundary integrity
                cur_e = max(cur_e, int(r['end']))
            else:
                merged.append((chrom, cur_s, cur_e))
                cur_s, cur_e = int(r['start']), int(r['end'])
        if cur_s is not None:
            merged.append((chrom, cur_s, cur_e))

    with open(out_prefix + "_candidates.bed", "w") as fh:
        for chrom, s, e in merged:
            fh.write(f"{chrom}\t{s}\t{e}\n")

    cand_df.to_csv(out_prefix + "_candidates_ranked.tsv", sep='\t', index=False)

    # ==============================
    # Plotting & Analysis Diagnostics
    # ==============================
    plotdir = out_prefix + "_plots"
    os.makedirs(plotdir, exist_ok=True)

    best_df = []
    gene_free_regions = [] 
    
    for chrom, g in df.groupby('chrom'):
        x = (g['start'] + g['end']) / 2 / 1000  
        internal = g[~g['is_excluded']].copy()
        
        if internal.empty:
            print(f"Skipping {chrom}: no internal windows after exclusion.")
            continue

        internal['accumulated_score'] = internal['centro_score'].rolling(
            window=args.accum_windows, center=True, min_periods=1
        ).sum()

        idx_glob = internal['accumulated_score'].idxmax()
        best_window = internal.loc[idx_glob]
        best_x = (best_window['start'] + best_window['end']) / 2 / 1000
        best_df.append(best_window)

        no_gene_region = None
        
        # Display the strictly gene-free centromere blocks in output plotting pipeline
        chrom_merged = [m for m in merged if m[0] == chrom]
        if chrom_merged:
            primary_m = max(chrom_merged, key=lambda item: item[2] - item[1])
            no_gene_region = (primary_m[1] / 1000, primary_m[2] / 1000)
            gene_free_regions.append((chrom, primary_m[1], primary_m[2], "best_candidate"))

        # Matplotlib visualization pipeline
        COLORS = {"TRF": "#0072B2", "TE": "#E69F00", "Meth": "#009E73", "Gene": "#56B4E9", "GC": "#CC79A7", "Coverage": "#D55E00", "Score": "#000000", "Centromere": "red"}
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 7), sharex=True)
        
        ax1.plot(x, g['trf_cov_n'], lw=0.75, label='TRF', color=COLORS["TRF"])
        ax1.plot(x, g['te_cov_n'], lw=0.75, label='TE', color=COLORS["TE"])
        ax1.plot(x, g['meth_diff_n'], lw=0.75, label='Meth', color=COLORS["Meth"])
        ax1.plot(x, g['gene_count_n'], lw=0.75, label='Gene', color=COLORS["Gene"])
        ax1.plot(x, g['gc_low_n'], lw=0.75, label='GC', color=COLORS["GC"])
        ax1.plot(x, g['cov_anom_n'], lw=0.75, label='Coverage', color=COLORS["Coverage"])
        ax1.set_ylabel("Normalized Feature Values")
        ax1.set_title(f"{chrom} - Features & Centromere Score")

        ax2.plot(x, g['centro_score'], lw=0.75, color=COLORS["Score"], label="Centromere Score")
        ax2.set_xlabel(f"{chrom} position (kb)")
        ax2.set_ylabel("Centromere Score")

        if no_gene_region is not None:
            ax2.axvspan(no_gene_region[0], no_gene_region[1], color=COLORS["Centromere"], alpha=0.3, label='Centromere Region')

        handles, labels = ax1.get_legend_handles_labels()
        h2, l2 = ax2.get_legend_handles_labels()
        handles.extend(h2); labels.extend(l2)

        fig.legend(handles, labels, loc='center right', fontsize='small')
        for ax in [ax1, ax2]: shade_excluded(ax, g)
        plt.tight_layout(rect=[0, 0, 0.85, 1])
        plt.savefig(f"{plotdir}/{chrom}_cen.pdf", dpi=150)
        plt.close()

    if best_df:
        best_per_chr = pd.DataFrame(best_df)[['chrom', 'start', 'end']]
        best_per_chr.to_csv(out_prefix + "_best_windows_marked.tsv", sep='\t', index=False, header=True)

    if gene_free_regions:
        gene_free_df = pd.DataFrame(gene_free_regions, columns=['chrom', 'start', 'end', 'name'])
        gene_free_df.to_csv(out_prefix + "_best_candidates.bed", sep='\t', index=False, header=False)

if __name__ == "__main__":
    main()
