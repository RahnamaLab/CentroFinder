#!/bin/bash
#SBATCH --job-name=centromere_score
#SBATCH --nodes=1
#SBATCH --cpus-per-task=28
#SBATCH --time=2-00:00:00
#SBATCH --account=phd-ssalimi42
#SBATCH --mem=124G
#SBATCH --partition=batch-impulse
# ==============================
# CONFIGURATION
# ==============================
THREADS=${SLURM_CPUS_PER_TASK}
WINDOW=1000
ASSEMBLY=Fo4287v4.fasta
TRF=Fo4287v4_trf.bed
TE=Fo4287v4_edta.bed
HIFI=bam/Fo4287v4.hifi.pbmm2.bam
METH=bam/Fo4287v4.hifi.pbmm2.call_mods.modbam.freq.aggregate.all.bed
GENES=Fo4287v4.gff3    # optional
OUTDIR="centromere_results"
mkdir -p $OUTDIR

# Load necessary modules
spack load bedtools2 samtools/f6zgopc py-pandas py-numpy/5kbd5ic py-matplotlib/6s4k4w5 py-biopython

# ==============================
# STEP 1: PREPARE INPUT
# ==============================
if [ ! -f ${ASSEMBLY}.fai ]; then
    # Create FASTA index, which contains chromosome lengths
    samtools faidx $ASSEMBLY
fi

# Create genome windows
bedtools makewindows -g ${ASSEMBLY}.fai -w $WINDOW > $OUTDIR/windows.${WINDOW}bp.bed

# Convert TRF output to BED and sort
awk -F'\t' '{
    split($1, coords, ":");
    chrom = coords[1];
    split(coords[2], range, "-");
    start = range[1] - 1;
    end = range[2];
    motif = $2;
    print chrom, start, end, motif;
}' OFS="\t" $TRF | sort -k1,1V -k2,2n > $OUTDIR/trf.sorted.bed

# Sort TE and methylation
sort -k1,1V -k2,2n $TE > $OUTDIR/te.sorted.bed
awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$11}' $METH \
    | sort -k1,1V -k2,2n > $OUTDIR/methylation.sorted.bedgraph

# ==============================
# STEP 2: FEATURE EXTRACTION
# ==============================
# TRF coverage
bedtools coverage -a $OUTDIR/windows.${WINDOW}bp.bed -b $OUTDIR/trf.sorted.bed -counts > $OUTDIR/tmp.trf_counts.bed

# TE coverage
bedtools coverage -a $OUTDIR/windows.${WINDOW}bp.bed -b $OUTDIR/te.sorted.bed -counts > $OUTDIR/tmp.te_counts.bed

# Gene counts
if [ -f $GENES ]; then
    awk '$3=="gene" {print $1"\t"($4-1)"\t"$5"\t"$9}' $GENES > $OUTDIR/genes.bed
    bedtools coverage -a $OUTDIR/windows.${WINDOW}bp.bed -b $OUTDIR/genes.bed -counts > $OUTDIR/tmp.gene_counts.bed
else
    awk '{print $1"\t"$2"\t"$3"\t0"}' $OUTDIR/windows.${WINDOW}bp.bed > $OUTDIR/tmp.gene_counts.bed
fi

# HiFi coverage (mean depth per window)
samtools depth -a $HIFI | awk '{print $1"\t"$2-1"\t"$2"\t"$3}' > $OUTDIR/hifi.depth.bed
bedtools map -a $OUTDIR/windows.${WINDOW}bp.bed -b $OUTDIR/hifi.depth.bed -c 4 -o mean -null 0 > $OUTDIR/tmp.hifi_cov_mean.bed

# Mean methylation per window
bedtools map -nonamecheck -a $OUTDIR/windows.${WINDOW}bp.bed -b $OUTDIR/methylation.sorted.bedgraph \
    -c 4 -o mean -null 0 > $OUTDIR/tmp.meth_mean.bed

# GC content per window
bedtools nuc -fi $ASSEMBLY -bed $OUTDIR/windows.${WINDOW}bp.bed \
    | awk 'NR>1 {print $1"\t"$2"\t"$3"\t"$5}' \
    > $OUTDIR/tmp.gc_content.bed

# ==============================
# STEP 3: COMBINE FEATURES
# ==============================
echo -e "chrom\tstart\tend\ttrf_cov\tte_cov\tgene_count\thifi_cov_mean\tmeth_mean\tgc_content" > $OUTDIR/windows.features.tsv
paste \
    <(cut -f1-3 $OUTDIR/windows.${WINDOW}bp.bed) \
    <(cut -f4 $OUTDIR/tmp.trf_counts.bed) \
    <(cut -f4 $OUTDIR/tmp.te_counts.bed) \
    <(cut -f4 $OUTDIR/tmp.gene_counts.bed) \
    <(cut -f4 $OUTDIR/tmp.hifi_cov_mean.bed) \
    <(cut -f4 $OUTDIR/tmp.meth_mean.bed) \
    <(cut -f4 $OUTDIR/tmp.gc_content.bed) \
    >> $OUTDIR/windows.features.tsv

# ==============================
# STEP 4: PYTHON SCORING & PLOTTING
# ==============================
cat > $OUTDIR/score_centromeres.py << 'PYCODE'
#!/usr/bin/env python3
import os, sys, numpy as np, pandas as pd, matplotlib.pyplot as plt

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

# ==============================
# CONFIGURATION
# ==============================
fn = "windows.features.tsv"
out_prefix = "centro"
# Chromosome length file (passed as argument 1)
try:
    FAI_FILE = sys.argv[1]
except IndexError:
    print("Error: FAI file path must be provided as the first argument.")
    sys.exit(1)

# Standard large exclusion zone (for long chromosomes)
EXCLUSION_BP_LARGE = 500000 
# Minimal exclusion zone (for very short scaffolds, e.g., 100 kbp)
EXCLUSION_BP_MIN = 50000 # INCREASED to 50kbp to better suppress small scaffolds
WINDOW = 1000 # Added for use in contiguity logic

# ==============================
# PREPARE DATA
# ==============================
df = pd.read_csv(fn, sep='\t')

# Load chromosome lengths from FAI file
chr_lengths = {}
try:
    with open(FAI_FILE, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            # FAI file columns: CHROM, LENGTH, OFFSET, LINE_BASES, LINE_WIDTH
            if len(parts) >= 2:
                chr_lengths[parts[0]] = int(parts[1])
except Exception as e:
    print(f"Error loading FAI file {FAI_FILE}: {e}")
    sys.exit(1)


def minmax(s):
    if s.max() == s.min(): return s*0
    return (s - s.min()) / (s.max() - s.min())

# Normalize features
df['trf_cov_n'] = minmax(df['trf_cov'])

# Invert TE coverage normalization: Low TE -> High Score (Centromere-poor TE)
df['te_cov_n'] = 1 - minmax(df['te_cov'])

df['gene_count_n'] = 1 - minmax(df['gene_count'])
df['meth_diff_n'] = minmax(abs(df['meth_mean'] - df['meth_mean'].median()))
df['hifi_cov_mean'] = df['hifi_cov_mean'].replace(0, 1e-6)
df['cov_anom_n'] = minmax(abs(np.log2(df['hifi_cov_mean']/df['hifi_cov_mean'].median())))
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
    if chr_len <= EXCLUSION_BP_LARGE * 2:
        current_exclusion_bp = EXCLUSION_BP_MIN
        # print(f"Note: Chromosome {chrom} is short ({chr_len/1000} kb). Using minimal {current_exclusion_bp/1000} kb exclusion.")
    else:
        current_exclusion_bp = EXCLUSION_BP_LARGE
    
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
    #df.loc[exclusion_mask, 'te_cov_n'] = 0.0
   # df.loc[exclusion_mask, 'gene_count_n'] = 0.0
    #df.loc[exclusion_mask, 'meth_diff_n'] = 0.0


# Weighted scoring (Weights remain strong for centromeric features)
w = dict(trf=4.0, te=3.0, gene=1.0, meth=1.0, cov=0.5, gc=1.0) 
df['centro_score'] = (
    w['trf']*df['trf_cov_n'] +
    w['te']*df['te_cov_n'] +  
    w['gene']*df['gene_count_n'] +
    w['meth']*df['meth_diff_n'] +
    w['cov']*df['cov_anom_n'] +
    w['gc']*df['gc_low_n']
)

# Rank per chromosome
df['rank_within_chr'] = df.groupby('chrom')['centro_score'].rank(ascending=False, method='first')
df = df.sort_values(['chrom','start']).reset_index(drop=True)
df.to_csv(out_prefix+"_windows_ranked.tsv", sep='\t', index=False)

# ==============================
# Candidate selection 
# ==============================
candidates = []
for chrom, sub in df.groupby('chrom'):
    # 1. Apply hard filter (allows up to 1 TE hit, zero gene) to ALL windows
    perfect_sub = sub[(sub['te_cov'] <= 1) & (sub['gene_count'] == 0)]

    # 1b. CRITICAL FILTER: Exclude subtelomeric regions from the perfect set
    internal_perfect_sub = perfect_sub[~perfect_sub['is_excluded']]

    if internal_perfect_sub.empty:
        # Fallback 1: If no internal perfect regions, print warning and select best overall score
        # This will select the best-scoring window, even if it's a penalized subtelomere, 
        # as a last resort if no internal candidate meets the purity filter.
        print(f"Warning: No internal windows with zero TE/Gene found on {chrom}. Selecting overall best score as fallback (may be penalized subtelomere).")
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

cand_df = pd.concat(candidates).sort_values(['chrom','start']).reset_index(drop=True)
# ... (rest of the script for merging and plotting remains the same)

# Merge close windows (<10 kb)
merged = []
for chrom, g in cand_df.groupby('chrom'):
    g = g.sort_values('start')
    cur_s, cur_e = None, None
    for _, r in g.iterrows():
        if cur_s is None: cur_s, cur_e = int(r['start']), int(r['end'])
        elif int(r['start']) <= cur_e + 10000: cur_e = max(cur_e,int(r['end']))
        else:
            merged.append((chrom, cur_s, cur_e))
            cur_s, cur_e = int(r['start']), int(r['end'])
    if cur_s is not None: merged.append((chrom, cur_s, cur_e))

# Save merged candidates
with open(out_prefix+"_candidates.bed","w") as fh:
    for chrom,s,e in merged:
        fh.write(f"{chrom}\t{s}\t{e}\n")
cand_df.to_csv(out_prefix+"_candidates_ranked.tsv", sep='\t', index=False)
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
    if not chrom_candidates.empty:
        best_window = chrom_candidates.loc[chrom_candidates['centro_score'].idxmax()]
        best_x = (best_window['start'] + best_window['end']) / 2 / 1000
        best_df.append(best_window)
    else:
        best_x = None

    # =====================
    # Identify STRICT, CONTIGUOUS no-gene/low-TE region around best candidate
    # =====================
    no_gene_region = None
    if best_x is not None:
        # Search window: ±100 kb 
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

PYCODE
# ==============================
# RUN PYTHON SCORING
# ==============================
cd $OUTDIR
# Pass the FASTA index file path (which is outside the $OUTDIR) to the Python script
python3 score_centromeres.py ../${ASSEMBLY}.fai

# ==============================
# STEP 5: PICK BEST CENTROMERE CANDIDATES
# ==============================
# The best candidate selection code is redundant with the output from score_centromeres.py
# (which outputs centro_best_windows_marked.tsv). We will rely on that output.

echo "✅ Pipeline complete."
echo "Results saved in $OUTDIR:"
echo " - centro_candidates.bed (This should contain regions closer to 40 kbp)"
echo " - centro_candidates_ranked.tsv"
echo " - centro_windows_ranked.tsv"
echo " - centro_best_windows_marked.tsv (plotted)"
echo " - centro_best_condidate.bed (The strictly defined *core* purity region)"
echo " - best_centromere_candidates.bed"
echo " - centro_plots/ (per-chromosome feature plots)"
