# Created By: Sharon Colson
# Creation Date: 12/01/2025
# Last Modified: 01/06/2026

# To run this on TN Tech Univ HPC:
#     spack load trf snakemake graphviz
# sample run line:
#     snakemake -np
#         OR
#     snakemake --cores 12
#         OR
#     snakemake --use-conda --cores 12 results/sample/METH_PACBIO/sample.hifi.pbmm2.bam
#     snakemake --use-conda --conda-frontend conda --cores 12 results/sample/CENTROMERE_SCORING/sample_1000/centro_candidates.bed
#     snakemake --use-conda --conda-frontend conda --cores 12 results/sample/CENTROMERE_SCORING/sample_1000/centro_best_windows_marked.tsv
# Create DAG: snakemake --dag results/Guy11_chr1/CENTROMERE_SCORING/sample.1000.te.sorted.bed | dot -Tsvg > centromere_pipeline_sample.svg

import os
import subprocess

configfile: "config.yaml"

MATCH     = config["trf_params"]["match"]
MISMATCH  = config["trf_params"]["mismatch"]
DELTA     = config["trf_params"]["delta"]
PM        = config["trf_params"]["pm"]
PI        = config["trf_params"]["pi"]
MINSCORE  = config["trf_params"]["minscore"]
MAXPERIOD = config["trf_params"]["maxperiod"]
OPTIONS   = config["trf_params"]["options"]

NANOPORE_DIR = config["nanopore_dir"]
PACBIO_DIR   = config["pacbio_dir"]

SAMPLES_BY_PLATFORM = config["samples"] or {}

NANOPORE_DICT = SAMPLES_BY_PLATFORM.get("nanopore") or {}
PACBIO_DICT   = SAMPLES_BY_PLATFORM.get("pacbio") or {}

# Fail if neither platform has samples
if not NANOPORE_DICT and not PACBIO_DICT:
    raise ValueError(
        "Config error: at least one of 'samples.nanopore' or "
        "'samples.pacbio' must be defined with ≥1 sample."
    )

# Extract sample lists
NANOPORE_SAMPLES = list(NANOPORE_DICT.keys())
PACBIO_SAMPLES   = list(PACBIO_DICT.keys())

SAMPLES_LIST = NANOPORE_SAMPLES + PACBIO_SAMPLES

WINDOW = config["window"]

def is_nanopore(sample):
    return sample in NANOPORE_SAMPLES

def is_pacbio(sample):
    return sample in PACBIO_SAMPLES

def get_platform(sample):
    if sample in NANOPORE_SAMPLES:
        return "nanopore"
    if sample in PACBIO_SAMPLES:
        return "pacbio"
    raise ValueError(f"Sample '{sample}' not found under config['samples']['nanopore'] or ['pacbio'].")

def get_base_dir(sample):
    platform = get_platform(sample)
    if platform == "nanopore":
        return NANOPORE_DIR
    return PACBIO_DIR

def get_path_with_ext(wildcards, ext):
    base = get_base_dir(wildcards.sample)
    return f"{base}/{wildcards.sample}/{wildcards.sample}.{ext}"

def get_fasta(wildcards):
    return get_path_with_ext(wildcards, "fasta")

def get_fai(wildcards):
    return get_path_with_ext(wildcards, "fasta.fai")

def get_gff3(wildcards):
    return get_path_with_ext(wildcards, "gff3")

def get_fastq(wildcards):
    return get_path_with_ext(wildcards, "fastq")

def get_bam(wildcards):
    return get_path_with_ext(wildcards, "subreads.bam")

# Order for TRF and filename suffix
TRF_NUMERIC_VALUES = [MATCH, MISMATCH, DELTA, PM, PI, MINSCORE, MAXPERIOD]

# Add options for the TRF parameter string
TRF_PARAM_STRING = " ".join(str(num) for num in TRF_NUMERIC_VALUES)
if OPTIONS:
    TRF_PARAM_STRING += f" {OPTIONS}"

# Build the .dat file name suffix, e.g. ".2.7.7.80.10.50.dat"
TRF_SUFFIX = "." + ".".join(str(num) for num in TRF_NUMERIC_VALUES) + ".dat"

rule all:
    input:
        expand("results/{sample}/CENTROMERE_SCORING/{sample}_{window}/centro_windows_ranked.tsv", sample=SAMPLES_LIST, window=WINDOW),
        expand("results/{sample}/CENTROMERE_SCORING/{sample}_{window}/centro_best_windows_marked.tsv", sample=SAMPLES_LIST, window=WINDOW),
        expand("results/{sample}/CENTROMERE_SCORING/{sample}_{window}/centro_candidates.bed", sample=SAMPLES_LIST, window=WINDOW),
        expand("results/{sample}/CENTROMERE_SCORING/{sample}_{window}/centro_candidates_ranked.tsv", sample=SAMPLES_LIST, window=WINDOW),
        expand("results/{sample}/CENTROMERE_SCORING/{sample}_{window}/centro_best_candidates.bed", sample=SAMPLES_LIST, window=WINDOW),

#### TRF ####
rule TRF_run:
    input:
        fasta = get_fasta
    output:
        dat = "results/{sample}/TRF/{sample}.fasta" + TRF_SUFFIX
    log:
        "results/{sample}/TRF/logs/run_trf_{sample}.log"
    conda:
        "envs/trf.yaml"
    params:
        trf_param_string=TRF_PARAM_STRING
    script:
        "scripts/run_trf.py"

rule TRF_convert_to_bed:
    input:
        rules.TRF_run.output
    output:
        "results/{sample}/TRF/{sample}_trf.bed"
    log:
        "results/{sample}/TRF/logs/convert_trf_to_bed_{sample}.log"
    conda:
        "envs/trf.yaml"
    shell:
        r"""
        mkdir -p $(dirname {output}) $(dirname {log})
        python3 trf2bed.py \
            --dat {input} \
            --bed {output} \
            --tool repeatseq &> {log}
        """

##### EDTA #####
rule EDTA_cds:
    input:
        fasta = get_fasta,
        gff   = get_gff3
    output:
        cds = "results/{sample}/EDTA/{sample}_cds.fasta"
    log:
        "results/{sample}/EDTA/logs/edta_cds_{sample}.log"
    conda:
        "envs/edta.yaml"
    shell:
        r"""
        mkdir -p $(dirname {output.cds}) $(dirname {log})

        gffread -x {output.cds} -g {input.fasta} {input.gff} &> {log}
        """

rule EDTA_gff2bed:
    input:
        gff = get_gff3
    output:
        bed = "results/{sample}/EDTA/{sample}.bed"
    log:
        "results/{sample}/EDTA/logs/edta_bed_{sample}.log"
    conda:
        "envs/edta.yaml"
    shell:
        r"""
        mkdir -p $(dirname {output.bed}) $(dirname {log})

        gff2bed < {input.gff} > {output.bed} 2> {log}
        """

rule EDTA_run:
    input:
        fasta = get_fasta,
        cds = rules.EDTA_cds.output.cds,
        bed = rules.EDTA_gff2bed.output.bed
    output:
        edta_gff3 = "results/{sample}/EDTA/{sample}.fasta.mod.EDTA.TEanno.gff3"
    log:
        "results/{sample}/EDTA/logs/edta_run_{sample}.log"
    conda:
        "envs/edta.yaml"
    params:
        container_bin   = config["container"]["binary"],
        container_binds = ",".join(config["container"]["binds"]),
        container_env   = " ".join(
            f'{k}={v}' for k, v in config["container"]["env"].items()
        ),
        edta_sif          = config["edta"]["sif"],
        edta_species      = config["edta"]["species"],
        edta_overwrite    = config["edta"]["overwrite"],
        edta_sensitive    = config["edta"]["sensitive"],
        edta_anno         = config["edta"]["anno"],
        edta_force        = config["edta"]["force"],
    threads: config["cpus_per_task"]
    shell:
        r"""
        workdir="results/{wildcards.sample}/EDTA/edta"
        mkdir -p "$(dirname {output.edta_gff3})" "$(dirname {log})" "$workdir"

        sample={wildcards.sample}
        logfile="{log}"

        cp {input.fasta} "$workdir/{wildcards.sample}.fasta"
        cp {input.cds}   "$workdir/{wildcards.sample}_cds.fasta"
        cp {input.bed}   "$workdir/{wildcards.sample}.bed"

        (
        cd "$workdir"

          {params.container_bin} run \
            --bind {params.container_binds} \
            --env {params.container_env} \
            {params.edta_sif} EDTA.pl \
            --genome {wildcards.sample}.fasta \
            --cds {wildcards.sample}_cds.fasta \
            --exclude {wildcards.sample}.bed \
            --species {params.edta_species} \
            --overwrite {params.edta_overwrite} \
            --sensitive {params.edta_sensitive} \
            --anno {params.edta_anno} \
            --threads {threads} \
            --force {params.edta_force}

        ) &> "$logfile"

        cp "$workdir/{wildcards.sample}.fasta.mod.EDTA.TEanno.gff3" {output.edta_gff3}
        """

rule EDTA_bed:
    input:
        edta_gff = rules.EDTA_run.output.edta_gff3
    output:
        "results/{sample}/EDTA/{sample}_edta.bed"
    log:
        "results/{sample}/EDTA/logs/edta_bed_final_{sample}.log"
    conda:
        "envs/edta.yaml"
    shell:
        r"""
        mkdir -p $(dirname {log})

        awk 'BEGIN{{OFS="\t"}}
             !/^#/ {{
                match($9,/classification=([^;]+)/,a);
                 class=a[1];
                if (class=="") class="NA";
                print $1, $4-1, $5, $3, class, $6, $7
             }}' {input.edta_gff} > {output} 2> {log}
        """

###### Meth Nanopore #####
rule MN_minimap2:
    input:
        fasta = get_fasta,
        fastq = get_fastq
    output:
        sam = "results/{sample}/METH_NANOPORE/{sample}.sam"
    log:
        "results/{sample}/METH_NANOPORE/logs/minimap2_{sample}.log"
    conda:
        "envs/minimap.yaml"
    threads: config["cpus_per_task"]
    shell:
        r"""
        mkdir -p "$(dirname {log})"

        minimap2 -t {threads} -ax map-ont -Y {input.fasta} {input.fastq} > {output.sam} 2> {log}
        """

rule MN_samtools_sort:
    input:
        sam = rules.MN_minimap2.output.sam
    output:
        bam = "results/{sample}/METH_NANOPORE/{sample}_sorted.bam"
    log:
        "results/{sample}/METH_NANOPORE/logs/mn_samtools_sort_{sample}.log"
    conda:
        "envs/minimap.yaml"
    shell:
        r"""
        mkdir -p "$(dirname {log})"

        samtools sort -@ {threads} -o {output.bam} {input.sam} 2> {log}
        """


rule MN_samtools_index:
    input:
        bam = rules.MN_samtools_sort.output.bam
    output:
        bai = "results/{sample}/METH_NANOPORE/{sample}_sorted.bam.bai"
    log:
        "results/{sample}/METH_NANOPORE/logs/mn_samtools_index_{sample}.log"
    conda:
        "envs/minimap.yaml"
    shell:
        r"""
        mkdir -p "$(dirname {log})"

        samtools index {input.bam} 2> {log}
        """

rule MN_modbam2bed:
    input:
        fasta = get_fasta,
        bam   = rules.MN_samtools_sort.output.bam,
        bai   = rules.MN_samtools_index.output.bai
    output:
        bed = "results/{sample}/METH_NANOPORE/{sample}_methyl.bed"
    log:
        "results/{sample}/METH_NANOPORE/logs/mn_modbam2bed_{sample}.log"
    conda:
        "envs/minimap.yaml"
    shell:
        r"""
        mkdir -p "$(dirname {log})"

        modbam2bed {input.fasta} {input.bam} > {output.bed} 2> {log}
        """


#### Meth Pacbio ####
rule MP_ccsmeth_call_hifi:
    input:
        bam   = get_bam
    output:
        bam = "results/{sample}/METH_PACBIO/{sample}.hifi.bam"
    log:
        "results/{sample}/METH_PACBIO/logs/ccsmeth_call_hifi_{sample}.log"
    conda:
        "envs/ccsmeth.yaml"
    threads: config["cpus_per_task"]
    shell:
        r"""
        mkdir -p "$(dirname {log})"

        ccsmeth call_hifi \
           --subreads {input.bam} \
           --threads {threads} \
           --output {output.bam} &> {log}
        """

rule MP_ccsmeth_align_reads:
    input:
        fasta = get_fasta,
        bam   = rules.MP_ccsmeth_call_hifi.output.bam
    output:
        bam   = "results/{sample}/METH_PACBIO/{sample}.hifi.pbmm2.bam"
    log:
        "results/{sample}/METH_PACBIO/logs/ccsmeth_align_reads_{sample}.log"
    conda:
        "envs/ccsmeth.yaml"
    threads: config["cpus_per_task"]
    shell:
        r"""
        mkdir -p "$(dirname {log})"

        ccsmeth align_hifi \
           --hifireads {input.bam} \
           --ref {input.fasta} \
           --output {output.bam} \
           --threads {threads} &> {log}
        """

rule MP_ccsmeth_call_mods:
    input:
        fasta = get_fasta,
        bam   = rules.MP_ccsmeth_align_reads.output.bam
    output:
        "results/{sample}/METH_PACBIO/{sample}.hifi.pbmm2.call_mods.modbam.bam"
    log:
        "results/{sample}/METH_PACBIO/logs/ccsmeth_call_mods_{sample}.log"
    conda:
        "envs/ccsmeth.yaml"
    params:
        model_file = config["ccsmeth"]["call_mod"]["model_file"],
        threads_call = config["ccsmeth"]["call_mod"]["threads_call"],
        model_type = config["ccsmeth"]["call_mod"]["model_type"],
        mode = config["ccsmeth"]["call_mod"]["mode"],
        out_prefix = "results/{sample}/METH_PACBIO/{sample}.hifi.pbmm2.call_mods"
    threads: config["cpus_per_task"]
    shell:
        r"""
        mkdir -p "$(dirname {log})"

        ccsmeth call_mods --input {input.bam} \
            --ref {input.fasta} \
            --model_file {params.model_file} \
            --output {params.out_prefix} \
            --threads {threads} \
            --threads_call {params.threads_call} \
            --model_type {params.model_type} \
            --mode {params.mode} &> {log}
        """

rule MP_ccsmeth_call_freqb:
    input:
        fasta = get_fasta,
        bam   = rules.MP_ccsmeth_call_mods.output
    output:
        "results/{sample}/METH_PACBIO/{sample}.hifi.pbmm2.call_mods.modbam.freq.aggregate.all.bed"
    log:
        "results/{sample}/METH_PACBIO/logs/ccsmeth_call_freqb_{sample}.log"
    conda:
        "envs/ccsmeth.yaml"
    params:
        model_file = config["ccsmeth"]["call_freqb"]["model_file"],
        call_mode = config["ccsmeth"]["call_freqb"]["call_mode"],
        out_prefix = "results/{sample}/METH_PACBIO/{sample}.hifi.pbmm2.call_mods.modbam.freq"
    threads: config["cpus_per_task"]
    shell:
        r"""
        mkdir -p "$(dirname {log})"

        ccsmeth call_freqb \
            --input_bam {input.bam} \
            --ref {input.fasta} \
            --output {params.out_prefix} \
            --threads {threads} \
            --sort --bed \
            --call_mode {params.call_mode} \
            --aggre_model {params.model_file} &> {log}
        """

############# Centromere Scoring ##############
rule CENTROMERE_SCORING_index_fai:
    input:
        fasta = get_fasta
    output:
        fai = "results/{sample}/CENTROMERE_SCORING/{sample}.fasta.fai"
    log:
        "results/{sample}/CENTROMERE_SCORING/logs/index_fai_{sample}.log"
    conda:
        "envs/centro.yaml"
    shell:
        r"""
        mkdir -p "$(dirname {log})"

        samtools faidx {input.fasta}

        cp {input.fasta}.fai {output.fai} &> {log}
        """

rule CENTROMERE_SCORING_make_windows:
    input:
        fai = rules.CENTROMERE_SCORING_index_fai.output.fai
    output:
        bed = "results/{sample}/CENTROMERE_SCORING/{sample}.windows.{window}bp.bed"
    log:
        "results/{sample}/CENTROMERE_SCORING/logs/{sample}_window_{window}.log"
    conda:
        "envs/centro.yaml"
    params:
        window = config["window"],
        do_sort = lambda wildcard: "true" if is_nanopore(wildcard.sample) else "false"
    shell:
        r"""
        mkdir -p "$(dirname {log})"

        if [ "{params.do_sort}" = "true" ]; then
            {{ bedtools makewindows -g {input.fai} -w {params.window} \
              | sort -k1,1V -k2,2n > {output.bed}; }} &> {log}
        else
            {{ bedtools makewindows -g {input.fai} -w {params.window} > {output.bed}; }} &> {log}
        fi

        # fail fast if output is empty
        test -s {output.bed} || {{ echo "ERROR: {output.bed} is empty" >&2; exit 1; }}
        """

rule CENTROMERE_SCORING_trf2bed_sort:
    input:
        trf_bed = rules.TRF_convert_to_bed.output
    output:
        sorted_bed = "results/{sample}/CENTROMERE_SCORING/{sample}.trf.sorted.bed"
    log:
        "results/{sample}/CENTROMERE_SCORING/logs/trf2bed_sort_{sample}.log"
    conda:
        "envs/centro.yaml"
    shell:
        r"""
        mkdir -p "$(dirname {log})"

        {{ awk -F'\t' '{{
        split($1, coords, ":");
        chrom = coords[1];
        split(coords[2], range, "-");
        start = range[1] - 1;
        end = range[2];
        motif = $2;
        print chrom, start, end, motif;
        }}' OFS="\t" {input.trf_bed} | sort -k1,1V -k2,2n > {output.sorted_bed}; }} &> {log}
        """

rule CENTROMERE_SCORING_sort_TE:
    input:
        edta = rules.EDTA_bed.output
    output:
        sorted_bed = "results/{sample}/CENTROMERE_SCORING/{sample}.te.sorted.bed"
    log:
        "results/{sample}/CENTROMERE_SCORING/logs/sorted_TE_{sample}.log"
    conda:
        "envs/centro.yaml"
    shell:
        r"""
        mkdir -p "$(dirname {log})"

        {{ sort -k1,1V -k2,2n {input.edta} > {output.sorted_bed}; }} &> {log}
        """

rule CENTROMERE_SCORING_sort_methylation:
    input:
        methyl=lambda wildcard: (
            rules.MN_modbam2bed.output.bed
            if is_nanopore(wildcard.sample)
            else rules.MP_ccsmeth_call_freqb.output
        )
    output:
        bedgraph = "results/{sample}/CENTROMERE_SCORING/{sample}.methylation.sorted.bedgraph"
    log:
        "results/{sample}/CENTROMERE_SCORING/logs/sorted_methylation_{sample}.log"
    conda:
        "envs/centro.yaml"
    params:
        do_sort = lambda wildcard: "true" if is_nanopore(wildcard.sample) else "false"
    shell:
        r"""
        mkdir -p "$(dirname {log})"

        if [ "{params.do_sort}" = "true" ]; then
            {{ awk 'BEGIN{{OFS="\t"}} !/^#/ && $5!="nan" && $5!="NA" && $5!="" {{
                val=$5
                if (val < 0.0001 && val > -0.0001) next
                if (val <= 1) val = val * 100
                print $1, $2, $3, val
            }}' {input.methyl} | sort -k1,1 -k2,2n > {output.bedgraph}; }} &> {log}
        else
            {{ awk 'BEGIN{{OFS="\t"}} {{print $1,$2,$3,$11}}' {input.methyl} \
                | sort -k1,1V -k2,2n > {output.bedgraph}; }} &> {log}
        fi
        """

rule CENTROMERE_SCORING_TRF_coverage:
    input:
        bed = rules.CENTROMERE_SCORING_make_windows.output.bed,
        trf_bed = rules.CENTROMERE_SCORING_trf2bed_sort.output.sorted_bed
    output:
        tmp_bed = "results/{sample}/CENTROMERE_SCORING/{sample}.{window}.tmp.trf_counts.bed"
    log:
        "results/{sample}/CENTROMERE_SCORING/logs/TRF_coverage_{sample}.{window}.log"
    conda:
        "envs/centro.yaml"
    params:
        window = config["window"]
    shell:
        r"""
        mkdir -p "$(dirname {log})"

        {{ bedtools coverage -a {input.bed} -b {input.trf_bed} -counts > {output.tmp_bed}; }} &> {log}
        """

rule CENTROMERE_SCORING_TE_coverage:
    input:
        bed = rules.CENTROMERE_SCORING_make_windows.output.bed,
        te_bed = rules.CENTROMERE_SCORING_sort_TE.output.sorted_bed
    output:
        tmp_bed = "results/{sample}/CENTROMERE_SCORING/{sample}.{window}.tmp.te_counts.bed"
    log:
        "results/{sample}/CENTROMERE_SCORING/logs/TE_coverage_{sample}.{window}.log"
    conda:
        "envs/centro.yaml"
    params:
        window = config["window"]
    shell:
        r"""
        mkdir -p "$(dirname {log})"

        {{ bedtools coverage -a {input.bed} -b {input.te_bed} -counts > {output.tmp_bed}; }} &> {log}
        """

rule CENTROMERE_SCORING_gene_counts:
    input:
        gff3 = get_gff3
    output:
        genes_bed = "results/{sample}/CENTROMERE_SCORING/{sample}.genes.bed"
    log:
        "results/{sample}/CENTROMERE_SCORING/logs/{sample}.genes.bed.log"
    conda:
        "envs/centro.yaml"
    shell:
        r"""
        mkdir -p "$(dirname {log})"

        {{ awk '$3=="gene" {{print $1"\t"($4-1)"\t"$5"\t"$9}}' {input.gff3} > {output.genes_bed}; }} &> {log}
        """

rule CENTROMERE_SCORING_gene_counts_bedtools_coverage:
    input:
        genes_bed = rules.CENTROMERE_SCORING_gene_counts.output.genes_bed,
        window_bed = rules.CENTROMERE_SCORING_make_windows.output.bed
    output:
        genes_bed = "results/{sample}/CENTROMERE_SCORING/{sample}.{window}.tmp.gene_counts.bed"
    log:
        "results/{sample}/CENTROMERE_SCORING/logs/{sample}.{window}.gene_counts.bedtools_coverage.log"
    conda:
        "envs/centro.yaml"
    shell:
        r"""
        mkdir -p "$(dirname {log})"

        {{ bedtools coverage -a {input.window_bed} -b {input.genes_bed} -counts > {output.genes_bed}; }} &> {log}
        """

rule CENTROMERE_SCORING_hifi_coverage:
    input:
        hifi = lambda wildcard: (
            rules.MN_samtools_sort.output.bam
            if is_nanopore(wildcard.sample)
            else rules.MP_ccsmeth_align_reads.output.bam
        )
    output:
        bed = "results/{sample}/CENTROMERE_SCORING/{sample}.hifi.depth.bed"
    log:
        "results/{sample}/CENTROMERE_SCORING/logs/hifi_coverage_{sample}.log"
    conda:
        "envs/centro.yaml"
    params:
        do_sort = lambda wildcard: "true" if is_nanopore(wildcard.sample) else "false"
    shell:
        r"""
        mkdir -p "$(dirname {log})"

        if [ "{params.do_sort}" = "true" ]; then
            {{ samtools depth -a {input.hifi} | awk '{{print $1"\t"$2-1"\t"$2"\t"$3}}' \
            | sort -k1,1V -k2,2n > {output.bed}; }} &> {log}
        else
            {{ samtools depth -a {input.hifi} | awk '{{print $1"\t"$2-1"\t"$2"\t"$3}}' > {output.bed}; }} &> {log}
        fi
        """

rule CENTROMERE_SCORING_hifi_coverage_bedtools_map:
    input:
        window = rules.CENTROMERE_SCORING_make_windows.output.bed,
        bed = rules.CENTROMERE_SCORING_hifi_coverage.output.bed
    output:
        bed = "results/{sample}/CENTROMERE_SCORING/{sample}.{window}.tmp.hifi_cov_mean.bed"
    log:
        "results/{sample}/CENTROMERE_SCORING/logs/hifi_coverage_bedtools_map_{sample}_{window}.log"
    conda:
        "envs/centro.yaml"
    shell:
        r"""
        mkdir -p "$(dirname {log})"

        {{ bedtools map -a {input.window} -b {input.bed} -c 4 -o mean -null 0 > {output.bed}; }} &> {log}
        """

rule CENTROMERE_SCORING_mean_methylation_per_window:
    input:
        window = rules.CENTROMERE_SCORING_make_windows.output.bed,
        bedgraph = rules.CENTROMERE_SCORING_sort_methylation.output.bedgraph
    output:
        bed = "results/{sample}/CENTROMERE_SCORING/{sample}.{window}.tmp.meth_mean.bed"
    log:
        "results/{sample}/CENTROMERE_SCORING/logs/mean_methylation_per_window_{sample}_{window}.log"
    conda:
        "envs/centro.yaml"
    params:
        do_sort = lambda wildcard: "true" if is_nanopore(wildcard.sample) else "false"
    shell:
        r"""
        mkdir -p "$(dirname {log})"

        if [ "{params.do_sort}" = "true" ]; then
            {{ bedtools map -a {input.window} -b {input.bedgraph} -c 4 -o mean -null 0 > {output.bed}; }} &> {log}
        else
            {{ bedtools map -nonamecheck -a {input.window} -b {input.bedgraph} -c 4 -o mean -null 0 > {output.bed}; }} &> {log}
        fi
        """

rule CENTROMERE_SCORING_calculate_gc_content_per_window:
    input:
        fasta = get_fasta,
        window = rules.CENTROMERE_SCORING_make_windows.output.bed
    output:
        bed = "results/{sample}/CENTROMERE_SCORING/{sample}.{window}.tmp.gc_content.bed"
    log:
        "results/{sample}/CENTROMERE_SCORING/logs/calculate_gc_content_per_window_{sample}_{window}.log"
    conda:
        "envs/centro.yaml"
    shell:
        r"""
        mkdir -p "$(dirname {log})"

        {{ bedtools nuc -fi {input.fasta} -bed {input.window} | awk 'NR>1 {{print $1"\t"$2"\t"$3"\t"$5}}' > {output.bed}; }} &> {log}
        """

rule CENTROMERE_SCORING_combine_features:
    input:
        windows = rules.CENTROMERE_SCORING_make_windows.output.bed,
        trf     = rules.CENTROMERE_SCORING_TRF_coverage.output.tmp_bed,
        te      = rules.CENTROMERE_SCORING_TE_coverage.output.tmp_bed,
        gene    = rules.CENTROMERE_SCORING_gene_counts_bedtools_coverage.output.genes_bed,
        hifi    = rules.CENTROMERE_SCORING_hifi_coverage_bedtools_map.output.bed,
        meth    = rules.CENTROMERE_SCORING_mean_methylation_per_window.output.bed,
        gc      = rules.CENTROMERE_SCORING_calculate_gc_content_per_window.output.bed
    output:
        tsv = "results/{sample}/CENTROMERE_SCORING/{sample}.{window}.windows.features.tsv"
    log:
        "results/{sample}/CENTROMERE_SCORING/logs/combine_features_{sample}.{window}.log"
    conda:
        "envs/centro.yaml"
    shell:
        r"""
        mkdir -p "$(dirname {log})"

        {{ echo -e "chrom\tstart\tend\ttrf_cov\tte_cov\tgene_count\thifi_cov_mean\tmeth_mean\tgc_content" > {output.tsv}; }} &> {log}

        {{ paste \
            <(cut -f1-3 {input.windows}) \
            <(cut -f4 {input.trf}) \
            <(cut -f4 {input.te}) \
            <(cut -f4 {input.gene}) \
            <(cut -f4 {input.hifi}) \
            <(cut -f4 {input.meth}) \
            <(cut -f4 {input.gc}) \
            >> {output.tsv}; }} &>> {log}
        """

rule CENTROMERE_SCORING_python:
    input:
        features = rules.CENTROMERE_SCORING_combine_features.output.tsv,
        fai = rules.CENTROMERE_SCORING_index_fai.output.fai
    output:
        windows_ranked_tsv = "results/{sample}/CENTROMERE_SCORING/{sample}_{window}/centro_windows_ranked.tsv",
        best_windows_tsv   = "results/{sample}/CENTROMERE_SCORING/{sample}_{window}/centro_best_windows_marked.tsv",
        candidates_bed    = "results/{sample}/CENTROMERE_SCORING/{sample}_{window}/centro_candidates.bed",
        candidates_ranked_tsv = "results/{sample}/CENTROMERE_SCORING/{sample}_{window}/centro_candidates_ranked.tsv",
        best_candidates_bed = "results/{sample}/CENTROMERE_SCORING/{sample}_{window}/centro_best_candidates.bed"
    log:
        "results/{sample}/CENTROMERE_SCORING/{sample}_{window}/logs/centromere_scoring_final_{sample}.log"
    conda:
        "envs/centro.yaml"
    params:
        outdir = "results/{sample}/CENTROMERE_SCORING/{sample}_{window}",
        platform = lambda wc: get_platform(wc.sample),
        exclusion_bp_large = config["exclusion_bp_large"],
        exclusion_bp_min = config["exclusion_bp_min"],
        window = config["window"],
        trf = config["trf"],
        te = config["te"],
        gene = config["gene"],
        meth = config["meth"],
        cov = config["cov"],
        gc = config["gc"]
    shell:
        r"""
        mkdir -p "{params.outdir}/logs"

        python3 score_centromeres.py \
          --features "{input.features}" \
          --fai "{input.fai}" \
          --outdir "{params.outdir}" \
          --platform "{params.platform}" \
          --trf "{params.trf}" \
          --te "{params.te}" \
          --gene "{params.gene}" \
          --meth "{params.meth}" \
          --cov "{params.cov}" \
          --gc "{params.gc}" \
          --exclusion-bp-large "{params.exclusion_bp_large}" \
          --exclusion-bp-min "{params.exclusion_bp_min}" \
          --window "{wildcards.window}" \
          &> "{log}"
        """

