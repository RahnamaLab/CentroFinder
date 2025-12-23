# Created By: Sharon Colson
# Creation Date: 12/01/2025
# Last Modified: 12/09/2025

# To run this on TTU HPC:
#     spack load trf snakemake graphviz
# sample run line:
#     snakemake -np
#         OR
#     snakemake --cores 12
#         OR
#     snakemake --use-conda --cores 12 results/Fo4287v4/METH_PACBIO/Fo4287v4.hifi.pbmm2.bam
# Create DAG: snakemake --dag results/Guy11_chr1/CENTROMERE_SCORING/Guy11_chr1.1000.te.sorted.bed | dot -Tsvg > centromere_pipeline_Guy11_chr1.svg

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

SAMPLES_DICT = config["samples"]
SAMPLES_LIST = list(SAMPLES_DICT.keys())

WINDOW = config["window"]

def get_base_dir(sample):
    platform = SAMPLES_DICT[sample]["platform"].lower()
    if platform == "nanopore":
        return NANOPORE_DIR
    elif platform == "pacbio":
        return PACBIO_DIR
    else:
        raise ValueError(f"Unknown platform: '{platform}' for sample: '{sample}'. "
                          "Please check the platform specified in your config.yaml for this sample."
                         )

def get_path_with_ext(wildcards, ext):
    base = get_base_dir(wildcards.sample)
    return f"{base}/{wildcards.sample}/{wildcards.sample}.{ext}"

def get_fasta(wildcards):
    return get_path_with_ext(wildcards, "fasta")

def get_gff3(wildcards):
    return get_path_with_ext(wildcards, "gff3")

def get_fastq(wildcards):
    return get_path_with_ext(wildcards, "fastq")

def get_bam(wildcards):
    return get_path_with_ext(wildcards, "subreads.bam")

def is_nanopore(sample):
    return SAMPLES_DICT[sample]["platform"].lower() == "nanopore"

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
        expand("results/{sample}/TRF/{sample}_trf.bed", sample=SAMPLES_LIST),
        expand("results/{sample}/EDTA/{sample}_edta.bed", sample=SAMPLES_LIST),
        expand("results/{sample}/METH_NANOPORE/{sample}_methyl.bed", sample=SAMPLES_LIST),
        expand("results/{sample}/METH_PACBIO/{sample}.hifi.pbmm2.call_mods.modbam.freq.aggregate.all.bed", sample=SAMPLES_LIST),
        # Centromere Scoring Notes
#        expand("results/{sample}/CENTROMERE_SCORING/{sample}.fasta.fai", sample=SAMPLES_LIST)
        expand("results/{sample}/CENTROMERE_SCORING/windows.{sample}.{window}bp.bed", sample=SAMPLES_LIST, window=WINDOW),
        expand("results/{sample}/CENTROMERE_SCORING/{sample}.trf.sorted.bed", sample=SAMPLES_LIST),
        expand("results/{sample}/CENTROMERE_SCORING/{sample}.te.sorted.bed", sample=SAMPLES_LIST),
        expand("results/{sample}/CENTROMERE_SCORING/{sample}.methylation.sorted.bedgraph", sample=SAMPLES_LIST)

#### TRF ####
rule run_trf:
    input:
        get_fasta
    output:
        "results/{sample}/TRF/{sample}.fasta" + TRF_SUFFIX
    log:
        "results/{sample}/TRF/logs/run_trf_{sample}.log"
    run:
        # Ensure directories exist
        os.makedirs(os.path.dirname(output[0]), exist_ok=True)
        os.makedirs(os.path.dirname(log[0]), exist_ok=True)

        # Work directory for TRF is the results dir for this file
        results_dir = os.path.dirname(output[0])

        # Make input path relative to results_dir
        input_rel = os.path.relpath(input[0], results_dir)

        trf_command = f"trf {input_rel} {TRF_PARAM_STRING}"

        # This code modified from https://stackoverflow.com/questions/45613881/what-would-be-an-elegant-way-of-preventing-snakemake-from-failing-upon-shell-r-e
        try:
            # Run TRF; this will raise CalledProcessError on non-zero exit codes
            proc_output = subprocess.check_output(
                trf_command,
                shell=True,
                cwd=results_dir,
                stderr=subprocess.STDOUT,
            )

            # Log normal output
            with open(log[0], "wb") as lf:
                lf.write(proc_output)
                lf.write(b"\nTRF exit code: 0\n")

        except subprocess.CalledProcessError as exc:
            # Log TRF output and exit code even on non-zero exit
            with open(log[0], "wb") as lf:
                if exc.output:
                    lf.write(exc.output)
                lf.write(f"\nTRF exit code: {exc.returncode}\n".encode())

            # If TRF did not produce the expected .dat file, THEN treat as failure
            if not os.path.exists(output[0]):
                raise

        # Final safety check: make sure the .dat file exists
        if not os.path.exists(output[0]):
            raise Exception(
                f"TRF failed to produce expected output file: {output[0]}"
            )

rule convert_trf_to_bed:
    input:
        "results/{sample}/TRF/{sample}.fasta" + TRF_SUFFIX
    output:
        "results/{sample}/TRF/{sample}_trf.bed"
    log:
        "results/{sample}/TRF/logs/convert_trf_to_bed_{sample}.log"
    shell:
        r"""
        mkdir -p $(dirname {output}) $(dirname {log})
        python3 trf2bed.py \
            --dat {input} \
            --bed {output} \
            --tool repeatseq &> {log}
        """

##### EDTA #####
rule edta_cds:
    input:
        fasta = get_fasta,
        gff   = get_gff3
    output:
        cds = "results/{sample}/EDTA/{sample}_cds.fasta"
    log:
        "results/{sample}/EDTA/logs/edta_cds_{sample}.log"
    shell:
        r"""
        mkdir -p $(dirname {output.cds}) $(dirname {log})

        gffread -x {output.cds} -g {input.fasta} {input.gff} &> {log}
        """

rule edta_gff2bed:
    input:
        gff = get_gff3
    output:
        bed = "results/{sample}/EDTA/{sample}.bed"
    log:
        "results/{sample}/EDTA/logs/edta_bed_{sample}.log"
    shell:
        r"""
        mkdir -p $(dirname {output.bed}) $(dirname {log})

        gff2bed < {input.gff} > {output.bed} 2> {log}
        """

rule edta_run:
    input:
        fasta = get_fasta,
        cds = "results/{sample}/EDTA/{sample}_cds.fasta",
        bed = "results/{sample}/EDTA/{sample}.bed"
    output:
        edta_gff3 = "results/{sample}/EDTA/{sample}.fasta.mod.EDTA.TEanno.gff3"
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
    log:
        "results/{sample}/EDTA/logs/edta_run_{sample}.log"
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

rule edta_bed:
    input:
        edta_gff = "results/{sample}/EDTA/{sample}.fasta.mod.EDTA.TEanno.gff3"
    output:
        "results/{sample}/EDTA/{sample}_edta.bed"
    log:
        "results/{sample}/EDTA/logs/edta_bed_final_{sample}.log"
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
rule mn_minimap2:
    input:
        fasta = get_fasta,
        fastq = get_fastq
    output:
        sam = "results/{sample}/METH_NANOPORE/{sample}.sam"
    log:
        "results/{sample}/METH_NANOPORE/logs/minimap2_{sample}.log"
    threads: config["cpus_per_task"]
    shell:
        r"""
        mkdir -p "$(dirname {log})"

        minimap2 -t {threads} -ax map-ont -Y {input.fasta} {input.fastq} > {output.sam} 2> {log}
        """

rule mn_samtools_sort:
    input:
        sam = "results/{sample}/METH_NANOPORE/{sample}.sam"
    output:
        bam = "results/{sample}/METH_NANOPORE/{sample}_sorted.bam"
    log:
        "results/{sample}/METH_NANOPORE/logs/mn_samtools_sort_{sample}.log"
    shell:
        r"""
        mkdir -p "$(dirname {log})"

        samtools sort -@ {threads} -o {output.bam} {input.sam} 2> {log}
        """


rule mn_samtools_index:
    input:
        bam = "results/{sample}/METH_NANOPORE/{sample}_sorted.bam"
    output:
        bai = "results/{sample}/METH_NANOPORE/{sample}_sorted.bam.bai"
    log:
        "results/{sample}/METH_NANOPORE/logs/mn_samtools_index_{sample}.log"
    shell:
        r"""
        mkdir -p "$(dirname {log})"

        samtools index {input.bam} 2> {log}
        """

rule mn_modbam2bed:
    input:
        fasta = get_fasta,
        bam   = "results/{sample}/METH_NANOPORE/{sample}_sorted.bam",
        bai   = "results/{sample}/METH_NANOPORE/{sample}_sorted.bam.bai"
    output:
        bed = "results/{sample}/METH_NANOPORE/{sample}_methyl.bed"
    log:
        "results/{sample}/METH_NANOPORE/logs/mn_modbam2bed_{sample}.log"
    params:
        modbam2bed = config["modbam2bed"]
    shell:
        r"""
        mkdir -p "$(dirname {log})"

        {params.modbam2bed} {input.fasta} {input.bam} > {output.bed} 2> {log}
        """


#### Meth Pacbio ####
rule ccsmeth_call_hifi:
    input:
        bam   = get_bam
    output:
        bam = "results/{sample}/METH_PACBIO/{sample}.hifi.bam"
    log:
        "results/{sample}/METH_PACBIO/logs/ccsmeth_call_hifi_{sample}.log"
    threads: config["cpus_per_task"]
    conda:
        "envs/ccsmeth.yaml"
    shell:
        r"""
        mkdir -p "$(dirname {log})"

        ccsmeth call_hifi \
           --subreads {input.bam} \
           --threads {threads} \
           --output {output.bam} &> {log}
        """

rule ccsmeth_align_reads:
    input:
        fasta = get_fasta,
        bam   = "results/{sample}/METH_PACBIO/{sample}.hifi.bam"
    output:
        bam   = "results/{sample}/METH_PACBIO/{sample}.hifi.pbmm2.bam"
    log:
        "results/{sample}/METH_PACBIO/logs/ccsmeth_align_reads_{sample}.log"
    threads: config["cpus_per_task"]
    conda:
        "envs/ccsmeth.yaml"
    shell:
        r"""
        mkdir -p "$(dirname {log})"

        ccsmeth align_hifi \
           --hifireads {input.bam} \
           --ref {input.fasta} \
           --output {output.bam} \
           --threads {threads} &> {log}
        """

rule ccsmeth_call_mods:
    input:
        fasta = get_fasta,
        bam   = "results/{sample}/METH_PACBIO/{sample}.hifi.pbmm2.bam"

    output:
        "results/{sample}/METH_PACBIO/{sample}.hifi.pbmm2.call_mods.modbam.bam"
    log:
        "results/{sample}/METH_PACBIO/logs/ccsmeth_call_mods_{sample}.log"
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

rule ccsmeth_call_freqb:
    input:
        fasta = get_fasta,
        bam   = "results/{sample}/METH_PACBIO/{sample}.hifi.pbmm2.call_mods.modbam.bam"
    output:
        "results/{sample}/METH_PACBIO/{sample}.hifi.pbmm2.call_mods.modbam.freq.aggregate.all.bed"
    log:
        "results/{sample}/METH_PACBIO/logs/ccsmeth_call_freqb_{sample}.log"
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
rule centromere_scoring_index_fai:
    input:
        fasta = get_fasta
    output:
        fai = "results/{sample}/CENTROMERE_SCORING/{sample}.fasta.fai"
    log:
        "results/{sample}/CENTROMERE_SCORING/logs/index_fai_{sample}.log"
    shell:
        r"""
        mkdir -p "$(dirname {log})"

        samtools faidx {input.fasta}

        cp {input.fasta}.fai {output.fai} &> {log}
        """

rule centromere_scoring_make_windows:
    input:
        fai = rules.centromere_scoring_index_fai.output.fai
    output:
        bed = "results/{sample}/CENTROMERE_SCORING/windows.{sample}.{window}bp.bed"
    log:
        "results/{sample}/CENTROMERE_SCORING/logs/{sample}_window_{window}.log"
    params:
        window = config["window"],
        do_sort = lambda wildcard: "true" if is_nanopore(wildcard.sample) else "false"
    shell:
        r"""
        mkdir -p "$(dirname {log})"

        if [ "{params.do_sort}" = "true" ]; then
            bedtools makewindows -g {input.fai} -w {params.window} \
              | sort -k1,1V -k2,2n > {output.bed} &> {log}
        else
            bedtools makewindows -g {input.fai} -w {params.window} > {output.bed} &> {log}
        fi
        """

rule centromere_scoring_trf2bed_sort:
    input:
#        bed = rules.centromere_scoring_make_windows.output.bed,
        trf_bed = rules.convert_trf_to_bed.output
    output:
        sorted_bed = "results/{sample}/CENTROMERE_SCORING/{sample}.trf.sorted.bed"
    log:
        "results/{sample}/CENTROMERE_SCORING/logs/trf2bed_sort_{sample}.log"
    shell:
        r"""
        mkdir -p "$(dirname {log})"

        awk -F'\t' '{{
        split($1, coords, ":");
        chrom = coords[1];
        split(coords[2], range, "-");
        start = range[1] - 1;
        end = range[2];
        motif = $2;
        print chrom, start, end, motif;
        }}' OFS="\t" {input.trf_bed} | sort -k1,1V -k2,2n > {output.sorted_bed} &> {log}
        """

rule centromere_scoring_sort_TE:
    input:
        edta = rules.edta_bed.output
    output:
        sorted_bed = "results/{sample}/CENTROMERE_SCORING/{sample}.te.sorted.bed"
    log:
        "results/{sample}/CENTROMERE_SCORING/logs/sorted_TE_{sample}.log"
    shell:
        r"""
        mkdir -p "$(dirname {log})"

        sort -k1,1V -k2,2n {input.edta} > {output.sorted_bed} &> {log}
        """

rule centromere_scoring_sort_methylation:
    input:
 #       fai = rules.centromere_scoring_index_fai.output.fai
        methyl=lambda wildcard: (
            rules.mn_modbam2bed.output.bed
            if is_nanopore(wildcard.sample)
            else rules.ccsmeth_call_freqb.output
        )
#        methyl = rules.mn_modbam2bed.output.bed
    output:
        bedgraph = "results/{sample}/CENTROMERE_SCORING/{sample}.methylation.sorted.bedgraph"
    log:
        "results/{sample}/CENTROMERE_SCORING/logs/sorted_methylation_{sample}.log"
    params:
#        window = config["window"],
        do_sort = lambda wildcard: "true" if is_nanopore(wildcard.sample) else "false"
    shell:
        r"""
        mkdir -p "$(dirname {log})"

        if [ "{params.do_sort}" = "true" ]; then
            awk 'BEGIN{{OFS="\t"}} !/^#/ && $5!="nan" && $5!="NA" && $5!="" {{
                val=$5
                if (val < 0.0001 && val > -0.0001) next
                if (val <= 1) val = val * 100
                print $1, $2, $3, val
            }}' {input.methyl} | sort -k1,1 -k2,2n > {output.bedgraph} &> {log}
        else
            awk 'BEGIN{{OFS="\t"}} {{print $1,$2,$3,$11}}' {input.methyl} \
                | sort -k1,1V -k2,2n > {output.bedgraph} &> {log}
        fi
        """
