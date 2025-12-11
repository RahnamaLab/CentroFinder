# Created By: Sharon Colson
# Creation Date: 12/01/2025
# Last Modified: 12/09/2025

# To run this on TTU HPC:
#     spack load trf snakemake graphviz
# sample run line:
#     snakemake -np
#         OR
#     snakemake --cores 12
# Create DAG: snakemake --dag | dot -Tsvg > test.svg

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
        expand("results/{sample}/EDTA/{sample}_edta.bed", sample=SAMPLES_LIST)
#        expand("results/{sample}/METHYL/{sample}_methyl.bed", sample=SAMPLES_LIST)

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
        mkdir -p "$(dirname {output.edta_gff3})" "$(dirname {log})"

        workdir="results/{wildcards.sample}/EDTA/edta"
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
        mkdir -p $(dirname {output}) $(dirname {log})

        awk 'BEGIN{{OFS="\t"}}
             !/^#/ {{
                match($9,/classification=([^;]+)/,a);
                 class=a[1];
                if (class=="") class="NA";
                print $1, $4-1, $5, $3, class, $6, $7
             }}' {input.edta_gff} > {output} 2> {log}
        """

###### Meth Nanopore #####
#rule mn_minimap2:
#    input:
#        fasta = get_fasta
#        fastq = get_fastq
#    output:
#        "results/{sample}/meth_nanopore/{sample}.sam"
#    log:
#        "logs/{sample}/minimap2_{sample}.log"
#    shell:
#
#rule mn_samtools_sort:
#    input:
#    output:
#    log:
#    shell:
#
#rule mn_samtools_index:
#    input:
#    output:
#    log:
#    shell:
#
#rule mn_modbam2bed:
#    input:
#    output:
#    log:
#    shell:


