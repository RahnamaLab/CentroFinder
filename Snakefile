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
PACBIO_DIR = config["pacbio_dir"]

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
    return f"{base}/{wildcards.sample}.{ext}"

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
        expand("results/{sample}_trf.bed", sample=SAMPLES_LIST),
        expand("results/{sample}_edta.bed", sample=SAMPLES_LIST)

#### TRF ####
rule run_trf:
    input:
        get_fasta
    output:
        "results/{sample}.fasta" + TRF_SUFFIX
    log:
        "logs/run_trf_{sample}.log"
    run:
        os.makedirs(os.path.dirname(output[0]), exist_ok=True)
        os.makedirs(os.path.dirname(log[0]), exist_ok=True)

        # Set the work directory for TRF to "results/"
        results_dir = os.path.dirname(output[0])

        # Show TRF where the input will be from "results/"
        input_rel = os.path.relpath(input[0], results_dir)

        trf_command = f"trf {input_rel} {TRF_PARAM_STRING}"

        # This code modified from https://stackoverflow.com/questions/45613881/what-would-be-an-elegant-way-of-preventing-snakemake-from-failing-upon-shell-r-e
        try:
            proc_output = subprocess.check_output(trf_command, shell=True, cwd=results_dir, stderr=subprocess.STDOUT)

            # Log what TRF is doing
            with open(log[0], "wb") as lf:
                lf.write(proc_output)

        # an exception is raised by check_output() for non-zero exit codes (usually returned to indicate failure)
        except subprocess.CalledProcessError as exc:
            # Capture errors in the log as well
            with open(log[0], "wb") as lf:
                if exc.output:
                    lf.write(exc.output)
                lf.write(f"\nTRF exit code: {exc.returncode}\n".encode())

            # Current releases of TRF (up to 4.09.1) return an exit code of the number of TRs processed
            # Anything less than 3 should be treated as an error code. Otherwise, this should pass
            if exc.returncode is not None and exc.returncode >= 3:
                pass
            else:
                raise

        # Make certain that output file is created
        if not os.path.exists(output[0]):
            raise Exception(f"TRF did not produce expected output file: {output[0]}")

rule convert_trf_to_bed:
    input:
        "results/{sample}.fasta" + TRF_SUFFIX
    output:
        "results/{sample}_trf.bed"
    log:
        "logs/convert_trf_to_bed_{sample}.log"
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
        gff = get_gff3
    output:
        cds = "results/edta/{sample}/{sample}_cds.fasta"
    log:
        "logs/edta_cds_{sample}.log"
    shell:
        r"""
        mkdir -p $(dirname {output.cds}) $(dirname {log})

        gffread -x {output.cds} -g {input.fasta} {input.gff} &> {log}
        """

rule edta_gff2bed:
    input:
        gff = get_gff3
    output:
        bed = "results/edta/{sample}/{sample}.bed"
    log:
        "logs/edta_bed_{sample}.log"
    shell:
        r"""
        mkdir -p $(dirname {output.bed}) $(dirname {log})

        gff2bed < {input.gff} > {output.bed} 2> {log}
        """

rule edta_run:
    input:
        fasta = get_fasta,
        cds = "results/edta/{sample}/{sample}_cds.fasta",
        bed = "results/edta/{sample}/{sample}.bed"
    output:
        edta_gff3 = "results/edta/{sample}/{sample}.fasta.mod.EDTA.TEanno.gff3"
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
        "logs/edta_run_{sample}.log"
    shell:
        r"""
        mkdir -p "$(dirname {output.edta_gff3})" "$(dirname {log})"

        workdir=$(dirname {output.edta_gff3})
        sample={wildcards.sample}
        logfile="{log}"

        cp {input.fasta} "${{workdir}}/"

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
        """

rule edta_bed:
    input:
        edta_gff = "results/edta/{sample}/{sample}.fasta.mod.EDTA.TEanno.gff3"
    output:
        "results/{sample}_edta.bed"
    log:
        "logs/edta_bed_final_{sample}.log"
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
