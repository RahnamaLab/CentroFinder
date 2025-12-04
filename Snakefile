# Created By: Sharon Colson
# Creation Date: 12/01/2025
# Last Modified:

# To run this on TTU HPC:
#     spack load trf snakemake graphviz
# sample run line:
#     snakemake -np --config sample=Guy11_Final_S4
#         OR
#     snakemake --cores 12 --config sample=Guy11_Final_S4
# Create DAG: snakemake --dag --config sample=Guy11_Final_S4 | dot -Tsvg > test.svg

#from pathlib import Path

configfile: "config.yaml"

rule all:
    input:
       expand("results/{sample}_trf.bed", sample=config["samples"])

rule trf1:
    input:
        "data/nanopore/{sample}.fasta"
    output:
        "results/{sample}.fasta.2.7.7.80.10.50.2000.dat"  ##### FIXME: hardcoded variables
    log:
        "logs/trf1_{sample}.log"
    run:
        os.makedirs(os.path.dirname(output[0]), exist_ok=True)
        os.makedirs(os.path.dirname(log[0]), exist_ok=True)

        # Set the work directory for TRF to "results/" 
        results_dir = os.path.dirname(output[0])

        # Show TRF where the input will be from "results/"
        input_rel = os.path.relpath(input[0], results_dir)

        # This code modified from https://stackoverflow.com/questions/45613881/what-would-be-an-elegant-way-of-preventing-snakemake-from-failing-upon-shell-r-e
        try:
            proc_output = subprocess.check_output(f"trf {input_rel} 2 7 7 80 10 50 2000 -h", shell=True, cwd=results_dir, stderr=subprocess.STDOUT)

            # Log what TRF is doing
            with open(log[0], "wb") as lf:
                lf.write(proc_output)

        # an exception is raised by check_output() for non-zero exit codes (usually returned to indicate failure)
        except subprocess.CalledProcessError as exc: 
            # Capture errors in the log as well
            with open(log[0], "wb") as lf: 
                if exc.output:
                    lf.write(exc.output)

            if exc.returncode == 7: #### FIXME: this is taken from the above hardcoded numbers
                # this exit code is OK
                pass
            else:
                # for all others, re-raise the exception
                raise

rule trf2:
    input:
        "results/{sample}.fasta.2.7.7.80.10.50.2000.dat"     ##### FIXME: hardcoded variables
    output:
        "results/{sample}_trf.bed"
    log:
        "logs/trf2_{sample}.log"
    shell:
        r"""
        python3 trf2bed.py \
            --dat {input} \
            --bed {output} \
            --tool repeatseq &> {log}
        """
