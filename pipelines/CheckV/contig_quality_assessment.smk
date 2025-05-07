import pandas as pd 

configfile: "config/CheckV.yaml"

rule all:
    input: 
        expand("results/{genera}/CheckV/complete_genomes.tsv", genera=config["genera"])

rule CheckV:
    """
    Run CheckV on input sequences
    """
    input: 
        seqs=config["seqs"]
    output: 
        "results/{genera}/CheckV/complete_genomes.tsv"
    params:
        genera=config["genera"],
        outdir = "results/{genera}/CheckV/"
    shell: 
        """
        module unload miniconda
        source activate /home/flg9/.conda/envs/CheckV
        checkv download_database ./
        export CHECKVDB=./
        checkv end_to_end -d checkv-db-v1.5/ {input.seqs} {params.outdir} -t 8
        """