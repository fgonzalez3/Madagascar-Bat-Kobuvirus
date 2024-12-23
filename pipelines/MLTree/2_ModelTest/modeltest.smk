import pandas as pd 

configfile: "config/modeltest.yaml"

rule all:
    input:
        expand("results/{genera}/modeltest/{genera}_MSA.fasta.out", genera=config["genera"])

rule modeltest:
    """
    Run ModelTest on MSA
    """
    input: 
        aln=config["aln"]
    output: 
        "results/{genera}/modeltest/{genera}_MSA.fasta.out"
    conda:
        "envs/modeltest.yaml"
    params:
        genera=config["genera"]
    shell: 
        """
        modeltest-ng -d nt -i {input.aln} > {output} -t ml -p 4
        """