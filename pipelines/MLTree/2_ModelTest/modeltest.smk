import pandas as pd 

configfile: "config/modeltest.yaml"

rule all:
    input:
        "results/modeltest/Kobuvirus_MSA.fasta.out"

rule modeltest:
    """
    Run ModelTest on Kobuvirus alignment
    """
    input: 
        aln=config["aln"]
    output: 
        "results/modeltest/Kobuvirus_MSA.fasta.out"
    conda:
        "envs/modeltest.yaml"
    shell: 
        """
        modeltest-ng -d nt -i {input.aln} > {output} -t ml -p 8
        """