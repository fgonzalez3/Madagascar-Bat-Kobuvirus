import pandas as pd 

configfile: "config/MSA.yaml"

rule all:
    input: 
        expand("results/{genera}/MAFFT/{genera}_MSA.fasta", genera=config["genera"])

rule mafft:
    """
    MAFFT alignment 
    """
    input: 
        seqs=config["seqs"]
    output: 
        "results/{genera}/MAFFT/{genera}_MSA.fasta"
    conda:
        "envs/mafft.yaml"
    params:
        genera=config["genera"]
    shell: 
        """
        mafft {input.seqs} > {output}
        """