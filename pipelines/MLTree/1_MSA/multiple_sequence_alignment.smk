import pandas as pd 

configfile: "config/MSA.yaml"

rule all:
    input: 
        "results/MAFFT/Kobuvirus_MSA.fasta"

rule mafft:
    """
    MAFFT alignment for whole genome NT Kobuvirus sequences 
    """
    input: 
        seqs=config["seqs"]
    output: 
        "results/MAFFT/Kobuvirus_MSA.fasta"
    conda:
        "envs/mafft.yaml"
    shell: 
        """
        mafft {input.seqs} > {output}
        """