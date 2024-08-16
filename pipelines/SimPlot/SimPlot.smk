configfile: "config/SimPlot.yaml"

rule all:
    input:
        "results/SimPlot/NT_identity.csv", 
        "results/mafft/NT_alignment.fasta",
        "results/SimPlot/AA_identity.csv", 
        "results/mafft/AA_alignment.fasta"
        
rule SimPlot_NT:
    """
    Create NT SimPlot for Kobuvirus reps
    """
    input:
        NT_seqs=config["NT_seqs"]
    output:
        NT_aln = "results/mafft/NT_alignment.fasta",
        NT_sim = "results/SimPlot/NT_identity.csv"
    conda:
        "envs/mafft.yaml"
    shell:
        """
        mafft {input.NT_seqs} > {output.NT_aln}
        pysimplot -i {output.NT_aln} -o {output.NT_sim}
        """

rule SimPlot_AA:
    """
    Create AA SimPlot for Kobuvirus reps
    """
    input:
        AA_seqs=config["AA_seqs"]
    output:
        AA_aln = "results/mafft/AA_alignment.fasta",
        AA_sim = "results/SimPlot/AA_identity.csv"
    conda:
        "envs/mafft.yaml"
    shell:
        """
        mafft {input.AA_seqs} > {output.AA_aln}
        pysimplot -i {output.AA_aln} -o {output.AA_sim}
        """