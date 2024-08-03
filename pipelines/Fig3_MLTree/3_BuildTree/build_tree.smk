configfile: "config/SP_parnas.yaml"

rule all:
    input:
        "results/RAxML/Kobuvirus.newick"
       
rule raxml:
    """
    Build Kobuvirus ML tree
    """
    input:
        aln=config["aln"]
    output:
        "results/RAxML/Kobuvirus.newick"
    envs:
        "envs/raxml.yaml"
    shell:
        """
        raxml-ng-mpi --all --msa {input.aln} --model --prefix T3 --seed 12 --threads 8 --bs-metric fbp, tbe
        """