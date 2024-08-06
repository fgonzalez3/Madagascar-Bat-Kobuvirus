configfile: "config/build_tree.yaml"

rule all:
    input:
        "results/RAxML/Kobuviru_nt.newick"
       
rule raxml:
    """
    Build Kobuvirus ML tree
    """
    input:
        aln=config["aln"]
    output:
        "results/RAxML/Kobuvirus.newick"
    params:
        model=config["model"]
    envs:
        "envs/raxml.yaml"
    shell:
        """
        raxml-ng-mpi --all --msa {input.aln} --model {params.model} --prefix T3 --seed 12 --threads 8 --bs-metric fbp, tbe
        """