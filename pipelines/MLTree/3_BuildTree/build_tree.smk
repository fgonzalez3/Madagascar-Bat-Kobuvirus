configfile: "config/build_tree.yaml"

rule all:
    input:
        expand("{genera}.raxml.supportFBP", genera=config["genera"])
       
rule raxml:
    """
    Build ML tree
    """
    input:
        aln=config["aln"]
    output:
        "{genera}.raxml.supportFBP"
    params:
        model=config["model"], 
        genera=config["genera"]
    conda:
        "envs/raxml.yaml"
    shell:
        """
        raxml-ng-mpi --all --msa {input.aln} --model {params.model} --prefix {params.genera} --seed 12 --threads 4 --bs-metric fbp,tbe
        """