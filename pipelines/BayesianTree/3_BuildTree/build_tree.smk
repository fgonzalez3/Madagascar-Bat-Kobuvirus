configfile: "config/build_tree.yaml"

rule all:
    input:
        "Kobuvirus_BayesianTree.log", 
        "Kobuvirus_BayesianTree.trees"
       
rule BEAST2:
    """
    Build Kobuvirus Bayesian Tree
    """
    input:
        xml=config["xml"]
    output:
        "Kobuvirus_BayesianTree.log", 
        "Kobuvirus_BayesianTree.trees"
    params:
        genera=config["genera"]
    shell:
        """
        module unload miniconda
        module load Beast/2.7.6-GCC-12.2.0-CUDA-12.0.0
        module load beagle-lib/4.0.1-GCC-12.2.0-CUDA-12.0.0
        beast -threads 9 -beagle_cpu -seed 777 {input.xml}
        """