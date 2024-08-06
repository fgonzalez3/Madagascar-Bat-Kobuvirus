configfile: "config/ref_free_BLAST.yaml"

rule all:
    input:
        "Fecal_Kobuvirus_sequences_DEDUP.fasta",
        "Fecal_Kobuvirus_sequences_DEDUP.fasta.clstr",
        "BLASTDB/Kobuvirus_NT.ndb",
        "BLASTDB/Kobuvirus_NT.nhr",
        "BLASTDB/Kobuvirus_NT.nin",
        "BLASTDB/Kobuvirus_NT.njs",
        "BLASTDB/Kobuvirus_NT.nog",
        "BLASTDB/Kobuvirus_NT.nos",
        "BLASTDB/Kobuvirus_NT.not",
        "BLASTDB/Kobuvirus_NT.nsq",
        "BLASTDB/Kobuvirus_NT.ntf",
        "BLASTDB/Kobuvirus_NT.nto",
        "BLASTDB/Kobuvirus_AA.pdb",
        "BLASTDB/Kobuvirus_AA.phr",
        "BLASTDB/Kobuvirus_AA.pin",
        "BLASTDB/Kobuvirus_AA.pjs",
        "BLASTDB/Kobuvirus_AA.pog",
        "BLASTDB/Kobuvirus_AA.pos",
        "BLASTDB/Kobuvirus_AA.pot",
        "BLASTDB/Kobuvirus_AA.psq",
        "BLASTDB/Kobuvirus_AA.ptf",
        "BLASTDB/Kobuvirus_AA.pto",
        "results/BLAST/Kobuvirus_blastn_feces_nt.txt",
        "results/BLAST/Kobuvirus_blastx_feces_aa.txt",
        "results/BLAST_parse/Kobuviruses_unique_fecal_contigs_nt.txt",
        "results/BLAST_parse/Kobuviruses_unique_fecal_sampleIDs_nt.txt",
        "results/BLAST_parse/Kobuviruses_fecal_100len5eval_nt.txt",
        "results/BLAST_parse/Kobuviruses_fecal_100len100bit_nt.txt",
        "results/BLAST_parse/Kobuviruses_unique_fecal_contigs_aa.txt",
        "results/BLAST_parse/Kobuviruses_unique_fecal_sampleIDs_aa.txt",
        "results/BLAST_parse/Kobuviruses_fecal_100len5eval_aa.txt",
        "results/BLAST_parse/Kobuviruses_fecal_100len100bit_aa.txt", 
        "results/BLAST_parse/Kobuviruses_fecal_100len5eval_nt_hiqual.txt",
        "results/BLAST_parse/Kobuviruses_fecal_100len100bit_nt_hiqual.txt",
        "results/BLAST_parse/Kobuviruses_fecal_100len5eval_aa_hiqual.txt",
        "results/BLAST_parse/Kobuviruses_fecal_100len100bit_aa_hiqual.txt" 

rule cdhit:
    """
    Deduplicate fecal contigs using CD-HIT
    """
    input:
        fec=config["fecseqs"]
    output:
        dedups="Fecal_Kobuvirus_sequences_DEDUP.fasta", 
        clstr="Fecal_Kobuvirus_sequences_DEDUP.fasta.clstr"
    conda:
        "envs/CD-HIT.yaml"
    shell:
        """
        cd-hit-est -i {input.fec} -o Fecal_Kobuvirus_sequences_DEDUP.fasta -T 8 -M 0 -c 0.95
        """

rule makeblastdb:
    """
    Make a Kobuvirus BLAST database
    """
    input:
        seqs=config["seqs"]
    output:
        "BLASTDB/Kobuvirus_NT.ndb",
        "BLASTDB/Kobuvirus_NT.nhr",
        "BLASTDB/Kobuvirus_NT.nin",
        "BLASTDB/Kobuvirus_NT.njs",
        "BLASTDB/Kobuvirus_NT.nog",
        "BLASTDB/Kobuvirus_NT.nos",
        "BLASTDB/Kobuvirus_NT.not",
        "BLASTDB/Kobuvirus_NT.nsq",
        "BLASTDB/Kobuvirus_NT.ntf",
        "BLASTDB/Kobuvirus_NT.nto",
        "BLASTDB/Kobuvirus_AA.pdb",
        "BLASTDB/Kobuvirus_AA.phr",
        "BLASTDB/Kobuvirus_AA.pin",
        "BLASTDB/Kobuvirus_AA.pjs",
        "BLASTDB/Kobuvirus_AA.pog",
        "BLASTDB/Kobuvirus_AA.pos",
        "BLASTDB/Kobuvirus_AA.pot",
        "BLASTDB/Kobuvirus_AA.psq",
        "BLASTDB/Kobuvirus_AA.ptf",
        "BLASTDB/Kobuvirus_AA.pto"
    shell:
        """
        module unload miniconda
        module load BLAST+/2.15.0-gompi-2022b
        makeblastdb -in {input.seqs} -dbtype nucl -parse_seqids -out BLASTDB/Kobuvirus_NT
        makeblastdb -in {input.seqs} -dbtype prot -parse_seqids -out BLASTDB/Kobuvirus_AA
        """

rule blast:
    """
    Run BLAST on Kobuvruses using previously curated database
    """
    input:
        dedupseqs=rules.cdhit.output.dedups
    output:
        BLASTn_raw="results/BLAST/Kobuvirus_blastn_feces_nt.txt",
        BLASTx_raw="results/BLAST/Kobuvirus_blastx_feces_aa.txt"
    shell:
        """
        module unload miniconda
        module load BLAST+/2.15.0-gompi-2022b
        blastn -word_size 10 -query {input.dedupseqs} -db BLASTDB/Kobuvirus_NT -outfmt '6 qseqid nident pident length evalue bitscore sgi sacc stitle' -max_target_seqs 10 -out {output.BLASTn_raw} -num_threads 8 -evalue 0.001
        blastx -word_size 3 -query {input.dedupseqs} -db BLASTDB/Kobuvirus_AA -outfmt '6 qseqid nident pident length evalue bitscore sgi sacc stitle' -max_target_seqs 10 -out {output.BLASTx_raw} -num_threads 8 -evalue 0.001
        """

rule blast_init_parse:
    """
    Create a few summaries of BLAST output data
    """
    input:
        BLASTn_parse=rules.blast.output.BLASTn_raw,
        BLASTx_parse=rules.blast.output.BLASTx_raw
    output:
        out_nt1="results/BLAST_parse/Kobuviruses_unique_fecal_contigs_nt.txt", 
        out_nt2="results/BLAST_parse/Kobuviruses_unique_fecal_sampleIDs_nt.txt", 
        out_nt3="results/BLAST_parse/Kobuviruses_fecal_100len5eval_nt.txt", 
        out_nt4="results/BLAST_parse/Kobuviruses_fecal_100len100bit_nt.txt",
        out_aa1="results/BLAST_parse/Kobuviruses_unique_fecal_contigs_aa.txt",
        out_aa2="results/BLAST_parse/Kobuviruses_unique_fecal_sampleIDs_aa.txt",
        out_aa3="results/BLAST_parse/Kobuviruses_fecal_100len5eval_aa.txt",
        out_aa4="results/BLAST_parse/Kobuviruses_fecal_100len100bit_aa.txt"
    shell:
        """
        cat {input.BLASTn_parse} | awk '{{print $1}}' | sort | uniq > {output.out_nt1}
        cat {input.BLASTn_parse} | awk -F\\_ '{{print $1"_"$2}}' | sort | uniq > {output.out_nt2}
        cat {input.BLASTn_parse} | awk -F\\t '($4>99 && $5<0.00001)' > {output.out_nt3}
        cat {input.BLASTn_parse} | awk -F\\t '($4>99 && $6>99)' > {output.out_nt4}
        cat {input.BLASTx_parse} | awk '{{print $1}}' | sort | uniq > {output.out_aa1}
        cat {input.BLASTx_parse} | awk -F\\_ '{{print $1"_"$2}}' | sort | uniq > {output.out_aa2}
        cat {input.BLASTx_parse} | awk -F\\t '($4>99 && $5<0.00001)' > {output.out_aa3}
        cat {input.BLASTx_parse} | awk -F\\t '($4>99 && $6>99)' > {output.out_aa4}
        """

rule blast_final_parse:
    """
    Create summaries of high quality hits from initial BLAST parse
    """
    input:
        in1=rules.blast_init_parse.output.out_nt3, 
        in2=rules.blast_init_parse.output.out_nt4, 
        in3=rules.blast_init_parse.output.out_aa3, 
        in4=rules.blast_init_parse.output.out_aa4
    output:
        out1="results/BLAST_parse/Kobuviruses_fecal_100len5eval_nt_hiqual.txt", 
        out2="results/BLAST_parse/Kobuviruses_fecal_100len100bit_nt_hiqual.txt",
        out3="results/BLAST_parse/Kobuviruses_fecal_100len5eval_aa_hiqual.txt",
        out4="results/BLAST_parse/Kobuviruses_fecal_100len100bit_aa_hiqual.txt"
    shell:
        """
        cat {input.in1} | awk '{{print $1}}' | sort | uniq > {output.out1}
        cat {input.in2} | awk '{{print $1}}' | sort | uniq > {output.out2}
        cat {input.in3} | awk '{{print $1}}' | sort | uniq > {output.out3}
        cat {input.in4} | awk '{{print $1}}' | sort | uniq > {output.out4}
        """