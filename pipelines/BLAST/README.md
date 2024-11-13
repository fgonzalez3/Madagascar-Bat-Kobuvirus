# BLAST Directory Description

Note - To get this workflow running, all you'll need to change is the config file being specified within the snakemake workflow

Two separate BLAST analyses were run, an NCBI-guided analysis using previously sequenced kobuviruses and another guided analysis using our novel kobuvirus (OP287812)

For the NCBI-guided analysis, I downloaded all nucleotide and protein full genome reference sequences under the kobuvirus ID within NCBI virus (last accessed: October 2021). These files are specified within the ref_free_BLAST.yaml config file in the config subdirectory

The OP287812-guided analysis was run following identifying the full-length kobuvirus genome from the NCBI-guided analysis to trim out low-quality hits. This file is specified within the ref_based_BLAST.yaml config file in the config sub directory

Hits were further parsed using 100 aln & 100 bit scores or 100 aln and 5 eval scores for both nucleotide and amino acid searches

Two contigs were common hits across the above four parses, our full genome which we would submit to NCBI as OP287812 and our partial sequence which we would submit to NCBI as OR082796. These were deemed true positives
