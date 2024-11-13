#preprep kobuvirus phylogeny

rm(list = ls())

library(plyr)
library(dplyr)

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("msa")

source("https://bioconductor.org/biocLite.R")
biocLite("msa")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.14")

#set wd
setwd('/Users/mfv2446/Desktop/kobuvirus/trees/offline_trees/kobu_tree/')
homewd=getwd()

#load the datasets and query 

kobuclass <- read.csv(file = 'kobuviruses_taxid_194960.csv', header = T, stringsAsFactors = F)
kobunonclass <- read.csv(file = 'unclassified_ kobuviruses_tax id 655314.csv', header = T, stringsAsFactors = F)

head(kobuclass)
head(kobunonclass)

#look at the unique hosts 
sort(unique(kobuclass$Host))
sort(unique(kobuclass$Species[kobuclass$Host==""])) 
kobuclass.dat = subset(kobuclass, Host=="Chiroptera" | Host =="Aselliscus stoliczkanus" | Host =="Chaerephon plicatus"|
                   Host=="Chaerephon plicatus" | Host == "Eonycteris spelaea" | Host =="Hipposideros pratti" | Host =="Hypsugo pulveratus"|
                   Host== "Hypsugo savii" | Host == "Laephotis capensis" | Host =="Pipistrellus abramus" | Host =="Pipistrellus kuhlii" |
                   Host == "Rhinolophus" | Host == "Rhinolophus affinis" | Host == "Rhinolophus blasii" | Host == "Rhinolophus cornutus" |
                   Host == "Rhinolophus ferrumequinum" | Host == "Rhinolophus hipposideros" | Host == "Rhinolophus macrotis" |
                   Host == "Rhinolophus malayanus" | Host == "Rhinolophus pusillus"  | Host == "Rhinolophus sinicus" | Host == "Rhinolophus stheno" |
                   Host == "Rousettus leschenaultii" | Host == "Rousettus sp."  | Host == "Tylonycteris pachypus" | Host == "Vespertilio sinensis")

sort(unique(kobunonclass$Host))
sort(unique(kobunonclass$Species[kobunonclass$Host=='']))
kobunonclass.dat = subset(kobunonclass, Host=="Chiroptera" | Host =="Aselliscus stoliczkanus" | Host =="Chaerephon plicatus"|
                   Host=="Chaerephon plicatus" | Host == "Eonycteris spelaea" | Host =="Hipposideros pratti" | Host =="Hypsugo pulveratus"|
                   Host== "Hypsugo savii" | Host == "Laephotis capensis" | Host =="Pipistrellus abramus" | Host =="Pipistrellus kuhlii" |
                   Host == "Rhinolophus" | Host == "Rhinolophus affinis" | Host == "Rhinolophus blasii" | Host == "Rhinolophus cornutus" |
                   Host == "Rhinolophus ferrumequinum" | Host == "Rhinolophus hipposideros" | Host == "Rhinolophus macrotis" |
                   Host == "Rhinolophus malayanus" | Host == "Rhinolophus pusillus"  | Host == "Rhinolophus sinicus" | Host == "Rhinolophus stheno" |
                   Host == "Rousettus leschenaultii" | Host == "Rousettus sp."  | Host == "Tylonycteris pachypus" | Host == "Vespertilio sinensis")

#no bat hosts found, combine the two 
kobu.all <- rbind(kobuclass,kobunonclass) #120 genomes
kobu.all.ref = subset(kobu.all, Sequence_Type == 'RefSeq')
sort(unique(kobu.all$Host))
kobu.all <- kobu.all[!duplicated(kobu.all),] #still 120 genomes

#get the text to download from NCBI
sort(unique(kobu.all$Accession))
kobu.all <- kobu.all[!duplicated(kobu.all$Accession),] #drops down to 101 genomes

accession_num <- paste(c(kobu.all$Accession), collapse = ',')

#now put this into your web browser to download 
text.for.NCBI <- paste0('http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&rettype=fasta&retmode=text&id=', accession_num)

#once downloaded, send to MAFFT for alignment
rm(list = ls())

#then, after alignment is ready, prepare the names for RAxML (no space, semicolon, colon, parenthesis, dash, slash, comma, quote allowed in name (should all just be underscore))
library(seqinr)
#library(msa)


alignment_kobu <- read.alignment(file = '/Users/mfv2446/Desktop/2C_peptide_translation_min.fasta', format='fasta', forceToLower = F)

tmp <- as.list(alignment_kobu$nam)

change.spacing <- function(df){
  df_new <- sapply(strsplit(df, '-'), function(x) x[[1]])
  return(df_new)
}


names_new.1 = c(unlist(lapply(tmp, change.spacing)))
names_new.1 <- lapply(names_new.1, gsub, pattern = ".1", replacement = "", fixed = TRUE)
names_new.1 <- lapply(names_new.1, gsub, pattern = ".2", replacement = "", fixed = TRUE)
names_new.1 <- lapply(names_new.1, gsub, pattern = ".3", replacement = "", fixed = TRUE)
names_new.1 <- lapply(names_new.1, gsub, pattern = " - 2C peptide", replacement = "", fixed = TRUE)
names_new.1 <- lapply(names_new.1, gsub, pattern = "_", replacement = "", fixed = TRUE)



class(alignment_kobu$seq)
write.fasta(sequences = as.list(alignment_kobu$seq), names = names_new.1, file.out = '/Users/mfv2446/Desktop/kobuvirus/trees/offline_trees/kobu_tree/picorna_align_2.fasta', as.string = T, open = 'w')
write.fasta(sequences = as.list(alignment_kobu$seq), names = names_new.1, file.out = '/Users/mfv2446/Desktop/picorna_align.fasta', as.string = T, open = 'w')
write.fasta(sequences = as.list(alignment_kobu$seq), names = names_new.1, file.out = '/Users/mfv2446/Desktop/extra_kobu.fasta', as.string = T, open = 'w')

