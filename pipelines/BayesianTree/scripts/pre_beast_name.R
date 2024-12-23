rm(list=ls())

library(plyr)
library(dplyr)
library(seqinr)

#load the dataset and query
dat <- read.csv(file = "BEAST_Metadata_LuSequences18DEC2024.csv", header = TRUE, stringsAsFactors = FALSE)
dat$collection_date <- as.Date(dat$collection_date, format = "%m/%d/%y")
head(dat)

head(dat)

#now, get the fasta file 
fasta.dat <- read.fasta(file= paste0("BEASTMSALuSequences18DEC2024Extraction.fasta"), forceDNAtolower = F, as.string = T)

names(fasta.dat)

fasta.meta <-  cbind.data.frame(tip_label = names(fasta.dat))

#add beast name to the main dataset
dat$collection_date <- as.Date(dat$collection_date)
#dat$collection_date <- as.character(dat$collection_date)
#dat$collection_date[is.na(dat$collection_date)] <- paste0(dat$collection_year[is.na(dat$collection_date)], "07-31")
#dat$collection_date <- as.Date(dat$collection_date)

dat$beast_name <- paste0(dat$accession_num, "_", dat$collection_date)

dat.merge <- dplyr::select(dat, tip_label, beast_name)


setdiff(dat.merge$tip_label, fasta.meta$tip_label)
setdiff( fasta.meta$tip_label,dat.merge$tip_label)

#fasta.meta$tip_label[fasta.meta$tip_label=="NODE_4"] <- "NODE_4_eidolon_dupreanum_2018-07-27"

#and add to data
fasta.meta <- merge(fasta.meta, dat.merge, all.x = T, by="tip_label", sort = F)
head(fasta.meta)
fasta.meta$beast_name
#and write over:
write.fasta(fasta.dat, names=fasta.meta$beast_name, file =paste0("BEASTMSALuSequencesRenamed18DEC2024.fasta"))


