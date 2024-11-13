rm(list=ls())

library(ggplot2)
library(ggtree)
library(ape)
library(plyr)
library(phytools)
library(phangorn)
library(ggnewscale)
library(phylobase)
library(stringr)
library(tidyr)
library(ggpubr)
library(grid)
library(gridExtra)
library(dplyr)
library(aptheme)

# tree data --------------------------------------------------------------------

# load tree
final.kobu <-  read.tree(file = paste0("data/T3.raxml.supportFBP"))
final.kobu$tip.label <- gsub("NC(\\d+)", "NC_\\1", final.kobu$tip.label)

# root it using a rabovirus as an outgroup
final.rooted.kobu <- root(final.kobu, which(final.kobu$tip.label=='NC_026314'))

# take a quick look in base R
p <- ggtree(final.rooted.kobu) + 
  geom_nodelab(aes(label=label), size=1, nudge_x=-0.01, nudge_y=0.25) +
  geom_tiplab(align= FALSE, linetype="dotted", linesize = 0.1, size = 2.5) + geom_point(colour='red')

final.rooted.kobu <- drop.tip(final.rooted.kobu, c('MN602325', 'MF352427'))

# shorten outgroup branch length 
outgroup_index <- which(final.rooted.kobu$tip.label == 'NC_026314')
branch_index <- which(final.rooted.kobu$edge[,2] == outgroup_index)
final.rooted.kobu$edge.length[branch_index] <- 8  

# check that length looks ok 
p <- ggtree(final.rooted.kobu) + 
  geom_nodelab(aes(label=label), size=1, nudge_x=-0.01, nudge_y=0.25) +
  geom_tiplab(align= FALSE, linetype="dotted", linesize = 0.1, size = 2.5) + geom_point(colour='red')

final.rooted.kobu <- drop.tip(final.rooted.kobu, c('MN602325', 'MF352427'))

p

# tree metadata ----------------------------------------------------------------

# load manual tree csv for kobuvirus 
# and fill in for NA values 
kobu.manual <- read.csv(file=paste0('data/kobuvirus_manual.csv'), 
                        header=T, stringsAsFactors = F)

kobu.manual <- read.csv(file=paste0('data/kobuvirus_manual.csv'), 
                        header=T, stringsAsFactors = F, na = "")

kobu.manual[is.na(kobu.manual)] = "NA"

# use dplyr to only include columns that will be in the final tip label

colnames(kobu.manual)

kobu.manual <- kobu.manual %>%
  dplyr::select(Accession_Number, Species, Geo_Location, host, Year)

# tree plot --------------------------------------------------------------------

# check unique hosts that will be used to color tip labels 
# and assign colors to them
unique(kobu.manual$host)

# pick order for the labels

C =c("Avian Kobuvirus","Bat Kobuvirus","Bovine Kobuvirus","Canid Kobuvirus",
     "Caprine Kobuvirus","Feline Kobuvirus","Human Kobuvirus","Ovine Kobuvirus",
     "Porcine Kobuvirus","Rabbit Kobuvirus",
     #"Rabovirus", 
     "Rodent Kobuvirus",
     "Rodent Rabovirus")

# visualize again 
p <- ggtree(final.rooted.kobu) %<+% kobu.manual + 
  geom_tippoint(aes(color=host), size=2) + 
  geom_tiplab(size=3) + 
  #geom_nodelab(size=1) +
  #scale_color_manual(values=colz, breaks=C) + 
  theme(legend.position = c(.70, .70), legend.title = element_blank())
p

# categorize node labels based on a threshold
final.rooted.kobu$node.label[final.rooted.kobu$node.label < 90] <- 0
final.rooted.kobu$node.label[final.rooted.kobu$node.label >= 90] <- 1
final.rooted.kobu$node.label <- as.character(final.rooted.kobu$node.label)
final.rooted.kobu$node.label[final.rooted.kobu$node.label == "0"] <- "<90"
final.rooted.kobu$node.label[final.rooted.kobu$node.label == "1"] <- ">=90"
final.rooted.kobu$node.label <- factor(final.rooted.kobu$node.label, levels = c("<90", ">=90"))

# create a color mapping
posfilz <- c('<90' = "white", '>=90' = "black")

# extract the node labels and ensure they are correctly aligned with the nodes in the tree
node_labels <- data.frame(node = (Ntip(final.rooted.kobu) + 1):(Ntip(final.rooted.kobu) + Nnode(final.rooted.kobu)),
                          node_label = final.rooted.kobu$node.label)

# convert the tree to a tibble and join with node labels
tree_data <- as_tibble(final.rooted.kobu) %>%
  left_join(node_labels, by = "node")

# add a "novel" category
kobu.manual$novel <- 0
kobu.manual$novel[kobu.manual$Accession_Number == "OP287812"] <- 1
kobu.manual$novel[kobu.manual$Accession_Number == "OR082796"] <- 1
kobu.manual$novel[kobu.manual$Accession_Number == "KJ641691"] <- 2
kobu.manual$novel[kobu.manual$Accession_Number == "KJ641686"] <- 2
kobu.manual$novel <- as.factor(kobu.manual$novel)

# add a bat host category
kobu.manual$bat_host <- "non-bat-host"
kobu.manual$bat_host[kobu.manual$Accession_Number %in% c("OP287812", "OR082796", "KJ641686", "KJ641691")] <- "bat-host"
kobu.manual$bat_host <- as.factor(kobu.manual$bat_host)

# tip labels on manual csv and raxml tree do not match, so fix that 
# first make a new df with the labels numbered by root tip on the raxml tree

kobu.dat <- data.frame(Accession_Number=final.rooted.kobu$tip.label, 
                       num =1:length(final.rooted.kobu$tip.label))

# then create a 3rd df and right join the original manual and the 2nd df 
# by accession_number onto this 3rd df 

kobux <- join(kobu.manual, kobu.dat, by = "Accession_Number", match = "all", 
              type = "right")


kobux$new_label <- kobux$Accession_Number

# create new tip labels

kobux$new_label[!is.na(kobux$new_label)] <- paste(kobux$Accession_Number[!is.na(kobux$Accession_Number)], " | ", 
                                                        kobux$Species[!is.na(kobux$Species)], " | ",
                                                        kobux$host[!is.na(kobux$host)], " | ",
                                                        kobux$Geo_Location[!is.na(kobux$Geo_Location)], " | ",
                                                        kobux$Year[!is.na(kobux$Year)])

kobux$new_label <- paste(kobux$Accession_Number, " | ", 
                                                  kobux$Species, " | ",
                                                  kobux$host, " | ",
                                                  kobux$Geo_Location, " | ",
                                                  kobux$Year)

kobux$Accession_Number <- paste(kobux$Accession_Number, " | ", 
                         kobux$Species, " | ",
                         kobux$host, " | ",
                         kobux$Geo_Location, " | ",
                         kobux$Year)


kobux$Accession_Number <- kobux$new_label

final.rooted.kobu$tip.label <- kobux$Accession_Number

final.rooted.kobu$tip.label

# final tree 

# plot the tree
p <- ggtree(final.rooted.kobu) %<+% tree_data %<+% kobux +
  geom_tippoint(aes(color=host, fill=host), size=2) +
  new_scale_fill() +
  geom_nodepoint(aes(fill = node_label), shape = 21, size = 1, color = "black") + 
  scale_fill_manual(values=posfilz) + 
  new_scale_fill() + 
  geom_treescale(y = 5, x = 20, fontsize = 4, offset = 1, color = "black", width = 0.5) + 
  geom_tiplab(geom="label", label.size = 0, alpha=.3, size=2.8, show.legend=F) +
  xlim(c(0,30)) +
  theme(legend.position = c(0.7, 0.5), legend.title = element_blank())

p

# now export for visualization in adobe 
ggsave(file = paste0("SuppFig1.png"),
       units="mm",  
       width=400, 
       height=250, 
       #limitsize = F,
       scale=1)#, 
