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

kobuvirus_treefile <- "./data/LuSequencesPlusNewMadaSequences/KobuvirusLuPlusNewMadaSeqs_05MAY2025_FLG.raxml.supportFBP"
kobuvirus_metadata <- "./data/LuSequencesPlusNewMadaSequences/kobuvirus_manual_LuSequencesPlusNewMadaSeqs_05MAY2025.csv"

outpng <- "./SuppFig1.png"

# tree plot --------------------------------------------------------------------

get_tree_plot <- function(treefile, metadata) {
  
  # load in tree
  tree <- read.tree(file = paste0(treefile))

  # make sure reference genomes are lableled correctly
  tree$tip.label <- gsub("NC(\\d+)", "NC_\\1", tree$tip.label)
  
  # change tip names to match those assigned by GenBank
  # iterate through rows in tree data
  for (i in 1:length(tree$tip.label)) {
    
    # access data for the current row 
    row_data <- tree$tip.label[i]
    
    # match the tip labels to the GenBank accessions
    tree$tip.label[i] <-gsub("NODE1RR034B051", "PV833573", tree$tip.label[i])
    tree$tip.label[i] <-gsub("NODE2RR034B053", "PV833576", tree$tip.label[i])
    tree$tip.label[i] <-gsub("NODE4RR034B056", "PV833581", tree$tip.label[i])
    tree$tip.label[i] <-gsub("NODE1RR034B050", "PV833572", tree$tip.label[i])
    tree$tip.label[i] <-gsub("NODE1RR034B073", "PV833571", tree$tip.label[i])
    tree$tip.label[i] <-gsub("NODE156RR034B165", "PV833570", tree$tip.label[i])
    tree$tip.label[i] <-gsub("NODE2RR034B089", "PV833578", tree$tip.label[i])
    tree$tip.label[i] <-gsub("NODE3RR034B086", "PV833579", tree$tip.label[i])
    tree$tip.label[i] <-gsub("NODE2RR034B057", "PV833577", tree$tip.label[i])
    
  }
  
  # root tree with outgroup
  tree <- root(tree, outgroup = 'NC_026314')
  
  # modify length of outgroup for visualization
  outgroup_index <- which(tree$tip.label == 'NC_026314')
  branch_index <- which(tree$edge[,2] == outgroup_index)
  tree$edge.length[branch_index] <- 8 
  
  # normalize the node labels to a range between 0 and 1 for bootstrap support measure
  tree$node.label <- as.numeric(tree$node.label)
  
  min_val <- min(tree$node.label, na.rm = TRUE)
  max_val <- max(tree$node.label, na.rm = TRUE)
  tree$node.label <- (tree$node.label - min_val) / (max_val - min_val) * 100
  
  # load in metadata and fill for NA's
  treedata <- read.csv(file=paste0(metadata), 
                       header=T, stringsAsFactors = F)
  
  treedata <- read.csv(file=paste0(metadata), 
                       header=T, stringsAsFactors = F, na = "")
  
  treedata[is.na(treedata)] = "NA"
  
  
  # use dplyr to only include columns that will be in the final tip label
  treedata <- treedata %>%
    dplyr::select(Accession_Number, Species, Geo_Location, host, Year)
  
  print(head(treedata))
  
  
  # order the tree tip labels in the final plot 
  C =c("Avian Kobuvirus","Bat Kobuvirus","Bovine Kobuvirus","Canid Kobuvirus",
       "Caprine Kobuvirus","Feline Kobuvirus","Human Kobuvirus","Ovine Kobuvirus",
       "Porcine Kobuvirus","Rabbit Kobuvirus",
       "Rodent Kobuvirus",
       "Rodent Rabovirus")
  
  # extract the node labels and ensure they are correctly aligned with the nodes in the tree
  node_labels <- data.frame(node = (Ntip(tree) + 1):(Ntip(tree) + Nnode(tree)),
                            node_label = tree$node.label)
  
  # convert the tree to a tibble and join with node labels
  tree_data <- as_tibble(tree) %>%
    left_join(node_labels, by = "node")
  
  # tip labels on manual csv and raxml tree do not match, so fix that 
  # first make a new df with the labels numbered by root tip on the raxml tree
  kobu.dat <- data.frame(Accession_Number=tree$tip.label, 
                         num =1:length(tree$tip.label))
  
  # then create a 3rd df and right join the original manual and the 2nd df 
  # by accession_number onto this 3rd df 
  kobux <- join(treedata, kobu.dat, by = "Accession_Number", match = "all", 
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
  
  print(head(kobux))
  
  tree$tip.label <- kobux$Accession_Number
  
  tree$tip.label
  
  # visualize tree
  
  treeplot <- ggtree(tree) %<+% tree_data %<+% kobux +
    geom_tippoint(aes(color=host, fill=host), size=2) +
    new_scale_fill() +
    geom_nodepoint(aes(fill = node_label), shape = 21, size = 1, color = "black") + 
    scale_fill_gradient(low = "white", high = "black", breaks = seq(0, 100, by = 20), labels = seq(0, 100, by = 20),
                        guide = guide_colorbar(direction = "horizontal")) + 
    new_scale_fill() + 
    geom_treescale(y = 5, x = 20, fontsize = 4, offset = 1, color = "black", width = 0.5) + 
    geom_tiplab(geom="label", label.size = 0, alpha=.3, size=2.8, show.legend=F) +
    xlim(c(0,30)) +
    theme(legend.position = c(0.7, 0.5), legend.title = element_blank())
  
  return(treeplot)
  
}


# assemble plot ----------------------------------------------------------------

get_tree_plot(kobuvirus_treefile, kobuvirus_metadata)
ggsave(outpng, units="mm", width=400, height=250, scale=1)
