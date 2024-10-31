rm(list=ls())

# make fig 1b - nt full genome kobuvirus

# load in necessary packages 

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

# set working directory 

homewd= '/Users/flg9/Desktop/Developer/brook_lab/sort this mess out/kobu_tree/'

setwd(paste0(homewd))

# load fig 1b tree 
final.kobu <-  read.tree(file = paste0("07JUL2023/T3.raxml.supportFBP"))
final.kobu$tip.label <- gsub("NC(\\d+)", "NC_\\1", final.kobu$tip.label)

# root it using a rabovirus as an outgroup
final.rooted.kobu <- root(final.kobu, which(final.kobu$tip.label=='NC_026314'))

# take a quick look in base R
ggtree(final.rooted.kobu) + 
  geom_nodelab(aes(label=label), size=1, nudge_x=-0.01, nudge_y=0.25) +
  geom_tiplab(align= FALSE, linetype="dotted", linesize = 0.1, size = 2.5) + geom_point(colour='red')

final.rooted.kobu <- drop.tip(final.rooted.kobu, c('MN602325', 'MF352427'))

# load manual tree csv for kobuvirus 
# and fill in for NA values 

kobu.manual <- read.csv(file=paste0(homewd, 'kobuvirus_manual2.csv'), header=T, stringsAsFactors = F)

kobu.manual <- read.csv(file=paste0('kobuvirus_manual2.csv'), 
                        header=T, stringsAsFactors = F, na = "")

kobu.manual[is.na(kobu.manual)] = "NA"

# use dplyr to only include columns that will be in the final tip label

colnames(kobu.manual)

kobu.manual <- kobu.manual %>%
  dplyr::select(Accession_Number, Species, Geo_Location, host, Year)


# check unique hosts that will be used to color tip labels 
# and assign colors to them

unique(kobu.manual$host)

colz = c("Human Kobuvirus" = "red", "Bovine Kobuvirus" = "tomato", 
         "Porcine Kobuvirus" = "khaki", "Ovine Kobuvirus" = "yellow", 
         "Canid Kobuvirus" = "blue", "Rodent Kobuvirus" = "plum", 
         "Sewage Kobuvirus" = "pink", "Caprine Kobuvirus" = "darkolivegreen1",
         "Feline Kobuvirus" = "purple", "Avian Kobuvirus" = "darksalmon",  
         #"Rabovirus" = "brown", 
         "Rabbit Kobuvirus" = "darkgreen", 
         "Bat Kobuvirus" = "cyan")


# pick order for the labels
C =c("Human Kobuvirus","Bovine Kobuvirus","Porcine Kobuvirus","Ovine Kobuvirus",
     "Canid Kobuvirus","Rodent Kobuvirus","Sewage Kobuvirus","Caprine Kobuvirus",
     "Feline Kobuvirus","Avian Kobuvirus",
     #"Rabovirus", 
     "Rabbit Kobuvirus",
     "Bat Kobuvirus", "Bat Picornavirus")


# visualize again 
p <- ggtree(final.rooted.kobu) %<+% kobu.manual + 
  geom_tippoint(aes(color=host), size=2) + 
  geom_tiplab(size=3) + 
  geom_nodelab(size=1) +
  scale_color_manual(values=colz, breaks=C) + 
  theme(legend.position = c(.70, .70), legend.title = element_blank())
p

# Convert the phylo object to a ggtree object
ggtree_obj <- ggtree(final.rooted.kobu)

# Define the tips for which you want to find the MRCA
Aiv_to_collapse <- c("NC_027918",
                     "KT325852")

# Get the MRCA node
mrca_node <- getMRCA(final.rooted.kobu, Aiv_to_collapse)

# Add a "novel" category
kobu.manual$novel <- 0
kobu.manual$novel[kobu.manual$Accession_Number == "OP287812"] <- 1
kobu.manual$novel[kobu.manual$Accession_Number == "OR082796"] <- 1
kobu.manual$novel[kobu.manual$Accession_Number == "KJ641691"] <- 2
kobu.manual$novel[kobu.manual$Accession_Number == "KJ641686"] <- 2
kobu.manual$novel <- as.factor(kobu.manual$novel)

# Add a bat host category
kobu.manual$bat_host <- "non-bat-host"
kobu.manual$bat_host[kobu.manual$Accession_Number %in% c("OP287812", "OR082796", "KJ641686", "KJ641691")] <- "bat-host"
kobu.manual$bat_host <- as.factor(kobu.manual$bat_host)

# Assign shapes and colors for bat hosts
shapez <- c("bat-host" = 24, "non-bat-host" = 21)
colz2 <- c('2' = "pink", '1' = "yellow", '0' = "white")

# tip labels on manual csv and raxml tree do not match, so fix that 
# first make a new df with the labels numbered by root tip on the raxml tree

kobu.dat <- data.frame(Accession_Number=final.rooted.kobu$tip.label, 
                       num =1:length(final.rooted.kobu$tip.label))

# then create a 3rd df and right join the original manual and the 2nd df 
# by accession_number onto this 3rd df 

kobux <- join(kobu.manual, kobu.dat, by = "Accession_Number", match = "all", 
              type = "right")


kobux$new_label <- kobux$Accession_Number


# match tip labels from tree with csv file 
#kobux <- kobu.manual[match(kobu.dat$new_label, kobu.manual$new_label),]
#kobux <- kobux[!is.na(kobux$Accession_Number), ]


#kobux$Accession <- kobux$new_label
#final.rooted.kobu$tip.label <- kobux$Accession


# original 
# create new tip labels

#kobux$Accession_Number<- str_replace_all(kobux$Accession_Number, "NA", "")


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


#kobux$Accession_Number <- gsub("NA", "", kobux$Accession_Number)
#kobux$Accession_Number <- gsub("| |", "", kobux$Accession_Number)


final.rooted.kobu$tip.label <- kobux$Accession_Number

final.rooted.kobu$tip.label

# final tree 

# Plot the tree
p <- ggtree(final.rooted.kobu) %<+% kobux +
  geom_tippoint(aes(color=host, fill=host, shape=bat_host), size=2) +
  scale_color_manual(values=colz, breaks = C) +
  scale_fill_manual(values=colz) + 
  scale_shape_manual(values=shapez) + 
  new_scale_fill() + geom_treescale(y = -2) +
  geom_tiplab(geom="label", label.size = 0, alpha=.3, size=2.5, show.legend=F) +
  #scale_fill_manual(values=colz2) + 
  ggtitle("Kobuvirus Tree") + xlim(c(-15,32)) +
  theme(legend.position = c(0.7, 0.5), legend.title = element_blank())

# Print the plot
print(p)

final.rooted.kobu$node.label

# collapse clades
p <- collapse(p, node=152, mode="max", fill="khaki", color="black", size=0.1)
p <- collapse(p, node=141, mode="max", fill="tomato", color="black", size=0.1)
p <- collapse(p, node=114, mode="max", fill="blue", color="black", size=0.1)
p <- collapse(p, node=125, mode="max", fill="purple", color="black", size=0.1)
p <- collapse(p, node=134, mode="max", fill="red", color="black", size=0.1)
p <- collapse(p, node=199, mode="max", fill="darkolivegreen1", color="black", size=0.1)
p <- collapse(p, node=132, mode="max", fill="plum", color="black", size=0.1)

# Add the clade label
p <- p + geom_cladelab(node=152, label="Porcine Kobuvirus", bar=1, offset=.43, fontsize=3, color="black") +
  geom_cladelab(node=141, label="Bovine Kobuvirus", bar=1, offset=.28, fontsize=3, color="black") + 
  geom_cladelab(node=114, label="Canid Kobuvirus", bar=1, offset=.1, fontsize=3, color="black") +
  geom_cladelab(node=125, label="Feline Kobuvirus", bar=1, offset=.05, fontsize=3, color="black") + 
  geom_cladelab(node=134, label="Human Kobuvirus", bar=1, offset=.1, fontsize=3, color="black") +
  geom_cladelab(node=199, label="Caprine Kobuvirus", bar=1, offset=.1, fontsize=3, color="black") +
  geom_cladelab(node=132, label="Rodent Kobuvirus", bar=1, offset=.18, fontsize=3, color="black") 

p

# now export 

homewd= '/Users/flg9/Desktop/Developer/brook_lab/Madagascar-Bat-Kobuvirus/figures/2_kobuvirus_phylogeny/'

setwd(paste0(homewd))

ggsave(file = paste0(homewd, "Fig2b.png"),
       units="mm",  
       width=400, 
       height=250, 
       #limitsize = F,
       scale=1)#, 




# Assuming final.rooted.kobu is your tree object and kobux is your data frame

# Identify the outgroup label
outgroup_label <- "NC_026314  |  Rabovirus A  |  Rabovirus  |  Germany  |  2011"

# Ensure the outgroup label is in the tree
if (outgroup_label %in% final.rooted.kobu$tip.label) {
  # Find the index of the outgroup tip
  outgroup_index <- which(final.rooted.kobu$tip.label == outgroup_label)
  
  # Find the edge corresponding to the outgroup
  outgroup_edge <- which(final.rooted.kobu$edge[, 2] == outgroup_index)
  
  # Modify the branch length of the outgroup
  final.rooted.kobu$edge.length[outgroup_edge] <- final.rooted.kobu$edge.length[outgroup_edge] / 10  # Adjust the factor as needed
} else {
  stop("Outgroup label not found in the tree.")
}

# Identify the node
node_id <- 105

# Find the descendants of the node
# Find the descendants of the node
descendants <- getDescendants(final.rooted.kobu, node_id)

# Modify the branch lengths of all descendants (including internal nodes)
for (descendant in descendants) {
  edge_index <- which(final.rooted.kobu$edge[, 2] == descendant)
  final.rooted.kobu$edge.length[edge_index] <- final.rooted.kobu$edge.length[edge_index] / 5
}






ggtree(final.rooted.kobu) %<+% kobux + geom_tippoint(aes(color=host, fill=host, shape=bat_host)) +
  geom_nodelab(size=.5,nudge_x = -.04, nudge_y = .7) +
  scale_color_manual(values=colz) + 
  scale_fill_manual(values=colz) +
  scale_shape_manual(values=shapez) + 
  new_scale_fill() +
  geom_tiplab( aes(fill = novel), geom = "label", label.size = 0, alpha=.3, size=1.8, show.legend=F) +
  scale_fill_manual(values=colz2) + 
  theme(legend.position = c(.2,.85), legend.title = element_blank()) +
  xlim(0,70) + ggtitle("Full Genome Kobuvirus Tree") + theme_ap(family="") + 
  theme(legend.position = "bottom", legend.title = ) + geom_treescale()




########
p1 <- 
  
  ggtree(final.rooted.kobu) %<+% kobux + geom_tippoint(aes(fill= Host, shape=bat_host), shape=21) +
  geom_nodelab(size=.5,nudge_x = -.04, nudge_y = .7) + scale_fill_manual(values = colznc) + scale_shape_manual(values=shapezk) +
  new_scale_fill() + geom_tiplab( aes(fill = novel), geom = "label", label.size = 0, alpha=.3, size=1.8, show.legend=F) + 
  scale_fill_manual(values = colz2k) + geom_tiplab(aes(fill=novel), geom="label",label.size = 0, alpha=.3, size=1.8, show.legend=F) +
  new_scale_fill() + ggtitle("Kobuvirus Full Genome") + theme(legend.position = "bottom", legend.title = ) 



hostlegend <- get_legend(p1 + theme(legend.position = "bottom"))

p2 <- ggtree(final.rooted.kobu) %<+% kobux + geom_tippoint(aes(shape=bat_host)) +
  geom_nodelab(size=.5,nudge_x = -.04, nudge_y = .7) + scale_fill_manual(values = colznc) + scale_shape_manual(values=shapezk) +
  new_scale_fill() + geom_tiplab( aes(fill = novel), geom = "label", label.size = 0, alpha=.3, size=1.8, show.legend=F) + 
  scale_fill_manual(values = colz2k) + geom_tiplab(aes(fill=novel), geom="label",label.size = 0, alpha=.3, size=1.8, show.legend=F) +
  new_scale_fill() + ggtitle("Kobuvirus Phylogenetic Tree") + theme(legend.position = "bottom", legend.title =) 
 

batshapelegend <- get_legend(p2 + theme(legend.position = "bottom"))

ggtree(final.rooted.kobu) %<+% kobux + geom_tippoint(aes(fill=new_class,shape=bat_host), show.legend = F) +
  geom_nodelab(size=.5,nudge_x = -.04, nudge_y = .7) + scale_fill_manual(values = colznc) + scale_shape_manual(values=shapezk) +
  new_scale_fill() + geom_tiplab( aes(fill = novel), geom = "label", label.size = 0, alpha=.3, size=1.8, show.legend=F) + 
  scale_fill_manual(values = colz2k) + geom_tiplab(aes(fill=novel), geom="label",label.size = 0, alpha=.3, size=1.8, show.legend=F) +
  new_scale_fill() + ggtitle("Kobuvirus Phylogenetic Tree") + cowplot::get_legend(hostlegend) + cowplot::get_legend(batshapelegend) +
  theme(legend.position = "none", legend.title =)

ggtree(final.rooted.kobu) %<+% kobux + geom_tippoint(aes(fill=host,shape=bat_host), show.legend = T) +
  geom_nodelab(size=.5,nudge_x = -.04, nudge_y = .7) + scale_fill_manual(values = colznc) + scale_shape_manual(values=shapezk) +
  new_scale_fill() + geom_tiplab( aes(fill = novel), geom = "label", label.size = 0, alpha=.3, size=1.8, show.legend=F) + 
  scale_fill_manual(values = colz2k) + geom_tiplab(aes(fill=novel), geom="label",label.size = 0, alpha=.3, size=1.8, show.legend=F) +
  new_scale_fill() + ggtitle("Kobuvirus Phylogenetic Tree") + theme(legend.position = "bottom")


+ get_legend(hostlegend +theme(legend.position = "bottom")) + 
  get_legend(batshapelegend + theme(legend.position = "bottom")) 
  
