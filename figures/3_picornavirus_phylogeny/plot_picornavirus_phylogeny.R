rm(list = ls())

#load in libs 
libraries <- c("ggplot2", "ggtree", "ape", "plyr", "phytools", "phangorn", 
               "ggnewscale", "phylobase", "stringr", "tidyr", "ggpubr", 
               "grid", "gridExtra", "caper", "dplyr", "stringr")

lapply(libraries, library, character.only = TRUE)

picornavirus_treefile <- "./data/PicornavirusNorovirusOutgroup20DEC2024.raxml.supportFBP"
picornavirus_metadata <- "./data/picorna_manual.csv"

outpng <- "./Fig2.png"

# tree plot --------------------------------------------------------------------

get_tree_plot <- function(treefile, metadata) {
  
  # load in tree 
  tree <- read.tree(file=paste0(treefile))
  
  # root tree with outgroup and drop tip for cleaner visualization
  tree <- root(tree, outgroup = 'NC_044856')
  tree <- drop.tip(tree, 'NC_024766')
  
  # load in metadata
  treedata <- read.csv(file=paste0(metadata),
                       header=T, stringsAsFactors = F, na.strings = "")
  
  # use dplyr to only include columns that will be in the final tip label
  treedata <- treedata %>%
    dplyr::select(Accession_Number, Species, Genus, Geo_Location, Family, Year)
  
  # replace virus families and genus that are NA with unspecified terminology
  treedata$Family[is.na(treedata$Family)] = "NA"
  treedata$Family <- gsub("NA", "Unspecified", treedata$Family)
  
  treedata$Genus[is.na(treedata$Genus)] = "NA"
  treedata$Genus <- gsub("NA", "Unspecified", treedata$Genus)
  
  #check unique hosts that will be used to color tip labels 
  unique_nam <- unique(treedata$Genus)
  sorted_nam <- sort(unique_nam)
  
  # using genus as the legend labels
  C <- c("Aalivirus", "Ailurivirus", "Anativirus", "Aphthovirus", "Avihepatovirus", "Avisivirus", "Cardiovirus",
         "Cosavirus", "Crohivirus", "Dicipivirus", "Enterovirus", "Erbovirus", "Gallivirus", "Hepatovirus",
         "Hunnivirus", "Kobuvirus", "Kunsagivirus", "Limnipivirus", "Livupivirus", "Malagsivirus", "Megrivirus", "Mischivirus",
         "Mosavirus", "Oscivirus", "Parechovirus", "Passerivirus", "Poecivirus", "Potamipivirus", "Rabovirus",
         "Rafivirus", "Rosavirus", "Sakobuvirus", "Salivirus", "Sapelovirus", "Senecavirus", "Shanbavirus",
         "Tremovirus", "Unspecified")
  
  # create a df for the legend labels 
  df <- data.frame(Category = C, Value = 1:length(C))
  
  # create plot of assigned genus colors to match them to summarized clades 
  p <- ggplot(df, aes(x = Category, y = Value, color = Category)) +
    geom_point(size = 5) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  
  # extract the colors used in the plot
  ggcolors <- scales::hue_pal()(length(C))
  names(ggcolors) <- C
  
  # print the colors and their assigned genus to match 
  print(ggcolors)
  
  # convert the phylo object to a ggtree object
  ggtree_obj <- ggtree(tree)
  
  # define the tips for which you want to find the MRCA
  # not a great system to doing this, but it works in defining clades I wanted to summarize
  clades_to_collapse <- c("NC_026315",
                          "NC_026316")
  
  # get the MRCA node
  mrca_node <- getMRCA(tree, clades_to_collapse)
  
  # tip labels on manual csv and raxml tree do not match, so fix that 
  # first make a new df with the labels numbered by root tip on the raxml tree
  picorna.dat<- data.frame(Accession_Number=tree$tip.label, 
                           num =1:length(tree$tip.label))
  
  # then create a 3rd df and right join the original manual and the 2nd df 
  # by accession_number onto this 3rd df 
  picornax <- join(treedata, picorna.dat, by = "Accession_Number", match = "all", 
                   type = "right")
  
  # mask NA's and create new label
  picornax$new_label <- paste(picornax$Accession_Number,
                              ifelse(is.na(picornax$Species), "", paste(" | ", picornax$Species)),
                              ifelse(is.na(picornax$Family), "", paste(" | ", picornax$Family)),
                              ifelse(is.na(picornax$Geo_Location), "", paste(" | ", picornax$Geo_Location)),
                              ifelse(is.na(picornax$Year), "", paste(" | ", picornax$Year)),
                              sep = "")
  
  picornax$Accession_Number <- picornax$new_label
  
  tree$tip.label <- picornax$Accession_Number
  
  tree$tip.label
  
  # label nodes based on posterior values 
  # convert node labels to numeric values
  tree$node.label <- as.numeric(tree$node.label)
  
  # categorize node labels based on a threshold 
  tree$node.label[tree$node.label < 90] <- 0
  tree$node.label[tree$node.label >= 90] <- 1
  tree$node.label <- as.character(tree$node.label)
  tree$node.label[tree$node.label == "0"] <- "<90"
  tree$node.label[tree$node.label == "1"] <- ">=90"
  tree$node.label <- factor(tree$node.label, levels = c("<90", ">=90"))
  
  # create a color mapping
  posfilz <- c('<90' = "white", '>=90' = "black")
  
  # normalize the node labels to a range between 0 and 1
  tree$node.label <- as.numeric(tree$node.label)
  
  min_val <- min(tree$node.label, na.rm = TRUE)
  max_val <- max(tree$node.label, na.rm = TRUE)
  tree$node.label <- (tree$node.label - min_val) / (max_val - min_val) * 100
  
  
  # extract the node labels and ensure they are correctly aligned with the nodes in the tree
  node_labels <- data.frame(node = (Ntip(tree) + 1):(Ntip(tree) + Nnode(tree)),
                            node_label = tree$node.label)
  
  # convert the tree to a tibble and join with node labels
  tree_data <- as_tibble(tree) %>%
    left_join(node_labels, by = "node")
  
  # visualize the final tree
  treeplot <- ggtree(tree) %<+% tree_data  %<+% picornax +  
    geom_tippoint(aes(color = Genus, fill = Genus), size = 1.5) +
    new_scale_fill() +
    geom_nodepoint(aes(fill = node_label), shape = 21, size = 1, color = "black") +  # Fill the nodes
    scale_fill_gradient(low = "white", high = "black", breaks = seq(0, 100, by = 20), labels = seq(0, 100, by = 20),
                        guide = guide_colorbar(direction = "horizontal")) + 
    #geom_nodelab(aes(label = node_label), size = 2, nudge_x = -.1, nudge_y = .7) +  # Label the nodes
    #scale_fill_manual(values = posfilz) + 
    new_scale_fill() +
    #scale_shape_manual(values = shapez) + 
    new_scale_fill() +
    geom_tiplab(geom = "label", label.size = 0, alpha = .3, size = 2.7, show.legend = F) +
    #scale_fill_manual(values = colz2) + 
    theme(legend.position = c(0.8, 0.4), legend.title = element_blank()) +  
    xlim(c(0, 10)) + 
    #ggtitle("Full Genome Picornavirus Tree") + 
    geom_treescale(y = 5, x = 5, fontsize = 4, offset = 1, color = "black", width = 0.1) + 
    theme(plot.margin = margin(1, 1, 1, 4, "cm"))
  
  treeplot
  
  # collapse clades 
  treeplot <- treeplot %>%
    collapse(node = 183) + 
    geom_point2(aes(subset = (node == 183)), shape = 22, size = 1.5, fill = '#6FB000') 
  treeplot <- treeplot %>%
    collapse(node = 158) + 
    geom_point2(aes(subset = (node == 158)), shape = 22, size = 1.5, fill = '#00BB44')
  treeplot <- treeplot %>%
    collapse(node = 277) + 
    geom_point2(aes(subset = (node == 277)), shape = 22, size = 1.5, fill = '#ADA200')
  treeplot <- treeplot %>%
    collapse(node = 269) + 
    geom_point2(aes(subset = (node == 269)), shape = 22, size = 1.5, fill = '#BC9D00')
  treeplot <- treeplot %>%
    collapse(node = 257) + 
    geom_point2(aes(subset = (node == 257)), shape = 22, size = 1.5, fill = '#00B813')
  treeplot <- treeplot %>%
    collapse(node = 231) + 
    geom_point2(aes(subset = (node == 231)), shape = 22, size = 1.5, fill = '#E26EF7')
  treeplot <- treeplot %>%
    collapse(node = 265) + 
    geom_point2(aes(subset = (node == 265)), shape = 22, size = 1.5, fill = '#D177FF')
  treeplot <- treeplot %>%
    collapse(node = 252) + 
    geom_point2(aes(subset = (node == 252)), shape = 22, size = 1.5, fill = '#F763DF')
  treeplot <- treeplot %>%
    collapse(node = 176) + 
    geom_point2(aes(subset = (node == 176)), shape = 22, size = 1.5, fill = '#00C1A2')
  treeplot <- treeplot %>%
    collapse(node = 144) + 
    geom_point2(aes(subset = (node == 144)), shape = 22, size = 1.5, fill = '#E08B00')
  treeplot <- treeplot %>%
    collapse(node = 263) + 
    geom_point2(aes(subset = (node == 263)), shape = 22, size = 1.5, fill = '#00AFF8')
  treeplot <- treeplot %>%
    collapse(node = 151) + 
    geom_point2(aes(subset = (node == 151)), shape = 22, size = 1.5, fill = '#00BFC4')
  treeplot <- treeplot %>%
    collapse(node = 150) + 
    geom_point2(aes(subset = (node == 150)), shape = 22, size = 1.5, fill = '#00BD61')
  treeplot <- treeplot %>%
    collapse(node = 276) + 
    geom_point2(aes(subset = (node == 276)), shape = 22, size = 1.5, fill = '#F8766D')
  
  
  # rename summarized clades 
  treeplot <- treeplot + 
    geom_text2(aes(subset = (node == 183), label = "Enterovirus - Collapsed Clade"), 
               nudge_x = 0.65, nudge_y = 0.01, size = 2.7, color = "black") +
    geom_text2(aes(subset = (node == 158), label = "Hepatovirus - Collapsed Clade"), 
               nudge_x = 0.65, nudge_y = 0.01, size = 2.7, color = "black") + 
    geom_text2(aes(subset = (node == 277), label = "Cosavirus - Collapsed Clade"), 
               nudge_x = 0.6, nudge_y = 0.01, size = 2.7, color = "black") +
    geom_text2(aes(subset = (node == 269), label = "Cardiovirus & Senecavirus - Collapsed Clade"), 
               nudge_x = 0.95, nudge_y = 0.01, size = 2.7, color = "black") + 
    geom_text2(aes(subset = (node == 257), label = "Gallivirus - Collapsed Clade"), 
               nudge_x = 0.6, nudge_y = 0.01, size = 2.7, color = "black") +
    geom_text2(aes(subset = (node == 231), label = "Rosavirus & Kunsagivirus - Collapsed Clade"), 
               nudge_x = 0.92, nudge_y = 0.01, size = 2.7, color = "black") +
    geom_text2(aes(subset = (node == 265), label = "Rafivirus - Collapsed Clade"), 
               nudge_x = 0.6, nudge_y = 0.01, size = 2.7, color = "black") +  
    geom_text2(aes(subset = (node == 252), label = "Salivirus - Collapsed Clade"), 
               nudge_x = 0.6, nudge_y = 0.01, size = 2.7, color = "black") +
    geom_text2(aes(subset = (node == 176), label = "Limnipivirus - Collapsed Clade"), 
               nudge_x = 0.65, nudge_y = 0.01, size = 2.7, color = "black") +
    geom_text2(aes(subset = (node == 144), label = "Apthovirus - Collapsed Clade"), 
               nudge_x = 0.65, nudge_y = 0.01, size = 2.7, color = "black") + 
    geom_text2(aes(subset = (node == 263), label = "Oscivirus - Collapsed Clade"), 
               nudge_x = 0.62, nudge_y = 0.01, size = 2.7, color = "black") + 
    geom_text2(aes(subset = (node == 151), label = "Malagsivirus - Collapsed Clade"), 
               nudge_x = 0.65, nudge_y = 0.01, size = 2.7, color = "black") +
    geom_text2(aes(subset = (node == 150), label = "Hunnivirus - Collapsed Clade"), 
               nudge_x = 0.65, nudge_y = 0.01, size = 2.7, color = "black") + 
    geom_text2(aes(subset = (node == 276), label = "Ailurivirus - Collapsed Clade"), 
               nudge_x = 0.65, nudge_y = 0.01, size = 2.7, color = "black")
  
  return(treeplot)
    
}

# assemble plot ----------------------------------------------------------------

get_tree_plot(picornavirus_treefile, picornavirus_metadata)
ggsave(outpng, plot = p2, width = 15, height = 11, units = "in", dpi = 300)


# misc code --------------------------------------------------------------------

# non-utilized code
p2 <- ggtree(tree) %<+% tree_data %<+% picornax +  
  geom_tippoint(aes(color = Genus, fill = Genus), size = 1.5) +
  new_scale_fill() +
  geom_nodepoint(aes(fill = node_label), shape = 21, size = 1, color = "black") +  # Fill the nodes
  scale_fill_gradient(low = "white", high = "black", breaks = seq(0, 100, by = 20), labels = seq(0, 100, by = 20),
                      guide = guide_colorbar(direction = "horizontal")) + 
  geom_tiplab(geom = "label", label.size = 0, alpha = .3, size = 2.7, show.legend = F) +
  theme(legend.position = c(0.8, 0.4), legend.title = element_blank()) +  
  xlim(c(0, 10)) + 
  geom_treescale(y = 5, x = 6, fontsize = 4, offset = 1, color = "black", width = 0.1) + 
  theme(plot.margin = margin(1, 1, 1, 4, "cm"))
