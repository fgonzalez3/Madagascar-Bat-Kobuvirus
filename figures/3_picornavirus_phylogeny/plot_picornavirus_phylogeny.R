rm(list = ls())

#load in libs 
libraries <- c("ggplot2", "ggtree", "ape", "plyr", "phytools", "phangorn", 
               "ggnewscale", "phylobase", "stringr", "tidyr", "ggpubr", 
               "grid", "gridExtra", "caper", "dplyr", "stringr")

lapply(libraries, library, character.only = TRUE)

# tree data --------------------------------------------------------------------

# load the Fig 1a tree
final.picorna <-  read.tree(file = paste0("data/T3.raxml.supportFBP"))

# change Node_4 to Genbank name
final.picorna$tip.label <- gsub("NODE_4", "OP287812", final.picorna$tip.label)

# root it using cov
final.rooted.picorna <- root(final.picorna, which(final.picorna$tip.label=='NC_048212'))


# take a quick look
ggtree(final.rooted.picorna) + 
  geom_nodelab(aes(label=label), size=1, nudge_x=-0.01, nudge_y=0.25) +
  geom_tiplab(align= FALSE, linetype="dotted", linesize = 0.1, size = 2.5) + geom_point(colour='red')


# tree metadata ----------------------------------------------------------------

# load tree data for picorna nt tree
# and fill in for any NA values
picorna.manual <- read.csv(file = paste0('data/picorna_manual.csv'),
                           header = TRUE,
                           stringsAsFactors = FALSE,
                           na.strings = "")

# use dplyr to only include columns that will be in the final tip label
# mass package also loaded, so need to explicitly call dplyr 
colnames(picorna.manual)

picorna.manual <- picorna.manual %>%
  dplyr::select(Accession_Number, Species, Genus, Geo_Location, Family, Year)

picorna.manual$Family[is.na(picorna.manual$Family)] = "NA"
picorna.manual$Family <- gsub("NA", "Unspecified", picorna.manual$Family)

picorna.manual$Genus[is.na(picorna.manual$Genus)] = "NA"
picorna.manual$Genus <- gsub("NA", "Unspecified", picorna.manual$Genus)

#check unique hosts that will be used to color tip labels 
unique_nam <- unique(picorna.manual$Genus)

sorted_nam <- sort(unique_nam)

# using genus 

C <- c("Aalivirus", "Ailurivirus", "Anativirus", "Aphthovirus", "Avihepatovirus", "Avisivirus", "Cardiovirus",
            "Cosavirus", "Crohivirus", "Dicipivirus", "Enterovirus", "Erbovirus", "Gallivirus", "Hepatovirus",
            "Hunnivirus", "Kobuvirus", "Kunsagivirus", "Limnipivirus", "Livupivirus", "Malagsivirus", "Megrivirus", "Mischivirus",
            "Mosavirus", "Oscivirus", "Parechovirus", "Passerivirus", "Poecivirus", "Potamipivirus", "Rabovirus",
            "Rafivirus", "Rosavirus", "Sakobuvirus", "Salivirus", "Sapelovirus", "Senecavirus", "Shanbavirus",
            "Tremovirus", "Unspecified")


# tree plot ------------------------------------------------------------------

df <- data.frame(Category = C, Value = 1:length(C))

# create plot of assigned genus colors to match them to summarized clades 
p <- ggplot(df, aes(x = Category, y = Value, color = Category)) +
  geom_point(size = 5) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

p

# extract the colors used in the plot
ggcolors <- scales::hue_pal()(length(C))
names(ggcolors) <- C

# print the colors and their assigned genus to match 
print(ggcolors)

# convert the phylo object to a ggtree object
ggtree_obj <- ggtree(final.rooted.picorna)

# define the tips for which you want to find the MRCA
# not a great system to doing this, but it works in defining clades I wanted to summarize
clades_to_collapse <- c("NC_026316",
                     "NC_026315")

# get the MRCA node
mrca_node <- getMRCA(final.rooted.picorna, clades_to_collapse)
                      
# and add a "novel" category
picorna.manual$novel = 0
picorna.manual$novel[picorna.manual$Accession_Number=="OP287812"] <- 1

picorna.manual$novel <- as.factor(picorna.manual$novel)

#tip shapes for bat hosts (picorna) 
picorna.manual$bat_host <- 0

picorna.manual$bat_host[picorna.manual$Accession_Number=="OP287812"] <- "bat-host"
picorna.manual$bat_host[picorna.manual$Accession_Number=="NC_038313"] <- "bat-host"
picorna.manual$bat_host[picorna.manual$Accession_Number=="NC_038316"] <- "bat-host"
picorna.manual$bat_host[picorna.manual$Accession_Number=="NC_038961"] <- "bat-host"
picorna.manual$bat_host[picorna.manual$Accession_Number=="NC_034381"] <- "bat-host"
picorna.manual$bat_host[picorna.manual$Accession_Number=="NC_033820"] <- "bat-host"
picorna.manual$bat_host[picorna.manual$Accession_Number=="NC_030843"] <- "bat-host"
picorna.manual$bat_host[picorna.manual$Accession_Number=="NC_028366"] <- "bat-host"
picorna.manual$bat_host[picorna.manual$Accession_Number=="NC_026470"] <- "bat-host"
picorna.manual$bat_host[picorna.manual$Accession_Number=="NC_015934"] <- "bat-host"
picorna.manual$bat_host[picorna.manual$Accession_Number=="NC_015940"] <- "bat-host"
picorna.manual$bat_host[picorna.manual$Accession_Number=="NC_015941"] <- "bat-host"
picorna.manual$bat_host[picorna.manual$Accession_Number=="MN602325"] <- "bat-host"
picorna.manual$bat_host[picorna.manual$Accession_Number=="MF352419"] <- "bat-host"
picorna.manual$bat_host[picorna.manual$Accession_Number=="MF352423"] <- "bat-host"
picorna.manual$bat_host[picorna.manual$Accession_Number=="MF352427"] <- "bat-host"
picorna.manual$bat_host[picorna.manual$Accession_Number=="KJ641686"] <- "bat-host"
picorna.manual$bat_host[picorna.manual$Accession_Number=="KJ641687"] <- "bat-host"
picorna.manual$bat_host[picorna.manual$Accession_Number=="KJ641691"] <- "bat-host"
picorna.manual$bat_host[picorna.manual$Accession_Number=="KJ641693"] <- "bat-host"
picorna.manual$bat_host[picorna.manual$Accession_Number=="KJ641694"] <- "bat-host"
picorna.manual$bat_host[picorna.manual$Accession_Number=="KJ641696"] <- "bat-host"
picorna.manual$bat_host[picorna.manual$Accession_Number=="HQ595341"] <- "bat-host"
picorna.manual$bat_host[picorna.manual$Accession_Number=="HQ595343"] <- "bat-host"
picorna.manual$bat_host[picorna.manual$Accession_Number=="HQ595345"] <- "bat-host"
picorna.manual$bat_host[picorna.manual$Accession_Number=="NC_048212"] <- "bat-host"

picorna.manual$bat_host[picorna.manual$bat_host==0] <- "non-bat-host"
picorna.manual$bat_host[picorna.manual$bat_host==1] <- "bat-host"

picorna.manual$bat_host <- as.factor(picorna.manual$bat_host)

# tip labels on manual csv and raxml tree do not match, so fix that 
# first make a new df with the labels numbered by root tip on the raxml tree
picorna.dat<- data.frame(Accession_Number=final.rooted.picorna$tip.label, 
                      num =1:length(final.rooted.picorna$tip.label))

# then create a 3rd df and right join the original manual and the 2nd df 
# by accession_number onto this 3rd df 
picornax <- join(picorna.manual, picorna.dat, by = "Accession_Number", match = "all", 
              type = "right")

# mask NA's and create new label
picornax$new_label <- paste(picornax$Accession_Number,
                            ifelse(is.na(picornax$Species), "", paste(" | ", picornax$Species)),
                            ifelse(is.na(picornax$Family), "", paste(" | ", picornax$Family)),
                            ifelse(is.na(picornax$Geo_Location), "", paste(" | ", picornax$Geo_Location)),
                            ifelse(is.na(picornax$Year), "", paste(" | ", picornax$Year)),
                            sep = "")

picornax$Accession_Number <- picornax$new_label

final.rooted.picorna$tip.label <- picornax$Accession_Number

final.rooted.picorna$tip.label

# label nodes based on posterior values 
# convert node labels to numeric values
final.rooted.picorna$node.label <- as.numeric(final.rooted.picorna$node.label)

# categorize node labels based on a threshold 
final.rooted.picorna$node.label[final.rooted.picorna$node.label < 90] <- 0
final.rooted.picorna$node.label[final.rooted.picorna$node.label >= 90] <- 1
final.rooted.picorna$node.label <- as.character(final.rooted.picorna$node.label)
final.rooted.picorna$node.label[final.rooted.picorna$node.label == "0"] <- "<90"
final.rooted.picorna$node.label[final.rooted.picorna$node.label == "1"] <- ">=90"
final.rooted.picorna$node.label <- factor(final.rooted.picorna$node.label, levels = c("<90", ">=90"))

# create a color mapping
posfilz <- c('<90' = "white", '>=90' = "black")

# extract the node labels and ensure they are correctly aligned with the nodes in the tree
node_labels <- data.frame(node = (Ntip(final.rooted.picorna) + 1):(Ntip(final.rooted.picorna) + Nnode(final.rooted.picorna)),
                          node_label = final.rooted.picorna$node.label)

# convert the tree to a tibble and join with node labels
tree_data <- as_tibble(final.rooted.picorna) %>%
  left_join(node_labels, by = "node")

# visualize the final tree
p2 <- ggtree(final.rooted.picorna) %<+% tree_data  %<+% picornax +  
  geom_tippoint(aes(color = Genus, fill = Genus), size = 1.5) +
  new_scale_fill() +
  geom_nodepoint(aes(fill = node_label), shape = 21, size = 1, color = "black") +  # Fill the nodes
  #geom_nodelab(aes(label = node_label), size = 2, nudge_x = -.1, nudge_y = .7) +  # Label the nodes
  scale_fill_manual(values = posfilz) + 
  new_scale_fill() +
  #scale_shape_manual(values = shapez) + 
  new_scale_fill() +
  geom_tiplab(geom = "label", label.size = 0, alpha = .3, size = 2.7, show.legend = F) +
  #scale_fill_manual(values = colz2) + 
  theme(legend.position = c(0.8, 0.8), legend.title = element_blank()) +  
  xlim(c(0, 10)) + 
  #ggtitle("Full Genome Picornavirus Tree") + 
  geom_treescale(y = 5, x = 5, fontsize = 4, offset = 1, color = "black", width = 0.1) + 
  theme(plot.margin = margin(1, 1, 1, 4, "cm"))

p2

# collapse clades 
p2 <- p2 %>%
  collapse(node = 165) + 
  geom_point2(aes(subset = (node == 165)), shape = 22, size = 1.5, fill = '#6FB000') 
p2 <- p2 %>%
  collapse(node = 205) + 
  geom_point2(aes(subset = (node == 205)), shape = 22, size = 1.5, fill = '#00BB44')
p2 <- p2 %>%
  collapse(node = 219) + 
  geom_point2(aes(subset = (node == 219)), shape = 22, size = 1.5, fill = '#ADA200')
p2 <- p2 %>%
  collapse(node = 227) + 
  geom_point2(aes(subset = (node == 227)), shape = 22, size = 1.5, fill = '#BC9D00')
p2 <- p2 %>%
  collapse(node = 252) + 
  geom_point2(aes(subset = (node == 252)), shape = 22, size = 1.5, fill = '#00B813')
p2 <- p2 %>%
  collapse(node = 273) + 
  geom_point2(aes(subset = (node == 273)), shape = 22, size = 1.5, fill = '#E26EF7')
p2 <- p2 %>%
  collapse(node = 280) + 
  geom_point2(aes(subset = (node == 280)), shape = 22, size = 1.5, fill = '#D177FF')
p2 <- p2 %>%
  collapse(node = 275) + 
  geom_point2(aes(subset = (node == 275)), shape = 22, size = 1.5, fill = '#00C08E')
p2 <- p2 %>%
  collapse(node = 239) + 
  geom_point2(aes(subset = (node == 239)), shape = 22, size = 1.5, fill = '#00BD61')
p2 <- p2 %>%
  collapse(node = 238) + 
  geom_point2(aes(subset = (node == 238)), shape = 22, size = 1.5, fill = '#00BFC4')

# rename summarized clades 
p2 <- p2 + 
  geom_text2(aes(subset = (node == 165), label = "Enterovirus - Collapsed Clade"), 
             nudge_x = 0.65, nudge_y = 0.01, size = 2.7, color = "black") +
  geom_text2(aes(subset = (node == 205), label = "Hepatovirus - Collapsed Clade"), 
             nudge_x = 0.65, nudge_y = 0.01, size = 2.7, color = "black") + 
  geom_text2(aes(subset = (node == 219), label = "Cosavirus - Collapsed Clade"), 
             nudge_x = 0.6, nudge_y = 0.01, size = 2.7, color = "black") +
  geom_text2(aes(subset = (node == 227), label = "Cardiovirus & Senecavirus - Collapsed Clade"), 
             nudge_x = 0.95, nudge_y = 0.01, size = 2.7, color = "black") + 
  geom_text2(aes(subset = (node == 252), label = "Gallivirus - Collapsed Clade"), 
             nudge_x = 0.6, nudge_y = 0.01, size = 2.7, color = "black") +
  geom_text2(aes(subset = (node == 273), label = "Rosavirus - Collapsed Clade"), 
           nudge_x = 0.6, nudge_y = 0.01, size = 2.7, color = "black") +
  geom_text2(aes(subset = (node == 280), label = "Rafivirus - Collapsed Clade"), 
             nudge_x = 0.6, nudge_y = 0.01, size = 2.7, color = "black") +  
  geom_text2(aes(subset = (node == 275), label = "Kunsagavirus - Collapsed Clade"), 
             nudge_x = 0.7, nudge_y = 0.01, size = 2.7, color = "black") +
  geom_text2(aes(subset = (node == 239), label = "Hunnivirus - Collapsed Clade"), 
             nudge_x = 0.65, nudge_y = 0.01, size = 2.7, color = "black") +
  geom_text2(aes(subset = (node == 238), label = "Malagsivirus - Collapsed Clade"), 
             nudge_x = 0.65, nudge_y = 0.01, size = 2.7, color = "black")
  
p2

ggsave("Fig2.png", plot = p2, width = 15, height = 10, units = "in", dpi = 300)

