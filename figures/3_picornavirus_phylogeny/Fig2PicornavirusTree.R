####### Figure 2 ########
###### Picornavirus Full Genome Tree ##########

# this r snippet contains code for figure 2 creation 
# figure 2 is our picornavirus full genome tree 

# rm any remaining environments
rm(list = ls())

#load in libs 
libraries <- c("ggplot2", "ggtree", "ape", "plyr", "phytools", "phangorn", 
               "ggnewscale", "phylobase", "stringr", "tidyr", "ggpubr", 
               "grid", "gridExtra", "caper", "dplyr", "stringr")

lapply(libraries, library, character.only = TRUE)

# set wd
homewd= '/Users/flg9/Desktop/Developer/brook_lab/Madagascar-Bat-Kobuvirus/figures/3_picornavirus_phylogeny/'
setwd(paste0(homewd))

# load the Fig 1a tree
final.picorna <-  read.tree(file = paste0(homewd, "T3.raxml.supportFBP"))

# change Node_4 to Genbank name
final.picorna$tip.label <- gsub("NODE_4", "OP287812", final.picorna$tip.label)

# root it using cov
final.rooted.picorna <- root(final.picorna, which(final.picorna$tip.label=='NC_048212'))


# take a quick look
ggtree(final.rooted.picorna) + 
  geom_nodelab(aes(label=label), size=1, nudge_x=-0.01, nudge_y=0.25) +
  geom_tiplab(align= FALSE, linetype="dotted", linesize = 0.1, size = 2.5) + geom_point(colour='red')


# load tree data for picorna nt tree
# and fill in for any NA values
picorna.manual <- read.csv(file = paste0(homewd, 'picorna_manual.csv'),
                           header = TRUE,
                           stringsAsFactors = FALSE,
                           na.strings = "")


#picorna.manual[is.na(picorna.manual)] = "NA"

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

colz = c("Human" = "red", "Bovine" = "tomato", "Porcine" = "khaki", 
         "Canine " = "blue", "Rodent" = "plum", "Caprine" = "darkolivegreen1", 
         "Feline" = "purple", "Avian" = "darksalmon", "Bat" = "cyan", 
         "Ursid" = "cadetblue1", "Reptile" = "mediumvioletred", 
         "Erinaceidae" = "darkseagreen1", "Amphibian" = "yellow", 
         "Primate" = "pink", "Camelid" = "brown", "Equine" = "darkgreen", 
         "Fish" = "thistle", "Shrew" = "sienna", "Seal" = "wheat",
         "Marsupial" = "peru", "NA" = "black")

#pick order for the labels
C = c("Human", "Bovine", "Porcine", "Canine", "Rodent", "Caprine", "Feline", 
      "Avian", "Bat", "Ursid", "Reptile", "Erinaceidae", "Amphibian", "Primate", 
      "Camelid", "Equine", "Fish", "Shrew", "Seal", "Marsupial", "NA")

# using genus 

C <- c("Aalivirus", "Ailurivirus", "Anativirus", "Aphthovirus", "Avihepatovirus", "Avisivirus", "Cardiovirus",
            "Cosavirus", "Crohivirus", "Dicipivirus", "Enterovirus", "Erbovirus", "Gallivirus", "Hepatovirus",
            "Hunnivirus", "Kobuvirus", "Kunsagivirus", "Limnipivirus", "Livupivirus", "Malagsivirus", "Megrivirus", "Mischivirus",
            "Mosavirus", "Oscivirus", "Parechovirus", "Passerivirus", "Poecivirus", "Potamipivirus", "Rabovirus",
            "Rafivirus", "Rosavirus", "Sakobuvirus", "Salivirus", "Sapelovirus", "Senecavirus", "Shanbavirus",
            "Tremovirus", "Unspecified")

df <- data.frame(Category = C, Value = 1:length(C))

# Create a plot with the categories
p <- ggplot(df, aes(x = Category, y = Value, color = Category)) +
  geom_point(size = 5) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# Print the plot
print(p)

# Extract the colors used in the plot
ggcolors <- scales::hue_pal()(length(C))
names(ggcolors) <- C

# Print the colors
print(ggcolors)


# Convert the phylo object to a ggtree object
ggtree_obj <- ggtree(final.rooted.picorna)

# Define the tips for which you want to find the MRCA
clades_to_collapse <- c("NC_026316",
                     "NC_026315")

# Get the MRCA node
mrca_node <- getMRCA(final.rooted.picorna, clades_to_collapse)
                      
#and add a "novel" category
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
shapez = c("bat-host" =  24, "non-bat-host" = 21)
colz2 = c('1' =  "yellow", '0' = "white")


# clade labels
#AichiA <- getMRCA(final.rooted.kobu, c("OP287812", "KJ934637"))
#AichiB <- getMRCA(final.rooted.kobu, c("MN336260", "GU245693"))
#AichiC <- getMRCA(final.rooted.kobu, c("NC_023422", "KY234500"))
#AichiD <- getMRCA(final.rooted.kobu, c("NC_027918", "NC_027919"))
#AichiEF <- getMRCA(final.rooted.kobu, c("KJ641686", "KT325852"))


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

#picornax$new_label<- str_replace_all(picornax$new_label, "NA", "")
#picornax$new_label<- str_replace_all(picornax$new_label, "| ", "")

# label nodes based on posterior values 
# Assuming final.rooted.picorna is your ML tree
# Convert node labels to numeric values
final.rooted.picorna$node.label <- as.numeric(final.rooted.picorna$node.label)

# Categorize node labels based on a threshold (e.g., 90)
final.rooted.picorna$node.label[final.rooted.picorna$node.label < 90] <- 0
final.rooted.picorna$node.label[final.rooted.picorna$node.label >= 90] <- 1
final.rooted.picorna$node.label <- as.character(final.rooted.picorna$node.label)
final.rooted.picorna$node.label[final.rooted.picorna$node.label == "0"] <- "<90"
final.rooted.picorna$node.label[final.rooted.picorna$node.label == "1"] <- ">=90"
final.rooted.picorna$node.label <- factor(final.rooted.picorna$node.label, levels = c("<90", ">=90"))

# Create a color mapping
posfilz <- c('<90' = "white", '>=90' = "black")

# Extract the node labels and ensure they are correctly aligned with the nodes in the tree
node_labels <- data.frame(node = (Ntip(final.rooted.picorna) + 1):(Ntip(final.rooted.picorna) + Nnode(final.rooted.picorna)),
                          node_label = final.rooted.picorna$node.label)

# Convert the tree to a tibble and join with node labels
tree_data <- as_tibble(final.rooted.picorna) %>%
  left_join(node_labels, by = "node")

# Visualize the final tree
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
  geom_treescale(y = -5, fontsize = 4, offset = 1, color = "black", width = 0.1) + 
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
  
  
#p2 <- collapse(p2, node=205, mode="min", fill="seagreen3", color="black", size=0.1) 

# collapse clades
#p2 <- p2 + geom_cladelab(node=165, label="Enterovirus", bar=1, offset=.48, fontsize=3, color="black") +
 # geom_cladelab(node=205, label="Hepatovirus", bar=1, offset=.48, fontsize=3, color="black") 
p2


ggsave("large_tree_plot.png", plot = p2, width = 15, height = 10, units = "in", dpi = 300)


# take a closer look at kobuvirus clade
# to better visualize relationships for reader
getMRCA(final.rooted.picorna, c("NC_016769 |  Aichivirus C |  Porcine |  China |  2010", 
                                "KJ641686 |  Aichivirus F |  Bat |  China |  2010"))

clade_tree <- extract.clade(final.rooted.picorna, node = 254)
fig2b <- ggtree(clade_tree) %<+% picornax + 
  geom_tippoint(aes(color=Genus, fill=Genus, shape=bat_host), size=1) +
  geom_nodelab(size=2, nudge_x = -.1, nudge_y = .7) + 
  scale_shape_manual(values=shapez) + 
  new_scale_fill() +
  geom_tiplab(aes(fill = novel), geom = "label", label.size = 0, 
              alpha=.3, size=1.8, show.legend=F) +
  scale_fill_manual(values=colz2) + 
  theme(legend.position = c(0.6, 0.5), legend.title = element_blank()) +  
  xlim(c(0,10)) + 
  ggtitle("Kobuvirus Clade") + 
 #geom_treescale(y = -5, fontsize = 5, offset = 1, color = "black", width = 0.5) + 
  theme(plot.margin = margin(1, 1, 1, 4, "cm")) 
fig2b

# merge the two trees 
fig2z <- cowplot::plot_grid(fig2a, fig2b, ncol=2, 
                            nrow=1, labels = c("(A)", "(B)"), 
                            label_size = 22, label_x = .03, label_y = .98)
fig2z


fig2 <- cowplot::plot_grid(p2, p, 
                           nrow=1, ncol=2)

fig2

homewd="/Users/flg9/Desktop/Developer/brook_lab/Madagascar-Bat-Kobuvirus/figures/3_picornavirus_phylogeny/"
setwd(paste0(homewd))
 
ggsave(file = paste0(homewd, "Fig2_Picornavirus_Kobuvirus_Tree.png"),
        plot = fig2,
        units="mm",  
        width=170, 
        height=100,
       dpi = 320,
        #limitsize = F,
        scale=3)#,
 

