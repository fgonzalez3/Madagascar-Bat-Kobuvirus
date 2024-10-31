###### Beast Tree ########

# rm any remaining environments
rm(list = ls())

library(ggplot2)
library(ggtree)
library(ape)
library(ggnewscale)
library(lubridate)
library(treeio)
library(rBt)
library(stringr)
library(ape)
library(ggmsa)
library(seqinr)
library(dplyr)

# !remove path when uploading to git!


# load data --------------------------------------------------------------------

# first, read in the tree
homewd= "/Users/flg9/Desktop/Developer/brook_lab/Madagascar-Bat-Kobuvirus/results/BayesianTree/" 
setwd(paste0(homewd))

beast_tree <- read.beast(file = "Kobuvirus.MCC.tree")

# visualize initial tree structure with variance
ggtree(beast_tree) + 
  geom_nodelab(aes(label=label), size=1, nudge_x=-0.01, nudge_y=0.25) +
  geom_range("height_0.95_HPD", color='red', size=2, alpha=.5) + 
  geom_tiplab(align= FALSE, linetype="dotted", linesize = 0.1, size = 2.5) + geom_point(colour='red')

# now create a df that contains the accessions fromt the tree
# use cbind to do this, which pulls those accessions from tip label data 
# then create a new column, which we will modify
treedat <- cbind.data.frame(tip_name = beast_tree@phylo$tip.label)
treedat$beast_name <-treedat$tip_name

# now remove the capture data from the accessions 
treedat$accession_num <- sapply(strsplit(treedat$tip_name, "_"), function(x) x[[1]])

# and rename some of them
treedat$accession_num[treedat$accession_num=="NC"] <- c("NC_015936", "NC_027054", 
                                                        "NC_004421", "NC_027918", 
                                                        "NC_027919", "NC_026314")

# load the manual tree data that i created 
# this will be used to create new tip labels for the tree
dat <- read.csv(file = "BEAST_Metadata.csv", 
                header = T, stringsAsFactors = F)

# and omit any na 
dat<-na.omit(dat)

# change collection date to a more readable format 
dat$collection_date <- as.Date(dat$collection_date)

# visualize with some new elements 
mrsd.dat <- max(dat$collection_date) # our most recent sequence

p1 <- ggtree(beast_tree, mrsd=mrsd.dat) + theme_tree2() #+ geom_range("height_0.95_HPD", color='red', size=2, alpha=.5)
p1

tree.dat <- p1$data
node.sub <- dplyr::select(tree.dat, node, x)
names(node.sub) <- c("node", "nodetime")
print(node.sub)

# determine the range of node times
min_time <- min(node.sub$nodetime)
max_time <- max(node.sub$nodetime)

# set x-axis breaks and label dynamically based on above range
#breaks <- seq(from = floor(min_time), to = ceiling(max_time), by = 2000)
#labels <- abs(breaks - max_time)
#labels <- round(labels)

# label nodes based on posterior values 
beast_tree@data$posterior[beast_tree@data$posterior<0.9] <- 0
beast_tree@data$posterior[beast_tree@data$posterior>=0.9] <- 1
beast_tree@data$posterior <- round(beast_tree@data$posterior, 2)
beast_tree@data$posterior <- as.character(beast_tree@data$posterior)
beast_tree@data$posterior[beast_tree@data$posterior=="0"] <- "<0.9"
beast_tree@data$posterior[beast_tree@data$posterior=="1"] <- "=>0.9"
beast_tree@data$posterior <- factor(beast_tree@data$posterior, levels= c("<0.9", "=>0.9"))
posfilz <- c('<0.9'="white", '=>0.9'="black")


# date nodes of interest -------------------------------------------------------

# define node of interest
nodeOP287812 <- MRCA(beast_tree, which(beast_tree@phylo$tip.label == "OP287812_2018-07-27"), 
                     which(beast_tree@phylo$tip.label == "FJ890523_2008-07-15"))

nodeKobuVs <- MRCA(beast_tree, which(beast_tree@phylo$tip.label == "OP287812_2018-07-27"), 
                     which(beast_tree@phylo$tip.label == "KJ641691_2012-12-15"))

# find date that OP87812 branches off 
op287812.date <- round(node.sub$nodetime[nodeOP287812],0)
recent.age <- year(mrsd.dat) + yday(mrsd.dat)/365
op287812.mean <- round(recent.age - beast_tree@data$height[51], 0)

# upper and lower confidence intervals for this date 
op287812.uci <- round(recent.age-beast_tree@data$height_0.95_HPD[51][[1]][1],0)
op287812.lci <- round(recent.age-beast_tree@data$height_0.95_HPD[51][[1]][2],0)

# format the date for the node
OP287812.date <- paste0("~", op287812.mean, "\n[", op287812.lci, "-", op287812.uci, "]")

new.nodel.lab <- rep(NA, nrow(node.sub))
#new.nodel.lab[nodeall] <- paste0("~",orig.date)
new.nodel.lab[nodeOP287812] <- paste0("~", op287812.date)


# output the calculated dates
list(
  op287812_mean = op287812.mean,
  op287812_uci = op287812.uci,
  op287812_lci = op287812.lci,
  op287812_date = op287812.date
)

# find date that OP87812 branches off 
KobuClade.date <- round(node.sub$nodetime[nodeKobuVs],0)
recent.age <- year(mrsd.dat) + yday(mrsd.dat)/365
KobuClade.mean <- round(recent.age - beast_tree@data$height[47L], 0)

# upper and lower confidence intervals for this date 
KobuClade.uci <- round(recent.age-beast_tree@data$height_0.95_HPD[47L][[1]][1],0)
KobuClade.lci <- round(recent.age-beast_tree@data$height_0.95_HPD[47L][[1]][2],0)

# format the date for the node
OP287812.date <- paste0("~", op287812.mean, "\n[", op287812.lci, "-", op287812.uci, "]")

new.nodel.lab <- rep(NA, nrow(node.sub))
#new.nodel.lab[nodeall] <- paste0("~",orig.date)
new.nodel.lab[nodeOP287812] <- paste0("~", op287812.date)

# rename tips and plot ---------------------------------------------------------

# head tree data and make new clade column  
head(dat)

dat$Clade <- dat$type

# change col name in tree.dat from label to accession_num
# and merge treedat and dat by accession number 

# colnames takes df and changes name by column position within df 
colnames(tree.dat)[4] <- "accession_num"

dat.plot <- merge(treedat, dat, by="accession_num", all.x = T, sort=F)

# now make new labels for tree tips using paste function
head(dat.plot)
dat.plot$new_label = ""

dat.plot$new_label[!is.na(dat.plot$type)] <- paste(dat.plot$accession_num[!is.na(dat.plot$type)], " | ", 
                                                     #dat.plot$type[!is.na(dat.plot$type)], " | ", 
                                                     #dat.plot$host[!is.na(dat.plot$type)], " | ",
                                                     dat.plot$title[!is.na(dat.plot$title)], " | ",
                                                     dat.plot$country[!is.na(dat.plot$type)], " | ",
                                                     dat.plot$collection_year[!is.na(dat.plot$type)])

# check to make sure labels are fine 
dat.plot$new_label

# hard code tip label on original tree file to be new label on dat.plot
beast_tree@phylo$tip.label <- dat.plot$new_label

# create new df that only contains columns we need for final tree
dat.sub <- dplyr::select(dat.plot, new_label, collection_date, country, Clade, title)
head(dat.sub) # check

# change dat.sub$clade from character to factor
dat.sub$title <- as.factor(dat.sub$title)

# categorize bat kobuvirus sequences and highlight
dat.sub$novel = "no"
dat.sub$novel[dat.sub$country=="Madagascar"] <- "yes"
dat.sub$novel[dat.sub$new_label=="KJ641686  |  Bat Kobuvirus  |  China  |  2010"] <- "no_but_bat"
dat.sub$novel[dat.sub$new_label=="KJ641691  |  Bat Kobuvirus  |  China  |  2012"] <- "no_but_bat"

colz2 = c('yes' =  "yellow", 'no' = "white", 'no_but_bat' = "purple")



ggtree(beast_tree, mrsd=mrsd.dat, size=.8) %<+% dat.sub +
  geom_tippoint(aes(color=title), size=3.8) +
  theme_tree2() +
  theme(legend.position = c(.01,.7), plot.margin = unit(c(.2,24,3,3), "lines"), legend.title = element_blank()) +
  coord_cartesian(clip = "off") + 
  geom_nodepoint(aes(fill=posterior), shape=21, color="black", size=3, stroke=.1) +
  guides(fill_continuous = guide_legend(order = 2),col = guide_legend(order = 1)) +
  scale_x_continuous(breaks=c(1620, 1720, 1400, 0, 2020),
                     labels=c(400, 300, 200, 100, 0)) +
  xlab("years to MRCA") +
  ggnewscale::new_scale_fill() +
  geom_tiplab(aes(fill = novel), geom = "label", label.size = 0, alpha=.3,  show.legend=F, size=5) + scale_fill_manual(values=colz2) 
#xlim(c(0,6))


p3 <- ggtree(beast_tree, mrsd=mrsd.dat, size=.8) %<+% dat.sub +
  #geom_range(range = 'height_0.95_HPD', color = "red", alpha = .6, size = 0.5) + 
  geom_tippoint(aes(color=title), size=4) + 
  #geom_text2(aes(label=round(as.numeric(posterior), 2), 
  #  subset=as.numeric(posterior)> 0.9, 
  # x=branch), vjust=0) + 
  geom_nodelab(aes(label=new.nodel.lab), size=4, nudge_x = -750, nudge_y = 1, color="firebrick", fontface=2, geom="label", fill="white") +
  scale_fill_manual(values = posfilz) + 
  geom_nodepoint(aes(fill = posterior), shape = 21, color = "black", size = 2, stroke = .1) +
  theme_tree2() +
  theme(axis.text.x = element_text(size = 10),  
        legend.position = c(.04, .81), 
        plot.margin = unit(c(2, 20, 2, 3), "lines"),  # Increase plot margins
        legend.title = element_blank()) +
  scale_x_continuous(breaks = breaks, labels = labels) + 
  ggnewscale::new_scale_fill() +
  xlab("Years to MRCA") +
  geom_tiplab(aes(fill = novel), geom = "label", label.size = 0, 
              alpha=.3, show.legend=F, size=3) + 
  scale_fill_manual(values=colz2) +
                       coord_cartesian(clip = "off") #+ ggtitle("Kobuvirus Types")

p3

# now save plot 
homewd="/Users/flg9/Desktop/Developer/brook_lab/Madagascar-Bat-Kobuvirus/figures/4_bayesian_phylogeny/"
setwd(paste0(homewd))

ggsave(file = paste0(homewd, "Fig3BayesianTree.png"),
       plot = p3, 
       units="mm",  
       width=150, 
       height=100, 
       #limitsize = F,
       scale=3)#, 
