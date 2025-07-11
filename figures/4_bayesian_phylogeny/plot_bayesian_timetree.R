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

kobuvirus_beastfile <- "./data/KobuvirusLuSequencesBeast.MCMC.tree"
kobuvirus_beastmetadata <- "./data/BEAST_Metadata_LuSequences18DEC2024.csv"

outpng <- "./Fig3.png"

# tree plot --------------------------------------------------------------------

get_beast_tree <- function(treedata, metadata) {
  
  # read in beast file
  beast_tree <- read.beast(treedata)
  
  # now create a df that contains the accessions from the tree
  # use cbind to do this, which pulls those accessions from tip label data 
  # then create a new column, which we will modify
  treedat <- cbind.data.frame(tip_name = beast_tree@phylo$tip.label)
  treedat$beast_name <-treedat$tip_name
  
  # now remove the capture data from the accessions 
  treedat$accession_num <- sapply(strsplit(treedat$tip_name, "_"), function(x) x[[1]])
  
  # and rename some of them
  treedat$accession_num[treedat$accession_num=="NC"] <- c("NC_015936", "NC_016769", 
                                                          "NC_023422", 
                                                          "NC_026314", 
                                                          "NC_027054",
                                                          "NC_027918", "NC_027919")
  
  # tree metadata ----------------------------------------------------------------
  
  # load the manual tree data that i created 
  # this will be used to create new tip labels for the tree
  dat <- read.csv(file = metadata, 
                  header = T, stringsAsFactors = F)
  
  # and omit any na 
  dat<-na.omit(dat)
  
  # change collection date to a more readable format 
  dat$collection_date <- as.Date(dat$collection_date, format = "%m/%d/%y")
  #dat$collection_date <- as.Date(dat$collection_date)
  
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
  
  # normalize the posterior values to a range between 0 and 1
  beast_tree@data$posterior <- as.numeric(beast_tree@data$posterior)
  min_posterior <- min(beast_tree@data$posterior, na.rm = TRUE)
  max_posterior <- max(beast_tree@data$posterior, na.rm = TRUE)
  beast_tree@data$posterior <- (beast_tree@data$posterior - min_posterior) / (max_posterior - min_posterior)
  
  # create a grayscale color mapping
  posfilz <- gray(seq(0, 1, by = 0.1))
  
  # date nodes of interest -------------------------------------------------------
  
  # define node of interest
  nodeOP287812 <- MRCA(beast_tree, which(beast_tree@phylo$tip.label == "OP287812_2018-07-27"), 
                       which(beast_tree@phylo$tip.label == "FJ890523_2008-07-15"))
  
  nodeKobuVs <- MRCA(beast_tree, which(beast_tree@phylo$tip.label == "OP287812_2018-07-27"), 
                     which(beast_tree@phylo$tip.label == "KJ641691_2012-12-15"))
  nodeKJ934637 <- MRCA(beast_tree, which(beast_tree@phylo$tip.label == "OP287812_2018-07-27"), 
                       which(beast_tree@phylo$tip.label == "KJ934637_2011-07-15"))
  
  nodeLuSeqs <- MRCA(beast_tree, which(beast_tree@phylo$tip.label == "MF947435_2014-06-11"), 
                     which(beast_tree@phylo$tip.label == "MN116647_2018-07-15"))
  
  
  # check we are identifying the right nodes 
  p1 <- p1 + geom_point(data=subset(p1$data, node == nodeOP287812), aes(x=x, y=y), color='blue', size=3)
  p1 <- p1 + geom_point(data=subset(p1$data, node == nodeKobuVs), aes(x=x, y=y), color='red', size=3)
  p1 <- p1 + geom_point(data=subset(p1$data, node == nodeKJ934637), aes(x=x, y=y), color='green', size=3)
  p1
  
  # find date that OP87812 branches off 
  op287812.date <- round(node.sub$nodetime[nodeOP287812],0) # works fine
  recent.age <- year(mrsd.dat) + yday(mrsd.dat)/365 # works fine 
  op287812.mean <- round(recent.age - tree.dat$height[52],0) # works after adjustment
  
  # upper and lower confidence intervals for this date 
  op287812.uci <- round(recent.age - tree.dat$height_0.95_HPD[52][[1]][1],0)
  op287812.lci <- round(recent.age - tree.dat$height_0.95_HPD[52][[1]][2],0)
  
  # format the date for the node
  op287812.date <- paste0("~", op287812.mean, "\n[", op287812.lci, "-", op287812.uci, "]")
  
  # create a new node label vector with NA values
  new.nodel.lab <- rep(NA, nrow(node.sub))
  
  # assign the formatted date to the specific node
  new.nodel.lab[nodeOP287812] <- op287812.date
  
  
  # output the calculated dates
  list(
    op287812_mean = op287812.mean,
    op287812_uci = op287812.uci,
    op287812_lci = op287812.lci,
    op287812_date = op287812.date
  )
  
  
  # find date that kobuvirus clade originates 
  KobuClade.date <- round(node.sub$nodetime[nodeKobuVs],0) # works fine
  recent.age <- year(mrsd.dat) + yday(mrsd.dat)/365 # works fine
  KobuClade.mean <- round(recent.age - tree.dat$height[48],0) # works fine after adjustment
  
  # upper and lower confidence intervals for this date 
  KobuClade.uci <- round(recent.age - tree.dat$height_0.95_HPD[48][[1]][1],0)
  KobuClade.lci <- round(recent.age - tree.dat$height_0.95_HPD[48][[1]][2],0)
  
  # format the date for the node
  KobuClade.date <- paste0("~", KobuClade.mean, "\n[", KobuClade.lci, "-", KobuClade.uci, "]")
  
  # assign the formatted date to the specific node
  new.nodel.lab[nodeKobuVs] <- KobuClade.date
  
  
  
  # lastly find date that avian kobuvirus branches off
  nodeKJ934637.date <- round(node.sub$nodetime[nodeKJ934637],0) # works fine
  recent.age <- year(mrsd.dat) + yday(mrsd.dat)/365 # works fine
  KJ934637.mean <- round(recent.age - tree.dat$height[51],0) # works fine after adjustment
  
  # upper and lower confidence intervals for this date 
  KJ934637.uci <- round(recent.age - tree.dat$height_0.95_HPD[51][[1]][1],0)
  KJ934637.lci <- round(recent.age - tree.dat$height_0.95_HPD[51][[1]][2],0)
  
  # format the date for the node
  KJ934637.date <- paste0("~", KJ934637.mean, "\n[", KJ934637.lci, "-", KJ934637.uci, "]")
  
  # assign the formatted date to the specific node
  new.nodel.lab[nodeKJ934637] <- KJ934637.date
  
  
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
  
  p <- ggtree(beast_tree, mrsd=mrsd.dat, size=.8) %<+% dat.sub +
    geom_tippoint(aes(color=title), size=4) + 
    geom_nodepoint(aes(fill = posterior), shape = 21, color = "black", size = 2, stroke = .1) +
    scale_fill_gradient(low = "white", high = "black", breaks = seq(0, 1, by = 0.2), labels = seq(0, 1, by = 0.2),
                        guide = guide_colorbar(direction = "vertical")) + 
    theme_tree2() +
    theme(axis.text.x = element_text(size = 14),
          axis.title.x = element_text(size=14),  
          legend.position = c(.6, .6), 
          plot.margin = unit(c(2, 22, 2, 3), "lines"),  
          legend.title = element_blank(), 
          legend.text = element_text(size=14)) +
    scale_x_continuous(breaks = c(-611, 1200, 1400, 1600, 1800, 2000), 
                       labels = c("611 BCE","1200","1400", "1600", "1800", "2000" 
                       )) +
    ggnewscale::new_scale_fill() +
    xlab("Divergence Times") +
    geom_tiplab(geom = "label", label.size = 0, alpha=.3, show.legend=F, size=4.5) + 
    coord_cartesian(clip = "off")
  
  return(p)
  
}

# assemble plot ----------------------------------------------------------------

get_beast_tree(kobuvirus_beastfile, kobuvirus_beastmetadata)

# now save plot and export to adobe 
ggsave(outpng,
       plot = p3, 
       units="mm",  
       width=150, 
       height=100, 
       #limitsize = F,
       scale=3)#, 
