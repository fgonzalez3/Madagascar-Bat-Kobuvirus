rm(list=ls())

library(plyr)
library(dplyr)
library(ggplot2)
library(reshape2)
library(gggenes)
library(wesanderson)
library(LaCroixColoR)
library(devtools)
library(cowplot)

colors <- wes_palette("Darjeeling2")
#print(colors)

aa_sim_data <- "./data/orf/AA_identity.csv"
nt_sim_data <- "./data/nucleotide/NT_identity.csv"
coverage_data <- "./data/coverage/Coverage.csv"
position_data <- "./data/coverage/position.csv"

outpng <- "./Fig1.png"

# aa similarity plot -----------------------------------------------------------

get_aa_similarity <- function(aa_dat) {
  
  # load in data 
  aa_sim <- read.csv(file=paste0(aa_dat), header = T, stringsAsFactors = F)
  
  # change orientation of the csv file for easier parsing
  id.plot <- melt(aa_sim, id.vars = c("pointer"), measure.vars = c("NC_001918",  "KJ934637"))
  
  # since our variable name is our accession number, change to character
  # then head to make sure structure is fine 
  id.plot$variable <- as.character(id.plot$variable)
  print(head(id.plot))
  
  # let's then change the name of our variable to Accesssion
  # this will make it easier for people to look up the sequences 
  names(id.plot)[names(id.plot)=="variable"] <- "Accession"
  
  #quick check
  print(head(id.plot))
  
  # define polyprotein coordinates 
  # previous characterization of these genome positions will help with this 
  genome.df <- data.frame(position = c(1, 185, 
                                       186, 556, 
                                       557, 779,
                                       780, 1021,
                                       1022, 1157,
                                       1158, 1322,
                                       1323, 1657, 
                                       1658, 1750, 
                                       1751, 1776, 
                                       1777, 1966, 
                                       1967, 2430), 
                          Peptide = rep(c("L", "VP0", "VP3", "VP1", "2A", "2B", "2C", "3A", "3B", "3C", "3D"), each=2))
  
  genome.df$Peptide <- factor(genome.df$Peptide, levels = unique(genome.df$Peptide))
  
  # get max & average similarity among strains
  max(id.plot$value)
  
  mean(id.plot$value) #sim among both strains to query (0.7435408) # draw mean line at this val
  
  id.plot_kj <- id.plot %>% filter(Accession == "KJ934637")
  mean(id.plot_kj$value) # sim to bird kov (0.7323649)
  
  id.plot_nc <- id.plot %>% filter(Accession == "NC_001918")
  mean(id.plot_nc$value) # sim to human kov (0.7547168)
  
  mean(id.plot_nc$value) 
  
  # plot 
  p1 <- ggplot(id.plot) + 
    geom_line(aes(x=pointer, y=value, color=Accession), size=1) + 
    geom_ribbon(data=genome.df, aes(x=position, ymin=-5, ymax=-0.5, fill=Peptide), color="black") + 
    geom_hline(yintercept = mean(id.plot$value), linetype ="dashed", color="black") + 
    facet_grid() + 
    theme_bw() + 
    xlab("Genome position") + 
    ylab("Amino acid similarity (%)") +
    theme(
      panel.grid = element_blank(), 
      strip.text = element_text(face="italic", size=14), 
      strip.background = element_rect(fill="white"), 
      legend.position = "bottom",  # Set legend position to bottom
      legend.direction = "horizontal",  # Set legend direction to horizontal
      legend.text = element_text(face="italic", size = 12),
      axis.text = element_text(size=12), 
      axis.title = element_text(size=14, hjust=0.47)
    ) + 
    ggtitle("Amino acid similarity to OP287812") + 
    scale_x_continuous(breaks=c(0, 800, 1600, 2400), labels = c(0, 800, 1600, 2400)) + 
    scale_fill_manual(values = wes_palette("Darjeeling2", 11, type = "continuous")) +
    guides(fill = FALSE)
  
  return(p1)
  
}

# nt similarity plot -----------------------------------------------------------

get_nt_similarity <- function(nt_dat) {
  
  # read in data
  nt_sim <- read.csv(file=paste0(nt_dat), header = T, stringsAsFactors = F)
  
  # change orientation of the csv file for easier parsing
  id.plot2 <- melt(nt_sim, id.vars = c("pointer"), measure.vars = c("KJ934637",  "NC_001918"))
  
  # since our variable name is our accession number, change to character
  # then head to make sure structure is fine 
  id.plot2$variable <- as.character(id.plot2$variable)
  print(head(id.plot2))
  
  # let's then change the name of our variable to strain
  # this is essential for visualization at the end
  names(id.plot2)[names(id.plot2)=="variable"] <- "Accession"
  
  # define nt genome positions excluding the UTRs
  genome.df2 <- data.frame(position = c(1, 555, # L
                                        556, 1668, #VP0
                                        1669, 2337, #VP3
                                        2338, 3063, #VP1
                                        3064, 3471, #2A
                                        3472, 3966, #2B
                                        3967, 4971, #2C
                                        4972, 5250, #3A
                                        5251, 5328, #3B
                                        5329, 5898, #3C
                                        5899, 7401), #3D
                           Peptide = rep(c("L", "VP0", "VP3", "VP1", "2A", "2B", "2C", "3A", "3B", "3C", "3D"), each=2))
  
  genome.df2$Peptide <- factor(genome.df2$Peptide, levels = unique(genome.df2$Peptide))
  
  # get max & avg similarity 
  max(id.plot2$value) 
  
  mean(id.plot2$value) #sim among both strains to query (0.5288817)
  
  id.plot_kj <- id.plot2 %>% filter(Accession == "KJ934637")
  mean(id.plot_kj$value) # sim to bird kov (0.5009596)
  
  id.plot_nc <- id.plot2 %>% filter(Accession == "NC_001918")
  mean(id.plot_nc$value) # sim to human kov (0.5568037)
  
  # and plot w legend 
  p2 <- ggplot(id.plot2) + geom_line(aes(x=pointer, y=value, color=Accession), size=1) + 
    geom_ribbon(data = genome.df2, aes(x=position, ymin=-5, ymax=-0.5, fill=Peptide), color="black") +
    geom_hline(yintercept = mean(id.plot2$value), linetype ="dashed", color="black") + 
    scale_y_continuous(limits = c(-5, 100)) + 
    facet_grid() + theme_bw() + xlab("Genome position") + ylab("Nucleotide similarity (%)") +
    theme(panel.grid = element_blank(), strip.text = element_text(face="italic", size=14), 
          strip.background = element_rect(fill="white"), 
          #legend.position = "bottom", legend.direction = "horizontal", #legend.box = "vertical",
          legend.position = "none",
          legend.text = element_text(face="italic", size = 12),
          axis.text = element_text(size=12), axis.title = element_text(size=14, hjust=0.47)) + 
    ggtitle("Nucleotide similarity to OP287812") +
    scale_fill_manual(values = wes_palette("Darjeeling2", 11, type = "continuous"))
  
  return(p2)
  
}

# coverage plot ----------------------------------------------------------------

get_coverage_plot <- function(coverage, position) {
  
  # read in data
  datcovg <- read.csv(file = paste0(coverage_data), header = T, stringsAsFactors = F)
  position <- read.csv(file = paste0(position_data), header = T, stringsAsFactors = F)
  
  # get mean depth
  mean(datcovg$Coverage)
  
  #datcovg$Coverage <- as.character(datcovg$Coverage)
  names(datcovg)[names(datcovg)=="Coverage"] <- "Coverage"
  
  # this will give raw coverage
  #datcovg$Coverage <- datcovg$Coverage/100
  
  # this will give rpm
  datcovg$Coverage <- datcovg$Coverage/31.282
  
  #colnames(genome.df3) <- c("Position", "Peptide")
  
  datcovg2 <- right_join(position, datcovg, by = "Position")
  
  #define genome positions
  genome.df3 <- data.frame(Position = c(1, 555,
                                        556, 1668,
                                        1669, 2337,
                                        2338, 3063,
                                        3064, 3471,
                                        3472, 3966,
                                        3967, 4971, 
                                        4972, 5250, 
                                        5251, 5328, 
                                        5329, 5898, 
                                        5899, 7305),
                           Peptide = rep(c("L", "VP0", "VP3", "VP1", "2A", "2B", "2C", "3A", "3B", "3C", "3D"), each=2))
  
  
  p3 <- ggplot(datcovg2) + geom_area(aes(x=Position, y=Coverage, fill = Peptide), size=1, show.legend = F) +
    geom_ribbon(data=genome.df3, aes(x = Position, ymin=0, ymax=0,fill = Peptide), color="black", show.legend = F) + 
    # geom_hline(yintercept = mean(datcovg2$Coverage), linetype = "dashed", color="black") + 
    facet_grid() + theme_bw() + xlab("Genome position") + ylab("Read Depth (rpm)") + 
    theme(panel.grid = element_blank(), strip.text = element_text(face="italic", size=14), 
          strip.background = element_rect(fill="white"), 
          legend.position = "bottom", legend.direction = "horizontal",legend.box = "horizontal",
          legend.text = element_text(face="italic", size = 12),
          axis.text = element_text(size=12), axis.title = element_text(size=14, hjust=0.47)) + 
    ggtitle("OP287812 Read Depth") + scale_x_continuous(breaks=c(0,2000/1,4000/1,6000/1, 8000/1), 
                                                        labels = c(0,2000, 4000, 6000, 8000)) + 
    scale_fill_manual(breaks = c("L", "VP0", "VP3", "VP1", "2A", "2B", "2C", "3A", "3B", "3C", "3D"), 
                      values = c(wes_palette("Darjeeling2", 11, type = "continuous")))
  
  return(p3)
  
  
}

# final plot -------------------------------------------------------------------

Kobu_AASimPlot <- get_aa_similarity(aa_sim_data)
Kobu_NTSimPlot <- get_nt_similarity(nt_sim_data)
Kobu_CoveragePlot <- get_coverage_plot(coverage_data,position_data)
  
Fig1 <-  cowplot::plot_grid(Kobu_AASimPlot, Kobu_NTSimPlot, Kobu_CoveragePlot,
                            nrow=3, ncol=1)
Fig1

ggsave(outpng,
       plot=fig5,
       units="mm",  
       width=95, 
       height=70, 
       #limitsize = F,
       scale=4)#, 

