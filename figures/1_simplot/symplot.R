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
print(colors)

# start out with aa plot 
homewd="/Users/flg9/Desktop/Developer/brook_lab/Madagascar-Bat-Kobuvirus/figures/1_simplot/orf/"
setwd(paste0(homewd))

# load in aa simplot data 
aa_sim <- read.csv(file = paste0(homewd, "AA_identity.csv"), header = T, stringsAsFactors = F)

# move to long 
# this will change orientation of csv file for easier visualization
id.plot <- melt(aa_sim, id.vars = c("pointer"), measure.vars = c("NC_001918",  "KJ934637"))

# since our variable name is our accession number, change to character
# then head to make sure structure is fine 
id.plot$variable <- as.character(id.plot$variable)
head(id.plot)

# let's then change the name of our variable to Accesssion
# this will make it easier for people to look up the sequences 
names(id.plot)[names(id.plot)=="variable"] <- "Accession"

#id.plot$strain[id.plot$strain=="Eidolon_helvum"] <- "E. helvum Kobuvirus"

#id.plot$strain <- factor(id.plot$strain, levels = c("E. helvum Kobuvirus",  "R. madagascariensis Kobuvirus"))

#quick check
head(id.plot)

# ignore this

#id.plot$value[id.plot$value<0] <- 0
#id.plot$Query[id.plot$Query=="Eidolon_helvum"] <- "Eidolon helvum"
#id.plot$Query[id.plot$Query=="Rousettus_madagascariensis"] <- "Rousettus madagascariensis"
#id.plot$Query <- factor(id.plot$Query, labels = c("Pteropus rufus", "Rousettus madagascariensis"))

# 
# genome.df <- data.frame(position = c(1, 4316,
#                                      4316, 7002,
#                                      7002, 8324,
#                                      8324, 8580,
#                                      8580, 8659,
#                                      8659, 8881,
#                                      8881, 9378,
#                                      9378, 9820), 
#                         gene = rep(c("ORF1a", "ORF1b", "S", "NS3", "E", "M", "N", "NS7"), each=2))

# define polyprotein coordinates 
# previous characterization of these genome positions will help with this 
# i am using the minimum and max positions out of the three sequences to define these positions 
genome.df <- data.frame(position = c(1, 185, # max is position 185 for OP287812, this is carried throughout gene definitions
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

# normalize similarity 
# id.plot$value <- id.plot$value/100
genome.df$Peptide <- factor(genome.df$Peptide, levels = unique(genome.df$Peptide))

# get max & avg similarity among strains

max(id.plot$value)

mean(id.plot$value) #sim among both strains to query (0.7435408) # draw mean line at this val

id.plot_kj <- id.plot %>% filter(Accession == "KJ934637")
mean(id.plot_kj$value) # sim to bird kov (0.7323649)

id.plot_nc <- id.plot %>% filter(Accession == "NC_001918")
mean(id.plot_nc$value) # sim to human kov (0.7547168)

mean(id.plot_nc$value) 
  
# give the sequences that are being plotted a color

colz= c("NC_001918" = "firebrick3", "KJ934637" = "forestgreen")

# plot 

p3 <- ggplot(id.plot) + 
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

p3
 
 p3 + p2


## next, let's do our full genome nucleotide similarity

homewd="/Users/flg9/Desktop/Developer/brook_lab/Madagascar-Bat-Kobuvirus/figures/1_simplot/nucleotide/"
setwd(paste0(homewd))

nt_sim <- read.csv(file = paste0(homewd, "NT_identity.csv"), header = T, stringsAsFactors = F)

# move to long
# this will change orientation of csv file for easier visualization

id.plot2 <- melt(nt_sim, id.vars = c("pointer"), measure.vars = c("KJ934637",  "NC_001918"))

# since our variable name is our accession number, change to character
# then head to make sure structure is fine 

id.plot2$variable <- as.character(id.plot2$variable)
head(id.plot2)


# let's then change the name of our variable to strain
# this is essential for visualization at the end

names(id.plot2)[names(id.plot2)=="variable"] <- "Accession"

#id.plot2$strain[id.plot2$strain=="Eidolon_helvum"] <- "E. helvum Kobuvirus"
#id.plot2$strain[id.plot2$strain=="Rousettus_madagascariensis"] <- "R. madagascariensis Kobuvirus"

#id.plot$strain <- factor(id.plot$strain, levels = c("HKU9", "GCCDC1", "GX2018.BtCoV92", "E. helvum bat coronavirus", "P. rufus Nobecovirus", "R. madagascariensis Nobecovirus"))
#id.plot$strain <- factor(id.plot$strain, levels = c("E. helvum Kobuvirus",  "R. madagascariensis Kobuvirus"))
#head(id.plot)


#and plot

#id.plot$value[id.plot$value<0] <- 0
#id.plot$Query[id.plot$Query=="Eidolon_helvum"] <- "Eidolon helvum"
#id.plot$Query[id.plot$Query=="Rousettus_madagascariensis"] <- "Rousettus madagascariensis"
#id.plot$Query <- factor(id.plot$Query, labels = c("Pteropus rufus", "Rousettus madagascariensis"))

# 
# genome.df <- data.frame(position = c(1, 4316,
#                                      4316, 7002,
#                                      7002, 8324,
#                                      8324, 8580,
#                                      8580, 8659,
#                                      8659, 8881,
#                                      8881, 9378,
#                                      9378, 9820), 
#                         gene = rep(c("ORF1a", "ORF1b", "S", "NS3", "E", "M", "N", "NS7"), each=2))

# define genome positions 
# these are the old genome positions including the UTRs
# opting to remove these for easier visualization in paper 
#genome.df2 <- data.frame(position = c(1, 679,
#                                     680, 1234,
#                                     1235, 2347,
#                                     2348, 3016,
#                                     3017, 3742,
#                                     3743, 4150,
#                                     4151, 4645, 
#                                     4646, 5650, 
#                                     5651, 5929, 
#                                     5930, 6007, 
#                                     6008, 6577,
#                                     6578, 7984,
#                                     7985, 8505), 
#                        Peptide = rep(c("5UTR", "L", "VP0", "VP3", "VP1", "2A", "2B", "2C", "3A", "3B", "3C", "3D", "3UTR"), each=2))

# new genome positions excluding the UTRs
genome.df2 <- data.frame(position = c(1, 555, # again use the max bp position to define gene positions # L
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



# normalize similarity 
# id.plot2$value <- id.plot2$value/100
genome.df2$Peptide <- factor(genome.df2$Peptide, levels = unique(genome.df2$Peptide))

# get max & avg similarity 

max(id.plot2$value) 

mean(id.plot2$value) #sim among both strains to query (0.5288817)

id.plot_kj <- id.plot2 %>% filter(Accession == "KJ934637")
mean(id.plot_kj$value) # sim to bird kov (0.5009596)

id.plot_nc <- id.plot2 %>% filter(Accession == "NC_001918")
mean(id.plot_nc$value) # sim to human kov (0.5568037)


#colz= c("HKU9"="firebrick3", "GCCDC1"="magenta", "GX2018.BtCoV92" = "purple", "E. helvum bat coronavirus" = "royalblue", "P. rufus Nobecovirus" = "goldenrod", "R. madagascariensis Nobecovirus" = "forestgreen")
#colz= c("HKU9"="firebrick3", "E. helvum bat coronavirus" = "royalblue", "P. rufus Nobecovirus" = "goldenrod", "R. madagascariensis Nobecovirus" = "forestgreen")

#colz= c("E. helvum Kobuvirus" = "firebrick3", "R. madagascariensis Kobuvirus" = "forestgreen")

#colz= c("NC_001918" = "firebrick3", "KJ934637" = "forestgreen")


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

p2

p2 + p3


## now let's plot coverage 

homewd="/Users/flg9/Desktop/Developer/brook_lab/Madagascar-Bat-Kobuvirus/figures/1_simplot/coverage/"
setwd(paste0(homewd))

datcovg <- read.csv(file = paste0(homewd, "Coverage copy.csv"), header = T, stringsAsFactors = F)
position <- read.csv(file = paste0(homewd, "position copy.csv"), header = T, stringsAsFactors = F)

# get mean depth
mean(datcovg$Coverage)

#datcovg$Coverage <- as.character(datcovg$Coverage)
names(datcovg)[names(datcovg)=="Coverage"] <- "Coverage"

# this will give raw coverage
datcovg$Coverage <- datcovg$Coverage/100

# this will give rpm
datcovg$Coverage <- datcovg$Coverage/31.282
#genome.df2$Peptide <- factor(genome.df2$Peptide, levels = unique(genome.df2$Peptide))
colnames(genome.df2) <- c("Position", "Peptide")

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


p1 <- ggplot(datcovg2) + geom_area(aes(x=Position, y=Coverage, fill = Peptide), size=1, show.legend = F) +
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

p1
                    
p1 + p2 + p3

  #lacroix_palette("PeachPear", n = 13, type = "continuous") 

#Darjeeling2
#Moonrise2

ggplot(datcovg2) + geom_area(aes(x=Position, y=Coverage, fill = Peptide), size=1) +
  geom_ribbon(data=genome.df3, aes(x = Position, ymin=0, ymax=0,fill = Peptide), color="black") + 
  facet_grid() + theme_bw() + xlab("Genome position") + ylab("Coverage (rpm)") + 
  theme(panel.grid = element_blank(), strip.text = element_text(face="italic", size=14), 
        strip.background = element_rect(fill="white"), 
        legend.position = "bottom", legend.direction = "horizontal",legend.box = "horizontal",
        legend.text = element_text(face="italic", size = 12),
        axis.text = element_text(size=12), axis.title = element_text(size=14)) + 
  ggtitle("OP287812 Coverage") + scale_x_continuous(breaks=c(0,2000/1,4000/1,6000/1, 8000/1), 
                                                    labels = c(0,2000, 4000, 6000, 8000)) + 
  scale_fill_manual(values = wes_palette("Darjeeling2", 13, type = "continuous"))
 
ggplot(datcovg2) + geom_area(aes(x=Position, y=Coverage, fill = Peptide), size=1) +
  geom_ribbon(data=genome.df3, aes(x = Position, ymin=0, ymax=0,fill = Peptide), color="black") + 
  facet_grid() + theme_bw() + xlab("Genome position") + ylab("Coverage (rpm)") + 
  theme(panel.grid = element_blank(), strip.text = element_text(face="italic", size=14), 
        strip.background = element_rect(fill="white"), 
        legend.position = "bottom", legend.direction = "horizontal",legend.box = "horizontal",
        legend.text = element_text(face="italic", size = 12),
        axis.text = element_text(size=12), axis.title = element_text(size=14)) + 
  ggtitle("OP287812 Coverage") + scale_x_continuous(breaks=c(0,2000/1,4000/1,6000/1, 8000/1), 
                                                    labels = c(0,2000, 4000, 6000, 8000)) + 
  scale_fill_manual(values = wes_palette("Moonrise2", 13, type = "continuous"))

ggplot(datcovg2) + geom_area(aes(x=Position, y=Coverage, fill = Peptide), size=1) +
  geom_ribbon(data=genome.df3, aes(x = Position, ymin=0, ymax=0,fill = Peptide), color="black") + 
  facet_grid() + theme_bw() + xlab("Genome position") + ylab("Coverage (rpm)") + 
  theme(panel.grid = element_blank(), strip.text = element_text(face="italic", size=14), 
        strip.background = element_rect(fill="white"), 
        legend.position = "bottom", legend.direction = "horizontal",legend.box = "horizontal",
        legend.text = element_text(face="italic", size = 12),
        axis.text = element_text(size=12), axis.title = element_text(size=14)) + 
  ggtitle("OP287812 Coverage") + scale_x_continuous(breaks=c(0,2000/1,4000/1,6000/1, 8000/1), 
                                                    labels = c(0,2000, 4000, 6000, 8000)) + 
  scale_fill_manual(values = wes_palette("GrandBudapest1", 13, type = "continuous"))

## plot all together
fig5 <- cowplot::plot_grid(p1,p2,p3, labels = "AUTO", 
                           label_size = 16, nrow = 3, ncol = 1, 
                           label_colour = "#D32D41")

fig5 <- cowplot::plot_grid(p1, p2, p3, 
                          nrow=3, ncol=1)

fig5

# and save
homewd="/Users/flg9/Desktop/Developer/brook_lab/Madagascar-Bat-Kobuvirus/figures/1_simplot/"
setwd(paste0(homewd))

ggsave(file = paste0(homewd, "Fig1_SimPlot.png"),
       plot=fig5,
       units="mm",  
       width=95, 
       height=70, 
       #limitsize = F,
       scale=4)#, 


