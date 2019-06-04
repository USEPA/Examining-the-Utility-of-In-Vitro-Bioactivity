#-----------------------------------------------------------------------------------#
# APCRA case study: Generate all Main Text and Supplemental Figures
# 
# Katie Paul Friedman, paul-friedman.katie@epa.gov
# Updated February 2019
#-----------------------------------------------------------------------------------#

#-----------------------------------------------------------------------------------#
# loading libraries and sources
#-----------------------------------------------------------------------------------#

rm(list = ls())

library(cowplot)
library(data.table)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(ggthemes)
library(grid)
library(gridExtra)
library(httk)
library(plyr)
library(RColorBrewer)
library(reshape2)
library(RMySQL)
library(scales)
library(stringr)
library(tcpl)
library(xlsx)

#setwd("/APCRA_final_retrospective/R") # set working directory and all others will be relative
load('../pod_ratio/pod_ratio_master_final.RData')
load('../in_vitro/source_invitro_data.RData')

#-----------------------------------------------------------------------------------#
# Figure 2 Heatmap of Substances Available, by AcTOR Use Category
#-----------------------------------------------------------------------------------#
# specifically for Figure 2 heatmap
source("https://bioconductor.org/biocLite.R")
biocLite("ComplexHeatmap")
library(ComplexHeatmap)
library(circlize)

# load source files

load('../expocast/source/UseData-012015.RData') # uses from AcTOR
casrn.list <- pod.ratio.master[,CASRN]

# using AcTOR use categories

colnames(NewUseTable)
case.use.actor <- data.table(NewUseTable)
case.use.actor <- case.use.actor[CASRN %in% casrn.list]

str(case.use.actor)
case.use.actor<-as.data.frame(case.use.actor)
case.use.actor[,'ANTIMICROBIAL'] <- as.numeric(as.character(case.use.actor[,'ANTIMICROBIAL']))
case.use.actor[,'COLORANT'] <- as.numeric(as.character(case.use.actor[,'COLORANT']))
case.use.actor[,'FERTILIZER'] <- as.numeric(as.character(case.use.actor[,'FERTILIZER']))
case.use.actor[,'FOOD.ADDITIVE'] <- as.numeric(as.character(case.use.actor[,'FOOD.ADDITIVE']))
case.use.actor[,'FRAGRANCE'] <- as.numeric(as.character(case.use.actor[,'FRAGRANCE']))
case.use.actor[,'HERBICIDE'] <- as.numeric(as.character(case.use.actor[,'HERBICIDE']))
case.use.actor[,'PERSONAL.CARE.PRODUCT'] <- as.numeric(as.character(case.use.actor[,'PERSONAL.CARE.PRODUCT']))
case.use.actor[,'PETROCHEMICAL'] <- as.numeric(as.character(case.use.actor[,'PETROCHEMICAL']))
case.use.actor[,'PESTICIDE.INERT'] <- as.numeric(as.character(case.use.actor[,'PESTICIDE.INERT']))
case.use.actor[,'FLAME.RETARDANT'] <- as.numeric(as.character(case.use.actor[,'FLAME.RETARDANT']))

case.use.actor <- as.data.table(case.use.actor)
colnames(case.use.actor)
colSums(case.use.actor[,c(4:19)])
rowSums(case.use.actor[,c(4:19)]) # more than one use category can be assigned
case.use.actor[,pesticide := sum(PESTICIDE.ACTIVE.NO.CONSUMER,
                                 PESTICIDE.ACTIVE.AND.CONSUMER,
                                 HERBICIDE,
                                 ANTIMICROBIAL), by=list(CASRN)]
case.use.pest <- case.use.actor[pesticide > 0]
length(unique(case.use.pest$CASRN)) #314 chemicals
(314/448)*100 #70% pesticide uses

# Use ComplexHeatMap from BioConductor to make heatmap of chemical use type from AcTOR

colors <- brewer.pal(6, "OrRd")
colnames(case.use.actor)

heat_mat <- as.matrix(case.use.actor[,c(4:19)])
heat_mat[is.na(heat_mat)] <- 0
rownames(heat_mat) <- case.use.actor[, COMMON.NAME]

colors = structure(c("#2166ac", "#f7f7f7"), names = c("1", "0"))

heatmap_fig <- Heatmap(matrix = heat_mat, 
                       cluster_columns = TRUE, 
                       col= colors,
                       #col = colorRamp2(breaks = c(1, 0),
                       #                  colors = c("#2166ac", "#f7f7f7")),
                       #col = list(type=c("1" = "#2166ac","0"= "#f7f7f7")),
                       show_row_names = FALSE, 
                       #row_dend_width = unit(3, "cm"),
                       column_names_max_height = unit(8, "cm"),
                       column_names_gp = gpar(fontsize = 12),
                       clustering_method_columns = "ward.D",
                       clustering_method_rows = "ward.D",
                       clustering_distance_rows = "euclidean",
                       clustering_distance_columns = "euclidean",
                       show_row_dend = FALSE,
                       show_column_dend = FALSE,
                       heatmap_legend_param = list(title = ""))

file.dir <- paste("Figures/", sep="")
file.name <- paste("/Figure_2_chemical_use_type_", Sys.Date(), ".png", sep="")
file.path <- paste(file.dir, file.name, sep="")
dir.create(path=file.dir, showWarnings = FALSE, recursive = TRUE)
png(file.path, width = 4000, height = 4000, res=480)
draw(heatmap_fig)
dev.off()

#-----------------------------------------------------------------------------------#
# Figure 3 Comparison of ExpoCast, POD-NAM, and POD traditional
#-----------------------------------------------------------------------------------#

# Make long data for plotting and visualization

colnames(pod.ratio.master)

pod.ratio.master.long <- melt.data.table(pod.ratio.master, 
                                         id.vars = c('DTXSID',
                                                     'CASRN',
                                                     'Name',
                                                     'pod.ratio.95',
                                                     'pod.nam.50',
                                                     'pod.nam.95'),
                                         measure.vars = c('log.expo.95',
                                                          'pod.nam.95',
                                                          'log.max.AED',
                                                          'log.p5.POD',
                                                          'ber.50',
                                                          'ber.95',
                                                          'ber.50.pod.50',
                                                          'ber.95.pod.50'),
                                         variable.name = c('Dose'))

predictivity_sort <- sort(unique(pod.ratio.master.long[,pod.ratio.95]))
pod.ratio.master.long <- pod.ratio.master.long[order(-pod.ratio.95)]
pod.ratio.master <- pod.ratio.master[order(-pod.ratio.95)]

pod.ratio.plot <- pod.ratio.master.long[Dose %in% c('log.expo.95',
                                                    'pod.nam.95',
                                                    'log.max.AED',
                                                    'log.p5.POD')]

# Set up annotation boxes and sorting

predictivity_0_cutoff <- which(predictivity_sort>0)[1]
predictivity_2_cutoff<-which(predictivity_sort> 2.5)[1]

manu_1_arrow_x <- 4.5

size_box <- 24
top_first_box <- 448
bottom_first_box <- top_first_box - size_box
top_second_box <- predictivity_2_cutoff
bottom_second_box <- top_second_box - size_box
top_third_box <- predictivity_0_cutoff
bottom_third_box <- top_third_box - 48

# Image 1 ----------  

manuscript_1  <- ggplot(data = pod.ratio.plot,
                        aes(x = value,
                            y = reorder(Name,pod.ratio.95), 
                            group = factor(Dose),
                            shape = factor(Dose), 
                            color = factor(Dose)))+
  geom_point() +
  geom_segment(aes(x=pod.nam.50,
                   xend=pod.nam.95,
                   y=reorder(Name,pod.ratio.95),
                   yend=reorder(Name,pod.ratio.95)),
               color = "#009E73",
               show.legend=FALSE)+
  xlab("log10 mg/kg-bw/day") +
  ylab('Chemical')+
  annotate("rect", 
           xmin = -9, 
           xmax = 4.5, 
           ymin = bottom_first_box, 
           ymax = top_first_box, 
           fill = "red", 
           alpha = 0.1, 
           size = 1, 
           color = "red") +
  annotate("rect", 
           xmin = -9, 
           xmax = 4.5, 
           ymin = bottom_second_box,  
           ymax = top_second_box, 
           fill = "red", 
           alpha = 0.1, 
           size = 1, 
           color = "red") +
  annotate("rect", 
           xmin = -9, 
           xmax = 4.5, 
           ymin = bottom_third_box,   
           ymax = top_third_box,  
           fill = "red", 
           alpha = 0.1, 
           size = 1, 
           color = "red") +
  geom_segment(aes(x = manu_1_arrow_x, 
                   y = top_first_box, 
                   xend = manu_1_arrow_x, 
                   yend = top_third_box),
               color = "black",
               arrow = arrow(type = "closed", 
                             angle = 10, 
                             ends="both")) +
  geom_segment(aes(x = manu_1_arrow_x, 
                   y = top_third_box, 
                   xend = manu_1_arrow_x, 
                   yend = bottom_third_box),
               color = "black",
               arrow = arrow(type = "closed", 
                             angle = 10, 
                             ends="both")) +
  annotate(geom="text", 
           x = 5, 
           y = mean(c(top_first_box, top_third_box)), 
           label = "POD ratio > 0", 
           color = "black",
           size = 3, 
           fontface="bold", 
           angle=90, 
           vjust=0.5) +
  annotate(geom="text", 
           x = 5.25, 
           y = mean(c(top_third_box, bottom_third_box)), 
           #y = mean(c(top_third_box)),
           label = "POD ratio \n < 0", 
           color = "black",
           size = 3, 
           fontface="bold", 
           angle=90, 
           vjust=0.5) +
  theme_bw() +
  theme(panel.border = element_blank(), 
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text = element_text(size=10),
        axis.title = element_text(size=12, face='bold'))+
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank())+
  theme(legend.position="bottom", legend.title=element_blank())+
  scale_x_continuous(breaks=seq(-8,5,1)) +
  coord_cartesian(xlim = c(-8, 5)) +
  scale_colour_manual(breaks=c('log.expo.95',
                               'pod.nam.95',
                               'log.max.AED',
                               'log.p5.POD'),
                      values=c("#999999",
                               "#009E73",
                               "#000000",
                               "#56B4E9"),
                      labels=c("ExpoCast","POD-NAM", "max AED", "POD-traditional")) +
  scale_shape_manual(breaks=c('log.expo.95',
                              'pod.nam.95',
                              'log.max.AED',
                              'log.p5.POD'),
                     values=c(19,19,17,15),
                     labels=c("ExpoCast","POD-NAM", "max AED", "POD-traditional"))


# Image 2a ----------   

casrn_a <- pod.ratio.master[c(1:24), CASRN]
pod.ratio.plot.a <- pod.ratio.plot[CASRN %in% casrn_a]

manuscript_2a <- ggplot(data = pod.ratio.plot.a, 
                        aes(x = value, 
                            y = reorder(strtrim(Name, 30), pod.ratio.95),
                            group = factor(Dose), 
                            color = factor(Dose), 
                            shape = factor(Dose))) + 
  coord_cartesian(xlim = c(-8, 5))+
  geom_point() +
  geom_segment(aes(x=pod.nam.50,
                   xend=pod.nam.95,
                   y=reorder(Name,pod.ratio.95),
                   yend=reorder(Name,pod.ratio.95)),
               color = "#009E73",
               show.legend=FALSE)+
  xlab("log10 mg/kg-bw/day") +
  theme_bw() +
  theme(legend.position="none") +
  theme(axis.title.y = element_blank()) +
  theme(axis.text.y = element_text(family = "sans", 
                                   face = "bold"))+
  scale_x_continuous(breaks=seq(-8,5,1)) +
  scale_colour_manual(breaks=c('log.expo.95',
                               'pod.nam.95',
                               'log.max.AED',
                               'log.p5.POD'),
                      values=c("#999999",
                                "#009E73",
                               "#000000",
                               "#56B4E9"),
                      labels=c("ExpoCast","POD-NAM", "max AED", "POD-traditional")) +
  scale_shape_manual(breaks=c('log.expo.95',
                              'pod.nam.95',
                              'log.max.AED',
                              'log.p5.POD'),
                     values=c(19,19,17,15),
                     labels=c("ExpoCast","POD-NAM", "max AED", "POD-traditional"))+ 
  annotate(geom="text", 
           x = 4.25, 
           y = 24, 
           label = "A", 
           color = "black", 
           hjust = 0, 
           vjust = 1, 
           size = 10)

# Image 2b ----------   

casrn_b <- pod.ratio.master[c(263:287), CASRN]
pod.ratio.plot.b <- pod.ratio.plot[CASRN %in% casrn_b]

manuscript_2b <- ggplot(data = pod.ratio.plot.b, 
                        aes(x = value, 
                            y = reorder(strtrim(Name,30), pod.ratio.95),
                            group = factor(Dose), 
                            color = factor(Dose), 
                            shape = factor(Dose))) + 
  coord_cartesian(xlim = c(-8,5))+
  geom_point() +
  geom_segment(aes(x=pod.nam.50,
                   xend=pod.nam.95,
                   y=reorder(Name,pod.ratio.95),
                   yend=reorder(Name,pod.ratio.95)),
               color = "#009E73",
               show.legend=FALSE)+
  xlab("log10 mg/kg-bw/day") +
  theme_bw() +
  theme(legend.position="none") +
  theme(axis.title.y = element_blank()) +
  theme(axis.text.y = element_text(family = "sans", 
                                   face = "bold"))+
  scale_x_continuous(breaks=seq(-8,5,1)) +
  scale_colour_manual(breaks=c('log.expo.95',
                               'pod.nam.95',
                               'log.max.AED',
                               'log.p5.POD'),
                      values=c("#999999",
                              "#009E73",
                               "#000000",
                               "#56B4E9"),
                      labels=c("ExpoCast","POD-NAM", "max AED", "POD-traditional")) +
  scale_shape_manual(breaks=c('log.expo.95',
                              'pod.nam.95',
                              'log.max.AED',
                              'log.p5.POD'),
                     values=c(19,19,17,15),
                     labels=c("ExpoCast","POD-NAM", "max AED", "POD-traditional"))+ 
  annotate(geom="text", 
           x = 4.25, 
           y = 24, 
           label = "B", 
           color = "black", 
           hjust = 0, 
           vjust = 1, 
           size = 10)

# Image 2c ----------   

casrn_c <- pod.ratio.master[c(401:448), CASRN]
pod.ratio.plot.c <- pod.ratio.plot[CASRN %in% casrn_c]

manuscript_2c <- ggplot(data = pod.ratio.plot.c, 
                        aes(x = value, 
                            y = reorder(strtrim(Name,30), pod.ratio.95),
                            group = factor(Dose), 
                            color = factor(Dose), 
                            shape = factor(Dose))) + 
  coord_cartesian(xlim = c(-8,5))+
  geom_point() +
  geom_segment(aes(x=pod.nam.50,
                   xend=pod.nam.95,
                   y=reorder(Name,pod.ratio.95),
                   yend=reorder(Name,pod.ratio.95)),
               color = "#009E73",
               show.legend=FALSE)+
  xlab("log10 mg/kg-bw/day") +
  theme_bw() +
  theme(legend.position="none") +
  theme(axis.title.y = element_blank()) +
  theme(axis.text.y = element_text(family = "sans", 
                                   face = "bold"))+
  scale_x_continuous(breaks=seq(-8,5,1)) +
  scale_colour_manual(breaks=c('log.expo.95',
                               'pod.nam.95',
                               'log.max.AED',
                               'log.p5.POD'),
                      values=c("#999999",
                               "#009E73",
                               "#000000",
                               "#56B4E9"),
                      labels=c("ExpoCast","POD-NAM", "max AED", "POD-traditional")) +
  scale_shape_manual(breaks=c('log.expo.95',
                              'pod.nam.95',
                              'log.max.AED',
                              'log.p5.POD'),
                     values=c(19,19,17,15),
                     labels=c("ExpoCast","POD-NAM", "max AED", "POD-traditional"))+ 
  annotate(geom="text", 
           x = 4.25, 
           y = 48, 
           label = "C", 
           color = "black", 
           hjust = 0, 
           vjust = 1, 
           size = 10)


layout <- matrix(c(1,1,2,2,2,
                   1,1,3,3,3,
                   1,1,4,4,4,
                   1,1,4,4,4), nrow=4, byrow=TRUE)

# multiplot was obtained from
# http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_%28ggplot2%29/
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

file.dir <- paste("Figures/", sep="")
file.name <- paste("/Figure_3_POD_ratio_rank_range_7Feb19", Sys.Date(), ".png", sep="")
file.path <- paste(file.dir, file.name, sep="")
dir.create(path=file.dir, showWarnings = FALSE, recursive = TRUE)
png(file.path, width = 6000, height = 7000, res=480)
multiplot(manuscript_1, manuscript_2a, manuscript_2b, manuscript_2c, layout = layout)
dev.off()

#-----------------------------------------------------------------------------------#
# Figure 4 Bioactivity:Exposure Ratio (BER) cumulative frequency
#-----------------------------------------------------------------------------------#

# create inset table of chems with BER < 0

colnames(pod.ratio.master)
ber.table <- pod.ratio.master[ber.95 < 0]
colnames(ber.table)
ber.table <- ber.table[,c('Name',
                          'ber.95',
                          'pod.nam.95',
                          'log.expo.95')] # need to round these values
ber.table <- ber.table[order(ber.95)]
ber.table <- ber.table[,Name,round(ber.table[,2:4], 2)]
setnames(ber.table, c("Name",
                      "ber.95",
                      "pod.nam.95",
                      "log.expo.95"), c("Chemical Name",
                                   "BER, 95%-ile",
                                   "log10(POD-NAM, 95th-quantile)",
                                   "log10(ExpoCast 95%-ile)"))

cols <- matrix(c("#000000"),nrow(ber.table), ncol(ber.table))
tt <- ttheme_minimal(
  core=list(fg_params = list(col = cols, cex=1.0),
            bg_params = list(col=NA)),
  rowhead=list(bg_params = list(col=NA), fg_params = list(fontface = "bold",
                                                          cex=1.0)),
  colhead=list(bg_params = list(col=NA), fg_params = list(cex=1.0)))
t.ber <- tableGrob(ber.table, theme = tt)


# cumulative distribution plot of BERs and list top 50 chemicals in a table grob

ber.plot <- ggplot(pod.ratio.master.long[Dose %in% c('ber.50', 
                                                     'ber.95', 
                                                     'ber.50.pod.50', 
                                                     'ber.95.pod.50')], aes(value, color=Dose))+
  stat_ecdf(geom='step')+
  scale_y_continuous(trans = 'log10',
                     breaks= c(0.01, 0.1,0.2,0.3,0.4,0.5,0.75,1))+
  ylab("Cumulative Frequency") +
  xlab('Bioactivity-to-Exposure Ratio (BER)')+
  theme_bw() +
  theme(
    axis.line = element_line(colour = "black"),
    axis.text = element_text(size=10),
    axis.title = element_text(size=12, face='bold'))+
  theme(axis.text.y = element_text(family = "sans", face = "bold", size=12))+
  theme(legend.position="right", legend.title=element_blank())+
  scale_x_continuous(breaks=seq(-5,10,1)) +
  coord_cartesian(xlim = c(-5, 10)) +
  scale_colour_manual(breaks=c('ber.95.pod.50', 'ber.95', 'ber.50.pod.50', 'ber.50'),
                      values=c("#E69F00","#000000", "#0072B2", "#D55E00"),
                      labels=c("BER,95 50th%-ile","BER,95 95th%-ile", 'BER,50 50th%-ile', 'BER,50 95th%-ile'))+
  geom_vline(xintercept=0, lty='dashed', color='red')+
  geom_hline(yintercept=0.0245, lty='dashed', color='red')

fig.ber<-plot_grid(ber.plot, t.ber,
                   nrow=2,
                   rel_heights=c(1,1),
                   labels=c("A", "B"),
                   label_size=18,
                   label_fontface = 'bold')


file.dir <- paste("Figures/", sep="")
file.name <- paste("/Figure_4_BER_dist_list_7Feb19", Sys.Date(), ".png", sep="")
file.path <- paste(file.dir, file.name, sep="")
dir.create(path=file.dir, showWarnings = FALSE, recursive = TRUE)
png(file.path, width = 5000, height = 5000, res=480)
fig.ber
dev.off()

#-----------------------------------------------------------------------------------#
# Figure 5 ExpoCast and In vitro bioactivity for substances with BER <0
#-----------------------------------------------------------------------------------#

ber.cas <- pod.ratio.master[ber.95 < 0, CASRN]

invitro.dist$ber.95 <- pod.ratio.master$ber.95[match(pod.ratio.master$CASRN,
                                                          invitro.dist$chemcas)]
invitro.dist[,logAC50p5 := log10(AC50p5)]
invitro.dist <- invitro.dist[order(ber.95)]

toxcast.ber <- toxcast.master[casn %in% ber.cas]
toxcast.ber$ber.95 <- pod.ratio.master$ber.95[match(toxcast.ber$casn,
                                                    pod.ratio.master$CASRN)]
# make plot of distribution of hits for ber <0 chemicals
p <- ggplot(data=toxcast.ber,
            aes(y=modl_ga,
                x=reorder(chnm, -ber.95))) +
  geom_boxplot() +
  geom_point(data=invitro.dist[chemcas %in% ber.cas],
             aes(y=logAC50p5,
                 x=reorder(chnm, -ber.95)),
             shape=17,
             color="#E69F00", size=3)+
  scale_y_continuous(breaks=seq(-6,3,1))+
  ylab('log10(AC50; micromolar)')+
  #xlab('Frequency')+
  theme(axis.title.y = element_blank())+
  coord_flip()

hist_top <- ggplot() + geom_histogram(aes(invitro.dist$logAC50p5)) +
  theme(axis.title=element_blank()
        #axis.text.x = element_blank()
  )+
  scale_x_continuous(breaks=seq(-5,2,1))+
  coord_cartesian(xlim=c(-5,2))

ber.box <- ggarrange(hist_top,p,
                     heights=c(1,4),
                     ncol=1, nrow=2,
                     align='hv')

p2 <- ggplot(data=pod.ratio.master[CASRN %in% ber.cas],
             aes(x=log.expo.95,
                 y=reorder(Name, -ber.95))) +
  geom_point(size=2) +
  coord_cartesian(xlim=c(-7,-1))+
  scale_x_continuous(breaks=seq(-7,-1,1))+
  xlab('log10(ExpoCast,95th; mg/kg/day)')+
  theme(axis.title.y = element_blank())


hist_top2 <- ggplot() + geom_histogram(aes(pod.ratio.master$log.expo.95)) +
  theme(axis.title=element_blank())+
  scale_x_continuous(breaks=seq(-7,-1,1))+
  coord_cartesian(xlim=c(-7,-1))

fig5bc <- ggarrange(hist_top2,hist_top,p2,p,
                      heights=c(1,4),
                      ncol=2, nrow=2,
                      align='hv',
                      labels=c("B","C"),
                      font.label=list(size=18, face='bold'))

file.dir <- paste("Figures/", sep="")
file.name <- paste("/Figure_5_BER_expocast_dist_", Sys.Date(), ".png", sep="")
file.path <- paste(file.dir, file.name, sep="")
dir.create(path=file.dir, showWarnings = FALSE, recursive = TRUE)
png(file.path, width = 6000, height = 4000, res=480)
fig5bc
dev.off()

fig5a <- ggplot(data=pod.ratio.master,
                    aes(y=log.expo.95, x=pod.nam.95))+
  geom_point(colour = "black", size = 2, position='jitter')+
  geom_label_repel(data=pod.ratio.master[ber.95<0],
                   aes(y=log.expo.95, 
                       x=pod.nam.95), fill=NA, label=pod.ratio.master[ber.95<0]$Name, color="black", size=3)+
  
  theme_bw()+
  geom_vline(lty='dashed', color='red', xintercept=-1.229734)+
  geom_hline(lty='dashed', color='red', yintercept=-5.189432)+
  theme(axis.title.x = element_text(face='bold', size=14),
        axis.title.y = element_text(face='bold', size=14),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12))+
  scale_y_continuous(breaks=seq(-6,-2,1))+
  scale_x_continuous(breaks=seq(-7.5,3,1))+
  coord_cartesian(xlim=c(-7.5,3), ylim=c(-6,-2))+
  ylab('log10 ExpoCast SEEM2 95th-%ile')+
  xlab('log10 POD-NAM,95')+
  theme(legend.position='none')

fig5aplot <- ggarrange(fig5a,
                    heights=c(1),
                    ncol=1, nrow=1,
                    #align='hv',
                    labels=c("A"),
                    font.label=list(size=18, face='bold'))


file.dir <- paste("Figures/", sep="")
file.name <- paste("/Figure_5a_expo_v_podnam_", Sys.Date(), ".png", sep="")
file.path <- paste(file.dir, file.name, sep="")
dir.create(path=file.dir, showWarnings = FALSE, recursive = TRUE)
png(file.path, width = 3000, height = 3000, res=600)
fig5aplot
dev.off()




#-----------------------------------------------------------------------------------#
# Figure 6 ExpoCast v Health Canada Exposure Assessment Values
#-----------------------------------------------------------------------------------#

# source information needed from hc.exp
load("../other_data/hc_exposure_total.RData")
colnames(hc.exp)

expo.lm <- summary((lm(log.expocast.total.95 ~log.max.Consumer.or.Env.mgkgd, data = hc.exp)))
expo.lm
#Residuals:
#  Min      1Q  Median      3Q     Max 
#-2.1152 -0.7175  0.2035  0.7589  1.4268 


exp.plot1 <- ggplot(data=hc.exp,
                   aes(y=log.expocast.total.med, x=log.max.Consumer.or.Env.mgkgd))+
  geom_point(colour = "black", size = 2, position='jitter')+
  geom_label_repel(aes(y=log.expocast.total.med, 
                 x=log.max.Consumer.or.Env.mgkgd), fill=NA, label=hc.exp$CASRN, color="black", size=3)+
  
  geom_smooth(aes(y=log.expocast.total.med, x=log.max.Consumer.or.Env.mgkgd), method='lm', se=FALSE, color='#999999')+
  #geom_abline(size=1, color='#999999')+
  #geom_segment(aes(x=log.max.Consumer.or.Env.mgkgd, 
  #                 y=log.expocast.total.med, 
  #                 xend=log.max.Consumer.or.Env.mgkgd, 
  #                 yend=log.expocast.total.95), colour="black",size=1, linetype='dashed', alpha=0.5)+
  theme_bw()+
  theme(axis.title.x = element_text(face='bold', size=14),
        axis.title.y = element_text(face='bold', size=14),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12))+
  scale_y_continuous(breaks=seq(-7,-1,1))+
  scale_x_continuous(breaks=seq(-6,1.5,1))+
  coord_cartesian(xlim=c(-6,1.5), ylim=c(-7,-1))+
  ylab('log10 ExpoCast SEEM2 Median')+
  xlab('log10 Health Canada Max')+
  theme(legend.position='none')


exp.plot2 <- ggplot(data=hc.exp,
                    aes(y=log.expocast.total.95, x=log.max.Consumer.or.Env.mgkgd))+
  geom_point(colour = "black", size = 2, position='jitter')+
  geom_label_repel(aes(y=log.expocast.total.95, 
                       x=log.max.Consumer.or.Env.mgkgd), fill=NA, label=hc.exp$CASRN, color="black", size=3)+
  
  geom_smooth(aes(y=log.expocast.total.95, x=log.max.Consumer.or.Env.mgkgd), method='lm', se=FALSE, color='#999999')+
 
  #geom_segment(aes(x=log.max.Consumer.or.Env.mgkgd, 
  #                 y=log.expocast.total.med, 
  #                 xend=log.max.Consumer.or.Env.mgkgd, 
  #                 yend=log.expocast.total.95), colour="black",size=1, linetype='dashed', alpha=0.5)+
  theme_bw()+
  theme(axis.title.x = element_text(face='bold', size=14),
        axis.title.y = element_text(face='bold', size=14),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12))+
  scale_y_continuous(breaks=seq(-7,-1,1))+
  scale_x_continuous(breaks=seq(-6,1.5,1))+
  coord_cartesian(xlim=c(-6,1.5), ylim=c(-7,-1))+
  ylab('log10 ExpoCast SEEM2 95th%-ile')+
  xlab('log10 Health Canada Max')+
  theme(legend.position='none')

exp.comp <- ggarrange(exp.plot1,exp.plot2,
                      widths =c(1,1),
                      ncol=2,
                      align='hv',
                      labels=c("A","B"),
                      font.label=list(size=18, face='bold'))



file.dir <- paste("Figures/", sep="")
file.name <- paste("/Figure_6_HC_exposure_comp_Feb2019_", Sys.Date(), ".png", sep="")
file.path <- paste(file.dir, file.name, sep="")
dir.create(path=file.dir, showWarnings = FALSE, recursive = TRUE)
png(file.path, width = 4500, height = 3000, res=600)
exp.comp
dev.off()

#-----------------------------------------------------------------------------------#
# Figure 7 POD ratio distribution and max AED evaluation
#-----------------------------------------------------------------------------------#

# POD ratio distribution 

median(pod.ratio.master$pod.ratio.95) #2.085778
mean(pod.ratio.master$pod.ratio.95) #2.047868
mad(pod.ratio.master$pod.ratio.95) # 1.578271
sd(pod.ratio.master$pod.ratio.95) #1.596057
range(pod.ratio.master$pod.ratio.95) #-2.681086  7.518733

median(pod.ratio.master$pod.ratio.50) #1.19186
mean(pod.ratio.master$pod.ratio.50) #1.260497
mad(pod.ratio.master$pod.ratio.50) #1.494351
sd(pod.ratio.master$pod.ratio.50) #1.550665
range(pod.ratio.master$pod.ratio.50) #-2.921813  7.053226

# histogram

hist.POD.ratio <- ggplot(data=pod.ratio.master) +
  geom_histogram(aes(pod.ratio.master$pod.ratio.95,fill="#0072B2"),
                binwidth = 0.5, 
                alpha=0.25)+
  geom_histogram(aes(pod.ratio.master$pod.ratio.50,fill="#D55E00"),
                 binwidth = 0.5, 
                 alpha=0.25)+
  geom_vline(xintercept=0, linetype='solid', color='black', size=1) +
  geom_vline(xintercept=2.09, linetype='dashed', color="#0072B2", size=1)+
  geom_vline(xintercept=1.19, linetype='dashed', color="#D55E00", size=1)+
  theme(panel.border = element_blank(), 
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.title.x = element_text(face='bold', size=14),
        axis.title.y = element_text(face='bold', size=14),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12))+
  scale_y_continuous(breaks=seq(0,110,10))+
  scale_x_continuous(breaks=seq(-5,10,1))+
  xlab('log10 POD-traditional:POD-NAM Ratio')+
  ylab('Frequency')+
  scale_fill_manual("POD-NAM type", values=c("#0072B2","#D55E00"),labels=c("POD-NAM, 95","POD-NAM, 50"))

# Compare max ADE (using 100 uM) and minimum POD traditional

pod.ratio.master[,max.aed.ratio := log.p5.POD - log.max.AED, by=list(CASRN)]

pod.ratio.lt0.maxaed.lt.tox <- nrow(pod.ratio.master[pod.ratio.95 < 0 & log.max.AED < log.p5.POD]) #0
pod.ratio.gt0.maxaed.lt.tox <- nrow(pod.ratio.master[pod.ratio.95 > 0 & log.max.AED < log.p5.POD]) #240
pod.ratio.gt0 <- nrow(pod.ratio.master[pod.ratio.95 > 0]) #400
pod.ratio.gt0.maxaed.gt.tox <- pod.ratio.gt0 - pod.ratio.gt0.maxaed.lt.tox #160
pod.ratio.lt0 <- 448-pod.ratio.gt0

# percent of the time that the max AED is less than the 5th percentile POD when the pod ratio is greater than zero
(pod.ratio.gt0.maxaed.lt.tox/pod.ratio.gt0)*100 #60%

# percent of the time when the max AED is less than the 5th percentile POD regardless of POD ratio
(pod.ratio.gt0.maxaed.gt.tox/448)*100 #36%

160/448 #36%
240/448 #54%
48/448 #11%

max.aed.summary <- data.frame(Condition = c('log10 POD ratio, 95 > 0',
                                            'log10 POD ratio, 95 < 0'),
                              maxADE_lt_minPOD = c(pod.ratio.gt0.maxaed.lt.tox, pod.ratio.lt0.maxaed.lt.tox),
                              maxADE_gt_minPOD = c(pod.ratio.gt0.maxaed.gt.tox, pod.ratio.lt0))

names(max.aed.summary) <- c("Condition", "Max AED < 5th-%ile POD", "Max AED > 5th-%ile POD")

ttheme_default <- ttheme_default()

tbl <- tableGrob(max.aed.summary, 
                 rows = NULL,
                 theme=ttheme_default)
get.gpar()
print(max.aed.summary)


fig.aed <-plot_grid(hist.POD.ratio, tbl, 
               labels=c("A", "B"),
               label_size=18,
               label_fontface = 'bold')


file.dir <- paste("../figures/", sep="")
file.name <- paste("/Figure_7_POD_ratio_distribution", Sys.Date(), ".png", sep="")
file.path <- paste(file.dir, file.name, sep="")
dir.create(path=file.dir, showWarnings = FALSE, recursive = TRUE)
png(file.path, width = 7500, height = 2500, res=600)
fig.aed
dev.off()

#-----------------------------------------------------------------------------------#
# Figure 9 Study type enrichment in POD ratio < 0 set
#-----------------------------------------------------------------------------------#

# Hypergeometric test to see if there is enrichment of study type in pod.ratio < 0

load('../toxval/pod_master_ra_uniq.RData')

r.d.TP <- nrow(pod.ra.uniq[pod.ratio.binary ==0 & min.pod.ra.class=='developmental / reproductive']) #3
r.d.FP <- nrow(pod.ra.uniq[pod.ratio.binary ==1 & min.pod.ra.class=='developmental / reproductive']) #41
pod.ratio.zero <- nrow(pod.ra.uniq[pod.ratio.binary==0]) #48
pod.ratio.one <- nrow(pod.ra.uniq[pod.ratio.binary==1]) #400

# are repro_dev studies over-represented in the min.pod.ra.class? #no
## (success-in-sample, success-in-left-part, failure-in-sample, failure-in-left-part)
fisher.test(matrix(c(r.d.TP,r.d.FP,(pod.ratio.zero-r.d.TP),pod.ratio.one-r.d.FP), nrow=2, ncol=2),
            alternative='greater')

# 	Fisher's Exact Test for Count Data

#data:  
#  p-value = 0.8772
#alternative hypothesis: true odds ratio is greater than 1
#95 percent confidence interval:
#  0.1477221       Inf
#sample estimates:
#  odds ratio 
#0.5843341 

# are chronic studies over-represented in the min.pod.ra.class? #no

chr.TP <- nrow(pod.ra.uniq[pod.ratio.binary ==0 & min.pod.ra.class=='chronic']) #28
chr.FP <- nrow(pod.ra.uniq[pod.ratio.binary ==1 & min.pod.ra.class=='chronic']) #244

fisher.test(matrix(c(chr.TP,chr.FP,(pod.ratio.zero-chr.TP),pod.ratio.one-chr.FP),
                   nrow=2,ncol=2),alternative="greater")

#Fisher's Exact Test for Count Data

# data:  
#p-value = 0.6984
#alternative hypothesis: true odds ratio is greater than 1
#95 percent confidence interval:
#  0.5154177       Inf
#sample estimates:
#  odds ratio 
#0.8953074 

# make data.frame: dev/repro

devrep.summary <- data.frame(Condition = c('log10-POD ratio,95 < 0','log10-POD ratio,95 > 0'),
                             Dev_repro_E_minPOD = c(r.d.TP, r.d.FP),
                             Dev_repro_DNE_minPOD = c(pod.ratio.zero-r.d.TP, pod.ratio.one-r.d.FP))
names(devrep.summary) <- c("Condition", "Dev/Repro is min POD", "Dev/Repro is not min POD")
ttheme_default <- ttheme_default()

# make table: dev/repro

devrep.tbl <- tableGrob(devrep.summary, 
                        rows = NULL,
                        theme=ttheme_default)

# make data.frame: chronic

chr.summary <- data.frame(Condition = c('log10-POD ratio,95 < 0','log10-POD ratio,95 > 0'),
                          chr_E_minPOD = c(chr.TP, chr.FP),
                          chr_DNE_minPOD = c(pod.ratio.zero-chr.TP, pod.ratio.one-chr.FP))

names(chr.summary) <- c("Condition", "Chronic is min POD", "Chronic is not min POD")

# make table: chronic

chr.tbl <- tableGrob(chr.summary, 
                     rows = NULL,
                     theme=ttheme_default)

# arrange tables in a figure

fisher.tables <- ggarrange(devrep.tbl, chr.tbl,nrow=2,
                           align='v', heights=c(2,2))

file.dir <- paste("Figures/", sep="")
file.name <- paste("/Figure_9_study_type_tables_", Sys.Date(), ".png", sep="")
file.path <- paste(file.dir, file.name, sep="")
dir.create(path=file.dir, showWarnings = FALSE, recursive = TRUE)
png(file.path, width = 4000, height = 1500, res=600)
fisher.tables
dev.off()

#-----------------------------------------------------------------------------------#
# Figure 10 POD-NAM compared to TTC
#-----------------------------------------------------------------------------------#

# TTC ratio distribution 
# ttc ratio is defined as pod-nam - ttc

ttc.ratio.median <- median(pod.ratio.master$ttc.ratio.95, na.rm=TRUE) #2.25
mean(pod.ratio.master$ttc.ratio.95, na.rm=TRUE) #2.35
mad(pod.ratio.master$ttc.ratio.95, na.rm=TRUE) #1.93
sd(pod.ratio.master$ttc.ratio.95, na.rm=TRUE) #1.96
range(pod.ratio.master$ttc.ratio.95, na.rm=TRUE) #-4.215299  7.432732

pod.ratio.master.ttc <- pod.ratio.master[ttc.label %in% c('0.0025','0.3','1.5','30','9')]

pod.ratio.master.ttc[ttc.ratio.95<0]
nrow(pod.ratio.master.ttc) #433
nrow(pod.ratio.master.ttc[ttc.ratio.95<0]) #44
44/433 #10% of the time

hist.ttc <- ggplot(data=pod.ratio.master.ttc,
                   aes(pod.ratio.master.ttc$ttc.ratio.95)) +
  geom_histogram(binwidth = 0.5, 
                 #center=0, 
                 #boundary=0.5, 
                 #fill='dark gray',
                 alpha=0.8
                 #col='dark gray'
  )+
  geom_vline(xintercept=0, linetype='dashed', color='black', size=1) +
  geom_vline(xintercept=ttc.ratio.median, linetype='solid', color='black', size=1)+
  theme(panel.border = element_blank(), 
        #panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.title.x = element_text(face='bold', size=14),
        axis.title.y = element_text(face='bold', size=14),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12))+
  scale_y_continuous(breaks=seq(0,110,10))+
  scale_x_continuous(breaks=seq(-5,10,1))+
  xlab('log10 POD-NAM 95:TTC Ratio')+
  ylab('Frequency')

# TTC vs POD
colnames(pod.ratio.master)
pod.ratio.master[, log.ttc.mkd.label := as.factor(round(log.ttc.mkd,2))]
count(pod.ratio.master$ttc.label)
pod.ratio.master[ttc.label=='Exclusion', log.ttc.mkd.label := 'Exclusion']
pod.ratio.master[ttc.label=='No Structure', log.ttc.mkd.label := 'NS']
str(pod.ratio.master$log.ttc.mkd.label)
pod.ratio.master <- pod.ratio.master[order(log.ttc.mkd.label)]
count(pod.ratio.master$log.ttc.mkd.label)
#x freq
#1      -5.6  141
#2     -3.52   36
#3     -2.82  212
#4     -2.05    5
#5     -1.52   39
#6 Exclusion   12
#7        NS    3

# make this a violin plot
ttc.by.group <- ggplot(data=pod.ratio.master, 
                       aes(x=log.ttc.mkd.label, 
                           y=pod.nam.95, 
                           group=log.ttc.mkd.label,
                           color=log.ttc.mkd.label,
                           shape=log.ttc.mkd.label)) +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75), size=1, colour='#000000')+
  geom_jitter()+
  theme(panel.border = element_blank(), 
        #panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.title.x = element_text(face='bold', size=14),
        axis.title.y = element_text(face='bold', size=14),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12))+
  scale_colour_manual(values=c('#999999','#999999','#999999','#999999','#999999','#999999','#999999'),
                      labels=c('-5.60','-3.52','-2.82', '-2.05', '-1.52', 'Exclusion', 'NS')) +
  scale_shape_manual(values=c(19,19,19,19,19,19,19),
                     labels=c('-5.60','-3.52','-2.82', '-2.05', '-1.52', 'Exclusion', 'NS')) +
  scale_y_continuous(breaks=seq(-8,3,1))+
  #scale_x_continuous(breaks=seq(-5,10,1))+
  xlab('log10 TTC')+
  ylab('log10 POD-NAM, 95')+
  theme(legend.position='none')

# now put the two TTC plots together as A and B

ttc.fig<-plot_grid(hist.ttc, ttc.by.group, 
                   labels=c("A", "B"),
                   label_size=18,
                   label_fontface = 'bold')


file.dir <- paste("../figures/", sep="")
file.name <- paste("/Figure_10_TTC_comparison_", Sys.Date(), ".png", sep="")
file.path <- paste(file.dir, file.name, sep="")
dir.create(path=file.dir, showWarnings = FALSE, recursive = TRUE)
png(file.path, width = 6600, height = 4000, res=480)
ttc.fig
dev.off()

#-----------------------------------------------------------------------------------#
# Supplemental Figures
#-----------------------------------------------------------------------------------#

#---------------------------------------------------------------------#
# Supplemental Figure 1 on hitpercent and flags for ToxCast
#---------------------------------------------------------------------#

load(file="../toxcast/source_invitro_data.RData")
# relationship between hitpct and modl_ga_delta is not clear; correlation but mostly looks like a blob
# cannot necessarily infer that large modl_ga_delta means a low hitpct

summary(lm(hitpct ~ modl_ga_delta, data=invitrodb.hitc1.mc5))
cor.test(invitrodb.hitc1.mc5$hitpct, invitrodb.hitc1.mc5$modl_ga_delta, method='spearman')

hitpct <- ggplot(invitrodb.hitc1.mc5, aes(x=modl_ga_delta, y=hitpct)) + 
  geom_point() +
  theme_few()+geom_smooth(method=lm)
hitpct

# how many curves have low hitpct?

nrow(invitrodb.hitc1.mc5[hitpct>0.625]) 
length(unique(invitrodb.hitc1.mc5[hitpct>0.625]$casn)) 
nrow(invitrodb.hitc1.mc5[hitpct>0.5]) 
count(invitrodb.hitc1.mc5[hitpct<0.5]$mc6_flags)
# most prevalent patterns are: 17,16,11,6  (259 times) and 17,16,11,7  (130 times)
count(invitrodb.hitc1.mc5[hitpct<0.3]$mc6_flags)
nrow(invitrodb.hitc1.mc5[hitpct<0.5 & flag.length > 2]) 
range(invitrodb.hitc1.mc5$flag.length, na.rm=T) # 1 to 6 flags

flags.hitpct <- ggplot(data=toxcast.master, aes(x=factor(flag.length), y=hitpct))+
  geom_boxplot()+
  xlab('Number of Caution Flags on ToxCast Curve Fit')+
  ylab('Hit percent')

flags.hitpct.unfiltered <- ggplot(data=invitrodb.hitc1.mc5, aes(x=factor(flag.length), y=hitpct))+
  geom_boxplot()+
  xlab('Number of Caution Flags on ToxCast Curve Fit')+
  ylab('Hit percent')

double.flags.hitpct <- ggarrange(
  flags.hitpct.unfiltered,
  flags.hitpct,
  widths=c(1,1),
  ncol=2, nrow=1,
  align='h',
  labels=c("A","B"),
  font.label=list(size=18, face='bold'))

file.dir <- paste("Figures/", sep="")
file.name <- paste("/Figure_double_flags_hitpct_", Sys.Date(), ".png", sep="")
file.path <- paste(file.dir, file.name, sep="")
dir.create(path=file.dir, showWarnings = FALSE, recursive = TRUE)
png(file.path, width = 6000, height = 4000, res=480)
double.flags.hitpct
dev.off()

#---------------------------------------------------------------------#
# Supplemental Figure 2 on in vitro data and filtering
#---------------------------------------------------------------------#

#need to use the hitcsum to understand how many AC50s will be in the distribution

load(file='../toxcast/toxcast_master.RData')
hmc5 <- unique(toxcast.master[,list(casn,chnm,hitcsum)])
count(hmc5, 'hitcsum')
median(hmc5$hitcsum, na.rm=TRUE) #56
range(hmc5$hitcsum) #0 1351

hist.hitcsum <- ggplot(data=hmc5, aes(hmc5$hitcsum))+
  geom_histogram(binwidth=10)+
  theme(panel.border = element_blank(), 
        #panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.title = element_text(face='bold', size=20),
        axis.text = element_text(size=18))+
  xlab('Hitcall Sum')+
  ylab('Frequency')+
  scale_x_continuous(breaks=seq(0,1300,200))
hist.hitcsum

hmc5[hitcsum < 100]
hmc5[casn=='50-06-6']
file.dir <- paste("Figures/", sep="")
file.name <- paste("/Supp_Figure_2_hitcsum_", Sys.Date(), ".png", sep="")
file.path <- paste(file.dir, file.name, sep="")
dir.create(path=file.dir, showWarnings = FALSE, recursive = TRUE)
png(file.path, width = 4000, height = 4000, res=480)
hist.hitcsum
dev.off()

#---------------------------------------------------------------------#
# Supplemental Figure 3 on using 5th percentile AC50
#---------------------------------------------------------------------#

#-----------------------------------------------------------------------------------#
# Plot distribution of AC50s versus 5th-ile for all chemicals
#-----------------------------------------------------------------------------------#

# make long format data for ggplot2

invitro.dist.long <- melt.data.table(invitro.dist,
                                     id.vars = c('chemcas', 'chnm','hitcsum','avg.dist'),
                                     measure.vars = c('min_modl_ga_uM', 'AC50p5', 'min.ec10.um'),
                                     variable.name = c('percentile'))

invitro.dist.long <- invitro.dist.long[order(-avg.dist, chnm)]

# plot compares ToxCast min AC50, 5th percentile, and the HIPPTox minEC10 (only 50 chem)
invitro.dist.plot <- ggplot(data=invitro.dist.long,
                            aes(x=value,
                                y=reorder(strtrim(chnm,30), -avg.dist),
                                group=factor(percentile))) +
  geom_point(aes(col=factor(percentile)))+
  scale_colour_manual(values=c("#999999","#0072B2", "#E69F00"),
                      labels=c("min ToxCast", "5%", "min HIPPTox")) +
  scale_shape_manual(values=c(19,19, 19),
                     labels=c("min ToxCast", "5%", "min HIPPTox")) +
  scale_x_log10()+
  theme_bw() +
  coord_flip()+
  theme(panel.border = element_blank(), 
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text = element_text(size=10),
        axis.title = element_text(size=12, face='bold'))+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())+
  theme(legend.position="bottom", legend.title=element_blank())+
  xlab('log10 micromolar')+
  ylab('Chemical')

file.dir <- paste("Figures/", sep="")
file.name <- paste("/Supp_Figure_3_min_v_5th", Sys.Date(), ".png", sep="")
file.path <- paste(file.dir, file.name, sep="")
dir.create(path=file.dir, showWarnings = FALSE, recursive = TRUE)
png(file.path, width = 5000, height = 4000, res=480)
invitro.dist.plot
dev.off() 

#---------------------------------------------------------------------#
# Supplemental Figure 4 on Frequency of ToxCast assay driving min AC50
#---------------------------------------------------------------------#

#-----------------------------------------------------------------------------------#
# Plot distribution of AC50s versus 5th-ile for all chemicals
#-----------------------------------------------------------------------------------#

# make long format data for ggplot2

invitro.dist.long <- melt.data.table(invitro.dist,
                                     id.vars = c('chemcas', 'chnm','hitcsum','avg.dist'),
                                     measure.vars = c('min_modl_ga_uM', 'AC50p5', 'min.ec10.um'),
                                     variable.name = c('percentile'))

invitro.dist.long <- invitro.dist.long[order(-avg.dist, chnm)]

# plot compares ToxCast min AC50, 5th percentile, and the HIPPTox minEC10 (only 50 chem)
invitro.dist.plot <- ggplot(data=invitro.dist.long,
                            aes(x=value,
                                y=reorder(strtrim(chnm,30), -avg.dist),
                                group=factor(percentile))) +
  geom_point(aes(col=factor(percentile)))+
  scale_colour_manual(values=c("#999999","#0072B2", "#E69F00"),
                      labels=c("min ToxCast", "5%", "min HIPPTox")) +
  scale_shape_manual(values=c(19,19, 19),
                     labels=c("min ToxCast", "5%", "min HIPPTox")) +
  scale_x_log10()+
  theme_bw() +
  coord_flip()+
  theme(panel.border = element_blank(), 
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text = element_text(size=10),
        axis.title = element_text(size=12, face='bold'))+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())+
  theme(legend.position="bottom", legend.title=element_blank())+
  xlab('log10 micromolar')+
  ylab('Chemical')

file.dir <- paste("Figures/", sep="")
file.name <- paste("/Figure_min_v_5th", Sys.Date(), ".png", sep="")
file.path <- paste(file.dir, file.name, sep="")
dir.create(path=file.dir, showWarnings = FALSE, recursive = TRUE)
png(file.path, width = 5000, height = 4000, res=480)
invitro.dist.plot
dev.off() 

#---------------------------------------------------------------------#
# Supplemental information on HTTK
#---------------------------------------------------------------------#

# on httk
load('AED_from_httk_18Dec18.RData')
colnames(aed.df)
aed.df[,'AED_5'] <- as.numeric(as.character(aed.df[,'AED_5']))
aed.df[,'AED_50'] <- as.numeric(as.character(aed.df[,'AED_50']))
aed.df[,'AED_95'] <- as.numeric(as.character(aed.df[,'AED_95']))
aed.df[,'AED_hipptox_50'] <- as.numeric(as.character(aed.df[,'AED_hipptox_50']))
aed.df[,'AED_hipptox_95'] <- as.numeric(as.character(aed.df[,'AED_hipptox_95']))
aed.df[,'max_AED'] <- as.numeric(as.character(aed.df[,'max_AED']))

aed <- as.data.table(aed.df)
aed[,log.AED.5 := log10(AED_5)]
aed[,log.AED.50 := log10(AED_50)]
aed[,log.AED.95 := log10(AED_95)]
aed[,log.AED.hipptox.50 := log10(AED_hipptox_50)]
aed[,log.AED.hipptox.95 := log10(AED_hipptox_95)]
aed[,log.max.AED := log10(max_AED)]

aed[,ind.var := AED_5/AED_95]
aed[,med.95 := AED_50/AED_95]

median(aed$ind.var) #19.2
mean(aed$ind.var) #26
range(aed$ind.var) #4.483931 94.388409
sd(aed$ind.var) #186.562

nrow(aed[ind.var < 19.2]) #224
(224/448)*100 #50%
nrow(aed[ind.var < 26]) #265

median(aed$med.95) #6.3
mean(aed$med.95) #7.6
range(aed$med.95) #1.740716 18.603057
sd(aed$med.95) #4.6
nrow(aed[med.95 < 6.3]) #225


# create supplemental Figure 5 on the size of the AED5 to AED95 interval

fold.hist.50.95 <- ggplot(data=aed, aes(aed$med.95)) +
  geom_histogram(binwidth = 0.5, 
                 alpha=0.8)+
  theme(panel.border = element_blank(), 
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.title.x = element_text(face='bold', size=12),
        axis.title.y = element_text(face='bold', size=12),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12))+
  scale_y_continuous(breaks=seq(0,80,10))+
  scale_x_continuous(breaks=seq(0,25,2))+
  coord_cartesian(ylim=c(0,80))+
  xlab('AED50/AED95')+
  ylab('Frequency')

fold.hist.5.95 <-ggplot(data=aed, aes(aed$ind.var)) +
  geom_histogram(binwidth = 2, 
                 alpha=0.8)+
  theme(panel.border = element_blank(), 
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.title.x = element_text(face='bold', size=12),
        axis.title.y = element_text(face='bold', size=12),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12))+
  scale_y_continuous(breaks=seq(0,80,10))+
  scale_x_continuous(breaks=seq(0,100,10))+
  coord_cartesian(ylim=c(0,80))+
  xlab('AED5/AED95')+
  ylab('Frequency')

aed.int <- ggarrange(fold.hist.5.95,fold.hist.50.95,
                     #heights=c(1,4),
                     widths=c(1,1),
                     ncol=2, nrow=1,
                     align='hv',
                     labels=c("A","B"),
                     font.label=list(size=18, face='bold'))

file.dir <- paste("Figures/", sep="")
file.name <- paste("/Supp_Figure_5_httk_interindivid_var_ab_8feb19", Sys.Date(), ".png", sep="")
file.path <- paste(file.dir, file.name, sep="")
dir.create(path=file.dir, showWarnings = FALSE, recursive = TRUE)
png(file.path, width = 6000, height = 3000, res=600)
aed.int
dev.off()

# Supplemental Figure 6

# show the effects of restrictive v nonrestrictive clearance

# first calculate with and without restrictive clearance for AED95
invitro.dist.subset<-as.data.frame(subset(invitro.dist, chemcas %in% get_cheminfo(species='Human')))

length(unique(invitro.dist.subset$chemcas)) # 448
colnames(invitro.dist.subset)
aed.df2=data.frame()
for (casn in invitro.dist.subset[,'chemcas'])
{
  AED_95_rc<-(calc_mc_oral_equiv(conc=subset(invitro.dist.subset,chemcas==casn)[['AC50p5']], chem.cas=casn, 
                                 which.quantile=c(0.95), species='Human', method='dr', 
                                 well.stirred.correction=T, restrictive.clearance=T, output.units='mg'))
  AED_95_nrc<-(calc_mc_oral_equiv(conc=subset(invitro.dist.subset,chemcas==casn)[['AC50p5']], chem.cas=casn, 
                                  which.quantile=c(0.95), species='Human', method='dr', 
                                  well.stirred.correction=T, restrictive.clearance=F, output.units='mg'))
  aed.df2<-rbind(aed.df2,cbind(casn, AED_95_rc, AED_95_nrc))
}

save(aed.df2, file='AED_from_httk_rc_v_nrc.RData')

aed.df2[,'AED_95_rc'] <- as.numeric(as.character(aed.df2[,'AED_95_rc']))
aed.df2[,'AED_95_nrc'] <- as.numeric(as.character(aed.df2[,'AED_95_nrc']))

aed2 <- as.data.table(aed.df2)
aed2[,log.AED.95.rc := log10(AED_95_rc)]
aed2[,log.AED.95.nrc := log10(AED_95_nrc)]

nrc.v.rc <- ggplot(data=aed2,
                   aes(x=log.AED.95.rc, y=log.AED.95.nrc))+
  geom_point(colour = "black", size = 2, position='jitter')+
  geom_abline(size=1, color='#999999')+
  theme_bw()+
  theme(axis.title.x = element_text(face='bold', size=14),
        axis.title.y = element_text(face='bold', size=14),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12))+
  ylab('AED, 95%ile, nonrestrictive')+
  xlab('AED, 95%ile, restrictive')

file.dir <- paste("Figures/", sep="")
file.name <- paste("/Supp_Figure_6_httk_resCL_", Sys.Date(), ".png", sep="")
file.path <- paste(file.dir, file.name, sep="")
dir.create(path=file.dir, showWarnings = FALSE, recursive = TRUE)
png(file.path, width = 3000, height = 3000, res=600)
nrc.v.rc
dev.off()

#-----------------------------------------------------------------------------------#
# Supplemental Figure 7 Comparison of POD-NAM to POD-traditional using allometric scaling
#-----------------------------------------------------------------------------------#

# HED POD ratio distribution 

pod.ratio.hed.master <- pod.ratio.master[!is.na(hed.pod.ratio.95)]
pod.ratio.hed.median <- median(pod.ratio.hed.master$hed.pod.ratio.95) #1.33
mean(pod.ratio.hed.master$hed.pod.ratio.95) #1.35
mad(pod.ratio.hed.master$hed.pod.ratio.95) # 1.55
sd(pod.ratio.hed.master$hed.pod.ratio.95) # 1.58
range(pod.ratio.hed.master$hed.pod.ratio.95) #-3.281239  6.728248

# POD ratio distribution

hist.hed.POD.ratio <- ggplot(data=pod.ratio.hed.master, aes(pod.ratio.hed.master$hed.pod.ratio.95)) +
  geom_histogram(binwidth = 0.5, 
                 #center=0, 
                 #boundary=0.5, 
                 #fill='dark gray',
                 alpha=0.8
                 #col='dark gray'
  )+
  geom_vline(xintercept=-0.5, linetype='dashed', color='black', size=1) +
  geom_vline(xintercept=pod.ratio.hed.median, linetype='solid', color='black', size=1)+
  theme(panel.border = element_blank(), 
        #panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.title.x = element_text(face='bold', size=14),
        axis.title.y = element_text(face='bold', size=14),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12))+
  scale_y_continuous(breaks=seq(0,110,10))+
  scale_x_continuous(breaks=seq(-5,10,1))+
  xlab('log10 POD-traditional:POD-NAM 95 Ratio')+
  ylab('Frequency')

file.dir <- paste("Figures/", sep="")
file.name <- paste("/Figure_11_HED_POD_ratio_", Sys.Date(), ".png", sep="")
file.path <- paste(file.dir, file.name, sep="")
dir.create(path=file.dir, showWarnings = FALSE, recursive = TRUE)
png(file.path, width = 3000, height = 3000, res=600)
hist.hed.POD.ratio
dev.off()






