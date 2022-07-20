# script to plot the coverage
#### script for R version = r/3.6.1
# libraries
library(data.table)
library(zoo)
library(ggplot2)
library(tidyverse)
library(dplyr)

args=(commandArgs(TRUE))

if(length(args)==0){
    stop("No arguments supplied.")
}else{
    for(i in 1:length(args)){
      eval(parse(text=args[[i]]))
    }
}

if( getDTthreads() < 10 || getDTthreads() > 10){
  setDTthreads(threads = 10)
}else{print(paste("WARNING: data.table package is running with ", getDTthreads(), " threads.", sep=''))}
#STEP <- as.numeric(STEP)
STEP <- 1000
#SIZE <- as.numeric(SIZE)
SIZE <- 1000

#SAMPLE <- SAMPLE
SAMPLE <- paste("T20.LC.1", sep = "")


# 1) row coverage
#df <- fread(COVERAGE_FILE) # read the file
df <- fread("depth_alignments.depth") # read the file with the depth coverage for each base (fastq aligned to USMA 521 v2)

names(df) <- c("Chr","Position","depth") # rename the headers

sample.median <- median(df$depth) #estimate the global coverage

depth2plot <- df[, .(window.start = rollapply(Position, width=SIZE, by=STEP, FUN=min, align="left"),window.end = rollapply(Position, width=SIZE, by=STEP, FUN=max, align="left"),coverage.median = rollapply(depth, width=SIZE, by=STEP, FUN=median, align="left")), .(Chr)] # estimate the median coverage by windows, the windows size is defined in width.

depth2plot$Chr <- gsub("USMA_521_v2_", "", depth2plot$Chr)

depth2plot$Chr <- as.numeric(depth2plot$Chr)


plot <- ggplot(depth2plot, aes(x = window.end, y = coverage.median, colour = as.factor(Chr))) + 
  geom_point() + 
  facet_grid(.~ Chr, space = 'free_x', scales = "free_x" ) + 
  geom_hline(yintercept = (sample.median*2), linetype = "dashed", size = 1, color = "gray")+
  geom_hline(yintercept = (sample.median*3), linetype = "dashed", size = 1, color = "gray")+
  theme_bw() + 
  scale_y_continuous(limits = c(0,sample.median*3)) + 
  geom_hline(yintercept = sample.median, linetype = "dashed", color = "red", size = 2) + 
  labs(title = paste("Sample: ", SAMPLE, " (", SIZE, " non-overlapping bp)", sep = ""), 
       subtitle = "Reads aligned to USMA 521 v2 reference",
       x = "Chromosome", 
       y = "Median coverage depth") + 
  theme(legend.position = "none", panel.spacing.x = grid::unit(0, "cm"), 
        panel.border = element_rect(colour = "grey", size = 0.1), panel.ontop = FALSE, 
        axis.title.x = element_text(face = "bold", color = "black", size = 22), 
        axis.title.y = element_text(face = "bold", color = "black", size = 22), 
        axis.text.x = element_blank(), axis.text.y = element_text(color = "black", size = 18), 
        axis.ticks.x= element_blank(), 
        plot.title = element_text(face = "bold", color = "red", size = 24, hjust = 0.5), 
        plot.subtitle = element_text( color = "red", size = 20, hjust = 0.5), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        strip.background = element_blank(), strip.text = element_text(size = 14, face = "bold")) # make the plot
plot
#ggsave(OUTPUT_PLOT, plot = plot, width = 20, height = 12, units = "in", dpi = 300) # save the plot
ggsave("Plot1.bwa.USMA.521.v2.brut.png", plot = plot, width = 20, height = 12, units = "in", dpi = 300) # save the plot
ggsave("Plot2.bwa.USMA.521.v2.brut.png", plot = plot, width = 20, height = 12, units = "in", dpi = 300) # save the plot


## calculate the normalized coverage

# 1) median coverage by chr
depth4chr <- df %>% group_by(Chr) %>% summarise(Median.cov.chr = median(depth))
depth4chr$Norm.global.coverage <- depth4chr$Median.cov.chr/sample.median 


depth4chr$Chr <- as.numeric(gsub("USMA_521_v2_", "", depth4chr$Chr))

depth2plot <- depth2plot %>% left_join(select(depth4chr, Chr, Median.cov.chr), by = c("Chr" = "Chr"))

Norm.cov <- depth2plot %>% group_by(Chr) %>% mutate(Norm.cov = coverage.median/sample.median) %>% setDT()
Norm.cov <- setDT(Norm.cov)
Norm.cov.mitoch <- Norm.cov[Chr == 24]
Norm.cov <- Norm.cov[Chr != 24]
Norm.cov <- Norm.cov[Chr != 25]
Norm.cov <- Norm.cov[Chr != 26]
Norm.cov <- Norm.cov[Chr != 27]
Norm.cov <- Norm.cov[Chr != 28]

plot.2 <- ggplot(Norm.cov, aes(x = window.end, y = Norm.cov, color = as.factor(Chr))) +
  geom_area(aes(alpha = 0.2, fill = as.factor(Chr))) +
  scale_y_continuous(limits = c(0,5))+
  geom_hline(yintercept = 1) +
  #geom_hline(yintercept = 2, linetype = "dashed", size = 1, color = "gray")+
  #geom_hline(yintercept = 3, linetype = "dashed", size = 1, color = "gray")+
  facet_wrap(.~ Chr, scales = "free_x", strip.position = "bottom") +
  theme_bw() + 
  labs(title = paste("Sample: ", SAMPLE, " (", SIZE, " non-overlapping bp)", sep = ""), 
       subtitle = "Reads aligned to USMA 521 v2 reference",
       x = "Chromosome", 
       y = "Normalized coverage") + 
  theme(legend.position = "none", panel.spacing.x = grid::unit(0, "cm"), 
        panel.border = element_rect(colour = "white", size = 0.1), panel.ontop = FALSE, 
        axis.title.x = element_text(face = "bold", color = "black", size = 22), 
        axis.title.y = element_text(face = "bold", color = "black", size = 22), 
        axis.text.x = element_blank(), axis.text.y = element_text(color = "black", size = 18), 
        axis.ticks.x= element_blank(), 
        plot.title = element_text(face = "bold", color = "red", size = 24, hjust = 0.5), 
        plot.subtitle = element_text( color = "red", size = 20, hjust = 0.5), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        strip.background = element_blank(), strip.text = element_text(size = 14, face = "bold")); plot.2

#ggsave(OUTPUT_PLOT_2, plot = plot.2, width = 20, height = 12, units = "in", dpi = 300) # save the plot
ggsave("Plot3.bwa.USMA.521.v2.mediane.png", plot = plot.2, width = 20, height = 12, units = "in", dpi = 300) # save the plot
ggsave("Plot4.bwa.USMA.521.v2.mediane.png", plot = plot.2, width = 20, height = 12, units = "in", dpi = 300) # save the plot

# extract chr9
df.chr9 <- Norm.cov[Chr == "9"]
plot.3 <- ggplot(df.chr9, aes(x = window.end, y = Norm.cov, color = as.factor(Chr))) +
  geom_area(aes(alpha = 0.2, fill = as.factor(Chr))) +
  scale_y_continuous(limits = c(0,5))+
  geom_hline(yintercept = 1) +
  geom_hline(yintercept = 2, linetype = "dashed", size = 1, color = "gray")+
  geom_hline(yintercept = 3, linetype = "dashed", size = 1, color = "gray")+
  geom_vline(xintercept = 150000, linetype = "dashed", size = 1, color = "red", alpha = 0.4)+
  #facet_wrap(.~ Chr, scales = "free_x", strip.position = "bottom") +
  theme_bw() + 
  labs(title = paste("Sample: ", SAMPLE, " (", SIZE, " non-overlapping bp)", sep = ""), 
       subtitle = "Reads aligned to USMA 521 v2 reference",
       x = "Chromosome 9", 
       y = "Normalized coverage") + 
  theme(legend.position = "none", panel.spacing.x = grid::unit(0, "cm"), 
        panel.border = element_rect(colour = "white", size = 0.1), panel.ontop = FALSE, 
        axis.title.x = element_text(face = "bold", color = "black", size = 22), 
        axis.title.y = element_text(face = "bold", color = "black", size = 22), 
        axis.text.x = element_blank(), axis.text.y = element_text(color = "black", size = 18), 
        axis.ticks.x= element_blank(), 
        plot.title = element_text(face = "bold", color = "red", size = 24, hjust = 0.5), 
        plot.subtitle = element_text( color = "red", size = 20, hjust = 0.5), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        strip.background = element_blank(), strip.text = element_text(size = 14, face = "bold")); plot.3

ggsave("Plot5.bwa.USMA.521.v2.mediane.png", plot = plot.3, width = 20, height = 12, units = "in", dpi = 300) # save the plot

## mitoch

plot.3 <- ggplot(Norm.cov.mitoch, aes(x = window.end, y = Norm.cov, color = as.factor(Chr))) +
  geom_area(aes(alpha = 0.2, fill = as.factor(Chr))) +
  scale_y_continuous()+
  geom_hline(yintercept = 1) +
  facet_wrap(.~ Chr, scales = "free_x", strip.position = "bottom") +
  theme_bw() + 
  labs(title = paste("Sample: ", SAMPLE, " (", SIZE, " non-overlapping bp)", sep = ""), 
       x = "Chromosome", 
       y = "Median coverage depth") + 
  theme(legend.position = "none", panel.spacing.x = grid::unit(0, "cm"), 
        panel.border = element_rect(colour = "white", size = 0.1), panel.ontop = FALSE, 
        axis.title.x = element_text(face = "bold", color = "black", size = 22), 
        axis.title.y = element_text(face = "bold", color = "black", size = 22), 
        axis.text.x = element_blank(), axis.text.y = element_text(color = "black", size = 18), 
        axis.ticks.x= element_blank(), 
        plot.title = element_text(face = "bold", color = "red", size = 24, hjust = 0.5), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        strip.background = element_blank(), strip.text = element_text(size = 14, face = "bold")); plot.2
plot.3
ggsave("Plot6.bwa.USMA.521.v2.mediane.png", plot = plot.3, width = 20, height = 12, units = "in", dpi = 300)

rm(list = ls())
