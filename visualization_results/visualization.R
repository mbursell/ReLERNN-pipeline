setwd("~/Documents/BIO624/output")
library(dplyr)
library(ggplot2)

# Test with example
examp<-read.table("example.PREDICT.BSCORRECTED.txt",header = TRUE, sep = "", dec = ".")
examp_2L<-examp[examp$chrom=="b'2L'",]
ggplot(examp_2L, aes(start, recombRate)) +
  geom_point() +
  geom_line(data = examp_2L) +
  geom_ribbon(aes(ymin = CI95LO, ymax = CI95HI), alpha = 0.2) +
  labs(x = "Genomic Position (bp)", y = "Recombination Rate (c/bp)", title = "Example 2L") +
  theme_minimal()

# Load in data for all three populations
local_pred<-read.table("116k_WS/local_chicken_autosomes_biallelic.PREDICT.txt",header = TRUE, sep = "", dec = ".")
red_pred<-read.table("116k_WS/red_junglefowl_autosomes_biallelic.PREDICT.txt",header = TRUE, sep = "", dec = ".")
yunbao_pred<-read.table("116k_WS/yunbao_chicken_autosomes_biallelic.PREDICT.txt",header = TRUE, sep = "", dec = ".")

# Specify chromosome number
chr_list<-unique(local_pred$chrom) 
i<-28 # i in 1:28
chr_i<-chr_list[i] 
# number of chromosome
n_i<-gsub(".*?([0-9]+).*", "\\1", chr_i)
local_pred_chr_i<-local_pred[local_pred$chrom==chr_i,]
red_pred_chr_i<-red_pred[red_pred$chrom==chr_i,]
yunbao_pred_chr_i<-yunbao_pred[yunbao_pred$chrom==chr_i,]

# Plot chromosome region vs recombination rate
legend_colors <- c("Local Chicken"="blueviolet", "Red Junglefowl"="darksalmon", "Yuanbao Chicken"="lightgreen")

ggplot() +
  #geom_point(size = 0.25,color='blue') +
  ylim(1.8e-08, 2.4e-08) +
  geom_line(data=local_pred_chr_i, aes(start, recombRate, color="Local Chicken"), alpha = 0.75) +
  geom_line(data=red_pred_chr_i, aes(start, recombRate, color="Red Junglefowl"), alpha = 0.75) + 
  geom_line(data=yunbao_pred_chr_i, aes(start, recombRate, color="Yuanbao Chicken"), alpha = 0.75) +
  labs(x = "Genomic Position (bp)", 
       y = "Recombination Rate (c/bp)", 
       title = paste0("Chr ", n_i, ", Predict, 10 epochs, Window size 116kb"),
       colour = "") +
  scale_color_manual(values=legend_colors) +
  annotate("text",  x=Inf, y = Inf, label = paste0("Chr ", n_i), fontface = "bold", size = 14/.pt, vjust=2, hjust=2) +
  theme_bw()

# Boxplot
box_chr_i<-cbind(local_pred_chr_i$recombRate,red_pred_chr_i$recombRate,yunbao_pred_chr_i$recombRate)
colnames(box_chr_i)<-c("Local Chicken","Red Junglefowl","Yuanbao Chicken")
box_chr_i<-as.data.frame(box_chr_i)

require(reshape2)
ggplot(data = melt(box_chr_i), aes(x=variable, y=value, fill=variable)) +
  ylim(1.8e-08, 2.4e-08) +
  geom_boxplot() +
  labs(x = "", 
       y = "Recombination Rate (c/bp)", 
       title =  paste0("Chr ", n_i),
       colour = "") +
  #scale_color_manual(values = c("blueviolet","darksalmon","lightgreen")) +
  scale_fill_manual(values = c("blueviolet","darksalmon","lightgreen"))+
  guides(fill=guide_legend(title="")) +
  #annotate("text",  x=Inf, y = Inf, label = "Chr 1", fontface = "bold", size = 14/.pt, vjust=2, hjust=2) +
  theme_bw()