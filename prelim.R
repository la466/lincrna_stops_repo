####
# prelim 
####

library(ggplot2)

file <- read.csv('results/sim_orf_lengths.csv', head=T)file=read.csv()
cor <- cor.test(file$gc, file$z, method="spearman")
bin <- binom.test(nrow(file[file$z > 0,]), nrow(file), p=0.5, alternative = "g")

p <- ggplot(data=file, aes(x=gc, y=z)) +
  geom_point(color="RoyalBlue", size=0.6) +
  geom_smooth(method="lm", se=T, lwd=0.5, color="black") +
  geom_hline(yintercept=0, lty=2) +
  annotate("text", x=0.7, y=28, label=paste("rho = ", round(cor$estimate, 3), ", P = ", cor$p.value), cex=3) +
  annotate("text", x=0.6, y=26, label=paste("one tail binomial test: ns = ", bin$statistic, ", nt = ", bin$parameter, ", P = ", bin$p.value), cex=3)
p
ggsave("results/plots/orf_length_sim.pdf", plot = p)



z_plot <- function(fileData, alternative=NULL) {
  
  fileData$sigz <- ifelse(abs(fileData$z) > 1.96, "Significant", "Non-significant")
  cor <- cor.test(fileData$gc, fileData$z, method="spearman")
  sig <- fileData[fileData$sigz == "Significant",]
  
  bin <- binom.test(nrow(fileData[fileData$z > 0,]), nrow(fileData), p=0.5, alternative = alternative)
  bin_sig <- binom.test(nrow(sig[sig$z > 0,]), nrow(sig), alternative = alternative)
  
  p <- ggplot(data=fileData, aes(x=gc, y=z, col=sigz)) +
    geom_point(size=0.6) +
    geom_smooth(method="lm", se=T, lwd=0.5, color="black") +
    geom_hline(yintercept=0, lty=2) +
    scale_color_manual(values=c("grey", "RoyalBlue")) +
    annotate("text", x=0.7, y=27, label=paste("rho = ", round(cor$estimate, 3), ", P = ", cor$p.value), cex=3) +
    annotate("text", x=0.6, y=26, label=paste("one tail binomial test: n successes = ", bin$statistic, ", n trials = ", bin$parameter, ", P = ", bin$p.value), cex=3) +
    annotate("text", x=0.6, y=25, label=paste("if significant z: one tail binomial test: n successes = ", bin_sig$statistic, ", n trials = ", bin_sig$parameter, ", P = ", bin_sig$p.value), cex=3)
  p
  return(p)
}


lengths <- read.csv("results/sim_orf_lengths_zs.csv", head=T)
lengths_plot <- z_plot(lengths, alternative="g")
ggsave("results/plots/lincRNA_orf_length_sim.pdf", plot=lengths_plot, width=10, height=7)

lincRNA <- read.csv("results/sim_stop_count_zs.csv", head=T)
lincRNAs <- z_plot(lincRNA, alternative = "l")
ggsave("results/plots/lincRNA_stop_count_sim.pdf", plot=lincRNAs, width=10, height=7)

coding <- read.csv("results/sim_stop_count_coding_exons_zs.csv", head=T)
coding_exons <- z_plot(coding, alternative="l")
ggsave("results/plots/coding_exons_stop_count_sim.pdf", plot=coding_exons, width=10, height=7)

non_coding <- read.csv("results/sim_stop_count_non_coding_exons_zs.csv", head=T)
non_coding_exons <- z_plot(non_coding, alternative="l")
non_coding_exons
ggsave("results/plots/non_coding_exons_stop_count_sim.pdf", plot=non_coding_exons, width=10, height=7)

library(ggplot2)
library(gridExtra)

file <- read.csv("results/lincRNA_expression_zs.csv", head=T)
tissues <- colnames(file)[3:(length(colnames(file)))]

tissue_plot <- function(tissue) {
  exists <- file[file[[tissue]] >= 0,]
  qnt <- quantile(exists[[tissue]], 0.95)
  exists <- exists[exists[[tissue]] <= qnt,]
  cor <- cor.test(exists$z, exists[[tissue]], method="spearman")
  
  p <- ggplot(data=exists, aes(x=z, y=exists[[tissue]])) +
    geom_point(size=0.4) +
    geom_smooth(method="lm", se=T, lwd=0.5, color="red") + 
    geom_vline(xintercept=0, lty=2) +
    ylab(tissue) + 
    annotate("text", x=Inf, y=max(exists[[tissue]]), label=paste("rho: ", round(cor$estimate, digits=3)), hjust=1, vjust=1) +
    annotate("text", x=Inf, y=max(exists[[tissue]]), label=paste("p: ", round(cor$p.value, digits=3)), hjust=1, vjust=3)
  return(p)
}

plots = lapply(tissues, function(i) tissue_plot(i) )
all_plots <- grid.arrange(grobs = plots, ncol = 5) ## display plot
ggsave("results/plots/expression.pdf", plot = all_plots, width=12, height=10)

file <- read.csv("results/sim_stop_count_mm_zs.csv", head=T)
head(file)
plot(file$gc, file$z)
nrow(file[file$z <0,])
nrow(file)
binom.test(nrow(file[file$z <0,]), nrow(file), alternative = "g")




#######
# Motif
#######

library(ggplot2)
library(ggpubr)

int3 <- read.csv("results2/motif_simulations/sim_int3_stop_count.csv", head=T)
RESCUE <- read.csv("results2/motif_simulations/sim_RESCUE_stop_count.csv", head=T)
ke <- read.csv("results2/motif_simulations/sim_ke400_stop_count.csv", head=T)

ise <- read.csv("results2/motif_simulations/sim_ises_wang_stop_count.csv", head=T)
ess<- read.csv("results2/motif_simulations/sim_ess_fas_hex2_stop_count.csv", head=T)
ess2<- read.csv("results2/motif_simulations/sim_ess_fas_hex3_stop_count.csv", head=T)

rbp_motifs <- read.csv("results/sim_rbp_motifs_stop_count.csv", head=T)

hist_plot <- function(data, xlab, ylab) {
  sims <- data[data$sim_no != "real",]
  real <- data[data$sim_no == "real",]
  min <- min(sims$stop_count, real$stop_count) - 5
  max <- max(sims$stop_count, real$stop_count) + 5
  p <- ggplot(sims, aes(x=stop_count)) +
    geom_histogram(binwidth = 1, fill="RoyalBlue", color="black") + 
    scale_x_continuous(limits = c(max(min), max)) +
    geom_vline(xintercept = real$stop_count, lty=2, color="red") +
    labs(x=xlab, y=ylab)
  return(p)
}

emperical_p <- function(data, column, side="less") {
  sims <- data[data$sim_no != "real",]
  real <- data[data$sim_no == "real",]
  if(side == "less") {
    count <- nrow(sims[sims[[column]] <= real[[column]], ])
  } else {
    count <- nrow(sims[sims[[column]] >= real[[column]], ])
  }
  p <- (count + 1) / (nrow(sims) + 1)
  return(p)
}


int3_hist <- hist_plot(int3, "GC", "Stop codon count")
int3_plot <- ggarrange(int3_hist, ncol = 1, nrow = 1, labels = c("INT3"))
ggsave("results2/plots/motif_simulation_int3.pdf", plot=int3_plot, width=10, height=7)
emperical_p(int3, "stop_count")

RESCUE_hist <- hist_plot(RESCUE, "GC", "Stop codon count")
ke_hist <- hist_plot(ke, "GC", "Stop codon count")
other_motifs_plot <- ggarrange(RESCUE_hist, ke_hist, ncol = 2, nrow = 1, labels = c("RESCUE", "Ke400"))
ggsave("results2/plots/motif_simulation_other_eses.pdf", plot=other_motifs_plot, width=14, height=7)
emperical_p(RESCUE, "stop_count")
emperical_p(ke, "stop_count")

rbp_motifs_hist <- hist_plot(rbp_motifs, "GC", "Stop codon count")
rbp_motifs_hist_plot <- ggarrange(rbp_motifs_hist, ncol = 1, nrow = 1, labels = c("RBP Motifs"))
ggsave("results2/plots/motif_simulation_rbp_motifs.pdf", plot=rbp_motifs_hist_plot, width=10, height=7)
emperical_p(rbp_motifs, "stop_count")rbp_motifs_hist <- hist_plot(rbp_motifs, "GC", "Stop codon count")

ise_hist <- hist_plot(ise, "GC", "Stop codon count")
ise_hist_plot <- ggarrange(ise_hist, ncol = 1, nrow = 1, labels = c("ISEs"))
ggsave("results2/plots/motif_simulation_ise.pdf", plot=ise_hist_plot, width=10, height=7)
emperical_p(ise, "stop_count")

ess_hist <- hist_plot(ess, "GC", "Stop codon count")
ess2_hist <- hist_plot(ess2, "GC", "Stop codon count")
ess_plot <- ggarrange(RESCUE_hist, ke_hist, ncol = 2, nrow = 1, labels = c("ESS", "ESS2"))
ggsave("results2/plots/motif_simulation_ess.pdf", plot=ess_plot, width=14, height=7)
emperical_p(RESCUE, "stop_count")
emperical_p(ke, "stop_count")


######
# Motif density
######

library(ggplot2)
library(stringr)
library(stringi)

file <- read.csv("results3/motif_densities/int3.csv", head=T)


motif_density_plot <- function(data) {
  stops = file[file$motifs == "TAA_TAG_TGA",]
  others = file[file$motifs != "TAA_TAG_TGA",]
  p <- ggplot(data = others, aes(x=id, y=density)) +
    geom_point(size=0.6, color="#555555")
  p <- p + geom_hline(yintercept=stops$density, lty=2, col="#555555")
  p <- p + geom_point(aes(x=stops$id, y=stops$density), size=1.8, color="red")
  return(p)
}



ggsave("results3/plots/int3_motif_densities.pdf", width=10, height=7)




less_than <- nrow(others[others$density > stops$density,])
binom.test(less_than, nrow(others), alternative="g")


file <- read.csv("results3/int3_sim_motif_densities.csv", head=T)

density_boxplot_by_purine <- function(data) {
  colnames(file) <- c("id", "motifs", "density")
  file$purine = round((str_count(file$motifs, "A") + str_count(file$motifs, "G")) / (stri_length(file$motifs) - str_count(file$motifs, "_")), 4)
  stops_purine = file$purine[file$motifs=="TAA_TAG_TGA"]
  stops_density =  file$density[file$motifs=="TAA_TAG_TGA"]
  file$purine <- as.factor(file$purine)
  file$col <- ifelse(file$purine == stops_purine, "a", "b")
  p <- ggplot(file, aes(x=file$purine, y=file$density)) +
    stat_boxplot(geom ='errorbar') + 
    geom_boxplot(aes(fill=file$col)) +
    scale_fill_manual(values=c("RoyalBlue", "#dddddd")) + 
    geom_hline(yintercept = stops_density, lty=2) + 
    labs(x = "Query codons purine content", y="Query codons ND") + 
    annotate("text", x=min(as.numeric(file$purine)), hjust=0.2, y=stops_density + 0.05, label="TAA,TAG,TGA", cex=3) +
    annotate("text", x=min(as.numeric(file$purine)), hjust=0.1, y=stops_density - 0.05, label=round(stops_density,4), cex=3) + 
    theme(legend.position="none")
  return(p)
}

int3_file <- read.csv("results3/int3_sim_motif_densities.csv", head=T)
int3 <- density_boxplot_by_purine(int3_file)
int3
ggsave("results3/plots/int3_nd_purine_content.pdf", width=10, height=7)

colnames(file) <- c("id", "motifs", "density")
file$purine = round((str_count(file$motifs, "A") + str_count(file$motifs, "G")) / (stri_length(file$motifs) - str_count(file$motifs, "_")), 4)

stops = file[file$motifs=="TAA_TAG_TGA",]
purines = file[file$purine == stops$purine,]
binom.test(nrow(purines[purines$density <= stops$density,]), nrow(purines), alternative = "l")


plot <- motif_density_plot(file)
ggsave("results3/plots/int3_motif_nd.pdf", plot=plot, width=10, height=7)



binom.test(nrow(others[others$density <= stops$density,]), nrow(others), alternative = "l")


stops$py <- (str_count(stops$motifs, "A") + str_count(stops$motifs, "G")) / (stri_length(stops$motifs) - str_count(stops$motifs, "_"))
others$py = (str_count(others$motifs, "A") + str_count(others$motifs, "G")) / (stri_length(others$motifs) - str_count(others$motifs, "_"))
same <- others[others$py == stops$py,]
plot(same$id, same$density)
points(stops$id, stops$density)
nrow(same[same$density <= stops$density,])
nrow(same)


###### Motif Stuff 


library(ggplot2)
library(stringr)
library(stringi)

file <- read.csv("results3/motif_densities/int3.csv", head=T)


motif_density_plot <- function(data) {
  stops = file[file$motifs == "TAA_TAG_TGA",]
  others = file[file$motifs != "TAA_TAG_TGA",]
  p <- ggplot(data = others, aes(x=id, y=density)) +
    geom_point(size=0.6, color="#555555")
  p <- p + geom_hline(yintercept=stops$density, lty=2, col="#555555")
  p <- p + geom_point(aes(x=stops$id, y=stops$density), size=1.8, color="red")
  return(p)
}

ggsave("results3/plots/int3_motif_densities.pdf", width=10, height=7)




less_than <- nrow(others[others$density > stops$density,])
binom.test(less_than, nrow(others), alternative="g")

get_n <- function(x, vjust=0){
  data = data.frame(y = max(x)+vjust,label = paste("N = ", length(x), sep=""))
  return(data)
}

grouped_boxplot <- function(data, x, y, vjust, xlab, ylab) {
  p <- ggplot(data, aes(x=data[[x]], y=data[[y]])) + geom_boxplot(fill="RoyalBlue")
  p <- p + stat_summary(fun.data = get_n, fun.args = list("vjust" = vjust), geom = "text", aes(group=x))
  p <- p + geom_hline(yintercept=0, lty=2, color="red")
  p <- p + labs(x=xlab, y=ylab)
  return(p)
}




density_boxplot_by_purine <- function(data) {
  colnames(data) <- c("id", "motifs", "density")
  data$purine = round((str_count(data$motifs, "A") + str_count(data$motifs, "G")) / (stri_length(data$motifs) - str_count(data$motifs, "_")), 4)
  stops_purine = data$purine[data$motifs=="TAA_TAG_TGA"]
  stops_density =  data$density[data$motifs=="TAA_TAG_TGA"]
  data$purine <- as.factor(data$purine)
  data$col <- ifelse(data$purine == stops_purine, "a", "b")
  p <- ggplot(data, aes(x=data$purine, y=data$density)) +
    stat_boxplot(geom ='errorbar') + 
    geom_boxplot(aes(fill=data$col)) +
    scale_fill_manual(values=c("RoyalBlue", "#dddddd")) + 
    geom_hline(yintercept = stops_density, lty=2) + 
    labs(x = "Query codons purine content", y="Query codons ND") + 
    annotate("text", x=min(as.numeric(data$purine)), hjust=0.2, y=stops_density + 0.05, label="TAA,TAG,TGA", cex=3) +
    annotate("text", x=min(as.numeric(data$purine)), hjust=0.1, y=stops_density - 0.05, label=round(stops_density,4), cex=3) + 
    theme(legend.position="none")
  return(p)
}

purine_binom_test <- function(data) {
  colnames(data) <- c("id", "motifs", "density")
  data$purine = round((str_count(data$motifs, "A") + str_count(data$motifs, "G")) / (stri_length(data$motifs) - str_count(data$motifs, "_")), 4)
  stops_purine = data$purine[data$motifs=="TAA_TAG_TGA"]
  stops_density =  data$density[data$motifs=="TAA_TAG_TGA"]
  data$purine <- as.factor(data$purine)
  remove_stops<- data[data$motifs != "TAA_TAG_TGA",]
  matched_purine <- remove_stops[remove_stops$purine == stops_purine,]
  binom.test(nrow(matched_purine[matched_purine$density <= stops_density,]), nrow(matched_purine), alternative = "l")
}
purine_binom_test(int3)

int3 <- read.csv("results3/int3_sim_motif_densities.csv", head=T)
head(int3)
int3 <- density_boxplot_by_purine(int3)
int3

ggsave("results3/plots/int3_nd_purine_content.pdf", plot = int3, width=10, height=7)


RESCUE <- read.csv("results3/RESCUE_sim_motif_genome_gc_matched_densities.csv", head=T)
head(RESCUE)
RESCUE <- density_boxplot_by_purine(RESCUE)
RESCUE

ggsave("results3/plots/RESCUE_nd_purine_content.pdf", plot = RESCUE, width=10, height=7)


ises <- read.csv("results3/ises_wang_sim_motif_densities.csv", head=T)
head(ises)
ises <- density_boxplot_by_purine(ises)
ises
ggsave("results3/plots/ises_nd_purine_content.pdf", plot = ises, width=10, height=7)

ises <- read.csv("results3/ises_wang_sim_motif_densities.csv", head=T)
head(ises)
ises <- density_boxplot_by_purine(ises)
ises
ggsave("results3/plots/ises_nd_purine_content.pdf", plot = ises, width=10, height=7)


colnames(file) <- c("id", "motifs", "density")
file$purine = round((str_count(file$motifs, "A") + str_count(file$motifs, "G")) / (stri_length(file$motifs) - str_count(file$motifs, "_")), 4)

stops = file[file$motifs == "TAA_TAG_TGA",]
purines = file[file$purine == stops$purine,]
binom.test(nrow(purines[purines$density <= stops$density,]), nrow(purines), alternative = "l")


plot <- motif_density_plot(file)
ggsave("results3/plots/int3_motif_nd.pdf", plot=plot, width=10, height=7)



binom.test(nrow(others[others$density <= stops$density,]), nrow(others), alternative = "l")


stops$py <- (str_count(stops$motifs, "A") + str_count(stops$motifs, "G")) / (stri_length(stops$motifs) - str_count(stops$motifs, "_"))
others$py = (str_count(others$motifs, "A") + str_count(others$motifs, "G")) / (stri_length(others$motifs) - str_count(others$motifs, "_"))
same <- others[others$py == stops$py,]
plot(same$id, same$density)
points(stops$id, stops$density)
nrow(same[same$density <= stops$density,])
nrow(same)

