library(ggplot2)

get_n <- function(x, vjust=0){
  data = data.frame(y = max(x)+vjust,label = paste("n = ", length(x), sep=""))
  return(data)
}

density_scatterplot <- function(data, ycol = "density", xlab = "", ylab = "") {
  data$id = as.numeric(seq(1:nrow(data)))
  stops = data[data$codons == "TAA_TAG_TGA",]
  not_stops = data[data$codons != "TAA_TAG_TGA",]
  p <- ggplot() + 
    geom_point(aes(x = not_stops$id, y = not_stops$nd), cex = 0.5, pch = 16) + 
    geom_point(aes(x = stops$id, y = stops$nd), cex = 3, pch = 16, col = "red") + 
    # scale_x_continuous(labels = not_stops$codons, breaks = seq(1:nrow(not_stops))) +
    theme(axis.text.x = element_blank()) + 
    geom_hline(yintercept = 0, lty = 2) + 
    labs(x = xlab, y = ylab)
    
  return(p)
}

density_boxplot <- function(data, ycol = "density") {
  stops_purine = data$purine[data$codons=="TAA_TAG_TGA"]
  stops_density =  data[[ycol]][data$codons=="TAA_TAG_TGA"]
  data$purine <- as.factor(data$purine)
  data$col <- ifelse(data$purine == stops_purine, "a", "b")
  p <- ggplot(data, aes(x=data$purine, y=data[[ycol]])) +
    stat_boxplot(geom ='errorbar') + 
    geom_boxplot(aes(fill=data$col)) +
    scale_fill_manual(values=c("RoyalBlue", "#dddddd")) + 
    geom_hline(yintercept = stops_density, lty=2) + 
    labs(x = "Query codons purine content", y="Query codons ND") + 
    annotate("text", x=min(as.numeric(data$purine)), hjust=0.2, y=stops_density + 0.05, label="TAA,TAG,TGA", cex=3) +
    annotate("text", x=min(as.numeric(data$purine)), hjust=0.1, y=stops_density - 0.05, label=round(stops_density,4), cex=3) + 
    theme(legend.position="none") + 
    scale_x_discrete(labels = c("0-0.1", "0.1-0.2", "0.2-0.3", "0.3-0.4", "0.4-0.5", "0.5-0.6", "0.6-0.7", "0.7-0.8", "0.8-0.9", "0.9-1", "1")) + 
    stat_summary(fun.data = get_n, fun.args = list("vjust" = 0.1), geom = "text", aes(group = "purine"), size=3)
  return(p)
}

binom_test <- function(data, ycol = "density", group = NULL) {
  stops = data[data$codons == "TAA_TAG_TGA",]
  not_stops = data[data$codons != "TAA_TAG_TGA",]
  if (!is.null(group)) {
    not_stops = not_stops[not_stops[[group]] == stops[[group]],]
  }
  b = binom.test(nrow(not_stops[not_stops[[ycol]] <= stops[[ycol]], ]), nrow(not_stops), alternative = "l")
  return(b)
}

filepath = "clean_run/motif_tests/int3_densities.csv"
file = read.csv(filepath, head = T)

head(file)
file$gc = substr(file$gc_content, 0, 3)
file$purine = substr(file$purine_content, 0, 3)
stops =  file[file$codons == "TAA_TAG_TGA",]
gc_matched = file[file$gc == stops$gc,]



plot1 = density_scatterplot(gc_matched, ycol = "nd", xlab = "GC matched codon sets (sorted alphabetically)", ylab = "ND")
plot2 = density_boxplot(gc_matched, ycol = "nd")
plot1
plot2
binom_test(gc_matched, ycol = "nd")
binom_test(gc_matched, ycol = "nd", group = "purine")

cairo_pdf(filename = "clean_run/plots/gc_matched_codons_nd.pdf", width = 10)
plot1
dev.off()
cairo_pdf("clean_run/plots/gc_matched_codons_purine_nd.pdf", width = 10)
plot2
dev.off()


filepath = "clean_run/motif_tests/ises_wang_densities.csv"
file = read.csv(filepath, head = T)
file$gc = signif(file$gc_content, 1)
file$purine = signif(file$purine_content, 1)
stops =  file[file$codons == "TAA_TAG_TGA",]
gc_matched = file[file$gc == stops$gc,]


density_scatterplot(gc_matched, ycol = "nd", xlab = "", ylab = "")
binom_test(gc_matched, ycol = "nd")
density_boxplot(gc_matched, ycol = "nd")
binom_test(gc_matched, ycol = "nd", group = "purine")

ggsave("clean_run/plots/gc_matched_codons_purine_nd.pdf", plot = plot1)

not_stops_purine = gc_matched[gc_matched$codons != "TAA_TAG_TGA" & gc_matched$purine == stops$purine,]
nrow(not_stops_purine[not_stops_purine$nd <= -0.0426,])




density_histogram<- function(data, xlab = "", ylab = "", binwidth = 0.01) {
  real = data[data$sim_id == "real",]
  sims = data[data$sim_id != "real",]
  p <- ggplot() + 
    geom_histogram(aes(sims$stop_density), binwidth = binwidth, fill = "RoyalBlue", col = "black") +
    geom_vline(xintercept = real$stop_density, lty = 1, size = 1.5, col = "red") +
    labs(x = xlab, y = ylab)
  return(p)
}

emperical_p <- function(data, id_col, id_differentiator, test_col) {
  query_val = data[data[[id_col]] == id_differentiator,]
  test_vals = data[data[[id_col]] != id_differentiator,]
  p <- (nrow(test_vals[test_vals[[test_col]] <= query_val[[test_col]],]) + 1) / (nrow(test_vals) + 1)
  return(p)
}

real_density <- function(data, id_col, id_differentiator, test_col) {
  real = data[data[[id_col]] == id_differentiator,]
  return (real[[test_col]])
}

normalised_density <- function(data, id_col, id_differentiator, test_col) {
  real = data[data[[id_col]] == id_differentiator,]
  sims = data[data[[id_col]] != id_differentiator,]
  nd = (real[[test_col]] - mean(sims[[test_col]])) / mean(sims[[test_col]])
  return (nd)
}

filepath = "clean_run/motif_tests/int3_stop_codon_densities.csv"
data = read.csv(filepath, head = T)
head(data)
plot <- density_histogram(data, xlab = "Stop codon density", ylab = "Count", binwidth = 0.01)
ggsave("clean_run/plots/dinucleotide_matched_stop_densities_int3.pdf", plot = plot)
stop_density = real_density(data, "sim_id", "real", "stop_density")
normalised_density(data, "sim_id", "real", "stop_density")
emperical_p(data, "sim_id", "real", "stop_density")

filepath = "clean_run/motif_tests/RESCUE_stop_codon_densities.csv"
data = read.csv(filepath, head = T)
plot <- density_histogram(data, xlab = "Stop codon density", ylab = "Count", binwidth = 0.01)
ggsave("clean_run/plots/dinucleotide_matched_stop_densities_RESCUE.pdf", plot = plot)
stop_density = real_density(data, "sim_id", "real", "stop_density")
stop_density
normalised_density(data, "sim_id", "real", "stop_density")
emperical_p(data, "sim_id", "real", "stop_density")

filepath = "clean_run/motif_tests/ke400_stop_codon_densities.csv"
data = read.csv(filepath, head = T)
plot <- density_histogram(data, xlab = "Stop codon density", ylab = "Count", binwidth = 0.01)
ggsave("clean_run/plots/dinucleotide_matched_stop_densities_ke400.pdf", plot = plot)
stop_density = real_density(data, "sim_id", "real", "stop_density")
stop_density
normalised_density(data, "sim_id", "real", "stop_density")
emperical_p(data, "sim_id", "real", "stop_density")

filepath = "clean_run/motif_tests/PESE_stop_codon_densities.csv"
data = read.csv(filepath, head = T)
plot <- density_histogram(data, xlab = "Stop codon density", ylab = "Count", binwidth = 0.01)
ggsave("clean_run/plots/dinucleotide_matched_stop_densities_PESE.pdf", plot = plot)
stop_density = real_density(data, "sim_id", "real", "stop_density")
stop_density
normalised_density(data, "sim_id", "real", "stop_density")
emperical_p(data, "sim_id", "real", "stop_density")

filepath = "clean_run/motif_tests/ESR_stop_codon_densities.csv"
data = read.csv(filepath, head = T)
plot <- density_histogram(data, xlab = "Stop codon density", ylab = "Count", binwidth = 0.01)
ggsave("clean_run/plots/dinucleotide_matched_stop_densities_ESR.pdf", plot = plot)
stop_density = real_density(data, "sim_id", "real", "stop_density")
stop_density
normalised_density(data, "sim_id", "real", "stop_density")
emperical_p(data, "sim_id", "real", "stop_density")

filepath = "clean_run/motif_tests/ises_wang_stop_codon_densities.csv"
data = read.csv(filepath, head = T)
plot = density_histogram(data, xlab = "Stop codon density", ylab = "Count", binwidth = 0.005)
ggsave("clean_run/plots/dinucleotide_matched_stop_densities_ises_wang.pdf", plot = plot)
stop_density = real_density(data, "sim_id", "real", "stop_density")
stop_density
normalised_density(data, "sim_id", "real", "stop_density")
emperical_p(data, "sim_id", "real", "stop_density")

filepath = "clean_run/motif_tests/ess_fas_hex2_stop_codon_densities.csv"
data = read.csv(filepath, head = T)
plot = density_histogram(data, xlab = "Stop codon density", ylab = "Count", binwidth = 0.005)
ggsave("clean_run/plots/dinucleotide_matched_stop_densities_ess_hex2.pdf", plot = plot)
stop_density = real_density(data, "sim_id", "real", "stop_density")
stop_density
normalised_density(data, "sim_id", "real", "stop_density")
emperical_p(data, "sim_id", "real", "stop_density")

filepath = "clean_run/motif_tests/ess_fas_hex3_stop_codon_densities.csv"
data = read.csv(filepath, head = T)
plot = density_histogram(data, xlab = "Stop codon density", ylab = "Count", binwidth = 0.005)
ggsave("clean_run/plots/dinucleotide_matched_stop_densities_ess_hex3.pdf", plot = plot)
stop_density = real_density(data, "sim_id", "real", "stop_density")
stop_density
normalised_density(data, "sim_id", "real", "stop_density")
emperical_p(data, "sim_id", "real", "stop_density")

filepath = "clean_run/motif_tests/rbp_motifs_cds_stop_codon_densities.csv"
data = read.csv(filepath, head = T)
plot = density_histogram(data, xlab = "Stop codon density", ylab = "Count", binwidth = 0.001)
ggsave("clean_run/plots/dinucleotide_matched_stop_densities_rbp_motifs_cds.pdf", plot = plot)
stop_density = real_density(data, "sim_id", "real", "stop_density")
stop_density
normalised_density(data, "sim_id", "real", "stop_density")
emperical_p(data, "sim_id", "real", "stop_density")

filepath = "clean_run/motif_tests/rbp_motifs_non_cds_stop_codon_densities.csv"
data = read.csv(filepath, head = T)
plot = density_histogram(data, xlab = "Stop codon density", ylab = "Count", binwidth = 0.001)
ggsave("clean_run/plots/dinucleotide_matched_stop_densities_rbp_motifs_non_cds.pdf", plot = plot)
stop_density = real_density(data, "sim_id", "real", "stop_density")
stop_density
normalised_density(data, "sim_id", "real", "stop_density")
emperical_p(data, "sim_id", "real", "stop_density")


filepath = "clean_run/motif_tests/int3_densities.csv"
file = read.csv(filepath, head = T)
file$gc = signif(file$gc_content, 1)
file$purine = signif(file$purine_content, 1)
stops =  file[file$codons == "TAA_TAG_TGA",]
# gc_matched = file[file$gc == stops$gc,]
close = c('TCA_TTA_TTG', 'AAA_TTA_TTG', 'AAG_TTA_TTG', 'AGA_TTA_TTG', 'AAA_TCA_TTA', 'AAG_TCA_TTA', 'AGA_TCA_TTA', 'AAA_AAG_TTA', 'AAA_AGA_TTA', 'AAG_AGA_TTA', 'AAA_TCA_TTG', 'AAG_TCA_TTG', 'AGA_TCA_TTG', 'AAA_AAG_TTG', 'AAA_AGA_TTG', 'AAG_AGA_TTG', 'AAA_AAG_TCA', 'AAA_AGA_TCA', 'AAG_AGA_TCA', 'AAA_AAG_AGA')
close_cases = file[file$codons %in% close,]
nrow(close_cases)

boxplot(close_cases$nd)
abline(h = stops$nd)
nrow(close_cases[close_cases$nd <= stops$nd,]) / nrow(close_cases)

less = close_cases[close_cases$nd <= stops$nd,]
less


less$gc = substr(less$gc_content, 0, 3)
less



## exon dinucleotide controls
library(ggplot2)


motif_set_test = function(filepath) {
  data = read.csv(filepath, head = T)
  stop_density = real_density(data, "sim_id", "real", "stop_density")
  nd = normalised_density(data, "sim_id", "real", "stop_density")
  p = emperical_p(data, "sim_id", "real", "stop_density")
  print(stop_density)
  print(nd)
  print(p)
}
filepath = "clean_run/motif_tests/int3_stop_codon_densities_exon_dinucleotides.csv"
motif_set_test(filepath)
filepath = "clean_run/motif_tests/RESCUE_stop_codon_densities_exon_dinucleotides.csv"
filepath = "clean_run/motif_tests/ess_fas_hex2_stop_codon_densities_exon_dinucleotides.csv"
filepath = "clean_run/motif_tests/ess_fas_hex3_stop_codon_densities_exon_dinucleotides.csv"
filepath = "clean_run/motif_tests/ises_wang_stop_codon_densities_exon_dinucleotides.csv"
filepath = "clean_run/motif_tests/ke400_stop_codon_densities_exon_dinucleotides.csv"
filepath = "clean_run/motif_tests/ESR_stop_codon_densities_exon_dinucleotides.csv"
filepath = "clean_run/motif_tests/PESE_stop_codon_densities_exon_dinucleotides.csv"
filepath = "clean_run/motif_tests/rbp_motifs_non_cds_stop_codon_densities_exon_dinucleotides.csv"
filepath = "clean_run/motif_tests/rbp_motifs_cds_stop_codon_densities_exon_dinucleotides.csv"
motif_set_test(filepath)


