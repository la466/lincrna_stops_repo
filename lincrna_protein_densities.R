library(ggplot2)
library(gridExtra)
library(ggpubr)
library(reshape2)

intron_ese_plot <- function(data, title = NULL) {
  p = ggplot(data, aes(x = log10(data$intron_size), y = data$ese_density)) + 
    geom_point(size = 0.8) + 
    geom_smooth(method='lm', col = "black") + 
    scale_y_continuous(limits = c(0, 1)) +
    labs(x = "log10 intron size", y = "ESE density", title = title) +
    theme(plot.title = element_text(hjust = 0.5))
  return (p)
}

ese_sets = c("int3", "RESCUE", "PESE", "ESR",  "combined_eses")

cors_family = data.frame(motif_set = character(), lincrna_rho = double(), lincrna_p = double(), protein_coding_rho = double(), protein_coding_p = double(), mean_lincrna = double(), mean_pc = double(), wilcox_text_score = double(), wilcox_test_p = double())
cors_flanks = data.frame(motif_set = character(), lincrna_rho = double(), lincrna_p = double(), protein_coding_rho = double(), protein_coding_p = double(), mean_lincrna = double(), mean_pc = double(), wilcox_text_score = double(), wilcox_test_p = double())

for (set_name in ese_sets) {
  print(set_name)
  
  lincrna = read.csv(paste("clean_run/tests/ese_densities/", set_name, "_lincrna_ese_densities.csv", sep = ""), head = T)
  pc = read.csv(paste("clean_run/tests/ese_densities/", set_name, "_pc_ese_densities.csv", sep = ""), head = T)
  cor1 = cor.test(log10(lincrna$intron_size), lincrna$ese_density, method = "spearman")
  cor2 = cor.test(log10(pc$intron_size), pc$ese_density, method = "spearman")
  w_test = wilcox.test(lincrna$ese_density, pc$ese_density)
  cor_output = data.frame(motif_set = set_name, lincrna_rho = cor1$estimate, lincrna_p = cor1$p.value, protein_coding_rho = cor2$estimate, protein_coding_p = cor2$p.value, mean_lincrna = mean(lincrna$ese_density), mean_pc = mean(pc$ese_density), wilcox_text_score = w_test$statistic, wilcox_test_p = w_test$p.value)
  cors_family = rbind(cors_family, cor_output)

  lincrna = read.csv(paste("clean_run/tests/ese_densities/", set_name, "_lincrna_ese_densities_flanks.csv", sep = ""), head = T)
  pc = read.csv(paste("clean_run/tests/ese_densities/", set_name, "_pc_ese_densities_flanks.csv", sep = ""), head = T)
  cor1 = cor.test(log10(lincrna$intron_size), lincrna$ese_density, method = "spearman")
  cor2 = cor.test(log10(pc$intron_size), pc$ese_density, method = "spearman")
  w_test = wilcox.test(lincrna$ese_density, pc$ese_density)
  cor_output = data.frame(motif_set = paste(set_name, "_flanks", sep = ""), lincrna_rho = cor1$estimate, lincrna_p = cor1$p.value, protein_coding_rho = cor2$estimate, protein_coding_p = cor2$p.value, mean_lincrna = mean(lincrna$ese_density), mean_pc = mean(pc$ese_density), wilcox_text_score = w_test$statistic, wilcox_test_p = w_test$p.value)
  cors_flanks = rbind(cors_flanks, cor_output)
  
  plot = ggarrange(
    intron_ese_plot(lincrna, "LincRNA"),
    intron_ese_plot(pc, "Protein-coding coding exons"),
    ncol = 2,
    labels = c("A", "B")
  )
  ggsave(plot = plot, file = paste("clean_run/plots/intron_size_ese_density_flanks_", set_name, ".pdf", sep = ""), width = 10, height = 5)
  
}

cors_family$adjusted_p = p.adjust(cors_family$wilcox_test_p, method = "bonferroni")
cors_flanks$adjusted_p = p.adjust(cors_flanks$wilcox_test_p, method = "bonferroni")

cors_family[nrow(cors_family)+1,] <- NA
cors = rbind(cors_family, cors_flanks)
write.csv(cors, file = "clean_run/tests/ese_densities/intron_size_ese_density_correlations.csv", row.names = F, na = "")



stop_density_plot <- function(data, title = NULL) {
  data.melt = melt(data, id.vars = c("id", "intron_size"), measure.vars = c("stop_ese_density", "non_stop_ese_density"), value.name = "density")
  legend_labs = labels = c("Stop ESEs", "Non stop ESEs")
  ggplot(data.melt, aes(x = log10(intron_size), y = density, group = variable)) +
    geom_point(cex = 1, aes(colour = variable, shape = variable)) +
    geom_smooth(method = 'lm', aes(col = variable), show.legend = F) +
    scale_color_manual(name = "", labels = legend_labs, values = c("RoyalBlue", "black")) +
    scale_shape_manual(name = "", labels = legend_labs, values = c(4, 16)) +
    scale_y_continuous(limits = c(0, 0.6)) +
    labs(x = "log10 intron size", y = "ESE density", title = title) +
    theme(
      plot.title = element_text(hjust = 0.5),
      legend.position = c(0.85,0.9),
      legend.title = element_blank(),
      legend.background = element_blank(),
      legend.key = element_blank(),
    ) +
    guides(shape = guide_legend(override.aes = list(size = 2)))
}

fischers_z <- function(x, y1, y2) {
  # create linear models for both
  lm1 = lm(y1 ~ log10(x))
  lm2 = lm(y2 ~ log10(x))
  # get summaries
  summary1 <- summary(lm1)
  summary2 <- summary(lm2)
  # get coefficients
  coeffs1 <- summary1$coefficients
  coeffs2 <- summary2$coefficients
  # get slopes
  slope1 <- coeffs1[2,1]
  slope2 <- coeffs2[2,1]
  # get std error
  se1 <- coeffs1[2,2]
  se2 <- coeffs2[2,2]
  # calculate z
  z <- abs((slope1 - slope2)/(sqrt(se1^2 + se2^2)))
  return(z)
}


p_from_z <- function(z) {
   #  calculate p value
  pval <- 2*pnorm(z, lower.tail=FALSE)
  return(pval)
}



lincrna = read.csv(paste("clean_run/tests/ese_densities/", "INT3", "_lincrna_ese_densities.csv", sep = ""), head = T)
pc = read.csv(paste("clean_run/tests/ese_densities/", "INT3", "_pc_ese_densities.csv", sep = ""), head = T)


plot = ggarrange(
  intron_ese_plot(lincrna, title = "LincRNA"),
  intron_ese_plot(pc, title = "Protein-coding coding exons"),
  stop_density_plot(lincrna, title = "LincRNA"),
  stop_density_plot(pc, title = "Protein-coding coding exons"),
  ncol = 2,
  nrow = 2,
  labels = c("A", "B", "C", "D")
)
ggsave(plot = plot, file = "clean_run/plots/intron_size_ese_density_plot.pdf", width = 12, height = 10)



ese_sets = c("int3", "RESCUE", "PESE", "ESR", "combined_eses")
# ese_sets = c("int3", "RESCUE")

correlations = data.frame(
    motif_set = character(),
    lincrna_stop_ese_rho = double(),
    lincrna_stop_ese_p = double(),
    lincrna_non_stop_ese_rho = double(),
    lincrna_non_stop_ese_p = double(),
    lincrna_fischers_z = double(),
    lincrna_diff_p = double(),
    pc_stop_ese_rho = double(),
    pc_stop_ese_p = double(),
    pc_non_stop_ese_rho = double(),
    pc_non_stop_ese_p = double(),
    pc_fischers_z = double(),
    pc_diff_p = double()
  )
correlations_flanks = data.frame(
  motif_set = character(),
  lincrna_stop_ese_rho = double(),
  lincrna_stop_ese_p = double(),
  lincrna_non_stop_ese_rho = double(),
  lincrna_non_stop_ese_p = double(),
  fischers_z = double(),
  diff_p = double(),
  pc_stop_ese_rho = double(),
  pc_stop_ese_p = double(),
  pc_non_stop_ese_rho = double(),
  pc_non_stop_ese_p = double(),
  fischers_z = double(),
  diff_p = double()
  )

types = c("", "_flanks")

for (type in types) {
  for (set_name in ese_sets) {
    print(set_name)
    lincrna = read.csv(paste("clean_run/tests/ese_densities/", set_name, "_lincrna_ese_densities", type, ".csv", sep = ""), head = T)
    pc = read.csv(paste("clean_run/tests/ese_densities/", set_name, "_pc_ese_densities", type, ".csv", sep = ""), head = T)
    
    lincra = lincrna[log10(lincrna$intron_size) != 0,]
    pc = pc[log10(pc$intron_size) != 0,]
    
    cor1 = cor.test(log10(lincrna$intron_size), lincrna$stop_ese_density, method = "spearman")
    cor2 = cor.test(log10(lincrna$intron_size), lincrna$non_stop_ese_density, method = "spearman")
    cor3 = cor.test(log10(pc$intron_size), pc$stop_ese_density, method = "spearman")
    cor4 = cor.test(log10(pc$intron_size), pc$non_stop_ese_density, method = "spearman")
    
    z1 = fischers_z(log10(lincrna$intron_size), lincrna$stop_ese_density, lincrna$non_stop_ese_density)
    z2 = fischers_z(log10(pc$intron_size), pc$stop_ese_density, pc$non_stop_ese_density)

    cor_output = data.frame(
      motif_set = paste(set_name, type, sep = ""),
      lincrna_stop_ese_rho = cor1$estimate,
      lincrna_stop_ese_p = cor1$p.value,
      lincrna_non_stop_ese_rho = cor2$estimate,
      lincrna_non_stop_ese_p = cor2$p.value,
      lincrna_fischers_z = z1,
      lincrna_diff_p = p_from_z(z1),
      pc_stop_ese_rho = cor3$estimate,
      pc_stop_ese_p = cor3$p.value,
      pc_non_stop_ese_rho = cor4$estimate,
      pc_non_stop_ese_p = cor4$p.value,
      pc_fischers_z = z2,
      pc_diff_p = p_from_z(z2)
    )
    if (type == "") {
      correlations = rbind(correlations, cor_output)
    } else if (type == "_all_sequences") {
      correlations_all = rbind(correlations_all, cor_output)
    } else if (type == "_flanks") {
      correlations_flanks = rbind(correlations_flanks, cor_output)
    }
  }
}

correlations[nrow(correlations)+1,] <- NA
cors = rbind(correlations, correlations_flanks)
write.csv(cors, file = "clean_run/tests/ese_densities/intron_size_ese_density_stop_non_stop_correlations.csv", row.names = F, na = "")


cor_outputs = data.frame(
  motif_set = character(),
  ese_rho = double(),
  ese_p = double(),
  stop_ese_rho = double(),
  stop_ese_p = double(),
  non_stop_ese_rho = double(),
  non_stop_ese_p = double(),
  fischers_z = double(),
  diff_p = double()
)

for (type in types) {
  data = read.csv(paste("clean_run/tests/ese_densities/PESE_non_coding_ese_densities", type, ".csv", sep = ""), head = T)
  data = data[log10(data$intron_size) != 0,]
  cor = cor.test(log10(data$intron_size), data$ese_density, method = "spearman")
  cor1 = cor.test(log10(data$intron_size), data$stop_ese_density, method = "spearman")
  cor2 = cor.test(log10(data$intron_size), data$non_stop_ese_density, method = "spearman")

  z1 = fischers_z(log10(data$intron_size), data$stop_ese_density, data$non_stop_ese_density)
  
  cor_output = data.frame(
    motif_set = paste("PESE", type, sep = ""),
    ese_rho = cor$estimate,
    ese_p = cor$p.value,
    stop_ese_rho = cor1$estimate,
    stop_ese_p = cor1$p.value,
    non_stop_ese_rho = cor2$estimate,
    non_stop_ese_p = cor2$p.value,
    fischers_z = z1,
    diff_p = p_from_z(z1)
  )
  cor_outputs = rbind(cor_outputs, cor_output)
}

write.csv(cor_outputs, file = "clean_run/tests/ese_densities/intron_size_ese_density_stop_non_stop_correlations_non_coding.csv", row.names = F, na = "")

