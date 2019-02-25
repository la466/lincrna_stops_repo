library(ggplot2)
library(gridExtra)
library(ggpubr)

intron_ese_plot <- function(data, title = NULL) {
  p = ggplot(data, aes(x = log10(data$intron_size), y = data$ese_density)) + 
    geom_point(size = 0.8) + 
    geom_smooth(method='lm', col = "black") + 
    scale_y_continuous(limits = c(0, 1)) +
    labs(x = "log10 intron size", y = "ESE density", title = title) +
    theme(plot.title = element_text(hjust = 0.5))
  return (p)
}

# ese_sets = c("int3", "RESCUE", "ke400", "PESE", "ESR")
ese_sets = c("int3", "RESCUE", "PESE", "ESR", "ke400", "combined_eses")

cors_family = data.frame(motif_set = character(), lincrna_rho = double(), lincrna_p = double(), protein_coding_rho = double(), protein_coding_p = double(), mean_lincrna = double(), mean_pc = double(), wilcox_text_score = double(), wilcox_test_p = double())
cors_all = data.frame(motif_set = character(), lincrna_rho = double(), lincrna_p = double(), protein_coding_rho = double(), protein_coding_p = double(), mean_lincrna = double(), mean_pc = double(), wilcox_text_score = double(), wilcox_test_p = double())
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

  lincrna = read.csv(paste("clean_run/tests/ese_densities/", set_name, "_lincrna_ese_densities_all_sequences.csv", sep = ""), head = T)
  pc = read.csv(paste("clean_run/tests/ese_densities/", set_name, "_pc_ese_densities_all_sequences.csv", sep = ""), head = T)
  cor1 = cor.test(log10(lincrna$intron_size), lincrna$ese_density, method = "spearman")
  cor2 = cor.test(log10(pc$intron_size), pc$ese_density, method = "spearman")
  w_test = wilcox.test(lincrna$ese_density, pc$ese_density)
  cor_output = data.frame(motif_set = paste(set_name, "_all_sequences", sep = ""), lincrna_rho = cor1$estimate, lincrna_p = cor1$p.value, protein_coding_rho = cor2$estimate, protein_coding_p = cor2$p.value, mean_lincrna = mean(lincrna$ese_density), mean_pc = mean(pc$ese_density), wilcox_text_score = w_test$statistic, wilcox_test_p = w_test$p.value)
  cors_all = rbind(cors_all, cor_output)

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
cors_all$adjusted_p = p.adjust(cors_all$wilcox_test_p, method = "bonferroni")
cors_flanks$adjusted_p = p.adjust(cors_flanks$wilcox_test_p, method = "bonferroni")

cors_family[nrow(cors_family)+1,] <- NA
cors_all[nrow(cors_all)+1,] <- NA
cors = rbind(cors_family, cors_all, cors_flanks)
write.csv(cors, file = "clean_run/tests/ese_densities/intron_size_ese_density_correlations.csv", row.names = F, na = "")




lincrna = read.csv(paste("clean_run/tests/ese_densities/", "combined_eses", "_lincrna_ese_densities_flanks.csv", sep = ""), head = T)
pc = read.csv(paste("clean_run/tests/ese_densities/", "combined_eses", "_pc_ese_densities_flanks.csv", sep = ""), head = T)

lincrna$ratio = lincrna$stop_ese_density / lincrna$non_stop_ese_density
plot(log10(lincrna$intron_size), lincrna$ratio)
cor.test(log10(lincrna$intron_size), lincrna$ratio, method = "spearman")

pc$ratio = pc$stop_ese_density / pc$non_stop_ese_density
plot(log10(pc$intron_size), pc$ratio)
cor.test(log10(pc$intron_size), pc$ratio, method = "spearman")

w = wilcox.test(lincrna$ese_density, pc$ese_density)
w
