library(ggplot2)
library(ggpubr)
library(gridExtra)

density_histogram<- function(data, xlab = "", ylab = "", binwidth = 0.01, title = "") {
  real = data[data$sim_id == "real",]
  sims = data[data$sim_id != "real",]
  p <- ggplot() + 
    geom_histogram(aes(sims$stop_density), binwidth = binwidth, fill = "#d1d2d4", col = "#222222") +
    geom_vline(xintercept = real$stop_density, lty = 1, size = 2, col = "#e74b4f") +
    labs(x = xlab, y = ylab, title = title) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5),
      plot.margin = unit(c(0.2,0.1,0.2,0.1), units = "in")
    ) 
  return(p)
}

density_histogram(data, xlab = "Stop codon density", ylab = "Count", binwidth = 0.01, title = "INT3")

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
int3_plot <- density_histogram(data, xlab = "Stop codon density", ylab = "Count", binwidth = 0.01, title = "INT3")
ggsave("clean_run/plots/dinucleotide_matched_stop_densities_int3.pdf", plot = int3_plot)
stop_density = real_density(data, "sim_id", "real", "stop_density")
normalised_density(data, "sim_id", "real", "stop_density")
emperical_p(data, "sim_id", "real", "stop_density")

filepath = "clean_run/motif_tests/RESCUE_stop_codon_densities.csv"
data = read.csv(filepath, head = T)
rescue_plot <- density_histogram(data, xlab = "Stop codon density", ylab = "Count", binwidth = 0.008, title = "RESCUE")
ggsave("clean_run/plots/dinucleotide_matched_stop_densities_RESCUE.pdf", plot = rescue_plot)
stop_density = real_density(data, "sim_id", "real", "stop_density")
stop_density
normalised_density(data, "sim_id", "real", "stop_density")
emperical_p(data, "sim_id", "real", "stop_density")

filepath = "clean_run/motif_tests/ke400_stop_codon_densities.csv"
data = read.csv(filepath, head = T)
ke400_plot <- density_histogram(data, xlab = "Stop codon density", ylab = "Count", binwidth = 0.003, title = "Ke400")
ggsave("clean_run/plots/dinucleotide_matched_stop_densities_ke400.pdf", plot = ke400_plot)
stop_density = real_density(data, "sim_id", "real", "stop_density")
stop_density
normalised_density(data, "sim_id", "real", "stop_density")
emperical_p(data, "sim_id", "real", "stop_density")

filepath = "clean_run/motif_tests/PESE_stop_codon_densities.csv"
data = read.csv(filepath, head = T)
pese_plot <- density_histogram(data, xlab = "Stop codon density", ylab = "Count", binwidth = 0.001, title = "PESE")
ggsave("clean_run/plots/dinucleotide_matched_stop_densities_PESE.pdf", plot = pese_plot)
stop_density = real_density(data, "sim_id", "real", "stop_density")
stop_density
normalised_density(data, "sim_id", "real", "stop_density")
emperical_p(data, "sim_id", "real", "stop_density")

filepath = "clean_run/motif_tests/ESR_stop_codon_densities.csv"
data = read.csv(filepath, head = T)
esr_plot <- density_histogram(data, xlab = "Stop codon density", ylab = "Count", binwidth = 0.005, title = "ESR")
ggsave("clean_run/plots/dinucleotide_matched_stop_densities_ESR.pdf", plot = esr_plot)
stop_density = real_density(data, "sim_id", "real", "stop_density")
stop_density
normalised_density(data, "sim_id", "real", "stop_density")
emperical_p(data, "sim_id", "real", "stop_density")

filepath = "clean_run/motif_tests/ises_wang_stop_codon_densities.csv"
data = read.csv(filepath, head = T)
ise_plot = density_histogram(data, xlab = "Stop codon density", ylab = "Count", binwidth = 0.0075, title = "ISE")
ggsave("clean_run/plots/dinucleotide_matched_stop_densities_ises_wang.pdf", plot = ise_plot)
stop_density = real_density(data, "sim_id", "real", "stop_density")
stop_density
normalised_density(data, "sim_id", "real", "stop_density")
emperical_p(data, "sim_id", "real", "stop_density")

filepath = "clean_run/motif_tests/iss_stop_codon_densities.csv"
data = read.csv(filepath, head = T)
iss_plot = density_histogram(data, xlab = "Stop codon density", ylab = "Count", binwidth = 0.0075, title = "ISS")
ggsave("clean_run/plots/dinucleotide_matched_stop_densities_iss.pdf", plot = iss_plot)
stop_density = real_density(data, "sim_id", "real", "stop_density")
stop_density
normalised_density(data, "sim_id", "real", "stop_density")
emperical_p(data, "sim_id", "real", "stop_density")

filepath = "clean_run/motif_tests/ess_fas_hex2_stop_codon_densities.csv"
data = read.csv(filepath, head = T)
fas_hex2_plot = density_histogram(data, xlab = "Stop codon density", ylab = "Count", binwidth = 0.01, title = "FAS-hex2 ESS")
ggsave("clean_run/plots/dinucleotide_matched_stop_densities_ess_hex2.pdf", plot = fas_hex2_plot)
stop_density = real_density(data, "sim_id", "real", "stop_density")
stop_density
normalised_density(data, "sim_id", "real", "stop_density")
emperical_p(data, "sim_id", "real", "stop_density")

filepath = "clean_run/motif_tests/ess_fas_hex3_stop_codon_densities.csv"
data = read.csv(filepath, head = T)
fas_hex3_plot = density_histogram(data, xlab = "Stop codon density", ylab = "Count", binwidth = 0.005, title = "FAS-hex3 ESS")
ggsave("clean_run/plots/dinucleotide_matched_stop_densities_ess_hex3.pdf", plot = fas_hex3_plot)
stop_density = real_density(data, "sim_id", "real", "stop_density")
stop_density
normalised_density(data, "sim_id", "real", "stop_density")
emperical_p(data, "sim_id", "real", "stop_density")

filepath = "clean_run/motif_tests/rbp_motifs_cds_stop_codon_densities.csv"
data = read.csv(filepath, head = T)
rbp_cds_plot = density_histogram(data, xlab = "Stop codon density", ylab = "Count", binwidth = 0.005, title = "RBPs (CDS-binding)")
ggsave("clean_run/plots/dinucleotide_matched_stop_densities_rbp_motifs_cds.pdf", plot = rbp_cds_plot)
stop_density = real_density(data, "sim_id", "real", "stop_density")
stop_density
normalised_density(data, "sim_id", "real", "stop_density")
emperical_p(data, "sim_id", "real", "stop_density")

filepath = "clean_run/motif_tests/rbp_motifs_non_cds_stop_codon_densities.csv"
data = read.csv(filepath, head = T)
rbp_non_cds_plot = density_histogram(data, xlab = "Stop codon density", ylab = "Count", binwidth = 0.004, title = "RBPs (non CDS-binding)")
ggsave("clean_run/plots/dinucleotide_matched_stop_densities_rbp_motifs_non_cds.pdf", plot = rbp_non_cds_plot)
stop_density = real_density(data, "sim_id", "real", "stop_density")
stop_density
normalised_density(data, "sim_id", "real", "stop_density")
emperical_p(data, "sim_id", "real", "stop_density")


ese_plots = grid.arrange(
  int3_plot,
  rescue_plot,
  ke400_plot,
  pese_plot,
  esr_plot,
  ncol = 2
)

other_plots = grid.arrange(
  ise_plot,
  iss_plot,
  # fas_hex2_plot,
  # fas_hex3_plot,
  rbp_cds_plot,
  # rbp_non_cds_plot,
  ncol = 1
)

plot = ggarrange(
  ese_plots,
  other_plots,
  labels = c("A", "B"),
  widths = c(2, 1)
)
plot
ggsave("clean_run/plots/motif_sets_stop_densities.pdf", width = 12, height= 8, plot = plot)

