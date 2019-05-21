stop_density_simulations = read.csv("clean_run/tests/lincrna/stop_density_simulation_outputs.csv", head = T)
median(stop_density_simulations$stop_codon_density)
median(stop_density_simulations$median_simulated_stop_codon_density)



file = read.csv("clean_run/tests/lincrna/exon_region_densities.csv", head = T)

head(file)
data = file[, c("upstream_density", "core_density", "downstream_density")]
data = melt(data)
head(data)


plot = ggplot(data = data, aes(x = variable, y = value)) +
  stat_boxplot(geom ='errorbar') +
  geom_boxplot(fill = "RoyalBlue") +
  labs(x = "Exon region", y = "Stop codon density")
plot

cairo_pdf("clean_run/plots/lincrna_exon_regions.pdf", width = 10)  
plot
dev.off()


kruskal.test(file$upstream_density, file$core_density, file$downstream_density)
x = c(file$upstream_density, file$core_density, file$downstream_density)
g = factor(rep(1:3, c(length(file$upstream_density), length(file$core_density), length(file$downstream_density))), labels = c("ud", "cd", "dd"))
dunn.test(x, g, method = "bonferroni")


file = read.csv("clean_run/tests/lincrna/compare_exon_intron_stop_density.csv", head = T)

head(file)

wilcox.test(file$exon_density, file$intron_density, paired = T)

head(file)
data = file[, c("exon_density", "intron_density")]
data = melt(data)

head(data)


plot = ggplot(data = data, aes(x = variable, y = value)) +
  stat_boxplot(geom ='errorbar') +
  geom_boxplot(fill = "RoyalBlue") +
  labs(x = "Region", y = "Stop codon density")
plot

cairo_pdf("clean_run/plots/lincrna_exon_intron_densities.pdf", width = 10)  
plot
dev.off()

wilcox.test(file$exon_density, file$intron_density, paired = T)

median(file$exon_density)
median(file$intron_density)

