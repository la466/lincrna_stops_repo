library(ggplot2)
library(reshape2)

density_boxplot = function(data, xlab = "", ylab = "") {
  p <- ggplot(data = data, aes(x = data$variable, y = data$value)) + 
    stat_boxplot(geom ='errorbar') +
    geom_boxplot(fill = "RoyalBlue") +
    scale_x_discrete(labels = c("Exons", "Introns", "Intron scaled")) + 
    labs(x = xlab, y = ylab)
  return(p)
}

filepath  = "clean_run/tests/compare_exon_intron_stop_density.csv"
file = read.csv(filepath, head = T)

data = file[, c("exon_density", "intron_density", "intron_density_scaled")]
data = melt(data)


plot1 = density_boxplot(data, ylab = "Stop codon density")
cairo_pdf("clean_run/plots/exon_intron_stop_density.pdf", width = 10)
plot1
dev.off()

wilcox.test(file$exon_density, file$intron_density, paired = T)
wilcox.test(file$exon_density, file$intron_density_scaled, paired = T)

median(file$exon_density)
median(file$intron_density)
median(file$intron_density_scaled)

median(file$exon_gc)
median(file$intron_gc)


nc_file = read.csv("clean_run/tests/non_coding_exons/non_coding_exons.csv", head = T)
colnames(nc_file) = c("non-coding exons", "non-coding exons scaled")
nc_data = melt(nc_file)
nc_data = nc_data[, c("variable", "value")]

full_data = rbind(data, nc_data)

plot = density_boxplot(full_data, "Type", "Stop density")
cairo_pdf("clean_run/plots/non_coding_exon_comparison.pdf", width = 12)
plot
dev.off()

kruskal.test(full_data$value~full_data$variable)
x = full_data$value

# x = c(file$upstream_density, file$core_density, file$downstream_density)
g = factor(rep(1:5, c(length(file$exon_density), length(file$intron_density), length(file$intron_density_scaled), length(nc_file$density), length(nc_file$scaled_density))), labels = c("e", "i", "is", "nc", "ncs"))
dunn.test(x, g, method = "bonferroni")
