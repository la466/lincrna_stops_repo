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

