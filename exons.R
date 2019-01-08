library(ggplot2)
library(reshape2)
library(Cairo)
library(dunn.test)


filepath = "clean_run/tests/exon_regions/exon_regions_density.csv"
file = read.csv(filepath, head = T)

head(file)
data = file[, c("upstream_density", "core_density", "downstream_density")]
data = melt(data)
head(data)


plot = ggplot(data = data, aes(x = variable, y = value)) +
  stat_boxplot(geom ='errorbar') +
  geom_boxplot(fill = "RoyalBlue") +
  labs(x = "Exon region", y = "Stop codon density")
plot

cairo_pdf("clean_run/plots/exon_regions.pdf", width = 10)  
plot
dev.off()


kruskal.test(file$upstream_density, file$core_density, file$downstream_density)
x = c(file$upstream_density, file$core_density, file$downstream_density)
g = factor(rep(1:3, c(length(file$upstream_density), length(file$core_density), length(file$downstream_density))), labels = c("ud", "cd", "dd"))
dunn.test(x, g, method = "bonferroni")

wilcox.test(file$upstream_density, file$core_density, paired = T)
median(file$upstream_density)
median(file$core_density)


t.test(file$upstream_density, file$core_density, paired = T)




binom.test(nrow(file[file$upstream_density <= file$core_density,]), nrow(file), alternative = "g")
nrow(file)

filepath = "clean_run/tests/exon_regions/exon_regions_excess.csv"
file = read.csv(filepath, head = T)

head(file)

file$diff = file$upstream_counts - file$expected_upstream_counts
file$chi = (file$diff*file$diff) / file$expected_upstream_counts
file = file[file$expected_upstream_counts != 0,]
sum(file$chi)
mean(file$diff)
