library(ggplot2)
fill_colour = "#d1d2d4"
line_colour = "#222222"
red_colour = "#e74b4f"

file = read.csv("clean_run/tests/lincrna/sim_orf_lengths_zs.csv", head = T)

head(file)

cor.test(file$gc, file$normalised_length, method = "spearman")

head(file)

nrow(file[file$z_score > 0,])
nrow(file[file$z_score > 0,]) / nrow(file)

nrow(file[file$z_score > 1.96,])
file$p = 2*pnorm(-abs(file$z_score))

nrow(file[file$p < 0.05 & file$z_score > 0 & !is.na(file$z_score),])
nrow(file[file$p < 0.05 & file$z_score > 0 & !is.na(file$z_score),])

cor.test(file$gc, file$z_score, method = "spearman")


plot = ggplot(data = file, aes(x = file$gc, y = file$z_score)) + 
  geom_point(col = ifelse(file$p < 0.05, "blue", "black")) + 
  geom_smooth(method = "lm", se = FALSE, col = "black") + 
  geom_hline(yintercept = 0, lty = 2) + 
  labs(x = "GC", y = "Z score") + 
  theme_minimal()

plot

ggsave(plot = plot, "clean_run/plots/lincRNA_orf_length_z.pdf", width = 10, height = 7)


head(file)



cor.test(file$gc, file$real, method = "spearman")

plot(file$gc, file$real, cex = 0.6, col = "blue", pch = 16)
abline(lm(file$real~file$gc), col = "blue")
points(file$gc, file$mean_sims, cex = 0.6, col = "red", pch = 16)
abline(lm(file$mean_sims~file$gc), col = "red")

cor.test(file$gc, file$mean_sims, method = "spearman")

head(file)
nrow(file[file$z_score > 1.96,])
nrow(file)

plot(file$gc, file$real)
cor.test(file$gc, file$real, method = "spearman")

lengths_file = read.csv("clean_run/tests/lincrna/lincRNA_lengths.csv", head = T)
head(lengths_file)

plot(lengths_file$gc, lengths_file$length)
cor.test(lengths_file$gc, lengths_file$length, method = "spearman")

hist(lengths_file$gc)



file = read.csv("clean_run/tests/lincrna/cabili/orf_length_sim.csv", head = T)
head(file)



plot(file$gc, file$z_score)
nrow(file[file$z_score > 1.96,])
nrow(file)
binom.test(273, 1921, p = 0.05, alternative = "g")




file = read.csv("clean_run/tests/lincrna/orf_lengths/hangauer_sim_orf_lengths_zs_grouped.csv", head = T)
head(file)




plot1 = gc_zscore_plot(file)
plot1

nrow(file[file$nd > 0,])
binom.test(nrow(file[file$nd > 0,]), nrow(file), alternative = "g")
binom.test(nrow(file[file$nd > 0 & file$empirical_p < 0.05,]), nrow(file), p = 0.05,  alternative = "g")
binom.test(nrow(file[file$nd < 0 & file$empirical_p < 0.05,]), nrow(file), p = 0.05,  alternative = "g")
binom.test(nrow(file[file$z > 0,]), nrow(file), alternative = "g")
binom.test(nrow(file[file$z > 0 & file$p < 0.05,]), nrow(file), p = 0.05,  alternative = "g")
binom.test(nrow(file[file$z < 0 & file$p < 0.05,]), nrow(file), p = 0.05,  alternative = "g")

cor.test(file$gc, file$nd, method = "spearman")
cor.test(file$gc, file$z, method = "spearman")


# used this

file = read.csv("clean_run/tests/lincrna/orf_length_sim/cabili_sim_orf_lengths_zs_grouped.csv", head = T)
file = file[!is.nan(file$real),]

nrow(file[file$real >= 300,])
nrow(file[file$real >= 300 & file$z > 0,])

gc_zscore_plot = function(data) {
  plot = ggplot(data = data, aes(x = data$gc, y = data$z)) + 
    geom_point(col = ifelse(data$p < 0.05, "RoyalBlue", "black"), size = ifelse(data$p < 0.05, 1.5, 1), pch = ifelse(data$p < 0.05, 18, 16) ) + 
    geom_smooth(method = "lm", se = FALSE, col = "black", fullrange = TRUE) + 
    geom_hline(yintercept = 0, lty = 2) + 
    scale_y_continuous(limits = c(-5, 10)) + 
    labs(x = "GC", y = "Z score") + 
    theme_minimal()
  return(plot)
}


nrow(file[file$real >= 370 & file$z > 0,]) / nrow(file)



head(file)
median(file$real)
max(file$real)
median(file$median_sims)
max(file$median_sims)

zplot = gc_zscore_plot(file)
plot
ggsave(plot = plot, "clean_run/plots/lincrna_orf_length_sim.pdf", width = 7, height = 5)

max(file$z)

nrow(file[file$nd > 0,])
binom.test(nrow(file[file$nd > 0,]), nrow(file), alternative = "g")
binom.test(nrow(file[file$nd > 0 & file$empirical_p < 0.05,]), nrow(file), p = 0.05,  alternative = "g")
binom.test(nrow(file[file$nd < 0 & file$empirical_p < 0.05,]), nrow(file), p = 0.05,  alternative = "g")
binom.test(nrow(file[file$z > 0,]), nrow(file), alternative = "g")
binom.test(nrow(file[file$z > 0 & file$p < 0.05,]), nrow(file), p = 0.05,  alternative = "g")
binom.test(nrow(file[file$z < 0 & file$p < 0.05,]), nrow(file), p = 0.05,  alternative = "g")

file[file$nd < 0 & file $empirical_p < 0.05,]

head(file)


cor.test(file$gc, file$nd, method = "spearman")
cor.test(file$gc, file$z, method = "spearman")

library(ggpubr)
plot = ggarrange(
  plot1,
  plot2,
  labels = c("A", "B")
)
ggsave(plot = plot, "clean_run/plots/lincrna_orf_length_sims.pdf", width = 10, height = 5)


filepath = paste("clean_run/tests/lincrna/total_lengths1/", 350, ".csv", sep = "") 
file = read.csv(filepath, head = T)
head(file)
real = file[file$id == "real",]
sims = file[file$id != "real",]
mean(sims$greater)
real$greater
max(sims$max_orf)
nrow(sims[sims$greater >= real$greater,])
(real$greater - mean(sims$greater)) / real$total


new_data = data.frame("threshold" = double(), "greater" = double())
thresholds = seq(200,600,10)
thresholds
for (threshold in thresholds) {
  filepath = paste("clean_run/tests/lincrna/total_lengths1/", threshold, ".csv", sep = "") 
  file_data = read.csv(filepath, head = T)
  real = file_data[file_data$id == "real",]
  sims = file_data[file_data$id != "real",]
  threshold_excess = (real$greater - mean(sims$greater)) / real$total
  outline = c(threshold, threshold_excess)
  new_data = rbind(new_data, outline)
}
colnames(new_data) <- c("threshold", "greater")



model1 = loess(greater ~ threshold, new_data)
predicts = predict(model1)
line_points = predict(model1)
pass_threshold = approx(x=model1$fitted, y=model1$x, xout=0.05)$y
pass_threshold

# set up a new data frame with predicted values in intervals of 5
# need this to fill in under the smoothed line
new_thresholds = data.frame(threshold=seq(200,600,1))
vals1 = predict(model1, new_thresholds)
model_data = data.frame(threshold = new_thresholds, greater = vals1)

new_data

threshold_plot = ggplot(data = new_data, aes(x = threshold, y = greater)) + 
  geom_ribbon(data=subset(new_data, 200 <= threshold & threshold <= 300), aes(ymin = 0, ymax = predict(loess(greater ~ threshold))), fill = fill_colour) +
  geom_ribbon(data=subset(model_data, 300 <= threshold & threshold <= round(pass_threshold)), aes(ymin = 0, ymax = predict(loess(greater ~ threshold))), fill = red_colour) +
  geom_hline(yintercept = 0.05, lty = 2) + 
  # geom_vline(xintercept = pass_threshold) +
  geom_point() +
  # geom_line(aes(x = thresholds, y = predicts)) +
  geom_smooth(method = "loess", col = line_colour, se = T, size = 0.8) +
  scale_x_continuous(breaks = seq(0, 600, 25)) +
  scale_y_continuous(breaks = seq(0, 0.4, 0.05)) +
  labs(x = "ORF length threshold", y = "Excess proportion of sequences with maximum\nORF lengths greater expected") +
  theme_minimal() + 
  theme(
    panel.grid.minor.x  = element_blank(),
    panel.grid.minor.y = element_blank()
  )

plot = ggarrange(
  zplot, 
  threshold_plot,
  ncol = 2,
  labels = c("A", "B")
)

ggsave(plot = plot, filename = "clean_run/plots/orf_lengths_z_thresholds.pdf", width = 12, height = 5)






