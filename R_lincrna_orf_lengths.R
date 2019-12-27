library(ggplot2)
library(ggpubr)

grey = "#d1d2d4"
fill_colour = "#2678b2"
line_colour = "#222222"
red_colour = "#d4292f"


gc_zscore_plot = function(data) {
  print(head(data))
plot = ggplot(data = data, aes(x = data$gc, y = data$z)) + 
  geom_point(col = ifelse(data$p < 0.05, fill_colour, "black"), size = ifelse(data$p < 0.05, 1.5, 1), pch = ifelse(data$p < 0.05, 18, 16) ) + 
  geom_smooth(method = "lm", se = FALSE, col = "black", fullrange = TRUE) + 
  geom_hline(yintercept = 0, lty = 2) + 
  scale_y_continuous(limits = c(-5, 10)) + 
  labs(x = "GC", y = "Z score") + 
  theme_minimal()
return(plot)
}

predict_threshold = function(directory, alpha) {
  new_data = data.frame("threshold" = double(), "greater" = double())
  thresholds = seq(200,600,10)
  thresholds
  for (threshold in thresholds) {
    filepath = paste(directory, threshold, ".csv", sep = "") 
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
  pass_threshold = approx(x=model1$fitted, y=model1$x, xout=alpha)$y
  return(pass_threshold)
}

threshold_plot = function(directory, alpha) {
  # set up a new data frame with predicted values in intervals of 5
  # need this to fill in under the smoothed line
  
  new_data = data.frame("threshold" = double(), "greater" = double())
  thresholds = seq(200,600,10)
  thresholds
  for (threshold in thresholds) {
    filepath = paste(directory, threshold, ".csv", sep = "") 
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
  pass_threshold = approx(x=model1$fitted, y=model1$x, xout=alpha)$y
  
  new_thresholds = data.frame(threshold=seq(200,600,1))
  vals1 = predict(model1, new_thresholds)
  model_data = data.frame(threshold = new_thresholds, greater = vals1)
  

  
  plot = ggplot(data = new_data, aes(x = threshold, y = greater)) + 
    geom_ribbon(data=subset(new_data, 200 <= threshold & threshold <= 300), aes(ymin = 0, ymax = predict(loess(greater ~ threshold))), fill = grey) +
    geom_ribbon(data=subset(model_data, 300 <= threshold & threshold <= round(pass_threshold)), aes(ymin = 0, ymax = predict(loess(greater ~ threshold))), fill = red_colour) +
    geom_hline(yintercept = 0.05, lty = 2) + 
    # geom_vline(xintercept = pass_threshold) +
    geom_point() +
    # geom_line(aes(x = thresholds, y = predicts)) +
    geom_smooth(method = "loess", col = line_colour, se = T, size = 0.8) +
    scale_x_continuous(breaks = seq(0, 600, 25)) +
    scale_y_continuous(breaks = seq(0, 0.4, 0.05)) +
    labs(x = "ORF length threshold", y = "Excess proportion of sequences with longest\nORF lengths greater expected") +
    theme_minimal() + 
    theme(
      panel.grid.minor.x  = element_blank(),
      panel.grid.minor.y = element_blank()
    )
  return(plot)
}


###

file = read.csv("clean_run/tests/lincrna/orf_length_sim/cabili_sim_orf_lengths_zs_grouped.csv", head = T)
file = file[!is.nan(file$real),]

# number longer than sims
nrow(file[file$z > 0,])
nrow(file[file$z > 0,]) / nrow(file)
# significantly longer
nrow(file[file$z > 0 & file$p < 0.05, ,])
nrow(file[file$z > 0 & file$p < 0.05, ,]) / nrow(file)
# less
nrow(file[file$z < 0,])
nrow(file[file$z < 0 & file$p < 0.05,])
# correlation
cor.test(file$gc, file$z, method = "spearman")

#binomial tests
binom.test(nrow(file[file$z > 0,]), nrow(file), alternative = "g")
binom.test(nrow(file[file$z > 0 & file$p < 0.05,]), nrow(file), p = 0.05,  alternative = "g")
binom.test(nrow(file[file$z < 0 & file$p < 0.05,]), nrow(file), p = 0.05,  alternative = "g")

# plot
zplot = gc_zscore_plot(file)
ggsave(plot = zplot, "clean_run/plots/lincrna_orf_length_sim.pdf", width = 7, height = 5)


# some stats on orf lengths
median(file$real)
max(file$real)
median(file$median_sims)
max(file$median_sims)
max(file$z)


# now look at how many exceed thresholds using each of the threshold parameters
# longer than 300 and compare with simulants
filepath = paste("clean_run/tests/lincrna/total_lengths1/", "300", ".csv", sep = "") 
threshold_file = read.csv(filepath, head = T)

real = threshold_file[threshold_file$id == "real",]
real$greater
sims = threshold_file[threshold_file$id != "real",]
max(sims$max_orf)
mean(sims$greater)
sd(sims$greater)


# p value for number of sets with greater
(nrow(sims[sims$greater >= real$greater,])+1) / (nrow(sims) + 1)
# excess excess 
(real$greater - mean(sims$greater)) / real$total


alpha_threshold = predict_threshold("clean_run/tests/lincrna/total_lengths1/", 0.05)
alpha_threshold

alpha_plot = threshold_plot("clean_run/tests/lincrna/total_lengths1/", 0.05)
zplot

plot = ggarrange(
  zplot, 
  alpha_plot,
  ncol = 2,
  labels = c("A", "B"),
  widths = c(4, 6)
)


ggsave(plot = plot, filename = "clean_run/plots/orf_lengths_z_thresholds.pdf", width = 10, height = 5)
ggsave(plot = plot, filename = "clean_run/plots/orf_lengths_z_thresholds.eps", width = 10, height = 5)
dev.off()





