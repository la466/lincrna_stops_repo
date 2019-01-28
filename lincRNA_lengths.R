library(ggplot2)

file = read.csv("clean_run/tests/lincrna/sim_orf_lengths_zs.csv", head = T)

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
  labs(x = "GC", y = "Z score")

ggsave(plot = plot, "clean_run/plots/lincRNA_orf_length_sim.pdf", width = 10, height = 7)



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

