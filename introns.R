install.packages("ggExtra")
library(ggplot2)
library(ggpubr)
library(Cairo)
library(LSD)
library(ggExtra)

fill_colour = "#d1d2d4"
line_colour = "#222222"
red_colour = "#e74b4f"


file = read.csv("clean_run/tests/introns/intron_hexamers.csv", head = T)


real = file[file$id == "real",]
sims = file[file$id != "real",]

max_purine = max(real$purine_content, sims$purine_content)
min_purine = min(real$purine_content, sims$purine_content)
max_a = max(real$a_content, sims$a_content)
min_a = min(real$a_content, sims$a_content)
max_g = max(real$g_content, sims$g_content)
min_g = min(real$g_content, sims$g_content)

heatscatter(sims$g_content, sims$a_content, xlim = c(min_g, max_g), ylim = c(min_a, max_a))


hist = ggplot() + 
  geom_histogram(aes(sims$purine_content), binwidth = 0.01, fill = fill_colour, col = line_colour) +
  geom_vline(xintercept = real$purine_content, lty = 1, size = 1.5, col = red_colour) +
  labs(x = "Motif set purine content", y = "Count")

# scatter =  ggplot() + 
#   geom_point(aes(x = sims$g_content, y = sims$a_content)) +
#   geom_point(aes(x = real$g_content, y = real$a_content), col = "red", pch = 18, cex = 3) +
#   scale_x_continuous(limits = c(min_g, max_g)) +
#   scale_y_continuous(limits = c(min_a, max_a)) +
#   labs(x = "Motif set G proportion", y = "Motif set A proportion")

plot = ggplot(data=sims, aes(g_content,a_content)) + 
  geom_point(cex = 0.8, col = "RoyalBlue") +
  geom_point(aes(x = real$g_content, y = real$a_content), pch = 18, cex = 3, col = red_colour) +
  labs(x = "Motif set G proportion", y = "Motif set A proportion") +
  guides(alpha="none") +
  scale_x_continuous(limits = c(min_g, max_g)) +
  scale_y_continuous(limits = c(min_a, max_a)) + 
  theme_minimal() +
  theme(legend.position = "none")
plot = ggMarginal(plot, type="histogram", fill = fill_colour, col = line_colour)

plot

# final_plot = ggarrange(hist, scatter, ncol = 2, nrow = 1, labels = c("A", "B"))
final_plot = ggarrange(hist, plot, ncol = 2, nrow = 1, labels = c("A", "B"))

ggsave("clean_run/plots/intron_hexamers.eps", plot = final_plot, width = 12, height = 5)



empirical_p = function(real, sims, col, greater = TRUE) {
  if (greater) {
    p = (nrow(sims[sims[[col]] >= real[[col]],]) + 1) / (nrow(sims) + 1)
  } else {
    p = (nrow(sims[sims[[col]] <= real[[col]],]) + 1) / (nrow(sims) + 1)
  }
  return (p)
}

empirical_p(real, sims, "purine_content", greater = T)
empirical_p(real, sims, "a_content", greater = T)
empirical_p(real, sims, "g_content", greater = T)



