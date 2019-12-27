library(ggplot2)
library(reshape2)
library(ggpubr)
library(data.table)
library(png)

grey = "#d1d2d4"
fill_colour = "#2678b2"
line_colour = "#222222"
red_colour = "#d4292f"



data = read.csv("clean_run/tests/lincrna/stop_density/stop_density_regions.csv", head = T)

melt = melt(data, id.vars = c("region"), measure.vars = c("ese_density", "stop_density"))
melt$region <- factor(melt$region, levels = c("five", "core", "three"))
melt$variable <- factor(melt$variable, c("ese_density", "stop_density"))

head(melt)

region_plot = ggplot(melt, aes(fill=variable, y=value, x=region)) +
  geom_bar(position="dodge", stat="identity") +
  scale_fill_manual(values = c("RoyalBlue", "#e74b4f"), labels = c("ESE", "Stop codon")) +
  scale_x_discrete(labels = c("5' flank", "Core", "3' flank")) +
  scale_y_continuous(limits = c(0, 0.2)) +
  labs(x = "", y = "Density") +
  theme_minimal() +
  theme(
    legend.title = element_blank()
  )
region_plot

ggsave("clean_run/plots/flank_core_ese_stop_density.pdf", plot = region_plot, width = 7, height = 4)


## plot with the per seq densities

data = read.csv("clean_run/tests/lincrna/stop_density/stop_density_regions_per_seq.csv", head = F)

data.t = transpose(data, fill = 0)
head(data.t)
colnames(data.t) <- as.character(unlist(data.t[1,]))
data.t = as.data.frame(data.t[-1, ])
data.t <- sapply( data.t, as.numeric )

i = data.frame(data.t)

meltdata = melt(i)
meltdata$region = c(rep("5'", nrow(i)*2), rep("Core", nrow(i)*2), rep("3'", nrow(i)*2))
meltdata$region <- factor(meltdata$region, levels=unique(meltdata$region))
meltdata$colour = rep(c(rep("ESE", nrow(i)), rep("Stop codon", nrow(i))), 3)

head(meltdata)

img <- readPNG("./source_data/exon.png")
plot = ggplot(data = meltdata, aes(x = region, y = value, fill = colour)) +
  annotation_raster(img, xmin=0.25, xmax=3.75, ymin=0.9, ymax=1.1) +
  stat_boxplot(geom = "errorbar") +
  geom_boxplot() +
  geom_vline(xintercept = 1.5, lty = 2, col = "#aaaaaa") +
  geom_vline(xintercept = 2.5, lty = 2, col = "#aaaaaa") +
  scale_fill_manual(values = c(fill_colour, red_colour)) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Region", y = "Density") +
  theme_minimal() +
  theme(
    legend.title = element_blank()
  ) +
  annotate("text", x = 2, y = 1, label = "Exon")
plot
ggsave(plot, filename = "clean_run/plots/exon_regions_per_seq.pdf", width = 8, height = 6)
ggsave(plot, filename = "clean_run/plots/exon_regions_per_seq.eps", width = 8, height = 6)
ggsave(plot, filename = "clean_run/plots/exon_regions_per_seq.png", width = 8, height = 6)
