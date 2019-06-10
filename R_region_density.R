library(ggplot2)
library(reshape2)
library(ggpubr)

data = read.csv("clean_run/tests/lincrna/stop_density/stop_density_regions.csv", head = T)

head(data)



melt = melt(data)
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
