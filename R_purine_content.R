library(ggplot2)
library(reshape2)

fill_colour = "#d1d2d4"
line_colour = "#222222"
red_colour = "#e74b4f"


purine_boxplot = function(data) {
  data = purine
  colnames(data) = c("id", "Exons", "Introns")
  data = melt(data)
  plot = ggplot() + 
    geom_boxplot(aes(x = data$variable, y = data$value), fill = fill_colour) + 
    labs(x = "", y = "Purine content") + 
    geom_hline(yintercept = 0.784, col = red_colour, size = 1.1) +
    theme_minimal()
  return(plot)
}
 

# exon intron purine
purine = read.csv("clean_run/tests/purine_content/exon_intron_purine.csv", head = T)
nrow(purine)
wilcox.test(purine$exon_purine, purine$intron_purine, paired = T)
median(purine$exon_purine)
median(purine$intron_purine)

plot = purine_boxplot(purine)
plot
ggsave(plot = plot, file = "clean_run/plots/exon_intron_purine_boxplot.pdf", width = 4, height = 4)

# no eses
purine_no_eses = read.csv("clean_run/tests/purine_content/exon_intron_purine_no_eses.csv", head = T)
median(purine_no_eses$exon_purine)
median(purine_no_eses$intron_purine)
wilcox.test(purine_no_eses$exon_purine, purine_no_eses$intron_purine, paired= T)


# purine content of motifs
file = read.csv("clean_run/tests/purine_content/int3_purine_content.csv", head = T)
file = file[complete.cases(file), ]
nrow(file)
# all ese
median(file$ese_purine)
median(file$non_ese_purine)
wilcox.test(file$ese_purine, file$non_ese_purine, paired = T)
# just terminal 50 nts
median(file$flanking_50_ese_purine)
median(file$flanking_50_non_ese_purine)
wilcox.test(file$flanking_50_ese_purine, file$flanking_50_non_ese_purine, paired = T)





