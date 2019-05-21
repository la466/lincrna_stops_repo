library(ggplot2)
library(ggpubr)
library(gridExtra)

fill_colour = "#d1d2d4"
line_colour = "#222222"
red_colour = "#e74b4f"

get_n <- function(x, vjust=0){
  data = data.frame(y = max(x)+vjust,label = paste("n = ", length(x), sep=""))
  return(data)
}
normalised_density_plot = function(data, stops) {
  # plot of normalised densities #
  plot = ggplot() +
    geom_histogram(aes(data$nd), bins = 30, col = line_colour, fill = fill_colour) + 
    geom_vline(xintercept = stops$nd, lty = 1, cex = 2, col = red_colour) + 
    labs(x = "Codon set normalised density (ND)", y = "Count") +
    theme_minimal()
  return(plot)
}

purine_boxplot = function(data, stops) {
  data$colour <- ifelse(data$purine == stops$purine, "a", "b")
  plot = ggplot(data, aes(x = data$purine, y = data$nd)) + 
    geom_hline(yintercept = 0, lty=1) +
    stat_boxplot(geom ='errorbar') +
    geom_boxplot(aes(fill=data$colour)) +
    scale_fill_manual(values=c("RoyalBlue", fill_colour)) +
    geom_hline(yintercept = stops$nd, lty=2) +
    labs(x = "Codon set purine content", y="ND") + 
    annotate("text", x=min(as.numeric(data$purine)), hjust= -0.2, y= stops$nd + 0.05, label="TAA,TAG,TGA", cex=3) +
    annotate("text", x=min(as.numeric(data$purine)), hjust= -0.5, y= stops$nd - 0.05, label=round(stops$nd,3), cex=3) +
    theme(legend.position="none") +
    scale_x_discrete(labels = c("0-0.1", "0.1-0.2", "0.2-0.3", "0.3-0.4", "0.4-0.5", "0.5-0.6", "0.6-0.7", "0.7-0.8", "0.8-0.9", "0.9-1", "1")) +
    stat_summary(fun.data = get_n, fun.args = list("vjust" = 0.1), geom = "text", aes(group = "purine"), size = 3) + 
    theme_minimal() +
    theme(
      legend.position = "none",
    )
  return(plot)
  }


binom_test <- function(data, ycol = "density", group = NULL) {
  stops = data[data$codons == "TAA_TAG_TGA",]
  not_stops = data[data$codons != "TAA_TAG_TGA",]
  if (!is.null(group)) {
    not_stops = not_stops[not_stops[[group]] == stops[[group]],]
  }
  b = binom.test(nrow(not_stops[not_stops[[ycol]] <= stops[[ycol]], ]), nrow(not_stops), alternative = "l")
  return(b)
}

####

filepath = "clean_run/motif_tests/int3_densities.csv"
file = read.csv(filepath, head = T)
file$gc = substr(file$gc_content, 0, 3)
file$purine = substr(file$purine_content, 0, 3)

nrow(file)


stops =  file[file$codons == "TAA_TAG_TGA",]
gc_matched = file[file$gc == stops$gc & file$codons != "TAA_TAG_TGA",]

pm = file[file$purine_content == stops$purine_content,]

nrow(pm)
nrow(pm[pm$nd > stops$nd,]) / nrow(pm)

binom.test(nrow(pm[pm$nd > stops$nd,]), nrow(pn), alternative =)

hist(gc_matched$density)
abline(v = stops$density)
nrow(gc_matched[gc_matched$density > stops$density,])
nrow(gc_matched)

norm_density_gc_plot = normalised_density_plot(gc_matched, stops)
# ggsave(plot = norm_density_gc_plot, "clean_run/plots/codon_sets_densities_nds.pdf", width = 12, height= 5, plot = plot)

# normalised density binom test
binom.test(nrow(gc_matched[gc_matched$nd > stops$nd,]), nrow(gc_matched), alternative = "g")

# match by purine too
purine_matched = gc_matched[gc_matched$purine_content == stops$purine_content,]
binom.test(nrow(purine_matched[purine_matched$nd > stops$nd,]), nrow(purine_matched), alternative = "g")
# plot
purine_plot = purine_boxplot(gc_matched, stops)
purine_plot

plot = ggarrange(
  norm_density_gc_plot,
  purine_plot,
  labels = c("A", "B"),
  nrow = 2, 
  ncol = 1
)

ggsave(plot = plot, file = "clean_run/plots/codon_nd_histogram_boxplot.pdf", width = 8, height = 10)
