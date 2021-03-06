library(ggplot2)
library(ggpubr)
library(gridExtra)

grey = "#d1d2d4"
fill_colour = "#2678b2"
line_colour = "#222222"
red_colour = "#d4292f"

get_n <- function(x, vjust=0){
  data = data.frame(y = max(x)+vjust,label = paste("n = ", length(x), sep=""))
  return(data)
}
normalised_density_plot = function(data, stops) {
  # plot of normalised densities #
  plot = ggplot() +
    geom_histogram(aes(data$nd), bins = 30, col = line_colour, fill = fill_colour) + 
    geom_vline(xintercept = stops$nd, lty = 1, cex = 2, col = red_colour) + 
    labs(x = "Codon set (FE)", y = "Count") +
    theme_minimal()
  return(plot)
}

purine_boxplot = function(data, stops) {
  # data$colour <- ifelse(data$purine == stops$purine, "a", "b")
  plot = ggplot(data, aes(x = data$purine, y = data$nd)) + 
    geom_hline(yintercept = 0, lty=1) +
    stat_boxplot(geom ='errorbar') +
    geom_boxplot(fill=grey) +
    scale_fill_manual(values=c(fill_colour, grey)) +
    geom_hline(yintercept = stops$nd, lty=2) +
    labs(x = "Codon set purine content", y="FE") + 
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
# filepath = "clean_run/motif_tests/int3_densities_strictly_no_stops.csv"
# filepath = "clean_run/motif_tests/int3_densities_non_overlapping.csv"
file = read.csv(filepath, head = T)

file$gc = substr(file$gc_content, 0, 3)
file$purine = substr(file$purine_content, 0, 3)



stops =  file[file$codons == "TAA_TAG_TGA",]
gc_matched = file[file$gc == stops$gc & file$codons != "TAA_TAG_TGA",]
purine_matched = file[file$purine_content == stops$purine_content,]

nrow(purine_matched)

# gc matched
nrow(gc_matched)
nrow(gc_matched[gc_matched$nd > stops$nd,])
binom.test(nrow(gc_matched[gc_matched$nd > stops$nd,]), nrow(gc_matched), alternative = "g")

norm_density_gc_plot = normalised_density_plot(gc_matched, stops)
norm_density_gc_plot
# ggsave(plot = norm_density_gc_plot, "clean_run/plots/codon_sets_nd_gc_matched.pdf", width = 12, height= 5, plot = plot)

# purine matched
nrow(purine_matched)
nrow(purine_matched[purine_matched$nd > stops$nd,])
binom.test(nrow(purine_matched[purine_matched$nd > stops$nd,]), nrow(purine_matched), alternative = "g")

norm_density_purine_plot = normalised_density_plot(purine_matched, stops)
ggsave(plot = norm_density_purine_plot, "clean_run/plots/codon_sets_nd_purine_matched.pdf", width = 6, height= 5)


# match by gc and purine
gc_purine_matched = gc_matched[gc_matched$purine_content == stops$purine_content,]
binom.test(nrow(gc_purine_matched[gc_purine_matched$nd > stops$nd,]), nrow(gc_purine_matched), alternative = "g")
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
plot
ggsave(plot = plot, file = "clean_run/plots/codon_histogram_boxplot_fe.pdf", width = 8, height = 10)
ggsave(plot = plot, file = "clean_run/plots/codon_histogram_boxplot_fe.eps", width = 8, height = 10)

