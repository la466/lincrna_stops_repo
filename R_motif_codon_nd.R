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

density_scatterplot <- function(data, ycol = "density", xlab = "", ylab = "") {
  data$id = as.numeric(seq(1:nrow(data)))
  stops = data[data$codons == "TAA_TAG_TGA",]
  not_stops = data[data$codons != "TAA_TAG_TGA",]
  p <- ggplot() +
    geom_point(aes(x = not_stops$id, y = not_stops$nd), cex = 0.5, pch = 16) +
    geom_point(aes(x = stops$id, y = stops$nd), cex = 3, pch = 16, col = "red") +
    # scale_x_continuous(labels = not_stops$codons, breaks = seq(1:nrow(not_stops))) +
    theme(axis.text.x = element_blank()) +
    geom_hline(yintercept = 0, lty = 2) +
    labs(x = xlab, y = ylab)

  return(p)
}

density_boxplot <- function(data, ycol = "density") {
  stops_purine = data$purine[data$codons=="TAA_TAG_TGA"]
  stops_density =  data[[ycol]][data$codons=="TAA_TAG_TGA"]
  data$purine <- as.factor(data$purine)
  data$col <- ifelse(data$purine == stops_purine, "a", "b")
  p <- ggplot(data, aes(x=data$purine, y=data[[ycol]])) +
    stat_boxplot(geom ='errorbar') +
    geom_boxplot(aes(fill=data$col)) +
    scale_fill_manual(values=c("RoyalBlue", "#dddddd")) +
    geom_hline(yintercept = stops_density, lty=2) +
    labs(x = "Query codons purine content", y="Query codons ND") +
    annotate("text", x=min(as.numeric(data$purine)), hjust=0.2, y=stops_density + 0.05, label="TAA,TAG,TGA", cex=3) +
    annotate("text", x=min(as.numeric(data$purine)), hjust=0.1, y=stops_density - 0.05, label=round(stops_density,4), cex=3) +
    theme(legend.position="none") +
    scale_x_discrete(labels = c("0-0.1", "0.1-0.2", "0.2-0.3", "0.3-0.4", "0.4-0.5", "0.5-0.6", "0.6-0.7", "0.7-0.8", "0.8-0.9", "0.9-1", "1")) +
    stat_summary(fun.data = get_n, fun.args = list("vjust" = 0.1), geom = "text", aes(group = "purine"), size=3)
  return(p)
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

stops =  file[file$codons == "TAA_TAG_TGA",]
gc_matched = file[file$gc == stops$gc & file$codons != "TAA_TAG_TGA",]


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

