library(ggplot2)
library(gridExtra)
library(ggpubr)
library(reshape2)

all_ese_plot <- function(data) {
  p = ggplot(data, aes(x = log10(data$intron_size), y = data$ese_density)) + 
    geom_point(size = 0.8) + 
    geom_smooth(method='lm', col = "black", se = F) + 
    scale_y_continuous(limits = c(0, 1)) +
    labs(x = "log10 intron size", y = "Density") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))
  return (p)
}

split_ese_plot <- function(data, title = NULL) {
  data.melt = melt(data, id.vars = c("id", "intron_size"), measure.vars = c("stop_ese_density", "non_stop_ese_density"), value.name = "density")
  legend_labs = labels = c("Stop ESEs", "Non stop ESEs")
  ggplot(data.melt, aes(x = log10(intron_size), y = density, group = variable)) +
    geom_point(cex = 1, aes(colour = variable, shape = variable)) +
    geom_smooth(method = 'lm', aes(col = variable), show.legend = F, se = FALSE) +
    scale_color_manual(name = "", labels = legend_labs, values = c("RoyalBlue", "black")) +
    scale_shape_manual(name = "", labels = legend_labs, values = c(4, 16)) +
    scale_y_continuous(limits = c(0, 0.6)) +
    labs(x = "log10 intron size", y = "Density", title = title) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5),
      legend.position = c(0.85,0.9),
      legend.title = element_blank(),
      legend.background = element_blank(),
      legend.key = element_blank(),
    ) +
    guides(shape = guide_legend(override.aes = list(size = 2)))
}


fischers_z <- function(x, y1, y2) {
  # create linear models for both
  lm1 = lm(y1 ~ log10(x))
  lm2 = lm(y2 ~ log10(x))
  # get summaries
  summary1 <- summary(lm1)
  summary2 <- summary(lm2)
  # get coefficients
  coeffs1 <- summary1$coefficients
  coeffs2 <- summary2$coefficients
  # get slopes
  slope1 <- coeffs1[2,1]
  slope2 <- coeffs2[2,1]
  # get std error
  se1 <- coeffs1[2,2]
  se2 <- coeffs2[2,2]
  # calculate z
  z <- abs((slope1 - slope2)/(sqrt(se1^2 + se2^2)))
  return(z)
}


p_from_z <- function(z) {
  #  calculate p value
  pval <- 2*pnorm(z, lower.tail=FALSE)
  return(pval)
}

# for the stop / non stop eses separately

stop_nonstop_correlations = function(set_types, types, ese_sets, output_name) {
  blank_dataframe = data.frame(
    motif_set = character(),
    ese_rho = double(),
    ese_p = double(),
    stop_ese_rho = double(),
    stop_ese_p = double(),
    non_stop_ese_rho = double(),
    non_stop_ese_p = double(),
    fischers_z = double(),
    diff_p = double()
  )
  
  for (set_type in set_types) {
    cors = c()
    for (type in types) {
      correlations = blank_dataframe
      for (set_name in ese_sets) {
        print(set_name)
        data = read.csv(paste("clean_run/tests/ese_densities/", set_name, "_", set_type, "_ese_densities", type, ".csv", sep = ""), head = T)
        data = data[log10(data$intron_size) != 0,]
        
        cor = cor.test(log10(data$intron_size), data$ese_density, method = "spearman")
        cor1 = cor.test(log10(data$intron_size), data$stop_ese_density, method = "spearman")
        cor2 = cor.test(log10(data$intron_size), data$non_stop_ese_density, method = "spearman")
        z = fischers_z(log10(data$intron_size), data$stop_ese_density, data$non_stop_ese_density)
        
        cor_output = data.frame(
          motif_set = paste(set_name, type, sep = ""),
          ese_rho = cor$estimate,
          ese_p = cor$p.value,
          stop_ese_rho = cor1$estimate,
          stop_ese_p = cor1$p.value,
          non_stop_ese_rho = cor2$estimate,
          non_stop_ese_p = cor2$p.value,
          fischers_z = z,
          diff_p = p_from_z(z)
        )
        correlations = rbind(correlations, cor_output)
      }
      cors = rbind(cors, correlations)
      cors[nrow(cors)+1,] <- NA
    }
    output_file = paste("clean_run/tests/ese_densities/", output_name, "_", set_type, "_intron_size_ese_density_stop_non_stop_correlations.csv", sep = "")
    print(output_file)
    write.csv(cors, file = output_file, row.names = F, na = "")
    
  }
}
####

set_types = c("pc", "lincrna")
ese_sets = c("int3")
types = c("", "_flanks")
stop_nonstop_correlations(set_types, types, ese_sets, "processed_output")

pc = read.csv(paste("clean_run/tests/ese_densities/int3_pc_ese_densities.csv", sep = ""), head = T)
pc_plot = ggarrange(
  all_ese_plot(pc),
  split_ese_plot(pc),
  ncol = 2,
  labels = c("A", "B")
)
ggsave(plot = pc_plot, file = "clean_run/plots/pc_intron_size_ese_density_plot_flanks.pdf", width = 9, height = 4)

lincrna = read.csv(paste("clean_run/tests/ese_densities/int3_linc_ese_densities_all_seq.csv", sep = ""), head = T)
lincrna_plot = ggarrange(
  all_ese_plot(lincrna),
  split_ese_plot(lincrna),
  ncol = 2,
  labels = c("A", "B")
)
ggsave(plot = lincrna_plot, file = "clean_run/plots/lincrna_intron_size_ese_density_plot.pdf", width = 9, height = 4)
