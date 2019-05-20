file = read.csv("clean_run/intron_density/motif_sets_intron_density.csv", head = T)

nrow(file)

sim_ids = colnames(file)[3:length(colnames(file))]
real = file$real
length(sim_ids)

diff_props = c()
p_values = c()
for (sim_id in sim_ids) {
  print(sim_id)
  sim = file[[sim_id]]  
  diff = real - sim
  # get the proportion of genes where the simulated density is greater than the real density, i.e 
  # the difference is negative
  diff_prop = sum(diff <= 0) / length(diff)
  test = wilcox.test(real, sim, paired = T)
  
  diff_props = c(diff_props, diff_prop)
  p_values = c(p_values, test$p.value)
}

sum(diff_props > 0.5)
binom.test(sum(diff_props > 0.5), length(diff_props), alternative = "g")

adjusted_ps = p.adjust(p_values, method = "bonferroni")

significantly_greater = c()
for (i in seq(1:length(diff_props))) {
  diff = diff_props[i]
  p = p_values[i]
  if(diff > 0.5 & p < 0.05) {
    significantly_greater = c(significantly_greater, 1)
  }
}

sum(significantly_greater)
binom.test(sum(significantly_greater), length(adjusted_ps), p = 0.05, alternative = "g")


library(ggplot2)
plot = ggplot() +
  geom_histogram(aes(x = diff_props), bins = 20, fill = "RoyalBlue", col = "black") + 
  labs(x = "Proportion of sequences with real ESE density less than simulant set density", y = "Frequency")
ggsave("clean_run/plots/intron_ese_density_simulant_proportion.pdf", width = 10, height = 7)

empirical_p = function(data) {
  real = data[2]
  sims = data[3:length(data)]
  p = (sum(sims <= real) + 1) / (length(sims) + 1)
  return(p)
}

within_gene = apply(file, 1, empirical_p)
sum(within_gene < 0.05)
binom.test(sum(within_gene < 0.05), length(within_gene), p = 0.05, alternative = "g")
within_gene = p.adjust(within_gene, method = "bonferroni")
sum(within_gene < 0.05)

