file = read.csv("clean_run/intron_density/motif_sets_intron_density.csv", head = T)
sim_ids = colnames(file)[3:length(colnames(file))]

head(file)
real = file$real

sim_ids = sim_ids[1:500]

diffs = c()
ps = c()
for (sim_id in sim_ids) {
  sim = file[[sim_id]]  
  diff = real - sim
  diff = sum(diff < 0) / length(diff)
  diffs = c(diffs, diff)
  
  test = wilcox.test(real, sim, paired = T)
  ps = c(ps, test$p.value)
}


values = c()
for (i in seq(1, length(ps))) {
  diff = diffs[i]
  p = ps[i]
  print(diff)
  if (diff < 0.5 & p < 0.05) {
    values = c(values, p)
  } else {
    # values = c(values, 1-p)
  }
}
print(length(values))

hist(values)

p = length(diffs)


binom.test(sum(diffs > 0.5), length(diffs), alternative = "g", p = 0.5)



