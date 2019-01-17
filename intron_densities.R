file = read.csv("clean_run/intron_density/motif_sets_intron_density.csv", head = T)
sim_ids = colnames(file)[3:length(colnames(file))]

head(file)
real = file$real

diffs = c()
for (sim_id in sim_ids) {
  sim = file[[sim_id]]  
  diff = real - sim
  diff = sum(diff < 0) / length(diff)
  diffs = c(diffs, diff)
}

p = length(diffs)
binom.test(sum(diffs > 0.5), length(diffs), alternative = "g", p = 0.5)

