ds_file = read.csv("clean_run/tests/ese_ds/ese_ds1.csv", head = T)

empirical_p = function(data, column, alternative = "less") {
  real = data[data$id == "real",]
  sims = data[data$id != "real",]
  
  if (alternative == "less") {
    p = (nrow(sims[sims[[column]] <= real[[column]],]) + 1) / (nrow(sims) + 1)
  } else {
    p = (nrow(sims[sims[[column]] >= real[[column]],]) + 1) / (nrow(sims) + 1)
  }
  return(p)
}

ds_file$ese_ds[ds_file$id == "real"]
median(ds_file$ese_ds[ds_file$id != "real"])
ese_p = empirical_p(ds_file, "ese_ds")
ese_p
non_ese_p = empirical_p(ds_file, "non_ese_ds")

ds_file$ese_stop_ds[ds_file$id == "real"]
ds_file$ese_non_stop_ds[ds_file$id == "real"]

empirical_p(ds_file, "ese_non_stop_ds")
empirical_p(ds_file, "ese_stop_ds")

ds_mutation_file = read.csv("clean_run/tests/ese_ds/ese_ds_mutation.csv", head = T)
head(ds_mutation_file)

ds_mutation_file[ds_mutation_file$id == "real",]
empirical_p(ds_mutation_file, "one_away_ds")
empirical_p(ds_mutation_file, "others_ds")

real = ds_mutation_file[ds_mutation_file$id == "real",]
real$one_away_ds
real$others_ds


file = read.csv("clean_run/tests/ese_ds/codon_ds_stats.csv", head=T)

head(file)

file$diff = file$stops_ese_ds = file$non_stops_ese_ds
real = file[file$id == "real",]
real
sims = file[file$id != "real",]

(nrow(sims[sims$all_ese_ds <= real$all_ese_ds,]) + 1) / (nrow(sims) + 1)
(nrow(sims[sims$non_ese_ds <= real$non_ese_ds,]) + 1) / (nrow(sims) + 1)
(nrow(sims[sims$stops_ese_ds <= real$stops_ese_ds,]) + 1) / (nrow(sims) + 1)
(nrow(sims[sims$non_stops_ese_ds <= real$non_stops_ese_ds,]) + 1) / (nrow(sims) + 1)
(nrow(sims[sims$diff <= real$diff,]) + 1) / (nrow(sims) + 1)

real
head(sims)
real


nrow(sims)

nrow(sims[sims$all_ese_ds <= real$all_ese_ds,])
nrow(sims[sims$non_ese_ds <= real$non_ese_ds,])
nrow(sims[sims$stops_ese_ds <= real$stops_ese_ds,])
nrow(sims[sims$non_stops_ese_ds <= real$non_stops_ese_ds,])


sims
nd = (real$all_ese_ds - mean(sims$all_ese_ds)) / mean(sims$all_ese_ds)

real$diff = (real$stops_ese_ds / real$non_stops_ese_ds)
sims$diff = (sims$stops_ese_ds / sims$non_stops_ese_ds) * (1 - abs(nd))

nrow(sims[sims$diff > 0 & sims$diff >= real$diff,])
real
sims

h



emperical_p <- function(data, id_col, id_differentiator, test_col) {
  query_val = data[data[[id_col]] == id_differentiator,]
  test_vals = data[data[[id_col]] != id_differentiator,]
  p <- (nrow(test_vals[test_vals[[test_col]] <= query_val[[test_col]],]) + 1) / (nrow(test_vals) + 1)
  return(p)
}

nd <- function(data, id_col, id_differentiator, test_col) {
  real = data[data[[id_col]] == id_differentiator,]
  sims = data[data[[id_col]] != id_differentiator,]
  nd = (real[[test_col]] - mean(real[[test_col]])) / mean(real[[test_col]])
  return (nd)
}

head(file)

emperical_p(file, "sim_id", "real", "all_ese_ds")

median(sims$all_ese_ds)
real$all_ese_ds

median(sims$stops_ese_ds)
median(sims$non_stops_ese_ds)

real$all_ese_stop_hits_ds
median(sims$all_ese_stop_hits_ds)
emperical_p(file, "sim_id", "real", "all_ese_stop_hits_ds")
nd(file, "sim_id", "real", "all_ese_stop_hits_ds")

(real$all_ese_stop_hits_ds - mean(sims$all_ese_stop_hits_ds)) / mean(sims$all_ese_stop_hits_ds)
(real$all_ese_non_stop_hits_ds - mean(sims$all_ese_non_stop_hits_ds)) / mean(sims$all_ese_non_stop_hits_ds)



real$all_ese_non_stop_hits_ds
median(sims$all_ese_non_stop_hits_ds)
emperical_p(file, "sim_id", "real", "all_ese_non_stop_hits_ds")

real$ratio
median(sims$ratio)
head(sims)



file = read.csv("clean_run/tests/ese_ds/codon_ds_stats_mutation.csv", head=T)
head(file)
file$diff = file$stop_mutation_ds - file$non_stops_mutation_ds
real = file[file$id == "real",]
sims = file[file$id != "real",]

nrow(sims[sims$diff < 0 & sims$diff <= real$diff,])
nrow(sims)
plot(sims$stop_mutation_ds, sims$non_stops_mutation_ds)
abline(lm(sims$non_stops_mutation_ds ~ sims$stop_mutation_ds), col="red")
lm = lm(sims$non_stops_mutation_ds ~ sims$stop_mutation_ds)

expected = lm$coefficients[1] + (lm$coefficients[1] * real$stop_mutation_ds)
expected
real$non_stops_mutation_ds



file = read.csv("clean_run/tests/ese_ds/ese_ds_mutation1.csv", head = T)
head(file)
nrow(file)

file$relative = file$one_away_ds / file$others_ds
real = file[file$id == "real",]
sims = file[file$id != "real",]

real$one_away_ds

mean(real$relative)

nrow(sims[sims$one_away_ds <= real$one_away_ds,]) / nrow(sims)
nrow(sims[sims$others_ds <= real$others_ds,]) / nrow(sims)

hist(sims$others_ds)
abline(v = real$others_ds)

head(file)
