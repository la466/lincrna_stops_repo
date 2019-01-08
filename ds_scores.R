file = read.csv("clean_run/tests/ese_ds/codon_ds_stats.csv", head=T)

file$ratio = file$all_ese_stop_hits_ds / file$all_ese_non_stop_hits_ds
real = file[file$sim_id == "real",]
sims = file[file$sim_id != "real",]

emperical_p <- function(data, id_col, id_differentiator, test_col) {
  query_val = data[data[[id_col]] == id_differentiator,]
  test_vals = data[data[[id_col]] != id_differentiator,]
  p <- (nrow(test_vals[test_vals[[test_col]] <= query_val[[test_col]],]) + 1) / (nrow(test_vals) + 1)
  return(p)
}


emperical_p(file, "sim_id", "real", "all_ese_ds")

median(sims$all_ese_ds)
real$all_ese_ds

median(sims$stops_ese_ds)
median(sims$non_stops_ese_ds)

real$all_ese_stop_hits_ds
median(sims$all_ese_stop_hits_ds)
emperical_p(file, "sim_id", "real", "all_ese_stop_hits_ds")

real$all_ese_non_stop_hits_ds
median(sims$all_ese_non_stop_hits_ds)
emperical_p(file, "sim_id", "real", "all_ese_non_stop_hits_ds")

real$ratio
median(sims$ratio)
head(sims)
