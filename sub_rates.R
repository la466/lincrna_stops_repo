
run = function() {
  
  file = read.csv("clean_run/tests/lincrna/substitution_rates/lincrna_int3_substitution_rates_motif.csv", head = T)
  
  real = file[file$id == "real",]
  sims = file[file$id != "real",]

  # ese rate
  print("ESE")
  print(real$ese_rate)
  print(median(sims$ese_rate))
  print((nrow(sims[sims$ese_rate <= real$ese_rate,]) +1) / (nrow(sims) + 1))
  
  # non ese rate
  print("Non ESE")
  print(real$non_ese_rate)
  print(median(sims$non_ese_rate))
  print((nrow(sims[sims$non_ese_rate <= real$non_ese_rate,]) +1) / (nrow(sims) + 1))
  # ese stop 
  print("ESE stop")
  print(real$ese_stop_rate)
  print(median(sims$ese_stop_rate))
  print((nrow(sims[sims$ese_stop_rate <= real$ese_stop_rate,]) +1) / (nrow(sims) + 1))
  # ese non stop 
  print("ESE non-stop")
  print(real$ese_non_stop_rate)
  print(median(sims$ese_non_stop_rate))
  print((nrow(sims[sims$ese_non_stop_rate <= real$ese_non_stop_rate,]) +1) / (nrow(sims) + 1))
  # non ese stop 
  print("non-ESE stop")
  print(real$non_ese_stop_rate)
  print(median(sims$non_ese_stop_rate))
  print((nrow(sims[sims$non_ese_stop_rate <= real$non_ese_stop_rate,]) +1) / (nrow(sims) + 1))# non ese stop 
  print("non-ESE non stop")
  print(real$non_ese_non_stop_rate)
  print(median(sims$non_ese_non_stop_rate))
  print((nrow(sims[sims$non_ese_non_stop_rate <= real$non_ese_non_stop_rate,]) +1) / (nrow(sims) + 1))
  # ESE diff
  print("ESE diff")
  print(real$relative_ese_diff)
  print(median(sims$relative_ese_diff))
  print((nrow(sims[sims$relative_ese_diff <= real$relative_ese_diff,]) +1) / (nrow(sims) + 1))
  # non ESE diff
  print("non-ESE diff")
  print(real$relative_non_ese_diff)
  print(median(sims$relative_non_ese_diff))
  print((nrow(sims[sims$relative_non_ese_diff <= real$relative_non_ese_diff,]) +1) / (nrow(sims) + 1))
  # log diff ratio
  print("log diff ratio")
  print(real$log_diff_ratio)
  print(median(sims$log_diff_ratio))
  print((nrow(sims[sims$log_diff_ratio <= real$log_diff_ratio,]) +1) / (nrow(sims) + 1))
  
}
run() 



