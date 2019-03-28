
run = function(filepath, set) {
  
  file = read.csv(filepath, head = T)
  
  real = file[file$id == "real",]
  sims = file[file$id != "real",]
  
  output = c(set)

  output = c(output, real$ese_rate, median(sims$ese_rate), (nrow(sims[sims$ese_rate <= real$ese_rate,]) +1) / (nrow(sims) + 1))
  output = c(output, real$non_ese_rate, median(sims$non_ese_rate), (nrow(sims[sims$non_ese_rate <= real$non_ese_rate,]) +1) / (nrow(sims) + 1))
  output = c(output, real$ese_stop_rate, median(sims$ese_stop_rate), (nrow(sims[sims$ese_stop_rate <= real$ese_stop_rate,]) +1) / (nrow(sims) + 1))
  output = c(output, real$ese_non_stop_rate, median(sims$ese_non_stop_rate), (nrow(sims[sims$ese_non_stop_rate <= real$ese_non_stop_rate,]) +1) / (nrow(sims) + 1))
  output = c(output, real$non_ese_stop_rate, median(sims$non_ese_stop_rate), (nrow(sims[sims$non_ese_stop_rate <= real$non_ese_stop_rate,]) +1) / (nrow(sims) + 1))
  output = c(output, real$non_ese_non_stop_rate, median(sims$non_ese_non_stop_rate), (nrow(sims[sims$non_ese_non_stop_rate <= real$non_ese_non_stop_rate,]) +1) / (nrow(sims) + 1))
  output = c(output, real$relative_ese_diff, median(sims$relative_ese_diff), (nrow(sims[sims$relative_ese_diff >= real$relative_ese_diff,]) +1) / (nrow(sims) + 1))
  output = c(output, real$relative_non_ese_diff, median(sims$relative_non_ese_diff), (nrow(sims[sims$relative_non_ese_diff >= real$relative_non_ese_diff,]) +1) / (nrow(sims) + 1))
  output = c(output, real$log_diff_ratio, median(sims$log_diff_ratio), (nrow(sims[sims$log_diff_ratio <= real$log_diff_ratio,]) +1) / (nrow(sims) + 1))
  
  return(output)
}

output1 = run("clean_run/tests/lincrna/substitution_rates/lincrna_int3_substitution_rates_motif.csv", set = "lincrna_int3") 
output2 = run("clean_run/tests/lincrna/substitution_rates/lincrna_RESCUE_substitution_rates_motif.csv", set = "lincrna_RESCUE") 
output3 = run("clean_run/tests/lincrna/substitution_rates/pc_int3_substitution_rates_motif.csv", set = "protein_coding_int3") 

headers1 = c("",
  "ESE", "", "",
  "Non-ESE", "", "",
  "ESE stops", "", "",
  "ESE non-stops", "", "",
  "Non ESE stops", "", "",
  "Non ESE non-stops", "", "",
  "ESE diff", "", "",
  "Non ESE diff", "", "",
  "LogDiffRatio", "", ""
)
headers2 = c("set",
  "real", "median sim", "p",
  "real", "median sim", "p",
  "real", "median sim", "p",
  "real", "median sim", "p",
  "real", "median sim", "p",
  "real", "median sim", "p",
  "real", "median sim", "p",
  "real", "median sim", "p",
  "real", "median sim", "p"
)


data = matrix(0, ncol = length(headers1), nrow = 0)
colnames(data) = headers1
data = rbind(data, headers2)
data = rbind(data, output1)
data = rbind(data, output2)
data = rbind(data, output3)

write.table(data, file = "temp_files/_output1.csv",  sep = ",", row.names = F)