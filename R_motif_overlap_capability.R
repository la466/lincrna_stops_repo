# overlap capability
overlap_capability_results = function(file) {
  data = read.csv(file, head = T)
  real = data[data$id == "real",]
  sims = data[data$id != "real",]

  stop_nd = (real$stop_overlaps - mean(sims$stop_overlaps)) / mean(sims$stop_overlaps)
  non_stop_nd = (real$non_stop_overlaps - mean(sims$non_stop_overlaps)) / mean(sims$non_stop_overlaps)
  diff_nd = (real$diff - mean(sims$diff)) / mean(sims$diff)
  
  stop_p = (nrow(sims[sims$stop_overlaps <= real$stop_overlaps,]) + 1) / (nrow(sims) + 1)
  non_stop_p = (nrow(sims[sims$non_stop_overlaps >= real$non_stop_overlaps,]) + 1) / (nrow(sims) + 1)
  diff_p = (nrow(sims[sims$diff >= real$diff,]) + 1) / (nrow(sims) + 1)
  
  output = data.frame(
    stop_nd,
    non_stop_nd,
    diff_nd,
    stop_p,
    non_stop_p,
    diff_p
  )
  return(output)
}


####

overlap_capability_results("clean_run/tests/ese_overlaps/int3_motif_overlap_capability.csv")
overlap_capability_results("clean_run/tests/ese_overlaps/RESCUE_motif_overlap_capability.csv")
