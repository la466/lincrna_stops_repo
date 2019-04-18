file = read.csv("clean_run/tests/ese_overlaps/int3_motif_overlap_capability.csv", head = T)

real = file[file$id == "real",]
sims = file[file$id != "real",]

real

(nrow(sims[sims$stop_overlaps <= real$stop_overlaps,]) + 1) / (nrow(sims) + 1)
(nrow(sims[sims$non_stop_overlaps >= real$non_stop_overlaps,]) + 1) / (nrow(sims) + 1)
(nrow(sims[sims$diff >= real$diff,]) + 1) / (nrow(sims) + 1)


file = read.csv("clean_run/tests/ese_overlaps/RESCUE_motif_overlap_capability.csv", head = T)


real = file[file$id == "real",]
sims = file[file$id != "real",]

real

(nrow(sims[sims$stop_overlaps >= real$stop_overlaps,]) + 1) / (nrow(sims) + 1)
(nrow(sims[sims$non_stop_overlaps >= real$non_stop_overlaps,]) + 1) / (nrow(sims) + 1)
(nrow(sims[sims$diff >= real$diff,]) + 1) / (nrow(sims) + 1)
