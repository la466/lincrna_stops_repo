file = read.csv("clean_run/tests/ese_overlaps/int3_motif_overlap_capability.csv", head = T)

real = file[file$id == "real",]
sims = file[file$id != "real",]

real

(nrow(sims[sims$stop_overlaps <= real$stop_overlaps,]) + 1) / (nrow(sims) + 1)
(nrow(sims[sims$non_stop_overlaps >= real$non_stop_overlaps,]) + 1) / (nrow(sims) + 1)
(nrow(sims[sims$diff >= real$diff,]) + 1) / (nrow(sims) + 1)


o = c(24454, 72943)
total = sum(o)

e = c(9457/91688*total, 82231/91688*total)
e

c = ((o-e)^2)/e
c = sum(c)
c
