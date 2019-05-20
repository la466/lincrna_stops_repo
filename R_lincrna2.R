file = read.csv("clean_run/lincrna/tests/stop_densities.csv", head = T)

wilcox.test(file$exon, file$intron, paired = T)
wilcox.test(file$exon_gc, file$intron_gc, paired = T)
median(file$exon)
median(file$intron)
median(file$exon_gc)
median(file$intron_gc)


file = read.csv("clean_run/lincrna/tests/stop_densities_nd.csv", head = T)
wilcox.test(file$exon, file$intron, paired = T)
nrow(file)
