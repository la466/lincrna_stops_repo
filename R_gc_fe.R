# file = read.csv("clean_run/motif_tests/int3_densities.csv", head = T)
# gc_matched = file$codons[file$gc_content == 2/9 & file$codons != "TAA_TAG_TGA"]
# write.table(gc_matched, "clean_run/gc_matched_sets.bed", row.names = F, quote= F)

file = "clean_run/gc_fes.csv"
data = read.csv(file, head = T)

non_stops = data[data$codon_set != "TAA_TAG_TGA",]
stops = data[data$codon_set == "TAA_TAG_TGA",]

binom.test(nrow(non_stops[non_stops$fe >= stops$fe,]), nrow(non_stops), p = 0.5, alternative = "g")
binom.test(nrow(non_stops[non_stops$fe >= stops$fe & non_stops$included_stops == 0,]), nrow(non_stops[non_stops$included_stops == 0,]), p = 0.5, alternative = "g")