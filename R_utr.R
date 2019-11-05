utr_ese_densities = read.csv("clean_run/utr_ese_densities.csv", head = F)
utr_densities = transpose_data(utr_ese_densities)

wilcox.test(utr_densities$single, utr_densities$multi)
median(utr_densities$single, na.rm =  T)
median(utr_densities$multi, na.rm =  T)
