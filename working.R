library(data.table)

transpose_data = function(data) {
  data.t = transpose(data, fill = 0)
  head(data.t)
  colnames(data.t) <- as.character(unlist(data.t[1,]))
  data.t = as.data.frame(data.t[-1, ])
  data.t <- sapply( data.t, as.numeric )
  head(data.t)
  data = data.frame(data.t)
  return(data)
}

# file = "temp_data/compare_densities1.csv"
file = "temp_data/utr_densities.csv"

data = read.csv(file, head = F)
data = transpose_data(data)

head(data)
nrow(data[!is.na(data$single_density),])

boxplot(data$single_density, data$multi_density)
boxplot(data$single_nd, data$multi_nd, ylim = c(-1, 2))


head(data)

median(data$single_density, na.rm = T)
median(data$multi_density, na.rm = T)
median(data$single_nd, na.rm = T)
median(data$multi_nd, na.rm = T)

wilcox.test(data$single_nd, data$multi_nd)
wilcox.test(data$single_density, data$multi_density)

data

median(data$single_nd, na.rm = T)



file1 = "temp_data/utr_densities.csv"
data1 = read.csv(file1, head = T)
head(data1)

boxplot(data1$nd)

binom.test(nrow(data1[data1$nd < 0,]), nrow(data1), p = 0.5)
median(data1$nd, na.rm = T)

wilcox.test(data$multi_nd, data1$nd)
