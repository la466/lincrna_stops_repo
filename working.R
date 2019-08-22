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

file = "temp_data/compare_densities1.csv"



data = read.csv(file, head = F)
data = transpose_data(data)

head(data)

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


tail(data)

(median(data$single, na.rm = T) - median(data$multi, na.rm = T)) / median(data$single, na.rm = T)


wilcox.test(data$single_nd, data$multi_nd)


nrow(data[!is.na(data$single_density),])

data1 = data[!is.name(data$single_density),]
nrow(data1)

nrow(data)
tail(data)

data$multi_density
