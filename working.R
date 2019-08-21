library(data.table)

file = "temp_data/compare_densities.csv"

data = read.csv(file, head = F)

data.t = transpose(data, fill = 0)
head(data.t)
colnames(data.t) <- as.character(unlist(data.t[1,]))
data.t = as.data.frame(data.t[-1, ])
data.t <- sapply( data.t, as.numeric )
head(data.t)
data = data.frame(data.t)

head(data)

boxplot(data$single_density, data$multi_density)
boxplot(data$single_nd, data$multi_nd)

wilcox.test(data$single, data$multi)

(median(data$single, na.rm = T) - median(data$multi, na.rm = T)) / median(data$single, na.rm = T)


wilcox.test(data$single_nd, data$multi_nd)

data1 = data[!is.name(data$single_density),]
nrow(data1)

nrow(data)
tail(data)

data$multi_density
