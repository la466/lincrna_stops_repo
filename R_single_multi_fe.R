library(reshape2)

single = read.csv("temp_files/gc_output_single.csv", head = T)
multi = read.csv("temp_files/gc_output_multi.csv", head = T)

single$group = "single"
multi$group = "multi"

single = single[is.finite(single$fe),]
multi = multi[is.finite(multi$fe),]

nrow(single)
nrow(multi)

wilcox.test(single$density, multi$density)
median(single$density)
median(multi$density)

wilcox.test(single$gc, multi$gc)
median(single$gc)
median(multi$gc)

wilcox.test(single$fe, multi$fe, alternative = "g")
median(single$fe)
median(multi$fe)

cor.test(single$gc, single$fe, method = "spearman")
cor.test(multi$gc, multi$fe, method = "spearman")

# join data and perform loess
data = rbind(single, multi)
data = data[is.finite(data$fe),]
nrow(data)

median(data$fe[data$group == "single"])
median(data$fe[data$group == "multi"])

data.lo <- loess(data$fe~data$gc)
plot(data$gc, data$fe,pch=19,cex=0.1)

wilcox.test(data.lo$residuals[data$group == "single"], data.lo$residuals[data$group == "multi"], alternative = "g")
median(data.lo$residuals[data$group == "single"])
median(data.lo$residuals[data$group == "multi"])


