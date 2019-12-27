library(reshape2)

single = read.csv("clean_run/gc_output_single.csv", head = T)
multi = read.csv("clean_run/gc_output_multi.csv", head = T)

single = single[is.finite(single$fe),]
multi = multi[is.finite(multi$fe),]

nrow(single)
nrow(multi)

single$group = "single"
multi$group = "multi"
single$group <- as.factor(single$group)
multi$group <- as.factor(multi$group)

median(single$density)
median(multi$density)
wilcox.test(single$density, multi$density, alternative = "g")

median(single$gc)
median(multi$gc)
wilcox.test(single$gc, multi$gc)


lm1 <- lm(single$fe~single$gc)
summary(lm1)
lm2 <- lm(multi$fe~multi$gc)
summary(lm2)

co1 <- summary(lm1)$coefficients[2,1]
se1 <- summary(lm1)$coefficients[2,2]
co2 <- summary(lm2)$coefficients[2,1]
se2 <- summary(lm2)$coefficients[2,2]

# Z test: prop A second ≠ prop A codons
z <- (co2 - co1) / sqrt((se2**2) + (se1**2))
pvalue2sided <- 2*pnorm(-abs(z))
pvalue2sided

# join data and perform loess
data = rbind(single, multi)

lm_gc = lm(data$fe ~ data$gc)
summary(lm_gc)

median(single$fe)
median(multi$fe)
wilcox.test(single$fe, multi$fe, alternative = "g")

median(single$ese_density)
median(multi$ese_density)
wilcox.test(single$ese_density, multi$ese_density, alternative = "l")

# plot(single$gc, single$fe, xlim = c(0, 1), cex = 0.2, col = "red", xlab = "GC", ylab = "FE")
# abline(lm(single$fe~single$gc), col = "red")
# points(multi$gc, multi$fe, pch = 17, cex = 0.2, col = "blue")
# abline(lm(multi$fe~multi$gc), col = "blue")
# abline(lm(data$fe~data$gc), lwd = 3)
# 
# cor.test(single$gc, single$fe, method = 'spearman')
# cor.test(multi$gc, multi$fe, method = 'spearman')
