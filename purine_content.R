file = read.csv("clean_run/purine_content.csv", head=T)

head(file)

qqnorm(file$exon_purine_content)
qqnorm(file$intron_purine_content)

wilcox.test(file$exon_purine_content, file$intron_purine_content, paired = T)
median(file$exon_purine_content)
median(file$intron_purine_content)

cor.test(file$exon_purine_content, file$intron_purine_content, method = "spearman")

wilcox.test(file$exon_core_purine_content, file$intron_core_purine_content, paired = T)
median(file$exon_core_purine_content)
median(file$intron_core_purine_content)


file = read.csv("temp_files/random_purine_content.csv", head = T)
head(file)

real = file[file$id == "real",]
sims = file[file$id != "real",]

btest = binom.test(nrow(sims[sims$pvalue < 0.05,]), nrow(sims), alternative = "g")
btest

exons = (nrow(sims[sims$median_exon >= real$median_exon,]) + 1) / (nrow(sims) + 1)
intron = (nrow(sims[sims$median_intron <= real$median_intron,]) + 1) / (nrow(sims) + 1)


file = read.csv("temp_files/purine_content_no_eses.csv", head = T)
real = file[file$id == "real",]
sims = file[file$id != "real",]


