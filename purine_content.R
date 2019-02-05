file = read.csv("clean_run/purine_content.csv", head=T)
nrow(file)
head(file)

qqnorm(file$exon_purine_content)
qqnorm(file$intron_purine_content)

wilcox.test(file$exon_purine_content, file$intron_purine_content, paired = T)
median(file$exon_purine_content)
median(file$intron_purine_content)

nrow(file)

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


file = read.csv("clean_run/purine_content_no_eses.csv", head = T)
head(file)
nrow(file)
wilcox.test(file$exon, file$intron, paired = T)
median(file$exon)
median(file$intron)

nrow(file)

file = read.csv("clean_run/tests/purine_content/exon_intron_purine_content_with_ese.csv", head = T)

real = file[file$id == "real",]
sims = file[file$id != "real",]

head(file)
real
