file = read.csv("clean_run/purine_content.csv", head=T)

head(file)

shapiro.test(file$exon_purine_content)

qqnorm(file$exon_purine_content)
qqnorm(file$intron_purine_content)

wilcox.test(file$exon_purine_content, file$intron_purine_content, paired = T)
median(file$exon_purine_content)
median(file$intron_purine_content)

cor.test(file$exon_purine_content, file$intron_purine_content, method = "spearman")

tail(file)
