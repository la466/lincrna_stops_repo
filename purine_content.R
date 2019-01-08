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

wilcox.test(file$exon, file$intron, paired = T)
wilcox.test(file$sim_exon, file$sim_intron, paired = T)

