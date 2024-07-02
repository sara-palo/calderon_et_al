invisible(lapply(list(
  "stringr",
  "biomaRt",
  "ggplot2",
  "pheatmap",
  "DESeq2",
  "readxl"
), FUN = function(x) {
  suppressPackageStartupMessages(library(x, character.only = T))
}))

set.seed(12345)

external_data <- 'BALL-1988S-HTSeq'

# load genes of interest (GOI)
goi_up <- read_xlsx('human_data/Leuk signature cluster 5 human.xlsx', sheet = 1, col_names = 'gene')
goi_dn <- read_xlsx('human_data/Leuk signature cluster 5 human.xlsx', sheet = 2, col_names = 'gene')

goi <- c(goi_up$gene, goi_dn$gene)

sample_table <- read.table(paste0(external_data, '/subtypes.tsv'),
                           header = TRUE, dec = ",", sep = "\t", quote = '""')

colums_oi <- c("patient", "fusion", "age", "gender", "primary.subtype", "RNA.seq.library")
sample_table <- sample_table[colums_oi]

fusions_oi <- c("KMT2A-AFF1", "ETV6-RUNX1", "BCR-ABL1", "NoFusion") 
sample_table <- sample_table[sample_table$fusion %in% fusions_oi, ]

subtype_oi <- c("KMT2A", "ETV6-RUNX1", "Ph", "High hyperdiploid", "Low hypodiploid") # 
sample_table <- sample_table[sample_table$primary.subtype %in% subtype_oi, ]

sample_table <- subset.data.frame(sample_table, !age == "NA")

sample_table$age_group <- cut(sample_table$age, c(
  0, 1, 16, 40,
  max(sample_table$age)
))
levels(sample_table$age_group) <- c("0-1", "1-16", "16-40", ">40")

for (i in 1:dim(sample_table)[1]){
  if(sample_table[i, 'fusion'] == 'NoFusion'){
    sample_table[i, 'fusion'] <- sample_table[i, 'primary.subtype']
  }
}

# order sample table by type
level_order <-  c("KMT2A-AFF1", "ETV6-RUNX1", "BCR-ABL1", "High hyperdiploid", "Low hypodiploid")
sample_table$fusion <- factor(sample_table$fusion, levels = level_order)

# get file names
file_names <- c()
for (i in 1:dim(sample_table)[1]) {
  file_names[i] <- list.files(
    path = external_data,
    pattern = as.character(sample_table$patient[i]),
    ignore.case = TRUE, full.names = TRUE
  )[1]
  # hard-coded to only take the first sample,
  # works for this one but might not for others
}
sample_table$file_name <- file_names

# Load samples into DEseq2
de_sample_table <- data.frame(
  sampleNames = sample_table$patient,
  fileName = as.character(sample_table$file_name),
  age_group = sample_table$age_group,
  type = sample_table$fusion,
  age = sample_table$age, 
  gender = sample_table$gender,
  RNA_lib = sample_table$RNA.seq.library
)

counts <- DESeqDataSetFromHTSeqCount(
  sampleTable = de_sample_table, design = ~RNA_lib
)
counts <- DESeq(counts)

vsd <- vst(counts, blind = FALSE)
plotPCA(vsd, "RNA_lib") + theme_bw()

assay(vsd) <- limma::removeBatchEffect(assay(vsd), vsd$RNA_lib)
plotPCA(vsd, "RNA_lib") + theme_bw()

ensembl <- useMart(
  biomart = "ENSEMBL_MART_ENSEMBL", host = "grch37.ensembl.org",
  path = "/biomart/martservice", dataset = "hsapiens_gene_ensembl"
)
ids_up <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters = "hgnc_symbol",
  values = goi_up$gene,
  mart = ensembl
)

ids2 <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters = "hgnc_symbol",
  values = c('HMGA2', 'LIN28B'),
  mart = ensembl
)

intgroup <- c("type", "age", "gender", "RNA_lib", "age_group") #

select <- rownames(vsd) %in% ids_up$ensembl_gene_id

length(select)
table(select)

ensm_to_sym <- data.frame(
  ensblID = rownames(assay(vsd)[select, ]),
  gene_sym = getBM(
    attributes = "hgnc_symbol", filters = "ensembl_gene_id",
    values = rownames(assay(vsd)[select, ]), mart = ensembl
  )
)

intgroup_df <- as.data.frame(colData(vsd)[, intgroup, drop = FALSE])

vsd_df <- as.data.frame(t(assay(vsd)[select, ]))
vsd_df <- sapply(vsd_df, function(vsd_df) {
  (vsd_df -
     mean(vsd_df)) /
    sd(vsd_df)
})
means <- rowMeans(vsd_df)
vsd_df <- cbind.data.frame(intgroup_df, means)

p1 <- ggplot(data = vsd_df, aes_string(
  y = "means", x = "type",
  fill = "type"
)) +
  geom_boxplot() +
  coord_fixed() +
  theme_bw() +
  theme(aspect.ratio = 2) +
  ylab("Z-scored mean expression") +
  xlab("") + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
print(p1)
ggsave(plot = p1, 'figs/upreg_sign_subtypes.pdf')

res_aov <- aov(means ~ type, data = vsd_df)
res <- TukeyHSD(res_aov)
res
write.csv(res$type, 'figs/upreg_sign_subtypes_stats.csv')


select <- rownames(vsd) %in% ids2[ids2$hgnc_symbol == 'HMGA2', 'ensembl_gene_id']
table(select)

HMGA2 <- assay(vsd)[select, ]
HMGA2 <- (HMGA2 - mean(HMGA2)) / sd(HMGA2)
intgroup_df$HMGA2 <- HMGA2

select <- rownames(vsd) %in% ids2[ids2$hgnc_symbol == 'LIN28B', 'ensembl_gene_id']
table(select)

LIN28B <- assay(vsd)[select, ]
LIN28B <- (LIN28B - mean(LIN28B)) / sd(LIN28B)
intgroup_df$LIN28B <- LIN28B


ggplot(data = intgroup_df, aes_string(
  y = "HMGA2", x = "type",
  fill = "type"
)) +
  geom_boxplot() +
  coord_fixed() +
  theme_bw() +
  theme(aspect.ratio = 2) +
  ylab("Z-scored expression") +
  xlab("") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ggsave('figs/HMGA2.pdf')

res_aov <- aov(HMGA2 ~ type, data = intgroup_df)
res <- TukeyHSD(res_aov)
res
write.csv(res$type, 'figs/HMGA2_stats.csv')


ggplot(data = intgroup_df, aes_string(
  y = "LIN28B", x = "type",
  fill = "type"
)) +
  geom_boxplot() +
  coord_fixed() +
  theme_bw() +
  theme(aspect.ratio = 2) +
  ylab("Z-scored expression") +
  xlab("") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ggsave('figs/LIN28B.pdf')


res_aov <- aov(LIN28B ~ type, data = intgroup_df)
res <- TukeyHSD(res_aov)
res
write.csv(res$type, 'figs/LIN28B_stats.csv')

