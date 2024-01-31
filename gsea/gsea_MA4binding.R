library(fgsea)
library(data.table)
library(ggplot2)
library(readxl)
library(dplyr)
library(stringr)

set.seed(12345)

embryo_preproB <- read.table('gsea/Embryo_PreproB_Cre+_vs_Cre-_for_Preprob_Only_NoCut.txt', header = T)

kerry <- read_xls('gsea/kerry_2017_s1.xls', sheet = 1)
kerry_all <- kerry$`All targets, genes names (2597 unique)` %>% unique() %>% na.omit()

# manually fix the fact that the excel sheet was published with some gene names converted to dates 
kerry_all <- str_replace(kerry_all, pattern = '37500', replacement =  'SEPT2') %>% 
  str_replace(pattern = '38961', replacement = 'SEPT6') %>% 
  str_replace(pattern = '42248', replacement = 'SEPT15') %>%
  str_replace(pattern = '38412', replacement = 'MARCH5') %>%
  str_replace(pattern = '40787', replacement = 'SEPT11') %>%
  str_replace(pattern = '39142', replacement = 'MARCH7') %>%
  str_replace(pattern = '38777', replacement = 'MARCH6') %>%
  str_replace(pattern = '40057', replacement = 'SEPT9')

# the mouse_human_genes file was generated in December 2023 through:
# mouse_human_genes <- read.csv("http://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt", sep="\t")
# write.csv(mouse_human_genes, 'mouse_human_genes.csv')

mouse_human_genes <- read.csv('gsea/mouse_human_genes.csv', row.names = 1)

m2h <- data.frame()

#loop over the DEGs and retrieve corresponding human genes
for (gene in rownames(embryo_preproB)) {
  class_key <- (
    mouse_human_genes %>% filter(Symbol == gene &
                                   Common.Organism.Name == "mouse, laboratory")
  )[['DB.Class.Key']]
  if (!identical(class_key, integer(0))) {
    human_genes = (
      mouse_human_genes %>% filter(DB.Class.Key == class_key &
                                     Common.Organism.Name == "human")
    )[, "Symbol"]
    for (human_gene in human_genes) {
      m2h = rbind(m2h, c(gene, human_gene))
    }
  }
}

colnames(m2h) <- c('mouseGene', 'humanGene')

duplicated <- m2h$mouseGene[duplicated(m2h$mouseGene)]
m2h[m2h$mouseGene %in% duplicated, ]

#for the mouse genes with 2+ human IDs, if one is the same as the mouse ID but capitalized, remove the other IDs
for(dup in duplicated){
  tmp <- m2h[m2h$mouseGene %in% dup, ]
  if(toupper(dup) %in% tmp$humanGene){
    m2h <- subset(m2h, ! (mouseGene == dup & ! humanGene == toupper(dup)))
  }
}

duplicated <- m2h$mouseGene[duplicated(m2h$mouseGene)]
m2h[m2h$mouseGene %in% duplicated, ]

m2h[duplicated(m2h$mouseGene), 'humanGene'][m2h[duplicated(m2h$mouseGene), 'humanGene'] %in% kerry_all]

# For the ones where one of the gene IDs is in the Kerry gene set, keep that one and discard the rest
for(dup in duplicated){
  tmp <- m2h[m2h$mouseGene %in% dup, ]
  
  if(length(kerry_all[(kerry_all %in% tmp$humanGene)]) > 0){
    use <- data.frame(mouseGene = dup, humanGene = kerry_all[(kerry_all %in% tmp$humanGene)])
    m2h <- subset(m2h, ! (mouseGene == use$mouseGene & ! humanGene == use$humanGene))
  }
}

duplicated <- m2h$mouseGene[duplicated(m2h$mouseGene)]
table(m2h[m2h$mouseGene %in% duplicated, 'humanGene'] %in% kerry_all) #none of the remaining are in the Kerry gene set

embryo_preproB$humanGene <- m2h[match(rownames(embryo_preproB), m2h$mouseGene), 'humanGene']

# for the genes with no human counterpart, use the caps version of the mouse name
embryo_preproB[is.na(embryo_preproB$humanGene), 'humanGene'] <- toupper(rownames(embryo_preproB[is.na(embryo_preproB$humanGene),]))

# there are still some duplicated human gene names
embryo_preproB$humanGene[duplicated(embryo_preproB$humanGene)]
table(embryo_preproB$humanGene[duplicated(embryo_preproB$humanGene)] %in% kerry_all) # only 7 are part of the Kerry gene set

ranks_cre_pos <- embryo_preproB$avg_log2FC
names(ranks_cre_pos) <- embryo_preproB$humanGene

table(names(ranks_cre_pos) %in% kerry_all)

plotEnrichment(kerry_all,
               ranks_cre_pos) + labs(title = "Kerry SEM all MLL-AF4 binding")
ggsave('kerry_SEM_GSEA.pdf')

fgseaRes_kerry <- fgsea(pathways = list('kerry SEM' = kerry_all), 
                        stats = ranks_cre_pos)
write.csv(data.frame(fgseaRes_kerry)[,1:7], 'kerry_SEM_GSEA.csv')

write.csv(embryo_preproB, 'preproB_cre+_humanized_genes.csv')


