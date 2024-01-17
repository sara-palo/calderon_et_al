library(tidyverse)
library(ggrepel)
library(ggplot2)
library(ggbreak) 
library(patchwork)

#EMBRYO VOLCANO PLOT

NewPPB_EM <- read.table("volcano_plots/Embryo_PreproB_Cre-_vs_Cre+.txt", header = TRUE)
NewPPB_EM

NewPPB_EM$NEWLOG <- NewPPB_EM$avg_log2FC*(-1)

NewPPB_EM$diffexpressed <- "Unchanged"
NewPPB_EM$diffexpressed[NewPPB_EM$NEWLOG >= 0.5 & NewPPB_EM$p_val_adj <= 0.05] <- "Up-regulated"
NewPPB_EM$diffexpressed[NewPPB_EM$NEWLOG <= -0.5 & NewPPB_EM$p_val_adj <= 0.05] <- "Down-regulated"

# find the top 10 positive and negative DEGs
tmp_p <- subset(NewPPB_EM, p_val_adj <= 0.05 & NEWLOG >= 0.5)
tmp_n <- subset(NewPPB_EM, p_val_adj <= 0.05 & NEWLOG <= -0.5)

top_pos_EM <- head(order(tmp_p[, "NEWLOG"], decreasing = T), 10)
top_pos_EM <- rownames(tmp_p[top_pos_EM, ])

top_neg_EM <- head(order(tmp_n[, "NEWLOG"], decreasing = F), 10)
top_neg_EM <- rownames(tmp_n[top_neg_EM, ])

# label the top genes
NewPPB_EM$label <- NA
NewPPB_EM[c(top_pos_EM, top_neg_EM), 'label'] <- c(top_pos_EM, top_neg_EM)

NewPPB_EM %>%
  ggplot(aes(
    x = NEWLOG,
    y = -log10(p_val_adj),
    col = diffexpressed,
    label = label
  )) +
  geom_point(size = 1) +
  xlab(expression("log"[2] * "FC")) +
  ylab(expression("-log"[10] * "Pvalue_adj")) +
  theme_minimal() +
  scale_color_manual(values = c("dodgerblue3", "gray50", "firebrick3")) +
  guides(colour = guide_legend(override.aes = list(size = 0.5))) +
  geom_vline(xintercept = c(-0.5, 0.5), col = "black") +
  geom_hline(yintercept = -log10(0.05), col = "black") +
  geom_label_repel(
    data = NewPPB_EM[c(top_pos_EM, top_neg_EM),],
    box.padding = 0.4,
    fill = 'white',
    max.overlaps = 20,
    show.legend = F,
    fontface = "italic"
  ) +
  labs(colour = NULL) +
  xlim(-3, 3)
ggsave("embryo_preproB_volcano.pdf")


#Adult VOLCANO PLOT

NewPPB_AD <- read.table("volcano_plots/Adult_PreproB_Cre-_vs_Cre+.txt", header = TRUE)
NewPPB_AD

NewPPB_AD$NEWLOG <- NewPPB_AD$avg_log2FC*(-1)

NewPPB_AD$diffexpressed <- "Unchanged"
NewPPB_AD$diffexpressed[NewPPB_AD$NEWLOG >= 0.5 & NewPPB_AD$p_val_adj <= 0.05] <- "Up-regulated"
NewPPB_AD$diffexpressed[NewPPB_AD$NEWLOG <= -0.5 & NewPPB_AD$p_val_adj <= 0.05] <- "Down-regulated"

tmp_p <- subset(NewPPB_AD, p_val_adj <= 0.05 & NEWLOG >= 0.5)
tmp_n <- subset(NewPPB_AD, p_val_adj <= 0.05 & NEWLOG <= -0.5)

top_pos_AD <- head(order(tmp_p[, "NEWLOG"], decreasing = T), 10)
top_pos_AD <- rownames(tmp_p[top_pos_AD, ])

top_neg_AD <- head(order(tmp_n[, "NEWLOG"], decreasing = F), 10)
top_neg_AD <- rownames(tmp_n[top_neg_AD, ])

NewPPB_AD$label <- NA
NewPPB_AD[c(top_pos_AD, top_neg_AD), 'label'] <- c(top_pos_AD, top_neg_AD)

NewPPB_AD %>%
  ggplot(aes(
    x = NEWLOG,
    y = -log10(p_val_adj),
    col = diffexpressed,
    label = label
  )) +
  geom_point(size = 1) +
  xlab(expression("log"[2] * "FC")) +
  ylab(expression("-log"[10] * "Pvalue_adj")) +
  theme_minimal() +
  scale_color_manual(values = c("dodgerblue3", "gray50", "firebrick3")) +
  guides(colour = guide_legend(override.aes = list(size = 0.5))) +
  geom_vline(xintercept = c(-0.5, 0.5), col = "black") +
  geom_hline(yintercept = -log10(0.05), col = "black") +
  geom_label_repel(
    data = NewPPB_AD[c(top_pos_AD, top_neg_AD),],
    box.padding = 0.4,
    fill = 'white',
    max.overlaps = 20,
    show.legend = F,
    fontface = "italic"
  ) +
  labs(colour = NULL) +
  xlim(-3, 3)
ggsave("adult_preproB_volcano.pdf")
