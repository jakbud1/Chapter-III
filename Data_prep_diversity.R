source("ASV_post_processing.R")

library(phyloseq)
library(microViz)
library(vegan)
library(tidyr)

df_sample_meta <- read_xlsx("Input/meta_samples.xlsx") %>% 
  column_to_rownames("ID") 

### V4 - Depth and diversity ----------------------------------------------
## Phyloseq's objects
ASV_V4 <- otu_table(df_V4_asv, taxa_are_rows = TRUE)
TAX_V4 <- tax_table(df_V4_tax)
samples <- sample_data(df_sample_meta)

V4_ps <- phyloseq(ASV_V4, TAX_V4, samples); rm(ASV_V4, TAX_V4, df_V4_tax)

sample_names(V4_ps); rank_names(V4_ps)

## Rarefraction curve
tiff("Output/V4_rarecurve.tiff", units = "px", width = 2000, height = 1100, pointsize = 16)
rarecurve(t(df_V4_asv), step = 200, cex = 1.75, cex.axis = 2.2, cex.lab = 1.5 ,col = "darkblue", ylab = "16S V4 rDNA ASVs", labelsize = 10)
dev.off()

## Export diversity indexes 
V4_out_rich <- estimate_richness(V4_ps, measure = c("Observed"))

### V9 - Depth and diversity ----------------------------------------------
## Phyloseq's objects 
ASV_V9 <- otu_table(df_V9_asv, taxa_are_rows = TRUE)
TAX_V9 <- tax_table(df_V9_tax)
samples <- sample_data(df_sample_meta)

V9_ps <- phyloseq(ASV_V9, TAX_V9, samples); rm(ASV_V9, TAX_V9, df_V9_tax)

sample_names(V9_ps); rank_names(V9_ps)

## Rarefraction curve
tiff("Output/V9_rarecurve.tiff", units = "px", width = 2000, height = 1100, pointsize = 16)
rarecurve(t(df_V9_asv),  step = 200, cex = 1.75, cex.axis = 2.2, cex.lab = 1.5, col = "darkgreen", ylab = "18S V9 rDNA ASVs")
dev.off()

## Export diversity indexes
V9_out_rich <- estimate_richness(V9_ps, measure = c("Observed"))

### Calculate relative abundance averaged for each Phylum (V4) and Class (V9) -----------------------
## V4
V4_rlt_in <- df_V4_s
V4_rlt_in$Phylum <- as.character(ifelse(is.na(df_V4_s$Phylum), "NA", df_V4_s$Phylum))

Phyla_V4 <- unique(V4_rlt_in$Phylum)

V4_rlt <- data.frame(matrix(NA, nrow = 157, ncol = length(Phyla_V4)+1))

for (i in 1:length(Phyla_V4)) {
  V4_rlt[,i] <- (colSums(subset(V4_rlt_in, Phylum == Phyla_V4[i])[,c(3:159)])/colSums(V4_rlt_in[,c(3:159)]))*100
  colnames(V4_rlt)[i] <- Phyla_V4[i]
}

V4_rlt[max(length(V4_rlt))] <- colnames(V4_rlt_in[c(3:159)])
colnames(V4_rlt)[max(length(V4_rlt))] <- "ID"

# Tests whether this loops work properly:
table(rowSums(V4_rlt[,c(1:39)]))

## V9
V9_rlt_in <- df_V9_s
V9_rlt_in$Class <- as.character(ifelse(is.na(df_V9_s$Class), "NA", df_V9_s$Class))

Class_V9 <- unique(V9_rlt_in$Class)

V9_rlt <- data.frame(matrix(NA, nrow = 155, ncol = length(Class_V9) + 1))

for (i in 1:length(Class_V9)) {
  V9_rlt[,i] <- (colSums(subset(V9_rlt_in, Class == Class_V9[i])[,c(3:157)])/colSums(V9_rlt_in[,c(3:157)]))*100
  colnames(V9_rlt)[i] <- Class_V9[i] 
}

V9_rlt[max(length(V9_rlt))] <- colnames(V9_rlt_in[c(3:157)])
colnames(V9_rlt)[max(length(V9_rlt))] <- "ID"

# Tests whether this loops work properly:
table((colSums(subset(df_V9_s, Class == "Metazoa")[,c(3:157)])/colSums(df_V9_s[,c(3:157)]))*100 == V9_rlt$Metazoa)
table(rowSums(V9_rlt[,c(1:27)]))

### Clear -----------------------------------------------------------------
rm(V4_ps, V9_ps, samples, 
   i, Class_V9, Phyla_V4, V9_rlt_in, V4_rlt_in, V4_t)

### Save data -------------------------------------------------------------
write.csv(df_track, "Output/Track_controls_drop.csv")

