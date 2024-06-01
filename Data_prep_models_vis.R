### Global ---------------------------------------------------------------------
source("Data_prep_diversity.R")

df_glacier_meta <- read_xlsx("Input/meta_glaciers.xlsx")

df_sample_meta$ID <- rownames(df_sample_meta)

### df - diversity analysis --------------------------------------------
## Cryoconite level 
V4_out_rich$ID <- sub('X', '', rownames(V4_out_rich))
V9_out_rich$ID <- sub('X', '', rownames(V9_out_rich))

df_all <- df_sample_meta

df_all <- merge(df_all, V4_out_rich[,c(1:2)], by = "ID"); rm(V4_out_rich)

names(df_all)[13] <- c("V4_observed")

df_cryoconite_dv <- merge(df_all, V9_out_rich[,c(1:2)], by = "ID"); rm(V9_out_rich, df_all)
names(df_cryoconite_dv)[14] <- c("V9_observed")

## Glacier level 
df_glac_dv <- aggregate(df_cryoconite_dv[,c(2,9,10,13,14)], .~ glacier, FUN = "mean")
df_glac_dv <- left_join(df_glac_dv, df_glacier_meta[,c(2:4,7,9)], by = "glacier")

### df - Tardigrada ----------------------------------------------------
df_glac_tar <- df_sample_meta

## Glacier level 
df_glac_tar <- aggregate(df_sample_meta[,c(1,2,4:6)], .~ glacier, FUN = "mean")

### df - Diversity across glaciers ----------------------------------------
## V4
# Transfer to long format 
V4_rlt_long <- gather(merge(V4_rlt, df_sample_meta[c(1:70, 81:170),c(1,12)], by = "ID"), 
                      Phylum, Relative_share, -c(ID,glacier), factor_key = TRUE)

# Calculate mean by glacier
V4_rlt_long <- aggregate(V4_rlt_long, Relative_share ~ glacier + Phylum, FUN = "mean")

# Drop Phyla with mean relative share < 0.1
keep_V4 <- aggregate(V4_rlt_long, Relative_share ~ Phylum, FUN = "mean")[round(aggregate(V4_rlt_long, Relative_share ~ Phylum, FUN = "mean")$Relative_share) > 0.1,]

V4_rlt_long <- V4_rlt_long[V4_rlt_long$Phylum %in% keep_V4$Phylum,]

# Add "Others" group
Others_V4 <- aggregate(V4_rlt_long, Relative_share ~ glacier, FUN = "sum") 
Others_V4$Relative_share <- 100 - Others_V4$Relative_share
Others_V4$Phylum <- "Others"

V4_rlt_long <- rbind(V4_rlt_long, Others_V4); rm(Others_V4)
V4_rlt_wide <- spread(V4_rlt_long, Phylum, Relative_share); rm(V4_rlt)
V4_rlt_wide <- merge(V4_rlt_wide, df_glac_dv[,c(1:3,8,9)], by = "glacier")

## V9
# Transfer to long format 
V9_rlt_long <- gather(merge(V9_rlt, df_sample_meta[c(1:70, 81:170),c(1,12)], by = "ID"), 
                      Class, Relative_share, -c(ID,glacier), factor_key = TRUE)

# Calculate mean by glacier
V9_rlt_long <- aggregate(V9_rlt_long, Relative_share ~ glacier + Class, FUN = "mean")

# Drop Class with mean relative share < 0.1
keep_V9 <- aggregate(V9_rlt_long, Relative_share ~ Class, FUN = "mean")[round(aggregate(V9_rlt_long, Relative_share ~ Class, FUN = "mean")$Relative_share) > 0.1,]

V9_rlt_long <- V9_rlt_long[V9_rlt_long$Class %in% keep_V9$Class,]

# Add "Others" group
Others_V9 <- aggregate(V9_rlt_long, Relative_share ~ glacier, FUN = "sum") 
Others_V9$Relative_share <- 100 - Others_V9$Relative_share
Others_V9$Class <- "Others"

V9_rlt_long <- rbind(V9_rlt_long, Others_V9)
V9_rlt_wide <- spread(V9_rlt_long, Class, Relative_share); rm(V9_rlt)
V9_rlt_wide <- merge(V9_rlt_wide, df_glac_dv[,c(1:3,8,9)], by = "glacier")

### Clear -----------------------------------------------------------------
rm(Others_V9, keep_V4, keep_V9, df_V9_asv, df_V4_asv)

