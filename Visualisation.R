### Global ----------------------------------------------------------------
source("Models.R")

library(RColorBrewer)
library(ggplot2)
library(data.table)
library(corrplot)

### Add Tardigrada data of ASV data 
df_animals_gl <- aggregate(. ~ glacier, df_sample_meta[,c(1, 4, 5)], FUN = sum)
df_animals_gl <- merge(df_animals_gl, df_glac_dv[,c(1:5,6,10,11)], by = "glacier")

df_animals_gl$Tard_dens <- df_animals_gl$tard_count/df_animals_gl$volume_ml

V4_cor <- merge(V4_rlt_wide[,c(1,3:13, 15:18)], df_animals_gl[,c(1,11)], by = "glacier" )
V9_cor <- merge(V9_rlt_wide[,c(1,3:9, 11:14)], df_animals_gl[,c(1,11)], by = "glacier" )

### Color palletes
brewer.pal(n = 11, name = "PuOr")
brewer.pal(n = 11, name = "BrBG")

### Descriptive statistics ------------------------------------------------
## Radioactivity 
min(df_glac_dv$Pb210); max(df_glac_dv$Pb210)
mean(df_glac_dv$Pb210); sd(df_glac_dv$Pb210)
median(df_glac_dv$Pb210); IQR(df_glac_dv$Pb210)

min(df_glac_dv$Cs137); max(df_glac_dv$Cs137)
mean(df_glac_dv$Cs137); sd(df_glac_dv$Cs137)
median(df_glac_dv$Cs137); IQR(df_glac_dv$Cs137)

## Diversity
cor(df_animals_gl$V4_observed, df_animals_gl$V9_observed)

nrow(df_V4_s); ncol(df_V4_s)-11 # Number V4 ASVs and samples at the end 
nrow(df_V9_s); ncol(df_V9_s)-13 # Number V4 ASVs and samples at the end 

table(nchar(df_V4_s$ASV_seq)); median(nchar(df_V4_s$ASV_seq)) # V4 ASV seq length distribution
table(nchar(df_V9_s$ASV_seq)); median(nchar(df_V9_s$ASV_seq)) # V9 ASV seq length distribution

## Tardigrada
hist(df_glac_tar$tard_dens)
min(df_glac_tar$tard_dens); max(df_glac_tar$tard_dens)
median(df_glac_tar$tard_dens); IQR(df_glac_tar$tard_dens)

### Radioactivity ---------------------------------------------------------
df_total <- df_glac_tar[,c(1,4,5)]
df_total$total_rad <- df_total$Pb210 + df_total$Cs137

df_box <- melt(as.data.table(df_sample_meta[,c(1,3,8,9)]), id.vars = c("code", "glacier"), variable.name = "Isotope")
df_box$glacier <- factor(df_box$glacier, levels = df_total[order(df_total$total_rad), ]$glacier)

levels(df_box$Isotope) <- c("Pb-210", "Cs-137")

g1 <- ggplot(df_box, aes(glacier, value, color = Isotope, fill = Isotope))

g1.o <- g1 + geom_boxplot(alpha = 0.5) +
  theme_classic(base_size = 18) + ylab(expression(Activity~concentration~(Bq~kg^{-1}))) + 
  scale_color_manual(values = c("#8a6587", "#8394b4")) + scale_fill_manual(values = c("#8a6587", "#8394b4")) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "right") + xlab("Glacier");g1.o
ggsave("Output/Fig. 1.tiff", dpi = 300, scale = 1, width = 2900, height = 1500, units = "px", bg = "white")

### Diversity - Relative share variation between glaciers ------------------
# V4 - relative share 
V4_rlt_long_gg <- V4_rlt_long
V4_rlt_long_gg$glacier <- factor(V4_rlt_long_gg$glacier, levels = df_glac_dv[order(df_glac_dv$rad_total), ]$glacier)


gg_div_g_V4 <- ggplot(V4_rlt_long_gg, aes(x = glacier, y = Relative_share, fill = Phylum))
gg_div_g_V4 + geom_bar(stat = "identity", color = "black") + theme_classic(base_size = 17) + 
  ylab("Relative taxa abundance (%)") + xlab("Glacier") + scale_fill_brewer(palette = "BrBG") + 
  guides(fill=guide_legend(ncol = 1, byrow = TRUE)) +
  theme(legend.text = element_text(size = 14, colour = "black"), legend.position = "right", legend.key.size = unit(0.7,"line"), axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("Output/Fig. 2a.tiff", dpi = 300, scale = 1, width = 3000, height = 1600, units = "px", bg = "white")

# V9 - relative share
V9_rlt_long_gg <- V9_rlt_long
V9_rlt_long_gg$glacier <- factor(V9_rlt_long_gg$glacier, levels = df_glac_dv[order(df_glac_dv$rad_total), ]$glacier)

gg_div_g_V9 <- ggplot(V9_rlt_long_gg, aes(x = glacier, y = Relative_share, fill = Class))
gg_div_g_V9 + geom_bar(stat = "identity", color = "black") + theme_classic(base_size = 17) + 
  ylab("Relative taxa abundance (%)") + xlab("Glacier") + scale_fill_brewer(palette = "PuOr") + 
  guides(fill=guide_legend(ncol = 1, byrow = TRUE)) +
  theme(legend.text = element_text(size = 14, colour = "black"), legend.position = "right", legend.margin = margin(), legend.key.size = unit(0.7,"line"), axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("Output/Fig. 2b.tiff", dpi = 300, scale = 1, width = 3000, height = 1600, units = "px", bg = "white")

# V4 - ASV 
ggplot(df_cryoconite_dv, aes(x = glacier, y = V4_observed)) + 
  geom_boxplot()

### Diversity - correlations ----------------------------------------------
## V4 
names(V4_cor)[13:17] <- c("Pb-210", "Cs-137", 
                          "Glacier altitude", "Glacier size", "Tardigrada density")

V4_cor_m <- cor(V4_cor[,-1])

png("output/S1.a.png", width = 1500, height = 1600, pointsize = 40)
corrplot(V4_cor_m, method = 'circle', type = 'lower', insig='blank',
         number.cex = 0.8, diag = FALSE, tl.col = 'black', tl.srt = 45)
dev.off()

## V9 
names(V9_cor)[9:13] <- c("Pb-210", "Cs-137", 
                          "Glacier altitude", "Glacier size", "Tardigrada density")

V9_cor_m <- cor(V9_cor[,-1])

png("output/S1.b.png", width = 1500, height = 1600, pointsize = 40)
corrplot(V9_cor_m, method = 'circle', type = 'lower', insig='blank',
         number.cex = 0.8, diag = FALSE, tl.col = 'black', tl.srt = 45)
dev.off()

### Models results --------------------------------------------------------
## V4 
vis.V4.1 <- visreg(m_V4_R.1, "rad_total", scale = "response",  
                   gg = TRUE, plot = FALSE)
ggplot(vis.V4.1$fit, aes(y = visregFit, x = rad_total)) + geom_line(col = "#542788") + 
  geom_ribbon(aes(ymin = visregLwr, ymax = visregUpr), alpha = 0.3, col = "#542788") + 
  geom_point(data = data.frame(m_V4_R.2$frame), aes(x = rad_total, y = V4_observed), shape = 21, col = "#542788") +
  theme_classic() + xlab(expression(paste("", "Total radioactivity (Bq kg" ^-1,")"))) + 
  ylab("Average number of ASV on glacier")
ggsave("Output/Fig. 3a.png", scale = 0.5, width = 2000, height = 1800, units = "px", bg = "white")

vis.V4.2 <- visreg2d(m_V4_R.1, "rad_total", "rad_ratio", scale = "response", plot.type = "gg")
vis.V4.2 + theme_classic() + xlab(expression(paste("", "Total radioactivity (Bq kg" ^-1,")"))) + 
  ylab(expression(paste(""^"210", "Pb", "/"^"137", "Cs (activity ratio)"))) + 
  theme(legend.position = "bottom")
ggsave("Output/Fig. 3c.png", scale = 0.5, width = 2000, height = 1800, units = "px", bg = "white")

## V9 
vis.V9.1 <- visreg(m_V9_R.2, "rad_total", scale = "response",  
                   gg = TRUE, plot = FALSE)
ggplot(vis.V9.1$fit, aes(y = visregFit, x = rad_total)) + geom_line(col = "#01665E") + 
  geom_ribbon(aes(ymin = visregLwr, ymax = visregUpr), alpha = 0.3, col = "#01665E") + 
  geom_point(data = data.frame(df_glac_dv), 
             aes(x = rad_total, y = V9_observed), color = "#01665E", shape = 21) +
  theme_classic() + xlab(expression(paste("", "Total radioactivity (Bq kg" ^-1,")"))) + 
  ylab("Average number of ASV on glacier")
ggsave("Output/Fig. 3b.png", scale = 0.5, width = 2000, height = 1800, units = "px", bg = "white")

vis.V9.2 <- visreg2d(m_V9_R.1, "rad_total", "rad_ratio", scale = "response", plot.type = "gg")
vis.V9.2 + theme_classic() + xlab(expression(paste("", "Total radioactivity (Bq kg" ^-1,")"))) + 
  ylab(expression(paste(""^"210", "Pb", "/"^"137", "Cs (activity ratio)"))) + 
  theme(legend.position = "bottom")
ggsave("Output/Fig. 3d.png", scale = 0.5, width = 2000, height = 1800, units = "px", bg = "white")


## Tardigrada 
# seq_data
df_V9_tar <- subset(df_V9_s, Order == "Tardigrada")
df_V9_tar$sum.ASV[which.max(df_V9_tar$sum.ASV)]/sum(df_V9_tar$sum.ASV)
df_V9_tar$sum.ASV[20]/sum(df_V9_tar$sum.ASV)

vis.tar.1 <- visreg(mT.1, "rad_total", scale = "response",  
                   gg = TRUE, plot = FALSE)

ggplot(vis.tar.1$fit, aes(y = visregFit, x = rad_total)) + geom_line(col = "#d65917") + 
  geom_ribbon(aes(ymin = visregLwr, ymax = visregUpr), alpha = 0.3, col = "#d65917") + 
  geom_point(data = data.frame(mT.1$frame), aes(x = rad_total, y = tard_count), shape = 21, col = "#d65917") +
  theme_classic() + xlab(expression(paste("", "Total radioactivity (Bq kg" ^-1,")"))) + 
  ylab("Tardigrada abundance") + xlim(1500,16500)
ggsave("Output/Fig. 4.png", scale = 0.5, width = 2000, height = 1300, units = "px", bg = "white")