setwd("/n/holylabs/LABS/dwang_lab/Lab/T2D_meta/3.analysis/dropbox/fig3/")

library(tidyverse)
load("/n/holylabs/LABS/dwang_lab/Lab/T2D_meta/3.analysis/0323/species_trend.RData")
meta_1trend <- meta_1
meta_3trend <- meta_3

gdata::keep(meta_1trend, meta_3trend, sure=T)

load("/n/holylabs/LABS/dwang_lab/Lab/T2D_meta/3.analysis/0323/species_T2D_main.RData")
gdata::keep(meta_1trend, meta_3trend, 
            meta_1, meta_3, taxa, sure=T)

spe_all <- Reduce(union, list(
  meta_1trend[meta_1trend$qval.fdr<0.25,]$feature,
  meta_3trend[meta_3trend$qval.fdr<0.25,]$feature,
  meta_1[meta_1$qval.fdr<0.25,]$feature,
  meta_3[meta_3$qval.fdr<0.25,]$feature
))

# exclude "s__Desulfovibrio_piger"
spe_all <- spe_all[-32]

data1trend <- meta_1trend[meta_1trend$feature %in% spe_all, ]
data3trend <- meta_3trend[meta_1trend$feature %in% spe_all, ]
data1 <- meta_1[meta_1$feature %in% spe_all, ]
data3 <- meta_3[meta_3$feature %in% spe_all, ]

data1trend$model <- "Ordinal: base model"
data3trend$model <- "Ordinal: base + metformin"
data1$model <- "Binary: base model"
data3$model <- "Binary: base + metformin"

heatmap_trend <- rbind(
  data1trend[, c("feature", "model", "coef", "pval", "qval.fdr")],
  data3trend[, c("feature", "model", "coef", "pval", "qval.fdr")]
)

heatmap <- rbind(
  data1[, c("feature", "model", "coef", "pval", "qval.fdr")],
  data3[, c("feature", "model", "coef", "pval", "qval.fdr")]
)

heatmap_trend$star <- ifelse(heatmap_trend$qval.fdr<0.25, "*", "")
heatmap_trend$species <- gsub("s__", "", heatmap_trend$feature)
heatmap_trend$species <- gsub("_", " ", heatmap_trend$species)

heatmap$star <- ifelse(heatmap$qval.fdr<0.25, "*", "")
heatmap$species <- gsub("s__", "", heatmap$feature)
heatmap$species <- gsub("_", " ", heatmap$species)

colnames(taxa)[8] <- "feature"
heatmap <- merge(taxa[, c("feature", "X2")], heatmap, by="feature", all.y=T)
heatmap$X2 <- gsub("p__", "", heatmap$X2)

summary(heatmap$coef)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# -1.1302 -0.3275  0.1307  0.0741  0.4151  1.4483

summary(heatmap_trend$coef)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# -0.47354 -0.16133  0.07361  0.04414  0.21554  0.71066

sum(heatmap$coef>1.2) # 2
heatmap$coef[heatmap$coef>=1.2] <- 1.2
sum(heatmap_trend$coef>0.6) # 2
heatmap_trend$coef[heatmap_trend$coef>=0.6] <- 0.6

sum(data1trend$coef>0) # 53
data1trend <- merge(taxa[, c("feature", "X2")], data1trend, by="feature", all.y=T)
data1trend$X2 <- gsub("p__", "", data1trend$X2)
data1trend$species <- gsub("s__","",data1trend$feature)
data1trend$species <- gsub("_"," ",data1trend$species)

data1trend <- data1trend %>% arrange(desc(X2), desc(species))

heatmap$sign <- ifelse((heatmap$feature %in% data1trend[data1trend$coef>0,]$feature),1,0)
heatmap$model <- factor(heatmap$model, levels=c("Binary: base model", 
                                                "Binary: base + metformin"))
heatmap_trend$sign <- ifelse((heatmap_trend$feature %in% data1trend[data1trend$coef>0,]$feature),1,0)
heatmap_trend$model <- factor(heatmap_trend$model, levels=c("Ordinal: base model", 
                                                            "Ordinal: base + metformin"))

heatmap1 <- heatmap %>% filter(sign==1)
heatmap1$species <- factor(heatmap1$species, levels=data1trend[data1trend$coef>0,]$species)
heatmap2 <- heatmap %>% filter(sign==0)
heatmap2$species <- factor(heatmap2$species, levels=data1trend[data1trend$coef<0,]$species)

heatmap_trend1 <- heatmap_trend %>% filter(sign==1)
heatmap_trend1$species <- factor(heatmap_trend1$species, levels=data1trend[data1trend$coef>0,]$species)
heatmap_trend2 <- heatmap_trend %>% filter(sign==0)
heatmap_trend2$species <- factor(heatmap_trend2$species, levels=data1trend[data1trend$coef<0,]$species)
  
table(heatmap[heatmap$qval.fdr<0.25,]$model)
# Binary: base model       Binary: base + metformin 
# 75                       30 

table(heatmap_trend[heatmap_trend$qval.fdr<0.25,]$model)
# Ordinal: base model       Ordinal: base + metformin 
# 89                        36

beta <- data.frame(mod1=meta_1$coef,
                   mod3=meta_3$coef,
                   mod1trend=meta_1trend$coef,
                   mod3trend=meta_3trend$coef)
cor(beta)
#                mod1      mod3 mod1trend mod3trend
# mod1      1.0000000 0.8664462 0.9393763 0.8134984
# mod3      0.8664462 1.0000000 0.8028757 0.8687571
# mod1trend 0.9393763 0.8028757 1.0000000 0.9071312
# mod3trend 0.8134984 0.8687571 0.9071312 1.0000000

phycol <- c("#1BAD4C", "#A7AC36", "#4A6EB6", "#282274", "#F03B20", "#1AAAAA")
names(phycol) <- c("Actinobacteria", "Bacteroidetes", "Euryarchaeota", "Firmicutes",
                   "Fusobacteria", "Proteobacteria")

library(ggplot2)
heat1_0 <- ggplot(heatmap1, aes(x=1, y=species)) +
  geom_tile(aes(fill=X2)) +
  labs(fill="Phylum") +
  scale_fill_manual(values=phycol) +
  xlab("") +
  theme_bw() +
  theme(panel.border = element_rect(size=1.2), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.ticks.y = element_line(colour = "black", size = 0.8),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(face = "italic", color = "black", size = 14),
        axis.text.x = element_blank(),
        axis.title =element_blank(),
        legend.position = "none",
        legend.text = element_text(colour = "black", size = 14),
        legend.title = element_text(colour = "black", size = 14)
  )

heat1_trend <- ggplot(heatmap_trend1, aes(x=model, y=species)) +
  geom_tile(aes(fill=coef)) +
  scale_fill_gradient2(low="#08519c",
                       mid="white",
                       high="#bd0026",
                       midpoint = 0,
                       guide = "colourbar",
                       breaks=c(-0.5,0,0.7),
                       limits=c(-0.5,0.7)) +
  scale_x_discrete(position = "top") +
  geom_text(aes(label=star), color="black", size=7, nudge_y = -0.2) +
  labs(color=expression(beta * " coefficient")) +
  xlab("") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_line(colour = "black", size = 0.8),
        axis.title.x = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_line(colour = "black", size = 0.8),
        axis.text.x = element_text(angle = 45, vjust = 0, hjust = 0, size = 14, colour = "black"),
        # axis.ticks.length=unit(0.07,"inch"),
        panel.border = element_rect(color = "black", fill = NA, size = 1.5)) +
  theme(legend.position = "none")

heat1 <- ggplot(heatmap1, aes(x=model, y=species)) +
  geom_tile(aes(fill=coef)) +
  scale_fill_gradient2(low="#5e3c99",
                       mid="white",
                       high="#e66101",
                       midpoint = 0,
                       guide = "colourbar",
                       breaks=c(-1.2,0,1.2),
                       limits=c(-1.2,1.2)) +
  scale_x_discrete(position = "top") +
  geom_text(aes(label=star), color="black", size=7, nudge_y = -0.2) +
  labs(color=expression(beta * " coefficient")) +
  xlab("") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_line(colour = "black", size = 0.8),
        axis.text.x = element_text(angle = 45, vjust = 0, hjust = 0, size = 14, colour = "black"),
        # axis.ticks.length=unit(0.07,"inch"),
        panel.border = element_rect(color = "black", fill = NA, size = 1.5)) +
  theme(legend.position = "none")

heat2_0 <- ggplot(heatmap2, aes(x=1, y=species)) +
  geom_tile(aes(fill=X2)) +
  scale_fill_manual(values=phycol) +
  scale_y_discrete(position = "right") +
  xlab("") +
  theme_bw() +
  theme(panel.border = element_rect(size=1.2), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.ticks.y = element_line(colour = "black", size = 0.8),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(face = "italic", color = "black", size = 14),
        axis.text.x = element_blank(),
        axis.title =element_blank(),
        legend.position = "none")

heat2_trend <- ggplot(heatmap_trend2, aes(x=model, y=species)) +
  geom_tile(aes(fill=coef)) +
  scale_fill_gradient2(low="#08519c",
                       mid="white",
                       high="#bd0026",
                       midpoint = 0,
                       guide = "colourbar",
                       breaks=c(-0.5,0,0.7),
                       limits=c(-0.5,0.7)) +
  scale_x_discrete(position = "top") +
  scale_y_discrete(position = "right") +
  geom_text(aes(label=star), color="black", size=7, nudge_y = -0.2) +
  labs(fill=expression(beta * " coefficient")) +
  xlab("") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_line(colour = "black", size = 0.8),
        axis.text.x = element_text(angle = 45, vjust = 0, hjust = 0, size = 14, colour = "black"),
        axis.ticks.length=unit(0.04,"inch"),
        panel.border = element_rect(color = "black", fill = NA, size = 1.5),
        legend.position = "none",
        legend.text = element_text(colour = "black", size = 14),
        legend.title = element_text(colour = "black", size = 14),
        # legend.justification = c(2, 0), 
        legend.box.margin=margin(c(0,-100,0,0))) +
  guides(fill = guide_colourbar(title.position="top",
                                title.hjust = 0.5,
                                barheight = 1.5,
                                barwidth = 10))

heat2 <- ggplot(heatmap2, aes(x=model, y=species)) +
  geom_tile(aes(fill=coef)) +
  scale_fill_gradient2(low="#5e3c99",
                       mid="white",
                       high="#e66101",
                       midpoint = 0,
                       guide = "colourbar",
                       breaks=c(-1.2,0,1.2),
                       limits=c(-1.2,1.2)) +
  scale_x_discrete(position = "top") +
  scale_y_discrete(position = "right") +
  geom_text(aes(label=star), color="black", size=7, nudge_y = -0.2) +
  labs(fill=expression(beta * " coefficient")) +
  xlab("") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_line(colour = "black", size = 0.8),
        axis.title.x = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_line(colour = "black", size = 0.8),
        axis.text.x = element_text(angle = 45, vjust = 0, hjust = 0, size = 14, colour = "black"),
        axis.ticks.length=unit(0.04,"inch"),
        panel.border = element_rect(color = "black", fill = NA, size = 1.5),
        legend.position = "none",
        legend.text = element_text(colour = "black", size = 14),
        legend.title = element_text(colour = "black", size = 14),
        # legend.justification = c(2, 0), 
        legend.box.margin=margin(c(0,-100,0,0))) +
  guides(fill = guide_colourbar(title.position="top",
                                title.hjust = 0.5,
                                barheight = 1.5,
                                barwidth = 10))

pdf("/n/holylabs/LABS/dwang_lab/Lab/T2D_meta/3.analysis/dropbox/fig3/figure3av3.pdf",
    width = 8.3, height = 35, onefile = F) # Open a new pdf file
egg::ggarrange(heat1_trend, heat1, heat1_0, 
               heat2_trend, heat2, heat2_0, nrow=2, ncol=3,
               widths = c(1,1,0.2), heights=c(54,45))
dev.off() # Close the file

# legend
trendlegend <- ggplot(heatmap_trend2, aes(x=model, y=species)) +
  geom_tile(aes(fill=coef)) +
  scale_fill_gradient2(low="#08519c",
                       mid="white",
                       high="#bd0026",
                       midpoint = 0,
                       guide = "colourbar",
                       breaks=c(-0.5,0,0.7),
                       limits=c(-0.5,0.7)) +
  scale_x_discrete(position = "bottom") +
  scale_y_discrete(position = "right") +
  geom_text(aes(label=star), color="black", size=7, nudge_y = -0.2) +
  labs(fill=expression(beta * " coefficient (ordinal)")) +
  xlab("") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_line(colour = "black", size = 0.8),
        axis.title.x = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_line(colour = "black", size = 0.8),
        axis.text.x = element_text(angle = 45, vjust = 0, hjust = 0, size = 14, colour = "black"),
        axis.ticks.length=unit(0.04,"inch"),
        panel.border = element_rect(color = "black", fill = NA, size = 1.5),
        legend.position = "top",
        legend.text = element_text(colour = "black", size = 14),
        legend.title = element_text(colour = "black", size = 14)) +
  guides(fill = guide_colourbar(title.position="top",
                                title.hjust = 0.5,
                                barheight = 1.5,
                                barwidth = 10))

legend <- ggplot(heatmap2, aes(x=model, y=species)) +
  geom_tile(aes(fill=coef)) +
  scale_fill_gradient2(low="#5e3c99",
                       mid="white",
                       high="#e66101",
                       midpoint = 0,
                       guide = "colourbar",
                       breaks=c(-1.2,0,1.2),
                       limits=c(-1.2,1.2)) +
  scale_x_discrete(position = "top") +
  scale_y_discrete(position = "right") +
  geom_text(aes(label=star), color="black", size=7, nudge_y = -0.2) +
  labs(fill=expression(beta * " coefficient (binary)")) +
  xlab("") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_line(colour = "black", size = 0.8),
        axis.title.x = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_line(colour = "black", size = 0.8),
        axis.text.x = element_text(angle = 45, vjust = 0, hjust = 0, size = 14, colour = "black"),
        axis.ticks.length=unit(0.04,"inch"),
        panel.border = element_rect(color = "black", fill = NA, size = 1.5),
        legend.position = "bottom",
        legend.text = element_text(colour = "black", size = 14),
        legend.title = element_text(colour = "black", size = 14)) +
  guides(fill = guide_colourbar(title.position="top",
                                title.hjust = 0.5,
                                barheight = 1.5,
                                barwidth = 10))

pdf("/n/holylabs/LABS/dwang_lab/Lab/T2D_meta/3.analysis/dropbox/fig3/figure3a_legend1.pdf",
    width = 3, height = 2, onefile = F) # Open a new pdf file
grid::grid.draw(ggpubr::get_legend(trendlegend))
dev.off() # Close the file

pdf("/n/holylabs/LABS/dwang_lab/Lab/T2D_meta/3.analysis/dropbox/fig3/figure3a_legend2.pdf",
    width = 3, height = 2, onefile = F) # Open a new pdf file
grid::grid.draw(ggpubr::get_legend(legend))
dev.off() # Close the file

save.image(file="/n/holylabs/LABS/dwang_lab/Lab/T2D_meta/3.analysis/dropbox/fig3/figure3a.RData")

# mean/se for trend
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -3.4213 -2.4149  2.3084  0.6684  3.0459  5.1359 

# mean/se for T2DCon
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -3.6626 -2.4319  2.3659  0.4945  2.9224  4.3925
