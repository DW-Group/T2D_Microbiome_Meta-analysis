setwd("/n/holylabs/LABS/dwang_lab/Lab/T2D_meta/3.analysis/dropbox/fig2/")

library(tidyverse)
load("/n/holylabs/LABS/dwang_lab/Lab/T2D_meta/3.analysis/0323/species_0.RData")
gdata::keep(metadata, spe_all, spe_relab,
            species, species_zero, species_final,
            species_zero_adj0, species_zero_adj, sure=T)

studycol <- c("#466983", "#339900", "#837B8D", 
              "#809900", "#5DB1DD", "#802268",
              "#CC9900", "#CE3D32","#C75127", "#996600")
names(studycol) <- c("DIRECT-PLUS (ISR)", "Fromentin_2022 (DEU/DNK/FRA)", "HPFS (USA)", 
                     "Karlsson_2013 (SWE)", "NHSII (USA)", "Qin_2012 (CHN)",
                     "SOL (USA)", "Wu_2020a (SWE)", "Wu_2020b (SWE)", "Zhong_2019 (CHN)")

metadata$study <- as.character(metadata$study)
metadata$study <- ifelse(metadata$study=="Fromentin_2022 (DEU_DNK_FRA)",
                         "Fromentin_2022 (DEU/DNK/FRA)",
                         metadata$study)

####################################################################
#               figure 1b (distribution for species)
####################################################################
table(spe_all$present)
# 1  2  3  4  5  6  7  8  9 10 
# 34 30 11 12 11 18 20  7 25 71

spe_all <- merge(spe_all, spe_relab, by="species")
universal <- spe_all %>% filter(present==10) %>% arrange(desc(relab))
universal <- universal[1:25, ]
universal$group <- "Universal"
overlap <- spe_all %>% filter(present!=10 & present!=1) %>% arrange(desc(relab))
overlap <- overlap[1:25, ]
overlap$group <- "Overlapping"
solo <- spe_all %>% filter(present==1) %>% arrange(desc(relab))
solo <- solo[1:25, ]
solo$group <- "Singular"

spe_plot <- rbind(universal,overlap,solo)
spe_plot2 <- spe_plot %>% select(-"relab") %>%
  gather(study, relab, str_subset(colnames(spe_plot), "_relab")[1:10]) %>%
  filter(relab>0)
spe_plot2 <- spe_plot2[, c(1,13:15)] # species, group, relab
spe_plot2$group <- factor(spe_plot2$group, levels=c("Universal", "Overlapping", "Singular"))
spe_plot2$species <- gsub("s__", "", spe_plot2$species)
spe_plot2$species <- gsub("_", " ", spe_plot2$species)

spe_order <- spe_plot2 %>% group_by(group, species) %>%
  summarise(mean=mean(relab)) %>%
  ungroup(species) %>% arrange(group, desc(mean))
spe_plot2$species <- factor(spe_plot2$species, levels=spe_order$species)

spe_plot2 <- spe_plot2 %>% mutate(study=case_when(study=="karlsson_relab" ~ "Karlsson_2013 (SWE)",
                                                  study=="mbs_relab" ~ "NHSII (USA)",
                                                  study=="mlvs_relab" ~ "HPFS (USA)",
                                                  study=="qin_relab" ~ "Qin_2012 (CHN)",
                                                  study=="zhong_relab" ~ "Zhong_2019 (CHN)",
                                                  study=="direct_relab" ~ "DIRECT-PLUS (ISR)",
                                                  study=="wu_relab1" ~ "Wu_2020a (SWE)",
                                                  study=="wu_relab2" ~ "Wu_2020b (SWE)",
                                                  study=="sol_relab" ~ "SOL (USA)",
                                                  study=="pedersen_relab" ~ "Fromentin_2022 (DEU/DNK/FRA)"))
spe_plot2$study <- factor(spe_plot2$study, 
                          levels=c("DIRECT-PLUS (ISR)", "Fromentin_2022 (DEU/DNK/FRA)", "HPFS (USA)", 
                                   "Karlsson_2013 (SWE)", "NHSII (USA)", "Qin_2012 (CHN)", 
                                   "SOL (USA)", "Wu_2020a (SWE)", "Wu_2020b (SWE)", "Zhong_2019 (CHN)"))

plain <- function(x,...) {
  format(x, ..., scientific = FALSE, drop0trailing = TRUE)
}

library(ggsci)
library(ggplot2)
spe_dist <- ggplot(spe_plot2, aes(x=species, y=relab)) +
  geom_boxplot() +
  geom_point(aes(fill=study), shape=21, size=3) +
  scale_fill_manual(values=studycol) +
  ylab("Relative abundance (%)") +
  xlab("") +
  scale_y_log10(breaks=c(10, 1, 0.1, 0.01), labels=plain) +
  facet_wrap( ~ group, scales = "free_x") +
  theme(strip.text.x = element_text(hjust = 0.5, color="black", size=14)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,
                                   size = 11, face = "italic", color="black"),
        axis.text.y = element_text(size = 14, color="black"),
        axis.title.y = element_text(size = 14, color="black")) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black")) +
  theme(legend.position = c(.85, .8),
        legend.title = element_blank(),
        legend.text = element_text(size = 10, color="black"),
        legend.background=element_blank(),
        legend.key.height= unit(.3, 'cm'),
        legend.key.width= unit(.3, 'cm'),
        legend.key=element_rect(fill="white")) +
  guides(fill = guide_legend(nrow=5, byrow=TRUE))

ggsave(spe_dist, filename="/n/holylabs/LABS/dwang_lab/Lab/T2D_meta/3.analysis/dropbox/fig2/figure2a.pdf",
       width =14, height = 6.4)

####################################################################
#                  figure s2 (MMUPHin compare)
####################################################################
# manual TSS normalization
# https://rdrr.io/cran/hilldiv/src/R/tss.R
species_zero_adj0 <- sweep(species_zero_adj0, 2, colSums(species_zero_adj0), FUN="/")
summary(colSums(species_zero_adj0)) # should be all 1

species_zero <- sweep(species_zero, 2, colSums(species_zero), FUN="/")
summary(colSums(species_zero)) # should be all 1

# PCoA data before batch correction
library(vegan)
identical(colnames(species_zero), metadata$id) # TRUE
distance.matrix <- vegdist(t(species_zero), method="bray")

## do the MDS math (this is basically eigen value decomposition)
mds.stuff <- cmdscale(distance.matrix, eig=TRUE, x.ret=TRUE, k=20)

## calculate the percentage of variation that each MDS axis accounts for...
mds.var.per <- round(mds.stuff$eig/sum(mds.stuff$eig)*100, 1)

## now make a fancy looking plot that shows the MDS axes and the variation:
mds.values <- mds.stuff$points
mds.data <- data.frame(id=rownames(mds.values),
                       X=mds.values[,1],
                       Y=mds.values[,2])

identical(mds.data$id, metadata$id) # TRUE
mds.data <- cbind(mds.data, metadata[,-1])

#############################################
#     PCoA data after batch correction
#############################################
identical(colnames(species_zero_adj0), metadata$id) # TRUE
distance.matrix2 <- vegdist(t(species_zero_adj0), method="bray")

## do the MDS math (this is basically eigen value decomposition)
mds.stuff2 <- cmdscale(distance.matrix2, eig=TRUE, x.ret=TRUE, k=20)

## calculate the percentage of variation that each MDS axis accounts for...
mds.var.per2 <- round(mds.stuff2$eig/sum(mds.stuff2$eig)*100, 1)

## now make a fancy looking plot that shows the MDS axes and the variation:
mds.values2 <- mds.stuff2$points
mds.data2 <- data.frame(id=rownames(mds.values2),
                        X=mds.values2[,1],
                        Y=mds.values2[,2])

identical(mds.data2$id, metadata$id) # TRUE
mds.data2 <- cbind(mds.data2, metadata[,-1])

set.seed(1234)
R2_before <- adonis(distance.matrix ~ study, data=metadata) # R2=0.08456
set.seed(1234)
R2_after <- adonis(distance.matrix2 ~ study, data=metadata) # R2=0.04026

# draw PCoA plot
mds.data$study <- factor(mds.data$study, levels=c("DIRECT-PLUS (ISR)", "Fromentin_2022 (DEU/DNK/FRA)", "HPFS (USA)", 
                                                  "Karlsson_2013 (SWE)", "NHSII (USA)", "Qin_2012 (CHN)",
                                                  "SOL (USA)", "Wu_2020a (SWE)", "Wu_2020b (SWE)", "Zhong_2019 (CHN)"))
mds.data2$study <- factor(mds.data2$study, levels=c("DIRECT-PLUS (ISR)", "Fromentin_2022 (DEU/DNK/FRA)", "HPFS (USA)", 
                                                    "Karlsson_2013 (SWE)", "NHSII (USA)", "Qin_2012 (CHN)",
                                                    "SOL (USA)", "Wu_2020a (SWE)", "Wu_2020b (SWE)", "Zhong_2019 (CHN)"))
p_study1 <- ggplot(data=mds.data, aes(x=-X, y=-Y)) +
  # geom_point(size=3.5, aes(color=study), shape=19, alpha=0.3) +
  geom_point(size=3.5, aes(fill=study), shape=21) +
  scale_fill_manual(values=studycol) +
  xlab(paste("PCo1 - ", mds.var.per[1], "%", sep="")) +
  ylab(paste("PCo2 - ", mds.var.per[2], "%", sep="")) +
  labs(title="") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 14, colour = "black"),
        axis.text.y = element_text(size = 14, colour = "black"),
        axis.title.x = element_text(size = 14, colour = "black"),
        axis.title.y = element_text(size = 14, colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        plot.title = element_text(hjust = 0.5, size = 14, colour = "black"),
        legend.position="bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 11, colour = "black"),
        legend.spacing.x = unit(0.1, 'cm')) +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE, override.aes = list(size = 4))) +
  annotate(geom="text", x=0.42, y=0.36, col="black", size=4.5,
           label="Before correction, R2=8.5%")

p_study2 <- ggplot(data=mds.data2, aes(x=X, y=Y)) +
  # geom_point(size=3.5, aes(color=study), shape=19, alpha=0.3) +
  geom_point(size=3.5, aes(fill=study), shape=21) +
  scale_fill_manual(values=studycol) +
  xlab(paste("PCo1 - ", mds.var.per2[1], "%", sep="")) +
  ylab(paste("PCo2 - ", mds.var.per2[2], "%", sep="")) +
  labs(title="") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 14, colour = "black"),
        axis.text.y = element_text(size = 14, colour = "black"),
        axis.title.x = element_text(size = 14, colour = "black"),
        axis.title.y = element_text(size = 14, colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        plot.title = element_text(hjust = 0.5, size = 14, colour = "black"),
        legend.position = "none") +
  annotate(geom="text", x=0.42, y=0.32, col="black", size=4.5,
           label="After correction, R2=4.0%")

pdf("/n/holylabs/LABS/dwang_lab/Lab/T2D_meta/3.analysis/dropbox/fig2/figure_s1.pdf",
    width = 10, height = 5.8, onefile = F) # Open a new pdf file
ggpubr::ggarrange(p_study1, p_study2, 
                  ncol=2, legend = "bottom", common.legend = T)
dev.off() # Close the file

####################################################################
#                       figure 1c (PCoA)
####################################################################
statuscol <- c("#5DB1DD", "#CC9900", "#C75127")
names(statuscol) <- c("Con", "Pre", "T2D")
p_study3 <- ggplot(data=mds.data2, aes(x=X, y=Y)) +
  geom_point(size=4, aes(fill=status_new), shape=21) +
  scale_color_manual(values = statuscol) +
  xlab(paste("PCo1 - ", mds.var.per2[1], "%", sep="")) +
  ylab(paste("PCo2 - ", mds.var.per2[2], "%", sep="")) +
  labs(title="PCoA plot by diabetes status") +
  scale_x_continuous(breaks = c(-0.3,0,0.3,0.6), limits=c(-0.32,0.66)) +
  scale_y_continuous(breaks = c(-0.4,-0.2,0,0.2), limits=c(-0.51,0.38)) +
  scale_fill_manual(values = statuscol) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 16, colour = "black"),
        axis.text.y = element_text(size = 16, colour = "black"),
        axis.title.x = element_text(size = 16, colour = "black"),
        axis.title.y = element_text(size = 16, colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        plot.title = element_text(size = 16, hjust = 0.5)) +
  theme(legend.position = c(.88, .88),
        legend.title = element_blank(),
        legend.key.height= unit(.7, 'cm'),
        legend.key.width= unit(.3, 'cm'),
        legend.text = element_text(size = 16, colour = "black")) +
  annotate(geom="text", x=-0.08, y=0.37, col="black", size=5.2,
           label="PERMANOVA P<0.001")

p_bacter <- ggplot(data=mds.data2, aes(x=X, y=Y)) +
  geom_point(size=3.6, aes(fill=Bacteroide), shape=21) +
  scale_fill_gradient(low = "#56B1F7", high = "#132B43") +
  xlab(paste("PCo1 - ", mds.var.per2[1], "%", sep="")) +
  ylab(paste("PCo2 - ", mds.var.per2[2], "%", sep="")) +
  labs(title="PCoA plot by Bacteroidetes") +
  scale_x_continuous(breaks = c(-0.3,0,0.3,0.6), limits=c(-0.32,0.66)) +
  scale_y_continuous(breaks = c(-0.4,-0.2,0,0.2), limits=c(-0.51,0.36)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 14, colour = "black"),
        axis.text.y = element_text(size = 14, colour = "black"),
        axis.title.x = element_text(size = 14, colour = "black"),
        axis.title.y = element_text(size = 14, colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.position="none")

p_firm <- ggplot(data=mds.data2, aes(x=X, y=Y)) +
  geom_point(size=3.6, aes(fill=Firmicu), shape=21) +
  scale_fill_gradient(low = "#56B1F7", high = "#132B43") +
  xlab(paste("PCo1 - ", mds.var.per2[1], "%", sep="")) +
  ylab(paste("PCo2 - ", mds.var.per2[2], "%", sep="")) +
  labs(title="PCoA plot by Firmicutes",
       fill="Relative \n abundance \n (%)") +
  scale_x_continuous(breaks = c(-0.3,0,0.3,0.6), limits=c(-0.32,0.66)) +
  scale_y_continuous(breaks = c(-0.4,-0.2,0,0.2), limits=c(-0.51,0.36)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 14, colour = "black"),
        axis.text.y = element_text(size = 14, colour = "black"),
        axis.title.x = element_text(size = 14, colour = "black"),
        axis.title.y = element_text(size = 14, colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title.align=0.5,
        legend.direction="vertical")

library(patchwork)
pdf("/n/holylabs/LABS/dwang_lab/Lab/T2D_meta/3.analysis/dropbox/fig2/figure2b.pdf",
    width = 10, height = 6, onefile = F) # Open a new pdf file
(p_study3 | ((p_bacter / p_firm) + plot_layout(guides = 'collect'))) + 
  plot_layout(ncol = 2, widths = c(2.2, 1))
dev.off() # Close the file

save.image(file="/n/holylabs/LABS/dwang_lab/Lab/T2D_meta/3.analysis/dropbox/fig2/figure2.RData")

