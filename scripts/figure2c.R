setwd("/n/holylabs/LABS/dwang_lab/Lab/T2D_meta/3.analysis/R1/fig2")
library(tidyverse)
library(ggplot2)

load("/n/holylabs/LABS/dwang_lab/Lab/T2D_meta/3.analysis/dropbox/fig3/species_resid_callout_all.RData")

metadata$study[metadata$study=="Fromentin_2022 (DEU_DNK_FRA)"] <- "Fromentin_2022 (DEU/DNK/FRA)"
metadata$study[metadata$study=="SOL (USA)"] <- "HCHS/SOL (USA)"
studycol <- c("#466983", "#339900", "#837B8D", 
              "#809900", "#5DB1DD", "#802268",
              "#CC9900", "#CE3D32","#C75127", "#996600")
names(studycol) <- c("DIRECT-PLUS (ISR)", "Fromentin_2022 (DEU/DNK/FRA)", "HPFS (USA)", 
                     "Karlsson_2013 (SWE)", "NHSII (USA)", "Qin_2012 (CHN)",
                     "HCHS/SOL (USA)", "Wu_2020a (SWE)", "Wu_2020b (SWE)", "Zhong_2019 (CHN)")

metadata$status_new <- ifelse(metadata$status_new=="Con", "Control", metadata$status_new)
metadata$status_new <- ifelse(metadata$status_new=="Pre", "Prediabetes", metadata$status_new)

# "s__Roseburia_intestinalis", "s__Clostridium_bolteae", "s__Escherichia_coli"
# "s__Bacteroides_fragilis", "s__Flavonifractor_plautii", "s__Alistipes_putredinis"

# "s__Roseburia_intestinalis"
spe00 <- data.frame(id=colnames(spe_resid[,2:8118]),
                    spe=as.numeric(t(spe_resid[spe_resid$spe=="s__Roseburia_intestinalis", -1])) + spe_means[spe_means$spe=="s__Roseburia_intestinalis",]$mean)
spe00 <- merge(metadata[, c("id", "study", "status_new")], spe00, by="id") %>% filter(!is.na(spe))
spe00plot <- spe00 %>% group_by(study, status_new) %>% summarise(mean=mean(spe), sd=sd(spe))

summary(spe00plot$mean)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -4.625  -3.749  -3.545  -3.482  -3.323  -2.177

ymin <- min(spe00plot$mean)-(max(spe00plot$mean)-min(spe00plot$mean))/5
ymax <- max(spe00plot$mean)+(max(spe00plot$mean)-min(spe00plot$mean))/5
p0 <- ggplot(aes(x=status_new, y=mean), data=spe00plot) +
  labs(y=expression(log[10]*" (abundance)")) +
  geom_boxplot(width=0.6, outlier.shape = NA) +
  geom_jitter(aes(fill=study), size=3, width = 0.07, shape=21) +
  scale_fill_manual(values=studycol) +
  scale_y_continuous(limits=c(ymin,ymax),
                     labels = scales::number_format(accuracy = 0.1)) +
  labs(title="Roseburia intestinalis") +
  ylab("") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1.2),
        axis.ticks.x = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "italic", size=11, colour="black"),
        axis.title.x = element_blank(),
        axis.text = element_text(color="black", size=10),
        axis.title.y = element_text(color="black", size=12),
        legend.position = "none") +
  annotate(geom="text", x=1.2, y=ymin+(ymax-ymin)/18, size=3.7, col="black", label="FDR=0.17")

# "s__Escherichia_coli"
spe22 <- data.frame(id=colnames(spe_resid[,2:8118]),
                    spe=as.numeric(t(spe_resid[spe_resid$spe=="s__Escherichia_coli", -1])) + spe_means[spe_means$spe=="s__Escherichia_coli",]$mean)
spe22 <- merge(metadata[, c("id", "study", "status_new")], spe22, by="id") %>% filter(!is.na(spe))
spe22plot <- spe22 %>% group_by(study, status_new) %>% summarise(mean=mean(spe), sd=sd(spe))

summary(spe22plot$mean)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# -5.242  -4.735  -4.522  -4.491  -4.292  -3.768

ymin <- min(spe22plot$mean)-(max(spe22plot$mean)-min(spe22plot$mean))/5
ymax <- max(spe22plot$mean)+(max(spe22plot$mean)-min(spe22plot$mean))/5
p2 <- ggplot(aes(x=status_new, y=mean), data=spe22plot) +
  labs(y=expression(log[10]*" (abundance)")) +
  geom_boxplot(width=0.6, outlier.shape = NA) +
  geom_jitter(aes(fill=study), size=3, width = 0.07, shape=21) +
  scale_fill_manual(values=studycol) +
  scale_y_continuous(limits=c(ymin,ymax)) +
  labs(title="Escherichia coli", fill="Study") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1.2),
        axis.ticks.x = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "italic", size=11, colour="black"),
        axis.title.x = element_blank(),
        axis.text = element_text(color="black", size=10),
        axis.title.y = element_text(color="black", size=12),
        legend.position="none") +
  annotate(geom="text", x=1.3, y=ymax-(ymax-ymin)/18, size=3.7, col="black", label="FDR<0.001")

# "s__Flavonifractor_plautii"
spe55 <- data.frame(id=colnames(spe_resid[,2:8118]),
                    spe=as.numeric(t(spe_resid[spe_resid$spe=="s__Flavonifractor_plautii", -1])) + spe_means[spe_means$spe=="s__Flavonifractor_plautii",]$mean)
spe55 <- merge(metadata[, c("id", "study", "status_new")], spe55, by="id") %>% filter(!is.na(spe))
spe55plot <- spe55 %>% group_by(study, status_new) %>% summarise(mean=mean(spe), sd=sd(spe))

summary(spe55plot$mean)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# -4.938  -4.325  -4.184  -4.176  -4.083  -3.096

ymin <- min(spe55plot$mean)-(max(spe55plot$mean)-min(spe55plot$mean))/5
ymax <- max(spe55plot$mean)+(max(spe55plot$mean)-min(spe55plot$mean))/5
p5 <- ggplot(aes(x=status_new, y=mean), data=spe55plot) +
  labs(y=expression(log[10]*" (abundance)")) +
  geom_boxplot(width=0.6, outlier.shape = NA) +
  geom_jitter(aes(fill=study), size=3, width = 0.07, shape=21) +
  scale_fill_manual(values=studycol) +
  scale_y_continuous(limits=c(ymin,ymax)) +
  labs(title="Flavonifractor plautii") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1.2),
        axis.ticks.x = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "italic", size=12, colour="black"),
        axis.title.x = element_blank(),
        axis.text = element_text(color="black", size=10),
        axis.title.y = element_text(color="black", size=12),
        legend.position = "none") +
  annotate(geom="text", x=1.3, y=ymax-(ymax-ymin)/18, size=3.7, col="black", label="FDR=0.004")

# Parasutterella_excrementihominis
spe66 <- data.frame(id=colnames(spe_resid[,2:8118]),
                    spe=as.numeric(t(spe_resid[spe_resid$spe=="s__Parasutterella_excrementihominis", -1])) + spe_means[spe_means$spe=="s__Parasutterella_excrementihominis",]$mean)
spe66 <- merge(metadata[, c("id", "study", "status_new")], spe66, by="id") %>% filter(!is.na(spe))
spe66plot <- spe66 %>% group_by(study, status_new) %>% summarise(mean=mean(spe), sd=sd(spe))

summary(spe66plot$mean)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# -5.037  -4.626  -4.518  -4.499  -4.394  -3.794

ymin <- min(spe66plot$mean)-(max(spe66plot$mean)-min(spe66plot$mean))/5
ymax <- max(spe66plot$mean)+(max(spe66plot$mean)-min(spe66plot$mean))/5
p6 <- ggplot(aes(x=status_new, y=mean), data=spe66plot) +
  labs(y=expression(log[10]*" (abundance)")) +
  geom_boxplot(width=0.6, outlier.shape = NA) +
  geom_jitter(aes(fill=study), size=3, width = 0.07, shape=21) +
  scale_fill_manual(values=studycol) +
  scale_y_continuous(limits=c(ymin,ymax)) +
  ylab("") +
  labs(title="Parasutterella nexcrementihominis") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1.2),
        axis.ticks.x = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "italic", size=11, colour="black"),
        axis.title.x = element_blank(),
        axis.text = element_text(color="black", size=10),
        axis.title.y = element_text(color="black", size=12),
        legend.position = "none") +
  annotate(geom="text", x=1.2, y=ymin+(ymax-ymin)/18, size=3.7, col="black", label="FDR=0.17")

pdf("/n/holylabs/LABS/dwang_lab/Lab/T2D_meta/3.analysis/R1/fig2/figure2c.pdf",
    width = 5, height = 4.8, onefile = F) # Open a new pdf file
ggpubr::ggarrange(p2,p0,p5,p6, nrow=2, ncol=2)
dev.off() # Close the file

legend <- ggplot(aes(x=status_new, y=mean), data=spe22plot) +
  labs(y=expression(log[10]*" (abundance)")) +
  geom_boxplot(width=0.6, outlier.shape = NA) +
  geom_jitter(aes(fill=study), size=3, width = 0.07, shape=21) +
  scale_fill_manual(values=studycol) +
  scale_y_continuous(limits=c(ymin,ymax)) +
  labs(title="Escherichia coli", fill="Study") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1.2),
        axis.ticks.x = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "italic", size=11, colour="black"),
        axis.title.x = element_blank(),
        axis.text = element_text(color="black", size=10),
        axis.title.y = element_text(color="black", size=12),
        legend.position="right",
        legend.title = element_text(color="black", size=10),
        legend.key.height= unit(.5, 'cm'),
        legend.key.width= unit(.5, 'cm'),
        legend.spacing.x = unit(0.03, 'cm'),
        legend.spacing.y = unit(0.01, 'cm')) +
  guides(fill = guide_legend(nrow = 10, byrow = T,
                             override.aes = list(size = 3),
                             title.position = "top"))

pdf("/n/holylabs/LABS/dwang_lab/Lab/T2D_meta/3.analysis/R1/fig2/figure2c_legend.pdf",
    width = 3, height = 3, onefile = F) # Open a new pdf file
grid::grid.draw(ggpubr::get_legend(legend))
dev.off() # Close the file

