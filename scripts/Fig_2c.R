setwd("/n/holylabs/LABS/dwang_lab/Lab/T2D_meta/3.analysis/dropbox/fig3/")

load("/n/holylabs/LABS/dwang_lab/Lab/T2D_meta/3.analysis/0323/species_0.RData")
gdata::keep(metadata, species_zero_adj, sure=T)

identical(colnames(species_zero_adj), rownames(metadata)) # TRUE
table(metadata$study)
# DIRECT-PLUS (ISR) Fromentin_2022 (DEU_DNK_FRA)           HPFS (USA)          Karlsson_2013 (SWE) 
# 280                         1005                          925                          145 
# NHSII (USA)               Qin_2012 (CHN)                 SOL (USA)               Wu_2020a (SWE) 
# 814                          344                         2881                          991 
# Wu_2020b (SWE)             Zhong_2019 (CHN) 
# 484                          248 

metadata$study <- as.character(metadata$study)
metadata[is.na(metadata$randid), ]$randid <- 1:6378

species_zero_adj <- species_zero_adj[, metadata$id]
identical(colnames(species_zero_adj), rownames(metadata)) # TRUE

metadata$metfYes <- ifelse(metadata$metf=="Yes",1,0)
metadata$metfNA <- ifelse(metadata$metf=="Missing",1,0)

species_zero_adj_raw <- species_zero_adj
species_zero_adj[is.na(species_zero_adj)] <- 0

library(MMUPHin)
meta <- lm_meta(feature_abd = as.matrix(species_zero_adj),
                batch = "study",
                exposure = "age",
                covariates = c("sex", "bmi", "metfNA", "metfYes"),
                covariates_random = 'randid',
                data = metadata,
                control = list(verbose = T,
                               rma_method = "FE",
                               transform = "LOG",
                               output = "species_resid_all"))

direct <- as.data.frame(readRDS("/n/holylabs/LABS/dwang_lab/Lab/T2D_meta/3.analysis/dropbox/fig3/species_resid_all/DIRECT-PLUS (ISR)/residuals.rds"))
# check number of non-NA features
# species_direct <- species_zero_adj_raw[, metadata[metadata$study=="DIRECT-PLUS (ISR)", ]$id]
# species_direct <- species_direct[!is.na(species_direct[,1]), ]
colnames(direct) <- metadata[metadata$study=="DIRECT-PLUS (ISR)", ]$id
direct$spe <- rownames(direct)

mbs <- as.data.frame(readRDS("/n/holylabs/LABS/dwang_lab/Lab/T2D_meta/3.analysis/dropbox/fig3/species_resid_all/NHSII (USA)/residuals.rds"))
colnames(mbs) <- metadata[metadata$study=="NHSII (USA)", ]$id
mbs$spe <- rownames(mbs)

fromentin <- as.data.frame(readRDS("/n/holylabs/LABS/dwang_lab/Lab/T2D_meta/3.analysis/dropbox/fig3/species_resid_all/Fromentin_2022 (DEU_DNK_FRA)/residuals.rds"))
colnames(fromentin) <- metadata[metadata$study=="Fromentin_2022 (DEU_DNK_FRA)", ]$id
fromentin$spe <- rownames(fromentin)

karlsson <- as.data.frame(readRDS("/n/holylabs/LABS/dwang_lab/Lab/T2D_meta/3.analysis/dropbox/fig3/species_resid_all/Karlsson_2013 (SWE)/residuals.rds"))
colnames(karlsson) <- metadata[metadata$study=="Karlsson_2013 (SWE)", ]$id
karlsson$spe <- rownames(karlsson)

mlvs <- as.data.frame(readRDS("/n/holylabs/LABS/dwang_lab/Lab/T2D_meta/3.analysis/dropbox/fig3/species_resid_all/HPFS (USA)/residuals.rds"))
colnames(mlvs) <- metadata[metadata$study=="HPFS (USA)", ]$id
mlvs$spe <- rownames(mlvs)

shenzhen <- as.data.frame(readRDS("/n/holylabs/LABS/dwang_lab/Lab/T2D_meta/3.analysis/dropbox/fig3/species_resid_all/Qin_2012 (CHN)/residuals.rds"))
colnames(shenzhen) <- metadata[metadata$study=="Qin_2012 (CHN)", ]$id
shenzhen$spe <- rownames(shenzhen)

sol <- as.data.frame(readRDS("/n/holylabs/LABS/dwang_lab/Lab/T2D_meta/3.analysis/dropbox/fig3/species_resid_all/SOL (USA)/residuals.rds"))
colnames(sol) <- metadata[metadata$study=="SOL (USA)", ]$id
sol$spe <- rownames(sol)

suzhou <- as.data.frame(readRDS("/n/holylabs/LABS/dwang_lab/Lab/T2D_meta/3.analysis/dropbox/fig3/species_resid_all/Zhong_2019 (CHN)/residuals.rds"))
colnames(suzhou) <- metadata[metadata$study=="Zhong_2019 (CHN)", ]$id
suzhou$spe <- rownames(suzhou)

wua <- as.data.frame(readRDS("/n/holylabs/LABS/dwang_lab/Lab/T2D_meta/3.analysis/dropbox/fig3/species_resid_all/Wu_2020a (SWE)/residuals.rds"))
colnames(wua) <- metadata[metadata$study=="Wu_2020a (SWE)", ]$id
wua$spe <- rownames(wua)

wub <- as.data.frame(readRDS("/n/holylabs/LABS/dwang_lab/Lab/T2D_meta/3.analysis/dropbox/fig3/species_resid_all/Wu_2020b (SWE)/residuals.rds"))
colnames(wub) <- metadata[metadata$study=="Wu_2020b (SWE)", ]$id
wub$spe <- rownames(wub)

spe_resid <- Reduce(function(x, y) merge(x=x, y=y, by=c("spe"), all=T),
                    list(direct, mbs, karlsson, fromentin, mlvs,
                         sol, shenzhen, suzhou, wua, wub))

spe_resid$spe <- as.numeric(gsub("T", "", spe_resid$spe))
spe_resid <- spe_resid %>% arrange(spe)
spe_resid$spe <- rownames(species_zero_adj)

species_zero_adj_raw$spe <- rownames(species_zero_adj)
species_zero_adj_raw <- species_zero_adj_raw[, c(8118, 1:8117)]

species_zero_adj_raw_v2 <- species_zero_adj_raw
for (i in 1:dim(species_zero_adj_raw_v2)[1]) {
  aa <- as.numeric(t(species_zero_adj_raw_v2[i, -1]))
  species_zero_adj_raw_v2[i, ][species_zero_adj_raw_v2[i, ]==0] <- min(aa[aa>0], na.rm=T)/2
}

species_zero_adj_raw_v2[,-1] <- log10(species_zero_adj_raw_v2[,-1])
spe_means <- data.frame(spe=species_zero_adj_raw_v2$spe,
                        mean=rowMeans(species_zero_adj_raw_v2[,-1], na.rm=T))

save.image(file="/n/holylabs/LABS/dwang_lab/Lab/T2D_meta/3.analysis/dropbox/fig3/species_resid_callout_all.RData")

metadata$study[metadata$study=="Fromentin_2022 (DEU_DNK_FRA)"] <- "Fromentin_2022 (DEU/DNK/FRA)"
studycol <- c("#466983", "#339900", "#837B8D", 
              "#809900", "#5DB1DD", "#802268",
              "#CC9900", "#CE3D32","#C75127", "#996600")
names(studycol) <- c("DIRECT-PLUS (ISR)", "Fromentin_2022 (DEU/DNK/FRA)", "HPFS (USA)", 
                     "Karlsson_2013 (SWE)", "NHSII (USA)", "Qin_2012 (CHN)",
                     "SOL (USA)", "Wu_2020a (SWE)", "Wu_2020b (SWE)", "Zhong_2019 (CHN)")

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
  labs(y=expression("Adjusted "*log[10]*" (abundance)")) +
  geom_boxplot(width=0.6, outlier.shape = NA) +
  geom_jitter(aes(fill=study), size=3, width = 0.07, shape=21) +
  scale_fill_manual(values=studycol) +
  scale_y_continuous(limits=c(ymin,ymax)) +
  labs(title="Roseburia intestinalis") +
  ylab("") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1.2),
        axis.ticks = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "italic", size=11, colour="black"),
        axis.title.x = element_blank(),
        axis.text.x = element_text(color="black", size=11),
        axis.title.y = element_text(color="black", size=11),
        axis.text.y = element_blank(),
        legend.position = "none") +
  annotate(geom="text", x=1.2, y=ymin+(ymax-ymin)/18, size=3.7, col="black", label="FDR=0.17")

# "s__Clostridium_bolteae"
spe11 <- data.frame(id=colnames(spe_resid[,2:8118]),
                    spe=as.numeric(t(spe_resid[spe_resid$spe=="s__Clostridium_bolteae", -1])) + spe_means[spe_means$spe=="s__Clostridium_bolteae",]$mean)
spe11 <- merge(metadata[, c("id", "study", "status_new")], spe11, by="id") %>% filter(!is.na(spe))
spe11plot <- spe11 %>% group_by(study, status_new) %>% summarise(mean=mean(spe), sd=sd(spe))

summary(spe11plot$mean)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -5.784  -5.352  -5.186  -5.164  -5.059  -3.849

ymin <- min(spe11plot$mean)-(max(spe11plot$mean)-min(spe11plot$mean))/5
ymax <- max(spe11plot$mean)+(max(spe11plot$mean)-min(spe11plot$mean))/5
p1 <- ggplot(aes(x=status_new, y=mean), data=spe11plot) +
  geom_boxplot(width=0.6, outlier.shape = NA) +
  geom_jitter(aes(fill=study), size=3, width = 0.07, shape=21) +
  scale_fill_manual(values=studycol) +
  scale_y_continuous(limits=c(ymin,ymax)) +
  labs(title="Clostridium bolteae") +
  labs(y=expression("Adjusted "*log[10]*" (abundance)")) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1.2),
        axis.ticks = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "italic", size=11, colour="black"),
        axis.title.x = element_blank(),
        axis.text.x = element_text(color="black", size=11),
        axis.title.y = element_text(color="black", size=11),
        axis.text.y = element_blank(),
        legend.position="bottom",
        legend.title = element_blank(),
        legend.key.height= unit(.5, 'cm'),
        legend.key.width= unit(.5, 'cm'),
        legend.spacing.x = unit(0.03, 'cm'),
        legend.spacing.y = unit(0.01, 'cm')) +
  guides(fill = guide_legend(nrow = 3, byrow = T,
                             override.aes = list(size = 3))) +
  annotate(geom="text", x=1.2, y=ymax-(ymax-ymin)/18, size=3.7, col="black", label="FDR=0.07")

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
  labs(y=expression("Adjusted "*log[10]*" (abundance)")) +
  geom_boxplot(width=0.6, outlier.shape = NA) +
  geom_jitter(aes(fill=study), size=3, width = 0.07, shape=21) +
  scale_fill_manual(values=studycol) +
  scale_y_continuous(limits=c(ymin,ymax)) +
  labs(title="Escherichia coli") +
  ylab("") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1.2),
        axis.ticks = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "italic", size=11, colour="black"),
        axis.title.x = element_blank(),
        axis.text.x = element_text(color="black", size=11),
        axis.title.y = element_text(color="black", size=11),
        axis.text.y = element_blank(),
        legend.position = "none") +
  annotate(geom="text", x=1.2, y=ymax-(ymax-ymin)/18, size=3.7, col="black", label="FDR<0.001")

# "s__Bacteroides_fragilis"
spe44 <- data.frame(id=colnames(spe_resid[,2:8118]),
                    spe=as.numeric(t(spe_resid[spe_resid$spe=="s__Bacteroides_fragilis", -1])) + spe_means[spe_means$spe=="s__Bacteroides_fragilis",]$mean)
spe44 <- merge(metadata[, c("id", "study", "status_new")], spe44, by="id") %>% filter(!is.na(spe))
spe44plot <- spe44 %>% group_by(study, status_new) %>% summarise(mean=mean(spe), sd=sd(spe))

summary(spe44plot$mean)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -5.622  -4.849  -4.576  -4.623  -4.325  -3.085

ymin <- min(spe44plot$mean)-(max(spe44plot$mean)-min(spe44plot$mean))/5
ymax <- max(spe44plot$mean)+(max(spe44plot$mean)-min(spe44plot$mean))/5
p4 <- ggplot(aes(x=status_new, y=mean), data=spe44plot) +
  geom_boxplot(width=0.6, outlier.shape = NA) +
  geom_jitter(aes(fill=study), size=3, width = 0.07, shape=21) +
  scale_fill_manual(values=studycol) +
  scale_y_continuous(limits=c(ymin,ymax)) +
  theme_bw() +
  labs(y=expression("Adjusted "*log[10]*" (abundance)")) +
  labs(title="Bacteroides fragilis") +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1.2),
        axis.ticks = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "italic", size=11, colour="black"),
        axis.title.x = element_blank(),
        axis.text.x = element_text(color="black", size=11),
        axis.title.y = element_text(color="black", size=11),
        axis.text.y = element_blank(),
        legend.position = "none") +
  annotate(geom="text", x=1.2, y=ymax-(ymax-ymin)/18, size=3.7, col="black", label="FDR=0.07")

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
  labs(y=expression("Adjusted "*log[10]*" (abundance)")) +
  geom_boxplot(width=0.6, outlier.shape = NA) +
  geom_jitter(aes(fill=study), size=3, width = 0.07, shape=21) +
  scale_fill_manual(values=studycol) +
  scale_y_continuous(limits=c(ymin,ymax)) +
  ylab("") +
  labs(title="Flavonifractor plautii") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1.2),
        axis.ticks = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "italic", size=12, colour="black"),
        axis.title.x = element_blank(),
        axis.text.x = element_text(color="black", size=12),
        axis.title.y = element_text(color="black", size=12),
        axis.text.y = element_blank(),
        legend.position = "none") +
  annotate(geom="text", x=1.2, y=ymax-(ymax-ymin)/18, size=3.7, col="black", label="FDR=0.004")

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
  labs(y=expression("Adjusted "*log[10]*" (abundance)")) +
  geom_boxplot(width=0.6, outlier.shape = NA) +
  geom_jitter(aes(fill=study), size=3, width = 0.07, shape=21) +
  scale_fill_manual(values=studycol) +
  scale_y_continuous(limits=c(ymin,ymax)) +
  ylab("") +
  labs(title="Parasutterella nexcrementihominis") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1.2),
        axis.ticks = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "italic", size=11, colour="black"),
        axis.title.x = element_blank(),
        axis.text.x = element_text(color="black", size=11),
        axis.title.y = element_text(color="black", size=11),
        axis.text.y = element_blank(),
        legend.position = "none") +
  annotate(geom="text", x=1.2, y=ymin+(ymax-ymin)/18, size=3.7, col="black", label="FDR=0.17")

pdf("/n/holylabs/LABS/dwang_lab/Lab/T2D_meta/3.analysis/dropbox/fig3/species_resid_callout_all.pdf",
    width = 6.75, height = 5.7, onefile = F) # Open a new pdf file
ggpubr::ggarrange(p1,p2,p0,p4,p5,p6, nrow=2, ncol=3, 
                  legend = "bottom", common.legend = T)
dev.off() # Close the file

