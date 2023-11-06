rm(list = ls())
library(tidyverse)
library(vegan)
library(RColorBrewer)
library(diverse)
library(grid)
library(ggpubr)
library(viridis)
library(gridExtra)
library(cowplot)
library(haven)
library(data.table)
library(sjmisc)

(date_log <- Sys.Date())
setwd("/n/holylabs/LABS/dwang_lab/Lab/T2D_meta/3.analysis/dropbox/fig_s3")

t2d_col <- c(
  Con = "#5DB1DD",
  Pre = "#CC9900",
  T2D = "#C75127"
)

# rank by sample size
study_col <- c(
  `SOL (USA)`= "#CC9900",
  # `Fromentin_2022 (DEU/DNK/FRA)` = "#339900",
  `Fromentin_2022 (DEU_DNK_FRA)` = "#339900",
  `HPFS (USA)` = "#837B8D", 
  `NHSII (USA)` = "#5DB1DD",
  `Wu_2020a (SWE)`= "#CE3D32",
  `Wu_2020b (SWE)`= "#C75127",
  `Qin_2012 (CHN)`= "#802268",
  `DIRECT-PLUS (ISR)` = "#466983",
  `Zhong_2019 (CHN)`= "#996600",
  `Karlsson_2013 (SWE)` = "#809900"
)
# study_col <- c(
#   `DIRECT-PLUS (ISR)` = "#466983",
#   `Fromentin_2022 (DEU/DNK/FRA)` = "#339900",
#   `HPFS (USA)` = "#837B8D", 
#   `Karlsson_2013 (SWE)` = "#809900", 
#   `NHSII (USA)` = "#5DB1DD",
#   `Qin_2012 (CHN)`= "#802268",
#   `SOL (USA)`= "#CC9900",
#   `Wu_2020a (SWE)`= "#CE3D32",
#   `Wu_2020b (SWE)`= "#C75127",
#   `Zhong_2019 (CHN)`= "#996600"
# )


load("data_for_fig_s3.RData")

head(path_all);dim(path_all) #292/22
table(path_all$present)
# 1  2  3  4  5  6  7  8  9 10 
# 42 33 23 21 19 17 22 22 28 65 
path_all <- filter(path_all,present==10)
dim(path_all) #99

head(enzy_all);dim(enzy_all) #1473/20
table(enzy_all$present) 
#  1   2   3   4   5   6   7   8   9  10 
# 223 114  93  72  92  97 114 130 183 356 
enzy_all <- filter(enzy_all,present==10)
dim(enzy_all) #365/20

# Log Transformation
# LOG <- function(x) {
#   y <- replace(x, x == 0, min(x[x>0],na.rm = T) / 2)
#   return(log10(y))
# }

# Species Distributions --------------
dim(species) #238/8123; 239/8117
head(colnames(species))
spe_ave_relab <- 
  data.frame(species=rownames(species),
             ave_relab=rowSums(species,na.rm = T)/ncol(species))
spe_ave_relab <- arrange(spe_ave_relab,desc(ave_relab))
head(spe_ave_relab)
(top10_spe <- head(spe_ave_relab$species,10))
head(colnames(species))
head(rownames(species))
head(colSums(species,na.rm = T))
species <- species + 0.1
# species <- species/100
species <- species %>%
  mutate(across(everything(),log10)) #log10
class(species)
species_top10 <- species[top10_spe,] %>% 
  rownames_to_column(var = "species") 
(spe_order <- factor(species_top10$species,levels = species_top10$species))
head(colnames(species_top10))
tail(colnames(species_top10)) 
species_top10_long <- 
  gather(species_top10,key = "ID",value = "reabun",S112:x10MCx1541)
head(species_top10_long)

## order samples by reabundance of P.copri
head(spe_metadata)
table(spe_metadata$status_new)
spe_metadata$status_new <- 
  factor(spe_metadata$status_new,levels = c("Con","Pre","T2D"))

table(spe_metadata$study)
spe_metadata$study <- 
  factor(spe_metadata$study, 
         levels=c("SOL (USA)","Fromentin_2022 (DEU_DNK_FRA)","HPFS (USA)",
                  "NHSII (USA)","Wu_2020a (SWE)","Wu_2020b (SWE)",
                  "Qin_2012 (CHN)","DIRECT-PLUS (ISR)",
                  "Zhong_2019 (CHN)", "Karlsson_2013 (SWE)"))

p.copri <- species[c("s__Prevotella_copri","s__Bacteroides_vulgatus"),] %>% 
  rotate_df() %>% 
  rownames_to_column(var = "id") %>% 
  inner_join(spe_metadata,by=c("id")) %>% 
  arrange(desc(s__Prevotella_copri)) %>% 
  arrange(status_new) %>% 
  arrange(study)
head(p.copri);dim(p.copri) 
p.copri$id <- factor(p.copri$id,levels = p.copri$id)

sample_order <- factor(p.copri$id,levels = p.copri$id)

species_top10_long$ID <- factor(species_top10_long$ID,levels = sample_order)
species_top10_long$species <- 
  factor(species_top10_long$species,levels = spe_order)

head(species_top10_long)

## species heatmap ------
p_spe <- ggplot(species_top10_long, aes(ID,species)) +
  geom_tile(aes(fill = reabun)) +
  scale_fill_viridis_c(name ="Relative abundance\nlog10 scale", na.value = "white",
                       direction = -1, option = "D",
                       breaks=c(-1,0,1,2),limits=c(-1,2))+
  scale_y_discrete(limits = rev(levels(spe_order)),
                   labels = c("Roseburia faecis","Ruminococcus bromii","Bacteroides dorei",
                              "Bacteroides stercoris","Alistipes putredinis","Eubacterium rectale",
                              "Faecalibacterium prausnitzii","Bacteroides uniformis",
                              "Bacteroides vulgatus","Prevotella copri")) +
  theme(legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.position = "none",
        plot.margin = unit(c(1,0,1,5.5),"point"),
        plot.title = element_blank(),
        axis.title=element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y=element_text(size=10, face = "italic"),
        axis.text.x = element_blank())
p_spe

scales::show_col(t2d_col)
p_t2d <- ggplot(p.copri,aes(x=id,y=1,fill=factor(status_new))) +
  geom_col(width = 1)+
  scale_fill_manual(values = t2d_col)+
  labs(y="T2D status")+
  theme_void()+
  theme(axis.title.y = element_text(size=11),
        plot.margin = unit(c(1,0,1,5.5),"point"),
        legend.position = "none")
p_t2d


p_study <- ggplot(p.copri,aes(x=id,y=1,fill=study)) +
  geom_col(width = 1)+
  scale_fill_manual(values = study_col)+
  labs(y="Study")+
  theme_void()+
  theme(axis.title.y = element_text(size=11),
        plot.margin = unit(c(5.5,0,1,5.5),"point"),
        legend.position = "none")
p_study

blank <- ggplot()+theme_void()
pdf("fig_s3_species.pdf",width = 12, height = 5, onefile = F)
egg::ggarrange(p_study,p_t2d,p_spe,ncol = 1,heights = c(0.1,0.1,2))
dev.off()

## legend -----
## heatmap
p_spe_lgd <- ggplot(species_top10_long, aes(ID,species)) +
  geom_tile(aes(fill = reabun)) +
  scale_fill_viridis_c(name ="Relative abundance\nlog10 scale", 
                       na.value = "white",direction = -1, option = "D",
                       breaks=c(-1,0,1,2),#labels=c("Minimum",0.5,"Maximum"),
                       limits=c(-1,2))+
  scale_y_discrete(limits = rev(levels(spe_order)),
                   labels = c("Roseburia faecis","Ruminococcus bromii","Bacteroides dorei",
                              "Bacteroides stercoris","Alistipes putredinis","Eubacterium rectale",
                              "Faecalibacterium prausnitzii","Bacteroides uniformis",
                              "Bacteroides vulgatus","Prevotella copri")) +
  theme(legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.position = "right",
        plot.title = element_blank(),
        axis.title=element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y=element_text(size=10, face = "italic"),
        axis.text.x = element_blank())+
  guides(fill = guide_colourbar(title.position="top", 
                                title.hjust = 0.5,
                                barheight = 6,
                                barwidth = 1))
p_spe_lgd
spe_lgd <- ggpubr::get_legend(p_spe_lgd)

pdf("species_plot_legend.pdf", height = 3,width = 3)
grid.draw(spe_lgd)
dev.off()

t2d_col
## t2d
p_t2d_lgd <- ggplot(p.copri,aes(x=id,y=1,fill=status_new)) +
  geom_col(width = 1)+
  scale_fill_manual(breaks = c("Con","Pre","T2D"),
                    values = c("#5DB1DD", "#CC9900", "#C75127"),
                    labels = c("Control","Prediabetes","T2D"),
                    name=NULL)+
  labs(y="T2D status")+
  theme_void()+
  theme(axis.title.y = element_text(size=11))
p_t2d_lgd
t2d_lgd <- ggpubr::get_legend(p_t2d_lgd)

pdf("t2d_status_legend.pdf", height = 3,width = 3)
grid.draw(t2d_lgd)
dev.off()

## study
p_study_lgd <- ggplot(p.copri,aes(x=id,y=1,fill=study)) +
  geom_col(width = 1)+
  scale_fill_manual(values = study_col,name=NULL,
                    labels = c("SOL (USA)",
                               "Fromentin_2022 (DEU/DNK/FRA)",
                               "HPFS (USA)", 
                               "NHSII (USA)",
                               "Wu_2020a (SWE)",
                               "Wu_2020b (SWE)",
                               "Qin_2012 (CHN)",
                               "DIRECT-PLUS (ISR)",
                               "Zhong_2019 (CHN)",
                               "Karlsson_2013 (SWE)"))+
  labs(y="Study")+
  theme_void()+
  theme(axis.title.y = element_text(size=11))
p_study_lgd
study_lgd <- ggpubr::get_legend(p_study_lgd)

pdf("study_legend.pdf", height = 3,width = 3)
grid.draw(study_lgd)
dev.off()

#  Pathway Distribution -------------
dim(pathway) #293/8115
head(colnames(pathway),30)
head(rownames(pathway),30)
sum(pathway$S153,na.rm = T)

pwy_ave_relab <- data.frame(pathway=rownames(pathway),
                            ave_relab=rowSums(pathway,na.rm = T)/ncol(pathway))
pwy_ave_relab <- pwy_ave_relab %>% 
  filter(pathway %in% path_all$path) %>% 
  arrange(desc(ave_relab))
head(pwy_ave_relab)
(top10_pathway <- head(pwy_ave_relab$pathway,10))

pathway <- pathway*100
summary(colSums(pathway,na.rm = T))
pathway <- pathway + 0.1
# enzyme <- enzyme/100
pathway <- pathway %>%
  mutate(across(everything(),log10)) #log10
class(pathway)
pathway <- filter(pathway,rownames(pathway) %in% pwy_meta$feature)

summary(pathway$S112)
table(colSums(is.na(pathway)))

## top 10 pathway
# order pathway by relative abundance

pathway_top10 <- pathway %>% 
  rownames_to_column(var = "feature") %>%
  filter(feature %in% top10_pathway) %>% 
  left_join(select(pwy_relab,c("feature","pathname","pathway")),by="feature") %>% 
  select(c("feature","pathname","pathway"),everything())
dim(pathway_top10)
(pathway_top10[1:6])
pathway_top10$pathname <- Hmisc::capitalize(pathway_top10$pathname)

pathway_top10$feature <- factor(pathway_top10$feature,top10_pathway)
pathway_top10 <- pathway_top10 %>% 
  arrange(desc(feature))
pathway_top10[1:6]
summary(colSums(pathway_top10[,4:ncol(pathway_top10)]))
colnames(pathway_top10[,4:ncol(pathway_top10)])[colSums(pathway_top10[,4:ncol(pathway_top10)])== -30]

pathway_top10_tf <- pathway_top10 %>% 
  select(-c("feature","pathname")) %>% 
  column_to_rownames(var = "pathway") %>% 
  sjmisc::rotate_df()
colSums(is.na(pathway_top10_tf))

pathway_top10_long <- 
  gather(pathway_top10,key = "samples",value = "re_abundance",4:ncol(pathway_top10))
head(pathway_top10_long)
pathway_top10_long$pathname <- 
  factor(pathway_top10_long$pathname,levels = pathway_top10$pathname)

pathway_top10_long$samples <- 
  factor(pathway_top10_long$samples,sample_order)

## pathway heatmap -----
# plot pathway distribution 
pwy_plot <- ggplot(pathway_top10_long, aes(samples, pathname)) +
  geom_tile(aes(fill = re_abundance)) +
  scale_fill_viridis_c(name ="Relative abundance\nlog10 scale", 
                       na.value = "grey80",#"white",
                       direction = -1, option = "G",
                       breaks=c(-1,0,1,2),limits=c(-1,2))+
  theme(legend.position = "none",
        legend.title = element_text(size=12),
        plot.margin = unit(c(1,0,1,5.5),"point"), #top, right, bottom, left
        plot.title = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title=element_blank(),
        axis.text.y=element_text(size = 10),
        axis.text.x = element_blank())+
  guides(fill = guide_colourbar(title.position="top", 
                                title.hjust = 0.5,
                                barheight = 6,
                                barwidth = 1))
print(pwy_plot)

## Stack_area  -------
# top 10 pathway with the top 5 contributing species 
load("pwy_spe_full.RData")
head(colnames(pwy_all_spe))

top10_pwy_all_spe <- pwy_all_spe %>% filter(path %in% top10_pathway)
dim(top10_pwy_all_spe) 
head(colnames(top10_pwy_all_spe))

top10_pwy_all_spe$ave_percentage <- 
  rowMeans(top10_pwy_all_spe[,4:ncol(top10_pwy_all_spe)],na.rm = T)
top10_pwy_all_spe_v2 <- 
  top10_pwy_all_spe[,c("path","genus","species","ave_percentage")]
head(top10_pwy_all_spe_v2)

top10_pwy_n <- top10_pwy_all_spe_v2 %>% 
  group_by(path) %>% 
  summarize(count=n())
top10_pwy_n

top10_pwy_top5_spe <- top10_pwy_all_spe_v2 %>% 
  group_by(path) %>% 
  mutate(rank_per=rank(desc(ave_percentage))) %>% 
  filter(rank_per<=5) %>% 
  arrange(path,rank_per)
top10_pwy_top5_spe
length(unique(top10_pwy_top5_spe$path))
length(unique(top10_pwy_top5_spe$species)) #18

head(top10_pwy_top5_spe)
top10_pwy_top5_spe <- top10_pwy_top5_spe %>% 
  left_join(pwy_relab,by=c("path"="feature"))
top10_pwy_top5_spe$pathname <- Hmisc::capitalize(top10_pwy_top5_spe$pathname)

dim(top10_pwy_top5_spe)
length(unique(top10_pwy_top5_spe$species)) # 18
summary(top10_pwy_top5_spe$ave_percentage)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 3.709   5.210   7.446   8.338  11.042  17.133 
top10_pwy_top5_spe$ave_percentage <- top10_pwy_top5_spe$ave_percentage * 100

top10_pwy_top5_spe$pathname <- factor(top10_pwy_top5_spe$pathname,
                                 levels = pathway_top10$pathname)
top10_pwy_top5_spe$species <- gsub("_"," ", 
                                   gsub("s__","",top10_pwy_top5_spe$species)) 
head(top10_pwy_top5_spe$species)
sort(unique(top10_pwy_top5_spe$species))

library(paletteer) # https://r-charts.com/color-palettes/#discrete
# paletteer_d("ggsci::category20b_d3")
paletteer_d("ggsci::default_igv")

bug_color_value <-   
  c("Alistipes putredinis"="#CC9900FF",
    "Bacteroides dorei"="#FFC20AFF",
    "Bacteroides eggerthii" = "#5050FFFF",
    "Bacteroides fragilis" = "#14FFB1FF", 
    "Bacteroides massiliensis" = "#749B58FF",
    "Bacteroides stercoris"="#99CC00FF",
    "Bacteroides uniformis"="#AE1F63FF",
    "Bacteroides vulgatus"="#CE3D32FF",
    "Bifidobacterium adolescentis" = "#F0E685FF",
    "Bifidobacterium longum" = "#466983FF",
    "Escherichia coli" = "#003399FF", 
    "Eubacterium rectale" = "#008099FF", 
    "Faecalibacterium prausnitzii" = "#7A65A5FF",
    "Klebsiella pneumoniae" = "#5DB1DDFF", 
    "Oscillibacter sp CAG 241" = "#802268FF", 
    "Parabacteroides distasonis" = "#CDDEB7FF", 
    "Prevotella copri"="#339900FF",
    "Ruminococcus torques"="#BA6338FF")
scales::show_col(bug_color_value)
bug_color_value[duplicated(bug_color_value)]

head(top10_pwy_top5_spe)
table(top10_pwy_top5_spe$species %in% names(bug_color_value))
colSums(is.na(top10_pwy_top5_spe[,c("pathname","ave_percentage","species")]))
write.xlsx(top10_pwy_top5_spe,"top10_pwy_top5_spe_tmp.xlsx")

pwy_bar <- ggplot(top10_pwy_top5_spe, aes(x=pathname, y=ave_percentage)) + 
  geom_col(aes(fill=species), position = position_stack(reverse = TRUE)) +
  scale_fill_manual(values=bug_color_value) +
  scale_y_continuous(breaks=c(0,25,50), 
                     limits=c(0,50),
                     labels=c("0%","25%","50%")) +
  guides(fill=guide_legend(ncol =1)) +
  coord_flip() +
  labs(y="") +
  theme_bw() +
  theme(panel.border = element_rect(size=0.2), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_line(size = 0.3),
        axis.text.x = element_blank(),
        axis.title.x = element_text(hjust=0.5, size=10, colour = "black"),
        axis.text.y=element_blank(),
        axis.title.y=element_blank(),
        plot.margin = unit(c(1,5.5,1,1),"point"), #top, right, bottom, left
        legend.position = "none",
        legend.text = element_text(colour = "black", size = 14),
        legend.title = element_blank())
print(pwy_bar)

blank <- ggplot()+theme_void()
pdf("fig_s3_pathway.pdf", width = 15, height = 5, onefile = F) 
egg::ggarrange(pwy_plot, pwy_bar, nrow = 1, ncol =  2, widths = c(5, 1))
dev.off()

## legend --------
pwy_plot_lgd <- ggplot(pathway_top10_long, aes(samples, pathname)) +
  geom_tile(aes(fill = re_abundance)) +
  scale_fill_viridis_c(name ="Relative abundance\nlog10 scale", na.value = "white",
                       direction = -1, option = "G",
                       breaks=c(-1,0,1,2),limits=c(-1,2))+
  theme(legend.position = "left",
        legend.text = element_text(size=12),
        legend.title = element_text(size=12), 
        plot.title = element_blank(),
        axis.title=element_blank(),
        axis.text.y=element_text(size = 16),
        axis.text.x = element_blank())+
  guides(fill = guide_colourbar(title.position="top", 
                                title.hjust = 0.5,
                                barheight = 6,
                                barwidth = 1))
pwy_plot_lgd
pwy_plot_legend <- get_legend(pwy_plot_lgd)

pdf("pwy_plot_legend.pdf", width = 3, height = 3, onefile = F) 
grid.draw(pwy_plot_legend)
dev.off() 

# pwy bar legend
pwy_bar_tmp <- ggplot(top10_pwy_top5_spe, aes(x=pathname, y=ave_percentage)) + 
  geom_hline(yintercept = 0, linetype="dashed", size=1.2) +
  geom_col(aes(fill=species), position = position_stack(reverse = TRUE)) +
  scale_fill_manual(values=bug_color_value) +
  guides(fill=guide_legend(ncol =1)) +
  coord_flip() +
  labs(y="") +
  theme_bw() +
  theme(panel.border = element_rect(size=1.2), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.ticks = element_line(colour = "black", size = 1.2),
        axis.text.x = element_text(size=16, colour = "black"),
        axis.title.x = element_text(hjust=0.5, size=16, colour = "black"),
        axis.text.y=element_blank(),
        axis.title.y=element_blank(),
        legend.position = "right",
        legend.text = element_text(colour = "black", size = 14,face = "italic"),
        legend.title = element_blank())
pwy_bar_legend <- get_legend(pwy_bar_tmp)

pdf("pwy_bar_legend.pdf", width = 5, height = 5, onefile = F) # Open a new pdf file
grid.draw(pwy_bar_legend)
dev.off()

# EC Distribution ------------
head(ecs_meta$feature)
dim(enzyme) #1473/8122
head(colnames(enzyme))
head(rownames(enzyme))
summary(colSums(enzyme,na.rm = T))
head(enzy_all)
ecs_ave_relab <- data.frame(ec=rownames(enzyme),
                            ave_relab=rowSums(enzyme,na.rm=T)/ncol(enzyme))
ecs_ave_relab <- ecs_ave_relab %>% 
  filter(ec %in% enzy_all$enzyme) %>% 
  arrange(desc(ave_relab))
(top10_ec <- head(ecs_ave_relab$ec,10))

enzyme <- filter(enzyme,rownames(enzyme) %in% ecs_meta$feature)

head(colSums(enzyme,na.rm = T))
enzyme <- enzyme*100
enzyme <- enzyme + 0.1
# enzyme <- enzyme/100
enzyme <- enzyme %>%
  mutate(across(everything(),log10))
class(enzyme)

enzyme_top10 <- enzyme[top10_ec,] %>% 
  rownames_to_column(var = "enzyme") 
(ec_order <- factor(enzyme_top10$enzyme,levels = enzyme_top10$enzyme))
head(colnames(enzyme_top10))
tail(colnames(enzyme_top10)) 
enzyme_top10_long <- 
  gather(enzyme_top10,key = "ID",value = "reabun",S112:HRS049423)
head(enzyme_top10_long)

enzyme_top10_long$ID <- factor(enzyme_top10_long$ID,levels = sample_order)
enzyme_top10_long$enzyme <- 
  factor(enzyme_top10_long$enzyme,levels = ec_order)
levels(enzyme_top10_long$enzyme)

colSums(is.na(enzyme_top10_long))
summary(enzyme_top10_long$reabun)
table(enzyme_top10_long$reabun=="-Inf")
table(enzyme_top10_long$reabun == -3)
table(grepl("[a-z]",enzyme_top10_long$reabun))
enzyme_top10_long$reabun[grepl("[a-z]",enzyme_top10_long$reabun)]
min(enzyme_top10_long$reabun)

ec_tmp <- enzyme_top10_long[enzyme_top10_long$enzyme=="EC3.2.1.23: Beta-galactosidase",]

## EC heatmap ------
# double-check include EC
unique(enzyme_top10_long$enzyme)

ecs_plot <- ggplot(enzyme_top10_long, aes(ID,enzyme)) +
  geom_tile(aes(fill = reabun)) +
  scale_fill_viridis_c(name ="Relative abundance\nlog10 scale", na.value = "grey80",
                       direction = -1, option = "F",
                       breaks=c(-1,0,1,2),limits=c(-1,2))+
  scale_y_discrete(limits = rev(levels(ec_order)),
                   labels = c( 
                     "EC6.3.5.5: Carbamoyl-phosphate synthase (glutamine-hydrolyzing)",
                     "EC3.1.11.6: Exodeoxyribonuclease VII",
                     "EC7.1.2.2: H+-transporting two-sector ATPase",
                     "EC1.6.5.11: NADH dehydrogenase (quinone)",
                     "EC2.7.7.6: DNA-directed RNA polymerase",
                     "EC3.2.1.23: Beta-galactosidase",
                     "EC2.7.7.7: DNA-directed DNA polymerase",
                     "EC5.2.1.8: Peptidylprolyl isomerase",
                     "EC3.6.4.12: DNA helicase",
                     "EC2.7.13.3: Histidine kinase")) +
  theme(legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.position = "none",
        plot.margin = unit(c(1,0,1,5.5),"point"), #top, right, bottom, left
        plot.title = element_blank(),
        axis.title=element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y=element_text(size=10), 
        axis.text.x = element_blank())
print(ecs_plot)

## Stack_area  -------
# top 10 pathway with the top 5 contributing species 
load("ecs_spe_full.RData")
head(colnames(ecs_all_spe))
top10_ecs_all_spe <- ecs_all_spe %>% filter(enzy %in% top10_ec)
dim(top10_ecs_all_spe) 
head(colnames(top10_ecs_all_spe))

top10_ecs_all_spe$ave_percentage <- 
  rowMeans(top10_ecs_all_spe[,4:ncol(top10_ecs_all_spe)],na.rm = T)
top10_ecs_all_spe_v2 <- 
  top10_ecs_all_spe[,c("enzy","genus","species","ave_percentage")]
head(top10_ecs_all_spe_v2)
table(top10_ecs_all_spe_v2$enzy)

top10_ecs_top5_spe <- top10_ecs_all_spe_v2 %>% 
  group_by(enzy) %>% 
  mutate(rank_per=rank(desc(ave_percentage))) %>% 
  filter(rank_per<=5) %>% 
  arrange(enzy,rank_per)
top10_ecs_top5_spe
length(unique(top10_ecs_top5_spe$enzy)) 
length(unique(top10_ecs_top5_spe$species)) #9
top10_ecs_top5_spe$species <- 
  gsub("_"," ",gsub("s__","",top10_ecs_top5_spe$species)) 

length(unique(top10_pwy_top5_spe$species)) # 18
length(unique(c(top10_ecs_top5_spe$species,top10_pwy_top5_spe$species)))
summary(top10_ecs_top5_spe$ave_percentage)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.03134 0.04986 0.06703 0.07378 0.08584 0.16818
top10_ecs_top5_spe$ave_percentage <- top10_ecs_top5_spe$ave_percentage * 100

top10_ecs_top5_spe$enzy <- 
  factor(top10_ecs_top5_spe$enzy, levels = ec_order)
unique(top10_ecs_top5_spe$species)

(species <- sort(union(unique(top10_pwy_top5_spe$species),
                      unique(top10_ecs_top5_spe$species))))
bug_color_value <-   
  c("Alistipes putredinis"="#CC9900FF",
    "Bacteroides dorei"="#FFC20AFF",
    "Bacteroides eggerthii" = "#5050FFFF",
    "Bacteroides fragilis" = "#14FFB1FF", 
    "Bacteroides massiliensis" = "#749B58FF",
    "Bacteroides stercoris"="#99CC00FF",
    "Bacteroides uniformis"="#AE1F63FF",
    "Bacteroides vulgatus"="#CE3D32FF",
    "Bifidobacterium adolescentis" = "#F0E685FF", 
    "Bifidobacterium longum" = "#466983FF", 
    "Escherichia coli" = "#003399FF",
    "Eubacterium rectale" = "#008099FF", 
    "Faecalibacterium prausnitzii" = "#7A65A5FF",
    "Klebsiella pneumoniae" = "#5DB1DDFF", 
    "Oscillibacter sp CAG 241" = "#802268FF", 
    "Parabacteroides distasonis" = "#CDDEB7FF", 
    "Prevotella copri"="#339900FF",
    "Ruminococcus torques"="#BA6338FF")
scales::show_col(bug_color_value)
bug_color_value[duplicated(bug_color_value)]

ecs_bar <- ggplot(top10_ecs_top5_spe, aes(x=enzy, y=ave_percentage)) + 
  geom_col(aes(fill=species), position = position_stack(reverse = TRUE)) +
  scale_fill_manual(values=bug_color_value) +
  scale_y_continuous(breaks=c(0,25,50), 
                     limits=c(0,50),
                     labels=c("0%","25%","50%")) +
  guides(fill=guide_legend(ncol =1)) +
  coord_flip() +
  labs(y=expression("\n\n Mean difference\n in contribution %")) +
  theme_bw() +
  theme(panel.border = element_rect(size=0.2), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.ticks.x = element_line(size = 0.3),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(size=8,colour = "black"), 
        axis.title.x = element_text(hjust=0.5, vjust=(-21), size=16, colour = "black"),
        axis.text.y=element_blank(),
        axis.title.y=element_blank(),
        plot.margin = unit(c(1,5.5,1,1),"point"), 
        legend.position = "none",
        legend.text = element_text(colour = "black", size = 14),
        legend.title = element_blank())
print(ecs_bar)

pdf("fig_s3_ec.pdf", width = 15, height = 5, onefile = F) 
egg::ggarrange(ecs_plot, ecs_bar, nrow=1, ncol=2,widths = c(5,1))
dev.off() 

## legend --------
ecs_plot_lgd <- ggplot(enzyme_top10_long, aes(ID,enzyme)) +
  geom_tile(aes(fill = reabun)) +
  scale_fill_viridis_c(name ="Relative abundance\nlog10 scale", na.value = "white",
                       direction = -1, option = "F",
                       breaks=c(-1,0,1,2),limits=c(-1,2))+
  scale_y_discrete(limits = rev(levels(ec_order)),
                   labels = c( "EC6.3.5.5: Carbamoyl-phosphate synthase (glutamine-hydrolyzing)",
                               "EC3.1.11.6: Exodeoxyribonuclease VII",
                               "EC7.1.2.2: H+-transporting two-sector ATPase",
                               "EC1.6.5.11: NADH dehydrogenase (quinone)",
                               "EC2.7.7.6: DNA-directed RNA polymerase",
                               "EC3.2.1.23: Beta-galactosidase",
                               "EC2.7.7.7: DNA-directed DNA polymerase",
                               "EC5.2.1.8: Peptidylprolyl isomerase",
                               "EC3.6.4.12: DNA helicase",
                               "EC2.7.13.3: Histidine kinase")) +
  theme(legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.position = "left",
        plot.title = element_blank(),
        axis.title=element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y=element_text(size=10), #, face = "italic"
        axis.text.x = element_blank())+
  guides(fill = guide_colourbar(title.position="top", 
                                title.hjust = 0.5,
                                barheight = 6,
                                barwidth = 1))
ecs_plot_lgd
ecs_plot_legend <- get_legend(ecs_plot_lgd)

pdf("ecs_plot_legend.pdf", width = 3, height = 3, onefile = F) 
grid.draw(ecs_plot_legend)
dev.off() 

# ec bar legend
ecs_bar_lgd <- ggplot(top10_ecs_top5_spe, aes(x=enzy, y=ave_percentage)) + 
  geom_hline(yintercept = 0, linetype="dashed", size=1.2) +
  geom_col(aes(fill=species), position = position_stack(reverse = TRUE)) +
  scale_fill_manual(values=bug_color_value) +
  guides(fill=guide_legend(ncol =1)) +
  coord_flip() +
  labs(y=expression("\n\n Mean difference\n in contribution %")) +
  theme_bw() +
  theme(panel.border = element_rect(size=1.2), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.ticks = element_line(colour = "black", size = 1.2),
        axis.text.x = element_text(size=16, colour = "black"),
        axis.title.x = element_text(hjust=0.5, vjust=(-21), size=16, colour = "black"),
        axis.text.y=element_blank(),
        axis.title.y=element_blank(),
        legend.position = "right",
        legend.text = element_text(colour = "black", size = 14, face = "italic"),
        legend.title = element_blank())
ecs_bar_lgd
ecs_bar_legend <- get_legend(ecs_bar_lgd)

pdf("ecs_bar_legend.pdf", width = 5, height = 5, onefile = F) 
grid.draw(ecs_bar_legend)
dev.off()

# Combine all plots -------
blank <- ggplot()+theme_void()

pdf(paste0("fig_s3_",date_log,".pdf"), 
    width = 12, height = 8, onefile = F)
egg::ggarrange(p_study,blank,
               p_t2d,blank,
               p_spe,blank,
               pwy_plot,pwy_bar,
               ecs_plot,ecs_bar,
               nrow = 5,ncol = 2,
               widths = c(5,1),
               heights = c(1,1,10,10,10))

dev.off()

