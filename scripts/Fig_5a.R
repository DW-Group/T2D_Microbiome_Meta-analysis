rm(list = ls())
(date_log <- Sys.Date())
setwd("/n/holylabs/LABS/dwang_lab/Lab/T2D_meta/3.analysis/anpan/")
library(tidyverse)
library(readr)
library(openxlsx)
library(data.table)
library(viridis)
library(hrbrthemes)
library(tidyr)
library(RColorBrewer)
select <- dplyr::select

fig6_path <- "/n/holylabs/LABS/dwang_lab/Lab/T2D_meta/3.analysis/dropbox/fig6/"

# T2D status species----
load("/n/holylabs/LABS/dwang_lab/Lab/T2D_meta/3.analysis/0323/species_T2D_main.RData")
head(meta3$meta_fits)
head(meta_3)
meta_3$feature <- gsub("s__","",meta_3$feature)
meta_3 <- meta_3 %>% 
  mutate(t2d_status=case_when(qval.fdr<0.25&coef>0 ~ "T2D enriched",
                              qval.fdr<0.25&coef<0 ~ "T2D depleted",
                              TRUE ~ "Non-Sig"))
table(meta_3$t2d_status)
spe_t2d_status <- select(meta_3,c("feature","t2d_status"))
spe_t2d_status$t2d_status <- 
  factor(spe_t2d_status$t2d_status,
         levels = c("T2D enriched","T2D depleted","Non-Sig"))

# species list ------
species_list <- 
  read.table("species_list/species_for_strain.txt",quote="\"", comment.char="")
(top_species <- species_list$V1)
length(top_species)
top_species2 <- gsub("s__","",top_species)

## species
species_allgene <- 
  fread("gsea/gsea20_tvalue_bmi/species_sig_allgene_20230510.txt",header = F)
head(species_allgene);dim(species_allgene) #40/1

# GO term category ----
go_info <- fread("gsea/sig_go_term_Daniel_v4.csv", header = T, drop = "V1")
head(go_info);dim(go_info)

# *Gene model results -----
re_gene_com0 <- 
  fread(paste0("anpan_batch_bmi_minmax5_2023-04-04/all_bug_gene_terms.tsv.gz"))
head(re_gene_com0)
length(unique(re_gene_com0$bug_name)) #69

table(top_species2 %in% re_gene_com0$bug_name)
top_species2[!top_species2 %in% re_gene_com0$bug_name]
(spe_gene_model <- top_species2[top_species2 %in% re_gene_com0$bug_name])

re_gene_com0[re_gene_com0$bug_name=="Streptococcus_gordonii",]
re_gene_com0$Direction <- 
  ifelse(re_gene_com0$estimate>0,"T2D_enriched","T2D_depleted")
re_gene_com0 <- filter(re_gene_com0,q_global<0.25) 
dim(re_gene_com0)
length(unique(re_gene_com0$bug_name)) #64
write.xlsx(
  re_gene_com0,
  "anpan_batch_bmi_minmax5_2023-04-04/anpan_sig_genes_all_species.xlsx")

species_info <- re_gene_com0 %>% 
  group_by(bug_name) %>% 
  summarise(sig_gene_n=n(),
            pos_gene_n=sum(estimate>0),
            neg_gene_n=sum(estimate<0),
            pos_gene_percent=pos_gene_n*100/sig_gene_n) %>% 
  as.data.frame()
species_info$species_info <-
  paste0(species_info$bug_name," (n=",species_info$sig_gene_n,", %",
         round(species_info$pos_gene_percent,0)," +)")
head(species_info)
summary(species_info$sig_gene_n)

re_gene_com <- re_gene_com0 %>% 
  left_join(species_info,by="bug_name") %>% 
  left_join(spe_t2d_status,by=c("bug_name"="feature")) %>%
  filter(sig_gene_n >= 50) %>%
  arrange(desc(sig_gene_n)) %>% 
  arrange(t2d_status)
head(re_gene_com);length(unique(re_gene_com$bug_name)) #14

re_gene_com$species_info <- factor(re_gene_com$species_info,
                                   levels = unique(re_gene_com$species_info))
length(unique(re_gene_com$bug_name)) #64 #39
re_gene_com$species_update <- gsub("_"," ",re_gene_com$bug_name)
re_gene_com$species_update <-
  factor(re_gene_com$species_update,
         levels = rev(unique(re_gene_com$species_update)))

# *GSEA results --------
re_allgene_tmp <- 
  read.xlsx("gsea/gsea20_tvalue_bmi/gsea_allgene_q.25_Akkermansia_muciniphila.xlsx")
head(re_allgene_tmp)
re_allgene_tmp$species <- "Akkermansia_muciniphila"

re_allgene0 <- re_allgene_tmp[0,]
for (i in species_allgene$V1){
  re_temp <- read.xlsx(paste0("gsea/gsea20_tvalue_bmi/gsea_allgene_q.25_",i,".xlsx"))
  re_temp$species <- i
  re_allgene0 <- bind_rows(re_allgene0,re_temp)
}
dim(re_allgene0) #205/10; 186/10; 175/10; 184/10
head(re_allgene0)
# write.xlsx(re_allgene0,
#            "gsea/gsea20_tvalue_bmi/gsea_results_all_species_with_sig_diff_genes.xlsx")

re_allgene <- re_allgene0 %>% 
  left_join(go_info,by=c("pathway","go_description")) %>%
  left_join(spe_t2d_status,by=c("species"="feature")) %>%
  filter(abs(ES)>=0.5) %>% # remove ES <0.5
  filter(size>=10) %>%  # remove gene size <10
  arrange(desc(size)) %>% 
  arrange(t2d_status) %>%
  arrange(category)
# write.xlsx(re_allgene,
#            "gsea/gsea20_tvalue_bmi/gsea_results_all_species_with_sig_diff_genes_ES0.5_size10.xlsx")

head(re_allgene);dim(re_allgene) #70/11; 62/12; 65/12
length(unique(re_allgene$go_description)) #43; 40
length(unique(re_allgene$species)) #26; 25
re_allgene$species_update <- gsub("_"," ",re_allgene$species)
re_allgene$go_description <- gsub("\\[BP] ","",re_allgene$go_description)
re_allgene$go_description <- 
  ifelse(re_allgene$go_description=="dTDP-rhamnose biosynthetic process",
         "dTDP-rhamnose biosynthetic process",
         Hmisc::capitalize(re_allgene$go_description))
head(re_allgene$go_description)
re_allgene$go_description[re_allgene$go_description=="Phosphoenolpyruvate-dependent sugar phosphotransferase system"] <- 
  "Sugar phosphotransferase system"  
re_allgene$go_description[re_allgene$go_description=="Lipopolysaccharide core region biosynthetic process"] <- 
  "LPS core region biosynthetic process"

## manully remove GO term
unique(re_allgene$go_description)
re_allgene <-
  filter(re_allgene,
         !go_description %in% c("Molybdate ion transport",
                                "DNA topological change",
                                "Nitrogen fixation",
                                "Peptidyl-L-beta-methylthioaspartic acid biosynthetic process from peptidyl-aspartic acid",
                                "Seryl−tRNA aminoacylation",
                                "Thiamine biosynthetic process",
                                "Cell redox homeostasis",
                                "Histidine biosynthetic process"))

re_allgene$go_id_description <- 
  paste0(re_allgene$pathway,", ",re_allgene$go_description)
re_allgene$go_id_description


re_allgene$go_id_description <- 
  factor(re_allgene$go_id_description,levels = unique(re_allgene$go_id_description))
re_allgene$species_update <- 
  factor(re_allgene$species_update,levels = rev(unique(re_allgene$species_update))) #rev
head(re_allgene)
re_allgene$category <- 
  factor(re_allgene$category,levels = (unique(re_allgene$category))) #rev

length(unique(re_allgene$go_id_description)) #35
length(unique(re_allgene$species_update)) #24

# Overlapping species -----
table(unique(re_allgene$species_update) %in% unique(re_gene_com$species_update))
(spe_overlap <- intersect(re_allgene$species_update,
                          re_gene_com$species_update))

re_gene_com <- filter(re_gene_com,species_update %in% spe_overlap)
dim(re_gene_com);length(unique(re_gene_com$species_update))

re_allgene <- filter(re_allgene,species_update %in% spe_overlap)
dim(re_allgene);length(unique(re_allgene$species_update))

re_allgene$species_update <- 
  factor(re_allgene$species_update,levels = rev(unique(re_allgene$species_update)))
levels(re_allgene$species_update)

levels(re_gene_com$species_update)
re_gene_com$species_update <- 
  factor(re_gene_com$species_update,levels = levels(re_allgene$species_update))
re_gene_com <- arrange(re_gene_com,species_update)
unique(re_gene_com$species_update)

## 1.T2D status gene model -----
re_t2d_plot_gene <- spe_t2d_status %>% 
  filter(feature %in% re_gene_com$bug_name)
head(re_t2d_plot_gene);dim(re_t2d_plot_gene) #26; 67
re_t2d_plot_gene$species_update <- gsub("_"," ",re_t2d_plot_gene$feature)
re_t2d_plot_gene$species_update <- 
  factor(re_t2d_plot_gene$species_update,
         levels = (levels(re_gene_com$species_update))) #rev
p_t2d_status <- 
  ggplot(re_t2d_plot_gene,aes(x=1,y=species_update,fill=t2d_status))+
  geom_col(width = 1)+
  scale_fill_manual(
    breaks = c("Non-Sig","T2D depleted","T2D enriched"),
    values = c("grey","#5050FFFF","#CE3D32FF")
    )+
  coord_cartesian(expand = FALSE)+
  theme_bw()+
  theme(axis.title = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 14,face = "italic",color = "black"), 
        axis.ticks.x = element_blank(),
        plot.margin = unit(c(5.5,0,5.5,1),"point"),
        legend.position = "none")+
  theme(panel.grid.major = element_blank())
p_t2d_status

## legend 
p_t2d_lgd <- ggplot(re_t2d_plot_gene,aes(x=1,y=species_update,fill=t2d_status))+
  geom_col(width = 1)+
  scale_fill_manual(
    breaks = c("T2D enriched","T2D depleted","Non-Sig"),
    values = c("#CE3D32FF","#5050FFFF","grey"),
    labels = c("T2D enriched","T2D depleted","Non-significant"),
    name = "Species")+
  coord_cartesian(expand = FALSE)+
  theme_bw()+
  theme(axis.title = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 12,face = "italic",color="black"),
        axis.ticks.x = element_blank(),
        legend.position = "right")+
  theme(panel.grid.major = element_blank())
p_t2d_lgd
p_t2d_legend <- ggpubr::get_legend(p_t2d_lgd)

pdf(paste0(fig6_path,"species_t2d_status_legend.pdf"),
    width = 5,height = 5,onefile = F)
grid::grid.draw(p_t2d_legend)
dev.off()

## 2.gene numbers barplot -----
head(species_info)
species_info_long <- species_info %>% 
  filter(bug_name %in% re_gene_com$bug_name) %>% 
  gather(key = gene_ind, value = gene_num, pos_gene_n:neg_gene_n)
head(species_info_long)
length(unique(species_info_long$bug_name))

dat_gene_num <- species_info_long
dat_gene_num$bug_name <- 
  factor(dat_gene_num$bug_name,levels = (unique(re_gene_com$bug_name)))

## cannot be 0 before log transformation
dat_gene_num$gene_num <- 
  ifelse(dat_gene_num$gene_num==0,1,dat_gene_num$gene_num)

p_gene_num <- ggplot(dat_gene_num,aes(x=log10(gene_num),y=bug_name,fill=gene_ind))+
  geom_bar(stat = "identity",position = "dodge",width = 0.7)+ #position_dodge(0.8)
  scale_fill_manual(values = c("#1F4690","#B73E3E"))+
  labs(x="T2D-associated genes")+
  theme_bw()+
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 8),
        axis.ticks.y = element_blank(),
        plot.margin = unit(c(5.5,1,5.5,0),"point"),
        legend.position = "none"
        )
p_gene_num

## 3.boxplot -----
p_box <- 
  ggplot(re_gene_com,aes(y=species_update,x=abs(statistic),color=Direction))+
  geom_boxplot(position = position_dodge(0.7),alpha=0.3,width = 0.7)+
  scale_color_manual(values = c("#1F4690","#B73E3E"))+
  labs(x="Absolute t statistic")+ #t statistics; effect size
  theme_bw()+
  theme(axis.title.y = element_blank(),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_blank(),
        plot.margin = unit(c(5.5,1,5.5,1),"point"),
        legend.position = "none")
print(p_box)

## legend
## boxplot
p_box_lgd <- 
  ggplot(re_gene_com,aes(y=species_update,x=abs(statistic),color=Direction))+
  geom_boxplot(position = position_dodge(0.8),alpha=0.3)+
  scale_color_manual(values = c("#1F4690","#B73E3E"))+
  labs(x="Absolute t statistic \nfrom anpan gene model")+
  theme_bw()+
  theme(axis.title.y = element_blank(),
        axis.text.x = element_text(size = 11),
        axis.text.y = element_text(size=11,face = "italic"),
        legend.position = "left") 
p_box_lgd
p_box_legend <- ggpubr::get_legend(p_box_lgd)

pdf(paste0(fig6_path,"gene_model_boxplot_legend.pdf"),
    width = 5,height = 5,onefile = F)
grid::grid.draw(p_box_legend)
dev.off()

## 4.sample size barplot -----
(i=top_species2[2])
re_tmp <- fread("anpan_batch_bmi_minmax5_2023-04-04/filter_stats/filtered_Akkermansia_muciniphila.tsv.gz")
dim(re_tmp)
head(colnames(re_tmp))
table(re_tmp$T2D)
re_tmp$species <- i
re_tmp <- select(re_tmp,c("species","T2D"))
re_tmp_n <- re_tmp %>% 
  group_by(species,T2D) %>% 
  summarise(n=n(),.groups = "keep")
head(re_tmp_n)

(re_samplesize <- re_tmp_n[0,])
for (i in spe_gene_model) {
  re_tmp <- 
    fread(paste0("anpan_batch_bmi_minmax5_2023-04-04/filter_stats/filtered_",i,".tsv.gz"))
  re_tmp$species <- i
  re_tmp <- select(re_tmp,c("species","T2D"))
  re_tmp_n <- re_tmp %>% 
    group_by(species,T2D) %>% 
    summarise(n=n(),.groups = "keep")
  re_samplesize <- bind_rows(re_samplesize,re_tmp_n)
}
head(re_samplesize);length(unique(re_samplesize$species))
re_samplesize$T2D <- ifelse(re_samplesize$T2D==TRUE,"T2D","Control")

re_samplesize <- filter(re_samplesize,species %in% re_gene_com$bug_name)
head(re_samplesize);length(unique(re_samplesize$species)) #39

re_samplesize$species <- 
  factor(re_samplesize$species,levels = (unique(re_gene_com$bug_name))) #rev

scales::show_col(c("#0E8388","#DF7857","white",
                   "#749B58", "#CC9900", "#C75127",
                   "#1F4690","#B73E3E"))
p_sample_size <- ggplot(re_samplesize,aes(y=species,x=n,fill=T2D))+
  geom_bar(stat = "identity",position = "dodge",width = 0.7)+
  scale_fill_manual(values = c("#72b0da","#C75127"))+  
  scale_x_continuous(breaks = c(0,1000,2000))+ 
  labs(x="Sample size",y="Species")+
  theme_bw()+
  theme(#axis.text.y = element_blank(),
        axis.text.x = element_text(size = 8),
        axis.title.y = element_blank(),
        plot.margin = unit(c(5.5,1,5.5,1),"point"),
        legend.position = "none"
        )
p_sample_size

## legend
p_sample_size_lgd <- ggplot(re_samplesize,aes(y=species,x=n,fill=T2D))+
  geom_bar(stat = "identity",position = "stack")+
  scale_fill_manual(values = c("#72b0da","#C75127"))+ 
  scale_x_continuous(breaks = c(0,1000,2000),position = "top")+
  labs(x="Sample size",y="Species")+
  theme_bw()+
  theme(#axis.text.y = element_blank(),
        axis.text.x = element_text(size = 11),
        axis.title.y = element_blank(),
        plot.margin = unit(c(5.5,2,5.5,0),"point")#,
        # legend.position = "none"
  )
p_sample_size_lgd

p_sample_size_legend <- ggpubr::get_legend(p_sample_size_lgd)

pdf(paste0(fig6_path,"legend_sample_size.pdf"),
    width = 5,height = 5,onefile = F)
grid::grid.draw(p_sample_size_legend)
dev.off()


## 5.GSEA bubble plot ------
breaklist <- seq(-1,1,by=0.001)
red_blue <- rev(brewer.pal(n=11,name="RdBu"))
scales::show_col(red_blue)
col_red_blue <- colorRampPalette(red_blue)(length(breaklist))

purple_green <- rev(brewer.pal(n=11,name="PRGn"))
scales::show_col(purple_green)
col_purple_green <- colorRampPalette(purple_green)(length(breaklist))

p_bub <- ggplot(re_allgene,
                aes(x=go_id_description,y=species_update,color=ES,size=size))+
  
  geom_point(alpha=0.7)+
  scale_size(range = c(3,10),
             name = "Num. of Genes")+
  labs(x="Species",y="GO term")+
  scale_color_gradientn(colours = col_red_blue,name="Enrichment score")+
  theme_bw()+
  theme(axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title = element_text(size = 13),
        axis.title.x = element_blank(),
        plot.margin = unit(c(5.5,1,5.5,1),"point"),
        legend.position = "none")
print(p_bub)

## legend
p_bub_lgd <- ggplot(re_allgene,
                    aes(x=species_update,y=go_id_description,
                        color=ES,size=size))+
  geom_point(alpha=0.8)+
  scale_size(range = c(3,8),
             name = "Num. of Genes")+
  labs(x="Species",y="GO term")+
  scale_color_gradientn(colours = col_red_blue,name="Enrichment score")+
  theme(#axis.text.y = element_blank(),
    #axis.text.x = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(size = 13, angle = 45,face = "italic",
                               hjust = 1,vjust = 1),
    axis.title = element_text(size = 13),
    axis.title.x = element_blank()#,
    # legend.position = "none"
  )
print(p_bub_lgd)
p_bub_legend <- ggpubr::get_legend(p_bub_lgd)

pdf(paste0(fig6_path,"legend_bubble_plot.pdf"), 
    width = 8, height = 5, onefile = F) # Open a new pdf file
grid::grid.draw(p_bub_legend)
dev.off()


## 6.GO category annotation ------
table(re_allgene$category)
length(unique(re_allgene$category)) #13
unique(re_allgene$category)
paletteer::paletteer_d("ggsci::default_igv")
(col_panel <- c(paletteer::paletteer_d("ggsci::default_igv")[1:14],
                "grey70","grey90"))
scales::show_col(col_panel)
names(col_panel) <- c(
  "Amino acid metabolism",
  "Bacterial structural components",
  "Cell motility",
  # "Damaged DNA repair",
  "DNA replication & transcription",
  "Fatty acid metabolism",
  "Genetic Rearrangement", #E.coli heatmap
  "Glucose metabolism",
  "Phage and HGT",
  "Proteolysis", #E.coli heatmap
  "Quorum sensing",
  "Signal peptide processing",
  "Signal transduction",
  "Stress response",
  "Virulence and antibiotic resistance",
  "Other", #E.coli heatmap
  "Unknown" #E.coli heatmap
  )

(breaks_go <- as.character(unique(re_allgene$go_id_description)))
(labels_go <- ifelse(str_length(breaks_go) < 35,
                     breaks_go,
                     paste0(str_sub(breaks_go,1,25),"...",
                            str_sub(breaks_go,-10,-1))))

p_cat <- ggplot(re_allgene,aes(y=1,x=go_id_description))+
  geom_tile(aes(fill=category),width=1)+ 
  labs(fill="Category")+ 
  scale_fill_manual(values=col_panel) +
  coord_cartesian(expand = FALSE)+
  scale_x_discrete(breaks=breaks_go,labels=labels_go)+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(size = 12,angle = 30,vjust = 1,hjust = 1), #,colour = "black"
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title = element_blank(),
        plot.margin = unit(c(0,5.5,5.5,5.5),"point"),
        legend.position = "none")
print(p_cat)

## legend
p_cat_lgd <- ggplot(re_allgene,aes(x=1,y=(go_id_description)))+
  geom_tile(aes(fill=category),width=1)+ #color=category,
  labs(fill="Category")+ #,color="Category"
  scale_fill_manual(values=col_panel) +
  guides(fill = guide_legend(reverse = FALSE))+  
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.y = element_text(size = 12,colour = "black"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title = element_text(size = 13),
        axis.title.x = element_blank())
print(p_cat_lgd)
cat_legend <- cowplot::get_legend(p_cat_lgd)
pdf(paste0(fig6_path,"go_category_legend_",date_log,".pdf"), 
    width = 5, height = 5, onefile = F) # Open a new pdf file
grid::grid.draw(cat_legend)
dev.off()

## T2D status GSEA -----
re_t2d_plot <- spe_t2d_status %>% 
  filter(feature %in% re_allgene$species)
head(re_t2d_plot);dim(re_t2d_plot) #26; 23
re_t2d_plot$species_update <- gsub("_"," ",re_t2d_plot$feature)
re_t2d_plot$species_update <- 
  factor(re_t2d_plot$species_update,
         levels = (unique(re_allgene$species_update))) #rev
p_t2d_gsea <- 
  ggplot(re_t2d_plot,aes(y=1,x=species_update,fill=t2d_status))+
  geom_col(width = 1)+
  scale_fill_manual(
    breaks = c("Non-Sig","T2D depleted","T2D enriched"),
    values = c("grey","#1F4690","#B73E3E"))+
  coord_cartesian(expand = FALSE)+
  theme_bw()+
  theme(axis.title = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 12,face = "italic",
                                   angle = 45,hjust = 1,vjust = 1),
        axis.ticks.y = element_blank(),
        plot.margin = unit(c(0,5.5,5.5,5.5),"point"),
        legend.position = "none")+
  theme(panel.grid.major = element_blank())
p_t2d_gsea

## legend 
p_t2d_lgd <- ggplot(re_t2d_plot,aes(x=1,y=species_update,fill=t2d_status))+
  geom_col(width = 1)+
  scale_fill_manual(
    breaks = c("T2D enriched","T2D depleted","Non-Sig"),
    values = c("#B73E3E","#1F4690","grey"),
    labels = c("T2D enriched","T2D depleted","Non-significant"),
    name = "T2D status")+
  coord_cartesian(expand = FALSE)+
  theme_bw()+
  theme(axis.title = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 11,face = "italic"),
        axis.ticks.x = element_blank(),
        legend.position = "right")+
  theme(panel.grid.major = element_blank())
p_t2d_lgd
p_t2d_legend <- ggpubr::get_legend(p_t2d_lgd)

pdf(paste0(fig6_path,"species_t2d_status_legend.pdf"),
    width = 5,height = 5,onefile = F)
grid::grid.draw(p_t2d_legend)
dev.off()


# Combine figure ------
blank <- ggplot()+theme_void()
pdf(paste0(fig6_path,"fig_6ac_",date_log,".pdf"),
    width = 16,height = 8,onefile = F)
egg::ggarrange(
  # blank,blank,blank,blank,blank,blank,
  p_t2d_status,blank,p_gene_num,p_box,p_bub, #p_sample_size,
  blank,blank,blank,blank,p_cat, #blank,
  nrow = 2,ncol = 5,
  heights =  c(5,0.1),
  widths = c(0.06,-0.01,1,1,6)) #0.6,
dev.off()




