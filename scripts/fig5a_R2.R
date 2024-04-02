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

fig5_path <- "/n/holylabs/LABS/dwang_lab/Lab/T2D_meta/3.analysis/R2/fig5/"

# T2D status species----
t2d_depleted_spe <- c("s__Bacteroides_plebeius",
                      "s__Butyrivibrio_crossotus",
                      "s__Coprococcus_eutactus",
                      "s__Firmicutes_bacterium_CAG_238",
                      "s__Roseburia_sp_CAG_182",
                      "s__Ruminococcus_lactaris",
                      "s__Clostridium_sp_CAG_167",
                      "s__Desulfovibrio_piger",
                      "s__Oscillibacter_sp_57_20",
                      "s__Turicibacter_sanguinis")
length(t2d_depleted_spe) #10
t2d_enriched_spe <- c("s__Bacteroides_fragilis",
                      "s__Bifidobacterium_bifidum",
                      "s__Clostridium_citroniae",
                      "s__Dorea_sp_CAG_317",
                      "s__Escherichia_coli",
                      "s__Flavonifractor_plautii",
                      "s__Streptococcus_parasanguinis",
                      "s__Streptococcus_salivarius",
                      "s__Clostridium_bolteae")
length(t2d_enriched_spe) #9
spe_t2d_status <- data.table(
  species=c(t2d_depleted_spe,t2d_enriched_spe),
  t2d_status=c(rep("T2D depleted",10),
               rep("T2D enriched",9))
)
spe_t2d_status
spe_t2d_status$species <- gsub("s__","",spe_t2d_status$species)

# species list ------
species_list <- 
  read.table("species_list/species_for_strain.txt",quote="\"", comment.char="")
(top_species <- species_list$V1)
length(top_species)
top_species2 <- gsub("s__","",top_species)

# GO term category ----
go_info <- fread("gsea/sig_go_term_category.csv", header = T, drop = "V1")
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
re_gene_com1 <- filter(re_gene_com0,q_global<0.25) 
dim(re_gene_com1) #9648
length(unique(re_gene_com1$bug_name)) #64

species_info <- re_gene_com1 %>% 
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

re_gene_com <- re_gene_com1 %>% 
  left_join(species_info,by="bug_name") %>% 
  left_join(spe_t2d_status,by=c("bug_name"="species")) %>%
  filter(sig_gene_n >= 50) %>%
  arrange(desc(sig_gene_n)) %>% 
  arrange(t2d_status)
head(re_gene_com);length(unique(re_gene_com$bug_name)) #14
table(is.na(re_gene_com$t2d_status))
re_gene_com$t2d_status[is.na(re_gene_com$t2d_status)] <- "Non-significant"


re_gene_com$species_info <- factor(re_gene_com$species_info,
                                   levels = unique(re_gene_com$species_info))
length(unique(re_gene_com$bug_name)) #64 #39
re_gene_com$species_update <- gsub("_"," ",re_gene_com$bug_name)
re_gene_com$species_update <-
  factor(re_gene_com$species_update,
         levels = rev(unique(re_gene_com$species_update)))

# *GSEA results --------
re_gsea_tmp <- 
  read.xlsx("gsea/gsea20_tval_r2/gsea_allgene_q.25_Klebsiella_pneumoniae.xlsx")
head(re_gsea_tmp)
re_gsea_tmp$species <- "Klebsiella_pneumoniae"

re_gsea0 <- re_gsea_tmp[0,]
for (i in top_species2){ #species_allgene$V1
  file_tmp <- paste0("gsea/gsea20_tval_r2/gsea_allgene_q.25_",i,".xlsx")
  if(file.exists(file_tmp)){
    re_temp <- read.xlsx(file_tmp)
    re_temp$species <- i
    re_gsea0 <- bind_rows(re_gsea0,re_temp)
  }
}
dim(re_gsea0) #171/10
head(re_gsea0)

# output for supplemental table
summary(re_gsea0$padj)
dim(re_gsea0) #184/10
write.xlsx(re_gsea0,paste0(fig5_path,"/re_gsea_sig_go_term_for_stable.xlsx"))

# 
re_gsea <- re_gsea0 %>% 
  left_join(go_info,by=c("pathway","go_description")) %>% 
  left_join(spe_t2d_status,by=c("species")) %>% 
  filter(size>=5) %>%  # remove gene size <5
  filter(padj < 0.1) 
table(is.na(re_gsea$t2d_status))
re_gsea$t2d_status[is.na(re_gsea$t2d_status)] <- "Non-significant"

re_gsea <- re_gsea %>% 
  arrange(desc(size)) %>% 
  arrange(t2d_status) %>% 
  arrange(category)
summary(re_gsea$size)

re_gsea <- re_gsea %>% 
  mutate(genenum_cat=case_when(
    size<30 ~ 1,
    size>=30&size<60 ~ 2,
    size>=60 ~ 3
  ))
table(re_gsea$genenum_cat)
write.xlsx(re_gsea,paste0(fig5_path,"gsea_re_all.xlsx"))

head(re_gsea);dim(re_gsea) # 76/12
length(unique(re_gsea$go_description)) #43; 40
length(unique(re_gsea$species)) #26; 25
re_gsea$species_update <- gsub("_"," ",re_gsea$species)
re_gsea$go_description <- gsub("\\[BP] ","",re_gsea$go_description)
re_gsea$go_description <- 
  ifelse(re_gsea$go_description=="dTDP-rhamnose biosynthetic process",
         "dTDP-rhamnose biosynthetic process",
         Hmisc::capitalize(re_gsea$go_description))
head(re_gsea$go_description)
re_gsea$go_description[re_gsea$go_description=="Phosphoenolpyruvate-dependent sugar phosphotransferase system"] <- 
  "Sugar phosphotransferase system"  
re_gsea$go_description[re_gsea$go_description=="Lipopolysaccharide core region biosynthetic process"] <- 
  "LPS core region biosynthetic process"

## remove specific GO terms
unique(re_gsea$go_description)
re_gsea <-
  filter(re_gsea,
         !go_description %in% c("Terpenoid biosynthetic process",
                                "Thiamine biosynthetic process",
                                "Protoporphyrinogen IX biosynthetic process",
                                "NAD biosynthetic process",
                                "Coenzyme A biosynthetic process"))

re_gsea$go_id_description <- 
  paste0(re_gsea$pathway,", ",re_gsea$go_description)
re_gsea$go_id_description
unique(re_gsea$go_id_description) #35

re_gsea$go_id_description <- 
  factor(re_gsea$go_id_description,levels = unique(re_gsea$go_id_description))
re_gsea$species_update <- 
  factor(re_gsea$species_update,levels = rev(unique(re_gsea$species_update))) #rev
head(re_gsea)
re_gsea$category <- 
  factor(re_gsea$category,levels = (unique(re_gsea$category))) #rev
length(unique(re_gsea$go_id_description)) #35
length(unique(re_gsea$species_update)) #29

# Overlapping species -----
table(unique(re_gsea$species_update) %in% unique(re_gene_com$species_update))
(spe_overlap <- 
    intersect(re_gsea$species_update, re_gene_com$species_update)) #15

re_gene_com <- filter(re_gene_com,species_update %in% spe_overlap)
dim(re_gene_com);length(unique(re_gene_com$species_update))

re_gsea <- filter(re_gsea,species_update %in% spe_overlap)
dim(re_gsea);length(unique(re_gsea$species_update))
re_gsea$species_update <- 
  factor(re_gsea$species_update,levels = rev(unique(re_gsea$species_update)))
levels(re_gsea$species_update)

levels(re_gene_com$species_update)
re_gene_com$species_update <- 
  factor(re_gene_com$species_update,levels = levels(re_gsea$species_update))
re_gene_com <- arrange(re_gene_com,species_update)
unique(re_gene_com$species_update)

## 1.gene numbers barplot -----
head(species_info)
species_info_long <- species_info %>% 
  filter(bug_name %in% re_gene_com$bug_name) %>% 
  gather(key = gene_ind, value = gene_num, pos_gene_n:neg_gene_n)
head(species_info_long)
length(unique(species_info_long$bug_name))

dat_gene_num <- species_info_long
dat_gene_num$bug_name <- gsub("_"," ",dat_gene_num$bug_name)
dat_gene_num$bug_name <- 
  factor(dat_gene_num$bug_name,levels = (unique(re_gene_com$species_update)))

## cannot be 0 before log transformation
dat_gene_num$gene_num <- 
  ifelse(dat_gene_num$gene_num==0,1,dat_gene_num$gene_num)

p_gene_num <- ggplot(dat_gene_num,aes(x=log10(gene_num),y=bug_name,fill=gene_ind))+
  geom_bar(stat = "identity",position = "dodge",width = 0.7)+ 
  scale_fill_manual(values = c("#1F4690","#B73E3E"))+
  labs(x="T2D-associated genes")+
  theme_bw()+
  theme(axis.title.y = element_blank(),
        axis.text.y = element_text(size = 14,face = "italic",color = "black"), #11
        # axis.text.y = element_blank(),
        axis.text.x = element_text(size = 8),
        # axis.ticks.y = element_blank(),
        plot.margin = unit(c(5.5,1,5.5,0),"point"),
        legend.position = "none"
        )
p_gene_num

## 2.boxplot -----
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

## legend boxplot
p_box_lgd <- 
  ggplot(re_gene_com,aes(y=species_update,x=abs(statistic),color=Direction))+
  geom_boxplot(position = position_dodge(0.8),alpha=0.3)+
  scale_color_manual(values = c("#1F4690","#B73E3E"))+
  labs(x="Absolute t statistic \nfrom anpan gene model")+
  theme_bw()+
  theme(axis.title.y = element_blank(),
        axis.text.x = element_text(size = 11),
        # axis.text.y = element_blank(),
        axis.text.y = element_text(size=11,face = "italic"),
        legend.position = "left") #vjust = 0, 
p_box_lgd
p_box_legend <- ggpubr::get_legend(p_box_lgd)

pdf(paste0(fig5_path,"gene_model_boxplot_legend.pdf"),
    width = 5,height = 5,onefile = F)
grid::grid.draw(p_box_legend)
dev.off()

## 3.GSEA bubble plot ------
length(unique(re_gsea$go_description)) #31
breaklist <- seq(-1,1,by=0.001)
red_blue <- rev(brewer.pal(n=11,name="RdBu"))
scales::show_col(red_blue)
col_red_blue <- colorRampPalette(red_blue)(length(breaklist))

purple_green <- rev(brewer.pal(n=11,name="PRGn"))
scales::show_col(purple_green)
col_purple_green <- colorRampPalette(purple_green)(length(breaklist))

length(unique(re_gsea$go_description))
re_gsea <- re_gsea %>% 
  mutate(padj_annot=case_when(
    padj>=0.05&padj<0.1 ~ "*",
    padj<0.05 ~ "#" ),
    padj_annot1=case_when(
      padj>=0.05&padj<0.1 ~ "*",
      TRUE ~ ""),
    padj_annot2=case_when(
      padj<0.05 ~ "#" , 
      TRUE ~ "")
    )
table(re_gsea$padj_annot) 
table(re_gsea$padj_annot1)
table(re_gsea$genenum_cat)
p_bub <- ggplot(re_gsea,
                aes(x=go_id_description,y=species_update,color=ES,size=genenum_cat))+
  geom_point(alpha=0.7)+
  geom_text(aes(label=padj_annot1),color="black",size=7,nudge_y = -0.2)+
  geom_text(aes(label=padj_annot2),color="black",size=4,nudge_y = 0)+
  scale_size(range = c(8,13),
             breaks = c(1,2,3),
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
p_bub_lgd <- ggplot(re_gsea,
                    aes(x=go_id_description,y=species_update,
                        color=ES,size=genenum_cat))+
  geom_point(alpha=0.8)+
  scale_size(range = c(8,13),
             breaks = c(1,2,3),
             name = "Num. of Genes")+
  geom_text(aes(label=padj_annot1),color="black",size=7,nudge_y = -0.2)+
  geom_text(aes(label=padj_annot2),color="black",size=4,nudge_y = 0)+
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

pdf(paste0(fig5_path,"legend_bubble_plot.pdf"), 
    width = 8, height = 5, onefile = F) # Open a new pdf file
grid::grid.draw(p_bub_legend)
dev.off()

## 4.GO category annotation ------
table(re_gsea$category)
length(unique(re_gsea$category)) #13;14; 12
unique(re_gsea$category)
paletteer::paletteer_d("ggsci::default_igv")
(col_panel <- c(paletteer::paletteer_d("ggsci::default_igv")[1:15],
                "grey70","grey90"))
scales::show_col(col_panel)
length(col_panel)
names(col_panel) <- c(
  "Amino acid metabolism",
  "Bacterial structural components",
  "Cell motility",
  "DNA replication & transcription",
  "Fatty acid metabolism",
  "Genetic Rearrangement", #E.coli heatmap
  "Glucose metabolism",
  "Phage and HGT",
  "Proteolysis", #E.coli heatmap
  "Quorum sensing",
  "DNA methylation", # newly added
  "Signal transduction",
  "Stress response",
  "Virulence and antibiotic resistance",
  "Damaged DNA repair", # newly added
  "Other", #E.coli heatmap
  "Unknown" #E.coli heatmap
  )

(breaks_go <- as.character(unique(re_gsea$go_id_description)))
(labels_go <- ifelse(str_length(breaks_go) < 35,
                     breaks_go,
                     paste0(str_sub(breaks_go,1,30),"...",
                            str_sub(breaks_go,-15,-1))))

p_cat <- ggplot(re_gsea,aes(y=1,x=go_id_description))+
  geom_tile(aes(fill=category),width=1)+ 
  labs(fill="Category")+ 
  scale_fill_manual(values=col_panel) +
  coord_cartesian(expand = FALSE)+
  scale_x_discrete(breaks=breaks_go,labels=labels_go)+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(size = 11,angle = 30,vjust = 1,hjust = 1,
                                   colour = "black"), 
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title = element_blank(),
        plot.margin = unit(c(0,5.5,5.5,5.5),"point"),
        legend.position = "none")
print(p_cat)


## legend
p_cat_lgd <- ggplot(re_gsea,aes(x=1,y=(go_id_description)))+
  geom_tile(aes(fill=category),width=1)+ 
  labs(fill="Category")+ 
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
pdf(paste0(fig5_path,"go_category_legend_",date_log,".pdf"), 
    width = 5, height = 5, onefile = F) # Open a new pdf file
grid::grid.draw(cat_legend)
dev.off()

# Combine fig5a ------
blank <- ggplot()+theme_void()
pdf(paste0(fig5_path,"fig5a_",date_log,".pdf"),
    width = 15.5,height = 8,onefile = F)
egg::ggarrange(
  p_gene_num,p_box,p_bub, 
  blank,blank,p_cat,
  nrow = 2,ncol = 3,
  heights =  c(5,0.1),
  widths = c(1,1,6)
  ) 
dev.off()




