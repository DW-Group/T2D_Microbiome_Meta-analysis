rm(list = ls())
setwd("/n/holylabs/LABS/dwang_lab/Lab/T2D_meta/3.analysis/anpan/")

library(tidyverse)
library(readr)
library(anpan)
library(ape)
library(vegan)
library(openxlsx)
library(data.table)
select <- dplyr::select
(date_log <- Sys.Date())
fig4_path <- "/n/holylabs/LABS/dwang_lab/Lab/T2D_meta/3.analysis/dropbox/fig5/"

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
gdata::keep(spe_t2d_status,select,date_log,fig4_path,sure = T)

## species list -----
species_list <- 
  read.table("species_list/species_for_strain.txt",quote="\"", comment.char="")
(top_species <- species_list$V1)
length(top_species)
top_species2 <- gsub("s__","",top_species)

# GSEA species list
gsea_allgene <- fread("gsea_allgene_results_2023-03-20.csv")
(gsea_species <- gsea_allgene$species_update)


# PGLMM results -----
species_pglmm <- fread("PGLMM_results/species_pglmm_20230406.txt",header = F)
head(species_pglmm)
length(species_pglmm$V1) #68
top_species[!top_species2 %in% species_pglmm$V1]

re_pglmm_t <- c()
for (i in species_pglmm$V1){
  speName <- i
  re_pglmm_tmp <- fread(paste0("PGLMM_results/PGLMM_re_filtered_",speName,"_sig.3_std.2.csv"))
  re_pglmm_t <- plyr::rbind.fill(re_pglmm_t,re_pglmm_tmp)
}
head(re_pglmm_t);dim(re_pglmm_t) 
re_pglmm_t$elpd_2se <- abs(re_pglmm_t$elpd_diff) - 2*re_pglmm_t$se_diff

## supplemental table ----
spe_sample_size <- read.xlsx("spe_gene_num_check.xlsx") %>% 
  select(bug_name,sample_size,t2d_cases,controls) %>% 
  setnames(c("species","sample_size","t2d_cases","controls"))
head(spe_sample_size)
spe_sample_size[,2:4] <- map_df(spe_sample_size[,2:4],as.numeric)
re_pglmm_output <- re_pglmm_t %>% 
  left_join(spe_sample_size,by="species") %>% 
  select("species",everything())

# write.xlsx(re_pglmm_output,
#            paste0("PGLMM_results/combine_results/PGLMM_results_all_species_",
#                   date_log,".xlsx"))

## sig species base+metf model ------
re_pglmm_metf <- re_pglmm_t %>% 
  filter(model=="base_fit") %>% 
  filter(model_metf=="base+metf") %>%
  filter(elpd_diff< -2) %>%
  # filter(elpd_2se<0) %>% 
  mutate(sig=case_when(elpd_2se>0 ~ "*", TRUE ~ "")) %>% 
  left_join(spe_t2d_status,by=c("species"="feature")) %>% 
  arrange(desc(elpd_diff)) 
(re_pglmm_metf)


re_pglmm_metf$species_update <-
  gsub("_"," ",re_pglmm_metf$species)
re_pglmm_metf$species_update <- 
  factor(re_pglmm_metf$species_update,levels = re_pglmm_metf$species_update)
dim(re_pglmm_metf) #27

p_pglmm_metf <- 
  ggplot(re_pglmm_metf,aes(x= -elpd_diff,y=species_update))+
  geom_errorbar(aes(xmin= -(elpd_diff-se_diff),xmax= -(elpd_diff+se_diff)),
                width=0.2,color="#1F4690")+
  geom_point(color="#1F4690")+
  theme_bw()+
  labs(x="Difference in ELPD")+
  theme(axis.title.x = element_text(size = 13),
        axis.text.x = element_text(size = 11),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        plot.margin = unit(c(5.5,0,5.5,2),"point"))
p_pglmm_metf

## 1.T2D status -----
re_t2d_plot <- spe_t2d_status %>% 
  filter(feature %in% re_pglmm_metf$species)
head(re_t2d_plot);dim(re_t2d_plot) #27
re_t2d_plot$species_update <- gsub("_"," ",re_t2d_plot$feature)
re_t2d_plot$species_update <- 
  factor(re_t2d_plot$species_update,
         levels = re_pglmm_metf$species_update)
p_t2d_status <- ggplot(re_t2d_plot,aes(x=1,y=species_update,fill=t2d_status))+
  geom_col(width = 1)+
  scale_fill_manual(
    breaks = c("Non-Sig","T2D depleted","T2D enriched"),
    values = c("grey","#5050FFFF","#CE3D32FF"))+ 
  coord_cartesian(expand = FALSE)+
  theme_bw()+
  theme(axis.title = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 13,face = "italic",color = "black"),
        axis.ticks.x = element_blank(),
        plot.margin = unit(c(5.5, 2, 5.5, 0), "point"),
        legend.position = "none")
p_t2d_status

## 2.sample size barplot -----
(i=top_species2[2])
re_tmp <- 
  fread(paste0("anpan_batch_bmi_minmax5_2023-04-04/filter_stats/filtered_",
               i,".tsv.gz"))
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
for (i in re_pglmm_metf$species) {
  re_tmp <- 
    fread(paste0("anpan_batch_bmi_minmax5_2023-04-04/filter_stats/filtered_",i,".tsv.gz"))
  re_tmp$species <- i
  re_tmp <- select(re_tmp,c("species","T2D"))
  re_tmp_n <- re_tmp %>% 
    group_by(species,T2D) %>% 
    summarise(n=n(),.groups = "keep")
  re_samplesize <- bind_rows(re_samplesize,re_tmp_n)
}
head(re_samplesize)
re_samplesize$T2D <- ifelse(re_samplesize$T2D==TRUE,"T2D","Control")

re_samplesize$species <- 
  factor(re_samplesize$species,levels = re_pglmm_metf$species)

scales::show_col(c("#0E8388","#DF7857","white",
                   "#749B58", "#CC9900", "#C75127",
                   "#1F4690","#B73E3E"))
theme_bw()$plot.margin
p_sample_size <- ggplot(re_samplesize,aes(y=species,x=n,fill=T2D))+
  geom_bar(stat = "identity",position = "stack")+
  scale_fill_manual(values = c("#72b0da","#C75127"))+ 
  scale_x_continuous(breaks = c(0,1000,2000))+
  labs(x="Sample size",y="Species")+
  theme_bw()+
  theme(axis.text.y = element_blank(),
        axis.text.x = element_text(size = 11,color = "black"),
        axis.title.x = element_text(size = 13,colour = "black"),
        legend.position = "none",
        plot.margin = unit(c(2, 10, 5.5, 5.5), "points"),
        axis.title.y = element_blank())
p_sample_size

## combine plots -----
blank <- ggplot()+theme_void()
pdf(paste0(fig4_path,"fig4a_",date_log,".pdf"), 
    width = 8,height = 5,onefile = F) 
egg::ggarrange(
  p_t2d_status,
  blank,
  p_pglmm_metf,#blank,
  p_sample_size,
  nrow=1, widths = c(0.15,-0.15,3,1.5)) #-0.05,
dev.off()

## legends ----
## t2d status
p_t2d_lgd <- ggplot(re_t2d_plot,aes(x=1,y=species_update,fill=t2d_status))+
  geom_col(width = 1)+
  scale_fill_manual(
    breaks = c("T2D enriched","T2D depleted","Non-Sig"),
    values = c("#CE3D32FF","#5050FFFF","grey"),
    labels = c("T2D enriched","T2D depleted","Non-significant"),
    name = "Species status")+
  theme_bw()+
  theme(axis.title = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 11,face = "italic"),
        axis.ticks.x = element_blank(),
        legend.position = "right")
p_t2d_lgd
p_t2d_legend <- ggpubr::get_legend(p_t2d_lgd)

pdf(paste0(fig4_path,"species_t2d_status_legend.pdf"),
    width = 5,height = 5,onefile = F)
grid::grid.draw(p_t2d_legend)
dev.off()

## sample size
p_sample_size_lgd <- 
  ggplot(re_samplesize,aes(y=species,x=n,fill=T2D))+
  geom_bar(stat = "identity",position = "stack")+
  scale_fill_manual(values = c("#72b0da","#C75127"))+ 
  scale_x_continuous(breaks = c(0,1000,2000))+
  labs(x="Sample size",y="Species")+
  theme_bw()+
  theme(axis.text.y = element_blank(),
        axis.text.x = element_text(size = 10),
        # legend.position = "none",
        axis.title.y = element_blank())
p_sample_size_lgd
p_sample_size_legend <- ggpubr::get_legend(p_sample_size_lgd)

pdf(paste0(fig4_path,"sample_size_legend.pdf"),
    width = 5,height = 5,onefile = F)
grid::grid.draw(p_sample_size_legend)
dev.off()


