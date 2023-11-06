rm(list = ls())
setwd("/n/holylabs/LABS/dwang_lab/Lab/T2D_meta/3.analysis/anpan/")
(date_log <- Sys.Date())
library(tidyverse)
library(readr)
library(anpan)
library(ape)
library(vegan)
library(openxlsx)
library(data.table)
library(fgsea)
library(viridis)
library(hrbrthemes)
library(tidyr)

select <- dplyr::select

species_list <- fread("species_list/species_for_strain.txt",header = F)
(top_species <- species_list$V1)
length(top_species) #70
top_species2 <- gsub("s__","",top_species)

anpan_batch_list <- 
  fread("anpan_batch_bmi_minmax5_2023-04-04/filter_stats/species_filtered_list.txt",header = F)
head(anpan_batch_list);dim(anpan_batch_list) #70
top_species2[!top_species2 %in% anpan_batch_list$V1]

## GO after annotation -----
go_match_uniref <- fread("~/database/go_output_long_20.txt",header = F)

table(str_length(go_match_uniref$V2))
head(go_match_uniref)
dim(go_match_uniref) #77245
colnames(go_match_uniref) <- c("gene_id","go_term")
go_match_uniref <- go_match_uniref %>% 
  separate(go_term,c("go_id","go_description"),":_")
go_match_uniref$go_description <- 
  str_replace_all(go_match_uniref$go_description,"_"," ")

table(duplicated(go_match_uniref$gene_id))
table(duplicated(go_match_uniref))

go_info <- go_match_uniref %>% 
  select(c("go_id","go_description")) %>% 
  distinct()
dim(go_info) #450

## anpan results -----
## for sbatch
speName <- str_replace("VARNAME","s__","")

anpan_re <- 
  fread(paste0("anpan_batch_bmi_minmax5_2023-04-04/model_stats/",speName,"_gene_terms.tsv.gz"))

## result clean -------
dim(anpan_re) #6399
head(anpan_re)
anpan_re <- anpan_re %>% 
  separate(gene,c("gene_id","gene_name"),": ",remove = F)
anpan_re$uniprot_id <- str_replace(anpan_re$gene_id,"UniRef90_","")
head(anpan_re)
anpan_re <- anpan_re %>% 
  left_join(go_match_uniref,by=c("gene_id"))
dim(anpan_re) #6476
head(anpan_re)
table(duplicated(anpan_re$gene_id))
table(is.na(anpan_re$gene_id))
table(is.na(anpan_re$go_id))
length(unique(anpan_re$go_id)) #224

anpan_re_gsea <- filter(anpan_re,!is.na(go_id))
dim(anpan_re_gsea) #14856
table(duplicated(anpan_re_gsea$gene_id)) #77

table(go_match_uniref$gene_id %in% anpan_re$gene_id)
go_match_uniref_inc <- filter(go_match_uniref,gene_id %in% anpan_re$gene_id)
dim(go_match_uniref_inc) #531

## GSEA ------
## format data
anpan_re_gsea_rank <- anpan_re_gsea[order(desc(anpan_re_gsea$statistic)),]
head(anpan_re_gsea_rank)
table(duplicated(anpan_re_gsea_rank$gene_id))

gene_rank <- anpan_re_gsea_rank$statistic
names(gene_rank) <- anpan_re_gsea_rank$gene_id

anpan_re_gsea_genesize_rank <- anpan_re_gsea_rank %>% 
  group_by(go_id) %>% 
  summarize(count=n()) %>% 
  arrange(count)
head(anpan_re_gsea_genesize_rank)
table(anpan_re_gsea_genesize_rank$count)
  
go <- unique(anpan_re_gsea_genesize_rank$go_id)
head(go)
length(go) #223
table(duplicated(go))

pathway <- list()
for (i in seq(length(go))){
  anpan_re_tmp <- filter(anpan_re_gsea_rank,go_id == go[i])
  pathway_tmp <- anpan_re_tmp$gene_id
  pathway[[i]] <- pathway_tmp
}
names(pathway) <- go
head(pathway)

set.seed(123)
gsea_mod <- fgsea(
  pathways = pathway,
  stats = gene_rank,
  minSize = 5,
  maxSize = Inf,
  eps = 1e-10,
  scoreType = c("std", "pos", "neg"),
  nproc = 0,
  gseaParam = 1,
  BPPARAM = NULL,
  nPermSimple = 1000,
  absEps = NULL
)
dim(gsea_mod)


gsea_mod_re <- gsea_mod %>% 
  left_join(go_info, by=c("pathway"="go_id")) %>% 
  as.data.frame()
head(gsea_mod_re)
dim(gsea_mod_re)

## GSEA plot ---
source("~/r_script/t2d_meta/anpan/6.gsea/plotEnrichment2.R")
table(gsea_mod_re$padj<0.25)
summary(gsea_mod_re$padj)

gsea_mod_re_sig <- filter(gsea_mod_re,padj<0.25)
gsea_mod_re_sig
if(nrow(gsea_mod_re_sig)>0){
  write.xlsx(gsea_mod_re_sig,
             paste0("gsea/gsea20_tvalue_bmi/gsea_allgene_q.25_",speName,".xlsx"))
  ggplot(gsea_mod_re_sig,
         aes(x=ES,y=reorder(go_description,desc(ES)),color=padj,size=size))+
    geom_point(alpha=0.5)+
    scale_size(range = c(3,10),name = "Num. of Genes")+
    labs(x="Enrichment Score",y="GO")+
    theme(axis.text.y = element_text(size = 12),
          axis.title = element_text(size = 13))
  ggsave(paste0("gsea/gsea20_tvalue_bmi/gsea_allgene_bmi_q.25_",speName,".png"),
         width = 25,height = 10,units = "cm",dpi = 300)
  
  ## GSEA plot
  gsea_mod_re_sig$pathway
  gsea_mod_re_sig$go_description
  nrow(gsea_mod_re_sig) #5
  for (i in seq(nrow(gsea_mod_re_sig))) {
    i=16
    go_id <- gsea_mod_re_sig$pathway[i]
    (go_name <- gsea_mod_re_sig$go_description[i])
    go_id2 <- gsub(":","_",go_id)
    pdf(paste0("gsea/gsea_plot/gsea_plot_",speName,"_",go_id2,".pdf"),
        width = 8,height = 4)
    print(plotEnrichment2(pathway = pathway[[go_id]],stats = gene_rank,
                         ticksSize = 0.5)+ labs(title=go_name))
    dev.off()
  }
}
head(gsea_mod)
