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
library(sjmisc)

library(ggtree)
library(ggtreeExtra)
library(phyloseq)
library(dplyr)
library(tidytree)
library(ggnewscale)
library(ggtreeExtra)
library(ggtree)
library(treeio)
library(tidytree)
library(ggstar)
library(ggplot2)
library(ggnewscale)
library(TDbook)

select <- dplyr::select

pca = function(mat) {
  centered_mat = scale(mat, scale = FALSE)
  svd_res = svd(centered_mat)
  eigs = svd_res$d^2
  tenth_max_k = sum(eigs > (eigs[1]*.1)) + 1
  message(paste0("k = ", tenth_max_k))
  d = diag(svd_res$d)[,1:tenth_max_k]
  svd_res$u %*% d
}

## for RStudio ------
cmdstanr::set_cmdstan_path(path="/n/home05/zmei/.cmdstan/cmdstan-2.30.1")

top_species <-
  read.table("species_list/species_for_strain.txt",quote="\"", comment.char="")
(top_species <- top_species$V1)
length(top_species)

head(top_species)
top_species2 <- str_replace(top_species,"s__","")
(top_species2 <- top_species2[top_species2!="Olsenella_scatoligenes"])
(speName <- top_species2[30]) #39; 30;

## filtered gene table -----
dat_filtered <- 
  fread(paste0("anpan_batch_bmi_minmax5_2023-04-04/filter_stats/filtered_",
  speName,".tsv.gz"))
head(colnames(dat_filtered));dim(dat_filtered) #1873/10285

## Metadata ------
load("metadata_T2D.RData")
head(metadata_T2D);dim(metadata_T2D) #4117/27
table(metadata_T2D$study)
metafile <- filter(metadata_T2D, sample_id %in% dat_filtered$sample_id)
head(metafile);dim(metafile) #1873/27

# phylo_effect posterior
phylo_effect <- 
  read.xlsx(paste0("pglmm_phylo_effect_posterior/phylo_effect_metf_",
                   speName,".xlsx"))
head(phylo_effect)
rownames(metafile) <- metafile$sample_id

## gene pre/abs table -----
genefile <- dat_filtered %>% 
  select(-c("age","sex","study","T2D")) %>% 
  column_to_rownames(var = "sample_id")
head(genefile[1:6])
genefile[genefile=="TRUE"] <- 1
genefile[genefile=="FALSE"] <- 0
summary(colSums(genefile))

## Build tree  -------
pca_df <- pca(genefile)
dim(pca_df) #1201/4
head(rownames(pca_df))
rownames(pca_df) <- rownames(genefile)

distfile <- vegdist(pca_df, method="euclidean")
class(distfile)

treefile <- nj(distfile) |>
  ape::ladderize()
tree_df <- as_tibble(treefile)
tree_df <- tree_df %>% 
  left_join(metafile,by=c("label"="sample_id")) 
head(tree_df)
summary(tree_df$node)
treefile2 <- as.treedata(tree_df)
treefile2

## plot tree -----
study_col <- c(
  `SOL (USA)`= "#CC9900",
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

p1 <- ggtree(treefile2,
             layout = "fan", 
             size = 0.15, 
             open.angle = 10) +
  xlim(-15,NA)+
  geom_tippoint(mapping = aes(color=T2D),
                size=0.5,
                show.legend = T)+
  scale_color_manual(breaks = c("FALSE","TRUE"),
                     values = c("#336fac","#b82929"))
p1

t2d_lgd <- ggpubr::get_legend(p1)
pdf("tree_plot/tree_legend_t2d.pdf", width=6, height=6, onefile = TRUE)
grid::grid.draw(t2d_lgd)
dev.off()

(p2 <- p1 + geom_fruit(
  geom = geom_tile,
  aes(x=0,y=label,fill=study),
  width = 4,
  offset = 0.08,
  show.legend = T,
  inherit.aes = F)+
  scale_fill_manual(values = study_col))

study_lgd <- ggpubr::get_legend(p2)
pdf("tree_plot/tree_legend_study.pdf", width=6, height=6, onefile = TRUE)
grid::grid.draw(study_lgd)
dev.off()

# add phylo_effect posterior
head(metafile)
head(phylo_effect)
summary(phylo_effect$mean)

(p3 <- p2 + new_scale_fill()+
  geom_fruit(data = phylo_effect,
  geom = geom_tile,
  aes(x=0,
      y=sample_id,
      fill=mean),
  width = 4,
  offset = 0.12,
  show.legend = T,
  inherit.aes = F)+
  scale_fill_viridis_c(option = "G",direction = -1,
                       name="Posterior mean phylogenetic effects"))
post_lgd <- ggpubr::get_legend(p3)
pdf(paste0("tree_plot/tree_legend_posterior_",speName,".pdf"),
    width=4, height=8, onefile = TRUE)
grid::grid.draw(post_lgd)
dev.off()
##find  node and parent--------
col_sec <- c(pal_d3('category20')(20),'gray')
(p_tree <- ggtree(treefile2,
                  layout = "fan", 
                  open.angle = 10, 
                  size = 0.15
                  ) +
    xlim(-15,NA)+
    geom_hilight(node=1410,type = "gradient",
                 extendto=40, alpha=0.15, fill=col_sec[1], color="grey40",
                 size=0.15)+
    geom_hilight(node=1406,type = "gradient",
                 extendto=40, alpha=0.15, fill=col_sec[12], color="grey40",
                 size=0.15)+
    geom_hilight(node=1407,type = "gradient",
                 extendto=40, alpha=0.15, fill=col_sec[3], color="grey40",
                 size=0.15)+
    geom_hilight(node=1485,type = "gradient",
                 extendto=40, alpha=0.2, fill=col_sec[5], color="grey40",
                 size=0.15)+
    geom_tippoint(mapping = aes(color=T2D),
                  size=0.4,
                  show.legend = F)+
    scale_color_manual(breaks = c("FALSE","TRUE"),
                       values = c("#336fac","#b82929"))+
    geom_fruit(
      geom = geom_tile,
      aes(x=0,y=label,fill=study),
      width = 4,
      offset = 0.12,
      show.legend = F,
      inherit.aes = F)+
    scale_fill_manual(values = study_col)+
    new_scale_fill()+
    geom_fruit(data = phylo_effect,
               geom = geom_tile,
               aes(x=0,
                   y=sample_id,
                   fill=mean),
               width = 4,
               offset = 0.12,
               show.legend = F,
               inherit.aes = F)+ 
    scale_fill_viridis_c(option = "G",direction = -1,
                         name="Posterior mean phylogenetic effects"))
pdf(paste0("tree_plot/circular_tree_total_filtered_",
           speName,"_V3_",date_log,".pdf"),
    width=4, height=4, onefile = TRUE)
p_tree
dev.off()

