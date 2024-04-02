# 20230403 
# update the study variable (10 studies instead of 9) 
# and further adjusted for BMI
rm(list = ls())
setwd("/n/holylabs/LABS/dwang_lab/Lab/T2D_meta/3.analysis/anpan/")
(date_log <- Sys.Date())
(lbl <- "_sig.3_std.2")
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

## for sbatch ------
cmdstanr::set_cmdstan_path(path="/n/home05/zmei/apps/CmdStan/cmdstan-2.30.1")
speName <- "VARNAME"; speName <- gsub("s__","",speName)

## for RStudio ------
# cmdstanr::set_cmdstan_path(path="/n/home05/zmei/.cmdstan/cmdstan-2.30.1")
# top_species <-
#   read.table("species_list/species_for_strain.txt",quote="\"", comment.char="")
# (top_species <- top_species$V1)
# length(top_species)
# 
# head(top_species)
# top_species2 <- str_replace(top_species,"s__","")
# head(top_species2)
# (top_species2 <- top_species2[top_species2!="Olsenella_scatoligenes"])
# (speName <- top_species2[39])
# 30; 39; 60


## filtered gene table -----
dat_filtered <- 
  fread(paste0("anpan_batch_bmi_minmax5_2023-04-04/filter_stats/filtered_",
               speName,".tsv.gz"))
head(colnames(dat_filtered));dim(dat_filtered) #1201/7465

## Metadata ------
load("metadata_T2D.RData")
head(metadata_T2D);dim(metadata_T2D) #4117/27

table(metadata_T2D$study,metadata_T2D$status_new)
table(metadata_T2D$metformin,metadata_T2D$status_new)
t1 <- chisq.test(metadata_T2D$study,metadata_T2D$status_new)
t2 <- chisq.test(metadata_T2D$metformin,metadata_T2D$status_new)
# metafile <- select(dat_filtered,
#                    c("sample_id","age","sex","study","T2D"))
metafile <- metadata_T2D %>% 
  filter(sample_id %in% dat_filtered$sample_id) #%>% 
  # filter(!is.na(bmi))
metafile <- metafile %>% 
  mutate(metf=case_when(metformin==1 ~ "Yes",
                        metformin==0 ~ "No",
                        TRUE ~ "Missing"))

head(metafile);dim(metafile)
colSums(is.na(metafile[,c("age", "sex", "bmi", "study","metf")]))

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
# distfile <- vegdist(pca_df, method="jaccard")
class(distfile)

treefile <- nj(distfile) |>
  ape::ladderize()
head(treefile$edge.length)
summary(treefile$edge.length)
treefile$edge.length[which(treefile$edge.length<0)]
head(treefile$edge)
head(treefile$tip.label)
head(treefile$Nnode)
table(treefile$tip.label %in% metafile$sample_id)

tree_df <- as_tibble(treefile)
tree_df <- tree_df %>% 
  left_join(metafile,by=c("label"="sample_id"))
head(tree_df)

treefile2 <- as.treedata(tree_df)
treefile2


# PGLMM -----------
## base model 
# treefile$edge.length <- treefile$edge.length + 0.00001
re_pglmm <- anpan_pglmm(meta_file = metafile,
                        tree_file = treefile,
                        outcome = "T2D",
                        covariates = c("age", "sex", "bmi", "study"),
                        family = "binomial",
                        out_dir = paste0("pglmm_filtered/",speName,lbl),
                        bug_name = speName,
                        reg_noise = FALSE,
                        loo_comparison = TRUE,
                        omit_na = TRUE,
                        show_plot_tree = FALSE,
                        beta_sd = c(0.1,0.2,0.167,0.2),
                        sigma_phylo_scale = 0.333)

source("~/r_script/t2d_meta/anpan/3.anpan.gene.model/visualization2.R")

p = plot_tree_with_post2(tree_file = treefile,
                        meta_file = metafile,
                        fit = re_pglmm$pglmm_fit,
                        covariates = c("age", "sex", "bmi", "study"),
                        outcome = "T2D",
                        color_bars = TRUE,
                        omit_na = TRUE,
                        labels = re_pglmm$model_input$sample_id
                        # labels = NULL
                        )
pdf(paste0("tree_plot/tree_posterior_",speName,lbl,".pdf"),
    width = 16,height = 9)
p
dev.off()

if (is.null(re_pglmm$error)) {
  print(speName)
  print(re_pglmm$loo$comparison)
}


## base+metf model
re_pglmm_metf <- anpan_pglmm(meta_file = metafile,
                             tree_file = treefile,
                             outcome = "T2D",
                             covariates = c("age", "sex", "bmi", "study", "metf"),
                             family = "binomial",
                             bug_name = speName,
                             out_dir = paste0("pglmm_filtered_metformin/",speName,lbl),
                             reg_noise = FALSE,
                             loo_comparison = TRUE,
                             omit_na = TRUE,
                             run_diagnostics = TRUE,
                             show_plot_tree = FALSE,
                             beta_sd = c(0.1,0.2,0.167,0.2,0.3),
                             sigma_phylo_scale = 0.333)


if (is.null(re_pglmm_metf$error)) {
  print(speName)
  print(re_pglmm_metf$loo$comparison)
}


p_metf = plot_tree_with_post2(tree_file = treefile,
                        meta_file = metafile,
                        fit = re_pglmm_metf$pglmm_fit,
                        covariates = c("age", "sex", "bmi", "study", "metf"),
                        outcome = "T2D",
                        color_bars = TRUE,
                        omit_na = TRUE,
                        labels = re_pglmm_metf$model_input$sample_id
                        # labels = "Sample"
                        )
pdf(paste0("tree_plot/tree_posterior_metformin_",speName,lbl,".pdf"),
    width = 16,height = 9)
p_metf
dev.off()
