setwd("/n/holylabs/LABS/dwang_lab/Lab/T2D_meta/3.analysis/dropbox/fig_s11/")
library(tidyverse)

#############################
#         pathway 
#############################
load("/n/holylabs/LABS/dwang_lab/Lab/T2D_meta/3.analysis/0323/pwy_ranef/pwy_ranef0.RData")
load("/n/holylabs/LABS/dwang_lab/Lab/T2D_meta/3.analysis/0323/metaphlan4/metadata_wthclade.RData")

pcopri <- bug_pwy_dat %>% filter(species=='s__Prevotella_copri' &
                                   pathway %in% c("BRANCHED-CHAIN-AA-SYN-PWY: superpathway of branched amino acid biosynthesis",
                                                  "ILEUSYN-PWY: L-isoleucine biosynthesis I (from threonine)",
                                                  "VALSYN-PWY: L-valine biosynthesis",
                                                  "PWY-7111: pyruvate fermentation to isobutanol (engineered)"))

# clade3cat2: clade A >75%
pcopri <- merge(pcopri, metadata_wthclade[,c("id","clade3cat2")], by="id", all.x=T)

pcopri <- pcopri %>% 
  mutate(status_new=case_when(status_new=="T2D" ~ 1,
                              status_new=="Con" ~ 0,
                              status_new=="Pre" ~ 9)) %>% filter(clade3cat2 != "No Pcopri")

# pwy_ranef by clade, among all participants
pcopriA <- pcopri %>% select("species", "pathway", "log10_species_abd", "log10_pwy_abd", "clade3cat2")
colnames(pcopriA) <- c("bug", "pwy", "log10_species_abd", "log10_pwy_abd", "clade")

library(anpan)
pcopriA$clade <- ifelse(pcopriA$clade=="Clade A dominant", 1, 0)
pA <- anpan_pwy_ranef(bug_pwy_dat = pcopriA,
                      group_ind = "clade")

aa <- pA$summary_df[[1]]

gdata::keep(pcopri, pcopriA, pA, aa, sure=T)
save.image("/n/holylabs/LABS/dwang_lab/Lab/T2D_meta/3.analysis/dropbox/fig_s11/fig_s11b_pwy.RData")

#############################
#         enzyme 
#############################
load("/n/holylabs/LABS/dwang_lab/Lab/T2D_meta/3.analysis/0323/pwy_ranef/ecs_ranef0.RData")
load("/n/holylabs/LABS/dwang_lab/Lab/T2D_meta/3.analysis/0323/metaphlan4/metadata_wthclade.RData")

pcopri <- bug_ecs_dat %>% filter(species=='s__Prevotella_copri' &
                                   enzy=='EC2.6.1.42: Branched-chain-amino-acid transaminase')
pcopri <- merge(pcopri, metadata_wthclade[,c("id","clade3cat2")], by="id", all.x=T)

pcopri <- pcopri %>% 
  mutate(status_new=case_when(status_new=="T2D" ~ 1,
                              status_new=="Con" ~ 0,
                              status_new=="Pre" ~ 9)) %>% filter(clade3cat2 != "No Pcopri")

gdata::keep(pcopri, sure=T)

# pwy_ranef by clade, among all participants
pcopriA <- pcopri %>% select("species", "enzy", "log10_species_abd", "log10_ecs_abd", "clade3cat2")
colnames(pcopriA) <- c("bug", "pwy", "log10_species_abd", "log10_pwy_abd", "clade")

pcopriA$clade <- ifelse(pcopriA$clade=="Clade A dominant", 1, 0)
pA <- anpan_pwy_ranef(bug_pwy_dat = pcopriA,
                      group_ind = "clade")

aa <- pA$summary_df[[1]]

gdata::keep(pcopri, pcopriA, pA, aa, sure=T)
save.image("/n/holylabs/LABS/dwang_lab/Lab/T2D_meta/3.analysis/dropbox/fig_s11/fig_s11b_ec.RData")

####################################
#     plot by clades for pwy
####################################
load("/n/holylabs/LABS/dwang_lab/Lab/T2D_meta/3.analysis/dropbox/fig_s11/fig_s11b_pwy.RData")

library(ggplot2)
pwy <- unique(pcopriA$pwy)
pwy <- pwy[c(1,2,4,3)]
pcopriA$clade <- ifelse(pcopriA$clade==1, "Clade A dominant", "Other clades")

colclade <- c("#C75127","#5DB1DD")
names(colclade) <- c("Clade A dominant", "Other clades")

boxp_pwy <- list()
for (i in c(4)){
  temp <- pcopriA[pcopriA$pwy==pwy[i], ]
  mod <- lm(log10_pwy_abd ~ log10_species_abd, data=temp)
  temp$pwy_abd <- 10^(resid(mod)+mean(temp$log10_pwy_abd))
  temp <- temp %>% group_by(clade) %>% mutate(mean=mean(pwy_abd),
                                              sd=sd(pwy_abd)) %>% filter(pwy_abd<=(mean+3*sd)) 
  
  boxp_pwy[[i]] <- ggplot(aes(x=clade, y=pwy_abd, fill=clade), data=temp) +
    geom_jitter(aes(color=clade), width=0.1, alpha=0.7) +
    geom_boxplot(width=0.3, outlier.color = NA) +
    scale_fill_manual(values=colclade) +
    scale_color_manual(values=colclade) +
    scale_y_continuous(breaks = c(0,0.5,1,1.5,2), limits=c(0,2.15)) +
    labs(
         # title="BRANCHED-CHAIN-AA-SYN-PWY:\nsuperpathway of branched\namino acid biosynthesis",
         # title="ILEUSYN-PWY: L-isoleucine\nbiosynthesis I (from threonine)",
         # title="VALSYN-PWY: L-valine biosynthesis",
         # title="PWY-7111: pyruvate fermentation\nto isobutanol (engineered)",
         x="",
         y="pathway_abd (%)",
         color="") +
    annotate(geom="text", x=0.8, y=2, size=5, col="black", label="P<0.001") + 
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.text = element_text(size = 14, colour = "black"),
          axis.title = element_text(size = 14, colour = "black"),
          # axis.title.y = element_blank(),
          # axis.text.y = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA, size=1),
          axis.ticks = element_line(colour = "black", size = 1),
          axis.ticks.length = unit(.15, "cm"),
          legend.position = "none",
          plot.title = element_text(size = 12, hjust = 0.5))
}

gdata::keep(boxp_pwy, colclade, sure=T)

# plot by clades for ec
load("/n/holylabs/LABS/dwang_lab/Lab/T2D_meta/3.analysis/dropbox/fig_s11/fig_s11b_ec.RData")

pcopriA$clade <- ifelse(pcopriA$clade==1, "Clade A dominant", "Other clades")

mod <- lm(log10_pwy_abd ~ log10_species_abd, data=pcopriA)
pcopriA$pwy_abd <- exp(resid(mod)+mean(pcopriA$log10_pwy_abd))

pcopriA <- pcopriA %>% group_by(clade) %>% mutate(pp=quantile(pwy_abd, 0.99))

boxp_ec <- ggplot(aes(x=clade, y=pwy_abd, fill=clade), data=pcopriA[pcopriA$pwy_abd<pcopriA$pp,]) +
  geom_jitter(aes(color=clade), width=0.1, alpha=0.7) +
  geom_boxplot(width=0.3, outlier.color = NA) +
  scale_y_continuous(breaks = c(0,0.3,0.6,0.9,1.2), limits=c(0,1.3)) +
  scale_fill_manual(values=colclade) +
  scale_color_manual(values=colclade) +
  labs(title="EC2.6.1.42: Branched-chain-\namino-acid transaminase",
       x="",
       y="enzyme_abd (%)",
       color="") +
  annotate(geom="text", x=0.8, y=1.2, size=5, col="black", label="P<0.001") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 14, colour = "black"),
        # axis.text.y = element_blank(),
        axis.title = element_text(size = 14, colour = "black"),
        # axis.title.y = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        axis.ticks = element_line(colour = "black", size = 1),
        axis.ticks.length = unit(.15, "cm"),
        legend.position = "none",
        plot.title = element_text(size = 12, hjust = 0.5))

pdf("/n/holylabs/LABS/dwang_lab/Lab/T2D_meta/3.analysis/dropbox/fig_s11/fig_s11b.pdf",
    width = 15, height = 4.2, onefile = F) # Open a new pdf file
egg::ggarrange(boxp_pwy[[1]], boxp_pwy[[2]], boxp_pwy[[3]], boxp_ec, nrow = 1)
dev.off() # Close the file

save.image(file="/n/holylabs/LABS/dwang_lab/Lab/T2D_meta/3.analysis/dropbox/fig_s11/fig_s11b.RData")

