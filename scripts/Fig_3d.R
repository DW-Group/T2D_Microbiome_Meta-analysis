setwd("/n/holylabs/LABS/dwang_lab/Lab/T2D_meta/3.analysis/R1/fig3")
load("/n/holylabs/LABS/dwang_lab/Lab/T2D_meta/3.analysis/0323/pwy_ranef/pwy_ranef2_22.RData")

library(tidyverse)
# plot for pcopri
aa <- fit_22[[1]]$summary_df[[1]]
aa$lab <- substr(aa$variable,1,11)
bbb <- aa %>% filter(lab=="pwy_effects") %>% arrange(desc(mean))
pcopri <- bug_pwy_dat2 %>% filter(bug==spe_22[1] & pwy %in% bbb[1:4,]$pwy)
pcopri$pwy <- factor(pcopri$pwy, levels=bbb[1:4,]$pwy)

bbb$col <- ifelse(bbb$mean<0.2, "Non-hit", "Hit")
ptcol <- c("red", "black")
names(ptcol) <- c("Hit", "Non-hit")
bbb$pwy <- factor(bbb$pwy, levels=rev(bbb$pwy))

library(ggplot2)
# scatter plot
pp_pwy <- list()
status_color <- c("#5DB1DD", "#C75127")
names(status_color) <- c("Con", "T2D")

for (i in c(1)) {
  bugdata <- bug_pwy_dat2 %>% filter(bug==spe_22[1] & pwy %in% bbb[i,]$pwy)
  bugdata$status_new <- ifelse(bugdata$status_new==1, "T2D", "Con")
  
  distr <- bugdata %>% summarise(species99=quantile(log10_species_abd, probs=0.99),
                                 pwy99=quantile(log10_pwy_abd, probs=0.99))
  bugdata99 <- bugdata %>% filter(log10_species_abd<=distr$species99 & log10_pwy_abd<=distr$pwy99)
  
  pp_pwy[[i]] <- ggplot(aes(x=10^(log10_species_abd), y=10^(log10_pwy_abd), color=status_new), data=bugdata99) +
    geom_point(size=1.8) +
    scale_color_manual(values = status_color) +
    geom_smooth(method = lm, formula=y~x+0, se = FALSE, size=1.5) +
    scale_x_continuous(breaks=c(0,20,40,60,80), limits=c(0,85)) +
    scale_y_continuous(breaks=c(0,0.5,1.0,1.5,2.0,2.5), limits=c(0,2.6)) + # 2.2
    labs(
         title="BRANCHED-CHAIN-AA-SYN-PWY:\nsuperpathway of branched\namino acid biosynthesis",
         # title="ILEUSYN-PWY: L-isoleucine\nbiosynthesis I (from threonine)",
         # title="VALSYN-PWY:\nL-valine biosynthesis",
         x="",
         # x="P.copri relative abundance (%)",
         # title="PWY-7111: pyruvate fermentation\nto isobutanol (engineered)", x="species_abd (%)",
         y="Pathway relative abundance (%)",
         color="") +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.text = element_text(size = 14, colour = "black"),
          axis.title = element_text(size = 14, colour = "black"),
          # axis.text.y = element_blank(),
          # axis.title.y = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA, size=1),
          axis.ticks = element_line(colour = "black", size = 1),
          axis.ticks.length = unit(.15, "cm"),
          legend.position ="none",
          # legend.position =c(0.15,0.9),
          # legend.title = element_blank(),
          # legend.text = element_text(size = 14, color="black"),
          plot.title = element_text(size = 12, hjust = 0.5, vjust=1))
}

save.image("/n/holylabs/LABS/dwang_lab/Lab/T2D_meta/3.analysis/R1/fig3/figure3d_pwy.RData")

gdata::keep(pp_pwy, status_color, sure=T)

# for BCAA enzyme
load("/n/holylabs/LABS/dwang_lab/Lab/T2D_meta/3.analysis/0323/pwy_ranef/ecs_ranef2_pcopri.RData")

aa <- pcopri$summary_df[[1]]
aa$lab <- substr(aa$variable,1,11)
bb <- aa %>% filter(lab=="pwy_effects") %>% arrange(desc(mean))

# EC2.6.1.42: Branched-chain-amino-acid transaminase
bugdata <- bug_ecs_dat2 %>% filter(bug==spe_22[1] & pwy %in% bb[5,]$pwy)
bugdata$status_new <- ifelse(bugdata$status_new==1, "T2D", "Con")

distr <- bugdata %>% summarise(species99=quantile(log10_species_abd, probs=0.99),
                               pwy99=quantile(log10_pwy_abd, probs=0.99))
bugdata99 <- bugdata %>% filter(log10_species_abd<=distr$species99 & log10_pwy_abd<=distr$pwy99)

pp_ec <- ggplot(aes(x=10^(log10_species_abd), y=10^(log10_pwy_abd), color=status_new), data=bugdata) +
  geom_point(size=1.8) +
  scale_color_manual(values = status_color) +
  scale_x_continuous(breaks=c(0,20,40,60,80), limits=c(0,91)) +
  scale_y_continuous(breaks=c(0,0.1,0.2,0.3,0.4,0.5), limits=c(0,0.55)) +
  geom_smooth(method = lm, formula=y~x+0, se = FALSE, size=1.5) +
  labs(title=stringr::str_wrap(bb[5,]$pwy, width = 27),
       x="",
       y="Enzyme relative abundance (%)",
       color="") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 14, colour = "black"),
        # axis.text.x = element_blank(),
        axis.title = element_text(size = 14, colour = "black"),
        # axis.title.x = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        axis.ticks = element_line(colour = "black", size = 1),
        axis.ticks.length = unit(.15, "cm"),
        legend.position = "none",
        legend.title = element_blank(),
        legend.text = element_text(size = 14, color="black"),
        plot.title = element_text(size = 12, hjust = 0.5))

pdf("/n/holylabs/LABS/dwang_lab/Lab/T2D_meta/3.analysis/R1/fig3/figure3d.pdf",
    width = 13, height = 4, onefile = F) # Open a new pdf file
egg::ggarrange(pp_pwy[[1]], pp_pwy[[3]], 
               pp_pwy[[4]], pp_ec, nrow = 1)
dev.off() # Close the file

