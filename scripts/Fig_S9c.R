load("/n/holylabs/LABS/dwang_lab/Lab/T2D_meta/3.analysis/dropbox/fig_s11/fig_s11b_pwy.RData")
gdata::keep(pcopri, sure=T)

# pwy_ranef by T2D within clade
pcopriA_T2D <- pcopri %>% filter(status_new != 9 & clade3cat2 =="Clade A dominant") %>%
  select("species", "pathway", "log10_species_abd", "log10_pwy_abd", "status_new")
colnames(pcopriA_T2D) <- c("bug", "pwy", "log10_species_abd", "log10_pwy_abd", "status_new")
pcopriBCD_T2D <- pcopri %>% filter(status_new != 9 & clade3cat2 =="Other clades") %>%
  select("species", "pathway", "log10_species_abd", "log10_pwy_abd", "status_new")
colnames(pcopriBCD_T2D) <- c("bug", "pwy", "log10_species_abd", "log10_pwy_abd", "status_new")

library(anpan)
pA_T2D <- anpan_pwy_ranef(bug_pwy_dat = pcopriA_T2D,
                          group_ind = "status_new")
pBCD_T2D <- anpan_pwy_ranef(bug_pwy_dat = pcopriBCD_T2D,
                            group_ind = "status_new")

# scatter plot among clade A
aa <- pA_T2D$summary_df[[1]]
aa$lab <- substr(aa$variable,1,11)
bbb <- aa %>% filter(lab=="pwy_effects")

status_color <- c("#5DB1DD", "#C75127")
names(status_color) <- c("Con", "T2D")

pp_pwy_cladeA <- list()
for (i in c(3)) {
  bugdata <- pcopriA_T2D %>% filter(pwy==bbb[i,]$pwy)
  bugdata$status_new <- ifelse(bugdata$status_new==1, "T2D", "Con")
  
  distr <- bugdata %>% summarise(species99=quantile(log10_species_abd, probs=0.99),
                                 pwy99=quantile(log10_pwy_abd, probs=0.99))
  bugdata99 <- bugdata %>% filter(log10_species_abd<=distr$species99 & log10_pwy_abd<=distr$pwy99)
  
  slope <- aa[aa$variable=="species_beta[1]",]$mean
  int_con <- aa[aa$variable=="global_intercept",]$mean +
    aa[aa$variable==paste0("pwy_intercepts[",i,"]"),]$mean
  int_t2d <- int_con + aa[aa$variable==paste0("pwy_effects[",i,"]"),]$mean
  
  pp_pwy_cladeA[[i]] <- ggplot(aes(x=log10_species_abd, y=log10_pwy_abd), data=bugdata99) +
    geom_point(aes(color=status_new), size=1.7) +
    scale_color_manual(values = status_color) +
    geom_abline(slope = slope,
                intercept = int_con,
                color="#5DB1DD",
                size=1.2) +
    geom_abline(slope = slope,
                intercept = int_t2d,
                color="#C75127",
                size=1.2) +
    scale_x_continuous(breaks=c(-1,0,1,2), limits=c(-1.6,2)) +
    scale_y_continuous(breaks=c(-4,-3,-2,-1,0), limits=c(-4.2,0.5)) +
    labs(
         # title="BRANCHED-CHAIN-AA-SYN-PWY:\nsuperpathway of branched\namino acid biosynthesis",
         # title="ILEUSYN-PWY: L-isoleucine\nbiosynthesis I (from threonine)",
         # title="PWY-7111: pyruvate fermentation\nto isobutanol (engineered)",
         # title="VALSYN-PWY: L-valine biosynthesis",
         # x="",
         # y="",
         x=expression(Log[10]*" (species_abd)"),
         y=expression(Log[10]*" (pathway_abd)"),
         color="") +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.text = element_text(size = 16, colour = "black"),
          axis.title = element_text(size = 16, colour = "black"),
          panel.border = element_rect(colour = "black", fill=NA, size=1),
          axis.ticks = element_line(colour = "black", size = 1),
          axis.ticks.length = unit(.15, "cm"),
          legend.position ="none",
          plot.title = element_text(size = 14, hjust = 0.5))
}

# scatter plot among clade BCD
aa <- pBCD_T2D$summary_df[[1]]
aa$lab <- substr(aa$variable,1,11)
bbb <- aa %>% filter(lab=="pwy_effects")

status_color <- c("#5DB1DD", "#C75127")
names(status_color) <- c("Con", "T2D")

pp_pwy_cladeBCD <- list()
for (i in c(3)) {
  bugdata <- pcopriBCD_T2D %>% filter(pwy==bbb[i,]$pwy)
  bugdata$status_new <- ifelse(bugdata$status_new==1, "T2D", "Con")
  
  distr <- bugdata %>% summarise(species99=quantile(log10_species_abd, probs=0.99),
                                 pwy99=quantile(log10_pwy_abd, probs=0.99))
  bugdata99 <- bugdata %>% filter(log10_species_abd<=distr$species99 & log10_pwy_abd<=distr$pwy99)
  
  slope <- aa[aa$variable=="species_beta[1]",]$mean
  int_con <- aa[aa$variable=="global_intercept",]$mean +
    aa[aa$variable==paste0("pwy_intercepts[",i,"]"),]$mean
  int_t2d <- int_con + aa[aa$variable==paste0("pwy_effects[",i,"]"),]$mean
  
  pp_pwy_cladeBCD[[i]] <- ggplot(aes(x=log10_species_abd, y=log10_pwy_abd), data=bugdata99) +
    geom_point(aes(color=status_new), size=1.7) +
    scale_color_manual(values = status_color) +
    geom_abline(slope = slope,
                intercept = int_con,
                color="#5DB1DD",
                size=1.2) +
    geom_abline(slope = slope,
                intercept = int_t2d,
                color="#C75127",
                size=1.2) +
    scale_x_continuous(breaks=c(-2,-1,0,1,2), limits=c(-2.4,2)) +
    scale_y_continuous(breaks=c(-3,-2,-1,0), limits=c(-3.5,0.4)) +
    labs(
         title="",
         # x="",
         # y="",
         x=expression(Log[10]*" (species_abd)"),
         y=expression(Log[10]*" (pathway_abd)"),
         color="") +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.text = element_text(size = 16, colour = "black"),
          axis.title = element_text(size = 16, colour = "black"),
          panel.border = element_rect(colour = "black", fill=NA, size=1),
          axis.ticks = element_line(colour = "black", size = 1),
          axis.ticks.length = unit(.15, "cm"),
          legend.position ="none",
          plot.title = element_text(size = 14, hjust = 0.5))
}

save.image(file="/n/holylabs/LABS/dwang_lab/Lab/T2D_meta/3.analysis/dropbox/fig_s11/fig_s11c_pwy.RData")

gdata::keep(pp_pwy_cladeA, pp_pwy_cladeBCD, sure=T)

#######################################
#               ECs
#######################################
load("/n/holylabs/LABS/dwang_lab/Lab/T2D_meta/3.analysis/dropbox/fig_s11/fig_s11b_ec.RData")
gdata::keep(pcopri, sure=T)

# pwy_ranef by T2D within clade
pcopriA_T2D <- pcopri %>% filter(status_new != 9 & clade3cat2 =="Clade A dominant") %>%
  select("species", "enzy", "log10_species_abd", "log10_ecs_abd", "status_new")
colnames(pcopriA_T2D) <- c("bug", "pwy", "log10_species_abd", "log10_pwy_abd", "status_new")
pcopriBCD_T2D <- pcopri %>% filter(status_new != 9 & clade3cat2 =="Other clades") %>%
  select("species", "enzy", "log10_species_abd", "log10_ecs_abd", "status_new")
colnames(pcopriBCD_T2D) <- c("bug", "pwy", "log10_species_abd", "log10_pwy_abd", "status_new")

pA_T2D <- anpan_pwy_ranef(bug_pwy_dat = pcopriA_T2D,
                          group_ind = "status_new")
pBCD_T2D <- anpan_pwy_ranef(bug_pwy_dat = pcopriBCD_T2D,
                            group_ind = "status_new")

# scatter plot among clade A
status_color <- c("#5DB1DD", "#C75127")
names(status_color) <- c("Con", "T2D")

aa <- pA_T2D$summary_df[[1]]
aa$lab <- substr(aa$variable,1,11)
bbb <- aa %>% filter(lab=="pwy_effects")

bugdata <- pcopriA_T2D %>% filter(pwy==bbb[1,]$pwy)
bugdata$status_new <- ifelse(bugdata$status_new==1, "T2D", "Con")

distr <- bugdata %>% summarise(species99=quantile(log10_species_abd, probs=0.99),
                               pwy99=quantile(log10_pwy_abd, probs=0.99))
bugdata99 <- bugdata %>% filter(log10_species_abd<=distr$species99 & log10_pwy_abd<=distr$pwy99)

slope <- aa[aa$variable=="species_beta[1]",]$mean
int_con <- aa[aa$variable=="global_intercept",]$mean +
  aa[aa$variable==paste0("pwy_intercepts[",1,"]"),]$mean
int_t2d <- int_con + aa[aa$variable==paste0("pwy_effects[",1,"]"),]$mean

pp_ecs_cladeA <- ggplot(aes(x=log10_species_abd, y=log10_pwy_abd), data=bugdata99) +
  geom_point(aes(color=status_new), size=1.7) +
  scale_color_manual(values = status_color) +
  geom_abline(slope = slope,
              intercept = int_con,
              color="#5DB1DD",
              size=1.2) +
  geom_abline(slope = slope,
              intercept = int_t2d,
              color="#C75127",
              size=1.2) +
  scale_x_continuous(breaks=c(-2,-1,0,1,2), limits=c(-2,2)) +
  scale_y_continuous(breaks=c(-6,-5,-4,-3,-2,-1), limits=c(-6,-0.5)) +
  labs(title="EC2.6.1.42: Branched-chain-\namino-acid transaminase",
       x="",
       y="",
       color="") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 16, colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        axis.ticks = element_line(colour = "black", size = 1),
        axis.ticks.length = unit(.15, "cm"),
        legend.position ="none",
        plot.title = element_text(size = 14, hjust = 0.5))

# scatter plot among clade B
status_color <- c("#5DB1DD", "#C75127")
names(status_color) <- c("Con", "T2D")

aa <- pBCD_T2D$summary_df[[1]]
aa$lab <- substr(aa$variable,1,11)
bbb <- aa %>% filter(lab=="pwy_effects")

bugdata <- pcopriBCD_T2D %>% filter(pwy==bbb[1,]$pwy)
bugdata$status_new <- ifelse(bugdata$status_new==1, "T2D", "Con")

distr <- bugdata %>% summarise(species99=quantile(log10_species_abd, probs=0.99),
                               pwy99=quantile(log10_pwy_abd, probs=0.99))
bugdata99 <- bugdata %>% filter(log10_species_abd<=distr$species99 & log10_pwy_abd<=distr$pwy99)

slope <- aa[aa$variable=="species_beta[1]",]$mean
int_con <- aa[aa$variable=="global_intercept",]$mean +
  aa[aa$variable==paste0("pwy_intercepts[",1,"]"),]$mean
int_t2d <- int_con + aa[aa$variable==paste0("pwy_effects[",1,"]"),]$mean

pp_ecs_cladeBCD <- ggplot(aes(x=log10_species_abd, y=log10_pwy_abd), data=bugdata99) +
  geom_point(aes(color=status_new), size=1.7) +
  scale_color_manual(values = status_color) +
  geom_abline(slope = slope,
              intercept = int_con,
              color="#5DB1DD",
              size=1.2) +
  geom_abline(slope = slope,
              intercept = int_t2d,
              color="#C75127",
              size=1.2) +
  scale_x_continuous(breaks=c(-2,-1,0,1,2), limits=c(-2.4,2)) +
  scale_y_continuous(breaks=c(-5,-4,-3,-2,-1), limits=c(-5.2,-0.5)) +
  labs(title="",
       x="",
       y="",
       color="") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 16, colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        axis.ticks = element_line(colour = "black", size = 1),
        axis.ticks.length = unit(.15, "cm"),
        legend.position ="none",
        plot.title = element_text(size = 14, hjust = 0.5))

save.image(file="/n/holylabs/LABS/dwang_lab/Lab/T2D_meta/3.analysis/dropbox/fig_s11/fig_s11c_ec.RData")

gdata::keep(pp_pwy_cladeA, pp_pwy_cladeBCD,
            pp_ecs_cladeA, pp_ecs_cladeBCD, sure=T)

blank <- ggplot()+theme_void()
pdf("/n/holylabs/LABS/dwang_lab/Lab/T2D_meta/3.analysis/dropbox/fig_s11/fig_s11c.pdf",
    width = 16, height = 8.1, onefile = F) # Open a new pdf file
egg::ggarrange(pp_pwy_cladeA[[1]], pp_pwy_cladeA[[2]], pp_pwy_cladeA[[4]], pp_ecs_cladeA,
               blank, blank, blank, blank,
               pp_pwy_cladeBCD[[1]], pp_pwy_cladeBCD[[2]], pp_pwy_cladeBCD[[4]], pp_ecs_cladeBCD,
               nrow = 3, heights = c(1,-0.15,1))
dev.off() # Close the file
