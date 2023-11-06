# dbRDA plot for metf/BMI, 8 studies

load("/n/holylabs/LABS/dwang_lab/Lab/T2D_meta/3.analysis/0323/species_T2D_main.RData")
gdata::keep(metadata_T2D, species_T2D_zero, sure=T)
table(metadata_T2D$status_new, metadata_T2D$metf)
#     Missing   No  Yes
# Con       0 3134    6
# T2D     102 1255  617

library(tidyverse)
metadata_T2D_new <- metadata_T2D %>%
  filter(!is.na(bmi)) %>%
  filter((status_new=="T2D" & metf != "Missing") | (status_new=="Con" & metf == "No"))
table(metadata_T2D_new$status_new, metadata_T2D_new$metf)
#       No  Yes
# Con 3127    0
# T2D 1252  615

species_T2D_zero_new <- species_T2D_zero[, metadata_T2D_new$id]
metadata_T2D_new <- metadata_T2D_new %>%
  mutate(T2Dmetf=case_when(status_new=="Con" ~ "Con",
                           status_new=="T2D" & metf=="No" ~ "T2D, metf-",
                           status_new=="T2D" & metf=="Yes" ~ "T2D, metf+"))

######### dbRDA by metformin #########
library(vegan)
library(ggvegan)
distance.matrix <- vegdist(t(species_T2D_zero_new), method="bray")

dbRDA <- capscale(distance.matrix ~ age+sex+bmi+metf,
                  comm=t(species_T2D_zero_new),
                  data=metadata_T2D_new)
anova(dbRDA, by="terms")
# Permutation test for capscale under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 999
# 
# Model: capscale(formula = distance.matrix ~ age + sex + bmi + metf, 
# data = metadata_T2D_new, comm = t(species_T2D_zero_new))
# Df SumOfSqs       F Pr(>F)    
# age         1     5.14 12.6173  0.001 ***
# sex         1     6.65 16.3312  0.001 ***
# bmi         1     5.76 14.1426  0.001 ***
# metf        1     2.83  6.9531  0.001 ***
# Residual 4989  2031.19  

ss <- summary(dbRDA)
ss_species <- data.frame(ss$species)
ss_sites <- data.frame(ss$sites)
identical(rownames(ss_sites), metadata_T2D_new$id) # TRUE
ss_sites$group <- factor(metadata_T2D_new$T2Dmetf,
                         levels=c("Con", "T2D, metf-", "T2D, metf+"))

dotcol <- c("#5DB1DD", "#CC9900", "#C75127")
names(dotcol) <- c("Con", "T2D, metf-", "T2D, metf+")

plot_metf <- ggplot(aes(x=CAP1, y=-CAP2), data=ss_sites) +
  geom_point(aes(fill=group), shape=21, size=4) +
  scale_fill_manual(values=dotcol) +
  stat_ellipse(type = "t", aes(color=group), lwd=1.2) +
  scale_color_manual(values=dotcol) +
  scale_x_continuous(breaks = c(-4,0,4,8), limits = c(-6,9)) +
  scale_y_continuous(breaks = c(-8,-4,0,4,8), limits = c(-11,8)) +
  xlab("dbRDA1") +
  ylab("dbRDA2") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(size=2)) +
  theme(axis.title = element_blank(),
        axis.text = element_text(size=20),
        legend.title = element_blank(),
        legend.text = element_text(size=13),
        legend.position = c(0.80, 0.13))

x_box <- ggplot(aes(y=group, x=-CAP1), data=ss_sites) +
  geom_boxplot(aes(fill=group), width=.4, outlier.shape = NA) +
  scale_fill_manual(values=dotcol) +
  scale_x_continuous(breaks = c(-4,0,4,8), limits = c(-6,9)) +
  xlab("dbRDA1") +
  ylab("") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_blank()) +
  theme(axis.title = element_text(size=20),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none")

y_box <- ggplot(aes(x=group, y=-CAP2), data=ss_sites) +
  geom_boxplot(aes(fill=group), width=.4, outlier.shape = NA) +
  scale_fill_manual(values=dotcol) +
  scale_y_continuous(breaks = c(-8,-4,0,4,8), limits = c(-11,8)) +
  xlab("") +
  ylab("dbRDA2") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_blank()) +
  theme(axis.title = element_text(size=20),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none")

blank <- ggplot() +
  theme(panel.background = element_rect(fill = 'white', colour = 'white'))

pdf("/n/holylabs/LABS/dwang_lab/Lab/T2D_meta/3.analysis/dropbox/fig_s6/figure_s6b.pdf",
    width = 5.5, height = 5.5, onefile = F)
egg::ggarrange(y_box, plot_metf,
               blank, x_box, nrow=2, ncol=2,
               widths = c(1,9),
               heights = c(9,1))
dev.off()

save.image(file="/n/holylabs/LABS/dwang_lab/Lab/T2D_meta/3.analysis/dropbox/fig_s6/figure_s6b.RData")
