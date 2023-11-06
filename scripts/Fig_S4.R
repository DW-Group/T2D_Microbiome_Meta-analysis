load("/n/holylabs/LABS/dwang_lab/Lab/T2D_meta/3.analysis/dropbox/fig2/figure2.RData")

gdata::keep(mds.data2, mds.var.per2, species_zero_adj0, spe_all, sure=TRUE)
identical(colnames(species_zero_adj0), mds.data2$id) # TRUE

library(tidyverse)
spe_all <- spe_all %>% arrange(desc(relab))
spe_9 <- spe_all[1:9,]$species
spe_lab <- gsub("s__", "", spe_9)

for (i in 1:9) {
  print(quantile(t(species_zero_adj0[rownames(species_zero_adj0)==spe_9[i],])*100,
                 probs=c(0.9,0.95,0.98)))
}

# 90%      95%      98% 
# 42.67521 58.92360 72.01703 
# 90%      95%      98% 
# 19.82257 26.19556 36.01521 
# 90%      95%      98% 
# 18.11825 23.46260 30.28440 
# 90%      95%      98% 
# 14.41105 17.76241 21.95716 
# 90%      95%      98% 
# 13.56730 18.18783 24.22623 
# 90%       95%       98% 
# 9.954138 13.033533 17.852758 
# 90%       95%       98% 
# 9.676695 16.437705 26.377816 
# 90%      95%      98% 
# 8.66281 14.27247 21.88481 
# 90%       95%       98% 
# 7.385362 10.510069 15.309381

library(ggplot2)
pcoaplot <- list()
for (i in c(1,2)) {
  identical(colnames(species_zero_adj0), mds.data2$id) # TRUE
  mds.data2$species <- t(species_zero_adj0[rownames(species_zero_adj0)==spe_9[i],])*100
  mds.data2$species[mds.data2$species>60] <- 60
  pcoaplot[[i]] <- ggplot(data=mds.data2, aes(x=X, y=Y)) +
    geom_point(size=3.6, aes(fill=species), shape=21) +
    scale_fill_viridis_c(direction=-1, limits=c(0,60), breaks=c(0,20,40,60)) +
    xlab(paste("PCo1 - ", mds.var.per2[1], "%", sep="")) +
    ylab(paste("PCo2 - ", mds.var.per2[2], "%", sep="")) +
    labs(title=spe_lab[i]) +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.x = element_text(size = 14, colour = "black"),
          axis.text.y = element_text(size = 14, colour = "black"),
          axis.title.x = element_text(size = 14, colour = "black"),
          axis.title.y = element_text(size = 14, colour = "black"),
          panel.border = element_rect(colour = "black", fill=NA, size=1),
          plot.title = element_text(size = 14, face = "italic", hjust = 0.5),
          legend.position="none")
}

i=3
mds.data2$species <- t(species_zero_adj0[rownames(species_zero_adj0)==spe_9[i],])*100
mds.data2$species[mds.data2$species>60] <- 60
pcoaplot[[i]] <- ggplot(data=mds.data2, aes(x=X, y=Y)) +
  geom_point(size=3.6, aes(fill=species), shape=21) +
  scale_fill_viridis_c(direction=-1, limits=c(0,60), breaks=c(0,20,40,60)) +
  xlab(paste("PCo1 - ", mds.var.per2[1], "%", sep="")) +
  ylab(paste("PCo2 - ", mds.var.per2[2], "%", sep="")) +
  labs(title=spe_lab[i],
       fill="Relative \n abundance \n (%)") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 14, colour = "black"),
        axis.text.y = element_text(size = 14, colour = "black"),
        axis.title.x = element_text(size = 14, colour = "black"),
        axis.title.y = element_text(size = 14, colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        plot.title = element_text(size = 14, face = "italic", hjust = 0.5),
        legend.title.align=0.5,
        legend.direction="vertical")

for (i in c(4,5,7,8)) {
  identical(colnames(species_zero_adj0), mds.data2$id) # TRUE
  mds.data2$species <- t(species_zero_adj0[rownames(species_zero_adj0)==spe_9[i],])*100
  mds.data2$species[mds.data2$species>20] <- 20
  pcoaplot[[i]] <- ggplot(data=mds.data2, aes(x=X, y=Y)) +
    geom_point(size=3.6, aes(fill=species), shape=21) +
    scale_fill_viridis_c(direction=-1, limits=c(0,20), breaks=c(0,5,10,15,20)) +
    xlab(paste("PCo1 - ", mds.var.per2[1], "%", sep="")) +
    ylab(paste("PCo2 - ", mds.var.per2[2], "%", sep="")) +
    labs(title=spe_lab[i]) +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.x = element_text(size = 14, colour = "black"),
          axis.text.y = element_text(size = 14, colour = "black"),
          axis.title.x = element_text(size = 14, colour = "black"),
          axis.title.y = element_text(size = 14, colour = "black"),
          panel.border = element_rect(colour = "black", fill=NA, size=1),
          plot.title = element_text(size = 14, face = "italic", hjust = 0.5),
          legend.position="none")
}

for (i in c(6,9)) {
  mds.data2$species <- t(species_zero_adj0[rownames(species_zero_adj0)==spe_9[i],])*100
  mds.data2$species[mds.data2$species>20] <- 20
  pcoaplot[[i]] <- ggplot(data=mds.data2, aes(x=X, y=Y)) +
    geom_point(size=3.6, aes(fill=species), shape=21) +
    scale_fill_viridis_c(direction=-1, limits=c(0,20), breaks=c(0,5,10,15,20)) +
    xlab(paste("PCo1 - ", mds.var.per2[1], "%", sep="")) +
    ylab(paste("PCo2 - ", mds.var.per2[2], "%", sep="")) +
    labs(title=spe_lab[i],
         fill="Relative \n abundance \n (%)") +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.x = element_text(size = 14, colour = "black"),
          axis.text.y = element_text(size = 14, colour = "black"),
          axis.title.x = element_text(size = 14, colour = "black"),
          axis.title.y = element_text(size = 14, colour = "black"),
          panel.border = element_rect(colour = "black", fill=NA, size=1),
          plot.title = element_text(size = 14, face = "italic", hjust = 0.5),
          legend.title.align=0.5,
          legend.direction="vertical")
}

library(patchwork)
pdf("/n/holylabs/LABS/dwang_lab/Lab/T2D_meta/3.analysis/dropbox/fig_s2/figure_9pcoa.pdf",
    width = 13, height = 12, onefile = F) # Open a new pdf file
(pcoaplot[[1]] | pcoaplot[[2]] | pcoaplot[[3]])/(pcoaplot[[4]] | pcoaplot[[5]] | pcoaplot[[6]])/(pcoaplot[[7]] | pcoaplot[[8]] | pcoaplot[[9]])
dev.off() # Close the file
