setwd("/n/holylabs/LABS/dwang_lab/Lab/T2D_meta/3.analysis/dropbox/fig_s11/")
library(tidyverse)

load("/n/holylabs/LABS/dwang_lab/Lab/T2D_meta/3.analysis/0323/metaphlan4/metaphlan4v2.RData")

pcopri2 <- merge(metadata[,c("id","study","status_new")], pcopri)
pcopri2 <- gather(pcopri2, clade, relab, c(clade_A,clade_B,clade_C,clade_D))
pcopri2 <- pcopri2 %>% mutate(clade=case_when(clade=="clade_A" ~ "Clade A",
                                              clade=="clade_B" ~ "Clade B",
                                              clade=="clade_C" ~ "Clade C",
                                              clade=="clade_D" ~ "Clade D"))
pcopri2$abs <- ifelse(pcopri2$relab>0, "Present", "Absent")
pcopri2$clade <- factor(pcopri2$clade, levels=c("Clade D", "Clade C", "Clade B", "Clade A"))

pcopri$study <- as.character(pcopri$study)
pcopri$study <- ifelse(pcopri$study=="Fromentin_2022 (DEU_DNK_FRA)", 
                       "Fromentin_2022 (DEU/DNK/FRA)", pcopri$study)
pcopri$study <- factor(pcopri$study,
                        levels=c("DIRECT-PLUS (ISR)", "Qin_2012 (CHN)", "Zhong_2019 (CHN)", "SOL (USA)", 
                                 "Fromentin_2022 (DEU/DNK/FRA)", "HPFS (USA)", "Karlsson_2013 (SWE)", 
                                 "NHSII (USA)", "Wu_2020a (SWE)", "Wu_2020b (SWE)"))
pcopri$clade_A <- ifelse(pcopri$clade_A>0,1,0)
pcopri$clade_B <- ifelse(pcopri$clade_B>0,1,0)
pcopri$clade_C <- ifelse(pcopri$clade_C>0,1,0)
pcopri$clade_D <- ifelse(pcopri$clade_D>0,1,0)
pcopri <- pcopri %>% arrange(study, desc(clade_A), desc(clade_B), desc(clade_C), desc(clade_D))
pcopri2$id <- factor(pcopri2$id, levels=pcopri$id)
pcopri2$study <- factor(pcopri2$study,
                       levels=c("DIRECT-PLUS (ISR)", "Qin_2012 (CHN)", "Zhong_2019 (CHN)", "SOL (USA)", 
                                "Fromentin_2022 (DEU/DNK/FRA)", "HPFS (USA)", "Karlsson_2013 (SWE)", 
                                "NHSII (USA)", "Wu_2020a (SWE)", "Wu_2020b (SWE)"))

colr <- c("#fb6a4a", "#f0f0f0")
names(colr) <- c("Present", "Absent")

library(ggplot2)
p1 <- ggplot(aes(x=id, y=clade), data=pcopri2) +
  geom_tile(aes(fill=abs)) +
  scale_fill_manual(values=colr) +
  theme_bw() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(color = "black", size = 14),
        axis.title =element_blank(),
        legend.position = "none")

studycol <- c("#466983", "#339900", "#837B8D", 
              "#809900", "#5DB1DD", "#802268",
              "#CC9900", "#CE3D32","#C75127", "#996600")
names(studycol) <- c("DIRECT-PLUS (ISR)", "Fromentin_2022 (DEU/DNK/FRA)", "HPFS (USA)", 
                     "Karlsson_2013 (SWE)", "NHSII (USA)", "Qin_2012 (CHN)",
                     "SOL (USA)", "Wu_2020a (SWE)", "Wu_2020b (SWE)", "Zhong_2019 (CHN)")
studycol <- studycol[c(1,6,10,7,2:5,8,9)]

p2 <- ggplot(aes(x=id, y=clade), data=pcopri2) +
  geom_tile(aes(fill=study)) +
  scale_fill_manual(values=studycol) +  xlab("") +
  theme_bw() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title =element_blank(),
        legend.position = "none")

blank <- ggplot()+theme_void()

pdf("/n/holylabs/LABS/dwang_lab/Lab/T2D_meta/3.analysis/dropbox/fig_s11/fig_s11a.pdf", width = 15, height = 2.5, onefile = F) # Open a new pdf file
egg::ggarrange(p2, blank, p1,
               ncol = 1, heights = c(0.1,-0.1,1))
dev.off()

pdf("/n/holylabs/LABS/dwang_lab/Lab/T2D_meta/3.analysis/dropbox/fig_s11/fig_s11a_legend1.pdf", width = 5, height = 10, onefile = F) # Open a new pdf file
ggplot(aes(x=id, y=clade), data=pcopri2) +
  geom_tile(aes(fill=abs)) +
  scale_fill_manual(values=colr) +
  theme_bw() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title =element_blank(),
        legend.position = "right",
        legend.title = element_blank(),
        legend.text = element_text(color = "black", size = 14))
dev.off()

pdf("/n/holylabs/LABS/dwang_lab/Lab/T2D_meta/3.analysis/dropbox/fig_s11/fig_s11a_legend2.pdf", width = 5, height = 10, onefile = F) # Open a new pdf file
ggplot(aes(y=id, x=clade), data=pcopri2) +
  geom_tile(aes(fill=study)) +
  scale_fill_manual(values=studycol) +  xlab("") +
  theme_bw() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title =element_blank(),
        legend.position = "right",
        legend.title = element_blank(),
        legend.text = element_text(color = "black", size = 14))
dev.off()

save.image("/n/holylabs/LABS/dwang_lab/Lab/T2D_meta/3.analysis/dropbox/fig_s11/fig_s11a.RData")
