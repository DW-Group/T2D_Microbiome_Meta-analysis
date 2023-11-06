setwd("/n/holylabs/LABS/dwang_lab/Lab/T2D_meta/3.analysis/0323/RF8/")
library(tidyverse)

datapath <- dir(path="/n/holylabs/LABS/dwang_lab/Lab/T2D_meta/3.analysis/0323/RF8/", pattern = "_estimations.txt")

datapathNo <- datapath[1:200]
groupNo <- gsub("_estimations.txt","",datapathNo)
groupnameNo <- cbind(path=datapathNo,
                     group=groupNo,
                     data.frame(str_split(groupNo, "_", simplify = T)))
RFdataNo <- list()
for (j in 1:length(groupNo)) {
  temp1 <- read.delim(paste0("/n/holylabs/LABS/dwang_lab/Lab/T2D_meta/3.analysis/0323/RF8/",
                             datapathNo[j]),
                      header=FALSE)
  temp2 <- temp1[2:1001,]
  
  temp3 <- data.frame()
  for (i in 1:20) {
    temp11 <- temp2[(50*(i-1)+1):(50*i), ]
    temp22 <- data.frame(index=as.numeric(t(temp11[temp11$V1=="sample index", -1])),
                         true=as.numeric(t(temp11[temp11$V1=="true labels", -1])),
                         est=as.numeric(t(temp11[temp11$V1=="estimated labels", -1])),
                         pred=as.numeric(t(temp11[temp11$V1=="estimated probabilities", -1])),
                         rep=i)
    temp3 <- rbind(temp3, temp22)
  }
  
  temp3 <- temp3 %>% filter(!is.na(index)) %>%
    group_by(index) %>% summarize(true=mean(true),
                                  pred=mean(pred),
                                  est=case_when(pred>0.5~1,
                                                TRUE~0))
  RFdataNo[[j]] <- temp3
}

datapathYes <- datapath[201:400]
groupYes <- gsub("_estimations.txt","",datapathYes)
groupnameYes <- cbind(path=datapathYes,
                      group=groupYes,
                      data.frame(str_split(groupYes, "_", simplify = T)))
RFdataYes <- list()
for (j in 1:length(groupYes)) {
  temp1 <- read.delim(paste0("/n/holylabs/LABS/dwang_lab/Lab/T2D_meta/3.analysis/0323/RF8/",
                             datapathYes[j]),
                      header=FALSE)
  temp2 <- temp1[2:1001,]
  
  temp3 <- data.frame()
  for (i in 1:20) {
    temp11 <- temp2[(50*(i-1)+1):(50*i), ]
    temp22 <- data.frame(index=as.numeric(t(temp11[temp11$V1=="sample index", -1])),
                         true=as.numeric(t(temp11[temp11$V1=="true labels", -1])),
                         est=as.numeric(t(temp11[temp11$V1=="estimated labels", -1])),
                         pred=as.numeric(t(temp11[temp11$V1=="estimated probabilities", -1])),
                         rep=i)
    temp3 <- rbind(temp3, temp22)
  }
  
  temp3 <- temp3 %>% filter(!is.na(index)) %>%
    group_by(index) %>% summarize(true=mean(true),
                                  pred=mean(pred),
                                  est=case_when(pred>0.5~1,
                                                TRUE~0))
  RFdataYes[[j]] <- temp3
}

gdata::keep(RFdataNo, groupNo, RFdataYes, groupYes, sure=T)
save.image(file="/n/holylabs/LABS/dwang_lab/Lab/T2D_meta/3.analysis/0323/RF8/ROC_v8_1.RData")

library(pROC)
auc_No <- data.frame(group=groupNo,
                     auc=NA,
                     lower=NA,
                     upper=NA)
rocdataNo <- list()
for (i in 1:length(groupNo)) {
  temp <- roc(as.numeric(RFdataNo[[i]]$true), as.numeric(RFdataNo[[i]]$pred))
  
  auc_No[i,]$auc <- ci.auc(temp)[2]
  auc_No[i,]$lower <- ci.auc(temp)[1]
  auc_No[i,]$upper <- ci.auc(temp)[3]
  
  temp_ci <- ci.se(temp, specificities=temp$specificities)
  rocdataNo[[i]] <- data.frame(FPR=1-temp$specificities,
                               TPR=temp$sensitivities,
                               lower=temp_ci[,1],
                               upper=temp_ci[,3],
                               group=groupNo[i])
}

auc_Yes <- data.frame(group=groupYes,
                      auc=NA,
                      lower=NA,
                      upper=NA)
rocdataYes <- list()
for (i in 1:length(groupYes)) {
  temp <- roc(as.numeric(RFdataYes[[i]]$true), as.numeric(RFdataYes[[i]]$pred))
  
  auc_Yes[i,]$auc <- ci.auc(temp)[2]
  auc_Yes[i,]$lower <- ci.auc(temp)[1]
  auc_Yes[i,]$upper <- ci.auc(temp)[3]
  
  temp_ci <- ci.se(temp, specificities=temp$specificities)
  rocdataYes[[i]] <- data.frame(FPR=1-temp$specificities,
                                TPR=temp$sensitivities,
                                lower=temp_ci[,1],
                                upper=temp_ci[,3],
                                group=groupYes[i])
}

save.image(file="/n/holylabs/LABS/dwang_lab/Lab/T2D_meta/3.analysis/0323/RF8/ROC_v8_2.RData")

summary(auc_No[1:100,]$auc) # 0.69 (0.67, 0.70)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.6728  0.6806  0.6869  0.6859  0.6905  0.7010
quantile(auc_No[1:100,]$auc, prob=c(0.025, 0.975))
# 2.5%     97.5% 
# 0.6736679 0.6987122 

summary(auc_No[101:200,]$auc) # 0.74 (0.73, 0.75), 7.68%
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.7275  0.7345  0.7387  0.7386  0.7421  0.7516
quantile(auc_No[101:200,]$auc, prob=c(0.025, 0.975))
# 2.5%     97.5% 
# 0.7288677 0.7473124 

summary(auc_Yes[1:100,]$auc) # 0.77 (0.75, 0.79)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.7442  0.7618  0.7684  0.7686  0.7769  0.7964
quantile(auc_Yes[1:100,]$auc, prob=c(0.025, 0.975))
# 2.5%     97.5% 
# 0.7468889 0.7861486

summary(auc_Yes[101:200,]$auc) # 0.86 (0.85, 0.87), 12.24%
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.8488  0.8582  0.8627  0.8627  0.8671  0.8744 
quantile(auc_Yes[101:200,]$auc, prob=c(0.025, 0.975))
# 2.5%     97.5% 
# 0.8515227 0.8729832 

# kappa value comparing two models
library(psych)
kappa_No <- data.frame(type=1:100,
                       kappa=NA,
                       lower=NA,
                       upper=NA)
for (i in 1:100) {
  temp <- cohen.kappa(x=cbind(RFdataNo[[i]]$est, RFdataNo[[i+100]]$est))
  kappa_No[i,]$kappa <- temp$confid[1,2]
  kappa_No[i,]$lower <- temp$confid[1,1]
  kappa_No[i,]$upper <- temp$confid[1,3]
}

kappa_Yes <- data.frame(type=1:100,
                        kappa=NA,
                        lower=NA,
                        upper=NA)
for (i in 1:100) {
  temp <- cohen.kappa(x=cbind(RFdataYes[[i]]$est, RFdataYes[[i+100]]$est))
  kappa_Yes[i,]$kappa <- temp$confid[1,2]
  kappa_Yes[i,]$lower <- temp$confid[1,1]
  kappa_Yes[i,]$upper <- temp$confid[1,3]
}

summary(kappa_No$kappa)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.3970  0.4279  0.4384  0.4385  0.4481  0.4755
summary(kappa_Yes$kappa)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.3953  0.4480  0.4689  0.4664  0.4865  0.5264 

quantile(kappa_No$kappa, prob=c(0.025, 0.975))
# 2.5%    97.5% 
# 0.406940 0.470716
quantile(kappa_Yes$kappa, prob=c(0.025, 0.975))
# 2.5%     97.5% 
# 0.4169977 0.5077364 

t_No <- mean(kappa_No$kappa)/sd(kappa_No$kappa)
t_Yes <- mean(kappa_Yes$kappa)/sd(kappa_Yes$kappa)

save.image(file="/n/holylabs/LABS/dwang_lab/Lab/T2D_meta/3.analysis/0323/RF8/ROC_v8_3.RData")

library(ggplot2)
xs = c(0, 1)
beta = c(0, 1)
ys = cbind(1, xs) %*% beta

No_0 <- bind_rows(rocdataNo[c(1:100)])
No_1 <- bind_rows(rocdataNo[c(101:200)])
Yes_0 <- bind_rows(rocdataYes[c(1:100)])
Yes_1 <- bind_rows(rocdataYes[c(101:200)])

No_00 <- No_0 %>% group_by(FPR) %>% 
  summarize(TPR=mean(TPR)) %>%
  arrange(FPR) %>%
  mutate(rankFPR=1:1068) %>%
  arrange(TPR) %>%
  mutate(rankTPR=1:1068) %>% filter(rankFPR==rankTPR)
No_11 <- No_1 %>% group_by(FPR) %>% 
  summarize(TPR=mean(TPR)) %>%
  arrange(FPR) %>%
  mutate(rankFPR=1:1068) %>%
  arrange(TPR) %>%
  mutate(rankTPR=1:1068) %>% filter(rankFPR==rankTPR)

Yes_00 <- Yes_0 %>% group_by(FPR) %>% 
  summarize(TPR=mean(TPR)) %>%
  arrange(FPR) %>%
  mutate(rankFPR=1:616) %>%
  arrange(TPR) %>%
  mutate(rankTPR=1:616) %>% filter(rankFPR==rankTPR)
Yes_11 <- Yes_1 %>% group_by(FPR) %>% 
  summarize(TPR=mean(TPR)) %>%
  arrange(FPR) %>%
  mutate(rankFPR=1:616) %>%
  arrange(TPR) %>%
  mutate(rankTPR=1:616) %>% filter(rankFPR==rankTPR)

No_00$type0 <- "Basic model"
No_11$type0 <- "Basic model + species"
No_00$type <- "AUC: 0.69 (0.67, 0.70)"
No_11$type <- "AUC: 0.74 (0.73, 0.75)"

Yes_00$type0 <- "Basic model"
Yes_11$type0 <- "Basic model + species"
Yes_00$type <- "AUC: 0.77 (0.75, 0.79)"
Yes_11$type <- "AUC: 0.86 (0.85, 0.87)"

col0 <- c("#CC9900", "#C75127")
names(col0) <- c("Basic model", "Basic model + species")
colNo <- c("#C75127", "#CC9900")
names(colNo) <- c("AUC: 0.74 (0.73, 0.75)", "AUC: 0.69 (0.67, 0.70)")
colYes <- c("#C75127", "#CC9900")
names(colYes) <- c("AUC: 0.86 (0.85, 0.87)", "AUC: 0.77 (0.75, 0.79)")

Nodata <- rbind(No_00, No_11)
Yesdata <- rbind(Yes_00, Yes_11)

No_ROC <- ggplot() +
  geom_segment(aes(x = xs[1], xend = xs[2], y = ys[1], yend = ys[2]),
               colour = "grey",lty = 2, size=1) +
  geom_line(aes(x=FPR, y=TPR, colour = type), data=Nodata, size=1.6) +
  scale_color_manual(values = colNo) +
  scale_x_continuous(breaks = c(0,0.2,0.4,0.6,0.8,1), limits=c(0,1)) +
  scale_y_continuous(breaks = c(0,0.2,0.4,0.6,0.8,1), limits=c(0,1)) +
  xlab("") +
  ylab("True positive rate") +
  ggtitle("Con vs. T2D, metf-") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(size=2)) +
  theme(axis.title.y = element_text(size=20, colour="black"),
        axis.text.y = element_text(size=20, colour="black"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        plot.title = element_text(size=20, hjust=0.5, colour="black"),
        legend.title = element_blank(),
        legend.text = element_text(size=18),
        legend.background = element_blank(),
        legend.position = c(0.58, 0.12)) +
  guides(colour = guide_legend(keyheight = 1.5))

Yes_ROC <- ggplot() +
  geom_segment(aes(x = xs[1], xend = xs[2], y = ys[1], yend = ys[2]),
               colour = "grey",lty = 2, size=1) +
  geom_line(aes(x=FPR, y=TPR, colour = type), data=Yesdata, size=1.6) +
  scale_x_continuous(breaks = c(0,0.2,0.4,0.6,0.8,1), limits=c(0,1)) +
  scale_y_continuous(breaks = c(0,0.2,0.4,0.6,0.8,1), limits=c(0,1)) +
  scale_color_manual(values = colYes) +
  xlab("False positive rate") +
  ylab("True positive rate") +
  ggtitle("Con vs. T2D, metf+") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(size=2)) +
  theme(axis.title = element_text(size=20, colour="black"),
        axis.text = element_text(size=20, colour="black"),
        plot.title = element_text(size=20, hjust=0.5, colour="black"),
        legend.title = element_blank(),
        legend.text = element_text(size=18),
        legend.background = element_blank(),
        legend.position = c(0.58, 0.12)) +
  guides(colour = guide_legend(keyheight = 1.5))

pdf("/n/holylabs/LABS/dwang_lab/Lab/T2D_meta/3.analysis/0323/RF8/ROC_v8.pdf",
    width = 4.6, height = 9, onefile = F)
egg::ggarrange(No_ROC, Yes_ROC, nrow = 2)
dev.off()

legend <- ggplot() +
  geom_segment(aes(x = xs[1], xend = xs[2], y = ys[1], yend = ys[2]),
               colour = "grey",lty = 2, size=1) +
  geom_line(aes(x=FPR, y=TPR, colour = type0), data=Nodata, size=1.6) +
  scale_color_manual(values = col0) +
  scale_x_continuous(breaks = c(0,0.2,0.4,0.6,0.8,1), limits=c(0,1)) +
  scale_y_continuous(breaks = c(0,0.2,0.4,0.6,0.8,1), limits=c(0,1)) +
  xlab("") +
  ylab("True positive rate") +
  ggtitle("Con vs. T2D, metf-") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(size=2)) +
  theme(axis.title.y = element_text(size=20, colour="black"),
        axis.text.y = element_text(size=20, colour="black"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        plot.title = element_text(size=20, hjust=0.5, colour="black"),
        legend.title = element_blank(),
        legend.text = element_text(size=18),
        legend.background = element_blank(),
        legend.position = "right")

pdf("/n/holylabs/LABS/dwang_lab/Lab/T2D_meta/3.analysis/0323/RF8/ROC_v8_legend.pdf",
    width = 4, height = 1, onefile = F) # Open a new pdf file
grid::grid.draw(ggpubr::get_legend(legend))
dev.off() # Close the file

