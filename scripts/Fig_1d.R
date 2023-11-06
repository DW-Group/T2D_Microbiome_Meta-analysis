# based on correct batch effect overall
setwd("/n/holylabs/LABS/dwang_lab/Lab/T2D_meta/3.analysis/dropbox/fig2/")

library(tidyverse)
load("/n/holylabs/LABS/dwang_lab/Lab/T2D_meta/3.analysis/0323/species_0.RData")
summary(colSums(species_zero_adj, na.rm = T)) # should be all 1
species_zero_adj[is.na(species_zero_adj)] <- 0

gdata::keep(metadata, species_zero_adj, sure=T)

# study, region, age, sex, bmi, status_new, metformin, insulin
# tc, ldlc, hdlc, tg, tg_hdlc, fbg, finsulin, hba1c, HOMA-IR, HOMA-B
# crp, leptin, adiponectin

# HOMA-IR= (glucose in mmol/L x insulin in mIU/mL)/22.5
# HOMA-B= (20 x insulin in mIU/mL)/(glucose in mmol/L - 3.5)
metadata$homaIR <- metadata$fbg*metadata$finsulin/22.5
metadata$homaB <- ifelse(metadata$fbg<=3.5,
                         NA,
                         ((20*metadata$finsulin)/(metadata$fbg-3.5)))
metadata$tg_hdlc <- metadata$tg/metadata$hdlc

covar <- c("study", "age", "sex", "bmi", "status_new", "metf", "insul", "fbg", "finsulin", 
           "homaIR", "homaB", "hba1c", "tc", "ldlc", "hdlc", "tg",  "crp")

#####################################
#    permanova among all
#####################################
library(tidyverse)
library(vegan)
adon=list()
for (i in 1:length(covar)){
  tempvar <- covar[i]
  tempmeta <- metadata %>% select(c("id", all_of(tempvar)))
  tempmeta <- tempmeta[!is.na(tempmeta[,2]), ]
  tempid <- tempmeta$id
  tempspec <- species_zero_adj[, tempid]
  # identical(colnames(tempspec), tempmeta$id) # TRUE
  tempchk <- data.frame(spec=rownames(tempspec), sum=rowSums(tempspec))
  tt <- tempchk[tempchk$sum!=0, ]$spec
  tempspec <- tempspec[(rownames(tempspec) %in% tt),]
  tempdist <- vegdist(t(tempspec), method="bray")
  set.seed(1234)
  adon[[i]] <- adonis(as.formula(paste0("tempdist ~ ", tempvar)), data = tempmeta)
}

save.image(file="/n/holylabs/LABS/dwang_lab/Lab/T2D_meta/3.analysis/dropbox/fig2/permanova_1.RData")

#############################################
#     permanova for individual studies
#############################################
study <- as.character(unique(metadata$study))
adon_karlsson=list()
adon_mbs=list()
adon_mlvs=list()
adon_qin=list()
adon_zhong=list()
adon_direct=list()
adon_sol=list()
adon_wu1=list()
adon_wu2=list()
adon_pederson=list()

for (i in 1:length(study)) {
  studyname <- study[i]
  studymeta <- metadata[metadata$study == studyname, ]
  for (j in 1:length(covar)) {
    tempvar <- covar[j]
    tempmeta <- studymeta %>% select(c("id", all_of(tempvar)))
    tempmeta <- tempmeta[!is.na(tempmeta[,2]), ]
    if (i==1) {
      if (dim(tempmeta)[1]==0 | length(unique(tempmeta[,2]))==1) {
        adon_karlsson[[j]]=0
      } else {
        tempid <- tempmeta$id
        tempspec <- species_zero_adj[, tempid]
        # identical(colnames(tempspec), tempmeta$id) # TRUE
        tempchk=data.frame(spec=rownames(tempspec), sum=rowSums(tempspec))
        tt=tempchk[tempchk$sum!=0, ]$spec
        tempspec=tempspec[(rownames(tempspec) %in% tt),]
        tempdist=vegdist(t(tempspec))
        set.seed(1234)
        adon_karlsson[[j]]=adonis(as.formula(paste0("tempdist ~ ", tempvar)), data = tempmeta)
      }
    } else if (i==2) {
      if (dim(tempmeta)[1]==0 | length(unique(tempmeta[,2]))==1) {
        adon_mbs[[j]]=0
      } else {
        tempid <- tempmeta$id
        tempspec <- species_zero_adj[, tempid]
        # identical(colnames(tempspec), tempmeta$id) # TRUE
        tempchk=data.frame(spec=rownames(tempspec), sum=rowSums(tempspec))
        tt=tempchk[tempchk$sum!=0, ]$spec
        tempspec=tempspec[(rownames(tempspec) %in% tt),]
        tempdist=vegdist(t(tempspec))
        set.seed(1234)
        adon_mbs[[j]]=adonis(as.formula(paste0("tempdist ~ ", tempvar)), data = tempmeta)
      }
    } else if (i==3) {
      if (dim(tempmeta)[1]==0 | length(unique(tempmeta[,2]))==1) {
        adon_mlvs[[j]]=0
      } else {
        tempid <- tempmeta$id
        tempspec <- species_zero_adj[, tempid]
        # identical(colnames(tempspec), tempmeta$id) # TRUE
        tempchk=data.frame(spec=rownames(tempspec), sum=rowSums(tempspec))
        tt=tempchk[tempchk$sum!=0, ]$spec
        tempspec=tempspec[(rownames(tempspec) %in% tt),]
        tempdist=vegdist(t(tempspec))
        set.seed(1234)
        adon_mlvs[[j]]=adonis(as.formula(paste0("tempdist ~ ", tempvar)), data = tempmeta)
      }
    } else if (i==4) {
      if (dim(tempmeta)[1]==0 | length(unique(tempmeta[,2]))==1) {
        adon_qin[[j]]=0
      } else {
        tempid <- tempmeta$id
        tempspec <- species_zero_adj[, tempid]
        # identical(colnames(tempspec), tempmeta$id) # TRUE
        tempchk=data.frame(spec=rownames(tempspec), sum=rowSums(tempspec))
        tt=tempchk[tempchk$sum!=0, ]$spec
        tempspec=tempspec[(rownames(tempspec) %in% tt),]
        tempdist=vegdist(t(tempspec))
        set.seed(1234)
        adon_qin[[j]]=adonis(as.formula(paste0("tempdist ~ ", tempvar)), data = tempmeta)
      }
    } else if (i==5) {
      if (dim(tempmeta)[1]==0 | length(unique(tempmeta[,2]))==1) {
        adon_zhong[[j]]=0
      } else {
        tempid <- tempmeta$id
        tempspec <- species_zero_adj[, tempid]
        # identical(colnames(tempspec), tempmeta$id) # TRUE
        tempchk=data.frame(spec=rownames(tempspec), sum=rowSums(tempspec))
        tt=tempchk[tempchk$sum!=0, ]$spec
        tempspec=tempspec[(rownames(tempspec) %in% tt),]
        tempdist=vegdist(t(tempspec))
        set.seed(1234)
        adon_zhong[[j]]=adonis(as.formula(paste0("tempdist ~ ", tempvar)), data = tempmeta)
      }
    }  else if (i==6) {
      if (dim(tempmeta)[1]==0 | length(unique(tempmeta[,2]))==1) {
        adon_direct[[j]]=0
      } else {
        tempid <- tempmeta$id
        tempspec <- species_zero_adj[, tempid]
        # identical(colnames(tempspec), tempmeta$id) # TRUE
        tempchk=data.frame(spec=rownames(tempspec), sum=rowSums(tempspec))
        tt=tempchk[tempchk$sum!=0, ]$spec
        tempspec=tempspec[(rownames(tempspec) %in% tt),]
        tempdist=vegdist(t(tempspec))
        set.seed(1234)
        adon_direct[[j]]=adonis(as.formula(paste0("tempdist ~ ", tempvar)), data = tempmeta)
      }
    }  else if (i==7) {
      if (dim(tempmeta)[1]==0 | length(unique(tempmeta[,2]))==1) {
        adon_sol[[j]]=0
      } else {
        tempid <- tempmeta$id
        tempspec <- species_zero_adj[, tempid]
        # identical(colnames(tempspec), tempmeta$id) # TRUE
        tempchk=data.frame(spec=rownames(tempspec), sum=rowSums(tempspec))
        tt=tempchk[tempchk$sum!=0, ]$spec
        tempspec=tempspec[(rownames(tempspec) %in% tt),]
        tempdist=vegdist(t(tempspec))
        set.seed(1234)
        adon_sol[[j]]=adonis(as.formula(paste0("tempdist ~ ", tempvar)), data = tempmeta)
      }
    }  else if (i==8) {
      if (dim(tempmeta)[1]==0 | length(unique(tempmeta[,2]))==1) {
        adon_wu1[[j]]=0
      } else {
        tempid <- tempmeta$id
        tempspec <- species_zero_adj[, tempid]
        # identical(colnames(tempspec), tempmeta$id) # TRUE
        tempchk=data.frame(spec=rownames(tempspec), sum=rowSums(tempspec))
        tt=tempchk[tempchk$sum!=0, ]$spec
        tempspec=tempspec[(rownames(tempspec) %in% tt),]
        tempdist=vegdist(t(tempspec))
        set.seed(1234)
        adon_wu1[[j]]=adonis(as.formula(paste0("tempdist ~ ", tempvar)), data = tempmeta)
      }
    } else if (i==9) {
      if (dim(tempmeta)[1]==0 | length(unique(tempmeta[,2]))==1) {
        adon_wu2[[j]]=0
      } else {
        tempid <- tempmeta$id
        tempspec <- species_zero_adj[, tempid]
        # identical(colnames(tempspec), tempmeta$id) # TRUE
        tempchk=data.frame(spec=rownames(tempspec), sum=rowSums(tempspec))
        tt=tempchk[tempchk$sum!=0, ]$spec
        tempspec=tempspec[(rownames(tempspec) %in% tt),]
        tempdist=vegdist(t(tempspec))
        set.seed(1234)
        adon_wu2[[j]]=adonis(as.formula(paste0("tempdist ~ ", tempvar)), data = tempmeta)
      }
    } else if (i==10) {
      if (dim(tempmeta)[1]==0 | length(unique(tempmeta[,2]))==1) {
        adon_pederson[[j]]=0
      } else {
        tempid <- tempmeta$id
        tempspec <- species_zero_adj[, tempid]
        # identical(colnames(tempspec), tempmeta$id) # TRUE
        tempchk=data.frame(spec=rownames(tempspec), sum=rowSums(tempspec))
        tt=tempchk[tempchk$sum!=0, ]$spec
        tempspec=tempspec[(rownames(tempspec) %in% tt),]
        tempdist=vegdist(t(tempspec))
        set.seed(1234)
        adon_pederson[[j]]=adonis(as.formula(paste0("tempdist ~ ", tempvar)), data = tempmeta)
      }
    }
  }
}

save.image(file="/n/holylabs/LABS/dwang_lab/Lab/T2D_meta/3.analysis/dropbox/fig2/permanova_2.RData")

gdata::keep(covar, study, metadata, species_zero_adj,
            adon, adon_karlsson, adon_direct, adon_mbs, adon_mlvs,
            adon_qin, adon_zhong, adon_sol, adon_wu1, adon_wu2, adon_pederson, sure=T)

##########################
#     extract results
##########################
permanova <- data.frame(covar=covar,
                        lab=NA,
                        All=NA,
                        DIRECT_PLUS=NA,
                        Fromentin=NA,
                        HPFS=NA,
                        Karlsson=NA,
                        NHS2=NA,
                        Qin=NA,
                        SOL=NA,
                        Wu1=NA,
                        Wu2=NA,
                        Zhong=NA)

for (i in 1:length(covar)){
  if(class(adon[[i]]) == "adonis"){
    permanova[i,"All"]=round(adon[[i]]$aov.tab$R2[1],4)*100
  } else {
    permanova[i,"All"]=NA
  }
}

for (i in 1:length(covar)){
  if(class(adon_direct[[i]]) == "adonis"){
    permanova[i,"DIRECT_PLUS"]=round(adon_direct[[i]]$aov.tab$R2[1],4)*100
  } else {
    permanova[i,"DIRECT_PLUS"]=NA
  }
}

for (i in 1:length(covar)){
  if(class(adon_karlsson[[i]]) == "adonis"){
    permanova[i,"Karlsson"]=round(adon_karlsson[[i]]$aov.tab$R2[1],4)*100
  } else {
    permanova[i,"Karlsson"]=NA
  }
}

for (i in 1:length(covar)){
  if(class(adon_mbs[[i]]) == "adonis"){
    permanova[i,"NHS2"]=round(adon_mbs[[i]]$aov.tab$R2[1],4)*100
  } else {
    permanova[i,"NHS2"]=NA
  }
}

for (i in 1:length(covar)){
  if(class(adon_mlvs[[i]]) == "adonis"){
    permanova[i,"HPFS"]=round(adon_mlvs[[i]]$aov.tab$R2[1],4)*100
  } else {
    permanova[i,"HPFS"]=NA
  }
}

for (i in 1:length(covar)){
  if(class(adon_qin[[i]]) == "adonis"){
    permanova[i,"Qin"]=round(adon_qin[[i]]$aov.tab$R2[1],4)*100
  } else {
    permanova[i,"Qin"]=NA
  }
}

for (i in 1:length(covar)){
  if(class(adon_zhong[[i]]) == "adonis"){
    permanova[i,"Zhong"]=round(adon_zhong[[i]]$aov.tab$R2[1],4)*100
  } else {
    permanova[i,"Zhong"]=NA
  }
}

for (i in 1:length(covar)){
  if(class(adon_sol[[i]]) == "adonis"){
    permanova[i,"SOL"]=round(adon_sol[[i]]$aov.tab$R2[1],4)*100
  } else {
    permanova[i,"SOL"]=NA
  }
}

for (i in 1:length(covar)){
  if(class(adon_pederson[[i]]) == "adonis"){
    permanova[i,"Fromentin"]=round(adon_pederson[[i]]$aov.tab$R2[1],4)*100
  } else {
    permanova[i,"Fromentin"]=NA
  }
}

for (i in 1:length(covar)){
  if(class(adon_wu1[[i]]) == "adonis"){
    permanova[i,"Wu1"]=round(adon_wu1[[i]]$aov.tab$R2[1],4)*100
  } else {
    permanova[i,"Wu1"]=NA
  }
}

for (i in 1:length(covar)){
  if(class(adon_wu2[[i]]) == "adonis"){
    permanova[i,"Wu2"]=round(adon_wu2[[i]]$aov.tab$R2[1],4)*100
  } else {
    permanova[i,"Wu2"]=NA
  }
}

permanova$lab <- c("Study", "Age", "Sex", "BMI", "Diabetes status",
                   "Metformin use", "Insulin use", "Fasting glucose", 
                   "Fasting insulin", "HOMA-IR", "HOMA-B", "HbA1c", 
                   "Total cholesterol", "LDL cholesterol",
                   "HDL cholesterol", "Triacylglyceride", "hs-CRP")

##########################
#     extract p
##########################
permanova_p <- data.frame(covar=covar,
                          lab=NA,
                          All=NA,
                          DIRECT_PLUS=NA,
                          Fromentin=NA,
                          HPFS=NA,
                          Karlsson=NA,
                          NHS2=NA,
                          Qin=NA,
                          SOL=NA,
                          Wu1=NA,
                          Wu2=NA,
                          Zhong=NA)

for (i in 1:length(covar)){
  if(class(adon[[i]]) == "adonis"){
    permanova_p[i,"All"]=adon[[i]]$aov.tab$`Pr(>F)`[1]
  } else {
    permanova_p[i,"All"]=NA
  }
}

for (i in 1:length(covar)){
  if(class(adon_direct[[i]]) == "adonis"){
    permanova_p[i,"DIRECT_PLUS"]=adon_direct[[i]]$aov.tab$`Pr(>F)`[1]
  } else {
    permanova_p[i,"DIRECT_PLUS"]=NA
  }
}

for (i in 1:length(covar)){
  if(class(adon_karlsson[[i]]) == "adonis"){
    permanova_p[i,"Karlsson"]=adon_karlsson[[i]]$aov.tab$`Pr(>F)`[1]
  } else {
    permanova_p[i,"Karlsson"]=NA
  }
}

for (i in 1:length(covar)){
  if(class(adon_mbs[[i]]) == "adonis"){
    permanova_p[i,"NHS2"]=adon_mbs[[i]]$aov.tab$`Pr(>F)`[1]
  } else {
    permanova_p[i,"NHS2"]=NA
  }
}

for (i in 1:length(covar)){
  if(class(adon_mlvs[[i]]) == "adonis"){
    permanova_p[i,"HPFS"]=adon_mlvs[[i]]$aov.tab$`Pr(>F)`[1]
  } else {
    permanova_p[i,"HPFS"]=NA
  }
}

for (i in 1:length(covar)){
  if(class(adon_qin[[i]]) == "adonis"){
    permanova_p[i,"Qin"]=adon_qin[[i]]$aov.tab$`Pr(>F)`[1]
  } else {
    permanova_p[i,"Qin"]=NA
  }
}

for (i in 1:length(covar)){
  if(class(adon_zhong[[i]]) == "adonis"){
    permanova_p[i,"Zhong"]=adon_zhong[[i]]$aov.tab$`Pr(>F)`[1]
  } else {
    permanova_p[i,"Zhong"]=NA
  }
}

for (i in 1:length(covar)){
  if(class(adon_sol[[i]]) == "adonis"){
    permanova_p[i,"SOL"]=adon_sol[[i]]$aov.tab$`Pr(>F)`[1]
  } else {
    permanova_p[i,"SOL"]=NA
  }
}

for (i in 1:length(covar)){
  if(class(adon_pederson[[i]]) == "adonis"){
    permanova_p[i,"Fromentin"]=adon_pederson[[i]]$aov.tab$`Pr(>F)`[1]
  } else {
    permanova_p[i,"Fromentin"]=NA
  }
}

for (i in 1:length(covar)){
  if(class(adon_wu1[[i]]) == "adonis"){
    permanova_p[i,"Wu1"]=adon_wu1[[i]]$aov.tab$`Pr(>F)`[1]
  } else {
    permanova_p[i,"Wu1"]=NA
  }
}

for (i in 1:length(covar)){
  if(class(adon_wu2[[i]]) == "adonis"){
    permanova_p[i,"Wu2"]=adon_wu2[[i]]$aov.tab$`Pr(>F)`[1]
  } else {
    permanova_p[i,"Wu2"]=NA
  }
}

permanova_p$lab <- c("Study", "Age", "Sex", "BMI", "Diabetes status",
                   "Metformin use", "Insulin use", "Fasting glucose", 
                   "Fasting insulin", "HOMA-IR", "HOMA-B", "HbA1c", 
                   "Total cholesterol", "LDL cholesterol",
                   "HDL cholesterol", "Triacylglyceride",  "hs-CRP")

library(openxlsx)
write.xlsx(list('r2'=permanova,
                'p'=permanova_p), 
           file="/n/holylabs/LABS/dwang_lab/Lab/T2D_meta/3.analysis/dropbox/fig2/permanova.xlsx")
save.image(file="/n/holylabs/LABS/dwang_lab/Lab/T2D_meta/3.analysis/dropbox/fig2/permanova_3.RData")

