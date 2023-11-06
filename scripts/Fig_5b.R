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
select <- dplyr::select


# manual heatmap E.coli-------
anpan_re_all <-
  fread("anpan_batch_bmi_minmax5_2023-04-04/all_bug_gene_terms.tsv.gz")
head(anpan_re_all);dim(anpan_re_all)

# anpan_sig <- anpan_re_all %>% 
#   filter(bug_name=="Escherichia_coli") %>% 
#   filter(q_global < 0.25)
# dim(anpan_sig) #1439/11
# 
# anpan_ecoli <- anpan_re_all %>% 
#   filter(bug_name=="Escherichia_coli") 
# dim(anpan_ecoli) #77703/11; 105510/11

top_spe <- fread("species_list/species_for_strain.txt",header = FALSE)
head(top_spe);dim(top_spe) #70
(top_species <- top_spe$V1)
top_species2 <- gsub("s__","",top_species)
(bugName=top_species2[36]) #36

## filtered gene data --------
## samples should have certain percent (like 1%) of genes
perct_samp <- 0.01

## gens should present at certain percent (like 5%) of samples
perct_gene <- 0.05  
beta_thresh <- 1.2   

anpan_re_sig <- anpan_re_all %>% 
  filter(bug_name == bugName) %>% 
  filter(abs(estimate) > beta_thresh)
dim(anpan_re_sig) 

(anpan_test <- filter(anpan_re_sig,abs(estimate) > 2))

model_input <-
  fread(paste0("anpan_batch_bmi_minmax5_2023-04-04/filter_stats/filtered_",
               bugName,".tsv.gz")) %>% as.data.frame()
head(model_input[,1:10]);
dim(model_input) 

## prevalence difference in case-control/total prevalence -------
# prevalence among all samples
pre_t <- data.frame(
  gene = colnames(model_input)[7:ncol(model_input)],
  pre_t = colSums(model_input[,7:ncol(model_input)])/nrow(model_input)
)
head(pre_t);dim(pre_t)

## prevalence among case/control
head(anpan_re_sig$gene)
gene_pre_t2d_con <- model_input %>%
  select(c("T2D",anpan_re_sig$gene)) %>% 
  group_by(T2D) %>% 
  summarise_all(mean) %>% 
  column_to_rownames(var = "T2D") %>% 
  sjmisc::rotate_df() %>% 
  rownames_to_column(var = "gene")
dim(gene_pre_t2d_con) #610/3
head(gene_pre_t2d_con)


## diff of prevalence
gene_pre_t2d_con <- gene_pre_t2d_con %>% 
  left_join(pre_t,by="gene")
head(gene_pre_t2d_con);dim(gene_pre_t2d_con) #610/4
colSums(is.na(gene_pre_t2d_con))
gene_pre_t2d_con$diff_perc <- 
  abs(gene_pre_t2d_con$`FALSE` - gene_pre_t2d_con$`TRUE`)*100/gene_pre_t2d_con$pre_t
summary(gene_pre_t2d_con$diff_perc)
table(gene_pre_t2d_con$diff_perc>50)
table(gene_pre_t2d_con$diff_perc>80)
table(gene_pre_t2d_con$diff_perc>100)
gene_pre_t2d_con_diff50 <- filter(gene_pre_t2d_con,diff_perc<50)
gene_pre_t2d_con_diff80 <- filter(gene_pre_t2d_con,diff_perc<80)
gene_pre_t2d_con_diff100 <- filter(gene_pre_t2d_con,diff_perc<100)


#### manually remove genes ------
genes_del_v0 <- c(
  "UniRef90_UPI000CE68AB0: NO_NAME",
  "UniRef90_A0A2U2TTM2: NO_NAME",
  "UniRef90_UPI000D67ECB0: NO_NAME",
  "UniRef90_UPI000C9AC6AF: NO_NAME"
)

genes_del_v1 <- c(
  "UniRef90_Q1RFY5: NO_NAME",
  "UniRef90_A0A0E9NSK5: NO_NAME",
  "UniRef90_Q7WUC8: NO_NAME",
  "UniRef90_A0A376FSK0: DNA primase TraC",
  "UniRef90_H4UU75: NO_NAME",
  "UniRef90_A0A376VWG4: TraJ transfer ATPase",
  "UniRef90_A0A2X1LQD8: RNA polymerase factor sigma-32",
  "UniRef90_A0A1V0DNP4: IncI1 plasmid conjugative transfer pilus-tip adhesin protein PilV",
  "UniRef90_UPI000C7A76F2: NO_NAME",
  "UniRef90_A0A2X1NJZ4: YgaA",
  "UniRef90_A0A376UXS0: NO_NAME",
  "UniRef90_A0A376UXB0: ParB-like nuclease",
  "UniRef90_UPI000CD59623: phosphopantetheine--protein transferase",
  "UniRef90_A0A2Y8S5Q6: NO_NAME",
  "UniRef90_UPI000B3F5AC9: malate transporter",
  
  "UniRef90_A0A1X0YJ48: NO_NAME",
  "UniRef90_UPI0004D401BF: NO_NAME",
  "UniRef90_A0A376P4T9: Transcription regulatory protein",
  "UniRef90_A0A2X7GFW0: Flagellin",
  "UniRef90_A0A376M4X1: Multiphosphoryl transfer protein 1 [includes phosphoenolpyruvate-protein phosphotransferasephosphocarrier protein Hp fructose-like phosphotransferase enzyme IIA component]",
  "UniRef90_A0A376VFB7: Sulfite reductase (NADPH) hemoprotein beta subunit",
  "UniRef90_A0A376K377: Adhesin for cattle intestine colonization",
  "UniRef90_A0A368IY81: NO_NAME",
  "UniRef90_A0A073HVW9: RHS Repeat family protein (Fragment)",
  "UniRef90_A0A073HS85: RHS Repeat family protein (Fragment)",
  "UniRef90_A0A377LIL4: Putative adhesin/invasin-like protein",
  "UniRef90_A0A2S8JX43: Fimbrial protein",
  "UniRef90_A0A1X0YHX1: NO_NAME",
  "UniRef90_A0A1X0YUQ4: NO_NAME",
  "UniRef90_I6FIA7: NO_NAME",
  "UniRef90_A0A377LL65: RhsD protein",
  "UniRef90_A0A1X0YI09: Outer membrane usher protein htrE",
  "UniRef90_A0A1X0YIJ1: DsORF-h1",
  "UniRef90_UPI0003F95708: fimbrial-like adhesin",
  "UniRef90_A0A368IJP4: DUF2544 domain-containing protein",
  "UniRef90_A0A376VQZ5: Putative type III effector protein",
  "UniRef90_A0A1X0YM65: NO_NAME",
  "UniRef90_A0A377LK00: Outer membrane usher protein",
  "UniRef90_A0A377LJQ9: Outer membrane usher protein",
  "UniRef90_A0A1X0YW45: Type III effector (Fragment)",
  "UniRef90_A0A377LE58: Protein",
  "UniRef90_A0A0I2ABB5: NO_NAME",
  "UniRef90_A0A244B9E1: NO_NAME",
  "UniRef90_A0A376VU91: Putative type III effector protein",
  "UniRef90_A0A377LNX1: Putative lipoprotein",
  "UniRef90_A0A1X0YIP5: Type II secretion system protein GspD (Fragment)",
  "UniRef90_A0A377CZP7: Putative peptidoglycan-binding protein",
  "UniRef90_A0A376JYI7: Low specificity L-threonine aldolase",
  "UniRef90_A0A376K6T5: Protein involved in detoxification of methylglyoxal",
  "UniRef90_A0A377DK98: Putrescine ABC transporter membrane protein",
  "UniRef90_UPI000A38FC0C: autotransporter outer membrane beta-barrel domain-containing protein",
  "UniRef90_A0A376U2K1: Thiamine ABC transport system, permease protein",
  "UniRef90_A0A376VLE8: Outer membrane protein",
  "UniRef90_A0A0K3NRK6: Deubiquitinase",
  "UniRef90_A0A369F042: Inverse autotransporter adhesin-like protein YeeJ (Fragment)",
  "UniRef90_A0A376VF15: Hydrogenase maturation protein hypF",
  "UniRef90_A0A0P7NFQ6: NO_NAME",
  "UniRef90_A0A376P0L1: ATP-dependent chaperone protein ClpB",
  "UniRef90_A0A399BT75: Tail assembly protein",
  "UniRef90_A0A376ZE81: Host specificity protein",
  "UniRef90_A0A377EF00: Portal protein",
  "UniRef90_A0A0G3KF41: NO_NAME",
  "UniRef90_A0A376KCC7: Type I restriction-modification system,specificity subunit S",
  "UniRef90_UPI00067D5D33: NO_NAME",
  "UniRef90_A0A2X1KVU7: Anthranilate synthase component II [includes glutamine amidotransferase anthranilate phosphoribosyltransferase]",
  "UniRef90_A0A243V6Z4: NO_NAME",
  "UniRef90_A0A376LB53: Inner membrane lipoprotein",
  "UniRef90_A0A376PS53: Protein of uncharacterized function (DUF2724)",
  "UniRef90_A0A1V2STD3: Cro/Cl family transcriptional regulator",
  "UniRef90_A0A0T5XGF1: NO_NAME",
  "UniRef90_A0A0T5XFE1: Adhesin",
  "UniRef90_A0A0T5XRL3: NO_NAME",
  "UniRef90_L0M6S3: NO_NAME",
  "UniRef90_A0A0A7RLM0: DNA polymerase V subunit UmuC",
  "UniRef90_A0A0T5XMV5: LysR family transcriptional regulator",
  "UniRef90_A0A1D3U0I6: NO_NAME",
  "UniRef90_A0A1X9TNY9: Phage portal protein",
  "UniRef90_A0A1D3U0G9: NO_NAME",
  "UniRef90_A0A2N5AES8: NO_NAME",
  "UniRef90_A0A1X9TPE6: DNA-invertase",
  "UniRef90_A0A2X3H607: Short-chain dehydrogenase/reductase SDR",
  "UniRef90_A0A0K4B412: NO_NAME",
  "UniRef90_UPI0003EDF9E3: NO_NAME",
  "UniRef90_A0A080JDU4: NO_NAME",
  "UniRef90_A0A271R8F9: NO_NAME",
  "UniRef90_A0A0B1GTH9: NO_NAME",
  "UniRef90_A0A1X3KCR7: NO_NAME",
  "UniRef90_A0A2H5C3Z9: DNA polymerase III subunit epsilon",
  
  "UniRef90_A0A073GZP3: NO_NAME",
  "UniRef90_A0A377DGU6: Putative flagellin structural protein",
  "UniRef90_A0A073UV16: NO_NAME",
  "UniRef90_A0A0D0LT16: NO_NAME",
  "UniRef90_A0A0A1AHV5: HK97 gp10 family phage protein",
  "UniRef90_A0A2S7GG14: Anaerobic C4-dicarboxylate transporter",
  "UniRef90_A0A376SHL1: Integral membrane protein-component of typeIII secretion apparatus",
  "UniRef90_A0A2X1KVF2: Putative oxidoreductase",
  "UniRef90_A0A073UHD2: NO_NAME",
  "UniRef90_A0A1D7Q682: Uncharacterized conserved protein",
  "UniRef90_A0A181WU75: NO_NAME",
  "UniRef90_A0A377LD74: Type-1 fimbrial protein",
  "UniRef90_A0A376LG26: Putative side tail phage protein",
  "UniRef90_D8A928: NO_NAME",
  "UniRef90_A0A0A1A518: NO_NAME",
  "UniRef90_A0A073H4P9: NO_NAME",
  "UniRef90_A0A1X0PCQ6: Type III effector protein",
  "UniRef90_A0A1D7PVX0: Fimbrial protein",
  "UniRef90_A0A0A1A6U8: Glycosyl hydrolases family 43",
  "UniRef90_A0A0A1A9B0: AAA domain protein",
  "UniRef90_A0A171D2Z3: RhsD protein",
  "UniRef90_A0A0A1ACG3: Baseplate protein",
  "UniRef90_A0A369G4Y5: Type VI secretion system protein TssA",
  "UniRef90_A0A0A1A357: NO_NAME",
  "UniRef90_A0A1D7PZ26: NO_NAME",
  "UniRef90_A0A3A3E6I0: NO_NAME",
  "UniRef90_A0A376D3H2: Terminase, ATPase subunit",
  "UniRef90_UPI000B5886B6: NO_NAME",
  "UniRef90_A0A1M1N6E4: Phage tail protein (Fragment)",
  "UniRef90_A0A0K4T2L4: Glyceraldehyde-3-phosphate dehydrogenase",
  "UniRef90_A0A2X2HVB8: Transcriptional antiterminator BglG",
  "UniRef90_E9TND0: NO_NAME",
  "UniRef90_A0A1B0QY72: NO_NAME",
  "UniRef90_A0A192CCA8: NO_NAME",
  "UniRef90_A0A376LAJ1: Type III secretion system protein",
  
  "UniRef90_UPI000B7AF189: colanic acid biosynthesis acetyltransferase wcaB domain protein",
  "UniRef90_A0A080GCT2: NO_NAME",
  "UniRef90_A0A3A6RRJ9: NO_NAME",
  "UniRef90_A0A2X1LDS3: NO_NAME",
  "UniRef90_I2X925: NO_NAME",
  "UniRef90_A0A2Y0APK8: NO_NAME",
  "UniRef90_A0A0A1AI49: NO_NAME",
  "UniRef90_A0A369DRR9: NO_NAME",
  "UniRef90_A0A377BB94: NO_NAME",
  "UniRef90_A0A377BCL2: Phage tail core protein",
  "UniRef90_A0A0A3UJL8: Nuclease",
  "UniRef90_A0A0A3TNR2: ATP-binding sugar transporter from pro-phage family protein",
  "UniRef90_A0A0A1AI80: Phage portal protein",
  "UniRef90_A0A367DDB1: ATP-dependent Clp protease proteolytic subunit",
  "UniRef90_A0A377BEY1: Putative phage tail sheath protein",
  "UniRef90_A0A377BC22: Phage terminase large subunit (GpA)",
  "UniRef90_A0A0A1AEI5: NO_NAME",
  "UniRef90_A0A0A3TKD9: Phage protein U",
  "UniRef90_A0A367DAU2: NO_NAME",
  "UniRef90_A0A0A3TKL8: ATPase",
  "UniRef90_A0A2X1LCG3: Protein of uncharacterized function (DUF1064)",
  "UniRef90_A0A0A1AKA0: NO_NAME",
  "UniRef90_A0A2H5C0V6: NO_NAME",
  "UniRef90_A0A0A1AF34: Phage terminase large subunit (GpA)",
  "UniRef90_A0A0V9GPZ4: Nuclease",
  "UniRef90_A0A377AKI1: NO_NAME",
  "UniRef90_A0A0A3TNQ9: NO_NAME",
  "UniRef90_A0A0A1AF19: NO_NAME",
  "UniRef90_UPI000BAADF66: NO_NAME",
  "UniRef90_UPI000BAAC12D: NO_NAME",
  "UniRef90_A0A377AKI8: NO_NAME",
  "UniRef90_A0A0A1AK91: NO_NAME",
  "UniRef90_A0A2A2C269: NO_NAME",
  "UniRef90_A0A0A1AF39: Capsid protein",
  "UniRef90_I2XHI3: Phage portal protein, lambda-like family protein (Fragment)",
  "UniRef90_A0A377BDL4: Chromosome segregation protein SMC",
  "UniRef90_A0A2H5C0M2: NO_NAME",
  "UniRef90_A0A0A1AKB5: Baseplate assembly protein",
  "UniRef90_A0A1X3HWZ1: NO_NAME",
  "UniRef90_A0A0A3TNT0: NO_NAME",
  "UniRef90_A0A0A1AFQ4: Baseplate J protein",
  "UniRef90_A0A1T1J9S7: Phage tail tape measure protein",
  "UniRef90_A0A0V9GQT1: NO_NAME",
  "UniRef90_A0A377BB85: NO_NAME",
  "UniRef90_A0A2X1MND2: NO_NAME",
  "UniRef90_H4UU27: NO_NAME",
  "UniRef90_A0A377BCS7: NO_NAME",
  "UniRef90_A0A377BCB4: NO_NAME",
  "UniRef90_A0A377BFT5: NO_NAME",
  "UniRef90_UPI00044E60EF: NO_NAME",
  "UniRef90_A0A377AJZ2: NO_NAME",
  "UniRef90_A0A377AIY7: ParB-like nuclease",
  "UniRef90_H4UU55: NO_NAME",
  "UniRef90_A0A0A1AFN1: NO_NAME",
  "UniRef90_A0A0A1AF24: NO_NAME",
  "UniRef90_A0A377BBS5: NO_NAME",
  "UniRef90_A0A377BDR3: NO_NAME",
  "UniRef90_A0A0A1AEJ0: Late control D family protein",
  "UniRef90_A0A377AKJ8: NO_NAME",
  "UniRef90_I2XKC5: Phage terminase large subunit GpA-like protein",
  "UniRef90_A0A377BBQ8: Chromosome segregation protein SMC",
  "UniRef90_A0A0V9GQ27: NO_NAME",
  "UniRef90_A0A2X1LCH9: NO_NAME",
  "UniRef90_A0A377BF30: Phage-related baseplate assembly protein J",
  "UniRef90_A0A377AJP2: Phage late control D family protein",
  "UniRef90_A0A377BD84: Putative control protein for phage late genes expression",
  "UniRef90_UPI000DA568DE: late control D family protein",
  "UniRef90_A0A1V3UZT0: Helix-turn-helix domain-containing protein",
  "UniRef90_A0A377BBK2: Phage baseplate assembly protein V",
  "UniRef90_A0A2X1N8G5: Multidrug ABC transporter ATPase/permease"
)

genes_del_v2 <- c(
  "UniRef90_A0A369DCW7: NO_NAME",
  "UniRef90_A0A0D7C3B7: NO_NAME",
  "UniRef90_T2FLI5: Mrr",
  "UniRef90_S4WZ38: NO_NAME",
  "UniRef90_UPI000B3F0A92: NO_NAME",
  "UniRef90_A0A1S6GKB0: Transcriptional regulator",
  "UniRef90_UPI00079FEAD7: NO_NAME",
  "UniRef90_UPI000DDCCF9B: chromosome partitioning protein ParB",
  "UniRef90_A0A376LQA0: Para-aminobenzoate synthase component I",
  "UniRef90_A0A376JGC4: Putative flavoprotein",
  "UniRef90_A0A377D5L8: DeoR-family transcriptional regulator",
  "UniRef90_A0A243UU65: NO_NAME",
  "UniRef90_A0A376SPE5: Tail component of prophage CP-933R",
  "UniRef90_I2S6C7: Phage minor tail protein U-like protein",
  "UniRef90_A0A376Q3H9: Major head protein",
  "UniRef90_A0A3A6S559: NO_NAME",
  "UniRef90_UPI000DAE5CCF: NO_NAME",
  "UniRef90_UPI000C0BD5DD: NO_NAME",
  "UniRef90_UPI000CF4DC9E: NO_NAME",
  "UniRef90_A0A2X1PQR4: NO_NAME",
  "UniRef90_A0A376SLF4: Exonuclease family protein",
  "UniRef90_UPI000B426D7C: type II toxin-antitoxin system RelE/ParE family toxin",
  "UniRef90_G4VUP5: Transposase IS1, ORF A",
  "UniRef90_L3P026: Rhs element Vgr protein",
  "UniRef90_UPI000C1FCD3F: peptidase",
  "UniRef90_UPI000BE5F8DD: peptidase",
  "UniRef90_UPI0003A9B647: PerC family transcriptional regulator",
  "UniRef90_UPI0003EE6164: PerC family transcriptional regulator",
  "UniRef90_A0A1D7Q7V5: NO_NAME",
  "UniRef90_A0A2X6FP33: HNH endonuclease family protein",
  "UniRef90_A0A2T1LN63: DUF2496 domain-containing protein",
  "UniRef90_A0A2T1LF01: Transcriptional regulator",
  "UniRef90_A0A376R1A4: Fimbrillin",
  "UniRef90_A0A0P7NZF5: NO_NAME",
  "UniRef90_A0A142EB17: NO_NAME",
  "UniRef90_A0A142EAW6: NO_NAME",
  "UniRef90_A0A142EB00: Antirestriction protein",
  "UniRef90_A0A2X9AC13: NO_NAME",
  "UniRef90_A0A2H5BZN6: NO_NAME",
  "UniRef90_A0A2X6DIF7: NO_NAME",
  "UniRef90_A0A0Q3IGB9: Lactate permease",
  "UniRef90_A0A237NFB6: Phage tail assembly protein",
  "UniRef90_A0A2S8HUF1: AraC family transcriptional regulator",
  "UniRef90_A0A2T1LGG1: MFS transporter",
  "UniRef90_A0A2T1LGG2: Protein YhiD",
  "UniRef90_A0A3F3QPR3: NO_NAME",
  "UniRef90_A0A376QCT1: TrbC-like protein",
  "UniRef90_A0A383HTJ0: NO_NAME",
  "UniRef90_A0A1X3KKY1: NO_NAME",
  "UniRef90_A0A0L7AL35: APC family permease",
  "UniRef90_A0A2X7RDK8: IS orf",
  "UniRef90_B7LBM5: Putative terminase small subunit",
  "UniRef90_UPI00083D4F62: phage tail protein",
  "UniRef90_E2QDD9: NO_NAME",
  "UniRef90_A0A2Y0NAD2: NO_NAME",
  "UniRef90_A0A377EIY9: Protein",
  "UniRef90_A0A0V9S860: NO_NAME",
  "UniRef90_A0A2P6IED0: Transcriptional regulator",
  "UniRef90_V0VF20: NO_NAME",
  "UniRef90_UPI000BE507A4: NO_NAME",
  "UniRef90_A0A2U9TP66: Ash family protein",
  "UniRef90_A0A209NSH2: Phage tail protein",
  "UniRef90_UPI000D0BB317: phage tail tape measure protein",
  "UniRef90_UPI000D6DD024: phage tail tape measure protein",
  "UniRef90_A0A377DEC5: NO_NAME",
  "UniRef90_A0A2S5U4P5: Host cell division inhibitor Icd-like protein",
  "UniRef90_UPI00038FCE9B: NO_NAME",
  "UniRef90_W1AXM3: NO_NAME",
  "UniRef90_UPI0007AC5EAB: DUF1327 domain-containing protein",
  "UniRef90_A0A2T3T734: DUF1133 domain-containing protein",
  "UniRef90_A0A2W8A8J7: DUF1391 domain-containing protein",
  "UniRef90_H4IRI1: Phage portal, lambda family protein",
  "UniRef90_Q9XB21: DNA-binding protein HU",
  "UniRef90_A0A0K4C257: AlpA family phage regulatory protein",
  "UniRef90_A0A376LDT0: Prophage protein",
  "UniRef90_A0A166SSS0: Putative transposase",
  "UniRef90_UPI0009A49669: transposase",
  "UniRef90_UPI000BE7D43D: NO_NAME",
  "UniRef90_UPI00053BD19B: NO_NAME",
  "UniRef90_A0A1X3KAR2: Side tail fiber protein-like protein",
  "UniRef90_UPI0004D3E5F0: adhesin",
  "UniRef90_UPI00093396AE: NO_NAME",
  "UniRef90_A0A2X7VXZ3: H-NS histone family protein",
  "UniRef90_A0A2S5C5C4: NO_NAME",
  "UniRef90_A0A2Y8JU28: NO_NAME",
  "UniRef90_A0A2X8GC90: NO_NAME",
  "UniRef90_A0A377NPU1: Heat resistant agglutinin 1",
  "UniRef90_A0A1Q4PMN4: NO_NAME",
  "UniRef90_A0A376T2Z4: Invasin",
  "UniRef90_A0A376STP5: Pantoate--beta-alanine ligase",
  "UniRef90_A0A148I0J9: Transcriptional regulator",
  "UniRef90_A0A2A3V9V5: NO_NAME",
  "UniRef90_A0A0K4NTQ7: Putative partitioning protein B",
  "UniRef90_A0A1M0EA82: DNA topoisomerase III",
  "UniRef90_A0A1M0S5W4: Integrating conjugative element protein",
  "UniRef90_A0A1Q4PBI5: DNA primase",
  "UniRef90_A0A1M0EAY4: Integrating conjugative element membrane protein",
  "UniRef90_A0A0C2B619: SciB domain protein",
  "UniRef90_A0A1D3U3Z4: AraC family transcriptional regulator",
  "UniRef90_A0A2T1LC15: Helicase SNF2",
  "UniRef90_A0A2T1LC21: ATP-dependent DNA helicase RecQ",
  "UniRef90_A0A2J7KWF5: NO_NAME",
  "UniRef90_A0A2Y0T492: IS1414, transposase",
  "UniRef90_V0Z598: NO_NAME",
  "UniRef90_V0YGQ9: NO_NAME",
  "UniRef90_A0A0K3GQT4: Tellurite resistance protein",
  "UniRef90_A0A0F3V8E2: NO_NAME",
  "UniRef90_A0A2K3TUS6: NO_NAME",
  "UniRef90_A0A2X7RI44: Putative transposase",
  "UniRef90_A0A2J7KL40: IS66 Orf2 family protein",
  "UniRef90_A0A1D7PI01: P pili regulatory PapB protein",
  "UniRef90_A0A372HLP4: NO_NAME",
  "UniRef90_A0A0F3V941: Type I restriction endonuclease"
)

## green
genes_del_v3 <- c(
  "Q99XQ7: Elongation factor Ts",
  "A0A376WAT6: NO_NAME",
  "A0A2P9E689: Macrophage stimulating factor",
  
  "UPI000C0BD5DD: NO_NAME",
  
  "A0A209NSH2: Phage tail protein",
  "A0A0V9S860: NO_NAME",
  "A0A2P6IED0: Transcriptional regulator",
  "V0VF20: NO_NAME"
)
(genes_del_v3 <- 
    paste0("UniRef90_",genes_del_v3))

## blue rectangle
genes_del_v4 <- c(
  "A0A377D3Q1: Transient receptor potential locus",
  "A0A376K8W6: Domain of uncharacterized function (DUF1972)",
  "A0A376VBV7: Putative cytoplasmic protein",
  "A0A376MP41: YD repeat-containing protein",
  "A0A377DDJ9: RhsA",
  "D7ZHQ3: RHS repeat-associated core domain protein",
  "A0A377LHH4: VgrG protein, Encoded within repeats that are hotspots for chromosomal duplication formation Function of protein is uncharacterized",
  "A0A0G3JYI9: NO_NAME",
  "A0A377L7U4: RhsC element core protein RshC",
  "A0A376Q238: Magnesium chelatase",
  "A0A377LKN9: DNA helicase II",
  "S1IFT9: NO_NAME",
  "A0A2S8JY30: NO_NAME",
  "UPI0009890E1A: transposase",
  "A0A377LF13: Multidrug resistance protein",
  "A0A377DC91: Nitrogen regulation protein NR(I)"
)
genes_del_v4 <- 
  paste0("UniRef90_",genes_del_v4)
head(genes_del_v4)

## red rectangle
genes_del_v5 <- c(
  "K4V9P3: NO_NAME",
  "UPI000645980F: NO_NAME",
  "A0A368J033: Lysis protein (Fragment)",
  "A0A376VSW5: Protein of uncharacterized function (DUF826)",
  "A0A029JN07: NO_NAME",
  "Q93D68: PilL",
  "A0A376VT02: Minor pilin subunit",
  "A0A029J8I8: NO_NAME",
  "A0A0Q3B8W5: Protein of uncharacterized function (DUF1378)",
  "A0A0A0FR20: Late control gene D protein",
  "A0A0A0GQB7: NO_NAME",
  "A0A377D577: Putative tail protein from prophage putative tail length tape measure motif",
  "A0A377D6D4: NO_NAME",
  "A0A377D4K3: Peptidase S6, IgA endopeptidase from phage origin",
  "A0A377CX27: Peptidase S6, IgA endopeptidase from phage origin",
  "A0A376W2F5: Phage protein",
  "A0A377D3I9: NO_NAME",
  "A0A0Q3AAK3: NO_NAME",
  "A0A376W4X8: Putative endonuclease from prophage, replication protein A (GpA)",
  "A0A029JTM3: NO_NAME",
  "A0A2Y0ZKN9: Putative colicin lysis protein",
  "A0A193LS31: NO_NAME",
  "A0A0H0RKG4: NO_NAME",
  "A0A376W625: YjhS",
  "A0A376DKT9: NO_NAME",
  "A0A2X1PJU1: Phage protein YjhS encoded within prophage CP-933O",
  "A0A029JTZ1: NO_NAME",
  "A0A376MG55: Peptidase S6, IgA endopeptidase from phage origin",
  "A0A376DHT2: Peptidase S6, IgA endopeptidase from phage origin",
  "UPI000C7CE0AE: autotransporter outer membrane beta-barrel domain-containing protein",
  "A0A2X1NGZ7: YjhS",
  "A0A376W208: Peptidase S6, IgA endopeptidase from phage origin",
  "A0A0A0GN36: NO_NAME",
  "A0A376DN51: NO_NAME",
  "A0A2T3SZC8: Shufflon system plasmid conjugative transfer pilus tip adhesin PilV",
  "UPI0007A0A40D: NO_NAME",
  "A0A376W4H1: Putative portal protein",
  "A0A160RGS4: NO_NAME",
  "A0A0Q3BKP6: NO_NAME",
  "A0A368J247: NO_NAME",
  "A0A0Q3BCB8: NO_NAME",
  "A0A376W7M3: Phage anti-repressor protein",
  "A0A2W8A8J7: DUF1391 domain-containing protein",
  "A0A368J2Z0: NO_NAME",
  "A0A376M5V3: Minor tail protein G",
  "A0A085P3E8: Tail protein",
  "I7AQ17: HlyD secretion family protein",
  "UPI00098AF521: DNA-packaging protein",
  "A0A0A0GY59: NO_NAME",
  "L4UQ03: NO_NAME",
  "H4UTT5: Replication initiation protein",
  "UPI000452D5D2: NO_NAME",
  
  "UPI000CF4DC9E: NO_NAME",
  "A0A376TA49: Phage protein",
  "A0A376LYR5: Phage protein",
  "UPI00092F829A: resolvase",
  "A0A376WX64: Bacteriophage protein",
  "UPI0006A5559C: phage tail protein",
  "A0A0K4C548: NO_NAME",
  "UPI0003EFA011: carbohydrate-binding protein",
  "A0A383HQC1: Tail component encoded by prophage CP-933N",
  "A0A376W6R6: Tail component of prophage CP-933X",
  "UPI000C7DC529: NO_NAME",
  "A0A2W7UCI6: Phage tail tape measure protein (Fragment)",
  "A0A368IY01: DUF1983 domain-containing protein (Fragment)",
  "UPI0007077D42: phage tail tape measure protein",
  "A0A025BQ58: Homoserine acetyltransferase (Fragment)",
  "A0A376DKX3: Host specificity protein J from prophage",
  "A0A376W9W1: Host specificity protein J from prophage",
  "Q9LA62: ORF-401-like protein",
  "A0A2X1NGS8: Putative tail fiber protein from prophage",
  "A0A376GDS0: Putative tail fiber protein from prophage",
  "A0A377CWJ9: Host specificity protein J from prophage",
  "A0A376D6L9: Putative tail fiber protein from prophage",
  "L3BHT2: NO_NAME",
  "UPI0005CDE880: NO_NAME",
  "T9U2J4: NO_NAME",
  "A0A376W3N7: Tail fiber protein",
  "A0A2T3T734: DUF1133 domain-containing protein",
  "A0A376VWX3: NO_NAME",
  "A0A377D377: Putative tail fiber protein from prophage",
  "A0A0D6ZY29: DUF1133 domain-containing protein"
)

genes_del_v5 <- 
  paste0("UniRef90_",genes_del_v5)
head(genes_del_v5)

gene_pre_filter <- filter(gene_pre_t2d_con,diff_perc<50)
dim(gene_pre_filter) #74
head(genes_del_v6 <- gene_pre_filter$gene)

genes_del_v7 <- 
  c("A0A0F6F402: NO_NAME",
    "A0A2B7MJL5: Protein ren",
    "S1EJF9: NO_NAME",
    "A0A345ETE5: PTS galactitol transporter subunit IIB",
    "A0A2X7QWU6: YD repeat-containing protein",
    "A0A377K9C5: Protein of uncharacterized function (DUF2732)",
    "UPI000DFF11E8: integrase",
    "UPI000BA96172: integrase",
    "A0A2H9FGR3: Integrase",
    "A0A385G0I1: NO_NAME",
    "UPI000BA30171: integrase",
    "A0A152BWF2: NO_NAME",
    "E1IZE6: Aerobactin siderophore biosynthesis protein IucB family protein (Fragment)",
    "I2RBJ4: NO_NAME",
    "W0AT18: NO_NAME",
    "A0A376RIC0: Putative aldolase")
genes_del_v7 <- 
  paste0("UniRef90_",genes_del_v7)

length(genes_del_v1) #197
length(genes_del_v2) #113
length(genes_del_v3) #7
length(genes_del_v4) #16
length(genes_del_v5) #82
length(genes_del_v6) #74
length(genes_del_v7) #17


table(genes_del_v1 %in% anpan_re_sig$gene)
genes_del_v1[!genes_del_v1 %in% anpan_re_sig$gene]

table(genes_del_v2 %in% anpan_re_sig$gene)
genes_del_v2[!genes_del_v2 %in% anpan_re_sig$gene]

table(genes_del_v3 %in% anpan_re_sig$gene)
genes_del_v3[!genes_del_v3 %in% anpan_re_sig$gene]

table(genes_del_v4 %in% anpan_re_sig$gene)
genes_del_v4[!genes_del_v4 %in% anpan_re_sig$gene]

table(genes_del_v5 %in% anpan_re_sig$gene)
genes_del_v5[!genes_del_v5 %in% anpan_re_sig$gene]

table(genes_del_v6 %in% anpan_re_sig$gene)
table(genes_del_v7 %in% anpan_re_sig$gene)

genes_del1 <- c(genes_del_v1,
                genes_del_v2,
                genes_del_v3,
                genes_del_v4
)
length(genes_del1) #330

genes_del2 <- c(genes_del_v1,
                genes_del_v2,
                genes_del_v3,
                genes_del_v4,
                genes_del_v5
)
length(genes_del2) #412

genes_del3 <- c(genes_del_v1,
                genes_del_v2,
                genes_del_v3,
                # genes_del_v4,
                genes_del_v5
)
length(genes_del3) #396

# remove red and diff_pre
genes_del4 <- unique(c(
  genes_del_v1,
  genes_del_v2,
  genes_del_v3,
  # genes_del_v4,
  genes_del_v5,
  genes_del_v6
))
length(genes_del4) #404

# remove blue and diff_pre
genes_del5 <- unique(c(
  genes_del_v1,
  genes_del_v2,
  genes_del_v3,
  genes_del_v4,
  # genes_del_v5,
  genes_del_v6
))
length(genes_del5) #342

# remove blue, diff_pre, and inconsistant block
genes_del6 <- c(genes_del_v1,
                genes_del_v2,
                genes_del_v3,
                genes_del_v5,
                genes_del_v7
)
length(genes_del6) #412

#### filter samples -------
anpan_re_sig <- filter(anpan_re_sig,!gene %in% genes_del_v0)
anpan_re_sig2 <- filter(anpan_re_sig, !gene %in% genes_del3) # main results
dim(anpan_re_sig2) 


# based on sig genes; with 2% genes or less (or 90% genes or more?)
model_sig <- model_input[,c("study","sex","age","bmi","sample_id","T2D",
                            anpan_re_sig2$gene)]
dim(model_sig) #916/416

model_sig2 <- model_sig %>% 
  filter((rowSums(model_sig[,7:ncol(model_sig)]) >= 
            (ncol(model_sig)-6)*perct_samp))
dim(model_sig2) #592/419

#### filter genes ------
## based on 
model_input3 <- model_input %>%
  filter(sample_id %in% model_sig2$sample_id) %>% as.data.table()

dim(model_input3) #603/82568
head(colnames(model_input3),10)
colnames(model_input3) <- gsub("UniRef90_","",colnames(model_input3))
table(model_input3$`Q1RFY5: NO_NAME`)

## gene model results ------
anpan_re <- filter(anpan_re_all,bug_name == bugName)
head(anpan_re);dim(anpan_re) #105510/11
anpan_re$gene <- gsub("UniRef90_","",anpan_re$gene)
anpan_re <- anpan_re %>% 
  separate(gene,c("gene_id","gene_name"),": ",remove = F)

## beta * -log10(p)
anpan_re$beta_p <- anpan_re$estimate * (-log10(anpan_re$p.value))

genes_del0_name <- gsub("UniRef90_","",genes_del_v0)
genes_del1_name <- gsub("UniRef90_","",genes_del1)
genes_del2_name <- gsub("UniRef90_","",genes_del2)
genes_del3_name <- gsub("UniRef90_","",genes_del3)
genes_del4_name <- gsub("UniRef90_","",genes_del4)
genes_del5_name <- gsub("UniRef90_","",genes_del5)
genes_del6_name <- gsub("UniRef90_","",genes_del6)
head(genes_diff_name <- gsub("UniRef90_","",gene_pre_t2d_con_diff100$gene))

anpan_re2 <- anpan_re %>% 
  filter(!gene %in% genes_del0_name) %>%
  filter(!gene %in% genes_del3_name) %>% #main
  filter(gene %in% colnames(model_input3))
head(anpan_re2);dim(anpan_re2) #85508/13; 78201/13; 105313/13; 105200/13

## heatmap ------
q_threshold = 1
(beta_threshold = beta_thresh)
head(anpan_re2$q_bug_wise)
head(anpan_re2$q_global)
res <- anpan_re2[order(metarank_global)]
gene_level_df = res
dim(gene_level_df) #105118/13

gene_level_df = gene_level_df[gene_level_df$q_global < q_threshold]
gene_level_df = gene_level_df[abs(gene_level_df$estimate) >= beta_threshold]
gene_levels = gene_level_df$gene
dim(gene_level_df) #218/13

### Match categories --------
gene_cat <- read.xlsx("genes_heatmap_ddw2.xlsx")
gene_level_df <- gene_level_df %>% 
  inner_join(gene_cat,by="gene_id")

### Get the order of genes ------
gene_mat = model_input3 |>
  dplyr::select('sample_id', all_of(gene_levels)) |>
  tibble::column_to_rownames("sample_id") |>
  as.matrix()
dim(gene_mat) #592/218
head(gene_mat[1:3,1:3])
gene_mat <- 1*gene_mat #convert to numeric

g_clust = hclust(dist(t(gene_mat)))
gene_levels = colnames(gene_mat)[g_clust$order]
length(gene_levels)
head(gene_levels,10)

### remove genes keeping orders -----
gene_levels2 <-
  gene_levels[!gene_levels %in% gsub("UniRef90_","",genes_del_v7)]
length(gene_levels2)

gene_mat2 = model_input3 |>
  dplyr::select('sample_id', all_of(gene_levels2)) |>
  tibble::column_to_rownames("sample_id") |>
  as.matrix()
head(gene_mat2[1:3,1:3])
gene_mat2 <- 1*gene_mat2
dim(gene_mat2) #592/202

### Get the order of samples --------
covariates = c("study") #,"age","bmi","sex"
outcome = "T2D"
select_cols = c("sample_id",covariates, outcome)
color_bar_df = model_input3 |>
  dplyr::select(dplyr::all_of(select_cols)) |>
  unique()
sample_clust = gene_mat |>
  pca()  |>
  dist() |>
  hclust()

tree_list = ctl_case_trees(sample_clust,
                           model_input3,
                           outcome)
ctl_tree  = tree_list[[1]]
case_tree = tree_list[[2]]
s_levels  = tree_list[[3]]

ctl_tree =  ctl_tree$result
case_tree = case_tree$result

color_bar_df$sample_id = 
  factor(color_bar_df$sample_id,levels = s_levels)

bug_covariate = "present"
fill_scale = scale_fill_manual(values = c("FALSE"  = "dodgerblue4", "TRUE"  = "chartreuse"))
guide_obj = guides(fill = guide_legend(title.position = 'bottom', title.hjust = .5))

### plot data -------
model_input_long = data.table::melt(model_input3 |> 
                                  dplyr::select(all_of(select_cols),
                                                # all_of(gene_levels)
                                                all_of(gene_levels2)
                                                ),
                                id.vars = c(covariates, outcome, "sample_id"),
                                variable.name = "gene",
                                value.name = bug_covariate)
plot_data = model_input_long |>
  mutate(gene = factor(gene, 
                       levels = gene_levels2
                       # levels = gene_levels
                       ),
         sample_id = factor(sample_id,
                            levels = levels(color_bar_df$sample_id)))

(n_healthy = sum(color_bar_df[[outcome]] ==0)) #255
(n_case = sum(color_bar_df[[outcome]] ==1)) #337



### 1.annotation plot -------
anno_plot = plot_color_bars(color_bar_df = color_bar_df,
                            model_input = model_input_long,
                            covariates = covariates,
                            offset = NULL,
                            outcome = outcome,
                            binary_outcome = TRUE)
print(anno_plot)

color_bar_i = 1 # the color bar subplot in the patchwork
study_i = 3 # the scale used for study on the plot
t2d_i = 1 

(aes_name = anno_plot$scales$scales[[study_i]]$aesthetics)
(guide = anno_plot$scales$scales[[study_i]]$guide)

guide$available_aes = aes_name
anno_plot$scales$scales[[study_i]] <- 
  scale_fill_manual(aesthetics = aes_name,
                    guide = guide,
                    values = study_col) 

(aes_name_t2d = anno_plot$scales$scales[[t2d_i]]$aesthetics)
(guide_t2d = anno_plot$scales$scales[[t2d_i]]$guide)
anno_plot$scales$scales[[1]] <- 
  scale_fill_manual(aesthetics = aes_name_t2d,
                    guide = guide_t2d,
                    values = c("#72b0da","#C75127"))
anno_plot

### 2.heatmap tile ---------
plot_data = as.data.table(anpan_re2)[plot_data, on = 'gene']
plot_data$sample_id = factor(plot_data$sample_id,
                             levels = levels(color_bar_df$sample_id))
plot_data$gene = factor(plot_data$gene,levels = rev(gene_levels2))
lab_df = plot_data[,.(gene)] |> unique()
y_scale = scale_y_discrete(breaks = lab_df$gene,
                           position = 'left')
ns = dplyr::n_distinct(plot_data$sample_id)
ng = length(gene_levels2)
glab_frac = ifelse(ng > 50, (50/ng)^.62, 1)
continuous_genes = dplyr::n_distinct(gene_mat2 |> as.vector()) > 2
line_color = ifelse(continuous_genes, "grey", "black")
black_vline = geom_vline(lwd = .5,
                         color = line_color,
                         xintercept = n_healthy + .5)
heatmap_tile = plot_data |>
  mutate(present = as.logical(present)) |>
  ggplot(aes(y = gene, x = sample_id)) +
  geom_tile(aes(fill = present)) +
  black_vline +
  fill_scale +
  guide_obj +
  y_scale +
  labs(x = "samples",
       y = NULL) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        # axis.text.y = element_text(size = ggplot2::rel(glab_frac)),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.border = element_blank()) +
  coord_cartesian(expand = FALSE)
print(heatmap_tile)

### 3.gene category -----
## add gene category
dim(plot_data) #119584/29
plot_data <- plot_data %>% 
  left_join(gene_cat,by="gene_id")

table(plot_data$Category)
length(unique(plot_data$Category)) #10
unique(plot_data$Category)
paletteer::paletteer_d("ggsci::default_igv")
(col_panel <- c(paletteer::paletteer_d("ggsci::default_igv")[1:14],
                "grey70","grey90"))
scales::show_col(col_panel)
names(col_panel) <- c(
  "Amino acid metabolism",
  "Bacterial structural components",
  "Cell motility",
  # "Damaged DNA repair",
  "DNA replication & transcription",
  "Fatty acid metabolism",
  "Genetic Rearrangement", #E.coli heatmap
  "Glucose metabolism",
  "Phage and HGT",
  "Proteolysis", #E.coli heatmap
  "Quorum sensing",
  "Signal peptide processing",
  "Signal transduction",
  "Stress response",
  "Virulence and antibiotic resistance",
  "Other", #E.coli heatmap
  "Unknown" #E.coli heatmap
)

p_cat <- ggplot(plot_data,aes(x=1,y=gene))+
  geom_tile(aes(fill=Category),width=1)+ #color=category,
  labs(fill="Category")+ #,color="Category"
  scale_fill_manual(values=col_panel) +
  coord_cartesian(expand = FALSE)+
  # scale_x_discrete(breaks=breaks_go,labels=labels_go)+
  theme_void()+
  theme(#panel.grid = element_blank(),
        #axis.text.x = element_blank(),
        #axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = ggplot2::rel(glab_frac),hjust = 1),
        #axis.ticks.y = element_blank(),
        #axis.title = element_blank(),
        #plot.margin = unit(c(5.5,5.5,5.5,5.5),"point"),
        legend.position = "none")
print(p_cat)

## legend
p_cat_lgd <- ggplot(plot_data,aes(x=1,y=gene))+
  geom_tile(aes(fill=Category),width=1)+ #color=category,
  labs(fill="Category")+ #,color="Category"
  scale_fill_manual(values=col_panel) +
  coord_cartesian(expand = FALSE)+
  # scale_x_discrete(breaks=breaks_go,labels=labels_go)+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title = element_blank(),
        plot.margin = unit(c(5.5,5.5,5.5,5.5),"point"),
        legend.position = "right")
print(p_cat_lgd)

gene_cat_legend <- cowplot::get_legend(p_cat_lgd)
pdf(paste0("anpan_1_by_1/gene_category_legend_",date_log,".pdf"), 
    width = 5, height = 5, onefile = F) # Open a new pdf file
grid::grid.draw(gene_cat_legend)
dev.off()

### 4.p value bar -----
head(plot_data)
library(RColorBrewer)
breaklist <- seq(-1,1,by=0.001)
red_blue <- rev(brewer.pal(n=11,name="RdBu"))
scales::show_col(red_blue)
col_pal <- colorRampPalette(red_blue)(length(breaklist))

p_bar <- ggplot(plot_data,aes(x=1,y=gene))+
  geom_tile(aes(fill=beta_p),width=1)+ #color=category,
  labs(fill="-log10(P)*Beta")+ #,color="Category"
  scale_fill_gradientn(colors = (col_pal)) + #rev
  coord_cartesian(expand = FALSE)+
  theme_void()+
  theme(legend.position = "none")
print(p_bar)

## legend
p_bar_lgd <- ggplot(plot_data,aes(x=1,y=gene))+
  geom_tile(aes(fill=beta_p),width=1)+ #color=category,
  labs(fill="-log10(P)*Beta")+ #,color="Category"
  scale_fill_gradientn(colors = (col_pal)) + #rev
  coord_cartesian(expand = FALSE)+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title = element_blank(),
        plot.margin = unit(c(5.5,5.5,5.5,5.5),"point"),
        legend.position = "right")
print(p_bar_lgd)

p_bar_legend <- cowplot::get_legend(p_bar_lgd)
pdf(paste0("anpan_1_by_1/p_bar_legend_",date_log,".pdf"), 
    width = 5, height = 5, onefile = F) # Open a new pdf file
grid::grid.draw(p_bar_legend)
dev.off()

### *combine plots ---------
## show trees 
n = n_healthy + n_case
ctl_width =  5 * n_healthy/n
case_width = 5 * n_case/n
tree_plot = patchwork::wrap_plots(ctl_tree, case_tree) +
  patchwork::plot_layout(nrow = 1, widths = c(ctl_width, case_width))
(title_str = paste("E.coli", " (n = ", ns, ")", sep = "", collapse = ""))
subtitle_str = paste0(length(gene_levels), " genes with Q below ", q_threshold, 
                      " and abs(coefficient) above ", beta_threshold)

threshold_warning_string = "Abundance color shown on log scale."
design_str = "
    #AAAAAAAAAAAAAAAAAAAAAAAAAA#
    DBBBBBBBBBBBBBBBBBBBBBBBBBBC
    DBBBBBBBBBBBBBBBBBBBBBBBBBBC
    DBBBBBBBBBBBBBBBBBBBBBBBBBBC
    DBBBBBBBBBBBBBBBBBBBBBBBBBBC
    DBBBBBBBBBBBBBBBBBBBBBBBBBBC
    DBBBBBBBBBBBBBBBBBBBBBBBBBBC
    DBBBBBBBBBBBBBBBBBBBBBBBBBBC
    DBBBBBBBBBBBBBBBBBBBBBBBBBBC
    DBBBBBBBBBBBBBBBBBBBBBBBBBBC
    DBBBBBBBBBBBBBBBBBBBBBBBBBBC
    "
print(heatmap_tile)
p = patchwork::wrap_plots(anno_plot, heatmap_tile,p_bar,p_cat,
                          ncol = 3,
                          guides = 'collect',
                          design = design_str
) +
  patchwork::plot_annotation(title = title_str,
                             # caption = threshold_warning_string,
                             subtitle = subtitle_str)
print(p)

pdf(paste0("anpan_1_by_1/heatmap_ecoli_manual_",
           date_log, ".pdf"),
    width = 15,height = 5) # height = 20/8
p
dev.off()

# pdf(paste0("anpan_1_by_1/heatmap_ecoli_manual_long_",
#            date_log, ".pdf"),
#     width = 15,height = 12) # height = 20/8
# p
# dev.off()
