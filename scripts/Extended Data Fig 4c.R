load(file="/n/holylabs/LABS/dwang_lab/Lab/T2D_meta/3.analysis/R2/analysis/abund0020/metf_species/species_nometf_noinsul_ordinal.RData")

library(ggplot2)
dataplot <- data.frame(all=meta_3v2$coef, 
                       nometf=meta_nometf_results$coef)
cor(dataplot$all, dataplot$nometf, method="spearman") # 0.95

pdf("/n/holylabs/LABS/dwang_lab/Lab/T2D_meta/3.analysis/R2/figS6/fig_s6c_ordinal.pdf", 
    width = 3.8, height = 3.8, onefile = F) # Open a new pdf file
ggplot(aes(x=all, y=nometf), data=dataplot) +
  geom_point() +
  geom_hline(yintercept = 0, col="red") +
  geom_vline(xintercept = 0, col="red") +
  labs(x="Coefficients from analysis among\nall participants adjusting for metf",
       y="Coefficients from analysis among\nparticipants not using metf",
       title="") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid = element_blank(),
        axis.text = element_text(colour = "black", size=12),
        axis.title = element_text(colour = "black", size=12),
        panel.border = element_rect(colour = "black", size=1.2)) +
  # annotate(geom="text", x=-0.7, y=0.9, col="black", label="r=0.98")
  annotate(geom="text", x=-0.2, y=0.45, col="black", label="r=0.95")
dev.off() # Close the file
