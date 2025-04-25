#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)==1) {
  # default output file
  args[2] = "out.txt"
}

output_file = args[1]

# Trying on scamper - note that cannot install tidyverse, so installing the various used modules within it
library(dplyr)
library(ggplot2)
library(tidyr)
library(lme4)
library(stringr)
library(purrr)
library(broom.mixed)
library(patchwork)
library(viridis)
library(hexbin)

# Using R/4.4.2
# SessionInfo() would return:
# [1] hexbin_1.28.5       viridis_0.6.5       viridisLite_0.4.2  
# [4] patchwork_1.3.0     broom.mixed_0.2.9.6 purrr_1.0.2        
# [7] stringr_1.5.1       lme4_1.1-36         Matrix_1.7-1       
# [10] tidyr_1.3.1         ggplot2_3.5.1       dplyr_1.1.4       


load("IVM_vs_CTL_F0_results05Mar2025.obj")
load("MOX_vs_CTL_F0_results05Mar2025.obj")
load("SEL_vs_CTL_F0_results05Mar2025.obj")
load("IVM_vs_MOX_F0_results05Mar2025.obj")
load("IVM_vs_MOX_vs_CTL_F0_results05Mar2025.obj")


add_corrections <- function(df) { 
  df$pval_corr <- p.adjust(df$p,method="holm")
  df$log10pval_corr <- -1.0 * log10(df$pval_corr)
  df$is.significant <- df$pval_corr < 0.001
  df
}

new_p_values_IC <- add_corrections(p_values_IC)
new_p_values_MC <- add_corrections(p_values_MC)
new_p_values_SEL <- add_corrections(p_values_SEL)
new_p_values_IM <- add_corrections(p_values_IM)
new_p_values_IMC <- add_corrections(p_values_IMC)

# Calculate the Bonferroni corrected p-values, using an alpha value of 0.01

#p_values_MC_bon<-subset(p_values_MC, p_values_MC$neglogp >= 7.97234912)
#p_values_IC_bon<-subset(p_values_IC, p_values_IC$neglogp >= 7.97234912)
#p_values_SEL_bon<-subset(p_values_SEL, p_values_SEL$neglogp >= 7.97234912)
#p_values_IM_bon<-subset(p_values_IM, p_values_IM$neglogp >= 7.97234912)
#p_values_IMC_bon<-subset(p_values_IMC, p_values_IMC$neglogp >= 7.97234912)

MC_Holm_positions<-new_p_values_MC[,c("snp_ID","m1.coef","neglogp","pval_corr","log10pval_corr","is.significant")]
write.table(MC_Holm_positions, file="p_values_MC_Holm.txt", row.names=FALSE, sep="\t", quote= FALSE)

IC_Holm_positions<-new_p_values_IC[,c("snp_ID","m1.coef","neglogp","pval_corr","log10pval_corr","is.significant")]
write.table(IC_Holm_positions, file="p_values_IC_Holm.txt", row.names=FALSE, sep="\t", quote= FALSE)

SEL_Holm_positions<-new_p_values_SEL[,c("snp_ID","m1.coef","neglogp","pval_corr","log10pval_corr","is.significant")]
write.table(SEL_Holm_positions, file="p_values_SEL_Holm.txt", row.names=FALSE, sep="\t", quote= FALSE)

IM_Holm_positions<-new_p_values_IM[,c("snp_ID","m1.coef","neglogp","pval_corr","log10pval_corr","is.significant")]
write.table(IM_Holm_positions, file="p_values_IM_Holm.txt", row.names=FALSE, sep="\t", quote= FALSE)

IMC_Holm_positions<-new_p_values_IMC[,c("snp_ID","m1.coef","neglogp","pval_corr","log10pval_corr","is.significant")]
write.table(IMC_Holm_positions, file="p_values_IMC_Holm.txt", row.names=FALSE, sep="\t", quote= FALSE)



#plot P values_MC
MC<-ggplot(data=new_p_values_MC,aes(x=POS/1e6, y=m1.coef)) + 
labs(x="Genomic position (Mb)", y="Slope coefficient")+
geom_point(colour="grey", alpha=0.3, size=0.3)+
geom_point(data=subset(new_p_values_MC, new_p_values_MC$is.significant),  aes(x=POS/1e6, y=m1.coef), colour="red", alpha=0.5, size=0.3)+
geom_hline(yintercept=0, linetype=2, colour="white")+
scale_x_continuous(
    breaks = seq(0, 55, 5),
    minor_breaks = seq(1, 54, 1)
  )+
theme_light()+
theme(axis.title.x = element_text(size=10), axis.text.x = element_text(size=8), axis.title.y = element_text(size=10), axis.text.y = element_text(size=8), strip.text.x = element_text(size = 6, face = "bold", margin = margin(0.5,0,0.5,0, "mm")), strip.text.y = element_text(size = 6, face = "bold", angle=0), plot.title = element_text(size=10, hjust=0.5), legend.position = "top", legend.title = element_text(size=8), legend.text = element_text(size=6), legend.key.height=unit(0.2,"cm"), panel.spacing.x = unit(0.2, "mm"), panel.spacing.y = unit(0.75, "mm"))


#plot P values_IC
IC<-ggplot(data=new_p_values_IC,aes(x=POS/1e6, y=m1.coef)) + 
labs(x="Genomic position (Mb)", y="Slope coefficient")+
geom_point(colour="grey", alpha=0.3, size=0.3)+
geom_point(data=subset(new_p_values_IC, new_p_values_IC$is.significant),  aes(x=POS/1e6, y=m1.coef), colour="red", alpha=0.5, size=0.3)+
geom_hline(yintercept=0, linetype=2, colour="white")+
scale_x_continuous(
    breaks = seq(0, 55, 5),
    minor_breaks = seq(1, 54, 1)
  )+
theme_light()+
theme(axis.title.x = element_text(size=10), axis.text.x = element_text(size=8), axis.title.y = element_text(size=10), axis.text.y = element_text(size=8), strip.text.x = element_text(size = 6, face = "bold", margin = margin(0.5,0,0.5,0, "mm")), strip.text.y = element_text(size = 6, face = "bold", angle=0), plot.title = element_text(size=10, hjust=0.5), legend.position = "top", legend.title = element_text(size=8), legend.text = element_text(size=6), legend.key.height=unit(0.2,"cm"), panel.spacing.x = unit(0.2, "mm"), panel.spacing.y = unit(0.75, "mm"))

#plot P values_SEL
SEL<-ggplot(data=new_p_values_SEL,aes(x=POS/1e6, y=m1.coef)) + 
labs(x="Genomic position (Mb)", y="Slope coefficient")+
geom_point(colour="grey", alpha=0.3, size=0.3)+
geom_point(data=subset(new_p_values_SEL, new_p_values_SEL$is.significant), aes(x=POS/1e6, y=m1.coef), colour="red", alpha=0.5, size=0.3)+
geom_hline(yintercept=0, linetype=2, colour="white")+
scale_x_continuous(
    breaks = seq(0, 55, 5),
    minor_breaks = seq(1, 54, 1)
  )+
theme_light()+
theme(axis.title.x = element_text(size=10), axis.text.x = element_text(size=8), axis.title.y = element_text(size=10), axis.text.y = element_text(size=8), strip.text.x = element_text(size = 6, face = "bold", margin = margin(0.5,0,0.5,0, "mm")), strip.text.y = element_text(size = 6, face = "bold", angle=0), plot.title = element_text(size=10, hjust=0.5), legend.position = "top", legend.title = element_text(size=8), legend.text = element_text(size=6), legend.key.height=unit(0.2,"cm"), panel.spacing.x = unit(0.2, "mm"), panel.spacing.y = unit(0.75, "mm"))



#plot P values_IM
IM<-ggplot(data=new_p_values_IM,aes(x=POS/1e6, y=m1.coef)) + 
labs(x="Genomic position (Mb)", y="Slope coefficient")+
geom_point(colour="grey", alpha=0.3, size=0.3)+
geom_point(data=subset(new_p_values_IM, new_p_values_IM$is.significant), aes(x=POS/1e6, y=m1.coef), colour="red", alpha=0.5, size=0.3)+
geom_hline(yintercept=0, linetype=2, colour="white")+
scale_x_continuous(
    breaks = seq(0, 55, 5),
    minor_breaks = seq(1, 54, 1)
  )+
theme_light()+
theme(axis.title.x = element_text(size=10), axis.text.x = element_text(size=6), axis.title.y = element_text(size=10), axis.text.y = element_text(size=4), strip.text.x = element_text(size = 6, face = "bold", margin = margin(0.5,0,0.5,0, "mm")), strip.text.y = element_text(size = 6, face = "bold", angle=0), plot.title = element_text(size=10, hjust=0.5), legend.position = "top", legend.title = element_text(size=8), legend.text = element_text(size=6), legend.key.height=unit(0.2,"cm"), panel.spacing.x = unit(0.2, "mm"), panel.spacing.y = unit(0.75, "mm"))


#plot P values_IMC
IMC<-ggplot(data=new_p_values_IMC,aes(x=POS/1e6, y=m1.coef)) + 
labs(x="Genomic position (Mb)", y="Slope coefficient")+
geom_point(colour="grey", alpha=0.3, size=0.3)+
geom_point(data=subset(new_p_values_IMC, new_p_values_IMC$is.significant), aes(x=POS/1e6, y=m1.coef), colour="red", alpha=0.5, size=0.3)+
geom_hline(yintercept=0, linetype=2, colour="white")+
scale_x_continuous(
    breaks = seq(0, 55, 5),
    minor_breaks = seq(1, 54, 1)
  )+
theme_light()+
theme(axis.title.x = element_text(size=10), axis.text.x = element_text(size=6), axis.title.y = element_text(size=10), axis.text.y = element_text(size=4), strip.text.x = element_text(size = 6, face = "bold", margin = margin(0.5,0,0.5,0, "mm")), strip.text.y = element_text(size = 6, face = "bold", angle=0), plot.title = element_text(size=10, hjust=0.5), legend.position = "top", legend.title = element_text(size=8), legend.text = element_text(size=6), legend.key.height=unit(0.2,"cm"), panel.spacing.x = unit(0.2, "mm"), panel.spacing.y = unit(0.75, "mm"))


#patch<- IC / MC / SEL / IM / IMC + plot_annotation(tag_levels = 'A') + plot_layout (axis_titles = "collect")
#patch

patch<- IC / MC / SEL + plot_annotation(tag_levels = 'A') + plot_layout (axis_titles = "collect")
patch

ggsave(paste(output_file,"_Holm.tiff", sep=""), plot=patch, device=tiff, width=2250, height=3250, units="px", dpi=320) 
ggsave(paste(output_file,"_Holm.pdf", sep=""), plot=patch, device=pdf, width=2250, height=3250, units="px", dpi=320) 

ggsave("MOX_vs_CTL_F0toF3_Holm.tiff", plot=MC, device=tiff, width=2250, height=2250, units="px", dpi=320) 
ggsave("IVM_vs_CTL_F0toF3_Holm.tiff", plot=IC, device=tiff, width=2250, height=2250, units="px", dpi=320) 
ggsave("SEL_vs_CTL_F0toF3_Holm.tiff", plot=SEL, device=tiff, width=2250, height=2250, units="px", dpi=320) 



#plot P values_IC
ICmaha<-ggplot(data=new_p_values_IC,aes(x=POS/1e6, y=log10pval_corr)) + 
labs(x="Genomic position (Mb)", y="-log10(p-value)")+
geom_point(colour="grey", alpha=0.3, size=0.3)+
geom_point(data=subset(new_p_values_IC, new_p_values_IC$is.significant), aes(x=POS/1e6, y=log10pval_corr, colour=as.factor(m1.coef < 0)), alpha=0.5)+
geom_hline(yintercept=0, linetype=2, colour="white")+
scale_x_continuous(
    breaks = seq(0, 55, 5),
    minor_breaks = seq(1, 54, 1)
  )+
theme_light()+
theme(axis.title.x = element_text(size=10), axis.text.x = element_text(size=6), axis.title.y = element_text(size=10), axis.text.y = element_text(size=4), strip.text.x = element_text(size = 6, face = "bold", margin = margin(0.5,0,0.5,0, "mm")), strip.text.y = element_text(size = 6, face = "bold", angle=0), plot.title = element_text(size=10, hjust=0.5), legend.position = "top", legend.title = element_text(size=8), legend.text = element_text(size=6), legend.key.height=unit(0.2,"cm"), panel.spacing.x = unit(0.2, "mm"), panel.spacing.y = unit(0.75, "mm"))


#plot P values_MC
MCmaha<-ggplot(data=new_p_values_MC,aes(x=POS/1e6, y=log10pval_corr)) + 
labs(x="Genomic position (Mb)", y="-log10(p-value)")+
geom_point(colour="grey", alpha=0.3, size=0.3)+
geom_point(data=subset(new_p_values_MC, new_p_values_MC$is.significant), aes(x=POS/1e6, y=log10pval_corr, colour=as.factor(m1.coef < 0)), alpha=0.5)+
geom_hline(yintercept=0, linetype=2, colour="white")+
scale_x_continuous(
    breaks = seq(0, 55, 5),
    minor_breaks = seq(1, 54, 1)
  )+
theme_light()+
theme(axis.title.x = element_text(size=10), axis.text.x = element_text(size=6), axis.title.y = element_text(size=10), axis.text.y = element_text(size=4), strip.text.x = element_text(size = 6, face = "bold", margin = margin(0.5,0,0.5,0, "mm")), strip.text.y = element_text(size = 6, face = "bold", angle=0), plot.title = element_text(size=10, hjust=0.5), legend.position = "top", legend.title = element_text(size=8), legend.text = element_text(size=6), legend.key.height=unit(0.2,"cm"), panel.spacing.x = unit(0.2, "mm"), panel.spacing.y = unit(0.75, "mm"))

#plot P values_SEL
SELmaha<-ggplot(data=new_p_values_SEL,aes(x=POS/1e6, y=log10pval_corr)) + 
labs(x="Genomic position (Mb)", y="-log10(p-value)")+
geom_point(colour="grey", alpha=0.3, size=0.3)+
geom_point(data=subset(new_p_values_SEL, new_p_values_SEL$is.significant), aes(x=POS/1e6, y=log10pval_corr, colour=as.factor(m1.coef < 0)), alpha=0.5)+
geom_hline(yintercept=0, linetype=2, colour="white")+
scale_x_continuous(
    breaks = seq(0, 55, 5),
    minor_breaks = seq(1, 54, 1)
  )+
theme_light()+
theme(axis.title.x = element_text(size=10), axis.text.x = element_text(size=6), axis.title.y = element_text(size=10), axis.text.y = element_text(size=4), strip.text.x = element_text(size = 6, face = "bold", margin = margin(0.5,0,0.5,0, "mm")), strip.text.y = element_text(size = 6, face = "bold", angle=0), plot.title = element_text(size=10, hjust=0.5), legend.position = "top", legend.title = element_text(size=8), legend.text = element_text(size=6), legend.key.height=unit(0.2,"cm"), panel.spacing.x = unit(0.2, "mm"), panel.spacing.y = unit(0.75, "mm"))


#plot P values_IM
IMmaha<-ggplot(data=new_p_values_IM,aes(x=POS/1e6, y=log10pval_corr)) + 
labs(x="Genomic position (Mb)", y="-log10(p-value)")+
geom_point(colour="grey", alpha=0.3, size=0.3)+
geom_point(data=subset(new_p_values_IM, new_p_values_IM$is.significant), aes(x=POS/1e6, y=log10pval_corr, colour=as.factor(m1.coef < 0)), alpha=0.5)+
geom_hline(yintercept=0, linetype=2, colour="white")+
scale_x_continuous(
    breaks = seq(0, 55, 5),
    minor_breaks = seq(1, 54, 1)
  )+
theme_light()+
theme(axis.title.x = element_text(size=10), axis.text.x = element_text(size=6), axis.title.y = element_text(size=10), axis.text.y = element_text(size=4), strip.text.x = element_text(size = 6, face = "bold", margin = margin(0.5,0,0.5,0, "mm")), strip.text.y = element_text(size = 6, face = "bold", angle=0), plot.title = element_text(size=10, hjust=0.5), legend.position = "top", legend.title = element_text(size=8), legend.text = element_text(size=6), legend.key.height=unit(0.2,"cm"), panel.spacing.x = unit(0.2, "mm"), panel.spacing.y = unit(0.75, "mm"))


#plot P values_IMC
IMCmaha<-ggplot(data=new_p_values_IMC,aes(x=POS/1e6, y=log10pval_corr)) + 
labs(x="Genomic position (Mb)", y="-log10(p-value)")+
geom_point(colour="grey", alpha=0.3, size=0.3)+
geom_point(data=subset(new_p_values_IMC, new_p_values_IMC$is.significant), aes(x=POS/1e6, y=log10pval_corr, colour=as.factor(m1.coef < 0)), alpha=0.5)+
geom_hline(yintercept=0, linetype=2, colour="white")+
scale_x_continuous(
    breaks = seq(0, 55, 5),
    minor_breaks = seq(1, 54, 1)
  )+
theme_light()+
theme(axis.title.x = element_text(size=10), axis.text.x = element_text(size=6), axis.title.y = element_text(size=10), axis.text.y = element_text(size=4), strip.text.x = element_text(size = 6, face = "bold", margin = margin(0.5,0,0.5,0, "mm")), strip.text.y = element_text(size = 6, face = "bold", angle=0), plot.title = element_text(size=10, hjust=0.5), legend.position = "top", legend.title = element_text(size=8), legend.text = element_text(size=6), legend.key.height=unit(0.2,"cm"), panel.spacing.x = unit(0.2, "mm"), panel.spacing.y = unit(0.75, "mm"))



patch<- ICmaha / MCmaha / SELmaha / IMmaha / IMCmaha + plot_annotation(tag_levels = 'A') + plot_layout (axis_titles = "collect")
patch

ggsave(paste(output_file,"_Holm_Manhattan.tiff", sep=""), plot=patch, device=tiff, width=2250, height=2250, units="px", dpi=320) 

