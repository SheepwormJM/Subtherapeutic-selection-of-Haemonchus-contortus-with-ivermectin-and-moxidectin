# Diversity along the genome, calculated by GrENEDALF
# Using Grenedalf v0.3.0
# Using the ss57 file made with the Selection Line samples only

https://github.com/lczech/grenedalf/wiki/Subcommand:-diversity

Will be using the corrected calcuations, based on the original PoPoolation calculations

Will look at Pi and Tajima's D (The later might be of greater interest in the selection line experiment than in most field samples).

1. Calculate diversity measures by 10kb windows

```
module load grenedalf/0.3.0

grenedalf diversity \
--sync-path /home/jenni/data_folder/working_folder/selection_lines_poolseq/mpileup/Jan24_correct_order_syncs/SLps_only_co_Hcwbps18_q20Q30_noindels_ss57.java.sync \
--reference-genome-fasta-file /home/jenni/data_folder/wbps18/haemonchus_contortus.PRJEB506.WBPS18.genomic.fa \
--rename-samples-file /home/jenni/data_folder/working_folder/selection_lines_poolseq/pairwise_comparisons/LINES_IN_CORRECT_ORDER/grenedalf/Jan24_F0_vs_all_sel_lines_only/SL_ss57_sample_names.txt \
--filter-sample-min-count 2 \
--filter-sample-min-coverage 57 \
--filter-sample-max-coverage 57 \
--window-type sliding --window-sliding-width 10000 --window-sliding-stride 10000 \
--pool-sizes 400 \
--separator-char tab \
--out-dir /home/jenni/data_folder/working_folder/selection_lines_poolseq/pairwise_comparisons/LINES_IN_CORRECT_ORDER/grenedalf/Jan24_F0_vs_all_sel_lines_only/ss57/diversity \
--file-prefix 400_Gf_F0_vs_SL_only_Sel_lines_Hcwbps18_q20Q30_noindels_ss57_10kb_diversity \
--log-file grenedalf_10kb_diversity_ss57.log \
2>error_log_10kb_diversity &
```

2. Melt the table:
```
module load R/4.3.1

R

div<-read.table("400_Gf_F0_vs_SL_only_Sel_lines_Hcwbps18_q20Q30_noindels_ss57_10kb_diversitydiversity.csv", header=TRUE)
summary(div)

div$mid<-((div$start + div$end)/2)

library(reshape2)
div_melt<-(melt(div, id.vars=c("chrom","start","end","mid")))
head(div_melt)



write.table(div_melt, file="400_Gf_F0_vs_SL_only_Sel_lines_Hcwbps18_q20Q30_noindels_ss57_10kb_diversity_MELTED.txt", row.names=FALSE, sep="\t", quote= FALSE) 
```
Correlation is stronger between snps and theta_pi_rel than for snps and fst: 
```
> cor(div$F0.snp_count, div$F0.theta_pi_rel)
[1] 0.5200153
> cor(div$F0.snp_count, div$F3_M_L3.theta_pi_rel)
[1] 0.50509
> cor(div$F3_M_L3.snp_count, div$F3_M_L3.theta_pi_rel)
[1] 0.5473724
```
And for Tajima's D:
```
> cor(div$F3_M_L3.snp_count, div$F3_M_L3.tajimas_d)
[1] 0.2200548
> cor(div$F3_M_L3.snp_count, div$F3_M_L3.tajimas_d)
[1] 0.2200548
> 
> cor(div$F3_M_L3.theta_pi_rel, div$F3_M_L3.tajimas_d)
[1] 0.5711798
```

# Pi

Select just the Pi:

Note that there are two of these: theta_pi_abs and theta_pi_rel. From (https://www.biorxiv.org/content/biorxiv/early/2022/02/04/2022.02.02.477408.full.pdf) I think that they must represent the sample and an estimate of the population, but not certain. Theta_pi_rel is the one that I want to use.

## theta_pi_rel
```
# Filter the data to retain each LINE separately:

library(stringr)
library(dplyr)

div_melt_pi<- div_melt %>%
 filter(str_detect(variable, "theta_pi_rel")) 

div_melt_L2<- div_melt_pi %>%
  filter(str_detect(variable, "L2"))

div_melt_L3<- div_melt_pi %>%
  filter(str_detect(variable, "L3"))

div_melt_F0<- div_melt_pi %>%
  filter(!str_detect(variable, "F0"))


# Only want the following for those with IVM for F1 (as no F1 line 1 ivm sample)
# So, instead do it this way:

div_L1_pirel<-div[,c(1:3,14,35,56,91,112,133,154,175,193)]
names(div_L1_pirel)

# Now to melt the table. Need package reshape2

library(reshape2)
win_melt_L1 <- melt(div_L1_pirel, id.vars=c("chrom","start","end","mid"))
head(win_melt_L1)

# Then, add in the extra facet:
win_melt_L1$variable.new <- factor(win_melt_L1$variable, levels=(c("F1_C_L1.theta_pi_rel", "F2_C_L1.theta_pi_rel", "F3_C_L1.theta_pi_rel", "", "F2_I_L1.theta_pi_rel", "F3_I_L1.theta_pi_rel", "F1_M_L1.theta_pi_rel", "F2_M_L1.theta_pi_rel", "F3_M_L1.theta_pi_rel")))   


```
```
# Plot by line, faceted by generation and anthelmintic

library(ggplot2)
library(patchwork)

x<- ggplot(data=subset(div_melt_F0, div_melt_F0$chrom == "hcontortus_chr5_Celeg_TT_arrow_pilon"), mapping=aes(x=mid/1000000, y=log10(value))) +
geom_point(colour="#619CFF", alpha=0.3) +
labs(y="log10(Theta_pi_rel)", x="Basepair (Mb)")+
theme(axis.title.x = element_text(size=14), axis.text.x = element_text(size=12), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12), strip.text.x = element_text(size = 14))+
ggtitle('F0 MHco3(ISE)')

x2<- ggplot(data=subset(win_melt_L1, win_melt_L1$chrom == "hcontortus_chr5_Celeg_TT_arrow_pilon"), mapping=aes(x=mid/1000000, y=log10(value))) +
geom_point(colour="#619CFF", alpha=0.3) +
labs(y="log10(Theta_pi_rel)", x="Basepair (Mb)")+
theme(axis.title.x = element_text(size=14), axis.text.x = element_text(size=12), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12), strip.text.x = element_text(size = 14))+
facet_wrap(~variable.new, ncol=3, drop=FALSE, dir="v")+
ggtitle('Line 1 generations')

x3<- ggplot(data=subset(div_melt_L2, div_melt_L2$chrom == "hcontortus_chr5_Celeg_TT_arrow_pilon"), mapping=aes(x=mid/1000000, y=log10(value))) +
geom_point(colour="#619CFF", alpha=0.3) +
labs(y="log10(Theta_pi_rel)", x="Basepair (Mb)")+
theme(axis.title.x = element_text(size=14), axis.text.x = element_text(size=12), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12), strip.text.x = element_text(size = 14))+
facet_wrap(~variable, ncol=3, drop=TRUE, dir="v")+
ggtitle('Line 2 generations')


x4<- ggplot(data=subset(div_melt_L3, div_melt_L3$chrom == "hcontortus_chr5_Celeg_TT_arrow_pilon"), mapping=aes(x=mid/1000000, y=log10(value))) +
geom_point(colour="#619CFF", alpha=0.3) +
labs(y="log10(Theta_pi_rel)", x="Basepair (Mb)")+
theme(axis.title.x = element_text(size=14), axis.text.x = element_text(size=12), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12), strip.text.x = element_text(size = 14))+
facet_wrap(~variable, ncol=3, drop=TRUE, dir="v")+
ggtitle('Line 3 generations')


pdf("400_Gf_SLps_Hcwbps18_q20Q30_ss57_sellinesonly.10kb_theta_pi_rel_feb2024.pdf",width = 16,height = 12,useDingbats = FALSE)
patchwork <- x + x2 + x3 + x4 + plot_layout(guides = 'collect')

patchwork + plot_annotation(
  title = 'Haemonchus contortus Selection Lines experiment',
   caption = 'Subsampled depth = 57, Grenedalf v0.3.0, using corrected PoPoolation2 calculations',
tag_levels = 'A'
)
dev.off()

# And zoom to the 'region' under selection...

y<-x+
coord_cartesian(xlim=c(27.5, 48.868368))

y2<-x2 +
coord_cartesian(xlim=c(27.5, 48.868368))


y3<-x3+
coord_cartesian(xlim=c(27.5, 48.868368))


y4<-x4+
coord_cartesian(xlim=c(27.5, 48.868368))



pdf("chrVlocus_400_Gf_SLps_Hcwbps18_q20Q30_ss57_sellinesonly.10kb_theta_pi_rel_feb2024.pdf",width = 16,height = 12,useDingbats = FALSE)
patchwork <- y + y2 + y3 + y4 + plot_layout(guides = 'collect')

patchwork + plot_annotation(
  title = 'Haemonchus contortus Selection Lines experiment',
   caption = 'Subsampled depth = 57, Grenedalf v0.3.0, using corrected PoPoolation2 calculations',
tag_levels = 'A'
)
dev.off()

```
```
scp jenni@scamper.mvls.gla.ac.uk:/home/jenni/data_folder/working_folder/selection_lines_poolseq/pairwise_comparisons/LINES_IN_CORRECT_ORDER/grenedalf/Jan24_F0_vs_all_sel_lines_only/ss57/diversity/400_Gf_SLps_Hcwbps18_q20Q30_ss57_sellinesonly.10kb_theta_pi_rel_feb2024.pdf ./400_Gf_SLps_Hcwbps18_q20Q30_ss57_sellinesonly.10kb_theta_pi_rel_feb2024.pdf

scp jenni@scamper.mvls.gla.ac.uk:/home/jenni/data_folder/working_folder/selection_lines_poolseq/pairwise_comparisons/LINES_IN_CORRECT_ORDER/grenedalf/Jan24_F0_vs_all_sel_lines_only/ss57/diversity/chrVlocus_400_Gf_SLps_Hcwbps18_q20Q30_ss57_sellinesonly.10kb_theta_pi_rel_feb2024.pdf ./chrVlocus_400_Gf_SLps_Hcwbps18_q20Q30_ss57_sellinesonly.10kb_theta_pi_rel_feb2024.pdf

scp jenni@scamper.mvls.gla.ac.uk:/home/jenni/data_folder/working_folder/selection_lines_poolseq/pairwise_comparisons/LINES_IN_CORRECT_ORDER/grenedalf/Jan24_F0_vs_all_sel_lines_only/ss57/diversity/F0_chrVlocus_400_Gf_SLps_Hcwbps18_q20Q30_ss57_sellinesonly.10kb_theta_pi_rel_feb2024.pdf ./F0_chrVlocus_400_Gf_SLps_Hcwbps18_q20Q30_ss57_sellinesonly.10kb_theta_pi_rel_feb2024.pdf
# Note - not sure where the code is for this plot - can't see the title anywhere - and the pdf seems to be an empty plot. 
```
<img width="468" alt="image" src="https://github.com/user-attachments/assets/1f265827-5f9d-4002-beba-08641d1aa984">

<img width="468" alt="image" src="https://github.com/user-attachments/assets/2f642854-0852-4819-b309-5397409338e1">
