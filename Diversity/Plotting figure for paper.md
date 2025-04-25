Colouring by rank, and plotting ChrV so that can easily compare between IVM and MOX lines.

# Diversity along the genome, calculated by GrENEDALF

_Using Grenedalf v0.3.0_

Using the ss57 file made with the Selection Line samples only

https://github.com/lczech/grenedalf/wiki/Subcommand:-diversity

Will be using the corrected calcuations, based on the original PoPoolation calculations

Will look at Pi and Tajima's D (The later might be of greater interest in the selection line experiment than in most field samples).

1. Calculate diversity measures by **100kb** windows

```
module load grenedalf/0.3.0

grenedalf diversity \
--sync-path /home/jenni/data_folder/working_folder/selection_lines_poolseq/mpileup/Jan24_correct_order_syncs/SLps_only_co_Hcwbps18_q20Q30_noindels_ss57.java.sync \
--reference-genome-fasta-file /home/jenni/data_folder/wbps18/haemonchus_contortus.PRJEB506.WBPS18.genomic.fa \
--rename-samples-file /home/jenni/data_folder/working_folder/selection_lines_poolseq/pairwise_comparisons/LINES_IN_CORRECT_ORDER/grenedalf/Jan24_F0_vs_all_sel_lines_only/SL_ss57_sample_names.txt \
--filter-sample-min-count 2 \
--filter-sample-min-coverage 57 \
--filter-sample-max-coverage 57 \
--window-type sliding --window-sliding-width 100000 --window-sliding-stride 100000 \
--pool-sizes 400 \
--separator-char tab \
--out-dir /home/jenni/data_folder/working_folder/selection_lines_poolseq/pairwise_comparisons/LINES_IN_CORRECT_ORDER/grenedalf/Jan24_F0_vs_all_sel_lines_only/ss57/diversity \
--file-prefix 400_Gf_F0_vs_SL_only_Sel_lines_Hcwbps18_q20Q30_noindels_ss57_100kb_diversity \
--log-file grenedalf_10kb_diversity_ss57.log \
2>error_log_10kb_diversity &
```

Rename the chr so are shorter for plotting:
```
sed 's/hcontortus_chr//g' 400_Gf_F0_vs_SL_only_Sel_lines_Hcwbps18_q20Q30_noindels_ss57_100kb_diversitydiversity.csv | sed 's/_Celeg_TT_arrow_pilon//g' > 400_Gf_F0_vs_SL_only_Sel_lines_Hcwbps18_q20Q30_noindels_ss57_100kb.renamed.diversity.csv 
```

2. Melt the table:
```
module load R/4.3.1

R

#[1] stringr_1.5.0   patchwork_1.3.0 ggplot2_3.4.2   dplyr_1.1.2    
#[5] reshape2_1.4.4 

div<-read.table("400_Gf_F0_vs_SL_only_Sel_lines_Hcwbps18_q20Q30_noindels_ss57_100kb.renamed.diversity.csv", header=TRUE)
summary(div)

div$mid<-((div$start + div$end)/2)

library(reshape2)
div_melt<-(melt(div, id.vars=c("chrom","start","end","mid")))
head(div_melt)



write.table(div_melt, file="400_Gf_F0_vs_SL_only_Sel_lines_Hcwbps18_q20Q30_noindels_ss57_100kb.renamed.diversity_MELTED.txt", row.names=FALSE, sep="\t", quote= FALSE) 
```

# Pi

Select just the Pi:

Note that there are two of these: theta_pi_abs and theta_pi_rel. From (https://www.biorxiv.org/content/biorxiv/early/2022/02/04/2022.02.02.477408.full.pdf) I think that they must represent the sample and an estimate of the population, but not certain. Theta_pi_rel is the one that I want to use.
```
# Give R a list of all filenames, can use a wildcard:

fileNames <- Sys.glob("400*100kb*renamed*diversity*csv")

# Then loop through them all:
for (fileName in fileNames) {

  # read data:
  df <- read.table(fileName,
    header = TRUE,
    sep = "\t",
    na.strings = "NA")

summary<-summary(df)
write.table(summary, paste0(fileName,'.summary.txt'), row.names=FALSE, sep="\t", quote= FALSE)

  # extract theta_pi_rel columns, create a mid-window co-ordinate:
# get_theta_pi_rel_cols <- function (my_data, half_window)
theta_pi_df <- get_theta_pi_rel_cols(df, 50000)

summary<-summary(theta_pi_df)
write.table(summary, paste0(fileName,'theta_pi_rel.summary.txt'), row.names=FALSE, sep="\t", quote= FALSE)}
```
To get the noramlised pi ratio between each sample and the normalised pi for the F0 sample (Sample/F0), using the output table from the function get_theta_pi_rel_cols() :
```
# The header line:
chrom   mid     F0.theta_pi_rel F1_C_L1.theta_pi_rel    F1_C_L2.theta_pi_rel    F1_C_L3.theta_pi_rel    F2_C_L1.theta_pi_rel    F2_C_L2.theta_pi_rel    F2_C_L3.theta_pi_rel        F3_C_L1.theta_pi_rel    F3_C_L2.theta_pi_rel    F3_C_L3.theta_pi_rel    F1_I_L2.theta_pi_rel    F1_I_L3.theta_pi_rel    F2_I_L1.theta_pi_rel    F2_I_L2.theta_pi_rel        F2_I_L3.theta_pi_rel    F3_I_L1.theta_pi_rel    F3_I_L2.theta_pi_rel    F3_I_L3.theta_pi_rel    F1_M_L1.theta_pi_rel    F1_M_L2.theta_pi_rel    F1_M_L3.theta_pi_rel        F2_M_L1.theta_pi_rel    F2_M_L2.theta_pi_rel    F2_M_L3.theta_pi_rel    F3_M_L1.theta_pi_rel    F3_M_L2.theta_pi_rel    F3_M_L3.theta_pi_rel

# ss57
cat 400_Gf_F0_vs_SL_only_Sel_lines_Hcwbps18_q20Q30_noindels_ss57_100kb.renamed.diversity.csv.theta_pi_rel.txt | tail -n +2 | awk 'BEGIN {OFS=FS="\t"} {print $1, $2, $3,($3 / 0.02182 ), (($4 / 0.02181  ) / ($3 / 0.02182   )),(($5 / 0.02172   )/ ($3 / 0.02182   )),(($6 / 0.02183   ) / ($3 / 0.02182   )),(($7 / 0.02180   )/ ($3 / 0.02182   )),(($8 / 0.02184 )/ ($3 / 0.02182   )),(($9 / 0.02193 )/ ($3 / 0.02182   )),(($10 / 0.02192  ) / ($3 / 0.02182   )),(($11 / 0.02191 ) / ($3 / 0.02182   )),(($12 /0.02204  )/($3 / 0.02182   )),(($13 / 0.02219    ) / ($3 / 0.02182   )),(($14 / 0.02125    ) / ($3 / 0.02182   )),(($15 / 0.02001   )/ ($3 / 0.02182   )),(($16 / 0.02115     ) / ($3 / 0.02182   )),(($17 / 0.02143     )/ ($3 / 0.02182   )),(($18 / 0.02140 )/ ($3 / 0.02182   )),(($19 / 0.01956   )/ ($3 / 0.02182   )),(($20 / 0.02105     )/ ($3 / 0.02182   )),(($21 / 0.02533      ) / ($3 / 0.02182   )),(($22 / 0.02180     ) / ($3 / 0.02182   )),(($23 / 0.02158      ) / ($3 / 0.02182   )),(($24 / 0.02462      ) / ($3 / 0.02182   )),(($25 / 0.02120    )/ ($3 / 0.02182   )),(($26 / 0.02063     ) / ($3 / 0.02182   )),(($27 / 0.02536       )/ ($3 / 0.02182   )),(($28 / 0.02035   )/ ($3 / 0.02182   )),(($29 / 0.02054    )/ ($3 / 0.02182   ))}' > 400_Gf_F0_vs_SL_only_Sel_lines_Hcwbps18_q20Q30_noindels_ss57_100kb.renamed.diversity.csv.theta_pi_rel.RATIO_F0.txt

# header:
chrom   mid  F0_pi   F0.theta_pi_rel F1_C_L1.theta_pi_rel    F1_C_L2.theta_pi_rel    F1_C_L3.theta_pi_rel    F2_C_L1.theta_pi_rel    F2_C_L2.theta_pi_rel    F2_C_L3.theta_pi_rel        F3_C_L1.theta_pi_rel    F3_C_L2.theta_pi_rel    F3_C_L3.theta_pi_rel    F1_I_L2.theta_pi_rel    F1_I_L3.theta_pi_rel    F2_I_L1.theta_pi_rel    F2_I_L2.theta_pi_rel        F2_I_L3.theta_pi_rel    F3_I_L1.theta_pi_rel    F3_I_L2.theta_pi_rel    F3_I_L3.theta_pi_rel    F1_M_L1.theta_pi_rel    F1_M_L2.theta_pi_rel    F1_M_L3.theta_pi_rel        F2_M_L1.theta_pi_rel    F2_M_L2.theta_pi_rel    F2_M_L3.theta_pi_rel    F3_M_L1.theta_pi_rel    F3_M_L2.theta_pi_rel    F3_M_L3.theta_pi_rel
```
sed 's/\.theta_pi_rel//g' header | sed 's/ \+/\t/g' > header_ratio
```
cat header_ratio 400_Gf_F0_vs_SL_only_Sel_lines_Hcwbps18_q20Q30_noindels_ss57_100kb.renamed.diversity.csv.theta_pi_rel.RATIO_F0.txt > tmp
less tmp
mv tmp 400_Gf_F0_vs_SL_only_Sel_lines_Hcwbps18_q20Q30_noindels_ss57_100kb.renamed.diversity.csv.theta_pi_rel.RATIO_F0.txt
```

To sort and rank the ratio - so that the greatest reduction is #1 and the greatest increase (in post cf pre) is #n.
```
# Will try for each comparison separately:

# My file looks like:

chrom	mid	F0_pi	F0	F1_C_L1	F1_C_L2	F1_C_L3	F2_C_L1	F2_C_L2	F2_C_L3	F3_C_L1	F3_C_L2	F3_C_L3	F1_I_L2	F1_I_L3	F2_I_L1	F2_I_L2	F2_I_L3	F3_I_L1	F3_I_L2	F3_I_L3	F1_M_L1	F1_M_L2	F1_M_L3	F2_M_L1	F2_M_L2	F2_M_L3	F3_M_L1	F3_M_L2	F3_M_L3
4	50000	0.00565564217	0.991355	0.976985	0.988541	1.00496	1.00103	1.03259	1.0421	1.05755	0.99557	0.854002	0.918812	1.01392	0.871337	1.01273	0.991815	1.06994	0.900253	0.841686	0.958588	0.942451	0.795475	0.980723	1.017	0.711578	1.0022	1.01394
4	150000	0.00554552421	0.990221	0.985105	0.999087	0.993464	0.997465	1.0531	1.05038	1.03836	0.992669	0.872924	0.920394	0.993604	0.902009	0.992495	0.996655	1.05546	0.922339	0.845638	0.994342	0.942165	0.795437	0.985304	1.00631	0.723654	1.03124	0.995412


# Remove the header, so can sort the file:
cat 400_Gf_F0_vs_SL_only_Sel_lines_Hcwbps18_q20Q30_noindels_ss57_100kb.renamed.diversity.csv.theta_pi_rel.RATIO_F0.txt | tail -n +2 > tmp

wc -l tmp
# 2837 tmp

# Then, remove lines with nan - note that doing '-nan' which is what is actually written in the file won't work with grep. 
cat tmp | grep -v 'nan' > tmp2

wc -l tmp2
# 2828 tmp2

# Add line number to file (i.e. rank by ratio). Note that the column is always added to the START of the file. (nl is what will add the line number)
cat tmp2 | sort -k 4 -n | nl | sort -k 6 -n | nl | sort -k 8 -n | nl | sort -k 10 -n | nl | sort -k 12 -n | nl | sort -k 14 -n | nl | sort -k 16 -n | nl | sort -k 18 -n | nl | sort -k 20 -n | nl | sort -k 22 -n | nl | sort -k 24 -n | nl | sort -k 26 -n | nl | sort -k 28 -n | nl  | sort -k 30 -n | nl | sort -k 32 -n | nl | sort -k 34 -n | nl | sort -k 36 -n | nl | sort -k 38 -n | nl | sort -k 40 -n | nl | sort -k 42 -n | nl | sort -k 44 -n | nl | sort -k 46 -n | nl | sort -k 48 -n | nl | sort -k 50 -n | nl | sort -k 52 -n | nl | sort -k 54 -n | nl | sort -k 56 -n | nl > tmp_sorted


# F3_M_L3R F3_M_L2R F3_M_L1R F2_M_L3R F2_M_L2R F2_M_L1R F1_M_L3R F1_M_L2R F1_M_L1R F3_I_L3R F3_I_L2R F3_I_L1R F2_I_L3R F2_I_L2R F2_I_L1R F1_I_L3R F1_I_L2R F3_C_L3R F3_C_L2R F3_C_L1R F2_C_L3R F2_C_L2R F2_C_L1R F1_C_L3R F1_C_L2R F1_C_L1R F0R chrom	mid	F0_pi	F0	F1_C_L1	F1_C_L2	F1_C_L3	F2_C_L1	F2_C_L2	F2_C_L3	F3_C_L1	F3_C_L2	F3_C_L3	F1_I_L2	F1_I_L3	F2_I_L1	F2_I_L2	F2_I_L3	F3_I_L1	F3_I_L2	F3_I_L3	F1_M_L1	F1_M_L2	F1_M_L3	F2_M_L1	F2_M_L2	F2_M_L3	F3_M_L1	F3_M_L2	F3_M_L3

# Add the new header back on!
cat new_header_rank tmp_sorted > 400_Gf_F0_vs_SL_only_Sel_lines_Hcwbps18_q20Q30_noindels_ss57_100kb.renamed.diversity.csv.theta_pi_rel.RATIO_F0.txt_ranked

sed 's/ \+/\t/g' 400_Gf_F0_vs_SL_only_Sel_lines_Hcwbps18_q20Q30_noindels_ss57_100kb.renamed.diversity.csv.theta_pi_rel.RATIO_F0.txt_ranked | sed 's/\t\+/\t/g' | sed 's/^\t//g' > tmp

mv tmp 400_Gf_F0_vs_SL_only_Sel_lines_Hcwbps18_q20Q30_noindels_ss57_100kb.renamed.diversity.csv.theta_pi_rel.RATIO_F0.txt_ranked

wc -l 400_Gf_F0_vs_SL_only_Sel_lines_Hcwbps18_q20Q30_noindels_ss57_100kb.renamed.diversity.csv.theta_pi_rel.RATIO_F0.txt_ranked
# 2829 400_Gf_F0_vs_SL_only_Sel_lines_Hcwbps18_q20Q30_noindels_ss57_100kb.renamed.diversity.csv.theta_pi_rel.RATIO_F0.txt_ranked
```
Therefore, there are 2828 100kb windows in the ranked file (remove the header line from the total). 

Therefore, the top 1% with the greatest reduction relative to F0 would be ranked 28 or higher (1-28), and the lowest 1% (those with the greatest increase in diversity post-tx) would be 2801 or lower.

Let's split the file and then melt the sub-files, and recombine to be able to plot all samples together. 
```
awk '{print NF}' 400_Gf_F0_vs_SL_only_Sel_lines_Hcwbps18_q20Q30_noindels_ss57_100kb.renamed.diversity.csv.theta_pi_rel.RATIO_F0.txt_ranked | uniq > out
```
There are 57 columns. Of which, 3 are in the middle and separate from the rest. 

Therefore I should have 27 rank columns and 27 ratio columns - 1+9+8+9 = 27 so correct. :) 

```
cut -f 1-29 400_Gf_F0_vs_SL_only_Sel_lines_Hcwbps18_q20Q30_noindels_ss57_100kb.renamed.diversity.csv.theta_pi_rel.RATIO_F0.txt_ranked > tmp_rank_cols
cut -f 28-30 400_Gf_F0_vs_SL_only_Sel_lines_Hcwbps18_q20Q30_noindels_ss57_100kb.renamed.diversity.csv.theta_pi_rel.RATIO_F0.txt_ranked > tmp_pi_cols
cut -f 28-29,31- 400_Gf_F0_vs_SL_only_Sel_lines_Hcwbps18_q20Q30_noindels_ss57_100kb.renamed.diversity.csv.theta_pi_rel.RATIO_F0.txt_ranked > tmp_ratio_cols
```
Then, melt them in R, and combine together: 
```
module load R/4.3.1
R
```
```
library(reshape2)
library(dplyr)
library(ggplot2)
library(patchwork)

df<-read.table("tmp_rank_cols", header=TRUE, sep = "\t",  na.strings = "NA")
df2<-read.table("tmp_pi_cols", header=TRUE, sep = "\t",  na.strings = "NA")
df3<-read.table("tmp_ratio_cols", header=TRUE, sep = "\t",  na.strings = "NA")

df_melt<-melt(df, id.vars=c("chrom","mid"))
df3_melt<-melt(df3, id.vars=c("chrom","mid"))

df_melt$variable<-gsub("R","", df_melt$variable)

# Join tables, keeping only rows present in both tables
df4<-inner_join(df_melt, df3_melt, by=c("chrom"="chrom", "mid"="mid", "variable"="variable"))
```
```
# Rename chromosomes:
df4$chrom<-gsub("1","I", df4$chrom)
df4$chrom<-gsub("2","II", df4$chrom)
df4$chrom<-gsub("3","III", df4$chrom)
df4$chrom<-gsub("4","IV", df4$chrom)
df4$chrom<-gsub("5","V", df4$chrom)


df2$chrom<-gsub("1","I", df2$chrom)
df2$chrom<-gsub("2","II", df2$chrom)
df2$chrom<-gsub("3","III", df2$chrom)
df2$chrom<-gsub("4","IV", df2$chrom)
df2$chrom<-gsub("5","V", df2$chrom)

```
To try plotting chr V for all samples, colouring by the rank: 

```
# Filter the data to retain each TREATMENT separately:

library(stringr)
library(dplyr)

mox<- df4 %>%
  filter(str_detect(variable, "M"))

ctl<- df4 %>%
  filter(str_detect(variable, "C"))

ivm<- df4 %>%
  filter(str_detect(variable, "I"))


# Then, add in the extra facet:
ivm$variable.new <- factor(ivm$variable, levels=(c("", "F1_I_L2", "F1_I_L3", "F2_I_L1", "F2_I_L2", "F2_I_L3", "F3_I_L1", "F3_I_L2", "F3_I_L3")))   

head(ctl)
  chrom      mid variable value.x  value.y
1     3 16950000  F3_C_L3    2653 1.083350
2     5 37750000  F3_C_L3     961 0.980703
3     5 36750000  F3_C_L3    1114 0.986404
4     5 37150000  F3_C_L3    1800 1.010280
5     5 13250000  F3_C_L3    2784 1.152580
6     5 37950000  F3_C_L3     705 0.970568

```
And try plotting everything....
```
library(ggplot2)
library(patchwork)

IVM<- ggplot(data=ivm, mapping=aes(x=mid/1e6, y=log10(value.y))) +
geom_point(aes(colour = value.x), alpha=1, size=0.3) +
scale_colour_gradientn('Rank',
  colours = c("red","grey","blue"),
                         values = c(0,0.01,0.99,1)
)+
labs(y="Log10{Ratio of [(Sample.Pi/Sample.Median.Pi) /\n(F0.Pi/F0.Median.Pi)]} (100kb window)", x="Genomic position (Mb)")+
facet_grid(variable ~ chrom, scales="free_x", switch="x")+
scale_x_continuous(
    breaks = seq(0, 55, 20),
    minor_breaks = seq(5, 50, 5)
  )+
theme_light()+
theme(axis.title.x = element_text(size=10), axis.text.x = element_text(size=6), axis.title.y = element_text(size=10), axis.text.y = element_text(size=4), strip.text.x = element_text(size = 6, face = "bold", margin = margin(0.5,0,0.5,0, "mm")), strip.text.y = element_text(size = 6, face = "bold", angle=0), plot.title = element_text(size=10, hjust=0.5), legend.position = "top", legend.title = element_text(size=8), legend.text = element_text(size=6), legend.key.height=unit(0.2,"cm"), panel.spacing.x = unit(0.2, "mm"), panel.spacing.y = unit(0.75, "mm"))+
ggtitle('Ivermectin lines')

MOX<- ggplot(data=mox, mapping=aes(x=mid/1e6, y=log10(value.y))) +
geom_point(aes(colour = value.x), alpha=1, size=0.3) +
scale_colour_gradientn('Rank',
  colours = c("red","grey","blue"),
                         values = c(0,0.01,0.99,1)
)+
labs(y="Log10{Ratio of [(Sample.Pi/Sample.Median.Pi) /\n(F0.Pi/F0.Median.Pi)]} (100kb window)", x="Genomic position (Mb)")+
facet_grid(variable ~ chrom, scales="free_x", switch="x")+
scale_x_continuous(
    breaks = seq(0, 55, 20),
    minor_breaks = seq(5, 50, 5)
  )+
theme_light()+
theme(axis.title.x = element_text(size=10), axis.text.x = element_text(size=6), axis.title.y = element_text(size=10), axis.text.y = element_text(size=4), strip.text.x = element_text(size = 6, face = "bold", margin = margin(0.5,0,0.5,0, "mm")), strip.text.y = element_text(size = 6, face = "bold", angle=0), plot.title = element_text(size=10, hjust=0.5), legend.position = "top", legend.title = element_text(size=8), legend.text = element_text(size=6), legend.key.height=unit(0.2,"cm"), panel.spacing.x = unit(0.2, "mm"), panel.spacing.y = unit(0.75, "mm"))+
ggtitle('Moxidectin lines')


CTL<- ggplot(data=ctl, mapping=aes(x=mid/1e6, y=log10(value.y))) +
geom_point(aes(colour = value.x), alpha=1, size=0.3) +
scale_colour_gradientn('Rank',
  colours = c("red","grey","blue"),
                         values = c(0,0.01,0.99,1)
)+
labs(y="Log10{Ratio of [(Sample.Pi/Sample.Median.Pi) /\n(F0.Pi/F0.Median.Pi)]} (100kb window)", x="Genomic position (Mb)")+
facet_grid(variable ~ chrom, scales="free_x", switch="x")+
scale_x_continuous(
    breaks = seq(0, 55, 20),
    minor_breaks = seq(5, 50, 5)
  )+
theme_light()+
theme(axis.title.x = element_text(size=10), axis.text.x = element_text(size=6), axis.title.y = element_text(size=10), axis.text.y = element_text(size=4), strip.text.x = element_text(size = 6, face = "bold", margin = margin(0.5,0,0.5,0, "mm")), strip.text.y = element_text(size = 6, face = "bold", angle=0), plot.title = element_text(size=10, hjust=0.5), legend.position = "top", legend.title = element_text(size=8), legend.text = element_text(size=6), legend.key.height=unit(0.2,"cm"), panel.spacing.x = unit(0.2, "mm"), panel.spacing.y = unit(0.75, "mm"))+
ggtitle('Control lines')

F0<- ggplot(data=df2, mapping=aes(x=mid/1e6, y=log10(F0_pi))) +
geom_point(aes(colour = chrom), alpha=1, size=0.3) +
labs(y="Log10(Relative Theta Pi)\n(100kb window)", x="Genomic position (Mb)")+
facet_grid(chrom ~ ., scales="free_x", switch="x")+
scale_x_continuous(
    breaks = seq(0, 55, 20),
    minor_breaks = seq(5, 50, 5)
  )+
theme_light()+
theme(axis.title.x = element_text(size=10), axis.text.x = element_text(size=6), axis.title.y = element_text(size=10), axis.text.y = element_text(size=4), strip.text.x = element_text(size = 6, face = "bold", margin = margin(0.5,0,0.5,0, "mm")), strip.text.y = element_text(size = 6, face = "bold", angle=0), plot.title = element_text(size=10, hjust=0.5), legend.position = "none", panel.spacing.x = unit(0.2, "mm"), panel.spacing.y = unit(0.75, "mm"))+
guides(colour = "none") +
ggtitle('F0 Pi')

patch<- IVM + MOX + CTL + F0 + plot_annotation(tag_levels = 'A') + plot_layout (axis_titles = "collect")
patch

#ggsave("EmergingMLR_manuscript_FigX_light_titles_Pi_ratio.tiff", patch, device=tiff, width=2250, height=2250, units="px", dpi=320)
ggsave("EmergingMLR_manuscript_FigX_light_titles_Pi_ratio_withromannumerals.tiff", patch, device=tiff, width=2250, height=2250, units="px", dpi=320)
```
_Image with roman numerals:_
![image](https://github.com/user-attachments/assets/30273a7f-2337-41cc-81df-5f6b1edea855)

_Original image:_
![image](https://github.com/user-attachments/assets/7f8c422d-c6fb-4d64-b743-d10493472fea)

_With simpler y-axis:_
```
library(ggplot2)
library(patchwork)

IVM<- ggplot(data=ivm, mapping=aes(x=mid/1e6, y=log10(value.y))) +
geom_point(aes(colour = value.x), alpha=1, size=0.3) +
scale_colour_gradientn('Rank',
  colours = c("red","grey","blue"),
                         values = c(0,0.01,0.99,1)
)+
labs(y="Log10(Pi ratio)\n(100kb window)", x="Genomic position (Mb)")+
facet_grid(variable ~ chrom, scales="free_x", switch="x")+
scale_x_continuous(
    breaks = seq(0, 55, 20),
    minor_breaks = seq(5, 50, 5)
  )+
theme_light()+
theme(axis.title.x = element_text(size=10), axis.text.x = element_text(size=6), axis.title.y = element_text(size=10), axis.text.y = element_text(size=4), strip.text.x = element_text(size = 6, face = "bold", margin = margin(0.5,0,0.5,0, "mm")), strip.text.y = element_text(size = 6, face = "bold", angle=0), plot.title = element_text(size=10, hjust=0.5), legend.position = "top", legend.title = element_text(size=8), legend.text = element_text(size=6), legend.key.height=unit(0.2,"cm"), panel.spacing.x = unit(0.2, "mm"), panel.spacing.y = unit(0.75, "mm"))+
ggtitle('Ivermectin lines')

MOX<- ggplot(data=mox, mapping=aes(x=mid/1e6, y=log10(value.y))) +
geom_point(aes(colour = value.x), alpha=1, size=0.3) +
scale_colour_gradientn('Rank',
  colours = c("red","grey","blue"),
                         values = c(0,0.01,0.99,1)
)+
labs(y="Log10(Pi ratio)\n(100kb window)", x="Genomic position (Mb)")+
facet_grid(variable ~ chrom, scales="free_x", switch="x")+
scale_x_continuous(
    breaks = seq(0, 55, 20),
    minor_breaks = seq(5, 50, 5)
  )+
theme_light()+
theme(axis.title.x = element_text(size=10), axis.text.x = element_text(size=6), axis.title.y = element_text(size=10), axis.text.y = element_text(size=4), strip.text.x = element_text(size = 6, face = "bold", margin = margin(0.5,0,0.5,0, "mm")), strip.text.y = element_text(size = 6, face = "bold", angle=0), plot.title = element_text(size=10, hjust=0.5), legend.position = "top", legend.title = element_text(size=8), legend.text = element_text(size=6), legend.key.height=unit(0.2,"cm"), panel.spacing.x = unit(0.2, "mm"), panel.spacing.y = unit(0.75, "mm"))+
ggtitle('Moxidectin lines')


CTL<- ggplot(data=ctl, mapping=aes(x=mid/1e6, y=log10(value.y))) +
geom_point(aes(colour = value.x), alpha=1, size=0.3) +
scale_colour_gradientn('Rank',
  colours = c("red","grey","blue"),
                         values = c(0,0.01,0.99,1)
)+
labs(y="Log10(Pi ratio)\n(100kb window)", x="Genomic position (Mb)")+
facet_grid(variable ~ chrom, scales="free_x", switch="x")+
scale_x_continuous(
    breaks = seq(0, 55, 20),
    minor_breaks = seq(5, 50, 5)
  )+
theme_light()+
theme(axis.title.x = element_text(size=10), axis.text.x = element_text(size=6), axis.title.y = element_text(size=10), axis.text.y = element_text(size=4), strip.text.x = element_text(size = 6, face = "bold", margin = margin(0.5,0,0.5,0, "mm")), strip.text.y = element_text(size = 6, face = "bold", angle=0), plot.title = element_text(size=10, hjust=0.5), legend.position = "top", legend.title = element_text(size=8), legend.text = element_text(size=6), legend.key.height=unit(0.2,"cm"), panel.spacing.x = unit(0.2, "mm"), panel.spacing.y = unit(0.75, "mm"))+
ggtitle('Control lines')

F0<- ggplot(data=df2, mapping=aes(x=mid/1e6, y=log10(F0_pi))) +
geom_point(aes(colour = chrom), alpha=1, size=0.3) +
labs(y="Log10(Pi)\n(100kb window)", x="Genomic position (Mb)")+
facet_grid(chrom ~ ., scales="free_x", switch="x")+
scale_x_continuous(
    breaks = seq(0, 55, 20),
    minor_breaks = seq(5, 50, 5)
  )+
theme_light()+
theme(axis.title.x = element_text(size=10), axis.text.x = element_text(size=6), axis.title.y = element_text(size=10), axis.text.y = element_text(size=4), strip.text.x = element_text(size = 6, face = "bold", margin = margin(0.5,0,0.5,0, "mm")), strip.text.y = element_text(size = 6, face = "bold", angle=0), plot.title = element_text(size=10, hjust=0.5), legend.position = "none", panel.spacing.x = unit(0.2, "mm"), panel.spacing.y = unit(0.75, "mm"))+
guides(colour = "none") +
ggtitle('F0 Pi')

patch<- IVM + MOX + CTL + F0 + plot_annotation(tag_levels = 'A') + plot_layout (axis_titles = "collect")
patch

#ggsave("EmergingMLR_manuscript_FigX_light_titles_Pi_ratio.tiff", patch, device=tiff, width=2250, height=2250, units="px", dpi=320)
#ggsave("EmergingMLR_manuscript_FigX_light_titles_Pi_ratio_withromannumerals.tiff", patch, device=tiff, width=2250, height=2250, units="px", dpi=320)
ggsave("EmergingMLR_manuscript_FigX_light_titles_Pi_ratio_withromannumerals_simplerYaxis.tiff", patch, device=tiff, width=2250, height=2250, units="px", dpi=320)
```

And try plotting just ChrV....
```
library(ggplot2)
library(patchwork)

IVM<- ggplot(data=subset(ivm, ivm$chrom == "V"), mapping=aes(x=mid/1e6, y=log10(value.y))) +
geom_point(aes(colour = value.x), alpha=1, size=0.3) +
scale_colour_gradientn('Rank',
  colours = c("red","grey","blue"),
                         values = c(0,0.01,0.99,1)
)+
labs(y="Log10{Ratio of [(Sample.Pi/Sample.Median.Pi) /\n(F0.Pi/F0.Median.Pi)]} (100kb window)", x="Genomic position (Mb)")+
facet_grid(variable ~ chrom, scales="free_x", switch="x")+
scale_x_continuous(
    breaks = seq(0, 55, 5),
    minor_breaks = seq(1, 54, 1)
  )+
theme_light()+
theme(axis.title.x = element_text(size=10), axis.text.x = element_text(size=6), axis.title.y = element_text(size=10), axis.text.y = element_text(size=4), strip.text.x = element_text(size = 6, face = "bold", margin = margin(0.5,0,0.5,0, "mm")), strip.text.y = element_text(size = 6, face = "bold", angle=0), plot.title = element_text(size=10, hjust=0.5), legend.position = "top", legend.title = element_text(size=8), legend.text = element_text(size=6), legend.key.height=unit(0.2,"cm"), panel.spacing.x = unit(0.2, "mm"), panel.spacing.y = unit(0.75, "mm"))+
ggtitle('Ivermectin lines')

MOX<- ggplot(data=subset(mox, mox$chrom == "V"), mapping=aes(x=mid/1e6, y=log10(value.y))) +
geom_point(aes(colour = value.x), alpha=1, size=0.3) +
scale_colour_gradientn('Rank',
  colours = c("red","grey","blue"),
                         values = c(0,0.01,0.99,1)
)+
labs(y="Log10{Ratio of [(Sample.Pi/Sample.Median.Pi) /\n(F0.Pi/F0.Median.Pi)]} (100kb window)", x="Genomic position (Mb)")+
facet_grid(variable ~ chrom, scales="free_x", switch="x")+
scale_x_continuous(
    breaks = seq(0, 55, 5),
    minor_breaks = seq(1, 54, 1)
  )+
theme_light()+
theme(axis.title.x = element_text(size=10), axis.text.x = element_text(size=6), axis.title.y = element_text(size=10), axis.text.y = element_text(size=4), strip.text.x = element_text(size = 6, face = "bold", margin = margin(0.5,0,0.5,0, "mm")), strip.text.y = element_text(size = 6, face = "bold", angle=0), plot.title = element_text(size=10, hjust=0.5), legend.position = "top", legend.title = element_text(size=8), legend.text = element_text(size=6), legend.key.height=unit(0.2,"cm"), panel.spacing.x = unit(0.2, "mm"), panel.spacing.y = unit(0.75, "mm"))+
ggtitle('Moxidectin lines')


CTL<- ggplot(data=subset(ctl, ctl$chrom == "V"), mapping=aes(x=mid/1e6, y=log10(value.y))) +
geom_point(aes(colour = value.x), alpha=1, size=0.3) +
scale_colour_gradientn('Rank',
  colours = c("red","grey","blue"),
                         values = c(0,0.01,0.99,1)
)+
labs(y="Log10{Ratio of [(Sample.Pi/Sample.Median.Pi) /\n(F0.Pi/F0.Median.Pi)]} (100kb window)", x="Genomic position (Mb)")+
facet_grid(variable ~ chrom, scales="free_x", switch="x")+
scale_x_continuous(
    breaks = seq(0, 55, 5),
    minor_breaks = seq(1, 54, 1)
  )+
theme_light()+
theme(axis.title.x = element_text(size=10), axis.text.x = element_text(size=6), axis.title.y = element_text(size=10), axis.text.y = element_text(size=4), strip.text.x = element_text(size = 6, face = "bold", margin = margin(0.5,0,0.5,0, "mm")), strip.text.y = element_text(size = 6, face = "bold", angle=0), plot.title = element_text(size=10, hjust=0.5), legend.position = "top", legend.title = element_text(size=8), legend.text = element_text(size=6), legend.key.height=unit(0.2,"cm"), panel.spacing.x = unit(0.2, "mm"), panel.spacing.y = unit(0.75, "mm"))+
ggtitle('Control lines')

F0<- ggplot(data=subset(df2, df2$chrom == "V"), mapping=aes(x=mid/1e6, y=log10(F0_pi))) +
geom_point(colour = "#619CFF", alpha=1, size=0.3) +
labs(y="Log10(Relative Theta Pi)\n(100kb window)", x="Genomic position (Mb)")+
facet_grid(. ~ chrom, scales="free_x", switch="x")+
scale_x_continuous(
    breaks = seq(0, 55, 5),
    minor_breaks = seq(1, 54, 1)
  )+
theme_light()+
theme(axis.title.x = element_text(size=10), axis.text.x = element_text(size=6), axis.title.y = element_text(size=10), axis.text.y = element_text(size=4), strip.text.x = element_text(size = 6, face = "bold", margin = margin(0.5,0,0.5,0, "mm")), strip.text.y = element_text(size = 6, face = "bold", angle=0), plot.title = element_text(size=10, hjust=0.5), legend.position = "none", panel.spacing.x = unit(0.2, "mm"), panel.spacing.y = unit(0.75, "mm"))+
guides(colour = "none") +
ggtitle('F0 Pi')

patch<- IVM + MOX + CTL + F0 + plot_annotation(tag_levels = 'A') + plot_layout (axis_titles = "collect")
patch

#ggsave("EmergingMLR_manuscript_FigX_light_titles_Pi_ratio.tiff", patch, device=tiff, width=2250, height=2250, units="px", dpi=320)
ggsave("EmergingMLR_manuscript_FigX_light_titles_Pi_ratio_ChrV_withromannumerals.tiff", patch, device=tiff, width=2250, height=2250, units="px", dpi=320)

```
_Image with roman numerals:_
![image](https://github.com/user-attachments/assets/b2d8bdd0-f336-4539-b2fc-9e93f04145a9)

_Original image:_
![image](https://github.com/user-attachments/assets/7e257ce1-1273-49f6-ab62-2526cb41fddf)

_With simpler y-axis:_
```
library(ggplot2)
library(patchwork)

IVM<- ggplot(data=subset(ivm, ivm$chrom == "V"), mapping=aes(x=mid/1e6, y=log10(value.y))) +
geom_point(aes(colour = value.x), alpha=1, size=0.3) +
scale_colour_gradientn('Rank',
  colours = c("red","grey","blue"),
                         values = c(0,0.01,0.99,1)
)+
labs(y="Log10(Pi ratio)\n(100kb window)", x="Genomic position (Mb)")+
facet_grid(variable ~ chrom, scales="free_x", switch="x")+
scale_x_continuous(
    breaks = seq(0, 55, 5),
    minor_breaks = seq(1, 54, 1)
  )+
theme_light()+
theme(axis.title.x = element_text(size=10), axis.text.x = element_text(size=6), axis.title.y = element_text(size=10), axis.text.y = element_text(size=4), strip.text.x = element_text(size = 6, face = "bold", margin = margin(0.5,0,0.5,0, "mm")), strip.text.y = element_text(size = 6, face = "bold", angle=0), plot.title = element_text(size=10, hjust=0.5), legend.position = "top", legend.title = element_text(size=8), legend.text = element_text(size=6), legend.key.height=unit(0.2,"cm"), panel.spacing.x = unit(0.2, "mm"), panel.spacing.y = unit(0.75, "mm"))+
ggtitle('Ivermectin lines')

MOX<- ggplot(data=subset(mox, mox$chrom == "V"), mapping=aes(x=mid/1e6, y=log10(value.y))) +
geom_point(aes(colour = value.x), alpha=1, size=0.3) +
scale_colour_gradientn('Rank',
  colours = c("red","grey","blue"),
                         values = c(0,0.01,0.99,1)
)+
labs(y="Log10(Pi ratio)\n(100kb window)", x="Genomic position (Mb)")+
facet_grid(variable ~ chrom, scales="free_x", switch="x")+
scale_x_continuous(
    breaks = seq(0, 55, 5),
    minor_breaks = seq(1, 54, 1)
  )+
theme_light()+
theme(axis.title.x = element_text(size=10), axis.text.x = element_text(size=6), axis.title.y = element_text(size=10), axis.text.y = element_text(size=4), strip.text.x = element_text(size = 6, face = "bold", margin = margin(0.5,0,0.5,0, "mm")), strip.text.y = element_text(size = 6, face = "bold", angle=0), plot.title = element_text(size=10, hjust=0.5), legend.position = "top", legend.title = element_text(size=8), legend.text = element_text(size=6), legend.key.height=unit(0.2,"cm"), panel.spacing.x = unit(0.2, "mm"), panel.spacing.y = unit(0.75, "mm"))+
ggtitle('Moxidectin lines')


CTL<- ggplot(data=subset(ctl, ctl$chrom == "V"), mapping=aes(x=mid/1e6, y=log10(value.y))) +
geom_point(aes(colour = value.x), alpha=1, size=0.3) +
scale_colour_gradientn('Rank',
  colours = c("red","grey","blue"),
                         values = c(0,0.01,0.99,1)
)+
labs(y="Log10(Pi ratio)\n(100kb window)", x="Genomic position (Mb)")+
facet_grid(variable ~ chrom, scales="free_x", switch="x")+
scale_x_continuous(
    breaks = seq(0, 55, 5),
    minor_breaks = seq(1, 54, 1)
  )+
theme_light()+
theme(axis.title.x = element_text(size=10), axis.text.x = element_text(size=6), axis.title.y = element_text(size=10), axis.text.y = element_text(size=4), strip.text.x = element_text(size = 6, face = "bold", margin = margin(0.5,0,0.5,0, "mm")), strip.text.y = element_text(size = 6, face = "bold", angle=0), plot.title = element_text(size=10, hjust=0.5), legend.position = "top", legend.title = element_text(size=8), legend.text = element_text(size=6), legend.key.height=unit(0.2,"cm"), panel.spacing.x = unit(0.2, "mm"), panel.spacing.y = unit(0.75, "mm"))+
ggtitle('Control lines')

F0<- ggplot(data=subset(df2, df2$chrom == "V"), mapping=aes(x=mid/1e6, y=log10(F0_pi))) +
geom_point(colour = "#619CFF", alpha=1, size=0.3) +
labs(y="Log10(Pi)\n(100kb window)", x="Genomic position (Mb)")+
facet_grid(. ~ chrom, scales="free_x", switch="x")+
scale_x_continuous(
    breaks = seq(0, 55, 5),
    minor_breaks = seq(1, 54, 1)
  )+
theme_light()+
theme(axis.title.x = element_text(size=10), axis.text.x = element_text(size=6), axis.title.y = element_text(size=10), axis.text.y = element_text(size=4), strip.text.x = element_text(size = 6, face = "bold", margin = margin(0.5,0,0.5,0, "mm")), strip.text.y = element_text(size = 6, face = "bold", angle=0), plot.title = element_text(size=10, hjust=0.5), legend.position = "none", panel.spacing.x = unit(0.2, "mm"), panel.spacing.y = unit(0.75, "mm"))+
guides(colour = "none") +
ggtitle('F0 Pi')

patch<- IVM + MOX + CTL + F0 + plot_annotation(tag_levels = 'A') + plot_layout (axis_titles = "collect")
patch

#ggsave("EmergingMLR_manuscript_FigX_light_titles_Pi_ratio.tiff", patch, device=tiff, width=2250, height=2250, units="px", dpi=320)
#ggsave("EmergingMLR_manuscript_FigX_light_titles_Pi_ratio_ChrV_withromannumerals.tiff", patch, device=tiff, width=2250, height=2250, units="px", dpi=320)
ggsave("EmergingMLR_manuscript_FigX_light_titles_Pi_ratio_ChrV_withromannumerals_simpler_y-axis.tiff", patch, device=tiff, width=2250, height=2250, units="px", dpi=320)

```


## Tajima's D
Select just the Tajima's D:

Note that I am using the 'renamed' grenedalf output file which has the chromosomes as simply a number (1, 2 etc) rather than the long name.
```
#R v4.3.1
#[1] stringr_1.5.0   patchwork_1.3.0 ggplot2_3.4.2   dplyr_1.1.2    
#[5] reshape2_1.4.4

# Give R a list of all filenames, can use a wildcard:

fileNames <- Sys.glob("400*100kb*renamed*diversity*csv")

# Then loop through them all:
for (fileName in fileNames) {

  # read data:
  df <- read.table(fileName,
    header = TRUE,
    sep = "\t",
    na.strings = "NA")

summary<-summary(df)
write.table(summary, paste0(fileName,'.summary.txt'), row.names=FALSE, sep="\t", quote= FALSE)

  # extract tajima's D columns, create a mid-window co-ordinate:
# get_tajimasD_cols <- function (my_data, half_window)
tajimasD_df <- get_tajimasD_cols(df, 50000)

summary<-summary(tajimasD_df)
write.table(summary, paste0(fileName,'tajimasD.summary.txt'), row.names=FALSE, sep="\t", quote= FALSE)}
```

```
# Melt the data:
library(reshape2)

dfmelt<-melt(tajimasD_df, id.vars = c("chrom", "mid"))

# Remove the .tajimas_d from each variable
dfmelt$variable<-gsub(".tajimas_d","", dfmelt$variable)
```
```
# Rename chromosomes:
dfmelt$chrom<-gsub("1","I", dfmelt$chrom)
dfmelt$chrom<-gsub("2","II", dfmelt$chrom)
dfmelt$chrom<-gsub("3","III", dfmelt$chrom)
dfmelt$chrom<-gsub("4","IV", dfmelt$chrom)
dfmelt$chrom<-gsub("5","V", dfmelt$chrom)
```
```
# Filter the data to retain each TREATMENT separately:

library(stringr)
library(dplyr)

f0<- dfmelt %>%
  filter(str_detect(variable, "F0"))

mox<- dfmelt %>%
  filter(str_detect(variable, "M"))

ctl<- dfmelt %>%
  filter(str_detect(variable, "C"))

ivm<- dfmelt %>%
  filter(str_detect(variable, "I"))


# Then, add in the extra facet (not convinced this is working, not sure why but not really important for these plots):
ivm$variable.new <- factor(ivm$variable, levels=c("", "F1_I_L2", "F1_I_L3", "F2_I_L1", "F2_I_L2", "F2_I_L3", "F3_I_L1", "F3_I_L2", "F3_I_L3"))   

head(ctl)
  chrom    mid          variable      value
1     4  50000 F1_C_L1.tajimas_d 0.60367862
2     4 150000 F1_C_L1.tajimas_d 0.32450977
3     4 250000 F1_C_L1.tajimas_d 0.54889517
4     4 350000 F1_C_L1.tajimas_d 0.92030134
5     4 450000 F1_C_L1.tajimas_d 0.06628736
6     4 550000 F1_C_L1.tajimas_d 0.37454241
```
And try plotting everything
```
library(ggplot2)
library(patchwork)

IVM<- ggplot(data=ivm, mapping=aes(x=mid/1e6, y=value)) +
geom_point(aes(colour = value), alpha=1, size=0.3) +
scale_colour_steps2('Tajimas D')+
labs(y="Tajima's D (100kb window)", x="Genomic position (Mb)")+
facet_grid(variable ~ chrom, scales="free_x", switch="x")+
scale_x_continuous(
    breaks = seq(0, 55, 20),
    minor_breaks = seq(5, 50, 5)
  )+
theme_light()+
theme(axis.title.x = element_text(size=10), axis.text.x = element_text(size=6), axis.title.y = element_text(size=10), axis.text.y = element_text(size=4), strip.text.x = element_text(size = 6, face = "bold", margin = margin(0.5,0,0.5,0, "mm")), strip.text.y = element_text(size = 6, face = "bold", angle=0), plot.title = element_text(size=10, hjust=0.5), legend.position = "top", legend.title = element_text(size=8), legend.text = element_text(size=6), legend.key.height=unit(0.2,"cm"), panel.spacing.x = unit(0.2, "mm"), panel.spacing.y = unit(0.75, "mm"))+
ggtitle('Ivermectin lines')

MOX<- ggplot(data=mox, mapping=aes(x=mid/1e6, y=value)) +
geom_point(aes(colour = value), alpha=1, size=0.3) +
scale_colour_steps2('Tajimas D')+
labs(y="Tajima's D (100kb window)", x="Genomic position (Mb)")+
facet_grid(variable ~ chrom, scales="free_x", switch="x")+
scale_x_continuous(
    breaks = seq(0, 55, 20),
    minor_breaks = seq(5, 50, 5)
  )+
theme_light()+
theme(axis.title.x = element_text(size=10), axis.text.x = element_text(size=6), axis.title.y = element_text(size=10), axis.text.y = element_text(size=4), strip.text.x = element_text(size = 6, face = "bold", margin = margin(0.5,0,0.5,0, "mm")), strip.text.y = element_text(size = 6, face = "bold", angle=0), plot.title = element_text(size=10, hjust=0.5), legend.position = "top", legend.title = element_text(size=8), legend.text = element_text(size=6), legend.key.height=unit(0.2,"cm"), panel.spacing.x = unit(0.2, "mm"), panel.spacing.y = unit(0.75, "mm"))+
ggtitle('Moxidectin lines')


CTL<- ggplot(data=ctl, mapping=aes(x=mid/1e6, y=value)) +
geom_point(aes(colour = value), alpha=1, size=0.3) +
scale_colour_steps2('Tajimas D')+
labs(y="Tajima's D (100kb window)", x="Genomic position (Mb)")+
facet_grid(variable ~ chrom, scales="free_x", switch="x")+
scale_x_continuous(
    breaks = seq(0, 55, 20),
    minor_breaks = seq(5, 50, 5)
  )+
theme_light()+
theme(axis.title.x = element_text(size=10), axis.text.x = element_text(size=6), axis.title.y = element_text(size=10), axis.text.y = element_text(size=4), strip.text.x = element_text(size = 6, face = "bold", margin = margin(0.5,0,0.5,0, "mm")), strip.text.y = element_text(size = 6, face = "bold", angle=0), plot.title = element_text(size=10, hjust=0.5), legend.position = "top", legend.title = element_text(size=8), legend.text = element_text(size=6), legend.key.height=unit(0.2,"cm"), panel.spacing.x = unit(0.2, "mm"), panel.spacing.y = unit(0.75, "mm"))+
ggtitle('Control lines')

F0<- ggplot(data=f0, mapping=aes(x=mid/1e6, y=value)) +
geom_point(aes(colour = value), alpha=1, size=0.3) +
scale_colour_steps2('Tajimas D')+
labs(y="Tajima's D (100kb window)", x="Genomic position (Mb)")+
facet_grid(chrom ~ ., scales="free_x", switch="x")+
scale_x_continuous(
    breaks = seq(0, 55, 20),
    minor_breaks = seq(5, 50, 5)
  )+
theme_light()+
theme(axis.title.x = element_text(size=10), axis.text.x = element_text(size=6), axis.title.y = element_text(size=10), axis.text.y = element_text(size=4), strip.text.x = element_text(size = 6, face = "bold", margin = margin(0.5,0,0.5,0, "mm")), strip.text.y = element_text(size = 6, face = "bold", angle=0), plot.title = element_text(size=10, hjust=0.5), legend.position = "top", legend.title = element_text(size=8), legend.text = element_text(size=6), legend.key.height=unit(0.2,"cm"), panel.spacing.x = unit(0.2, "mm"), panel.spacing.y = unit(0.75, "mm"))+
ggtitle('F0 Tajimas D')

patch<- IVM + MOX + CTL + F0 + plot_annotation(tag_levels = 'A') + plot_layout (axis_titles = "collect")
patch

#ggsave("EmergingMLR_manuscript_FigX_light_titles_TajimasD.tiff", patch, device=tiff, width=2250, height=2250, units="px", dpi=320)
ggsave("EmergingMLR_manuscript_FigX_light_titles_TajimasD_withromannumerals.tiff", patch, device=tiff, width=2250, height=2250, units="px", dpi=320)

```
_Image with roman numerals:_
![image](https://github.com/user-attachments/assets/90881d8f-e852-4a93-9f22-96f8227e376b)

_Original image:_
![image](https://github.com/user-attachments/assets/e1b9c947-f1ec-431e-9ab5-39969019ea5f)




And to plot just Chr5:
```
library(ggplot2)
library(patchwork)

IVM<- ggplot(data=subset(ivm, ivm$chrom == "V"), mapping=aes(x=mid/1e6, y=value)) +
geom_point(aes(colour = value), alpha=1, size=0.3) +
scale_colour_steps2('Tajimas D')+
labs(y="Tajima's D (100kb window)", x="Genomic position (Mb)")+
facet_grid(variable ~ chrom, scales="free_x", switch="x")+
scale_x_continuous(
    breaks = seq(0, 55, 5),
    minor_breaks = seq(1, 54, 1)
  )+
theme_light()+
theme(axis.title.x = element_text(size=10), axis.text.x = element_text(size=6), axis.title.y = element_text(size=10), axis.text.y = element_text(size=4), strip.text.x = element_text(size = 6, face = "bold", margin = margin(0.5,0,0.5,0, "mm")), strip.text.y = element_text(size = 6, face = "bold", angle=0), plot.title = element_text(size=10, hjust=0.5), legend.position = "top", legend.title = element_text(size=8), legend.text = element_text(size=6), legend.key.height=unit(0.2,"cm"), panel.spacing.x = unit(0.2, "mm"), panel.spacing.y = unit(0.75, "mm"))+
ggtitle('Ivermectin lines')

MOX<- ggplot(data=subset(mox, mox$chrom == "V"), mapping=aes(x=mid/1e6, y=value)) +
geom_point(aes(colour = value), alpha=1, size=0.3) +
scale_colour_steps2('Tajimas D')+
labs(y="Tajima's D (100kb window)", x="Genomic position (Mb)")+
facet_grid(variable ~ chrom, scales="free_x", switch="x")+
scale_x_continuous(
    breaks = seq(0, 55, 5),
    minor_breaks = seq(1, 54, 1)
  )+
theme_light()+
theme(axis.title.x = element_text(size=10), axis.text.x = element_text(size=6), axis.title.y = element_text(size=10), axis.text.y = element_text(size=4), strip.text.x = element_text(size = 6, face = "bold", margin = margin(0.5,0,0.5,0, "mm")), strip.text.y = element_text(size = 6, face = "bold", angle=0), plot.title = element_text(size=10, hjust=0.5), legend.position = "top", legend.title = element_text(size=8), legend.text = element_text(size=6), legend.key.height=unit(0.2,"cm"), panel.spacing.x = unit(0.2, "mm"), panel.spacing.y = unit(0.75, "mm"))+
ggtitle('Moxidectin lines')


CTL<- ggplot(data=subset(ctl, ctl$chrom == "V"), mapping=aes(x=mid/1e6, y=value)) +
geom_point(aes(colour = value), alpha=1, size=0.3) +
scale_colour_steps2('Tajimas D')+
labs(y="Tajima's D (100kb window)", x="Genomic position (Mb)")+
facet_grid(variable ~ chrom, scales="free_x", switch="x")+
scale_x_continuous(
    breaks = seq(0, 55, 5),
    minor_breaks = seq(1, 54, 1)
  )+
theme_light()+
theme(axis.title.x = element_text(size=10), axis.text.x = element_text(size=6), axis.title.y = element_text(size=10), axis.text.y = element_text(size=4), strip.text.x = element_text(size = 6, face = "bold", margin = margin(0.5,0,0.5,0, "mm")), strip.text.y = element_text(size = 6, face = "bold", angle=0), plot.title = element_text(size=10, hjust=0.5), legend.position = "top", legend.title = element_text(size=8), legend.text = element_text(size=6), legend.key.height=unit(0.2,"cm"), panel.spacing.x = unit(0.2, "mm"), panel.spacing.y = unit(0.75, "mm"))+
ggtitle('Control lines')

F0<- ggplot(data=subset(f0, f0$chrom == "V"), mapping=aes(x=mid/1e6, y=value)) +
geom_point(aes(colour = value), alpha=1, size=0.3) +
scale_colour_steps2('Tajimas D')+
labs(y="Tajima's D (100kb window)", x="Genomic position (Mb)")+
facet_grid(chrom ~ ., scales="free_x", switch="x")+
scale_x_continuous(
    breaks = seq(0, 55, 5),
    minor_breaks = seq(1, 54, 1)
  )+
theme_light()+
theme(axis.title.x = element_text(size=10), axis.text.x = element_text(size=6), axis.title.y = element_text(size=10), axis.text.y = element_text(size=4), strip.text.x = element_text(size = 6, face = "bold", margin = margin(0.5,0,0.5,0, "mm")), strip.text.y = element_text(size = 6, face = "bold", angle=0), plot.title = element_text(size=10, hjust=0.5), legend.position = "top", legend.title = element_text(size=8), legend.text = element_text(size=6), legend.key.height=unit(0.2,"cm"), panel.spacing.x = unit(0.2, "mm"), panel.spacing.y = unit(0.75, "mm"))+
ggtitle('F0 Tajimas D')

patch<- IVM + MOX + CTL + F0 + plot_annotation(tag_levels = 'A') + plot_layout (axis_titles = "collect")
patch

#ggsave("EmergingMLR_manuscript_FigX_light_titles_TajimasD_CHR5.tiff", patch, device=tiff, width=2250, height=2250, units="px", dpi=320)
ggsave("EmergingMLR_manuscript_FigX_light_titles_TajimasD_CHR5_withromannumerals.tiff", patch, device=tiff, width=2250, height=2250, units="px", dpi=320)

```
_Image with roman numerals:_
![image](https://github.com/user-attachments/assets/d83e5200-9e74-4c40-a3f8-ae36315ea4bb)


_Original image_
![image](https://github.com/user-attachments/assets/65b20964-b440-44d6-8dc2-424893a3c279)
