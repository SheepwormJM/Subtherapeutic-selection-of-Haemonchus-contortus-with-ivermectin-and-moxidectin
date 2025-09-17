Was asked by a reviewer about the reference allele for F0 - was it normally the major allele? 

Note that from the AF model we have snps which increase in the reference allele as well as those that reduce.

Note that the AF model in the submitted manuscript used only those SNPs between 30-45 Mb that had a reference allele frequency of <0.9 in at least one of the F3 selected samples. 

```
# Working from the file: ChrV_locus_AC_depth_allsites_NINE_REPLICATES_refAFless0.9.txt
# pwd = /home/jenni/data_folder/working_folder/selection_lines_poolseq/pairwise_comparisons/LINES_IN_CORRECT_ORDER/grenedalf/Jan24_F0_vs_all_sel_lines_only/all_sites/AF/wholegenome_depth

# This file has the reference allele count, the 'alt' allele count and the depth (calculated from the mpileup) for all SNPs within the ChrV:30-45 Mb region that have a reference allele frequency of <0.9 in at least one F3 selected sample. Therefore, any SNPs that have a reference allele frequency of >0.9 in all F3 selected samples, or have NA in all of the selected F3 samples have already been discarded.

cut -f 1-2,5,59 ChrV_locus_AC_depth_allsites_NINE_REPLICATES_refAFless0.9.txt > F0_V30-45.txt

awk 'BEGIN {OFS=FS="\t"} {print $0, ($3/$4)}' F0_V30-45.txt > F0_V30-45_AF.txt

# Get those with the ref AF  >0.5 and the total number of SNPs we are considering:

wc -l F0_V30-45_AF.txt > TOTAL_SNPS_CONSIDERING
awk 'BEGIN {OFS=FS="\t"} $5 >0.5 {print $0}' F0_V30-45_AF.txt > F0_V30-45_AF_RefMajor.txt
wc -l F0_V30-45_AF_RefMajor.txt > TOTAL_SNPS_REFMAJOR_F0

# 938316 F0_V30-45_AF.txt
# 815427 F0_V30-45_AF_RefNotMajor.txt

Therefore 86.9% of SNPs we have included in the model have the major allele as the reference allele in the F0 sample. 
```
To show where these are all found.... 
```
> df<-read.table("F0_V30-45_AF.txt", header=T)
> head(df)
  CHROM      POS F0.REF_CNT F0.DEPTH    X.nan
1  chr5 30000139         10       20 0.500000
2  chr5 30001031         19       27 0.703704
3  chr5 30001033         20       29 0.689655
4  chr5 30001056         17       27 0.629630
5  chr5 30001060         20       30 0.666667
6  chr5 30001092         21       21 1.000000
> library(ggplot2)
> plot<-ggplot(data=df, aes(x=POS/1e6, y=X.nan))+
+ geom_point(alpha=0.5, size=0.3)+
+ geom_hline(yintercept=0.5)
> ggsave("REF_AF_F0.tiff", plot)
```
<img width="2100" height="2100" alt="image" src="https://github.com/user-attachments/assets/f76c8bf0-616c-40bd-a8f3-927619b4f2e4" />

```
library(ggplot2)
library(viridis)

df<-read.table("F0_V30-45_AF.txt", header=T)
head(df)
plot<-ggplot(data=df, aes(x=POS/1e6, y=X.nan))+
labs(y="Reference allele frequency", x="Genomic position (Mb)")+
geom_hex(binwidth=c(0.1,0.01)) +
scale_fill_viridis()+
scale_x_continuous(
    breaks = seq(0, 55, 1),
    minor_breaks = seq(1, 55, 1)
  )+
geom_hline(yintercept=0.5)
ggsave("REF_AF_F0_geomhex.tiff", plot)
```
<img width="2100" height="2100" alt="image" src="https://github.com/user-attachments/assets/b46210a4-ff6d-4386-b940-fddc232425c9" />


## Determining the reference allele frequency for F0 for those SNPs that were excluded from the model because they had >0.9 AF for the reference allele in all F3 IVM and F3 MOX samples.

This is in case these are 'near fixation' in the selected lines and have an alternate allele as the major allele in F0. 

```
# First, need to calculate the AF for the reference allele for each of the samples:

module load R/4.4.2
R
df<-read.table("Chr5_30-45Mb_AC_depth_allsites.txt", header=T)

df$F3_I_L1_RAF<-(df$F3_I_L1.REF_CNT / df$F3_I_L1.DEPTH)
df$F3_I_L2_RAF<-(df$F3_I_L2.REF_CNT / df$F3_I_L2.DEPTH)
df$F3_I_L3_RAF<-(df$F3_I_L3.REF_CNT / df$F3_I_L3.DEPTH)

df$F3_M_L1_RAF<-(df$F3_M_L1.REF_CNT / df$F3_M_L1.DEPTH)
df$F3_M_L2_RAF<-(df$F3_M_L2.REF_CNT / df$F3_M_L2.DEPTH)
df$F3_M_L3_RAF<-(df$F3_M_L3.REF_CNT / df$F3_M_L3.DEPTH)

df$F0_RAF<-(df$F0.REF_CNT / df$F0.DEPTH)

summary(df)
# Note, have NA where have a depth of zero for a site

write.table(df, file="F0_AND_F3refAF_Chr5_30-45Mb_AC_depth_allsites.txt", row.names=FALSE, sep="\t", quote= FALSE)
q()
n


# Then will use python to identify which lines to keep and which to filter:

# First, make the input file a csv file:
sed 's/\t/,/g' F0_AND_F3refAF_Chr5_30-45Mb_AC_depth_allsites.txt > F0_AND_F3refAF_Chr5_30-45Mb_AC_depth_allsites.csv


python3
# Opens up a mini python console within scamper

# Importing Libraries
import pandas as pd
import numpy as np

df=pd.read_csv("F0_AND_F3refAF_Chr5_30-45Mb_AC_depth_allsites.csv") # Add ,header=None if no header. Although then found hard to manipulate later.
df.head


# Use the 'Where' method to compare the values
# The values were stored in the new column
# Note that this didn't work when I had no header (and so headers were 0,1 and 2 etc).

# Identify sites where in all F3 selected samples, the reference allele frequency is < 0.9. Label these "keep" - these were the SNPs used in the AF model in the paper
df['AF_F3'] = np.where(((df['F3_I_L1_RAF'] < 0.9) & (df['F3_I_L1_RAF'] != 'NA')) | ((df['F3_I_L2_RAF'] < 0.9) & (df['F3_I_L2_RAF'] != 'NA')) | ((df['F3_I_L3_RAF'] < 0.9) & (df['F3_I_L3_RAF'] != 'NA')) | ((df['F3_M_L1_RAF'] < 0.9) & (df['F3_M_L1_RAF'] != 'NA')) | ((df['F3_M_L2_RAF'] < 0.9) & (df['F3_M_L2_RAF'] != 'NA')) | ((df['F3_M_L3_RAF'] < 0.9) & (df['F3_M_L3_RAF'] != 'NA')), 'keep', np.nan)

# Create four new columns identifying whether the F0 reference allele frequency is <0.9, 0.75, 0.5 or <0.25. 
df['AF_F0_9'] = np.where(((df['F0_RAF'] < 0.9) & (df['F0_RAF'] != 'NA')), 'F0_l0.9_REF', np.nan)
df['AF_F0_75'] = np.where(((df['F0_RAF'] < 0.75) & (df['F0_RAF'] != 'NA')), 'F0_l0.75_REF', np.nan)
df['AF_F0_5'] = np.where(((df['F0_RAF'] < 0.5) & (df['F0_RAF'] != 'NA')), 'F0_l0.5_REF', np.nan)
df['AF_F0_25'] = np.where(((df['F0_RAF'] < 0.25) & (df['F0_RAF'] != 'NA')), 'F0_l0.25_REF', np.nan)


# Save as new table:

df.to_csv('F0_marked_Chr5_30-45Mb_AC_depth_allsites.csv', sep='\t',na_rep='NaN', index=False)
# Checked that it will keep any site with at ref AF<0.9 in at least one F3 selected sample.

# Then, want to keep only those rows withOUT 'keep' and look at what proportion of these have low reference allele frequencies in the F0 (i.e. the reference allele could be resistant)

head -n 1 F0_marked_Chr5_30-45Mb_AC_depth_allsites.csv > header
grep -e 'keep' F0_marked_Chr5_30-45Mb_AC_depth_allsites.csv > F0_keep_Chr5_30-45Mb_AC_depth_allsites.csv
grep -v 'keep' F0_marked_Chr5_30-45Mb_AC_depth_allsites.csv > NOT_keep_Chr5_30-45Mb_AC_depth_allsites.csv

grep -e 'F0_l0.25_REF' NOT_keep_Chr5_30-45Mb_AC_depth_allsites.csv > F0_l0.25_REF_NOT_keep_Chr5_30-45Mb_AC_depth_allsites.csv
grep -e 'F0_l0.5_REF' NOT_keep_Chr5_30-45Mb_AC_depth_allsites.csv > F0_l0.5_REF_NOT_keep_Chr5_30-45Mb_AC_depth_allsites.csv
grep -e 'F0_l0.75_REF' NOT_keep_Chr5_30-45Mb_AC_depth_allsites.csv > F0_l0.75_REF_NOT_keep_Chr5_30-45Mb_AC_depth_allsites.csv
grep -e 'F0_l0.9_REF' NOT_keep_Chr5_30-45Mb_AC_depth_allsites.csv > F0_l0.9_REF_NOT_keep_Chr5_30-45Mb_AC_depth_allsites.csv

wc -l F0_l*
wc -l NOT_*
wc -l F0_keep*
wc -l Chr5_30-45Mb_AC_depth_allsites.txt
     177 F0_l0.25_REF_NOT_keep_Chr5_30-45Mb_AC_depth_allsites.csv
     934 F0_l0.5_REF_NOT_keep_Chr5_30-45Mb_AC_depth_allsites.csv
   18205 F0_l0.75_REF_NOT_keep_Chr5_30-45Mb_AC_depth_allsites.csv
   93920 F0_l0.9_REF_NOT_keep_Chr5_30-45Mb_AC_depth_allsites.csv
 3684953 NOT_keep_Chr5_30-45Mb_AC_depth_allsites.csv   # Those SNPs not used in the model 


  938315 F0_keep_Chr5_30-45Mb_AC_depth_allsites.csv  # Those used for the model (REF AF <0.9 and not NA in all F3 selected samples)
4623268 Chr5_30-45Mb_AC_depth_allsites.txt  # All SNPs in the locus

```
