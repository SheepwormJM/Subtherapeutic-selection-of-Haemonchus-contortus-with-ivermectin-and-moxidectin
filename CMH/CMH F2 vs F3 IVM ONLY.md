Comparing F2 selected with F3 selected for IVM selection lines.

Notes:
- If we have the same snps/backgrounds under selection for both ivm and mox then could see a large number of snps for the ctl vs selected.
- If we do not have strong selection continuing between 1/16th and 1/10th dose of IVM (Fst peaks don't look that much bigger) then may find that we lose a lot of snps that are actually under selection by both drugs with the F2 vs F3 comparison. But equally, if there are enough genetic backgounds under selection around the causal snp(s) AND there is continued selection pressure with the increasing dose, or even if there is strong enough continued selection, then we may narrow the locus for the 1/10th doses.
- Note that this may not mean that lower doses do not also select for other snps.

# Calculate CMH statistics between pairwise populations with PoPoolation2

Using an updated version of the cmh.pl script found on sourceforge https://sourceforge.net/p/popoolation2/code/HEAD/tree/trunk/.

Obtained using by downloading and then using scp to transfer to Scamper.

Using Popoolation2 to run the CMH on the pool-seq as no CMH test yet available in Grenedalf.

Note that in PoPoolation2, unlike in GrENEDALF/0.2.0, the min-count applies considering all populations simultaneously. Note that in Grenedalf this is per sample. The PoPoolation2 says:
```
--min-count
        the minimum count of the minor allele. used for SNP identification.
        SNPs will be identified considering all populations simultanously.
        default=2
```
In my Fst and diversity calculations for GrENEDALF I used a per sample min count of 2, however for the CMH I will ask for a min count across all populations of 5. This may therefore include more SNPs than were included by GreneDALF. Note also that this is considering all my 29 populations in the sync file - not just those I am comparing.

Remember that the sync file now has the selection lines in the correct order.

### F2 generation selected sample paired with each of IVM and MOX F3 samples:
```
Generation F2:F3_IVM
13	16
14	17
15	18
```
Using all potential sites (i.e. not using the subsampled mpileup):

```
#!/bin/bash -l

############# SLURM SETTINGS #############
#SBATCH --job-name=cmh_F3_IVM_all_F2sel       # some descriptive job name of your choice
#SBATCH --output=%x-%J.out      # output file name will contain job name + job ID
#SBATCH --error=%x-%J.err       # error file name will contain job name + job ID
#SBATCH --time=0-10:00:00       # time limit for the whole run, in the form of d-hh:mm:ss, also accepts mm, mm:ss, hh:mm:ss, d-hh, d-hh:mm
#SBATCH --mem=10G                # memory required per node, in the form of [num][M|G|T]
#SBATCH --ntasks=1              # number of Slurm tasks to be launched, increase for multi-process runs ex. MPI
#SBATCH --cpus-per-task=1       # number of processor cores to be assigned for each task, default is 1, increase for multi-threaded runs
#SBATCH --ntasks-per-node=1     # number of tasks to be launched on each allocated node

############# LOADING MODULES (optional) #############
#module load apps/xxx
#module load libs/xxx

############# MY CODE #############

# Comparing F2 selected samples with F3 samples of selected ivm lines:

perl /home/jenni/popoolation2-code-r205-trunk/cmh-test.pl \
--input /home/jenni/data_folder/working_folder/selection_lines_poolseq/mpileup/Jan24_correct_order_syncs/SL_ps_only_Hcwbps18_q20Q30_noindel_lines_in_correct_order_Jan24.java.sync \
--output F2sel_vs_F3_IVM_ONLY_SLPS_Hcwbps18_q20Q30_noindels_popn2_updated_oct24.cmh \
--min-count 5 --min-coverage 10 --max-coverage 2% --population 13-16,14-17,15-18 --remove-temp
```

Notes on the output format using the new script - note it is no longer outputting a pvalue!
```
Input is a single tab delimited file which contains a lightwight representation of every pileup file.
Every pileup file represents a population and will be parsed into a list of A-count:T-count:C-count:G-count:N-count:*-count

 2L     5002    G       0:0:0:17:0:0    0:0:0:28:0:0    0:0:0:31:0:0    0:0:0:35:0:0    0:1:0:33:0:0    0:3:0:31:0:0
 2L     5009    A       16:0:0:0:0:0    26:0:0:0:0:0    29:0:1:0:0:0    36:0:0:0:0:0    34:0:0:0:0:0    32:0:1:0:0:0
 2L     5233    G       0:0:5:46:0:0    0:0:0:43:0:0    0:0:0:60:0:0    0:0:3:61:0:0    0:0:0:56:0:0    0:0:0:48:0:0
 
 col 1: reference chromosome
 col 2: position in the reference chromosome
 col 3: reference genome base
 col 4: population 1
 col 5: population 2
 col n: population n-3
 
 population data are in the form
 A:T:C:G:N:*
 A: count of character A
 T: count of character T
 C: count of character C
 G: count of character G
 N: count of character N
 *: deletion, count of deletion
 
=head2 Output

 2L     5002    G       0:0:0:17:0:0    0:0:0:28:0:0    0:0:0:31:0:0    0:0:0:35:0:0    0:1:0:33:0:0    0:3:0:31:0:0    0.609   0.1
 2L     5009    A       16:0:0:0:0:0    26:0:0:0:0:0    29:0:1:0:0:0    36:0:0:0:0:0    34:0:0:0:0:0    32:0:1:0:0:0    0.957   0.0
 2L     5233    G       0:0:5:46:0:0    0:0:0:43:0:0    0:0:0:60:0:0    0:0:3:61:0:0    0:0:0:56:0:0    0:0:0:48:0:0    0.8088  0.2


 col 1: reference chromosome
 col 2: position in the reference chromosome
 col 3: reference genome base
 col 4: population 1
 col 5: population 2
 col n: population n-3
 col n+1: cmh -log10(p-value)
 col n+2: log-odds ratio
 Note: If user gives --population  1-13,2-6,3-7 and --select-population 1,13,2,6,3,7 then SNP calling and p-value will be calculated only for selected populations but still all population will be printed in output file just to keep all sync file information.
```
Note therefore, that numbers < 1 are p values of >p=0.1, while numbers > 2 are p values of >p=0.01, and numbers > 1 are p values of >p=0.001 etc.

Ok, so now that I have the CMH file, what to do with it?

I will have some sites which have Nan for a pvalue. Remove these.
```
# First, remove the sites with Nan.

grep -v 'NaN' F2sel_vs_F3_IVM_ONLY_SLPS_Hcwbps18_q20Q30_noindels_popn2_updated_oct24.cmh > tmp.cmh
```
From now on use only those sites with a computed CMH value - i.e. for totals use those in the tmp.cmh file!
Adjust for multiple comparisons (i.e. loads of snps!)

1. Do an FDR - try 5%
```
wc -l *.cmh

# all sites allowed:
   18141828 F2sel_vs_F3_IVM_ONLY_SLPS_Hcwbps18_q20Q30_noindels_popn2_updated_oct24.cmh
   15602150 tmp.cmh
 
# So some sites had a SNP only in one comparison and not both. 

# Get the number of columns in the tmp file (to find out which column has the pvalue in it - the last one):
awk '{print NF}' tmp.cmh | uniq > nf_tmp

# Sort by the pvalue, in reverse
# sort -k 10 -nr tmp.cmh > tmp_sort # THIS HAS NOT WORKED.

sort -k 31 -r -g tmp.cmh > tmp_sort # This does work - it can use scientific numbers as input (ie. 9.9e-18 is known to be a larger number than 1.2e-30). If you also use -r then it goes from largest down to smallest...

# Add numbers to the file, increasing by one per line.
# %.0f means to print as a number without any numbers after the decimal point. ie. to print an integer.
# See here for more on printing numbers: https://bl831.als.lbl.gov/~gmeigs/scripting_help/printf_awk_notes.txt
# See here for more on doing this: https://www.gnu.org/software/gawk/manual/html_node/Format-Modifiers.html
# Note that the column is printed to the left.
awk '{printf("%.0f %s\n", (NR), $0)}' tmp_sort > tmp_sort_numbers
```
### Note - have not done an FDR for this dataset - just the plot towards the end.
Then to calculate the qvalue, using Benjamini-Hochberg proceedure:

https://www.statisticshowto.com/benjamini-hochberg-procedure/

First, as we have the -log10(pvalue) I will make a column containing the pvalue...
```
awk '{printf("%e %s\n", (10^-($32)), $0)}' tmp_sort_numbers > tmp_sort_pvalue
```
1. Put the individual p-values in ascending order.
2. Assign ranks to the p-values. For example, the smallest has a rank of 1, the second smallest has a rank of 2.
3. Calculate each individual p-value’s Benjamini-Hochberg critical value, using the formula (i/m)Q, where: i = the individual p-value’s rank, m = total number of tests, Q = the false discovery rate (a percentage, chosen by you, for example, if you wanted to allow a FDR of 25% then multiply by 0.25).
4. Compare your original p-values to the critical B-H from Step 3; find the largest p value that is smaller than the critical value.

# Allow an FDR of 5% 
```
# Get the qvalue to be in scientific notation (%e)

# For all sites:
awk '{printf("%e %s\n", ($2/18022431)*0.05, $0)}' tmp_sort_pvalue > tmp_sort_qvalue
```
Let's try using python to identify the largest pvalue that is smaller than the qvalue:

You want the file to be a csv (not tab delimited) for pandas to read it as a dataframe (https://pandas.pydata.org/docs/user_guide/io.html#io)
```
awk '{OFS="\t"}{ for (i=4; i<=NF; i++)   printf $i "\t"}{print $1, "\t", $2,"\t", $3}' tmp_sort_qvalue > tmp_new_sort_qvalue
# Add header - ensure use tabs
# chrom snp ref.allele F0 F1_CTL1 F1_CTL2 F1_CTL3 F2_CTL1 F2_CTL2 F2_CTL3 F3_CTL1 F3_CTL2 F3_CTL3 F1_IVM2 F1_IVM3 F2_IVM1 F2_IVM2 F2_IVM3 F3_IVM1 F3_IVM2 F3_IVM3 F1_MOX1 F1_MOX2 F1_MOX3 F2_MOX1 F2_MOX2 F2_MOX3 F3_MOX1 F3_MOX2 F3_MOX3 -log10p logoddsratio qvalue pvalue rank 
cat header tmp_new_sort_qvalue > tmp
sed 's/ /\t/g' tmp | sed 's/\t\t\t/\t/g' > tmp2 # The addition of the qvalue and rank to the cmh file was added using spaces rather than tabs to separate.
mv tmp2 tmp_sort_qvalue
sed 's/\t/,/g' tmp_sort_qvalue > tmp_sort_qvalue.csv

# First, install a module/package that I want (that is not already there):
python3 -m pip install pandas

# THEN...

python3
# Opens up a mini python console within scamper

# Importing Libraries
import pandas as pd
import numpy as np

df=pd.read_csv("tmp_sort_qvalue.csv") # Add ,header=None if no header. Although then found hard to manipulate later.
df.head

# Where method to compare the values
# The values were stored in the new column
# Note that this didn't work when I had no header (and so headers were 0,1 and 2).
df['qvaluecut'] = np.where((df['pvalue'] < df['qvalue']), df['rank'], np.nan)

# Ok, so it has worked. How to export it?
df.to_csv('F2F3_cmh_fdr_5pc.csv', sep='\t',na_rep='NaN', index=False)

# index=FALSE removes the rownames
```
Then to extract rows without NaN
```
# 'all sites':
grep -v 'NaN' F2F3_cmh_fdr_5pc.csv > grep_not_nan_F2F3_cmh_fdr_5pc.csv
less grep_not_nan_F2F3_cmh_fdr_5pc.csv

wc -l grep_not_nan_F2F3_cmh_fdr_5pc.csv
# 227087 grep_not_nan_F2F3_cmh_fdr_5pc.csv
```
So we have 227k SNPs, about 1.3% of all SNPs.

2. Do a Bonferroni correction

https://www.statology.org/bonferroni-correction/

If I use an **alpha value of 0.01** (i.e. pvalue of 0.01), then would do:

# Bonferroni correction to obtain new p value cut-off for multiple tests:

original a/total tests n = new a
 
0.01/18022431 = 5.54864102407e-10

And convert this to a -log10 value:

-log10(5.54864102407e-10)= 9.25581337

So, only accept SNPs as significant which have -log10(pvalue) of **9.26 or higher**.
```
awk 'BEGIN {FS=OFS="\t"} $31 > 9.26 {print $0}' F2F3_cmh_fdr_5pc.csv > BF_0.01_noNan_F2F3_cmh_fdr_5pc.csv

wc -l  BF_0.01_noNan_F2F3_cmh_fdr_5pc.csv
# 21735 BF_0.01_noNan_F2F3_cmh_fdr_5pc.csv
```
So **0.12% of SNPs** in the CMH noNAN file.

If I use an **alpha value of 0.05** (i.e. pvalue of 0.05), then would do:
```
0.05/18022431 = 2.77432051203e-9

And convert this to a -log10 value:

-log10(2.77432051203e-9)= 8.55684337

awk 'BEGIN {FS=OFS="\t"} $31 > 8.56 {print $0}' F2F3_cmh_fdr_5pc.csv > BF_0.05_noNan_F2F3_cmh_fdr_5pc.csv

wc -l  BF_0.05_noNan_F2F3_cmh_fdr_5pc.csv
# 27674 BF_0.05_noNan_F2F3_cmh_fdr_5pc.csv
```
So still approximately **0.15% of SNPs** in the CMH noNaN file. 


### I have done the below... 

# Ok, but how many of these SNPs are really interesting? I.e. are moving away from the reference allele following IVM treatment? 

Hypothesis: greater number of SNPs in ChrV locus that are in the top 1% of ranked SNPs than elsewhere in genome

To test this, we need to have the proportion of SNPs in the ChrV locus that are ranked between 1-82892, and ideally are changing in the correct direction.

And the proportion elsewhere etc.
```
Get just the chromosomal scaffolds individually:
for i in {1..5} ;
do
grep -e chr${i} tmp_sort_numbers > tmp_sort_numbers_chr_${i} ;
done &

grep -e 'chrX' tmp_sort_numbers > tmp_sort_numbers_chr_X &


# And get the SNPs for each 500 kb window.
# awk -v assigns a value to a variable before awk starts
for FILE in tmp_sort_numbers_chr_* ; 
do
for i in $(seq 1 500000 52000000); do
    cat ${FILE} |\
    awk -v pos=$i 'BEGIN {FS=OFS="\t"} {if ($2>=pos && $2<pos+500000) print $0, "locus_"pos"_"pos+499999}' > tmp_snps_${FILE}_${i} ;
    done ;
done &
```
2. Calculate the number of SNPs in each Mb. And split the filename to get chr and start and end points of windows:
```
wc -l tmp_snps* | sed 's/tmp_snps_tmp_sort_numbers_//g' | sed 's/_/\t/g' >> SNPS_SUM_500KB.txt
grep -v 'total' SNPS_SUM_500KB.txt > tmp
mv tmp SNPS_SUM_500KB.txt
# Note, had checked the above without the grep -v part and the total snps on the chromosomes was the same as the total in the noNaN cmh file. 

# Get rid of the empty windows:
awk 'BEGIN {OFS="\t"} $1 > 0 {print $0}' SNPS_SUM_500KB.txt | sort -n -k 1 -r > SNPS_SUM_500KB.txt_tmp
```
3. Find out how many of the SNPs are in the top 1pc rank AND have a log odds ratio suggesting a reduction in the reference allele in the post-treatment populations relative to the control populations (i.e. logoddsratio < 1).
```
for i in tmp_snps* ;
do
awk 'BEGIN {OFS="\t"} $1 < 156021 && $33 < 1 {print $0, FILENAME}' ${i} >> SNPS_500KB_lorless1_top1pc.txt ;
done
```
How many SNPs are in the top 1pc and have a log odds ratio < 1 on the chromosomes?
```
# all sites:
wc -l SNPS_500KB_lorless1_top1pc.txt
#  63836 SNPS_500KB_lorless1_top1pc.txt
```
Then you want to find the number of SNPs are in the top 1pc and have a log odds ratio < 1 per 500Kb window.
```
# Print the last two columns in the file ('locus_24000001_24500000' and 'tmp_snps_tmp_sort_numbers_chr_5_24000001' i.e. the locus and the filename inserted above - you want both as only one has the chromosome info and the other has the locus start and end!)
cat SNPS_500KB_lorless1_top1pc.txt | awk 'BEGIN {OFS="\t"} {print $34, $35}' | uniq -c > uniq_count_SNPS_500KB_lorless1_top1pc.txt
# Worked fine
```
Get the tables ready to combine together in R:
```
# 122 locus_1_1000000 tmp_snps_tmp_sort_numbers_chr_1_1
cat uniq_count_SNPS_500KB_lorless1_top1pc.txt | sed 's/tmp_snps_tmp_sort_numbers_//g' | sed 's/locus_//g' | sed 's/_/\t/g' | sed 's/ \+ /\t/g' | sed 's/ /\t/g' | sed 's/^\t//g'> uniq_count_SNPS_500KB_lorless1_top1pc.txt_tmp

#      70139 chr  X       31000001
cat SNPS_SUM_500KB.txt_tmp | sed 's/ \+ /\t/g' | sed 's/ /\t/g' | sed 's/^\t//g' | awk 'BEGIN {OFS="\t"} {print $0, ($4 + 499999)}' > SNPS_SUM_500KB.txt_tmp2
```
Note, I then removed the last 500kb window from each chromosome as it would not be a complete 500kb and could produce odd proportions/results.
```
# First, I reverse sorted the file, then I used nano to manually edit it:
sort -r -n -k 6 uniq_count_SNPS_500KB_lorless1_top1pc.txt_tmp > uniq_count_SNPS_500KB_lorless1_top1pc.txt_tmp_sorted

sort -r -n -k 5 SNPS_SUM_500KB.txt_tmp2 > SNPS_SUM_500KB.txt_tmp2_sorted
```
```
hcontortus_chr4_Celeg_TT_arrow_pilon    51826579        54      60      61
hcontortus_chr5_Celeg_TT_arrow_pilon    48868368        52690464        60      61
hcontortus_chr2_Celeg_TT_arrow_pilon    47382676        102373359       60      61
hcontortus_chrX_Celeg_TT_arrow_pilon    46012677        150545801       60      61
hcontortus_chr1_Celeg_TT_arrow_pilon    45774638        197325410       60      61
hcontortus_chr3_Celeg_TT_arrow_pilon    43560352      
```
4. Combine the tables in R to then calculate the proportion of SNPs within the window that are in the top 1pc and with the probability that the refallele is lower in the post-tx than the pre-tx sample:

```
# R version 4.3.1

library(dplyr)
library(ggplot2)
library(patchwork)

df<-read.table("SNPS_SUM_500KB.txt_tmp2_sorted", header=FALSE, sep="\t")
df2<-read.table("uniq_count_SNPS_500KB_lorless1_top1pc.txt_tmp_sorted", header=FALSE, sep="\t")

head(df)
head(df2)

# Need to make the chromosome number in the uniq_count file as a character - no ChrX SNPs unlike in the SNPS file - and inner_join() needs them to be the same type to join.
df2$V5 <- as.character(df2$V5)
# Join tables, keeping only rows present in both tables (this will account for removal of windows with no SNPs
df3<-inner_join(df, df2, by=c("V3"="V5", "V4"="V2", "V5"="V3"))

# Calculate the proportion of SNPs in each window that are in the top 1pc by rank, and have log odds ratio < 1:
df3$prop<-(df3$V1.y/df3$V1.x)

summary(df3)
```
```
      V1.x            V2                 V3                  V4          
 Min.   : 5631   Length:564         Length:564         Min.   :       1  
 1st Qu.:22578   Class :character   Class :character   1st Qu.:11500001  
 Median :28932   Mode  :character   Mode  :character   Median :23250001  
 Mean   :27586                                         Mean   :23317377  
 3rd Qu.:33636                                         3rd Qu.:35000001  
 Max.   :44157                                         Max.   :51000001  
       V5                V1.y            V4.y                 V6          
 Min.   :  500000   Min.   :   1.0   Length:564         Min.   :       1  
 1st Qu.:12000000   1st Qu.:  19.0   Class :character   1st Qu.:11500001  
 Median :23750000   Median :  42.0   Mode  :character   Median :23250001  
 Mean   :23817376   Mean   : 113.1                      Mean   :23317377  
 3rd Qu.:35500000   3rd Qu.:  87.0                      3rd Qu.:35000001  
 Max.   :51500000   Max.   :2079.0                      Max.   :51000001  
      prop          
 Min.   :3.985e-05  
 1st Qu.:7.752e-04  
 Median :1.427e-03  
 Mean   :4.169e-03  
 3rd Qu.:2.780e-03  
 Max.   :7.420e-02  
```
```
# Plot the data:

rect<-data.frame(xmin = c(NA, NA, NA, NA, 30, NA), xmax = c(NA,NA,NA,NA,45,NA), 
                 ymin = c(NA, NA, NA, NA, 0, NA), ymax = c(NA, NA, NA, NA, 0.15, NA), 
                 alpha = c(NA, NA, NA, NA, 0.2, NA),
                 fill = c(NA, NA, NA, NA, "red", NA))

x<-ggplot(data=df3, aes(x=((V4 + 249999)/1e6), y=(prop)))+
  geom_line(linewidth=0.8)+
  #geom_point(size=0.5)+
geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, alpha=alpha, fill=fill),
            data = transform(rect, V3 = as.character(5)),
            inherit.aes = FALSE)+
facet_grid(.~V3, scales = "free_x", switch="x")+
labs(y="Proportion of SNPs per 500 Kb window that are in top 1%\nby rank of CMH p-values, with log odds ratio < 1", x= "Basepair (Mb)", title="Proportion of SNPs in 500 Kb windows along nuclear chromosomes that are in the top 1% CMH SNPs by rank,\nand have log odds ratio < 1 (probability that ref allele reduced in F3 sample cf F2 sample)")+
   scale_x_continuous(
 breaks = seq(0, 55, 10),
         minor_breaks = seq(5, 50, 5)
     )+
theme_light()+
theme(axis.title.x = element_text(size=10), axis.text.x = element_text(size=6), axis.title.y = element_text(size=10), axis.text.y = element_text(size=6), strip.text.x = element_text(size = 6, face = "bold", margin = margin(0.5,0,0.5,0, "mm")), strip.text.y = element_text(size = 6, face = "bold", angle=0), plot.title = element_text(size=10, hjust=0.5), legend.position = "none", panel.spacing.x = unit(0.2, "mm"), panel.spacing.y = unit(0.75, "mm"))

x

# Plot just ChrV:
y<-ggplot(data=subset(df3, df3$V3 == 5), aes(x=((V4 + 249999)/1e6), y=(prop)))+
  #geom_line(linewidth=0.8)+
  geom_point(size=2)+
geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, alpha=alpha, fill=fill),
            data = transform(rect, V3 = as.character(5)),
            inherit.aes = FALSE)+
facet_grid(.~V3, scales = "free_x", switch="x")+
labs(y="Proportion of SNPs per 500 Kb window that are in top 1%\nby rank of CMH p-values, with log odds ratio < 1", x= "Basepair (Mb)", title="Proportion of SNPs in 500 Kb windows along nuclear chromosomes that are in the top 1% CMH SNPs by rank,\nand have log odds ratio < 1 (probability that ref allele reduced in F3 sample cf F2 sample)")+
   scale_x_continuous(
 breaks = seq(0, 55, 5),
         minor_breaks = seq(1, 54, 1)
     )+
theme_light()+
theme(axis.title.x = element_text(size=10), axis.text.x = element_text(size=6), axis.title.y = element_text(size=10), axis.text.y = element_text(size=6), strip.text.x = element_text(size = 6, face = "bold", margin = margin(0.5,0,0.5,0, "mm")), strip.text.y = element_text(size = 6, face = "bold", angle=0), plot.title = element_text(size=10, hjust=0.5), legend.position = "none", panel.spacing.x = unit(0.2, "mm"), panel.spacing.y = unit(0.75, "mm"))

y

patch <- x/y
patch2 <- patch + plot_annotation(tag_levels = 'A') + plot_layout(guides = 'collect', axis_titles = 'collect', axes = 'collect') & theme(plot.tag = element_text(size = 20))
patch2

ggsave("HcSLPS_SNP_PROP_CMH_Chr_500kb_point_patch_F2_vs_F3_IVMonly.tiff", width=17, height=10)
```
![image](https://github.com/user-attachments/assets/6933b775-89da-417b-85dd-239ca1e1715a)

_Above: F2 vs F3 generation (only IVM)_
