Needs updated still... 

Using Grenedalf, calculate allele counts for the reference and the most common alternate allele for each site. 

Do this for all sites, any depth. 

Note that the allele frequency is calculated from the sum of the ref and alt alleles by Grenedalf. So there could be an elevated ref allele frequency if a different alt allele was present in the data, as Grenedalf selects the most common alternate allele across samples to report. An alternative would be to get the ref allele count and the total depth count for each site and then divide by that. 

## Using all sites, all depths

1. Run Grenedalf on the mpileup file to obtain allele frequencies for all sites.

Extract just the sites I am interested in - those from V:30-45Mb from the sync file:

(pwd:Scamper,  /home/jenni/data_folder/working_folder/selection_lines_poolseq/mpileup/Jan24_correct_order_syncs)
```
grep -e 'chr5' SL_ps_only_Hcwbps18_q20Q30_noindel_lines_in_correct_order_Jan24.java.sync | awk 'BEGIN {FS=OFS="\t"} $2 > 30000000 && $2 < 45000000 {print $0}' > Chr5_locus_sellinesonly_Hcwbps18_q20Q30_noindel_lines_in_correct_order_Jan24.java.sync
```
Make the sample names file: 
```
sed 's/SLps_only_co_Hcwbps18_q20Q30_noindels_ss57/Chr5_locus_sellinesonly_Hcwbps18_q20Q30_noindel_lines_in_correct_order_Jan24/g' SL_ss57_sample_names.txt > SL_sample_names.txt
```

```
#!/bin/bash -l

############# SLURM SETTINGS #############
#SBATCH --job-name=AF_allsites       # some descriptive job name of your choice
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

module load grenedalf/0.3.0

############# MY CODE #############
# pwd = /home/jenni/data_folder/working_folder/selection_lines_poolseq/pairwise_comparisons/LINES_IN_CORRECT_ORDER/grenedalf/Jan24_F0_vs_all_sel_lines_only/all_sites/AF

# Calculate the allele counts (ref and main alt) for all sites from V:30-45 Mb

grenedalf frequency \
--sync-path /home/jenni/data_folder/working_folder/selection_lines_poolseq/mpileup/Jan24_correct_order_syncs/Chr5_locus_sellinesonly_Hcwbps18_q20Q30_noindel_lines_in_correct_order_Jan24.java.sync \
--reference-genome-fasta-file /home/jenni/data_folder/wbps18/haemonchus_contortus.PRJEB506.WBPS18.genomic.fa \
--rename-samples-file /home/jenni/data_folder/working_folder/selection_lines_poolseq/pairwise_comparisons/LINES_IN_CORRECT_ORDER/grenedalf/Jan24_F0_vs_all_sel_lines_only/all_sites/AF/SL_sample_names.txt \
--write-sample-counts \
--separator-char tab \
--out-dir /home/jenni/data_folder/working_folder/selection_lines_poolseq/pairwise_comparisons/LINES_IN_CORRECT_ORDER/grenedalf/Jan24_F0_vs_all_sel_lines_only/all_sites/AF \
--file-prefix Chr5_30-45Mb_400_Gf_F0_vs_SL_only_Hcwbps18_q20Q30_noindels_AF_sample_counts \
--log-file grenedalf_AF_all.log \
2>error_log_all_af &
```

2. Calculate the total depth per site.

Note that using samtools depth to do this outputs a different total depth from the mpileup sync file that Grenedalf uses as input. Therefore, will calculate the depth from the mpileup sync file.

This script will do the following: 
- It will loop through the columns 4 to 30 (in the sync file there are 30 columns, the first three are concerned with the position and the SNP allele, the rest have the per sample data
- It will then print each column, print the first 4 bases, and the deletions from that column and finally it will sum the counts of those 4 bases and deletions (so it ignores N's).
- NICE! :D 

```
# awk -v assigns a value to a variable before awk starts
for i in {4..30}; 
do
    cat /home/jenni/data_folder/working_folder/selection_lines_poolseq/mpileup/Jan24_correct_order_syncs/Chr5_locus_sellinesonly_Hcwbps18_q20Q30_noindel_lines_in_correct_order_Jan24.java.sync | \
    awk -v col=$i 'BEGIN {FS=OFS="\t"} {print $col }' | awk 'BEGIN {FS=OFS=":"} {print $1,$2,$3,$4,$6}' | sed 's/\:/+/g' | bc > sum_col${i} ;
done &
```
Right. Next - use paste to put the columns back together. 
```
cut -f 1-3 /home/jenni/data_folder/working_folder/selection_lines_poolseq/mpileup/Jan24_correct_order_syncs/Chr5_locus_sellinesonly_Hcwbps18_q20Q30_noindel_lines_in_correct_order_Jan24.java.sync > positions

paste positions sum_col4 sum_col5 sum_col6 sum_col7 sum_col8 sum_col9 sum_col10 sum_col11 sum_col12 sum_col13 sum_col14 sum_col15 sum_col16 sum_col17 sum_col18 sum_col19 sum_col20 sum_col21 sum_col22 sum_col23 sum_col24 sum_col25 sum_col26 sum_col27 sum_col28 sum_col29 sum_col30 > from_mpileupfile_SL_only_Dec_24_Chr5_30-45Mb.ACTGdepth
```

Next, remove the second reference allele column from the file - this would mess up the R script later.
```
cut -f 4- from_mpileupfile_SL_only_Dec_24_Chr5_30-45Mb.ACTGdepth > tmp
cut -f 1-2 from_mpileupfile_SL_only_Dec_24_Chr5_30-45Mb.ACTGdepth > tmp2

paste tmp2 tmp > tmp3

mv tmp3 from_mpileupfile_SL_only_Dec_24_Chr5_30-45Mb.ACTGdepth
```

3. Combine the grenedalf allele count file and the depth per site file.

```
module load R/4.3.1
R
```
Then, in R: 
```
library(dplyr)

df<-read.table("Chr5_30-45Mb_400_Gf_F0_vs_SL_only_Hcwbps18_q20Q30_noindels_AF_sample_countsfrequency.csv", header=TRUE, sep="\t")
df2<-read.table("from_mpileupfile_SL_only_Dec_24_Chr5_30-45Mb.ACTGdepth", header=FALSE, sep="\t")

head(df)
head(df2)

# Join tables, keeping only rows present in both tables (this will account for removal of windows with no SNPs
df3<-inner_join(df, df2, by=c("CHROM"="V1", "POS"="V2"))

head(df3)

write.table(df3, file="Chr5_30-45Mb_AC_depth_allsites.txt", row.names=FALSE, sep="\t", quote= FALSE) 
```
Then, attach the header.

Note that to get it to work well with the R script I am adding .DEPTH to each of the depth column headers
```
ln -s ~/data_folder/selection_lines_poolseq_bams/header_depth_sellinesonly ./header_depth_sellinesonly

sed 's/\t/.DEPTH\t/g' header_depth_sellinesonly | sed 's/$/.DEPTH/g' > header_depth_sellinesonly_DEPTH

cut -f 3- header_depth_sellinesonly_DEPTH > tmp

paste header_allele_count tmp > tmp2

cat tmp2 Chr5_30-45Mb_AC_depth_allsites.txt > tmp3

# Then I removed the second old header line using nano.

mv tmp3 Chr5_30-45Mb_AC_depth_allsites.txt

wc -l Chr5_30-45Mb_AC_depth_allsites.txt
# 4623268 Chr5_30-45Mb_AC_depth_allsites.txt
```

Am struggling to get it to run due to size on my laptop and installation of libraries on scamper and access on MARS. 

So, will try removing all snps that have a reference allele frequency of > 90% in all three F3 IVM And F3 MOX lines. These will essentially have only a low alternate allele frequency (whether it is the major AC or a different AC and are therefore SNPs that are unlikely to be under selection. As the model runs independently for each SNP removing them should not matter - but just remember they have been thinned when looking at data. 

First, need to calculate the AF for the reference allele for each of the samples:
```
module load R/4.4.2
R
df<-read.table("Chr5_30-45Mb_AC_depth_allsites.txt", header=T)

df$F3_I_L1_RAF<-(df$F3_I_L1.REF_CNT / df$F3_I_L1.DEPTH)
df$F3_I_L2_RAF<-(df$F3_I_L2.REF_CNT / df$F3_I_L2.DEPTH)
df$F3_I_L3_RAF<-(df$F3_I_L3.REF_CNT / df$F3_I_L3.DEPTH)

df$F3_M_L1_RAF<-(df$F3_M_L1.REF_CNT / df$F3_M_L1.DEPTH)
df$F3_M_L2_RAF<-(df$F3_M_L2.REF_CNT / df$F3_M_L2.DEPTH)
df$F3_M_L3_RAF<-(df$F3_M_L3.REF_CNT / df$F3_M_L3.DEPTH)

summary(df)
# Note, have NA where have a depth of zero for a site

write.table(df, file="F3refAF_Chr5_30-45Mb_AC_depth_allsites.txt", row.names=FALSE, sep="\t", quote= FALSE)
```

Then will use python to identify which lines to keep and which to filter:
```
# First, make the input file a csv file:
sed 's/\t/,/g' F3refAF_Chr5_30-45Mb_AC_depth_allsites.txt > F3refAF_Chr5_30-45Mb_AC_depth_allsites.csv


python3
# Opens up a mini python console within scamper

# Importing Libraries
import pandas as pd
import numpy as np

df=pd.read_csv("F3refAF_Chr5_30-45Mb_AC_depth_allsites.csv") # Add ,header=None if no header. Although then found hard to manipulate later.
df.head


# Use the 'Where' method to compare the values
# The values were stored in the new column
# Note that this didn't work when I had no header (and so headers were 0,1 and 2 etc).
df['AF_F3'] = np.where(((df['F3_I_L1_RAF'] < 0.9) & (df['F3_I_L1_RAF'] != 'NA')) | ((df['F3_I_L2_RAF'] < 0.9) & (df['F3_I_L2_RAF'] != 'NA')) | ((df['F3_I_L3_RAF'] < 0.9) & (df['F3_I_L3_RAF'] != 'NA')) | ((df['F3_M_L1_RAF'] < 0.9) & (df['F3_M_L1_RAF'] != 'NA')) | ((df['F3_M_L2_RAF'] < 0.9) & (df['F3_M_L2_RAF'] != 'NA')) | ((df['F3_M_L3_RAF'] < 0.9) & (df['F3_M_L3_RAF'] != 'NA')), 'keep', np.nan)

# Save as new table:

df.to_csv('marked_Chr5_30-45Mb_AC_depth_allsites.csv', sep='\t',na_rep='NaN', index=False)
```
Checked that it will keep any site with at ref AF<0.9 in at least one F3 selected sample. 

Then, want to keep only those rows with 'keep' 
```
head -n 1 marked_Chr5_30-45Mb_AC_depth_allsites.csv > header
grep -e 'keep' marked_Chr5_30-45Mb_AC_depth_allsites.csv > keep_Chr5_30-45Mb_AC_depth_allsites.csv
cat header keep_Chr5_30-45Mb_AC_depth_allsites.csv > header_keep_Chr5_30-45Mb_AC_depth_allsites.csv
```
Have removed a lot of rows!
```
4623268 Chr5_30-45Mb_AC_depth_allsites.txt
938315 keep_Chr5_30-45Mb_AC_depth_allsites.csv
```
So, to make the file smaller still... remove what I don't need - F0 and RAF columns:
```
awk '{print NF}' header_keep_Chr5_30-45Mb_AC_depth_allsites.csv | uniq > nf_file
# Have 92 columns in total. Therefore...
cut -f 1-4 header_keep_Chr5_30-45Mb_AC_depth_allsites.csv > first
cut -f 7-85 header_keep_Chr5_30-45Mb_AC_depth_allsites.csv > second
paste first second > tmp
mv tmp noF0_keep_Chr5_30-45Mb_AC_depth_allsites.csv
```
Ok, the file size is now only 247 Mb. However, still have a fair number of SNPs (<1M, but not by much). To see if it will work in R on my laptop. If not then will split into groups of snps. 

First, remove the F0.depth column! I forgot about this: 
```
cut -f 1-56 noF0_keep_Chr5_30-45Mb_AC_depth_allsites.csv > first
cut -f 58-83 noF0_keep_Chr5_30-45Mb_AC_depth_allsites.csv > second
paste first second > noF0_keep_Chr5_30-45Mb_AC_depth_allsites.csv
```
```
scp -r jenni@scamper.mvls.gla.ac.uk:/home/jenni/data_folder/working_folder/selection_lines_poolseq/pairwise_comparisons/LINES_IN_CORRECT_ORDER/grenedalf/Jan24_F0_vs_all_sel_lines_only/all_sites/AF/noF0_keep_Chr5_30-45Mb_AC_depth_allsites.csv ./noF0_keep_Chr5_30-45Mb_AC_depth_allsites.csv
```

Debugging.. 
```
# in python3 try to extract snps which have missing data for at least one sample - i.e. no depth of coverage. Using the 1-99999 file of data with F0 columns removed.
df['DEPTH_ZERO'] = np.where(((df['F3_I_L1.DEPTH'] == 0)) | ((df['F3_I_L2.DEPTH'] == 0)) | ((df['F3_I_L3.DEPTH'] == 0 )) | ((df['F3_M_L1.DEPTH'] == 0)) | ((df['F3_M_L2.DEPTH'] == 0 )) | ((df['F3_M_L3.DEPTH'] == 0 )), 'keep', np.nan)
df.to_csv('marked_tmp.csv', sep='\t',na_rep='NaN', index=False)
```

Rename the lines so that have nine replicates for the model (each line a different replicate) rather than grouping by treatment: 
```
# Changed the header line using sed 

sed 's/I_L1/I_L4/g' noF0_keep_Chr5_30-45Mb_AC_depth_allsites.csv | sed 's/I_L2/I_L5/g' | sed 's/I_L3/I_L6/g' | sed 's/M_L1/M_L7/g' | sed 's/M_L2/M_L8/g' | sed 's/M_L3/M_L9/g' > noF0_keep_Chr5_30-45Mb_AC_depth_allsites_NINE_REPLICATES.csv
```

# To check the number of SNPs over the entire genome/chrV and be able to compare with the locus: 
```
#!/bin/bash -l

############# SLURM SETTINGS #############
#SBATCH --job-name=AF_allsites_wholegenome       # some descriptive job name of your choice
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

module load grenedalf/0.3.0

############# MY CODE #############
# pwd = /home/jenni/data_folder/working_folder/selection_lines_poolseq/pairwise_comparisons/LINES_IN_CORRECT_ORDER/grenedalf/Jan24_F0_vs_all_sel_lines_only/all_sites/AF

# Calculate the allele counts (ref and main alt) for all sites from V:30-45 Mb

grenedalf frequency \
--sync-path /home/jenni/data_folder/working_folder/selection_lines_poolseq/mpileup/Jan24_correct_order_syncs/SL_ps_only_Hcwbps18_q20Q30_noindel_lines_in_correct_order_Jan24.java.sync \
--reference-genome-fasta-file /home/jenni/data_folder/wbps18/haemonchus_contortus.PRJEB506.WBPS18.genomic.fa \
--rename-samples-file /home/jenni/data_folder/working_folder/selection_lines_poolseq/pairwise_comparisons/LINES_IN_CORRECT_ORDER/grenedalf/Jan24_F0_vs_all_sel_lines_only/all_sites/AF/SL_sample_names.txt \
--write-sample-counts \
--separator-char tab \
--out-dir /home/jenni/data_folder/working_folder/selection_lines_poolseq/pairwise_comparisons/LINES_IN_CORRECT_ORDER/grenedalf/Jan24_F0_vs_all_sel_lines_only/all_sites/AF \
--file-prefix Wholegenome_400_Gf_F0_vs_SL_only_Hcwbps18_q20Q30_noindels_AF_sample_counts \
--log-file grenedalf_AF_all_wholegenome.log
```
```
wc -l Wholegenome_400_Gf_F0_vs_SL_only_Hcwbps18_q20Q30_noindels_AF_sample_countsfrequency.csv
# 101150381 Wholegenome_400_Gf_F0_vs_SL_only_Hcwbps18_q20Q30_noindels_AF_sample_countsfrequency.csv
```
Therefore have 101,150,380 SNPs across the genome

Next, check the proportion of SNPs with ref AF < 0.9. 

First, need to calculate the depth for all sites using the mpileup file:

```
# awk -v assigns a value to a variable before awk starts
for i in {4..30}; 
do
    cat /home/jenni/data_folder/working_folder/selection_lines_poolseq/mpileup/Jan24_correct_order_syncs/SL_ps_only_Hcwbps18_q20Q30_noindel_lines_in_correct_order_Jan24.java.sync | \
    awk -v col=$i 'BEGIN {FS=OFS="\t"} {print $col }' | awk 'BEGIN {FS=OFS=":"} {print $1,$2,$3,$4,$6}' | sed 's/\:/+/g' | bc > sum_col${i} ;
done &

```

Right. Next - use paste to put the columns back together. 
```
cut -f 1-3 /home/jenni/data_folder/working_folder/selection_lines_poolseq/mpileup/Jan24_correct_order_syncs/SL_ps_only_Hcwbps18_q20Q30_noindel_lines_in_correct_order_Jan24.java.sync > positions

paste positions sum_col4 sum_col5 sum_col6 sum_col7 sum_col8 sum_col9 sum_col10 sum_col11 sum_col12 sum_col13 sum_col14 sum_col15 sum_col16 sum_col17 sum_col18 sum_col19 sum_col20 sum_col21 sum_col22 sum_col23 sum_col24 sum_col25 sum_col26 sum_col27 sum_col28 sum_col29 sum_col30 > from_mpileupfile_SL_only_Dec_24_Wholegenome.ACTGdepth
```

Next, remove the second reference allele column from the file - this would mess up the R script later.
```
cut -f 4- from_mpileupfile_SL_only_Dec_24_Wholegenome.ACTGdepth > tmp
cut -f 1-2 from_mpileupfile_SL_only_Dec_24_Wholegenome.ACTGdepth > tmp2

paste tmp2 tmp > tmp3

mv tmp3 from_mpileupfile_SL_only_Dec_24_Wholegenome.ACTGdepth
```

```
#!/bin/bash -l

############# SLURM SETTINGS #############
#SBATCH --job-name=join_R_script       # some descriptive job name of your choice
#SBATCH --output=%x-%J.out      # output file name will contain job name + job ID
#SBATCH --error=%x-%J.err       # error file name will contain job name + job ID
#SBATCH --time=0-3:00:00       # time limit for the whole run, in the form of d-hh:mm:ss, also accepts mm, mm:ss, hh:mm:ss, d-hh, d-hh:mm
#SBATCH --mem=25G                # memory required per node, in the form of [num][M|G|T]
#SBATCH --ntasks=1              # number of Slurm tasks to be launched, increase for multi-process runs ex. MPI
#SBATCH --cpus-per-task=1       # number of processor cores to be assigned for each task, default is 1, increase for multi-threaded runs
#SBATCH --ntasks-per-node=1     # number of tasks to be launched on each allocated node

############# LOADING MODULES (optional) #############
#module load apps/xxx
#module load libs/xxx

module load R/4.3.1

############# MY CODE #############
# pwd = /home/jenni/data_folder/working_folder/selection_lines_poolseq/pairwise_comparisons/LINES_IN_CORRECT_ORDER/grenedalf/Jan24_F0_vs_all_sel_lines_only/all_sites/AF/wholegenome_depth

Rscript /home/jenni/data_folder/working_folder/selection_lines_poolseq/pairwise_comparisons/LINES_IN_CORRECT_ORDER/grenedalf/Jan24_F0_vs_all_sel_lines_only/all_sites/AF/wholegenome_depth/join.R 
```
join.R script below:
```
#!/usr/bin/env Rscript

library(dplyr)
# dplyr_1.1.2

df<-read.table("Wholegenome_400_Gf_F0_vs_SL_only_Hcwbps18_q20Q30_noindels_AF_sample_countsfrequency.csv", header=TRUE, sep="\t")
df2<-read.table("from_mpileupfile_SL_only_Dec_24_Wholegenome.ACTGdepth", header=FALSE, sep="\t")

head(df)
head(df2)

# Join tables, keeping only rows present in both tables (this will account for removal of windows with no SNPs
df3<-inner_join(df, df2, by=c("CHROM"="V1", "POS"="V2"))

head(df3)

write.table(df3, file="Wholegenome_AC_depth_allsites.txt", row.names=FALSE, sep="\t", quote= FALSE) 
```

Note, I then checked the file and changed the header to provide the sample names.
```
# Change the chr names:
sed 's/hcontortus_//g' Wholegenome_AC_depth_allsites.txt | sed 's/_Celeg_TT_arrow_pilon//g' > chrrenamed_Wholegenome_AC_depth_allsites.txt &
```

To obtain only those sites which have ref AF <0.9 merge with the previous file... 
```
awk 'NR == 1 ; $1 == "chr5" && $2 > 30000000 && $2 < 49000000 {print $0}' Wholegenome_AC_depth_allsites_NINE_REPLICATES.txt > ChrV_locus_AC_depth_allsites_NINE_REPLICATES.txt

```
Then, in R: 

```
#!/usr/bin/env Rscript
# Using R v 4.3.1

library(dplyr)
# dplyr_1.1.2

df<-read.table("ChrV_locus_AC_depth_allsites_NINE_REPLICATES.txt", header=TRUE, sep="\t") # New table
df2<-read.table("keep_Chr5_30-45Mb_AC_depth_allsites.csv", header=TRUE, sep="\t") # Previous table filtered to retain only locus sites with ref AF <0.9

head(df)
head(df2)

# Join tables, keeping only rows present in both tables (this will account for removal of windows with no SNPs
df3<-inner_join(df, df2, by=c("POS"="POS"))

head(df3)

write.table(df3, file="ChrV_locus_AC_depth_allsites_NINE_REPLICATES_refAFless0.9.txt", row.names=FALSE, sep="\t", quote= FALSE) 
```


Also, I want to calculate the AF for the reference allele for each of the samples and compare the AF < 0.9 in the locus vs the rest of the genome:
```
module load R/4.4.2
```
calculate_AF_rscript.r below:
```
#!/usr/bin/env Rscript
# Using R v4.4.2

df<-read.table("Wholegenome_AC_depth_allsites.txt", header=T)

df$F3_I_L1_RAF<-(df$F3_I_L1.REF_CNT / df$F3_I_L1.DEPTH)
df$F3_I_L2_RAF<-(df$F3_I_L2.REF_CNT / df$F3_I_L2.DEPTH)
df$F3_I_L3_RAF<-(df$F3_I_L3.REF_CNT / df$F3_I_L3.DEPTH)

df$F3_M_L1_RAF<-(df$F3_M_L1.REF_CNT / df$F3_M_L1.DEPTH)
df$F3_M_L2_RAF<-(df$F3_M_L2.REF_CNT / df$F3_M_L2.DEPTH)
df$F3_M_L3_RAF<-(df$F3_M_L3.REF_CNT / df$F3_M_L3.DEPTH)

summary(df)
# Note, have NA where have a depth of zero for a site

write.table(df, file="F3refAF_Wholegenome_AC_depth_allsites.txt", row.names=FALSE, sep="\t", quote= FALSE)
