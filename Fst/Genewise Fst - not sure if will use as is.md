# To run the genewise Fst:

1. Get the gff3 file from WBPS18.
```
wget https://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS18/species/haemonchus_contortus/PRJEB506/haemonchus_contortus.PRJEB506.WBPS18.annotations.gff3.gz

gunzip haemonchus_contortus.PRJEB506.WBPS18.annotations.gff3.gz
```
Convert it to a GTF file using AGAT (https://github.com/NBISweden/AGAT/tree/master, kindly installed by Oliver at path help on the Sanger).

Note, the first time I ran this it added 133 L2 to link L1 and L3 and then deleted a load of L1 - pseudogenes, and the same ones. So, after raising this on AGAT, NBISweden/AGAT#385 I adjusted the feature file ```agat --help agat levels --expose```
```
agat_convert_sp_gff2gtf.pl --gff haemonchus_contortus.PRJEB506.WBPS18.annotations.gff3 --o haemonchus_contortus.PRJEB506.WBPS18.annotations.gtf
```
Transfer this to scamper:
```
C:
cd "C:/Users/jmi45g/OneDrive - University of Glasgow/BIOINFORMATICS/Roz's Hcon sequencing"
pscp -P 2227 jm62@localhost:/nfs/users/nfs_j/jm62/scratch/genomes/HcWBPS18/haemonchus_contortus.PRJEB506.WBPS18.annotations.gtf ./haemonchus_contortus.PRJEB506.WBPS18.annotations.gtf

scp ./haemonchus_contortus.PRJEB506.WBPS18.annotations.gtf jenni@scamper.mvls.gla.ac.uk:/home/jenni/haemonchus_contortus.PRJEB506.WBPS18.annotations.gtf
```

Using awk, filter the gtf so that only 'gene' in column three is retained.
```
cat /home/jenni/data_folder/working_folder/selection_lines_poolseq/pairwise_comparisons/haemonchus_contortus.PRJEB506.WBPS18.annotations.gtf | awk '$3 == "gene"' > haemonchus_contortus.PRJEB506.WBPS18.genes.gtf &
```
And run Grenedalf v0.3.0 for the genewise Fst:
```
grenedalf fst --sync-path /home/jenni/data_folder/working_folder/selection_lines_poolseq/mpileup/Jan24_correct_order_syncs/SLps_only_co_Hcwbps18_q20Q30_noindels_ss57.java.sync --reference-genome-fasta-file /home/jenni/data_folder/wbps18/haemonchus_contortus.PRJEB506.WBPS18.genomic.fa --rename-samples-file /home/jenni/data_folder/working_folder/selection_lines_poolseq/pairwise_comparisons/LINES_IN_CORRECT_ORDER/grenedalf/Jan24_F0_vs_all_sel_lines_only/SL_ss57_sample_names.txt --filter-sample-min-count 2 --filter-sample-min-coverage 57 --filter-sample-max-coverage 57 --window-type regions --window-region-gff /home/jenni/data_folder/wbps18/haemonchus_contortus.PRJEB506.WBPS18.genes.gtf --method unbiased-nei --pool-sizes 400 --omit-na-windows --comparand-list /home/jenni/data_folder/working_folder/selection_lines_poolseq/pairwise_comparisons/LINES_IN_CORRECT_ORDER/grenedalf/Jan24_F0_vs_all_sel_lines_only/SL_pairwise_comparisons.txt --separator-char tab --out-dir /home/jenni/data_folder/working_folder/selection_lines_poolseq/pairwise_comparisons/LINES_IN_CORRECT_ORDER/grenedalf/Jan24_F0_vs_all_sel_lines_only/ss57/genewise_fst --file-prefix 400_Gf_F0_vs_SL_only_Sel_lines_Hcwbps18_q20Q30_noindels_ss57_genewise_fst --log-file grenedalf_gtf_gene_fst_ss57.log 2>error_log_gtf_gene_fst &
```
Ok, now what to do? 

Q: Are the genes with variants the SAME in all lines? Are they the same between IVM and MOX? 
 How to answer? - Could sort, rank and then get the average of all the ranking per gene and then select the top 1%? And plot this with IVM on x and MOX on y and see if a correlation? Could then label genes that are within ChrV, or just, say the top 10, or those of interest (e.g. cky-1)

Q: Which genes are the most important? How to handle low diversity with high Fst vs higher diversity with lower Fst? 
 Could I rank them? And take the top 1%? Could I plot and colour by snps/gene? Could I look at, say, the top 1% in snps/gene < 10 and the top 1% of snps/gene >10? Would 10 really be the ideal cut-off? 

Q: Which genes also have a reduction in ref AF? Which genes also have a variant of potential interest by snpeff?
 To need to intersect with the AF data and the snpeff data


Q: Could I do a heatmap of the genes and Fst too? 


```
cd /home/jenni/data_folder/working_folder/selection_lines_poolseq/pairwise_comparisons/LINES_IN_CORRECT_ORDER/grenedalf/Jan24_F0_vs_all_sel_lines_only/ss57/genewise_fst

wc -l 400*
# 17890 -1 = genes in the file. Therefore the top 1% would be equivalent to the 178 top genes.
```
Checked for a correlation between the Fst and the number of SNPs in a gene, excluding NA results. Doesn't look to be particularly strong. So that is positive.
```
> cor(df$F0.F1_C_L1, df$snps, use="complete.obs")
[1] -0.01234874
> cor(df$F0.F1_M_L1, df$snps, use="complete.obs")
[1] 0.001875482
> cor(df$F0.F3_M_L1, df$snps, use="complete.obs")
[1] 0.01396289
> cor(df$F0.F3_M_L2, df$snps, use="complete.obs")
[1] -0.005571892
> cor(df$F0.F3_M_L3, df$snps, use="complete.obs")
[1] -0.006968888
> cor(df$F0.F3_I_L3, df$snps, use="complete.obs")
[1] 0.001156476
> cor(df$F0.F3_I_L2, df$snps, use="complete.obs")
[1] 0.007487804
> cor(df$F0.F3_I_L1, df$snps, use="complete.obs")
[1] -0.003762768
> cor(df$F3_I_L2.F3_I_L3, df$snps, use="complete.obs")
[1] -0.00513386
> cor(df$F3_M_L2.F3_M_L3, df$snps, use="complete.obs")
[1] 0.02409203
> cor(df$F3_M_L1.F3_M_L3, df$snps, use="complete.obs")
[1] 0.02869724
```
```
# To rank each gene by the sample do: 

# Get the header:
head -n 1 400_Gf_F0_vs_SL_only_Sel_lines_Hcwbps18_q20Q30_noindels_ss57_genewise_fstfst.csv > header
# Get the rest of the file:
tail -n +2 400_Gf_F0_vs_SL_only_Sel_lines_Hcwbps18_q20Q30_noindels_ss57_genewise_fstfst.csv > tmp
# Sort using -g to allow for scientific notation, nan and negative numbers. Use -r to reverse, so that the largest is at the top (nan at the bottom).
# nl will add line numbers at the start of the file.
# Increase by 2 each time. This will progress through the IVM F3 samples, then the MOX F3 samples cf F0.
sort -g -k25 -r tmp | nl  > tmp1
sort -g -k27 -r tmp1 | nl > tmp2
sort -g -k29 -r tmp2 | nl > tmp3
sort -g -k31 -r tmp3 | nl > tmp4
sort -g -k33 -r tmp4 | nl > tmp5
sort -g -k35 -r tmp5 | nl > tmp6

# Checked was sane. Good. 
# Adjust the header with nano and cat the two together:

cat header tmp6 > tmp7
```
Then get the average using awk:
```
awk 'BEGIN {OFS="\t"} {print (($1 + $2 + $3)/3)}' tmp7 > avg_rank_mox
# adjust header (from 0 to avg_rank_mox)
awk 'BEGIN {OFS="\t"} {print (($4 + $5 + $6)/3)}' tmp7 > avg_rank_ivm
# adjust header (from 0 to avg_rank_ivm)
paste avg_rank_mox avg_rank_ivm tmp7 > tmp8
```
In total I have 17889 genes (wc -l tmp8, takeaway 1 for the header line)
```
module load R/4.3.1

df<-read.table("tmp8", header=T)
summary(df)

x<-ggplot(data=df, aes(x=-log10(avg_rank_ivm/17889), y=-log10(avg_rank_mox/17889)))+
geom_point(aes(colour=chrom))+
geom_hline(yintercept = (-log10(178/17889)), linetype = 2)+
geom_vline(xintercept = (-log10(178/17889)), linetype = 2)+
geom_abline(intercept = 0, slope = 1)+
geom_smooth()

pdf("HcPS_SLonly_genewise_avgrankivmmox_Feb2024.pdf",width = 16,height = 12,useDingbats = FALSE)
x
dev.off()

> cor(df$avg_rank_mox, df$avg_rank_ivm)
[1] 0.6777554
```
```
scp jenni@scamper.mvls.gla.ac.uk:/home/jenni/data_folder/working_folder/selection_lines_poolseq/pairwise_comparisons/LINES_IN_CORRECT_ORDER/grenedalf/Jan24_F0_vs_all_sel_lines_only/ss57/genewise_fst/HcPS_SLonly_genewise_avgrankivmmox_Feb2024.pdf ./HcPS_SLonly_genewise_avgrankivmmox_Feb2024.pdf 
```
_Average rank for ivm and mox for genes, with geom_smooth plotted. hline and vline showing the 178/178889 ranking cut-off - 178 would be 1% of genes_
![image](https://github.com/SheepwormJM/Hitting-A-Moving-Target-PoolSeq-Generations/assets/55552826/b5efe00a-1f0a-40a6-8809-fbdc41e39887)

Try plotting the genewise Fst as a heatmap - want the data in long format (so similar to how I usually have it) - need an X a Y and a Z (the value to colour) - e.g. sample, genename, Fst value

Try this blog: https://www.royfrancis.com/a-guide-to-elegant-tiled-heatmaps-in-r-2019/ 
```
#basic ggplot
p <- ggplot(m3, aes(x=year, y=state, fill=count))+
      geom_tile()

#modified ggplot
p <- ggplot(m3, aes(x=year, y=state, fill=count))+
  #add border white colour of line thickness 0.25
  geom_tile(colour="white", size=0.25)+
  #remove x and y axis labels
  labs(x="", y="")+
  #remove extra space
  scale_y_discrete(expand=c(0, 0))+
  #define new breaks on x-axis
  scale_x_discrete(expand=c(0, 0),
                    breaks=c("1930", "1940", "1950", "1960", "1970", "1980", "1990", "2000"))+
  #set a base size for all fonts
  theme_grey(base_size=8)+
  #theme options
  theme(
    #bold font for legend text
    legend.text=element_text(face="bold"),
    #set thickness of axis ticks
    axis.ticks=element_line(size=0.4),
    #remove plot background
    plot.background=element_blank(),
    #remove plot border
    panel.border=element_blank()
  )

Want to get the genenames too. 

```
```
# Get the chromosome, start, end and gene name from the genes.gtf file made earlier
# Note that I have had to initially specifcy that the input FS was the same as the OFS to get the right colums
# But then I deliberately let anything be a field separator (it would then accept a space) to select just the gene name itself
awk 'BEGIN {FS=OFS="\t"} {print $1, $4, $5, $9}' haemonchus_contortus.PRJEB506.WBPS18.genes.gtf | awk 'BEGIN {OFS="\t"} {print $1, $2, $3, $5}' | sed 's/"gene://g' | sed 's/";//g' > tmp2

mv tmp2 genenames.gtf.bed
```
```
awk 'BEGIN {OFS="\t"} {print $9, $10, $11}' tmp8 > tmp8_coords

paste tmp8_coords tmp8 > tmp9.bed
```

```
module load bedtools/2.31

bedtools intersect -a /home/jenni/data_folder/wbps18/genenames.gtf.bed -b tmp9.bed -wa -wb -f 2E-9 > tmp10 | cut -f 4- > tmp10
```

Gene numbers don't make sense:
```
19778 /home/jenni/data_folder/wbps18/genenames.gtf.bed
17890 tmp9.bed (includes a header line)
22040 tmp10
```
Ok, so Grenedalf has omitted some genes. Perhaps they had no variant sites. 

Bedtools has repeated some genes: 
```
hcontortus_chr5_Celeg_TT_arrow_pilon    48822312        48835932        HCON_00163820   hcontortus_chr5_Celeg_TT_arrow_pilon    48835871        48857893           1563.33 1813.67   hcontortus_chr5_Celeg_TT_arrow_pilon    48822312        48835932        HCON_00163820   hcontortus_chr5_Celeg_TT_arrow_pilon    48822312        48835932           1322.67 1613

hcontortus_chr5_Celeg_TT_arrow_pilon    48835871        48857893        HCON_00163840   hcontortus_chr5_Celeg_TT_arrow_pilon    48835871        48857893           1563.33 1813.67   hcontortus_chr5_Celeg_TT_arrow_pilon    48835871        48857893        HCON_00163840   hcontortus_chr5_Celeg_TT_arrow_pilon    48822312        48835932           1322.67 1613
```

_The UTRs overlap_
![image](https://github.com/SheepwormJM/Hitting-A-Moving-Target-PoolSeq-Generations/assets/55552826/038847cd-2d94-42da-a145-571d52fbf476)

Try dpplyr package instead! https://github.com/rstudio/cheatsheets/blob/main/data-transformation.pdf

_Note, got an error when loading file - had an nbis gene in it. Checked - some weird non-translated gene. Decided to try leaving in, but removed the " before it with sed._

```
sed 's/"//g' /home/jenni/data_folder/wbps18/genenames.gtf.bed > /home/jenni/data_folder/wbps18/newgenenames.gtf.bed
```

```
module load R/4.3.1
R

library(dplyr)

df<-read.table("/home/jenni/data_folder/wbps18/newgenenames.gtf.bed", header=F)
df2<-read.table("tmp9.bed", header=T)
head(df)
names(df)
names(df2)

# Join the two tables, using the chrom, start and end positions.
# Note that for col1 = col2 it will assume col1 is in the first df and col2 is in the second
# By using all three it will ensure that only if they all match will it join, so it will not duplicate
# It will output the remaining columns in each df, but only the matching columns once

new<-inner_join(df,  df2, by=c("V1"="chrom","V2"="start","V3"="end"))

# Get just the top 10% Fst ranked genes in either MOX or IVM F3 generation (supposedly, but using the average, so will in fact be less)
tenpc<-subset(new,new$avg_rank_mox < 1790 | new$avg_rank_ivm < 1790)
# Have just 1361 genes remaining (by using summary, had 17889 for new)

df3<-tenpc[,c(4,17:51)]
names(df3)

# melt the table in order to plot a heatmap
library(reshape2)
meltdf3<-(melt(df3, id.vars=c("V4")))

# plot a heatmap of the Fst for each gene by each sample comparison
library(ggplot2)
#basic ggplot
p <- ggplot(meltdf3, aes(x=variable, y=V4, fill=value))+
      geom_tile()


pdf("HcPS_SLonly_genewise_heatmap_top10pcFst_Feb2024.pdf",width = 16,height = 12,useDingbats = FALSE)
p
dev.off()

# Get just the top 1% Fst ranked genes in either MOX or IVM F3 generation (supposedly, but using the average, so will in fact be less)
onepc<-subset(new,new$avg_rank_mox < 179 | new$avg_rank_ivm < 179)
# Have just 159 genes remaining (by using summary, had 17889 for new)

df4<-onepc[,c(4,17:51)]
names(df4)

# melt the table in order to plot a heatmap
library(reshape2)
meltdf4<-(melt(df4, id.vars=c("V4")))

# plot a heatmap of the Fst for each gene by each sample comparison
library(ggplot2)
#basic ggplot
p <- ggplot(meltdf4, aes(x=variable, y=V4, fill=value))+
      geom_tile()


pdf("HcPS_SLonly_genewise_heatmap_onepcFst_Feb2024.pdf",width = 16,height = 12,useDingbats = FALSE)
p
dev.off()


# Get just the top 0.1% Fst ranked genes in either MOX or IVM F3 generation (supposedly, but using the average, so will in fact be less)
pointonepc<-subset(new,new$avg_rank_mox < 18 | new$avg_rank_ivm < 18)
# Have just 9 genes remaining (by using summary, had 17889 for new)

# Note that there is some variation in the top 9 genes!
      V1                  V2                 V3                V4           
 Length:9           Min.   :37239956   Min.   :37240258   Length:9          
 Class :character   1st Qu.:37934935   1st Qu.:37941953   Class :character  
 Mode  :character   Median :41124535   Median :41142589   Mode  :character  
                    Mean   :39968904   Mean   :39978298                     
                    3rd Qu.:41274614   3rd Qu.:41283249                     
                    Max.   :41309214   Max.   :41318844                     
  avg_rank_mox     avg_rank_ivm    rank_F0.F3_M_L3  rank_F0.F3_M_L2 
 Min.   : 2.667   Min.   : 1.333   Min.   : 2.000   Min.   : 1.000  
 1st Qu.: 4.667   1st Qu.: 8.000   1st Qu.: 4.000   1st Qu.: 3.000  
 Median : 7.333   Median : 9.333   Median : 8.000   Median : 5.000  
 Mean   :13.407   Mean   :20.630   Mean   : 9.444   Mean   : 6.778  
 3rd Qu.:13.667   3rd Qu.:29.667   3rd Qu.:13.000   3rd Qu.: 9.000  
 Max.   :57.333   Max.   :59.667   Max.   :21.000   Max.   :15.000

df5<-pointonepc[,c(4,17:51)]
names(df5)

# melt the table in order to plot a heatmap
library(reshape2)
meltdf5<-(melt(df5, id.vars=c("V4")))

# plot a heatmap of the Fst for each gene by each sample comparison
library(ggplot2)
#basic ggplot
p <- ggplot(meltdf5, aes(x=variable, y=V4, fill=value))+
      geom_tile()


pdf("HcPS_SLonly_genewise_heatmap_pointonepcFst_Feb2024.pdf",width = 16,height = 12,useDingbats = FALSE)
p
dev.off()


write.table(pointonepc, file="HcPS_SLonly_genewise_Fst_pointonepc_ivm_or_moc.txt", row.names=FALSE, sep="\t", quote= FALSE)
write.table(onepc, file="HcPS_SLonly_genewise_Fst_onepc_ivm_or_moc.txt", row.names=FALSE, sep="\t", quote= FALSE)
write.table(tenpc, file="HcPS_SLonly_genewise_Fst_tenpc_ivm_or_moc.txt", row.names=FALSE, sep="\t", quote= FALSE)
write.table(new, file="HcPS_SLonly_genewise_Fst_with_gene_names.txt", row.names=FALSE, sep="\t", quote= FALSE)
```
```
scp jenni@scamper.mvls.gla.ac.uk:/home/jenni/data_folder/working_folder/selection_lines_poolseq/pairwise_comparisons/LINES_IN_CORRECT_ORDER/grenedalf/Jan24_F0_vs_all_sel_lines_only/ss57/genewise_fst/HcPS_SLonly_genewise_heatmap_top10pcFst_Feb2024.pdf ./HcPS_SLonly_genewise_heatmap_top10pcFst_Feb2024.pdf

scp jenni@scamper.mvls.gla.ac.uk:/home/jenni/data_folder/working_folder/selection_lines_poolseq/pairwise_comparisons/LINES_IN_CORRECT_ORDER/grenedalf/Jan24_F0_vs_all_sel_lines_only/ss57/genewise_fst/HcPS_SLonly_genewise_heatmap_onepcFst_Feb2024.pdf ./HcPS_SLonly_genewise_heatmap_onepcFst_Feb2024.pdf

scp -r jenni@scamper.mvls.gla.ac.uk:/home/jenni/data_folder/working_folder/selection_lines_poolseq/pairwise_comparisons/LINES_IN_CORRECT_ORDER/grenedalf/Jan24_F0_vs_all_sel_lines_only/ss57/genewise_fst/HcPS_SLonly_genewise_FST_outputs ./HcPS_SLonly_genewise_FST_outputs
```
_Heatmap of 1361 top ranked (avg mox or ivm) genes, Fst value for colour_
![image](https://github.com/SheepwormJM/Hitting-A-Moving-Target-PoolSeq-Generations/assets/55552826/c0f9da4f-4347-4062-9da4-b2593455a57c)


![image](https://github.com/SheepwormJM/Hitting-A-Moving-Target-PoolSeq-Generations/assets/55552826/7367452d-29c6-469d-8ef3-a7c31867c214)

![image](https://github.com/SheepwormJM/Hitting-A-Moving-Target-PoolSeq-Generations/assets/55552826/758bdd67-9ae5-4701-a566-63ed11925f3d)
