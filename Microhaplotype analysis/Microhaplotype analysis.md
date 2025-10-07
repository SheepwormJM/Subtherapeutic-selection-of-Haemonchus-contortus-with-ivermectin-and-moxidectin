To generate microhaplotypes from read data, rather than using individual SNPs to investigate:
  1. Whether this helps in resolving the locus
  2. Whether multiple common haplotypes are under selection by IVM or MOX, or only one

Trying, first, the **mhFromLowDepSeq** by Thomas Delomas and Stuart Willis:
https://github.com/delomast/mhFromLowDepSeq 

Requires as input:
- Bam files of aligned reads and their corresponding population.... The bam files must be sorted by position and indexed.
- A file of 'substitution SNPs' - i.e. not indels or complex variants, just basic SNPs

Note that the python script will not exclude duplicate reads. I still have these in my bam files (marked but not removed) - therefore, will first need to remove these. I will keep in secondary/split hits as these could be true alignments and if removed could remove particular microhaplotypes from the data which are associated with TEs. 

I am also going to remove MAPQ<20, in line with the mpileup I had made before.

```
# This was what was run before - to remove unmapped reads - either the read itself or with it's mate unmapped
# samtools view --threads 4 -F 12 -b tmp.merged.sorted.marked.bam -o $"{sample_name}".merged.sorted.marked.bam
```
Get the F3 population bam files ready to loop:
```
cat Populations_F3only_bams.txt | awk '{print $1}' > BAMS.txt

sed 's/.*\///g' BAMS.txt | sed 's/\.bam//g' > BAMS_out.txt  # Remove everything up to and including the last / in each line in the file, then remove the .bam.

paste -d "\t" BAMS.txt BAMS_out.txt | sed 's/^/samtools view --threads 4 -F 1036 -q 20 -b /g' | sed 's/\t/ -o /g' | sed 's/$/\.dupsremoved\.bam/g' > SAMTOOLS_BAMS.txt
```
```
(This did not work.... tried loads of other ideas. Wanted to loop and change the output file. Should be removing the front of the file, which it does but only the very first \:

do
name=${NAME##/users}
echo $NAME
echo $name ;
done < BAMS.txt

name=$"{sed 's/$.*$\///g' $NAME | sed 's/.bam//g'}"\




while read NAME ;
do
awk '{print $1}' > input_file
input_file={awk '{print $1}'}
output_file=${NAME}[2]
echo $input_file
echo $output_file ;
done < NEW_BAMS.txt
```
```
#!/bin/bash -l

############# SLURM SETTINGS #############
#SBATCH --job-name=filter_bams       # some descriptive job name of your choice
#SBATCH --output=%x-%J.out      # output file name will contain job name + job ID
#SBATCH --error=%x-%J.err       # error file name will contain job name + job ID
#SBATCH --time=0-12:00:00       # time limit for the whole run, in the form of d-hh:mm:ss, also accepts mm, mm:ss, hh:mm:ss, d-hh, d-hh:mm
#SBATCH --mem=20G                # memory required per node, in the form of [num][M|G|T]
#SBATCH --ntasks=1              # number of Slurm tasks to be launched, increase for multi-process runs ex. MPI
#SBATCH --cpus-per-task=1       # number of processor cores to be assigned for each task, default is 1, increase for multi-threaded runs
#SBATCH --ntasks-per-node=1     # number of tasks to be launched on each allocated node

############# LOADING MODULES (optional) #############
#module load apps/xxx
#module load libs/xxx

module load apps/miniforge/24.7.1/bin
conda activate bwa_splitter

############# MY CODE #############

samtools view --threads 4 -F 1036 -q 20 -b /users/jmi45g/project0005/HCON/DATA/SLPS_Illumina/Aligned_bam_files/bam_files/3-14949_postCTL3.merged.sorted.marked.bam -o 3-14949_postCTL3.merged.sorted.marked.dupsremoved.bam
samtools view --threads 4 -F 1036 -q 20 -b /users/jmi45g/project0005/HCON/DATA/SLPS_Illumina/Aligned_bam_files/bam_files/1-14031_CTL3.merged.sorted.marked.bam -o 1-14031_CTL3.merged.sorted.marked.dupsremoved.bam
samtools view --threads 4 -F 1036 -q 20 -b /users/jmi45g/project0005/HCON/DATA/SLPS_Illumina/Aligned_bam_files/bam_files/2-14036_postCTL3.merged.sorted.marked.bam -o 2-14036_postCTL3.merged.sorted.marked.dupsremoved.bam
samtools view --threads 4 -F 1036 -q 20 -b /users/jmi45g/project0005/HCON/DATA/SLPS_Illumina/Aligned_bam_files/bam_files/6-14917_postIVM3.merged.sorted.marked.bam -o 6-14917_postIVM3.merged.sorted.marked.dupsremoved.bam
samtools view --threads 4 -F 1036 -q 20 -b /users/jmi45g/project0005/HCON/DATA/SLPS_Illumina/Aligned_bam_files/bam_files/5-14916_postIVM3.merged.sorted.marked.bam -o 5-14916_postIVM3.merged.sorted.marked.dupsremoved.bam
samtools view --threads 4 -F 1036 -q 20 -b /users/jmi45g/project0005/HCON/DATA/SLPS_Illumina/Aligned_bam_files/bam_files/4-14907_postIVM3.merged.sorted.marked.bam -o 4-14907_postIVM3.merged.sorted.marked.dupsremoved.bam
samtools view --threads 4 -F 1036 -q 20 -b /users/jmi45g/project0005/HCON/DATA/SLPS_Illumina/Aligned_bam_files/bam_files/9-14954_postMOX3.merged.sorted.marked.bam -o 9-14954_postMOX3.merged.sorted.marked.dupsremoved.bam
samtools view --threads 4 -F 1036 -q 20 -b /users/jmi45g/project0005/HCON/DATA/SLPS_Illumina/Aligned_bam_files/bam_files/7-14029_postMOX3.merged.sorted.marked.bam -o 7-14029_postMOX3.merged.sorted.marked.dupsremoved.bam
samtools view --threads 4 -F 1036 -q 20 -b /users/jmi45g/project0005/HCON/DATA/SLPS_Illumina/Aligned_bam_files/bam_files/8-14953_postMOX3.merged.sorted.marked.bam -o 8-14953_postMOX3.merged.sorted.marked.dupsremoved.bam
```

```
#!/bin/bash -l

############# SLURM SETTINGS #############
#SBATCH --job-name=index_bams       # some descriptive job name of your choice
#SBATCH --output=%x-%J.out      # output file name will contain job name + job ID
#SBATCH --error=%x-%J.err       # error file name will contain job name + job ID
#SBATCH --time=0-5:00:00       # time limit for the whole run, in the form of d-hh:mm:ss, also accepts mm, mm:ss, hh:mm:ss, d-hh, d-hh:mm
#SBATCH --mem=20G                # memory required per node, in the form of [num][M|G|T]
#SBATCH --ntasks=1              # number of Slurm tasks to be launched, increase for multi-process runs ex. MPI
#SBATCH --cpus-per-task=1       # number of processor cores to be assigned for each task, default is 1, increase for multi-threaded runs
#SBATCH --ntasks-per-node=1     # number of tasks to be launched on each allocated node

############# LOADING MODULES (optional) #############
#module load apps/xxx
#module load libs/xxx

module load apps/miniforge/24.7.1/bin
conda activate bwa_splitter

############# MY CODE #############

for i in *dupsremoved.bam ;
do
samtools index -b ${i} ;
done

```
```
# Python code available at the Github - calc_mh_freq.py
# The python script utilized the packages Pysam (https://www.github.com/pysam-developers/pysam), numpy, and numba

python3 calc_mh_freq.py -m ./example/filePop.txt -s ./example/snpLoc.txt -o outputPool.txt -af -w 65 -pool
```

Hmm... thoughts: 
1. Max window size - default is 125 bp. Obviously wouldn't want it to be longer than the read length, unless reads overlap but equally would I want it to be the full read length if the quality reduced over the read? How far along do I want it?
2. What will I do with the microhaplotype information? I will have heterozygosity, can get allele frequencies, and read depth for the locus for each populations.
2.a.  I could identify low coverage areas (<20) and filter these/treat these separately, as get increased bias with coverage lower than 20.
2.b.  I could plot the heterozygosity over the chromosome, and plot the allele frequency over the chromosome too...
   2.c. I could look for areas which have very different heterozygosity between control/parental and selected? Or where have change in either. But unsure whether looking for a reduction or increase in allele frequency/heterozygosity!
   2.d. compare the heterozygosity/allele frequency between IVM and MOX re inheritance - could expect IVM to be higher than MOX...
   2.e. Not sure how I can actually 'see' the microhaplotype alleles and their overlap.... I think that this will basically output the co-ordinates of SNPs, but not the actual sequences, or their individual frequencies. So can't tell what mh-alleles are being selected.
4. How many SNPs do I want max in a locus? They have a default of 25/125 bp. Do I have more or less in my nematodes?? I could look at all, with a cut-off of SNPs the same as the window length, say, and see what my proportions are - would need to calculate the number of snps of course. Or could look at the density of snps as a starting value and go with 3 or 5 SD either side of this? Or look at areas of high depth (collapsed) and calculate the number of snps there and then cut off below this?? Actually - the number of SNPs allowed in a locus is NOT the same as the number of alleles (set much much much higher). The number of SNPs in a locus could be more realistically related to the nucleotide diversity. So, for us.... this is pretty darn high.** Go with the upper 90% value for the pi? **

```
I have 150 bp PE reads, so a locus of 125 bp (default) would be very reasonable.
```
```
# Pi results from 100kb Grenedalf window analysis:

chrom        mid        F0.theta_pi_rel F1_C_L1.theta_pi_rel    F1_C_L2.theta_pi_rel    F1_C_L3.theta_pi_rel    F2_C_L1.theta_pi_rel    F2_C_L2.theta_pi_rel    F2_C_L3.theta_pi_rel        F3_C_L1.theta_pi_rel    F3_C_L2.theta_pi_rel    F3_C_L3.theta_pi_rel    F1_I_L2.theta_pi_rel    F1_I_L3.theta_pi_rel    F2_I_L1.theta_pi_rel        F2_I_L2.theta_pi_rel    F2_I_L3.theta_pi_rel    F3_I_L1.theta_pi_rel    F3_I_L2.theta_pi_rel    F3_I_L3.theta_pi_rel    F1_M_L1.theta_pi_rel    F1_M_L2.theta_pi_rel        F1_M_L3.theta_pi_rel    F2_M_L1.theta_pi_rel    F2_M_L2.theta_pi_rel    F2_M_L3.theta_pi_rel    F3_M_L1.theta_pi_rel    F3_M_L2.theta_pi_rel    F3_M_L3.theta_pi_rel
1:458   Min.   :   50000        Min.   :0.00000         Min.   :0.00000         Min.   :0.00000         Min.   :0.00000         Min.   :0.00000         Min.   :0.00000     Min.   :0.00000         Min.   :0.00000         Min.   :0.00000         Min.   :0.00000         Min.   :0.00000         Min.   :0.00000         Min.   :0.00000     Min.   :0.00000         Min.   :0.00000         Min.   :0.00000         Min.   :0.00000         Min.   :0.00000         Min.   :0.00000         Min.   :0.00000     Min.   :0.00000         Min.   :0.00000         Min.   :0.00000         Min.   :0.00000         Min.   :0.00000         Min.   :0.00000         Min.   :0.00000  
2:474   1st Qu.:11850000        1st Qu.:0.01252         1st Qu.:0.01241         1st Qu.:0.01236         1st Qu.:0.01254         1st Qu.:0.01228         1st Qu.:0.01245     1st Qu.:0.01252         1st Qu.:0.01240         1st Qu.:0.01264         1st Qu.:0.01247         1st Qu.:0.01258         1st Qu.:0.01141         1st Qu.:0.01068     1st Qu.:0.01217         1st Qu.:0.01202         1st Qu.:0.01182         1st Qu.:0.01074         1st Qu.:0.01217         1st Qu.:0.01551         1st Qu.:0.01231     1st Qu.:0.01238         1st Qu.:0.01495         1st Qu.:0.01165         1st Qu.:0.01159         1st Qu.:0.01563         1st Qu.:0.01129         1st Qu.:0.01056  
3:436   Median :23650000        Median :0.02182         Median :0.02181         Median :0.02172         Median :0.02183         Median :0.02180         Median :0.02184     Median :0.02193         Median :0.02192         Median :0.02191         Median :0.02204         Median :0.02219         Median :0.02125         Median :0.02001     Median :0.02115         Median :0.02143         Median :0.02140         Median :0.01956         Median :0.02105         Median :0.02533         Median :0.02180     Median :0.02158         Median :0.02462         Median :0.02120         Median :0.02063         Median :0.02536         Median :0.02035         Median :0.02054  
4:519   Mean   :23714117        Mean   :0.02087         Mean   :0.02090         Mean   :0.02078         Mean   :0.02090         Mean   :0.02085         Mean   :0.02093     Mean   :0.02101         Mean   :0.02096         Mean   :0.02104         Mean   :0.02100         Mean   :0.02119         Mean   :0.02023         Mean   :0.01961     Mean   :0.02058         Mean   :0.02067         Mean   :0.02054         Mean   :0.01915         Mean   :0.02045         Mean   :0.02396         Mean   :0.02107     Mean   :0.02064         Mean   :0.02333         Mean   :0.02049         Mean   :0.02035         Mean   :0.02374         Mean   :0.02011         Mean   :0.02010  
5:489   3rd Qu.:35450000        3rd Qu.:0.02948         3rd Qu.:0.02966         3rd Qu.:0.02931         3rd Qu.:0.02949         3rd Qu.:0.02942         3rd Qu.:0.02961     3rd Qu.:0.02956         3rd Qu.:0.02947         3rd Qu.:0.02965         3rd Qu.:0.02967         3rd Qu.:0.02999         3rd Qu.:0.02886         3rd Qu.:0.02775     3rd Qu.:0.02885         3rd Qu.:0.02899         3rd Qu.:0.02874         3rd Qu.:0.02700         3rd Qu.:0.02855         3rd Qu.:0.03314         3rd Qu.:0.03012     3rd Qu.:0.02897         3rd Qu.:0.03233         3rd Qu.:0.02884         3rd Qu.:0.02895         3rd Qu.:0.03266         3rd Qu.:0.02849         3rd Qu.:0.02882  
X:461   Max.   :51850000        Max.   :0.08707         Max.   :0.07819         Max.   :0.09276         Max.   :0.06018         Max.   :0.06579         Max.   :0.08691     Max.   :0.06183         Max.   :0.07018         Max.   :0.06904         Max.   :0.06644         Max.   :0.05870         Max.   :0.07895         Max.   :0.07426     Max.   :0.07482         Max.   :0.06336         Max.   :0.06141         Max.   :0.06220         Max.   :0.06354         Max.   :0.06489         Max.   :0.08159     Max.   :0.07034         Max.   :0.06256         Max.   :0.07018         Max.   :0.08503         Max.   :0.06292         Max.   :0.07018         Max.   :0.05911 


# So, Pi would suggest that for every 100 bp of read sequenced these Hcon populations have a mean of 2 bp as SNPs, but 75% reads could have 3 bp as SNPs, and some have up to almost 9bp in some populations.

# However, this is per read. The SNPs in the locus in the paper, and SNPs overall are much higher - maybe 1 SNP every 3 bp:
"Despite the use of an inbred parental isolate, and seemingly severe bottlenecking of the population at the first round of selection, over 100 million SNPs were retained across the genome in selected lines, averaging 1 SNP every 2.8 bp.  Even with further filtering based on depth and a minimum allele count, between 15.5 to 18 million SNPs were retained genome-wide depending on the samples used for the Cochran-Mantel-Haenzel test (Table 5). Without any depth filtering, 4.6 million SNPs were identified within the 30-45 Mb region of Chromosome V, of which 938,315 SNPs had an alternate allele frequency of at least 0.1."
```
```
# So, to ask for a max number of SNPs of about 44 SNPs in a 125 bp locus as max would fit with genome wide number of SNPs.
# However, for common SNPs (Alt AF >= 0.1) - it would be less than 1 SNP per 125 bp locus, which seems too low cf the Pi (the latter calculated with ss57 data)

# Therefore, to try:
- max 25 SNPs 
- max 40 SNPs (based on 90% of max pi)
```

5. Minimum allele frequency:
```
-minAF FLOAT minimum allele frequency to keep an allele in the model when pruning. Note that pruning is not performed after the last cycle, so the final result may have alleles with frequencies below this minimum. default 0.001
```
Hum.... Could try the default, but also try upping it to, say, 0.01 given the 'common' allele frequency was 0.1 in the paper and we are thinking along the lines of being selected. However, also this could be influenced by the number - we have 200 individuals in a pool, 400 chromosomes. And a variable coverage, but never high enough to sample all chromosomes, and indeed, may commonly sample 12.5% (50x). Therefore, an AF of 0.001 is very low for this coverage, surely... ? If increased so that have at least 2 reads.... but then, it could be real and rare and only be there once. If only there once, the AF would be, for 50X coverage... 0.02. Maybe raise it to this, or to 0.04 - two reads required for 50 X. 
```
# Therefore, to try:
- minAF of 0.02
- minAf of 0.04
```
Actually - the minAF is not really dependent on the read depth as such, I don't think, although it must have some influence upon this. And it can be very very small - see the Appendix 1 in the paper. So, for want of better understanding of their model, will stick with their default value of 0.001. 

First things first!
- The Bam files are on MARS.
- MARS has python3, but it doesn't have all the required packages.
```
python3 -m pip install pysam
python3 -m pip install numba
```

The bam files list (with paths to bams) is here: ```/users/jmi45g/project0005/HCON/DATA/SLPS_Illumina/Aligned_bam_files/bam_files/SL_only_Dec24_bam.list``` - checked correct order and all good. Used ```paste``` and a sample names file to put together to make the ```Populations_bams.txt```

Where to get the SNP position file??

I could get it from Grenedalf, calculating allele frequencies... 
```
sed 's/Chr5_locus_sellinesonly_Hcwbps18_q20Q30_noindel_lines_in_correct_order_Jan24/SL_ps_only_Hcwbps18_q20Q30_noindel_lines_in_correct_order_Jan24/g' /home/jenni/data_folder/working_folder/selection_lines_poolseq/pairwise_comparisons/LINES_IN_CORRECT_ORDER/grenedalf/Jan24_F0_vs_all_sel_lines_only/all_sites/AF/SL_sample_names.txt > /home/jenni/data_folder/working_folder/selection_lines_poolseq/pairwise_comparisons/LINES_IN_CORRECT_ORDER/grenedalf/Jan24_F0_vs_all_sel_lines_only/all_sites/AF/SL_sample_names_for_full_genome.txt
```

```
#!/bin/bash -l

############# SLURM SETTINGS #############
#SBATCH --job-name=AF_allsites_entiregenome       # some descriptive job name of your choice
#SBATCH --output=%x-%J.out      # output file name will contain job name + job ID
#SBATCH --error=%x-%J.err       # error file name will contain job name + job ID
#SBATCH --time=0-24:00:00       # time limit for the whole run, in the form of d-hh:mm:ss, also accepts mm, mm:ss, hh:mm:ss, d-hh, d-hh:mm
#SBATCH --mem=10G                # memory required per node, in the form of [num][M|G|T]
#SBATCH --ntasks=1              # number of Slurm tasks to be launched, increase for multi-process runs ex. MPI
#SBATCH --cpus-per-task=1       # number of processor cores to be assigned for each task, default is 1, increase for multi-threaded runs
#SBATCH --ntasks-per-node=1     # number of tasks to be launched on each allocated node

############# LOADING MODULES (optional) #############
#module load apps/xxx
#module load libs/xxx

module load grenedalf/0.3.0

############# MY CODE #############

# Calculate the allele counts (ref and main alt) for all sites for the entire genome

grenedalf frequency \
--sync-path /home/jenni/data_folder/working_folder/selection_lines_poolseq/mpileup/Jan24_correct_order_syncs/SL_ps_only_Hcwbps18_q20Q30_noindel_lines_in_correct_order_Jan24.java.sync \
--reference-genome-fasta-file /home/jenni/data_folder/wbps18/haemonchus_contortus.PRJEB506.WBPS18.genomic.fa \
--rename-samples-file /home/jenni/data_folder/working_folder/selection_lines_poolseq/pairwise_comparisons/LINES_IN_CORRECT_ORDER/grenedalf/Jan24_F0_vs_all_sel_lines_only/all_sites/AF/SL_sample_names_for_full_genome.txt \
--write-sample-counts \
--separator-char tab \
--out-dir /home/jenni/data_folder/working_folder/selection_lines_poolseq/pairwise_comparisons/LINES_IN_CORRECT_ORDER/grenedalf/Jan24_F0_vs_all_sel_lines_only/all_sites/AF \
--file-prefix Full_genome_400_Gf_F0_vs_SL_only_Hcwbps18_q20Q30_noindels_AF_sample_counts \
--log-file grenedalf_AF_all_fullgenome.log
```
Set the above running - 4th Aug. 

Make the SNP file, without a header, ```Chr\tPos```: 

```
awk 'BEGIN {OFS=FS="\t"} {print $1,$2}' Full_genome_400_Gf_F0_vs_SL_only_Hcwbps18_q20Q30_noindels_AF_sample_countsfrequency.csv | tail -n +2 > SLPS_snps.txt
```
Transfer it to MARS. 

Sort the SNP location file by chromosome and position.
Hmmm.... Not at all sure exactly how well this will work. To try the file as is and see what it does.
```
sed 's/hcontortus_chr//g' SLPS_snps.txt | sed 's/_Celeg_TT_arrow_pilon//g' | sort -n -k1,1, -k2,2 > SLPS_snps_sorted.txt
```

Copy and put the python script in the analysis folder.
```
#!/bin/bash -l

############# SLURM SETTINGS #############
#SBATCH --job-name=mhap_allsites_entiregenome      # some descriptive job name of your choice
#SBATCH --output=%x-%J.out      # output file name will contain job name + job ID
#SBATCH --error=%x-%J.err       # error file name will contain job name + job ID
#SBATCH --time=0-48:00:00       # time limit for the whole run, in the form of d-hh:mm:ss, also accepts mm, mm:ss, hh:mm:ss, d-hh, d-hh:mm
#SBATCH --mem=20G                # memory required per node, in the form of [num][M|G|T]
#SBATCH --ntasks=1              # number of Slurm tasks to be launched, increase for multi-process runs ex. MPI
#SBATCH --cpus-per-task=1       # number of processor cores to be assigned for each task, default is 1, increase for multi-threaded runs
#SBATCH --ntasks-per-node=1     # number of tasks to be launched on each allocated node

############# LOADING MODULES (optional) #############
#module load apps/xxx
#module load libs/xxx

############# MY CODE #############
# pwd = /users/jmi45g/project0005/HCON/DATA/SLPS_Illumina/mhFromLowDepSeq_analysis

python3 calc_mh_freq.py -m ./Populations_bams.txt -s ./SLPS_snps.txt -o SLPS_mhap_25SNPs.txt -af -w 125 -pool -ms 25

python3 calc_mh_freq.py -m ./Populations_bams.txt -s ./SLPS_snps.txt -o SLPS_mhap_40SNPs.txt -af -w 125 -pool -ms 40
```

Was taking a very long time - 1 day and 4 h in had done the first 30 Mb of Chr4. Therefore, decided to:
1. Split the snp files by chromosomes and set running separate jobs
2. Prevent output of allele frequencies (-af) for all but Chr5 (not interested and makes files huge)

The 25 SNPs/window completed for some.... mito, Chr1, 2, 3 and X.  Chr5 looked near to completion. 
The 40 SNPS/window failed for all but the mito genome. Also seemed to have got through only a very tiny proportion of the genome. Not sure if this is due to having 40 SNPS/window or having perhaps reduced cpus??? 

However, got the output frequencies of SNPs per window for those that had completed and plotted these - using just 25 or just 40 SNPs is a huge underrepresentation of windows!! Decided to try using 76 SNPs per window, aware that this could take far far too long. There are no 'threaded' options in the github page that I can see, so have set running again with a single one. Have only started Chr1 and ChrV. Have decided to run ChrV without asking for individual AF to be output as file size is simply massive with them and would not be using them all anyway.

```
python3 calc_mh_freq.py -m ./Populations_bams.txt -s ./Chr5_SLPS_snps.txt -o Chr5_SLPS_mhap_76SNPs_noaf.txt -w 125 -pool -ms 76
```
Will also run just the locus on ChrV to try to speed things up and get answers for it by the end of the week! Will make a new SNPs file.
```
python3 calc_mh_freq.py -m ./Populations_bams.txt -s ./Ch5_locus_30-45Mb_SLPS_snps.txt -o Chr5_SLPS_mhap_76SNPs_noaf.txt -w 125 -pool -ms 76
```

I also then tried submitting the Chr5 locus for only the F3 generation (CTL, IVM and MOX): Populations_F3only_bams.txt

```
python3 calc_mh_freq.py -m ./Populations_F3only_bams.txt -s ./Ch5_locus_30-45Mb_SLPS_snps.txt -o Chr5_SLPS_mhap_76SNPs_locus_noaf_F3only.txt -w 125 -pool -ms 76

```



```
# To get snps:
grep -e 'chr1' SLPS_snps.txt > Chr1_SLPS_snps.txt
grep -e 'chr2' SLPS_snps.txt > Chr2_SLPS_snps.txt
grep -e 'chr3' SLPS_snps.txt > Chr3_SLPS_snps.txt
grep -e 'chr4' SLPS_snps.txt > Chr4_SLPS_snps.txt
grep -e 'chr5' SLPS_snps.txt > Chr5_SLPS_snps.txt
grep -e 'chrX' SLPS_snps.txt > ChrX_SLPS_snps.txt
grep -e 'mito' SLPS_snps.txt > Mito_SLPS_snps.txt
```


Ok, it is still taking a long time. I now have filtered bam files for the F3 generations (see above). I'm going to set running 15 jobs, each for a 1 Mb chunk of the locus.

To get the SNPs for each 1 Mb window:
```
# And get the SNPs for each 500 kb window.
# awk -v assigns a value to a variable before awk starts
for FILE in Ch5_locus_30-45Mb_SLPS_snps.txt ; 
do
for i in $(seq 30000000 1000000 45000000); do
    cat ${FILE} |\
    awk -v pos=$i 'BEGIN {FS=OFS="\t"} {if ($2>=pos && $2<pos+1000000) print $0}' > tmp_snps_${FILE}_${i} ;
    done ;
done &
```
Adjust the population bam list file:
```
sed 's/marked.bam/marked.dupsremoved.bam/g' Populations_F3only_bams.txt | sed 's/Aligned_bam_files\/bam_files/mhFromLowDepSeq_analysis/g' > Populations_F3only_FILTEREDbams.txt
```

And submit: 
```
#!/bin/bash -l

############# SLURM SETTINGS #############
#SBATCH --job-name=mhap_allsites_entiregenome      # some descriptive job name of your choice
#SBATCH --output=%x-%J.out      # output file name will contain job name + job ID
#SBATCH --error=%x-%J.err       # error file name will contain job name + job ID
#SBATCH --time=3-00:00:00       # time limit for the whole run, in the form of d-hh:mm:ss, also accepts mm, mm:ss, hh:mm:ss, d-hh, d-hh:mm
#SBATCH --mem=7G                # memory required per node, in the form of [num][M|G|T]
#SBATCH --ntasks=1              # number of Slurm tasks to be launched, increase for multi-process runs ex. MPI
#SBATCH --cpus-per-task=1       # number of processor cores to be assigned for each task, default is 1, increase for multi-threaded runs
#SBATCH --ntasks-per-node=1     # number of tasks to be launched on each allocated node

############# LOADING MODULES (optional) #############
#module load apps/xxx
#module load libs/xxx

############# MY CODE #############
# pwd = /users/jmi45g/project0005/HCON/DATA/SLPS_Illumina/mhFromLowDepSeq_analysis


# Ran a separate script for each of the following:



python3 calc_mh_freq.py -m ./Populations_F3only_FILTEREDbams.txt -s ./tmp_snps_Ch5_locus_30-45Mb_SLPS_snps.txt_30000000 -o Chr5_SLPS_mhap_76SNPs_locus_noaf_F3only_30000000.txt -w 125 -pool -ms 76

python3 calc_mh_freq.py -m ./Populations_F3only_FILTEREDbams.txt -s ./tmp_snps_Ch5_locus_30-45Mb_SLPS_snps.txt_31000000 -o Chr5_SLPS_mhap_76SNPs_locus_noaf_F3only_31000000.txt -w 125 -pool -ms 76

python3 calc_mh_freq.py -m ./Populations_F3only_FILTEREDbams.txt -s ./tmp_snps_Ch5_locus_30-45Mb_SLPS_snps.txt_32000000 -o Chr5_SLPS_mhap_76SNPs_locus_noaf_F3only_32000000.txt -w 125 -pool -ms 76

python3 calc_mh_freq.py -m ./Populations_F3only_FILTEREDbams.txt -s ./tmp_snps_Ch5_locus_30-45Mb_SLPS_snps.txt_33000000 -o Chr5_SLPS_mhap_76SNPs_locus_noaf_F3only_33000000.txt -w 125 -pool -ms 76

python3 calc_mh_freq.py -m ./Populations_F3only_FILTEREDbams.txt -s ./tmp_snps_Ch5_locus_30-45Mb_SLPS_snps.txt_34000000 -o Chr5_SLPS_mhap_76SNPs_locus_noaf_F3only_34000000.txt -w 125 -pool -ms 76

python3 calc_mh_freq.py -m ./Populations_F3only_FILTEREDbams.txt -s ./tmp_snps_Ch5_locus_30-45Mb_SLPS_snps.txt_35000000 -o Chr5_SLPS_mhap_76SNPs_locus_noaf_F3only_35000000.txt -w 125 -pool -ms 76

python3 calc_mh_freq.py -m ./Populations_F3only_FILTEREDbams.txt -s ./tmp_snps_Ch5_locus_30-45Mb_SLPS_snps.txt_36000000 -o Chr5_SLPS_mhap_76SNPs_locus_noaf_F3only_36000000.txt -w 125 -pool -ms 76

python3 calc_mh_freq.py -m ./Populations_F3only_FILTEREDbams.txt -s ./tmp_snps_Ch5_locus_30-45Mb_SLPS_snps.txt_37000000 -o Chr5_SLPS_mhap_76SNPs_locus_noaf_F3only_37000000.txt -w 125 -pool -ms 76

python3 calc_mh_freq.py -m ./Populations_F3only_FILTEREDbams.txt -s ./tmp_snps_Ch5_locus_30-45Mb_SLPS_snps.txt_38000000 -o Chr5_SLPS_mhap_76SNPs_locus_noaf_F3only_38000000.txt -w 125 -pool -ms 76

python3 calc_mh_freq.py -m ./Populations_F3only_FILTEREDbams.txt -s ./tmp_snps_Ch5_locus_30-45Mb_SLPS_snps.txt_39000000 -o Chr5_SLPS_mhap_76SNPs_locus_noaf_F3only_39000000.txt -w 125 -pool -ms 76

python3 calc_mh_freq.py -m ./Populations_F3only_FILTEREDbams.txt -s ./tmp_snps_Ch5_locus_30-45Mb_SLPS_snps.txt_40000000 -o Chr5_SLPS_mhap_76SNPs_locus_noaf_F3only_40000000.txt -w 125 -pool -ms 76

python3 calc_mh_freq.py -m ./Populations_F3only_FILTEREDbams.txt -s ./tmp_snps_Ch5_locus_30-45Mb_SLPS_snps.txt_41000000 -o Chr5_SLPS_mhap_76SNPs_locus_noaf_F3only_41000000.txt -w 125 -pool -ms 76

python3 calc_mh_freq.py -m ./Populations_F3only_FILTEREDbams.txt -s ./tmp_snps_Ch5_locus_30-45Mb_SLPS_snps.txt_42000000 -o Chr5_SLPS_mhap_76SNPs_locus_noaf_F3only_42000000.txt -w 125 -pool -ms 76

python3 calc_mh_freq.py -m ./Populations_F3only_FILTEREDbams.txt -s ./tmp_snps_Ch5_locus_30-45Mb_SLPS_snps.txt_43000000 -o Chr5_SLPS_mhap_76SNPs_locus_noaf_F3only_43000000.txt -w 125 -pool -ms 76

python3 calc_mh_freq.py -m ./Populations_F3only_FILTEREDbams.txt -s ./tmp_snps_Ch5_locus_30-45Mb_SLPS_snps.txt_44000000 -o Chr5_SLPS_mhap_76SNPs_locus_noaf_F3only_44000000.txt -w 125 -pool -ms 76
```


# Plotting the results. 

(Originally I had the Expected Heterozygosity for each population and the Allele frequencies for each haplotype for each population as so: 
```
Chr     Pos     F0      F1_C_L1 F1_C_L2 F1_C_L3 F2_C_L1 F2_C_L2 F2_C_L3 F3_C_L1 F3_C_L2 F3_C_L3 F1_I_L2 F1_I_L3 F2_I_L1 F2_I_L2 F2_I_L3 F3_I_L1 F3_I_L2 F3_I_L3 F1_M_L1 F1_M_L2 F1_M_L3 F2_M_L1 F2_M_L2 F2_M_L3 F3_M_L1 F3_M_L2 F3_M_L3 NumReads_F0     NumReads_F1_C_L1        NumReads_F1_C_L2        NumReads_F1_C_L3        NumReads_F2_C_L1        NumReads_F2_C_L2        NumReads_F2_C_L3        NumReads_F3_C_L1        NumReads_F3_C_L2        NumReads_F3_C_L3        NumReads_F1_I_L2        NumReads_F1_I_L3        NumReads_F2_I_L1        NumReads_F2_I_L2        NumReads_F2_I_L3        NumReads_F3_I_L1        NumReads_F3_I_L2        NumReads_F3_I_L3        NumReads_F1_M_L1        NumReads_F1_M_L2        NumReads_F1_M_L3        NumReads_F2_M_L1        NumReads_F2_M_L2        NumReads_F2_M_L3        NumReads_F3_M_L1        NumReads_F3_M_L2        NumReads_F3_M_L3        AlleleFreq_F0   AlleleFreq_F1_C_L1      AlleleFreq_F1_C_L2      AlleleFreq_F1_C_L3      AlleleFreq_F2_C_L1      AlleleFreq_F2_C_L2      AlleleFreq_F2_C_L3      AlleleFreq_F3_C_L1      AlleleFreq_F3_C_L2      AlleleFreq_F3_C_L3      AlleleFreq_F1_I_L2      AlleleFreq_F1_I_L3      AlleleFreq_F2_I_L1      AlleleFreq_F2_I_L2      AlleleFreq_F2_I_L3      AlleleFreq_F3_I_L1      AlleleFreq_F3_I_L2      AlleleFreq_F3_I_L3      AlleleFreq_F1_M_L1      AlleleFreq_F1_M_L2      AlleleFreq_F1_M_L3      AlleleFreq_F2_M_L1      AlleleFreq_F2_M_L2      AlleleFreq_F2_M_L3      AlleleFreq_F3_M_L1      AlleleFreq_F3_M_L2      AlleleFreq_F3_M_L3
hcontortus_chr4_Celeg_TT_arrow_pilon    1197,1198,1208,1209,1220,1223,1226,1251,1298,1299,1307,1311,1316        0.0     NA      0.014318697048328222    0.11729925176094358     0.15649458094692692     0.060581588379464146    0.1321781340771515      NA      0.5172383437903809      0.04167332891407638     0.5468683625080613      0.0
     1.8504219314596781e-06  0.08040206366937341     0.32723705193303587     0.0     0.07048245546554721     0.15665126635450521     1.3974857144827268e-06  0.19635149807751362     0.01432499365317641     0.0     NA      0.042546002369478364    0.05647253071577485     0.16142175428358885     0.15579383320396667     100     100     100     100     100     100     100     100     100     100     100     100     100     100     100     83      100     100     100     100     100     100     100     100     100     100     100     CTCCCACCCTATC:1.0       NA      CTCCCACCCTATC:0.9927886478763851,CTCCCACCCTATT:0.007211352123616066     CTCCCACCCTATC:0.9374361371897939,CTCCCACGCTATC:0.06256386281020658      CTCCCACCCTATC:0.9151475264705856,CTCCCTCCCTATC:0.07714264810387596,CTCCCACCCTATT:0.007709454664803822,CTCCCTCCCTATT:2.259961600891298e-19,CTCCCACCCTATG:3.671098700160269e-07,CTCCCTCCCTATG:3.6508651601503047e-09      CTCCCACCCTATC:0.9687315125938991,AGCCCACCAGATC:1.2163368888573622e-10,CTCCCACCCTATG:8.022495979733337e-07,AGCCCACCAGATG:7.104127715804271e-19,CTCCCACCCTATT:1.7911846166697315e-13,AGCCCACCAGATT:0.03126768503468921    CTCCCACCCTATC:0.929434527004487,CTCCCACTCTATC:0.06251874124628001,CTCCCACCCTATT:0.008045677514566162,CTCCCACTCTATT:1.5142086435053437e-19,CTCCCACCCTATA:1.051791198185782e-06,CTCCCACTCTATA:2.443468466785587e-09       NA      CTCCCACCCTATC:0.6513471973808659,CACCCACCCTATC:0.18663180046308672,CTCTCACCCTATC:0.1537481600022613,CTCCCACCCTATT:0.005609354714489579,CACCCACCCTATT:0.0026634874175004287,CTCTCACCCTATT:7.668333491916033e-24,CTCCCACCCTATA:8.24814604300612e-12,CACCCACCCTATA:6.645702304660402e-12,CTCTCACCCTATA:6.903176812199657e-12       CTCCCACCCTATC:0.9787100782294496,AGCCCACCAGATC:1.608800681059387e-07,CTCCCACCCTATT:1.427246697915508e-09,AGCCCACCAGATT:0.021289759463235483     CTCCCACCCTATC:0.6535995905933943,CTGCCACCCTATC:0.09970142028578727,AGGCCACCCTATC:0.03244369935299324,CTCTCACCCTATC:0.08495515428953163,AGCTCACCCTATC:0.027645417884106994,CTCCCTCCCTATC:0.08080867141393198,AGCCCACCAGATC:1.3286266038389324e-08,CTCCCACCCTATT:1.7231847750534239e-12,CTGCCACCCTATT:2.1669700951023653e-20,AGGCCACCCTATT:6.567137725768135e-21,CTCTCACCCTATT:1.8487607859462827e-20,AGCTCACCCTATT:5.6028260258215995e-21,CTCCCTCCCTATT:1.7620715412279514e-20,AGCCCACCAGATT:0.020846032892265685        CTCCCACCCTATC:1.0       CTCCCACCCTATC:0.9999990747881783,CTCCCACCCTATA:9.252118223431084e-07    CTCCCACCCTATC:0.9583570807878364,CTACCACCCTATC:0.032837654619267,CTCCCACCCTATT:0.0003675281296872293,CTACCACCCTATT:0.008437736463210391 CTCCCACCCTATC:0.8088213587241124,CGCCCACCCTATC:0.10776144655327441,TTTCAACCCTATC:0.08341719472261501    CTCCCACCCTATC:1.0       CTCCCACCCTATC:0.963420729216148,CTCCCCCCCTATC:0.03657927078385398       CTCCCACCCTATC:0.9150322820345642,CTCCCATCCTATC:0.0775188890468351,CTCCCACCCTATT:0.007448377062376928,CTCCCATCCTATT:3.2083366391072567e-18,CTCCCACCCTATA:4.399120209091651e-07,CTCCCATCCTATA:1.194420212946444e-08       CTCCCACCCTATC:0.9999993012566545,CTCCCACCCTATA:6.987433438996381e-07    CTCCCACCCTATC:0.8919394542724974,CTCCCATCCTATC:0.08760092140374966,CTCCCACCAGATC:1.422679520309345e-07,CTCCCACCCTATT:5.454501732018073e-11,CTCCCATCCTATT:5.24741031345241e-18,CTCCCACCAGATT:0.02045948200125612 CTCCCACCCTATC:0.9927854534920971,CTCCCACCCTATT:0.0072145465079041125    CTCCCACCCTATC:1.0       NA      CTCCCACCCTATC:0.9782541210525302,CTCCCACCAGATC:2.3747710533114315e-07,CTCCCACCCTATT:1.4916807552913437e-08,CTCCCACCAGATT:0.021745626553555684
   CTCCCACCCTATC:0.9710900355395998,CTCCCAGCCTATC:0.02130066268499601,CTCCCACCCTATA:4.954650356115765e-07,CTCCCAGCCTATA:3.0295428891074934e-10,CTCCCACCCTATT:0.0076088060074146055,CTCCCAGCCTATT:3.051429077893874e-24     CTTCCACCCTATC:0.0885524056964166,CTCCCACCCTATC:0.9114475943035828       CTCCCACCCTATC:0.9164534497048215,CTCCTACCCTATC:0.06213570864616992,CTCCCACCAGATC:4.0160923330069296e-10,CTCCCACCCTATA:1.9890369601086233e-07,CTCCTACCCTATA:4.804746805726271e-07,CTCCCACCAGATA:5.998488281086648e-16,CTCCCACCCTATT:5.3687439155764244e-15,CTCCTACCCTATT:7.565892153312685e-12,CTCCCACCAGATT:0.02141016186145171
```

Now I do not have the allele frequencies:
```
Chr     Pos     F3_C_L1 F3_C_L2 F3_C_L3 F3_I_L1 F3_I_L2 F3_I_L3 F3_M_L1 F3_M_L2 F3_M_L3 NumReads_F3_C_L1        NumReads_F3_C_L2        NumReads_F3_C_L3   NumReads_F3_I_L1        NumReads_F3_I_L2        NumReads_F3_I_L3        NumReads_F3_M_L1        NumReads_F3_M_L2        NumReads_F3_M_L3
hcontortus_chr5_Celeg_TT_arrow_pilon    37000001,37000002,37000003,37000004,37000006,37000007,37000009,37000012,37000013,37000014,37000015,37000017,37000019,37000034,37000039,37000042,37000049,37000052,37000054,37000057,37000059,37000062,37000063,37000064,37000066,37000067,37000068,37000069,37000072,37000073,37000085,37000086,37000087,37000088,37000090,37000091   0.4599850300197168      0.5440728436992743      0.5212307680406698 0.4288569471534275      0.3371671285986553      0.4336674278973107      0.44357301798541227     0.3980835694304813      0.5042734231122581 100     100     100     100     100     100     100     100     100

```
Extract just the Chr, the Positions and the Expected Heterozygosity for each sample:
```
awk '{print NF}' file
# 20, therefore have (20-2)/2 = 9 samples

# Select the Expected heterozygosity columns and rename the chromosomes to make shorter for plotting:
for i in Chr5*txt ;
do
cut -f1-11 $i | sed 's/hcontortus_chr//g' | sed 's/_Celeg_TT_arrow_pilon//g' > HE_${i}  ;
done

# The program calculates microhaplotypes in overlapping sliding windows. To select every 125th line in the files so as to select each microhaplotype "region" at most once. This will lose some entirely, but that is fine.

for i in HE_* ;
do
awk 'NR % 125 == 0' $i > q125_${i} ;
done


# Get the first and last snp for each microhaplotype, to enable plotting later, and count the number of SNPs in each microhaplotype:
# Extract the snp positions:
for i in q125* ;
do
awk 'BEGIN {FS=OFS="\t"} {print $2}' $i | awk 'BEGIN {FS=","}{OFS="\t"} {print $0}' > SNPs_in_mhs_${i} ;
done

# Count the number of snps in each microhaplotype:
for i in SNPs_in_mhs* ;
do
awk 'BEGIN {FS=","}{print NF}' $i > N${i} ;
done

# Extract the first and last SNP and calculate a mid-point for plotting:
for i in SNPs_in_mhs* ;
do
cat $i | awk 'BEGIN {FS=","}{OFS="\t"}{print $1, $NF, (($NF + $1)/2) }' > FnL_${i} ;
done

# Edit the header to read Start\tEnd\tMid

for i in FnL* ;
do
cat header $i > head_${i} ;
done

# Add the original header back on the q125 files:
# (Chr     Pos     F3_C_L1 F3_C_L2 F3_C_L3 F3_I_L1 F3_I_L2 F3_I_L3 F3_M_L1 F3_M_L2 F3_M_L3)
for i in q125* ;
do
cat HE_header $i > head_${i} ;
done

# And add a header to the NSNPs files:
for i in NSNPs* ;
do
cat NSNPs_header $i > head_${i} ;
done



# Paste the new columns back onto the old:
for i in {0..9} ;
do 
paste head_NSNPs_in_mhs_q125_HE_Chr5_SLPS_mhap_76SNPs_locus_noaf_F3only_3${i}* head_FnL_SNPs_in_mhs_q125_HE_Chr5_SLPS_mhap_76SNPs_locus_noaf_F3only_3${i}* head_q125_HE_Chr5_SLPS_mhap_76SNPs_locus_noaf_F3only_3${i}* | cut -f1-5,7-> ForRplot_HE_SLPS_mhap_76SNPs_3${i}Mb.txt ;
done

for i in {0..4} ;
do 
paste head_NSNPs_in_mhs_q125_HE_Chr5_SLPS_mhap_76SNPs_locus_noaf_F3only_4${i}* head_FnL_SNPs_in_mhs_q125_HE_Chr5_SLPS_mhap_76SNPs_locus_noaf_F3only_4${i}* head_q125_HE_Chr5_SLPS_mhap_76SNPs_locus_noaf_F3only_4${i}* | cut -f1-5,7-> ForRplot_HE_SLPS_mhap_76SNPs_4${i}Mb.txt ;
done

# Combine the files:
for i in ForR* ;
do 
cat $i | tail -n +2 >> ForRplot_Chr5_30-45Mb_HE_SLPS_mhap_76SNPs.txt ;
done

head -n 1 ForRplot_HE_SLPS_mhap_76SNPs_30Mb.txt > final_header

cat final_header ForRplot_Chr5_30-45Mb_HE_SLPS_mhap_76SNPs.txt > FINAL_ForRplot_Chr5_30-45Mb_HE_SLPS_mhap_76SNPs.txt

```


1. Plot the Expected Heterozygosity:

```
module load apps/R/4.4.1/gcc-8.5.0+openblas-0.3.28

R

library(stringr)
library(patchwork)
library(ggplot2)
library(dplyr)
library(reshape2)
library(viridis)

#[1] reshape2_1.4.4  dplyr_1.1.4     ggplot2_3.5.2   patchwork_1.3.1
#[5] stringr_1.5.1  viridis_0.6.5

df<-read.table("FINAL_ForRplot_Chr5_30-45Mb_HE_SLPS_mhap_76SNPs.txt", header=TRUE, sep="\t", na.strings = "NA")
summary(df)
```
```
    NSNPs           Start               End                Mid          
 Min.   : 3.00   Min.   :30000800   Min.   :30000923   Min.   :30000900  
 1st Qu.:36.00   1st Qu.:33122240   1st Qu.:33122364   1st Qu.:33122302  
 Median :47.00   Median :36496540   Median :36496663   Median :36496600  
 Mean   :46.06   Mean   :36942424   Mean   :36942547   Mean   :36942485  
 3rd Qu.:57.00   3rd Qu.:40819230   3rd Qu.:40819353   3rd Qu.:40819300  
 Max.   :76.00   Max.   :44999872   Max.   :44999987   Max.   :44999900  
                                                                         
      Chr       F3_C_L1          F3_C_L2          F3_C_L3      
 Min.   :5   Min.   :0.0000   Min.   :0.0000   Min.   :0.0000  
 1st Qu.:5   1st Qu.:0.2265   1st Qu.:0.2341   1st Qu.:0.2273  
 Median :5   Median :0.4580   Median :0.4664   Median :0.4643  
 Mean   :5   Mean   :0.4249   Mean   :0.4320   Mean   :0.4283  
 3rd Qu.:5   3rd Qu.:0.6158   3rd Qu.:0.6197   3rd Qu.:0.6145  
 Max.   :5   Max.   :0.9745   Max.   :0.9537   Max.   :0.9580  
             NA's   :30       NA's   :41       NA's   :44      
    F3_I_L1          F3_I_L2          F3_I_L3          F3_M_L1      
 Min.   :0.0000   Min.   :0.0000   Min.   :0.0000   Min.   :0.0000  
 1st Qu.:0.2454   1st Qu.:0.2249   1st Qu.:0.2354   1st Qu.:0.4192  
 Median :0.5134   Median :0.4298   Median :0.4727   Median :0.6250  
 Mean   :0.4502   Mean   :0.4041   Mean   :0.4290   Mean   :0.5599  
 3rd Qu.:0.6287   3rd Qu.:0.5542   3rd Qu.:0.5907   3rd Qu.:0.7410  
 Max.   :0.9811   Max.   :0.9640   Max.   :0.9832   Max.   :0.9881  
 NA's   :337      NA's   :369      NA's   :221      NA's   :289     
    F3_M_L2          F3_M_L3      
 Min.   :0.0000   Min.   :0.0000  
 1st Qu.:0.1712   1st Qu.:0.1199  
 Median :0.3462   Median :0.2782  
 Mean   :0.3616   Mean   :0.3200  
 3rd Qu.:0.5362   3rd Qu.:0.5062  
 Max.   :0.9853   Max.   :0.9815  
 NA's   :389      NA's   :2500  
```
```
dfmelt<-(melt(df, id.vars=c("NSNPs","Start","End","Mid","Chr")))
head(dfmelt)

write.table(dfmelt, file="FINAL_ForRplot_Chr5_30-45Mb_HE_SLPS_mhap_76SNPs_MELTED.txt", row.names=FALSE, sep="\t", quote= FALSE) 

# Rename chromosomes:
dfmelt$Chr<-gsub("1","I", dfmelt$Chr)
dfmelt$Chr<-gsub("2","II", dfmelt$Chr)
dfmelt$Chr<-gsub("3","III", dfmelt$Chr)
dfmelt$Chr<-gsub("4","IV", dfmelt$Chr)
dfmelt$Chr<-gsub("5","V", dfmelt$Chr)


# Plot all data:
# binwidth=c(x-axis-units,y-axis-units)
ALL<- ggplot(data=dfmelt, mapping=aes(x=Mid/1e6, y=value)) +
geom_hex(binwidth=c(0.1,0.1)) +
scale_fill_viridis()+
labs(y="Expected Heterozygosity", x="Genomic position (Mb)")+
facet_grid(variable ~ Chr, scales="free_x", switch="x")+
scale_x_continuous(
    breaks = seq(0, 55, 1),
    minor_breaks = seq(1, 55, 1)
  )+
theme_light()+
theme(axis.title.x = element_text(size=10), axis.text.x = element_text(size=6), axis.title.y = element_text(size=10), axis.text.y = element_text(size=4), strip.text.x = element_text(size = 6, face = "bold", margin = margin(0.5,0,0.5,0, "mm")), strip.text.y = element_text(size = 6, face = "bold", angle=0), plot.title = element_text(size=10, hjust=0.5), legend.position = "top", legend.title = element_text(size=8), legend.text = element_text(size=6), legend.key.height=unit(0.2,"cm"), panel.spacing.x = unit(0.2, "mm"), panel.spacing.y = unit(0.75, "mm"))+
ggtitle('Chromosome V:30-45 Mb Microhaplotypes Expected Heterozygosity')

patch<- ALL + plot_annotation(tag_levels = 'A') + plot_layout (axis_titles = "collect")
patch

ggsave("Chr5_30-45Mb_HE_SLPS_mhap_76SNPs_p1p1.tiff", patch, device=tiff, width=2250, height=2250, units="px", dpi=320)

```
Note that the MOX_L3 had a lot more NAs than the others - NAs can be due to missing data or too many haplotypes: https://github.com/delomast/mhFromLowDepSeq/issues/10. 

<img width="2250" height="2250" alt="image" src="https://github.com/user-attachments/assets/469346e9-bfaa-47c3-a6a6-1f795910d298" />

<img width="2250" height="2250" alt="image" src="https://github.com/user-attachments/assets/8b4faa69-a679-4f5b-86b1-f3c007b415cd" />

Ok, it reflects what we already know, although not quite sure exactly how the increased het towards the left ties in with Tajima's D. 

There are clearly far too many micro haps to accurately identify whether something is being selected across all. However, to try looking at single copy genes (remember have filtered out MAPQ=0 reads, so protein kinases, repetitive regions etc may not be as accurate). 

1. Cky-1: HCON_00155390: Scaffold hcontortus_chr5_Celeg_TT_arrow_pilon: 37,487,982-37,497,398
```
# Filter the SNPs file to include only snps at this locus.

awk 'BEGIN {FS=OFS="\t"} {if ($2>=37487981 && $2<37497399) print $0}' Ch5_locus_30-45Mb_SLPS_snps.txt > Cky-1_SLPS_snps.txt

# Have 2402 SNPs in total. 
```
Note that I will allow all haplotypes to be included in the model as above, as I think that otherwise they are taken OUT of the model itself, rather than just the output file. But clearly I am looking for the most frequent ones. 
```
#!/bin/bash -l

############# SLURM SETTINGS #############
#SBATCH --job-name=mhap_allsites_cky1      # some descriptive job name of your choice
#SBATCH --output=%x-%J.out      # output file name will contain job name + job ID
#SBATCH --error=%x-%J.err       # error file name will contain job name + job ID
#SBATCH --time=0-12:00:00       # time limit for the whole run, in the form of d-hh:mm:ss, also accepts mm, mm:ss, hh:mm:ss, d-hh, d-hh:mm
#SBATCH --mem=7G                # memory required per node, in the form of [num][M|G|T]
#SBATCH --ntasks=1              # number of Slurm tasks to be launched, increase for multi-process runs ex. MPI
#SBATCH --cpus-per-task=1       # number of processor cores to be assigned for each task, default is 1, increase for multi-threaded runs
#SBATCH --ntasks-per-node=1     # number of tasks to be launched on each allocated node

############# LOADING MODULES (optional) #############
#module load apps/xxx
#module load libs/xxx

############# MY CODE #############
# pwd = /users/jmi45g/project0005/HCON/DATA/SLPS_Illumina/mhFromLowDepSeq_analysis

python3 calc_mh_freq.py -m ./Populations_F3only_FILTEREDbams.txt -s ./Cky-1_SLPS_snps.txt -o Cky1_SLPS_mhap_76SNPs_F3only.txt -w 125 -pool -ms 76 -af
```

2. HCON_00158050, a serpentine receptor (class E). Has orthologues but no paralogues. Scaffold hcontortus_chr5_Celeg_TT_arrow_pilon: 41,309,214-41,318,844

It sits about 41.3 Mb so at the region of slightly increased heterozygosity/increased haplotype count possibly, with apparent reduction in heterozygosity in the controls.

```
# Filter the SNPs file to include only snps at this locus.

awk 'BEGIN {FS=OFS="\t"} {if ($2>=41309214 && $2<41318844) print $0}' Ch5_locus_30-45Mb_SLPS_snps.txt > HCON_00158050_SLPS_snps.txt
# 3729 SNPs

python3 calc_mh_freq.py -m ./Populations_F3only_FILTEREDbams.txt -s ./HCON_00158050_SLPS_snps.txt -o HCON_00158050_SLPS_mhap_76SNPs_F3only.txt -w 125 -pool -ms 76 -af

```
```
# The program calculates microhaplotypes in overlapping sliding windows. To select every 125th line in the files so as to select each microhaplotype "region" at most once. This will lose some entirely, but that is fine.

awk 'NR % 125 == 0' Cky1_SLPS_mhap_76SNPs_F3only.txt > q125_Cky1_SLPS_mhap_76SNPs_F3only.txt
# Gives 9 microhaplotypes

awk 'NR % 125 == 0' HCON_00158050_SLPS_mhap_76SNPs_F3only.txt > q125_HCON_00158050_SLPS_mhap_76SNPs_F3only.txt
# Gives 15 microhaplotypes

```
Ok, I want to figure out whether:
1. There is a shift in the allele frequency for the selected lines relative to the control lines.
2. Whether the SAME allele is being selected across all lines

To do this I need to get the file usuable. 

```
# Get the first and last snp for each microhaplotype, to enable plotting later, and count the number of SNPs in each microhaplotype:
# Extract the snp positions:
for i in q125* ;
do
awk 'BEGIN {FS=OFS="\t"} {print $2}' $i | awk 'BEGIN {FS=","}{OFS="\t"} {print $0}' > SNPs_in_mhs_${i} ;
done

# Count the number of snps in each microhaplotype:
for i in SNPs_in_mhs* ;
do
awk 'BEGIN {FS=","}{print NF}' $i > N${i} ;
done

# Extract the first and last SNP and calculate a mid-point for plotting:
for i in SNPs_in_mhs* ;
do
cat $i | awk 'BEGIN {FS=","}{OFS="\t"}{print $1, $NF, (($NF + $1)/2) }' > FnL_${i} ;
done
```
## For cky1:
```
paste FnL_SNPs_in_mhs_q125_Cky1_SLPS_mhap_76SNPs_F3only.txt NSNPs_in_mhs_q125_Cky1_SLPS_mhap_76SNPs_F3only.txt q125_Cky1_SLPS_mhap_76SNPs_F3only.txt | cut -f1-5,25- > ForRplot_q125_Cky1_SLPS_mhap_76SNPs_F3only.txt

sed 's/hcontortus_chr//g' ForRplot_q125_Cky1_SLPS_mhap_76SNPs_F3only.txt | sed 's/_Celeg_TT_arrow_pilon//g' > new_ForRplot_q125_Cky1_SLPS_mhap_76SNPs_F3only.txt

cat header_af new_ForRplot_q125_Cky1_SLPS_mhap_76SNPs_F3only.txt > FINAL_ForRplot_q125_Cky1_SLPS_mhap_76SNPs_F3only.txt

```

```

module load apps/R/4.4.1/gcc-8.5.0+openblas-0.3.28

R

library(stringr)
library(patchwork)
library(ggplot2)
library(dplyr)
library(reshape2)
library(viridis)

#[1] reshape2_1.4.4  dplyr_1.1.4     ggplot2_3.5.2   patchwork_1.3.1
#[5] stringr_1.5.1  viridis_0.6.5

df<-read.table("FINAL_ForRplot_q125_Cky1_SLPS_mhap_76SNPs_F3only.txt", header=TRUE, sep="\t", na.strings = "NA")
summary(df)

dfmelt<-melt(df, id.vars=c("NSNPs","Start","End","Mid","Chr"))
head(dfmelt)

write.table(dfmelt, file="FINAL_ForRplot_q125_Cky1_SLPS_mhap_76SNPs_F3only_MELTED.txt", row.names=FALSE, sep="\t", quote= FALSE)
```

Then, expand the allele column and sort by the start position (as output by sample)
```
sed 's/,/\t/g' FINAL_ForRplot_q125_Cky1_SLPS_mhap_76SNPs_F3only_MELTED.txt | sort -k 2 > Sorted_FINAL_ForRplot_q125_Cky1_SLPS_mhap_76SNPs_F3only_MELTED.txt

# Very confusing to look at. Let's get the first allele for each (Note - not necessarily the most frequent allele for a sample! and compare.

cut -f4,6,7 Sorted_FINAL_ForRplot_q125_Cky1_SLPS_mhap_76SNPs_F3only_MELTED.txt > First_allele_cky1.txt

sed 's/:/\t/g' First_allele_cky1.txt | sort -k 1 | cut -f1-3 | uniq -c -f2 > out
```

Do - get all the alleles and their freq in separate columns. then melt and then plot the freq as a dot for each allele, coloured by the sample, faceted by the position along the x-axis..... 
Could get to be a very wide plot....... BUT IT WOULD WORK! :D 

So, going back to the ```FINAL_ForRplot_q125_Cky1_SLPS_mhap_76SNPs_F3only.txt``` need to select each sample individually, then spread out and replicate the header. 

```
for i in {6..14} ;
do
sed 's/AlleleFreq_//g' FINAL_ForRplot_q125_Cky1_SLPS_mhap_76SNPs_F3only.txt | cut -f$i | sed 's/,/\t/g' > Col${i} ;
done

for i in Col* ;
do
awk '{print FILENAME "\t" NF}' $i >> ncols ;
done

# I then manually changed the headers to have the max number of columns each sample required

cut -f1-5 FINAL_ForRplot_q125_Cky1_SLPS_mhap_76SNPs_F3only.txt > start

paste start Col* > Updated_ForRplot_q125_Cky1_SLPS_mhap_76SNPs_F3only.txt

```
Struggling with having empty columns on rows but with nothing in them... 
Try this instead. 

```
for i in {6..14} ;
do
awk 'BEGIN {OFS=FS="\t"} 
       NR==1 {max=NF} 
             {for(i=NF+1;i<=max;i++) $i="null"}1' Col$i > Col${i}_new ;
done
```
PERFECT! :D 

Let's rejoin them all.... 

```
paste start Col*new > Updated_ForRplot_q125_Cky1_SLPS_mhap_76SNPs_F3only.txt

```

```
module load apps/R/4.4.1/gcc-8.5.0+openblas-0.3.28

R

library(stringr)
library(patchwork)
library(ggplot2)
library(dplyr)
library(reshape2)
library(viridis)

#[1] reshape2_1.4.4  dplyr_1.1.4     ggplot2_3.5.2   patchwork_1.3.1
#[5] stringr_1.5.1  viridis_0.6.5

df<-read.table("Updated_ForRplot_q125_Cky1_SLPS_mhap_76SNPs_F3only.txt", header=TRUE, sep="\t", na.strings = "null")
summary(df)

dfmelt<-melt(df, id.vars=c("NSNPs","Start","End","Mid","Chr"))
head(dfmelt)

write.table(dfmelt, file="Updated_ForRplot_q125_Cky1_SLPS_mhap_76SNPs_F3only_MELTED.txt", row.names=FALSE, sep="\t", quote= FALSE)

```
Then, replace the : with a tab: 
```
sed 's/:/\t/g' Updated_ForRplot_q125_Cky1_SLPS_mhap_76SNPs_F3only_MELTED.txt > New_ForRplot_q125_Cky1_SLPS_mhap_76SNPs_F3only_MELTED.txt

# Manually updated the header line.

# Remove any lines with NA as these will now be shorter than the others and I cannot fix it for each haplotype - the data is simply missing. Also remove any haplotype for a sample with <0.01 frequency just to try to make the graph more plottable. 

awk 'BEGIN {OFS=FS="\t"} $7!="NA" {print $0}'  New_ForRplot_q125_Cky1_SLPS_mhap_76SNPs_F3only_MELTED.txt | awk 'BEGIN {OFS=FS="\t"} $8>=0.01 {print $0}' > Finally_New_ForRplot_q125_Cky1_SLPS_mhap_76SNPs_F3only_MELTED.txt

# Note that removing freq <0.01 takes it from 439 to 295 lines.
```


```
module load apps/R/4.4.1/gcc-8.5.0+openblas-0.3.28

R

library(stringr)
library(patchwork)
library(ggplot2)
library(dplyr)
library(reshape2)
library(viridis)

#[1] reshape2_1.4.4  dplyr_1.1.4     ggplot2_3.5.2   patchwork_1.3.1
#[5] stringr_1.5.1  viridis_0.6.5

df<-read.table("Finally_New_ForRplot_q125_Cky1_SLPS_mhap_76SNPs_F3only_MELTED.txt", header=TRUE, sep="\t", na.strings = "null")
summary(df)

# Then, replace the .2 to .xx with gsub():

df$variable.new<-gsub("\\...","", df$variable) # Replace everything that is .10 upwards
df$variable.new<-gsub("\\..","", df$variable.new) # Replace everything that is .2 to .9. Have to do this AFTER The one above or it doesn't work. 
df$variable.new<-as.factor(df$variable.new)
df$Haplotype<-as.factor(df$Haplotype)

plot<-ggplot(data=df, aes(x=Haplotype, y=Frequency))+
geom_col(aes(fill=variable.new))+
scale_fill_viridis_d()+
labs(y="Haplotype frequency", x="Haplotype")+
theme_light()+
theme(axis.title.x = element_text(size=10), axis.text.x = element_blank(), axis.title.y = element_text(size=10), axis.text.y = element_text(size=6), strip.text.x = element_text(size = 6, face = "bold", margin = margin(0.5,0,0.5,0, "mm")), strip.text.y = element_text(size = 6, face = "bold", angle=0), plot.title = element_text(size=10, hjust=0.5), legend.position = "none", panel.spacing.x = unit(0.2, "mm"), panel.spacing.y = unit(0.75, "mm"))+
facet_grid(variable.new~Mid, scales = "free_x", switch="x")+
ggtitle('F3 generation - a subset of haplotypes across Cky-1')

patch<-plot + plot_annotation(tag_levels = 'A') + plot_layout(axis_titles = "collect_y")
patch
ggsave("Subset_Cky1_mhap_freqs_F3gen.tiff", patch, device=tiff, width=2250, height=2250, units="px", dpi=320)
```
<img width="2250" height="2250" alt="image" src="https://github.com/user-attachments/assets/52858ea8-7d3f-4bab-ac28-4a4f00664fe0" />


## For HCON_00158050:
```
paste FnL_SNPs_in_mhs_q125_HCON_00158050_SLPS_mhap_76SNPs_F3only.txt NSNPs_in_mhs_q125_HCON_00158050_SLPS_mhap_76SNPs_F3only.txt q125_HCON_00158050_SLPS_mhap_76SNPs_F3only.txt | cut -f1-5,25- > ForRplot_q125_HCON_00158050_SLPS_mhap_76SNPs_F3only.txt

sed 's/hcontortus_chr//g' ForRplot_q125_HCON_00158050_SLPS_mhap_76SNPs_F3only.txt | sed 's/_Celeg_TT_arrow_pilon//g' > new_ForRplot_q125_HCON_00158050_SLPS_mhap_76SNPs_F3only.txt

cat header_af new_ForRplot_q125_HCON_00158050_SLPS_mhap_76SNPs_F3only.txt > FINAL_ForRplot_q125_HCON_00158050_SLPS_mhap_76SNPs_F3only.txt

```

```

module load apps/R/4.4.1/gcc-8.5.0+openblas-0.3.28

R

library(stringr)
library(patchwork)
library(ggplot2)
library(dplyr)
library(reshape2)
library(viridis)

#[1] reshape2_1.4.4  dplyr_1.1.4     ggplot2_3.5.2   patchwork_1.3.1
#[5] stringr_1.5.1  viridis_0.6.5

df<-read.table("FINAL_ForRplot_q125_HCON_00158050_SLPS_mhap_76SNPs_F3only.txt", header=TRUE, sep="\t", na.strings = "NA")
summary(df)

dfmelt<-melt(df, id.vars=c("NSNPs","Start","End","Mid","Chr"))
head(dfmelt)

write.table(dfmelt, file="FINAL_ForRplot_q125_HCON_00158050_SLPS_mhap_76SNPs_F3only_MELTED.txt", row.names=FALSE, sep="\t", quote= FALSE)
```
```
for i in {6..14} ;
do
sed 's/AlleleFreq_//g' FINAL_ForRplot_q125_HCON_00158050_SLPS_mhap_76SNPs_F3only.txt | cut -f$i | sed 's/,/\t/g' > Col${i} ;
done

for i in Col* ;
do
awk '{print FILENAME "\t" NF}' $i >> ncols ;
done

# I then manually changed the headers to have the max number of columns each sample required

cut -f1-5 FINAL_ForRplot_q125_HCON_00158050_SLPS_mhap_76SNPs_F3only.txt > start
```
```
# first scan to find the max number of columns and fill the missing entries in the second round.
# file{,} is expanded by bash to file file, a neat trick not to repeat the filename (and eliminates possible typos).
# Note that it will produce empty files unless have file{,}.
for NUMBER in {6..14} ;
do
awk 'BEGIN   {OFS=FS="\t"} 
       NR==FNR {max=max<NF?NF:max; next} 
               {for(i=NF+1;i<=max;i++) $i="null"}1' Col$NUMBER{,} > Col${NUMBER}_new ;
done
```
PERFECT! :D 

Just need to change the header line
```
for NUMBER in {6..14} ;
do
head -n 1 Col${NUMBER}_new > header_${NUMBER} ;
done

# Rename the nulls
sed 's/null/F3_C_L1/g' header_6 > new_header_6
sed 's/null/F3_C_L2/g' header_7 > new_header_7
sed 's/null/F3_C_L3/g' header_8 > new_header_8
sed 's/null/F3_I_L1/g' header_9 > new_header_9
sed 's/null/F3_I_L2/g' header_10 > new_header_10
sed 's/null/F3_I_L3/g' header_11 > new_header_11
sed 's/null/F3_M_L1/g' header_12 > new_header_12
sed 's/null/F3_M_L2/g' header_13 > new_header_13
sed 's/null/F3_M_L3/g' header_14 > new_header_14

for i in {6..14} ;
do
tail -n+2 Col${i}_new | cat new_header_$i - > final_col${i} ;
done
```
Let's rejoin them all.... 

```
paste start final_col* > Updated_ForRplot_q125_HCON_00158050_SLPS_mhap_76SNPs_F3only.txt

```

```
module load apps/R/4.4.1/gcc-8.5.0+openblas-0.3.28

R

library(stringr)
library(patchwork)
library(ggplot2)
library(dplyr)
library(reshape2)
library(viridis)

#[1] reshape2_1.4.4  dplyr_1.1.4     ggplot2_3.5.2   patchwork_1.3.1
#[5] stringr_1.5.1  viridis_0.6.5

df<-read.table("Updated_ForRplot_q125_HCON_00158050_SLPS_mhap_76SNPs_F3only.txt", header=TRUE, sep="\t", na.strings = "null")
summary(df)

dfmelt<-melt(df, id.vars=c("NSNPs","Start","End","Mid","Chr"))
head(dfmelt)

write.table(dfmelt, file="Updated_ForRplot_q125_HCON_00158050_SLPS_mhap_76SNPs_F3only_MELTED.txt", row.names=FALSE, sep="\t", quote= FALSE)

```
Then, replace the : with a tab: 
```
sed 's/:/\t/g' Updated_ForRplot_q125_HCON_00158050_SLPS_mhap_76SNPs_F3only_MELTED.txt > New_ForRplot_q125_HCON_00158050_SLPS_mhap_76SNPs_F3only_MELTED.txt

# Manually updated the header line.

# Remove any lines with NA as these will now be shorter than the others and I cannot fix it for each haplotype - the data is simply missing. Also remove any haplotype for a sample with <0.01 frequency just to try to make the graph more plottable. 

awk 'BEGIN {OFS=FS="\t"} $7!="NA" {print $0}'  New_ForRplot_q125_HCON_00158050_SLPS_mhap_76SNPs_F3only_MELTED.txt | awk 'BEGIN {OFS=FS="\t"} $8>=0.01 {print $0}' > Finally_New_ForRplot_q125_HCON_00158050_SLPS_mhap_76SNPs_F3only_MELTED.txt

# Note that removing freq <0.01 takes it from 1111 to 698 lines.
```


```
module load apps/R/4.4.1/gcc-8.5.0+openblas-0.3.28

R

library(stringr)
library(patchwork)
library(ggplot2)
library(dplyr)
library(reshape2)
library(viridis)

#[1] reshape2_1.4.4  dplyr_1.1.4     ggplot2_3.5.2   patchwork_1.3.1
#[5] stringr_1.5.1  viridis_0.6.5

df<-read.table("Finally_New_ForRplot_q125_HCON_00158050_SLPS_mhap_76SNPs_F3only_MELTED.txt", header=TRUE, sep="\t", na.strings = "null")
summary(df)

# Then, replace the .2 to .xx with gsub():

df$variable.new<-gsub("\\...","", df$variable) # Replace everything that is .10 upwards
df$variable.new<-gsub("\\..","", df$variable.new) # Replace everything that is .2 to .9. Have to do this AFTER The one above or it doesn't work. 
df$variable.new<-as.factor(df$variable.new)
df$Haplotype<-as.factor(df$Haplotype)

plot<-ggplot(data=df, aes(x=Haplotype, y=Frequency))+
geom_col(aes(fill=variable.new))+
scale_fill_viridis_d()+
labs(y="Haplotype frequency", x="Haplotype")+
theme_light()+
theme(axis.title.x = element_text(size=10), axis.text.x = element_blank(), axis.title.y = element_text(size=10), axis.text.y = element_text(size=6), strip.text.x = element_text(size = 6, face = "bold", margin = margin(0.5,0,0.5,0, "mm")), strip.text.y = element_text(size = 6, face = "bold", angle=0), plot.title = element_text(size=10, hjust=0.5), legend.position = "none", panel.spacing.x = unit(0.2, "mm"), panel.spacing.y = unit(0.75, "mm"))+
facet_grid(variable.new~Mid, scales = "free_x", switch="x")+
ggtitle('F3 generation - a subset of haplotypes across HCON_00158050')

patch<-plot + plot_annotation(tag_levels = 'A') + plot_layout(axis_titles = "collect_y")
patch
ggsave("Subset_HCON_00158050_mhap_freqs_F3gen.tiff", patch, device=tiff, width=2250, height=2250, units="px", dpi=320)
```
<img width="2250" height="2250" alt="image" src="https://github.com/user-attachments/assets/fa8906f2-e026-4836-aff1-5338ab79d2e1" />


To plot them both together for the supplementary:
```
module load apps/R/4.4.1/gcc-8.5.0+openblas-0.3.28

R

library(stringr)
library(patchwork)
library(ggplot2)
library(dplyr)
library(reshape2)
library(viridis)

#[1] reshape2_1.4.4  dplyr_1.1.4     ggplot2_3.5.2   patchwork_1.3.1
#[5] stringr_1.5.1  viridis_0.6.5

df<-read.table("Finally_New_ForRplot_q125_HCON_00158050_SLPS_mhap_76SNPs_F3only_MELTED.txt", header=TRUE, sep="\t", na.strings = "null")
summary(df)

# Then, replace the .2 to .xx with gsub():

df$variable.new<-gsub("\\...","", df$variable) # Replace everything that is .10 upwards
df$variable.new<-gsub("\\..","", df$variable.new) # Replace everything that is .2 to .9. Have to do this AFTER The one above or it doesn't work. 
df$variable.new<-as.factor(df$variable.new)
df$Haplotype<-as.factor(df$Haplotype)

plot_HCON_00158050<-ggplot(data=df, aes(x=Haplotype, y=Frequency))+
geom_col(aes(fill=variable.new))+
scale_fill_viridis_d()+
labs(y="Haplotype frequency", x="Haplotype")+
theme_light()+
theme(axis.title.x = element_text(size=10), axis.text.x = element_blank(), axis.title.y = element_text(size=8), axis.text.y = element_text(size=6), strip.text.x = element_text(size = 6, face = "bold", margin = margin(0.5,0,0.5,0, "mm")), strip.text.y = element_text(size = 6, face = "bold", angle=0), plot.title = element_text(size=10, hjust=0.5), legend.position = "none", panel.spacing.x = unit(0.2, "mm"), panel.spacing.y = unit(0.75, "mm"))+
facet_grid(variable.new~Mid, scales = "free_x", switch="x")


df2<-read.table("Finally_New_ForRplot_q125_Cky1_SLPS_mhap_76SNPs_F3only_MELTED.txt", header=TRUE, sep="\t", na.strings = "null")
summary(df2)

# Then, replace the .2 to .xx with gsub():

df2$variable.new<-gsub("\\...","", df2$variable) # Replace everything that is .10 upwards
df2$variable.new<-gsub("\\..","", df2$variable.new) # Replace everything that is .2 to .9. Have to do this AFTER The one above or it doesn't work. 
df2$variable.new<-as.factor(df2$variable.new)
df2$Haplotype<-as.factor(df2$Haplotype)

plot_Cky1<-ggplot(data=df2, aes(x=Haplotype, y=Frequency))+
geom_col(aes(fill=variable.new))+
scale_fill_viridis_d()+
labs(y="Haplotype frequency", x="Haplotype")+
theme_light()+
theme(axis.title.x = element_text(size=10), axis.text.x = element_blank(), axis.title.y = element_text(size=8), axis.text.y = element_text(size=6), strip.text.x = element_text(size = 6, face = "bold", margin = margin(0.5,0,0.5,0, "mm")), strip.text.y = element_text(size = 6, face = "bold", angle=0), plot.title = element_text(size=10, hjust=0.5), legend.position = "none", panel.spacing.x = unit(0.2, "mm"), panel.spacing.y = unit(0.75, "mm"))+
facet_grid(variable.new~Mid, scales = "free_x", switch="x")


patch<- plot_Cky1 / plot_HCON_00158050 + plot_annotation(tag_levels = 'A')
patch
ggsave("Subset_Cky1_A_HCON_00158050_B_mhap_freqs_F3gen.tiff", patch, device=tiff, width=2250, height=2250, units="px", dpi=320)

```
<img width="2250" height="2250" alt="image" src="https://github.com/user-attachments/assets/cf557735-7286-419c-a589-9c9214fdce97" />
