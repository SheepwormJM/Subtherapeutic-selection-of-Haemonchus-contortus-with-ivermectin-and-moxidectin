_Selection-lines-only analysis with the new ss57 mpileup file made in Jan24_

Ok, let's assume that the ISE genome assembly is SENSITIVE in genotype.

Then, let's assume that the control samples are almost entirely (if not entirely) SENSITIVE in genotype.

Then, let's assume that the RESISTANT genotype will increase in allele frequency overtime, and that it is Non-reference. 

1. If we further assume that it will be a SNP, then let's try identifying sites where the REFERENCE allele steadily reduces from the F0->I1->I2->I3 and/or from F0->M1->M2->M3, but not in any control lines. 

2. Let's also try identifying sites where the reference allele is zero in any F3 selected sample, but >=20 in all F0 and F3 control samples (there are a lot of sites throughout ChrV which have a reference allele count of <20 in these samples, likely due to diversity to the reference genome sample). 

Note that SNPs which reduce only in early generations, but never reach zero will not be identified. This method may identify important SNPs, but is more robust to detect regions under selection and strength of selection than indivdual SNPs of importance.

# So, make new AF count file for all, using ss57:

```
module load grenedalf/0.3.0

grenedalf frequency --sync-path /home/jenni/data_folder/working_folder/selection_lines_poolseq/mpileup/Jan24_correct_order_syncs/SLps_only_co_Hcwbps18_q20Q30_noindels_ss57.java.sync --reference-genome-fasta-file /home/jenni/data_folder/wbps18/haemonchus_contortus.PRJEB506.WBPS18.genomic.fa --rename-samples-file /home/jenni/data_folder/working_folder/selection_lines_poolseq/pairwise_comparisons/LINES_IN_CORRECT_ORDER/grenedalf/Jan24_F0_vs_all_sel_lines_only/SL_ss57_sample_names.txt --write-sample-counts --separator-char tab --out-dir /home/jenni/data_folder/working_folder/selection_lines_poolseq/pairwise_comparisons/LINES_IN_CORRECT_ORDER/grenedalf/Jan24_F0_vs_all_sel_lines_only/ss57/AF --file-prefix 400_Gf_F0_vs_SL_only_Sel_lines_Hcwbps18_q20Q30_noindels_ss57_AF_sample_counts --log-file grenedalf_AF_all_ss57.log 2>error_log_all_57_af &
```
_Processed 89952741 genome positions of the input file, from 6 chromosomes, and thereof skipped 55456798 due to being invariant sites._

Code to determine reducing ref allele frequency: 
1. Sites where reference allele reduces from F0 --> F3

```
sed 's/\t/,/g' 400_Gf_F0_vs_SL_only_Sel_lines_Hcwbps18_q20Q30_noindels_ss57_AF_sample_countsfrequency.csv | sed 's/_Celeg_TT_arrow_pilon//g' | sed 's/hcontortus_//g' > 400_Gf_F0_vs_SL_only_Sel_lines_Hcwbps18_q20Q30_noindels_ss57_AF_sample_counts.frequency.csv

python3

import pandas as pd
import numpy as np

df=pd.read_csv("400_Gf_F0_vs_SL_only_Sel_lines_Hcwbps18_q20Q30_noindels_ss57_AF_sample_counts.frequency.csv") # Add ,header=None if no header. Although then found hard to manipulate later.
df.head

# Note that as I am putting in text, the np.nan now puts it in as nan, rather than when I put in a value from another column when it became NaN. At least, that is what I see. Maybe it is also that the redMOXL1 is in ''?

# First, identify all SNPs which have the reference allele count steadily reducing from F0>F1>F2>F3 by line and by treatment. Also find them for the control lines too - these SNPS, which reduce in the CTL lines can be subsequently filtered.

df['ref_reducesMOX_L1'] = np.where((df['F3_M_L1.REF_CNT'] < df['F2_M_L1.REF_CNT']) & (df['F2_M_L1.REF_CNT'] < df['F1_M_L1.REF_CNT']) & (df['F1_M_L1.REF_CNT'] < df['F0.REF_CNT']), 'redMOXL1', np.nan)
df['ref_reducesMOX_L2'] = np.where((df['F3_M_L2.REF_CNT'] < df['F2_M_L2.REF_CNT']) & (df['F2_M_L2.REF_CNT'] < df['F1_M_L2.REF_CNT']) & (df['F1_M_L2.REF_CNT'] < df['F0.REF_CNT']), 'redMOXL2', np.nan)
df['ref_reducesMOX_L3'] = np.where((df['F3_M_L3.REF_CNT'] < df['F2_M_L3.REF_CNT']) & (df['F2_M_L3.REF_CNT'] < df['F1_M_L3.REF_CNT']) & (df['F1_M_L3.REF_CNT'] < df['F0.REF_CNT']), 'redMOXL3', np.nan)

# Remember no F1_I_L1 so to compare the F2 to F0 instead
df['ref_reducesIVM_L1'] = np.where((df['F3_I_L1.REF_CNT'] < df['F2_I_L1.REF_CNT']) & (df['F2_I_L1.REF_CNT'] < df['F0.REF_CNT']), 'redIVML1', np.nan)
df['ref_reducesIVM_L2'] = np.where((df['F3_I_L2.REF_CNT'] < df['F2_I_L2.REF_CNT']) & (df['F2_I_L2.REF_CNT'] < df['F1_I_L2.REF_CNT']) & (df['F1_I_L2.REF_CNT'] < df['F0.REF_CNT']), 'redIVML2', np.nan)
df['ref_reducesIVM_L3'] = np.where((df['F3_I_L3.REF_CNT'] < df['F2_I_L3.REF_CNT']) & (df['F2_I_L3.REF_CNT'] < df['F1_I_L3.REF_CNT']) & (df['F1_I_L3.REF_CNT'] < df['F0.REF_CNT']), 'redIVML3', np.nan)

df['ref_reducesCTL_L1'] = np.where((df['F3_C_L1.REF_CNT'] < df['F2_C_L1.REF_CNT']) & (df['F2_C_L1.REF_CNT'] < df['F1_C_L1.REF_CNT']) & (df['F1_C_L1.REF_CNT'] < df['F0.REF_CNT']), 'redCTL1', np.nan)
df['ref_reducesCTL_L2'] = np.where((df['F3_C_L2.REF_CNT'] < df['F2_C_L2.REF_CNT']) & (df['F2_C_L2.REF_CNT'] < df['F1_C_L2.REF_CNT']) & (df['F1_C_L2.REF_CNT'] < df['F0.REF_CNT']), 'redCTL2', np.nan)
df['ref_reducesCTL_L3'] = np.where((df['F3_C_L3.REF_CNT'] < df['F2_C_L3.REF_CNT']) & (df['F2_C_L3.REF_CNT'] < df['F1_C_L3.REF_CNT']) & (df['F1_C_L3.REF_CNT'] < df['F0.REF_CNT']), 'redCTL3', np.nan)

# Assuming a single allele is responsible for MOX-R in all three lines (which might not be the case) find the sites which reduce the reference allele count from F0>F1>F2>F3 in all three lines:
df['ref_reduces_allMOX'] = np.where((df['ref_reducesMOX_L1'] != "nan" ) & (df['ref_reducesMOX_L2'] != "nan") & (df['ref_reducesMOX_L3'] !=  "nan"), 'redallMOX', np.nan)
# Assuming a single allele is responsible for IVM-R in all three lines (which might not be the case):
df['ref_reduces_allIVM'] = np.where((df['ref_reducesIVM_L1'] != "nan" ) & (df['ref_reducesIVM_L2'] != "nan") & (df['ref_reducesIVM_L3'] !=  "nan"), 'redallIVM', np.nan)
# Assuming a single allele is responsible for both IVM-R AND MOX-R in all three lines (which might not be the case) find those which overlap above (so all 6 lines):
df['ref_reduces_bothIVMandMOX'] = np.where((df['ref_reduces_allIVM'] != "nan") & (df['ref_reduces_allMOX'] != "nan"), 'redboth', np.nan)

# Assuming the same allele causing resistance does NOT, by chance, reduce between any single pair of ctl lines... filter those which have snps that also reduce in ANY control line:
df['ref_reducesctl'] = np.where((df['ref_reducesCTL_L1'] != "nan") | (df['ref_reducesCTL_L2'] != "nan") | (df['ref_reducesCTL_L3'] != "nan"), 'redinaCTL', np.nan)

df['ref_reduces_mox_not_ctl'] = np.where((df['ref_reducesctl'] == "nan") & (df['ref_reduces_allMOX'] == "redallMOX"), 'redMOXnotCTL', np.nan)
df['ref_reduces_ivm_not_ctl'] = np.where((df['ref_reducesctl'] == "nan") & (df['ref_reduces_allIVM'] == "redallIVM"), 'redIVMnoCTL', np.nan)

df['ref_reduces_bothsel_not_ctl'] = np.where((df['ref_reducesctl'] == "nan") & (df['ref_reduces_bothIVMandMOX'] == "redboth"), 'redSELnotCTL', np.nan)


#Hmm. Let's add a column where the F3 gen has a fixed zero reference in at least one line, while none of the F3_CTL, or the F0 does. This will NOT necessarily mean that another CTL population does not have zero ref. But it is a good skim. 

#Remember that the reference allele count could still be very low in the F0 or control population! It only has to be 1 to allow a site to be included
#**AMENDMENT**: I have since plotted it and found many sites with ref allele =0 in at least one sample,likely due to diversity relative to the reference sequence/poor alignment, and then updated, so that the F0/CTL population has to have a reference allele count of >= 20 based on the allele counts in the mox F3 populations. 

#Also note that what I still will NOT be showing are any snps where the ref allele reduces from F0 to F1 or F2 or F3, but does not steadily reduce, and is not zero in an F3 sample. But let's give this a go.

# Let's identify any where reference allele is 0 in at least one F3 generation, but >0 in F0, and all F3 control lines
# NOTE - you want to put the or etc in separate brackets to the ands, otherwise it will pull out ones where it is zero in an F3 sample, but also zero in an F0 sample.

df['MOXref_is_zero_F3_but_not_F0_CTL'] = np.where((df['F0.REF_CNT'] >= 20) & (df['F3_C_L1.REF_CNT'] >= 20) & (df['F3_C_L2.REF_CNT'] >= 20) & (df['F3_C_L3.REF_CNT'] >= 20) & ((df['F3_M_L1.REF_CNT'] == 0) | (df['F3_M_L2.REF_CNT'] == 0) | (df['F3_M_L3.REF_CNT'] == 0)), 'noREF_F3MOX', np.nan)
df['IVMref_is_zero_F3_but_not_F0_CTL'] = np.where((df['F0.REF_CNT'] >= 20) & (df['F3_C_L1.REF_CNT'] >= 20) & (df['F3_C_L2.REF_CNT'] >= 20) & (df['F3_C_L3.REF_CNT'] >= 20) & ((df['F3_I_L1.REF_CNT'] == 0) | (df['F3_I_L2.REF_CNT'] == 0) | (df['F3_I_L3.REF_CNT'] == 0)), 'noREF_F3IVM', np.nan)
df['bothref_is_zero_F3sel'] = np.where ((df['MOXref_is_zero_F3_but_not_F0_CTL'] != "nan") & (df['IVMref_is_zero_F3_but_not_F0_CTL'] != "nan"), 'bothrefzeroF3sel', np.nan)

# Save as new table:

df.to_csv('Mar25_marked_400_Gf_F0_vs_SL_only_Sel_lines_Hcwbps18_q20Q30_noindels_ss57_AF_sample_counts.frequency.csv', sep='\t',na_rep='NaN', index=False)
# index=FALSE removes the rownames

quit()
```



Ok, now can filter for those snps you want to plot. Then plot them as follows perhaps:
- subset chrV
- plot as geom_point with alpha=0.3
- plot F3 faceted by line, colour by red_mox_line (so binary - grey vs blue?)
- then plot over the top a subset of data (subset by red_both) - colour these all something else - orange
- then plot over the top a subset of data (subset by red_both_not_ctl) - could these all something else. - red

_Note that there are more data files below than were used for the final plot - some were not sensible to include_

```
grep -E 'redMOXnotCTL' Mar25_marked_400_Gf_F0_vs_SL_only_Sel_lines_Hcwbps18_q20Q30_noindels_ss57_AF_sample_counts.frequency.csv > redallMOX_Mar25_marked_400_Gf_F0_vs_SL_only_Sel_lines_Hcwbps18_q20Q30_noindels_ss57_AF_sample_counts.frequency

grep -E 'redIVMnoCTL' Mar25_marked_400_Gf_F0_vs_SL_only_Sel_lines_Hcwbps18_q20Q30_noindels_ss57_AF_sample_counts.frequency.csv > redallIVM_Mar25_marked_400_Gf_F0_vs_SL_only_Sel_lines_Hcwbps18_q20Q30_noindels_ss57_AF_sample_counts.frequency

grep -E 'redSELnotCTL' Mar25_marked_400_Gf_F0_vs_SL_only_Sel_lines_Hcwbps18_q20Q30_noindels_ss57_AF_sample_counts.frequency.csv > redboth_Mar25_marked_400_Gf_F0_vs_SL_only_Sel_lines_Hcwbps18_q20Q30_noindels_ss57_AF_sample_counts.frequency

grep -E 'noREF_F3MOX' Mar25_marked_400_Gf_F0_vs_SL_only_Sel_lines_Hcwbps18_q20Q30_noindels_ss57_AF_sample_counts.frequency.csv > F3_MOX_Mar25_refzero_marked_400_Gf_F0_vs_SL_only_Sel_lines_Hcwbps18_q20Q30_noindels_ss57_AF_sample_counts.frequency

grep -E 'noREF_F3IVM' Mar25_marked_400_Gf_F0_vs_SL_only_Sel_lines_Hcwbps18_q20Q30_noindels_ss57_AF_sample_counts.frequency.csv > F3_IVM_Mar25_refzero_marked_400_Gf_F0_vs_SL_only_Sel_lines_Hcwbps18_q20Q30_noindels_ss57_AF_sample_counts.frequency

grep -E 'bothrefzeroF3sel' Mar25_marked_400_Gf_F0_vs_SL_only_Sel_lines_Hcwbps18_q20Q30_noindels_ss57_AF_sample_counts.frequency.csv > F3_bothsel_Mar25_refzero_marked_400_Gf_F0_vs_SL_only_Sel_lines_Hcwbps18_q20Q30_noindels_ss57_AF_sample_counts.frequency

wc -l *frequency
       0 F3_bothsel_Mar25_refzero_marked_400_Gf_F0_vs_SL_only_Sel_lines_Hcwbps18_q20Q30_noindels_ss57_AF_sample_counts.frequency
      14 F3_IVM_Mar25_refzero_marked_400_Gf_F0_vs_SL_only_Sel_lines_Hcwbps18_q20Q30_noindels_ss57_AF_sample_counts.frequency
    1158 F3_MOX_Mar25_refzero_marked_400_Gf_F0_vs_SL_only_Sel_lines_Hcwbps18_q20Q30_noindels_ss57_AF_sample_counts.frequency
   11234 redallIVM_Mar25_marked_400_Gf_F0_vs_SL_only_Sel_lines_Hcwbps18_q20Q30_noindels_ss57_AF_sample_counts.frequency
   37437 redallMOX_Mar25_marked_400_Gf_F0_vs_SL_only_Sel_lines_Hcwbps18_q20Q30_noindels_ss57_AF_sample_counts.frequency
    4071 redboth_Mar25_marked_400_Gf_F0_vs_SL_only_Sel_lines_Hcwbps18_q20Q30_noindels_ss57_AF_sample_counts.frequency
  
for i in *frequency ; do grep -E 'chr5' ${i} > ${i}_chrV ; done

wc -l *chrV
       0 F3_bothsel_Mar25_refzero_marked_400_Gf_F0_vs_SL_only_Sel_lines_Hcwbps18_q20Q30_noindels_ss57_AF_sample_counts.frequency_chrV
      12 F3_IVM_Mar25_refzero_marked_400_Gf_F0_vs_SL_only_Sel_lines_Hcwbps18_q20Q30_noindels_ss57_AF_sample_counts.frequency_chrV
    1158 F3_MOX_Mar25_refzero_marked_400_Gf_F0_vs_SL_only_Sel_lines_Hcwbps18_q20Q30_noindels_ss57_AF_sample_counts.frequency_chrV
   10031 redallIVM_Mar25_marked_400_Gf_F0_vs_SL_only_Sel_lines_Hcwbps18_q20Q30_noindels_ss57_AF_sample_counts.frequency_chrV
   36228 redallMOX_Mar25_marked_400_Gf_F0_vs_SL_only_Sel_lines_Hcwbps18_q20Q30_noindels_ss57_AF_sample_counts.frequency_chrV
    4065 redboth_Mar25_marked_400_Gf_F0_vs_SL_only_Sel_lines_Hcwbps18_q20Q30_noindels_ss57_AF_sample_counts.frequency_chrV
```

Transfer to OneDrive those tables that have the refzero SNPs for future reference:
```
C:
cd "Users\jmi45g\OneDrive - University of Glasgow\BIOINFORMATICS\Roz's Hcon sequencing\Selection lines poolseq\AF"

scp -r jenni@scamper.mvls.gla.ac.uk:/home/jenni/data_folder/working_folder/selection_lines_poolseq/pairwise_comparisons/LINES_IN_CORRECT_ORDER/grenedalf/Jan24_F0_vs_all_sel_lines_only/ss57/AF/new_analysis_March25 ./SLPS_ss57_AF_redrefAC_Mar25
```

Transfer the header to the smaller ChrV files:
```
head -n1 Mar25_marked_400_Gf_F0_vs_SL_only_Sel_lines_Hcwbps18_q20Q30_noindels_ss57_AF_sample_counts.frequency.csv > header_AF_counts_marked

cat header_AF_counts_marked redallMOX_Mar25_marked_400_Gf_F0_vs_SL_only_Sel_lines_Hcwbps18_q20Q30_noindels_ss57_AF_sample_counts.frequency_chrV > tmp
mv tmp redallMOX_Mar25_marked_400_Gf_F0_vs_SL_only_Sel_lines_Hcwbps18_q20Q30_noindels_ss57_AF_sample_counts.frequency_chrV

cat header_AF_counts_marked redallIVM_Mar25_marked_400_Gf_F0_vs_SL_only_Sel_lines_Hcwbps18_q20Q30_noindels_ss57_AF_sample_counts.frequency_chrV > tmp
mv tmp redallIVM_Mar25_marked_400_Gf_F0_vs_SL_only_Sel_lines_Hcwbps18_q20Q30_noindels_ss57_AF_sample_counts.frequency_chrV

cat header_AF_counts_marked redboth_Mar25_marked_400_Gf_F0_vs_SL_only_Sel_lines_Hcwbps18_q20Q30_noindels_ss57_AF_sample_counts.frequency_chrV > tmp
mv tmp redboth_Mar25_marked_400_Gf_F0_vs_SL_only_Sel_lines_Hcwbps18_q20Q30_noindels_ss57_AF_sample_counts.frequency_chrV


cat header_AF_counts_marked F3_MOX_Mar25_refzero_marked_400_Gf_F0_vs_SL_only_Sel_lines_Hcwbps18_q20Q30_noindels_ss57_AF_sample_counts.frequency_chrV > tmp
mv tmp F3_MOX_Mar25_refzero_marked_400_Gf_F0_vs_SL_only_Sel_lines_Hcwbps18_q20Q30_noindels_ss57_AF_sample_counts.frequency_chrV

cat header_AF_counts_marked F3_IVM_Mar25_refzero_marked_400_Gf_F0_vs_SL_only_Sel_lines_Hcwbps18_q20Q30_noindels_ss57_AF_sample_counts.frequency_chrV > tmp
mv tmp F3_IVM_Mar25_refzero_marked_400_Gf_F0_vs_SL_only_Sel_lines_Hcwbps18_q20Q30_noindels_ss57_AF_sample_counts.frequency_chrV

cat header_AF_counts_marked F3_bothsel_Mar25_refzero_marked_400_Gf_F0_vs_SL_only_Sel_lines_Hcwbps18_q20Q30_noindels_ss57_AF_sample_counts.frequency_chrV > tmp
mv tmp F3_bothsel_Mar25_refzero_marked_400_Gf_F0_vs_SL_only_Sel_lines_Hcwbps18_q20Q30_noindels_ss57_AF_sample_counts.frequency_chrV

```

Get just ChrV from the main counts file: 
```
grep -E 'chr5' Mar25_marked_400_Gf_F0_vs_SL_only_Sel_lines_Hcwbps18_q20Q30_noindels_ss57_AF_sample_counts.frequency.csv > chr5_Mar25_marked_400_Gf_F0_vs_SL_only_Sel_lines_Hcwbps18_q20Q30_noindels_ss57_AF_sample_counts.frequency

wc -l chr5_Mar25_marked_400_Gf_F0_vs_SL_only_Sel_lines_Hcwbps18_q20Q30_noindels_ss57_AF_sample_counts.frequency
# 5302376 chr5_Mar25_marked_400_Gf_F0_vs_SL_only_Sel_lines_Hcwbps18_q20Q30_noindels_ss57_AF_sample_counts.frequency
```
So, there are clearly lots of SNPs which are not reducing in the reference allele frequency. However, this information is captured in another figure (AF model) and so I will not include these points in this, but simply show those that are, focusing the attention on regions with change on chrV. 

Select just those columns for CHROMOSOME V files that allow me to plot MOXIDECTIN selected F3 samples:

```
# Get POS, F3_M_L1.REF_CNT, F3_M_L2.REF_CNT, F3_M_L3.REF_CNT, ref_reduces_mox_not_ctl - JUST those that are reducing F0->F3 in all MOX, but not in a control line.
# Will colour these in BLUE
awk 'BEGIN {OFS="\t"} {print $2,$53,$55,$57,$72}' redallMOX_Mar25_marked_400_Gf_F0_vs_SL_only_Sel_lines_Hcwbps18_q20Q30_noindels_ss57_AF_sample_counts.frequency_chrV  > F3_MOX_redallMOX_Mar25_marked_400_Gf_F0_vs_SL_only_Sel_lines_Hcwbps18_q20Q30_noindels_ss57_AF_sample_counts.frequency_chrV

#PPOS     F3_M_L1.REF_CNT F3_M_L2.REF_CNT F3_M_L3.REF_CNT ref_reduces_mox_not_ctl
#271169  20      11      14      redMOXnotCTL
#1089253 32      39      46      redMOXnotCTL
#1599137 42      45      48      redMOXnotCTL
#1626077 27      25      25      redMOXnotCTL
# ...


# Get POS, F3_M_L1.REF_CNT, F3_M_L2.REF_CNT, F3_M_L3.REF_CNT, ref_reduces_bothsel_not_ctl - JUST those that are reducing F0->F3 in all MOX and all IVM, but not in a control line.
# Will colour these in orange
awk 'BEGIN {OFS="\t"} {print $2,$53,$55,$57,$74}' redboth_Mar25_marked_400_Gf_F0_vs_SL_only_Sel_lines_Hcwbps18_q20Q30_noindels_ss57_AF_sample_counts.frequency_chrV > F3_MOX_redboth_Mar25_marked_400_Gf_F0_vs_SL_only_Sel_lines_Hcwbps18_q20Q30_noindels_ss57_AF_sample_counts.frequency_chrV

#POS     F3_M_L1.REF_CNT F3_M_L2.REF_CNT F3_M_L3.REF_CNT ref_reduces_bothsel_not_ctl
#18804210        11      6       6       redSELnotCTL
#19093210        4       3       4       redSELnotCTL
#20698324        19      20      25      redSELnotCTL
#20704132        21      18      29      redSELnotCTL
# ...


# Get POS, F3_M_L1.REF_CNT, F3_M_L2.REF_CNT, F3_M_L3.REF_CNT, MOXref_is_zero_F3_but_not_F0_CTL
# Will colour these in red
awk 'BEGIN {OFS="\t"} {print $2,$53,$55,$57,$75}' F3_MOX_Mar25_refzero_marked_400_Gf_F0_vs_SL_only_Sel_lines_Hcwbps18_q20Q30_noindels_ss57_AF_sample_counts.frequency_chrV > JUST_F3_MOX_Mar25_refzero_marked_400_Gf_F0_vs_SL_only_Sel_lines_Hcwbps18_q20Q30_noindels_ss57_AF_sample_counts.frequency_chrV

#POS     F3_M_L1.REF_CNT F3_M_L2.REF_CNT F3_M_L3.REF_CNT MOXref_is_zero_F3_but_not_F0_CTL
#8969223 13      0       9       noREF_F3MOX
#20199739        33      10      0       noREF_F3MOX
#25925841        16      11      0       noREF_F3MOX
#34453508        17      0       4       noREF_F3MOX
# ...
```

And for IVM plots:
```
# Get POS, F3_I_L1.REF_CNT, F3_I_L2.REF_CNT, F3_I_L3.REF_CNT, ref_reduces_ivm_not_ctl - JUST those that are reducing F0->F3 in all IVM, but not in a control line.
# Will colour these in BLUE
awk 'BEGIN {OFS="\t"} {print $2,$35,$37,$39,$73}' redallIVM_Mar25_marked_400_Gf_F0_vs_SL_only_Sel_lines_Hcwbps18_q20Q30_noindels_ss57_AF_sample_counts.frequency_chrV  > F3_redallIVM_Mar25_marked_400_Gf_F0_vs_SL_only_Sel_lines_Hcwbps18_q20Q30_noindels_ss57_AF_sample_counts.frequency_chrV

#POS     F3_I_L1.REF_CNT F3_I_L2.REF_CNT F3_I_L3.REF_CNT ref_reduces_ivm_not_ctl
#20489   12      13      13      redIVMnoCTL
#28080   20      13      16      redIVMnoCTL
#28129   20      13      18      redIVMnoCTL
# ...

# Get POS, F3_I_L1.REF_CNT, F3_I_L2.REF_CNT, F3_I_L3.REF_CNT, ref_reduces_bothsel_not_ctl - JUST those that are reducing F0->F3 in all MOX and all IVM, but not in a control line.
# Will colour these in orange
awk 'BEGIN {OFS="\t"} {print $2,$35,$37,$39,$74}' redboth_Mar25_marked_400_Gf_F0_vs_SL_only_Sel_lines_Hcwbps18_q20Q30_noindels_ss57_AF_sample_counts.frequency_chrV > F3_IVM_redboth_Mar25_marked_400_Gf_F0_vs_SL_only_Sel_lines_Hcwbps18_q20Q30_noindels_ss57_AF_sample_counts.frequency_chrV

#POS     F3_I_L1.REF_CNT F3_I_L2.REF_CNT F3_I_L3.REF_CNT ref_reduces_bothsel_not_ctl
#18804210        16      5       9       redSELnotCTL
#19093210        7       8       6       redSELnotCTL
#20698324        30      15      28      redSELnotCTL
#20704132        24      14      24      redSELnotCTL
# ...


# Get POS, F3_I_L1.REF_CNT, F3_I_L2.REF_CNT, F3_I_L3.REF_CNT, IVMref_is_zero_F3_but_not_F0_CTL
awk 'BEGIN {OFS="\t"} {print $2,$35,$37,$39,$76}' F3_IVM_Mar25_refzero_marked_400_Gf_F0_vs_SL_only_Sel_lines_Hcwbps18_q20Q30_noindels_ss57_AF_sample_counts.frequency_chrV > JUST_F3_IVM_Mar25_refzero_marked_400_Gf_F0_vs_SL_only_Sel_lines_Hcwbps18_q20Q30_noindels_ss57_AF_sample_counts.frequency_chrV

#POS     F3_I_L1.REF_CNT F3_I_L2.REF_CNT F3_I_L3.REF_CNT IVMref_is_zero_F3_but_not_F0_CTL
#16506512        12      0       8       noREF_F3IVM
#19464379        29      0       18      noREF_F3IVM
#19771818        27      0       8       noREF_F3IVM
# ...

```
Checked had the right columns and numbers of SNPs for all files. :-) Note the huge difference in ref counts of zero between IVM and MOX lines! 

```
module load R/4.4.2
R

# Load libraries required
library(ggplot2)
library(reshape2)
library(patchwork)
# patchwork_1.3.0 reshape2_1.4.4  ggplot2_3.5.1 

# Load all tables

# For moxidectin plots:
redmox<-read.table("F3_MOX_redallMOX_Mar25_marked_400_Gf_F0_vs_SL_only_Sel_lines_Hcwbps18_q20Q30_noindels_ss57_AF_sample_counts.frequency_chrV", header=T, sep="\t", na.strings = "nan")

redselmox<-read.table("F3_MOX_redboth_Mar25_marked_400_Gf_F0_vs_SL_only_Sel_lines_Hcwbps18_q20Q30_noindels_ss57_AF_sample_counts.frequency_chrV", header=T, sep="\t", na.strings = "nan")

refzeromox<-read.table("JUST_F3_MOX_Mar25_refzero_marked_400_Gf_F0_vs_SL_only_Sel_lines_Hcwbps18_q20Q30_noindels_ss57_AF_sample_counts.frequency_chrV", header=T, sep="\t", na.strings = "nan")


# For ivermectin plots:
redivm<-read.table("F3_redallIVM_Mar25_marked_400_Gf_F0_vs_SL_only_Sel_lines_Hcwbps18_q20Q30_noindels_ss57_AF_sample_counts.frequency_chrV", header=T, sep="\t", na.strings = "nan")

redselivm<-read.table("F3_IVM_redboth_Mar25_marked_400_Gf_F0_vs_SL_only_Sel_lines_Hcwbps18_q20Q30_noindels_ss57_AF_sample_counts.frequency_chrV", header=T, sep="\t", na.strings = "nan")

refzeroivm<-read.table("JUST_F3_IVM_Mar25_refzero_marked_400_Gf_F0_vs_SL_only_Sel_lines_Hcwbps18_q20Q30_noindels_ss57_AF_sample_counts.frequency_chrV", header=T, sep="\t", na.strings = "nan")

gsub

# Melt all tables
redmox_melt<- melt(redmox, id.vars=c("POS","ref_reduces_mox_not_ctl"))
head(redmox_melt)

redsel_melt_mox <- melt(redselmox, id.vars=c("POS","ref_reduces_bothsel_not_ctl"))
head(redsel_melt_mox)

refzeroM_melt_mox <- melt(refzeromox, id.vars=c("POS","MOXref_is_zero_F3_but_not_F0_CTL"))
head(refzeroM_melt_mox)

summary(redmox_melt)
summary(redsel_melt_mox)
summary(refzeroM_melt_mox)


redivm_melt<- melt(redivm, id.vars=c("POS","ref_reduces_ivm_not_ctl"))
head(redivm_melt)

redsel_melt_ivm <- melt(redselivm, id.vars=c("POS","ref_reduces_bothsel_not_ctl"))
head(redsel_melt_ivm)

refzeroI_melt_ivm <- melt(refzeroivm, id.vars=c("POS","IVMref_is_zero_F3_but_not_F0_CTL"))
head(refzeroI_melt_ivm)

summary(redivm_melt)
summary(redsel_melt_ivm)
summary(refzeroI_melt_ivm)

# Rename the samples:
redivm_melt$variable.new<-gsub("F3_I_L","F3 IVM", redivm_melt$variable)
redivm_melt$variable.new<-gsub(".REF_CNT","", redivm_melt$variable.new)

redsel_melt_ivm$variable.new<-gsub("F3_I_L","F3 IVM", redsel_melt_ivm$variable)
redsel_melt_ivm$variable.new<-gsub(".REF_CNT","", redsel_melt_ivm$variable.new)

refzeroI_melt_ivm$variable.new<-gsub("F3_I_L","F3 IVM", refzeroI_melt_ivm$variable)
refzeroI_melt_ivm$variable.new<-gsub(".REF_CNT","", refzeroI_melt_ivm$variable.new)

redmox_melt$variable.new<-gsub("F3_M_L","F3 MOX", redmox_melt$variable)
redmox_melt$variable.new<-gsub(".REF_CNT","", redmox_melt$variable.new)

redsel_melt_mox$variable.new<-gsub("F3_M_L","F3 MOX", redsel_melt_mox$variable)
redsel_melt_mox$variable.new<-gsub(".REF_CNT","", redsel_melt_mox$variable.new)

refzeroM_melt_mox$variable.new<-gsub("F3_M_L","F3 MOX", refzeroM_melt_mox$variable)
refzeroM_melt_mox$variable.new<-gsub(".REF_CNT","", refzeroM_melt_mox$variable.new)

```


```
# Plot the MOX samples F3 gen:
# 
# Note, have 0,57 or 57,0 for the coord_cartesian to reverse the y-axis. 
xM<-ggplot(data=redmox_melt, aes(x=POS/1e6, y=value))+
coord_cartesian(ylim=c(0,57), xlim=c(0,49))+
scale_x_continuous(
    breaks = seq(0, 50, 5),
    minor_breaks = seq(0, 50, 1)
  )+
scale_y_continuous(
    breaks = seq(0, 57, 10),
    minor_breaks = seq(0, 57, 5)
  )+
facet_grid(variable.new~.)+
labs(y="Reference allele count", x="Genomic position (Mb)")+
theme_light()+
theme(axis.title.x = element_text(size=14), axis.text.x = element_text(size=10), axis.title.y = element_text(size=14), axis.text.y = element_text(size=10), strip.text.x = element_text(size = 12, face = "bold", margin = margin(0.5,0,0.5,0, "mm")), strip.text.y = element_text(size = 12, face = "bold", angle=0), plot.title = element_text(size=10, hjust=0.5), legend.position = "none", panel.spacing.x = unit(0.2, "mm"), panel.spacing.y = unit(0.75, "mm"))

yM<-xM+
geom_point(data=redmox_melt, colour="#619CFF", alpha=1)+ # the chrV blue 
geom_point(data=redsel_melt_mox, colour="orange", alpha=1)+
geom_point(data=refzeroM_melt_mox, colour="#D55E00", alpha=1) # a red colour from a colour blind palette


# Plot the IVM samples F3 gen:

xI<-ggplot(data=redivm_melt, aes(x=POS/1e6, y=value))+
coord_cartesian(ylim=c(0,57), xlim=c(0,49))+
scale_x_continuous(
    breaks = seq(0, 50, 5),
    minor_breaks = seq(0, 50, 1)
  )+
scale_y_continuous(
    breaks = seq(0, 57, 10),
    minor_breaks = seq(0, 57, 5)
  )+
facet_grid(variable.new~.)+
labs(y="Reference allele count", x="Genomic position (Mb)")+
theme_light()+
theme(axis.title.x = element_text(size=14), axis.text.x = element_text(size=10), axis.title.y = element_text(size=14), axis.text.y = element_text(size=10), strip.text.x = element_text(size = 12, face = "bold", margin = margin(0.5,0,0.5,0, "mm")), strip.text.y = element_text(size = 12, face = "bold", angle=0), plot.title = element_text(size=10, hjust=0.5), legend.position = "none", panel.spacing.x = unit(0.2, "mm"), panel.spacing.y = unit(0.75, "mm"))

yI<-xI+
geom_point(data=redivm_melt, colour="#619CFF", alpha=1)+ # the chrV blue colour
geom_point(data=redsel_melt_ivm, colour="orange", alpha=1)+
geom_point(data=refzeroI_melt_ivm, colour="#D55E00", alpha=1) # a red colour from a colour blind palette

patchwork3 <- yI / yM
plotme2<- patchwork3 + plot_annotation(tag_levels = 'A') + plot_layout(axis_titles = "collect")

ggsave("ChrV_400_Gf_SLps_Hcwbps18_q20Q30_ss57_sellinesonly.AF_ref_count_MARCH2025_normY.tiff", plotme2, device=tiff, width=2250, height=3550, units="px", dpi=320)
```

_Chr V reference allele count, plotted March 2025_
![image](https://github.com/user-attachments/assets/940bae5a-abb0-4f97-aead-eeb3bcf86e36)
