To extract those SNPs which show extreme change and are considered significant following Holm correction in the AF time models.

```
# Total SNPs used in model = 938316
```
```
# MOX vs CTL:
awk 'BEGIN {OFS=FS="\t"} $6=="TRUE" {print $0}' p_values_MC_Holm.txt | wc -l > MC_sig_snps
# 9092 SNPs considered significant

awk 'BEGIN {OFS=FS="\t"} $2<-0.5 && $6=="TRUE" {print $0}' p_values_MC_Holm.txt > MC_neg_and_sig

chr5_33387923   -1.11803174057752       15.1853525763689        5.8979882519714e-10     9.22929609684171        TRUE
chr5_33387958   -1.11803174057752       15.1853525763689        5.8979882519714e-10     9.22929609684171        TRUE
chr5_35523631   -0.738944344455043      15.3918610461831        3.66604537553629e-10    9.43580216397198        TRUE
chr5_37381476   -0.559017006786758      285.997310438273        9.09436687479573e-281   280.041227529954        TRUE
chr5_39135563   -0.559017006786758      79.7908456683747        1.46297451044863e-74    73.834763240565 TRUE
chr5_39135577   -0.559017006786758      79.7908456683747        1.46297451044863e-74    73.834763240565 TRUE
chr5_39559227   -0.636095768839916      14.7608712776078        1.56738639085384e-09    8.80482392840052        TRUE
chr5_44017946   -0.558898646751354      15.9794503614772        9.47581754244903e-11    10.0233833102404        TRUE
```

```
# IVM vs CTL:
awk 'BEGIN {OFS=FS="\t"} $6=="TRUE" {print $0}' p_values_IC_Holm.txt | wc -l > IC_sig_snps
# 11871 SNPs considered significant

awk 'BEGIN {OFS=FS="\t"} $2<-0.5 && $6=="TRUE" {print $0}' p_values_IC_Holm.txt > IC_neg_and_sig

chr5_33856861   -0.568552049286365      16.4874028985853        2.52315326585629e-11    10.5980563679945        TRUE
chr5_44581709   -0.568552049286365      18.2160305273473        4.71327442785336e-13    12.326677272943 TRUE
```

```
# SEL vs CTL:
awk 'BEGIN {OFS=FS="\t"} $6=="TRUE" {print $0}' p_values_SEL_Holm.txt | wc -l > SEL_sig_snps
# 13818 SNPs considered significant

awk 'BEGIN {OFS=FS="\t"} $2<-0.5 && $6=="TRUE" {print $0}' p_values_SEL_Holm.txt > SEL_neg_and_sig

chr5_32689143   -0.565324549675811      14.8088865561541        1.39124002950598e-09    8.85659793507129        TRUE
chr5_39559441   -0.52286854067844       15.9470936076492        1.01203943862654e-10    9.99480256294667        TRUE
chr5_44581709   -0.565324549675811      16.8587238911175        1.24041818070356e-11    10.9064318769709        TRUE
```
