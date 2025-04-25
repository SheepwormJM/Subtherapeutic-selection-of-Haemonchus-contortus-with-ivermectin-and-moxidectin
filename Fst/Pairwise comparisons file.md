Pairwise comparisons file:
```
F0      F1_C_L1
F0      F1_C_L2
F0      F1_C_L3
F0      F1_I_L2
F0      F1_I_L3
F0      F1_M_L1
F0      F1_M_L2
F0      F1_M_L3
F0      F2_C_L1
F0      F2_C_L2
F0      F2_C_L3
F0      F2_I_L1
F0      F2_I_L2
F0      F2_I_L3
F0      F2_M_L1
F0      F2_M_L2
F0      F2_M_L3
F0      F3_C_L1
F0      F3_C_L2
F0      F3_C_L3
F0      F3_I_L1
F0      F3_I_L2
F0      F3_I_L3
F0      F3_M_L1
F0      F3_M_L2
F0      F3_M_L3
F3_M_L1 F3_M_L2
F3_M_L1 F3_M_L3
F3_M_L2 F3_M_L3
F3_I_L1 F3_I_L2
F3_I_L1 F3_I_L3
F3_I_L2 F3_I_L3
F3_I_L1 F3_M_L1
F3_I_L2 F3_M_L2
F3_I_L3 F3_M_L3
```
Note, when copy and paste, get multiple spaces, not a tab (check with cat -T). To replace multiple spaces with a single space and then all single spaces with a tab do the following: 

```
sed 's/  */ /g' SL_pairwise_comparisons.txt | sed 's/ /\t/g' > tmp
```
