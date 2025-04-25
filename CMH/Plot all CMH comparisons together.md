Produce a plot with all the CMH tests overlaid - the top 1% SNPs by pvalue rank, further filtered to retain only those with a reduction in the reference allele in the F3 selected sample relative to the other (i.e. a log odds ratio <1).

1. Get all files in a new folder, and rename the sym links:

```
# CTL vs IVM:
ln -s /home/jenni/data_folder/cmh_oct24_slps/ctl_ivm_F3/SNPS_SUM_500KB.txt_tmp2_sorted ./CTL_vs_IVM_SNPS_SUM_500KB.txt_tmp2_sorted
ln -s /home/jenni/data_folder/cmh_oct24_slps/ctl_ivm_F3/uniq_count_SNPS_500KB_lorless1_top1pc.txt_tmp_sorted ./CTL_vs_IVM_uniq_count_SNPS_500KB_lorless1_top1pc.txt_tmp_sorted

# CTL vs MOX:
ln -s /home/jenni/data_folder/cmh_oct24_slps/ctl_mox_F3/SNPS_SUM_500KB.txt_tmp2_sorted ./CTL_vs_MOX_SNPS_SUM_500KB.txt_tmp2_sorted
ln -s /home/jenni/data_folder/cmh_oct24_slps/ctl_mox_F3/uniq_count_SNPS_500KB_lorless1_top1pc.txt_tmp_sorted ./CTL_vs_MOX_uniq_count_SNPS_500KB_lorless1_top1pc.txt_tmp_sorted


# random CTL vs IVM or MOX:
ln -s /home/jenni/data_folder/cmh_oct24_slps/ctl_ivm_mox/SNPS_SUM_500KB.txt_tmp2_sorted ./ranCTL_vs_IVM_or_MOX_SNPS_SUM_500KB.txt_tmp2_sorted
ln -s /home/jenni/data_folder/cmh_oct24_slps/ctl_ivm_mox/uniq_count_SNPS_500KB_lorless1_top1pc.txt_tmp_sorted ./ranCTL_vs_IVM_or_MOX_uniq_count_SNPS_500KB_lorless1_top1pc.txt_tmp_sorted


# F2 vs F3 for IVM or MOX:
ln -s /home/jenni/data_folder/cmh_oct24_slps/F2_ivm_mox/SNPS_SUM_500KB.txt_tmp2_sorted ./F2_vs_F3_SNPS_SUM_500KB.txt_tmp2_sorted
ln -s /home/jenni/data_folder/cmh_oct24_slps/F2_ivm_mox/uniq_count_SNPS_500KB_lorless1_top1pc.txt_tmp_sorted ./F2_vs_F3_uniq_count_SNPS_500KB_lorless1_top1pc.txt_tmp_sorted

# F2 vs F3 IVM only:
ln -s /home/jenni/data_folder/cmh_oct24_slps/F2_IVM/SNPS_SUM_500KB.txt_tmp2_sorted ./F2_IVM_SNPS_SUM_500KB.txt_tmp2_sorted
ln -s /home/jenni/data_folder/cmh_oct24_slps/F2_IVM/uniq_count_SNPS_500KB_lorless1_top1pc.txt_tmp_sorted ./F2_IVM_uniq_count_SNPS_500KB_lorless1_top1pc.txt_tmp_sorted

# F2 vs F3 MOX only:
ln -s /home/jenni/data_folder/cmh_oct24_slps/F2_MOX/SNPS_SUM_500KB.txt_tmp2_sorted ./F2_MOX_SNPS_SUM_500KB.txt_tmp2_sorted
ln -s /home/jenni/data_folder/cmh_oct24_slps/F2_MOX/uniq_count_SNPS_500KB_lorless1_top1pc.txt_tmp_sorted ./F2_MOX_uniq_count_SNPS_500KB_lorless1_top1pc.txt_tmp_sorted

```
2. In R v4.4.2:
```
library(dplyr)
library(ggplot2)
library(patchwork)
library(viridis)

[1] viridis_0.6.5     viridisLite_0.4.2 patchwork_1.3.0   ggplot2_3.5.1    
[5] dplyr_1.1.4 

df<-read.table("CTL_vs_IVM_SNPS_SUM_500KB.txt_tmp2_sorted", header=FALSE, sep="\t")
df2<-read.table("CTL_vs_IVM_uniq_count_SNPS_500KB_lorless1_top1pc.txt_tmp_sorted", header=FALSE, sep="\t")

df3<-read.table("CTL_vs_MOX_SNPS_SUM_500KB.txt_tmp2_sorted", header=FALSE, sep="\t")
df4<-read.table("CTL_vs_MOX_uniq_count_SNPS_500KB_lorless1_top1pc.txt_tmp_sorted", header=FALSE, sep="\t")

df5<-read.table("ranCTL_vs_IVM_or_MOX_SNPS_SUM_500KB.txt_tmp2_sorted", header=FALSE, sep="\t")
df6<-read.table("ranCTL_vs_IVM_or_MOX_uniq_count_SNPS_500KB_lorless1_top1pc.txt_tmp_sorted", header=FALSE, sep="\t")

df7<-read.table("F2_vs_F3_SNPS_SUM_500KB.txt_tmp2_sorted", header=FALSE, sep="\t")
df8<-read.table("F2_vs_F3_uniq_count_SNPS_500KB_lorless1_top1pc.txt_tmp_sorted", header=FALSE, sep="\t")

df9<-read.table("F2_IVM_SNPS_SUM_500KB.txt_tmp2_sorted", header=FALSE, sep="\t")
df10<-read.table("F2_IVM_uniq_count_SNPS_500KB_lorless1_top1pc.txt_tmp_sorted", header=FALSE, sep="\t")

df11<-read.table("F2_MOX_SNPS_SUM_500KB.txt_tmp2_sorted", header=FALSE, sep="\t")
df12<-read.table("F2_MOX_uniq_count_SNPS_500KB_lorless1_top1pc.txt_tmp_sorted", header=FALSE, sep="\t")


# Need to make the chromosome number in the uniq_count file as a character - no ChrX SNPs unlike in the SNPS file - and inner_join() needs them to be the same type to join.
df2$V5 <- as.character(df2$V5)
df4$V5 <- as.character(df4$V5)
df6$V5 <- as.character(df6$V5)
df8$V5 <- as.character(df8$V5)
df10$V5 <- as.character(df10$V5)
df12$V5 <- as.character(df12$V5)

# Join tables, keeping only rows present in both tables (this will account for removal of windows with no SNPs
df13<-inner_join(df, df2, by=c("V3"="V5", "V4"="V2", "V5"="V3"))
df14<-inner_join(df3, df4, by=c("V3"="V5", "V4"="V2", "V5"="V3"))
df15<-inner_join(df5, df6, by=c("V3"="V5", "V4"="V2", "V5"="V3"))
df16<-inner_join(df7, df8, by=c("V3"="V5", "V4"="V2", "V5"="V3"))
df17<-inner_join(df9, df10, by=c("V3"="V5", "V4"="V2", "V5"="V3"))
df18<-inner_join(df11, df12, by=c("V3"="V5", "V4"="V2", "V5"="V3"))

# Calculate the proportion of SNPs in each window that are in the top 1pc by rank, and have log odds ratio < 1:
df13$prop<-(df13$V1.y/df13$V1.x)
df14$prop<-(df14$V1.y/df14$V1.x)
df15$prop<-(df15$V1.y/df15$V1.x)
df16$prop<-(df16$V1.y/df16$V1.x)
df17$prop<-(df17$V1.y/df17$V1.x)
df18$prop<-(df18$V1.y/df18$V1.x)

summary(df13)
summary(df14)
summary(df15)
summary(df16)
summary(df17)
summary(df18)

# Then, to allow the data to be easily plotted on the same plot:

df13$test <- paste0("CTL:IVM", df13$test)
df14$test <- paste0("CTL:MOX", df14$test)
df15$test <- paste0("ranCTL:SEL", df15$test)
df16$test <- paste0("F2:F3_SEL", df16$test)
df17$test <- paste0("F2:F3_IVM", df17$test)
df18$test <- paste0("F2:F3_MOX", df18$test)

# Stack the tables one on top of the other (i.e. concatenate without the headers)
df19<-bind_rows(df13,df14,df15,df16,df17,df18)

# Chromosome column is now V3 - I want to rename the chromsomes to be in roman numerals
df19$V3<-gsub("1","I", df19$V3)
df19$V3<-gsub("2","II", df19$V3)
df19$V3<-gsub("3","III", df19$V3)
df19$V3<-gsub("4","IV", df19$V3)
df19$V3<-gsub("5","V", df19$V3)


# Plot the data:

rect<-data.frame(xmin = c(NA, NA, NA, NA, 30, NA), xmax = c(NA,NA,NA,NA,45,NA), 
                 ymin = c(NA, NA, NA, NA, 0, NA), ymax = c(NA, NA, NA, NA, 0.15, NA), 
                 alpha = c(NA, NA, NA, NA, 0.2, NA),
                 fill = c(NA, NA, NA, NA, "red", NA))

x<-ggplot(data=df19, aes(x=((V4 + 249999)/1e6), y=(prop)))+
geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, alpha=alpha, fill=fill),
            data = transform(rect, V3 = as.character("V")),
            inherit.aes = FALSE)+
  geom_line(aes(colour=test), linewidth=0.8)+
scale_colour_viridis_d()+
facet_grid(.~V3, scales = "free_x", switch="x")+
labs(y="Proportion of SNPs per 500 Kb window that are in top 1%\nby rank of CMH p-values, with log-odds ratio < 1", x= "Genomic position (Mb)")+
   scale_x_continuous(
 breaks = seq(0, 55, 10),
         minor_breaks = seq(5, 50, 5)
     )+
theme_light()+
theme(legend.position = "none")+
guides(colour = guide_legend(override.aes = list(alpha = 1, size=1), title = "CMH Test"), alpha = "none", fill = "none")

#theme(axis.title.x = element_text(size=11), axis.text.x = element_text(size=10), axis.title.y = element_text(size=12), axis.text.y = element_text(size=10), strip.text.x = element_text(size = 12, face = "bold", margin = margin(0.5,0,0.5,0, "mm")), strip.text.y = element_text(size = 12, face = "bold", angle=0), plot.title = element_text(size=10, hjust=0.5), legend.position = "none", panel.spacing.x = unit(0.2, "mm"), panel.spacing.y = unit(0.75, "mm"))+


y<-ggplot(data=subset(df19, df19$V3 == "V"), aes(x=((V4 + 249999)/1e6), y=(prop)))+
geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, alpha=alpha, fill=fill),
            data = transform(rect, V3 = as.character("V")),
            inherit.aes = FALSE)+
 geom_line(aes(colour=test), linewidth=0.8)+
scale_colour_viridis_d()+
facet_grid(.~V3, scales = "free_x", switch="x")+
labs(y="Proportion of SNPs per 500 Kb window that are in top 1%\nby rank of CMH p-values, with log-odds ratio < 1", x= "Genomic position (Mb)")+
   scale_x_continuous(
 breaks = seq(0, 55, 5),
         minor_breaks = seq(1, 54, 1)
     )+
theme_light()+
theme(legend.position = "bottom")+
guides(colour = guide_legend(override.aes = list(alpha = 1, linewidth=3), title = "CMH Test"), alpha = "none", fill = "none")

#theme(axis.title.x = element_text(size=11), axis.text.x = element_text(size=10), axis.title.y = element_text(size=12), axis.text.y = element_text(size=10), strip.text.x = element_text(size = 12, face = "bold", margin = margin(0.5,0,0.5,0, "mm")), strip.text.y = element_text(size = 12, face = "bold", angle=0), plot.title = element_text(size=10, hjust=0.5), legend.position = "bottom", panel.spacing.x = unit(0.2, "mm"), panel.spacing.y = unit(0.75, "mm"))+
guides(colour = guide_legend(override.aes = list(alpha = 1, linewidth=1.5), title = "CMH Test"), alpha = "none", fill = "none")


patch <- x/y
patch2 <- patch + plot_annotation(tag_levels = 'A') + plot_layout(axis_titles = 'collect', axes = 'collect') & theme(plot.tag = element_text(size = 20))
patch2

ggsave("HcSLPS_SNP_PROP_CMH_Chr_500kb_point_patch_allcomparisons_romnumerals.tiff", device=tiff, width=2250, height=2250, units="px", dpi=320)
```
_Below: Updated image with roman numerals (without specifying all of theme text sizes etc, and with Genomic position as x-axis title): _
![image](https://github.com/user-attachments/assets/b1272b64-cf92-4677-9f67-ef4dc3843f82)

_Below: All CMH tests overlaid, original image_
![image](https://github.com/user-attachments/assets/121844fa-657d-4876-af92-979f06213bac)


_With simpler y-axis:_
```
library(dplyr)
library(ggplot2)
library(patchwork)
library(viridis)

[1] viridis_0.6.5     viridisLite_0.4.2 patchwork_1.3.0   ggplot2_3.5.1    
[5] dplyr_1.1.4 

df<-read.table("CTL_vs_IVM_SNPS_SUM_500KB.txt_tmp2_sorted", header=FALSE, sep="\t")
df2<-read.table("CTL_vs_IVM_uniq_count_SNPS_500KB_lorless1_top1pc.txt_tmp_sorted", header=FALSE, sep="\t")

df3<-read.table("CTL_vs_MOX_SNPS_SUM_500KB.txt_tmp2_sorted", header=FALSE, sep="\t")
df4<-read.table("CTL_vs_MOX_uniq_count_SNPS_500KB_lorless1_top1pc.txt_tmp_sorted", header=FALSE, sep="\t")

df5<-read.table("ranCTL_vs_IVM_or_MOX_SNPS_SUM_500KB.txt_tmp2_sorted", header=FALSE, sep="\t")
df6<-read.table("ranCTL_vs_IVM_or_MOX_uniq_count_SNPS_500KB_lorless1_top1pc.txt_tmp_sorted", header=FALSE, sep="\t")

df7<-read.table("F2_vs_F3_SNPS_SUM_500KB.txt_tmp2_sorted", header=FALSE, sep="\t")
df8<-read.table("F2_vs_F3_uniq_count_SNPS_500KB_lorless1_top1pc.txt_tmp_sorted", header=FALSE, sep="\t")

df9<-read.table("F2_IVM_SNPS_SUM_500KB.txt_tmp2_sorted", header=FALSE, sep="\t")
df10<-read.table("F2_IVM_uniq_count_SNPS_500KB_lorless1_top1pc.txt_tmp_sorted", header=FALSE, sep="\t")

df11<-read.table("F2_MOX_SNPS_SUM_500KB.txt_tmp2_sorted", header=FALSE, sep="\t")
df12<-read.table("F2_MOX_uniq_count_SNPS_500KB_lorless1_top1pc.txt_tmp_sorted", header=FALSE, sep="\t")


# Need to make the chromosome number in the uniq_count file as a character - no ChrX SNPs unlike in the SNPS file - and inner_join() needs them to be the same type to join.
df2$V5 <- as.character(df2$V5)
df4$V5 <- as.character(df4$V5)
df6$V5 <- as.character(df6$V5)
df8$V5 <- as.character(df8$V5)
df10$V5 <- as.character(df10$V5)
df12$V5 <- as.character(df12$V5)

# Join tables, keeping only rows present in both tables (this will account for removal of windows with no SNPs
df13<-inner_join(df, df2, by=c("V3"="V5", "V4"="V2", "V5"="V3"))
df14<-inner_join(df3, df4, by=c("V3"="V5", "V4"="V2", "V5"="V3"))
df15<-inner_join(df5, df6, by=c("V3"="V5", "V4"="V2", "V5"="V3"))
df16<-inner_join(df7, df8, by=c("V3"="V5", "V4"="V2", "V5"="V3"))
df17<-inner_join(df9, df10, by=c("V3"="V5", "V4"="V2", "V5"="V3"))
df18<-inner_join(df11, df12, by=c("V3"="V5", "V4"="V2", "V5"="V3"))

# Calculate the proportion of SNPs in each window that are in the top 1pc by rank, and have log odds ratio < 1:
df13$prop<-(df13$V1.y/df13$V1.x)
df14$prop<-(df14$V1.y/df14$V1.x)
df15$prop<-(df15$V1.y/df15$V1.x)
df16$prop<-(df16$V1.y/df16$V1.x)
df17$prop<-(df17$V1.y/df17$V1.x)
df18$prop<-(df18$V1.y/df18$V1.x)

summary(df13)
summary(df14)
summary(df15)
summary(df16)
summary(df17)
summary(df18)

# Then, to allow the data to be easily plotted on the same plot:

df13$test <- paste0("CTL:IVM", df13$test)
df14$test <- paste0("CTL:MOX", df14$test)
df15$test <- paste0("ranCTL:SEL", df15$test)
df16$test <- paste0("F2:F3_SEL", df16$test)
df17$test <- paste0("F2:F3_IVM", df17$test)
df18$test <- paste0("F2:F3_MOX", df18$test)

# Stack the tables one on top of the other (i.e. concatenate without the headers)
df19<-bind_rows(df13,df14,df15,df16,df17,df18)

# Chromosome column is now V3 - I want to rename the chromsomes to be in roman numerals
df19$V3<-gsub("1","I", df19$V3)
df19$V3<-gsub("2","II", df19$V3)
df19$V3<-gsub("3","III", df19$V3)
df19$V3<-gsub("4","IV", df19$V3)
df19$V3<-gsub("5","V", df19$V3)


# Plot the data:

rect<-data.frame(xmin = c(NA, NA, NA, NA, 30, NA), xmax = c(NA,NA,NA,NA,45,NA), 
                 ymin = c(NA, NA, NA, NA, 0, NA), ymax = c(NA, NA, NA, NA, 0.15, NA), 
                 alpha = c(NA, NA, NA, NA, 0.2, NA),
                 fill = c(NA, NA, NA, NA, "red", NA))

x<-ggplot(data=df19, aes(x=((V4 + 249999)/1e6), y=(prop)))+
geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, alpha=alpha, fill=fill),
            data = transform(rect, V3 = as.character("V")),
            inherit.aes = FALSE)+
  geom_line(aes(colour=test), linewidth=0.8)+
scale_colour_viridis_d()+
facet_grid(.~V3, scales = "free_x", switch="x")+
labs(y="SNP Proportion", x= "Genomic position (Mb)")+
   scale_x_continuous(
 breaks = seq(0, 55, 10),
         minor_breaks = seq(5, 50, 5)
     )+
theme_light()+
theme(legend.position = "none")+
guides(colour = guide_legend(override.aes = list(alpha = 1, size=1), title = "CMH Test"), alpha = "none", fill = "none")

#theme(axis.title.x = element_text(size=11), axis.text.x = element_text(size=10), axis.title.y = element_text(size=12), axis.text.y = element_text(size=10), strip.text.x = element_text(size = 12, face = "bold", margin = margin(0.5,0,0.5,0, "mm")), strip.text.y = element_text(size = 12, face = "bold", angle=0), plot.title = element_text(size=10, hjust=0.5), legend.position = "none", panel.spacing.x = unit(0.2, "mm"), panel.spacing.y = unit(0.75, "mm"))+


y<-ggplot(data=subset(df19, df19$V3 == "V"), aes(x=((V4 + 249999)/1e6), y=(prop)))+
geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, alpha=alpha, fill=fill),
            data = transform(rect, V3 = as.character("V")),
            inherit.aes = FALSE)+
 geom_line(aes(colour=test), linewidth=0.8)+
scale_colour_viridis_d()+
facet_grid(.~V3, scales = "free_x", switch="x")+
labs(y="SNP Proportion", x= "Genomic position (Mb)")+
   scale_x_continuous(
 breaks = seq(0, 55, 5),
         minor_breaks = seq(1, 54, 1)
     )+
theme_light()+
theme(legend.position = "bottom")+
guides(colour = guide_legend(override.aes = list(alpha = 1, linewidth=3), title = "CMH Test"), alpha = "none", fill = "none")

#theme(axis.title.x = element_text(size=11), axis.text.x = element_text(size=10), axis.title.y = element_text(size=12), axis.text.y = element_text(size=10), strip.text.x = element_text(size = 12, face = "bold", margin = margin(0.5,0,0.5,0, "mm")), strip.text.y = element_text(size = 12, face = "bold", angle=0), plot.title = element_text(size=10, hjust=0.5), legend.position = "bottom", panel.spacing.x = unit(0.2, "mm"), panel.spacing.y = unit(0.75, "mm"))+
guides(colour = guide_legend(override.aes = list(alpha = 1, linewidth=1.5), title = "CMH Test"), alpha = "none", fill = "none")


patch <- x/y
patch2 <- patch + plot_annotation(tag_levels = 'A') + plot_layout(axis_titles = 'collect', axes = 'collect') & theme(plot.tag = element_text(size = 20))
patch2

ggsave("HcSLPS_SNP_PROP_CMH_Chr_500kb_point_patch_allcomparisons_romnumerals_simplerYaxis.tiff", device=tiff, width=2250, height=2250, units="px", dpi=320)
```
![image](https://github.com/user-attachments/assets/2de0ddf1-f7c6-4fb6-b889-6540a4cab51e)
