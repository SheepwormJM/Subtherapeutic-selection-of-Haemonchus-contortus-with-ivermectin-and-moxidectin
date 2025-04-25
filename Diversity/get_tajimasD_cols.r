get_tajimasD_cols <- function (my_data, half_window) {

library(dplyr)
  
# Let's make Chromosome a factor not a character. 
df$chrom<-as.factor(df$chrom)
summary(df)
# Perfect!

# And calculate the mid-point of the window for plotting:
df$mid<-(df$end - half_window)

# Select only those columns with tajimas D:
df_new<-select(df,chrom,mid,contains("tajimas_d"))

# Try removing windows with 'Na' in at least one pwc. The idea being (from previous work on this by me) that this would remove some of the lowest SNP windows. The data was indicating that having more SNPs reduced Fst, and was significant, but there was still a fair amount of variation, and hard to define an appropriate cut-off. However, if have only a few SNPs then likely to have some pwc with no SNPs and hence Na (from grenedalf wiki website: '--na-entry TEXT=nan Set the text to use in the output for n/a and NaN entries (e.g., resulting from positions with no counts, or windows with no variants). This is useful to match formatting expectations of downstream software. Output'.)

# Filter out any windows which have 'NaN' in at least one pwc:

df2<-na.omit(df_new)
summary(df2)
  write.table(df2, paste0(fileName,'.tajimas_d.txt'), row.names=FALSE, sep="\t", quote= FALSE)
  return(df2)

}


