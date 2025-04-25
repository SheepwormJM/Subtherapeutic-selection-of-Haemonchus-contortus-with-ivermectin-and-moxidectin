#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)==1) {
  # default output file
  args[2] = "out.txt"
}

input_path = args[1]
input_file = args[2]
output_file = args[3]
my_date = args[4]


# Trying on scamper - note that cannot install tidyverse, so installing the various used modules within it

#module load R/4.4.2 

#then proceeding in R
library(dplyr)
library(ggplot2)
library(tidyr)
library(lme4)
library(stringr)
library(purrr)
library(broom.mixed)

# sessionInfo()
# [1] patchwork_1.3.0     broom.mixed_0.2.9.6 purrr_1.0.2        
# [4] stringr_1.5.1       tidyr_1.3.1         ggplot2_3.5.1      
# [7] dplyr_1.1.4         lme4_1.1-36         Matrix_1.7-1       

# Load data - note this will already have replicates as L1 -> L9 and contains a single F0.REF_CNT, F0.ALT_CNT and F0.DEPTH column in addition to F1 to F3 sample columns
all_data <- read.delim(paste(input_path, input_file, sep=""), header=TRUE)

head(all_data)

#try with F0 but downweight
#make some new column names
f0_colnames <- paste0(rep("F0",9),"_",rep(c("C","I","M"),3),"_",rep(c("L1","L2","L3","L4","L5","L6","L7","L8","L9"),each=1))
f0_ref_colnames <- paste0(f0_colnames,".REF_CNT")
f0_alt_colnames <- paste0(f0_colnames,".ALT_CNT")
f0_depth_colnames <- paste0(f0_colnames,".DEPTH")

#make 'copied' data for F0 for each of the above (C,I,M, L1-L9):
for (i in f0_ref_colnames)  {
    all_data[[i]] <- all_data$F0.REF_CNT
    print(ncol(all_data))
}

for (i in f0_alt_colnames)  {
    all_data[[i]] <- all_data$F0.ALT_CNT
    print(ncol(all_data))
}

for (i in f0_depth_colnames) { 
    all_data[[i]] <- all_data$F0.DEPTH
}

#remove original F0 columns
all_data <- all_data %>% select(!c("F0.REF_CNT","F0.ALT_CNT","F0.DEPTH"))

head(all_data)

#split sample names into timepoint, #replicate and treatment columns
#if you change anything above, make sure you check that the 'gather' command is specifying the right columns
reshaped_data <- all_data %>% gather(sample,count,F1_C_L1.REF_CNT:F0_M_L9.DEPTH) %>% separate(sample,c("sample","allele"),"[.]")  %>% separate(sample,c("timepoint","treatment","replicate"),"[_]") %>%  spread(allele,count)

head(reshaped_data)

# Make a copy of the reshaped_data so don't have to repeat the above
reshaped_data_original<-reshaped_data

###############################################################
#################### MOXIDECTIN VS CONTROL ####################

#for now select just moxidectin and control
#tmp<-subset(reshaped_data_original, reshaped_data_original$treatment == "M" | reshaped_data_original$treatment == "C")
#reshaped_data<-tmp

#add allele frequency column
#head(reshaped_data)
#reshaped_data <- reshaped_data %>% mutate(ref_freq=(REF_CNT / DEPTH))

#reshaped_data <- reshaped_data %>% mutate(timenum = as.numeric(str_remove(timepoint,"F")),snp_ID = str_c(CHROM,"_",POS))
#reshaped_data <- reshaped_data %>% mutate(scaled_timenum = scale(timenum))

#reshaped_data$treatment <- factor(reshaped_data$treatment)
#reshaped_data$replicate<- factor(reshaped_data$replicate)

#add data / likelihood weights
#reshaped_data$weight <- 1.0
#reshaped_data$weight[reshaped_data$timenum == 0] <- 1.0/6.0
#########


#key to fitting many models  is to nest the dataframe. We want to run one model for each SNP in the data:

#by_SNP <- reshaped_data %>% group_by( CHROM,POS,snp_ID) %>% nest()

#include replicate as a random effect using the lme4 library
#You also need the broom.mixed library for the next step too

# get_Mcoeff = Get the information regarding the SLOPE of the treated group relative to the control group (note, we are not including the intercept for either the control or treated here)
# This will using the tidy() function to put the model output for m1 (the better fit) into a table, from which the coefficients are then extracted and summed together to get the slope
# The slope will tell us whether the reference allele is reducing or increasing in frequency from F0 to F3 generations with treatment. 
# The pvalue (extracted below) tells us whether there is a significant difference between the control and the treated lines, but not whether it is important

# Testing larger dataset:

#library(broom.mixed)

#get_Mcoeff <- function(m1) {
#if (is.infinite(logLik(m1))) { 
#return(NaN)} else {
#coef <- tidy(m1)
#return (coef[coef$term == "scaled_timenum",]$estimate + coef[coef$term == "treatmentM:scaled_timenum",]$estimate)
#}
#}

# This next function runs two different models, m1 and m2. The replicate (line) is a random factor variable, and we assume that all lines have equal variance (they do have very similar variance), and all are tested independently - they're not grouped by treatment
# From the anova used to compare the two models we extract the p value and the moxidectin slope coefficient (calculated using the function above)


#do_both_fits <- function(df) {
#            m1 <- lmer(ref_freq ~ treatment*scaled_timenum+(1|replicate),data=df,weights=weight)
#            m2 <- lmer(ref_freq ~ treatment+scaled_timenum+(1|replicate),data=df,weights=weight)
#            a <- anova(m2,m1,test="LRT")    
#            return(list(p=a$`Pr(>Chisq)`[2],m1.coef=get_Mcoeff(m1)))
#}

#this also catches warnings - and lets the function continue past SNPs which have insufficient depth coverage - i.e. some samples have zero coverage. Otherwise it fails and nothing is output
#do_both_fits_safely <- function(df) {
#tryCatch(do_both_fits(df), error = function(err){NA}, warning=function(w){})
#}

# This actually runs the three functions above and outputs the results:
#all_SNP <- SNP_results <- by_SNP %>% 
#  mutate(fit_results = map(data, do_both_fits_safely))

# Got a warning with test data not noticed before when running MvsC (without F0 lines, but could have been there, who knows):
# fixed-effect model matrix is rank deficient so dropping 1 column / coefficient
# https://stackoverflow.com/questions/37090722/lme4lmer-reports-fixed-effect-model-matrix-is-rank-deficient-do-i-need-a-fi had some helpful info, although didn't understand it all, it is possible that trying a model without REML hasn't worked due to nested data?? Or missing data? Could be where have a depth of zero or something in some. 

#reshape a bit to get coef and p-values from list
#p_values_MC <- all_SNP %>% hoist(fit_results,"p","m1.coef") %>% mutate(neglogp=-1*log10(p))

# Make a binary Robject file that could be easily reloaded if desired to replot the plot!
#save(p_values_MC, file=(paste("MOX_vs_CTL_F0_results",my_date,".obj", sep="")))



###############################################################
#################### IVERMECTIN VS CONTROL ####################

#select just ivermectin and control 
tmp<-subset(reshaped_data_original, reshaped_data_original$treatment == "I" | reshaped_data_original$treatment == "C")
reshaped_data<-tmp

#add allele frequency column
reshaped_data <- reshaped_data %>% mutate(ref_freq=(REF_CNT / DEPTH))

# Convert the generation data to a scaled numerical number, and create a new column called snp_ID
reshaped_data <- reshaped_data %>% mutate(timenum = as.numeric(str_remove(timepoint,"F")),snp_ID = str_c(CHROM,"_",POS))
reshaped_data <- reshaped_data %>% mutate(scaled_timenum = scale(timenum))

# Make treatment and replicates (i.e. lines) as factors
reshaped_data$treatment <- factor(reshaped_data$treatment)
reshaped_data$replicate<- factor(reshaped_data$replicate)

#add data / likelihood weights
reshaped_data$weight <- 1.0
reshaped_data$weight[reshaped_data$timenum == 0] <- 1.0/6.0

#########

#key to fitting many models is to nest the dataframe. We want to run one model for each SNP in the data:

by_SNP <- reshaped_data %>% group_by( CHROM,POS,snp_ID) %>% nest()

#include replicate as a random effect using the lme4 library
#You also need the broom.mixed library for the next step too

# get_Icoeff = Get the information regarding the SLOPE of the treated group relative to the control group (note, we are not including the intercept for either the control or treated here)
# This will using the tidy() function to put the model output for m1 (the better fit) into a table, from which the coefficients are then extracted and summed together to get the slope
# The slope will tell us whether the reference allele is reducing or increasing in frequency from F1 to F3 generations with treatment. 
# The pvalue (extracted below) tells us whether there is a significant difference between the control and the treated lines, but not whether it is important

get_Icoeff <- function(m1) {
if (is.infinite(logLik(m1))) { 
return(NaN)} else {
coef <- tidy(m1)
return (coef[coef$term == "scaled_timenum",]$estimate + coef[coef$term == "treatmentI:scaled_timenum",]$estimate)
}
}

# This next function runs two different models, m1 and m2. The replicate (line) is a random factor variable, and we assume that all lines have equal variance (they do have very similar variance), and all are tested independently - they're not grouped by treatment
# From the anova used to compare the two models we extract the p value and the moxidectin slope coefficient (calculated using the function above)

do_both_fits <- function(df) {
m1 <- lmer(ref_freq ~ treatment*scaled_timenum+(1|replicate),data=df)
m2 <- lmer(ref_freq ~ treatment+scaled_timenum+(1|replicate),data=df)
a <- anova(m2,m1,test="LRT") 
return(list(p=a$`Pr(>Chisq)`[2],m1.coef=get_Icoeff(m1)))
}

#this also catches warnings - and lets the function continue past SNPs which have insufficient depth coverage - i.e. some samples have zero coverage. Otherwise it fails and nothing is output
do_both_fits_safely <- function(df) {
tryCatch(do_both_fits(df), error = function(err){NA}, warning=function(w){})
}

# This actually runs the three functions above and outputs the results:
SNP_results <- by_SNP %>% 
mutate(fit_results = map(data, do_both_fits_safely))

#reshape a bit to get coef and p-values from list, and get the neg log10 of the pvalue
p_values_IC <- SNP_results %>% hoist(fit_results,"p","m1.coef") %>% mutate(neglogp=-1*log10(p))

# Make an binary Robject file that could be easily reloaded if desired to replot the plot!
save(p_values_IC, file=(paste("IVM_vs_CTL_F0_results",my_date,".obj", sep="")))




###############################################################
#################### SELECTION (I & M) VS CONTROL #############


# Keep all three treatments: don't subselect
reshaped_data<-reshaped_data_original

#change both I and M to S for selected
reshaped_data$treatment<-gsub("I","S", reshaped_data$treatment)
reshaped_data$treatment<-gsub("M","S", reshaped_data$treatment)

#add allele frequency column
reshaped_data <- reshaped_data %>% mutate(ref_freq=(REF_CNT / DEPTH))

# Convert the generation data to a scaled numerical number, and create a new column called snp_ID
reshaped_data <- reshaped_data %>% mutate(timenum = as.numeric(str_remove(timepoint,"F")),snp_ID = str_c(CHROM,"_",POS))
reshaped_data <- reshaped_data %>% mutate(scaled_timenum = scale(timenum))

# Make treatment and replicates (i.e. lines) as factors
reshaped_data$treatment <- factor(reshaped_data$treatment)
reshaped_data$replicate<- factor(reshaped_data$replicate)

#add data / likelihood weights
reshaped_data$weight <- 1.0
reshaped_data$weight[reshaped_data$timenum == 0] <- 1.0/9.0

#########

#key to fitting many models is to nest the dataframe. We want to run one model for each SNP in the data:

by_SNP <- reshaped_data %>% group_by( CHROM,POS,snp_ID) %>% nest()

#include replicate as a random effect using the lme4 library
#You also need the broom.mixed library for the next step too

# get_Scoeff = Get the information regarding the SLOPE of the treated group relative to the control group (note, we are not including the intercept for either the control or treated here)
# This will using the tidy() function to put the model output for m1 (the better fit) into a table, from which the coefficients are then extracted and summed together to get the slope
# The slope will tell us whether the reference allele is reducing or increasing in frequency from F1 to F3 generations with treatment. 
# The pvalue (extracted below) tells us whether there is a significant difference between the control and the treated lines, but not whether it is important

get_Scoeff <- function(m1) {
if (is.infinite(logLik(m1))) { 
return(NaN)} else {
coef <- tidy(m1)
return (coef[coef$term == "scaled_timenum",]$estimate + coef[coef$term == "treatmentS:scaled_timenum",]$estimate)
}
}

# This next function runs two different models, m1 and m2. The replicate (line) is a random factor variable, and we assume that all lines have equal variance (they do have very similar variance), and all are tested independently - they're not grouped by treatment
# From the anova used to compare the two models we extract the p value and the moxidectin slope coefficient (calculated using the function above)

do_both_fits <- function(df) {
m1 <- lmer(ref_freq ~ treatment*scaled_timenum+(1|replicate),data=df)
m2 <- lmer(ref_freq ~ treatment+scaled_timenum+(1|replicate),data=df)
a <- anova(m2,m1,test="LRT") 
return(list(p=a$`Pr(>Chisq)`[2],m1.coef=get_Scoeff(m1)))
}

#this also catches warnings - and lets the function continue past SNPs which have insufficient depth coverage - i.e. some samples have zero coverage. Otherwise it fails and nothing is output
do_both_fits_safely <- function(df) {
tryCatch(do_both_fits(df), error = function(err){NA}, warning=function(w){})
}

# This actually runs the three functions above and outputs the results:
SNP_results <- by_SNP %>% 
mutate(fit_results = map(data, do_both_fits_safely))

#reshape a bit to get coef and p-values from list, and get the neg log10 of the pvalue
p_values_SEL <- SNP_results %>% hoist(fit_results,"p","m1.coef") %>% mutate(neglogp=-1*log10(p))

# Make an binary Robject file that could be easily reloaded if desired to replot the plot!
save(p_values_SEL, file=(paste("SEL_vs_CTL_F0_results",my_date,".obj", sep="")))




###############################################################
#################### IVERMECTIN VS MOXIDECTIN ####################

#for now select just ivermectin and moxidectin
tmp<-subset(reshaped_data_original, reshaped_data_original$treatment == "I" | reshaped_data_original$treatment == "M")
reshaped_data<-tmp

#add allele frequency column
reshaped_data <- reshaped_data %>% mutate(ref_freq=(REF_CNT / DEPTH))

reshaped_data <- reshaped_data %>% mutate(timenum = as.numeric(str_remove(timepoint,"F")),snp_ID = str_c(CHROM,"_",POS))
reshaped_data <- reshaped_data %>% mutate(scaled_timenum = scale(timenum))

reshaped_data$treatment <- factor(reshaped_data$treatment)
reshaped_data$replicate<- factor(reshaped_data$replicate)

#add data / likelihood weights
reshaped_data$weight <- 1.0
reshaped_data$weight[reshaped_data$timenum == 0] <- 1.0/6.0
#########


#key to fitting many models  is to nest the dataframe. We want to run one model for each SNP in the data:

by_SNP <- reshaped_data %>% group_by( CHROM,POS,snp_ID) %>% nest()

#include replicate as a random effect using the lme4 library
#You also need the broom.mixed library for the next step too

# get_Mcoeff = Get the information regarding the SLOPE of the treated group relative to the control group (note, we are not including the intercept for either the control or treated here)
# This will using the tidy() function to put the model output for m1 (the better fit) into a table, from which the coefficients are then extracted and summed together to get the slope
# The slope will tell us whether the reference allele is reducing or increasing in frequency from F0 to F3 generations with treatment. 
# The pvalue (extracted below) tells us whether there is a significant difference between the control and the treated lines, but not whether it is important

# Testing larger dataset:

library(broom.mixed)

get_Mcoeff <- function(m1) {
if (is.infinite(logLik(m1))) { 
return(NaN)} else {
coef <- tidy(m1)
return (coef[coef$term == "scaled_timenum",]$estimate + coef[coef$term == "treatmentM:scaled_timenum",]$estimate)
}
}

# This next function runs two different models, m1 and m2. The replicate (line) is a random factor variable, and we assume that all lines have equal variance (they do have very similar variance), and all are tested independently - they're not grouped by treatment
# From the anova used to compare the two models we extract the p value and the moxidectin slope coefficient (calculated using the function above)


do_both_fits <- function(df) {
            m1 <- lmer(ref_freq ~ treatment*scaled_timenum+(1|replicate),data=df,weights=weight)
            m2 <- lmer(ref_freq ~ treatment+scaled_timenum+(1|replicate),data=df,weights=weight)
            a <- anova(m2,m1,test="LRT")    
            return(list(p=a$`Pr(>Chisq)`[2],m1.coef=get_Mcoeff(m1)))
}

#this also catches warnings - and lets the function continue past SNPs which have insufficient depth coverage - i.e. some samples have zero coverage. Otherwise it fails and nothing is output
do_both_fits_safely <- function(df) {
tryCatch(do_both_fits(df), error = function(err){NA}, warning=function(w){})
}

# This actually runs the three functions above and outputs the results:
all_SNP <- SNP_results <- by_SNP %>% 
  mutate(fit_results = map(data, do_both_fits_safely))

#reshape a bit to get coef and p-values from list
p_values_IM <- all_SNP %>% hoist(fit_results,"p","m1.coef") %>% mutate(neglogp=-1*log10(p))

# Make a binary Robject file that could be easily reloaded if desired to replot the plot!
save(p_values_IM, file=(paste("IVM_vs_MOX_F0_results",my_date,".obj", sep="")))





###############################################################
#################### IVERMECTIN VS MOXIDECTIN VS CONTROL #############


# Keep all three treatments: don't subselect
reshaped_data<-reshaped_data_original

#add allele frequency column
reshaped_data <- reshaped_data %>% mutate(ref_freq=(REF_CNT / DEPTH))

# Convert the generation data to a scaled numerical number, and create a new column called snp_ID
reshaped_data <- reshaped_data %>% mutate(timenum = as.numeric(str_remove(timepoint,"F")),snp_ID = str_c(CHROM,"_",POS))
reshaped_data <- reshaped_data %>% mutate(scaled_timenum = scale(timenum))

# Make treatment and replicates (i.e. lines) as factors
reshaped_data$treatment <- factor(reshaped_data$treatment)
reshaped_data$replicate<- factor(reshaped_data$replicate)

#add data / likelihood weights
reshaped_data$weight <- 1.0
reshaped_data$weight[reshaped_data$timenum == 0] <- 1.0/9.0

#########

#key to fitting many models is to nest the dataframe. We want to run one model for each SNP in the data:

by_SNP <- reshaped_data %>% group_by( CHROM,POS,snp_ID) %>% nest()

#include replicate as a random effect using the lme4 library
#You also need the broom.mixed library for the next step too

# get_Scoeff = Get the information regarding the SLOPE of the treated group relative to the control group (note, we are not including the intercept for either the control or treated here)
# This will using the tidy() function to put the model output for m1 (the better fit) into a table, from which the coefficients are then extracted and summed together to get the slope
# The slope will tell us whether the reference allele is reducing or increasing in frequency from F1 to F3 generations with treatment. 
# The pvalue (extracted below) tells us whether there is a significant difference between the control and the treated lines, but not whether it is important

get_Scoeff <- function(m1) {
if (is.infinite(logLik(m1))) { 
return(NaN)} else {
coef <- tidy(m1)
return (coef[coef$term == "scaled_timenum",]$estimate + coef[coef$term == "treatmentS:scaled_timenum",]$estimate)
}
}

# This next function runs two different models, m1 and m2. The replicate (line) is a random factor variable, and we assume that all lines have equal variance (they do have very similar variance), and all are tested independently - they're not grouped by treatment
# From the anova used to compare the two models we extract the p value and the moxidectin slope coefficient (calculated using the function above)

do_both_fits <- function(df) {
m1 <- lmer(ref_freq ~ treatment*scaled_timenum+(1|replicate),data=df)
m2 <- lmer(ref_freq ~ treatment+scaled_timenum+(1|replicate),data=df)
a <- anova(m2,m1,test="LRT") 
return(list(p=a$`Pr(>Chisq)`[2],m1.coef=get_Scoeff(m1)))
}

#this also catches warnings - and lets the function continue past SNPs which have insufficient depth coverage - i.e. some samples have zero coverage. Otherwise it fails and nothing is output
do_both_fits_safely <- function(df) {
tryCatch(do_both_fits(df), error = function(err){NA}, warning=function(w){})
}

# This actually runs the three functions above and outputs the results:
SNP_results <- by_SNP %>% 
mutate(fit_results = map(data, do_both_fits_safely))

#reshape a bit to get coef and p-values from list, and get the neg log10 of the pvalue
p_values_IMC <- SNP_results %>% hoist(fit_results,"p","m1.coef") %>% mutate(neglogp=-1*log10(p))

# Make an binary Robject file that could be easily reloaded if desired to replot the plot!
save(p_values_IMC, file=(paste("IVM_vs_MOX_vs_CTL_F0_results",my_date,".obj", sep="")))







############
# MOX vs CTL = p_values_MC, file="MOX_vs_CTL_F0_results_date.obj"
# IVM vs CTL = p_values_IC, file="IVM_vs_CTL_F0_results_date.obj"
# IVM AND MOX ("SEL") vs CTL = p_values_SEL, file="SELECTED_vs_CTL_F0_results_date.obj"
# IVM vs MOX = p_values_IM, file="IVM_vs_MOX_F0_results_date.obj"
# IVM vs MOX vs CTL = p_values_IMC, file="IVM_vs_MOX_vs_CTL_F0_results_date.obj"
