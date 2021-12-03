# CoDa analysis of fractions spent in each of 8 behavioral classes

require(data.table)
require(compositions)
require(npmv)
library(plyr)
library(rstatix)

# for visualization
require(ggplot2)
library(ggpubr)

# load data set
dat <- read.table("/Users/bergeler/Documents/Paper_Ugne/Data/MMM_8_data_set.csv",header = TRUE, sep = ",")

dat <- within(dat,{
  group <- factor(group)
  mouse <- factor(mouse)
  day <- factor(day)
  class <- factor(class)
})

# there are some NaNs in the data set!(C57bl-FEMAL, mouse 16)
dat <- dat[!is.na(dat$prob),]

# get number of zero counts
sum(dat$prob < 10^(-15)) # --> no zeros

# bring data into right format for later use
tab <- dcast(data = setDT(dat),formula = group + mouse + day ~ class,fun.aggregate = sum,value.var = "prob")
tab <- as.data.frame(tab)

# Make data compositional
dataset.compositional = acomp(tab[,4:11])

# calculate ilr coordinates
dataset.ilr <- ilr(dataset.compositional)

# make data sets with ilr coordinates
dat_ilr = cbind(tab[,1:3],dataset.ilr)

#--- MANOVA ---#
# 
ggboxplot(
  dat_ilr, x = "group", y = c("V1","V2","V3","V4","V5","V6","V7"), 
  merge = TRUE, palette = "jco"
) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
  ggtitle('All days') + 
  theme(legend.position="right") +
  theme(plot.title = element_text(hjust = 0.5))

dat1_ilr <- dat_ilr[dat_ilr$day == 1,]
dat2_ilr <- dat_ilr[dat_ilr$day == 2,]
dat3_ilr <- dat_ilr[dat_ilr$day == 3,]
dat4_ilr <- dat_ilr[dat_ilr$day == 4,]

pdf("/Users/bergeler/Documents/Paper_Ugne/Figures/ilr_day1.pdf", width = 7, height = 4)
ggboxplot(
  dat1_ilr, x = "group", y = c("V1","V2","V3","V4","V5","V6","V7"), 
  merge = TRUE, palette = "jco"
) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1)) + 
  ggtitle('Day 1') + 
  theme(legend.position="right") +
  theme(plot.title = element_text(hjust = 0.5))+
  ylab('ilr values') +
  xlab('')
dev.off()

pdf("/Users/bergeler/Documents/Paper_Ugne/Figures/ilr_day2.pdf", width = 7, height = 4)
ggboxplot(
  dat2_ilr, x = "group", y = c("V1","V2","V3","V4","V5","V6","V7"), 
  merge = TRUE, palette = "jco"
) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1)) + 
  ggtitle('Day 2') + 
  theme(legend.position="right") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylab('ilr values') +
  xlab('')
dev.off()

pdf("/Users/bergeler/Documents/Paper_Ugne/Figures/ilr_day3.pdf", width = 7, height = 4)
ggboxplot(
  dat3_ilr, x = "group", y = c("V1","V2","V3","V4","V5","V6","V7"), 
  merge = TRUE, palette = "jco"
) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1)) + 
  ggtitle('Day 3') + 
  theme(legend.position="right") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylab('ilr values') +
  xlab('')
dev.off()

pdf("/Users/bergeler/Documents/Paper_Ugne/Figures/ilr_day4.pdf", width = 7, height = 4)
ggboxplot(
  dat4_ilr, x = "group", y = c("V1","V2","V3","V4","V5","V6","V7"), 
  merge = TRUE, palette = "jco"
) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1)) + 
  ggtitle('Day 4') + 
  theme(legend.position="right") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylab('ilr values') +
  xlab('')
dev.off()

# test multivariate normality
dat1_ilr %>%
  select(setdiff(colnames(dat1_ilr),c("group","mouse","day"))) %>%
  mshapiro_test()
# significant --> assumptions NOT satisfied --> use non-parametric test

dat2_ilr %>%
  select(setdiff(colnames(dat2_ilr),c("group","mouse","day"))) %>%
  mshapiro_test()
# significant --> assumptions NOT satisfied --> use non-parametric test

dat3_ilr %>%
  select(setdiff(colnames(dat3_ilr),c("group","mouse","day"))) %>%
  mshapiro_test()
# significant --> assumptions NOT satisfied --> use non-parametric test

dat4_ilr %>%
  select(setdiff(colnames(dat4_ilr),c("group","mouse","day"))) %>%
  mshapiro_test()
# significant --> assumptions NOT satisfied --> use non-parametric test

# non-parametric test 
# nonpartest(V1 | V2 | V3 | V4 | V5 | V6 | V7 ~ group, data = dat1_ilr) 
ssnonpartest(V1 | V2 | V3 | V4 | V5 | V6 | V7 ~ group, data = dat1_ilr, alpha = 0.05, factors.and.variables = TRUE)
# Rejected hypothesis of equality between two factor levels:
# L7-Tsc1-HOM L7-Tsc1-NEG (x)
# L7-Tsc1-HET L7-Tsc1-HOM

# Cntnap2-NEG L7-Tsc1-NEG 
# Cntnap2-NEG L7-Tsc1-HET
# Cntnap2-HOMO L7-Tsc1-NEG
# Cntnap2-HOMO L7-Tsc1-HOM
# Cntnap2-HOMO L7-Tsc1-HET
# Cntnap2-HET L7-Tsc1-NEG
# Cntnap2-HET L7-Tsc1-HOM 
# Cntnap2-HET L7-Tsc1-HET 

# C57bl-MALE L7-Tsc1-NEG (x)
# C57bl-MALE L7-Tsc1-HOM
# C57bl-MALE L7-Tsc1-HET 
# C57bl-MALE Cntnap2-NEG
# C57bl-MALE Cntnap2-HOMO
# C57bl-MALE Cntnap2-HET 

# C57bl-FEMALE L7-Tsc1-NEG (x)
# C57bl-FEMALE L7-Tsc1-HOM
# C57bl-FEMALE L7-Tsc1-HET 
# C57bl-FEMALE Cntnap2-NEG  
# C57bl-FEMALE Cntnap2-HET  

# C57bl-FEMALE C57bl-MALE (x)

# --> all groups are significantly different except (6): 
# C57bl-FEMALE vs Cntnap2-HOMO (x)

# Cntnap2-HET vs Cntnap2-HOMO
# Cntnap2-HET vs Cntnap2-NEG
# Cntnap2-HOMO vs Cntnap2-NEG

# L7-Tsc1-HET vs L7-Tsc1-NEG (x)

# Cntnap2-NEG vs L7-Tsc1-HOM

# nonpartest(V1 | V2 | V3 | V4 | V5 | V6 | V7 ~ group, data = dat2_ilr)
ssnonpartest(V1 | V2 | V3 | V4 | V5 | V6 | V7 ~ group, data = dat2_ilr, alpha = 0.05, factors.and.variables = TRUE)
# Rejected hypothesis of equality between two factor levels:
# L7-Tsc1-HET L7-Tsc1-HOM is rejected 
# Cntnap2-NEG L7-Tsc1-NEG is rejected 
# Cntnap2-NEG L7-Tsc1-HOM is rejected 
# Cntnap2-NEG L7-Tsc1-HET is rejected 
# Cntnap2-HET L7-Tsc1-NEG is rejected 
# Cntnap2-HET L7-Tsc1-HET is rejected 
# C57bl-MALE L7-Tsc1-NEG is rejected 
# C57bl-MALE L7-Tsc1-HOM is rejected 
# C57bl-MALE L7-Tsc1-HET is rejected 
# C57bl-MALE Cntnap2-NEG is rejected 
# C57bl-MALE Cntnap2-HOMO is rejected 
# C57bl-MALE Cntnap2-HET is rejected 
# C57bl-FEMALE L7-Tsc1-NEG is rejected 
# C57bl-FEMALE L7-Tsc1-HOM is rejected 
# C57bl-FEMALE L7-Tsc1-HET is rejected 
# C57bl-FEMALE C57bl-MALE is rejected 

# --> all groups are significantly different except (12): 
# L7-Tsc1-HET vs L7-Tsc1-NEG 
# L7-Tsc1-NEG vs L7-Tsc1-HOM

# L7-Tsc1-HET vs Cntnap2-HOMO
# L7-Tsc1-NEG vs Cntnap2-HOMO
# L7-Tsc1-HOM vs Cntnap2-HOMO
# L7-Tsc1-HOM vs Cntnap2-HET

# Cntnap2-NEG vs Cntnap2-HOMO
# Cntnap2-NEG vs Cntnap2-HET
# Cntnap2-HET vs Cntnap2-HOMO

# Cntnap2-NEG vs C57bl-FEMALE
# Cntnap2-HET vs C57bl-FEMALE
# Cntnap2-HOMO vs C57bl-FEMALE 

ssnonpartest(V1 | V2 | V3 | V4 | V5 | V6 | V7 ~ group, data = dat3_ilr, alpha = 0.05, factors.and.variables = TRUE)
# Rejected hypothesis of equality between two factor levels:
# L7-Tsc1-HOM L7-Tsc1-NEG is rejected 
# Cntnap2-NEG L7-Tsc1-NEG is rejected 
# Cntnap2-NEG L7-Tsc1-HOM is rejected 
# Cntnap2-NEG L7-Tsc1-HET is rejected 
# Cntnap2-HOMO L7-Tsc1-NEG is rejected 
# Cntnap2-HOMO L7-Tsc1-HOM is rejected 
# Cntnap2-HOMO L7-Tsc1-HET is rejected 
# Cntnap2-HET L7-Tsc1-NEG is rejected 
# Cntnap2-HET L7-Tsc1-HOM is rejected 
# Cntnap2-HET L7-Tsc1-HET is rejected 
# C57bl-MALE L7-Tsc1-NEG is rejected 
# C57bl-MALE L7-Tsc1-HOM is rejected 
# C57bl-MALE L7-Tsc1-HET is rejected 
# C57bl-MALE Cntnap2-NEG is rejected 
# C57bl-MALE Cntnap2-HOMO is rejected 
# C57bl-MALE Cntnap2-HET is rejected 
# C57bl-FEMALE L7-Tsc1-NEG is rejected 
# C57bl-FEMALE L7-Tsc1-HOM is rejected 
# C57bl-FEMALE L7-Tsc1-HET is rejected 
# C57bl-FEMALE C57bl-MALE is rejected 

# --> all groups are significantly different except (8): 
# L7-Tsc1-HOM vs L7-Tsc1-HET
# L7-Tsc1-NEG vs L7-Tsc1-HET

# Cntnap2-NEG vs Cntnap2-HOMO
# Cntnap2-NEG vs Cntnap2-HET
# Cntnap2-HOMO vs Cntnap2-HET

# Cntnap2-NEG vs C57bl-FEMALE
# Cntnap2-HOMO vs C57bl-FEMALE
# Cntnap2-HET vs C57bl-FEMALE

ssnonpartest(V1 | V2 | V3 | V4 | V5 | V6 | V7 ~ group, data = dat4_ilr, alpha = 0.05, factors.and.variables = TRUE)
# Rejected hypothesis of equality between two factor levels:
# Cntnap2-NEG L7-Tsc1-NEG is rejected 
# Cntnap2-NEG L7-Tsc1-HET is rejected 
# Cntnap2-HOMO L7-Tsc1-NEG is rejected 
# Cntnap2-HOMO L7-Tsc1-HET is rejected 
# Cntnap2-HET L7-Tsc1-NEG is rejected 
# Cntnap2-HET L7-Tsc1-HET is rejected 
# C57bl-MALE L7-Tsc1-NEG is rejected 
# C57bl-MALE L7-Tsc1-HOM is rejected 
# C57bl-MALE L7-Tsc1-HET is rejected 
# C57bl-MALE Cntnap2-NEG is rejected 
# C57bl-MALE Cntnap2-HOMO is rejected 
# C57bl-MALE Cntnap2-HET is rejected 
# C57bl-FEMALE L7-Tsc1-NEG is rejected 
# C57bl-FEMALE L7-Tsc1-HOM is rejected 
# C57bl-FEMALE L7-Tsc1-HET is rejected 
# C57bl-FEMALE C57bl-MALE is rejected 

# --> all groups are significantly different except (12): 
# Cntnap2-NEG vs Cntnap2-HOMO
# Cntnap2-NEG vs Cntnap2-HET
# Cntnap2-HOMO vs Cntnap2-HET 

# L7-Tsc1-NEG vs L7-Tsc1-HET
# L7-Tsc1-NEG vs L7-Tsc1-HOM
# L7-Tsc1-HOM vs L7-Tsc1-HET

# Cntnap2-NEG vs L7-Tsc1-HOM
# Cntnap2-HOMO vs L7-Tsc1-HOM 
# Cntnap2-HET vs L7-Tsc1-HOM

# Cntnap2-NEG vs C57bl-FEMALE
# Cntnap2-HOMO vs C57bl-FEMALE 
# Cntnap2-HET vs C57bl-FEMALE 
