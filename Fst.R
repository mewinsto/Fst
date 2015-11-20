##IMPORTS AND FORMATS SNP MATRIX ("GENO_MATRIX") FROM .GENO FILE (pyRAD OUTPUT)
#GENO_MATRIX HAS SAMPLES AS ROWS AND SNPS FOR COLUMNS
setwd("~/Desktop/FST/lucanoides//")
geno_data = scan("all_lanes_lucanoides_10k.usnps.geno",what = "character")
L = length(geno_data)
S = nchar(geno_data[[1]])
geno_matrix = matrix(data = NA, nrow = S, ncol = L)
for (i in 1:L){
   genostring = strsplit(geno_data[[i]], split="")
   for (j in 1:S){
     geno_matrix[j,i] = as.numeric(genostring[[1]][j])
   }
}
geno_matrix[geno_matrix == 9] = NA

##IMPORTS SAME SET OF SAMPLES FROM .SAMPLES FILE TO CONSTRUCT POPULATION VECTORS
samples = read.table("all_lanes_lucanoides_10k.sample",header=FALSE)
samples = samples$V1

#READS IN BURCH POP TO MAKE A TEST VECTOR, THEN CONSTRUCTS NULL VECTOR BASED ON POPS
testpop = read.table("all_lanes_lucanoides_10k_pop.txt",header=FALSE)
test_vector = {}
for (i in 1:S){
  test_vector[i] = testpop$V2[match(samples[i],testpop$V1)]
}

test_sub = {}
for (i in 1:S){
  if (test_vector[i] == 1){
    test_sub = c(test_sub,i)
  }
}
  
null_vector = sample(c(rep(1,sum(test_vector)),rep(0,S-sum(test_vector))),size=S)
null_sub = {}
for (i in 1:S){
  if (null_vector[i] == 1){
    null_sub = c(null_sub,i)
  }
}

##THIS MODULE CALCULATES PI (NUCLEOTIDE DIVERSITY)
pi = {}
for (i in 1:L){
  locus = {}
  locus = geno_matrix[,i]
  locus = locus[!is.na(locus)]
  locus_sum = 0
  for (j in 2:length(locus)){
    jcount = j - 1
    for (k in 1:jcount){
      diff = 0
      diff = abs(locus[j]-locus[k])
      locus_sum = locus_sum + diff
    }
  }
  locus_sum = (2*locus_sum)/(length(locus)^2) 
  pi[i] = locus_sum
}

write(pi,file="pi_lucanoides",ncol = 1)

##THIS MODULE CALCULATES FREQUENCIES AND FST FOR TEST
#Counts and mean frequencies of alleles
geno_means = vector(mode = "numeric", length = L)
geno_cover = vector(mode = "numeric", length = L)
geno_cover_p1 = vector(mode = "numeric", length = length(test_sub))
geno_cover_p2 = vector(mode = "numeric", length = L - length(test_sub))
for (i in 1:L){
  geno_means[i] = mean(geno_matrix[,i],na.rm=TRUE)/2
}
for (i in 1:L){
  geno_cover[i] = sum(!is.na(geno_matrix[,i]))
  geno_cover_p1[i] = sum(!is.na(geno_matrix[test_sub,i]))
  geno_cover_p2[i] = sum(!is.na(geno_matrix[-test_sub,i]))
}

#Use test_sub to look at Fst values
pop_1 = {}
pop_2 = {}
pop_mean = {}
FST = {}
for (i in 1:L){
  pop_1[i] = mean(geno_matrix[test_sub,i],na.rm=TRUE)/2
  pop_2[i] = mean(geno_matrix[-test_sub,i],na.rm=TRUE)/2
#  pop_mean[i] = (pop_1[i] + pop_2[i])/2
  FST[i] = (geno_cover_p1[i]*(abs(geno_means[i] - pop_1[i])^2) + geno_cover_p2[i]*abs(geno_means[i] - pop_2[i])^2)/(geno_cover[i]*geno_means[i]*(1 - geno_means[i]))
}

FST = abs(FST)

##WRITE FST TO FILE
write(FST,file="Fst_test_lucanoides.txt",ncolumns=1)
write(geno_cover,file="geno_cover_lucanoides.txt",ncolumns=1)

##THIS MODULE CALCULATES FREQUENCIES AND FST FOR NULL
#Counts and mean frequencies of alleles
geno_cover_p1_N = vector(mode = "numeric", length = length(null_sub))
geno_cover_p2_N = vector(mode = "numeric", length = L - length(null_sub))

for (i in 1:L){
  geno_cover[i] = sum(!is.na(geno_matrix[,i]))
  geno_cover_p1_N[i] = sum(!is.na(geno_matrix[null_sub,i]))
  geno_cover_p2_N[i] = sum(!is.na(geno_matrix[-null_sub,i]))
}

#Use null_sub to look at Fst values
pop_1_N = {}
pop_2_N = {}
pop_mean = {}
FST_N = {}
for (i in 1:L){
  pop_1_N[i] = mean(geno_matrix[null_sub,i],na.rm=TRUE)/2
  pop_2_N[i] = mean(geno_matrix[-null_sub,i],na.rm=TRUE)/2
  #  pop_mean[i] = (pop_1[i] + pop_2[i])/2
  FST_N[i] = (geno_cover_p1_N[i]*(abs(geno_means[i] - pop_1_N[i])^2) + geno_cover_p2_N[i]*abs(geno_means[i] - pop_2_N[i])^2)/(geno_cover[i]*geno_means[i]*(1 - geno_means[i]))
}

FST_N = abs(FST_N)

#WRITE FST TO FILE
write(FST_N,file="Fst_null_lucanoides.txt",ncolumns=1)

##THIS MODULE CALCULATES MAF FOR EACH LOCUS
af = {}
for (i in 1:L){
  L_maf = sum(!is.na(geno_matrix[,i]))
  L_sum = sum(geno_matrix[,i],na.rm=TRUE)
  af[i] = L_sum/(2*L_maf)
}

#MAKE MAF MINOR
maf = {}
for (i in 1:L){
  if (af[i] > 0.5){
    maf[i] = 1 - af[i]
  } else{
    maf[i] = af[i]
  }
}

#WRITE MAF OR AF TO FILE
write(maf,file="MAF_lucanoides.txt", ncolumns=1)
write(af, file="AF_lucanoides.txt", ncolumns=1)

##THIS MODULE CALCULATES H-W EQUILIBRIUM FOR EACH LOCUS
HW_chi = {}
for (i in 1:L){
  L_maf = sum(!is.na(geno_matrix[,i]))
  p = af[i]
  q = 1 - p
  exp_pp = L_maf*(p^2)
  exp_pq = L_maf*(2*p*q)
  exp_qq = L_maf*(q^2)
  act_pp = sum(geno_matrix[,i] == 2,na.rm=TRUE)
  act_pq = sum(geno_matrix[,i] == 1,na.rm=TRUE)
  act_qq = sum(geno_matrix[,i] == 0,na.rm=TRUE)
  HW_chi[i] = ((act_pp - exp_pp)^2/exp_pp) + ((act_pq - exp_pq)^2/exp_pq) + ((act_qq - exp_qq)^2/exp_qq)
}

#WRITE HARDY-WEINBERG CHI-SQUARED STAT (CAN DERIVE p-values FROM THIS)
write(HW_chi,file="HW_chi_lucanoides.txt",ncolumns=1)

#NUMBER OF INDIVIDUALS WITH SNP (FOR REMOVING SINGLETONS)
allele_class = {}
for (i in 1:L){
  allele_class[i] = maf[i]*geno_cover[i]
}

write(allele_class,file="AC_lucanoides.txt",ncolumns=1)
##READ IN FILES: NULL AND TEST AND PLOT AGAINST ONE ANOTHER, ALSO NUC FREQ
library(ggplot2)

pi = read.table("~/Desktop/FST/lucanoides/pi_lucanoides")
f_null = read.table("~/Desktop/FST/lucanoides/Fst_null_lucanoides.txt")
f_test = read.table("~/Desktop/FST/lucanoides/Fst_test_lucanoides.txt")
cover_test = read.table("~/Desktop/FST/lucanoides/geno_cover_lucanoides.txt") 
af = read.table("~/Desktop/FST/lucanoides/AF_lucanoides.txt")
HW_chi = read.table("~/Desktop/FST/lucanoides/HW_chi_lucanoides.txt")
allele_class = read.table("~/Desktop/FST/lucanoides/AC_lucanoides.txt")

FST_data = cbind(f_null,f_test,cover_test,pi,af,HW_chi,allele_class)
names(FST_data) = c("F_null","F_test","Cover","Pi","AF","HW","AC")

##SUBSETS FST DATA INTO COVERAGE GROUPS
for (i in min(cover_test):max(cover_test)){
  name = paste("FST_data_cover_",i, sep = "")
  assign(name,subset(FST_data, cover_test == i))
}

#SUBSETS FST DATA INTO GROUPS WITH FST = 1.0
FST_data_HT = subset(FST_data, F_test == 1.0)
FST_data_HN = subset(FST_data, F_null == 1.0)

#SUBSETS FST DATA INTO GROUPS WITH FST < 0.1
FST_data_LT = subset(FST_data, F_test < 0.1)
FST_data_LN = subset(FST_data, F_null < 0.1)

#SUBSETS FST DATA INTO GROUPS WITH Pi > 0.8
Hi_Pi = subset(FST_data, Pi > 0.8)
Lo_Pi = subset(FST_data, Pi < 0.2)

Pi_cat = cut(FST_data$Pi,breaks=10)

#SUBSETS FST DATA TO EXCLUDE SINGLETONS (AC == 1)
FST_data_AC = subset(FST_data, allele_class > 1)
FST_data_AC_HT = subset(FST_data_AC, F_test == 1)

##############
###PLOTTING###
##############

#Plot densities of these values
plot(density(FST_data$Pi,bw=0.01),xlim=c(0,1),main="Distribution of Major Allele Frequencies across Loci", xlab="Major Allele Frequency", sub="Species = E. lucanoides (G) -- N = 36,287 - bw = .01")
plot(density(FST_data$Cover,bw=0.5),xlim=c(4,16),main="Density Distribution of Individuals per Locus: E. lucanoides", xlab="Number of Individuals (Max = 19)",sub="N = 63,087, bw = 0.5")

##PLOT FST AGAINST GENO_COVER
plot(density(FST_data$F_test,na.rm=TRUE,bw=.004),ylim=c(0,60),main="Density Distribution for Fst Values: E. lucanoides (G)",xlab="Fst Values", sub = "N = 63,087   bw = 0.004")
lines(density(FST_data$F_null,na.rm=TRUE,bw=.004),col=2)

plot(FST_data$F_test,FST_data$Cover,col=2,main="Fst for all Loci by Number of Samples",xlab="Fst",ylab="Number of Individuals with Locus (Max: 161)",sub="Black = Null, Red = lucanoides/burchellii -- N = 253,297")
points(FST_data$F_null,FST_data$Cover)

#PLOTS NUCLEOTIDE DIVERSITY FOR HIGH VS. LOW FST LOCI
plot(density(FST_data_LT$Pi),main="High vs. Low Fst Loci Nucleotide Diversity Distributions (E. burchellii_BC)",xlim = c(0,1), xlab="Pi (Normalized)",sub="Black = Low (Fst < 0.1) [N = 55,574] -- Red = High (Fst = 1.0) [N = 4,866]")
lines(density(FST_data_HT$Pi),col=2)

plot(density(FST_data_HN$Cover,bw=1),xlim=c(0,161))
lines(density(FST_data_HT$Cover),col=2)

##Plots LOW DIVERSITY AND HIGH DIVERSITY LOCI
HP = ggplot(Hi_Pi, aes(x = Cover, y = F_test))
HP + geom_point(alpha = 0.05) + ggtitle(expression(atop("Fst and Individual Coverage for High Nucleotide Diversity Loci (>0.8)", atop(italic("E. drepanophorum (N = 15, L = 12,335)", "")))))

HP = ggplot(Lo_Pi, aes(x = Cover, y = F_test))
HP + geom_point(alpha = 0.05) + ggtitle(expression(atop("Fst and Individual Coverage for Low Nucleotide Diversity Loci (<0.2)", atop(italic("E. drepanophorum (N = 15, L = 7950)", "")))))


##PLOTS DENSITY DISTRIBUTION OF NUCLEOTIDE DIVERSITY FOR FST == 1.0
TOP_DENS = ggplot(FST_data_HT, aes(x = Pi))
TOP_DENS + geom_density(fill = "black") + ggtitle(expression(atop("Fst Distributions for Nucleotide Diversity Classes", atop(italic("E. burchellii & E. lucanoides (N = 85, L = 40,722)", ""))))) + labs(x = "Nucleotide Diversity", y = "Density")

##GGPLOT NUCLEOTIDE DIVERSITY CLASSES AGAINST FST FOR NULL AND TEST
p = ggplot(FST_data, aes(factor(Pi_cat), F_null))

p + geom_violin() + ggtitle(expression(atop("Fst Distributions for Nucleotide Diversity Classes", atop(italic("E. burchellii (N = 75, L = 168,518)", ""))))) + labs(x = "Nucleotide Diversity Class", y = "Fst Value")
p + geom_boxplot() + ggtitle(expression(atop("Fst Distributions for Nucleotide Diversity Classes: Null Hypothesis", atop(italic("E. lucanoides (N = 15, L = 63,087", ""))))) + labs(x = "Nucleotide Diversity Class", y = "Fst Value")

##PLOT DENSITY DISTRIBUTION OF NUCLEOTIDE DIVERSITY
DENS = ggplot(FST_data, aes(x = Pi))
DENS + geom_density(fill = "black") + labs(x = "Nucleotide Diversity", y = "Density")

hist(FST_data_HT$Cover,xlim=c(0,19),breaks=15,col=2,main="Histogram of Loci with FST = 1.0: E. lucanoides",xlab="Number of Individuals Sampled",ylab="Count",sub="Red = E. lucanoides : N = 32,879; Blue = Null : N = 3,791")
hist(FST_data_HN$Cover,xlim=c(0,37),breaks=9,add=T,col=4)

m = ggplot(FST_data_HT, aes(x=Cover))
m + geom_histogram(binwidth = 1, fill = "darkslategray4") + geom_histogram(data = FST_data_HN, fill = "firebrick4", binwidth = 1)

d = ggplot(Hi_Pi, aes(x=FST))

#write(FST_data_HT$Cover,file="FST_HT_burchellii_SCA.txt",ncol=1)
#write(FST_data_HN$Cover,file="FST_HN_burchellii_SCA.txt",ncol=1)

norm_FST = table(FST_data_HT$Cover)/(table(FST_data$Cover)[1:73])
plot(norm_FST,ylim=c(0,0.8),main="Percentage of Loci with Fst = 1.0: E. lucanoides (G)",xlab="Individual Coverage",ylab="Percentage of Loci",sub="N = 19, L = 36,287: Dotted line = 0.02 ")
abline(h=0.5,lty=3)

norm_df = as.data.frame(norm_FST)
q = ggplot(norm_df, aes(x = Freq))
qplot(norm_df$Var1, norm_df$Freq, ylim=c(0,1),size=3)
q + geom_histogram(binwidth = 0.01)
#write(norm_FST,file="norm_FST_burchellii",ncol=dim(norm_FST))

#plot(f_test$V1,cover_test$V1,col=2,main="Fst for all Loci by Number of Samples",xlab="Fst",ylab="Number of Individuals with Locus (Max: 75)",sub="Black = Null, Red = mexicanum")
#points(f_null$V1,cover_test$V1)

#plot(density(f_test$V1,na.rm=TRUE),col=2,xlim=c(0,1),ylim=c(0,30),main="Fst for all Loci: Density Distribution",xlab="Fst",ylab="Density",sub="Black = Null, Red = Test")
#lines((density(f_null$V1,na.rm=TRUE)))

###PLOTS WORKING WITH NUCLEOTIDE DIVERSITY AND FST
#Plots relationship between normalized Pi and Fst
plot(FST_data$Pi,FST_data$F_test,xlab="Pi (Normalized)",ylab="Fst",main="Nucleotide Diversity vs. Fst: E. burchellii",sub="N = 168,518")
plot(density(FST_data$Pi),xlab="Pi (Normalized)",main="Density Distribution for Nucleotide Diversity: E. burchellii",sub = "N = 168,518")
plot(density(FST_data$F_test,na.rm=TRUE),xlab="Fst",main="Density Distribution for Fst: E. burchellii",sub = "N = 168,518")

pi_fst_plot = ggplot(FST_data, aes(Pi, y = value))
pi_fst_plot + geom_point(alpha = 0.05, aes(y = F_null))
pi_fst_plot + geom_point(alpha = 0.05, aes(y = F_test))

#Plots Pi against Coverage of Individuals
plot(FST_data$Cover,FST_data$Pi)
