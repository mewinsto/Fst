#FST

library(hierfstat)

##IMPORTS AND FORMATS SNP MATRIX ("GENO_MATRIX") FROM .GENO FILE (pyRAD OUTPUT)
#GENO_MATRIX HAS SAMPLES AS ROWS AND SNPS FOR COLUMNS
setwd("~/Desktop/")
geno_data = scan("all_lanes_mex_10k.usnps.geno",what = "character")
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
samples = read.table("all_lanes_mex_10k.samples",header=FALSE)
samples = samples$V1

#READS IN BURCH POP TO MAKE A TEST VECTOR, THEN CONSTRUCTS NULL VECTOR BASED ON POPS
testpop = read.table("mexicanum_pop.txt",header=FALSE)
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
write(FST,file="Fst_test_mexicanum.txt",ncolumns=1)
write(geno_cover,file="geno_cover_mexicanum.txt",ncolumns=1)

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
write(FST_N,file="Fst_null_mexicanum.txt",ncolumns=1)

###PLOTTING 
#Plot densities of these values
plot(density(geno_means,bw=0.01),xlim=c(0,1),main="Distribution of Major Allele Frequencies across Loci", xlab="Major Allele Frequency", sub="Species = E. mexicanum -- N = 163,813 - bw = .01")
plot(density(geno_cover,bw=0.5),xlim=c(4,S),main="Density Distribution of Individuals per Locus: E. mexicanum", xlab="Number of Individuals (Max = 29)",sub="N = 163,813, bw = 0.5")

##PLOT FST AGAINST GENO_COVER
plot(density(FST_N,na.rm=TRUE,bw=.004),ylim=c(0,15),main="Density Distribution for Fst Values: E. mexicanum",xlab="Fst Values", sub = "N = 163,813   bw = 0.004")
lines(density(FST,na.rm=TRUE,bw=.004),col=2)

plot(FST,geno_cover,col=2,main="Fst for all Loci by Number of Samples",xlab="Fst",ylab="Number of Individuals with Locus (Max: 29)",sub="Black = Null, Red = mexicanum -- N = 163,813")
points(FST_N,geno_cover)

##READ IN FILES: NULL AND TEST AND PLOT AGAINST ONE ANOTHER
#f_null = read.table("Fst_null.txt")
#f_test = read.table("Fst_test.txt")
#cover_test = read.table("geno_cover.txt") 


#plot(f_test$V1,cover_test$V1,col=2,main="Fst for all Loci by Number of Samples",xlab="Fst",ylab="Number of Individuals with Locus (Max: 75)",sub="Black = Null, Red = mexicanum")
#points(f_null$V1,cover_test$V1)

#plot(density(f_test$V1,na.rm=TRUE),col=2,xlim=c(0,1),ylim=c(0,30),main="Fst for all Loci: Density Distribution",xlab="Fst",ylab="Density",sub="Black = Null, Red = Test")
#lines((density(f_null$V1,na.rm=TRUE)))
