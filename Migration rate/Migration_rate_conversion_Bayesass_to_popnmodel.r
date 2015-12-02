#This script converts the Bayesass 'migration rates' to migration rates as defined by population modellers

inpMig <- "W:/seq_genetics/analysis/collaboration_paper/1_genetics/Bayesass_ClaireR/output/run100/migration_rates_best.csv"
apps <- list("bris", "loga", "redl")
inpLeslie <- "W:/seq_genetics/analysis/collaboration_paper/2_population models/Model_parameters/Growth_rates/output_lambdas/LeslieMatrix_"
outDir <- "W:/seq_genetics/analysis/collaboration_paper/2_population models/Model_parameters/Migration_rates/"

######################
#Load parameters
######################
#Migration rates (Bayesass)
migBayes <- as.matrix(read.csv(inpMig, header=T)[,2:4])

#Initial population size (2008)
#N1=Bris, N2=Logan, N3=Redlands
N2008 <- c(bris=327,loga=450, redl=1502)
N2005 <- c(bris=721,loga=951, redl=2939)
N1996 <- c(bris=906,loga=1286, redl=4053)

#Generation time (years)
genTime <- 6

#####################
#Calculate per capita birth rate and stable age structure
#####################

#Make a list of the leslie matrices
leslieL <- lapply(apps, function(x) as.matrix(read.csv(paste0(inpLeslie, x, ".txt"))))


#Find the age structured fecundity
FecL <- lapply(leslieL, function(x){
		#stable age distribution (the right eigenvector of the dominant eigenvalue of the Leslie matrix)
		currEigen <-eigen(x)
		stable_age <- matrix(Re(currEigen$vectors[,1]/sum(currEigen$vectors[,1])), ncol=1)
		#fecundity
		F1S1 = x[1,2]
		F2S2 = x[1,3]
		F3S3 = x[1,4]
		F <- matrix(c(0,F1S1,F2S2,F3S3), nrow=1)
		return(stable_age %*% F)
		})

#Calculate the age structured survival rates
SurvL <- lapply(leslieL, function(x){
		x[1,]<-0 #assume offspring don't breed
		St <- as.matrix(c(1,0,0,0), ncol=1) #start with offspring
		St6 <- for (i in 1:genTime){
			St <- x%*%St
			} #how many survive to year 6
			#return(colSums(St))
			return(St)
		})

#per capita birthrate per generation
B <- lapply(c(1:3), function(i) {
		y <-genTime*FecL[[i]] %*% SurvL[[i]]
		return(colSums(y))
		})


#####################
#Calculate the migration rates per generation
#####################

z <- matrix(0, nrow=3, ncol=3) 

for (i in 1:3){
	for (j in 1:3){
		z[i,j] <- (N2008[j] * migBayes[j,i]) / (N1996[i] * B[[j]] + N2005[i])
	}
}
z[1,1]<-NA
z[2,2]<-NA
z[3,3]<-NA

dimnames(z)<- list(unlist(apps), unlist(apps))
write.csv(z, paste0(outDir, "Converted_migration_rates_per_generation.csv"))

#####################
#Calculate the annual migration rates
#####################

zAnnual <- z/genTime

write.csv(zAnnual, paste0(outDir, "Converted_migration_rates_annual.csv"))

######################

 stableagefun <- function(i){
  x<-leslieL[[i]]
  currEigen <-eigen(x)
  stable_age <- matrix(Re(currEigen$vectors[,1]/sum(currEigen$vectors[,1])), ncol=1)
  print(stable_age)
  }
 stableagefun(1)
 stableagefun(2)
 stableagefun(3)
 