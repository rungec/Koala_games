#KOALA SEQ GENETICS COLLABORATION PROJECT
#This script calculates the investment benefits and optimal investment spread
#Using theoretical values for koala dispersal in SEQ
#comparing different levels of assymetry

#setwd("C:/Users/uqgiacon/Dropbox/UQ Comp/Research/Game Theory Koalas/data")

#Setup inputs
inpMig <- "E:/Projects/seq_genetics/analysis/collaboration_paper/2_population models/Model_parameters/Migration_rates/Converted_migration_rates_annual.csv"
inpVitals <- "E:/Projects/seq_genetics/analysis/collaboration_paper/2_population models/Model_parameters/Growth_rates/output_vitalrates/"
inpLeslie <- "E:/Projects/seq_genetics/analysis/collaboration_paper/2_population models/Model_parameters/Growth_rates/output_lambdas/"
outDir <- "E:/Projects/seq_genetics/analysis/collaboration_paper/3_game theory/outputs/Random dispersal rates/Run 1/"
# inpMig <- "C:/Claire/GPEM_Postdoc/2_SEQ_Koala_Collaboration/2_population models/Model_parameters/Migration_rates/Converted_migration_rates_annual.csv"
# inpVitals <- "C:/Claire/GPEM_Postdoc/2_SEQ_Koala_Collaboration/2_population models/Model_parameters/Growth_rates/output_vitalrates/"
# inpLeslie <- "C:/Claire/GPEM_Postdoc/2_SEQ_Koala_Collaboration/2_population models/Model_parameters/Growth_rates/output_lambdas/"
# outDir <- "C:/Claire/GPEM_Postdoc/2_SEQ_Koala_Collaboration/3_game theory/outputs/Random dispersal Rates/"
apps <- list("bris", "loga", "redl")
library(plyr) 
library(dplyr) 
library(deSolve) # for ode solving 
library(nloptr)
library(snow)
library(doSNOW)
library(compiler) #compile code to speed it up

ncore <- 12 #Number of CPU for parallel

##################
#MODEL PARAMETERS
##################
#Sex ratio
Sex_ratio <- 0.5
#Generation time (years)
genTime <- 6
#Council investment ($m annual)
initInvest <- c(bris=40.6, loga=9.0, redl=10.8)
#Area of each LGA
area <- c(bris=8389/37377, loga=7728/37377, redl=21260/37377)

#Initial population size (2008), females only
#N1=Bris, N2=Logan, N3=Redlands
Popn2008 <- c(N1=327,N2=450, N3=1502)*Sex_ratio
Popn2005 <- c(N1=721,N2=951, N3=2939)*Sex_ratio
Popn1996 <- c(N1=906,N2=1286, N3=4053)*Sex_ratio
N1 <- Popn2008["N1"]
N2 <- Popn2008["N2"]
N3 <- Popn2008["N3"]

#50% decline over 3 generations
declper <- c(N1, N2, N3, N1+N2+N3)*0.5
names(declper) <- c(unlist(apps), "overall")

#Initial growth rates (2008) #females log(r) to convert lambda to r
r1 <- log(read.csv(paste0(inpLeslie, "LambdaF_bris.txt"))["Rates_ModAvFemales","mean"])
r2 <- log(read.csv(paste0(inpLeslie, "LambdaF_loga.txt"))["Rates_ModAvFemales","mean"])
r3 <- log(read.csv(paste0(inpLeslie, "LambdaF_redl.txt"))["Rates_ModAvFemales","mean"])

#Carrying capacity, females only
initK <- c(K1=805, K2=1098, K3=3855)*Sex_ratio
K1 <- initK["K1"]
K2 <- initK["K2"]
K3 <- initK["K3"]

#Annual migration rates (converted from Bayesass)
#Z <- as.matrix(read.csv(inpMig, header=T)[,2:4])

#Fertility and mortality rates
Vitals <- lapply(apps, function(i) {
		read.csv(paste0(inpVitals, "modav_fem", i, ".txt")) %>% colMeans()
		})

####################
Times <- seq(0,18,1) #number of years to model
#State <- seq(0,45,0.1) #budget to map
StateIncr <- 0.1
State <- seq(0,10,StateIncr) #budget to map

#Pars <- as.list(c(r1=r1,r2=r2, r3=r3, K1, K2, K3, Z=Z)) #other parameters
Yini <- c(N1, N2, N3) #the initial state

###################
#Dispersal rates Z
	#coefficient of asymmetry 0=asymmetric 1=symmetric #from Esposito et al 2014 Measuring Symmetry, Assymetry and Randomness in Neural Network Connectivity
	#generate numbers from beta distribution (0-50% migration over a generation) 
	#betadist <- sapply(1:6, function(x) rbeta(20000, 0.5,0.5)*50/600) #beta distribution #%of popn
	betadist <- matrix(sample(rbeta(60000, 0.5,0.5)*0.15/6,60000),ncol=6) #beta distribution #15%of popn moving over a generation
	#numKoalaDispMax <- rep(Popn2008, 2)*0.3/6 # 30% of popn moving over 6 years, as a proportion (ie divided by 100)
	numKoalaDisp <- betadist*sum(Popn2008)
	#betadist <- sapply(numKoalaDispMax, function(x) rbeta(20000, 0.5,0.5)*x) #beta distribution, number of koala in population
	#calculate the symmetric coef
	# symmcoef <- apply(betadist, 1, function(x) {
				# return(1-((2/6)*sum(c(abs(x[1]-x[4])/(x[1]+x[4]), abs(x[2]-x[5])/(x[2]+x[5]), abs(x[3]-x[6])/(x[3]+x[6])))))
				# })
	#symmcoef <- apply(numKoalaDisp, 1, function(x) {
			#	return(1-((2/6)*sum(c(abs(x[1]-x[2])/(x[1]+x[2]), abs(x[3]-x[4])/(x[3]+x[4]), abs(x[5]-x[6])/(x[5]+x[6])))))
			#	})	
	# betadist <- data.frame(betadist, symmcoef)
	# betadist$bin <- sapply(symmcoef, function(x) trunc(x*10))
	# bindata <- ddply(betadist, "bin", function(x) x[sample(1:nrow(x), 10),])
	symmcoef <- apply(numKoalaDisp, 1, function(x) {
				return(1-((2/6)*sum(c(abs(x[1]-x[2]), abs(x[3]-x[4]), abs(x[5]-x[6])))/(sum(Popn2008)*0.15/6))) #the difference between immigration and emigration between any two councils, divided by the maximum number of koalas moving in any direction (15% of total popn/6)
				})		
								
	#randomly sample migration rates that fit into equal bins of assymetry
	#10 bins, 10 samples from each
	numKoalaDisp <- data.frame(numKoalaDisp, symmcoef)
	numKoalaDisp$bin <- sapply(symmcoef, function(x) trunc(x*10))
	bindata <- ddply(numKoalaDisp, "bin", function(x) x[sample(1:nrow(x), 10),])
	#bindata <- as.matrix(bindata)
	bindata <- data.frame(bindata, sweep(bindata[,1:6],2, rep(Popn2008,2), '/')) #calculate the proportional dispersal rate
	names(bindata) <- c("N12","N21","N31","N13","N23","N32","symmcoef","bin","P12","P21","P31","P13","P23","P32")
	outPath <- paste0(outDir, "randomlygeneratedZvalues.csv")
	write.csv(bindata, outPath, row.names=FALSE)

####################
#Set up functions for population models
####################

#Continous model with migration
ThreePatchModC <- function(Time,State, Pars){
	with(as.list(c(State, Pars)),{
    	dN1 <- N1 * r1 * (1 - (N1 / K1)) - N1*(Z[1,2][[1]] + Z[1,3][[1]]) + (Z[2,1][[1]]*N2+Z[3,1][[1]]*N3)
    	dN2 <- N2 * r2 * (1 - (N2 / K2)) - N2*(Z[2,1][[1]] + Z[2,3][[1]]) + (Z[1,2][[1]]*N1+Z[3,2][[1]]*N3)
    	dN3 <- N3 * r3 * (1 - (N3 / K3)) - N3*(Z[3,1][[1]] + Z[3,2][[1]]) + (Z[1,3][[1]]*N1+Z[2,3][[1]]*N2)
    return(list(c(dN1,dN2,dN3)))
	})
}

#########################
#Set up function for return on investment
#########################

ReturnonInvest <- cmpfun(function(invest,i){
	#invest <- initInvest[as.character(apps[i])]
	#bris i=1 loga i=2 redl i=3
	#88 % of budget spent on cars, 12% of budget spent on dogs (Maxwell et al 2015); 
	#assuming impediments to management success (Ng et al 2014)
	#Define the cost function for dollars spent to reduce car mortality:
	prM_car = 0.23 + 0.77*exp((-invest*0.88/area[i])/5.93)
	# Define the cost function for dollars spent to reduce dog mortality:
	prM_dog = 0.56 + 0.44*exp((-invest*0.12/area[i])/0.81)

	#read in initial mortality rates
	currVitals <- Vitals[[i]]
	#The survival rate for koalas in age class 0 is unaffected by management: 
	S0 = 1 - currVitals["M0"]
	#The survival rates for koalas in age classes 1-3 is: 
	S1 = 1 - currVitals["MNAT1"] - currVitals["MDIS1"] - currVitals["MDOG1"]*prM_dog - currVitals["MCAR1"]*prM_car
	S2 = 1 - currVitals["MNAT2"] - currVitals["MDIS2"] - currVitals["MDOG2"]*prM_dog - currVitals["MCAR2"]*prM_car
	S3 = 1 - currVitals["MNAT3"] - currVitals["MDIS3"] - currVitals["MDOG3"]*prM_dog - currVitals["MCAR3"]*prM_car
	#Leslie matrix
	Leslie <- matrix(0,ncol=4,nrow=4)
	Leslie[1,2] <- currVitals["F1"] * S1
	Leslie[1,3] <- currVitals["F2"] * S2
	Leslie[1,4] <- currVitals["F3"] * S3
	Leslie[2,1] <- S0
	Leslie[3,2] <- S1
	Leslie[4,3] <- S2
	Leslie[4,4] <- S3
	
	rnew <- log(Re(eigen(Leslie)$values[1]))
	return(rnew)
})

######################
#Set up functions for scenarios 3 & 4 - Discrete optimisation
######################

#Create an array of the investments sets, value in each cell is the number of councils that meet 50% decline (0-3)
NashPossFun <- cmpfun(function(State){
	TargetDecline = 0.50*(Popn2008)
	#idealpop <- sum(Popn2008)*TargetDecline
	#Set up blank arrays
	NashPoss = array(NA,dim =c(length(State),length(State),length(State)))
	OverallBenefit = array(NA,dim =c(length(State),length(State),length(State)))
	#loop through arrays
	for (i in seq_along(State)){
		for (j in seq_along(State)){
			for (k in seq_along(State)){
				Pars <- as.list(c(r1=rMatrix[as.character(State[i]),1], r2=rMatrix[as.character(State[j]),2], r3=rMatrix[as.character(State[k]),3], K1, K2, K3, Z=Z))#update the Pars with the post-investment r values
				NewN <- ode(Yini,Times,ThreePatchModC,Pars)#recalculate the population size
				NashPoss[i,j,k] <- {NewN[length(Times),"N1"]>=TargetDecline["N1"]} +{NewN[length(Times),"N2"]>=TargetDecline["N2"]} + {NewN[length(Times),"N3"]>=TargetDecline["N3"]} 
				#If we don't need to keep the benefit arrays this is very mildly faster
				OverallBenefit[i,j,k] <- sum(NewN[nrow(NewN),2:4])/1139.5 #sum(Popn2008)
				}}}
				
			CostMatr <- which(OverallBenefit >=0.5, TRUE)
			if (length(CostMatr)==0){
			optimalcost <- c(NA, NA, NA)
			} else {
			optimalcost <- State[CostMatr[which.min(rowSums(CostMatr)),]]
			names(optimalcost)<-list("bris", "loga", "redl") #The optimal investment strategy to meet target of 50% overall
			}
	return(list(NashPoss=NashPoss, OverallOptim=optimalcost))
})

NashFun <-cmpfun(function(TargetMet,n){
#print(Sys.time())
Nash<-matrix(NA,nrow=nrow(TargetMet), byrow=TRUE)
for(i in 1:n){
	currSpend <- TargetMet[sample(nrow(TargetMet),1),] # random choice for start point 
	GoOrder <- sample(c(1:3)) # random order of council spend
	
	#find the best that the first player could do given the other two choices are held constant
	FirstChoices <- as.matrix(TargetMet[which(TargetMet[,GoOrder[2]]== currSpend[GoOrder[2]] & TargetMet[,GoOrder[3]]== currSpend[GoOrder[3]]),])
	if (dim(FirstChoices)[2]==1){ FirstChoices <- t(FirstChoices)} # to prevent format issues if there is only one option for the choice sets
	if (dim(FirstChoices)[1]==1){ newChoices <- FirstChoices } else { # if there is only one option for choice sets (if there is only one row call first choices new choices_?)
	newChoices <- FirstChoices[order(FirstChoices[,1]),][1,]}

	#find the best that the second player could do given the other choices are held constant
	SecondChoices <- as.matrix(TargetMet[which(TargetMet[,GoOrder[1]]== newChoices[GoOrder[1]] &TargetMet[,GoOrder[3]]== newChoices[GoOrder[3]]),])
	if (dim(SecondChoices)[2]==1){ SecondChoices <- t(SecondChoices)} 
	if (dim(SecondChoices)[1] ==1){ newChoices <- SecondChoices } else {
	newChoices <- SecondChoices[order(SecondChoices[,2]),][1,]	}

	#find the best that the third player could do given the other choices are held constant
	ThirdChoices <- as.matrix(TargetMet[which(TargetMet[,GoOrder[1]]== newChoices[GoOrder[1]] & TargetMet[,GoOrder[2]]== newChoices[GoOrder[2]]),])
	if (dim(ThirdChoices)[2]==1){ ThirdChoices <- t(ThirdChoices)}
	if (dim(ThirdChoices)[1]==1){ newChoices <- ThirdChoices } else {
	newChoices <- ThirdChoices[order(ThirdChoices[,3]),][1,]	}
	
	PrintVals <- rbind(PrintVals, c(currSpend, newChoices[1], newChoices[2], newChoices[3]))
}

return(PrintVals)
})


###########################
#Calculate r-values post investment
#in parallel
###########################
rMatrix <- matrix(nrow=length(State), ncol=3, dimnames=list(State, c("r1", "r2", "r3")))
councils <- list(1,2,3)

clust <- makeCluster(ncore, type="SOCK")
#export variables/functions
clusterEvalQ(clust, {library(plyr);library(dplyr);library(deSolve);library(nloptr);library(compiler)})
clusterExport(clust, c("area", "Vitals", "State", "ReturnonInvest", "rMatrix", "councils"))
registerDoSNOW(clust)

#Set up a matrix where each row represents r for each investment value
rMatrix[,1:3] <- t(laply(councils, .parallel=TRUE, function(i){
	ret=laply(State, function(invest) ReturnonInvest(invest,i))
	return(ret)
	}))

clust=stopCluster(clust)

outPath <-paste0(outDir, "GrowthRate_r_vsInvestment_by_council_10mby100k.csv")
write.csv(rMatrix, outPath)



######################
#Do in Parallel
######################

clust <- makeCluster(ncore, type="SOCK")
#export variables/functions
clusterEvalQ(clust, {library(plyr);library(dplyr);library(deSolve);library(nloptr);library(compiler)})

clusterExport(clust, c("area", "Popn2008", "N1", "N2", "N3", "r1", "r2", "r3", "K1", "K2", "K3", "Vitals", "Times", "StateIncr", "State", "Yini", "ThreePatchModC", "ReturnonInvest", "bindata", "NashPossFun","NashFun", "rMatrix"))
registerDoSNOW(clust)

#seq_len(nrow(bindata[1:5,]) for first 5 rows for testing
assymROI <- llply(seq_len(nrow(bindata)), .parallel=TRUE, function(currRow) {
	#new Z (assign to global so all functions can see it)
	# Z<<-matrix(c(NA, bindata[currRow,1], bindata[currRow,2], bindata[currRow,4], NA, bindata[currRow,3], bindata[currRow,5], bindata[currRow,6], NA), ncol=3, dimnames=NULL)
	#Z<<-matrix(c(NA, bindata[currRow,2], bindata[currRow,3], bindata[currRow,1], NA, bindata[currRow,6], bindata[currRow,4], bindata[currRow,5], NA), ncol=3, dimnames=NULL)
	Z<<-matrix(c(NA, bindata[currRow,"P21"], bindata[currRow,"P31"], bindata[currRow,"P12"], NA, bindata[currRow,"P32"], bindata[currRow,"P13"], bindata[currRow,"P23"], NA), ncol=3, dimnames=NULL)
	
	#Calculate best investment strategy for 50% overll decline, plus draw up the array where value = c(0:3) representing number of councils meeting 50% local decline
	NashList <- NashPossFun(State)
	
	#Table of all investments strategies meeting 50% local decline
	currNashPoss <<- NashList$NashPoss
	currTargetMet <<- which(currNashPoss == max(currNashPoss),TRUE)
	maxN <<- max(currNashPoss)
	dimnames(currTargetMet)<-list(NULL, c("bris", "loga", "redl"))

	#Sets of unstable Nash investment strategies 50% local
	PrintVals <<- matrix(NA,1,6)
	PrintVals <<- NashFun(currTargetMet, 10000)
	# the end sets of spend choices 
	EquilibriumSpendOptions <- aaply(na.omit(unique(PrintVals[,4:6])),c(1,2),function(x) State[x])
	#frequency of those sets of choices
	freq<-table(as.numeric(paste(na.omit(PrintVals[,4]),na.omit(PrintVals[,5]),na.omit(PrintVals[,6]), sep = ""
	)))
	if(length(EquilibriumSpendOptions)<=3){ 
		EquilibriumSpendOptions <-append(EquilibriumSpendOptions, freq)
	} else {EquilibriumSpendOptions <- cbind(EquilibriumSpendOptions, freq)
	}
	return(list(LocalOptim=EquilibriumSpendOptions, OverallOptim=NashList$OverallOptim, councilsmettargets=maxN, currZ=Z))
})
 
clust=stopCluster(clust) 

#if you forget to stopCluster run this in cmd 
#taskkill -im Rscript.exe -f

#Save data
OverallInvStrat<- laply(seq_len(nrow(bindata)), function(x) assymROI[[x]]$OverallOptim)
OverallInvStrat <- cbind(OverallInvStrat, bindata[,"symmcoef"], seq_len(nrow(bindata)))
dimnames(OverallInvStrat) <- list(NULL, c("bris", "loga", "redl", "symmcoef", "run"))
outPath <- paste0(outDir, "InvestmentStrategies_Overall50decline_100randomZ.csv")
write.csv(OverallInvStrat, outPath)

#Save data
LocalInvStrat <- ldply(seq_len(nrow(bindata)), function(x) {					
				currROI <-matrix(assymROI[[x]]$LocalOptim, ncol=4)
				currMaxN <- assymROI[[x]]$councilsmettargets
				return(cbind(currROI, rep(bindata[x,"symmcoef"], nrow(currROI)), rep(x, nrow(currROI)), rep(currMaxN, nrow(currROI))))
				})
names(LocalInvStrat) <- c("bris", "loga", "redl", "freq", "symmcoef", "run", "councilsmettargets")
outPath <- paste0(outDir, "InvestmentStrategies_Local50decline_100randomZ.csv")
write.csv(LocalInvStrat, outPath)

#Calculate how many koalas each budget choice results in (Overall scenario)
ROIs <- ldply(seq_len(nrow(bindata)), function(i) {
			newr <<- Map(ReturnonInvest, assymROI[[i]]$OverallOptim, 1:3)
			names(newr)<-NULL
			Z<<-assymROI[[i]]$currZ
			Pars <<- as.list(c(r1=newr[1], r2=newr[2], r3=newr[3], K1, K2, K3, Z=Z))
			NewN <- ode(Yini,Times,ThreePatchModC,Pars)
			numKoala <- NewN[nrow(NewN),2:4]
			declKoala <- 1-numKoala/Popn2008
			numKoalattl <- sum(numKoala)
			declKoalattl <- 1-numKoalattl/sum(Popn2008)
			return(append(assymROI[[i]]$OverallOptim, c(numKoala, numKoalattl, declKoala, declKoalattl)))
			})
names(ROIs) <- c("bris_invest", "loga_invest", "redl_invest", "numKoala_bris", "numKoala_loga", "numKoala_redl", "numKoala_overall", "decline_bris", "decline_loga", "decline_redl", "decline_overall")
outPath <- paste0(outDir, "ROI_overall50decline_100randomZ.csv")
write.csv(ROIs, outPath)

#Calculate how many koalas each budget choice results in (maximum possible investment)
ROIMAXs <- ldply(seq_len(nrow(bindata)), function(i) {
			print(i)
			#print(assymROI[[i]]$OverallOptim)
			newr <- Map(ReturnonInvest, c(10,10,10), 1:3)
			names(newr)<-NULL
			Z<<-assymROI[[i]]$currZ
			print(Z)
			Pars <<- as.list(c(r1=newr[1], r2=newr[2], r3=newr[3], K1, K2, K3, Z=Z))
			print(Pars)
			NewN <- ode(Yini,Times,ThreePatchModC,Pars)
			numKoala <- NewN[nrow(NewN),2:4]
			declKoala <- 1-numKoala/Popn2008
			numKoalattl <- sum(numKoala)
			declKoalattl <- 1-numKoalattl/sum(Popn2008)
			return(append(c(10,10,10), c(numKoala, numKoalattl, declKoala, declKoalattl)))
			})
names(ROIMAXs) <- c("bris_invest", "loga_invest", "redl_invest", "numKoala_bris", "numKoala_loga", "numKoala_redl", "numKoala_overall", "decline_bris", "decline_loga", "decline_redl", "decline_overall")
outPath <- paste0(outDir, "ROI_maxspend_100randomZ.csv")
write.csv(ROIMAXs, outPath)

ROILocal <- ldply(seq_len(nrow(LocalInvStrat)), function(i) {
			currBudgetLcl <-as.numeric(LocalInvStrat[i,c("bris", "loga", "redl")])
			newr <- Map(ReturnonInvest, currBudgetLcl, 1:3)
			names(newr)<-NULL
			currRun <- LocalInvStrat[i,"run"]
			Z<<-matrix(c(NA, bindata[currRun,1], bindata[currRun,2], bindata[currRun,4], NA, bindata[currRun,3], bindata[currRun,5], bindata[currRun,6], NA), ncol=3, dimnames=NULL)
			#print(Z)
			Pars <<- as.list(c(r1=newr[1], r2=newr[2], r3=newr[3], K1, K2, K3, Z=Z))
			#print(Pars)
			NewN <- ode(Yini,Times,ThreePatchModC,Pars)
			numKoala <- NewN[nrow(NewN),2:4]
			declKoala <- 1-numKoala/Popn2008
			numKoalattl <- sum(numKoala)
			declKoalattl <- 1-numKoalattl/sum(Popn2008)
			return(append(currBudgetLcl, c(numKoala, numKoalattl, declKoala, declKoalattl, currRun)))
			})
names(ROILocal) <- c("bris_invest", "loga_invest", "redl_invest", "numKoala_bris", "numKoala_loga", "numKoala_redl", "numKoala_overall", "decline_bris", "decline_loga", "decline_redl", "decline_overall", "run")
outPath <- paste0(outDir, "ROI_local50decline_100randomZ.csv")
write.csv(ROILocal, outPath)


#the end
##################
##################
