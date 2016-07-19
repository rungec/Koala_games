#KOALA SEQ GENETICS COLLABORATION PROJECT
#This script calculates the investment benefits and optimal investment spread under 5 scenarios 
#Using actual values of koala dispersal in SEQ
# see script Game_Investment_benefit_eval_dispersal_assymetry_150907.r for parallel implementation and calculation of investment under randomly generated migration rates

#setwd("C:/Users/uqgiacon/Dropbox/UQ Comp/Research/Game Theory Koalas/data")

#Setup inputs
#inpMig <- "C:/Claire/GPEM_Postdoc/2_SEQ_Koala_Collaboration/2_population models/Model_parameters/Migration_rates/Converted_migration_rates_annual.csv"
 inpMig <- "W:/seq_genetics/analysis/collaboration_paper/2_population models/Model_parameters/Migration_rates/Converted_migration_rates_annual.csv"
apps <- list("bris", "loga", "redl")
 inpVitals <- "W:/seq_genetics/analysis/collaboration_paper/2_population models/Model_parameters/Growth_rates/output_vitalrates/"
#inpVitals <- "C:/Claire/GPEM_Postdoc/2_SEQ_Koala_Collaboration/2_population models/Model_parameters/Growth_rates/output_vitalrates/"
 inpLeslie <- "W:/seq_genetics/analysis/collaboration_paper/2_population models/Model_parameters/Growth_rates/output_lambdas/"
#inpLeslie <- "C:/Claire/GPEM_Postdoc/2_SEQ_Koala_Collaboration/2_population models/Model_parameters/Growth_rates/output_lambdas/"
outDir <- "W:/seq_genetics/analysis/collaboration_paper/3_game theory/outputs/"
#outDir <- "C:/Claire/GPEM_Postdoc/2_SEQ_Koala_Collaboration/3_game theory/outputs/"

#setwd(outDir)

library(plyr) 
library(dplyr) 
library(deSolve) # for ode solving 
#library(microbenchmark)
library(nloptr)
#library(snow)
library(compiler)


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

#the best possible decline for this system (28.4% decline)
best <- 1-0.284

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
Z <- as.matrix(read.csv(inpMig, header=T)[,2:4])

#Fertility and mortality rates
Vitals <- lapply(apps, function(i) {
		read.csv(paste0(inpVitals, "modav_fem", i, ".txt")) %>% colMeans()
		})

####################
Times <- seq(0,18,1) #number of years to model
State <- seq(0,10,0.1) #budget to map
#StateIncr <- 0.001
#State <- seq(0,45,StateIncr) #budget to map

Pars <- as.list(c(r1=r1,r2=r2, r3=r3, K1, K2, K3, Z=Z)) #other parameters
Yini <- c(N1, N2, N3) #the initial state

####################
#Set up functions for population models
####################

#Continous model with migration
ThreePatchModC <- function(Time,State, Pars){
	with(as.list(c(State, Pars)),{
    	dN1 <- N1 * r1 * (1 - (N1 / K1)) - N1*(Z[1,2] + Z[1,3]) + (Z[2,1]*N2+Z[3,1]*N3)
    	dN2 <- N2 * r2 * (1 - (N2 / K2)) - N2*(Z[2,1] + Z[2,3]) + (Z[1,2]*N1+Z[3,2]*N3)
    	dN3 <- N3 * r3 * (1 - (N3 / K3)) - N3*(Z[3,1] + Z[3,2]) + (Z[1,3]*N1+Z[2,3]*N2)
    return(list(c(dN1,dN2,dN3)))
	})
}

#Continuous model ignoring migration
ThreePatchModCNoMig <- function(Time,State,Pars){
	with(as.list(c(State,Pars)),{
    	dN1 <- N1 * r1 * (1 - (N1 / K1)) #- N1*(Z[1,2] + Z[1,3]) + (Z[2,1]*N2+Z[3,1]*N3)
    	dN2 <- N2 * r2 * (1 - (N2 / K2)) #- N2*(Z[2,1] + Z[2,3]) + (Z[1,2]*N1+Z[3,2]*N3)
    	dN3 <- N3 * r3 * (1 - (N3 / K3)) #- N3*(Z[3,1] + Z[3,2]) + (Z[1,3]*N1+Z[2,3]*N2)
    return(list(c(dN1,dN2,dN3)))
	})
}


##############################
#Discrete model with migration
ThreePatchModD <- cmpfun(function(Times, Pars, Yini){
	popnMatrix <- matrix(ncol=4, nrow=length(Times))#set up a matrix
	popnMatrix[,1]<-Times #fill the first column with time
	popnMatrix[1,2:4] <- Yini #fill the first row with the starting population
	with(as.list(c(Yini,Pars)),{ #this line just tells it where to get the z and r values
	for (i in 2:length(Times)){
		#first calculate population size at t1 after immigration/emigration
			N1m <- N1 - N1*(Z[1,2] + Z[1,3]) + (Z[2,1]*N2+Z[3,1]*N3)
			N2m <- N2 - N2*(Z[2,1] + Z[2,3]) + (Z[1,2]*N1+Z[3,2]*N3)
			N3m <- N3 - N3*(Z[3,1] + Z[3,2]) + (Z[1,3]*N1+Z[2,3]*N2)
		#then calculate population size in t+1
			N1 <- N1m + N1m * r1 * (1 - (N1m / K1)) 
			N2 <- N2m + N2m * r2 * (1 - (N2m / K2))
			N3 <- N3m + N3m * r3 * (1 - (N3m / K3))
		popnMatrix[i,2:4]<-c(N1,N2,N3)
	}
	return(popnMatrix)
	})
})

#Discrete model no migration
ThreePatchModDNoMig <- cmpfun(function(Times, Pars, Yini){
	popnMatrix <- matrix(ncol=4, nrow=length(Times))#set up a matrix
	popnMatrix[,1]<-Times #fill the first column with time
	popnMatrix[1,2:4] <- Yini #fill the first row with the starting population
	with(as.list(c(Yini,Pars)),{ #this line just tells it where to get the z and r values
	for (i in 2:length(Times)){
		#then calculate population size in t+1
			N1 <- N1 + N1 * r1 * (1 - (N1 / K1)) 
			N2 <- N2 + N2 * r2 * (1 - (N2 / K2))
			N3 <- N3 + N3 * r3 * (1 - (N3 / K3))
		popnMatrix[i,2:4]<-c(N1,N2,N3)
	}
	return(popnMatrix)
	})
})

####To run models
#NewN returns a matrix of the number of koalas in each LGA in each year
#NewN <- ode(Yini,Times,ThreePatchModC,Pars)

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

###########################
#Calculate koala numbers post investment
###########################

################
#Scenario 1: Baseline
################
#no dispersal, investment resulting in 50% local decline
#Set up a 2D matrix where each row represents r for each investment value
	rMatrix <- matrix(nrow=length(State), ncol=3, dimnames=list(State, c("r1", "r2", "r3")))
	councils <- list(1,2,3)
	rMatrix[,1:3] <- sapply(councils, function(i){sapply(State, function(invest) ReturnonInvest(invest,i))})
	
	#outPath <-paste0(outDir, "GrowthRate_r_vsInvestment_by_council_to30m.csv")
	#write.csv(rMatrix, outPath)

#loop ThreePatchMod over each row of this matrix to get the benefit arrays
	BenefitS1 <- matrix(ncol=5, nrow=length(State), dimnames=list(State, c("State", unlist(apps), "overall"))) #2D
	BenefitS1[,"State"] <- State
	for (i in seq_along(State)){
		Pars <- as.list(c(r1=rMatrix[as.character(State[i]),1], r2=rMatrix[as.character(State[i]),2], r3=rMatrix[as.character(State[i]),3], K1, K2, K3))#update the Pars with the post-investment r values
		NewN <- ode(Yini,Times,ThreePatchModCNoMig,Pars)#recalculate the population size
	#the benefit is the population at the end of the time period
		BenefitS1[i,c("bris", "loga", "redl")] <- NewN[length(Times),c("N1", "N2", "N3")]
	}
	BenefitS1[,"overall"] <- rowSums(BenefitS1[,c("bris", "loga", "redl")])
	BenefitS1 <- data.frame(BenefitS1)
	#convert benefits from number of koalas to % decline
	BenefitS1[,c("brisD", "logaD", "redlD")] <- ddply(BenefitS1, 1, function(x) 1-(x[c("bris","loga", "redl")]/Popn2008[1:3]))[,2:4]
	BenefitS1[,"overallD"] <- 1-BenefitS1[,"overall"]/sum(Popn2008)

#############
#Optimal spend for local 50% decline	
	S1local <- sapply(apps, function(x){ 
					currDecline <- declper[x[1]]
					return(BenefitS1[min(which(BenefitS1[,x[1]] >= currDecline)), "State"])
	})

#############
#Optimal spend for overall 50% decline	
	nreps <- 10
	idealpop <- sum(Popn2008)*0.5 #50% decline
	#invest0 <-c(10,10,10)#fixed starting values for search
	#invest0 <-runif(3,0,45)#random starting values for search
	investList <- sapply(c(1:nreps), function(i) list(runif(3,min(State),max(State)))) #a list of random start values for cobyla fn
	fn<- function(invest) sum(invest) #thing to minimise
	#subject to inequality constraint hin(x)>=0 
	hin <-function(invest) {
			investReturn <- Map(ReturnonInvest, invest, 1:3) # the new growth rate due to investment
			Pars2 <- as.list(c(r1= investReturn[1],r2=investReturn[2], r3=investReturn[3], K1=K1,K2=K2, K3=K3))
			NewN <- ode(Yini,Times,ThreePatchModCNoMig,Pars2)
			return(sum(NewN[nrow(NewN),2:4])-idealpop) 
			}
	lower <-c(min(State),min(State),min(State)) 
	upper <-c(max(State),max(State),max(State))
	
	S1overall <- ldply(investList, function(i) {
				S1optimCOBYLA <- cobyla(i, fn, lower=lower, upper=upper, hin = hin, control = list(xtol_rel=1e-6, maxeval=500))
				#S1optimCOBYLA <- cobyla(i, fn, lower=lower, upper=upper, hin = hin, control = list(xtol_abs=rep(StateIncr, 3), maxeval=500)) #stops when increments by less than the bins in State
				return(t(unlist(S1optimCOBYLA[1:4])))
				}) #input list return df

	names(S1overall)[1:4] <- c("bris_invest", "loga_invest", "redl_invest", "overall_invest")

#############	
#save
outPath <- paste0(outDir, "InvestmentBenefit_Scenario1.csv")
write.csv(BenefitS1, outPath, row.names=FALSE)
outPath <- paste0(outDir, "OptimalInvestment_Scenario1_local50decline.csv")
write.csv(S1local, outPath, row.names=FALSE)
outPath <- paste0(outDir, "OptimalInvestment_Scenario1_overall50decline.csv")
write.csv(S1overall, outPath, row.names=FALSE)

##################
#Benefits from investing only in one council (and 0 in other councils) (with migration)
##################
for (j in 1:3){
		BenefitB <- matrix(ncol=5, nrow=length(State), dimnames=list(State, c("State", unlist(apps), "overall"))) #2D
		BenefitB[,"State"] <- State
		r1=r1<-rMatrix[as.character(State[1]),1]
		r2=r1<-rMatrix[as.character(State[1]),2]
		r3=r1<-rMatrix[as.character(State[1]),3]
	for (i in seq_along(State)){
		if (j==1) {
		r1<-rMatrix[as.character(State[i]),1]
		} else if (j==2) {
		r2<-rMatrix[as.character(State[i]),2]
		} else { 
		r3<-rMatrix[as.character(State[i]),3]
		}

		Pars <- as.list(c(r1=r1, r2=r2, r3=r3, K1, K2, K3))#update the Pars with the post-investment r values
		NewN <- ode(Yini,Times,ThreePatchModC,Pars)#recalculate the population size
	#the benefit is the population at the end of the time period
		BenefitB[i,c("bris", "loga", "redl")] <- NewN[length(Times),c("N1", "N2", "N3")]
	}
	BenefitB[,"overall"] <- rowSums(BenefitB[,c("bris", "loga", "redl")])
	BenefitB <- data.frame(BenefitB)
	#convert benefits from number of koalas to % decline
	BenefitB[,c("brisD", "logaD", "redlD")] <- ddply(BenefitB, 1, function(x) 1-(x[c("bris","loga", "redl")]/Popn2008[1:3]))[,2:4]
	BenefitB[,"overallD"] <- 1-BenefitB[,"overall"]/sum(Popn2008)
	
	outPath <- paste0(outDir, "InvestmentBenefit_", apps[[j]], "_only_Mig.csv")
	write.csv(BenefitB, outPath, row.names=FALSE)
	}


################
#Scenario 2: Baseline dispersal
################
##Investment benefit 2D matrix - calculates return but ASSUMES equal spend across each LGA
	BenefitS2 <- matrix(ncol=5, nrow=length(State), dimnames=list(State, c("State", unlist(apps), "overall"))) #2D
	BenefitS2[,"State"] <- State
	for (i in seq_along(State)){
		Pars <- as.list(c(r1=rMatrix[as.character(State[i]),1], r2=rMatrix[as.character(State[i]),2], r3=rMatrix[as.character(State[i]),3], K1, K2, K3, Z))#update the Pars with the post-investment r values
		NewN <- ode(Yini,Times,ThreePatchModC,Pars)#recalculate the population size
	#the benefit is the population at the end of the time period
		BenefitS2[i,c("bris", "loga", "redl")] <- NewN[length(Times),c("N1", "N2", "N3")]
	}
	BenefitS2[,"overall"] <- rowSums(BenefitS2[,c("bris", "loga", "redl")])
	BenefitS2 <- data.frame(BenefitS2)
	#convert benefits from number of koalas to % decline
	BenefitS2[,c("brisD", "logaD", "redlD")] <- ddply(BenefitS2, 1, function(x) 1-(x[c("bris","loga", "redl")]/Popn2008[1:3]))[,2:4]
	BenefitS2[,"overallD"] <- 1-BenefitS2[,"overall"]/sum(Popn2008)
	
#save
outPath <- paste0(outDir, "InvestmentBenefit_Scenario2.csv")
write.csv(BenefitS2, outPath, row.names=FALSE)

################
#Scenario 4: Benevolent social planner
################
#this calculates how much you would have to spend in each LGA to get an overall decline of no more than 50%
idealpop <- sum(Popn2008)*0.5 #50% decline

###SET PARAMETERS
#set up vars for cobyla -- Constrained Optimization BY Linear Approximations
nreps <- 10
#invest0 <-c(10,10,10)#fixed starting values for search
#invest0 <-runif(3,0,45)#random starting values for search
investList <- sapply(c(1:nreps), function(i) list(runif(3,min(State),max(State)))) #a list of random start values for cobyla fn
fn<- function(invest) sum(invest) #thing to minimise
#subject to inequality constraint hin(x)>=0 
hin <-function(invest) {
		investReturn <- Map(ReturnonInvest, invest, 1:3) # the new growth rate due to investment
		Pars2 <- as.list(c(r1= investReturn[1],r2=investReturn[2], r3=investReturn[3], K1=K1,K2=K2, K3=K3, Z=Z))
		NewN <- ode(Yini,Times,ThreePatchModC,Pars2)
		return(sum(NewN[nrow(NewN),2:4])-idealpop) 
		}
	lower <-c(min(State),min(State),min(State)) 
	upper <-c(max(State),max(State),max(State))

###OPTIMISE
S4optimDF <- ldply(investList, function(i) {
				S4optimCOBYLA <- cobyla(i, fn, lower=lower, upper=upper, hin = hin, control = list(xtol_rel=1e-6, maxeval=500))
				#S4optimCOBYLA <- cobyla(i, fn, lower=lower, upper=upper, hin = hin, control = list(xtol_abs=rep(StateIncr, 3), maxeval=500)) #stops when increments by less than the bins in State
				return(t(unlist(S4optimCOBYLA[1:4])))
				}) #input list return df

names(S4optimDF)[1:4] <- c("bris_invest", "loga_invest", "redl_invest", "overall_invest")
				
#Save
outPath <- paste0(outDir, "OptimalInvestment_Scenario4.csv")
write.csv(S4optimDF, outPath, row.names=FALSE)

################
#Scenario 4b: Benevolent social planner
################
#This calculates how much you would need to spend to get the best possible outcome (overall outcome)
idealpop <- sum(Popn2008)*best #best possible decline
#other parameters are the same as scenario 4 above
###OPTIMISE
S4BoptimDF <- ldply(investList, function(i) {
				S4optimCOBYLA <- cobyla(i, fn, lower=lower, upper=upper, hin = hin, control = list(xtol_rel=1e-6, maxeval=500))
				return(t(unlist(S4optimCOBYLA[1:4])))
				}) #input list return df

names(S4BoptimDF)[1:4] <- c("bris_invest", "loga_invest", "redl_invest", "overall_invest")
				
#Save
outPath <- paste0(outDir, "OptimalInvestment_Scenario4b.csv")
write.csv(S4BoptimDF, outPath, row.names=FALSE)

###CHECK solution
# solncheck <- function (s4optim) {
			# newr <- sapply(councils, function(i) ReturnonInvest(s4optim$par[i],i))
			# Pars <- as.list(c(r1=newr[1], r2=newr[2], r3=newr[3], K1, K2, K3, Z))
			# NewN <- ode(Yini,Times,ThreePatchModC,Pars)
			# return(sum(NewN[nrow(NewN),2:4])-idealpop)
			# }
			
#############
#SCENARIO 2 + ALL - return on investment across scenarios
#############
#for the same budget as scenario 1, how many koala do they get?
#for the min and max possible budget how many koala do they get?
#for the optimal spend (S4) how many koalas do they get? 

	S1overallbudget <- as.numeric(sapply(S1overall[which(S1overall$convergence==4),1:3], median))
	S4optimbudget <- as.numeric(sapply(S4optimDF[which(S4optimDF$convergence==4),1:3], median))
	S4Boptimbudget <- as.numeric(sapply(S4BoptimDF[which(S4BoptimDF$convergence==4),1:3], median))
	S2spend <- list(rep(min(State), 3), S1local, S1overallbudget, S4optimbudget, S4Boptimbudget, rep(max(State), 3))
		
	S2benefit <- function (budget, Z) {
			newr <<- Map(ReturnonInvest, budget, 1:3) # the new growth rate due to investment
			Pars <<- as.list(c(r1=newr[1], r2=newr[2], r3=newr[3], K1, K2, K3, Z))
			NewN <- ode(Yini,Times,ThreePatchModC,Pars)
			numKoala <- NewN[nrow(NewN),2:4]
			declKoala <- 1-numKoala/Popn2008
			numKoalattl <- sum(numKoala)
			declKoalattl <- 1-numKoalattl/sum(Popn2008)
			return(append(budget, c(numKoala, numKoalattl, declKoala, declKoalattl)))
	}

	S2benDF <- ldply(S2spend, S2benefit)
	names(S2benDF) <- c("bris_invest", "loga_invest", "redl_invest", "numKoala_bris", "numKoala_loga", "numKoala_redl", "numKoala_overall", "decline_bris", "decline_loga", "decline_redl", "decline_overall")
	
outPath <- paste0(outDir, "ReturnonInvestment_Benefits_AllScenarios.csv")
write.csv(S2benDF, outPath, row.names=FALSE)		


################
#Scenario 3: Non-collaborative simple game & discrete calculation of Scenario 4
################
#Set up functions
NashPossFun <- cmpfun(function(State){
TargetDecline = 0.50*Popn2008
NashPoss = array(NA,dim =c(length(State),length(State),length(State)))
OverallBenefit = array(NA,dim =c(length(State),length(State),length(State)))
for (i in seq_along(State)){
	for (j in seq_along(State)){
		for (k in seq_along(State)){
			Pars <- as.list(c(r1=rMatrix[as.character(State[i]),1], r2=rMatrix[as.character(State[j]),2], r3=rMatrix[as.character(State[k]),3], K1, K2, K3, Z=Z))#update the Pars with the post-investment r values
			NewN <- ode(Yini,Times,ThreePatchModC,Pars)#recalculate the population size
#the benefit is the percent overall decline in population by the time period
			NashPoss[i,j,k] <- {NewN[length(Times),"N1"]>=TargetDecline["N1"]} + {NewN[length(Times),"N2"]>=TargetDecline["N2"]} + {NewN[length(Times),"N3"]>=TargetDecline["N3"]}
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
	if (dim(SecondChoices)[1]==1){ newChoices <- SecondChoices } else {
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

#Run functions
NashPoss <- NashPossFun(State)
TargetMet <- which(NashPoss$NashPoss == max(NashPoss$NashPoss),TRUE)
print(max(NashPoss$NashPoss))
dimnames(TargetMet)<-list(NULL, c("bris", "loga", "redl"))

PrintVals <- matrix(NA,1,6)
PrintVals<-NashFun(TargetMet, 10000)

outPath <- paste0(outDir, "PrintVals_10000subsamples.csv")
write.csv(PrintVals, outPath)

# the end set of spend choices 
EquilibriumSpendOptions <- aaply(na.omit(unique(PrintVals[,4:6])),c(1,2),function(x) State[x])
outPath <- paste0(outDir, "Nash_10000subsamples.csv")
write.csv(EquilibriumSpendOptions, outPath)

# hack to get a frequency distribution of the different choices.
table(as.numeric(paste(na.omit(PrintVals[,4]),na.omit(PrintVals[,5]),na.omit(PrintVals[,6]), sep = "")))

# have to look at the choice index values to match it back up to the table of equilibrium spend options
unique(PrintVals[,4:6])

