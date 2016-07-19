### edited Claire's code to add a global search for the 
## 23/6/15

#setwd("C:/Users/uqgiacon/Dropbox/UQ Comp/Research/Game Theory Koalas/data")
rm(list = ls())

library(plyr) # don't know if this is actually used for anything?
library(dplyr) # don't know if this is actually used for anything?
library(deSolve) # for ode solving 

#Setup inputs
inpMig <- "C:/Users/uqgiacon/Dropbox/UQ Comp/Research/Game Theory Koalas/data/Converted_migration_rates_annual.csv"
apps <- list("bris", "loga", "redl")
inpVitals <- "C:/Users/uqgiacon/Dropbox/UQ Comp/Research/Game Theory Koalas/data/output_vitalrates/"
inpLeslie <- "C:/Users/uqgiacon/Dropbox/UQ Comp/Research/Game Theory Koalas/data/output_lambdas/"


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

####################
Times <- seq(0,18,1) #number of years to model
State <- seq(0,10,0.1)#c(1,10,45) #seq(0,45,0.1) #budget to map
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


#########################
#Set up function for return on investment
#########################

ReturnonInvest <- function(invest,i){
	#invest <- initInvest[as.character(apps[i])]
	#bris i=1 loga i=2 redl i=3
	#88 % of budget spent on cars, 12% of budget spent on dogs (Maxwell et al 2015); 
	#assuming impediments to management success (Ng et al 2014)
	#Define the cost function for dollars spent to reduce car mortality:
	prM_car = 0.23 + 0.77*exp((-invest*0.88/area[i])/5.93)
	# Define the cost function for dollars spent to reduce dog mortality:
	prM_dog = 0.56 + 0.44*exp((-invest*0.12/area[i])/0.81)

	#read in initial mortality rates
	Vitals <- read.csv(paste0(inpVitals, "modav_fem", apps[i], ".txt")) %>% colMeans()
	#The survival rate for koalas in age class 0 is unaffected by management: 
	S0 = 1 - Vitals["M0"]
	#The survival rates for koalas in age classes 1-3 is: 
	S1 = 1 - Vitals["MNAT1"] - Vitals["MDIS1"] - Vitals["MDOG1"]*prM_dog - Vitals["MCAR1"]*prM_car
	S2 = 1 - Vitals["MNAT2"] - Vitals["MDIS2"] - Vitals["MDOG2"]*prM_dog - Vitals["MCAR2"]*prM_car
	S3 = 1 - Vitals["MNAT3"] - Vitals["MDIS3"] - Vitals["MDOG3"]*prM_dog - Vitals["MCAR3"]*prM_car
	#Leslie matrix
	Leslie <- matrix(0,ncol=4,nrow=4)
	Leslie[1,2] <- Vitals["F1"] * S1
	Leslie[1,3] <- Vitals["F2"] * S2
	Leslie[1,4] <- Vitals["F3"] * S3
	Leslie[2,1] <- S0
	Leslie[3,2] <- S1
	Leslie[4,3] <- S2
	Leslie[4,4] <- S3
	
	rnew <- log(Re(eigen(Leslie)$values[1]))
	return(rnew)
}
###########################
#Calculate koala numbers post investment
###########################

#Set up a matrix where each row represents r for each investment value
rMatrix <- matrix(nrow=length(State), ncol=3, dimnames=list(State, c("r1", "r2", "r3")))
councils <- list(1,2,3)
rMatrix[,1:3] <- sapply(councils, function(i){sapply(State, function(invest) ReturnonInvest(invest,i))})

#You then loop ThreePatchMod over each row of this matrix to get the benefit arrays
#ie
#Benefit <- matrix(ncol=3, nrow=length(State), dimnames=list(State, apps)) 2D
BrisBenefit = array(NA,dim =c(length(State),length(State),length(State)))
LogaBenefit = array(NA,dim =c(length(State),length(State),length(State)))
RedlBenefit = array(NA,dim =c(length(State),length(State),length(State)))
for (i in seq_along(State)){
	for (j in seq_along(State)){
		for (k in seq_along(State)){
			Pars <- as.list(c(r1=rMatrix[as.character(State[i]),1], r2=rMatrix[as.character(State[j]),2], r3=rMatrix[as.character(State[k]),3], K1, K2, K3, Z=Z))#update the Pars with the post-investment r values
			NewN <- ode(Yini,Times,ThreePatchModC,Pars)#recalculate the population size
#the benefit is the population at the end of the time period but we could change this
#			BrisBenefit[i,j,k] <- NewN[length(Times),"N1"]#the benefit is the population at the end of the time period but we could change this
#			LogaBenefit[i,j,k] <- NewN[length(Times),"N2"]
#			RedlBenefit[i,j,k] <- NewN[length(Times),"N3"]

#the benefit is the percent overall decline in population by the time period
			BrisBenefit[i,j,k] <- abs(N1 - NewN[length(Times),"N1"])/N1
			LogaBenefit[i,j,k] <- abs(N2 - NewN[length(Times),"N2"])/N2
			RedlBenefit[i,j,k] <- abs(N3 - NewN[length(Times),"N3"])/N3
			#If we don't need to keep the benefit arrays
			NashPoss <- {BrisBenefit<=TargetDecline} + {LogaBenefit<=TargetDecline} + {RedlBenefit<=TargetDecline}

}}}

write.table(BrisBenefit, "C:/Users/uqgiacon/Dropbox/UQ Comp/Research/Game Theory Koalas/data/BrisBenefit25")
write.table(LogaBenefit, "C:/Users/uqgiacon/Dropbox/UQ Comp/Research/Game Theory Koalas/data/LogaBenefit25")
write.table(RedlBenefit, "C:/Users/uqgiacon/Dropbox/UQ Comp/Research/Game Theory Koalas/data/RedlBenefit25")

### which choices will the players make if their goal is to attain a decline in koala population size that is no higher than some target
###### can I get to those choices non iteratively

TargetDecline = 0.50 # 12 sets of choices that will attain this goal
#TargetDecline = 0.49 # 8 sets of choices that will attain this goal
#TargetDecline = 0.48 # 8 sets of choices
#TargetDecline = 0.47 # not possible for all three agencies
# NashLocatorBris = array(0,dim =c(length(State),length(State),length(State)))
# NashLocatorBris[which(BrisBenefit <= TargetDecline)] <- 1

# NashLocatorLoga = array(0,dim =c(length(State),length(State),length(State)))
# NashLocatorLoga[which(LogaBenefit <= TargetDecline)] <- 1

# NashLocatorRedl = array(0,dim =c(length(State),length(State),length(State)))
# NashLocatorRedl[which(RedlBenefit <= TargetDecline)] <- 1

## value indicates the number of agencies for which the target decline is met given that choice set (I need to label the dimensions for choices)

#NashPoss <- NashLocatorBris + NashLocatorLoga + NashLocatorRedl
NashPoss <- {BrisBenefit<=TargetDecline} + {LogaBenefit<=TargetDecline} + {RedlBenefit<=TargetDecline}


#### this approach didn't work ###############
#### need to fix the following statements so that they are always true. How to get the minimum budget amount. 

#NashOps <- which(NashPoss == max(NashPoss),TRUE)
#head(which(NashPoss == max(NashPoss),TRUE))
#head(NashOps[order(NashOps[,1]),])[1:10,]
#head(NashOps[order(NashOps[,2]),])
#head(NashOps[order(NashOps[,3]),])

#NashOps[which(NashOps[,1]==17),]

### if the sorting was correctly identifying the first row as the one that held the best (lowest) budget amount for all three
#TargetChoice <- which(NashPoss == max(NashPoss),TRUE)[1,] # the first instance (lowest budget amount for all three agencies) where the target is met

#BrisBenefit[TargetChoice[1],TargetChoice[2],TargetChoice[3]] # brisbane rate of decline at equilibrium
#State[TargetChoice[1]]# how much brisbane has to spend to get that benefit

#LogaBenefit[TargetChoice[1],TargetChoice[2],TargetChoice[3]] # logan rate of decline at equilibrium
#State[TargetChoice[2]]# how much Logan has to spend to get that benefit

#RedlBenefit[TargetChoice[1],TargetChoice[2],TargetChoice[3]] # Redlands rate of decline at equilibrium
#State[TargetChoice[3]]# how much Redlands has to spend to get that benefit

######################### Try iterative search for equilibrium 

#dimnames(NashPoss)<-list(State,State,State)
#In which places (values of State) is each council meeting their targets
TargetMet <- which(NashPoss == max(NashPoss),TRUE)
print(max(NashPoss))
dimnames(TargetMet)<-list(NULL, c("bris", "loga", "redl"))

PrintVals <- matrix(NA,1,6)

#Now we need to iterate over this to find the Nash(es)
### trying to shuffle order of council might make the output uninterpretable???

Nash<-matrix(NA,nrow=nrow(TargetMet), byrow=TRUE)
for(i in 1:100){
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

# the end set of spend choices 
EquilibriumSpendOptions <- aaply(na.omit(unique(PrintVals[,4:6])),c(1,2),function(x) State[x])
# hack to get a frequency distribution of the different choices.
freq<-table(as.numeric(paste(na.omit(PrintVals[,4]),na.omit(PrintVals[,5]),na.omit(PrintVals[,6]), sep = "")))

EquilibriumSpendOptions <- cbind(EquilibriumSpendOptions, freq)
outPath <- paste0(outDir, "Nash_10000subsamples2.csv")
outPath <- paste0(outDir, "Nash_10000_nomatch.csv")
write.csv(EquilibriumSpendOptions, outPath)

# have to look at the choice index values to match it back up to the table of equilibrium spend options
#unique(PrintVals[,4:6])


###############################
Not working yet
##############################################################################################################################

## LGA Choices when there is the option to contribute to a "Koala Fund"
## modifies "State" so that it incorporates the budget choices
## calculate final budget for each LGA given different allocation choices (but they all start with the same amount of money)

StartBudget <- 45 # total budget that each LGA has available to spend 
State <- c(0, 0.5, 1)  #c((0:10)/10) #choice of proportion of budget to spend
prop <- c(1/3, 1/3, 1/3) # proportion allocated from koala fund to each LGA

BrisBudget = array(NA,dim =c(length(State),length(State),length(State)))
LogBudget = array(NA,dim =c(length(State),length(State),length(State)))
RedBudget = array(NA,dim =c(length(State),length(State),length(State)))

for (i in seq_along(State)){
	for (j in seq_along(State)){
		for (k in seq_along(State)){
BrisBudget[i,j,k] <- State[i]*StartBudget + prop[1]*((1-State[i])*StartBudget + (1-State[j])*StartBudget + (1-State[k])*StartBudget)
LogBudget[i,j,k] <- State[j]*StartBudget + prop[2]*((1-State[i])*StartBudget + (1-State[j])*StartBudget + (1-State[k])*StartBudget)
RedBudget[i,j,k] <- State[k]*StartBudget + prop[3]*((1-State[i])*StartBudget + (1-State[j])*StartBudget + (1-State[k])*StartBudget)
}}} 


## NOW NEED TO CALCULATE r (return on investment) given each budget

Bris_r = array(NA,dim =c(length(State),length(State),length(State)))
Log_r = array(NA,dim =c(length(State),length(State),length(State)))
Red_r = array(NA,dim =c(length(State),length(State),length(State)))

for (i in seq_along(State)){
	for (j in seq_along(State)){
		for (k in seq_along(State)){
			Bris_r[i,j,k] <- ReturnonInvest(BrisBudget[i,j,k],1)
			Log_r[i,j,k] <- ReturnonInvest(LogBudget[i,j,k],2)
			Red_r[i,j,k] <- ReturnonInvest(RedBudget[i,j,k],3)
}}}		







###########################
#Calculate koala numbers post investment with koala fund accounted for
###########################

#Now have to apply ThreePatchMod to each cell of the arrays to get the benefit arrays

BrisBenefit = array(NA,dim =c(length(State),length(State),length(State)))
LogaBenefit = array(NA,dim =c(length(State),length(State),length(State)))
RedlBenefit = array(NA,dim =c(length(State),length(State),length(State)))
for (i in seq_along(State)){
	for (j in seq_along(State)){
		for (k in seq_along(State)){
			Pars <- as.list(c(r1=Bris_r[i,j,k], r2=Log_r[i,j,k], r3=Red_r[i,j,k], K1, K2, K3, Z=Z))#update the Pars with the post-investment r values
			NewN <- ode(Yini,Times,ThreePatchModC,Pars)#recalculate the population size
#the benefit is the population at the end of the time period but we could change this
			BrisBenefit[i,j,k] <- NewN[length(Times),"N1"]#the benefit is the population at the end of the time period but we could change this
			LogaBenefit[i,j,k] <- NewN[length(Times),"N2"]
			RedlBenefit[i,j,k] <- NewN[length(Times),"N3"]

#the benefit is the percent overall decline in population by the time period
#			BrisBenefit[i,j,k] <- abs(N1 - NewN[length(Times),"N1"])/N1
#			LogaBenefit[i,j,k] <- abs(N2 - NewN[length(Times),"N2"])/N2
#			RedlBenefit[i,j,k] <- abs(N3 - NewN[length(Times),"N3"])/N3
}}}

## Iteratively identify Nash equilibrium 
## represents the set of budget choices for the three regions that maximizes the number of koalas while considering what the others are doing

### each region makes a choice of how much it wants to spend
## initial choices
BrisChoice = 1
LogChoice = 1
RedChoice = 1

test <-vector() # 

while(sum(test)<3){

		#each player takes a turn to find the best choice they can make given the other's choices (order matters, but how much?)
		
 # maximize over Benefit surface
		newBrisChoice <- which.max(BrisBenefit[,LogChoice,RedChoice])		#find the best that Bris could do given the other choices
		newLogChoice <- which.max(LogaBenefit[newBrisChoice,,RedChoice])		#find the best that Logan could do given the other choices
		newRedChoice <- which.max(RedlBenefit[newBrisChoice,newLogChoice,])	#find the best that Logan could do given the other choices
 
	test <- c(BrisChoice==newBrisChoice, LogChoice==newLogChoice,RedChoice==newRedChoice) # stop trying when you don't move choices any more

		# if the new choices are better than the last, update the choice
		# if the new choices are better than the last, update the choice
		if (
			BrisBenefit[newBrisChoice,newLogChoice,newRedChoice] > BrisBenefit[BrisChoice,LogChoice,RedChoice]
			){BrisChoice <- newBrisChoice}

		if (
			LogaBenefit[newBrisChoice,newLogChoice,newRedChoice] > LogaBenefit[BrisChoice,LogChoice,RedChoice]
		){LogChoice <- newLogChoice}

		if (
			RedlBenefit[newBrisChoice,newLogChoice,newRedChoice] > RedlBenefit[BrisChoice,LogChoice,RedChoice]
		){RedChoice <- newRedChoice}	
		c(BrisChoice,LogChoice,RedChoice)
}		

NashChoices <- c(BrisChoice,LogChoice,RedChoice)

NashChoices

c(BrisBenefit[NashChoices[1],NashChoices[2],NashChoices[3]],LogaBenefit[NashChoices[1],NashChoices[2],NashChoices[3]],RedlBenefit[NashChoices[1],NashChoices[2],NashChoices[3]])





################################## old code ###############################################################
## Iteratively identify Nash equilibrium 
## represents the set of budget choices for the three regions that maximizes the number of koalas while considering what the others are doing
## seems like the benefit should really be koalas per spend?

### each region makes a choice of how much it wants to spend
## initial choices
## have to start with minimum values?

BrisChoice = 1
LogChoice = 1
RedChoice = 1

TargetDecline = 0.48

test <-vector() # have to zero it out or it stops before it does anything

while(sum(test)<3){

		#each player takes a turn to find the best choice they can make given the other's choices (order matters, but how much?)
		
 # maximize over Benefit surface
		#newBrisChoice <- which.max(BrisBenefit[,LogChoice,RedChoice])		#find the best that Bris could do given the other choices
		#newLogChoice <- which.max(LogaBenefit[newBrisChoice,,RedChoice])		#find the best that Logan could do given the other choices
		#newRedChoice <- which.max(RedlBenefit[newBrisChoice,newLogChoice,])	#find the best that Logan could do given the other choices

 # minimize cost within constraint

		if ( length(which(BrisBenefit[,LogChoice,RedChoice]<= TargetDecline))> 0
			){			
			newBrisChoice <- min(which(BrisBenefit[,LogChoice,RedChoice]<= TargetDecline))
		}else{
			newBrisChoice <- BrisChoice
		}
		if ( length(which(LogaBenefit[newBrisChoice,,RedChoice]<= TargetDecline)) >0
			){
			newLogChoice <- min(which(LogaBenefit[newBrisChoice,,RedChoice]<= TargetDecline))
		}else{
			newLogChoice <- LogChoice
		}
		if ( length(which(RedlBenefit[newBrisChoice,newLogChoice,]<= TargetDecline)) >0
			){
			newRedChoice <- min(which(RedlBenefit[newBrisChoice,newLogChoice,]<= TargetDecline))
		}else{
			newRedChoice <- RedChoice
		}

# this test might not work now because it would stop it without allowing anyone to move if they had no options... maybe that is ok.
	test <- c(BrisChoice==newBrisChoice, LogChoice==newLogChoice,RedChoice==newRedChoice) # stop trying when you don't move choices any more

# switched the direction of the choice decision  because now i am minimizing
# is this step actually selecting better options?

		# if the new choices are better than the last, update the choice
		if (
			BrisBenefit[newBrisChoice,newLogChoice,newRedChoice] < BrisBenefit[BrisChoice,LogChoice,RedChoice]
			){BrisChoice <- newBrisChoice}

		if (
			LogaBenefit[newBrisChoice,newLogChoice,newRedChoice] < LogaBenefit[BrisChoice,LogChoice,RedChoice]
		){LogChoice <- newLogChoice}

		if (
			RedlBenefit[newBrisChoice,newLogChoice,newRedChoice] < RedlBenefit[BrisChoice,LogChoice,RedChoice]
		){RedChoice <- newRedChoice}	
		c(BrisChoice,LogChoice,RedChoice)
}		

NashChoices <- c(BrisChoice,LogChoice,RedChoice)
NashChoices # first of the non-iteratively identified choices

# these are the choices that each player will take given what the other players are doing to attempt to each obtain their objective
c(BrisBenefit[NashChoices[1],NashChoices[2],NashChoices[3]],LogaBenefit[NashChoices[1],NashChoices[2],NashChoices[3]],RedlBenefit[NashChoices[1],NashChoices[2],NashChoices[3]])



