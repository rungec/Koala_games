###TEST script for Game_Investment_benefit_eval_dispersal_assymetry_150907.r
#Run only rows 1:5


clust <- makeCluster(ncore, type="SOCK")
#export variables/functions
clusterEvalQ(clust, {library(plyr);library(dplyr);library(deSolve);library(nloptr);library(compiler)})

clusterExport(clust, c("area", "Popn2008", "N1", "N2", "N3", "r1", "r2", "r3", "K1", "K2", "K3", "Vitals", "Times", "StateIncr", "State", "Yini", "ThreePatchModC", "ReturnonInvest", "bindata", "NashPossFun","NashFun", "rMatrix"))
registerDoSNOW(clust)

#seq_len(nrow(bindata[1:5,]) for first 5 rows for testing
assymROI <- llply(seq_len(nrow(bindata[1:5,])), .parallel=TRUE, function(currRow) {
	#new Z (assign to global so all functions can see it)
	Z<<-matrix(c(NA, bindata[currRow,1], bindata[currRow,2], bindata[currRow,4], NA, bindata[currRow,3], bindata[currRow,5], bindata[currRow,6], NA), ncol=3, dimnames=NULL)
	
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


#Save data
OverallInvStrat<- laply(seq_len(nrow(bindata[1:5,])), function(x) assymROI[[x]]$OverallOptim)
OverallInvStrat <- cbind(OverallInvStrat, bindata[1:5,"symmcoef"], seq_len(nrow(bindata[1:5,])))
dimnames(OverallInvStrat) <- list(NULL, c("bris", "loga", "redl", "symmcoef", "run"))
outPath <- paste0(outDir, "InvestmentStrategies_Overall50decline_100randomZ.csv")
write.csv(OverallInvStrat, outPath)

#Save data
LocalInvStrat <- ldply(seq_len(nrow(bindata[1:5,])), function(x) {					
				currROI <-matrix(assymROI[[x]]$LocalOptim, ncol=4)
				currMaxN <- assymROI[[x]]$councilsmettargets
				return(cbind(currROI, rep(bindata[x,"symmcoef"], nrow(currROI)), rep(x, nrow(currROI)), rep(currMaxN, nrow(currROI))))
				})
names(LocalInvStrat) <- c("bris", "loga", "redl", "freq", "symmcoef", "run", "councilsmettargets")
outPath <- paste0(outDir, "InvestmentStrategies_Local50decline_100randomZ.csv")
write.csv(LocalInvStrat, outPath)


ROIs <- ldply(1:5, function(i) {
			#print(assymROI[[i]]$OverallOptim)
			newr <<- Map(ReturnonInvest, assymROI[[i]]$OverallOptim, 1:3)
			names(newr)<-NULL
			Z<<-assymROI[[i]]$currZ
			#print(Z)
			Pars <<- as.list(c(r1=newr[1], r2=newr[2], r3=newr[3], K1, K2, K3, Z=Z))
			#print(Pars)
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









ROIMAXs <- ldply(1:5, function(i) {
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
outPath <- paste0(outDir, "ROI_overall50decline_100randomZ_max.csv")
write.csv(ROIMAXs, outPath)


#Calculate how many koalas each budget choice results in (Local scenario)
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
outPath <- paste0(outDir, "ROI_local50decline_100randomZ_max.csv")
write.csv(ROILocal, outPath)


#Extra test scripts

# ################
# #Scenario 4: Benevolent social planner
# ################

# S4fun <- cmpfun(function(Popn2008, K1, K2, K3, Z, State, StateIncr){
	# idealpop <- sum(Popn2008)*0.5 #50% decline
	# investList <- sapply(c(1:10), function(i) list(runif(3,min(State),max(State)))) #a list of random start values for cobyla fn
	# # fn<- function(invest) sum(invest) #thing to minimise
	# #subject to inequality  constraint hin(x)>=0 
	# hin <-function(invest) {
			# investReturn <- Map(ReturnonInvest, invest, 1:3) # the new growth rate due to investment
			# Pars2 <- as.list(c(r1= investReturn[1],r2=investReturn[2], r3=investReturn[3], K1=K1,K2=K2, K3=K3, Z=Z))
			# NewN <- ode(Yini,Times,ThreePatchModC,Pars2)
			# return(sum(NewN[nrow(NewN),2:4])-idealpop) 
			# }
		# lower <-c(min(State),min(State),min(State)) 
		# upper <-c(max(State),max(State),max(State))

	# ###OPTIMISE
	# S4optimDF <- ldply(investList, function(i) {
					# S4optimCOBYLA <- cobyla(i, sum, lower=lower, upper=upper, hin = hin, control = list(xtol_rel=1e-6, maxeval=500))
					# #S4optimCOBYLA <- cobyla(i, fn, lower=lower, upper=upper, hin = hin, control = list(ftol_abs=rep(StateIncr, 3), maxeval=500)) #stops when increments by less than the bins in State
					# return(t(unlist(S4optimCOBYLA[1:4])))
					# }) #input list return df

	# names(S4optimDF)[1:4] <- c("bris_invest", "loga_invest", "redl_invest", "overall_invest")
	# return(S4optimDF)
	# })




# S4discretefun <- cmpfun(function(State, rMatrix, K1, K2, K3, Z, Popn2008){
	# TargetDecline <-0.5
	# idealpop <- sum(Popn2008)*TargetDecline
	# #Set up 3D array of ROI for each council (Benefit is % decline)
	# OverallBenefit = array(NA,dim =c(length(State),length(State),length(State)))
	# BrisBenefit = array(NA,dim =c(length(State),length(State),length(State)))
	# LogaBenefit = array(NA,dim =c(length(State),length(State),length(State)))
	# RedlBenefit = array(NA,dim =c(length(State),length(State),length(State)))
	# for (i in seq_along(State)){ #try Map instead then mdapply
		# for (j in seq_along(State)){
			# for (k in seq_along(State)){
				# Pars <- as.list(c(r1=rMatrix[as.character(State[i]),1], r2=rMatrix[as.character(State[j]),2], r3=rMatrix[as.character(State[k]),3], K1, K2, K3, Z=Z))#update the Pars with the post-investment r values
				# NewN <- ode(Yini,Times,ThreePatchModC,Pars)#recalculate the population size
				# #the benefit is the percent overall decline in population by the time period
			# OverallBenefit[i,j,k] <- sum(NewN[nrow(NewN),2:4])-idealpop
			# BrisBenefit[i,j,k] <- abs(N1 - NewN[length(Times),"N1"])/N1
			# LogaBenefit[i,j,k] <- abs(N2 - NewN[length(Times),"N2"])/N2
			# RedlBenefit[i,j,k] <- abs(N3 - NewN[length(Times),"N3"])/N3
				# }}}
# #Calculate overall best spend (Scenario 4 - benevolent planner)
# OptimLocator = array(0,dim =c(length(State),length(State),length(State)))
# OptimLocator[which(OverallBenefit >= 0)] <- 1
# dimnames(OptimLocator)<-list(State,State,State)
# CostMatr <- which(OptimLocator == 1,TRUE)
# optimalcost <- State[CostMatr[which.min(rowSums(CostMatr)),]]
# names(optimalcost)<-c("bris", "loga", "redl")

# #Calculate Nash (Scenario 3 - non-cooperative game)
# NashLocatorBris = array(0,dim =c(length(State),length(State),length(State)))
# NashLocatorBris[which(BrisBenefit <= TargetDecline)] <- 1
# NashLocatorLoga = array(0,dim =c(length(State),length(State),length(State)))
# NashLocatorLoga[which(LogaBenefit <= TargetDecline)] <- 1
# NashLocatorRedl = array(0,dim =c(length(State),length(State),length(State)))
# NashLocatorRedl[which(RedlBenefit <= TargetDecline)] <- 1

# NashPoss <- NashLocatorBris + NashLocatorLoga + NashLocatorRedl
# dimnames(NashPoss)<-list(State,State,State)
# #In which places (values of State) is each council meeting their targets
# TargetMet <- which(NashPoss == max(NashPoss),TRUE)
# dimnames(TargetMet)<-list(NULL, c("bris", "loga", "redl"))
# #Now we need to iterate over this to find the Nash(es)

# Nash<-matrix(NA,nrow=nrow(TargetMet), byrow=TRUE)
# for(i in 1:nrow(TargetMet)){
	# currSpend <- TargetMet[i,]
	
	# #find the best that Bris could do given the other choices
	# newBrisChoice <- min(TargetMet[which(TargetMet[,"loga"]==currSpend["loga"]&TargetMet[,"redl"]==currSpend["redl"]),"bris"])
	# #find the best that Logan could do given the other choices
	# newLogChoice <- min(TargetMet[which(TargetMet[,"bris"]==currSpend["bris"]&TargetMet[,"redl"]==currSpend["redl"]),"loga"])	
	# #find the best that Logan could do given the other choices
	# newRedChoice <- min(TargetMet[which(TargetMet[,"bris"]==currSpend["bris"]&TargetMet[,"loga"]==currSpend["loga"]),"redl"])
	
	# if (
	# sum(newBrisChoice==currSpend["loga"], newLogChoice==currSpend["loga"],newRedChoice==currSpend["loga"])==3 )
	# {
	# Nash[i]<-TRUE
	# } else {
	# Nash[i]<-FALSE
	# }
	# }
# NashSpend <-State[TargetMet[which(Nash==TRUE),]]
# return(list(Optimal=optimalcost,NashSpend=NashSpend))
# })


# ######################
# #Do in Parallel
# ######################
# ncore<-3

# clust <- makeCluster(ncore, type="SOCK")
# #export variables/functions
# clusterEvalQ(clust, {library(plyr);library(dplyr);library(deSolve);library(nloptr);library(compiler)})

# clusterExport(clust, c("area", "Popn2008", "N1", "N2", "N3", "r1", "r2", "r3", "K1", "K2", "K3", "Vitals", "Times", "StateIncr", "State", "Pars", "Yini", "ThreePatchModC", "ReturnonInvest", "bindata", "S4fun","S4discretefun"))
# registerDoSNOW(clust)

# assymROI <- llply(seq_len(nrow(bindata[1:5,])), .parallel=TRUE, function(currRow) {
	# #new Z (assign to global so all functions can see it)
	# Z<<-matrix(c(NA, bindata[currRow,1], bindata[currRow,2], bindata[currRow,4], NA, bindata[currRow,3], bindata[currRow,5], bindata[currRow,6], NA), ncol=3, dimnames=NULL)
	
	# ret=S4discretefun(State, rMatrix, K1, K2, K3, Z, Popn2008)
	# return(list(data=ret))
	# #ret=S4fun(Popn2008, K1, K2, K3, Z, State, StateIncr)
	# #return(list(data=ret, value=as.numeric(sapply(ret[which(ret$convergence==4),1:3], median))))
# })
 
# clust=stopCluster(clust) 
 # #taskkill -im Rscript.exe -f










