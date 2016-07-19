#R SCRIPT TO GET PARAMETER ESTIMATES

inpDirRasters <- "E:/Projects/seq_genetics/analysis/collaboration_paper/2_population models/Growth_Rates/input_rasters/"
inpDirModels <- "E:/Projects/seq_genetics/analysis/collaboration_paper/2_population models/Growth_Rates/input_coda/model_rt_abu_MOD3/output/"
rastApp <- c("bris", "redl", "loga")
outDir <- "E:/Projects/seq_genetics/analysis/collaboration_paper/2_population models/Growth_rates/"

library(maptools)
library(coda)
library(plyr)

########################
#MAIN PROCESSING #load in functions first (they are below)
########################
#LOOP THROUGH LGAs
for (i in 1:3){

	#GET SPATIAL DATA
	#Abund08 <- as.matrix(readAsciiGrid(paste0(inpDirRasters, "abund08", rastApp[i], ".txt")))
	Abund08 <- as.matrix(readAsciiGrid(paste0(inpDirRasters, "density08", rastApp[i], ".txt")))
	F08100 <- as.matrix(readAsciiGrid(paste0(inpDirRasters, "f08100", rastApp[i], ".txt")))
	F08500 <- as.matrix(readAsciiGrid(paste0(inpDirRasters, "f08500", rastApp[i], ".txt")))
	F081000 <- as.matrix(readAsciiGrid(paste0(inpDirRasters, "f081000", rastApp[i], ".txt")))
	F082000 <- as.matrix(readAsciiGrid(paste0(inpDirRasters, "f082000", rastApp[i], ".txt")))
	Strata08 <- as.matrix(readAsciiGrid(paste0(inpDirRasters, "strata08", rastApp[i], ".txt")))

	#GET CODA OUTPUT
	Mod1_Output <- get.coda.list(paste0(inpDirModels, "mod1_run"))
	Mod2_Output <- get.coda.list(paste0(inpDirModels, "mod2_run"))
	Mod3_Output <- get.coda.list(paste0(inpDirModels, "mod3_run"))
	Mod4_Output <- get.coda.list(paste0(inpDirModels, "mod4_run"))

	#GET SUMMARY DATA
	Summary_08100<-get.summary.data(Strata08,Abund08,F08100)
	gc()
	Summary_08500<-get.summary.data(Strata08,Abund08,F08500)
	gc()
	Summary_081000<-get.summary.data(Strata08,Abund08,F081000)
	gc()
	Summary_082000<-get.summary.data(Strata08,Abund08,F082000)
	gc()

	#GET THE VITAL RATES
	Rates_Mod1 <- get.vitals(Mod1_Output,Summary_08100[,1],Summary_08100[,2])
	Rates_Mod2 <- get.vitals(Mod2_Output,Summary_08500[,1],Summary_08500[,2])
	Rates_Mod3 <- get.vitals(Mod3_Output,Summary_081000[,1],Summary_081000[,2])
	Rates_Mod4 <- get.vitals(Mod4_Output,Summary_082000[,1],Summary_082000[,2])
	Rates_ModAv <- list(Rates_Mod1[[1]] * 0.159827 + Rates_Mod2[[1]] * 0.066659 + Rates_Mod3[[1]] * 0.306159 + Rates_Mod4[[1]] * 0.467357,Rates_Mod1[[2]] * 0.159827 + Rates_Mod2[[2]] * 0.066659 + Rates_Mod3[[2]] * 0.306159 + Rates_Mod4[[2]] * 0.467357)

	#sf parameters
	Sf<-rbind(Mod1_Output[[1]],Mod1_Output[[2]],Mod1_Output[[3]])[,"sf"]*0.159827 + rbind(Mod2_Output[[1]],Mod2_Output[[2]],Mod2_Output[[3]])[,"sf"]*0.066659 + rbind(Mod3_Output[[1]],Mod3_Output[[2]],Mod3_Output[[3]])[,"sf"]* 0.306159+rbind(Mod4_Output[[1]],Mod4_Output[[2]],Mod4_Output[[3]])[,"sf"]*0.467357

	#cf parameters
	Cf2<-rbind(Mod1_Output[[1]],Mod1_Output[[2]],Mod1_Output[[3]])[,"cf[2]"] * 0.159827 + rbind(Mod2_Output[[1]],Mod2_Output[[2]],Mod2_Output[[3]])[,"cf[2]"]*0.066659 + rbind(Mod3_Output[[1]],Mod3_Output[[2]],Mod3_Output[[3]])[,"cf[2]"]* 0.306159+rbind(Mod4_Output[[1]],Mod4_Output[[2]],Mod4_Output[[3]])[,"cf[2]"]*0.467357

	#output vital rates
	write.table(Rates_Mod1[[1]],paste0(outDir, "output_vitalrates/", "mod1_fem", rastApp[i], ".txt"),sep=",",row.names=F)
	write.table(Rates_Mod1[[2]],paste0(outDir,"output_vitalrates/", "mod1_mal", rastApp[i], ".txt"),sep=",",row.names=F)
	write.table(Rates_Mod2[[1]],paste0(outDir,"output_vitalrates/", "mod2_fem", rastApp[i], ".txt"),sep=",",row.names=F)
	write.table(Rates_Mod2[[2]],paste0(outDir,"output_vitalrates/", "mod2_mal", rastApp[i], ".txt"),sep=",",row.names=F)
	write.table(Rates_Mod3[[1]],paste0(outDir,"output_vitalrates/", "mod3_fem", rastApp[i], ".txt"),sep=",",row.names=F)
	write.table(Rates_Mod3[[2]],paste0(outDir,"output_vitalrates/", "mod3_mal", rastApp[i], ".txt"),sep=",",row.names=F)
	write.table(Rates_Mod4[[1]],paste0(outDir,"output_vitalrates/", "mod4_fem", rastApp[i], ".txt"),sep=",",row.names=F)
	write.table(Rates_Mod4[[2]],paste0(outDir,"output_vitalrates/", "mod4_mal", rastApp[i], ".txt"),sep=",",row.names=F)
	write.table(Rates_ModAv[[1]],paste0(outDir,"output_vitalrates/", "modav_fem", rastApp[i], ".txt"),sep=",",row.names=F)
	write.table(Rates_ModAv[[2]],paste0(outDir,"output_vitalrates/", "modav_mal", rastApp[i], ".txt"),sep=",",row.names=F)

	#GET THE DENSITY PREDICTIONS
	Dens_Mod1 <- get.pred.dens(Mod1_Output)
	Dens_Mod2 <- get.pred.dens(Mod2_Output)
	Dens_Mod3 <- get.pred.dens(Mod3_Output)
	Dens_Mod4 <- get.pred.dens(Mod4_Output)
	Dens_ModAv <- Dens_Mod1 * 0.159827 + Dens_Mod2 * 0.066659 + Dens_Mod3 * 0.306159 + Dens_Mod4 * 0.467357

	#output density predictions
	write.table(Dens_Mod1,paste0(outDir,"output_density/", "mod1_dens", rastApp[i], ".txt"),sep=",",row.names=F)
	write.table(Dens_Mod2,paste0(outDir,"output_density/","mod2_dens", rastApp[i], ".txt"),sep=",",row.names=F)
	write.table(Dens_Mod3,paste0(outDir,"output_density/","mod3_dens", rastApp[i], ".txt"),sep=",",row.names=F)
	write.table(Dens_Mod4,paste0(outDir,"output_density/","mod4_dens", rastApp[i], ".txt"),sep=",",row.names=F)
	write.table(Dens_ModAv,paste0(outDir,"output_density/","modav_dens", rastApp[i], ".txt"),sep=",",row.names=F)
	
	modelList <- list(Rates_Mod1, Rates_Mod2, Rates_Mod3, Rates_Mod4, Rates_ModAv)
	
	#GET THE LESLIE MATRIX #for females only
	meanModAvRates <- apply(Rates_ModAv[[1]], 2, mean)
	LeslieMA <- get.leslie(meanModAvRates)
	write.table(LeslieMA, paste0(outDir,"output_lambdas/","LeslieMatrix_", rastApp[i], ".txt"),sep=",",row.names=F)
	
	#GET THE LAMBDA VALUES #for females only
	Lambda_All <- laply(modelList, get.lambda)
	rownames(Lambda_All) <- c("Rates_Mod1", "Rates_Mod2", "Rates_Mod3", "Rates_Mod4", "Rates_ModAvFemales")
	write.table(Lambda_All, paste0(outDir,"output_lambdas/","LambdaF_", rastApp[i], ".txt"),sep=",",row.names=T)
	
}


#######################
#FUNCTIONS
#######################
#gets the sumnmary (aggreagted) of the abundance and habitat data 
get.summary.data<-function(Strata,Abund,Hab){
	Output <- matrix(NA,nrow=(max(Strata,na.rm=T) - min(Strata,na.rm=T) + 1),ncol=2)
	dimnames(Output) <- list(NULL,c("ABUND","HAB"))
	for (i in min(Strata,na.rm=T):max(Strata,na.rm=T)){
		Selection <- ifelse(Strata == i,1,NA)
		Output[i,1] <- sum(Selection * Abund * 0.0625 ,na.rm=T) 
		Output[i,2] <- mean(Selection * Hab,na.rm=T)
	}
	return(Output)
}

get.coda.list<-function(Model.Path){
	Ch1<-read.coda(paste(Model.Path,"1.txt",sep=""),paste(Model.Path,"Index.txt",sep=""))
	Ch2<-read.coda(paste(Model.Path,"2.txt",sep=""),paste(Model.Path,"Index.txt",sep=""))
	Ch3<-read.coda(paste(Model.Path,"3.txt",sep=""),paste(Model.Path,"Index.txt",sep=""))
	All<-mcmc.list(Ch1,Ch2,Ch3)
	return(All)
}

#calculates average juvenile mortality for sex l
#jslp.l. = a vector of the posterior parameter values
fJM <- function(jslp.l.){
	return(apply(as.matrix(jslp.l.),MARGIN=1,function(X){return(1 - (exp(X) / (1 + exp(X))))}))
}

#calculates the average mortality rate for cause m, sex l, and age class k  
#m = index for cause, slp.l.k. = vector of posterior parameter values
#clp.l.m. = matrix of posterior parameter values for all four causes
#cf.m. = matrix of parameter values for all four causes
#A = abundance raster layer (matrix)
#H = proportion forest raster layer (matrix)
fM <- function(m,slp.l.k.,sf,clp.l.m.,cf.m.,A,H){
		mortalities <- function(X,Abund,Hab){ #set up the mortalities function
		Mort <- 1 - ((exp(X[1] + X[2] * Hab) / (1 + exp(X[1] + X[2] * Hab)))^365)
		if (m == 1)	{
			Cause <- 1 / (1 + exp(X[2 + 1] + X[5 + 1] * Hab) + exp(X[2 + 2] + X[5 + 2] * Hab) + exp(X[2 + 3] + X[5 + 3] * Hab))
		} else {
			Cause <- exp(X[2 + m - 1] + X[5 + m - 1] * Hab) / (1 + exp(X[2 + 1] + X[5 + 1] * Hab) + exp(X[2 + 2] + X[5 + 2] * Hab) + exp(X[2 + 3] + X[5 + 3] * Hab))
		}
		Total <- Abund * Mort * Cause
		Mean <- sum(Total,na.rm=T) / sum(Abund,na.rm=T)
		return(Mean)
	}
	Params <- as.matrix(cbind(slp.l.k.,sf,clp.l.m.,cf.m.)) #create a matrix of the parameter values
	Mortalities <- apply(Params,MARGIN=1,FUN=mortalities,Abund=A,Hab=H) #get the annual mortality rates for each parameter combination	
	return(Mortalities)#return the mortality rates
}

#calculates the average fecundity rates in terms of production of sex l  
#l = index of sex (1 = female, 2 = male), blp.k. = vector of posterior parameter values
#bslp = vector of posterior parameter values
#A = abundance raster layer (matrix)
#H = proportion forest raster layer (matrix)
fF <- function(blp.k.,bslp,l,A,H){
	#set up the fecundity function
	fecundities <- function(X,Sex,Abund,Hab){
		if (l == 1) {#females
		Fec <- (exp(X[1]) / (1 + exp(X[1]))) * (exp(X[2]) / (1 + exp(X[2])))
		} else if (l == 2) { #males
		Fec <- (exp(X[1]) / (1 + exp(X[1]))) * (1 - (exp(X[2]) / (1 + exp(X[2]))))
		}
		Total <- Abund * Fec
		Mean <- sum(Total,na.rm=T) / sum(Abund,na.rm=T)
		return(Mean)
	}
	Params <- as.matrix(cbind(blp.k.,bslp)) #create a matrix of the parameter values
	Fecundities <- apply(Params,MARGIN=1,FUN=fecundities,Sex=l,Abund=A,Hab=H) #get the annual fecundity rates for each parameter combination
	return(Fecundities)#return the fecundity rates
}

#code for calculating the vital rates
get.vitals <- function(Coda,A,H){
	Params <- Coda[[1]]
	if(length(Params) > 1){
		for (i in 2:length(Coda)){
			Params <- rbind(Params,Coda[[i]])
		}
	}
	#set up output matrices
	Output_Females <- matrix(NA,nrow=nrow(Params),ncol=16)
	dimnames(Output_Females) <- list(NULL,c("F1","F2","F3","M0","MNAT1","MCAR1","MDOG1","MDIS1","MNAT2","MCAR2","MDOG2","MDIS2","MNAT3","MCAR3","MDOG3","MDIS3"))
	Output_Males <- matrix(NA,nrow=nrow(Params),ncol=16)
	dimnames(Output_Males) <- list(NULL,c("F1","F2","F3","M0","MNAT1","MCAR1","MDOG1","MDIS1","MNAT2","MCAR2","MDOG2","MDIS2","MNAT3","MCAR3","MDOG3","MDIS3"))
		
	#get values - females
	Output_Females[,"F1"] <- fF(Params[,"blp[1]"],Params[,"bslp"],1,A,H)
	Output_Females[,"F2"] <- fF(Params[,"blp[2]"],Params[,"bslp"],1,A,H)
	Output_Females[,"F3"] <- fF(Params[,"blp[3]"],Params[,"bslp"],1,A,H)
	Output_Females[,"M0"] <-fJM(Params[,"jslp[1]"])
	Output_Females[,"MNAT1"] <- fM(1,Params[,"slp[1,1]"],Params[,"sf"],Params[,c("clp[1,2]","clp[1,3]","clp[1,4]")],cbind(Params[,c("cf[2]")],matrix(0,nrow=nrow(Params),ncol=1),matrix(0,nrow=nrow(Params),ncol=1)),A,H)
	Output_Females[,"MCAR1"] <- fM(2,Params[,"slp[1,1]"],Params[,"sf"],Params[,c("clp[1,2]","clp[1,3]","clp[1,4]")],cbind(Params[,c("cf[2]")],matrix(0,nrow=nrow(Params),ncol=1),matrix(0,nrow=nrow(Params),ncol=1)),A,H)
	Output_Females[,"MDOG1"] <- fM(3,Params[,"slp[1,1]"],Params[,"sf"],Params[,c("clp[1,2]","clp[1,3]","clp[1,4]")],cbind(Params[,c("cf[2]")],matrix(0,nrow=nrow(Params),ncol=1),matrix(0,nrow=nrow(Params),ncol=1)),A,H)
	Output_Females[,"MDIS1"] <- fM(4,Params[,"slp[1,1]"],Params[,"sf"],Params[,c("clp[1,2]","clp[1,3]","clp[1,4]")],cbind(Params[,c("cf[2]")],matrix(0,nrow=nrow(Params),ncol=1),matrix(0,nrow=nrow(Params),ncol=1)),A,H)
	Output_Females[,"MNAT2"] <- fM(1,Params[,"slp[1,2]"],Params[,"sf"],Params[,c("clp[1,2]","clp[1,3]","clp[1,4]")],cbind(Params[,c("cf[2]")],matrix(0,nrow=nrow(Params),ncol=1),matrix(0,nrow=nrow(Params),ncol=1)),A,H)
	Output_Females[,"MCAR2"] <- fM(2,Params[,"slp[1,2]"],Params[,"sf"],Params[,c("clp[1,2]","clp[1,3]","clp[1,4]")],cbind(Params[,c("cf[2]")],matrix(0,nrow=nrow(Params),ncol=1),matrix(0,nrow=nrow(Params),ncol=1)),A,H)
	Output_Females[,"MDOG2"] <- fM(3,Params[,"slp[1,2]"],Params[,"sf"],Params[,c("clp[1,2]","clp[1,3]","clp[1,4]")],cbind(Params[,c("cf[2]")],matrix(0,nrow=nrow(Params),ncol=1),matrix(0,nrow=nrow(Params),ncol=1)),A,H)
	Output_Females[,"MDIS2"] <- fM(4,Params[,"slp[1,2]"],Params[,"sf"],Params[,c("clp[1,2]","clp[1,3]","clp[1,4]")],cbind(Params[,c("cf[2]")],matrix(0,nrow=nrow(Params),ncol=1),matrix(0,nrow=nrow(Params),ncol=1)),A,H)
	Output_Females[,"MNAT3"] <- fM(1,Params[,"slp[1,3]"],Params[,"sf"],Params[,c("clp[1,2]","clp[1,3]","clp[1,4]")],cbind(Params[,c("cf[2]")],matrix(0,nrow=nrow(Params),ncol=1),matrix(0,nrow=nrow(Params),ncol=1)),A,H)
	Output_Females[,"MCAR3"] <- fM(2,Params[,"slp[1,3]"],Params[,"sf"],Params[,c("clp[1,2]","clp[1,3]","clp[1,4]")],cbind(Params[,c("cf[2]")],matrix(0,nrow=nrow(Params),ncol=1),matrix(0,nrow=nrow(Params),ncol=1)),A,H)
	Output_Females[,"MDOG3"] <- fM(3,Params[,"slp[1,3]"],Params[,"sf"],Params[,c("clp[1,2]","clp[1,3]","clp[1,4]")],cbind(Params[,c("cf[2]")],matrix(0,nrow=nrow(Params),ncol=1),matrix(0,nrow=nrow(Params),ncol=1)),A,H)
	Output_Females[,"MDIS3"] <- fM(4,Params[,"slp[1,3]"],Params[,"sf"],Params[,c("clp[1,2]","clp[1,3]","clp[1,4]")],cbind(Params[,c("cf[2]")],matrix(0,nrow=nrow(Params),ncol=1),matrix(0,nrow=nrow(Params),ncol=1)),A,H)
		
	#get values - males
	Output_Males[,"F1"] <- fF(Params[,"blp[1]"],Params[,"bslp"],2,A,H)
	Output_Males[,"F2"] <- fF(Params[,"blp[2]"],Params[,"bslp"],2,A,H)
	Output_Males[,"F3"] <- fF(Params[,"blp[3]"],Params[,"bslp"],2,A,H)
	Output_Males[,"M0"] <- fJM(Params[,"jslp[2]"])
	Output_Males[,"MNAT1"] <- fM(1,Params[,"slp[2,1]"],Params[,"sf"],Params[,c("clp[2,2]","clp[2,3]","clp[2,4]")],cbind(Params[,c("cf[2]")],matrix(0,nrow=nrow(Params),ncol=1),matrix(0,nrow=nrow(Params),ncol=1)),A,H)
	Output_Males[,"MCAR1"] <- fM(2,Params[,"slp[2,1]"],Params[,"sf"],Params[,c("clp[2,2]","clp[2,3]","clp[2,4]")],cbind(Params[,c("cf[2]")],matrix(0,nrow=nrow(Params),ncol=1),matrix(0,nrow=nrow(Params),ncol=1)),A,H)
	Output_Males[,"MDOG1"] <- fM(3,Params[,"slp[2,1]"],Params[,"sf"],Params[,c("clp[2,2]","clp[2,3]","clp[2,4]")],cbind(Params[,c("cf[2]")],matrix(0,nrow=nrow(Params),ncol=1),matrix(0,nrow=nrow(Params),ncol=1)),A,H)
	Output_Males[,"MDIS1"] <- fM(4,Params[,"slp[2,1]"],Params[,"sf"],Params[,c("clp[2,2]","clp[2,3]","clp[2,4]")],cbind(Params[,c("cf[2]")],matrix(0,nrow=nrow(Params),ncol=1),matrix(0,nrow=nrow(Params),ncol=1)),A,H)
	Output_Males[,"MNAT2"] <- fM(1,Params[,"slp[2,2]"],Params[,"sf"],Params[,c("clp[2,2]","clp[2,3]","clp[2,4]")],cbind(Params[,c("cf[2]")],matrix(0,nrow=nrow(Params),ncol=1),matrix(0,nrow=nrow(Params),ncol=1)),A,H)
	Output_Males[,"MCAR2"] <- fM(2,Params[,"slp[2,2]"],Params[,"sf"],Params[,c("clp[2,2]","clp[2,3]","clp[2,4]")],cbind(Params[,c("cf[2]")],matrix(0,nrow=nrow(Params),ncol=1),matrix(0,nrow=nrow(Params),ncol=1)),A,H)
	Output_Males[,"MDOG2"] <- fM(3,Params[,"slp[2,2]"],Params[,"sf"],Params[,c("clp[2,2]","clp[2,3]","clp[2,4]")],cbind(Params[,c("cf[2]")],matrix(0,nrow=nrow(Params),ncol=1),matrix(0,nrow=nrow(Params),ncol=1)),A,H)
	Output_Males[,"MDIS2"] <- fM(4,Params[,"slp[2,2]"],Params[,"sf"],Params[,c("clp[2,2]","clp[2,3]","clp[2,4]")],cbind(Params[,c("cf[2]")],matrix(0,nrow=nrow(Params),ncol=1),matrix(0,nrow=nrow(Params),ncol=1)),A,H)
	Output_Males[,"MNAT3"] <- fM(1,Params[,"slp[2,3]"],Params[,"sf"],Params[,c("clp[2,2]","clp[2,3]","clp[2,4]")],cbind(Params[,c("cf[2]")],matrix(0,nrow=nrow(Params),ncol=1),matrix(0,nrow=nrow(Params),ncol=1)),A,H)
	Output_Males[,"MCAR3"] <- fM(2,Params[,"slp[2,3]"],Params[,"sf"],Params[,c("clp[2,2]","clp[2,3]","clp[2,4]")],cbind(Params[,c("cf[2]")],matrix(0,nrow=nrow(Params),ncol=1),matrix(0,nrow=nrow(Params),ncol=1)),A,H)
	Output_Males[,"MDOG3"] <- fM(3,Params[,"slp[2,3]"],Params[,"sf"],Params[,c("clp[2,2]","clp[2,3]","clp[2,4]")],cbind(Params[,c("cf[2]")],matrix(0,nrow=nrow(Params),ncol=1),matrix(0,nrow=nrow(Params),ncol=1)),A,H)
	Output_Males[,"MDIS3"] <- fM(4,Params[,"slp[2,3]"],Params[,"sf"],Params[,c("clp[2,2]","clp[2,3]","clp[2,4]")],cbind(Params[,c("cf[2]")],matrix(0,nrow=nrow(Params),ncol=1),matrix(0,nrow=nrow(Params),ncol=1)),A,H)

	#return output
	return(list(Output_Females,Output_Males))
}

get.pred.dens <- function(Coda){
	Pred <- Coda[[1]]
	if(length(Pred) > 1){
		for (i in 2:length(Coda)){
			Pred <- rbind(Pred,Coda[[i]])
		}
	}
	#get output
	Output <- Pred[,c("meanD97[1]","meanD97[2]","meanD97[3]","meanD97[4]","meanD97[5]","meanD05[1]","meanD05[2]","meanD05[3]","meanD05[4]","meanD05[5]","meanD08[1]","meanD08[2]","meanD08[3]","meanD08[4]","meanD08[5]")]
	return(Output)
}
#########################
#function for setting up Leslie matrix
#input is a row from the output of get.vitals
#the female or male part
get.leslie <- function(Params){
	Leslie <- matrix(0,ncol=4,nrow=4)
	#get the aggregate survival rates
	S1 <- 1 - Params[5] - Params[6] - Params[7] - Params[8]
	S2 <- 1 - Params[9] - Params[10] - Params[11] - Params[12]
	S3 <- 1 - Params[13] - Params[14] - Params[15] - Params[16]
	Leslie[1,2] <- Params[1] * S1
	Leslie[1,3] <- Params[2] * S2
	Leslie[1,4] <- Params[3] * S3
	Leslie[2,1] <- 1 - Params[4]
	Leslie[3,2] <- S1
	Leslie[4,3] <- S2
	Leslie[4,4] <- S3
	
	return(Leslie)
}

#code for calculating expected lambda list
#argument is the output from get.vitals
#and is for females only
get.lambda.list <- function(Vitals){
	#get the aggregate mortality rates
	Lambdas <- apply(Vitals[[1]],MARGIN=1,FUN=function(X){eigen(get.leslie(X))$values[1]})
	return(Lambdas)#return the expected lambda
}

#code for calculating expected lambda
#argument is the output from get.vitals
#and is for females only
get.lambda <- function(Vitals){
	#get the aggregate mortality rates
	Lambdas <- apply(Vitals[[1]],MARGIN=1,FUN=function(X){eigen(get.leslie(X))$values[1]})
	LambdaStats <- data.frame(Re(mean(Lambdas)), matrix(Re(quantile(Lambdas,c(0.025,0.975))), nrow=1))
	names(LambdaStats)<- c("mean", names(quantile(Lambdas,c(0.025,0.975))))
	return(LambdaStats) #return the expected lambda
	}

#GETS MEAN MODEL AVERAGE LAMBDA #across males and females
get.mod.av.lambda <- function(Rates1,Rates2,Rates3,Rates4,W1,W2,W3,W4){
	MARates <- Rates1
	MARates[[1]] <- W1 * Rates1[[1]] + W2 * Rates2[[1]] + W3 * Rates3[[1]] + W4 * Rates4[[1]]
	MARates[[2]] <- W1 * Rates1[[2]] + W2 * Rates2[[2]] + W3 * Rates3[[2]] + W4 * Rates4[[2]]
	MALambda <- get.lambda(MARates)
	return(as.numeric(MALambda))
}

#GETS A LIST OF THE MEAN MODEL AVERAGE LAMBDA
Rates_Mod1[[1]] * 0.159827 + Rates_Mod2[[1]] * 0.066659 + Rates_Mod3[[1]] * 0.306159 + Rates_Mod4[[1]] * 0.467357
get.mod.av.lambda.list <- function(Rates1,Rates2,Rates3,Rates4,W1,W2,W3,W4){
	MARates <- Rates1
	MARates[[1]] <- W1 * Rates1[[1]] + W2 * Rates2[[1]] + W3 * Rates3[[1]] + W4 * Rates4[[1]]
	MARates[[2]] <- W1 * Rates1[[2]] + W2 * Rates2[[2]] + W3 * Rates3[[2]] + W4 * Rates4[[2]]
	MALambda <- get.lambda.list(MARates)
	return(as.numeric(MALambda))
}

