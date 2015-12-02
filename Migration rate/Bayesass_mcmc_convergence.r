#This script analyses mcmc convergence on Bayesass 3.0 outputs
#it assumes each chain (run) is stored in a separate folder listed in 'runs'

#install.packages('coda')
library(coda)
nruns <- 10 #the number of chains
runs <- paste0(rep("run0", nruns), 80:89)
runDir <- "W:/seq_genetics/analysis/collaboration_paper/Bayesass_ClaireR/output/"
CouncilIndex <- c("Brisbane", "Gold Coast", "Ipswich", "Logan", "Moreton Bay", 
"Redland", "Scenic Rim", "Somerset")

#sapply(runs, makeplot)

#makeplot("run005")


# Change this value to the actual burnin used for your MCMC run
burnin = 30000000
# Change this value to the actual sampling interval used for your MCMC run
sampling.interval = 10000
makeplot <- function(run){

	inpDir <- paste0(runDir, run, "/BA3trace.txt")
	x<-read.table(inpDir, sep="\t", header=T, stringsAsFactors=FALSE)
	#x1<-mcmc(x[,1:(ncol(x)-1)])
	x1<-mcmc(x$LogProb, thin=sampling.interval)

	pdf(paste0(runDir, run, "/traceplot.pdf"), width=10, height=5)
	# Plotting the likelihoods
	par(mfrow=c(1,2))
	plot(x$State,x$LogProb, xlab="State", ylab="LogProb", type="n", lwd=2, font.lab=2, main=run)
	lines(x$State[1:which(x$State==burnin)], x$LogProb[1:which(x$State==burnin)], col="grey70")
	lines(x$State[which(x$State>=burnin)], x$LogProb[which(x$State>=burnin)], col="black")
	densplot(x1)
	dev.off()
	pdf(paste0(runDir, run, "/autocorrplot.pdf"), width=5, height=5)
	autocorr.plot(x1)
	dev.off()
	#Calculate Bayesian deviance
	range = (x$State > burnin & x$State %% sampling.interval == 0)
	D = -2*mean(x$LogProb[range])
	return(D)
	}
	
BayesDev <- sapply(runs, makeplot)
write.csv(BayesDev, paste0(runDir, runs[1], "/BayesDeviance.csv"))
	
##############
#calculate convergence between runs
mcmclist <- function(run){
	inpDir <- paste0(runDir, run, "/BA3trace.txt")
	x<-read.table(inpDir, sep="\t", header=T, stringsAsFactors=FALSE)
	#x1<-mcmc(x[,1:(ncol(x)-1)])
	x1<-mcmc(x[c("State", "LogProb")], thin=sampling.interval, start=burnin) #thin=the sampling interval start=the burnin
	}
mh.list <- lapply(runs, mcmclist)
sink(paste0(runDir, runs[1], "/BayesStats2.csv"))
#Gelman & Rubin diagnostic
cat("Gelman & Rubin diagnostic", "\n")
gelman.diag(mh.list, autoburnin=FALSE)
cat("\n", "Summary", "\n")
lapply(mh.list, summary)
sink()

##############
#plot convergence

library(plyr)
library(dplyr)
library(stringr)

allmigRates <- c()
allinbRates <- c()
allalleFreq <- c()
nloci <- 6

for (i in seq_along(runs)){
run <- runs[i]
inpDir <- paste0(runDir, run, "/", run, ".txt")
x<-readLines(inpDir)
popns <- x[5]%>% strsplit(" ")%>%"[["(1)%>%"["(4)%>%as.numeric()
migIndex <- grep("Migration", x)
migLines <- x[(migIndex+2):(migIndex+1+popns)]
migRates <- migLines%>% str_trim()%>%strsplit(" ")%>%unlist(use.names=FALSE, recursive=FALSE)
migRates2 <- migRates[seq(2, length(migRates), 2)]%>%strsplit("\\(")%>%laply("[[", 1)%>%as.numeric() 
allmigRates <- cbind(allmigRates, migRates2)
migNames <- migRates[seq(1, length(migRates), 2)]%>%strsplit(":")%>%laply("[[", 1)

inbrIndex <- grep("Inbreeding", x)
inbrLines <- x[(inbrIndex+2):(inbrIndex+1+popns)]
inbrRates <- inbrLines%>%strsplit(" ")%>%unlist(use.names=FALSE, recursive=FALSE)
inbrRates <- inbrRates[seq(4, length(inbrRates), 4)]%>%strsplit("\\(")%>%laply("[[", 1)%>%as.numeric()
allinbRates <- cbind(allinbRates, inbrRates)

alleIndex <- grep("Allele Freq", x)
alleLines <- x[(alleIndex+2):length(x)]
alleLinesub <- alleLines[unlist(sapply(seq(1, by=14, length.out=popns), function(x) seq(x+2, by=2, length.out=nloci)))]
allefreqs <- alleLinesub%>%str_trim()%>%strsplit(" ")%>%unlist(use.names=FALSE, recursive=FALSE)%>%strsplit(":")%>%laply("[[", 2)%>%strsplit("\\(")%>%laply("[[", 1)%>%as.numeric()
alleles <- alleLinesub%>%str_trim()%>%strsplit(" ")%>%unlist(use.names=FALSE, recursive=FALSE)%>%strsplit(":")%>%laply("[[", 1)%>%as.factor()
allalleFreq <- cbind(allalleFreq, allefreqs)
}

allmigRates <- data.frame(allmigRates)
allinbRates <- data.frame(allinbRates)
allalleFreq <- data.frame(allalleFreq)
names(allmigRates)<- runs
names(allinbRates) <- runs
names(allalleFreq) <- runs

pdf(paste0(runDir, runs[1], "/corrplot_Inbreeding.pdf"), width=10, height=10)
pairs(allinbRates, pch=20, main="Inbreeding Coefficients", cex=0.5)
dev.off()
pdf(paste0(runDir, runs[1], "/corrplot_Migration.pdf"), width=10, height=10)
pairs(allmigRates, pch=20, main="Migration Rates", cex=0.5)
dev.off()
pdf(paste0(runDir, runs[1], "/corrplot_AlleleFreq.pdf"), width=10, height=10)
pairs(allalleFreq, pch=20, main="Allele Frequencies", cex=0.5)
dev.off()

migDF <- data.frame(migNames, allmigRates)
write.csv(migDF, paste0(runDir, runs[1], "/migration_rates.csv"), row.names=FALSE)

bestrun <- matrix(migDF[which.min(BayesDev)+1][[1]], nrow=popns,ncol=popns, byrow=TRUE)
bayespopnindex <- grep("Population Index", x)
bayespopnindex<-x[(bayespopnindex+2)]
BI<- bayespopnindex%>%str_trim()%>%strsplit(" ")%>%unlist(use.names=FALSE, recursive=FALSE)%>%strsplit("->")%>%laply("[", 2)%>%as.numeric()
bestrunnames <- CouncilIndex[BI]
bestrun<-data.frame(bestrun, row.names=bestrunnames)
names(bestrun)<-bestrunnames
write.csv(bestrun, paste0(runDir, runs[1], "/migration_rates_best.csv"), row.names=TRUE)


#############################
# calculateDeviance.R; Patrick G. Meirmans (2013) Non-convergence in Bayesian estimation of migration rates, Molecular Ecology Resources; Supplementary Material

# This script will calculate the Bayesian Deviance from the output of BayesAss (Wilson & Rannala 2003)
# For more information on the use of the deviance, see Faubet et al. (2007)
