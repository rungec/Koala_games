#This script makes figures from the output of the Koala Games project - Random dispersal rates
#100 sets of random dispersal rates are generated for the three LGAs from a beta(0.5, 0.5) distribution
#these have different assymetry coefficients (0 assymetric -1 symmetric migration between 3LGAs)
#The optimal response of each LGA is generated for each of these sets of random dispersal rates
#for two scenarios: Benevolent social planner (overall target), and non-collaborative game (individual target for each LGA)
#The maximum number of koalas is also calculated for each set of random dispersal rates for a spend of c(10, 10, 10) $m

options(stringsAsFactors=FALSE)
library(ggplot2)
library(gridExtra) #for grid.arrange()

wd <- "C:/Claire/GPEM_Postdoc/2_SEQ_Koala_Collaboration/3_game theory/outputs/Random dispersal rates/Run 1/"
setwd(wd)
outDir <- "C:/Claire/GPEM_Postdoc/2_SEQ_Koala_Collaboration/figures/"


maxspend <- read.csv("ROI_maxspend_100randomZ.csv")
maxspend$symmcoefBin <- as.factor(rep(1:10, each=10))
Z <- read.csv("randomlygeneratedZvalues.csv")

maxspendLong <- data.frame(symmcoef=rep(maxspend$symmcoefBin, 3), LGA=rep(c("bris", "loga", "redl"), each=100), decline=c(maxspend$decline_bris, maxspend$decline_loga, maxspend$decline_redl)*100)
#, invest=c(maxspend$bris_invest, maxspend$loga_invest, maxspend$redl_invest))

#plot variation in population decline for different assymetries and maximim spend of c(10, 10, 10) $m
p1 <- ggplot(maxspendLong, aes(symmcoef, decline, color=LGA))+
	geom_point()
	
p2 <- ggplot(maxspendLong, aes(symmcoef, decline, color=LGA))+
	geom_boxplot()
	
outPath <- paste0(outDir, "Koala_games_Randomdispersal_boxplotofdeclines_byLGA.png") 
png(filename=outPath, width=960, height=480)
grid.arrange(p1, p2, ncol=2)
dev.off()
	
#plot variation in overall outcome under maximum spend for different assymetries
p3 <- ggplot(maxspend, aes(symmcoefBin, decline_overall))+
		geom_point()

p4 <- ggplot(maxspend, aes(symmcoefBin, decline_overall))+
		geom_boxplot()

outPath <- paste0(outDir, "Koala_games_Randomdispersal_boxplotofdeclines_overall.png") 
png(filename=outPath, width=960, height=480)
grid.arrange(p3, p4, ncol=2)
dev.off()