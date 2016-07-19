#Make figures for koala games project
library(reshape)
library(ggplot2)
library(gridExtra)

inpDir <- "W:/seq_genetics/analysis/collaboration_paper/3_game theory/outputs_temp/"
outDir <- "W:/seq_genetics/analysis/collaboration_paper/figures/"

S1 <- read.csv(paste0(inpDir, "InvestmentBenefit_Scenario1.csv"), stringsAsFactors=FALSE)
S2 <- read.csv(paste0(inpDir, "InvestmentBenefit_Scenario2.csv"), stringsAsFactors=FALSE)

################
#Preliminary processing
################

dat <- rbind(melt(S2, id.vars="State", measure.vars=c("brisD", "logaD", "redlD")), melt(S1, id.vars="State", measure.vars=c("brisD", "logaD", "redlD")))
	dat$Scenario <- as.factor(rep(c("Dispersal", "No dispersal"), each=3*nrow(S1)))
	names(dat)[1:3] <- c("Investment", "Council", "Decline")

dat2 <- rbind(melt(S2, id.vars="State", measure.vars=c("overallD")), melt(S1, id.vars="State", measure.vars=c("overallD")))
	dat2$Scenario <- as.factor(rep(c("Dispersal", "No dispersal"), each=1*nrow(S1)))
	names(dat2)[1:3] <- c("Investment", "Council", "Decline")


###############################
#Plot of overall ROI
p1 <- ggplot(dat2, aes(x=Investment, y=Decline*100, linetype=Scenario)) + 
		geom_line(lwd=1) +
		theme_classic(17) + #remove grey #remove grids
		#theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
		labs(x="Investment ($m)", y="Percent decline")+
		theme(axis.title.x = element_text(vjust=-0.6),axis.title.y = element_text(vjust=1))+
		theme(legend.key = element_blank()) + #remove boxes from around legend items
		theme(legend.justification=c(1,1), legend.position=c(.95, 1))

#Plot of ROI by council
p2 <- ggplot(dat, aes(x=Investment, y=Decline*100, col=Council, linetype=Scenario)) + 
		geom_line(lwd=1) +
		theme_classic(17) + #remove grey #remove grids
		#theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
		labs(x="Investment ($m)", y="Percent decline")+
		theme(axis.title.x = element_text(vjust=-0.6),axis.title.y = element_text(vjust=1))+
		theme(legend.key = element_blank()) + #remove boxes from around legend items
		theme(legend.justification=c(1,1), legend.position=c(.95, 1))+
		scale_color_discrete(labels=c("Brisbane", "Logan", "Redlands"))+
		guides(linetype=FALSE)#turn off the legend for line type

outPath <- paste0(outDir, "Koala_games_ROI_DispersalvsNoDispersal.png") 
png(filename=outPath, width=960, height=480)
grid.arrange(p1+coord_cartesian(ylim=c(15,75)), p2+coord_cartesian(ylim=c(15,75)), ncol=2)
dev.off()

#plot above with zoom
outPath <- paste0(outDir, "Koala_games_ROI_DispersalvsNoDispersal_zoom.png") 
png(filename=outPath, width=960, height=480)
grid.arrange(p1+coord_cartesian(xlim=c(0,10), ylim=c(15,75)), p2+coord_cartesian(xlim=c(0,10), ylim=c(15,75)), ncol=2)
dev.off()



#outPath <- paste0(outDir, "Koala_games_ROI.png")
#	ggsave(filename=outPath)

