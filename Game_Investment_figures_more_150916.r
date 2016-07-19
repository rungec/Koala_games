

inpDir <- "W:/seq_genetics/analysis/collaboration_paper/3_game theory/outputs/Random dispersal rates/randomlygeneratedZvalues.csv"
outDir <- "W:/seq_genetics/analysis/collaboration_paper/figures/"

library(ggplot2)
library(reshape2)

a<-read.csv(inpDir, header=TRUE)
names(a)[1:6] <- c("LGA2 to LGA1","LGA3 to LGA1","LGA3 to LGA2","LGA1 to LGA2","LGA1 to LGA3","LGA2 to LGA3")
dataa <- melt(a[,1:6])

###################
#Make plots of the dispersal rates for the two-way interactions

p1 <- ggplot(data=dataa, aes(x=value, fill=variable))+
	theme_bw(17) + #get rid of grey bkg and gridlines
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
	coord_cartesian(xlim=c(0,max(dataa$value)), ylim=c(5,15))+ #set x and y limits
	labs(x="Annual dispersal rate", y="Frequency")+
	scale_fill_discrete(name="Dispersal from")+
	geom_density(alpha=0.25, colour=NA, adjust=0.8)
	#geom_density(alpha=0.125, colour="grey70", adjust=0.8, show_guide=FALSE)#get rid of the legend slashes

	
outPath <- paste0(outDir, "Density_plot_of_random_dispersal_rates.png")
ggsave(filename=outPath, width=10, height=7)
	
#for 3d figures:	
library(scatterplot3d)
library(rgl)
a<-read.csv(file.choose())
plot3d(a$bris,a$loga,a$redl)	
