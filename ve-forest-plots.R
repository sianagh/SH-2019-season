# ====================================================================
# Scripts to plot the VE estimates from southern hemisphere countries
# DATE: 20190904
# AUTHOR: Sheena Sullivan
# ====================================================================
# Load packages
if(!require("plyr")) install.packages("plyr", repos="http://cran.us.r-project.org")
if(!require("dplyr")) install.packages("dplyr", repos="http://cran.us.r-project.org")  #rbind.fill function
if(!require("metafor")) install.packages("metafor", repos="http://cran.us.r-project.org")

# ====================================================================
# plot functions
source("ve-plot-functions.R")  

# import data
load("shvedata.Rdata")
with(sh.data,table(study,typesubtype))
all.data <- subset(sh.data, typesubtype=="any")
ah1.data <- subset(sh.data, typesubtype=="ah1")
ah3.data <- subset(sh.data, typesubtype=="ah3")
b.data <- subset(sh.data, typesubtype=="b_")

# output vaccination coverage
write.csv(all.data[order(all.data$study,all.data$pop),c("study",'pop',"vncp")], file="coverage.csv")

# ====================================================================
# plot parms
# ====================================================================
a.lim=c(-20,100)           #limits for the x axis of the forest plots
a.lim.at=c(-20,0,50,100)   #values shown on the forest plots
x.lim <- c(-220,150)       #range of the x-axis
x.start <- x.lim[1]        #start point of the x-axis
case.lab <- c("Setting","V", "UV", "V", "UV") #sample size headers
case.pos <- c(-150,-120, -100, -70, -50)           #position of the study sample sizes and headers
vax.lab <- c("Postiive", "Negative")#vaccination status header
vax.lab.pos <- c(-110,-60) #position of the vaccination status headers
auth.lab <- "Network"      #Study group header
ve.lab.pos <- 150          #position of the VE estimates and header (right side of plot)
ve.lab <- "VE [95% CL]"    #header for the VE estimates
# ====================================================================

# ====================================================================
# Forest plots without summary estimates
# ====================================================================
par(mar=c(3,0,2,0))
# pdf("ve-forest-nosumm.pdf",7,10)
# windows(12,15)
# layout.show(
layout(matrix(c(1:3), 3,1), heights=c(5,5.6,5))
  # ) # set up layout for multiple plots  
plotme(ah1.data)
mtext("Influenza A(H1N1)pdm09", side=3, line=1, cex=0.8, adj=0, col="blue", font=2)
plotme(ah3.data)
mtext("Influenza A(H3N2)", side=3, line=1, cex=0.8, adj=0, col="blue", font=2)
plotme(b.data)
mtext("Influenza B", side=3, line=1, cex=0.8, adj=0, col="blue", font=2)
mtext("Notes: Estimates from REVELAC-I only include patients in a target group for vaccination; Australia uses adjuvanted TIV for elderly aged 65+",
      side=1, line=4, cex=0.5, adj=0)
dev.off()


# ====================================================================
# Forest plots with summary estimates
# ====================================================================
par(mar=c(3,0,2,0))
# windows(12,15)
pdf("ve-forest-summ.pdf",7,10)
# layout.show(
  layout(matrix(c(1:3), 3,1), heights=c(5,5.6,4))
  # ) # set up layout for multiple plots  
plotwithsumm(ah1.data, "H1")
mtext("Influenza A(H1N1)pdm09", side=3, line=1, cex=0.8, adj=0, col="blue", font=2)
par(mar=c(3,0,2,0))
plotwithsumm(ah3.data,"H3")
mtext("Influenza A(H3N2)", side=3, line=1, cex=0.8, adj=0, col="blue", font=2)
plotwithsumm(b.data,"B")
mtext("Influenza B", side=3, line=1, cex=0.8, adj=0, col="blue", font=2)
mtext("Notes: Estimates from REVELAC-I only include patients in a target group for vaccination; Australia uses adjuvanted TIV for elderly aged 65+",
      side=1, line=4, cex=0.5, adj=0)
dev.off()


# ====================================================================
# Funnel plots 
# ====================================================================
dat.names <- list(ah1.data,ah3.data,b.data)
viruses <- c("ah1","ah3","b_","any")
pdf("otherplots.pdf")
for(i in 1:3){
  dat <- as.data.frame(dat.names[i])
  dat <- subset(dat,pop=="All patients")
  dat$or <- 1-dat$ve/100
  dat$se <- getse(dat)
  # random effects
  res.re <- (rma(yi=or, sei=se, data=dat, measure="OR", method="DL"))
  par(mar=c(3,3,3,1))
  layout(matrix(c(1:3),3,1))
  forest(res.re, slab=paste(dat$study,dat$setting, dat$pop), main = paste(viruses[i],"- Random-Effects Model"))
  funnel(res.re, main="Funnel plot")
  qqnorm(res.re, main="QQ plot")
  # mixed effects
  res.me <- rma(yi=or, sei=se, mods = cbind(study), data = dat)
  forest(res.me, slab=paste(dat$study,dat$setting, dat$pop), main = paste(viruses[i],"- Mixed-Effects Model"))
  funnel(res.me, main="Funnel plot")
  qqnorm(res.me, main="QQ plot")
}
dev.off()
