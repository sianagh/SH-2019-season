# ----------------------------------------------------------------------------
# SUBJECT: Interim Southern Hemisphere VE paper - surveillance data plots
# AUTHOR: Sheena Sullivan
# DATE: 23/8/2019
# -----------------------------------------------------------------------------
library(mem)
library(plyr)
library(dplyr)

# ===================================================
# colours
rd <- rgb(1, 0, 0,0.4)
gn <- rgb(0, 1, 0,0.4)
bl <- rgb(0, 0, 1,0.4)
pp <- rgb(1, 0, 1,0.4)
gy1 <- rgb(0, 0, 0,0.05)
gy2 <- rgb(0, 0, 0,0.2)

# -----------------------------------------------------------------------------
# Functions
# -----------------------------------------------------------------------------
# plot the MEM curves
plotmem <- function(dat,y.lab="Weekly rate",this.year=2019,n.years=6, extraheight=0){  #not used in loop at end, but useful for plotting each jurisdiction's data
  this.year.col <- grep(this.year,colnames(dat))  #identifies this year's data - don't want this in the mem
  end.col <- this.year.col-1             #identifies the last column to be used by mem
  start.col <- this.year.col-n.years     #identifies the first column to be used by mem
  
  n.row <- nrow(dat)  #finds the number of rows of the data for the X axis
  
  dat <- dat[order(dat$week),]
  epi<-memmodel(dat[start.col:end.col])  #aligns the data for average of 5 years
  intensity <- memintensity(epi)   #intensity thresholds
  thresholds <- as.data.frame(intensity$intensity.thresholds)
  
  summary(epi)
  mean <- as.data.frame(epi$typ.curve)
  mean$V2[is.na(mean$V1)] <- NA
  mean$week <- as.numeric(dat$week)
  mean$thisyear <- dat[,this.year.col]
  # mean$V2[is.na(mean$V1) & is.na(mean$V3) & mean$V2 > lag(mean$V2)] <- NA
  ymax <- round_any(max(mean[,c("V2","thisyear")],na.rm=TRUE),5,ceiling)  #reduce size if planning to plot highest thresholds
  if(thresholds[3]>ymax){
    ymax <- as.numeric(thresholds[3])
  }

  # plot the mean and current year's data
  par(xpd=F)
  plot(NA, ylim=c(0,ymax), xlim=c(1,53), axes=F, 
       ylab = paste(y.lab), xlab="")
  axis(1, at=seq(1,53,2),cex.axis=1)
  axis(2, las=1, cex.axis=1)
  lines(mean$week,mean$V2, lty=1, lwd=2, col="orange")  #mean of past X years
  lines(mean$week,mean$thisyear, lty=2, lwd=2, col="brown")  #this year's data
  
  # baseline annotations
  points(min(row.names(mean)[mean$thisyear>=epi$epidemic.thresholds[1]],na.rm=T),
         epi$epidemic.thresholds[1],
         col="red",pch=1)
  points(max(row.names(mean)[mean$thisyear>=epi$epidemic.thresholds[2]],na.rm=T),
         epi$epidemic.thresholds[2],
         col="green",pch=1)
  
  # intensity annotations
  abline(h=thresholds[2],col=3,lty=2)
  text(x=46,y=thresholds[2]+extraheight,"Medium (40%)", cex=0.8)
  abline(h=thresholds[3],col=3,lty=2)
  text(x=46,y=thresholds[3]+extraheight,"High (90%)", cex=0.8)
  if(thresholds[4]<ymax){
    abline(h=thresholds[4],col=3,lty=2)
    text(x=46,y=thresholds[4]+extraheight,"Very high (97.5%)", cex=0.8)
  }
}
# plotmem(aust.ili)

plotswab <- function(swab.dat){
  ymax.swab <- round_any(max(swab.dat[,c(3:length(swab.dat))],na.rm=TRUE),5,ceiling)  #reduce size if planning to plot highest thresholds
  # swab.dat <- aust.ve
  plot(noncases~week, data=swab.dat, type="n", ylim=c(0,ymax.swab),  axe=F, 
       xlab = "", ylab = "No. cases")
  # rect(18,0,33,160, col=gy1,lty=0)
  polygon(c(0:53),c(NA,swab.dat$noncases,NA), col=gy1,lty=0)
  polygon(c(0:53),c(NA,swab.dat$fluahx,NA), col=rd,lty=0)
  polygon(c(0:53),c(NA,swab.dat$fluah1,NA), col=bl,lty=0)
  polygon(c(0:53),c(NA,swab.dat$fluah3,NA), col=pp,lty=0)
  polygon(c(0:53),c(NA,swab.dat$flub,NA), col=gn,lty=0)
  axis(1, at=seq(1,53,2),cex.axis=1)
  axis(side=2, las=1, cex.axis=1)
}

library(RColorBrewer)
plotswab2 <- function(swab.dat){
  ymax.swab <- round_any(max(swab.dat[,c(3:length(swab.dat))],na.rm=TRUE),5,ceiling)  #reduce size if planning to plot highest thresholds
  # swab.dat <- aust.ili.swab
  barplot(as.matrix(t(swab.dat[,5:2])), col=brewer.pal(4, "PuRd"),
          legend.text = NULL, axes = F)
  axis(1, at=seq(1,53,2),cex.axis=1)
  axis(side=2, las=1, cex.axis=1)
}

#==========================
# Australia
#==========================
# Aust ILI
aust.ili <- read.csv("aust-ili.csv")
plotmem(aust.ili)

aust.ili.swab <- read.csv("aust-ili-swab.csv")
plotswab2(aust.ili.swab)


# Aust SARI
aust.sari <- read.csv("aust-sari.csv")
plotmem(aust.sari)

aust.sari.swab <- read.csv("aust-sari-swab.csv")
plotswab2(aust.sari.swab)


# ============================================================================
# Chile
# ============================================================================
# Load packages
chile.sari <- read.csv("chile-sari.csv")
plotmem(chile.sari)

chile.swab <- read.csv("chile-swab.csv")
plotswab2(chile.swab)


# ============================================================================
# New Zealand
# ============================================================================
nz.ili <- read.csv("nz-ili.csv")
plotmem(nz.ili)

nz.ili.swab <- read.csv("nz-ili-swab.csv")
nz.sari.swab <- nz.sari.swab[order(nz.sari.swab$week),]
plotswab2(nz.ili.swab)

# SARI
nz.sari <- read.csv("nz-sari.csv")
plotmem(nz.sari)

nz.sari.swab <- read.csv("nz-sari-swab.csv")
colSums(nz.sari.swab)


# ============================================================================
# South Africa
# ============================================================================
rsa.ili <- read.csv("rsa-ili.csv")
plotmem(rsa.ili)

rsa.ili.swab <- read.csv("rsa-ili-swab.csv")
plotswab(rsa.ili.swab)

# ppos for rsa
plot(rsa.ili.swab$fluah3, type="l")
lines(rsa.ili.swab$noncases)
lines(with(rsa.ili.swab,fluah3/(fluah3+noncases))*100,col="red")
(with(rsa.ili.swab,fluah3/(fluah3+noncases)))

# ===================================================
# Plot all countries
# ===================================================
countries <- c("Australia \nPrimary care","Australia \nHospitals","Chile \nHospitals","New Zealand \nPrimary care",
               "New Zealand \nHospitals","South Africa \nPrimary care")
dat.names <- list(aust.ili,aust.sari,chile.sari,nz.ili,nz.sari,rsa.ili)
swab.dat.names <- list(aust.ili.swab,aust.sari.swab,chile.swab,nz.ili.swab,nz.sari.swab,rsa.ili.swab)
ytext <- c("ILI rate \n(per 1000)","Hospitalisations \n(per 1000 beds)","SARI rate \n(per 100)",
           "ILI rate \n(per 100,000)","SARI rate \n(per 100,000)","P&I presentations \n(per 100,000)")
num.years <- c(6,6,6,6,6,6)
# ===================================================
dev.off()
# windows(8,10)
pdf("surveillance.pdf", 8,10)
layout(matrix(c(1:21), 7, 3, byrow=T), widths=c(1,5,5), heights=rep(c(3),7)) # set up layout for multiple plots  
for(i in 1:length(dat.names)){
  par(mar=c(0,0,0,0))
  plot.new()
  mtext(paste0(countries[i]), line=-3, adj=0, cex=0.8, font=2)

  # Surveillance data
  dat <- as.data.frame(dat.names[i])
  dat <- subset(dat,week != 53)
  par(mar=c(1,6,2,1))
  plotmem(dat, y.lab=ytext[i],n.years=num.years[i])
  par(xpd=NA)

  # Swab data
  swab.dat <- as.data.frame(swab.dat.names[i])
  plotswab2(swab.dat)
}
par(mar=c(0,0,3,0))
plot.new()
plot.new()
mtext("Week", cex=0.7,line=-0.5,adj=0.5)
legend("left",legend=c("Pre-season threshold","Post-season threshold"), 
       col=c("red","green"), pch=1, pt.cex=1, horiz=T, border=NULL, box.lty=0, cex=1)
legend("bottomleft",legend=c("2013-2018 mean","2019"), 
       col=c("orange","brown"), lwd=c(2,3), lty=c(1,2), horiz=T, border=NULL, box.lty=0, cex=1)
plot.new()
mtext("Week", cex=0.7,line=-0.5,adj=0.5)
legend("bottom",c("A(H1N1)pdm09","A(H3N2)", "A(Hx)","B"), col=rev(brewer.pal(4, "PuRd")), pch=19, border=NULL, box.lty=0, cex=1)
# legend("bottomright",legend=c("A(Hx)","B","Negative"),
#          col=c(rd,gn,gy1), horiz=T, pch=19, border=NULL, box.lty=0, cex=1)
# legend("right",legend=c("A(H1N1)pdm09","A(H3N2)"),
#        col=c(bl,pp), horiz=T, pch=19, border=NULL, box.lty=0, cex=1)
dev.off()



