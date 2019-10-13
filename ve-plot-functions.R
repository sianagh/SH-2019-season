
# --------------------------------------------------
#  Functions and constants for the forest plots 
# --------------------------------------------------

library("plyr")
library("dplyr")
library("metafor")

# ====================================================================
# Y-axis functions
# these functions are used to identify the spacings and 
# break points for populations in the forest plots
# ====================================================================

# returns a vector of values for the y positions of each population group
get.y.pos <- function(dfn.in){
  dfn <- dfn.in
  for (i in unique(dfn.in$pop.order)) {  #add two rows for each population group
    dfn <- rbind.fill(dfn, 
                      data.frame(studygroup="",pop.order = i,stringsAsFactors = F),
                      data.frame(studygroup="",pop.order = i,stringsAsFactors = F))  #repeat to get two rows
  }  
  pops <- table(dfn$pop.order)           #number of estimates for each population
  y.pos <- cumsum(pops)                  #cumulative sum - to get break points
  return(y.pos)  
}
get.y.pos.extra <- function(dfn.in){
  dfn <- dfn.in
  for (i in unique(dfn.in$pop.order)) {  #add two rows for each population group
    dfn <- rbind.fill(dfn, 
                      data.frame(studygroup="",pop.order = i,stringsAsFactors = F),
                      data.frame(studygroup="",pop.order = i,stringsAsFactors = F),
                      data.frame(studygroup="",pop.order = i,stringsAsFactors = F),
                      data.frame(studygroup="",pop.order = i,stringsAsFactors = F))  #repeat to add rows between groups
  }  
  pops <- table(dfn$pop.order)           #number of estimates for each population
  y.pos <- cumsum(pops)                  #cumulative sum - to get break points
  return(y.pos)  
}
# returns a vector of numbers which represent the range for each group of estimates

get.rows <- function(dfn.in, y.pos){
  true.pops <- table(dfn.in$pop.order)   #get true number of estimates per population
  row.end <- y.pos-2                     #for each population, this is the end row, -2 because of the added rows
  row.start <- row.end                   #for each population, this is the start row (using loop below)
  for(i in 1:length(row.end)){                         
    row.start[i] <- row.end[i]-true.pops[i]+1  
  }
  rows<-vector()
  for (i in 1:length(row.end)){
    rows <- c(rows,row.start[i]:row.end[i])
  }
  return(rows)
}


# returns a vector of the population group names
get.pop.names <- function(dfn){
  x <- factor(dfn$pop)
  levels(x)
  pop.names <- levels(x)
  return(pop.names)
}


# returns a vector of the limits for the y-axis
get.y.lim <- function(y.pos){
  ystart <- -1
  yend <- max(y.pos)+2 #add 2 rows for labels at the top of the figure
  y.lim <- c(ystart,yend)
  return(y.lim)
}


# ====================================================================
# Forest plot function - no summary estimates
# ====================================================================
plotme <- function(dat){
  dat <- dat[with(dat, order(dat$study)),]
  y.pos<-get.y.pos.extra(dat)-2        # get y limits
  rows <- get.rows(dat, y.pos=y.pos) # get row positions
  pop.names <- levels(droplevels(dat$pop))
  y.lim <- get.y.lim(y.pos=y.pos) # y scale
  head.pos <- y.lim[2]-1          # y point for heading
  vax.lab.pos <- vax.lab.pos  #reset position of vaccination labels
  x.start <- x.start
  dat <- dat[with(dat, order(dat$population)),]
  dat <- dat[order(dat$pop.order,-as.numeric(dat$setting),-dat$group.order),]
  
  forest(dat$ve, ci.lb=dat$ll, ci.ub=dat$ul, 
         xlim=x.lim, alim=a.lim, at=a.lim.at, cex.axis=0.6, cex=0.6,
         refline=c(0,50), slab=c(as.character(dat$study)), 
         ilab=cbind(as.character(dat$setting),dat$vc,dat$uvc,dat$vnc,dat$uvnc), ilab.xpos=case.pos,
         ylim=y.lim,  
         rows=rows,
         xlab="Vaccine effectiveness",
         digits=0
  )
  
  # labels along top of plot
  op <- par(cex=.7, font=4)
  text(x.start, y.pos-1, pop.names, pos=4)  #population group names
  par(font=2)
  text(case.pos, head.pos, case.lab)        #"Flu+", "Flu-", "Flu+", "Flu-"
  text(vax.lab.pos, head.pos+1, vax.lab)  #"Vaccinated", "Unvaccinated"
  text(x.start, head.pos, auth.lab, pos=4)  #" Network"
  text(ve.lab.pos, head.pos, ve.lab, pos=2) #"VE [95% CL]"
  par(op)
}



# ============================================================
# Forst plot function with summary estimates
# ============================================================
changescale <- function(res){
  res.c<-res
  res.c$beta <- (1-(res$b))*100  
  res.c$ci.lb <- (1-(res$ci.ub))*100
  res.c$ci.ub <- (1-(res$ci.lb))*100
  return(res.c)
}

getse <- function(dat){
  dat$lower <- 1-dat$ul/100
  dat$upper <- 1-dat$ll/100
  dat$se <- (dat$upper-dat$lower)/3.92
  return(dat$se)
}

# dat<-all.data
plotwithsumm <- function(dat,virus){
  dat$or <- 1-dat$ve/100
  dat$se <- getse(dat)
  y.pos.res<-(get.y.pos.extra(dat))
  y.pos.res <- c(0,y.pos.res[-length(y.pos.res)])
  ages <- levels(droplevels(dat$pop))
  
  par(mar=c(5,1,2,0))
  plotme(dat)
  op <- par(font=3)
  for(i in 1:length(ages)){
    sb.dat <- subset(dat,pop==ages[i])
    if(nrow(sb.dat)>1){
      res <- changescale(rma(yi=or, sei=se, data=sb.dat, measure="OR", method="DL"))
      addpoly(res, row=y.pos.res[i], cex=0.7, col="blue",
              mlab=
                # expression(paste("RE Model (I"^"2","=",round(res$I2, digits=3),
                #                  tau^"2", "=",formatC(res$tau2, digits=2),")")))
                paste0("RE Model (I-squared: ", (formatC(res$I2, digits=3)),
                       "%; tau-squared: ",(formatC(res$tau2, digits=2)),")"))
      
    }else{
      text(x.lim[1],y.pos.res[i], paste0(""), cex=0.7, pos=4)
    }
  }
  # title(paste(virus))
  par(op)
}