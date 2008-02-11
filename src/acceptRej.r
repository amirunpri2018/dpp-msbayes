#8_24_07 still have to add the plots for new hyper parameter posteriors

#library and R code  and package stuff
source("make_pd2005.r")
source("loc2plot.r")
library(KernSmooth)
library(locfit)  # this need to be installed from CRAN 
# 1. To install locfit, download a package locfit_*.tar.gz from CRAN.
#    (replace * with the most recent version number of locfit, currently 1.5-3)
#    <http://cran.r-project.org/>
# 2. As a super user
#      R CMD INSTALL locfit_*.tar.gz
#    Or if you are using Mac OS-X, you can use sudo
#      sudo R CMD INSTALL locfit_*.tar.gz
# For more details, see <http://cran.r-project.org/doc/manuals/R-admin.html>
# There is a setion "Add-on packages" describing how to install packages.


#### Note that the following two vector variables are used in getData()
#    If the outputs of sumstatsvector (including ordering) is changed,
#    you can modify these two vectors accordingly, and everything should work.
#
# The first several columns contain the parameter values used for the
# simulations.  These values are sampled from the prior distributions.
# If you change the name here, you need to change a few statements
# in stdAnalysis()
# for unconstrained
params.from.priorDistn.orig <- c("Psi", "var.t", "E.t", "omega")
# constrained with 2 psi
params.from.priorDistn2 <- c("Tau1", "Tau2", "Psi1", "Psi2", "Psi", "var.t", "E.t", "omega")

# If there are 9 taxon pairs, there are 9 columns of "pi.b" stats for
# each taxon pairs, then 9 columns of "pi.w", ... etc.
# If you change the names here, don't forget to update the default value
# of used.stats in stdAnalysis
summary.stat.names <- c("pi.b", "pi.w", "pi", "wattTheta", "pi.net", "tajD", "tajD.denom", "pi.wPop2", "pi.wPop1", "wattTheta.Pop2", "wattTheta.Pop1", "tajD.denomPop2", "tajD.denomPop1", "ShannonsIndex.Between", "ShannonsIndex.Net", "ShannonsIndex.Pop1", "ShannonsIndex.Pop2"  )

#summary.stat.names <- c("pi.b", "pi.w", "pi", "wattTheta", "pi.net",
#                            "tajD", "tajD.denom", "pi.wPop2", "pi.wPop1", "wattTheta.Pop2", "wattTheta.Pop1", "tajD.denomPop2", "tajD.denomPop1")

# Just print out the stat names
printStatNames <- function(constrained=F) {
  cat("## params.from.priorDistn\n")
  if(! constrained) {
    cat(params.from.priorDistn.orig)
  } else {
    cat(params.from.priorDistn2)
  }
  cat("\n\n## summary.stat.names\n")
  cat(summary.stat.names, "\n")
}

##### The main function to do the standard analysis
# sim.infile and obs.infile are the file names of the input files
# usedStatNames: a list of stats used to summarize
#                These names should be from summary.stat.names
# return.res: If set to TRUE, the function returns a huge results in a list,
#             so you can analyze the results further.
# rejmethod: if True, it doesn't boether with the regression, and uses
#            simple rejection
stdAnalysis <- function(obs.infile, sim.infile, prior.infile, 
	                pdf.outfile="figs.pdf",
                        tol=0.002,
                        used.stats=c("pi","wattTheta","pi.net","tajD.denom"),
                        rejmethod=T, pre.rejected=F,
                        return.res=F, constrained=F
                        ) {
  if(constrained) {
    params.from.priorDistn <- params.from.priorDistn2
  } else {
    params.from.priorDistn <- params.from.priorDistn.orig  
  }
  
  simDat <- getData(sim.infile, params.from.priorDistn)
  if(pre.rejected) {
    prior.dat <- scan(prior.infile)
    prior.dat <- data.frame(matrix(prior.dat, ncol=length(params.from.priorDistn), byrow=T))
    names(prior.dat) <- params.from.priorDistn
  }
  nPairs <- simDat[["numTaxonPairs"]]
  simDat <- simDat[["dat"]]

  # construct the column names
  usedColNames <- as.vector(sapply(used.stats, paste,
                                   1:nPairs, sep="."))

  #load OBSERVED summary stat vector
  obsDat <-getData(obs.infile, params.from.priorDistn)
  if (obsDat[["numTaxonPairs"]] != nPairs) {
    cat("ERROR: The number of taxon pairs are not same between the\n      ",
        "observed data,", obsDat$numTaxonPairs, "pairs, and simulated data,",
        nPairs, "pairs.\n")
    return(NA)
  }
  obsDat <- obsDat[["dat"]]
  
  #acceptance/regression, .... ie  the meat
  print("Working on Step 1")
  result.omega <- makepdANY(as.vector(obsDat[1,usedColNames],mode="numeric"),
                        simDat[,"omega"], simDat[,usedColNames], tol,
                        rep(T,len=nrow(simDat)),rejmethod=rejmethod)
  print("Working on Step 2")
  result.Psi <- makepdANY(as.vector(obsDat[1,usedColNames],mode="numeric"),
                        simDat[,"Psi"], simDat[,usedColNames], tol,
                        rep(T,len=nrow(simDat)),rejmethod=rejmethod)
  print("Working on Step 3")
  result.E.t <- makepdANY(as.vector(obsDat[1,usedColNames],mode="numeric"),
                        simDat[,"E.t"], simDat[,usedColNames], tol,
                        rep(T,len=nrow(simDat)),rejmethod=rejmethod)
  
  if (constrained) {
    print("Working on Step 4")
    result.Tau1 <- makepdANY(as.vector(obsDat[1,usedColNames],mode="numeric"),
                             simDat[,"Tau1"], simDat[,usedColNames], tol,
                             rep(T,len=nrow(simDat)),rejmethod=rejmethod)
  
    print("Working on Step 5")
    result.Tau2 <- makepdANY(as.vector(obsDat[1,usedColNames],mode="numeric"),
                             simDat[,"Tau2"], simDat[,usedColNames], tol,
                             rep(T,len=nrow(simDat)),rejmethod=rejmethod)

    print("Working on Step 6")
    result.Psi1 <- makepdANY(as.vector(obsDat[1,usedColNames],mode="numeric"),
                             simDat[,"Psi1"], simDat[,usedColNames], tol,
                             rep(T,len=nrow(simDat)),rejmethod=rejmethod)
    
    print("Working on Step 7")
    result.Psi2 <- makepdANY(as.vector(obsDat[1,usedColNames],mode="numeric"),
                             simDat[,"Psi2"], simDat[,usedColNames], tol,
                             rep(T,len=nrow(simDat)),rejmethod=rejmethod)
  }
  
  # Transformations, absorbing boundary at 0.0
  result.omega$x[which(result.omega$x < 0)] <- 0
  result.Psi$x[which(result.Psi$x < 1)] <- 1
  result.E.t$x[which(result.E.t$x < 0)] <- 0
  result.Psi$x[which(result.Psi$x > nPairs)] <- nPairs

  if(constrained) {
    result.Tau1$x[which(result.Tau1$x < 0)] <- 0
    result.Tau2$x[which(result.Tau2$x < 0)] <- 0
    
    result.Psi1$x[which(result.Psi1$x < 1)] <- 1
    result.Psi2$x[which(result.Psi2$x < 1)] <- 1
    
    result.Psi1$x[which(result.Psi1$x > (nPairs-1))] <- (nPairs-1)
    result.Psi2$x[which(result.Psi2$x > (nPairs-1))] <- (nPairs-1)
  }

  # Print out the results
  cat("######### results #######\n")
  cat("### Mode estimation for Omega (=Var(t)/E(t))\n")
  res.mode <- loc1stats(result.omega$x,prob=0.95)[1]
  cat("# Mode\n")
  print(res.mode)
  cat("# 95 % quantile for Omega (=Var(t)/E(t))\n")  
  print(quantile(result.omega$x,prob=c(0.025,0.975)))
  
  cat("\n### Mode estimation for Psi (= # of possible divtimes)\n")
  res.mode <- loc1stats(result.Psi$x,prob=0.95)[1]
  cat("# Mode\n")
  print(res.mode)  
  cat("# 95 % quantile for Psi (= # of possible divtimes)\n")
  print(quantile(result.Psi$x,prob=c(0.025,0.975)))

  cat("\n### Mode estimation for E(t)\n")
  res.mode<-loc1stats(result.E.t$x,prob=0.95)[1]
  cat("# Mode\n")
  print(res.mode)  
  cat("# 95 % quantile for E(t)")
  print(quantile(result.E.t$x,prob=c(0.025,0.975)))

  if (constrained) {
    cat("\n### Mode estimation for Tau1\n")
    res.mode<-loc1stats(result.Tau1$x,prob=0.95)[1]
    cat("# Mode\n")
    print(res.mode)
    cat("# 95 % quantile for Tau1")
    print(quantile(result.Tau1$x,prob=c(0.025,0.975)))
    
    cat("\n### Mode estimation for Tau2\n")
    res.mode<-loc1stats(result.Tau2$x,prob=0.95)[1]
    cat("# Mode\n")
    print(res.mode)
    cat("# 95 % quantile for Tau2")
    print(quantile(result.Tau2$x,prob=c(0.025,0.975)))


    cat("\n### Mode estimation for Psi1 (= # of possible divtimes)\n")
    res.mode <- loc1stats(result.Psi1$x,prob=0.95)[1]
    cat("# Mode\n")
    print(res.mode)
    cat("# 95 % quantile for Psi1 (= # of possible divtimes in Tau1)\n")
    print(quantile(result.Psi1$x,prob=c(0.025,0.975)))
    
    cat("\n### Mode estimation for Psi2 (= # of possible divtimes)\n")
    res.mode <- loc1stats(result.Psi2$x,prob=0.95)[1]
    cat("# Mode\n")
    print(res.mode)
    cat("# 95 % quantile for Psi2 (= # of possible divtimes in Tau2)\n")
    print(quantile(result.Psi2$x,prob=c(0.025,0.975)))
  }


  # Print out figures
  pdf(pdf.outfile)

  if (pre.rejected) {
    if(constrained) {
      make.hist(prior.dat[,"Tau1"], result.Tau1, title="Most recent Tau")
      make.hist(prior.dat[,"Tau2"], result.Tau2, title="oldest Tau")
      make.hist(prior.dat[,"Psi1"], result.Psi1, title="#taxon pairs divergence Tau1")
      make.hist(prior.dat[,"Psi2"], result.Psi2, title="#taxon pairs divergence Tau2")
    }
    make.hist(prior.dat[,"omega"], result.omega, title="Omega (Var(t)/E(t))")
    make.hist(prior.dat[,"Psi"], result.Psi, title="Psi (# of possible divtimes)")
    make.hist(prior.dat[,"E.t"], result.E.t, title="E(t)")
  } else {
    if (constrained) {
      make.hist(simDat[,"Tau1"], result.Tau1, title="Most recent Tau")
      make.hist(simDat[,"Tau2"], result.Tau2, title="oldest Tau")
      make.hist(simDat[,"Psi1"], result.Psi1, title="#taxon pairs divergence Tau1")
      make.hist(simDat[,"Psi2"], result.Psi2, title="#taxon pairs divergence Tau2")
    }

    make.hist(simDat[,"omega"], result.omega, title="Omega (Var(t)/E(t))")
    make.hist(simDat[,"Psi"], result.Psi, title="Psi (# of possible divtimes)")
    make.hist(simDat[,"E.t"], result.E.t, title="E(t)")
  }
  # pdf("CrustCOITheta60Na30ssPibPi_jointdots7_18.pdf")
  #plot(result11$x,result22$x,xlim=c(0,0.7),ylim=c(0,2.0),lty=2,lwd=0.5)
  plot(result.omega$x,result.E.t$x,xlim=c(0,0.7),ylim=c(0,2.0),lty=2,lwd=0.5)
  # CHECK THIS

  plotKernDensity(result.omega,result.E.t,
                  title="Omega and E(t)")

  # this plot doesn't seem to work.
  ## pdf("Skink0.5Milltol0.002Na0.5.pdf") 
  #  loc2plot(result.omega$x,result.Psi$x,cprob=0.6,alpha=0.4,xlab="omega", ylab="Psi")
  #  points(result.omega$x,result.Psi$x,pch=1,cex=0.7,col="darkgray")
  ## dev.off()

  cat("######### end of results #######\n")
  dev.off()

  if (return.res)
    return (list(nPairs=nPairs, simDat=simDat, obsDat=obsDat, result.omega=result.omega, result.Psi=result.Psi, result.E.t=result.E.t))
  
}


# This function takes an file name as an argument, read in the data
# from the file, and assign appropriate names.
# Returns a list(dat, numTaxonPairs)
# dat is the data.frame, numTaxonPairs is an integer indicating the number
# of taxon pairs.
getData <- function (infile, params.from.priorDistn) {
#  dat <- read.table(infile, header=F)
  first.line <- scan(infile, nlines=1)
  dat <- scan(infile)
  dat <- data.frame(matrix(dat, ncol=length(first.line), byrow=T))
 # check this

  # number of taxon pairs can be calculated from
  nTaxPairs <-
    (ncol(dat) - length(params.from.priorDistn)) / (length(summary.stat.names))

  # Make the column Names
  # Now columnNames becomes a vector
  #   c("pi.b.1", "pi.b.2", "pi.b.3", "pi.w.1", "pi.w.2", "pi.w.3", ...)
  # for 3 taxon pair
  columnNames <- as.vector(sapply(summary.stat.names, paste,
                                  1:nTaxPairs, sep="."))

  columnNames <- c(params.from.priorDistn, columnNames)
  
  names(dat) <- columnNames  # assign the column names to the data.frame
  
  return (list(dat=dat, numTaxonPairs=nTaxPairs))
}

# -----------------------------------------------------------------------
# 2D Kernel density estimates: q1 X q2
# -----------------------------------------------------------------------
# horizontal view
plotKernDensity <- function (res1, res2, title="q1 and q2", xlab="q1", ylab="q2") {
  bwq1 <- dpik(res1$x)
  bwq2 <- dpik(res2$x)
  x <- cbind(res1$x, res2$x)
  est <- bkde2D(x, bandwidth=c(bwq1,bwq2))
  par(ask = FALSE)
  persp(est$x1, est$x2, est$fhat, theta = 145, phi = 25, col = "lightblue1",
        xlab=xlab, ylab =ylab,
        zlab = paste("Pr(", xlab, ", ", ylab, "| X)",sep=""),
        main = paste("Joint Density of", title), axes = TRUE, nticks = 5,
        ticktype = "detailed", ltheta = -135, lphi = 145, shade = TRUE)
}

make.hist <-function(vect, res.makepd, title="", xlim) {
  if(missing(xlim)) {
    hist(vect,br=15,col="black",prob=TRUE,main=paste(title, "Raw"),
       xlab=title)
  } else {
    hist(vect,br=15,col="black",xlim=xlim,prob=TRUE,main=paste(title, "Raw"),
       xlab=title)
  }
  lines(density(vect,bw=0.1),lty=2)
  if(missing(xlim)) {
    hist(res.makepd$x,br=15,col="black",prob=TRUE,
         main=paste(title, "Adj. "), xlab=title)
  } else {
    hist(res.makepd$x,br=15,col="black",xlim=xlim,prob=TRUE,
         main=paste(title, "Adj. "), xlab=title)
  }
  lines(density(res.makepd$x,bw=0.1))
}
