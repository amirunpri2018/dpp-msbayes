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


#### UPDATE HERE, this is not true any more
#### Note that the following two vector variables are used in getData()
#    If the outputs of sumstatsvector (including ordering) is changed,
#    you can modify these two vectors accordingly, and everything should work.
#
# The first several columns contain the parameter values used for the
# simulations.  These values are sampled from the prior distributions.
# If you change the name here, you need to change a few statements
# in stdAnalysis()
# for unconstrained
#params.from.priorDistn.orig <- c("Psi", "var.t", "E.t", "omega")
# constrained with 2 psi
#params.from.priorDistn2 <- c("Tau1", "Tau2", "Psi1", "Psi2", "Psi", "var.t", "E.t", "omega")

# If there are 9 taxon pairs, there are 9 columns of "pi.b" stats for
# each taxon pairs, then 9 columns of "pi.w", ... etc.
# If you change the names here, don't forget to update the default value
# of used.stats in stdAnalysis
#summary.stat.names <- c("pi.b", "pi.w", "pi", "wattTheta", "pi.net", "tajD", "tajD.denom", "pi.wPop2", "pi.wPop1", "wattTheta.Pop2", "wattTheta.Pop1", "tajD.denomPop2", "tajD.denomPop1", "ShannonsIndex.Between", "ShannonsIndex.Net", "ShannonsIndex.Pop1", "ShannonsIndex.Pop2"  )

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
                        return.res=F
                        ) {
#  if(constrained) {
#    params.from.priorDistn <- params.from.priorDistn2
#  } else {
#    params.from.priorDistn <- params.from.priorDistn.orig  
#  }
  
  simDat <- getData(sim.infile)
  if(pre.rejected) {
    priorHeader <- scan(prior.infile,what="character", nlines=1, quiet=T)
    if (length(priorHeader) != length(simDat$prior.names)) {
      cat("ERROR: prior file doesn't have the correct number of columns\n")
      return(NA)
    } else if (any(priorHeader != simDat$prior.names)) {
      cat("ERROR: prior file doesn't match simDat\n")
      return(NA)
    }
    prior.dat <- scan(prior.infile, skip=1, quiet=T)
    prior.dat <- data.frame(matrix(prior.dat, ncol=length(priorHeader), byrow=T))
    names(prior.dat) <- priorHeader
  }
  nPairs <- simDat[["numTaxonPairs"]]
  params.from.priorDistn <- simDat[["prior.names"]]
  summary.stat.names <- simDat[["summary.stats"]]
  simDat <- simDat[["dat"]]

  # construct the column names
  usedColNames <- as.vector(sapply(used.stats, paste,
                                   1:nPairs, sep="."))

  #load OBSERVED summary stat vector
  obsDat <-getData(obs.infile)
  if (obsDat[["numTaxonPairs"]] != nPairs) {
    cat("ERROR: The number of taxon pairs are not same between the\n      ",
        "observed data,", obsDat$numTaxonPairs, "pairs, and simulated data,",
        nPairs, "pairs.\n")
    return(NA)
  }
  obsDat <- obsDat[["dat"]]

  # acceptance/regression, .... ie  the meat
  # The 1st column is PRI.numTauClass, which should be removed from analysis
  result <- list(prior.names=
                 params.from.priorDistn[params.from.priorDistn != "PRI.numTauClass"])
  # The column tells the model of Psi
  if(simDat[1,"PRI.numTauClass"] != 0) { # constrained
    constrained <- T
    # get rid of PRI.Psi, which is always constant in constrained Psi
    result$prior.names <- result$prior.names[result$prior.names != "PRI.Psi"]
  } else {
    constrained <- F
  }
  
  # set min, max boundary of prior parameters, and verbose print message
  # PRI.omega >= 0, PRI.E.t >= 0, PRI.Psi > 0
  min.val <- list(PRI.omega = 0, PRI.E.t = 0, PRI.Psi = 1)
  max.val <- list(PRI.Psi = nPairs)
  verbose.print <- list(PRI.omega = "(=Var(t)/E(t))",
                        PRI.Psi="(= number of possible divtimes)",
                        PRI.E.t="(= E(t))")
  # PRI.Tau* >= 0
  tauNames <- params.from.priorDistn[grep("^PRI[.]Tau", params.from.priorDistn)]
  if(length(tauNames) > 0) {
    temp.val <- sapply(rep(0, length(tauNames)), list)
    names(temp.val) <- tauNames
    min.val <- c(min.val, temp.val)
  }
  # 1 <= PRI.Psi.* <= nPairs - (numTauClasss - 1)
  if(constrained) {
    psiNames <- params.from.priorDistn[grep("^PRI[.]Psi[.]", params.from.priorDistn)]
    temp.val <- sapply(rep(1, length(psiNames)), list)
    names(temp.val) <- psiNames
    min.val <- c(min.val, temp.val)
    # I could use length of psiNames below, instead of simDat[1,1] to get numTauClass
    temp.val <- sapply(rep(nPairs - (simDat[1,"PRI.numTauClass"]-1), length(psiNames)), list)
    names(temp.val) <- psiNames
    max.val <- c(max.val, temp.val)
    
    temp.val <- sapply(rep("(= number of possible divtimes)", length(psiNames)), list)
    names(temp.val) <- psiNames
    verbose.print <- c(verbose.print, temp.val)
  }
  
  # run makepdANY for each para
  for (i in 1:length(result$prior.names)) {
    thisPriorName <- result$prior.names[i]
    # might need to work on constrained vs unconstrained here
    temp <- list(makepdANY(as.vector(obsDat[1,usedColNames],mode="numeric"),
                           simDat[,thisPriorName], simDat[,usedColNames], tol,
                           rep(T,len=nrow(simDat)),rejmethod=rejmethod))
    names(temp) <- thisPriorName

    # absorbing boundary
    if ( !  is.null(min.val[[thisPriorName]])) {
      temp[[thisPriorName]]$x[which(temp[[thisPriorName]]$x<min.val[[thisPriorName]])] <- min.val[[thisPriorName]]
    }
    if ( !  is.null(max.val[[thisPriorName]])) {
      temp[[thisPriorName]]$x[which(temp[[thisPriorName]]$x>max.val[[thisPriorName]])] <- max.val[[thisPriorName]]
    }
    
    result <- c(result, temp)
  }

  cat("######### results #######\n")
  # make print loop for Mode and quantile
  for (i in 1:length(result$prior.names)) {
    thisPriorName <- result$prior.names[i]
    name.rm.PRI <- sub("PRI[.]", "", thisPriorName)
    if(! is.null(verbose.print[[thisPriorName]])) {
      additional.print <- verbose.print[[thisPriorName]]
    } else {
      additional.print <- ""
    }
    
    cat ("### Mode estimation for ", name.rm.PRI, " ", additional.print, "\n");
    res.mode <- loc1stats((result[[thisPriorName]])$x,prob=0.95)[1]
    cat ("# Mode\n")
    print(res.mode)
    cat ("# 95 % quantile for ", name.rm.PRI, " ", additional.print, "\n");
    print(quantile((result[[thisPriorName]])$x,prob=c(0.025,0.975)))
  }

  # Print out figures
  pdf(pdf.outfile, width=7.5, height=10, paper="letter")
  for (i in 1:length(result$prior.names)) {
    thisPriorName <- result$prior.names[i]
    name.rm.PRI <- sub("PRI[.]", "", thisPriorName)
    if(! is.null(verbose.print[[thisPriorName]])) {
      additional.print <- verbose.print[[thisPriorName]]
    } else {
      additional.print <- ""
    }

    this.title <- paste(name.rm.PRI, additional.print, sep=" ")
    old.mfcol <- par()$mfcol
    par(mfcol=c(3,1))  # 3 plots per page
    if (pre.rejected) {
      make.hist(prior.dat[,thisPriorName],result[[thisPriorName]], title=this.title, breaks=20)
      plot.bf(prior.dat[,thisPriorName],result[[thisPriorName]]$x, main=this.title)
    } else {
      make.hist(simDat[,thisPriorName],result[[thisPriorName]],title=this.title,breaks=20)
      plot.bf(simDat[,thisPriorName],result[[thisPriorName]]$x,main=this.title)
    }
    par(mfcol=old.mfcol)
  }
  
  # pdf("CrustCOITheta60Na30ssPibPi_jointdots7_18.pdf")
  #plot(result11$x,result22$x,xlim=c(0,0.7),ylim=c(0,2.0),lty=2,lwd=0.5)
  plot((result[["PRI.omega"]])$x,(result[["PRI.E.t"]])$x,lty=2,lwd=0.5)
  # CHECK THIS

  plotKernDensity(result[["PRI.omega"]],result[["PRI.E.t"]], xlab="Omega",
                  ylab="E(t)", title="Omega and E(t)")

  # this plot doesn't seem to work.
  ## pdf("Skink0.5Milltol0.002Na0.5.pdf") 
  #  loc2plot(result.omega$x,result.Psi$x,cprob=0.6,alpha=0.4,xlab="omega", ylab="Psi")
  #  points(result.omega$x,result.Psi$x,pch=1,cex=0.7,col="darkgray")
  ## dev.off()

  cat("######### end of results #######\n")
  dev.off()

  if (return.res)
    return (list(nPairs=nPairs, simDat=simDat, obsDat=obsDat, result=result))
  
}


# This function takes an file name as an argument, read in the data
# from the file, and assign appropriate names.
# Returns a list(dat, numTaxonPairs)
# dat is the data.frame, numTaxonPairs is an integer indicating the number
# of taxon pairs.

getData <- function (infile) {
  first.line <- scan(infile, what="character", nlines=1, quiet=T) #header
  dat <- scan(infile, skip=1, quiet=T)
  dat <- data.frame(matrix(dat, ncol=length(first.line), byrow=T))
  names(dat) <- first.line  # assign the column names to the data.frame

  prior.names <- first.line[grep("^PRI[.]", first.line)]
  num.prior <- length(prior.names)
  
  # sum stats column-names (header) have the following form
  #   c("pi.b.1", "pi.b.2", "pi.b.3", "pi.w.1", "pi.w.2", "pi.w.3", ...)
  # Here, I'm getting rid of .digits part and taking uniue names.
  sum.stat.names <- unique(sub("[.][0-9]+$", "", first.line[(num.prior+1):length(first.line)],fixed=F))  

  # number of taxon pairs can be calculated from
  nTaxPairs <-
    (ncol(dat) - num.prior) / (length(sum.stat.names))

  return (list(dat=dat, numTaxonPairs=nTaxPairs, prior.names=prior.names, summary.stats=sum.stat.names))
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

make.hist <-function(vect, res.makepd, title="", xlim, ...) {
  #old.mfcol <- par()$mfcol
  #par(mfcol=c(2,1))
  hist.col = "gray"
  if(missing(xlim)) {
    hist(vect,br=15,col=hist.col,prob=TRUE,main=paste(title, ": Prior Dist'n"),
       xlab=title, ...)
  } else {
    hist(vect,br=15,col=hist.col,xlim=xlim,prob=TRUE,main=paste(title, ": Prior Dist'n"),
       xlab=title, ...)
  }
  lines(density(vect,bw=0.1),lty=2)
  if(missing(xlim)) {
    hist(res.makepd$x,br=15,col=hist.col,prob=TRUE,
         main=paste(title, ": Posterior Dist'n "), xlab=title, ...)
  } else {
    hist(res.makepd$x,br=15,col=hist.col,xlim=xlim,prob=TRUE,
         main=paste(title, ": Posterior Dist'n "), xlab=title, ...)
  }
  lines(density(res.makepd$x,bw=0.1))
  #par(mfcol=old.mfcol)
}


plot.bf <- function(prior,posterior, ...) {
  bf.out <- make.bf.vect(prior,posterior)
  if (is.null(bf.out)) {
    return(NULL);
  }

  plot(bf.out$crit.vals, bf.out$bf, type="l",
       ylim=c(0,min(40, max(bf.out$bf) * 1.1)),
       xlab="Crit. val",ylab="Bayes Factor",...)
  line.type = 3
  abline (1, 0, lty=2)
  abline(3, 0, lty=line.type)
  abline(10, 0, lty=line.type)
  abline(1/3, 0, lty=line.type)
  abline(1/10, 0, lty=line.type)
}

# Compares two Model of x < crit.val and x >= crit.val, and
# calculate the bayes factor.
# It returns a data.frame with two columns
make.bf.vect <- function(prior, posterior, crit.vals = NULL, num.points=100) {

  prior.len <- length(prior)
  posterior.len <- length(posterior)
  if (prior.len * posterior.len == 0) {  # some error
    return(NULL);
  }

  if (is.null(crit.vals)) {
    temp.pos <- posterior[posterior != max(posterior) &
                          posterior != min(posterior)]
     post.range <- range(temp.pos)
     crit.vals <- seq(post.range[1], post.range[2], length.out=num.points)
  }
  
  bf.vect <- numeric (0)
  for(i in 1:length(crit.vals)) {
    prior.below <- length(which(prior < crit.vals[i]))
    posterior.below <- length(which(posterior < crit.vals[i]))
    if (prior.below == 0 || posterior.below == posterior.len) {
      this.bf <- Inf
    } else {
      this.bf <- posterior.below * (prior.len-prior.below)/
        ((posterior.len - posterior.below) * prior.below)
    }
    bf.vect <- c(bf.vect, this.bf)
  }

  return (data.frame(crit.vals=crit.vals, bf=bf.vect))
}
