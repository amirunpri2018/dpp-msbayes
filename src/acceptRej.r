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

# The first several columns of simDat contain the parameter values
# used for the simulations.  These values are sampled from the prior
# distributions.
# Example:
# When Psi constrainted to be 3
# PRI.numTauClass PRI.Tau.1       PRI.Tau.2       PRI.Tau.3       PRI.Psi.1       
# PRI.Psi.2       PRI.Psi.3       PRI.Psi PRI.var.t       PRI.E.t PRI.omega

# If there are 9 taxon pairs (or 3 taxon pairs with 3 genes per taxon
# pair), there are 9 columns of "pi.b" stats for each taxon pairs,
# then 9 columns of "pi.w", ... etc.
#
# Examples of summary stats (may not be up to date)
# c("pi.b", "pi.w", "pi", "wattTheta", "pi.net", "tajD", "tajD.denom", "pi.wPop2", "pi.wPop1", "wattTheta.Pop2", "wattTheta.Pop1", "tajD.denomPop2", "tajD.denomPop1", "ShannonsIndex.Between", "ShannonsIndex.Net", "ShannonsIndex.Pop1", "ShannonsIndex.Pop2"  )


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

    numTaxonPairs <- sum(simDat[1,paste("PRI.Psi.",
                                        1:simDat[1,"PRI.numTauClass"],
                                        sep="")])
    # This means there are 3 taxon pairs, and constrained to have 3 tau's.
    # No need to analyze Psi.1 etc because they are always set to 1.
    # If we don't remove these, the analysis chokes
    if(numTaxonPairs == simDat[1,"PRI.numTauClass"]) {
      noPsiAnalysis <- T
      result$prior.names <-
        result$prior.names[- grep("^PRI[.]Psi[.]", result$prior.names)]
    } else {
      noPsiAnalysis <- F
    }
  } else {
    constrained <- F
  }
  
  # set min, max boundary of prior parameters, and verbose print message
  # PRI.omega >= 0, PRI.E.t >= 0, PRI.Psi > 0
  min.val <- list(PRI.Psi = 1, PRI.var.t = 0, PRI.E.t = 0, PRI.omega = 0 )
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
    
    temp.val <- sapply(rep("(= number of taxon pairs that divergence at corresponding tau)", length(psiNames)), list)
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


if(is.null(result$PRI.E.t$vals)==FALSE){
	big_table<-c()
	for (i in 1:length(result$prior.names)){
	    thisPriorName <- result$prior.names[i]
	    big_table<-rbind(big_table,result[[thisPriorName]]$x,result[[thisPriorName]]$vals)
	    }
	rownames(big_table)<-rownames(big_table,do.NULL=FALSE) 
	    
	for (i in 1:length(result$prior.names)){
	    	thisPriorName <- result$prior.names[i]
		name.post<- sub("PRI[.]", "Pos.LLR.", thisPriorName)
		name.post.raw<-sub("PRI[.]", "Pos.wo.LLR.", thisPriorName)
		rownames(big_table)[(i*2-1)]<-name.post
		rownames(big_table)[i*2]<-name.post.raw
		}
		
		write.table(cbind(t(big_table),result$PRI.Psi$ss),file="posterior_table",row.names = FALSE)
}else{
	cat("Local linear regression is not used. Please use '-a' to get the accepted posteriors.\n\n")
}

  real.mode.mean.median <- NULL
  modeToPrint <- NULL
  fileToPrint <- paste("acceptedPriorSummary_")
  for(i in 1:length(used.stats))
    {
      fileToPrint <- paste(fileToPrint, used.stats[i], sep = "")
    }

  tempMatrix <- scan(file = obs.infile, what = double(0), skip = 1, nmax = 5)
  TempMatrix <- matrix(tempMatrix, ncol = 1, byrow = T)
  truePRI <- c(TempMatrix[2,1], TempMatrix[3,1], TempMatrix[4,1], TempMatrix[5,1])
  counter <- 0
 
  cat("######### results #######\n")
  # make print loop for Mode and quantile
  for (i in 1:length(result$prior.names)) {
    counter <- counter + 1
    thisPriorName <- result$prior.names[i]
    name.rm.PRI <- sub("PRI[.]", "", thisPriorName)
    if(! is.null(verbose.print[[thisPriorName]])) {
      additional.print <- verbose.print[[thisPriorName]]
    } else {
      additional.print <- ""
    }

    cat ("##### Summary of posterior distribution #####\n")
    cat ("#####",name.rm.PRI, additional.print, "#####\n")

    # When constarined with large numTauClasses, PRI.Psi.1 become
    # mostly 1, and locfit has following problem, see help(locfit.raw)
    # Warning: procv: density estimate, empty integration region
    # Error: newsplit: out of vertex space Error: Descend tree proble
    # So, simply printing the posterior mean, and mode from
    # accepted values
    if (length(grep("^PRI\\.Psi\\.[0-9]+$", thisPriorName)) == 1) {
      if (! rejmethod) 
        cat ("### With local-linear regression (by default)\n")
      else
        cat ("### With simple rejection method, NO LOCAL LINEAR REGRESSION\n")
      mean.median.vect <- c(mean((result[[thisPriorName]])$x),
                                 median((result[[thisPriorName]])$x))
      names(mean.median.vect) <- c("mean", "median")
      
      cat ("## Mean/Median\n")
      print(mean.median.vect)
      
      cat ("## 95 % quantile\n")
      print(quantile((result[[thisPriorName]])$x,prob=c(0.025,0.975)))
      if (! rejmethod)
        cat("\n### With simple rejection method, CAUTION: not using local-linear regression\n");
      cat ("## posterior distribution\n")
      if (! rejmethod) 
        post.distn.accRej <- table(result[[thisPriorName]]$vals)
      else
        post.distn.accRej <- table(result[[thisPriorName]]$x)
      temp.pd.ar <- c("frequencies", post.distn.accRej)
      names(temp.pd.ar)[1] <- name.rm.PRI
      print(temp.pd.ar)
      cat ("## Mode (from simple rejection method):\n")
      print(names(which(post.distn.accRej == max(post.distn.accRej))))
      
      real.mode.mean.median <- append(real.mode.mean.median, c(truePRI[counter],1,mean.median.vect), after = length(real.mode.mean.median))
      
    } else {  # Do regular summary, Not PRI.Psi.*
      cat ("### Mode\n")
      # With locfit 1.5_4, newsplit error of locfit() will stop the
      # analysis.  So I'm doing the error handling by myself with try().
      res.mode <- try(loc1stats((result[[thisPriorName]])$x,prob=0.95),silent=T)
      if(class(res.mode) == "try-error") {
        cat("NA\n")
        cat("** locfit failed to find mode, this occurs with an " ,
            "extreme\n** L-shaped posterior distribution, and the mode is ",
            "at the boundary (e.g. 0)\n", sep="")

        if(name.rm.PRI == "Psi")
          { modeToPrint <- 1.0 }
        else if((name.rm.PRI == "var.t")||(name.rm.PRI == "omega"))
          { modeToPrint <- 0.0 }
      } else {
        cat ("MODE:\n" )
        res.mode <- res.mode[1]
        print(res.mode)
        modeToPrint <- res.mode
      }
      cat ("### Mean/Median\n")
      mean.median.vect <- c(mean((result[[thisPriorName]])$x),
                                 median((result[[thisPriorName]])$x))
      names(mean.median.vect) <- c("mean", "median")
      print(mean.median.vect)
      
      cat ("### 95 % quantile\n")
      print(quantile((result[[thisPriorName]])$x,prob=c(0.025,0.975)))

      real.mode.mean.median <- append(real.mode.mean.median, c(truePRI[counter], modeToPrint, mean.median.vect),after = length(real.mode.mean.median))
    }
    cat("\n")
  }

  write(real.mode.mean.median, file = fileToPrint, ncol = 20, append = T)
  
  # Print out figures
  pdf(pdf.outfile, width=7.5, height=10, paper="letter")
layout(mat=matrix(1:2,2,1))
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
#    par(mfcol=c(2,1))  # 2 plots per page
    if (pre.rejected) {
      make.hist(prior.dat[,thisPriorName],result[[thisPriorName]], title=this.title, breaks=20)
      
    } else {
      make.hist(simDat[,thisPriorName],result[[thisPriorName]],title=this.title,breaks=20)
      plot.bf(simDat[,thisPriorName],result[[thisPriorName]]$x,main="Bayes Support for true Hyper-parameter value < threshold")
    }
#    par(mfcol=old.mfcol)
  }

#  print(thisPriorName)
  plot((result[["PRI.omega"]])$x,(result[["PRI.E.t"]])$x,lty=2,lwd=0.5,ylim=c(0,max(prior.dat[["PRI.E.t"]])),xlim=c(0,max(prior.dat[["PRI.omega"]])))

  rc <- try(plotKernDensity(result[["PRI.omega"]],result[["PRI.E.t"]],
                            xlab="Omega", ylab="E(t)", title="Omega and E(t)"))
  

if(class(rc) == "try-error") {
    cat("WARN: plotKernDensity failed for some reason, so the kernel density ",
        "plot was not created\n", file=stderr())
  }
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
  bwq1 <- try(dpik(res1$x),silent=T)
  bwq2 <- try(dpik(res2$x),silent=T)
  # I think bandwidth choice by dpik may fail if there aren't enough unique
  # values (e.g. mostly 0), I'm not sure following is ok or not, but
  # bw.nrd0 seems to be more robust
  if(class(bwq1) == "try-error") {
    cat("INFO: In plotKernDensity(), simpler bandwidth selection used for ",
            xlab, "\n", file=stderr())
    bwq1 <- bw.nrd0(res1$x)
  }
  if(class(bwq2) == "try-error") {
    cat("INFO: In plotKernDensity(), simpler bandwidth selection used for ",
            ylab, "\n", file=stderr())
    bwq2 <- bw.nrd0(res1$x)
  }  
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
#  par(mfcol=c(3,1))
bw_vect<-max(vect)/100
bw_res.makepd<-max(res.makepd$x)/50
  hist.col = "white"
  if(missing(xlim)) {
    hist(vect,col=hist.col,border ="white",xlim=c(0,max(vect)),ylim=c(0,max(density(vect,bw=bw_vect)$y,density(res.makepd$x,bw=bw_res.makepd)$y)),prob=TRUE,main=paste(title),xlab=title, ...)
  } else {
    hist(vect,col=hist.col,freq=F,xlim=xlim,prob=TRUE,main=paste(title),
       xlab=title, ...)
  }
  lines(density(vect,bw=bw_vect),lty=2,col="red",pch=3)
#  if(missing(xlim)) {
#    hist(res.makepd$x,col="blue",prob=TRUE,freq=F,xlim=c(0,max(vect)),
#         main=paste(title, ": Posterior Dist'n "), xlab=title, ...)
#  } else {
#    hist(res.makepd$x,col=hist.col,xlim=xlim,freq=F,prob=TRUE,
#         main=paste(title, ": Posterior Dist'n "), xlab=title, ...)
#  
  lines(density(res.makepd$x,bw=bw_res.makepd),lty=1,col="blue",pch=3)
  #par(mfcol=old.mfcol)
legend("topright",c("Prior","Posterior"),text.col=c("red","blue"),lty=c(2,1),col=c("red","blue"))

}

plot.bf <- function(prior,posterior,...) {
  bf.out <- make.bf.vect(prior,posterior)
  if (is.null(bf.out)) {
    return(NULL);
  }
  # finding some reasonable y.lim for the plot, all values here are arbitrary
  y.max <- 40
  if (max(bf.out$bf) * 1.1 < 40) {
    y.max <- max(bf.out$bf)* 1.1
  } else if (mean(bf.out$bf > 40) > 0.6) {
    y.max <- sort(bf.out$bf)[ceiling(length(bf.out$bf) * 0.4)] * 1.1
  }
  plot(bf.out$crit.vals, bf.out$bf, type="l",
       ylim=c(0,y.max), 
       xlab="Hyper Parameter Thresholds",ylab="Bayes Factor",...)
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

    post.range <- range(posterior)

    # A single value of posterior is observed (e.g. at the boundary)
    if(post.range[1] == post.range[2]) {
      if(post.range[1] == 0) {
        post.range <- c(-0.01, 0.01)
      } else  {
        post.range <- c(0.99, 1.01) * post.range[1]
      }
    }

    # following is probably not required
    pri.min <- min(prior)
    if (post.range[1] <= pri.min) {
      if(pri.min >= post.range[2]) {  # impossible, but better check
        cat("WARN: In make.bf.vect(), weird prior/posterior range encountered\n",
            file=stderr())
        return(NULL)
      }
      post.range[1] <- pri.min
    }
    crit.vals <- seq(post.range[1], post.range[2], length.out=num.points+2)
    # get rid of the ends to avoid bf = Inf
    crit.vals <- crit.vals[2:(length(crit.vals) - 1)]
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
