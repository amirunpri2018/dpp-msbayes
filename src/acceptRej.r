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
params.from.priorDistn <- c("Psi", "var.t", "E.t", "omega")

# If there are 9 taxon pairs, there are 9 columns of "pi.b" stats for
# each taxon pairs, then 9 columns of "pi.w", ... etc.
# If you change the names here, don't forget to update the default value
# of used.stats in stdAnalysis
summary.stat.names <- c("pi.b", "pi.w", "pi", "wattTheta", "pi.net",
                            "tajD", "tajD.denom")

# Just print out the stat names
printStatNames <- function() {
  cat("## params.from.priorDistn\n")
  cat(params.from.priorDistn)
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
stdAnalysis <- function(obs.infile, sim.infile, pdf.outfile="figs.pdf",
                        tol=0.002,
                        used.stats=c("pi","wattTheta","pi.net","tajD.denom"),
                        rejmethod=T,
                        return.res=F
                        ) {
  simDat <- getData(sim.infile)
  
  nPairs <- simDat[["numTaxonPairs"]]
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
  
  #acceptance/regression, .... ie  the meat
  print("here1")
  result.omega <- makepdANY(as.vector(obsDat[1,usedColNames],mode="numeric"),
                        simDat[,"omega"], simDat[,usedColNames], tol,
                        rep(T,len=nrow(simDat)),rejmethod=rejmethod)
  print("here2")
  result.Psi <- makepdANY(as.vector(obsDat[1,usedColNames],mode="numeric"),
                        simDat[,"Psi"], simDat[,usedColNames], tol,
                        rep(T,len=nrow(simDat)),rejmethod=rejmethod)
  print("here3")
  result.E.t <- makepdANY(as.vector(obsDat[1,usedColNames],mode="numeric"),
                        simDat[,"E.t"], simDat[,usedColNames], tol,
                        rep(T,len=nrow(simDat)),rejmethod=rejmethod)
  
  # Transformations, absorbing boundary at 0.0
  result.omega$x[which(result.omega$x < 0)] <- 0
  result.Psi$x[which(result.Psi$x < 0)] <- 0
  result.E.t$x[which(result.E.t$x < 0)] <- 0
  result.E.t$x[which(result.E.t$x > nPairs)] <- nPairs
  
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

  # Print out figures
  pdf(pdf.outfile)
  make.hist(simDat[,"omega"], result.omega, title="Omega (Var(t)/E(t))")
  make.hist(simDat[,"Psi"], result.Psi, title="Psi (# of possible divtimes)")
  make.hist(simDat[,"E.t"], result.E.t, title="E(t)")
  
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
getData <- function (infile) {
  dat <- read.table(infile, header=F)

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
