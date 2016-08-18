# v3.04 - Add adjust of p-values for two-tails tests
# v3.03 - Add TIE estimations (same event type) for pValue and statistic  
# v3.02 - Add estimation of all combinations
# v3.01 - Add valorate.plot.kaplan, valorate.risk, added kullback-leibler in valorate.plot.diff.empirical and correction of maximum
# v3.00 - Implements weights (should be included in the SUM OF V, see Peto, Tarone-Ware, Flemington weights in Kleinbaum : Survival Analysis - A self learning text)
# v2.31 - default for min.sampling.size raised to 1000
# v2.3 - Bug corrections for TIES estimations (also speed improvements)
# v2.2 - alpha function missing from scales package
# v2.1 - IMPLEMENTATION OF TIES
# v2.1 - Implementation of shades in valorate.plot.empirical and addition of valorate.plot.diff.empirical
# v2.0 - C implementation
# v2.0 - options to avoid (beta, weibull, etc) parameter estimations and estimated only if needed (to save processing time)
# v2.0 - vjx and vjcx matrices now only 1 matrix (only the division term)
# v1.2 - Removed "&& k.dens[k+1]" (not so sure since it will be futile)
# v1.1 - Added min.sampling.size
# v1.? - sampling.size is now used to express "total sampling" instead of sampling in any number of events.

# PENDINGS:

#Reverted Changed:
# v2.31 - pvalues are now 0-1 range, p-value in sampling is now estimated as V > sampling if V >= 0 or V < sampling when V < 0

require(methods)
require(survival)

setClass("valorate", representation(
	s="numeric", 
	n="numeric",
	events="numeric",
	parameters="list", 
	sampling.size="numeric",
	min.sampling.size="numeric",
	verbose="logical",
	save.sampling="logical",
	wcensored="numeric",
	wevents="numeric",
	order="numeric",
	subpop="environment",
	time="character",
	samplings="list",
	ties="list",
	sampling.ties="numeric",
	tiesame="list",
	tiesame.pos="numeric",
	tiesame.sampling="numeric",
	method="character",
	estimate.distribution.parameters="character",
	weights="numeric",
	weights.method="character",
	tails="numeric"
	))

valorate.libfile <- paste("valorate_samplings",.Platform$dynlib.ext,sep="")
if (file.exists(valorate.libfile)) {
	if (!is.loaded("valorate_sampling",PACKAGE="valorate_samplings")) {
		oklib <- NULL
		try(oklib <- dyn.load(valorate.libfile))
		if (is.null(oklib)) {
			cat("There are problems loadding C methods library [",valorate.libfile,"].\n")
		}			
	}
}


# Create valorate object 
new.valorate <- function(time, censored, rank, sampling.size=max(10000,200000/events), 
	min.sampling.size=1000, verbose=FALSE, save.sampling=TRUE, method="C", cPath="", 
	estimate.distribution.parameters=c("empirical","gaussian","beta","weibull")[1], 
	sampling.ties=30, weights=NULL, 
	weights.method=c("logrank","Wilcoxon","Tarone-Ware","Peto","Flemington-Harrington","Trevino","user")[1], 
	weights.parameters=list(p=1,q=1,t=3),
	tails=2) {
# Input: 
#	(1) time=Survival Times and censored=censored indicator (1=censored, 0=dead)
#	(2) time=Survival Times with censoring indicator "+" (censored=missing)
#	(3) rank=0/1 ranked subjects, 0 is censored, 1 is event

	# Create valorate object
	
	libfile <- paste(cPath, "valorate_samplings",.Platform$dynlib.ext,sep="")
	loaded <- is.loaded("valorate_sampling",PACKAGE="valorate_samplings")
	if (method == "C" && !loaded) {
		if (!file.exists(libfile)) {
			cat("C methods library have not been loaded and does not seem to be present on [",libfile,"]. Setting to R methods.\n")
			method <- "R"
		} else {
			oklib <- NULL
			try(oklib <- dyn.load(libfile))
			if (is.null(oklib)) {
				cat("There are problems loadding C methods library [",libfile,"], setting to R methods.\n")
				method <- "R"
			}			
		}
	}
	
	vro <- new("valorate")
	vro@time <- character(0)
	ties <- list()
	tiesame <- list()

	if (!missing(time) && !missing(censored) && missing(rank)) {
		if (is.numeric(time) && (is.numeric(censored) || is.logic(censored)) && length(time) == length(censored)) {
			o <- order(time)
			s <- 1-(1*censored[o])
			vro@time <- paste(time[o],ifelse(censored[o],"+",""),sep="")
			vro@order <- o
		} else {
			stop("time needs to be numeric, censored needs to be numeric (0/1) or logical.")
		}
	} else if (!missing(time) && missing(censored) && missing(rank)) {
		if (is.character(time)) {
			newtime <- gsub("[^0-9Ee\\.\\+\\-]", "", time)
			wc <- grep("\\+$",newtime)
			s <- rep(1,length(time))
			s[wc] <- 0
			thetimes <- as.numeric(gsub("\\+","",newtime))
			if (any(!is.finite(thetimes))) {
				stop(paste("problems with time parameter:",paste(newtime,collapse=", ")))
			}
			o <- order(thetimes)
			s <- s[o]
			vro@order <- o
			vro@time <- paste(thetimes[o],ifelse(s,"","+"),sep="")
		} else {
			stop("time needs to be character when not censored parameter is provided.")
		}
	} else if (missing(time) && missing(censored) && !missing(rank)) {
		if (is.numeric(rank) || is.logical(rank)) {
			s <- 1*(rank==1)
			vro@order <- 1:length(s)
			vro@time <- paste(1:length(s), ifelse(s,"","+"),sep="")
		} else {
			stop("rank needs to be numeric (0/1) or logical.")
		}
	} else {
		stop("Invalid input. Possible inputs: (1) time:numeric, censorded:numeric/logical, (2) time:character using '+' as censoring at the end, (3) rank:numeric/logical.")
	}

	weights.method <- match.arg(weights.method, c("logrank","Wilcoxon","Tarone-Ware","Peto","Flemington-Harrington","Trevino","user"))	
	if (!is.null(weights)) {
		if (!is.numeric(weights) || !all(is.finite(weights)) || length(weights) != length(s)) {
			stop("Weights needs to be valid numeric values and length equal to time or rank.")
		}
		#weights <- as.double(weights * length(time) / sum(weights))
		weights <- as.double(weights)
		weights.method <- "user"
	} else {
		# Survival Analysis - A Self Learning Text - Kleinbaum - Klein Springer 2005 pg 64
		if (weights.method == "Wilcoxon") {
			weights <- as.double(length(vro@time):1)
		} else if (weights.method == "Tarone-Ware") {
			weights <- as.double(sqrt(length(vro@time):1))
		} else if (weights.method == "Peto") {
			sf <- as.double(length(vro@time):1 / length(vro@time))
			weights <- sf
		} else if (weights.method == "Flemington-Harrington") {
			sf <- as.double(length(vro@time):1 / length(vro@time))
			weights <- as.double(sf^weights.parameters$p + (1-sf)^weights.parameters$q)
		} else if (weights.method == "Trevino") {
			# This method is mine (Victor Trevino).
			# It gives more weight to high density times
			if (missing(time)) {
				stop("Trevino method need time values instead of rank values.")
			}
			xh <- hist(time, breaks=seq(min(time),max(time),length=length(time)/weights.parameters$t), plot=FALSE)
			weights <- numeric(length(time))
			for (i in 1:length(time)) {
				wr <- which(xh$mids > time[i])
				if (length(wr) == 0) wr <- length(xh$mids)
				wr <- wr[1]
				wl <- rev(which(xh$mids < time[i]))
				if (length(wl) == 0) wl <- 1
				wl <- wl[1]
				if (wl == wr) {
					if (wl == 1) { wr <- 2; } else { wl = wr-1; }
				}
				m <- (xh$counts[wr]-xh$counts[wl]) / (xh$mids[wr]-xh$mids[wl])
				b <- xh$counts[wl]
				x <- (time[i] - xh$mids[wl])
				weights[i] <- m*x + b
			}
			weights <- as.double(weights[vro@order]*length(time)/sum(weights))
		} else {
			weights <- as.double(rep(1,len=length(vro@time)))
		}
	}

	# set parameters
	if (!missing(time)) vro@parameters$time <- time
	if (!missing(censored)) vro@parameters$censored <- censored
	if (!missing(rank)) vro@parameters$rank <- rank

	events <- sum(s)
	vro@s <- s 								# subjects ranked by time and event indicated by 1 (0 for censoring)
	vro@events <- as.integer(events) 		# the total number of events
	vro@n <- as.integer(length(s))			# the total number of subjects
	vro@verbose <- verbose					# verbose while running ?
	vro@save.sampling <- save.sampling 		# save sampling in object ? 
	vro@wcensored <- as.integer(which(s == 0))
	vro@wevents <- as.integer(which(s == 1))
	vro@subpop <- new.env()
	vro@sampling.size <- sampling.size		# the number of samples needed for sampling the distribution
	vro@min.sampling.size <- min.sampling.size
	vro@method <- ifelse(method=="C","C","R")
	vro@estimate.distribution.parameters <- estimate.distribution.parameters
	vro@tails <- ifelse(tails==1, 1, 2)
	if (length(vro@time)) {
		## Check for ties
		newtime <- gsub("[^0-9Ee\\.\\+\\-]", "", vro@time)
		thetimes <- as.numeric(gsub("\\+","",newtime))
		ttimes <- table(thetimes)
		if (any(ttimes > 1)) {
			rept <- as.numeric(names(ttimes[ttimes > 1]))
			for (i in rept) {
				w <- which(thetimes == i)
				if (length(unique(s[w])) > 1) {
					# There are events and censoring, create tie record.
					ties[[length(ties)+1]] <- w
				} else {
					tiesame[[length(tiesame)+1]] <- w
				}
			}
		}
		valorate.cat(vro, "There were",length(ties),"ties having events and censoring at same time point.\n")
		valorate.cat(vro, "There were",length(tiesame),"time points having ties of events or censoring.\n")
	}
	vro@ties <- ties
	vro@tiesame <- tiesame
	vro@tiesame.pos <- unlist(tiesame)
	vro@tiesame.sampling <- max(1,sampling.ties*(length(tiesame) > 0))
	vro@sampling.ties <- 0
	if (length(ties)) {
		lties <- unlist(lapply(ties, length))
		pties <- unlist(lapply(lties, function(x) valorate.perm(x,x)))
		# the permutations are:
		# elements in ties: 2, 3,  4,   5,   6
		# permutations    : 2, 6, 24, 120, 720
		ncombs <- max(cumprod(pties))
		ncomb <- min(sampling.ties, ncombs*4) ## maximum number of matrices		
		vro@sampling.ties <- ncomb
	}
	vro@weights <- weights
	vro@weights.method <- weights.method

	return ((vro))
}

# utility function
valorate.p.overlap <- function(balls, white, drawn, observed.white, lower.tail=FALSE) {
	black <- balls - white
	phyper(observed.white, white, black, drawn, lower.tail=lower.tail)
}

# utility function
valorate.hyper.density <- function(npop, nevents, nmut, kmut) {
	valorate.p.overlap(npop,nevents,nmut,kmut,lower.tail=TRUE)-valorate.p.overlap(npop,nevents,nmut,kmut-1,lower.tail=TRUE)
}

# utility function
valorate.comb <- function(n, x) { 
	return ( exp(lfactorial(n) - (lfactorial(x)+lfactorial(n-x))) )
}
# utility function
valorate.perm <- function(n,x) {
	return (exp(lfactorial(n)-lfactorial(n-x)))
}
# utility function
valorate.cat <- function(vro, ...) {
	if (vro@verbose) {
		cat(...)
		flush.console()
	}
}

#taken from http://druedin.com/2012/08/11/moving-averages-in-r/
# moving average utility function
valorate.mav <- function(x,n=5,avoid.na=FALSE){
	y <- stats::filter(x,rep(1/n,n), sides=2)
	if (avoid.na) {
		#VT:adapted to avoid NA
		wna <- which(is.na(y))
		if (length(wna)) {
			y[wna] <- x[wna]
		}	
	}
	y
}


# Prepare object to estimate n1 population
setGeneric("prepare.n1", function(vro, n1) standardGeneric("prepare.n1"))
setMethod("prepare.n1", signature("valorate", "numeric"), 
  function(vro, n1) {
	nxname <- paste("subpop", n1, sep="")
	if (nxname %in% ls(vro@subpop)) {
		# already exists, it doesn´t need to be recomputed
		return (invisible(vro))
	}

	##### Prepare object to add this n1 population

	#load parameters from object
	n <- vro@n
	s <- vro@s
	events <- vro@events
	nx <- n1 <- as.integer(n1)
	ne <- min(n1, events) # n events to be used
	sampling.size <- vro@sampling.size
	min.sampling.size <- vro@min.sampling.size
	k.dens <- valorate.hyper.density(n,events,nx,0:ne)
	valorate.cat(vro, "[[[ Estimating Log-Rank Distribution for a Risk Group of", n1, "]]]\n")
	valorate.cat(vro, "Densities:\n", k.dens,"\n")
	valorate.cat(vro, "Densities Log:\n", round(log10(k.dens),2),"\n")

	build.vjx <- function(s, n, nx) {
		#build V matrix,
		#rows (i) is the number of samples that still are in n1
		#cols (j) is the time ("n" different times = samples)
		vjx <- matrix(NA, ncol=n, nrow=nx) #array(NA, dim=c(nx,n,2))
		for (j in  1:n) { 
		    if (s[j]) {
		        den <- n - j + 1
		        vjx[,j] <- 1:nx/den
		        #for (i in 1:nx) {
		        #    num = i
		        #    div = num/den
		        #    vjx[i,j] = div;
		        #    ####if (i <= n-j && nx-i < j) {
		        #    #    vjx[i,j,1] = -div; # xj=0, sample is not included in n1
		        #    ####}
		        #    ####if (i <= n-j && nx-i < j) {
		        #    #    vjx[i,j,2] = 1-div; # xj=1, sample is included in n1
		        #    ####}
		        #}
		    }
		}
		vcjx <- vjx[,s[1:n]==1,drop=FALSE]		#vjx[,s[1:n]==1,,drop=FALSE]
		return (list(vjx=vjx, vcjx=vcjx))
	}
	vjm <- build.vjx(s, n, nx)
	#vjx <- vjm$vjx # Not Needed
	vcjx <- vjm$vcjx
	vcjx.n <- as.integer(1)

	ties <- vro@ties
	if (length(ties) > 0 && vro@sampling.ties > 0) {
		lties <- unlist(lapply(ties, length))
		pties <- unlist(lapply(lties, function(x) valorate.perm(x,x)))
		# the permutations are:
		# elements in ties: 2, 3,  4,   5,   6
		# permutations    : 2, 6, 24, 120, 720
		ncombs <- max(cumprod(pties))
		ncomb <- vro@sampling.ties #min(vro@sampling.ties, ncombs*4) ## maximum number of matrices
		valorate.cat(vro, "Estimating",ncomb,ifelse(is.finite(ncombs),paste("of",ncombs,sep=" "),""),"combinations for",length(ties),"ties in",ne,"events: ")
		for (i in 1:ncomb) {
			if (i %% 10 == 0) {
				valorate.cat(vro, ".")
				if (i %% 100 == 0) valorate.cat(vro, " ")
			}
			s <- vro@s
			for (j in 1:length(ties)) {
				s[ties[[j]]] <- sample(s[ties[[j]]])
				#vjx <- cbind(vjx, vjm$vjx) # Not needed
			}
			vjm <- build.vjx(s, n, nx)
			vcjx <- cbind(vcjx, vjm$vcjx)
		}
		valorate.cat(vro, "Done.\n")
		s <- vro@s
		vcjx.n <- as.integer(vcjx.n + ncomb)
	}


	#Estimate densities per each combination of events:censored for n1 subjects
	dens <- list()
	v <- numeric(min.sampling.size)
	wcensored <- vro@wcensored #which(s[1:n] == 0)
	wevents <- vro@wevents #which(s[1:n] == 1)
	ncensored <- as.integer(length(wcensored))
	nevents <- as.integer(length(wevents)) # should be equal to events
	combinations <- numeric(ne+1)
	verbose <- as.integer(vro@verbose)
	weightev <- as.double(vro@weights[wevents])
	inn1 <- integer(n)
	ldx  <- integer(n)
	valorate.cat(vro, "Simulating taking",nx,"samples of",n,"having",events,"events:\n")
	first <- as.integer(1)
	for (k in 0:ne) {
		allComb <- NULL
		allCombMatrix <- 0
		ncomb <- valorate.comb(n-events,nx-k) * valorate.comb(events,k)
		combinations[k+1] <- ncomb
		sim <- as.integer(round(max(min.sampling.size, min(ncomb*4, sampling.size*k.dens[k+1])))) #as.integer(min(ncomb*2, sampling.size)) #round(max(min.sampling.size, min(ncomb*2, sampling.size*k.dens[k+1])))
		if (sim*2 > ncomb && ncomb > 0) {
			mcens <- combn(wcensored, nx-k)
			mevnt <- combn(wevents, k)
			if (k > 0 && nx-k > 0) {
				allComb <- matrix(as.integer(0), nrow=nx, ncol=ncol(mcens)*ncol(mevnt))
				mk <- 0
				for (mi in 1:ncol(mcens)) {
					for (mj in 1:ncol(mevnt)) {
						mk <- mk + 1
						allComb[,mk] <- c(mcens[,mi],mevnt[,mj])
					}
				}
			} else if (k > 0) {
				allComb <- mevnt
			} else {
				allComb <- mcens
			}
			storage.mode(allComb) <- "integer"
			sim <- as.integer(round(ncomb,0))
		}
		if (sim != length(v)) {
			v <- numeric(sim)
		}
		valorate.cat(vro, "Sampling",k,"of",ne,"events, comb=",ncomb,"=",format(k.dens[k+1]*100,digits=3),"%, Size=")
		if (ncomb > 0) { ## && k.dens[k+1]
			valorate.cat(vro, sim,":")
			doAllCombinations <- (!is.null(allComb))
			if (doAllCombinations) {
				valorate.cat(vro,"[*All Combinations*]",dim(allComb)) #,class(allComb),mode(allComb),is.integer(allComb)
				allCombMatrix <- allComb
			}
			if (vro@method == "R") {
				for (i in 1:sim) {
					if (i %% 1000 == 0) {
						valorate.cat(vro, ".")
						if (i %% 10000 == 0)
							valorate.cat(vro, " ")
					}
					#Generate inn1
					inn1[] <- 0
					if (doAllCombinations) {
						inn1[allComb[,i]] <- 1
					} else {
						if (k < nx) {
							### censored
							inn1[sample(wcensored,nx-k)] <- 1
						}
						if (k > 0) {
							### events
							inn1[sample(wevents,k)] <- 1
						}
					}
					# Calculate the V statistic
					#vcjx <- vro@subpop[[nxname]]$vcjx
					einn1 <- inn1[wevents]#1+inn1[wevents]
					ldx <- nx - cumsum(c(0,inn1))[wevents]
					os <- if (vcjx.n == 1) 0 else round(runif(1)*(vcjx.n-1))*events # OffSet of the matrix of ties, each of n columns
					V <- 0
					for (j in 1:events) {
						if (ldx[j] == 0) break;
						V <- V + weightev[j] * (einn1[j]-vcjx[ldx[j],j+os]) #vcjx[ldx[j],j,einn1[j]]
					}
					v[i] <- V
					#v[i] <- sum(weightev*(einn1-vcjx[ldx,1:events+os]))
				}
				dens[[k+1]] <- v
			} else {
				v <- .C("valorate_sampling", v=v, sim, n, k, nx, wcensored, ncensored, wevents, nevents, weightev, vcjx, vcjx.n, inn1, ldx, first, verbose, allCombMatrix, PACKAGE="valorate_samplings")$v
				dens[[k+1]] <- v					
			}
			first <- as.integer(0)
		} else {
			valorate.cat(vro, 0,":")
			dens[[k+1]] <- 0
		}
		valorate.cat(vro, "\n")
	}

	# sumarize the samples in an empirical distribution considering the weights per number of events
	xmin <- unlist(lapply(dens,min,na.rm=TRUE))
	xmax <- unlist(lapply(dens,max,na.rm=TRUE))
	nmx <- max(xmax)
	nmn <- min(xmin)
	empirical <- 0
	empirical.breaks <- seq(nmn*ifelse(nmn<0,1.001,.999),nmx*ifelse(nmx>0,1.001,.999),len=1001)
	emp.hist <- list()
	for (k in 0:ne) {
		dens.hist <- hist(pmax(pmin(nmx,dens[[k+1]]),nmn),plot=FALSE,breaks=empirical.breaks)
		smv <- valorate.mav(dens.hist$count/sum(dens.hist$count), 10)
		smv[is.na(smv)] <- 0
		empirical <- empirical + k.dens[k+1] * smv
		#empirical <- empirical + k.dens[k+1] * dens.hist$count/sum(dens.hist$count)
		emp.hist[[k+1]] <- hist(dens[[k+1]],plot=FALSE,breaks=1001)
	}

	vro@subpop[[nxname]] <- list(
		nx=nx,
		ne=ne,
		n1=n1,
		#vjx=vjx, # Not Needed
		vcjx=vcjx, 
		vcjx.n=vcjx.n,
		m=unlist(lapply(dens,mean,na.rm=TRUE)),
		sd=unlist(lapply(dens,sd,na.rm=TRUE)),
		min=xmin,
		max=xmax,
		k.density=k.dens,
		combinations=combinations,
		empirical=empirical,
		empirical.breaks=empirical.breaks,
		emp.hist=emp.hist)
	if (vro@save.sampling) {
		vro@subpop[[nxname]]$sampling <- dens
	}


	# Gaussian fitting
	# This are always estimated since are anyway stored
	vro@subpop[[nxname]]$gaussian <- list(mean=vro@subpop[[nxname]]$m, sd=vro@subpop[[nxname]]$sd)
	if ("gaussian" %in% vro@estimate.distribution.parameters) {
		valorate.cat(vro, "Estimated Gaussian parameters for",nx,"samples:\n")
		valorate.cat(vro, "mean=",vro@subpop[[nxname]]$gaussian$mean,"\n")
		valorate.cat(vro, "sd=",vro@subpop[[nxname]]$gaussian$sd,"\n")
	}

	# Beta fitting
	if ("beta" %in% vro@estimate.distribution.parameters) {
		valorate.estimate.beta.parameters(vro, vro@subpop[[nxname]])
	}

	# Weibull fitting
	if ("weibull" %in% vro@estimate.distribution.parameters) {
		valorate.estimate.weibull.parameters(vro, vro@subpop[[nxname]])
	}

	return (invisible(vro))
})



valorate.p.value.normal <- function(vro, vrsubo, lrv, z) {
	if (length(z) > 1) {
		unlist(sapply(unlist(z),function(x) valorate.p.value.normal(vro, vrsubo, 0, x)))
	} else {
		p <- 1-pnorm(z)
		min(p,1-p) * vro@tails
		#p <- if (z >= 0) 1-pnorm(z) else pnorm(z)
		#p
		#both are the same, the first code is clearer
	}
}

valorate.p.value.chisq <- function(vro, vrsubo, lrv, z) {
	if (length(z) > 1) {
		unlist(sapply(unlist(z),function(x) valorate.p.value.chisq(vro, vrsubo, 0, x)))
	} else {
		#p <- 1-pchisq(z^2,df=1)
		#min(p,1-p)
		p <- 1-pchisq(z^2,df=1)
		p # here there is no adjust
	}
}

valorate.p.value.gaussian <- function(vro, vrsubo, lrv, z) {
	if (length(lrv) > 1) {
		unlist(sapply(unlist(lrv),function(x) valorate.p.value.gaussian(vro, vrsubo, x)))
	} else {
		ncomb <- valorate.comb(vro@n, vrsubo$nx)
		p <- sum((1-pnorm(lrv,mean=vrsubo$gaussian$mean,sd=vrsubo$gaussian$sd)) * vrsubo$k.density * ncomb, na.rm=TRUE) / ncomb
		min(p,1-p) * vro@tails
	}
}

valorate.p.value.weibull <- function(vro, vrsubo, lrv, z) {
	if (length(lrv) > 1) {
		unlist(sapply(unlist(lrv),function(x) valorate.p.value.weibull(vro, vrsubo, x)))
	} else {
		valorate.estimate.weibull.parameters(vro, vrsubo)
		ncomb <- valorate.comb(vro@n, vrsubo$nx)
		x <- (lrv-vrsubo$min)/(vrsubo$max-vrsubo$min)
		p <- sum((1-pweibull(x,vrsubo$weibull$k,vrsubo$weibull$l)) * vrsubo$k.density * ncomb, na.rm=TRUE) / ncomb
		min(p,1-p) * vro@tails
	}
}

valorate.p.value.beta <- function(vro, vrsubo, lrv, z) {
	if (length(lrv) > 1) {
		unlist(sapply(unlist(lrv),function(x) valorate.p.value.beta(vro, vrsubo, x)))
	} else {
		valorate.estimate.beta.parameters(vro, vrsubo)
		ncomb <- valorate.comb(vro@n, vrsubo$nx)
		x <- (lrv-vrsubo$min)/(vrsubo$max-vrsubo$min)
		p <- sum((1-pbeta(x,vrsubo$beta$alpha,vrsubo$beta$beta)) * vrsubo$k.density * ncomb, na.rm=TRUE) / ncomb
		min(p,1-p) * vro@tails
	}
}

valorate.estimate.beta.parameters <- function(vro, vrsubo) {
	if (is.null(vrsubo$beta)) {
		nx <- vrsubo$nx
		combinations <- vrsubo$combinations
		nxname <- paste("subpop", nx, sep="")
		valorate.cat(vro, "Estimating Beta parameters for",nx,"samples: ")
		estBetaParams <- function(mu, var) {
		  alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
		  beta <- alpha * (1 / mu - 1)
		  return(params = list(alpha = alpha, beta = beta))
		}		
		dens <- vrsubo$sampling
		e.b.a <- c()
		e.b.b <- c()
		for (i in 1:length(dens)) {
			## Scale dens[[i]] to 0-1 
			valorate.cat(vro, i, " ")
			e.b.a[i] <- NA
			e.b.b[i] <- NA			
			if (combinations[i] > 2) {
				x <- (dens[[i]]-min(dens[[i]]))/(max(dens[[i]])-min(dens[[i]]))
				xb <- estBetaParams(mean(x),var(x))
				e.b.a[i] <- xb$alpha
				e.b.b[i] <- xb$beta			
			}
		}
		valorate.cat(vro, "done\n")
		vro@subpop[[nxname]]$beta <- list(alpha=e.b.a, beta=e.b.b)
		valorate.cat(vro, "alpha=",vro@subpop[[nxname]]$beta$alpha,"\n")
		valorate.cat(vro, "beta=",vro@subpop[[nxname]]$beta$beta,"\n")	
	}
}


valorate.estimate.weibull.parameters <- function(vro, vrsubo) {
	if (is.null(vrsubo$weibull)) {
		nx <- vrsubo$nx
		nxname <- paste("subpop", nx, sep="")
		dens <- vrsubo$sampling
		combinations <- vrsubo$combinations
		valorate.cat(vro, "Estimating Weibull parameters for",nx,"samples: ")
		e.w.k <- c()
		e.w.l <- c()
		for (i in 1:length(dens)) {
			valorate.cat(vro, i, " ")
			e.w.k[i] <- NA
			e.w.l[i] <- NA
			if (combinations[i] > 2) {
				## Scale dens[[i]] to 0-1 
				x <- sort( (dens[[i]]-min(dens[[i]]))/(max(dens[[i]])-min(dens[[i]])) )
				#x <- x + x[2]/2 ##to avoid NA
				#first estimate
				#weibull plot: https://en.wikipedia.org/wiki/Weibull_distribution
				ly <- log(-log(1-cumsum(x)/(sum(x)))) #+x[1] to avoid NA
				lx <- log(x)
				wok <- which(!is.finite(lx) | !is.finite(ly))
				if (length(wok)) {
					lx <- lx[-wok]
					ly <- ly[-wok]
					x <- x[-wok]
				}
				if (length(x)) {
					xlm <- lm(ly~lx)
					k <- xlm$coefficients[2] ## aproximacion inicial
					if (is.na(k)) {
						## This is an empirical observation, when almost all x are the same
						## Almost never happen
						l <- median(x)
						k <- l*100+length(dens[[i]])*0.05
					} else {
						l <- exp(-xlm$coefficients[1]/k)						
						kweibullest <- function(x,k) 1 / (sum(x^k*log(x))/sum(x^k)-mean(log(x)))
						dk <- 1
						# aproximaciones siguientes
						j <- 1
						while(abs(dk) > 0.001 && j < 100) {
							knew = kweibullest(x,k)
							dk <- knew-k
							k = knew
							j <- j + 1
						}
						l <- mean(x^k)^(1/k)
					}
					e.w.k[i] <- k
					e.w.l[i] <- l
				}
			}
		}
		valorate.cat(vro, "done\n")
		vro@subpop[[nxname]]$weibull <- list(k=e.w.k, l=e.w.l)
		valorate.cat(vro, "k=",vro@subpop[[nxname]]$weibull$k,"\n")
		valorate.cat(vro, "l=",vro@subpop[[nxname]]$weibull$l,"\n")
	}
}

valorate.p.value.sampling <- function(vro, vrsubo, lrv, z) {
	if (length(lrv) > 1) {
		unlist(sapply(unlist(lrv),function(x) valorate.p.value.sampling(vro, vrsubo, x)))
	} else {
		ncomb <- valorate.comb(vro@n, vrsubo$nx)
		if (is.null(vrsubo$sampling)) {
			warning("Sampling is better for estimation but it was not saved. Estimated by histogram.")
			nplus <- unlist(lapply(vrsubo$emp.hist, function(x) {
				a <- sum(c(0,x$counts[x$mids > lrv]), na.rm=TRUE)
				b <- sum(c(0,x$counts[rev(which(x$mids <= lrv))[1]]/2), na.rm=TRUE) 
				(a + b*(a > 0))/ sum(x$counts)
				} ))
			nminus <- unlist(lapply(vrsubo$emp.hist, function(x) {
				a <- sum(c(0,x$counts[x$mids < lrv]), na.rm=TRUE)
				b <- sum(c(0,x$counts[which(x$mids >= lrv)[1]]/2), na.rm=TRUE) 
				(a + b*(a > 0))/ sum(x$counts)
				} ))

		} else {
			nminus <- pmin(1,unlist(lapply(vrsubo$sampling, function(x) (1+sum(x < lrv, na.rm=TRUE)) / length(x))))
			nplus <- pmin(1,unlist(lapply(vrsubo$sampling, function(x) (1+sum(x > lrv, na.rm=TRUE)) / length(x))))
		}
		pleft <- sum(nminus * vrsubo$k.density, na.rm=TRUE)
		pright <- sum(nplus * vrsubo$k.density, na.rm=TRUE)
		p <- min(pleft, pright)
		min(p, 1-p) * vro@tails
	}
}

valorate.p.value.all <- function(vro, vrsubo, lrv, z=NULL) {
	p <- list(
		normal=valorate.p.value.normal(vro, vrsubo, lrv, z),
		gaussian=valorate.p.value.gaussian(vro, vrsubo, lrv, z),
		weibull=valorate.p.value.weibull(vro, vrsubo, lrv, z),
		beta=valorate.p.value.beta(vro, vrsubo, lrv, z)
		)
	if (!is.null(z)) {
		p$chisq=valorate.p.value.chisq(vro, vrsubo, lrv, z)
	}
	if (!is.null(vrsubo$sampling)) {
		p$sampling = valorate.p.value.sampling(vro, vrsubo, lrv, z)
	}
	return (p)
}

#setGeneric("survdiff", function(vro, clusters, p.func) standardGeneric("survdiff"))
#setMethod("survdiff", signature("valorate", "numeric", "function"),
valorate.survdiff <- function(vro, clusters, p.func=valorate.p.value.sampling) {
	uc <- unique(clusters)
	if (length(uc) != 2) {
		warning("Clusters needs to be only 2.")
		p <- 1
		attr(p, c("Statistic(LR,Z,X2)")) <- c(LR=0,Z=0,X2=0)
		return (p)
	}
	if (length(clusters) != vro@n) {
		stop(paste("Cluster length=",length(cluster),"is different to n=",vro@n))
	}
	if (FALSE && any(table(clusters) == 1)) {
		warning("A cluster has only 1 sample.")
		p <- 1
		attr(p, c("Statistic(LR,Z,X2)")) <- c(LR=0,Z=0,X2=0)
		return (p)
	}
	clusters <- clusters[vro@order]
	c1 <- sum(clusters==uc[1])
	c2 <- length(clusters)-c1
	nx <- min(c1, c2)
	prepare.n1(vro, nx)
	nxname <- paste("subpop", nx, sep="")

	if (c1 < c2) {
		inn1 <- (clusters == uc[1])*1
	} else {
		inn1 <- (clusters == uc[2])*1
	}

	# Calculate the V statistic
	vcjx <- vro@subpop[[nxname]]$vcjx
	events <- vro@events
	wevents <- vro@wevents

	nsamp <- vro@tiesame.sampling
	inties <- which(inn1 %in% vro@tiesame.pos)
	if (length(inties) && length(unique(inn1[inties])) > 1) {
		# Look at tie positions to see which needs to be resampled
		tiepos <- list()
		for (tp in 1:length(vro@tiesame)) {
			w <- vro@tiesame[tp]
			if (length(unique(inn1[w])) > 1) {
				# these positions do need
				tiepos[[length(tiepos)+1]] <- w
			}
		}
	} else {
		nsamp <- 1
	}
	VR <- numeric(nsamp)
	VRZ <- numeric(nsamp)

	for (ivr in 1:nsamp) {
		if (ivr > 1) {
			# resample only active tie positions
			for (tp in 1:length(tiepos)) {
				w <- tiepos[tp]
				inn1[w] <- sample(inn1[w])
			}
		}
		V <- 0
		VZ <- 0
		einn1 <- inn1[wevents] #1+inn1[wevents]
		weightev <- vro@weights[wevents]
		ldx <- nx - cumsum(c(0,inn1))[wevents]
		for (j in 1:events) {
			if (ldx[j] == 0) break;
			V <- V + weightev[j]*(einn1[j] - vcjx[ldx[j],j]) #vcjx[ldx[j],j,einn1[j]]
			VZ<- VZ+             (einn1[j] - vcjx[ldx[j],j]) #vcjx[ldx[j],j,einn1[j]]
		}
		if (vro@sampling.ties > 0) {
			VV <- numeric(vro@sampling.ties+1)
			for (i in 1:vro@sampling.ties) {
				offset <- i*events
				v <- 0
				for (j in 1:events) {
					if (ldx[j] == 0) break;
					v <- v + weightev[j]*(einn1[j] - vcjx[ldx[j],j+offset]) #vcjx[ldx[j],j,einn1[j]]
				}
				VV[i] <- v
			}
			VV[vro@sampling.ties+1] <- V
			#V <- V / (vro@sampling.ties + 1)
			V <- mean(VV)		
		}
		VR[ivr] <- V
		VRZ[ivr] <- VZ
	}
	V <- mean(VR)
	VZ <- mean(VRZ)

	# Calculate variance if needed in normal test
	n   <- vro@n
	n1j <- (sum(inn1)-(c(0,cumsum(inn1))[1:n]))[wevents]
	nj  <- (n:1)[wevents]
	oj  <- (vro@s)[wevents]
	vj  <- oj*(n1j/nj)*(1-n1j/nj)*(nj-oj) / pmax(1,nj-1)

	names(V) <- NULL
	zv <- VZ/sqrt(sum(vj))

	# Estimate p-value
	p <- p.func(vro, vro@subpop[[nxname]], V, zv)
	attr(p, c("Statistic(LR,Z,X2,LRunw)")) <- c(LR=V,Z=zv,X2=zv^2,LR.unweigthed=VZ,pZ=1-pnorm(zv),pX2=1-pchisq(zv^2,df=1))
	if (nsamp == 1 && vro@sampling.ties > 0) {
		xp <- p.func(vro, vro@subpop[[nxname]], VV, zv)
		attr(p, c("diff.ties_CI(VRmin,VRmax,pMin,pMax,n)")) <- c(VRmin=min(VV),VRmax=max(VV),pMin=min(xp),pMax=max(xp),n=length(VV))
	}
	if (nsamp > 1) {
		xp <- p.func(vro, vro@subpop[[nxname]], VR, zv)
		attr(p, c("ties_CI(VRmin,VRmax,pMin,pMax,n)")) <- c(VRmin=min(VR),VRmax=max(VR),pMin=min(xp),pMax=max(xp),n=length(VR))
	}

	return (p)
}

valorate.psurvdiff <- function(vro, n1, v, z=NULL, p.func=valorate.p.value.sampling) {
	nx <- n1
	prepare.n1(vro, nx)
	nxname <- paste("subpop", nx, sep="")
	# Estimate p-value
	p <- p.func(vro, vro@subpop[[nxname]], v, z)
	attr(p, "statistic") <- v
	p
}

valorate.plot.empirical <- function(vro, n1, vstat=NULL, type="l", log="", add=FALSE, include=c("none","gaussian","beta","weibull","all")[1], xlab="valorate LR", ylab="density", main=paste("Empirical Density: Mutations=",n1,ifelse(is.null(vstat),"",paste0("\n(marked statistic at ",vstat,")"))), samp=vro@sampling.size, smooth=10, legends=FALSE, shades=c(6,8), transparency=0.25, lwd=2, xlim=range(c(mids,vstat))+(c(-0.05,+0.05)*abs(range(c(mids,vstat)))), minL=NA, minR=NA, ...) {
	if (length(n1) > 1) {
		if (all(n1 %in% c(0,1))) {
			##### n1 represents risk groups, estimate statistic and real n1
			vstat <- valorate.survdiff(vro, n1)
			vstat <- attributes(vstat)[[1]][1]
			n1 <- min(length(n1)-sum(n1), sum(n1))
		} else {
			n1 <- n1[1]
		}
	}
	nx <- n1
	prepare.n1(vro, nx)
	nxname <- paste("subpop", nx, sep="")
	sp <- vro@subpop[[nxname]]
	f <- if (add) points else plot
	mids <- sp$empirical.breaks[-1]-diff(sp$empirical.breaks)/2
	empv <- valorate.mav(sp$empirical,smooth)
	#if (any(empv == 0 | is.na(empv))) empv[empv == 0 | is.na(empv)] <- min(empv[empv > 0], na.rm=TRUE)
	if (is.finite(minL) & any((empv == 0 | is.na(empv)) & mids < 0)) empv[(empv == 0 | is.na(empv)) & mids < 0] <- minL
	if (is.finite(minR) & any((empv == 0 | is.na(empv)) & mids > 0)) empv[(empv == 0 | is.na(empv)) & mids > 0] <- minR
	f(mids, empv, type=type, log=log, xlab=xlab, ylab=ylab, main=main, lwd=lwd, xlim=xlim, ...)
	res <- NULL
	if (!is.null(vstat)) {
		if (length(vstat) > 1) {
			res <- numeric(length(vstat))
			for (i in 1:length(vstat)) {
				res[i] <- empv[which.min(abs(mids-vstat[i]))[1]]
			}
			rug(res)
		} else {
			avstat <- attributes(vstat)
			if (length(avstat)==1 && names(avstat)=="Statistic(LR,Z,X2)" && names(avstat[[1]])[1] == "LR") {
				vstat <- avstat[[1]][1]
			}
			v <- min(max(par("usr")[1],vstat), par("usr")[2])
			#rug(v, col=shades[1])
			miny <- min(empv[empv > 0],na.rm=TRUE)
			maxy <- max(empv,na.rm=TRUE)
			abline(v=v, lty=2, col=shades[1])
			left.x <- c(mids[mids < vstat],vstat,vstat,mids[1])
			left.y <- c(empv[mids < vstat],empv[mids >= vstat][1],miny,miny)
			left.y[!is.finite(left.y) | left.y == 0] <- miny
			right.x <- c(vstat, mids[mids > vstat], mids[length(mids)])
			right.y <- c(miny, empv[mids > vstat], miny)
			right.y[!is.finite(right.y) | right.y == 0] <- miny
			p <- valorate.p.value.sampling(vro, sp, vstat, 0)
			left <- sum(empv[mids < vstat], na.rm=TRUE)
			right <- sum(empv[mids > vstat], na.rm=TRUE)
			if (right > left) { right = 1-p; left = p; }
			else { right = p; left = 1-p; }
			left.col <- 2-1*(left < right)
			right.col <- 3-left.col
			# Taken from scales package 
			alpha <- function (colour, alpha = NA) 
				{
				    col <- grDevices::col2rgb(colour, TRUE)/255
				    if (length(colour) != length(alpha)) {
				        if (length(colour) > 1 && length(alpha) > 1) {
				            stop("Only one of colour and alpha can be vectorised")
				        }
				        if (length(colour) > 1) {
				            alpha <- rep(alpha, length.out = length(colour))
				        }
				        else if (length(alpha) > 1) {
				            col <- col[, rep(1, length(alpha)), drop = FALSE]
				        }
				    }
				    alpha[is.na(alpha)] <- col[4, ][is.na(alpha)]
				    new_col <- grDevices::rgb(col[1, ], col[2, ], col[3, ], alpha)
				    new_col[is.na(colour)] <- NA
				    new_col
				}
			polygon(left.x, left.y, col=alpha(shades[left.col],transparency))
			polygon(right.x, right.y, col=alpha(shades[right.col],transparency))
			l <- format(left,digits=4)
			r <- format(right,digits=4)
			mc <- max(nchar(c(l,r)))
			if (left < 1 && left > 0 && nchar(l) < mc && length(grep("e",l)) == 0 && length(grep("\\.",l)) > 0) l <- format(1-right, digits=nchar(r)-2) #l <- paste(l,paste(rep("0",mc-nchar(l)),collapse=""),sep="")
			if (right < 1 && right > 0 && nchar(r) < mc && length(grep("e",r)) == 0 && length(grep("\\.",r)) > 0) r <- format(1-left, digits=nchar(l)-2) #paste(r,paste(rep("0",mc-nchar(r)),collapse=""),sep="")
			srt <- ifelse(((vstat-strwidth(l)) < (par("usr")[1])) || ((vstat+strwidth(r)) > (par("usr")[2])),90,0)

			if (l=="1" && right > 1e-3) l <- format(1-right)
			if (r=="1" && left  > 1e-3) r <- format(1-left)

			if (l=="1") { l <- expression(phantom() %=~% 1) } 
			else { if (l=="0") { l <- expression(phantom() %=~% 0) } }

			if (r=="1") { r <- expression(phantom() %=~% 1) } 
			else { if (r=="0") { r <- expression(phantom() %=~% 0) } }

			text(vstat,maxy,l,adj=if(srt == 0) c(1.10,1.3) else c(1,-0.25) ,col=shades[left.col], cex=1, srt=srt)
			text(vstat,maxy,r,adj=if(srt == 0) c(-0.10,1.3) else c(1,1.25) ,col=shades[right.col], cex=1, srt=srt)
		}
	}
	if (is.null(vstat)) {
		res <- data.matrix(data.frame(x=mids, y=empv))
	}
	if (length(include) > 0 && any(include %in% c("gaussian","beta","weibull","all"))) {
		ns <- samp
		xc <- 2
		brks <- sp$empirical.breaks #c(-Inf,sp$empirical.breaks,Inf)
		mp <- NULL
		minmax <- function(x, min=0, max=1, quantiles=FALSE) {
			if (quantiles) {
				min = quantile(x, min, na.rm=TRUE)
				max = quantile(x, max, na.rm=TRUE)
			}
			x[x < min] <- min
			x[x > max] <- max
			x
		}
		if (any(include == "all")) include <- c("gaussian","beta","weibull")
		for (i in include) {
			rs <- list()
			ec <- 0
			for (j in 1:length(sp$k.density)) {
				if (i == "gaussian") {
					r <- rnorm(ns, m=sp$gaussian$mean[j], sd=sp$gaussian$sd[j])
				} else if (i == "beta") {
					valorate.estimate.beta.parameters(vro, sp)
					sp <- vro@subpop[[nxname]]	
					r <- rbeta(ns, sp$beta$alpha[j], sp$beta$beta[j])
				} else if (i == "weibull") {
					valorate.estimate.weibull.parameters(vro, sp)
					sp <- vro@subpop[[nxname]]	
					r <- rweibull(ns, sp$weibull$k[j], sp$weibull$l[j])
				} else r <- NULL
				if (!is.null(r)) {
					r <- if(i == "gaussian") r else (r*(sp$max[j]-sp$min[j]) + sp$min[j])
					h <- hist(pmax(pmin(r,max(sp$max)),min(sp$min)), breaks=brks, plot=FALSE)
					if (is.null(mp)) mp <- h$mids
					ec <- ec + h$counts * sp$k.density[j]
				}
			}
			ex <- ec*sum(empv,na.rm=TRUE)/sum(ec) # *max(empv, na.rm=TRUE)/max(exs, na.rm=TRUE)
			# minmax(ex,0.02,0.98,TRUE)
			exs <- valorate.mav(ex,smooth)
			points(mp, exs, type="l", col=xc, lwd=lwd[min(length(lwd),2)]) # ec*par("usr")[4]/max(ec) ... 
			xc <- xc + 1
		}
		if (legends) legend("topleft", legend=c("empirical",include), lty=1, col=0:length(include)+1)
	}
	invisible(res)
}



plot.kaplan.valorate <- function(data, time, status, logrank=NULL, plogrank=NULL, main= "", cluster, risk.groups=length(unique(cluster)), draw.main=FALSE, short.names=FALSE, mark.cex=1, mark=3, margins=TRUE, col=1:risk.groups+1, col.main=16, pName="p") {
	library(survival)
	ocox <- NULL
	try(ocox <- coxph(Surv(time, status) ~ ., data.frame(t(data))))
	if (is.null(ocox)) {
		shr <- list(conf.int=matrix(0,ncol=4,nrow=1),coefficients=matrix(0,ncol=10,nrow=1))
	} else {
		shr <- summary(ocox)  	
	}
	hr <- round(shr$conf.int[1,1],2)
	hrst <- paste(format(round(shr$conf.int[1,1],1),digits=4)," [",format(round(shr$conf.int[1,3],1),digits=4),"~",format(round(shr$conf.int[1,4],1),digits=4),"], p=",format(shr$coefficients[1,5], digits=4),sep="")
	censxrisk <- rep(0,risk.groups) ##cuanta cuantos censored hay en cada grupo 
	toview <- 1:(risk.groups+ifelse(draw.main,1,0)) ## las poblaciones que se van a mostrar en la leyenda 
	for(i in 1:risk.groups){
		w <- which(cluster == i)
		censxrisk[i] <- sum(status[w]) ##hace el conteo de censored por grupo 
	}

	if (margins) {
	  pp = par(mar=c(5.1+max(0,risk.groups-2), 2.5, 5.1, 1.6))
	  on.exit(par(pp))
	}
	lrLabel <- ifelse(short.names,", LR=",", Log-Rank=")
	hrLabel <- ifelse(short.names,"\nHR=","\nHazard Ratio= ")
	if (is.null(plogrank)) {
		xsd <- survdiff(Surv(time, status) ~ cluster)
		plogrank <- 1-pchisq(xsd$chisq,df=1)
	}
	plot(survfit(Surv(time, status) ~ cluster, 
		data=data.frame(time=time, status=status, cluster=cluster)), 
		col=col, 
		lty=1, lwd=2, 
		conf.int=FALSE,
		main=paste(main,lrLabel,round(logrank,1),", ",pName,"=",format(plogrank, digits=4),hrLabel,hrst,sep=""),
		ylim=c(0, 1.05),
		cex.main=1.25, cex=mark.cex, mark=mark)
	pos <- c("bottomleft", "bottomright", "top", "right")
	xt <- axTicks(1)
	xt <- sort(c(xt, xt[-length(xt)]+diff(xt)/2))
	if (draw.main) {  ##Plots the whole population line (gray)
	lines(survfit(Surv(time, status) ~ 1), mark.time=TRUE, col=col.main,
	  conf.int=FALSE, lty=1, lwd=1, mark=mark, cex=mark.cex, cex.main=1.25)
	lines(survfit(Surv(time, status) ~ cluster, 
	    data=data.frame(time=time, status=status, cluster=cluster)), 
	    col=col, 
	    lty=1, lwd=2, 
	    conf.int=FALSE,
	    cex.main=1.25, cex=mark.cex, mark=mark)
	}
	uk <- sort(unique(cluster))
	for (k in 1:length(uk)) { ##Genera y plotea el numero de paciente por cada grupo con relacion al tiempo. 
	w <- which(cluster == uk[k])
	axis(1, xt, labels=sapply(xt, function(x) sum(time[w] >= x)), tick=FALSE, line=k, col.axis=col[which(uk[k] == uk[order(uk)])])
	}
	if (draw.main) {
	axis(1, xt, labels=sapply(xt, function(x) sum(time >= x)), tick=FALSE, line=k+1, col.axis=col.main)  	
	}
}

valorate.plot.kaplan <- function(vro, data, p=valorate.survdiff(vro, data), main="", short.names=TRUE, draw.all=FALSE, 
	mark="|", mark.cex=0.75, margins=TRUE, col=2:3, col.all="skyblue", ...) {

	time <- as.numeric(sub("+","",vro@time,fixed=TRUE))
	status <- rep(1,length(time))
	status[grep("\\+",vro@time)] <- 0

	plot.kaplan.valorate(data=matrix(data[vro@order],nrow=1), plogrank=p, logrank=attributes(p)[[1]][1], 
		time=time, status=status, 
		main=main, short.names=short.names, draw.main=draw.all,
		cluster=1+as.vector(data[vro@order]), 
		mark=mark, mark.cex=mark.cex, margins=margins, col=col, col.main=col.all, pName="pValorate")

}


valorate.risk <- function(vro, data, ...) {

	time <- as.numeric(sub("+","",vro@time,fixed=TRUE))
	status <- rep(1,length(time))
	status[grep("\\+",vro@time)] <- 0

	data <- as.vector(data[vro@order])
	ocox <- NULL
	try(ocox <- coxph(Surv(time, status) ~ ., data.frame(data=data)))
	if (is.null(ocox)) {
		shr <- list(conf.int=matrix(0,ncol=4,nrow=1),coefficients=matrix(0,ncol=10,nrow=1))
	} else {
		shr <- summary(ocox)
	}
	hr <- shr$conf.int[1,1]
	attr(hr, "confidence.interval") <- shr$conf.int[1,3:4]
	attr(hr, "p.value") <- shr$coefficients[1,5]
	attr(hr, "hazard.ratio") <- shr
	hr
}



valorate.plot.diff.empirical <- function(vro, n1, type="l", log="", include=c("gaussian","beta","weibull","all")[4], xlab="valorate LR", ylab="density", main=paste("Differences of Densities: Mutations=",n1), samp=vro@sampling.size, smooth=10, ylim=c(min(miny,-maxy),max(-miny,maxy)), legends=TRUE, ...) {
	nx <- n1
	prepare.n1(vro, nx)
	nxname <- paste("subpop", nx, sep="")
	sp <- vro@subpop[[nxname]]
	ns <- samp
	#xc <- 2
	brks <- sp$empirical.breaks #c(-Inf,sp$empirical.breaks,Inf)
	mp <- NULL
	minmax <- function(x, min=0, max=1, quantiles=FALSE) {
		if (quantiles) {
			min = quantile(x, min, na.rm=TRUE)
			max = quantile(x, max, na.rm=TRUE)
		}
		x[x < min] <- min
		x[x > max] <- max
		x
	}
	if (any(include == "all")) include <- c("gaussian","beta","weibull")
	s <- c()
	xp <- list()
	miny <- 0
	maxy <- 0
	for (i in include) {
		rs <- list()
		ec <- 0
		for (j in 1:length(sp$k.density)) {
			if (i == "gaussian") {
				r <- rnorm(ns, m=sp$gaussian$mean[j], sd=sp$gaussian$sd[j])
			} else if (i == "beta") {
				valorate.estimate.beta.parameters(vro, sp)
				sp <- vro@subpop[[nxname]]	
				r <- rbeta(ns, sp$beta$alpha[j], sp$beta$beta[j])
			} else if (i == "weibull") {
				valorate.estimate.weibull.parameters(vro, sp)
				sp <- vro@subpop[[nxname]]	
				r <- rweibull(ns, sp$weibull$k[j], sp$weibull$l[j])
			} else r <- NULL
			if (!is.null(r)) {
				r <- if(i == "gaussian") r else (r*(sp$max[j]-sp$min[j]) + sp$min[j])
				h <- hist(pmax(pmin(r,max(sp$max)),min(sp$min)), breaks=brks, plot=FALSE)
				if (is.null(mp)) mp <- h$mids
				ec <- ec + h$counts * sp$k.density[j]
			}
		}
		ex <- ec/sum(ec)
		d <- sp$empirical - ex #(ex*max(sp$empirical)/max(minmax(ex,0.02,0.98,TRUE)))
		#points(mp, valorate.mav(d,smooth), type="l", col=xc) # ec*par("usr")[4]/max(ec) ... 
		ds <- valorate.mav(d,smooth)
		xp[[length(xp)+1]] <- ds
		s[i] <- sum(abs(d))
		KL.div <- function(p, q) {
			kld <- p*log(p/q)
			sum(kld[p > 0 & q > 0])
		}
		s[i] <- KL.div(sp$empirical, ex) #(ex*max(sp$empirical)/max(minmax(ex,0.02,0.98,TRUE))))
		miny <- min(miny, ds[is.finite(ds)])
		maxy <- max(maxy, ds[is.finite(ds)])
	}
	plot(0,0,xlim=range(sp$empirical.breaks), ylim=ylim, xlab=xlab, ylab=ylab, main=main, type="n", ...)
	abline(h=0, col=1)
	for (i in 1:length(xp)) {
		points(mp, xp[[i]], type="l", col=i+1)
	}
	if (legends) legend("topleft", legend=c("empirical",paste(include,round(s,5))), lty=1, col=0:length(include)+1, ncol=2)
}

valorate.plot.subpop.empirical <- function(vro, which=NULL, type="l", log="", xlim=NULL, smooth=10, legends=3, density=TRUE, ylim=c(ymin, ymax), ...) {
	xsp <- sort(as.numeric(gsub("subpop","",ls(vro@subpop))))
	if (!is.null(which)) {
		xsp <- xsp[xsp %in% which]
	}
	xmin <- Inf
	xmax <- -Inf
	ymax <- 0
	ymin <- Inf
	dmin <- Inf
	dmax <- -Inf
	for (nx in xsp) {
		nxname <- paste("subpop", nx, sep="")
		sp <- vro@subpop[[nxname]]
		dmin <- min(dmin, median(diff(sp$empirical.breaks)))
		dmax <- max(dmax, median(diff(sp$empirical.breaks)))
	}
	smoothing <- if (smooth > 0) function(x) valorate.mav(x, smooth) else function(x) x
	for (nx in xsp) {
		nxname <- paste("subpop", nx, sep="")
		sp <- vro@subpop[[nxname]]
		d <- median(diff(sp$empirical.breaks))
		y <- sp$empirical*(dmax/d)/sum(sp$empirical) #*(1+(dmax-d)/(dmax-dmin))
		xmin <- min(xmin, sp$min)
		xmax <- max(xmax, sp$max)
		ymax <- max(ymax, y)
		ymin <- min(ymin, y[y > 0])
	}
	if (is.null(xlim)) {
		xlim <- c(xmin, xmax)
	}
	f <- plot
	icol <- 1
	if (!density) {
		ymax <- 1
		femp <- function(y, factor) y / max(y)
	} else {
		femp <- function(y, factor) y * factor / sum(y)
	}
	for (nx in xsp) {
		nxname <- paste("subpop", nx, sep="")
		sp <- vro@subpop[[nxname]]
		d <- median(diff(sp$empirical.breaks))
		y <- femp(sp$empirical, dmax/d)
		f(sp$empirical.breaks[-1]-diff(sp$empirical.breaks)/2, smoothing(y), 
			type=type, log=log, xlim=xlim, ylim=ylim, col=icol, xlab="Log-Rank Statistic", ylab=paste(ifelse(density,"Normalized Densities","Scaled Densities")," (smooth=",smooth,")",sep=""),...)
		f <- points
		icol <- icol + 1
	}
	if (legends > 0) {
		legend("topleft", legend=xsp, lty=1, col=1:length(xsp), ncol=legends, cex=0.75)
	}
}


valorate.plot.subpop.empirical.to.0 <- function(vro, which=NULL, type="l", log="", xlim=NULL, smooth=10, legends=3, density=TRUE, ylim=c(ymin, ymax), ...) {
	xsp <- sort(as.numeric(gsub("subpop","",ls(vro@subpop))))
	if (!is.null(which)) {
		xsp <- xsp[xsp %in% which]
	}
	xmin <- Inf
	xmax <- -Inf
	ymax <- 0
	ymin <- Inf
	dmin <- Inf
	dmax <- -Inf
	for (nx in xsp) {
		nxname <- paste("subpop", nx, sep="")
		sp <- vro@subpop[[nxname]]
		dmin <- min(dmin, median(diff(sp$empirical.breaks)))
		dmax <- max(dmax, median(diff(sp$empirical.breaks)))
	}
	smoothing <- if (smooth > 0) function(x) valorate.mav(x, smooth) else function(x) x
	w0 <- c()
	wL <- c()
	for (nx in xsp) {
		nxname <- paste("subpop", nx, sep="")
		sp <- vro@subpop[[nxname]]
		d <- median(diff(sp$empirical.breaks))
		y <- sp$empirical*(dmax/d)/sum(sp$empirical) #*(1+(dmax-d)/(dmax-dmin))
		xmin <- min(xmin, sp$min)
		xmax <- max(xmax, sp$max)
		ymax <- max(ymax, y)
		ymin <- min(ymin, y[y > 0])
		w0[nx] <- which(sp$empirical.breaks >= 0)[1]
		wL[nx] <- length(sp$empirical.breaks)
	}
	if (is.null(xlim)) {
		xlim <- c(xmin, xmax)
	}
	f <- plot
	icol <- 1
	if (!density) {
		ymax <- 1
		femp <- function(y, factor) y / max(y)
	} else {
		femp <- function(y, factor) y * factor / sum(y)
	}
	neg <- max(w0-1, na.rm=TRUE)
	pos <- max(wL-w0, na.rm=TRUE)
	wLm <- median(w0, na.rm=TRUE)
	#cat("w0=",w0,"\n")
	#cat("neg=",neg,",pos=",pos,",wLm=",wLm,"\n")
	plot(0,0,type="n",xlim=c(wLm-neg,wLm*2+(pos-wLm)),ylim=ylim,xlab="'Scaled' Log-Rank Statistic", ylab=paste(ifelse(density,"Normalized Densities","Scaled Densities")," (smooth=",smooth,")",sep=""), xaxt="n", ...)
	for (nx in xsp) {
		nxname <- paste("subpop", nx, sep="")
		sp <- vro@subpop[[nxname]]
		d <- median(diff(sp$empirical.breaks))
		y <- femp(sp$empirical, dmax/d)
		sy <- smoothing(y)
		points(1:length(sy)+wLm-w0[nx], sy, 
			type=type, col=icol, ...)
		icol <- icol + 1
	}
	abline(v=wLm,col=1,lwd=0.5,lty=3)
	wlmed <- median(wL, na.rm=TRUE)
	ats <- wLm + seq(-1,1,len=9)*wlmed
	axis(side=1, at=ats, labels=round((ats-wLm)*2/wlmed,2))
	#axis(side=1)
	if (legends > 0) {
		legend("topleft", legend=xsp, lty=1, col=1:length(xsp), ncol=legends, cex=0.75)
	}
}

valorate.plot.sampling.densities <- function(vro, n1, type="l", log="", xlim=c(minx, maxx), ylim=c(miny, maxy), ncol=max(1,round(sqrt(n1)/2)), main="Sampling Densities Per Event and Weighted Sum", empirical=mean(mxy,na.rm=TRUE), rug=TRUE, add=FALSE, w.sum=TRUE, sampling=TRUE, ...) {
	nx <- n1
	prepare.n1(vro, nx)
	nxname <- paste("subpop", nx, sep="")
	sp <- vro@subpop[[nxname]]
	dens <- lapply(sp$sampling, density)
	#minx <- min(unlist(lapply(dens, function(x) min(x$x))))
	#maxx <- max(unlist(lapply(dens, function(x) max(x$x))))
	minx <- min(sp$min)
	maxx <- max(sp$max)	
	miny <- min(min(unlist(lapply(dens, function(x) min(x$y)))), min(sp$empirical))+1e-13
	my <- unlist(lapply(dens, function(x) max(x$y)))*1+1e-13
	maxy <- max(my)
	minxdens <- min(unlist(lapply(dens, function(x) min(x$x))))
	maxxdens <- max(unlist(lapply(dens, function(x) max(x$x))))
	dx <- maxxdens-minxdens
	L <- length(dens[[1]]$x)+1 #/2
	xaprx <- seq(minxdens,maxxdens,length=L)
	yaprx <- numeric(length(xaprx))
	mxy <- numeric(length(dens))
	if (!add) plot(1,1, type="n", log=log, xlim=c(minx, maxx), ylim=ylim, main=main, xlab="", ylab="", ...)
	for (i in (1:length(dens))[order(sp$k.density)]) {
		points(dens[[i]], type=type, col=(i+1))
		idx <- round(L * (dens[[i]]$x-minxdens) / dx + 1)
		w <- which(diff(idx) > 0)
		yaprx[idx[w]] <- yaprx[idx[w]] + dens[[i]]$y[w] * sp$k.density[i]
		if (rug) rug(sp$sampling[[i]], col=i+1, side=(i %% 2)*2+1)
		mxy[i] <- max(dens[[i]]$y, na.rm=TRUE)
	}
	is.logy <- (log=="y" || log=="xy")
	leg <- c(0:sp$ne,"W. Sum","Sampling")
	lok <- 1:length(leg)
	if (!sampling) lok <- lok[-length(lok)]
	if (!w.sum) lok <- lok[-grep("W. Sum",leg)]
	if (!add) legend(if(is.logy) "bottomleft" else "topleft",legend=leg[lok],lty=1,col=c(2+0:sp$ne,sp$ne+3,1)[lok], ncol=ncol, text.width=strwidth("0"),box.col=0)
	if (!add) text(sp$m, pmin(my,maxy*.75), paste("x ",format(sp$k.density,digits=2)), srt=90, col=2+0:sp$ne, adj=c(ifelse(is.logy,1,0),1))
	if (sampling) points(sp$empirical.breaks[-1]-diff(sp$empirical.breaks)/2, sp$empirical*mean(mxy)/max(sp$empirical, na.rm=TRUE), type=type, lwd=2, col=1, ...)
	if (w.sum) points(xaprx, yaprx*mean(mxy)/max(yaprx,na.rm=TRUE), type=type, col=sp$ne+3, lwd=3)
	axis(4,at=axTicks(2),round(max(sp$empirical)*empirical*axTicks(2)/max(axTicks(2)),5))
	invisible(list(x=xaprx, y=yaprx))
}


valorate.plot.sampling.densities.figure <- function(vro, n1, type="l", log="", xlim=c(minx, maxx), ylim=c(miny, maxy), main="", empirical=mean(mxy,na.rm=TRUE), rug=1, rug.size=1, w.sum=TRUE, sampling=TRUE, ncol=1, ...) {
	nx <- n1
	prepare.n1(vro, nx)
	nxname <- paste("subpop", nx, sep="")
	sp <- vro@subpop[[nxname]]
	dens <- lapply(sp$sampling, density)
	lens <- lapply(sp$sampling, length)
	#minx <- min(unlist(lapply(dens, function(x) min(x$x))))
	#maxx <- max(unlist(lapply(dens, function(x) max(x$x))))
	minx <- min(sp$min)
	maxx <- max(sp$max)	
	miny <- min(min(unlist(lapply(dens, function(x) min(x$y)))), min(sp$empirical))+1e-13
	my <- unlist(lapply(dens, function(x) max(x$y)))*1+1e-13
	maxy <- max(my)
	minxdens <- min(unlist(lapply(dens, function(x) min(x$x))))
	maxxdens <- max(unlist(lapply(dens, function(x) max(x$x))))
	dx <- maxxdens-minxdens
	L <- length(dens[[1]]$x)+1 #/2
	xaprx <- seq(minxdens,maxxdens,length=L)
	yaprx <- numeric(length(xaprx))
	mxy <- numeric(length(dens))
	par(mfrow=c(ceiling((length(dens)+sampling*1+w.sum*1)/ncol),ncol))
	mar <- par("mar")
	par(mar=c(0,mar[2],1.2,mar[4]))
	brks <- seq(minx,maxx,len=round(par("pin")[1]*96))
	mids <- brks[-1]+diff(brks)/2
	rside <- par("usr")[rug+2]
	rsign <- ifelse(rug==1,1,-1)
	for (i in (1:length(dens))[order(sp$k.density)]) {
		plot(1,1, type="n", log=log, xlim=c(minx, maxx), ylim=ylim, main=paste0(main,"(Events=",i-1," x ",sp$k.density[i],", size=",lens[[i]],")"), xlab="", ylab="", xaxt="n", ...)
		points(dens[[i]], type=type, col=(i+1))
		idx <- round(L * (dens[[i]]$x-minxdens) / dx + 1)
		w <- which(diff(idx) > 0)
		yaprx[idx[w]] <- yaprx[idx[w]] + dens[[i]]$y[w] * sp$k.density[i]
		#if (rug) rug(sp$sampling[[i]], col=i+1, side=rug)
		if (rug) {
			h <- hist(sp$sampling[[i]],breaks=brks,plot=FALSE)$counts
			segments(mids[h > 0], rside, mids[h > 0], rside+rsign*h[h > 0]*abs(par("usr")[3])*rug.size/max(h), col=i+1)
		}
		mxy[i] <- max(dens[[i]]$y, na.rm=TRUE)
	}
	is.logy <- (log=="y" || log=="xy")
	leg <- c(0:sp$ne,"W. Sum","Sampling")
	lok <- 1:length(leg)
	if (!sampling) lok <- lok[-length(lok)]
	if (!w.sum) lok <- lok[-grep("W. Sum",leg)]
	#if (!add) legend(if(is.logy) "bottomleft" else "topleft",legend=leg[lok],lty=1,col=c(2+0:sp$ne,sp$ne+3,1)[lok], ncol=ncol, text.width=strwidth("0"),box.col=0)
	#if (!add) text(sp$m, pmin(my,maxy*.75), paste("x ",format(sp$k.density,digits=2)), srt=90, col=2+0:sp$ne, adj=c(ifelse(is.logy,1,0),1))
	if (sampling) {
		par(mar=c(2,mar[2],2,mar[4]))
		plot(1,1, type="n", log=log, xlim=c(minx, maxx), ylim=ylim, main=paste0(main,"Sampling"), xlab="", ylab="", ...)
		points(sp$empirical.breaks[-1]-diff(sp$empirical.breaks)/2, sp$empirical*mean(mxy)/max(sp$empirical, na.rm=TRUE), type=type, lwd=2, col=1, ...)
	}
	if (w.sum) {
		par(mar=c(2,mar[2],2,mar[4]))
		plot(1,1, type="n", log=log, xlim=c(minx, maxx), ylim=ylim, main=paste0(main,"W.Sum"), xlab="", ylab="", ...)
		points(xaprx, yaprx*mean(mxy)/max(yaprx,na.rm=TRUE), type=type, col=sp$ne+3, lwd=3)
	}
	axis(4,at=axTicks(2),round(max(sp$empirical)*empirical*axTicks(2)/max(axTicks(2)),5))
	par(mar=mar)
	invisible(list(x=xaprx, y=yaprx))
}

