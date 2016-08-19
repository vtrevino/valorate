
# Victor Trevino, Emmanuel Martinez, and Jose Tamez

# VALORATE is a procedure to accurately estimate the p-value
# of the difference in two survival curves using the log-rank 
# test specially in the cases of largely unbalanced groups.
# Instead of using a normal or chi-squrare, VALORATE estimates
# the null-distribution by a weighted sum of conditional
# distributions over a co-occurrence parameters.
# VALORATE was designed for cancer genomics where the
# comparisons between survival groups are heavily unbalanced
# since the frequency of gene mutations is quite low.
# Nevertheless, VALORATE should work for standard log-rank tests.
#
# This script shows an example for using VALORATE from:
# (1) mutation matrix
# (2) time and status survival information
#
# Basically, this script follows:
# (1) Create a Valorate object
# (2) Initialize variables
# (3) Compute all p-values (valorate and standard survidiff)
# (4) Print top 20 genes
# (5) Print summary information for top gene

####### DATA PARAMETERS #######
# data is a numeric matrix where rows contains genes and 
# columns represent samples. data contains 1 for mutations 
# and 0 otherwise.
# time is a numeric vector of survival times 
# status is a numeric vector for event 
# status is 1=event/death/recurrence etc, and 0=censored
# The order of data columns, time, and status are assumed
# to correspond in order to samples 1, 2, 3, ...
#######

####### VALORATE PARAMETERS #######
sampling.size <- 100000		## Amount of sampling (larger=better=slower=more memory needed)
min.sampling.size <- 1000   ## Minimum number of sampling (per co-occurrence item)
verbose <- TRUE				## Display tracing information in output console file.
min.mut  <- 4				## minimum number of mutations
tie.sampling <- 30			## How ties are treated ? ==> 0 : Do not sample ties, > 0 : # of tie samplings
valorateMethod <- "C"		## "C" or "R" (R is slower), this imply generating the library (.DLL or .SO)
							## Use R CMD SHLIB valorate_samplings.c in UNIX-LIKE (Linux/Mac OS X) or equivalent in Windows

## Use this code to use the file mutations-BRCA.txt containing data from BRCA TCGA.
surv <- read.delim("mutations-BRCA.txt", comment.char="", nrow=1, as.is=TRUE, row.names=1)
survival <- as.character(surv[1,])
brca <- read.delim("mutations-BRCA.txt", comment.char="#", as.is=TRUE, row.names=1)
ok <- which(!is.na(survival) & survival != "NA" & survival != "" & survival != "0" & survival != "0+")
data <- data.matrix(brca[,ok])
time <- as.numeric(gsub("\\+$","",survival[ok]))
status <- numeric(length(time))
status[] <- 1
status[grep("\\+$", survival[ok])] <- 0


## STEP 1: CREATE VALORATE OBJECT
source("valorate.R", chdir=TRUE) #assumes valorate.R and valorate_samplings library in the same directory
v.obj <- new.valorate(
	time=time, censored=1-status, 
	sampling.size=sampling.size, 
	min.sampling.size=min.sampling.size, 
	sampling.ties=tie.sampling,
	method=valorateMethod, 
	verbose=verbose)

## STEP 2: INITIALIZE RESULT VECTORS
n <- nrow(data)
p.v <- numeric(n)
p.v[] <- 1 # default to 1
p.sd <- p.v
risk <- p.v
nmut <- apply(data != 0, 1, sum)


## STEP 3: CALCULATE VALORATE & Common Chi-Square test (using SURVDIFF) FOR EVERY GENE
timestart.valorate <- Sys.time() #format(Sys.time(), "%Y-%m-%d %H:%M:%S")
for (i in 1:n) {
	genei <- data[i,]
	if (nmut[i] >= min.mut) {
		# Only estimate valorate for genes having a minimum of mutations
		try(p.v[i] <- valorate.survdiff(v.obj, genei))
		try(p.sd[i] <- 1-pchisq(survdiff(Surv(time, status) ~ genei)$chisq, 1))
		try(risk[i] <- c(valorate.risk(v.obj, genei)[1],1)[1])
	}
}
timeend.valorate <- Sys.time() #format(Sys.time(), "%Y-%m-%d %H:%M:%S")


## STEP 4: SHOW TOP 20 genes
o <- order(p.v)[1:20]
print(data.frame(Row=o, Gene=rownames(data)[o], p.Valorate=p.v[o], p.SurvDiff=p.sd[o], Risk=risk[o]))

## STEP 5: SHOW USEFUL INFO for top gene
valorate.survdiff(v.obj, data[o[1], ])
valorate.plot.empirical(v.obj, data[o[1], ])



## Bonus: Comparing running time
v.obj.R <- new.valorate(
	time=time, censored=1-status, 
	sampling.size=sampling.size, 
	min.sampling.size=min.sampling.size, 
	sampling.ties=tie.sampling,
	method="R", 
	verbose=verbose)

v.obj.C <- new.valorate(
	time=time, censored=1-status, 
	sampling.size=sampling.size, 
	min.sampling.size=min.sampling.size, 
	sampling.ties=tie.sampling,
	method="C", 
	verbose=verbose)

w <- sample(which(nmut == 4), 1)
system.time(valorate.survdiff(v.obj.R, data[w, ]))
system.time(valorate.survdiff(v.obj.C, data[w, ]))

w <- sample(which(nmut > 50), 1)
system.time(print(valorate.survdiff(v.obj.R, data[w, ])))
system.time(print(valorate.survdiff(v.obj.C, data[w, ])))
