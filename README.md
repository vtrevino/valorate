# Valorate
VALORATE is a procedure to accurately estimate the p-value of the difference in two survival curves using the log-rank test specially in the cases of largely unbalanced groups. Instead of using a normal or chi-squrare, VALORATE estimates the null-distribution by a weighted sum of conditional distributions over a co-occurrence parameter. VALORATE was designed for cancer genomics where the comparisons between survival groups are heavily unbalanced since the frequency of gene mutations is quite low. Nevertheless, VALORATE should work for standard log-rank tests.

# Running Time
VALORATE estimations of null distribution are quite fast compared to other methods. Its running time is within the order of seconds or less (for 1 null distribution) depending on the "sampling.size", which is the main parameter. So, for non systematic or sporadic calculations, the running time is insignificant. Nevertheless, for systematic analyses such as in genomics, it may take few minutes. For example, for the BRCA dataset included and the default parameters (minimum mutations = 4), there are 5940 genes and 965 individuals. Running the 'run-valorate.R' in a Mac i7 2.5Ghz, the estimations take less than 8 minutes using the "core" calculations in C (which uses valorate_sampling.c). To estimate these 5940 p-values, VALORATE performed 48 rounds of null distribution estimations corresponding to genes mutated between 4 and 312 individuals. The running time using the "core" calculations in R (parameter method="R") is around 100 times slower (between 70 and 120 depending on the number of mutations -or lowest group size-), however this is useful for users not minded to build external C libraries.

# Building C libraries
To speed up calculations, it is recommended to use the C code. For this, it is needed to compile the code valorate_sampling.c and that the built library (.so or .dll) is located is the same directory than valorate.R (and using chdir=TRUE in the source calling). If having problems or require further information, please read the R writing extension help. In the following, it is assumed that R is already installed.

Mac OS X: 
Depending on your OS version, it may be needed the XCode and Command line utilities that can be obtained free from Apple.
For building use:
- Open a Terminal
- Change directory to the directory where valorate_sampling.c is
- type "R CMD SHLIB valorate_sampling.c <ENTER>"
- if there is a file valorate_sampling.so, everthing is ok.

Linux: 
It should work without additional installations. But here is more difficult to tell what would be needed because of all flavours of Linux. However, it should be straight forward for a common linux user.
- Open a Terminal
- Change directory to the directory where valorate_sampling.c is
- type "R CMD SHLIB valorate_sampling.c <ENTER>"
- if there is a file valorate_sampling.so, everthing is ok.

Windows: 
You need Rtools (for example see http://mcglinn.web.unc.edu/blog/linking-c-with-r-in-windows/).
- Open a Terminal
- Change directory to the directory where valorate_sampling.c is
- type "R CMD SHLIB valorate_sampling.c <ENTER>"
- if there is a file valorate_sampling.dll, everthing is ok.

# Example Running VALORATE
- Place all files in a local directory.
- Optional but highly recommended: Build C libraries
- Start R
- Run "run-valorate.R" (source("run-valorate.R"))

