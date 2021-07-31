# Valorate
VALORATE is a procedure to accurately estimate the p-value of the difference in two survival curves using the log-rank test specially in the cases of largely unbalanced groups. Instead of using a normal or chi-squrare, VALORATE estimates the null-distribution by a weighted sum of conditional distributions over a co-occurrence parameter. VALORATE was designed for cancer genomics where the comparisons between survival groups are heavily unbalanced since the frequency of gene mutations is quite low. Nevertheless, VALORATE should work for standard log-rank tests. VALORATE was developed in [Victor Trevino Lab] (http://bioinformatica.mty.itesm.mx/Valorate [broken]) at [Tecnológico de Monterrey] (http://www.itesm.mx) ¡Viva México!. 

[NOTE: The newest VALORATE R package is available through CRAN https://cran.r-project.org/web/packages/valorate/index.html].

# License
See License file.

# Running Time
VALORATE estimations of null distribution are quite fast compared to other methods. Its running time is within the order of seconds or less (for 1 null distribution) depending on the "sampling.size", which is the main parameter. So, for non systematic or sporadic calculations, the running time is insignificant. Nevertheless, for systematic analyses such as in genomics, it may take few minutes. For example, for the BRCA dataset included and the default parameters (minimum mutations = 4), there are 5940 genes and 965 individuals. Running the 'run-valorate.R' in a Mac i7 2.5Ghz, the estimations take less than 8 minutes using the "core" calculations in C (which uses valorate_sampling.c). To estimate these 5940 p-values, VALORATE performed 48 rounds of null distribution estimations corresponding to genes mutated between 4 and 312 individuals. The running time using the "core" calculations in R (parameter method="R") is around 100 times slower (between 70x and 120x depending on the number of mutations -or lowest group size-), however this is useful for users not minded to build external C libraries.

# Building C library
To speed up calculations, it is recommended to use the C code. For this, it is needed to compile the code valorate_sampling.c and that the built library (.so or .dll) is located is the same directory than valorate.R (and using chdir=TRUE in the source calling). If having problems or require further information, please read the R writing extension help. In the following, it is assumed that R is already installed.

Mac OS X: 
Depending on your OS version, it may be needed the XCode and Command line utilities that can be obtained free from Apple. First try without any installation, if fails, install XCode and Command line utilities.


Linux: 
Here is more difficult to tell what would be needed because of all flavours of Linux. However, it should work without additional installations. So, it should be straight forward for a common linux user.

Windows: 
You need Rtools (for example see http://mcglinn.web.unc.edu/blog/linking-c-with-r-in-windows/).

For building use:
- Open a Terminal
- Change directory to the directory where valorate_sampling.c is (don´t forget to replace the <..> for your actual folder)

    `cd <valorate directory>`

- Compile the C code within valorate (typing "R CMD SHLIB valorate_sampling.c <ENTER>" as below)

    `R CMD SHLIB valorate_sampling.c`

- Linux and Mac: if there is a file valorate_sampling.so, everthing is ok.

    `ls -l valorate_sampling.*`

- Windows:       if there is a file valorate_sampling.dll, everthing is ok.

    `dir valorate_sampling.*`

# Example Running VALORATE
- Place all files in a local directory (valorate.R, valorate-samplings.c, run-valorate.R, and the uncompressed mutations-BRCA.txt)
- Change default directory to above directory (don't forget to replace the <...>)

    `cd <dir>`

- Optional but highly recommended: Build C library (see steps above)

- Start R

    `R`

- Run "run-valorate.R" within R

    `source("run-valorate.R")`
    

# Known ISSUES
[October 9th 2016] The C code has a bias in the random selection of samples. The impact is very subtle and it has been corrected in the R package submitted to CRAN.
It is recommended to use the R package in pipelines.
See http://bioinformatica.mty.itesm.mx/valorateR
    
[NOTE: The newest VALORATE R package is available through CRAN https://cran.r-project.org/web/packages/valorate/index.html].
