// v1.3 - Use of all combinations array
// v1.2 - Weights implementation
// v1.1 - Ties implementation
// v1.0 - Initial development

/////////////////////////////
// COMPILATION IN MAC OS X

/////////////////////////////
// COMPILATION IN WINDOWS
// - Install Rtools
// - Include in PATH the CurrentRversion/bin directory
// - Compile using:
//		R CMD SHLIB valorate_samplings.c

#include <math.h>
#include <R.h>
#include <Rinternals.h>
//#include <Rinterface.h>
#include <R_ext/Lapack.h>
#include <Rmath.h>
#include <time.h>       /* time */


/*
SEXP valorate_sampling(SEXP s_v, SEXP s_sim, SEXP s_n, SEXP s_k, SEXP s_nx, SEXP s_wcensored, SEXP s_wevents, SEXP s_vcjx, SEXP s_inn1, SEXP s_ldx, SEXP s_randomize) {

	int k = INTEGER(s_k)[0];
	int nx = INTEGER(s_nx)[0];
	int sim = INTEGER(s_sim)[0];
	int n = INTEGER(s_n)[0];
	int* inn1 = INTEGER(s_inn1);
	int* ldx = INTEGER(s_ldx);
	double* v = REAL(s_v);
	int* wcensored = INTEGER(s_wcensored);
	int ncensored = LENGTH(s_wcensored);
	int* wevents = INTEGER(s_wevents);
	int nevents = LENGTH(s_wcensored);

	if (INTEGER(s_randomize)[0]) {
		srand(time(NULL));
	}
*/

void valorate_sampling(double* v, int* psim, int* pn, int* pk, int* pnx, int* wcensored, int* pncensored, int *wevents, int* pnevents, double* weightedevents, double* vcjx, int* pvcjx_n, int* inn1, int* ldx, int* prandomize, int* pdebug, int* allComb) {

	int i,j,m,r;
	double V;
	int sim = *psim;
	int n = *pn;
	int k = *pk;
	int nx = *pnx;
	int ncensored = *pncensored;
	int nevents = *pnevents;
	int debug = *pdebug;
	int vcjx_n = *pvcjx_n;
	int offset;
	//int matrixlen = nevents * nx;

	if (*prandomize) {
		srand(time(NULL));
	}

	int doAllCombinations = (allComb[0] > 0);
	int acK = 0;

	for (i=0; i < sim; i++) {

		if (debug && (i+1) % 1000 == 0) {
			Rprintf(".");
			if (debug && (i+1) % 10000 == 0) {
			   Rprintf(" ");
			}
		}

		// inn1[] <- 0

		for (j=0; j < n; j++) inn1[j] = 0;

		if (doAllCombinations) {
			for (j=0; j < nx; j++) {
				inn1[allComb[acK++]-1] = 1;
			}
		} else {
			if (k < nx) {			
				// Censored
				// inn1[sample(wcensored,nx-k)] <- 1
				m = nx-k;
				if (ncensored < m) m = ncensored; // FOR SECURITY and to avoid cycling
				for (j=0; j < m; j++) {
					for (r = rand() % ncensored; inn1[wcensored[r]-1] == 1; r = (r+1) % ncensored);
					inn1[wcensored[r]-1] = 1;
				}
			}
			if (k > 0) {
				// Events
				// inn1[sample(wevents,k)] <- 1
				m = k;
				if (nevents < k) m = nevents; // FOR SECURITY and to avoid cycling 
				for (j=0; j < m; j++) {
					for (r = rand() % nevents; inn1[wevents[r]-1] == 1; r = (r+1) % nevents);
					inn1[wevents[r]-1] = 1;
				}
			}			
		}
		// Calculate the V statistic
		// einn1 <- inn1[wevents] // einn1 is not necessary because array is indexed as 0 but it is needed in wevents
		// ldx <- nx - cumsum(c(0,inn1))[wevents] // ldx will be increased dinamically instead of saving an array
		// V <- 0
		// for (j in 1:events) {
		// 		if (ldx[j] > 0) V <- V + einn1[j] - vcjx[ldx[j],j]
		// }
		ldx[0] = nx - 1; // in 0 indexes
		for (j=1; j < n; j++) {
			ldx[j] = ldx[j-1] - inn1[j-1];
		}
		V = 0;
		offset = (vcjx_n == 1 ? 0 : (rand() % vcjx_n) * nevents);
		for (j=0; j < nevents && ldx[r=(wevents[j]-1)] >= 0; j++) {
			//m = ldx[r] + j*nx + inn1[r]*matrixlen;
			//if (m < 0 || m >=matrixlen*2) {
			//	Rprintf("****AGUAS***************************\n");
			//	Rprintf("*******************************\n m=%d ",m);
			//	Rprintf("j=%d,r=%d,wevents[j]=%d,ldr[r]=%d,inn1[r]=%d,matrixlen=%d\n",j,r,wevents[j],ldx[r],inn1[r],matrixlen);
			//}
			V += weightedevents[j] * (inn1[r]-vcjx[ldx[r] + (j+offset)*nx]); //vcjx[ldx[r] + j*nx + inn1[r]*matrixlen];
		}
		v[i] = V;
	}

	//return (s_v);
}
