// UNCMND_ak1.CPP - Program for minimizing a function of several variables (with no constraints.)
// Can also be used for non-linear Least Squares Data-fitting.
// Written in Microsoft Visual Studio Express 2013 for Windows Desktop
// 5 March 2017
//
// This program is a translation of the FORTRAN routine UNCMND
// written by Stephen Nash, George Mason University.
// From the book "Numerical Methods and Software"
// D. Kahaner, C. Moler, and S. Nash
// Prentice Hall, 1988
//
// To distinguish this program from others, an '_ak1' suffix has been appended to its name.
//
// In this program, data is input from a file to eliminate the need for a user to type data in via
// the console.

#include <iostream>
#include <fstream>
#include <cctype>
#include <cmath>
#include <vector>
#include <cfloat>

#define NUMCOEFF 20

using namespace std;

typedef vector<double> C1DArray;
typedef vector <vector<double>> C2DArray;

double fObj(int N, double* vec, C1DArray x);
double Euclid2Norm(int N, C1DArray inVec);
double sign(double y);
void FSTFDD_ak1(int N, double* vec, C1DArray& x, C1DArray& gradient_vec, double RNFVAL, double FVAL);
void FSTCDD_ak1(int N, double* vec, C1DArray& x, C1DArray& gradient_vec, double RNFVAL);
void LLTSLD_ak1(int N, C1DArray& pvec, C1DArray gvec, C2DArray A);
void LNSRCD_ak1_ak1(int N, C1DArray x0, C1DArray& x, C1DArray g, C1DArray& p, int *iretcd, bool *mxtake, double *tempfx, double STEPTL, double STEPMX, double* vec, double fxn0);
void OPTSTD_ak1(int N, C1DArray p, C1DArray x, C1DArray gtemp, int iretcd, int *itrmcd, int *icscmx, double tempfx, double GRADTL, double STEPTL, int MAXIT, bool mxtake, int itncnt);
void QRUPDD_ak1(int N, C1DArray& UVEC, C1DArray WVEC, C2DArray& A);
void uncmnd_ak1(int N, double* vec, C1DArray& x, int *info, double *fx, C1DArray x0, C1DArray p, C1DArray g, C1DArray gtemp, C2DArray& A, C1DArray YVEC, C1DArray UVEC, C1DArray WVEC);

double fObj(int N, double* vec, C1DArray x) {
	// The Objective Function: the Generalized Rosenbrock Function for N = 10.
	// N: the number of dimensions (i.e., number of independent variables.)
	// vec: the coefficient vector: not used in the present example
	// x: the present value of x

	double dummy = 1.0, temp1, dum1;

	for (int i = 1; i < N; ++i){
		temp1 = x[i] - x[i - 1] * x[i - 1];
		dum1 = 1.0 - x[i - 1];
		dummy += 100.0 * temp1 * temp1 + dum1 * dum1;
	} // End for i

	return dummy;
} //End fObj

double Euclid2Norm(int N, C1DArray inVec){

	double absxi, dummy, scale = 0.0, ssq = 1.0;

	for (int i = 0; i < N; ++i){
		if (inVec[i] != 0){
			absxi = fabs(inVec[i]);
			dummy = scale / absxi;
			ssq = 1.0 + ssq*dummy*dummy;
			scale = absxi;
		}//End if (inVec[i] != 0)
	}//End for i

	return scale*sqrt(ssq);

} //End Euclid2Norm

double sign(double y)	//If y < 0, return -1, else +1
{
	return ((y < 0) ? -1.0 : 1.0);
}

void FSTFDD_ak1(int N, double* vec, C1DArray& x, C1DArray& gradient_vec, double RNFVAL, double FVAL){
	// Computes Forward Difference gradient at position given by input vector x

	double dumfx, x_temp, stpsiz;

	for (int i = 0; i < N; ++i){
		x_temp = x[i];
		stpsiz = fabs(x_temp);
		if (stpsiz < 1.0) stpsiz = 1.0;
		stpsiz *= RNFVAL;
		x[i] += stpsiz;
		dumfx = fObj(N, vec, x);
		x[i] = x_temp;
		gradient_vec[i] = (dumfx - FVAL) / stpsiz;
	} // End for i

	return;
} // End FSTFDD_ak1

void FSTCDD_ak1(int N, double* vec, C1DArray& x, C1DArray& gradient_vec, double RNFVAL){
	// Computes Central Difference gradient at position given by input vector x

	double fxminus, fxplus, x_temp, stpsiz;

	for (int i = 0; i < N; ++i){
		x_temp = x[i];
		stpsiz = fabs(x_temp);
		if (stpsiz < 1.0) stpsiz = 1.0;
		stpsiz *= RNFVAL;
		x[i] += stpsiz;
		fxplus = fObj(N, vec, x);
		x[i] = x_temp - stpsiz;
		fxminus = fObj(N, vec, x);
		x[i] = x_temp;
		gradient_vec[i] = (fxplus - fxminus) / (2.0*stpsiz);
	} // End for i

	return;
} // End FSTCDD_ak1

void LLTSLD_ak1(int N, C1DArray& pvec, C1DArray gvec, C2DArray A){
	// Solve for Newton step: AP = -G

	double temp;
	int i, j;

	// FORSLD
	// Solve LX = B (Forward Solve); put result in pvec

	pvec[0] = -gvec[0] / A[0][0];
	for (i = 1; i < N; ++i) {
		temp = 0.0;
		for (j = 0; j < i; ++j)
			temp += A[i][j] * pvec[j];
		pvec[i] = -(gvec[i] + temp) / A[i][i];
	} // End for i

	// End of FORSLD

	// BAKSLD
	// Solve (L-Transpose)X = B (Back Solve); put result in pvec

	pvec[N - 1] /= A[N - 1][N - 1];
	for (i = (N - 2); i >= 0; --i) {
		temp = 0.0;
		for (j = (N - 1); j > i; --j)
			temp += A[j][i] * pvec[j];
		pvec[i] = (pvec[i] - temp) / A[i][i];
	} // End for i

	// End of BAKSLD

	return;
} // End LLTSLD_ak1

void LNSRCD_ak1(int N, C1DArray x0, C1DArray& x, C1DArray g, C1DArray& p, int *iretcd, bool *mxtake, double *tempfx, double STEPTL, double STEPMX, double* vec, double fxn0){
	// Find a next Newton iterate by line search

	bool fFlag = 1; // Flag to indicate first time through loop
	int i;
	double a = 0.0, ALMSTMX = 0.99*STEPMX, b, disc, almbda = 1.0, plmbda, pfpls, RMNLMB, SLN, SLP = 0.0, temp, THOSLP, tlmbda, t1, t2, t3;

	*mxtake = 0;
	*iretcd = 2;

	SLN = Euclid2Norm(N, p);

	if (SLN > STEPMX) { // Newton Step larger than maximum allowed
		t1 = STEPMX / SLN;
		for (i = 0; i < N; ++i)  p[i] *= t1;	// SCLMLD
		SLN = STEPMX;
	} // End if (SLN > STEPMX)

	// Calculate the dot product of vectors G and P
	for (i = 0; i < N; ++i)  SLP += g[i] * p[i];
	THOSLP = 0.0001*SLP;

	for (i = 0; i < N; ++i) {
		t1 = fabs(x0[i]);
		if (t1 < 1.0) t1 = 1.0;
		temp = fabs(p[i]) / t1;
		if (a < temp) a = temp;
	} // End for i

	RMNLMB = STEPTL / a;

	// Loop. Check if new iterate satisfactory. Generate new lambda if necessary
	do {

		for (i = 0; i < N; ++i)  x[i] = x0[i] + almbda*p[i];

		*tempfx = fObj(N, vec, x);

		if ((*tempfx) <= (fxn0 + THOSLP*almbda)) { // Solution found
			*iretcd = 0;
			if ((fFlag) && (SLN > ALMSTMX)) *mxtake = 1;
			break;
		} // End if (*tempfx <= (fxn0))

		// Solution not found yet

		if (almbda < RMNLMB) { // No satisfactory x found sufficiently distinct from x0
			*iretcd = 1;
			break;
		} // End if (almbda < RMNLMB)

		// Calculate new lambda
		if (fFlag) {  // First backtrack: quadratic fit
			tlmbda = *tempfx - (fxn0 + SLP);
			tlmbda = -0.5 * SLP / tlmbda;
			fFlag = 0;
		} // End if (fFlag)
		else { //else (!fFlag); all subsequent backtracks: cubic fit
			t1 = *tempfx - (fxn0 + almbda*SLP);
			t2 = pfpls - (fxn0 + plmbda*SLP);
			t3 = 1.0 / (almbda - plmbda);
			a = t3 * (t1 / (almbda*almbda) - t2 / (plmbda*plmbda));
			b = t3 * (t2*almbda / (plmbda*plmbda) - t1*plmbda / (almbda*almbda));
			disc = b*b - 3.0*a*SLP;

			if (disc > b*b)  // Only one positive critical point, must be a minimum
				tlmbda = sign(a)*sqrt(disc) - b;
			else  // else (disc <= b*b). Both critical points positive, first is minimum
				tlmbda = -(b + sign(a)*sqrt(disc));

			tlmbda /= (a*3.0);

			if (tlmbda > 0.5*almbda) tlmbda = 0.5*almbda;
		} // End else (!fFlag)

		plmbda = almbda;
		pfpls = *tempfx;

		if (tlmbda < 0.1*almbda) almbda *= 0.1;
		else almbda = tlmbda;

	} while (*iretcd >= 2); // End do-while (*iretcd >= 2)

	return;
} // End LNSRCD_ak1

void OPTSTD_ak1(int N, C1DArray p, C1DArray x, C1DArray gtemp, int iretcd, int *itrmcd, int *icscmx, double tempfx, double GRADTL, double STEPTL, int MAXIT, bool mxtake, int itncnt){
	// Check whether stopping criteria satisfied

	double grdtemp, rgx = 0.0, rsx = 0.0, stptemp, xtmpj = fabs(tempfx);

	*itrmcd = 0;

	if (iretcd == 1) { // Last global step failed to locate a point lower than x0
		*itrmcd = 3;
	}
	else { //else (iretcd != 1)
		// Find direction in which relative gradient maximum.
		// Find direction in which relative stepsize maximum
		// Check whether these quantities within tolerance

		if (xtmpj < 1.0) xtmpj = 1.0;

		for (int i = 0; i < N; ++i) {
			grdtemp = fabs(x[i]);
			if (grdtemp < 1.0) grdtemp = 1.0;
			stptemp = grdtemp;
			grdtemp *= fabs(gtemp[i]) / xtmpj;
			stptemp = fabs(p[i]) / stptemp;
			if (rgx < grdtemp) rgx = grdtemp;
			if (rsx < stptemp) rsx = stptemp;
		} // End for i

		if (rgx <= GRADTL){
			*itrmcd = 1;
		}
		else { // else (rgx > GRADTL)

			if (rsx <= STEPTL) {
				*itrmcd = 2;
			}
			else { // else rsx > STEPTL

				// Check iteration limit

				if (itncnt >= MAXIT) {
					*itrmcd = 4;
				}
				else { // else itncnt < MAXIT

					// Check number of consecutive steps/ STEPMX

					if (mxtake) {
						++(*icscmx);
						if ((*icscmx) >= 5) *itrmcd = 5;
					} // End if (mxtake)
					else
						*icscmx = 0;

				} // end else itncnt < MAXIT

			} // end else rsx > STEPTL

		} // End else (rgx > GRADTL)

	} // End else (iretcd != 1)

	return;
} // End of OPTSTD_ak1

void QRUPDD_ak1(int N, C1DArray& UVEC, C1DArray WVEC, C2DArray& A){

	int i, j, JJ, LLL;
	double QDUM, QC, QS, QY, QZ, T1, T2;

	// Determine last non-zero in U.
	for (i = (N - 1); i > 0; --i) if (UVEC[i] != 0.0) break;

	// (K - 1) Jacobi Rotations Transform

	for (j = (i - 1); j >= 0; --j){

		if (UVEC[j] != 0.0){

			// QRAX2D

			QDUM = sqrt(UVEC[j] * UVEC[j] + UVEC[j + 1] * UVEC[j + 1]);
			QC = UVEC[j] / QDUM;
			QS = -UVEC[j + 1] / QDUM;

			for (LLL = j; LLL < N; ++LLL){

				QY = A[j][LLL];
				QZ = A[j + 1][LLL];
				A[j][LLL] = QC*QY - QS*QZ;
				A[j + 1][LLL] = QS*QY + QC*QZ;

			} // End for LLL

			// End QRAX2D

			UVEC[j] = sqrt(UVEC[j] * UVEC[j] + UVEC[j + 1] * UVEC[j + 1]);

		} // End if (UVEC[j] != 0.0)
		else { // else (UVEC[j] == 0.0)

			// QRAX1D

			for (JJ = j; JJ < N; ++JJ){

				QDUM = A[j][JJ];
				A[j][JJ] = A[j + 1][JJ];
				A[j + 1][JJ] = QDUM;

			} // End for JJ

			// End QRAX1D

			UVEC[j] = UVEC[j + 1];

		} // End else (UVEC[j] == 0.0)

	} // End for j

	QDUM = UVEC[0];
	for (j = 0; j < N; ++j) A[0][j] += QDUM * WVEC[j];

	// (K - 1) Jacobi Rotations transform Upper Hessenberg R to Upper Triangular (R*)

	for (j = 0; j < i; ++j){

		if (A[j][j] != 0.0){

			T1 = A[j][j];
			T2 = -A[j + 1][j];

			// QRAX2D

			QDUM = sqrt(T1 * T1 + T2 * T2);
			QC = T1 / QDUM;
			QS = T2 / QDUM;

			for (LLL = j; LLL < N; ++LLL){

				QY = A[j][LLL];
				QZ = A[j + 1][LLL];
				A[j][LLL] = QC*QY - QS*QZ;
				A[j + 1][LLL] = QS*QY + QC*QZ;

			} // End for LLL

			// End QRAX2D

		} // End if (A[j][j] != 0.0)
		else { // else (A[j][j] == 0.0)

			// QRAX1D

			for (JJ = j; JJ < N; ++JJ){

				QDUM = A[j][JJ];
				A[j][JJ] = A[j + 1][JJ];
				A[j + 1][JJ] = QDUM;

			} // End for JJ

			// End QRAX1D

		} // End else (A[j][j] == 0.0)

	} // End for j

	return;
} // End of QRUPDD_ak1

void uncmnd_ak1(int N, double* vec, C1DArray& x, int *info, double *fx, C1DArray x0, C1DArray p, C1DArray g, C1DArray gtemp, C2DArray& A, C1DArray YVEC, C1DArray UVEC, C1DArray WVEC) {
	// This subroutine minimizes a smooth nonlinear function of N independent variables
	//
	// Termination Codes:
	// Info = 1: Terminated with gradient small, x is probably optimal
	// Info = 2: Terminated with stepsize small, x is probably optimal
	// Info = 3: Lower point cannot be found, x is probably optimal
	// Info = 4: Iteration limit (150) exceeded
	// Info = 5: Too many large steps, function may be unbounded

	// fx is the function value at x (the minimum that we sought)

	const int MAXIT = 150; // Maximum number of iterations that will be allowed.
	int i, icscmx, iretcd, itncnt = 0, j, k;
	bool fDifFLG = 1, mxtake;
	bool NOUPDT = 1, NSDIAG = 0, SKPUPD;		// Variables for SECFCD
	double ALP, DEN1, RELTOL, SNORM2, YNRM2;	// Variables for SECFCD
	double stpsiz, temp, tempfx = 0.0;

	const double EPSMCH = DBL_EPSILON; //Machine epsilon for type double
	const double GRADTL = pow(EPSMCH, 1.0 / 3.0);
	const double STEPTL = sqrt(EPSMCH);

	iretcd = -log10(EPSMCH);
	temp = pow(10.0, -iretcd);
	if (temp < EPSMCH) temp = EPSMCH;

	const double RNF = temp;
	const double RT2RNF = sqrt(RNF);
	const double RNFRT13 = pow(RNF, 1.0 / 3.0);

	// Check maximum step size
	stpsiz = Euclid2Norm(N, x0);
	if (stpsiz < 1.0) stpsiz = 1000.0;
	else stpsiz *= 1000.0;

	const double STEPMX = stpsiz;

	*fx = fObj(N, vec, x0);
	cout << "\nAt x0, the Objective Function = " << *fx << " \n\n";

	// Compute the gradients at x0; save these values in g array
	FSTFDD_ak1(N, vec, x0, g, RT2RNF, *fx);

	// OPTSTD

	// Find direction in which relative gradient maximum.
	// Check whether within tolerance

	DEN1 = fabs(*fx);
	if (DEN1 < 1.0) DEN1 = 1.0;

	for (i = 0; i < N; ++i) {
		temp = fabs(x0[i]);
		if (temp < 1.0) temp = 1.0;
		temp *= fabs(g[i]) / DEN1;
		if (tempfx < temp) tempfx = temp;
	} // End for i

	if (tempfx <= GRADTL) {
		*info = 1;
		for (i = 0; i < N; ++i) x[i] = x0[i];
		return;
	} // End if (tempfx <= GRADTL)

	// End of OPTSTD

	// Assume function is expensive to evaluate (IEXP = 1), so Hessian will be obtained by secant updates.

	// Main iterative loop
	for (k = 0; k <= MAXIT; ++k) {

		if (NOUPDT) for (i = 0; i < N; ++i) p[i] = -g[i];
		else LLTSLD_ak1(N, p, g, A); // Solve for Newton step: AP = -G

		// Decide whether to accept Newton Step X = X0 + P or to choose
		// X by a global strategy

		LNSRCD_ak1(N, x0, x, g, p, &iretcd, &mxtake, &tempfx, STEPTL, STEPMX, vec, *fx);

		if ((iretcd == 1) && (fDifFLG)) {

			// Could not find a satisfactory step using Forward Difference Gradient. Retry using Central Difference Gradient

			FSTCDD_ak1(N, vec, x0, g, RNFRT13);
			fDifFLG = 0;
			continue;	// Skip rest of k loop

		} // End if (iretcd == 1)

		++itncnt;

		// Calculate step for output
		for (i = 0; i < N; ++i) p[i] = x[i] - x0[i];

		// Calculate Gradient at x
		// Store this gradient in gtemp

		if (fDifFLG) FSTFDD_ak1(N, vec, x, gtemp, RT2RNF, tempfx);
		else FSTCDD_ak1(N, vec, x, gtemp, RNFRT13);

		// Check whether stopping criteria satisfied

		OPTSTD_ak1(N, p, x, gtemp, iretcd, info, &icscmx, tempfx, GRADTL, STEPTL, MAXIT, mxtake, itncnt);

		if (*info != 0) break;  // break out of for k loop

		// SECFCD	*************************************************************
		
		// Before A is changed, NOUPDT = 1; simplified calculations can be done.
		// After A updated, full routine must be done.

		for (i = 0; i < N; ++i) YVEC[i] = gtemp[i] - g[i];

		// Compute Dot Product of vectors p and Y
		DEN1 = 0.0;
		for (i = 0; i < N; ++i) DEN1 += p[i] * YVEC[i];

		SNORM2 = Euclid2Norm(N, p);
		YNRM2 = Euclid2Norm(N, YVEC);

		if (DEN1 >= (STEPTL*SNORM2*YNRM2)){

			ALP = sqrt(DEN1);

			if (NOUPDT){

				ALP /= SNORM2;

				for (i = 0; i < N; ++i){
					UVEC[i] = ALP*p[i];
					WVEC[i] = ALP*UVEC[i];
					A[i][i] = ALP;
				} // End for i

				ALP = 1.0;
				NOUPDT = 0;

			} // End if (NOUPDT)
			else { // else (!NOUPDT)

				// MVMLUD

				for (i = 0; i < N; ++i){
					temp = 0.0;
					for (j = i; j < N; ++j) temp += A[j][i] * p[j];
					UVEC[i] = temp;
				} // End for i

				// End of MVMLUD

				// MVMLLD

				for (i = 0; i < N; ++i){
					temp = 0.0;
					for (j = 0; j <= i; ++j) temp += A[i][j] * UVEC[j];
					WVEC[i] = temp;
				} // End for i

				// End of MVMLLD

				ALP /= Euclid2Norm(N, UVEC);

			} // End else (!NOUPDT)

			if (fDifFLG) RELTOL = RT2RNF;
			else RELTOL = RNF;

			SKPUPD = 1;

			i = 0;
			do { //  LINE 60 IN FORTRAN CODE
				temp = fabs(g[i]);
				if (temp < fabs(gtemp[i])) temp = fabs(gtemp[i]);
				if (fabs(YVEC[i] - WVEC[i++]) >= RELTOL*temp) SKPUPD = 0;
			} while ((SKPUPD) && (i < N));

			if (!SKPUPD){

				for (i = 0; i < N; ++i) WVEC[i] = YVEC[i] - ALP * WVEC[i];

				ALP /= DEN1;

				for (i = 0; i < N; ++i) UVEC[i] *= ALP;

				// Copy L into Upper Triangular part. Zero L
				// This step can be skipped the first time through--assuming A initialized to Identity Matrix

				if (NSDIAG){
					for (i = 1; i < N; ++i) {

						for (j = 0; j < i; ++j) {

							A[j][i] = A[i][j];
							A[i][j] = 0.0;

						} // End for j

					} // End for i
				} // End if NSDIAG
				else NSDIAG = 1;

				QRUPDD_ak1(N, UVEC, WVEC, A);

				// Upper Triangular part and diagonal of A now contain updated Cholesky Decomposition of Hessian. 
				// Copy back to Lower Triangular part.

				for (i = 1; i < N; ++i) {
					for (j = 0; j < i; ++j) A[i][j] = A[j][i];
				} // End for i

			} // End if (!SKPUPD)

		} // End if (DEN1 >= (STEPTL*SNORM2*YNRM2))

		// End of SECFCD	*************************************************************

		// Assign xpls to x0, gpls to g, and fpls to *fx
		*fx = tempfx;
		for (i = 0; i < N; ++i) {
			x0[i] = x[i];
			g[i] = gtemp[i];
		} // End for i

	} // End of for k loop; main iterative loop

	// Termination. Reset xpls, fpls, gpls if previous iterate solution

	if (*info == 3) {
		for (i = 0; i < N; ++i) x[i] = x0[i];
	} // End if (*info == 3)

	return;
} // End uncmnd

int main()
{
	char rflag = 0;	//Readiness flag 

	cout << "                     UNCMND_ak1   (5 March 2017)\n";
	cout << "=========================================================================== \n";
	cout << "This program minimizes the Generalized Rosenbrock Function for N = 10\n";
	cout << "(without constraints):\n\n";
	cout << "F = 1.0 + SUMMATION from i = 2 to i = N of the following:\n";
	cout << "[ 100 * (x_i - x_(i-1)^2)^2 + (1 - x_(i-1))^2 ]\n\n";
	
	cout << "Data should have been saved beforehand in a file named uncmnd_in.txt, which\n";
	cout << "should be in the same folder as the uncmnd_ak1 executable. The\n";
	cout << "first entry in this file should be N, the number of dimensions (10).\n";
	//cout << "The entries for the function should follow:\n";
	//cout << "A, B, C, ...\n\n";
	
	cout << "The program assumes an initial guess for the solution has been provided.\n";
	cout << "This data should come next. \n\n";

	cout << "\nThe data is assumed to be of type double. Variables used within this program\n";
	cout << "are type double.\n";
	cout << "\nThe solution is written to the file uncmnd_out.txt.\n";
	cout << "\nIs everything ready (are you ready to continue?)? If yes, Enter y. \n";
	cout << "Otherwise Enter any other key. \n";
	cin >> rflag;

	if (toupper(rflag) == 'Y') {

		double coeffVec[NUMCOEFF];			// Holds values of the function: A, B, C, D, . . .
		C1DArray g, gtemp, p, x0, xVec;		// g, gtemp, p, x0, and x vectors
		C1DArray YVEC, UVEC, WVEC;		// Y, U, and W vectors for SECFCD.
		C2DArray HessMat;					// Hessian Matrix
		int i, info, rDim;
		double fv;							// The function value

		ifstream in("uncmnd_in.txt", ios::in);

		if (!in) {
			cout << "Cannot open the input file.\n";
			return 0;
		}

		in >> rDim;  //Input rDim, the number of independent variables N, from the file

		cout << "\nThis problem is one of " << rDim << " dimensions.\n\n";

		if (rDim <= 0)  {
			cout << "Invalid N dimension entered. Program terminated. \n";
			in.close();
			return 0;
		}

		if (rDim == 1)  {
			cout << "For a single-dimensional problem, this utility is \n";
			cout << "not the best choice. It would be better to use a utility\n";
			cout << "that is specifically written for one-dimensional problems.\n";
			cout << "Program terminated. \n";
			in.close();
			return 0;
		}

		ofstream out("uncmnd_out.txt", ios::out);
		if (!out) {
			cout << "Cannot open the output file. Program terminated.\n";
			in.close(); //Close the input file before terminating
			return 0;
		}

		try { // Resize the g, gtemp, p, x0, and xVec arrays to their required sizes
			p.resize(rDim);
			g.resize(rDim);
			gtemp.resize(rDim);
			x0.resize(rDim);
			xVec.resize(rDim);

			YVEC.resize(rDim);
			UVEC.resize(rDim);
			WVEC.resize(rDim);

			HessMat.resize(rDim);
			for (i = 0; i < rDim; ++i) HessMat[i].resize(rDim);

		} // End of try block

		catch (bad_alloc& xa) { // Catch block, for exceptions
			in.close();
			out.close();
			cerr << "In catch block for resizing x0 and xVec arrays: " << xa.what() << "\n";
			cout << "\nEnter any key to continue. \n";
			cin >> rflag;
			return 0;
		} // End of catch block

		//for (i = 0; i < NUMCOEFF; ++i)  //Input the function coefficients to array
		//	in >> coeffVec[i];

		for (i = 0; i < rDim; ++i)  //Input the x0 vector from the file
			in >> x0[i];

		in.close();  //Close the input file

		out.precision(DBL_DIG);
/*
		cout << "The function coefficients follow:\n";
		cout << "A = " << coeffVec[0] << "\t\tB = " << coeffVec[1] << "\n";
		cout << "C = " << coeffVec[2] << "\t\tD = " << coeffVec[3] << "\n";
		cout << "E = " << coeffVec[4] << "\t\tF = " << coeffVec[5] << "\n";
		cout << "G = " << coeffVec[6] << "\t\tH = " << coeffVec[7] << "\n";
		cout << "I = " << coeffVec[8] << "\t\tJ = " << coeffVec[9] << "\n";
		cout << "K = " << coeffVec[10] << "\t\tL = " << coeffVec[11] << "\n";
		cout << "M = " << coeffVec[12] << "\t\tN = " << coeffVec[13] << "\n";
		cout << "O = " << coeffVec[14] << "\t\tP = " << coeffVec[15] << "\n";
		cout << "Q = " << coeffVec[16] << "\t\tR = " << coeffVec[17] << "\n";
		cout << "S = " << coeffVec[18] << "\t\tT = " << coeffVec[19] << "\n";
*/
		out << "\nThe values of x0 follow:\n";
		for (i = 0; i < rDim; ++i) out << "x0[" << i << "] = " << x0[i] << "\n";

		uncmnd_ak1(rDim, coeffVec, xVec, &info, &fv, x0, p, g, gtemp, HessMat, YVEC, UVEC, WVEC);

		out << "\ninfo = " << info << " \n";

		switch (info)
		{
		//case 0:	cout << "Optimal solution found.\n";
		//		break;
		case 1:	out << "Terminated with gradient small, x is probably optimal.\n";
				break;
		case 2:	out << "Terminated with stepsize small, x is probably optimal.\n";
				break;
		case 3:	out << "Lower point cannot be found, x is probably optimal.\n";
			break;
		case 4:	out << " Iteration limit (150) exceeded.\n";
				break;
		case 5:	out << "Too many large steps, function may be unbounded.\n";
				break;
		default:	out << "Invalid Error Code.\n";
		} //End switch

		out << "\nThe values of x follow:" << "  \n";
		for (i = 0; i < rDim; ++i) out << xVec[i] << "  \n";

		out << "\nAt this point, the Objective Function = " << fv << " \n";

	} //End if rflag = 'Y'
	else cout << "\nNot ready. Try again when ready with information. \n";

	cout << "\nProgram ran to completion. About to exit.\n";
	cout << "\nEnter any key to continue. \n";
	cin >> rflag;
	return 0;
}                      // End main program.