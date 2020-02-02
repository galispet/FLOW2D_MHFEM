
#include "exact_arithmetic.h"

#include <cmath>

/*****************************************************************************/
/*                                                                           */
/*  estimate()   Produce a one-word estimate of an expansion's value.        */
/*                                                                           */
/*  See either version of my paper for details.                              */
/*                                                                           */
/*****************************************************************************/
double Estimate(const int m, double * e) {

	double Q = e[0];

	for (int i = 1; i < m; i++)
		Q += e[i];

	return Q;

};


/*****************************************************************************/
/*                                                                           */
/*  fast_expansion_sum_zeroelim()   Sum two expansions, eliminating zero     */
/*                                  components from the output expansion.    */
/*                                                                           */
/*  Sets h = e + f.  See the long version of my paper for details.           */
/*                                                                           */
/*  If round-to-even is used (as with IEEE 754), maintains the strongly      */
/*  nonoverlapping property.  (That is, if e is strongly nonoverlapping, h   */
/*  will be also.)  Does NOT maintain the nonoverlapping or nonadjacent      */
/*  properties.                                                              */
/*                                                                           */
/*****************************************************************************/
int Fast_Expansion_Sum_Zero_Eliminiation(const int m, double * e, const int n, double * f, double *  h) {


	double Q = 0.0;
	double Qnew = 0.0;
	double hi = 0.0;

	int eindex = 0;
	int findex = 0;
	int hindex = 0;

	double enow = e[0];
	double fnow = f[0];

	if (abs(fnow) > abs(enow)) {

		Q = enow;
		enow = e[++eindex];

	}
	else {

		Q = fnow;
		fnow = f[++findex];

	}

	if (eindex < m && findex < n) {


		if (fnow > abs(enow)) {

			Fast_Two_Sum(enow, Q, Qnew, hi);
			enow = e[++eindex];

		}
		else {

			Fast_Two_Sum(fnow, Q, Qnew, hi);
			fnow = f[++findex];

		}

		Q = Qnew;

		if (hi != 0.0)
			h[hindex++] = hi;


		while (eindex < m && findex < n) {

			if (abs(fnow) > abs(enow)) {

				Two_Sum(Q, enow, Qnew, hi);
				enow = e[++eindex];

			}
			else {

				Two_Sum(Q, fnow, Qnew, hi);
				fnow = f[++findex];

			}

			Q = Qnew;

			if (hi != 0.0)
				h[hindex++] = hi;

		}
	}

	while (eindex < m) {

		Two_Sum(Q, enow, Qnew, hi);
		enow = e[++eindex];
		Q = Qnew;

		if (hi != 0.0)
			h[hindex++] = hi;

	}
	while (findex < n) {

		Two_Sum(Q, fnow, Qnew, hi);
		fnow = f[++findex];
		Q = Qnew;

		if (hi != 0.0)
			h[hindex++] = hi;

	}

	if (Q != 0.0 || hindex == 0)
		h[hindex++] = Q;

	return hindex;

};


/***************************************************************************** /
/*                                                                           */
/*  scale_expansion_zeroelim()   Multiply an expansion by a scalar,          */
/*                               eliminating zero components from the        */
/*                               output expansion.                           */
/*                                                                           */
/*  Sets h = be.											                 */
/*                                                                           */
/*****************************************************************************/
int Scale_Expansion_Zero_Eliminiation(const int m, double * e,  double b, double *  h) {


	double Q;
	double sum;
	double hh;
	double product1;
	double product0;
	int hindex;
	double enow;
	double bhi, blo;

	Split(b, bhi, blo);
	Two_Product_Presplit(e[0], b, bhi, blo, Q, hh);

	hindex = 0;

	if (hh != 0)
		h[hindex++] = hh;

	
	for (int i = 1; i < m; i++) {

		enow = e[i];

		Two_Product_Presplit(enow, b, bhi, blo, product1, product0);
		Two_Sum(Q, product0, sum, hh);


		if (hh != 0)
			h[hindex++] = hh;
		
		Fast_Two_Sum(product1, sum, Q, hh);

		if (hh != 0)
			h[hindex++] = hh;
		
	}


	if (Q != 0.0 || hindex == 0)
		h[hindex++] = Q;

	return hindex;

};