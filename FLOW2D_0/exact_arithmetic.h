#pragma once

#include <float.h>


/*****************************************************************************/
/*                                                                           */
/*  Routines for Arbitrary Precision Floating-point Arithmetic               */
/*  and Fast Robust Geometric Predicates                                     */
/*****************************************************************************/

/*****************************************************************************/
/*                                                                           */	
/*  Be sure to call exactinit() once, before calling any of the arithmetic   */
/*    functions or geometric predicates.  Also be sure to turn on the        */
/*    optimizer when compiling this file.                                    */
/*                                                                           */
/*                                                                           */
/*                                                                           */
/*                                                                           */
/*  An expansion is represented by an array of floating-point numbers,       */
/*    sorted from smallest to largest magnitude (possibly with interspersed  */
/*    zeros).  The length of each expansion is stored as a separate integer, */
/*    and each arithmetic function returns an integer which is the length    */
/*    of the expansion it created.                                           */
/*                                                                           */
/*                                                                           */
/*                                                                           */
/*  Several arithmetic functions are defined.  Their parameters are          */
/*                                                                           */
/*    e, f           Input expansions                                        */
/*    elen, flen     Lengths of input expansions (must be >= 1)              */
/*    h              Output expansion                                        */
/*    b              Input scalar                                            */
/*                                                                           */
/*                                                                           */
/*                                                                           */
/*  The arithmetic functions are                                             */
/*                                                                           */
/*    grow_expansion(elen, e, b, h)                                          */
/*    grow_expansion_zeroelim(elen, e, b, h)                                 */
/*    expansion_sum(elen, e, flen, f, h)                                     */
/*    expansion_sum_zeroelim1(elen, e, flen, f, h)                           */
/*    expansion_sum_zeroelim2(elen, e, flen, f, h)                           */
/*    fast_expansion_sum(elen, e, flen, f, h)                                */
/*    fast_expansion_sum_zeroelim(elen, e, flen, f, h)                       */
/*    linear_expansion_sum(elen, e, flen, f, h)                              */
/*    linear_expansion_sum_zeroelim(elen, e, flen, f, h)                     */
/*    scale_expansion(elen, e, b, h)                                         */
/*    scale_expansion_zeroelim(elen, e, b, h)                                */
/*    compress(elen, e, h)                                                   */
/*                                                                           */
/*                                                                           */
/*                                                                           */
/*                                                                           */
/*  All of these are described in the long version of the paper; some are    */
/*    described in the short version.  All return an integer that is the     */
/*    length of h.  Those with suffix _zeroelim perform zero elimination,    */
/*    and are recommended over their counterparts.  The procedure            */
/*    fast_expansion_sum_zeroelim() (or linear_expansion_sum_zeroelim() on   */
/*    processors that do not use the round-to-even tiebreaking rule) is      */
/*    recommended over expansion_sum_zeroelim().  Each procedure has a       */
/*    little note next to it (in the code below) that tells you whether or   */
/*    not the output expansion may be the same array as one of the input     */
/*    expansions.                                                            */
/*                                                                           */
/*                                                                           */
/*  If you look around below, you'll also find macros for a bunch of         */
/*    simple unrolled arithmetic operations, and procedures for printing     */
/*    expansions (commented out because they don't work with all C           */
/*    compilers) and for generating random floating-point numbers whose      */
/*    significand bits are all random.  Most of the macros have undocumented */
/*    requirements that certain of their parameters should not be the same   */
/*    variable; for safety, better to make sure all the parameters are       */
/*    distinct variables.  Feel free to send email to jrs@cs.cmu.edu if you  */
/*    have questions.                                                        */
/*                                                                           */
/*****************************************************************************/


/*****************************************************************************/
/*                                                                           */
/*                                                                           */
/* Many of the operations are broken up into two pieces, a main part that    */
/*   performs an approximate operation, and a "tail" that computes the       */
/*   roundoff error of that operation.                                       */
/*                                                                           */
/*                                                                           */
/* The operations :															 */
/*																			 */
/* Fast_Two_Sum(),Fast_Two_Diff(),Two_Sum(),Two_Diff(),Split(),Two_Product() */
/*																			 */
/* Each of these macros requires certain variables to be					 */
/*   defined in the calling routine.  The variables `bvirt', `c', `abig',    */
/*   `_i', `_j', `_k', `_l', `_m', and `_n' are declared `INEXACT' because   */
/*   they store the result of an operation that may incur roundoff error.    */
/*   The input parameter `x' (or the highest numbered `x_' parameter) must   */
/*   also be declared `INEXACT'.                                             */
/*                                                                           */
/*                                                                           */
/*****************************************************************************/



/*	splitter = ceil(p/2) , p is the mantisa of the floating-point number : double -> p = 52 + 1 */
const double SPLITTER = 134217729.0;

/*	Thresholds for adaptive 2D orientation test and 2D in circumcircle test						*/
const double resulterrbound = (3.0 + 8.0 * DBL_EPSILON) * DBL_EPSILON;

const double orientBoundA = (3.0 + 16.0 * DBL_EPSILON) * DBL_EPSILON;
const double orientBoundB = (2.0 + 12.0 * DBL_EPSILON) * DBL_EPSILON;
const double orientBoundC = (9.0 + 64.0 * DBL_EPSILON) * DBL_EPSILON * DBL_EPSILON;

const double incircleBoundA = (10.0 + 96.0 * DBL_EPSILON) * DBL_EPSILON;
const double incircleBoundB = (4.0 + 48.0 * DBL_EPSILON) * DBL_EPSILON;
const double incircleBoundC = (44.0 + 576.0 * DBL_EPSILON) * DBL_EPSILON * DBL_EPSILON;



double Estimate(const int m, double * e);
int Fast_Expansion_Sum_Zero_Eliminiation(const int m, double * e, const int n, double * f, double *  h);
int Scale_Expansion_Zero_Eliminiation(const int m, double * e, double b, double * h);





inline void Fast_Two_Sum_Tail(double a, double b, double &x, double &y) {

	const double b_virtual = x - a;
	y = b - b_virtual;

};
inline void Fast_Two_Sum(double a, double b, double &x, double &y) {

	x = a + b;
	Fast_Two_Sum_Tail(a, b, x, y);

};

inline void Fast_Two_Diff_Tail(double a, double b, double &x, double &y) {

	const double b_virtual = a - x;
	y = b_virtual - b;

};
inline void Fast_Two_Diff(double a, double b, double &x, double &y) {

	x = a - b;
	Fast_Two_Diff_Tail(a, b, x, y);

};

inline void Two_Sum_Tail(double a, double b, double &x, double &y) {

	const double b_virtual = x - a;
	const double a_virtual = x - b_virtual;

	const double b_roundoff = b - b_virtual;
	const double a_roundoff = a - a_virtual;

	y = a_roundoff + b_roundoff;

};
inline void Two_Sum(double a, double b, double &x, double &y) {

	x = a + b;
	Two_Sum_Tail(a, b, x, y);

};

inline void Two_Diff_Tail(double a, double b, double &x, double &y) {

	const double b_virtual = a - x;
	const double a_virtual = x + b_virtual;

	const double b_roundoff = b_virtual - b;
	const double a_roundoff = a - a_virtual;

	y = a_roundoff + b_roundoff;

};
inline void Two_Diff(double a, double b, double &x, double &y) {

	x = a - b;
	Two_Diff_Tail(a, b, x, y);

};

inline void Split(double a, double &a_hi, double &a_lo) {

	const double c = SPLITTER * a;
	const double a_big = c - a;

	a_hi = c - a_big;
	a_lo = c - a_hi;

};
inline void Two_Product_Tail(double a, double b, double &x, double &y) {

	double a_hi;
	double a_lo;

	double b_hi;
	double b_lo;

	Split(a, a_hi, a_lo);
	Split(b, b_hi, b_lo); 
	
	const double err1 = x - (a_hi * b_hi);
	const double err2 = err1 - (a_lo * b_hi);
	const double err3 = err2 - (a_hi * b_lo);

	y = a_lo * b_lo - err3;

};
inline void Two_Product(double a, double b, double &x, double &y) {

	x = a * b;
	Two_Product_Tail(a, b, x, y);

};

inline void Two_One_Diff(double a1, double a0, double b, double &x2, double &x1, double &x0) {

	double _i;

	Two_Diff(a0, b, _i, x0);
	Two_Sum(a1, _i, x2, x1);

};
inline void Two_Two_Diff(double a1, double a0, double b1, double b0, double &x3, double &x2, double &x1, double &x0) {

	double _j;
	double _0;

	Two_One_Diff(a1, a0, b0, _j, _0, x0);
	Two_One_Diff(_j, _0, b1, x3, x2, x1);

};

inline void Two_Product_Presplit(double a, double b, double bhi, double blo, double &x,double & y) {

	x = a * b;

	double ahi;
	double alo;

	Split(a, ahi, alo);

	const double err1 = x - (ahi * bhi);
	const double err2 = err1 - (alo * bhi);
	const double err3 = err2 - (ahi * blo);

	y = (alo * blo) - err3;

};


