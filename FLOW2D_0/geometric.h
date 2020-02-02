#pragma once

#include "shapes.h"
#include "exact_arithmetic.h"
#include "enumerators.h"


template<typename T> T sqr(T x) {

	return x * x;

};

/*****************************************************************************/
/*                                                                           */
/*    - Constants and thresholds                                             */
/*                                                                           */
/*                                                                           */
/*    EPSILON_EQUALITY		  : measures distance between two points         */
/*                                                                           */
/*                                                                           */
/*    - These are use in Fast geometric predicates                           */
/*                                                                           */
/*    EPSILON_ORIENTATION	  : if less than this EPS, orientation_fast		 */
/*							    returns '0.0'								 */
/*    EPSILON_IN_CIRCUMCIRCLE : if bigger than this EPS, in_circle_fast		 */
/*							    returns 'true'								 */
/*    EPSILON_ON_EDGE		  : if bigger than this (EPS x length of edge)   */
/*								on edge_fast  returns 'false'				 */
/*                            : Should be bigger than EPSILON_ORIENTATION    */
/*                                                                           */
/*                                                                           */
/*****************************************************************************/

double const EPSILON_EQUALITY =  1.e-5;
double const EPSILON_ORIENTATION = 1.e-5;
double const EPSILON_IN_CIRCUMCIRCLE = 1.e-5;
double const EPSILON_ON_EDGE = 1.e-3;		// Should be bigger than ORIENTATION -> Otherwise 'COLLINEAR' triangle could be formed => In circum-circle test can be non-symmetric => program failure
double const EPSILON_SEGMENT_ENCROACH = 1.e-4;

double const Pi = 3.14159265358979323846;




bool intersecting_edges(Edge * e, Edge * f);
bool intersecting_edges(Vertex * p1, Vertex * p2, Vertex * q1, Vertex * q2);

bool get_edge_intersection(Vertex * p1, Vertex * p2, Vertex * p3, Vertex * p4, double & x, double & y);


bool in_triangle_fast(Triangle * t, Vertex * v);
bool in_triangle_robust(Triangle * t, Vertex * v);

bool in_circle_fast(Triangle * t, Vertex * v);
bool in_circle_robust(Triangle * t, Vertex * v);

bool on_edge_fast(Vertex * a, Vertex * b, Vertex * v);
bool on_edge_fast(Edge * e, Vertex * v);

bool on_edge_robust(Vertex * a, Vertex * b, Vertex * v);
bool on_edge_robust(Edge * e, Vertex * v);

double orientation_fast(Vertex * a, Vertex * b, Vertex * v);
double orientation_robust(Vertex * a, Vertex * b, Vertex * v);


double orientationAdapt(Vertex * a, Vertex * b, Vertex * v, const double detsum);
bool in_circleAdapt(Triangle * t, Vertex * v, const double permanent);



bool is_encroached(Vertex * const a, Vertex * const b, Vertex * const v);
void get_mid_point(Edge * e, double &x_mid, double &y_mid);





enum PREDICATES { FAST, ROBUST };

template<GEOMETRIC_KERNEL GK>
class Predicates {


public:


	Predicates() {};
	~Predicates() {};


	PREDICATES const predicate = (GK == GEOMETRIC_KERNEL::INEXACT_PREDICATES_INEXACT_CONSTRUCT || GK == GEOMETRIC_KERNEL::INEXACT_PREDICATES_EXACT_CONSTRUCT) ? FAST : ROBUST;


	inline bool in_triangle(Triangle * const t, Vertex * const v) const {


		switch (predicate) {

		case FAST:
			return in_triangle_fast(t, v);
			break;
		case ROBUST:
			return in_triangle_robust(t, v);
			break;

		}

		return 0;

	};
	inline bool in_circle(Triangle * const t, Vertex * const v) const {


		switch (predicate) {

		case FAST:
			return in_circle_fast(t, v);
			break;
		case ROBUST:
			return in_circle_robust(t, v);
			break;

		}

		return 0;

	};
	inline bool on_edge(Vertex * const a, Vertex * const b, Vertex * const v) const {


		switch (predicate) {

		case FAST:
			return on_edge_fast(a, b, v);
			break;
		case ROBUST:
			return on_edge_robust(a, b, v);
			break;

		}

		return 0;

	};
	inline bool on_edge(Edge * const e, Vertex * const v) const {


		switch (predicate) {

		case FAST:
			return on_edge_fast(e, v);
			break;
		case ROBUST:
			return on_edge_robust(e, v);
			break;

		}

		return 0;

	};

	double orientation(Vertex * const a, Vertex * const b, Vertex * const c) const {


		switch (predicate) {

		case FAST:
			return orientation_fast(a, b, c);
			break;
		case ROBUST:
			return orientation_robust(a, b, c);
			break;

		}

		return 0;

	};



	bool intersecting_edges(Edge * const e, Edge * const f) const {


		Vertex * const p1 = e->a;
		Vertex * const p2 = e->b;

		Vertex * const q1 = f->a;
		Vertex * const q2 = f->b;

		const double s1 = (q1->x - p1->x)*(p2->y - p1->y) - (q1->y - p1->y)*(p2->x - p1->x);
		const double s2 = (q2->x - p1->x)*(p2->y - p1->y) - (q2->y - p1->y)*(p2->x - p1->x);

		const double s3 = (p1->x - q1->x)*(q2->y - q1->y) - (p1->y - q1->y)*(q2->x - q1->x);
		const double s4 = (p2->x - q1->x)*(q2->y - q1->y) - (p2->y - q1->y)*(q2->x - q1->x);

		const bool b1 = s1 * s2 <= 0.0;
		const bool b2 = s3 * s4 <= 0.0;

		return b1 && b2;

	};
	bool intersecting_edges(Vertex * const p1, Vertex * const p2, Vertex * const q1, Vertex * const q2) const {


		const double s1 = (q1->x - p1->x)*(p2->y - p1->y) - (q1->y - p1->y)*(p2->x - p1->x);
		const double s2 = (q2->x - p1->x)*(p2->y - p1->y) - (q2->y - p1->y)*(p2->x - p1->x);

		const double s3 = (p1->x - q1->x)*(q2->y - q1->y) - (p1->y - q1->y)*(q2->x - q1->x);
		const double s4 = (p2->x - q1->x)*(q2->y - q1->y) - (p2->y - q1->y)*(q2->x - q1->x);

		const bool b1 = s1 * s2 <= 0.0;
		const bool b2 = s3 * s4 <= 0.0;

		return b1 && b2;

	};


//private:


	bool in_triangle_fast(Triangle * const t, Vertex * const v) const {

		const bool cs1 = orientation_fast(t->vertices[0], t->vertices[1], v) >= 0.0;
		const bool cs2 = orientation_fast(t->vertices[1], t->vertices[2], v) >= 0.0;
		const bool cs3 = orientation_fast(t->vertices[2], t->vertices[0], v) >= 0.0;

		return cs1 && cs2 && cs3;

	};
	bool in_triangle_robust(Triangle * const t, Vertex * const v) const {

		const bool cs1 = orientation_robust(t->vertices[0], t->vertices[1], v) >= 0.0;
		const bool cs2 = orientation_robust(t->vertices[1], t->vertices[2], v) >= 0.0;
		const bool cs3 = orientation_robust(t->vertices[2], t->vertices[0], v) >= 0.0;

		return cs1 && cs2 && cs3;

	};

	bool in_circle_fast(Triangle * const t, Vertex * const v) const {


		const double px = v->x;
		const double py = v->y;

		const double avx = t->vertices[0]->x - px;
		const double avy = t->vertices[0]->y - py;
		const double bvx = t->vertices[1]->x - px;
		const double bvy = t->vertices[1]->y - py;
		const double cvx = t->vertices[2]->x - px;
		const double cvy = t->vertices[2]->y - py;

		const double abdet = avx * bvy - bvx * avy;
		const double bcdet = bvx * cvy - cvx * bvy;
		const double cadet = cvx * avy - avx * cvy;

		const double alift = avx * avx + avy * avy;
		const double blift = bvx * bvx + bvy * bvy;
		const double clift = cvx * cvx + cvy * cvy;

		return alift * bcdet + blift * cadet + clift * abdet >= EPSILON_IN_CIRCUMCIRCLE;

	};
	bool in_circle_robust(Triangle * const t, Vertex * const v) const {


		const double vx = v->x;
		const double vy = v->y;

		const double avx = t->vertices[0]->x - vx;
		const double bvx = t->vertices[1]->x - vx;
		const double cvx = t->vertices[2]->x - vx;

		const double avy = t->vertices[0]->y - vy;
		const double bvy = t->vertices[1]->y - vy;
		const double cvy = t->vertices[2]->y - vy;

		const double alift = avx * avx + avy * avy;
		const double blift = bvx * bvx + bvy * bvy;
		const double clift = cvx * cvx + cvy * cvy;

		const double bvx_cvy = bvx * cvy;
		const double cvx_bvy = cvx * bvy;

		const double cvx_avy = cvx * avy;
		const double avx_cvy = avx * cvy;

		const double avx_bvy = avx * bvy;
		const double bvx_avy = bvx * avy;


		const double det = alift * (bvx_cvy - cvx_bvy)
			+ blift * (cvx_avy - avx_cvy)
			+ clift * (avx_bvy - bvx_avy);

		const double permanent = (abs(bvx_cvy) + abs(cvx_bvy)) * alift
			+ (abs(cvx_avy) + abs(avx_cvy)) * blift
			+ (abs(avx_bvy) + abs(bvx_avy)) * clift;

		const double  errbound = incircleBoundA * permanent;


		if (abs(det) > errbound)
			return det > 0.0 ? true : false;

		// If the points are cocircular, we do not flip edge
		if (det == 0.0)
			return false;

		return in_circleAdapt(t, v, permanent);

	};

	bool on_edge_fast(Vertex * const a, Vertex * const b, Vertex * const v) const {


		const double ax = a->x;
		const double ay = a->y;

		const double bx = b->x;
		const double by = b->y;

		const double vx = v->x;
		const double vy = v->y;

		const double lenght_ab = sqrt(sqr(bx - ax) + sqr(by - ay));

		// Points are not collinear. Error is normalized to edge length
		if (abs(orientation_fast(a, b, v)) > EPSILON_ON_EDGE*lenght_ab)
			return false;

		const double dot = (vx - ax)*(bx - ax) + (vy - ay)*(by - ay);

		// point is under / above the segment
		if (dot < 0.0)
			return false;

		const double lenght_squared_ab = sqr(lenght_ab);

		// point is above / under the segment
		if (dot > lenght_squared_ab)
			return false;

		return true;

	};
	bool on_edge_fast(Edge * const e, Vertex * const v) const {


		Vertex * a = e->a;
		Vertex * b = e->b;

		const double ax = a->x;
		const double ay = a->y;

		const double bx = b->x;
		const double by = b->y;

		const double vx = v->x;
		const double vy = v->y;

		const double lenght_ab = sqrt(sqr(bx - ax) + sqr(by - ay));

		if (abs(orientation_fast(a, b, v)) > EPSILON_ON_EDGE*lenght_ab)
			return false;

		const double dot = (vx - ax)*(bx - ax) + (vy - ay)*(by - ay);

		if (dot < 0.0)
			return false;

		const double lenght_squared_ab = sqr(lenght_ab);

		if (dot > lenght_squared_ab)
			return false;

		return true;

	};

	bool on_edge_robust(Vertex * const a, Vertex * const b, Vertex * const v) const {


		const double ax = a->x;
		const double ay = a->y;

		const double bx = b->x;
		const double by = b->y;

		const double vx = v->x;
		const double vy = v->y;

		const double lenght_ab = sqrt(sqr(bx - ax) + sqr(by - ay));

		if (abs(orientation_robust(a, b, v)) > EPSILON_ON_EDGE*lenght_ab)
			return false;

		const double left = (vx - ax)*(bx - ax);
		const double right = (vy - ay)*(by - ay);

		if ((left < 0.0 && right <= 0.0) || (left <= 0.0 && right < 0.0))
			return false;


		const double dot = left + right;
		const double length_squared = sqr(lenght_ab);

		const double diff = dot - length_squared;
		const double diff_sum = dot + length_squared;

		const double errorBound = orientBoundA * diff_sum;

		// If the sign is trusted
		if (abs(diff) >= errorBound)
			return diff > 0.0 ? false : true;

		//return onEdgeAdapt(Edge * e, Vertex * v, diff_sum);

		if (diff == 0.0)
			return true;

		std::cout << "Implemenet 'onEdgeAdapat(Edge * e, Vertex * v, double diff_sum)'" << std::endl;
		system("pause");

		return true;

	};
	bool on_edge_robust(Edge * const e, Vertex * const v) const {

		Vertex * a = e->a;
		Vertex * b = e->b;

		const double ax = a->x;
		const double ay = a->y;

		const double bx = b->x;
		const double by = b->y;

		const double vx = v->x;
		const double vy = v->y;

		const double lenght_ab = sqrt(sqr(bx - ax) + sqr(by - ay));

		if (abs(orientation_robust(a, b, v)) > EPSILON_ON_EDGE*lenght_ab)
			return false;


		const double left = (vx - ax)*(bx - ax);
		const double right = (vy - ay)*(by - ay);

		if ((left < 0.0 && right <= 0.0) || (left <= 0.0 && right < 0.0))
			return false;


		const double dot = left + right;
		const double length_squared = sqr(lenght_ab);
		const double diff = dot - length_squared;
		const double diff_sum = dot + length_squared;


		const double errorBound = orientBoundA * diff_sum;

		// If the sign is trusted
		if (abs(diff) >= errorBound)
			return diff > 0.0 ? false : true;

		//return onEdgeAdapt(Edge * e, Vertex * v, diff_sum);


		std::cout << "Implemenet 'onEdgeAdapat(Edge * e, Vertex * v, double diff_sum)'" << std::endl;
		system("pause");

		return true;

	};

	double orientation_fast(Vertex * const a, Vertex * const b, Vertex * const c) const {

		const double det = (a->x - c->x) * (b->y - c->y) - (a->y - c->y) * (b->x - c->x);

		if (abs(det) <= EPSILON_ORIENTATION)
			return 0.0;

		return det;

	};
	double orientation_robust(Vertex * const a, Vertex * const b, Vertex * const c) const {


		const double left = (a->x - c->x) * (b->y - c->y);
		const double right = (a->y - c->y) * (b->x - c->x);

		const double det = left - right;

		double detsum = 0.0;

		if (left > 0.0) {

			if (right <= 0.0)
				return det;
			else
				detsum = left + right;

		}
		else if (left < 0.0) {

			if (right >= 0.0)
				return det;
			else
				detsum = -left - right;

		}
		else {

			return det;

		}

		const double errorBound = orientBoundA * detsum;

		if (abs(det) >= errorBound)
			return det;

		return orientationAdapt(a, b, c, detsum);

	};


	double orientationAdapt(Vertex * const a, Vertex * const b, Vertex * const v, const double detsum) const {


		double B[4];

		double av_x_tail;
		double av_y_tail;
		double bv_x_tail;
		double bv_y_tail;

		double left;
		double right;
		double left_tail;
		double right_tail;

		double av_x = a->x - v->x;
		double av_y = a->y - v->y;
		double bv_x = b->x - v->x;
		double bv_y = b->y - v->y;

		Two_Product(av_x, bv_y, left, left_tail);
		Two_Product(av_y, bv_x, right, right_tail);

		Two_Two_Diff(left, left_tail, right, right_tail, B[3], B[2], B[1], B[0]);


		double det = Estimate(4, B);
		double errorBound = orientBoundB * detsum;

		if (abs(det) >= errorBound)
			return det;


		Two_Diff_Tail(a->x, v->x, av_x, av_x_tail);
		Two_Diff_Tail(b->x, v->x, bv_x, bv_x_tail);
		Two_Diff_Tail(a->y, v->y, av_y, av_y_tail);
		Two_Diff_Tail(b->y, v->y, bv_y, bv_y_tail);

		if (av_x_tail == 0.0 && av_y_tail == 0.0 && bv_x_tail == 0.0 && bv_y_tail == 0.0)
			return det;


		errorBound = orientBoundC * detsum + resulterrbound * abs(det);
		det += (av_x * bv_y_tail + bv_y * av_x_tail) - (av_y * bv_x_tail + bv_x * av_y_tail);

		if (abs(det) >= errorBound)
			return det;



		// Exact evaluation
		double C1[8], C2[12], D[16];
		double u[4];

		double s1, t1;
		double s0, t0;


		Two_Product(av_x_tail, bv_y, s1, s0);
		Two_Product(av_y_tail, bv_x, t1, t0);

		Two_Two_Diff(s1, s0, t1, t0, u[3], u[2], u[1], u[0]);

		const int C1_length = Fast_Expansion_Sum_Zero_Eliminiation(4, B, 4, u, C1);


		Two_Product(av_x, bv_y_tail, s1, s0);
		Two_Product(av_y, bv_x_tail, t1, t0);

		Two_Two_Diff(s1, s0, t1, t0, u[3], u[2], u[1], u[0]);

		const int C2_length = Fast_Expansion_Sum_Zero_Eliminiation(C1_length, C1, 4, u, C2);


		Two_Product(av_x_tail, bv_y_tail, s1, s0);
		Two_Product(av_y_tail, bv_x_tail, t1, t0);

		Two_Two_Diff(s1, s0, t1, t0, u[3], u[2], u[1], u[0]);

		const int D_length = Fast_Expansion_Sum_Zero_Eliminiation(C2_length, C2, 4, u, D);


		return(D[D_length - 1]);

	};
	bool in_circleAdapt(Triangle * const t, Vertex * const v, const double permanent) const {


		double det;
		double errorBound;

		double adxtail, bdxtail, cdxtail, adytail, bdytail, cdytail;

		double bvx_cvy1, cvx_bvy1, cdxady1, adxcdy1, adxbdy1, bdxady1;
		double bvx_cvy0, cvx_bvy0, cdxady0, adxcdy0, adxbdy0, bdxady0;
		double bc[4], ca[4], ab[4];

		double axbc[8], axxbc[16], aybc[8], ayybc[16], adet[32];
		double bxca[8], bxxca[16], byca[8], byyca[16], bdet[32];
		double cxab[8], cxxab[16], cyab[8], cyyab[16], cdet[32];
		double abdet[64];
		double fin1[1152];

		int axbclen, axxbclen, aybclen, ayybclen, alen;
		int bxcalen, bxxcalen, bycalen, byycalen, blen;
		int cxablen, cxxablen, cyablen, cyyablen, clen;
		int ablen;
		int finlength;

		const double vx = v->x;
		const double vy = v->y;

		double avx = t->vertices[0]->x - vx;
		double bvx = t->vertices[1]->x - vx;
		double cvx = t->vertices[2]->x - vx;

		double avy = t->vertices[0]->y - vy;
		double bvy = t->vertices[1]->y - vy;
		double cvy = t->vertices[2]->y - vy;


		Two_Product(bvx, cvy, bvx_cvy1, bvx_cvy0);
		Two_Product(cvx, bvy, cvx_bvy1, cvx_bvy0);
		Two_Two_Diff(bvx_cvy1, bvx_cvy0, cvx_bvy1, cvx_bvy0, bc[3], bc[2], bc[1], bc[0]);

		axbclen = Scale_Expansion_Zero_Eliminiation(4, bc, avx, axbc);
		axxbclen = Scale_Expansion_Zero_Eliminiation(axbclen, axbc, avx, axxbc);
		aybclen = Scale_Expansion_Zero_Eliminiation(4, bc, avy, aybc);
		ayybclen = Scale_Expansion_Zero_Eliminiation(aybclen, aybc, avy, ayybc);

		alen = Fast_Expansion_Sum_Zero_Eliminiation(axxbclen, axxbc, ayybclen, ayybc, adet);

		Two_Product(cvx, avy, cdxady1, cdxady0);
		Two_Product(avx, cvy, adxcdy1, adxcdy0);
		Two_Two_Diff(cdxady1, cdxady0, adxcdy1, adxcdy0, ca[3], ca[2], ca[1], ca[0]);

		bxcalen = Scale_Expansion_Zero_Eliminiation(4, ca, bvx, bxca);
		bxxcalen = Scale_Expansion_Zero_Eliminiation(bxcalen, bxca, bvx, bxxca);
		bycalen = Scale_Expansion_Zero_Eliminiation(4, ca, bvy, byca);
		byycalen = Scale_Expansion_Zero_Eliminiation(bycalen, byca, bvy, byyca);

		blen = Fast_Expansion_Sum_Zero_Eliminiation(bxxcalen, bxxca, byycalen, byyca, bdet);

		Two_Product(avx, bvy, adxbdy1, adxbdy0);
		Two_Product(bvx, avy, bdxady1, bdxady0);
		Two_Two_Diff(adxbdy1, adxbdy0, bdxady1, bdxady0, ab[3], ab[2], ab[1], ab[0]);

		cxablen = Scale_Expansion_Zero_Eliminiation(4, ab, cvx, cxab);
		cxxablen = Scale_Expansion_Zero_Eliminiation(cxablen, cxab, cvx, cxxab);
		cyablen = Scale_Expansion_Zero_Eliminiation(4, ab, cvy, cyab);
		cyyablen = Scale_Expansion_Zero_Eliminiation(cyablen, cyab, cvy, cyyab);

		clen = Fast_Expansion_Sum_Zero_Eliminiation(cxxablen, cxxab, cyyablen, cyyab, cdet);

		ablen = Fast_Expansion_Sum_Zero_Eliminiation(alen, adet, blen, bdet, abdet);
		finlength = Fast_Expansion_Sum_Zero_Eliminiation(ablen, abdet, clen, cdet, fin1);



		det = Estimate(finlength, fin1);

		errorBound = incircleBoundB * permanent;

		if (abs(det) >= errorBound)
			return det > 0.0 ? true : false;


		Two_Diff_Tail(t->vertices[0]->x, vx, avx, adxtail);
		Two_Diff_Tail(t->vertices[0]->y, vy, avy, adytail);
		Two_Diff_Tail(t->vertices[1]->x, vx, bvx, bdxtail);
		Two_Diff_Tail(t->vertices[1]->y, vy, bvy, bdytail);
		Two_Diff_Tail(t->vertices[2]->x, vx, cvx, cdxtail);
		Two_Diff_Tail(t->vertices[2]->y, vy, cvy, cdytail);



		if (adxtail == 0.0 && bdxtail == 0.0 && cdxtail == 0.0 && adytail == 0.0 && bdytail == 0.0 && cdytail == 0.0)
			return det > 0.0 ? true : false;


		errorBound = incircleBoundC * permanent + resulterrbound * abs(det);

		det += ((avx * avx + avy * avy) * ((bvx * cdytail + cvy * bdxtail) - (bvy * cdxtail + cvx * bdytail))
			+ 2.0 * (avx * adxtail + avy * adytail) * (bvx * cvy - bvy * cvx))
			+ ((bvx * bvx + bvy * bvy) * ((cvx * adytail + avy * cdxtail) - (cvy * adxtail + avx * cdytail))
				+ 2.0 * (bvx * bdxtail + bvy * bdytail) * (cvx * avy - cvy * avx))
			+ ((cvx * cvx + cvy * cvy) * ((avx * bdytail + bvy * adxtail) - (avy * bdxtail + bvx * adytail))
				+ 2.0 * (cvx * cdxtail + cvy * cdytail) * (avx * bvy - avy * bvx));

		if (abs(det) >= errorBound)
			return det > 0.0 ? true : false;


		std::cout << "Implement the rest in 'in_circleAdapt(Triangle * t, Vertex * a, const double permanent)'" << std::endl;

		return false;

	};

};




/*

double barycentric0_of(Vertex * v) const {

	return orientation_robust(vertices[1], vertices[2], v) / orientation_robust(vertices[0], vertices[1], vertices[2]);

};
double barycentric1_of(Vertex * v) const {

	return orientation_robust(vertices[2], vertices[0], v) / orientation_robust(vertices[0], vertices[1], vertices[2]);

};
double barycentric2_of(Vertex * v) const {

	return orientation_robust(vertices[0], vertices[1], v) / orientation_robust(vertices[0], vertices[1], vertices[2]);

};

*/