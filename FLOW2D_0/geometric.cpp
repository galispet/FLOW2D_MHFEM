
#include "geometric.h"
#include "shapes.h"
#include "exact_arithmetic.h"

#include <cstdlib>
#include <iostream>



bool intersecting_edges(Edge * const e, Edge * const f) {


	Vertex * p1 = e->a;
	Vertex * p2 = e->b;

	Vertex * q1 = f->a;
	Vertex * q2 = f->b;

	const double s1 = (q1->x - p1->x)*(p2->y - p1->y) - (q1->y - p1->y)*(p2->x - p1->x);
	const double s2 = (q2->x - p1->x)*(p2->y - p1->y) - (q2->y - p1->y)*(p2->x - p1->x);

	const double s3 = (p1->x - q1->x)*(q2->y - q1->y) - (p1->y - q1->y)*(q2->x - q1->x);
	const double s4 = (p2->x - q1->x)*(q2->y - q1->y) - (p2->y - q1->y)*(q2->x - q1->x);

	const bool b1 = s1 * s2 <= 0.0;
	const bool b2 = s3 * s4 <= 0.0;

	return b1 && b2;

};
bool intersecting_edges(Vertex * const p1, Vertex * const p2, Vertex * const q1, Vertex * const q2) {


	const double s1 = (q1->x - p1->x)*(p2->y - p1->y) - (q1->y - p1->y)*(p2->x - p1->x);
	const double s2 = (q2->x - p1->x)*(p2->y - p1->y) - (q2->y - p1->y)*(p2->x - p1->x);

	const double s3 = (p1->x - q1->x)*(q2->y - q1->y) - (p1->y - q1->y)*(q2->x - q1->x);
	const double s4 = (p2->x - q1->x)*(q2->y - q1->y) - (p2->y - q1->y)*(q2->x - q1->x);

	const bool b1 = s1 * s2 <= 0.0;
	const bool b2 = s3 * s4 <= 0.0;


	//const double o1 = orientation_fast(p1, p2, q1);
	//const double o2 = orientation_fast(p1, p2, q2);

	//const double o3 = orientation_fast(q1, q2, p1);
	//const double o4 = orientation_fast(q1, q2, p2);

	//const bool h1 = o1 * o2 <= 0.0;
	//const bool h2 = o3 * o4 <= 0.0;
		

	return b1 && b2;

};

bool get_edge_intersection(Vertex * const p1, Vertex * const p2, Vertex * const p3, Vertex * const p4, double & x, double & y) {


	const double s1 = (p3->x - p1->x)*(p2->y - p1->y) - (p3->y - p1->y)*(p2->x - p1->x);
	const double s2 = (p4->x - p1->x)*(p2->y - p1->y) - (p4->y - p1->y)*(p2->x - p1->x);

	const double s3 = (p1->x - p3->x)*(p4->y - p3->y) - (p1->y - p3->y)*(p4->x - p3->x);
	const double s4 = (p2->x - p3->x)*(p4->y - p3->y) - (p2->y - p3->y)*(p4->x - p3->x);


	const double p12_x = p1->x - p2->x;
	const double p13_x = p1->x - p3->x;
	const double p14_x = p1->x - p4->x;
	const double p23_x = p2->x - p3->x;
	const double p34_x = p3->x - p4->x;

	const double p12_y = p1->y - p2->y;
	const double p13_y = p1->y - p3->y;
	const double p14_y = p1->y - p4->y;
	const double p23_y = p2->y - p3->y;
	const double p34_y = p3->y - p4->y;


	const double s11 = p13_x * p12_y;
	const double s12 = p13_y * p12_x;

	const double s21 = p14_x * p12_y;
	const double s22 = p14_y * p12_x;

	const double s31 = -p13_x * p34_y;
	const double s32 = -p13_y * p34_x;

	const double s41 = -p23_x * p34_y;
	const double s42 = -p23_y * p34_x;

	const bool b1 = (s11 - s12) * (s21 - s22) <= 0.0;
	const bool b2 = (s31 - s32) * (s41 - s42) <= 0.0;


	if (b1 && b2){


		double const denominator = -p34_x * p12_y + p12_x * p34_y;

		double const t1 = (p34_y*p13_x - p34_x * p13_y) / denominator;
		//double const t2 = (p12_y*p13_x - p12_x * p13_y) / denominator;


		// Collision detected
		x = p1->x - (t1 * p12_x);
		y = p1->y - (t1 * p12_y);

		return true;
	}

	return false; // No collision
}





bool in_triangle_fast(Triangle * const t, Vertex * const v) {

	const bool cs1 = orientation_fast(t->vertices[0], t->vertices[1], v) >= 0.0;
	const bool cs2 = orientation_fast(t->vertices[1], t->vertices[2], v) >= 0.0;
	const bool cs3 = orientation_fast(t->vertices[2], t->vertices[0], v) >= 0.0;

	return cs1 && cs2 && cs3;

};
bool in_triangle_robust(Triangle * const t, Vertex * const v) {

	const bool cs1 = orientation_robust(t->vertices[0], t->vertices[1], v) >= 0.0;
	const bool cs2 = orientation_robust(t->vertices[1], t->vertices[2], v) >= 0.0;
	const bool cs3 = orientation_robust(t->vertices[2], t->vertices[0], v) >= 0.0;

	return cs1 && cs2 && cs3;

};

bool in_circle_fast(Triangle * const t, Vertex * const v) {


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
bool in_circle_robust(Triangle * const t, Vertex * const v) {


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

bool on_edge_fast(Vertex * const a, Vertex * const b, Vertex * const v) {


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

	/*
	if (ax > bx && ay > by) {
		if (bx <= vx && vx <= ax && by <= vy && vy <= ay)
			return true;
	}
	else if (ax > bx && ay < by) {
		if (bx <= vx && vx <= ax && ay <= vy && vy <= by)
			return true;
	}
	else if (ax < bx && ay < by) {
		if (ax <= vx && vx <= bx && ay <= vy && vy <= by)
			return true;
	}
	else if (ax < bx && ay >by) {
		if (ax <= vx && vx <= bx && by <= vy && vy <= ay)
			return true;
	}

	return false;
	*/

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
bool on_edge_fast(Edge * const e, Vertex * const v) {


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

bool on_edge_robust(Vertex * const a, Vertex * const b, Vertex * const v) {


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
bool on_edge_robust(Edge * const e, Vertex * const v) {

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

double orientation_fast(Vertex * const a, Vertex * const b, Vertex * const v) {

	const double det = (a->x - v->x) * (b->y - v->y) - (a->y - v->y) * (b->x - v->x);

	if (abs(det) <= EPSILON_ORIENTATION)
		return 0.0;

	return det;

};
double orientation_robust(Vertex * const a, Vertex * const b, Vertex * const v) {


	const double left = (a->x - v->x) * (b->y - v->y);
	const double right = (a->y - v->y) * (b->x - v->x);

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

	return orientationAdapt(a, b, v, detsum);

};


double orientationAdapt(Vertex * const a, Vertex * const b, Vertex * const v, const double detsum) {


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

	if (abs(det) >= errorBound)	// PROBLEM OCCURS HERE !!! It is maybe computed wrong above
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
bool in_circleAdapt(Triangle * t, Vertex * v, const double permanent) {


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





bool is_encroached(Vertex * const a, Vertex * const b, Vertex * const v) {


	/*
	Vertex * const a = e->a;
	Vertex * const b = e->b;

	double const a_x = a->x;
	double const a_y = a->y;

	double const b_x = b->x;
	double const b_y = b->y;


	const double x_mid = 0.5*(a_x + b_x);
	const double y_mid = 0.5*(a_y + b_y);


	double const edge_radius = 0.5*sqrt(sqr(a_x - b_x) + sqr(a_y - b_y));
	double const vertex_from_edge_distance = sqrt(sqr(x_mid - v->x) + sqr(y_mid - v->y));

	return vertex_from_edge_distance < edge_radius;
	*/

	double const h1_x = v->x - a->x;
	double const h1_y = v->y - a->y;

	double const h2_x = v->x - b->x;
	double const h2_y = v->y - b->y;

	// If acute, then > 0 , else if obtuse, then < 0. If == 0 , then orthogonal. Encrochment happens with obtuse angle
	return h1_x * h2_x + h1_y * h2_y < -EPSILON_SEGMENT_ENCROACH;

};

void get_mid_point(Edge * e, double &x_mid, double &y_mid) {


	Vertex * const a = e->a;
	Vertex * const b = e->b;

	double const a_x = a->x;
	double const a_y = a->y;

	double const b_x = b->x;
	double const b_y = b->y;


	x_mid = 0.5*(a_x + b_x);
	y_mid = 0.5*(a_y + b_y);


};







