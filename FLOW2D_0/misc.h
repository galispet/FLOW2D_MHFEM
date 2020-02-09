#pragma once


#include "shapes.h"
#include "enumerators.h"
#include "geometric.h"
#include "coefficient_matrix.h"
#include "matrix.h"
#include "integration.h"

#include <iostream>



typedef double real;

double const TOL = DBL_EPSILON;

unsigned const MAX_IT = 30;
//unsigned const quadratureRule = 17;
//unsigned const error_quadRule = 17;

unsigned const quadrature_order = 6;


//bool const RT_ORDER = 0;

double const θ = 1.0;		// 1.0 = Backward Euler | 0.5 = Crank-Nicolson



typedef Vertex *	v_pointer;
typedef Edge *		e_pointer;
typedef Triangle *	t_pointer;



inline double sqr(double x) {

	return x * x;

};
inline double δij(unsigned const i, unsigned const j) {

	return i == j ? 1.0 : 0.0;

};

//double funnn(double x, double y) {
//	return pow(x + y,-0.5);
//};


inline void hard_copy(double *& copy_to, double *& copy_from, unsigned const n) {

	for (unsigned i = 0; i < n; i++)
		copy_to[i] = copy_from[i];

};



inline double barenblatt(double x, double y, double time) {

	double const norm_squared = x * x + y * y;

	//return (1.0 / pow(time, (1.0 / 3.0)))*fmax(1.0 - norm_squared / (12.0 * pow(time, (2.0 / 3.0))), 0.0);

	return (1.0 / sqrt(time))*fmax(1.0 - norm_squared / (16.0 * sqrt(time)), 0.0);

};

inline double barenblatt_dx(double x, double y, double time) {

	//if (abs(x) >= 2.0 *sqrt(3.0) * pow(t, 1.0 / 3.0))
	//	return 0.0;

	return x / (-6.0*time);

};
inline double barenblatt_dy(double x, double y, double time) {

	//if (abs(x) >= 2.0 *sqrt(3.0) * pow(t, 1.0 / 3.0))
	//	return 0.0;

	return y / (-6.0*time);

};


double integrate_edge(e_pointer const e, real time, double(*fun)(double, double, double)) {


	real const x0 = e->a->x;
	real const y0 = e->a->y;

	real const x1 = e->b->x;
	real const y1 = e->b->y;

	real const length = e->length();



	// Quadrature weights and points on reference segment [-1,1]
	gauss_quadrature_1D quad(quadrature_order);
	unsigned const num_quad_points = quad.NumberOfPoints;


	real integral = 0.0;

	for (unsigned n = 0; n < num_quad_points; n++) {


		// Quadrature point on reference segment [-1,1] and its weights
		real const ksi_ref		= quad.points[n];
		real const weight_ref	= quad.weights[n];

		real const x = x0 + 0.5*(1.0 + ksi_ref)*(x1 - x0);
		real const y = y0 + 0.5*(1.0 + ksi_ref)*(y1 - y0);

		integral += weight_ref * fun(x, y, time);

	}

	return 0.5 * length * integral;

};
double integrate_edge(e_pointer const e, double(*fun)(double, double)) {


	real const x0 = e->a->x;
	real const y0 = e->a->y;

	real const x1 = e->b->x;
	real const y1 = e->b->y;

	real const length = e->length();



	// Quadrature weights and points on reference segment [-1,1]
	gauss_quadrature_1D quad(quadrature_order);
	unsigned const num_quad_points = quad.NumberOfPoints;


	real integral = 0.0;

	for (unsigned n = 0; n < num_quad_points; n++) {


		// Quadrature point on reference segment [-1,1] and its weights
		real const ksi_ref		= quad.points[n];
		real const weight_ref	= quad.weights[n];

		real const x = x0 + 0.5*(1.0 + ksi_ref)*(x1 - x0);
		real const y = y0 + 0.5*(1.0 + ksi_ref)*(y1 - y0);

		integral += weight_ref * fun(x, y);

	}


	return 0.5 * length * integral;

};
double integrate_triangle(t_pointer const K, real time, double(*fun)(double, double, double)) {



	v_pointer const a = K->vertices[0];
	v_pointer const b = K->vertices[1];
	v_pointer const c = K->vertices[2];

	real const x0 = a->x;
	real const y0 = a->y;

	real const x1 = b->x;
	real const y1 = b->y;

	real const x2 = c->x;
	real const y2 = c->y;

	real const detJF = abs((x1 - x0)*(y2 - y0) - (x2 - x0)*(y1 - y0));


	// Quadrature weights and points on [-1,1]
	quadrature_triangle quad(quadrature_order);
	unsigned const num_quad_points = quad.NumberOfPoints;


	real integral = 0.0;

	for (unsigned n = 0; n < num_quad_points; n++) {


		real const s = quad.points_x[n];
		real const t = quad.points_y[n];

		real const w = quad.weights[n];


		// Gauss points on given triangle. This is the transformation : [-1,1] x [-1,1] square --> reference triangle ([0,0],[1,0],[0,1]) --> given triangle ([x0y,0],[x1,y1],[x2,y2]) 
		real const x = x0 + (x1 - x0)*s + (x2 - x0)*t;
		real const y = y0 + (y1 - y0)*s + (y2 - y0)*t;

		// Gauss quadrature
		integral += w * fun(x, y, time);


	}

	// Not quite sure, why there is the 1/2. See https://people.sc.fsu.edu/~jburkardt/datasets/quadrature_rules_tri/quadrature_rules_tri.html
	return 0.5 * detJF * integral;

};
double integrate_triangle(t_pointer const K, double(*fun)(double, double)) {


	v_pointer const a = K->vertices[0];
	v_pointer const b = K->vertices[1];
	v_pointer const c = K->vertices[2];

	real const x0 = a->x;
	real const y0 = a->y;

	real const x1 = b->x;
	real const y1 = b->y;

	real const x2 = c->x;
	real const y2 = c->y;

	real const detJF = abs((x1 - x0)*(y2 - y0) - (x2 - x0)*(y1 - y0));


	quadrature_triangle quad(quadrature_order);
	unsigned const num_quad_points = quad.NumberOfPoints;


	real integral = 0.0;

	for (unsigned n = 0; n < num_quad_points; n++) {


		real const s = quad.points_x[n];
		real const t = quad.points_y[n];

		real const w = quad.weights[n];

		// Gauss points on given triangle. This is the transformation : [-1,1] x [-1,1] square --> reference triangle ([0,0],[1,0],[0,1]) --> given triangle ([x0y,0],[x1,y1],[x2,y2]) 
		real const x = x0 + (x1 - x0)*s + (x2 - x0)*t;
		real const y = y0 + (y1 - y0)*s + (y2 - y0)*t;

		// Gauss quadrature
		integral += w * fun(x, y);

	}

	/*
	for (unsigned i = 0; i < num_quad_points; i++) {


		double const wi = quad.weights[i];
		double const ξi = quad.points_x[i];

		for (unsigned j = 0; j < num_quad_points; j++) {


			double const wj = quad.weights[j];
			double const ηj = quad.points_y[j];

			// Product weight
			double const w = (1.0 - ξi) / 8.0 * wi * wj;

			// Gauss points on reference triangle
			double const u = 0.5*(1.0 + ξi);
			double const v = 0.25*(1.0 - ξi)*(1.0 + ηj);

			// Gauss points on given triangle. This is the transformation : [-1,1] x [-1,1] square --> reference triangle ([0,0],[1,0],[0,1]) --> given triangle ([x0y,0],[x1,y1],[x2,y2])
			double const x = x0 + (x1 - x0)*u + (x2 - x0)*v;
			double const y = y0 + (y1 - y0)*u + (y2 - y0)*v;

			// Gauss quadrature
			s += w * fun(x, y);

		}
	}
	*/

	// Not quite sure, why there is the 1/2. See https://people.sc.fsu.edu/~jburkardt/datasets/quadrature_rules_tri/quadrature_rules_tri.html
	return 0.5 * detJF * integral;

};


inline double permeability(unsigned i, unsigned j, double x, double y) {


	double const k00 = 1.0;
	double const k01 = 0.0;
	double const k10 = 0.0;
	double const k11 = 1.0;

	if (i == 0 && j == 0) return k00;
	if (i == 0 && j == 1) return k01;
	if (i == 1 && j == 0) return k10;
	if (i == 1 && j == 1) return k11;

	return INFINITY;

};
inline void permeability(double x, double y, double(&out)[2][2]) {

	out[0][0] = 1.0;
	out[0][1] = 0.0;
	out[1][0] = 0.0;
	out[1][1] = 1.0;

};
inline double porosity(double x, double y) {

	return 1.0;

};
inline double viscosity(double x, double y) {

	return 0.5;

};
inline double source(double x, double y, double t) {

	return 0.0;

};
inline double iK(unsigned i, unsigned j, double x, double y) {

	double const k00 = permeability(0, 0, x, y);
	double const k01 = permeability(0, 1, x, y);
	double const k10 = permeability(1, 0, x, y);
	double const k11 = permeability(1, 1, x, y);

	double const idet = 1.0 / (k00 * k11 - k10 * k01);

	if (i == 0 && j == 0) return idet * k11;
	if (i == 0 && j == 1) return -idet * k01;
	if (i == 1 && j == 0) return -idet * k10;
	if (i == 1 && j == 1) return idet * k00;

	return INFINITY;

};
inline void iK(double x, double y, double(&out)[2][2]) {


	double K[2][2];
	permeability(x, y, K);

	double const k00 = K[0][0];
	double const k01 = K[0][1];
	double const k10 = K[1][0];
	double const k11 = K[1][1];

	double const idet = 1.0 / (k00 * k11 - k10 * k01);

	out[0][0] = idet * k11;
	out[0][1] = -idet * k01;
	out[1][0] = -idet * k10;
	out[1][1] = idet * k00;

};
template<typename real>
inline void permeability(real const x, real const y, Matrix<real> & K) {

	K(0, 0) = (real) 1.0;
	K(0, 1) = (real) 0.0;
	K(1, 0) = (real) 0.0;
	K(1, 1) = (real) 1.0;

};
template<typename real>
inline void iK(real const x, real const y, Matrix<real> & K) {


	real K[2][2];
	permeability(x, y, K);

	real const k00 = K[0][0];
	real const k01 = K[0][1];
	real const k10 = K[1][0];
	real const k11 = K[1][1];

	real const idet = 1.0 / (k00 * k11 - k10 * k01);

	K(0, 0) = (real)idet * k11;
	K(0, 1) = (real)-idet * k01;
	K(1, 0) = (real)-idet * k10;
	K(1, 1) = (real)idet * k00;

};

inline double eos(real const c) {

	return c;

};
inline double Deos(real const c) {

	return 1.0;

};



//inline double DIRICHLET_GAMMA_P_pressure(double x, double y, double time) {
//
//	return barenblatt(x, y, time);
//
//	//return x;
//
//};
//inline double DIRICHLET_GAMMA_P_concentration(double x, double y, double time) {
//
//	return barenblatt(x, y, time);
//
//	//return x;
//
//};
//inline double NEUMANN_GAMMA_Q_velocity(double x, double y, double time) {
//
//	return permeability(0, 0, x, y) * barenblatt_dx(x, y, time) / viscosity(x, y);
//
//};
//inline double NEUMANN_GAMMA_Q_concentration(double x, double y, double time) {
//
//	return barenblatt(x, y, time);
//
//};
//
//





inline unsigned LI(t_pointer const & K, e_pointer const & E, unsigned const & DOF) {

	return K->get_edge_index(E) + 3 * DOF;

};



inline real NEUMANN_GAMMA_Q_velocity(e_pointer const E, real const time) {

	return 0.0;

};
inline real DIRICHLET_GAMMA_Q_concentration(e_pointer const E, real const time) {

	return 0.0 * integrate_edge(E, time, barenblatt) / E->length();

};
inline real DIRICHLET_GAMMA_P_concentration(e_pointer const E, real const time) {

	return integrate_edge(E, time, barenblatt) / E->length();

};
inline real DIRICHLET_GAMMA_P_pressure(e_pointer const E, real const time) {

	return integrate_edge(E, time, barenblatt) / E->length();

};


inline real DIRICHLET_GAMMA_Q_concentration(real const s, real const t, real const time) {

	return 0.0 * barenblatt(s, t, time);

};
inline real DIRICHLET_GAMMA_P_concentration(real const s, real const t, real const time) {

	return barenblatt(s, t, time);

};
inline real DIRICHLET_GAMMA_P_pressure(real const s, real const t, real const time) {

	return barenblatt(s, t, time);

};


inline double phi1(real const s, real const t) {

	return 1.0;

};
inline double phi2(real const s, real const t) {

	return -1.0 + 2.0*s;

};
inline double phi3(real const s, real const t) {

	return -1.0 + 2.0*t;

};


// Sources integrals
double F1(t_pointer K, real const time) {


	v_pointer const a = K->vertices[0];
	v_pointer const b = K->vertices[1];
	v_pointer const c = K->vertices[2];

	double const x0 = a->x;
	double const y0 = a->y;

	double const x1 = b->x;
	double const y1 = b->y;

	double const x2 = c->x;
	double const y2 = c->y;


	// Quadrature weights and points on [-1,1]
	quadrature_triangle quad(quadrature_order);

	unsigned const num_quad_points = quad.NumberOfPoints;

	double integral = 0.0;

	for (unsigned i = 0; i < num_quad_points; i++) {


		double const w = quad.weights[i];
		double const s = quad.points_x[i];
		double const t = quad.points_y[i];

		// Product weight
		double const weight = 0.5 * w;

		// Gauss points on given triangle. This is the transformation : [-1,1] x [-1,1] square --> reference triangle ([0,0],[1,0],[0,1]) --> given triangle ([x0y,0],[x1,y1],[x2,y2]) 
		double const x = x0 + (x1 - x0)*s + (x2 - x0)*t;
		double const y = y0 + (y1 - y0)*s + (y2 - y0)*t;

		// Gauss quadrature
		integral += w * source(x, y, time)*phi1(x, y);

	}

	return integral;

};
double F2(t_pointer K, real const time) {


	v_pointer const a = K->vertices[0];
	v_pointer const b = K->vertices[1];
	v_pointer const c = K->vertices[2];

	double const x0 = a->x;
	double const y0 = a->y;

	double const x1 = b->x;
	double const y1 = b->y;

	double const x2 = c->x;
	double const y2 = c->y;


	// Quadrature weights and points on [-1,1]
	quadrature_triangle quad(quadrature_order);

	unsigned const num_quad_points = quad.NumberOfPoints;

	double integral = 0.0;

	for (unsigned i = 0; i < num_quad_points; i++) {


		double const w = quad.weights[i];
		double const s = quad.points_x[i];
		double const t = quad.points_y[i];

		// Product weight
		double const weight = 0.5 * w;

		// Gauss points on given triangle. This is the transformation : [-1,1] x [-1,1] square --> reference triangle ([0,0],[1,0],[0,1]) --> given triangle ([x0y,0],[x1,y1],[x2,y2]) 
		double const x = x0 + (x1 - x0)*s + (x2 - x0)*t;
		double const y = y0 + (y1 - y0)*s + (y2 - y0)*t;

		// Gauss quadrature
		integral += w * source(x, y, time)*phi2(x, y);

	}

	return integral;

};
double F3(t_pointer K, real const time) {


	v_pointer const a = K->vertices[0];
	v_pointer const b = K->vertices[1];
	v_pointer const c = K->vertices[2];

	double const x0 = a->x;
	double const y0 = a->y;

	double const x1 = b->x;
	double const y1 = b->y;

	double const x2 = c->x;
	double const y2 = c->y;


	// Quadrature weights and points on [-1,1]
	quadrature_triangle quad(quadrature_order);

	unsigned const num_quad_points = quad.NumberOfPoints;

	double integral = 0.0;

	for (unsigned i = 0; i < num_quad_points; i++) {


		double const w = quad.weights[i];
		double const s = quad.points_x[i];
		double const t = quad.points_y[i];

		// Product weight
		double const weight = 0.5 * w;

		// Gauss points on given triangle. This is the transformation : [-1,1] x [-1,1] square --> reference triangle ([0,0],[1,0],[0,1]) --> given triangle ([x0y,0],[x1,y1],[x2,y2]) 
		double const x = x0 + (x1 - x0)*s + (x2 - x0)*t;
		double const y = y0 + (y1 - y0)*s + (y2 - y0)*t;

		// Gauss quadrature
		integral += w * source(x, y, time)*phi3(x, y);

	}

	return integral;

};


/*****************************************************************************/
/*                                                                           */
/*    - My containers												         */
/*                                                                           */
/*****************************************************************************/

inline void evaluate_raviartthomas_basis(real const s, real const t, Matrix<real> & out) {

	out(0, 0) = (-3.0*s + 4.0*s*t + 4.0*s * s);
	out(1, 0) = (-3.0*t + 4.0*s*t + 4.0*t * t);

	out(0, 1) = (-1.0 + 5.0*s - 4.0*s * s);
	out(1, 1) = (t - 4.0*s * t);

	out(0, 2) = (s - 4.0*s * t);
	out(1, 2) = (-1.0 + 5.0*t - 4.0*t * t);

	out(0, 3) = (s + 4.0*t * s - 4.0*s * s);
	out(1, 3) = (-t - 4.0*s * t + 4.0*t * t);

	out(0, 4) = (-3.0 + 7.0*s + 6.0*t - 8.0*t * s - 4.0*s * s);
	out(1, 4) = (5.0*t - 4.0*s * t - 8.0*t * t);

	out(0, 5) = (-5.0*s + 4.0*t * s + 8.0*s * s);
	out(1, 5) = (3.0 - 6.0*s - 7.0*t + 8.0*s * t + 4.0*t * t);

	out(0, 6) = (16.0*s - 8.0*s * t - 16.0*s * s);
	out(1, 6) = (8.0*t - 16.0*s * t - 8.0*t * t);

	out(0, 7) = (8.0*s - 16.0*s * t - 8.0*s * s);
	out(1, 7) = (16.0*t - 8.0*s * t - 16.0*t * t);

};
inline void evaluate_raviartthomas_basis_divergence(real const s, real const t, Vector<real> & out) {

	out(0) = (-3.0 + 4.0*t + 8.0*s - 3.0 + 4.0*s + 8.0*t);
	out(1) = (5.0 - 8.0*s + 1.0 - 4.0*s);
	out(2) = (1.0 - 4.0*t + 5.0 - 8.0*t);
	out(3) = (1.0 + 4.0 * t - 8.0*s - 1.0 - 4.0*s + 8.0*t);
	out(4) = (7.0 - 8.0*t - 8.0*s + 5.0 - 4.0*s - 16.0*t);
	out(5) = (-5.0 + 16.0*s + 4.0*t - 7.0 + 8.0*s + 8.0*t);
	out(6) = (16.0 - 8.0*t - 32.0*s + 8.0 - 16.0*s - 16.0*t);
	out(7) = (8.0 - 16.0*s - 16.0*t + 16.0 - 8.0*s - 32.0*t);

};

inline void evaluate_polynomial_basis(real const s, real const t, Vector<real> & out) {

	out(0) = (+1.0);
	out(1) = (-1.0 + 2.0*s);
	out(2) = (-1.0 + 2.0*t);

};
inline void evaluate_polynomial_basis_gradient(real const s, real const t, Matrix<real> & out) {

	out(0, 0) = (0.0);
	out(1, 0) = (0.0);

	out(0, 1) = (2.0);
	out(1, 1) = (0.0);

	out(0, 2) = (0.0);
	out(1, 2) = (2.0);

}

inline void evaluate_edge_polynomial_basis(real const ksi, unsigned const El, Vector<real> & out) {

	switch (El) {

	case 0:
		out(0) = 1.0;
		out(1) = (2.0 / sqrt(2.0)) * (ksi - sqrt(2.0) / 2.0);
		return;
	case 1:
		out(0) = 1.0;
		out(1) = 2.0*(ksi - 0.5);
		return;
	case 2:
		out(0) = 1.0;
		out(1) = 2.0*(ksi - 0.5);
		return;

	}

};

inline void evaluate_edge_normal(Matrix<real> & out) {

	out(0, 0) = 1.0 / sqrt(2.0);
	out(1, 0) = 1.0 / sqrt(2.0);

	out(0, 1) = -1.0;
	out(1, 1) = 0.0;

	out(0, 2) = 0.0;
	out(1, 2) = -1.0;

};
inline void evaluate_edge_parametrization(real const ksi, unsigned const El, Vector<real> & out) {

	switch (El) {

	case 0:
		out(0) = (1.0 - ksi / sqrt(2.0));
		out(1) = (ksi / sqrt(2.0));
		return;
	case 1:
		out(0) = 0.0;
		out(1) = (1.0 - ksi);
		return;
	case 2:
		out(0) = ksi;
		out(1) = 0.0;
		return;

	}

};
inline void evaluate_edge_parametrization_derivative(real const ksi, unsigned const El, Vector<real> & out) {

	switch (El) {

	case 0:
		out(0) = -1.0 / sqrt(2.0);
		out(1) = 1.0 / sqrt(2.0);
		return;
	case 1:
		out(0) = 0.0;
		out(1) = -1.0;
		return;
	case 2:
		out(0) = 1.0;
		out(1) = 0.0;
		return;

	}

};


/*****************************************************************************/
/*                                                                           */
/*    - Eigen containers											         */
/*                                                                           */
/*****************************************************************************/

inline void evaluate_raviartthomas_basis(real const s, real const t, Eigen::MatrixXd & out) {


	out.coeffRef(0, 0) = (real)(-3.0*s + 4.0*s*t + 4.0*s * s);
	out.coeffRef(1, 0) = (real)(-3.0*t + 4.0*s*t + 4.0*t * t);

	out.coeffRef(0, 1) = (real)(-1.0 + 5.0*s - 4.0*s * s);
	out.coeffRef(1, 1) = (real)(t - 4.0*s * t);

	out.coeffRef(0, 2) = (real)(s - 4.0*s * t);
	out.coeffRef(1, 2) = (real)(-1.0 + 5.0*t - 4.0*t * t);

	out.coeffRef(0, 3) = (real)(s + 4.0*t * s - 4.0*s * s);
	out.coeffRef(1, 3) = (real)(-t - 4.0*s * t + 4.0*t * t);

	out.coeffRef(0, 4) = (real)(-3.0 + 7.0*s + 6.0*t - 8.0*t * s - 4.0*s * s);
	out.coeffRef(1, 4) = (real)(5.0*t - 4.0*s * t - 8.0*t * t);

	out.coeffRef(0, 5) = (real)(-5.0*s + 4.0*t * s + 8.0*s * s);
	out.coeffRef(1, 5) = (real)(3.0 - 6.0*s - 7.0*t + 8.0*s * t + 4.0*t * t);

	out.coeffRef(0, 6) = (real)(16.0*s - 8.0*s * t - 16.0*s * s);
	out.coeffRef(1, 6) = (real)(8.0*t - 16.0*s * t - 8.0*t * t);

	out.coeffRef(0, 7) = (real)(8.0*s - 16.0*s * t - 8.0*s * s);
	out.coeffRef(1, 7) = (real)(16.0*t - 8.0*s * t - 16.0*t * t);

};
inline void evaluate_raviartthomas_basis_divergence(real const s, real const t, Eigen::VectorXd & out) {


	out.coeffRef(0) = (real)(-3.0 + 4.0*t + 8.0*s - 3.0 + 4.0*s + 8.0*t);
	out.coeffRef(1) = (real)(5.0 - 8.0*s + 1.0 - 4.0*s);
	out.coeffRef(2) = (real)(1.0 - 4.0*t + 5.0 - 8.0*t);
	out.coeffRef(3) = (real)(1.0 + 4.0 * t - 8.0*s - 1.0 - 4.0*s + 8.0*t);
	out.coeffRef(4) = (real)(7.0 - 8.0*t - 8.0*s + 5.0 - 4.0*s - 16.0*t);
	out.coeffRef(5) = (real)(-5.0 + 16.0*s + 4.0*t - 7.0 + 8.0*s + 8.0*t);
	out.coeffRef(6) = (real)(16.0 - 8.0*t - 32.0*s + 8.0 - 16.0*s - 16.0*t);
	out.coeffRef(7) = (real)(8.0 - 16.0*s - 16.0*t + 16.0 - 8.0*s - 32.0*t);

};

inline void evaluate_polynomial_basis(real const s, real const t, Eigen::VectorXd & out) {

	out.coeffRef(0) = (+1.0);
	out.coeffRef(1) = (-1.0 + 2.0*s);
	out.coeffRef(2) = (-1.0 + 2.0*t);

};
inline void evaluate_polynomial_basis_gradient(real const s, real const t, Eigen::MatrixXd & out) {

	out.coeffRef(0, 0) = (0.0);
	out.coeffRef(1, 0) = (0.0);

	out.coeffRef(0, 1) = (2.0);
	out.coeffRef(1, 1) = (0.0);

	out.coeffRef(0, 2) = (0.0);
	out.coeffRef(1, 2) = (2.0);

};

inline void evaluate_edge_polynomial_basis(real const ksi, unsigned const El, Eigen::VectorXd & out) {


	switch (El) {

	case 0:
		out.coeffRef(0) = 1.0;
		out.coeffRef(1) = (2.0 / sqrt(2.0)) * (ksi - sqrt(2.0) / 2.0);
		return;
	case 1:
		out.coeffRef(0) = 1.0;
		out.coeffRef(1) = 2.0*(ksi - 0.5);
		return;
	case 2:
		out.coeffRef(0) = 1.0;
		out.coeffRef(1) = 2.0*(ksi - 0.5);
		return;

	}

};

inline void evaluate_edge_normal(Eigen::MatrixXd & out) {


	out.coeffRef(0, 0) = 1.0 / sqrt(2.0);
	out.coeffRef(1, 0) = 1.0 / sqrt(2.0);

	out.coeffRef(0, 1) = -1.0;
	out.coeffRef(1, 1) = 0.0;

	out.coeffRef(0, 2) = 0.0;
	out.coeffRef(1, 2) = -1.0;

};
inline void evaluate_edge_parametrization(real const ksi, unsigned const El, Eigen::VectorXd & out) {


	switch (El) {

	case 0:
		out.coeffRef(0) = (1.0 - ksi / sqrt(2.0));
		out.coeffRef(1) = (ksi / sqrt(2.0));
		return;
	case 1:
		out.coeffRef(0) = 0.0;
		out.coeffRef(1) = (1.0 - ksi);
		return;
	case 2:
		out.coeffRef(0) = ksi;
		out.coeffRef(1) = 0.0;
		return;

	}

};
inline void evaluate_edge_parametrization2(real const ksi, unsigned const El, Eigen::VectorXd & out) {


	switch (El) {

	case 0:
		out.coeffRef(0) = ksi / sqrt(2.0);
		out.coeffRef(1) = 1.0 - ksi / sqrt(2.0); 
		return;
	case 1:
		out.coeffRef(0) = 0.0;
		out.coeffRef(1) = ksi;
		return;
	case 2:
		out.coeffRef(0) = 1.0 - ksi;
		out.coeffRef(1) = 0.0;
		return;

	}

};
inline void evaluate_edge_parametrization_derivative(real const ksi, unsigned const El, Eigen::VectorXd & out) {


	switch (El) {

	case 0:
		out.coeffRef(0) = -1.0 / sqrt(2.0);
		out.coeffRef(1) = 1.0 / sqrt(2.0);
		return;
	case 1:
		out.coeffRef(0) = 0.0;
		out.coeffRef(1) = -1.0;
		return;
	case 2:
		out.coeffRef(0) = 1.0;
		out.coeffRef(1) = 0.0;
		return;

	}

};

inline Eigen::Vector2d const evaluate_edge_normal(unsigned const El) {


	Eigen::Vector2d normal;

	switch (El) {

	case 0:
		normal.coeffRef(0) = 1.0 / sqrt(2.0);
		normal.coeffRef(1) = 1.0 / sqrt(2.0);
		return normal;
	case 1:
		normal.coeffRef(0) = -1.0;
		normal.coeffRef(1) = 0.0;
		return normal;
	case 2:
		normal.coeffRef(0) = 0.0;
		normal.coeffRef(1) = -1.0;
		return normal;
	}

	return normal;
};



/*


real solver::velocity_in_normal_direction(t_pointer const K, e_pointer const E, real const x) {


	real orientations[3] = { 1.0, 1.0, 1.0 };
	real orientations2[3] = { 1.0, 1.0, 1.0 };
	Eigen::VectorXd parametrization(2);
	//Vector<real> parametrization(2);
	Eigen::MatrixXd normals(2, 3);
	Eigen::Vector2d normal;
	Eigen::MatrixXd JF(2, 2);
	Eigen::MatrixXd itJF(2, 2);
	Eigen::MatrixXd basisRaviartThomas(2, 8);


	unsigned const e_index_local = K->get_edge_index(E);

	//for (unsigned e = 0; e < 3; e++)
	//	orientations[e] = edgeOrientation(K->index, 0, e, 0);
	//for (unsigned e = 0; e < 3; e++)
	//	orientations2[e] = edgeOrientation2(K->index, 0, e, 0);

	//real const orient = edgeOrientation(K->index, 0, e_index_local, 0);
	//real const orient2 = edgeOrientation2(K->index, 0, e_index_local, 0);


	real const A = (real) 0.0;
	real const B = (real)(e_index_local != 0) ? 1.0 : sqrt(2.0);

	real const C = (real)(B - A) / 2.0;
	real const D = (real)(B + A) / 2.0;

	real const X = (real)x * C + D;


	v_pointer const a = K->vertices[0];
	v_pointer const b = K->vertices[1];
	v_pointer const c = K->vertices[2];

	real const x0 = (real)a->x;
	real const y0 = (real)a->y;

	real const x1 = (real)b->x;
	real const y1 = (real)b->y;

	real const x2 = (real)c->x;
	real const y2 = (real)c->y;

	real const detJF = (real)abs((x1 - x0)*(y2 - y0) - (x2 - x0)*(y1 - y0));

	JF(0, 0) = x1 - x0;
	JF(0, 1) = x2 - x0;
	JF(1, 0) = y1 - y0;
	JF(1, 1) = y2 - y0;


	//evaluate_edge_parametrization(X, e_index_local, parametrization);
	evaluate_edge_parametrization(X, e_index_local, parametrization);
	evaluate_edge_normal(normals);

	real const s = parametrization(0);
	real const t = parametrization(1);

	evaluate_raviartthomas_basis(s, t, basisRaviartThomas);

	

	itJF = (JF.inverse()).transpose();

	normal = itJF * (normals.col(e_index_local));

	real val = 0.0;

	for (unsigned dof = 0; dof < 2; dof++) {

		unsigned const j = LI(K, E, dof);

		real const dotProduct = normal.dot(JF * basisRaviartThomas.col(j)) / (detJF * (normal).norm());

		//real const dotProduct = normal.dot(basisRaviartThomas.col(j));

		val += velocities(K->index, 0, j, 0) * dotProduct;

	}

	return val;

};

real solver::continuity_of_normal_component(t_pointer const K, e_pointer const E, real const x) {


	//real orientations[3];
	Vector<real> parametrization(2);
	Matrix<real> normals(2, 3);
	Eigen::Vector2d normal;
	Eigen::MatrixXd JF(2, 2);
	Eigen::MatrixXd itJF(2, 2);
	Eigen::MatrixXd basisRaviartThomas(2, 8);


	unsigned const e_index_local = K->get_edge_index(E);

	//for (unsigned e = 0; e < 3; e++)
	//	orientations[e] = edgeOrientation(K->index, 0, e, 0);

	//real const orient = edgeOrientation(K->index, 0, e_index_local, 0);


	real const A = (real) 0.0;
	real const B = (real)(e_index_local != 0) ? 1.0 : sqrt(2.0);

	real const C = (real)(B - A) / 2.0;
	real const D = (real)(B + A) / 2.0;

	real const X = (real)x * C + D;


	v_pointer const a = K->vertices[0];
	v_pointer const b = K->vertices[1];
	v_pointer const c = K->vertices[2];

	real const x0 = (real)a->x;
	real const y0 = (real)a->y;

	real const x1 = (real)b->x;
	real const y1 = (real)b->y;

	real const x2 = (real)c->x;
	real const y2 = (real)c->y;

	real const detJF = (real)abs((x1 - x0)*(y2 - y0) - (x2 - x0)*(y1 - y0));

	JF(0, 0) = x1 - x0;
	JF(0, 1) = x2 - x0;
	JF(1, 0) = y1 - y0;
	JF(1, 1) = y2 - y0;


	evaluate_edge_parametrization(X, e_index_local, parametrization);
	evaluate_edge_normal(normals);

	real const s = parametrization(0);
	real const t = parametrization(1);

	evaluate_raviartthomas_basis(s, t, basisRaviartThomas);

	normal(0) = normals.getColumn(e_index_local)(0);
	normal(1) = normals.getColumn(e_index_local)(1);

	itJF = (JF.inverse()).transpose();



	real const constant = 1.0 / (detJF * (itJF * normal).norm());

	real val = 0.0;

	//for (unsigned dof = 0; dof < 2; dof++) {

	//	unsigned const j = LI(K, E, dof);

	//	val += orient * normal.dot(basisRaviartThomas.col(j));

	//}

	unsigned const j = LI(K, E, 0);
	val += normal.dot(basisRaviartThomas.col(j));

	val *= constant;

	return val;

};

*/
