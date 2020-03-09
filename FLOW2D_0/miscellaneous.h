#pragma once


#include "mesh_primitives.h"
#include "integration.h"

#include <Eigen/Dense>






/*****************************************************************************/
/*                                                                           */
/*    - Miscellaneous														 */
/*                                                                           */
/*****************************************************************************/
inline Real square(Real const x) {

	return x * x;

};
inline Real kroneckerDelta(unsigned const i, unsigned const j) {

	return i == j ? 1.0 : 0.0;

};
inline void hardCopy(Real * & copyTo, Real * & copyFrom, unsigned const n) {

	for (unsigned i = 0; i < n; i++)
		copyTo[i] = copyFrom[i];

};

inline Real barenblatt(Real const x, Real const y, Real const time) {

	return (1.0 / sqrt(time))*fmax(1.0 - (x*x + y * y) / (16.0 * sqrt(time)), 0.0);

};



/*****************************************************************************/
/*                                                                           */
/*    - Gauss numerical integration over triangle area						 */
/*                                                                           */
/*****************************************************************************/
template<unsigned QuadraturePrecision>
Real integrateTriangle(tm_pointer const K, Real time, Real(*fun)(Real, Real, Real)) {



	vm_pointer const a = K->vertices[0];
	vm_pointer const b = K->vertices[1];
	vm_pointer const c = K->vertices[2];

	Real const x0 = a->x;
	Real const y0 = a->y;

	Real const x1 = b->x;
	Real const y1 = b->y;

	Real const x2 = c->x;
	Real const y2 = c->y;

	Real const detJF = abs((x1 - x0)*(y2 - y0) - (x2 - x0)*(y1 - y0));


	quadrature_triangle<Real> quad(QuadraturePrecision);
	unsigned const num_quad_points = quad.NumberOfPoints;


	Real integral = 0.0;

	for (unsigned n = 0; n < num_quad_points; n++) {


		Real const s = quad.points_x[n];
		Real const t = quad.points_y[n];
		Real const w = quad.weights[n];

		// Gauss points on given triangle. This is the transformation : [-1,1] x [-1,1] square --> reference triangle ([0,0],[1,0],[0,1]) --> given triangle ([x0y,0],[x1,y1],[x2,y2]) 
		Real const x = x0 + (x1 - x0)*s + (x2 - x0)*t;
		Real const y = y0 + (y1 - y0)*s + (y2 - y0)*t;

		// Gauss quadrature
		integral += w * fun(x, y, time);


	}

	return  0.5 * detJF * integral;

};
template<unsigned QuadraturePrecision>
Real integrateTriangle(tm_pointer const K, Real(*fun)(Real, Real)) {


	vm_pointer const a = K->vertices[0];
	vm_pointer const b = K->vertices[1];
	vm_pointer const c = K->vertices[2];

	Real const x0 = a->x;
	Real const y0 = a->y;

	Real const x1 = b->x;
	Real const y1 = b->y;

	Real const x2 = c->x;
	Real const y2 = c->y;

	Real const detJF = abs((x1 - x0)*(y2 - y0) - (x2 - x0)*(y1 - y0));


	quadrature_triangle<Real> quad(QuadraturePrecision);
	unsigned const num_quad_points = quad.NumberOfPoints;


	Real integral = 0.0;

	for (unsigned n = 0; n < num_quad_points; n++) {


		Real const s = quad.points_x[n];
		Real const t = quad.points_y[n];

		Real const w = quad.weights[n];

		// Gauss points on given triangle. This is the transformation : [-1,1] x [-1,1] square --> reference triangle ([0,0],[1,0],[0,1]) --> given triangle ([x0y,0],[x1,y1],[x2,y2]) 
		Real const x = x0 + (x1 - x0)*s + (x2 - x0)*t;
		Real const y = y0 + (y1 - y0)*s + (y2 - y0)*t;

		// Gauss quadrature
		integral += w * fun(x, y);

	}

	return 0.5 * detJF * integral;

};



/*****************************************************************************/
/*                                                                           */
/*    - Gauss numerical integration over edge								 */
/*                                                                           */
/*****************************************************************************/
template<unsigned QuadraturePrecision>
Real integrateEdge(em_pointer const & E, Real time, Real(*fun)(Real, Real, Real)) {


	Real const x0 = E->a->x;
	Real const y0 = E->a->y;

	Real const x1 = E->b->x;
	Real const y1 = E->b->y;

	Real const length = E->length();


	gauss_quadrature_1D<Real> quad(QuadraturePrecision);
	unsigned const num_quad_points = quad.NumberOfPoints;


	Real integral = 0.0;

	for (unsigned n = 0; n < num_quad_points; n++) {


		Real const x = quad.points[n];
		Real const w = quad.weights[n];

		Real const X = x0 + 0.5*(1.0 + x)*(x1 - x0);
		Real const Y = y0 + 0.5*(1.0 + x)*(y1 - y0);

		integral += w * fun(X, Y, time);

	}

	return 0.5 * length * integral;

};
template<unsigned QuadraturePrecision>
Real integrateEdge(em_pointer const & E, Real(*fun)(Real, Real)) {


	Real const x0 = E->a->x;
	Real const y0 = E->a->y;

	Real const x1 = E->b->x;
	Real const y1 = E->b->y;

	Real const length = E->length();


	gauss_quadrature_1D<Real> quad(QuadraturePrecision);
	unsigned const num_quad_points = quad.NumberOfPoints;


	Real integral = 0.0;

	for (unsigned n = 0; n < num_quad_points; n++) {


		Real const x = quad.points[n];
		Real const w = quad.weights[n];

		Real const X = x0 + 0.5*(1.0 + x)*(x1 - x0);
		Real const Y = y0 + 0.5*(1.0 + x)*(y1 - y0);

		integral += w * fun(X, Y);

	}

	return 0.5 * length * integral;

};



/*****************************************************************************/
/*                                                                           */
/*    - Model properties													 */
/*                                                                           */
/*****************************************************************************/
inline Real equationOfState(Real const c) {

	return c;

};
inline Real equationOfStateDerivative(Real const c) {

	return 1.0;

};

inline Real porosity(Real const x, Real const y) {

	return 1.0;

};
inline Real viscosity(Real const x, Real const y) {

	return 0.5;

};
inline Real source(Real const x, Real const y, Real const time) {

	return 0.0;

};
inline void permeability(Real const x, Real const y, Eigen::Matrix<Real, 2, 2> & out) {

	out(0, 0) = 1.0;
	out(0, 1) = 0.0;
	out(1, 0) = 0.0;
	out(1, 1) = 1.0;

};



/*****************************************************************************/
/*                                                                           */
/*    - Boundary conditions (Neumann, Dirichlet)							 */
/*                                                                           */
/*****************************************************************************/
inline Real NEUMANN_GAMMA_Q_velocity(em_pointer const & E, Real const time) {

	return 0.0;

};
inline Real DIRICHLET_GAMMA_Q_concentration(em_pointer const & E, Real const time) {

	return 0.0;

};
template<unsigned QuadraturePrecision>
inline Real DIRICHLET_GAMMA_P_concentration(em_pointer const & E, Real const time) {

	return integrateEdge<QuadraturePrecision>(E, time, barenblatt) / E->length();

};
template<unsigned QuadraturePrecision>
inline Real DIRICHLET_GAMMA_P_pressure(em_pointer const & E, Real const time) {

	return integrateEdge<QuadraturePrecision>(E, time, barenblatt) / E->length();

};


inline Real NEUMANN_GAMMA_Q_velocity(Real const x, Real const y, Real const time) {

	return 0.0;

};
inline Real DIRICHLET_GAMMA_Q_concentration(Real const x, Real const y, Real const time) {

	return 0.0;

};
inline Real DIRICHLET_GAMMA_P_concentration(Real const x, Real const y, Real const time) {

	return barenblatt(x, y, time);

};
inline Real DIRICHLET_GAMMA_P_pressure(Real const x, Real const y, Real const time) {

	return barenblatt(x, y, time);

};





