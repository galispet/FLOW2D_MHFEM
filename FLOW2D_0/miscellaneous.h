#pragma once


#include <Eigen/Dense>

#include "mesh_primitives.h"
#include "integration.h"



/*****************************************************************************/
/*                                                                           */
/*    - Miscellaneous														 */
/*                                                                           */
/*****************************************************************************/
template<typename T>
inline T square(T const x) {

	return x * x;

};
template<typename T>
inline T kroneckerDelta(unsigned const i, unsigned const j) {

	return i == j ? 1.0 : 0.0;

};
template<typename T>
inline void hardCopy(T * & copyTo, T * & copyFrom, unsigned const n) {

	for (unsigned i = 0; i < n; i++)
		copyTo[i] = copyFrom[i];

};

template<typename T>
inline T barenblatt(T const x, T const y, T const time) {

	return (1.0 / sqrt(time))*fmax(1.0 - (x*x + y * y) / (16.0 * sqrt(time)), 0.0);

};



/*****************************************************************************/
/*                                                                           */
/*    - Gauss numerical integration over triangle area						 */
/*                                                                           */
/*****************************************************************************/
template<typename T>
T integrateTriangle(tm_pointer const K, T time, T(*fun)(T, T, T)) {



	vm_pointer const a = K->vertices[0];
	vm_pointer const b = K->vertices[1];
	vm_pointer const c = K->vertices[2];

	T const x0 = a->x;
	T const y0 = a->y;

	T const x1 = b->x;
	T const y1 = b->y;

	T const x2 = c->x;
	T const y2 = c->y;

	T const detJF = abs((x1 - x0)*(y2 - y0) - (x2 - x0)*(y1 - y0));


	quadrature_triangle quad(quadrature_order);
	unsigned const num_quad_points = quad.NumberOfPoints;


	T integral = 0.0;

	for (unsigned n = 0; n < num_quad_points; n++) {


		T const s = quad.points_x[n];
		T const t = quad.points_y[n];
		T const w = quad.weights[n];

		// Gauss points on given triangle. This is the transformation : [-1,1] x [-1,1] square --> reference triangle ([0,0],[1,0],[0,1]) --> given triangle ([x0y,0],[x1,y1],[x2,y2]) 
		T const x = x0 + (x1 - x0)*s + (x2 - x0)*t;
		T const y = y0 + (y1 - y0)*s + (y2 - y0)*t;

		// Gauss quadrature
		integral += w * fun(x, y, time);


	}

	return  0.5 * detJF * integral;

};
template<typename T>
T integrateTriangle(tm_pointer const K, T(*fun)(T, T)) {


	vm_pointer const a = K->vertices[0];
	vm_pointer const b = K->vertices[1];
	vm_pointer const c = K->vertices[2];

	T const x0 = a->x;
	T const y0 = a->y;

	T const x1 = b->x;
	T const y1 = b->y;

	T const x2 = c->x;
	T const y2 = c->y;

	T const detJF = abs((x1 - x0)*(y2 - y0) - (x2 - x0)*(y1 - y0));


	quadrature_triangle quad(quadrature_order);
	unsigned const num_quad_points = quad.NumberOfPoints;


	T integral = 0.0;

	for (unsigned n = 0; n < num_quad_points; n++) {


		T const s = quad.points_x[n];
		T const t = quad.points_y[n];

		T const w = quad.weights[n];

		// Gauss points on given triangle. This is the transformation : [-1,1] x [-1,1] square --> reference triangle ([0,0],[1,0],[0,1]) --> given triangle ([x0y,0],[x1,y1],[x2,y2]) 
		T const x = x0 + (x1 - x0)*s + (x2 - x0)*t;
		T const y = y0 + (y1 - y0)*s + (y2 - y0)*t;

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
template<typename T>
double integrateEdge(em_pointer const & E, T time, T(*fun)(T, T, T)) {


	T const x0 = E->a->x;
	T const y0 = E->a->y;

	T const x1 = E->b->x;
	T const y1 = E->b->y;

	T const length = E->length();


	gauss_quadrature_1D quad(quadrature_order);
	unsigned const num_quad_points = quad.NumberOfPoints;


	T integral = 0.0;

	for (unsigned n = 0; n < num_quad_points; n++) {


		T const x = quad.points[n];
		T const w = quad.weights[n];

		T const X = x0 + 0.5*(1.0 + x)*(x1 - x0);
		T const Y = y0 + 0.5*(1.0 + x)*(y1 - y0);

		integral += w * fun(X, Y, time);

	}

	return 0.5 * length * integral;

};
template<typename T>
double integrateEdge(em_pointer const & E, T(*fun)(T, T)) {


	T const x0 = E->a->x;
	T const y0 = E->a->y;

	T const x1 = E->b->x;
	T const y1 = E->b->y;

	T const length = E->length();


	gauss_quadrature_1D quad(quadrature_order);
	unsigned const num_quad_points = quad.NumberOfPoints;


	T integral = 0.0;

	for (unsigned n = 0; n < num_quad_points; n++) {


		T const x = quad.points[n];
		T const w = quad.weights[n];

		T const X = x0 + 0.5*(1.0 + x)*(x1 - x0);
		T const Y = y0 + 0.5*(1.0 + x)*(y1 - y0);

		integral += w * fun(X, Y);

	}

	return 0.5 * length * integral;

};








/*****************************************************************************/
/*                                                                           */
/*    - Model properties													 */
/*                                                                           */
/*****************************************************************************/
template<typename T>
inline T equationOfState(T const c) {

	return c;

};
template<typename T>
inline T equationOfStateDerivative(T const c) {

	return 1.0;

};

template<typename T>
inline T porosity(T const x, T const y) {

	return 1.0;

};
template<typename T>
inline T viscosity(T const x, T const y) {

	return 0.5;

};
template<typename T>
inline T source(T const x, T const y, T const time) {

	return 0.0;

};
template<typename T>
inline void permeability(T const x, T const y, Eigen::Matrix<T, 2, 2> & out) {

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
template<typename T>
inline T NEUMANN_GAMMA_Q_velocity(em_pointer const & E, T const time) {

	return 0.0;

};
template<typename T>
inline T DIRICHLET_GAMMA_Q_concentration(em_pointer const & E, real const time) {

	return 0.0;

};
template<typename T>
inline T DIRICHLET_GAMMA_P_concentration(em_pointer const & E, real const time) {

	return integrateEdge(E, time, barenblatt) / E->length();

};
template<typename T>
inline T DIRICHLET_GAMMA_P_pressure(em_pointer const & E, real const time) {

	return integrateEdge(E, time, barenblatt) / E->length();

};

template<typename T>
inline T NEUMANN_GAMMA_Q_velocity(T const x, T const y, T const time) {

	return 0.0;

};
template<typename T>
inline T DIRICHLET_GAMMA_Q_concentration(T const x, T const y, T const time) {

	return 0.0;

}; 
template<typename T>
inline T DIRICHLET_GAMMA_P_concentration(T const x, T const y, T const time) {

	return barenblatt(x, y, time);

};
template<typename T>
inline T DIRICHLET_GAMMA_P_pressure(T const x, T const y, T const time) {

	return barenblatt(x, y, time);

};






