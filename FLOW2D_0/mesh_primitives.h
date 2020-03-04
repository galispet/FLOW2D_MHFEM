#pragma once


#pragma once


#include "enumerators.h"

#include <cassert>
#include <algorithm>
#include <iostream>




/*
template<typename real>
class MESHVertex;
template<typename real>
class MESHEdge;
template<typename real>
class MESHTriangle;

template<typename real>
using vm_pointer = MESHVertex<real> *;
template<typename real>
using em_pointer = MESHEdge<real> *;
template<typename real>
using tm_pointer = MESHTriangle<real> *;
*/


class MESHVertex;
class MESHEdge;
class MESHTriangle;

typedef MESHVertex *	vm_pointer;
typedef MESHEdge *		em_pointer;
typedef MESHTriangle *	tm_pointer;



class MESHVertex {


	
	/*friend class Triangle;

	friend class Edge;



	friend class Mesh;

	*/
friend class PlanarStraightLineGraph;

	template<GEOMETRIC_KERNEL GK>
	friend class Triangulation;

	template <GEOMETRIC_KERNEL GK>
	friend class Predicates;
	friend class Mesh2;

	friend class MESHTriangle;
	friend class MESHEdge;

	//			Friend functions		 
	friend bool t_compare_x(tm_pointer const t1, tm_pointer const t2);
	friend bool t_compare_y(tm_pointer const t1, tm_pointer const t2);
	friend bool v_compare_x(vm_pointer const v1, vm_pointer const v2);
	friend bool v_compare_y(vm_pointer const v1, vm_pointer const v2);
	friend bool e_compare_x(em_pointer const e1, em_pointer const e2);
	friend bool e_compare_y(em_pointer const e1, em_pointer const e2);
	


public:


	/*****************************************************************************/
	/*                                                                           */
	/*    - Class members														 */
	/*                                                                           */
	/*****************************************************************************/
	double const x;
	double const y;

	unsigned index;


	/*****************************************************************************/
	/*                                                                           */
	/*    - Class constructor/destructor										 */
	/*                                                                           */
	/*****************************************************************************/
	MESHVertex(double const X, double const Y);
	~MESHVertex();


};

class MESHEdge {


	/*
	friend class Vertex;

	friend class Triangle;



	friend class Mesh;



	//				Friend functions	
	friend bool e_compare_x(Edge * const e1, Edge * const e2);
	friend bool e_compare_y(Edge * const e1, Edge * const e2);


	*/

	template<GEOMETRIC_KERNEL GK>
	friend class Triangulation;

	template <GEOMETRIC_KERNEL GK>
	friend class Predicates;
	friend class Mesh2;

	friend class MESHVertex;
	friend class MESHTriangle;
	friend class MESHEdge;

	friend bool t_compare_x(tm_pointer const t1, tm_pointer const t2);
	friend bool t_compare_y(tm_pointer const t1, tm_pointer const t2);
	friend bool v_compare_x(vm_pointer const v1, vm_pointer const v2);
	friend bool v_compare_y(vm_pointer const v1, vm_pointer const v2);
	friend bool e_compare_x(em_pointer const e1, em_pointer const e2);
	friend bool e_compare_y(em_pointer const e1, em_pointer const e2);

	friend bool marker_compare_neumann(em_pointer const e1, em_pointer const e2);
	friend bool marker_compare_dirichlet(em_pointer const e1, em_pointer const e2);
public:
	

	/*****************************************************************************/
	/*                                                                           */
	/*    - Class members														 */
	/*                                                                           */
	/*****************************************************************************/
	vm_pointer a;
	vm_pointer b;

	unsigned index;
	E_MARKER marker;

	tm_pointer neighbors[2] = { NULL,NULL };


	/*****************************************************************************/
	/*                                                                           */
	/*    - Class constructor/destructor										 */
	/*                                                                           */
	/*****************************************************************************/
	MESHEdge(vm_pointer const & va, vm_pointer const & vb);
	~MESHEdge();


	/*****************************************************************************/
	/*                                                                           */
	/*    - Class methods														 */
	/*                                                                           */
	/*****************************************************************************/
	double const length() const;

};

class MESHTriangle {


	/*
	friend class Vertex;

	friend class Edge;

	friend class Mesh;



	//				Friend functions		 

	friend bool t_compare_x(Triangle * const t1, Triangle * const t2);
	friend bool t_compare_y(Triangle * const t1, Triangle * const t2);

	*/
	template <GEOMETRIC_KERNEL GK>
	friend class Triangulation;

	template <GEOMETRIC_KERNEL GK>
	friend class Predicates;
	friend class Mesh2;

	friend class MESHVertex;
	friend class MESHTriangle;
	friend class MESHEdge;

	friend bool t_compare_x(tm_pointer const t1, tm_pointer const t2);
	friend bool t_compare_y(tm_pointer const t1, tm_pointer const t2);
	friend bool v_compare_x(vm_pointer const v1, vm_pointer const v2);
	friend bool v_compare_y(vm_pointer const v1, vm_pointer const v2);
	friend bool e_compare_x(em_pointer const e1, em_pointer const e2);
	friend bool e_compare_y(em_pointer const e1, em_pointer const e2);


public:


	/*****************************************************************************/
	/*                                                                           */
	/*    - Class members														 */
	/*                                                                           */
	/*****************************************************************************/
	vm_pointer vertices[3]	= { NULL,NULL,NULL };
	em_pointer edges[3]		= { NULL,NULL,NULL };
	tm_pointer neighbors[3] = { NULL,NULL,NULL };

	unsigned index;


	/*****************************************************************************/
	/*                                                                           */
	/*    - Class constructor/destructor										 */
	/*                                                                           */
	/*****************************************************************************/
	MESHTriangle(vm_pointer const & va, vm_pointer const & vb, vm_pointer const & vc);
	~MESHTriangle();



	/*****************************************************************************/
	/*                                                                           */
	/*    - Class methods														 */
	/*                                                                           */
	/*****************************************************************************/
	vm_pointer const & get_vertex(unsigned const i) const;
	vm_pointer const & get_vertex(em_pointer const & e) const;
	vm_pointer const & get_vertex_cw(vm_pointer const & v)  const;
	vm_pointer const & get_vertex_ccw(vm_pointer const & v) const;

	unsigned const get_edge_index(em_pointer const & e) const;

	double const area() const;


};



/*****************************************************************************/
/*																			 */
/*									VERTEX									 */
/*																			 */
/*****************************************************************************/
MESHVertex::MESHVertex(double const X, double const Y) : x(X), y(Y), index(0)  {

};
MESHVertex::~MESHVertex() {

};


/*****************************************************************************/
/*																			 */
/*									EDGE									 */
/*																			 */
/*****************************************************************************/
MESHEdge::MESHEdge(vm_pointer const & va, vm_pointer const & vb) : a(va), b(vb) {

};
MESHEdge::~MESHEdge() {

	a = NULL;
	b = NULL;

	neighbors[0] = NULL;
	neighbors[1] = NULL;

};

double const MESHEdge::length() const {

	double const x0 = a->x;
	double const y0 = a->y;

	double const x1 = b->x;
	double const y1 = b->y;

	return sqrt(sqr(x1 - x0) + sqr(y1 - y0));

};


/*****************************************************************************/
/*																			 */
/*									TRIANGLE								 */
/*																			 */
/*****************************************************************************/
MESHTriangle::MESHTriangle(vm_pointer const & va, vm_pointer const & vb, vm_pointer const & vc) : vertices{ va,vb,vc } {

};
MESHTriangle::~MESHTriangle() {

	vertices[0] = NULL;
	vertices[1] = NULL;
	vertices[2] = NULL;

	edges[0] = NULL;
	edges[1] = NULL;
	edges[2] = NULL;
	
	neighbors[0] = NULL;
	neighbors[1] = NULL;
	neighbors[2] = NULL;

};

vm_pointer const & MESHTriangle::get_vertex(unsigned const i) const {

	return vertices[i];

};
vm_pointer const & MESHTriangle::get_vertex(em_pointer const & e) const {

	if		(e == edges[0]) return vertices[0];
	else if (e == edges[1]) return vertices[1];
	else if (e == edges[2]) return vertices[2];

	assert(0 && "Triangle::get_vertex(Edge)");

	return vertices[0];

};
vm_pointer const & MESHTriangle::get_vertex_cw(vm_pointer const & v) const {

	if		(vertices[0] == v) return vertices[2];
	else if (vertices[1] == v) return vertices[0];
	else if (vertices[2] == v) return vertices[1];

	assert(0 && "Triangle::vertex_cw(Vertex)");

	return vertices[0];

};
vm_pointer const & MESHTriangle::get_vertex_ccw(vm_pointer const & v) const {

	if		(vertices[0] == v) return vertices[1];
	else if (vertices[1] == v) return vertices[2];
	else if (vertices[2] == v) return vertices[0];

	assert(0 && "Triangle::vertex_ccw(Vertex)");

	return vertices[0];

};

unsigned const MESHTriangle::get_edge_index(em_pointer const & e) const {

	if		(e == edges[0]) return 0;
	else if (e == edges[1]) return 1;
	else if (e == edges[2]) return 2;

	assert(0 && "Triangle::edge_index(Edge)");

	return -1;

};

double const MESHTriangle::area() const {

	vm_pointer const va = vertices[0];
	vm_pointer const vb = vertices[1];
	vm_pointer const vc = vertices[2];

	return 0.5*((va->x - vc->x) * (vb->y - vc->y) - (va->y - vc->y) * (vb->x - vc->x));

};
