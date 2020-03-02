#pragma once


#pragma once


#include "enumerators.h"

#include <vector>
#include "shapes.h"
#include "geometric.h"

#include <cassert>
#include <algorithm>
#include <iostream>






class MESHVertex;
class MESHEdge;
class MESHTriangle;


typedef MESHVertex *	vm_pointer;
typedef MESHEdge *		em_pointer;
typedef MESHTriangle *	tm_pointer;


class MESHVertex {


	/*
	friend class Triangle;

	friend class Edge;



	friend class Mesh;

	friend class PlanarStraightLineGraph;

	template<GEOMETRIC_KERNEL GK>
	friend class Triangulation;

	template <GEOMETRIC_KERNEL GK>
	friend class Predicates;


	//			Friend functions		 
	friend bool t_compare_x(Triangle * const t1, Triangle * const t2);
	friend bool t_compare_y(Triangle * const t1, Triangle * const t2);
	friend bool v_compare_x(Vertex * const v1, Vertex * const v2);
	friend bool v_compare_y(Vertex * const v1, Vertex * const v2);
	friend bool e_compare_x(Edge * const e1, Edge * const e2);
	friend bool e_compare_y(Edge * const e1, Edge * const e2);
	*/


public:


	/*****************************************************************************/
	/*                                                                           */
	/*    - Class members														 */
	/*                                                                           */
	/*****************************************************************************/
	double x;
	double y;

	unsigned index;


	/*****************************************************************************/
	/*                                                                           */
	/*    - Class constructor/destructor										 */
	/*                                                                           */
	/*****************************************************************************/
	MESHVertex(double const & const X, double const & const Y);
	~MESHVertex();


};

class MESHEdge {


	/*
	friend class Vertex;

	friend class Triangle;



	friend class Mesh;

	template<GEOMETRIC_KERNEL GK>
	friend class Triangulation;

	template <GEOMETRIC_KERNEL GK>
	friend class Predicates;


	//				Friend functions	
	friend bool e_compare_x(Edge * const e1, Edge * const e2);
	friend bool e_compare_y(Edge * const e1, Edge * const e2);
	friend bool marker_compare_neumann(Edge * const e1, Edge * const e2);
	friend bool marker_compare_dirichlet(Edge * const e1, Edge * const e2);

	*/

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
	MESHEdge(vm_pointer & const va, vm_pointer & const vb);
	~MESHEdge();


	/*****************************************************************************/
	/*                                                                           */
	/*    - Class methods														 */
	/*                                                                           */
	/*****************************************************************************/
	double const & length() const;

};

class MESHTriangle {


	/*
	friend class Vertex;

	friend class Edge;

	friend class Mesh;

	template <GEOMETRIC_KERNEL GK>
	friend class Triangulation;

	template <GEOMETRIC_KERNEL GK>
	friend class Predicates;

	//				Friend functions		 

	friend bool t_compare_x(Triangle * const t1, Triangle * const t2);
	friend bool t_compare_y(Triangle * const t1, Triangle * const t2);

	*/


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
	T_MARKER marker;


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
	vm_pointer const & get_vertex(em_pointer const & e) const;
	vm_pointer const & get_vertex_cw(vm_pointer const & v)  const;
	vm_pointer const & get_vertex_ccw(vm_pointer const & v) const;

	unsigned const & get_edge_index(em_pointer const & e) const;

	double const & area() const;


};



/*****************************************************************************/
/*																			 */
/*									VERTEX									 */
/*																			 */
/*****************************************************************************/
MESHVertex::MESHVertex(double const & const X, double const & const Y) : x(X), y(Y), index(0)  {

};
MESHVertex::~MESHVertex() {

};


/*****************************************************************************/
/*																			 */
/*									EDGE									 */
/*																			 */
/*****************************************************************************/
MESHEdge::MESHEdge(vm_pointer & const va, vm_pointer & const vb) : a(va), b(vb) {

};
MESHEdge::~MESHEdge() {

	a = NULL;
	b = NULL;

	neighbors[0] = NULL;
	neighbors[1] = NULL;

};

double const & MESHEdge::length() const {

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
MESHTriangle::MESHTriangle(vm_pointer const & va, vm_pointer const & vb, vm_pointer const & vc) : vertices{ va,vb,vc }, index(0), marker(T_MARKER::NONE) {};
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

vm_pointer const & MESHTriangle::get_vertex(em_pointer const & e) const {

	if		(e == edges[0]) return vertices[0];
	else if (e == edges[1]) return vertices[1];
	else if (e == edges[2]) return vertices[2];

	assert(0 && "Triangle::get_vertex(Edge)");

	return NULL;

};
vm_pointer const & MESHTriangle::get_vertex_cw(vm_pointer const & v) const {

	if		(vertices[0] == v) return vertices[2];
	else if (vertices[1] == v) return vertices[0];
	else if (vertices[2] == v) return vertices[1];

	assert(0 && "Triangle::vertex_cw(Vertex)");

	return NULL;

};
vm_pointer const & MESHTriangle::get_vertex_ccw(vm_pointer const & v) const {

	if		(vertices[0] == v) return vertices[1];
	else if (vertices[1] == v) return vertices[2];
	else if (vertices[2] == v) return vertices[0];

	assert(0 && "Triangle::vertex_ccw(Vertex)");

	return NULL;

};

unsigned const & MESHTriangle::get_edge_index(em_pointer const & e) const {

	if		(e == edges[0]) return 0;
	else if (e == edges[1]) return 1;
	else if (e == edges[2]) return 2;

	assert(0 && "Triangle::edge_index(Edge)");

	return -1;

};

double const & MESHTriangle::area() const {

	vm_pointer const va = vertices[0];
	vm_pointer const vb = vertices[1];
	vm_pointer const vc = vertices[2];

	return 0.5*((va->x - vc->x) * (vb->y - vc->y) - (va->y - vc->y) * (vb->x - vc->x));

};
