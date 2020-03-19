#pragma once



#include "enumerators.h"

#include <cassert>
#include <algorithm>
#include <iostream>


class MESHVertex;
class MESHEdge;
class MESHTriangle;

typedef MESHVertex *	vm_pointer;
typedef MESHEdge *		em_pointer;
typedef MESHTriangle *	tm_pointer;


typedef double Real;


/*****************************************************************************/
/*																			 */
/*									VERTEX									 */
/*																			 */
/*****************************************************************************/
class MESHVertex {


public:


	/*****************************************************************************/
	/*                                                                           */
	/*    - Class members														 */
	/*                                                                           */
	/*****************************************************************************/
	Real const x;
	Real const y;

	unsigned index;


	/*****************************************************************************/
	/*                                                                           */
	/*    - Class constructor/destructor										 */
	/*                                                                           */
	/*****************************************************************************/
	MESHVertex(Real const X, Real const Y);
	~MESHVertex();


};


MESHVertex::MESHVertex(Real const X, Real const Y) : x(X), y(Y), index(0) {

};
MESHVertex::~MESHVertex() {

};



/*****************************************************************************/
/*																			 */
/*									EDGE									 */
/*																			 */
/*****************************************************************************/
class MESHEdge {


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
	Real const length() const;

};


MESHEdge::MESHEdge(vm_pointer const & va, vm_pointer const & vb) : a(va), b(vb) {

};
MESHEdge::~MESHEdge() {

	a = NULL;
	b = NULL;

	neighbors[0] = NULL;
	neighbors[1] = NULL;

};

Real const MESHEdge::length() const {

	double const x0 = a->x;
	double const y0 = a->y;

	double const x1 = b->x;
	double const y1 = b->y;

	return sqrt((x1 - x0)*(x1 - x0) + (y1 - y0)*(y1 - y0));

};



/*****************************************************************************/
/*																			 */
/*									TRIANGLE								 */
/*																			 */
/*****************************************************************************/
class MESHTriangle {


public:


	/*****************************************************************************/
	/*                                                                           */
	/*    - Class members														 */
	/*                                                                           */
	/*****************************************************************************/
	vm_pointer vertices[3] = { NULL,NULL,NULL };
	em_pointer edges[3] = { NULL,NULL,NULL };
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

	Real const area() const;


};


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

	if (e == edges[0]) return vertices[0];
	else if (e == edges[1]) return vertices[1];
	else if (e == edges[2]) return vertices[2];

	assert(0 && "Triangle::get_vertex(Edge)");

	return vertices[0];

};
vm_pointer const & MESHTriangle::get_vertex_cw(vm_pointer const & v) const {

	if (vertices[0] == v) return vertices[2];
	else if (vertices[1] == v) return vertices[0];
	else if (vertices[2] == v) return vertices[1];

	assert(0 && "Triangle::vertex_cw(Vertex)");

	return vertices[0];

};
vm_pointer const & MESHTriangle::get_vertex_ccw(vm_pointer const & v) const {

	if (vertices[0] == v) return vertices[1];
	else if (vertices[1] == v) return vertices[2];
	else if (vertices[2] == v) return vertices[0];

	assert(0 && "Triangle::vertex_ccw(Vertex)");

	return vertices[0];

};

unsigned const MESHTriangle::get_edge_index(em_pointer const & e) const {

	if (e == edges[0]) return 0;
	else if (e == edges[1]) return 1;
	else if (e == edges[2]) return 2;

	assert(0 && "Triangle::edge_index(Edge)");

	return -1;

};

Real const MESHTriangle::area() const {

	vm_pointer const va = vertices[0];
	vm_pointer const vb = vertices[1];
	vm_pointer const vc = vertices[2];

	return 0.5*((va->x - vc->x) * (vb->y - vc->y) - (va->y - vc->y) * (vb->x - vc->x));

};