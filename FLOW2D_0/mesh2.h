#pragma once


#pragma once


#include "triangulation.h"
#include "enumerators.h"
#include "shapes.h"
#include "mesh_primitives.h"

#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>




bool tm_compare_x(tm_pointer const t1, tm_pointer const t2);
bool tm_compare_y(tm_pointer const t1, tm_pointer const t2);
bool vm_compare_x(vm_pointer const v1, vm_pointer const v2);
bool vm_compare_y(vm_pointer const v1, vm_pointer const v2);
bool em_compare_x(em_pointer const e1, em_pointer const e2);
bool em_compare_y(em_pointer const e1, em_pointer const e2);

bool mmarker_compare_neumann(em_pointer const e1, em_pointer const e2);
bool mmarker_compare_dirichlet(em_pointer const e1, em_pointer const e2);




class Mesh2 {


	std::vector<vm_pointer> vertices;
	std::vector<em_pointer> edges;
	std::vector<tm_pointer> triangles;


	unsigned NumberOf_Vertices = 0;
	unsigned NumberOf_Edges = 0;
	unsigned NumberOf_Triangles = 0;

	unsigned NumberOf_Dirichlet = 0;
	unsigned NumberOf_Neumann = 0;


	void clear();


	template < GEOMETRIC_KERNEL GK >
	void allocate_new_primitives(Triangulation<GK> & triangulation);
	template < GEOMETRIC_KERNEL GK >
	void set_adjacencies(Triangulation<GK> & triangulation);

public:



	template < GEOMETRIC_KERNEL GK>
	Mesh2(Triangulation<GK> & triangulation);
	~Mesh2();


	vm_pointer const & get_vertex(unsigned const i)		const;
	em_pointer const & get_edge(unsigned const i)		const;
	tm_pointer const & get_triangle(unsigned const i)	const;


	unsigned const get_number_of_vertices()		const;
	unsigned const get_number_of_edges()		const;
	unsigned const get_number_of_triangles()	const;

	unsigned const get_number_of_dirichlet_edges()	const;
	unsigned const get_number_of_neumann_edges()	const;

	void export_vertices(std::ofstream & stream)	const;
	void export_edges(std::ofstream & stream)		const;
	void export_triangles(std::ofstream & stream)	const;

};



template < GEOMETRIC_KERNEL GK>
Mesh2::Mesh2(Triangulation<GK> & triangulation) {


	NumberOf_Vertices = triangulation.get_number_of_vertices();
	NumberOf_Edges = triangulation.get_number_of_edges();
	NumberOf_Triangles = triangulation.get_number_of_triangles();

	NumberOf_Dirichlet = triangulation.get_num_dirichlet_edges();
	NumberOf_Neumann = triangulation.get_num_neumann_edges();


	/*****************************************************************************/
	/*                                                                           */
	/*    - Alocate new MESH vertices, edges, triangles							 */
	/*                                                                           */
	/*****************************************************************************/
	allocate_new_primitives(triangulation);


	/*****************************************************************************/
	/*                                                                           */
	/*    - Set neighboring triangles to triangles and edges					 */
	/*                                                                           */
	/*****************************************************************************/
	set_adjacencies(triangulation);


	/*****************************************************************************/
	/*                                                                           */
	/*    - Sort vertices and triangles by their increasing coordinates			 */
	/*                                                                           */
	/*****************************************************************************/
	std::stable_sort(vertices.begin(), vertices.end(), vm_compare_x);
	std::stable_sort(vertices.begin(), vertices.end(), vm_compare_y);

	std::stable_sort(triangles.begin(), triangles.end(), tm_compare_x);
	std::stable_sort(triangles.begin(), triangles.end(), tm_compare_y);


	/*****************************************************************************/
	/*                                                                           */
	/*    - Sort edges by their markers	and then by x,y coordinates				 */
	/*			: 1. Dirichlet													 */
	/*			: 2. Neumann													 */
	/*			: 3. None														 */
	/*                                                                           */
	/*    - Order of the edges is reversed, so that Dirichlet are last			 */
	/*		and Neumann edges are last-but-one (because of the matrix assembly)	 */
	/*                                                                           */
	/*****************************************************************************/
	std::stable_sort(edges.begin(), edges.end(), mmarker_compare_dirichlet);
	std::stable_sort(edges.begin(), edges.begin() + NumberOf_Dirichlet, em_compare_x);
	std::stable_sort(edges.begin(), edges.begin() + NumberOf_Dirichlet, em_compare_y);

	std::stable_sort(edges.begin() + NumberOf_Dirichlet, edges.end(), mmarker_compare_neumann);
	std::stable_sort(edges.begin() + NumberOf_Dirichlet, edges.begin() + NumberOf_Neumann + NumberOf_Dirichlet, em_compare_x);
	std::stable_sort(edges.begin() + NumberOf_Dirichlet, edges.begin() + NumberOf_Neumann + NumberOf_Dirichlet, em_compare_y);

	std::stable_sort(edges.begin() + NumberOf_Neumann + NumberOf_Dirichlet, edges.end(), em_compare_x);
	std::stable_sort(edges.begin() + NumberOf_Neumann + NumberOf_Dirichlet, edges.end(), em_compare_y);

	std::reverse(edges.begin(), edges.end());


	/*****************************************************************************/
	/*                                                                           */
	/*    - Re-indexing geometric primitives									 */
	/*                                                                           */
	/*****************************************************************************/
	for (size_t i = 0; i < edges.size(); i++)
		edges[i]->index = i;

	for (size_t i = 0; i < vertices.size(); i++)
		vertices[i]->index = i;

	for (size_t i = 0; i < triangles.size(); i++)
		triangles[i]->index = i;


	/*****************************************************************************/
	/*                                                                           */
	/*    - Sort edges end points with respect to vertices's index				 */
	/*		in ascending order													 */
	/*                                                                           */
	/*****************************************************************************/
	for (unsigned e = 0; e < NumberOf_Edges; e++) {

		vm_pointer const va = edges[e]->a;
		vm_pointer const vb = edges[e]->b;

		unsigned const vi = va->index;
		unsigned const vj = vb->index;

		if (vj < vi)
			std::swap(edges[e]->a, edges[e]->b);

	}



	std::cout << std::endl;
	std::cout << "/*****************/" << std::endl;
	std::cout << "/*               */" << std::endl;
	std::cout << "/*      MESH     */" << std::endl;
	std::cout << "/*               */" << std::endl;
	std::cout << "/*****************/" << std::endl;
	std::cout << std::endl;
	std::cout << "No. vertices	: ";
	std::cout << vertices.size() << endl;
	std::cout << "No. edges	: ";
	std::cout << edges.size() << endl;
	std::cout << "No. triangles	: ";
	std::cout << triangles.size() << endl << endl;
	std::cout << "No. Dirichlet edges	: ";
	std::cout << std::count_if(edges.begin(), edges.end(), [](auto e) {return e->marker == E_MARKER::DIRICHLET; }) << std::endl;
	std::cout << "No. Neumann edges	: ";
	std::cout << std::count_if(edges.begin(), edges.end(), [](auto e) {return e->marker == E_MARKER::NEUMANN; }) << std::endl;
	std::cout << std::endl;


};
Mesh2::~Mesh2() {

	clear();

};


template < GEOMETRIC_KERNEL GK >
void Mesh2::allocate_new_primitives(Triangulation<GK> & triangulation) {



	/*****************************************************************************/
	/*                                                                           */
	/*    - Deep Copy of vertices												 */
	/*                                                                           */
	/*****************************************************************************/
	for (size_t i = 0; i < NumberOf_Vertices; i++) {


		v_pointer const v = triangulation.vertices_tri[i];

		double const x = v->x;
		double const y = v->y;

		vm_pointer const new_vertex = new MESHVertex(x, y);

		new_vertex->index = v->index;

		vertices.push_back(new_vertex);

	}


	/*****************************************************************************/
	/*                                                                           */
	/*    - Deep Copy of edges													 */
	/*                                                                           */
	/*****************************************************************************/
	for (size_t i = 0; i < NumberOf_Edges; i++) {


		e_pointer const e = triangulation.edges_tri[i];

		vm_pointer const a = vertices[e->a->index];
		vm_pointer const b = vertices[e->b->index];

		em_pointer const new_edge = new MESHEdge(a, b);

		new_edge->index = e->index;
		new_edge->marker = e->marker;

		edges.push_back(new_edge);

	}


	/*****************************************************************************/
	/*                                                                           */
	/*    - Deep Copy of triangles												 */
	/*                                                                           */
	/*****************************************************************************/
	for (size_t i = 0; i < NumberOf_Triangles; i++) {


		t_pointer const t = triangulation.triangles_tri[i];

		v_pointer const A = t->vertices[0];
		v_pointer const B = t->vertices[1];
		v_pointer const C = t->vertices[2];

		vm_pointer const a = vertices[A->index];
		vm_pointer const b = vertices[B->index];
		vm_pointer const c = vertices[C->index];

		tm_pointer const new_triangle = new MESHTriangle(a, b, c);

		new_triangle->edges[0] = edges[t->edges[0]->index];
		new_triangle->edges[1] = edges[t->edges[1]->index];
		new_triangle->edges[2] = edges[t->edges[2]->index];

		new_triangle->index = t->index;

		triangles.push_back(new_triangle);

	}

};
template < GEOMETRIC_KERNEL GK >
void Mesh2::set_adjacencies(Triangulation<GK> & triangulation) {


	/*****************************************************************************/
	/*                                                                           */
	/*    - Set adjacent triangles to edges										 */
	/*                                                                           */
	/*****************************************************************************/
	for (size_t i = 0; i < NumberOf_Edges; i++) {


		e_pointer const e = triangulation.edges_tri[i];

		t_pointer const n0 = e->neighbors[0];
		t_pointer const n1 = e->neighbors[1];

		em_pointer const f = edges[i];

		if (n0)	f->neighbors[0] = triangles[n0->index];
		else	f->neighbors[0] = NULL;

		if (n1)	f->neighbors[1] = triangles[n1->index];
		else	f->neighbors[1] = NULL;

	}


	/*****************************************************************************/
	/*                                                                           */
	/*    - Set adjacent triangles to triangles									 */
	/*                                                                           */
	/*****************************************************************************/
	for (size_t i = 0; i < NumberOf_Triangles; i++) {


		t_pointer const t = triangulation.triangles_tri[i];

		t_pointer const n0 = t->neighbors[0];
		t_pointer const n1 = t->neighbors[1];
		t_pointer const n2 = t->neighbors[2];

		tm_pointer const k = triangles[i];

		if (n0)	k->neighbors[0] = triangles[n0->index];
		else	k->neighbors[0] = NULL;

		if (n1)	k->neighbors[1] = triangles[n1->index];
		else	k->neighbors[1] = NULL;

		if (n2)	k->neighbors[2] = triangles[n2->index];
		else	k->neighbors[2] = NULL;

	}

};



vm_pointer const & Mesh2::get_vertex(unsigned const i) const {

	return vertices[i];

}
em_pointer const & Mesh2::get_edge(unsigned const i) const {

	return edges[i];

}
tm_pointer const & Mesh2::get_triangle(unsigned const i) const {

	return triangles[i];

}



unsigned const Mesh2::get_number_of_vertices() const {

	return NumberOf_Vertices;

};
unsigned const Mesh2::get_number_of_edges() const {

	return NumberOf_Edges;

};
unsigned const Mesh2::get_number_of_triangles() const {

	return NumberOf_Triangles;

};


unsigned const Mesh2::get_number_of_dirichlet_edges() const {

	return NumberOf_Dirichlet;

};
unsigned const Mesh2::get_number_of_neumann_edges() const {

	return NumberOf_Neumann;

};


void Mesh2::export_vertices(std::ofstream & stream) const {


	for (size_t i = 0; i < vertices.size(); i++) {

		double const v_x = vertices[i]->x;
		double const v_y = vertices[i]->y;

		stream << v_x << "  " << v_y << std::endl;

	}

};
void Mesh2::export_edges(std::ofstream & stream) const {


	/*****************************************************************************/
	/*                                                                           */
	/*    - Gnuplot														         */
	/*                                                                           */
	/*****************************************************************************/
	/*for (size_t i = 0; i < edges.size(); i++) {

		double const v0_x = edges[i]->a->x;
		double const v0_y = edges[i]->a->y;

		double const v1_x = edges[i]->b->x;
		double const v1_y = edges[i]->b->y;

		stream << v0_x << "  " << v0_y << std::endl;
		stream << v1_x << "  " << v1_y << std::endl << std::endl;

	}*/

	/*****************************************************************************/
	/*                                                                           */
	/*    - Matlab														         */
	/*                                                                           */
	/*****************************************************************************/
	for (size_t i = 0; i < edges.size(); i++) {


		double const v0_x = edges[i]->a->x;
		double const v0_y = edges[i]->a->y;

		double const v1_x = edges[i]->b->x;
		double const v1_y = edges[i]->b->y;

		stream << v0_x << "  " << v0_y << "  " << v1_x << "  " << v1_y << std::endl;

	}

};
void Mesh2::export_triangles(std::ofstream & stream) const {


	for (size_t i = 0; i < triangles.size(); i++) {

		tm_pointer t = triangles[i];

		double const v0_x = t->vertices[0]->x;
		double const v0_y = t->vertices[0]->y;

		double const v1_x = t->vertices[1]->x;
		double const v1_y = t->vertices[1]->y;

		double const v2_x = t->vertices[2]->x;
		double const v2_y = t->vertices[2]->y;

		stream << v0_x << "  " << v0_y << std::endl;
		stream << v1_x << "  " << v1_y << std::endl;
		stream << v2_x << "  " << v2_y << std::endl;
		stream << v0_x << "  " << v0_y << std::endl << std::endl;

	}

};


void Mesh2::clear() {


	for (size_t i = 0; i < vertices.size(); i++)
		delete vertices[i];

	for (size_t i = 0; i < edges.size(); i++)
		delete edges[i];

	for (size_t i = 0; i < triangles.size(); i++)
		delete triangles[i];

	vertices.clear();
	edges.clear();
	triangles.clear();

};



bool tm_compare_x(tm_pointer const t1, tm_pointer const t2) {


	double const t1_min = std::fmin(std::fmin(t1->vertices[0]->x, t1->vertices[1]->x), t1->vertices[2]->x);
	double const t2_min = std::fmin(std::fmin(t2->vertices[0]->x, t2->vertices[1]->x), t2->vertices[2]->x);

	return t1_min < t2_min;

};
bool tm_compare_y(tm_pointer const t1, tm_pointer const t2) {


	double const t1_min = std::fmin(std::fmin(t1->vertices[0]->y, t1->vertices[1]->y), t1->vertices[2]->y);
	double const t2_min = std::fmin(std::fmin(t2->vertices[0]->y, t2->vertices[1]->y), t2->vertices[2]->y);

	return t1_min < t2_min;

};
bool vm_compare_x(vm_pointer const v1, vm_pointer const v2) {


	double const v1_min = v1->x;
	double const v2_min = v2->x;

	return v1_min < v2_min;

};
bool vm_compare_y(vm_pointer const v1, vm_pointer const v2) {


	double const v1_min = v1->y;
	double const v2_min = v2->y;

	return v1_min < v2_min;

};
bool em_compare_x(em_pointer const e1, em_pointer const e2) {


	double const e1_min = std::fmin(e1->a->x, e1->b->x);
	double const e2_min = std::fmin(e2->a->x, e2->b->x);

	return e1_min < e2_min;

};
bool em_compare_y(em_pointer const e1, em_pointer const e2) {


	double const e1_min = std::fmin(e1->a->y, e1->b->y);
	double const e2_min = std::fmin(e2->a->y, e2->b->y);

	return e1_min < e2_min;

};

bool mmarker_compare_neumann(em_pointer const e1, em_pointer const e2) {


	E_MARKER const e1_m = e1->marker;
	E_MARKER const e2_m = e2->marker;

	if (e1_m == E_MARKER::NEUMANN && e2_m != E_MARKER::NEUMANN)
		return true;

	return false;

};
bool mmarker_compare_dirichlet(em_pointer const e1, em_pointer const e2) {


	E_MARKER const e1_m = e1->marker;
	E_MARKER const e2_m = e2->marker;

	if (e1_m == E_MARKER::DIRICHLET && e2_m != E_MARKER::DIRICHLET)
		return true;

	return false;

};