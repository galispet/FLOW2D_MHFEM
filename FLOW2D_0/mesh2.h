#pragma once


#pragma once


#include "triangulation.h"
#include "enumerators.h"
#include "mesh_primitives.h"
#include "shapes.h"

#include <vector>
#include <iostream>
#include <fstream>

#include <cmath>


typedef Vertex *		v_pointer;
typedef Edge *			e_pointer;
typedef Triangle *		t_pointer;

typedef MESHVertex *	vm_pointer;
typedef MESHEdge *		em_pointer;
typedef MESHTriangle *	tm_pointer;



bool t_compare_x(Triangle * const t1, Triangle * const t2);
bool t_compare_y(Triangle * const t1, Triangle * const t2);
bool v_compare_x(Vertex * const v1, Vertex * const v2);
bool v_compare_y(Vertex * const v1, Vertex * const v2);
bool e_compare_x(Edge * const e1, Edge * const e2);
bool e_compare_y(Edge * const e1, Edge * const e2);


bool marker_compare_neumann(Edge * const e1, Edge * const e2);
bool marker_compare_dirichlet(Edge * const e1, Edge * const e2);



class Mesh2 {


	std::vector<vm_pointer> vertices;
	std::vector<em_pointer> edges;
	std::vector<tm_pointer> triangles;


	unsigned NumberOf_Vertices	= 0;
	unsigned NumberOf_Edges		= 0;
	unsigned NumberOf_Triangles	= 0;

	unsigned NumberOf_Dirichlet	= 0;
	unsigned NumberOf_Neumann	= 0;


	void clear();


	template < GEOMETRIC_KERNEL GK >
	void allocate_new_primitives(Triangulation<GK> & triangulation);
	template < GEOMETRIC_KERNEL GK >
	void set_adjacencies(Triangulation<GK> & triangulation);

public:



	template < GEOMETRIC_KERNEL GK>
	Mesh2(Triangulation<GK> & triangulation);
	~Mesh2();


	vm_pointer & get_vertex(unsigned const i)	const;
	em_pointer & get_edge(unsigned const i)		const;
	tm_pointer & get_triangle(unsigned const i) const;


	unsigned const get_number_of_vertices()		const;
	unsigned const get_number_of_edges()		const;
	unsigned const get_number_of_triangles()	const;

	unsigned const get_number_of_dirichlet_edges()	const;
	unsigned const get_number_of_neumann_edges()	const;

	void export_vertices(std::ofstream & stream)					const;
	void export_edges(std::ofstream & stream)						const;
	void export_triangles(std::ofstream & stream)					const;
	void export_triangle(std::ofstream & stream, t_pointer const t) const;

};



template < GEOMETRIC_KERNEL GK>
Mesh2::Mesh2(Triangulation<GK> & triangulation) {


	NumberOf_Vertices	= triangulation.get_number_of_vertices();
	NumberOf_Edges		= triangulation.get_number_of_edges();
	NumberOf_Triangles	= triangulation.get_number_of_triangles();

	NumberOf_Dirichlet	= triangulation.get_num_dirichlet_edges();
	NumberOf_Neumann	= triangulation.get_num_neumann_edges();


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
	std::stable_sort(vertices.begin(), vertices.end(), v_compare_x);
	std::stable_sort(vertices.begin(), vertices.end(), v_compare_y);

	std::stable_sort(triangles.begin(), triangles.end(), t_compare_x);
	std::stable_sort(triangles.begin(), triangles.end(), t_compare_y);


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
	std::stable_sort(edges.begin(), edges.end(), marker_compare_dirichlet);
	std::stable_sort(edges.begin(), edges.begin() + num_dirichlet, e_compare_x);
	std::stable_sort(edges.begin(), edges.begin() + num_dirichlet, e_compare_y);

	std::stable_sort(edges.begin() + num_dirichlet, edges.end(), marker_compare_neumann);
	std::stable_sort(edges.begin() + num_dirichlet, edges.begin() + num_neumann + num_dirichlet, e_compare_x);
	std::stable_sort(edges.begin() + num_dirichlet, edges.begin() + num_neumann + num_dirichlet, e_compare_y);
	
	std::stable_sort(edges.begin() + num_neumann + num_dirichlet, edges.end(), e_compare_x);
	std::stable_sort(edges.begin() + num_neumann + num_dirichlet, edges.end(), e_compare_y);

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
	for (unsigned e = 0; e < num_edges; e++) {

		v_pointer const va = edges[e]->a;
		v_pointer const vb = edges[e]->b;

		unsigned const vi = va->index;
		unsigned const vj = vb->index;

		if (vj < vi)
			std::swap(va, vb);

	}



	std::cout << std::endl;
	std::cout << "/*****************************************************************************/" << std::endl;
	std::cout << "/*																		   */" << std::endl;
	std::cout << "/*									 MESH								   */" << std::endl;
	std::cout << "/*																		   */" << std::endl;
	std::cout << "/*****************************************************************************/" << std::endl;
	std::cout << std::endl;
	std::cout << "\nNumber of vertices : ";
	std::cout << vertices.size() << endl;
	std::cout << "\nNumber of edges : ";
	std::cout << edges.size() << endl;
	std::cout << "\nNumber of triangles : ";
	std::cout << triangles.size() << endl;
	std::cout << "Number of Dirichlet edges : ";
	std::cout << std::count_if(edges.begin(), edges.end(), [](auto e) {return e->marker == E_MARKER::DIRICHLET; }) << std::endl;
	std::cout << "Number of Neumann edges : ";
	std::cout << std::count_if(edges.begin(), edges.end(), [](auto e) {return e->marker == E_MARKER::NEUMANN; }) << std::endl;
	std::cout << std::endl;
	std::cout << "/*****************************************************************************/" << std::endl;
	std::cout << "/*																		   */" << std::endl;
	std::cout << "/*									     								   */" << std::endl;
	std::cout << "/*																		   */" << std::endl;
	std::cout << "/*****************************************************************************/" << std::endl;


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

		new_vertex->index	= v->index;

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

		new_edge->index		= e->index;
		new_edge->marker	= e->marker;

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

		new_triangle->index		= t->index;
		new_triangle->marker	= t->marker;

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
	for (size_t i = 0; i < num_edges; i++) {


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
	for (size_t i = 0; i < num_triangles; i++) {


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


