#pragma once



#include "shapes.h"
#include "geometric.h"
#include "enumerators.h"
#include "PSLG.h"

#include <iostream>
#include <vector>
#include <list>

//#include <cmath>

#include <fstream>

using std::cout;
using std::endl;


//struct Kernel {
//
//	GEOMETRIC_KERNEL	GK;
//	CHECK_KERNEL		CK;
//
//	double EPSILON_EQUALITY = 1.e-8;
//	double EPSILON_ORIENTATION = 1.e-8;
//	double EPSILON_IN_CIRCUMCIRCLE = 1.e-7;
//	double EPSILON_ON_EDGE = 1.e-7;
//
//};



struct Handle {

	v_pointer v_handle = NULL;
	e_pointer e_handle = NULL;
	t_pointer t_handle = NULL;

};




template <GEOMETRIC_KERNEL GK>
class Triangulation{

	friend class Mesh;
	friend class Mesh2;

private:


	std::vector<v_pointer> vertices_tri;
	std::vector<e_pointer> edges_tri;
	std::vector<t_pointer> triangles_tri;

	v_pointer st_v0 = NULL;
	v_pointer st_v1 = NULL;
	v_pointer st_v2 = NULL;

	e_pointer st_e0 = NULL;
	e_pointer st_e1 = NULL;
	e_pointer st_e2 = NULL;


	// Angle and area bounds for refinement algorithms
	double Angle_Bound;
	double Area_Bound;
	double Minimal_Edge_Length;


	// Additonal container for constrainted edges (segments and subsegments)
	std::vector<e_pointer> segments_tri;

	// Container for additonal vertices resulting from inserting constraints (segments)	by intersection method. Try MAKE new vertices positioned in between two ending points on the edge -> better quality of triangles
	std::vector<v_pointer> constraints_new_vertices;

	// Triangles labled as Outside. They are to be deleted
	std::vector<t_pointer> triangles_outside;

	// Queue of bad triangles in refinement algorithm
	std::list<t_pointer> bad_triangles;

	// Queue of bad triangles in refinement algorithm
	std::list<e_pointer> encroached_segments;

	/*****************************************************************************/
	/*                                                                           */
	/*    - Each vertex of inputed PSLG is given a degree of how many            */
	/*      segments / subsegments is connected to this vertex                   */
	/*    - Used when handling small segments angles with the method of          */
	/*      concentric circular shells                                           */
	/*                                                                           */
	/*****************************************************************************/
	std::vector<unsigned> vertices_degrees_of_subsegments;

	// Geometric functions	:  Template parametr 'GK' provides 'EXACT' or 'INEXACT' computation
	Predicates<GK> const predicates;



	/*****************************************************************************/
	/*                                                                           */
	/*    - Construction and destruction of 'super triangle'                     */
	/*                                                                           */
	/*****************************************************************************/
	void super_triangle_construct(PlanarStraightLineGraph & input_pslg);
	void super_triangle_remove();



	/*****************************************************************************/
	/*                                                                           */
	/*    - marker triangles condemned to destruction !m! 666                    */
	/*                                                                           */
	/*****************************************************************************/
	void mark_triangle(t_pointer const t, int const i_prev);
	void remove_outside_triangles_and_edges();



	/*****************************************************************************/
	/*																			 */
	/*    - Inserts incrementaly single vertex into triangulation                */
	/*																			 */
	/*****************************************************************************/
	v_pointer insert_vertex(v_pointer const v, t_pointer t = NULL);



	/*****************************************************************************/
	/*                                                                           */
	/*    - Point location routine												 */
	/*	  - 'Handle' contains pointer to either: Vertex, Edge, Triangle			 */
	/*		 where the inserted point is located								 */
	/*																			 */
	/*****************************************************************************/
	LOCATION locate_vertex(Handle & location_result, v_pointer const p, t_pointer initial_triangle);



	/*****************************************************************************/
	/*                                                                           */
	/*    - Routines which create new triangles based on where the               */
	/*		inserted point is located											 */
	/*																			 */
	/*****************************************************************************/
	t_pointer insert_vertex_into_triangle(t_pointer const t, v_pointer const p);
	t_pointer insert_vertex_into_edge(e_pointer const e, v_pointer const p, e_pointer &subsegment1, e_pointer &subsegment2);



	/*****************************************************************************/
	/*                                                                           */
	/*    - Legalization of edges routine										 */
	/*	  - Rotate two neighboring triangles CW in order to 'flip' edge			 */
	/*	    which triangles share		  									     */
	/*																			 */
	/*****************************************************************************/
	void rotate_triangle_pair(t_pointer const t, v_pointer const v, t_pointer const ot, v_pointer const ov);
	bool legalize(t_pointer const t, v_pointer const v);



	/*****************************************************************************/
	/*                                                                           */
	/*    - Point localization algorithms                                        */
	/*																			 */
	/*****************************************************************************/
	t_pointer locate_triangle_straight_walk(t_pointer t, v_pointer const p);
	t_pointer locate_triangle_lawson_walk(t_pointer t, v_pointer const p);
	t_pointer locate_triangle_lawson_walk_remembering(t_pointer t, v_pointer const p);

	t_pointer locate_triangle_barycentric(t_pointer t, v_pointer const p);
	t_pointer locate_triangle_lawson_walk_stochastic(t_pointer t, v_pointer const p);
	t_pointer locate_triangle_lawson_walk_remembering_stochastic(t_pointer t, v_pointer const p);
	t_pointer locate_triangle_lawson_walk_fast_remembering(t_pointer t, v_pointer const p, unsigned const n);



	/*****************************************************************************/
	/*                                                                           */
	/*    - Finds an starting vertex 'a' from which the constrained edge is		 */
	/*		heading																 */
	/*																			 */
	/*****************************************************************************/
	t_pointer edge_location(v_pointer const a, v_pointer const b, E_MARKER marker);




	/*****************************************************************************/
	/*                                                                           */
	/*    - Input constrained edge												 */
	/*																			 */
	/*****************************************************************************/
	void insert_constraint_intersection(v_pointer const a, v_pointer const b, E_MARKER marker);



	/*****************************************************************************/
	/*                                                                           */
	/*	  - Frees allocated memory and erase from vector						 */
	/*                                                                           */
	/*    - 'remove' : removes triangle / edge / vertex							 */
	/*				   from triangulation data structure						 */
	/*    - 'clear' :  Clearing all data										 */
	/*                                                                           */
	/*****************************************************************************/
	template<typename T>
	void remove(std::vector<T> & vec, T t);
	template<typename T>
	void replace(std::vector<T> & vec, T t);
	void clear();





public:


	/*****************************************************************************/
	/*                                                                           */
	/*    - Constructor and destructor					                         */
	/*																			 */
	/*****************************************************************************/
	Triangulation(PlanarStraightLineGraph & input_pslg, std::vector<Vertex> seeds = NULL);
	~Triangulation();


	/*****************************************************************************/
	/*                                                                           */
	/*    - Ruppert's Delaunay refinement algorithm		                         */
	/*																			 */
	/*	  - Template parameter priority is order of triangles in				 */
	/*		bad triangle queue													 */
	/*			: BEST  = triangle with lowest ratio radius/edge goes first		 */
	/*			: WORST = triangle with biggest ratio radius/edge goes first	 */
	/*																			 */
	/*	  - 'split_encroached_segments' lookup and split all encroached			 */
	/*		subsegments.														 */
	/*			: Return true if any of subsegments has been split				 */
	/*			: Return false if none of subsegments has been split			 */
	/*																			 */
	/*																			 */
	/*****************************************************************************/
	template<REFINEMENT_PRIORITY priority>
	void refinement_ruppert(double const angle_bound, double const area_max = INFINITY);

	bool get_encroached_segments(v_pointer const p, std::vector<e_pointer> & encroached_segments);
	void split_encroached_segment(e_pointer const e, v_pointer const p);


	void split_encroached_segment(e_pointer const e, v_pointer const p, double const how_far);
	void concentric_shell_split(e_pointer const e, v_pointer const p);

	void get_bad_triangles_ruppert();
	void enqueue_bad_triangle(t_pointer const t);


	bool get_all_encroached_segments();



	/*****************************************************************************/
	/*                                                                           */
	/*	 - Chew's first Delaunay refinement algorithm							 */
	/*			: Will split all triangles, which circumradius is greater		 */
	/*			  then given parameter h_min (shortest edge in the mesh)		 */
	/*			: i.e. the bound is B = 1										 */
	/*			: For given parametr parameter h_min, defined segments			 */
	/*			  are firstly divided into subsegments of the lenght  			 */
	/*			  [ h , sqrt(3)*h ]												 */
	/*																			 */
	/*****************************************************************************/
	t_pointer get_skinny_triangle_chew1st(double &x_c, double &y_c);
	void refinement_chews1st(double const h_min);
	void divide_segments(double const h_min);



	/*****************************************************************************/
	/*                                                                           */
	/*    - Returns number of triangulation primitives							 */
	/*																			 */
	/*****************************************************************************/
	unsigned get_number_of_vertices() const;
	unsigned get_number_of_edges() const;
	unsigned get_number_of_triangles() const;



	/*****************************************************************************/
	/*                                                                           */
	/*    - Returns number of edges with 'NEUMANN' marker						 */
	/*    - Returns number of edges with 'DIRICHLET' marker						 */
	/*																			 */
	/*****************************************************************************/
	unsigned get_num_neumann_edges() const;
	unsigned get_num_dirichlet_edges() const;



	/*****************************************************************************/
	/*                                                                           */
	/*    - Write primitive's coordinates to a text-file						 */
	/*																			 */
	/*****************************************************************************/
	void export_vertices(std::ofstream & stream) const;
	void export_edges(std::ofstream & stream) const;
	void export_triangles(std::ofstream & stream) const;







public:


	//v_pointer insert_vertex(Vertex v) {
	//
	//	Vertex * const new_vertex = new Vertex(v.x, v.y);
	//
	//	input_vertex(new_vertex, NULL);
	//
	//	return new_vertex;
	//
	//};

	//v_pointer insert_vertex(v_pointer const v) {
	//
	//	return input_vertex(v, NULL);
	//
	//};

	/*****************************************************************************/
	/*																			 */
	/*    - Inserts vertices contained in vector into triangulation using		 */
	/*	    'insert_vertex'	on each vertex										 */
	/*																			 */
	/*****************************************************************************/

	//template<class input_iterator>
	//unsigned insert_vertex(input_iterator first, input_iterator last) {
	//
	//
	//	unsigned n = 0;
	//
	//	for (auto it = first; it != last; it++) {
	//
	//		Vertex const v = (*it);
	//
	//		double const x = v.x;
	//		double const y = v.y;
	//
	//		Vertex * const new_vertex = new Vertex(x, y);
	//
	//		input_vertex(new_vertex);
	//
	//		n++;
	//
	//	}
	//
	//	return n;
	//
	//};

	// Also implement 'insert_vertex(input_iterator_handles first, input_iterator_handles last)()', so we can do it with vertices already inserted -> Again !! make sure to implement safe Vertex_handle class, so user can't modifie or input anything

private:

	

public:




	/*****************************************************************************/
	/*                                                                           */
	/*    - Getters for the containers of objects of triangulation primitives	 */
	/*																			 */
	/*****************************************************************************/
	//std::vector<Triangle > get_triangles();
	//std::vector<Edge > get_edges();
	//std::vector<Vertex > get_new_vertices();


	/*****************************************************************************/
	/*                                                                           */
	/*    - Getters for the containers of handles to triangulation primitives	 */
	/*																			 */
	/*****************************************************************************/
	// Don't forget to secure access for the user from changing anything about these primitives !!!!
	//std::vector<Vertex_handle> get_vertices_handles();
	//std::vector<Edge_handle> get_edges_handles();
	//std::vector<Triangle_handle> get_triangles_handles();
	//std::vector<Vertex_handle> get_new_vertices_handles();


	/*****************************************************************************/
	/*                                                                           */
	/*    - Getters for triangulaion primitives									 */
	/*																			 */
	/*****************************************************************************/




public:

	unsigned num_deleted_vertices = 0;

	//void make_individual_hole(Vertex * const seed);
	//void make_individual_hole(Vertex seed);

	//template<class input_iterator>
	//void make_holes(input_iterator first_seed, input_iterator last_seed);

	//template<INSERT_CONSTRAINT algorithm, E_MARKER marker>
	//void insert_individual_constraint(Vertex * const a, Vertex * const b) {

	//	switch (algorithm) {
	//
	//	case INSERT_CONSTRAINT::INTERSECTION:
	//
	//		insert_constraint_intersection<marker>(a, b, a->adjacent_triangle);
	//		break;
	//
	//	case INSERT_CONSTRAINT::HALVING:
	//
	//		insert_constraint_halving<marker>(a, b, a->adjacent_triangle);
	//		break;
	//
	//	case INSERT_CONSTRAINT::FLIP:
	//
	//		insert_constraint_flip<marker>(a, b, a->adjacent_triangle);
	//		break;
	//
	//	case INSERT_CONSTRAINT::ENFORCE:
	//
	//		insert_constraint_enforce<marker>(a, b, a->adjacent_triangle);
	//		break;
	//
	//	}
	//
	//};

	//template<INSERT_CONSTRAINT algorithm, POLYGON_TYPE type, E_MARKER marker, class input_iterator>
	//unsigned insert_sequential_constraint(input_iterator first, input_iterator last) {
	//
	//
	//	// Need to create this temporary container, so there can be inserted constraints (vertices of constraints must be already in triangulation. That doesn't hold for first <-> last)
	//	std::vector<Vertex * > temp;
	//
	//
	//	unsigned n = 0;
	//
	//	for (auto it = first; it != last; it++) {
	//
	//
	//		Vertex const v = (*it);
	//
	//		double const x = v.x;
	//		double const y = v.y;
	//
	//		Vertex * const new_vertex = new Vertex(x, y);
	//
	//		Vertex * const possible_duplicate = insert_vertex(new_vertex);
	//
	//		// In the case, that a vertex of a vertices sequence is almost equal to some vertex already present in triangulation 
	//		if (possible_duplicate != new_vertex) {
	//
	//			temp.push_back(possible_duplicate);
	//			delete new_vertex;
	//
	//		}
	//
	//		else
	//			temp.push_back(new_vertex);
	//
	//		n++;
	//
	//	}
	//
	//
	//	unsigned const nv = n;
	//
	//	for (size_t i = 0; i < nv - 1; i++) {
	//
	//
	//		Vertex * const a = temp[i];
	//		Vertex * const b = temp[i + 1];
	//
	//		a->is_constrained = true;
	//		b->is_constrained = true;
	//
	//		insert_individual_constraint<algorithm, marker>(a, b);
	//
	//	}
	//
	//	//if (type == POLYGON_TYPE::CLOSED)
	//	if (type != POLYGON_TYPE::OPEN)
	//		insert_individual_constraint<algorithm, marker>(temp[nv - 1], temp[0]);
	//
	//
	//	return n;
	//
	//};


	/*****************************************************************************/
	/*                                                                           */
	/*    - Laplacian smoothing													 */
	/*																			 */
	/*****************************************************************************/

public:

	//template<SMOOTHING algorithm>
	//void apply_smoothing(unsigned const N) {


	//	switch (algorithm) {

	//	case SMOOTHING::LAPLACIAN:

	//		laplacian_smoothing(N);
	//		break;

	//	}

	//};


	/*****************************************************************************/
	/*                                                                           */
	/*    - Check for delaunay property of vertex's adjacent triangles 			 */
	/*		and legalize them if needed 										 */
	/*																			 */
	/*	  - Input is vector iterator. If not provided, then all vertices are	 */
	/*		checked																 */
	/*																			 */
	/*****************************************************************************/
	/*
	template<class input_iterator_handle>
	void legalize_vertices(input_iterator_handle first = vertices.begin(), input_iterator_handle last = vertices.end()) {


		unsigned n = 0;

		for (auto it = first; it != last; it++) {


			Vertex * const v = (*it);

			std::vector<Triangle *> adjacent_triangles;
			std::vector<Vertex *> adjacent_vertices;

			size_t const nv = get_adjacent_triangles_and_vertices(v, adjacent_triangles, adjacent_vertices);
			size_t const nt = adjacent_triangles.size();


			for (size_t i = 0; i < nt; i++) {

				Triangle_handle const t = adjacent_triangles[i];

				if (legalize(t, v))
					n++;

			}
		}

		cout << "\nNumber of additional flips: " << n << endl;

	};
	void legalize_vertices() {


		unsigned n = 0;

		size_t const N = vertices.size();


		for (size_t i = 0; i < N; i++) {


			Vertex * const v = vertices[i];

			std::vector<Triangle *> adjacent_triangles;
			std::vector<Vertex *> adjacent_vertices;

			if (!v->adjacent_triangle) {

				cout << "Possibly invalid triangulation (Solo constrained edges or vertices which are not part of any triangle) !" << endl;
				continue;

			}

			size_t const nv = get_adjacent_triangles_and_vertices(v, adjacent_triangles, adjacent_vertices);
			size_t const nt = adjacent_triangles.size();


			for (size_t i = 0; i < nt; i++) {

				Triangle_handle const t = adjacent_triangles[i];

				if (legalize(t, v))
					n++;

			}
		}

		//cout << "\nNumber of additional flips: " << n << endl;

	};
	*/

private:



	/*****************************************************************************/
	/*                                                                           */
	/*    - Additional data structures for easier implementation				 */
	/*																			 */
	/*	  : 'constrained_vertices' : If a vertex is included in constrained edge */
	/*	  : 'is_on_convex_hull' : If a vertex is on the convex hull				 */
	/*																			 */
	/*	  : 'adjacent_triangle_from_vertex' : Pointer to any vertex's adjacent   */
	/*										  triangle							 */
	/*																			 */
	/*****************************************************************************/


	//void laplacian_smoothing(unsigned const N);

	//unsigned get_adjacent_triangles_and_vertices(Vertex * const v, std::vector<Triangle *> & triangles, std::vector<Vertex *> & vertices);





	/*****************************************************************************/
	/*                                                                           */
	/*    - Delaunay's refinement algorithm										 */
	/*																			 */
	/*																			 */
	/*****************************************************************************/

//	double shortest_edge_length() {
//
//
//		double shortest_edge_squared = INFINITY;
//
//		for (size_t i = 0; i < triangles.size(); i++) {
//
//
//			Triangle * const t = triangles[i];
//
//			Vertex * const a = t->vertices[0];
//			Vertex * const b = t->vertices[1];
//			Vertex * const c = t->vertices[2];
//
//			double const a_x = a->x;
//			double const a_y = a->y;
//
//			double const b_x = b->x;
//			double const b_y = b->y;
//
//			double const c_x = c->x;
//			double const c_y = c->y;
//
//			double const A = b_x - a_x;
//			double const B = b_y - a_y;
//			double const C = c_x - a_x;
//			double const D = c_y - a_y;
//
//			double const H = c_x - b_x;
//			double const I = c_y - b_y;
//
//
//			double const AB = A * A + B * B;
//			double const BC = H * H + I * I;
//			double const CA = C * C + D * D;
//
//			double const this_min = std::fmin(AB, std::fmin(BC, CA));
//
//			if (this_min < shortest_edge_squared)
//				shortest_edge_squared = this_min;
//
//		}
//
//		return shortest_edge_squared;
//
//	}

};







/*****************************************************************************/
/*                                                                           */
/*    - Constructor and destructor					                         */
/*																			 */
/*****************************************************************************/
template <GEOMETRIC_KERNEL GK>
Triangulation<GK>::Triangulation(PlanarStraightLineGraph & input_pslg, std::vector<Vertex> seeds) {


	// I need at least on triangle in the triangulation. Super triangle is enclosing all vertices and segments given in PSLG
	super_triangle_construct(input_pslg);


	size_t const nv = input_pslg.vertices_pslg.size();
	size_t const ns = input_pslg.segments_pslg.size();


	// This will propably help after a change the write-in structures of primitives
	//vertices_tri.reserve(1000 * nv);
	//edges_tri.reserve(3 * 1000 * nv - 6);
	//triangles_tri.reserve(2 * 1000 * nv - 5);


	vertices_degrees_of_subsegments.resize(nv);


	// Deep copy of vertices(independent memory allocation)
	for (size_t i = 0; i < nv; i++) {


		v_pointer const v = input_pslg.vertices_pslg[i];

		v_pointer const new_vertex = new Vertex(v->x, v->y);

		new_vertex->marker = v->marker;

		insert_vertex(new_vertex);


		vertices_degrees_of_subsegments[new_vertex->index] = 0;

	}


	// Insert constrained segments defined by the user
	for (size_t i = 0; i < ns; i++) {


		// 0. index of first vertex 1. index of second vertex 2. edge_marker 3. constraint_method
		unsigned * const constraint = input_pslg.segments_pslg[i];

		v_pointer const a = vertices_tri[constraint[0]];
		v_pointer const b = vertices_tri[constraint[1]];

		// Constraint endpoints are automatically marked as constrained
		a->marker = V_MARKER::CONSTRAINED;
		b->marker = V_MARKER::CONSTRAINED;

		E_MARKER const marker = (E_MARKER)constraint[2];
		INSERT_CONSTRAINT const insert_method = (INSERT_CONSTRAINT)constraint[3];

		insert_constraint_intersection(a, b, marker);


		vertices_degrees_of_subsegments[a->index]++;
		vertices_degrees_of_subsegments[b->index]++;

	}


	// Mark triangles as outside in user-denoted hole using seed vertex
	for (size_t i = 0; i < seeds.size(); i++)
		mark_triangle(locate_triangle_straight_walk(triangles_tri.back(), &seeds[i]), -1);


	//Mark triangles outside, every triangle neighboring to super triangle and of course in concavities
	mark_triangle(st_v0->adjacent_triangle, -1);


	//Remove triangles and their edges which are marked outside
	remove_outside_triangles_and_edges();


	//Remove super_vertices and super_edges
	super_triangle_remove();


	// After deletetion of triangles, edges and vertices, we have to reindex it. Indexing is important for efficient for write-in / remove into data structers of primitives
	for (size_t i = 0; i < vertices_tri.size(); i++)
		vertices_tri[i]->index = i;

	for (size_t i = 0; i < edges_tri.size(); i++)
		edges_tri[i]->index = i;

	for (size_t i = 0; i < triangles_tri.size(); i++)
		triangles_tri[i]->index = i;


};
template <GEOMETRIC_KERNEL GK>
Triangulation<GK>::~Triangulation() {

	clear();

};




/*****************************************************************************/
/*                                                                           */
/*    - Construction and destruction of 'super triangle'                     */
/*                                                                           */
/*****************************************************************************/
template <GEOMETRIC_KERNEL GK>
void Triangulation<GK>::super_triangle_construct(PlanarStraightLineGraph & input_pslg) {


	// Get maximum and minimum coordinate on x,y axis
	double minX = input_pslg.vertices_pslg[0]->x;
	double minY = input_pslg.vertices_pslg[0]->y;
	double maxX = minX;
	double maxY = minY;

	size_t const nv = input_pslg.vertices_pslg.size();

	for (size_t i = 0; i < nv; i++) {

		v_pointer const v = input_pslg.vertices_pslg[i];

		if (v->x < minX) minX = v->x;
		if (v->y < minY) minY = v->y;
		if (v->x > maxX) maxX = v->x;
		if (v->y > maxY) maxY = v->y;

	}

	double const dx = maxX - minX;
	double const dy = maxY - minY;
	double const ox = (maxX + minX) / 2.0;
	double const oy = (maxY + minY) / 2.0;

	double const M = std::max(dx, dy);

	// Construct super triangle
	st_v0 = new Vertex(ox - 20 * M, oy - M);
	st_v1 = new Vertex(ox + 20 * M, oy - M);
	st_v2 = new Vertex(ox, oy + 20 * M);

	//st_v0 = new Vertex(-3 * M, -3 * M);
	//st_v1 = new Vertex(3 * M, 0.0);
	//st_v2 = new Vertex(0.0, 3 * M);

	st_v0->marker = V_MARKER::CONSTRAINED;
	st_v1->marker = V_MARKER::CONSTRAINED;
	st_v2->marker = V_MARKER::CONSTRAINED;

	Triangle * const super_triangle = new Triangle(st_v0, st_v1, st_v2);

	// Construct edges of super triangle
	st_e0 = new Edge(st_v1, st_v2);
	st_e1 = new Edge(st_v2, st_v0);
	st_e2 = new Edge(st_v0, st_v1);

	// Assign edges to super triangle
	super_triangle->set_edge(0) = st_e0;
	super_triangle->set_edge(1) = st_e1;
	super_triangle->set_edge(2) = st_e2;

	// Assign neighboring triangles to edges
	st_e0->set_neighbor(0) = super_triangle;
	st_e0->set_neighbor(1) = NULL;

	st_e1->set_neighbor(0) = super_triangle;
	st_e1->set_neighbor(1) = NULL;

	st_e2->set_neighbor(0) = super_triangle;
	st_e2->set_neighbor(1) = NULL;

	st_e0->is_constrained = true;
	st_e1->is_constrained = true;
	st_e2->is_constrained = true;


	// Add super triangle to initial triangulation
	triangles_tri.push_back(super_triangle);

	// Fisrt triangle in the triangulation
	super_triangle->index = 0;

};
template <GEOMETRIC_KERNEL GK>
void Triangulation<GK>::super_triangle_remove() {


	delete st_v0;
	delete st_v1;
	delete st_v2;

	delete st_e0;
	delete st_e1;
	delete st_e2;

	st_v0 = NULL;
	st_v1 = NULL;
	st_v2 = NULL;

	st_e0 = NULL;
	st_e1 = NULL;
	st_e2 = NULL;

};




/*****************************************************************************/
/*                                                                           */
/*    - marker triangles condemned to destruction !m! 666                    */
/*                                                                           */
/*****************************************************************************/
template <GEOMETRIC_KERNEL GK>
void Triangulation<GK>::mark_triangle(t_pointer const t, int const i_prev) {


	if (!t) {

		cout << "Skipping seed. It would lead to invalid triangulation." << endl;

		return;

	}

	t->marker = T_MARKER::OUTSIDE;

	for (unsigned i = 0; i < 3; i++) {


		// I came here from this triangle. It is already marked and processed
		if (i == i_prev)
			continue;

		// Set new direction of depth-search
		t_pointer const neighbor = t->neighbors[i];

		// If we are on the edge of super triangle, there are no neighbors
		if (!neighbor)
			continue;

		if (neighbor->marker == T_MARKER::OUTSIDE || t->edges[i]->is_constrained)
			continue;


		// Recurrent depth-search
		mark_triangle(neighbor, neighbor->get_neighbor_index(t));

	}

	// Insert this triangle into stack of triangles to be deleted (marked as outside)
	triangles_outside.push_back(t);

};
template <GEOMETRIC_KERNEL GK>
void Triangulation<GK>::remove_outside_triangles_and_edges() {


	std::vector<v_pointer> vertices_to_remove;

	for (size_t i = 0; i < triangles_outside.size(); i++) {


		t_pointer const t = triangles_outside[i];


		// Delete edges inside hole, which are not constrained
		for (unsigned j = 0; j < 3; j++) {


			e_pointer const e = t->edges[j];

			// If edge is already deleted from other iteration 
			if (!e)
				continue;

			t_pointer const neighbor = t->neighbors[j];

			v_pointer const a = e->a;
			v_pointer const b = e->b;


			if (!e->is_constrained) {


				// The edge 'e' is being destroyed. We need to set the pointer NULL to this edge
				neighbor->edges[neighbor->get_edge_index(e)] = NULL;

				remove(edges_tri, e);
				delete e;

				// Vertex inside a hole is pushed into 'vertices_to_remove' container, so that is uniquely present in the container
				if (a->marker != V_MARKER::CONSTRAINED) {


					bool already_in = false;

					for (size_t k = 0; k < vertices_to_remove.size(); k++) {

						if (a == vertices_to_remove[k]) {

							already_in = true;
							break;

						}
					}

					if (!already_in)
						vertices_to_remove.push_back(a);

				}

				if (b->marker != V_MARKER::CONSTRAINED) {

					bool already_in = false;

					for (size_t k = 0; k < vertices_to_remove.size(); k++) {

						if (b == vertices_to_remove[k]) {

							already_in = true;
							break;

						}
					}

					if (!already_in)
						vertices_to_remove.push_back(b);

				}

			}
			else if (neighbor) {

				neighbor->neighbors[neighbor->get_edge_index(e)] = NULL;

				e->neighbors[0] = neighbor;
				e->neighbors[1] = NULL;

				a->adjacent_triangle = neighbor;
				b->adjacent_triangle = neighbor;

			}

		}

		remove(triangles_tri, t);
		delete t;

	}

	// Test on some examples what is faster !!
	std::sort(vertices_to_remove.begin(), vertices_to_remove.end());
	std::unique(vertices_to_remove.begin(), vertices_to_remove.end());

	// Delete vertices inside hole
	for (size_t i = 0; i < vertices_to_remove.size(); i++) {

		v_pointer const v = vertices_to_remove[i];

		remove(vertices_tri, v);
		delete v;

		num_deleted_vertices++;

	}

	// Re-index vertices
	if (!vertices_to_remove.empty()) {

		for (size_t i = 0; i < vertices_tri.size(); i++)
			vertices_tri[i]->index = i;

	}


	vertices_to_remove.clear();
	triangles_outside.clear();

};


/*****************************************************************************/
/*																			 */
/*    - Inserts incrementaly single vertex into triangulation                */
/*																			 */
/*****************************************************************************/
template <GEOMETRIC_KERNEL GK>
v_pointer Triangulation<GK>::insert_vertex(v_pointer const v, t_pointer t) {


	Handle rh;

	switch (locate_vertex(rh, v, t)) {

	case LOCATION::IN_TRIANGLE: {

		insert_vertex_into_triangle(rh.t_handle, v);
		break;

	}
	case LOCATION::ON_EDGE: {

		e_pointer dummy = NULL;

		insert_vertex_into_edge(rh.e_handle, v, dummy, dummy);
		break;

	}
	case LOCATION::ON_VERTEX: {

		std::cout << "Duplicate vertex ignored. Returning pointer to existing vertex." << std::endl;
		return rh.v_handle;

	}
							  // Get rid of this somehow in the future. Also this is connected to point-location -> there are checks if the triangles is NULL : NONESENSE !! can't happen
	case LOCATION::NOT_VISIBLE: {

		std::cout << "NOT visible." << std::endl;
		return NULL;

	}
	}

	v->index = vertices_tri.size();

	vertices_tri.push_back(v);

	return v;

};




/*****************************************************************************/
/*                                                                           */
/*	  - Location of the vertex												 */
/*                                                                           */
/*	  - Handle contains pointer to either: Vertex, Edge, Triangle			 */
/*		 where the inserted point is located								 */
/*																			 */
/*****************************************************************************/
template <GEOMETRIC_KERNEL GK>
LOCATION Triangulation<GK>::locate_vertex(Handle & location_result, v_pointer const p, t_pointer initial_triangle) {


	t_pointer k = initial_triangle == NULL ? triangles_tri.back() : initial_triangle;


	t_pointer t = locate_triangle_straight_walk(k, p);
	//t_pointer t = locate_triangle_lawson_walk(k, p);
	//t_pointer t = locate_triangle_lawson_walk_remembering(k, p);

	//t_pointer t = locate_triangle_barycentric(k, p);
	//t_pointer t = locate_triangle_lawson_walk_stochastic(k, p);
	//t_pointer t = locate_triangle_lawson_walk_remembering_stochastic(k, p);
	//t_pointer t = locate_triangle_lawson_walk_fast_remembering(k, p, 50);

	if (!t)
		return LOCATION::NOT_VISIBLE;


	// For the validity :) -> that's why i need to make point location more robust
	if (!predicates.in_triangle_robust(t, p)) {

		std::cout << "Inserted vertex is not in the resulting Triangle. Using robust predicates ... ";

		for (unsigned i = 0; i < 3; i++) {

			Triangle * const k = t->neighbors[i];

			if (predicates.in_triangle_robust(k, p)) {

				t = k;

				cout << "Ok" << endl;
				break;

			}
			else if (i == 2) {

				cout << " : Not Ok. Trying to locate again." << endl;
				return locate_vertex(location_result, p, k);

			}
		}
	}


	CHECK_KERNEL CK = CHECK_KERNEL::CHECK_DUPLICATE;

	if (CK == CHECK_KERNEL::CHECK_DUPLICATE) {

		v_pointer const v0 = t->vertices[0];
		v_pointer const v1 = t->vertices[1];
		v_pointer const v2 = t->vertices[2];

		if (p->is_almost_equal(v0)) {

			location_result.v_handle = v0;
			return LOCATION::ON_VERTEX;

		}
		else if (p->is_almost_equal(v1)) {

			location_result.v_handle = v1;
			return LOCATION::ON_VERTEX;

		}
		else if (p->is_almost_equal(v2)) {

			location_result.v_handle = v2;
			return LOCATION::ON_VERTEX;

		}
	}

	// If used robust, then it might lead to creation of degenrate triangle
	bool const b0 = predicates.on_edge_fast(t->edges[0], p);
	bool const b1 = predicates.on_edge_fast(t->edges[1], p);
	bool const b2 = predicates.on_edge_fast(t->edges[2], p);


	if (b0) {

		location_result.e_handle = t->edges[0];
		return LOCATION::ON_EDGE;

	}
	else if (b1) {

		location_result.e_handle = t->edges[1];
		return LOCATION::ON_EDGE;

	}
	else if (b2) {

		location_result.e_handle = t->edges[2];
		return LOCATION::ON_EDGE;

	}

	location_result.t_handle = t;

	return LOCATION::IN_TRIANGLE;

};




/*****************************************************************************/
/*                                                                           */
/*    - Routines which create new triangles based on where the               */
/*		inserted point is located											 */
/*																			 */
/*****************************************************************************/
template<GEOMETRIC_KERNEL GK>
t_pointer Triangulation<GK>::insert_vertex_into_triangle(t_pointer const t, v_pointer const p) {


	Vertex * const v0 = t->vertices[0];
	Vertex * const v1 = t->vertices[1];
	Vertex * const v2 = t->vertices[2];

	Triangle * const n0 = t->neighbors[0];
	Triangle * const n1 = t->neighbors[1];
	Triangle * const n2 = t->neighbors[2];

	Triangle * const t0 = new Triangle(p, v1, v2);
	Triangle * const t1 = new Triangle(p, v2, v0);
	Triangle * const t2 = new Triangle(p, v0, v1);


	// If any of these triangles is / are (both possible in English) degenerate, insert it into appropriate edge instead
	if (predicates.orientation_fast(p, v1, v2) == 0.0) {


		std::cout << "It happened !" << std::endl;

		e_pointer dummy = NULL;

		insert_vertex_into_edge(t->edges[0], p, dummy, dummy);
		return NULL;

	}
	if (predicates.orientation_fast(p, v2, v0) == 0.0) {

		std::cout << "It happened !" << std::endl;

		e_pointer dummy = NULL;

		insert_vertex_into_edge(t->edges[1], p, dummy, dummy);
		return NULL;

	}
	if (predicates.orientation_fast(p, v0, v1) == 0.0) {

		std::cout << "It happened !" << std::endl;

		e_pointer dummy = NULL;

		insert_vertex_into_edge(t->edges[2], p, dummy, dummy);
		return NULL;

	}


	Edge * const e0 = new Edge(p, v0);
	Edge * const e1 = new Edge(p, v1);
	Edge * const e2 = new Edge(p, v2);

	// Assign edges to triangles
	t0->set_edge(0) = t->edges[0];
	t0->set_edge(1) = e2;
	t0->set_edge(2) = e1;

	t1->set_edge(0) = t->edges[1];
	t1->set_edge(1) = e0;
	t1->set_edge(2) = e2;

	t2->set_edge(0) = t->edges[2];
	t2->set_edge(1) = e1;
	t2->set_edge(2) = e0;

	// Assign neighboring triangles to edges
	t->edges[0]->set_neighbor(0) = t0;
	t->edges[0]->set_neighbor(1) = n0;

	t->edges[1]->set_neighbor(0) = t1;
	t->edges[1]->set_neighbor(1) = n1;

	t->edges[2]->set_neighbor(0) = t2;
	t->edges[2]->set_neighbor(1) = n2;


	// Set indeces of new triangles
	unsigned const nt = triangles_tri.size();		// see paper for fastering the remove() ->transformed to erase .... see the paper man

	t0->index = t->index;
	t1->index = nt;
	t2->index = nt + 1;

	// Set indeces of new edges
	unsigned const ne = edges_tri.size();
	e0->index = ne;
	e1->index = ne + 1;
	e2->index = ne + 2;


	// Set pointer to any vertex's 'p' adjacent triangle. Also there is a possibility that vertices of the located triangle 't' have adjacency to already deleted triangle 't', therefore we need to renew these adjacencies.
	p->adjacent_triangle = t0;

	v0->adjacent_triangle = t1;
	v1->adjacent_triangle = t2;
	v2->adjacent_triangle = t0;


	// Erase 't' from bad triangles if present
	bad_triangles.remove(t);

	// Remove and delete splitted triangle consisting of t0 , t1 , t2 from triangulation
	//triangles_tri[t->index] = t0;
	remove(triangles_tri, t);
	delete t;



	// Assign neighboring triangles to edges
	e0->set_neighbor(0) = t2;
	e0->set_neighbor(1) = t1;

	e1->set_neighbor(0) = t2;
	e1->set_neighbor(1) = t0;

	e2->set_neighbor(0) = t0;
	e2->set_neighbor(1) = t1;

	// Assign neighborhood between triangles
	t0->set_neighbor(t1);
	t0->set_neighbor(n0);

	t1->set_neighbor(t2);
	t1->set_neighbor(n1);

	t2->set_neighbor(t0);
	t2->set_neighbor(n2);

	// Add new triangles to triangulation
	triangles_tri.push_back(t0);
	triangles_tri.push_back(t1);
	triangles_tri.push_back(t2);

	// Add new edges to triangulation
	edges_tri.push_back(e0);
	edges_tri.push_back(e1);
	edges_tri.push_back(e2);

	// Perform edge-flip algorithm to legalize non-Delaunay edges
	legalize(t0, p);
	legalize(t1, p);
	legalize(t2, p);


	// Check if newly triangles are bad or not. Refinement algorithm use this container for bad triangles.
	if (!bad_triangles.empty() && t0->is_bad(Angle_Bound, Area_Bound)) enqueue_bad_triangle(t0);
	if (!bad_triangles.empty() && t1->is_bad(Angle_Bound, Area_Bound)) enqueue_bad_triangle(t1);
	if (!bad_triangles.empty() && t2->is_bad(Angle_Bound, Area_Bound)) enqueue_bad_triangle(t2);


	return t0;

};
template<GEOMETRIC_KERNEL GK>
t_pointer Triangulation<GK>::insert_vertex_into_edge(e_pointer e, v_pointer p, e_pointer &subsegment1, e_pointer &subsegment2) {


	t_pointer const t0 = e->neighbors[0];
	t_pointer const t1 = e->neighbors[1];

	bool const flag = e->is_constrained;
	E_MARKER const marker = e->marker;


	v_pointer v0 = NULL;
	v_pointer v1 = NULL;

	v_pointer w0 = NULL;
	v_pointer w1 = NULL;

	e_pointer e0 = NULL;
	e_pointer e1 = NULL;
	e_pointer e2 = NULL;
	e_pointer e3 = NULL;

	t_pointer n0 = NULL;
	t_pointer n1 = NULL;
	t_pointer n2 = NULL;
	t_pointer n3 = NULL;

	t_pointer k0 = NULL;
	t_pointer k1 = NULL;
	t_pointer k2 = NULL;
	t_pointer k3 = NULL;

	if (t0) {

		unsigned const v0_index = t0->get_edge_index(e);

		v0 = t0->vertices[v0_index];

		w0 = t0->get_vertex_ccw(v0);
		w1 = t0->get_vertex_cw(v0);

		n0 = t0->get_neighbor_ccw(v0);
		n1 = t0->get_neighbor_cw(v0);

		k0 = new Triangle(p, v0, w0);
		k1 = new Triangle(p, w1, v0);

		e1 = new Edge(p, v0);

	}
	if (t1) {

		unsigned const v1_index = t1->get_edge_index(e);

		v1 = t1->vertices[v1_index];

		w0 = t1->get_vertex_cw(v1);
		w1 = t1->get_vertex_ccw(v1);

		n2 = t1->get_neighbor_ccw(v1);
		n3 = t1->get_neighbor_cw(v1);

		k2 = new Triangle(p, v1, w1);
		k3 = new Triangle(p, w0, v1);

		e3 = new Edge(p, v1);

	}


	e0 = new Edge(p, w0);
	e2 = new Edge(p, w1);


	if (t0) {

		// Assign edges to triangles 'k0' and 'k1' created from triangle 't0'
		k0->set_edge(0) = t0->get_edge_ccw(v0);
		k0->set_edge(1) = e0;
		k0->set_edge(2) = e1;

		k1->set_edge(0) = t0->get_edge_cw(v0);
		k1->set_edge(1) = e1;
		k1->set_edge(2) = e2;

	}
	if (t1) {

		// Assign edges to triangles 'k2' and 'k3' created from triangle 't1'
		k2->set_edge(0) = t1->get_edge_ccw(v1);
		k2->set_edge(1) = e2;
		k2->set_edge(2) = e3;

		k3->set_edge(0) = t1->get_edge_cw(v1);
		k3->set_edge(1) = e3;
		k3->set_edge(2) = e0;

	}

	// Copy constraint flag of splitted edges to edges which are formed from it
	e0->is_constrained = flag;
	e2->is_constrained = flag;

	e0->marker = marker;
	e2->marker = marker;


	// Assign neighboring triangles to edges
	e0->set_neighbor(1) = k0;
	e2->set_neighbor(0) = k1;
	e2->set_neighbor(1) = k2;
	e0->set_neighbor(0) = k3;


	if (t0) {

		// Assign neighboring triangles to edges
		e1->set_neighbor(0) = k0;
		e1->set_neighbor(1) = k1;

		k0->edges[0]->set_neighbor(0) = k0;
		k0->edges[0]->set_neighbor(1) = n0;

		k1->edges[0]->set_neighbor(0) = k1;
		k1->edges[0]->set_neighbor(1) = n1;

		// Assign neighborhoods of triangles
		k0->set_neighbor(k1);
		k0->set_neighbor(n0);

		k1->set_neighbor(n1);

	}
	if (t1) {

		// Assign neighboring triangles to edges
		e3->set_neighbor(0) = k2;
		e3->set_neighbor(1) = k3;

		k2->edges[0]->set_neighbor(0) = k2;
		k2->edges[0]->set_neighbor(1) = n2;

		k3->edges[0]->set_neighbor(0) = k3;
		k3->edges[0]->set_neighbor(1) = n3;

		// Assign neighborhoods of triangles
		k2->set_neighbor(k3);
		k2->set_neighbor(n2);

		k3->set_neighbor(n3);

	}
	if (t0 && t1) {

		// Assign neighborhoods of triangles
		k1->set_neighbor(k2);
		k3->set_neighbor(k0);

	}


	// Set pointer to any vertex's 'p' adjacent triangle. Also there is a possibility that some vertices have pointers to deleted triangle 't'
	if (t0) p->adjacent_triangle = k0;
	else	p->adjacent_triangle = k3;

	if (t0) v0->adjacent_triangle = k1;
	if (t1) v1->adjacent_triangle = k2;

	if (t0) {

		w1->adjacent_triangle = k1;
		w0->adjacent_triangle = k0;

	}
	else {

		w1->adjacent_triangle = k2;
		w0->adjacent_triangle = k3;

	}


	// Set indeces of new edges
	size_t ne = edges_tri.size();

	e0->index = e->index;
	e2->index = ne;

	// Remove and delete splitted edge from the triangulation
	remove(edges_tri, e);
	delete e;

	// Add new edges to triangulation
	edges_tri.push_back(e0);
	edges_tri.push_back(e2);

	//ne++;

	subsegment1 = e0;
	subsegment2 = e2;


	if (t0) {


		// Set indeces of new triangles
		unsigned const nt = triangles_tri.size();		// see paper for fastering the remove() ->transformed to erase .... see the paper man

		k0->index = t0->index;
		k1->index = nt;


		// Erase 't0' from bad triangles if present
		bad_triangles.remove(t0);

		// Remove and delete splitted triangle from the triangulation
		//triangles_tri[t0->index] = k0;
		remove(triangles_tri, t0);
		delete t0;

		// Add new triangles to triangulation
		triangles_tri.push_back(k0);
		triangles_tri.push_back(k1);

		size_t const ne = edges_tri.size();
		e1->index = ne;
		//ne++;

		// Add new edges to triangulation
		edges_tri.push_back(e1);





	}
	if (t1) {

		// Set indeces of new triangles
		unsigned const nt = triangles_tri.size();		// see paper for fastering the remove() ->transformed to erase .... see the paper man

		k2->index = t1->index;
		k3->index = nt;


		// Erase 't1' from bad triangles if present
		bad_triangles.remove(t1);

		// Remove and delete splitted triangle from the triangulation
		//triangles_tri[t1->index] = k2;
		remove(triangles_tri, t1);
		delete t1;

		triangles_tri.push_back(k2);
		triangles_tri.push_back(k3);

		size_t const ne = edges_tri.size();
		e3->index = ne;
		//ne++;

		// Add new edges to triangulation
		edges_tri.push_back(e3);

	}



	// Perform edge-flip algorithm to legalize non-Delaunay edge
	if (t0) {

		legalize(k0, p);
		legalize(k1, p);

		bool const empty_queue = bad_triangles.empty();

		// Check if newly triangles are bad or not. Refinement algorithm use this container for bad triangles.
		if (!empty_queue && k0->is_bad(Angle_Bound, Area_Bound)) enqueue_bad_triangle(k0);
		if (!empty_queue && k1->is_bad(Angle_Bound, Area_Bound)) enqueue_bad_triangle(k1);

	}
	if (t1) {

		legalize(k2, p);
		legalize(k3, p);

		bool const empty_queue = bad_triangles.empty();

		// Check if newly triangles are bad or not. Refinement algorithm use this container for bad triangles.
		if (!empty_queue && k2->is_bad(Angle_Bound, Area_Bound)) enqueue_bad_triangle(k2);
		if (!empty_queue && k3->is_bad(Angle_Bound, Area_Bound)) enqueue_bad_triangle(k3);

	}


	if (t0)
		return k0;
	else
		return k3;


	/*
	Triangle * const t0 = e->neighbor(0);
	Triangle * const t1 = e->neighbor(1);

	if (t0) {


	}
	if(t1)
	int const v0_index = t0->edge_index(e);
	int const v1_index = t1->edge_index(e);

	Vertex * const v0 = t0->vertex(v0_index);
	Vertex * const v1 = t1->vertex(v1_index);

	Vertex * const w0 = t0->vertex_ccw(v0);		// Equivalently : t1->vertex_cw(v1)
	Vertex * const w1 = t0->vertex_cw(v0);		// Equivalently : t1->vertex_ccw(v1)

	Triangle * const n0 = t0->neighbor_ccw(v0);
	Triangle * const n1 = t0->neighbor_cw(v0);
	Triangle * const n2 = t1->neighbor_ccw(v1);
	Triangle * const n3 = t1->neighbor_cw(v1);

	bool const c_flag = e->is_constrained;
	E_MARKER const marker = e->marker;

	Triangle * const k0 = new Triangle(p, v0, w0);
	Triangle * const k1 = new Triangle(p, w1, v0);
	Triangle * const k2 = new Triangle(p, v1, w1);
	Triangle * const k3 = new Triangle(p, w0, v1);

	Edge * const e0 = new Edge(p, w0);
	Edge * const e1 = new Edge(p, v0);
	Edge * const e2 = new Edge(p, w1);
	Edge * const e3 = new Edge(p, v1);

	// Assign edges to triangles
	k0->set_edge(0) = t0->edge_ccw(v0);
	k0->set_edge(1) = e0;
	k0->set_edge(2) = e1;

	k1->set_edge(0) = t0->edge_cw(v0);
	k1->set_edge(1) = e1;
	k1->set_edge(2) = e2;

	k2->set_edge(0) = t1->edge_ccw(v1);
	k2->set_edge(1) = e2;
	k2->set_edge(2) = e3;

	k3->set_edge(0) = t1->edge_cw(v1);
	k3->set_edge(1) = e3;
	k3->set_edge(2) = e0;

	// Copy constraint flag of splitted edges to edges which are formed from it
	e0->is_constrained = c_flag;
	e2->is_constrained = c_flag;

	e0->marker = marker;
	e2->marker = marker;

	// Assign neighboring triangles to edges
	e0->set_neighbor(0) = k3;
	e0->set_neighbor(1) = k0;

	e1->set_neighbor(0) = k0;
	e1->set_neighbor(1) = k1;

	e2->set_neighbor(0) = k1;
	e2->set_neighbor(1) = k2;

	e3->set_neighbor(0) = k2;
	e3->set_neighbor(1) = k3;

	k0->edge(0)->set_neighbor(0) = k0;
	k0->edge(0)->set_neighbor(1) = n0;

	k1->edge(0)->set_neighbor(0) = k1;
	k1->edge(0)->set_neighbor(1) = n1;

	k2->edge(0)->set_neighbor(0) = k2;
	k2->edge(0)->set_neighbor(1) = n2;

	k3->edge(0)->set_neighbor(0) = k3;
	k3->edge(0)->set_neighbor(1) = n3;

	// Assign neighborhoods of triangles
	k0->set_neighbor(k1);
	k0->set_neighbor(n0);

	k1->set_neighbor(k2);
	k1->set_neighbor(n1);

	k2->set_neighbor(k3);
	k2->set_neighbor(n2);

	k3->set_neighbor(k0);
	k3->set_neighbor(n3);


	// Set pointer to any vertex's 'p' adjacent triangle. Also there is a possibility that some vertices have pointers to deleted triangle 't'
	p->adjacent_triangle = k0;

	v0->adjacent_triangle = k1;
	w1->adjacent_triangle = k2;
	v1->adjacent_triangle = k3;
	w0->adjacent_triangle = k0;

	// Remove and delete splitted triangles and edge from the triangulation
	remove(edges, e);
	delete e;

	remove(triangles, t0);
	delete t0;

	remove(triangles, t1);
	delete t1;

	// Add new triangles to triangulation
	triangles.push_back(k0);
	triangles.push_back(k1);
	triangles.push_back(k2);
	triangles.push_back(k3);

	// Add new edges to triangulation
	edges.push_back(e0);
	edges.push_back(e1);
	edges.push_back(e2);
	edges.push_back(e3);


	// Perform edge-flip algorithm to legalize non-Delaunay edge
	legalize(k0, p);
	legalize(k1, p);
	legalize(k2, p);
	legalize(k3, p);

	return k0;
	*/

};




/*****************************************************************************/
/*                                                                           */
/*    - Legalization of edges routines										 */
/*																			 */
/*****************************************************************************/
template<GEOMETRIC_KERNEL GK>
void Triangulation<GK>::rotate_triangle_pair(t_pointer const t, v_pointer const v, t_pointer const ot, v_pointer const ov) {


	// Save neighbors of triangles
	t_pointer const n0 = t->get_neighbor_ccw(v);
	t_pointer const n1 = t->get_neighbor_cw(v);
	t_pointer const n2 = ot->get_neighbor_ccw(ov);
	t_pointer const n3 = ot->get_neighbor_cw(ov);

	// Save edges of triangles
	e_pointer const flip_edge = t->edges[t->get_vertex_index(v)];

	e_pointer const e0 = t->get_edge_ccw(v);
	e_pointer const e1 = t->get_edge_cw(v);
	e_pointer const e2 = ot->get_edge_ccw(ov);
	e_pointer const e3 = ot->get_edge_cw(ov);

	// Flip edge by rotating neighboring triangles clock-wise, redefine edge
	t->rotate_triangle_cw(v, ov);
	ot->rotate_triangle_cw(ov, v);

	flip_edge->set(v, ov);

	// Assign edges to triangles
	t->set_edge(flip_edge);
	t->set_edge(e2);
	t->set_edge(e1);

	ot->set_edge(flip_edge);
	ot->set_edge(e0);
	ot->set_edge(e3);

	// Assign neighboring triangles to edges
	e0->set_neighbor(0) = ot;
	e0->set_neighbor(1) = n0;

	e1->set_neighbor(0) = t;
	e1->set_neighbor(1) = n1;

	e2->set_neighbor(0) = t;
	e2->set_neighbor(1) = n2;

	e3->set_neighbor(0) = ot;
	e3->set_neighbor(1) = n3;

	flip_edge->set_neighbor(0) = t;
	flip_edge->set_neighbor(1) = ot;

	// Assign neighborhood
	t->null_neighbors();
	ot->null_neighbors();

	if (n0) ot->set_neighbor(n0);
	if (n1) t->set_neighbor(n1);
	if (n2) t->set_neighbor(n2);
	if (n3) ot->set_neighbor(n3);

	t->set_neighbor(ot);

};
template<GEOMETRIC_KERNEL GK>
bool Triangulation<GK>::legalize(t_pointer const t, v_pointer const v) {


	unsigned const edgeToFlip = t->get_vertex_index(v);

	if (t->edges[edgeToFlip]->is_constrained)
		return false;

	t_pointer const ot = t->neighbors[edgeToFlip];

	// When hole or super triangle, 'ot' may be NULL
	if (!ot)
		return false;

	v_pointer const ov = ot->get_opposite_vertex(t, v);


	bool const in = predicates.in_circle(ot, v);

	if (in) {

		rotate_triangle_pair(t, v, ot, ov);

		if (!bad_triangles.empty()) {

			bad_triangles.remove(t);
			bad_triangles.remove(ot);

			if (t->is_bad(Angle_Bound, Area_Bound)) enqueue_bad_triangle(t);
			if (ot->is_bad(Angle_Bound, Area_Bound)) enqueue_bad_triangle(ot);

		}


		// Rotating triangles can make vertices's pointers of quadrilateral (except 'v'), point to different triangle, which does not contain respective vertices
		ov->adjacent_triangle = ot;
		t->get_vertex_cw(v)->adjacent_triangle = t;
		ot->get_vertex_cw(ov)->adjacent_triangle = ot;


		legalize(t, v);
		legalize(ot, v);

		return true;

	}

	return false;

};




/*****************************************************************************/
/*                                                                           */
/*    - Point localization algorithms                                        */
/*																			 */
/*****************************************************************************/
template<GEOMETRIC_KERNEL GK>
t_pointer Triangulation<GK>::locate_triangle_straight_walk(t_pointer t, v_pointer const p) {


	v_pointer q = t->vertices[0];
	v_pointer r = t->vertices[1];
	v_pointer l = t->vertices[2];

	v_pointer s = NULL;

	if (predicates.orientation(r, q, p) <= 0.0) {

		while (predicates.orientation(l, q, p) < 0.0) {

			r = l;
			t = t->get_neighbor(q, l);

			if (!t)
				return NULL;

			l = t->get_vertex_but(q, r);

		}
	}
	else {

		do {

			l = r;
			t = t->get_neighbor(q, r);

			if (!t)
				return NULL;

			r = t->get_vertex_but(q, l);

		} while (predicates.orientation(r, q, p) > 0.0);
	}

	while (predicates.orientation(p, r, l) <= 0.0) {

		t = t->get_neighbor(r, l);

		if (!t)
			return NULL;

		s = t->get_vertex_but(r, l);

		if (predicates.orientation(s, q, p) < 0.0)
			r = s;
		else
			l = s;

	}

	return t;

};
template<GEOMETRIC_KERNEL GK>
t_pointer Triangulation<GK>::locate_triangle_lawson_walk(t_pointer t, v_pointer const p) {


	bool found = false;

	v_pointer r = NULL;
	v_pointer l = NULL;

	unsigned const ccw_indeces[4] = { 0 , 1 , 2 , 0 };

	while (!found) {

		found = true;

		for (unsigned i = 0; i < 3; i++) {

			//r = t->vertices[i % 3];
			//l = t->vertices[(i + 1) % 3];

			// May be SLOWER -> check for large triangulation
			r = t->vertices[ccw_indeces[i]];
			l = t->vertices[ccw_indeces[i + 1]];

			if (predicates.orientation(r, l, p) < 0.0) {

				t = t->get_neighbor(l, r);
				found = false;

				break;

			}

		}

	}

	return t;

};
template<GEOMETRIC_KERNEL GK>
t_pointer Triangulation<GK>::locate_triangle_lawson_walk_remembering(t_pointer t, v_pointer const p) {


	t_pointer previous = t;

	bool found = false;

	unsigned const ccw_indeces[4] = { 0 , 1 , 2 , 0 };

	v_pointer r = NULL;
	v_pointer l = NULL;

	while (!found) {

		found = true;

		for (unsigned i = 0; i < 3; i++) {

			//r = t->vertex(i % 3);
			//l = t->vertex((i + 1) % 3);

			r = t->vertices[ccw_indeces[i]];
			l = t->vertices[ccw_indeces[i + 1]];

			if (previous != t->get_neighbor(l, r)) {

				if (predicates.orientation(r, l, p) < 0.0) {

					previous = t;
					t = t->get_neighbor(l, r);
					found = false;

					break;

				}

			}

		}

	}

	return t;

};


template<GEOMETRIC_KERNEL GK>
t_pointer Triangulation<GK>::locate_triangle_barycentric(t_pointer t, v_pointer const p) {


	v_pointer r = NULL;
	v_pointer l = NULL;

	// Computation of barycentric coordinates is optimized -> instead of dividing both 'b0','b1' by 'd' we multiply the 'b2' by 'd' -> much faster

	double c[3];
	unsigned const ccw_indeces[4] = { 0 , 1 , 2 , 0 };


	double d = predicates.orientation(t->vertices[0], t->vertices[1], t->vertices[2]);

	c[0] = predicates.orientation(t->vertices[1], t->vertices[2], p);
	c[1] = predicates.orientation(t->vertices[2], t->vertices[0], p);
	c[2] = d - c[0] - c[1];


	double min = c[0];

	r = t->vertices[1];
	l = t->vertices[2];

	if (c[1] < min) {

		r = t->vertices[2];
		l = t->vertices[0];

		min = c[1];

	}
	if (c[2] < min) {

		r = t->vertices[0];
		l = t->vertices[1];

		min = c[2];

	}

	t_pointer tau = NULL;

	while (min < 0.0) {

		tau = t->get_neighbor(r, l);

		unsigned int i = tau->get_neighbor_index(t);

		t = tau;

		c[ccw_indeces[i]] = -min;

		//c[ccw_indeces[i + 1]] = predicates.orientation(tau->vertices[ccw_indeces[i + 2]], tau->vertices[ccw_indeces[i]], p);	// Check the indeces
		c[(i + 1) % 3] = predicates.orientation(tau->vertices[(i + 2) % 3], tau->vertices[(i + 3) % 3], p);

		d = predicates.orientation(tau->vertices[0], tau->vertices[1], tau->vertices[2]);

		//c[ccw_indeces[i + 2]] = d - c[ccw_indeces[i]] - c[ccw_indeces[i + 1]];	// Check the indeces
		c[(i + 2) % 3] = d - c[i] - c[(i + 1) % 3];

		min = c[0];

		r = tau->vertices[1];
		l = tau->vertices[2];

		if (c[1] < min) {

			r = tau->vertices[2];
			l = tau->vertices[0];

			min = c[1];

		}
		if (c[2] < min) {

			r = tau->vertices[0];
			l = tau->vertices[1];

			min = c[2];

		}

	}

	return t;

};
template<GEOMETRIC_KERNEL GK>
t_pointer Triangulation<GK>::locate_triangle_lawson_walk_stochastic(t_pointer t, v_pointer const p) {


	t_pointer previous = t;

	bool found = false;

	//unsigned const ccw_indeces[5] = { 0 , 1 , 2 , 0 , 1 };

	v_pointer r = NULL;
	v_pointer l = NULL;

	while (!found) {

		found = true;

		unsigned const k = rand() % 3;

		for (unsigned i = k; i < k + 3; i++) {

			r = t->vertices[i % 3];
			l = t->vertices[(i + 1) % 3];

			//r = t->vertices[ccw_indeces[i]];
			//l = t->vertices[ccw_indeces[i + 1]];

			if (predicates.orientation(r, l, p) < 0.0) {

				previous = t;
				t = t->get_neighbor(l, r);
				found = false;

				break;

			}


		}

	}

	return t;

};
template<GEOMETRIC_KERNEL GK>
t_pointer Triangulation<GK>::locate_triangle_lawson_walk_remembering_stochastic(t_pointer t, v_pointer const p) {


	t_pointer previous = t;

	bool found = false;

	//unsigned const ccw_indeces[4] = { 0 , 1 , 2 , 0 };

	v_pointer r = NULL;
	v_pointer l = NULL;

	while (!found) {

		found = true;

		unsigned const k = rand() % 3;

		for (unsigned i = k; i < k + 3; i++) {

			r = t->vertices[i % 3];
			l = t->vertices[(i + 1) % 3];

			//r = t->vertices[ccw_indeces[i]];
			//l = t->vertices[ccw_indeces[i + 1]];

			if (previous != t->get_neighbor(l, r)) {

				if (predicates.orientation(r, l, p) < 0.0) {

					previous = t;
					t = t->get_neighbor(l, r);
					found = false;

					break;

				}


			}

		}

	}

	return t;

};
template<GEOMETRIC_KERNEL GK>
t_pointer Triangulation<GK>::locate_triangle_lawson_walk_fast_remembering(t_pointer t, v_pointer const p, unsigned const n) {


	t_pointer psi = t;

	t_pointer orig = t;

	v_pointer r = NULL;
	v_pointer l = NULL;

	unsigned const ccw_indeces[4] = { 0 , 1 , 2 , 0 };

	for (unsigned k = 0; k < n; k++) {

		unsigned j = 0;

		for (unsigned i = 0; i < 3; i++) {

			//l = t->vertex(i % 3);
			//r = t->vertex((i + 1) % 3);

			r = t->vertices[ccw_indeces[i]];
			l = t->vertices[ccw_indeces[i + 1]];

			j = i;

			if (psi != t->get_neighbor(r, l))
				break;

		}

		psi = t;

		if (predicates.orientation(l, r, p) < 0.0)
			t = t->get_neighbor(r, l);

		else {

			l = r;
			r = t->vertices[(j + 2) % 3];

			t = t->get_neighbor(r, l);

		}

		if (!t)
			return locate_triangle_lawson_walk(orig, p);

	}

	return locate_triangle_lawson_walk(t, p);

};




/*****************************************************************************/
/*                                                                           */
/*    - Input constrained edge												 */
/*																			 */
/*****************************************************************************/
template<GEOMETRIC_KERNEL GK>
void Triangulation<GK>::insert_constraint_intersection(v_pointer const a, v_pointer const b, E_MARKER marker) {



	// Search for the triangle 't' which contains 'a' and the constrained edge goes through the edge which is formed by the rest 2 vertices 'vertex_ccw(a)' <-> 'vertex_cw(a)'. If the edge is already there, NULL is returned
	t_pointer const t = edge_location(a, b, marker);


	if (!t)
		return;


	// Constrained edge is situated between these two points
	v_pointer const r = t->get_vertex_ccw(a);
	v_pointer const l = t->get_vertex_cw(a);


	if (on_edge_fast(a, b, r)) {


		e_pointer const f = t->get_edge(a, r);

		// Constrained edge goes through ccw vertex 'r' with respect to 'a', and so, part of constrained edge is already in there
		f->marker = marker;
		f->is_constrained = true;

		// 'f' already exists, therefore its adjacency is already set

		segments_tri.push_back(f);

		// Denote Vertices of constrained edge as 'constrained' -> When smoothing, these vertices won't be displaced
		a->marker = V_MARKER::CONSTRAINED;
		r->marker = V_MARKER::CONSTRAINED;

		// 'f' already exists, therefore its adjacency is already set

		// The rest of the c.edge must start at the vertex 'r' and end in 'b'. Therefore, it is best to lookup for the triangle containg 'r' from the opposite triangle 'ot'
		insert_constraint_intersection(r, b, marker);

		return;

	}
	else if (on_edge_fast(a, b, l)) {


		e_pointer const f = t->get_edge(a, l);

		// Constrained edge goes through cw vertex 'l' with respect to 'a', and so, part of constrained edge is already in there
		f->marker = marker;
		f->is_constrained = true;

		segments_tri.push_back(f);

		// Denote Vertices of constrained edge as 'constrained' -> When smoothing, these vertices won't be displaced
		a->marker = V_MARKER::CONSTRAINED;
		l->marker = V_MARKER::CONSTRAINED;

		// The rest of the c.edge must start at the vertex 'l' and end in 'b'. Therefore, it is best to lookup for the triangle containg 'l' from the opposite triangle 'ot'
		insert_constraint_intersection(l, b, marker);

		return;

	}


	// An edge between 'r' and 'l' vertices
	e_pointer const e = t->get_edge(r, l);


	// Opposite triangle
	t_pointer ot = t->get_opposite_triangle(a);
	v_pointer const s = ot->get_vertex_but(r, l);


	// Edge intersection coordinates
	double ix;
	double iy;


	if (s == b) {


		// Opposite vertex 's' in triangle 'ot' IS the ending vertex 'b' 
		get_edge_intersection(a, b, l, r, ix, iy);

		v_pointer const new_vertex = new Vertex(ix, iy);

		// Denote Vertices of constrained edge as 'constrained' -> When smoothing, these vertices won't be displaced
		new_vertex->index = vertices_tri.size();
		new_vertex->marker = V_MARKER::CONSTRAINED;
		a->marker = V_MARKER::CONSTRAINED;
		s->marker = V_MARKER::CONSTRAINED;

		vertices_tri.push_back(new_vertex);
		constraints_new_vertices.push_back(new_vertex);


		// Insert new_vertex defined by intersection of edges and Get any triangle containing 'new_vertex'

		e_pointer dummy = NULL;

		ot = insert_vertex_into_edge(e, new_vertex, dummy, dummy);


		// Rotary traverse through triangles around 'new_vertex', so that the resulting triangle contains 'a'. So we could denote new subsegments as constrained
		while (true) {

			if (!ot->contains(a))	ot = ot->get_neighbor_ccw(new_vertex);
			else					break;

		}


		e_pointer f = ot->get_edge(new_vertex, a);

		// From the resulting triangle 'ot' we are able to to set the type 'marker' of boundary and to set the part of constrained edge to be constrained
		f->is_constrained = true;
		f->marker = marker;

		//// Setting adjacency to this edge. We know, that 'ot' is on one side. The other side is simply neighbor through 'a' and 'new_vertex'
		//f->neighbors[0] = ot;
		//f->neighbors[1] = ot->get_neighbor(a, new_vertex);

		segments_tri.push_back(f);


		// The same here, but for the ending vertex 'b'
		while (true) {

			if (!ot->contains(s))	ot = ot->get_neighbor_ccw(new_vertex);
			else					break;

		}

		f = ot->get_edge(new_vertex, s);

		f->is_constrained = true;
		f->marker = marker;

		//// Setting adjacency to this edge. We know, that 'ot' is on one side. The other side is simply neighbor through 's' and 'new_vertex'
		//f->neighbors[0] = ot;
		//f->neighbors[1] = ot->get_neighbor(s, new_vertex);

		segments_tri.push_back(f);


		// Now we are done. No need for further introduction of additional vertices
		return;


	}
	else if (on_edge_fast(a, b, s)) {


		get_edge_intersection(a, b, l, r, ix, iy);

		v_pointer const new_vertex = new Vertex(ix, iy);

		new_vertex->index = vertices_tri.size();
		new_vertex->marker = V_MARKER::CONSTRAINED;
		// Denote Vertices of constrained edge as 'constrained' -> When smoothing, these vertices won't be displaced
		a->marker = V_MARKER::CONSTRAINED;
		s->marker = V_MARKER::CONSTRAINED;

		vertices_tri.push_back(new_vertex);
		constraints_new_vertices.push_back(new_vertex);


		e_pointer dummy = NULL;

		ot = insert_vertex_into_edge(e, new_vertex, dummy, dummy);


		while (true) {

			if (!ot->contains(a))	ot = ot->get_neighbor_ccw(new_vertex);
			else					break;

		}

		e_pointer f = ot->get_edge(new_vertex, a);

		f->is_constrained = true;
		f->marker = marker;

		//// Setting adjacency to this edge. We know, that 'ot' is on one side. The other side is simply neighbor through 'a' and 'new_vertex'
		//f->neighbors[0] = ot;
		//f->neighbors[1] = ot->get_neighbor(a, new_vertex);

		segments_tri.push_back(f);


		while (true) {

			if (!ot->contains(s))	ot = ot->get_neighbor_ccw(new_vertex);
			else					break;

		}

		f = ot->get_edge(new_vertex, s);

		f->is_constrained = true;
		f->marker = marker;

		//// Setting adjacency to this edge. We know, that 'ot' is on one side. The other side is simply neighbor through 'a' and 'new_vertex'
		//f->neighbors[0] = ot;
		//f->neighbors[1] = ot->get_neighbor(s, new_vertex);

		segments_tri.push_back(f);


		// Contrary to the last possibility 's == b', the c.edge passes through the opposite vertex 's', but vertex 's' is not an ending point. Therefore further subdivision is needed
		insert_constraint_intersection(s, b, marker);

	}
	else {

		get_edge_intersection(a, b, l, r, ix, iy);

		v_pointer const new_vertex = new Vertex(ix, iy);

		new_vertex->index = vertices_tri.size();
		new_vertex->marker = V_MARKER::CONSTRAINED;
		// Denote Vertices of constrained edge as 'constrained' -> When smoothing, these vertices won't be displaced
		a->marker = V_MARKER::CONSTRAINED;

		vertices_tri.push_back(new_vertex);
		constraints_new_vertices.push_back(new_vertex);


		e_pointer dummy = NULL;

		ot = insert_vertex_into_edge(e, new_vertex, dummy, dummy);


		while (true) {

			if (!ot->contains(a))	ot = ot->get_neighbor_ccw(new_vertex);
			else					break;

		}

		e_pointer const f = ot->get_edge(new_vertex, a);

		f->is_constrained = true;
		f->marker = marker;

		//// Setting adjacency to this edge. We know, that 'ot' is on one side. The other side is simply neighbor through 'a' and 'new_vertex'
		//f->neighbors[0] = ot;
		//f->neighbors[1] = ot->get_neighbor(a, new_vertex);

		segments_tri.push_back(f);

		// This is typical situation, when the c.edge doesn't pass through any vertex but the 'new_vertex'
		insert_constraint_intersection(new_vertex, b, marker);

	}

};




/*****************************************************************************/
/*                                                                           */
/*    - Finds an starting vertex 'a' from which the constrained edge is		 */
/*		heading																 */
/*																			 */
/*****************************************************************************/
template<GEOMETRIC_KERNEL GK>
t_pointer Triangulation<GK>::edge_location(v_pointer const a, v_pointer const b, E_MARKER marker) {


	// Locate triangle containg vertex 'a'  
	//t_pointer t = locate_triangle_straight_walk(triangles_tri.back(), a);
	//t_pointer t = locate_triangle_lawson_walk(triangles_tri.back(), a);
	//t_pointer t = locate_triangle_lawson_walk_remembering(triangles_tri.back(), a);

	// Locate triangle containg vertex 'a'  
	t_pointer t = a->adjacent_triangle;


	// For Debugging purpose :) 
	if (!predicates.in_triangle_robust(t, a))
		std::cout << "Edge location routine : 'a' not in Triangle 't'" << std::endl;


	// Circulate through triangles around the point 'a'	so we can check if the edge exists, or to get the direction to the vertice 'b'
	v_pointer l = t->get_vertex_cw(a);
	v_pointer r = t->get_vertex_ccw(a);


	if (orientation_robust(r, a, b) <= 0.0) {

		while (orientation_robust(l, a, b) < 0.0) {

			r = l;
			t = t->get_neighbor_cw(a);
			l = t->get_vertex_cw(a);

			//t = t->get_neighbor(a, l);
			//l = t->vertex_but(a, r);

		}
	}
	else {

		do {

			l = r;
			t = t->get_neighbor_ccw(a);
			r = t->get_vertex_ccw(a);

			//t = t->neighbor(a, r);
			//R = t->vertex_but(a, l);

		} while (orientation_robust(r, a, b) > 0.0);
	}


	// Now 't' is the starting triangle. Check its edges if there exists edge 'ab'
	if (l == b || r == b) {


		//std::cout << "Edge already exists. Skipping." << std::endl;

		e_pointer const e = t->get_edge(a, b);


		// If the segment is already created some time before, it is already marked as constrained. We have to skip or else this segment would be in the container twice 8-)
		if (e->is_constrained) {

			// If the user (by funny mistake) want two different marks on the same edge, laugh at him
			if (e->marker != marker)
				std::cout << "Ambiguous segment marker : [ " << e->a->index << " , " << e->b->index << " ]. Maintaing primal Mark." << std::endl;

			return NULL;

		}

		e->is_constrained = true;
		e->marker = marker;

		segments_tri.push_back(e);

		return NULL;

	}

	return t;

};




/*****************************************************************************/
/*                                                                           */
/*    - 'remove' : removes triangle / edge / vertex							 */
/*				   from triangulation data structure						 */
/*    - 'clear'  : completely clear all data								 */
/*                                                                           */
/*****************************************************************************/
template<GEOMETRIC_KERNEL GK> template<typename T>
void Triangulation<GK>::remove(std::vector<T> & vec, T t) {


	// It is much faster to lookup for an object 't' from the back
	for (auto it = vec.rbegin(); it != vec.rend(); it++) {


		// If an iterator is equal object 't' it is erased from vector and the search is finished
		if (*it == t) {

			vec.erase(it.base() - 1);
			break;

		}
	}

};
template<GEOMETRIC_KERNEL GK> template<typename T>
void Triangulation<GK>::replace(std::vector<T> & vec, T t) {

	vec[t->index] = t;

};
template<GEOMETRIC_KERNEL GK>
void Triangulation<GK>::clear() {


	v_pointer v = NULL;
	e_pointer e = NULL;
	t_pointer t = NULL;


	// First, free allocated memory
	for (size_t i = 0; i < vertices_tri.size(); i++) {

		v = vertices_tri[i];
		delete v;

	}
	for (size_t i = 0; i < edges_tri.size(); i++) {

		e = edges_tri[i];
		delete e;

	}
	for (size_t i = 0; i < triangles_tri.size(); i++) {

		t = triangles_tri[i];
		delete t;

	}

	// Secodnly, clear the vectors of the rubbish
	vertices_tri.clear();
	edges_tri.clear();
	triangles_tri.clear();

	// These containers contains already deleted primitives. No need for delete theme again (it would only lead to an program failure).
	segments_tri.clear();
	triangles_outside.clear();
	constraints_new_vertices.clear();

};




/*****************************************************************************/
/*                                                                           */
/*    - Returns number of triangulation primitives							 */
/*																			 */
/*****************************************************************************/
template<GEOMETRIC_KERNEL GK>
unsigned Triangulation<GK>::get_number_of_vertices() const {

	return (unsigned)vertices_tri.size();

};
template<GEOMETRIC_KERNEL GK>
unsigned Triangulation<GK>::get_number_of_edges() const {

	return (unsigned)edges_tri.size();

};
template<GEOMETRIC_KERNEL GK>
unsigned Triangulation<GK>::get_number_of_triangles() const {

	return (unsigned)triangles_tri.size();

};




/*****************************************************************************/
/*                                                                           */
/*    - Returns number of edges with 'NEUMANN' marker						 */
/*    - Returns number of edges with 'DIRICHLET' marker						 */
/*																			 */
/*****************************************************************************/
template<GEOMETRIC_KERNEL GK>
unsigned Triangulation<GK>::get_num_neumann_edges() const {

	return (unsigned)std::count_if(edges_tri.begin(), edges_tri.end(), [](e_pointer e) { return e->marker == E_MARKER::NEUMANN; });

};
template<GEOMETRIC_KERNEL GK>
unsigned Triangulation<GK>::get_num_dirichlet_edges() const {

	return (unsigned)std::count_if(edges_tri.begin(), edges_tri.end(), [](e_pointer e) {return e->marker == E_MARKER::DIRICHLET; });

};




/*****************************************************************************/
/*                                                                           */
/*    - Write primitive's coordinates to a text-file						 */
/*																			 */
/*****************************************************************************/
template<GEOMETRIC_KERNEL GK>
void Triangulation<GK>::export_vertices(std::ofstream & stream) const {


	size_t const nv = vertices_tri.size();

	for (size_t i = 0; i < nv; i++) {

		double const v_x = vertices_tri[i]->x;
		double const v_y = vertices_tri[i]->y;

		stream << v_x << "  " << v_y << std::endl;

	}

};
template<GEOMETRIC_KERNEL GK>
void Triangulation<GK>::export_edges(std::ofstream & stream) const {


	size_t const ne = edges_tri.size();


	//// Gnuplot
	for (size_t i = 0; i < ne; i++) {


		e_pointer const e = edges_tri[i];

		double const v0_x = e->a->x;
		double const v0_y = e->a->y;

		double const v1_x = e->b->x;
		double const v1_y = e->b->y;

		stream << v0_x << "  " << v0_y << std::endl;
		stream << v1_x << "  " << v1_y << std::endl << std::endl;

	}

	//// Matlab
	//for (size_t i = 0; i < ne; i++) {

	//	e_pointer const e = edges_tri[i];
	//
	//	double const v0_x = e->a->x;
	//	double const v0_y = e->a->y;
	//
	//	double const v1_x = e->b->x;
	//	double const v1_y = e->b->y;
	//
	//	stream << v0_x << "  " << v0_y << "  " << v1_x << "  " << v1_y << std::endl;
	//
	//}

};
template<GEOMETRIC_KERNEL GK>
void Triangulation<GK>::export_triangles(std::ofstream & stream) const {


	size_t const nt = triangles_tri.size();

	for (size_t i = 0; i < nt; i++) {


		t_pointer const t = triangles_tri[i];

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

	// Matlab
	//for (size_t i = 0; i < nt; i++) {
	//
	//
	//	t_pointer const t = triangles_tri[i];
	//
	//	double const v0_x = t->vertices[0]->x;
	//	double const v0_y = t->vertices[0]->y;
	//
	//	double const v1_x = t->vertices[1]->x;
	//	double const v1_y = t->vertices[1]->y;
	//
	//	double const v2_x = t->vertices[2]->x;
	//	double const v2_y = t->vertices[2]->y;
	//
	//	stream << v0_x << "  " << v0_y << "  " << v1_x << "  " << v1_y << "  " << v2_x << "  " << v2_y<< "  " << v0_x << "  " << v0_y <<  std::endl;
	//
	//}

};




/*****************************************************************************/
/*                                                                           */
/*    - Ruppert's Delaunay refinement algorithm		                         */
/*																			 */
/*	  - Template parameter priority is order of triangles in				 */
/*		bad triangle queue													 */
/*			: BEST  = triangle with lowest ratio radius/edge goes first		 */
/*			: WORST = triangle with biggest ratio radius/edge goes first	 */
/*																			 */
/*	  - 'split_encroached_segments' lookup and split all encroached			 */
/*		subsegments.														 */
/*			: Return true if any of subsegments has been split				 */
/*			: Return false if none of subsegments has been split			 */
/*																			 */
/*																			 */
/*****************************************************************************/

template<GEOMETRIC_KERNEL GK> template<REFINEMENT_PRIORITY priority>
void Triangulation<GK>::refinement_ruppert(double const angle_bound, double const area_max) {


	Angle_Bound = angle_bound;
	Area_Bound = area_max;


	get_bad_triangles_ruppert();

	bool continue_refining = true;
	bool split_segments = true;

	while (true) {


		continue_refining = false;


		// Split every encroached segments / subsegments there are		
		if (split_segments) {


			for (size_t i = 0; i < segments_tri.size(); i++) {


				e_pointer const e = segments_tri[i];

				t_pointer const n1 = e->neighbors[0];
				t_pointer const n2 = e->neighbors[1];

				if (n1) {

					v_pointer const apex = n1->get_vertex(e);


					if (is_encroached(e->a, e->b, apex)) {

						// Probably useless. Termination doesn't depend on that. Find out why
						continue_refining = true;

						split_encroached_segment(e, apex);

						i--;
						continue;

					}
				}

				if (n2) {

					v_pointer const apex = n2->get_vertex(e);

					if (is_encroached(e->a, e->b, apex)) {

						// Probably useless. Termination doesn't depend on that. Find out why
						continue_refining = true;

						split_encroached_segment(e, apex);

						i--;
						continue;

					}
				}
			}
		}




		// Get any skinny triangle and its circumcenter
		t_pointer const bad_triangle = !bad_triangles.empty() ? (priority == REFINEMENT_PRIORITY::BEST ? bad_triangles.front() : bad_triangles.back()) : NULL;


		// If there are no skinny triangles left, we are done
		if (!bad_triangle && !continue_refining)
			break;



		// Circum center coordinates
		double x_c;
		double y_c;

		bad_triangle->circum_center(x_c, y_c);


		// Circumcenter of skinny triangle
		Vertex circum_center(x_c, y_c);

		std::vector<e_pointer> encroached_segments2;


		if (get_encroached_segments(&circum_center, encroached_segments2)) {


			split_segments = true;

			size_t const k = encroached_segments2.size();

			//for (size_t i = 0; i < k; i++)
			//	split_encroached_segment_no_recursion(encroached_segments2[i]);

			for (size_t i = 0; i < encroached_segments2.size(); i++)
				split_encroached_segment(encroached_segments2[i], &circum_center);

		}
		else {

			split_segments = false;

			// We have to begin from given skinny triangle, so if the circumcenter is behind segments in a hole, it wont crash. Otherwise, we can go to new vertex through hole => there are no triangles (problem with walking algorithm)
			//v_pointer const possible_duplicate = insert_vertex(new Vertex(x_c, y_c), skinny_triangle);
			insert_vertex(new Vertex(x_c, y_c), bad_triangle);

		}

	}

	//std::stable_sort(vertices_tri.begin(), vertices_tri.end(), [](v_pointer v, v_pointer p) {return v->index < p->index; });
	//for (size_t i = 0; i < vertices_tri.size() - 1; i++) {

	//	unsigned const index = vertices_tri[i + 1]->index - vertices_tri[i]->index;
	//	//unsigned const index = vertices_tri[i]->index;

	//	cout << index << endl;

	//	if (index != 1)
	//		cout << "shit ---------------------------------------- \n----------------------------------------" << endl;

	//}


	//std::stable_sort(edges_tri.begin(), edges_tri.end(), [](e_pointer e, e_pointer f) {return e->index < f->index; });
	//for (size_t i = 0; i < edges_tri.size() - 1; i++) {

	//	unsigned const index = edges_tri[i + 1]->index - edges_tri[i]->index;
	//	//unsigned const index = edges_tri[i]->index;

	//	cout << index << endl;

	//	if (index != 1)
	//		cout << "shit ---------------------------------------- \n----------------------------------------" << endl;

	//}

	//std::stable_sort(triangles_tri.begin(), triangles_tri.end(), [](t_pointer t, t_pointer k) {return t->index < k->index; });
	//for (size_t i = 0; i < triangles_tri.size() - 1; i++) {

	//	unsigned const index = triangles_tri[i + 1]->index - triangles_tri[i]->index;
	//	//unsigned const index = triangles_tri[i]->index;

	//	cout << index << endl;

	//	if (index != 1)
	//		cout << "shit ---------------------------------------- \n----------------------------------------" << endl;

	//}


};


template<GEOMETRIC_KERNEL GK>
bool Triangulation<GK>::get_encroached_segments(v_pointer const p, std::vector<e_pointer> & encroached_segments) {


	for (size_t i = 0; i < segments_tri.size(); i++) {


		e_pointer const e = segments_tri[i];

		if (e->contains(p))
			continue;

		if (is_encroached(e->a, e->b, p))
			encroached_segments.push_back(e);

	}

	return !encroached_segments.empty();

};
template<GEOMETRIC_KERNEL GK>
void Triangulation<GK>::split_encroached_segment(e_pointer const e, v_pointer const p) {



	double x_mid;
	double y_mid;

	get_mid_point(e, x_mid, y_mid);


	v_pointer const new_vertex = new Vertex(x_mid, y_mid);

	new_vertex->index = vertices_tri.size();
	//new_vertex->marker = V_MARKER::FREE;		// By default, a vertex is FREE

	vertices_tri.push_back(new_vertex);


	remove(segments_tri, e);

	e_pointer subsegment1 = NULL;
	e_pointer subsegment2 = NULL;

	// Split segment. Also newly created vertex in the midpont of segment can encroache some other segment! Do not forget to include it to 'vertices_tri' . Done above!
	insert_vertex_into_edge(e, new_vertex, subsegment1, subsegment2);



	segments_tri.push_back(subsegment1);
	segments_tri.push_back(subsegment2);


	// Vertex that encroached older segment can still encroach newly created subsegments
	if (is_encroached(subsegment1->a, subsegment1->b, p))
		split_encroached_segment(subsegment1, p);

	if (is_encroached(subsegment2->a, subsegment2->b, p))
		split_encroached_segment(subsegment2, p);


};

template<GEOMETRIC_KERNEL GK>
void Triangulation<GK>::get_bad_triangles_ruppert() {


	double const Area_Max = Area_Bound;
	double const Bound = 0.5 / (sin(Pi * Angle_Bound / 180.0));

	for (size_t i = 0; i < triangles_tri.size(); i++) {


		t_pointer const t = triangles_tri[i];

		double const circum_radius_squared = t->circum_radius_squared();
		double const ratio_squared = t->ratio_squared();

		if (ratio_squared > sqr(Bound) || circum_radius_squared > sqr(Area_Max))
			bad_triangles.push_back(t);

	}

	bad_triangles.sort([](t_pointer const t, t_pointer const k) { return t->ratio_squared() < k->ratio_squared(); });

};
template<GEOMETRIC_KERNEL GK>
void Triangulation<GK>::enqueue_bad_triangle(t_pointer const t) {


	bool inserted = false;

	for (auto it = bad_triangles.begin(); it != bad_triangles.end(); it++) {

		t_pointer const k = *it;

		const double r1 = t->ratio_squared();
		const double r2 = k->ratio_squared();


		if (t->ratio_squared() < k->ratio_squared()) {

			inserted = true;

			bad_triangles.insert(it, t);

			break;

		}
	}

	if (!inserted)
		bad_triangles.push_back(t);

};







/*****************************************************************************/
/*																			 */
/*	 - Chew's fisrt Delaunay refinement algorithm							 */
/*			: Will split all triangles, which circumradius is greater		 */
/*			  then given parameter h_min (shortest edge in the mesh)		 */
/*			: i.e. the bound is B = 1										 */
/*			: For given parametr parameter h_min, defined segments			 */
/*			  are firstly divided into subsegments of the lenght  			 */
/*			  [ h , sqrt(3)*h ]												 */
/*																			 */
/*****************************************************************************/
template<GEOMETRIC_KERNEL GK>
t_pointer Triangulation<GK>::get_skinny_triangle_chew1st(double &x_c, double &y_c) {


	t_pointer worst_triangle = NULL;

	for (size_t i = 0; i < triangles_tri.size(); i++) {


		t_pointer const t = triangles_tri[i];

		double const circum_radius_squared = find_circumcenter(t->vertices[0], t->vertices[1], t->vertices[2], x_c, y_c);
		double const shortest_edge_squared = lenght_shortest_edge(t);

		if (circum_radius_squared > sqr(Minimal_Edge_Length))
			return t;

	}

	return NULL;

};


template<GEOMETRIC_KERNEL GK>
void Triangulation<GK>::divide_segments(double const h_min) {


	double h = h_min;

	// Check if there are shorter segments than h_min. If yes use the minimum 8-)
	for (size_t i = 0; i < segments_tri.size(); i++) {

		double const e_lenght = segments_tri[i]->length();

		if (e_lenght < h_min)
			h = e_lenght;

	}

	// Check if there are vertices closer to each other than h_min. Check if the distance of the vertices to any segments is closer than h_min


	// Divide segments into parts with the lenght of [ h , sqrt(3) * h ]


};

template<GEOMETRIC_KERNEL GK>
void Triangulation<GK>::refinement_chews1st(double const h_min) {


	bool continue_refining = true;
	bool split_segments = true;

	while (true) {

		continue_refining = false;

		// Split every encroached segments / subsegments there are
		if (split_segments) {


			// While there are encroached segments / subsegments split them
			for (size_t i = 0; i < vertices_tri.size(); i++) {

				// !!!! Check if the vertex is visible. If not, do not add the segment
				v_pointer const v = vertices_tri[i];

				std::vector<e_pointer> encroached_segments;

				if (get_encroached_segments(v, encroached_segments)) {

					continue_refining = true;

					for (size_t j = 0; j < encroached_segments.size(); j++)
						split_encroached_segment(encroached_segments[j], v);

				}
			}
		}


		// Circum center coordinates
		double x_c;
		double y_c;


		// Get any skinny triangle and its circumcenter
		t_pointer const skinny_triangle = get_skinny_triangle_chew1st(h_min, x_c, y_c);


		// If there are no skinny triangles left, we are done
		if (!skinny_triangle && !continue_refining)
			break;


		// Circumcenter of skinny triangle
		Vertex circum_center(x_c, y_c);

		std::vector<e_pointer> encroached_segments2;


		// !!!! Check if the circumcenter is visible. If not, do not add the segment
		if (get_encroached_segments(&circum_center, encroached_segments2)) {


			split_segments = true;

			size_t const k = encroached_segments2.size();

			for (size_t i = 0; i < k; i++)
				split_encroached_segment_no_recursion(encroached_segments2[i]);

			//for (size_t i = 0; i < encroached_segments2.size(); i++)
			//	split_encroached_segment(encroached_segments2[i], &circum_center);

		}
		else {

			split_segments = false;

			// We have to begin from given skinny triangle, so if the circumcenter is behind segments in a hole, it wont crash. Otherwise, we can go to new vertex through hole => there are no triangles (problem with walking algorithm)
			//v_pointer const possible_duplicate = insert_vertex(new Vertex(x_c, y_c), skinny_triangle);
			insert_vertex(new Vertex(x_c, y_c), skinny_triangle);

		}

	}

};































/*

template<GEOMETRIC_KERNEL GK>
bool Triangulation<GK>::get_all_encroached_segments() {


	for (size_t i = 0; i < segments_tri.size(); i++) {

		e_pointer const e = segments_tri[i];

		t_pointer const n1 = e->neighbors[0];
		t_pointer const n2 = e->neighbors[1];

		if (n1 && is_encroached(e->a, e->b, n1->get_vertex(e))) {

			encroached_segments.push_back(e);
			continue;

		}
		if (n2 && is_encroached(e->a, e->b, n2->get_vertex(e))) {

			encroached_segments.push_back(e);
			continue;

		}
	}

	return encroached_segments.empty();

};

template<GEOMETRIC_KERNEL GK>
void Triangulation<GK>::split_encroached_segment(e_pointer const e, v_pointer const p, double const how_far) {


	v_pointer const a = e->a;
	v_pointer const b = e->b;

	double const new_x = a->x + how_far * (b->x - a->x);
	double const new_y = a->y + how_far * (b->y - a->y);

	//double const new_x = (1.0 - how_far) * a->x + how_far * (b->x - a->x);
	//double const new_y = (1.0 - how_far) *a->y + how_far * (b->y - a->y);


	v_pointer const new_vertex = new Vertex(new_x, new_y);

	new_vertex->index = vertices_tri.size();

	vertices_tri.push_back(new_vertex);


	remove(segments_tri, e);

	encroached_segments.remove(e);

	e_pointer subsegment1 = NULL;
	e_pointer subsegment2 = NULL;

	// Split segment. Also newly created vertex in the midpont of segment can encroache some other segment! Do not forget to include it to 'vertices_tri' . Done above!
	insert_vertex_into_edge(e, new_vertex, subsegment1, subsegment2);



	segments_tri.push_back(subsegment1);
	segments_tri.push_back(subsegment2);


	t_pointer n1 = subsegment1->neighbors[0];
	t_pointer n2 = subsegment1->neighbors[1];

	//if (n1 && n1->area() < Area_Bound)
	//	bad_triangles.remove(n1);
	//if (n2 && n2->area() < Area_Bound)
	//	bad_triangles.remove(n2);

	if (n1 && is_encroached(subsegment1->a, subsegment1->b, n1->get_vertex(subsegment1)))
		encroached_segments.push_back(subsegment1);
	else if (n2 && is_encroached(subsegment1->a, subsegment1->b, n2->get_vertex(subsegment1)))
		encroached_segments.push_back(subsegment1);

	n1 = subsegment2->neighbors[0];
	n2 = subsegment2->neighbors[1];

	//if (n1 && n1->area() < Area_Bound)
	//	bad_triangles.remove(n1);
	//if (n2 && n2->area() < Area_Bound)
	//	bad_triangles.remove(n2);

	if (n1 && is_encroached(subsegment2->a, subsegment2->b, n1->get_vertex(subsegment2)))
		encroached_segments.push_back(subsegment2);
	else if (n2 && is_encroached(subsegment2->a, subsegment2->b, n2->get_vertex(subsegment2)))
		encroached_segments.push_back(subsegment2);

};




template<GEOMETRIC_KERNEL GK> template<REFINEMENT_PRIORITY priority>
void Triangulation<GK>::refinement_ruppert(double const angle_bound, double const area_max) {


	Angle_Bound = angle_bound;
	Area_Bound = area_max;


	get_bad_triangles_ruppert();

	bool continue_refining = true;
	bool split_segments = true;

	while (true) {


		continue_refining = false;


		// Split every encroached segments / subsegments there are		
		if (split_segments) {


			bool const no_encroached = get_all_encroached_segments();


			while (!encroached_segments.empty()) {


				std::ofstream edges_txt;
				edges_txt.open("C:\\Users\\pgali\\Desktop\\edgesPoints.txt");
				export_edges(edges_txt);
				edges_txt.close();

				e_pointer const e = encroached_segments.front();

				encroached_segments.remove(e);


				// Use concentric shells if only one endpoint is a vertex of PSLG. That is, it is not labeled as FREE vertex
				bool const is_pslg_vertex_a = e->a->marker == V_MARKER::CONSTRAINED ? true : false;
				bool const is_pslg_vertex_b = e->b->marker == V_MARKER::CONSTRAINED ? true : false;

				if (!is_pslg_vertex_a && !is_pslg_vertex_b) {

					t_pointer const n1 = e->neighbors[0];
					t_pointer const n2 = e->neighbors[1];

					if (n1) {

						v_pointer const apex = n1->get_vertex(e);

						if (is_encroached(e->a, e->b, apex)) {

							split_encroached_segment(e, apex, 0.5);
							continue;

						}
					}
					if (n2) {

						v_pointer const apex = n2->get_vertex(e);

						if (is_encroached(e->a, e->b, apex)) {

							split_encroached_segment(e, apex, 0.5);
							continue;

						}
					}

				}

				// is this encroached segment connected to other segment? That is if the ending point is PSLG vertex with a degree of segment >= 2 ?
				bool const use_concentric_a = is_pslg_vertex_a && vertices_degrees_of_subsegments[e->a->index] >= 2 ? true : false;
				bool const use_concentric_b = is_pslg_vertex_b && vertices_degrees_of_subsegments[e->b->index] >= 2 ? true : false;



				double const length = e->length();

				double nearestpoweroftwo = 1.0;

				while (length > 3.0 * nearestpoweroftwo)
					nearestpoweroftwo *= 2.0;
				while (length < 1.5 * nearestpoweroftwo)
					nearestpoweroftwo *= 0.5;

				double split_position = nearestpoweroftwo / length;

				if (!use_concentric_a && !use_concentric_b)
					split_position = 0.5;
				else if(use_concentric_b)
					split_position = 1.0 - split_position;


				v_pointer const a = e->a;
				v_pointer const b = e->b;

				double const _x = a->x + split_position * (b->x - a->x);
				double const _y = a->y + split_position * (b->y - a->y);

				//double const _x = (1.0 - split_position)* a->x + split_position * (b->x - a->x);
				//double const _y = (1.0 - split_position)* a->y + split_position * (b->y - a->y);


				v_pointer const new_vertex = new Vertex(_x, _y);

				new_vertex->index = vertices_tri.size();
				//new_vertex->marker = V_MARKER::FREE;		// By default, a vertex is FREE

				vertices_tri.push_back(new_vertex);

				remove(segments_tri, e);
				

				e_pointer subsegment1 = NULL;
				e_pointer subsegment2 = NULL;

				// Split segment. Also newly created vertex in the midpont of segment can encroache some other segment! Do not forget to include it to 'vertices_tri' . Done above!
				insert_vertex_into_edge(e, new_vertex, subsegment1, subsegment2);


				segments_tri.push_back(subsegment1);
				segments_tri.push_back(subsegment2);

				t_pointer n1 = subsegment1->neighbors[0];
				t_pointer n2 = subsegment1->neighbors[1];

				if (n1 && n1->area() < Area_Bound)
					bad_triangles.remove(n1);
				if (n2 && n2->area() < Area_Bound)
					bad_triangles.remove(n2);


				if ((n1 && is_encroached(subsegment1->a, subsegment1->b, n1->get_vertex(subsegment1))) || n2 && is_encroached(subsegment1->a, subsegment1->b, n2->get_vertex(subsegment1)))
					encroached_segments.push_back(subsegment1);

				n1 = subsegment2->neighbors[0];
				n2 = subsegment2->neighbors[1];

				if (n1 && n1->area() < Area_Bound)
					bad_triangles.remove(n1);
				if (n2 && n2->area() < Area_Bound)
					bad_triangles.remove(n2);

				if ((n1 && is_encroached(subsegment2->a, subsegment2->b, n1->get_vertex(subsegment2))) || n2 && is_encroached(subsegment2->a, subsegment2->b, n2->get_vertex(subsegment2)))
					encroached_segments.push_back(subsegment2);
				
			}
		}


		std::ofstream edges_txt;
		edges_txt.open("C:\\Users\\pgali\\Desktop\\edgesPoints.txt");
		export_edges(edges_txt);
		edges_txt.close();

		// Get any skinny triangle and its circumcenter
		t_pointer const bad_triangle = !bad_triangles.empty() ? (priority == REFINEMENT_PRIORITY::BEST ? bad_triangles.front() : bad_triangles.back()) : NULL;


		// If there are no skinny triangles left, we are done
		if (!bad_triangle && !continue_refining)
			break;



		// Circum center coordinates
		double x_c;
		double y_c;

		bad_triangle->circum_center(x_c, y_c);


		// Circumcenter of skinny triangle
		Vertex circum_center(x_c, y_c);

		std::vector<e_pointer> encroached_segments2;


		if (get_encroached_segments(&circum_center, encroached_segments2)) {

			while (!encroached_segments2.empty()) {

				split_segments = true;

				//size_t const k = encroached_segments2.size();



				std::ofstream edges_txt;
				edges_txt.open("C:\\Users\\pgali\\Desktop\\edgesPoints.txt");
				export_edges(edges_txt);
				edges_txt.close();

				e_pointer const e = encroached_segments2.front();

				remove(encroached_segments2, e);


				// Use concentric shells if only one endpoint is a vertex of PSLG. That is, it is not labeled as FREE vertex
				bool const is_pslg_vertex_a = e->a->marker == V_MARKER::CONSTRAINED ? true : false;
				bool const is_pslg_vertex_b = e->b->marker == V_MARKER::CONSTRAINED ? true : false;

				if (!is_pslg_vertex_a && !is_pslg_vertex_b) {

					split_encroached_segment(e, &circum_center);
					continue;


				}

				// is this encroached segment connected to other segment? That is if the ending point is PSLG vertex with a degree of segment >= 2 ?
				bool const use_concentric_a = is_pslg_vertex_a && vertices_degrees_of_subsegments[e->a->index] >= 2 ? true : false;
				bool const use_concentric_b = is_pslg_vertex_b && vertices_degrees_of_subsegments[e->b->index] >= 2 ? true : false;



				double const length = e->length();

				double nearestpoweroftwo = 1.0;

				while (length > 3.0 * nearestpoweroftwo)
					nearestpoweroftwo *= 2.0;
				while (length < 1.5 * nearestpoweroftwo)
					nearestpoweroftwo *= 0.5;

				double split_position = nearestpoweroftwo / length;

				if (!use_concentric_a && !use_concentric_b)
					split_position = 0.5;
				//else if (use_concentric_b)
				//	split_position = 1.0 - split_position;


				v_pointer const a = e->a;
				v_pointer const b = e->b;

				double const _x = a->x + split_position * (b->x - a->x);
				double const _y = a->y + split_position * (b->y - a->y);

				//double const _x = (1.0 - split_position)* a->x + split_position * (b->x - a->x);
				//double const _y = (1.0 - split_position)* a->y + split_position * (b->y - a->y);


				v_pointer const new_vertex = new Vertex(_x, _y);

				new_vertex->index = vertices_tri.size();
				//new_vertex->marker = V_MARKER::FREE;		// By default, a vertex is FREE

				vertices_tri.push_back(new_vertex);

				remove(segments_tri, e);


				e_pointer subsegment1 = NULL;
				e_pointer subsegment2 = NULL;

				// Split segment. Also newly created vertex in the midpont of segment can encroache some other segment! Do not forget to include it to 'vertices_tri' . Done above!
				insert_vertex_into_edge(e, new_vertex, subsegment1, subsegment2);


				segments_tri.push_back(subsegment1);
				segments_tri.push_back(subsegment2);


				if (is_encroached(subsegment1->a, subsegment1->b, &circum_center))
					encroached_segments2.push_back(subsegment1);

				if (is_encroached(subsegment2->a, subsegment2->b, &circum_center))
					encroached_segments2.push_back(subsegment2);

			}


			//for (size_t i = 0; i < k; i++)
			//	split_encroached_segment_no_recursion(encroached_segments2[i]);

			//for (size_t i = 0; i < encroached_segments2.size(); i++)
			//	split_encroached_segment(encroached_segments2[i], &circum_center);

		}
		else {

			split_segments = false;

			// We have to begin from given skinny triangle, so if the circumcenter is behind segments in a hole, it wont crash. Otherwise, we can go to new vertex through hole => there are no triangles (problem with walking algorithm)
			//v_pointer const possible_duplicate = insert_vertex(new Vertex(x_c, y_c), skinny_triangle);
			insert_vertex(new Vertex(x_c, y_c), bad_triangle);

		}

	}

};


template<GEOMETRIC_KERNEL GK>
bool Triangulation<GK>::get_encroached_segments(v_pointer const p, std::vector<e_pointer> & encroached_segments) {


	for (size_t i = 0; i < segments_tri.size(); i++) {


		e_pointer const e = segments_tri[i];

		if (e->contains(p))
			continue;

		if (is_encroached(e->a, e->b, p))
			encroached_segments.push_back(e);

	}

	return !encroached_segments.empty();

};
template<GEOMETRIC_KERNEL GK>
void Triangulation<GK>::split_encroached_segment(e_pointer const e, v_pointer const p) {



	double x_mid;
	double y_mid;

	get_mid_point(e, x_mid, y_mid);


	v_pointer const new_vertex = new Vertex(x_mid, y_mid);

	new_vertex->index = vertices_tri.size();
	//new_vertex->marker = V_MARKER::FREE;		// By default, a vertex is FREE

	vertices_tri.push_back(new_vertex);


	remove(segments_tri, e);

	encroached_segments.remove(e);

	e_pointer subsegment1 = NULL;
	e_pointer subsegment2 = NULL;

	// Split segment. Also newly created vertex in the midpont of segment can encroache some other segment! Do not forget to include it to 'vertices_tri' . Done above!
	insert_vertex_into_edge(e, new_vertex, subsegment1, subsegment2);



	segments_tri.push_back(subsegment1);
	segments_tri.push_back(subsegment2);


	t_pointer n1 = subsegment1->neighbors[0];
	t_pointer n2 = subsegment1->neighbors[1];

	if (n1 && n1->area() < Area_Bound)
		bad_triangles.remove(n1);
	if (n2 && n2->area() < Area_Bound)
		bad_triangles.remove(n2);

	if ((n1 && is_encroached(subsegment1->a, subsegment1->b, n1->get_vertex(subsegment1))) || n2 && is_encroached(subsegment1->a, subsegment1->b, n2->get_vertex(subsegment1)))
		encroached_segments.push_back(subsegment1);

	n1 = subsegment2->neighbors[0];
	n2 = subsegment2->neighbors[1];

	if (n1 && n1->area() < Area_Bound)
		bad_triangles.remove(n1);
	if (n2 && n2->area() < Area_Bound)
		bad_triangles.remove(n2);

	if ((n1 && is_encroached(subsegment2->a, subsegment2->b, n1->get_vertex(subsegment2))) || n2 && is_encroached(subsegment2->a, subsegment2->b, n2->get_vertex(subsegment2)))
		encroached_segments.push_back(subsegment2);



	// Vertex that encroached older segment can still encroach newly created subsegments
	//if (is_encroached(subsegment1->a, subsegment1->b, p)) {

	//	if (subsegment1->a->marker == V_MARKER::CONSTRAINED && subsegment1->b->marker != V_MARKER::CONSTRAINED)
	//		concentric_shell_split(subsegment1, p);
	//	else if (subsegment1->a->marker != V_MARKER::CONSTRAINED && subsegment1->b->marker == V_MARKER::CONSTRAINED)
	//		concentric_shell_split(subsegment1, p);
	//	else
	//		split_encroached_segment(subsegment1, p);

	//}


	//if (is_encroached(subsegment2->a, subsegment2->b, p)) {

	//	if (subsegment2->a->marker == V_MARKER::CONSTRAINED && subsegment2->b->marker != V_MARKER::CONSTRAINED)
	//		concentric_shell_split(subsegment2, p);
	//	else if (subsegment2->a->marker != V_MARKER::CONSTRAINED && subsegment2->b->marker == V_MARKER::CONSTRAINED)
	//		concentric_shell_split(subsegment2, p);
	//	else
	//		split_encroached_segment(subsegment2, p);

	//}


};

template<GEOMETRIC_KERNEL GK>
void Triangulation<GK>::concentric_shell_split(e_pointer const e, v_pointer const p) {



	

};

template<GEOMETRIC_KERNEL GK>
void Triangulation<GK>::get_bad_triangles_ruppert() {


	double const Area_Max = Area_Bound;
	double const Bound = 0.5 / (sin(Pi * Angle_Bound / 180.0));

	for (size_t i = 0; i < triangles_tri.size(); i++) {


		t_pointer const t = triangles_tri[i];

		double const circum_radius_squared = t->circum_radius_squared();
		double const ratio_squared = t->ratio_squared();

		if (ratio_squared > sqr(Bound) || circum_radius_squared > sqr(Area_Max))
			bad_triangles.push_back(t);

	}

	bad_triangles.sort([](t_pointer const t, t_pointer const k) { return t->ratio_squared() < k->ratio_squared(); });

};
template<GEOMETRIC_KERNEL GK>
void Triangulation<GK>::enqueue_bad_triangle(t_pointer const t) {


	bool inserted = false;

	for (auto it = bad_triangles.begin(); it != bad_triangles.end(); it++) {

		t_pointer const k = *it;

		const double r1 = t->ratio_squared();
		const double r2 = k->ratio_squared();


		if (t->ratio_squared() < k->ratio_squared()) {

			inserted = true;

			bad_triangles.insert(it, t);

			break;

		}
	}

	if (!inserted)
		bad_triangles.push_back(t);

};


*/