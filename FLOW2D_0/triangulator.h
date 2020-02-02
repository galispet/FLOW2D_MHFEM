#pragma once

#include "shapes.h"
#include "geometric.h"
#include "enumerators.h"

#include <iostream>
#include <vector>
#include <list>

#include <fstream>


using std::cout;
using std::endl;


typedef std::iterator<std::input_iterator_tag, Vertex>		input_iterator;
typedef std::iterator<std::input_iterator_tag, Vertex *>	input_iterator_handle;

typedef Vertex *										Vertex_handle;
typedef Edge *											Edge_handle;
typedef Triangle *										Triangle_handle;



struct Handle {

	Vertex_handle Vertex_Handle			= NULL;
	Edge_handle	Edge_Handle				= NULL;
	Triangle_handle	Triangle_Handle		= NULL;

};


// This must be implemnted, so in the function 'insert_individual_constraint' user can't pass 'Vertex *' as argument -> it will be allowed only for 'Vertex_handle'. Therefore it is sure, that the vertex is already in triangulation (compare with : insert_individual_constraint(Vertex * a, Vertex * b) )
//class Vertex_handle {
//
//
//	Vertex_handle() {};
//	~Vertex_handle() {};
//
//private:
//
//	Vertex * v_handle = NULL;
//
//	
//	Vertex_handle & operator=(Vertex_handle & v) { v_handle = v.v_handle; };
//
//	Vertex * operator->() { return v_handle; };
//	Vertex & operator*() { return *v_handle; };
//
//};



template <GEOMETRIC_KERNEL GK, CHECK_KERNEL CK>
class CDT {


	friend class Mesh;

public:


	/*****************************************************************************/
	/*                                                                           */
	/*    - Constructor and destructor					                         */
	/*																			 */
	/*****************************************************************************/
	CDT();
	~CDT();


public:

	/*****************************************************************************/
	/*																			 */
	/*    - Inserts incrementaly single vertex into triangulation                */
	/*																			 */
	/*****************************************************************************/
	Vertex_handle insert_vertex(Vertex v) {

		Vertex * const new_vertex = new Vertex(v.x, v.y);

		input_vertex(new_vertex, NULL);

		return new_vertex;

	};
	Vertex_handle insert_vertex(Vertex_handle const v) {

		return input_vertex(v, NULL);

	};


	/*****************************************************************************/
	/*																			 */
	/*    - Inserts vertices contained in vector into triangulation using		 */
	/*	    'insert_vertex'	on each vertex										 */
	/*																			 */
	/*****************************************************************************/
	template<class input_iterator>
	unsigned insert_vertex(input_iterator first, input_iterator last) {


		unsigned n = 0;

		for (auto it = first; it != last; it++) {

			Vertex const v = (*it);

			double const x = v.x;
			double const y = v.y;

			Vertex * const new_vertex = new Vertex(x, y);

			input_vertex(new_vertex);

			n++;

		}

		return n;

	};
	// Also implement 'insert_vertex(input_iterator_handles first, input_iterator_handles last)()', so we can do it with vertices already inserted -> Again !! make sure to implement safe Vertex_handle class, so user can't modifie or input anything

private:

	Vertex_handle input_vertex(Vertex_handle const v, Triangle * t = NULL) {



		// Resizing the super triangle if needed
		if (abs(v->x) > maximum_coordinate || abs(v->y) > maximum_coordinate){

			maximum_coordinate = fmax(abs(v->y), abs(v->x));

			double const mid = 0.0;									 // 0.5*(minimum_coordinate + maximum_coordinate);
			double const M = 2 * maximum_coordinate;					// maximum_coordinate - minimum_coordinate;

			st_v0->set(mid - 20 * M, mid - M);
			st_v1->set(mid + 20 * M, mid - M);
			st_v2->set(mid, mid + 20 * M);
		
		}


		Handle rh;

		switch (locate(rh, v, t)) {

		case LOCATION::IN_TRIANGLE:

			triangle_event(rh.Triangle_Handle, v);
			break;

		case LOCATION::ON_EDGE:

			edge_event(rh.Edge_Handle, v);
			break;

		case LOCATION::ON_VERTEX:

			std::cout << "Duplicate vertex ignored" << std::endl;
			return rh.Vertex_Handle;
			//return v;

		}


		v->index = num_vertices;

		vertices.push_back(v);

		num_vertices = vertices.size();
		num_edges = edges.size();
		num_triangles = triangles.size();

		return v;

	};

public:


	/*****************************************************************************/
	/*																			 */
	/*    - Inserts single vertex into provided edge 'e'						 */
	/*    - Returns pointer to any one of four created triangle, because		 */	
	/*      we need to denoted some of created edges as 'constrained'			 */	
	/*																			 */
	/*****************************************************************************/
	Triangle * insert_vertex_into_edge(Vertex * const v, Edge * const e) {


		v->index = num_vertices;

		Triangle * const k0 = edge_event(e, v);

		vertices.push_back(v);

		num_vertices = vertices.size();
		num_edges = edges.size();
		num_triangles = triangles.size();

		return k0;

	};


	/*****************************************************************************/
	/*                                                                           */
	/*    - Free the allocated memory. Destroying mesh etc.                      */
	/*																			 */
	/*****************************************************************************/
	void clear();



	/*****************************************************************************/
	/*                                                                           */
	/*    - Getters for the containers of objects of triangulation primitives	 */
	/*																			 */
	/*****************************************************************************/
	std::vector<Triangle > get_triangles();
	std::vector<Edge > get_edges();
	std::vector<Vertex > get_new_vertices();


	/*****************************************************************************/
	/*                                                                           */
	/*    - Getters for the containers of handles to triangulation primitives	 */
	/*																			 */
	/*****************************************************************************/
	// Don't forget to secure access for the user from changing anything about these primitives !!!!
	std::vector<Vertex_handle> get_vertices_handles();
	std::vector<Edge_handle> get_edges_handles();
	std::vector<Triangle_handle> get_triangles_handles();
	std::vector<Vertex_handle> get_new_vertices_handles();


	/*****************************************************************************/
	/*                                                                           */
	/*    - Getters for triangulaion primitives									 */
	/*																			 */
	/*****************************************************************************/
	Vertex_handle get_vertex(unsigned i) { 
		
		return vertices[i];
	
	}
	Edge_handle get_edge(unsigned i) {
		
		return edges[i]; 
	
	}
	Triangle_handle get_triangle(unsigned i) {
		
		return triangles[i]; 
	
	}


	/*****************************************************************************/
	/*                                                                           */
	/*    - Returns number of triangulation primitives							 */
	/*																			 */
	/*****************************************************************************/
	unsigned get_num_vertices() { return (unsigned)vertices.size(); };
	unsigned get_num_edges() { return (unsigned)edges.size(); };
	unsigned get_num_triangles() { return (unsigned)triangles.size(); };


	/*****************************************************************************/
	/*                                                                           */
	/*    - Returns number of edges with 'NEUMANN' marker						 */
	/*    - Returns number of edges with 'DIRICHLET' marker						 */
	/*																			 */
	/*****************************************************************************/
	unsigned get_num_neumann_edges() { return (unsigned)std::count_if(edges.begin(), edges.end(), [](Edge_handle e) {return e->marker == E_MARKER::NEUMANN; }); };
	unsigned get_num_dirichlet_edges() { return (unsigned)std::count_if(edges.begin(), edges.end(), [](Edge_handle e) {return e->marker == E_MARKER::DIRICHLET; }); };


	/*****************************************************************************/
	/*                                                                           */
	/*    - Write primitive's coordinates to a text-file						 */
	/*																			 */
	/*****************************************************************************/
	void export_vertices(std::ofstream & stream) const;
	void export_edges(std::ofstream & stream) const;
	void export_triangles(std::ofstream & stream) const;


private:


	/*****************************************************************************/
	/*                                                                           */
	/*    - Number of each triangulation primitive								 */
	/*    - Number of edges with boundary marks 'NEUMANN' and 'DIRICHLET'		 */
	/*																			 */
	/*****************************************************************************/
	public:

	unsigned num_vertices = 0;
	unsigned num_edges = 0;
	unsigned num_triangles = 0;

	//unsigned num_dirichlet = 0;
	//unsigned num_neumann = 0;

	private:


	/*****************************************************************************/
	/*                                                                           */
	/*    - Maximum magnitude in dimension of coordinates                        */
	/*    - It is used for 'blowing up' super triangle	                         */
	/*																			 */
	/*****************************************************************************/
	double maximum_coordinate;


	/*****************************************************************************/
	/*                                                                           */
	/*    - Containers for triangulation primitives			                     */
	/*																			 */
	/*****************************************************************************/
	std::vector<Vertex *> vertices;
	std::vector<Edge *> edges;
	std::vector<Triangle *> triangles;


	/*****************************************************************************/
	/*                                                                           */
	/*    - Container for additonal vertices resulting from inserting			 */
	/*		constraints															 */
	/*																			 */
	/*****************************************************************************/
	std::vector<Vertex *> constraints_new_vertices;


	//MAKE new vertices positioned in between two ending points on the edge -> better quality of triangles


	/*****************************************************************************/
	/*                                                                           */
	/*    - Containers for steiner vertices					                     */
	/*    - Used to maintain Delaunay property of constrained triangulation      */
	/*																			 */
	/*    - Not yet implemented								                     */
	/*																			 */
	/*****************************************************************************/
	std::vector<Vertex *> steiner_points;


	/*****************************************************************************/
	/*                                                                           */
	/*    - Vertices and edges of super triangle			                     */
	/*																			 */
	/*****************************************************************************/
	Vertex * st_v0;
	Vertex * st_v1;
	Vertex * st_v2;

	Edge * st_e0;
	Edge * st_e1;
	Edge * st_e2;


	/*****************************************************************************/
	/*                                                                           */
	/*    - Geometric functions													 */
	/*    - Template parametr 'GK' provides 'EXACT' or 'INEXACT' computation     */
	/*																			 */
	/*****************************************************************************/
	Predicates<GK> const predicates;





	/*****************************************************************************/
	/*                                                                           */
	/*    - Point location routine												 */
	/*	  - Return handles to either: Vertex, Edge, Triangle					 */
	/*	    where the inserted point is located									 */
	/*																			 */
	/*****************************************************************************/
	LOCATION locate(Handle & location_result, Vertex * const p, Triangle * initial_triangle = NULL) {


		Triangle * t = initial_triangle == NULL ? locate_triangle_straight_walk(triangles.back(), p) : initial_triangle;

		//Triangle * const t = locate_triangle_barycentric(triangles.back(), p);
		//Triangle * const t = locate_triangle_straight_walk(triangles.back(), p);

		//Triangle * const t = locate_triangle_lawson_walk(triangles.back(), p);
		//Triangle * const t = locate_triangle_lawson_walk_remembering(triangles.back(), p);
		//Triangle * const t = locate_triangle_lawson_walk_stochastic(triangles.back(), p);
		//Triangle * const t = locate_triangle_lawson_walk_remembering_stochastic(triangles.back(), p);
		//Triangle * const t = locate_triangle_lawson_walk_fast_remembering(triangles.back(), p, 10);


		// For the validity :) -> that's why i need to make point location more robust
		if (!predicates.in_triangle_robust(t, p)) {

			std::cout << "Inserted vertex is not in the resulting Triangle. Checking neighbors ... ";

			for (unsigned i = 0; i < 3; i++) {

				Triangle * const k = t->neighbors[i];

				if (predicates.in_triangle_robust(k, p)) {

					t = k;

					cout << "Ok" << endl;
					break;

				}
				else if (i == 2) {

					cout << " : Not Ok. Trying to locate again." << endl;
					return locate(location_result, p, k);

				}
			}
		}
			


		if (CK == CHECK_KERNEL::CHECK_DUPLICATE) {

			Vertex * const v0 = t->vertex(0);
			Vertex * const v1 = t->vertex(1);
			Vertex * const v2 = t->vertex(2);

			if (p->almost_equal(v0)) {

				location_result.Vertex_Handle = v0;
				return LOCATION::ON_VERTEX;

			}
			else if (p->almost_equal(v1)) {

				location_result.Vertex_Handle = v1;
				return LOCATION::ON_VERTEX;

			}
			else if (p->almost_equal(v2)) {

				location_result.Vertex_Handle = v2;
				return LOCATION::ON_VERTEX;

			}
		}


		bool const b0 = predicates.on_edge(t->edge(0), p);
		bool const b1 = predicates.on_edge(t->edge(1), p);
		bool const b2 = predicates.on_edge(t->edge(2), p);


		if (b0) {

			location_result.Edge_Handle = t->edge(0);
			return LOCATION::ON_EDGE;

		}
		else if (b1) {

			location_result.Edge_Handle = t->edge(1);
			return LOCATION::ON_EDGE;

		}
		else if (b2) {

			location_result.Edge_Handle = t->edge(2);
			return LOCATION::ON_EDGE;

		}

		location_result.Triangle_Handle = t;

		return LOCATION::IN_TRIANGLE;

	};


	/*****************************************************************************/
	/*                                                                           */
	/*    - Point localization algorithms                                        */
	/*																			 */
	/*****************************************************************************/
	Triangle * locate_triangle_barycentric(Triangle * t, Vertex * const p);
	Triangle * locate_triangle_straight_walk(Triangle * t, Vertex * const p);

	Triangle * locate_triangle_lawson_walk(Triangle * t, Vertex * const p);
	Triangle * locate_triangle_lawson_walk_remembering(Triangle * t, Vertex * const p);
	Triangle * locate_triangle_lawson_walk_stochastic(Triangle * t, Vertex * const p);
	Triangle * locate_triangle_lawson_walk_remembering_stochastic(Triangle * t, Vertex * const p);
	Triangle * locate_triangle_lawson_walk_fast_remembering(Triangle * t, Vertex * const p, unsigned const n);


	/*****************************************************************************/
	/*                                                                           */
	/*    - Legalization of edges routine										 */
	/*	  - Rotate two neighboring triangles CW in order to 'flip' edge			 */
	/*	    which triangles share		  									     */
	/*																			 */
	/*****************************************************************************/
	void rotate_triangle_pair(Triangle * t, Vertex * v, Triangle * ot, Vertex * ov);
	bool legalize(Triangle * t, Vertex * v);
	

	/*****************************************************************************/
	/*                                                                           */
	/*    - Routines which create new triangles based on where the               */
	/*		inserted point is located											 */
	/*																			 */
	/*****************************************************************************/
	void triangle_event(Triangle * t, Vertex * p);
	Triangle * edge_event(Edge * e, Vertex * p);


	/*****************************************************************************/
	/*                                                                           */
	/*    - Construction and destruction of 'super triangle'                     */
	/*	  - Delete all Triangles and edges connected to super triangle.			 */
	/*		Deletion of super triangle itself									 */
	/*																			 */
	/*****************************************************************************/
private:	
	
	void super_triangle_construct();
	void super_triangle_remove();
	

	/*****************************************************************************/
	/*                                                                           */
	/*    - Remove triangle or edge from respective vectors		                 */
	/*	  - Frees allocated memory and erase from vector						 */
	/*																			 */
	/*****************************************************************************/
	template<typename T>
	void remove(std::vector<T> & vec, T t);




	/*****************************************************************************/
	/*                                                                           */
	/*    - Constrained triangulation							                 */
	/*																			 */
	/*****************************************************************************/

	private:

		template<E_MARKER marker>
		Triangle * edge_location(Vertex * const a, Vertex * const b,  Triangle * const alpha = NULL);

		template<E_MARKER marker>
		void insert_constraint_intersection(Vertex * const a, Vertex * const b, Triangle * alpha = NULL);
		template<E_MARKER marker>
		void insert_constraint_flip(Vertex * const a, Vertex * const b, Triangle * alpha = NULL);
		template<E_MARKER marker>
		void insert_constraint_enforce(Vertex * const a, Vertex * const b, Triangle * alpha = NULL) {};

		// Doesn't work correctly ->try to fix it (Sometimes it tends to insert 'infinite' number of mid-vertices and loop wants to go 4ever. also some weird constrained edges shows up and some do not. try N1 = 100, N2 = 50, airfoil3)
		template<E_MARKER marker>
		void insert_constraint_halving(Vertex * const a, Vertex * const b, Triangle * alpha = NULL);


		void mark_triangle(Triangle * const t);
		void remove_triangles_in_hole(Vertex * const v, Triangle * alpha = NULL);

	

	public:

		unsigned num_deleted_vertices = 0;

		void make_individual_hole(Vertex * const seed);
		void make_individual_hole(Vertex seed);

		template<class input_iterator>
		void make_holes(input_iterator first_seed, input_iterator last_seed);

		template<INSERT_CONSTRAINT algorithm, E_MARKER marker>
		void insert_individual_constraint(Vertex * const a, Vertex * const b) {



			switch (algorithm) {

			case INSERT_CONSTRAINT::INTERSECTION:

				insert_constraint_intersection<marker>(a, b, a->adjacent_triangle);
				break;

			case INSERT_CONSTRAINT::HALVING:

				insert_constraint_halving<marker>(a, b, a->adjacent_triangle);
				break;

			case INSERT_CONSTRAINT::FLIP:

				insert_constraint_flip<marker>(a, b, a->adjacent_triangle);
				break;

			case INSERT_CONSTRAINT::ENFORCE:

				insert_constraint_enforce<marker>(a, b, a->adjacent_triangle);
				break;

			}

		};
		template<INSERT_CONSTRAINT algorithm, POLYGON_TYPE type, E_MARKER marker, class input_iterator>
		unsigned insert_sequential_constraint(input_iterator first, input_iterator last) {


			// Need to create this temporary container, so there can be inserted constraints (vertices of constraints must be already in triangulation. That doesn't hold for first <-> last)
			std::vector<Vertex * > temp;


			unsigned n = 0;

			for (auto it = first; it != last; it++) {


				Vertex const v = (*it);

				double const x = v.x;
				double const y = v.y;

				Vertex * const new_vertex = new Vertex(x, y);
				
				Vertex * const possible_duplicate = insert_vertex(new_vertex);

				// In the case, that a vertex of a vertices sequence is almost equal to some vertex already present in triangulation 
				if (possible_duplicate != new_vertex) {

					temp.push_back(possible_duplicate);
					delete new_vertex;

				}
					
				else
					temp.push_back(new_vertex);

				n++;

			}


			unsigned const nv = n;

			for (size_t i = 0; i < nv - 1; i++) {


				Vertex * const a = temp[i];
				Vertex * const b = temp[i + 1];

				a->is_constrained = true;
				b->is_constrained = true;

				insert_individual_constraint<algorithm, marker>(a, b);

			}

			//if (type == POLYGON_TYPE::CLOSED)
			if (type != POLYGON_TYPE::OPEN)
				insert_individual_constraint<algorithm, marker>(temp[nv - 1], temp[0]);


			return n;

		};
	

	/*****************************************************************************/
	/*                                                                           */
	/*    - Laplacian smoothing													 */
	/*																			 */
	/*****************************************************************************/

	public:

		template<SMOOTHING algorithm>
		void apply_smoothing(unsigned const N) {


			switch (algorithm) {

			case SMOOTHING::LAPLACIAN:

				laplacian_smoothing(N);
				break;

			}

		};


		/*****************************************************************************/
		/*                                                                           */
		/*    - Check for delaunay property of vertex's adjacent triangles 			 */
		/*		and legalize them if needed 										 */
		/*																			 */
		/*	  - Input is vector iterator. If not provided, then all vertices are	 */
		/*		checked																 */
		/*																			 */
		/*****************************************************************************/
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
		/*	  : Triangle quality metric : measures of qualities of triangles		 */
		/*																			 */
		/*****************************************************************************/
		std::vector<double> mesh_quality_angle;			// Minimum angle
		std::vector<double> mesh_quality_radius;		// Incircle radius
		std::vector<double> mesh_quality_diameter;		// Length of shortest edge

		void laplacian_smoothing(unsigned const N);

		void set_constrained_convex_hull(std::vector<bool> & is_on_convex_hull);
		
		unsigned get_adjacent_triangles_and_vertices(Vertex * const v, std::vector<Triangle *> & triangles, std::vector<Vertex *> & vertices);

		void compute_mesh_quality(MESH_QUALITY quality);




		/*****************************************************************************/
		/*                                                                           */
		/*    - Delaunay's refinement algorithm										 */
		/*																			 */
		/*																			 */
		/*****************************************************************************/

		double shortest_edge_length() {


			double shortest_edge_squared = INFINITY;

			for (size_t i = 0; i < triangles.size(); i++) {


				Triangle * const t = triangles[i];

				Vertex * const a = t->vertices[0];
				Vertex * const b = t->vertices[1];
				Vertex * const c = t->vertices[2];

				double const a_x = a->x;
				double const a_y = a->y;

				double const b_x = b->x;
				double const b_y = b->y;

				double const c_x = c->x;
				double const c_y = c->y;

				double const A = b_x - a_x;
				double const B = b_y - a_y;
				double const C = c_x - a_x;
				double const D = c_y - a_y;

				double const H = c_x - b_x;
				double const I = c_y - b_y;


				double const AB = A * A + B * B;
				double const BC = H * H + I * I;
				double const CA = C * C + D * D;

				double const this_min = std::fmin(AB, std::fmin(BC, CA));

				if (this_min < shortest_edge_squared)
					shortest_edge_squared = this_min;

			}

			return shortest_edge_squared;

		}

		public:

			void refinement_delaunay(double const min_angle);
			void refinement_ruppert(double const Bound);





};






/*****************************************************************************/
/*                                                                           */
/*    - Constructor and destructor					                         */
/*																			 */
/*****************************************************************************/
template<GEOMETRIC_KERNEL GK, CHECK_KERNEL CK>
CDT<GK, CK>::CDT() {

	super_triangle_construct();

};
template<GEOMETRIC_KERNEL GK, CHECK_KERNEL CK>
CDT<GK, CK>::~CDT() {

	clear();

	// These are not deleted in 'clear()'
	delete st_v0;
	delete st_v1;
	delete st_v2;

};



/*****************************************************************************/
/*                                                                           */
/*    - Point location algorithms                                            */
/*																			 */
/* 	locate_triangle_barycentric - quite fast								 */
/* 								- on structered unsorted meshes will loop	 */
/*                                                                           */
/* 	locate_triangle_straight_walk - very fast								 */
/* 								  - on structered unsorted meshes will fail	 */
/*                                                                           */
/*                                                                           */
/*                                                                           */
/* 	locate_triangle_lawson_walk - moderatery fast							 */
/* 								- didn't fail								 */
/*                                                                           */
/* 	locate_triangle_lawson_walk_remembering - faster than original			 */
/* 										    - didn't fail					 */
/*                                                                           */
/* 	locate_triangle_lawson_walk_stochastic - slower than original			 */
/* 										   - didn't fail					 */
/*                                                                           */
/* 	locate_triangle_lawson_walk_remembering_stochastic - slower than original*/
/* 													   - didn't fail		 */
/*                                                                           */
/*                                                                           */
/*                                                                           */
/*                                                                           */
/*   TODO - Take care of special cases in straight walk algorithm and		 */
/*		    barycentric coordinates										     */
/*        - Namely if some point happen to lie on the line 'pq' etc.		 */
/*        - If point on line, repeat 'initialization' step and twist around  */
/*		  - e.g. here are robust codes :									 */
/*																			 */
/* https://github.com/eppz/Triangle.NET/blob/master/Triangle.NET/Triangle/TriangleLocator.cs */
/* https://github.com/libigl/triangle/blob/master/triangle.c								 */
/*																	         */
/*                                                                           */
/*                                                                           */
/*                                                                           */
/*****************************************************************************/
template<GEOMETRIC_KERNEL GK, CHECK_KERNEL CK>
Triangle * CDT<GK, CK>::locate_triangle_barycentric(Triangle * t, Vertex * const p) {


	Vertex * r = NULL;
	Vertex * l = NULL;

	// Computation of barycentric coordinates is optimized -> instead of dividing both 'b0','b1' by 'd' we multiply the 'b2' by 'd' -> much faster

	double c[3];

	double d = predicates.orientation(t->vertex(0), t->vertex(1), t->vertex(2));

	c[0] = predicates.orientation(t->vertex(1), t->vertex(2), p);
	c[1] = predicates.orientation(t->vertex(2), t->vertex(0), p);
	c[2] = d - c[0] - c[1];


	double min = c[0];

	r = t->vertex(1);
	l = t->vertex(2);

	if (c[1] < min) {

		r = t->vertex(2);
		l = t->vertex(0);

		min = c[1];

	}
	if (c[2] < min) {

		r = t->vertex(0);
		l = t->vertex(1);

		min = c[2];

	}

	Triangle * tau = NULL;

	while (min < 0.0) {

		tau = t->neighbor(r, l);

		unsigned int i = tau->neighbor_index(t);

		t = tau;

		c[i] = -min;
		c[(i + 1) % 3] = predicates.orientation(tau->vertex((i + 2) % 3), tau->vertex((i + 3) % 3), p);

		d = predicates.orientation(tau->vertex(0), tau->vertex(1), tau->vertex(2));

		c[(i + 2) % 3] = d - c[i] - c[(i + 1) % 3];

		min = c[0];

		r = tau->vertex(1);
		l = tau->vertex(2);

		if (c[1] < min) {

			r = tau->vertex(2);
			l = tau->vertex(0);

			min = c[1];

		}
		if (c[2] < min) {

			r = tau->vertex(0);
			l = tau->vertex(1);

			min = c[2];

		}

	}

	return t;

};
template<GEOMETRIC_KERNEL GK, CHECK_KERNEL CK>
Triangle * CDT<GK, CK>::locate_triangle_straight_walk(Triangle * t, Vertex * const p) {


	Vertex * q = t->vertex(0);
	Vertex * r = t->vertex(1);
	Vertex * l = t->vertex(2);


	Vertex * s = NULL;

	if (predicates.orientation(r, q, p) <= 0.0) {

		while (predicates.orientation(l, q, p) < 0.0) {

			r = l;
			t = t->neighbor(q, l);
			l = t->vertex_but(q, r);

		}
	}
	else {

		do {

			l = r;
			t = t->neighbor(q, r);
			r = t->vertex_but(q, l);

		} while (predicates.orientation(r, q, p) > 0.0);
	}

	while (predicates.orientation(p, r, l) <= 0.0) {

		t = t->neighbor(r, l);
		s = t->vertex_but(r, l);

		if (predicates.orientation(s, q, p) < 0.0)
			r = s;
		else
			l = s;

	}

	return t;

};

template<GEOMETRIC_KERNEL GK, CHECK_KERNEL CK>
Triangle * CDT<GK, CK>::locate_triangle_lawson_walk(Triangle * t, Vertex * const p) {


	bool found = false;

	Vertex * r = NULL;
	Vertex * l = NULL;

	while (!found) {

		found = true;

		for (unsigned i = 0; i < 3; i++) {

			r = t->vertex(i % 3);
			l = t->vertex((i + 1) % 3);

			if (predicates.orientation(r, l, p) < 0.0) {

				t = t->neighbor(l, r);
				found = false;

				break;

			}

		}

	}

	return t;

};
template<GEOMETRIC_KERNEL GK, CHECK_KERNEL CK>
Triangle * CDT<GK, CK>::locate_triangle_lawson_walk_remembering(Triangle * t, Vertex * const p) {


	Triangle * previous = t;

	bool found = false;

	Vertex * r = NULL;
	Vertex * l = NULL;

	while (!found) {

		found = true;

		for (unsigned i = 0; i < 3; i++) {

			r = t->vertex(i % 3);
			l = t->vertex((i + 1) % 3);

			if (previous != t->neighbor(l, r)) {

				if (predicates.orientation(r, l, p) < 0.0) {

					previous = t;
					t = t->neighbor(l, r);
					found = false;

					break;

				}

			}

		}

	}

	return t;

};
template<GEOMETRIC_KERNEL GK, CHECK_KERNEL CK>
Triangle * CDT<GK, CK>::locate_triangle_lawson_walk_stochastic(Triangle * t, Vertex * const p) {


	Triangle * previous = t;

	bool found = false;

	Vertex * r = NULL;
	Vertex * l = NULL;

	while (!found) {

		found = true;

		unsigned const k = rand() % 3;

		for (unsigned i = k; i < k + 3; i++) {

			r = t->vertex(i % 3);
			l = t->vertex((i + 1) % 3);

			if (predicates.orientation(r, l, p) < 0.0) {

				previous = t;
				t = t->neighbor(l, r);
				found = false;

				break;

			}


		}

	}

	return t;

};
template<GEOMETRIC_KERNEL GK, CHECK_KERNEL CK>
Triangle * CDT<GK, CK>::locate_triangle_lawson_walk_remembering_stochastic(Triangle * t, Vertex * const p) {


	Triangle * previous = t;

	bool found = false;

	Vertex * r = NULL;
	Vertex * l = NULL;

	while (!found) {

		found = true;

		unsigned const k = rand() % 3;

		for (unsigned i = k; i < k + 3; i++) {

			r = t->vertex(i % 3);
			l = t->vertex((i + 1) % 3);

			if (previous != t->neighbor(l, r)) {

				if (predicates.orientation(r, l, p) < 0.0) {

					previous = t;
					t = t->neighbor(l, r);
					found = false;

					break;

				}


			}

		}

	}

	return t;

};
template<GEOMETRIC_KERNEL GK, CHECK_KERNEL CK>
Triangle * CDT<GK, CK>::locate_triangle_lawson_walk_fast_remembering(Triangle * t, Vertex * const p, unsigned const n) {


	Triangle * psi = t;

	Triangle * orig = t;

	Vertex * r = NULL;
	Vertex * l = NULL;

	for (unsigned k = 0; k < n; k++) {

		unsigned j = 0;

		for (unsigned i = 0; i < 3; i++) {

			l = t->vertex(i % 3);
			r = t->vertex((i + 1) % 3);

			j = i;

			if (psi != t->neighbor(r, l))
				break;

		}

		psi = t;

		if (predicates.orientation(l, r, p) < 0.0)
			t = t->neighbor(r, l);

		else {

			l = r;
			r = t->vertex((j + 2) % 3);

			t = t->neighbor(r, l);

		}

		if (!t)
			return locate_triangle_lawson_walk(orig, p);

	}

	return locate_triangle_lawson_walk(t, p);

};


/*****************************************************************************/
/*                                                                           */
/*    - Legalization of edges routines										 */
/*																			 */
/*****************************************************************************/
template<GEOMETRIC_KERNEL GK, CHECK_KERNEL CK>
void CDT<GK, CK>::rotate_triangle_pair(Triangle * t, Vertex * v, Triangle * ot, Vertex * ov) {


	// Save neighbors of triangles
	Triangle * const n0 = t->neighbor_ccw(v);
	Triangle * const n1 = t->neighbor_cw(v);
	Triangle * const n2 = ot->neighbor_ccw(ov);
	Triangle * const n3 = ot->neighbor_cw(ov);

	// Save edges of triangles
	Edge * const flip_edge = t->edge(t->vertex_index(v));

	Edge * const e0 = t->edge_ccw(v);
	Edge * const e1 = t->edge_cw(v);
	Edge * const e2 = ot->edge_ccw(ov);
	Edge * const e3 = ot->edge_cw(ov);

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
template<GEOMETRIC_KERNEL GK, CHECK_KERNEL CK>
bool CDT<GK, CK>::legalize(Triangle * t, Vertex * v) {


	int const edgeToFlip = t->vertex_index(v);

	if (t->edge(edgeToFlip)->is_constrained)
		return false;

	Triangle * const ot = t->neighbor(edgeToFlip);

	// When hole 'ot' may be NULL
	if (!ot)
		return false;

	Vertex * const ov = ot->opposite_vertex(t, v);

	// Think about this more ->it should be useless, because, when edges do not intersect, the 'v' vertex is not contained in the triangle 'ot'
	//if (!intersecting_edges(t->vertex_ccw(v), t->vertex_cw(v), ov, v) && predicates.in_circle(ot, v)) {
	//
	//	 std::cout << "Shit" << std::endl;
	//	return;
	//
	//}

	//const bool in = in_circle_fast(ot, v);
	//const bool in = in_circle_robust(ot, v);

	bool const in = predicates.in_circle(ot, v);

	if (in) {

		rotate_triangle_pair(t, v, ot, ov);

		// Rotating triangles can make vertices's pointers of quadrilateral (except 'v'), point to different triangle, which does not contain respective vertices
		ov					->adjacent_triangle = ot;
		t->vertex_cw(v)		->adjacent_triangle = t;
		ot->vertex_cw(ov)	->adjacent_triangle = ot;


		legalize(t, v);
		legalize(ot, v);

		return true;

	}

	return false;

};


/*****************************************************************************/
/*                                                                           */
/*    - Constructing new triangles routines					                 */
/*			: 'triangle_event' - Vertex is inserted into an triangle		 */
/*			: 'edge_event' - Vertex is inserted into an edge				 */
/*																			 */
/*****************************************************************************/
template<GEOMETRIC_KERNEL GK, CHECK_KERNEL CK>
void CDT<GK, CK>::triangle_event(Triangle * t, Vertex * p) {


	Vertex * const v0 = t->vertex(0);
	Vertex * const v1 = t->vertex(1);
	Vertex * const v2 = t->vertex(2);

	Triangle * const n0 = t->neighbor(0);
	Triangle * const n1 = t->neighbor(1);
	Triangle * const n2 = t->neighbor(2);

	Triangle * const t0 = new Triangle(p, v1, v2);
	Triangle * const t1 = new Triangle(p, v2, v0);
	Triangle * const t2 = new Triangle(p, v0, v1);

	Edge * const e0 = new Edge(p, v0);
	Edge * const e1 = new Edge(p, v1);
	Edge * const e2 = new Edge(p, v2);

	// Assign edges to triangles
	t0->set_edge(0) = t->edge(0);
	t0->set_edge(1) = e2;
	t0->set_edge(2) = e1;

	t1->set_edge(0) = t->edge(1);
	t1->set_edge(1) = e0;
	t1->set_edge(2) = e2;

	t2->set_edge(0) = t->edge(2);
	t2->set_edge(1) = e1;
	t2->set_edge(2) = e0;

	// Assign neighboring triangles to edges
	t->edge(0)->set_neighbor(0) = t0;
	t->edge(0)->set_neighbor(1) = n0;

	t->edge(1)->set_neighbor(0) = t1;
	t->edge(1)->set_neighbor(1) = n1;

	t->edge(2)->set_neighbor(0) = t2;
	t->edge(2)->set_neighbor(1) = n2;



	/******************************/
	//unsigned const nt = triangles.size();

	//t0->index = t->index;
	//t1->index = nt;
	//t2->index = nt + 1;
	/******************************/

	// Set pointer to any vertex's 'p' adjacent triangle. Also there is a possibility that some vertices have pointers to deleted triangle 't'
	p->adjacent_triangle = t0;

	v0->adjacent_triangle = t1;
	v1->adjacent_triangle = t2;
	v2->adjacent_triangle = t0;


	// Remove and delete splitted triangle consisting of t0 , t1 , t2 from triangulation
	remove(triangles, t);
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
	triangles.push_back(t0);
	triangles.push_back(t1);
	triangles.push_back(t2);

	// Add new edges to triangulation
	edges.push_back(e0);
	edges.push_back(e1);
	edges.push_back(e2);

	// Perform edge-flip algorithm to legalize non-Delaunay edge
	legalize(t0, p);
	legalize(t1, p);
	legalize(t2, p);

};
template<GEOMETRIC_KERNEL GK, CHECK_KERNEL CK>
Triangle * CDT<GK, CK>::edge_event(Edge * e, Vertex * p) {


	// Assumption : all edges are not on the super triangle => Thus every edge have both neighbors non-NULL
	Triangle * const t0 = e->neighbor(0);
	Triangle * const t1 = e->neighbor(1);

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

};


/*****************************************************************************/
/*                                                                           */
/*    - Construction and destruction of 'super triangle'                     */
/*                                                                           */
/*****************************************************************************/
template<GEOMETRIC_KERNEL GK, CHECK_KERNEL CK>
void CDT<GK, CK>::super_triangle_construct() {


	double const mid = 0.0;
	double const range = 10;

	maximum_coordinate = 2 * range;

	// Construct super triangle with initial size
	st_v0 = new Vertex(mid - 20 * maximum_coordinate, mid - maximum_coordinate);
	st_v1 = new Vertex(mid + 20 * maximum_coordinate, mid - maximum_coordinate);
	st_v2 = new Vertex(mid, mid + 20 * maximum_coordinate);

	Triangle * const t = new Triangle(st_v0, st_v1, st_v2);

	// Construct edges of super triangle
	st_e0 = new Edge(st_v1, st_v2);
	st_e1 = new Edge(st_v2, st_v0);
	st_e2 = new Edge(st_v0, st_v1);

	// Assign edges to super triangle
	t->set_edge(0) = st_e0;
	t->set_edge(1) = st_e1;
	t->set_edge(2) = st_e2;

	// Assign neighboring triangles to edges
	st_e0->set_neighbor(0) = t;
	st_e0->set_neighbor(1) = NULL;

	st_e1->set_neighbor(0) = t;
	st_e1->set_neighbor(1) = NULL;

	st_e2->set_neighbor(0) = t;
	st_e2->set_neighbor(1) = NULL;

	st_e0->is_constrained = true;
	st_e1->is_constrained = true;
	st_e2->is_constrained = true;


	// Set neighborhood of super tiangle to NULL
	t->null_neighbors();

	// Add super triangle to initial triangulation
	triangles.push_back(t);

	// Add edges initial triangulation		- DO I NEED TO ADD these EDGES to the edge vector?
	edges.push_back(st_e0);
	edges.push_back(st_e1);
	edges.push_back(st_e2);

	//t->index = 0;

};
template<GEOMETRIC_KERNEL GK, CHECK_KERNEL CK>
void CDT<GK, CK>::super_triangle_remove() {


	Vertex * const a = st_v0;
	Vertex * const b = st_v1;
	Vertex * const c = st_v2;

	// Delete triangles connected to vertices of super triangle
	for (size_t i = 0; i < triangles.size(); i++) {

		Triangle * const t = triangles[i];

		//t->index = i;

		if (t->contains(a) || t->contains(b) || t->contains(c)) {


			// Delete reference to these triangles (from the true triangulation). Also it is possible to denoted edges convex hull as constrained -> commented out for the time being
			if (t->neighbors[0]) {

				unsigned const index0 = t->neighbors[0]->edge_index(t->edges[0]);
				t->neighbors[0]->neighbors[index0] = NULL;

				//t->edges[0]->is_constrained = true;
				//t->edges[0]->marker = E_MARKER::OUTSIDE;

			}
			if (t->neighbors[1]) {

				unsigned const index1 = t->neighbors[1]->edge_index(t->edges[1]);
				t->neighbors[1]->neighbors[index1] = NULL;

				//t->edges[1]->is_constrained = true;
				//t->edges[1]->marker = E_MARKER::OUTSIDE;

			}
			if (t->neighbors[2]) {

				unsigned const index2 = t->neighbors[2]->edge_index(t->edges[2]);
				t->neighbors[2]->neighbors[index2] = NULL;

				//t->edges[2]->is_constrained = true;
				//t->edges[2]->marker = E_MARKER::OUTSIDE;

			}

			remove(triangles, t);
			delete t;

			i--;

		}
	}

	// Delete edges connected to vertices of super triangle
	for (size_t i = 0; i < edges.size(); i++) {

		Edge * const e = edges[i];

		if (e->contains(a) || e->contains(b) || e->contains(c)) {

			remove(edges, e);
			delete e;

			i--;

		}
	}


	delete st_v0;
	delete st_v1;
	delete st_v2;

	num_edges = edges.size();
	num_triangles = triangles.size();


	//triangles.erase(std::remove_if(begin(triangles), end(triangles),
	//							   [a, b, c](Triangle * t) {return t->contains(a) || t->contains(b) || t->contains(c); }),
	//			    end(triangles));
	//
	// Delete edges connected to vertices of super triangle
	//edges.erase(std::remove_if(begin(edges), end(edges),
	//						   [a, b, c](Edge * e) {return e->contains(a) || e->contains(b) || e->contains(c); }),
	//			end(edges));
	//
	//delete st_e0;
	//delete st_e1;
	//delete st_e2;

};


/*****************************************************************************/
/*                                                                           */
/*    - Getters for triangulation primitives			                     */
/*                                                                           */
/*****************************************************************************/
template<GEOMETRIC_KERNEL GK, CHECK_KERNEL CK>
std::vector<Triangle> CDT<GK, CK>::get_triangles() {

	std::vector<Triangle> temp;

	for (size_t i = 0; i < triangles.size(); i++)
		temp.push_back(*triangles[i]);

	return temp;

};
template<GEOMETRIC_KERNEL GK, CHECK_KERNEL CK>
std::vector<Edge> CDT<GK, CK>::get_edges() {

	std::vector<Edge> temp;

	for (size_t i = 0; i < edges.size(); i++)
		temp.push_back(*edges[i]);

	return temp;

};
template<GEOMETRIC_KERNEL GK, CHECK_KERNEL CK>
std::vector<Vertex> CDT<GK, CK>::get_new_vertices() {

	std::vector<Vertex> temp;

	for (size_t i = 0; i < constraints_new_vertices.size(); i++)
		temp.push_back(*constraints_new_vertices[i]);

	return temp;

};


template<GEOMETRIC_KERNEL GK, CHECK_KERNEL CK>
std::vector<Vertex_handle> CDT<GK, CK>::get_vertices_handles() {

	return vertices;

};
template<GEOMETRIC_KERNEL GK, CHECK_KERNEL CK>
std::vector<Triangle_handle> CDT<GK, CK>::get_triangles_handles() {

	return triangles;

};
template<GEOMETRIC_KERNEL GK, CHECK_KERNEL CK>
std::vector<Edge_handle> CDT<GK, CK>::get_edges_handles() {

	return edges;

};
template<GEOMETRIC_KERNEL GK, CHECK_KERNEL CK>
std::vector<Vertex_handle> CDT<GK, CK>::get_new_vertices_handles() {

	return constraints_new_vertices;

};


/*****************************************************************************/
/*                                                                           */
/*    - Writing vertices / triangles / edges into 'txt' files                */
/*                                                                           */
/*****************************************************************************/
template<GEOMETRIC_KERNEL GK, CHECK_KERNEL CK>
void CDT<GK, CK>::export_vertices(std::ofstream & stream) const {

	size_t const n = vertices.size();

	for (size_t i = 0; i < n; i++) {

		double const v_x = vertices[i]->x;
		double const v_y = vertices[i]->y;

		stream << v_x << "  " << v_y << std::endl;

	}

};
template<GEOMETRIC_KERNEL GK, CHECK_KERNEL CK>
void CDT<GK, CK>::export_edges(std::ofstream & stream) const {

	size_t const n = edges.size();


	// Gnuplot
	for (size_t i = 0; i < n; i++) {

		//if (edges[i]->is_constrained) {
			double const v0_x = edges[i]->getA()->x;
			double const v0_y = edges[i]->getA()->y;

			double const v1_x = edges[i]->getB()->x;
			double const v1_y = edges[i]->getB()->y;

			stream << v0_x << "  " << v0_y << std::endl;
			stream << v1_x << "  " << v1_y << std::endl << std::endl;
		//}

	}

	// Matlab
	//for (size i = 0; i < n; i++) {
	//
	//	double const v0_x = edges[i]->getA()->x;
	//	double const v0_y = edges[i]->getA()->y;
	//
	//	double const v1_x = edges[i]->getB()->x;
	//	double const v1_y = edges[i]->getB()->y;
	//
	//	stream << v0_x << "  " << v0_y << "  " << v1_x << "  " << v1_y << std::endl;
	//
	//}

};
template<GEOMETRIC_KERNEL GK, CHECK_KERNEL CK>
void CDT<GK, CK>::export_triangles(std::ofstream & stream) const {

	size_t const n = triangles.size();

	for (size_t i = 0; i < n; i++) {

		double const v0_x = triangles[i]->vertex(0)->x;
		double const v0_y = triangles[i]->vertex(0)->y;

		double const v1_x = triangles[i]->vertex(1)->x;
		double const v1_y = triangles[i]->vertex(1)->y;

		double const v2_x = triangles[i]->vertex(2)->x;
		double const v2_y = triangles[i]->vertex(2)->y;

		stream << v0_x << "  " << v0_y << std::endl;
		stream << v1_x << "  " << v1_y << std::endl;
		stream << v2_x << "  " << v2_y << std::endl;
		stream << v0_x << "  " << v0_y << std::endl << std::endl;

	}

};


/*****************************************************************************/
/*                                                                           */
/*    - 'remove' : removes triangle / edge / vertex							 */
/*				   from triangulation data structure						 */
/*    - 'clear'  : completely clear all data								 */
/*                                                                           */
/*****************************************************************************/
template<GEOMETRIC_KERNEL GK, CHECK_KERNEL CK>
template<typename T>
void CDT<GK, CK>::remove(std::vector<T> & vec, T t) {


	// It is much faster to lookup for an object 't' from the back
	for (auto it = vec.rbegin(); it != vec.rend(); it++) {


		// If an iterator is equal object 't' it is erased from vector and the search is finished
		if (*it == t) {

			vec.erase(it.base() - 1);
			break;

		}
	}

};
template<GEOMETRIC_KERNEL GK, CHECK_KERNEL CK>
void CDT<GK, CK>::clear() {


	Vertex	 * v = NULL;
	Edge	 * e = NULL;
	Triangle * t = NULL;


	// First, free allocated memory
	for (size_t i = 0; i < vertices.size(); i++) {

		v = vertices[i];
		delete v;

	}
	for (size_t i = 0; i < edges.size(); i++) {

		e = edges[i];
		delete e;

	}
	for (size_t i = 0; i < triangles.size(); i++) {

		t = triangles[i];
		delete t;

	}

	// Secodnly, clear the vectors of the rubbish
	vertices.clear();
	edges.clear();
	triangles.clear();

	//seeds.clear();
	constraints_new_vertices.clear();			// Vertices from here are contained in the vector 'vertices' also -> these vertices are therefore already deleted

	mesh_quality_angle.clear();
	mesh_quality_radius.clear();
	mesh_quality_diameter.clear();


	// So we can do new triangulation (super triangle is deleted in above lines)
	super_triangle_construct();

};


/*****************************************************************************/
/*                                                                           */
/*    - Input constrained edge using										 */
/*			: subdivision algorithm - if inserted edge intersects other edge */
/*									  intersections vertex is computed and	 */
/*									  is inserted into triangulation		 */
/*								    - iteratively we obtain constrained		 */
/*									  edge consisting of many smaller edges	 */
/*																			 */
/*****************************************************************************/
template<GEOMETRIC_KERNEL GK, CHECK_KERNEL CK> template<E_MARKER marker>
void CDT<GK, CK>::insert_constraint_intersection(Vertex * const a, Vertex * const b, Triangle * alpha) {



	// Search for the triangle 't' which contains 'a' and the constrained edge goes through the edge which is formed by the rest 2 vertices 'vertex_ccw(a)' <-> 'vertex_cw(a)'. If the edge is already there, NULL is returned
	Triangle * const t = edge_location<marker>(a, b, alpha);

	if (!t)
		return;


	// Constrained edge is situated between these two points
	Vertex * const r = t->vertex_ccw(a);
	Vertex * const l = t->vertex_cw(a);


	if (on_edge_fast(a, b, r)) {

		// Constrained edge goes through ccw vertex 'r' with respect to 'a', and so, part of constrained edge is already in there
		t->edge(a, r)->marker = marker;
		t->edge(a, r)->is_constrained = true;

		// Denote Vertices of constrained edge as 'constrained' -> When smoothing, these vertices won't be displaced
		a->is_constrained = true;
		r->is_constrained = true;

		// The rest of the c.edge must start at the vertex 'r' and end in 'b'. Therefore, it is best to lookup for the triangle containg 'r' from the opposite triangle 'ot'
		insert_constraint_intersection<marker>(r, b, t);

		return;

	}
	else if (on_edge_fast(a, b, l)) {

		// Constrained edge goes through cw vertex 'l' with respect to 'a', and so, part of constrained edge is already in there
		t->edge(a, l)->marker = marker;
		t->edge(a, l)->is_constrained = true;

		// Denote Vertices of constrained edge as 'constrained' -> When smoothing, these vertices won't be displaced
		a->is_constrained = true;
		l->is_constrained = true;

		// The rest of the c.edge must start at the vertex 'l' and end in 'b'. Therefore, it is best to lookup for the triangle containg 'l' from the opposite triangle 'ot'
		insert_constraint_intersection<marker>(l, b, t);

		return;

	}
	
	
	// An edge between 'r' and 'l' vertices
	Edge * const e = t->edge(r, l);

	// Opposite triangle
	Triangle * ot = t->opposite_triangle(a);// neighbor(r, l);
	Vertex * const s = ot->vertex_but(r, l);

	// Edge intersection coordinates
	double ix;
	double iy;

	
	if (s == b) {

		// Opposite vertex 's' in triangle 'ot' IS the ending vertex 'b' 
		get_edge_intersection(a, b, l, r, ix, iy);

		Vertex * const new_vertex = new Vertex(ix, iy);

		constraints_new_vertices.push_back(new_vertex);


		// Insert new_vertex defined by intersection of edges and Get any triangle containing 'new_vertex'
		ot = insert_vertex_into_edge(new_vertex, e);


		// 'new_vertex' obtains its index in the 'insert_vertex_into_edge' function. 
		new_vertex->is_constrained = true;


		// Rotary traverse through triangles around 'new_vertex', so that the resulting triangle contains 'a'. 
		while (true) {

			if (!ot->contains(a))	ot = ot->neighbor_ccw(new_vertex);
			else					break;

		}

		// From the resulting triangle 'ot' we are able to to set the type 'marker' of boundary and to set the part of constrained edge to be constrained
		ot->edge(new_vertex, a)->is_constrained = true;
		ot->edge(new_vertex, a)->marker = marker;

		// Denote Vertices of constrained edge as 'constrained' -> When smoothing, these vertices won't be displaced
		a->is_constrained = true;
		new_vertex->is_constrained = true;


		// The same here, but for the ending vertex 'b'
		while (true) {

			if (!ot->contains(s))	ot = ot->neighbor_ccw(new_vertex);
			else					break;

		}

		ot->edge(new_vertex, s)->is_constrained = true;
		ot->edge(new_vertex, s)->marker = marker;

		// Denote Vertices of constrained edge as 'constrained' -> When smoothing, these vertices won't be displaced
		s->is_constrained = true;

		// Now we are done. No need for further introduction of additional vertices
		return;

	}
	else if (on_edge_fast(a, b, s)) {


		get_edge_intersection(a, b, l, r, ix, iy);

		Vertex * const new_vertex = new Vertex(ix, iy);

		constraints_new_vertices.push_back(new_vertex);


		ot = insert_vertex_into_edge(new_vertex, e);

		new_vertex->is_constrained = true;


		while (true) {

			if (!ot->contains(a))	ot = ot->neighbor_ccw(new_vertex);
			else					break;

		}

		ot->edge(new_vertex, a)->is_constrained = true;
		ot->edge(new_vertex, a)->marker = marker;

		// Denote Vertices of constrained edge as 'constrained' -> When smoothing, these vertices won't be displaced
		a->is_constrained = true;
		new_vertex->is_constrained = true;


		while (true) {

			if (!ot->contains(s))	ot = ot->neighbor_ccw(new_vertex);
			else					break;

		}

		ot->edge(new_vertex, s)->is_constrained = true;
		ot->edge(new_vertex, s)->marker = marker;

		// Denote Vertices of constrained edge as 'constrained' -> When smoothing, these vertices won't be displaced
		s->is_constrained = true;


		// Contrary to the last possibility 's == b', the c.edge passes through the opposite vertex 's', but vertex 's' is not an ending point. Therefore further subdivision is needed
		insert_constraint_intersection<marker>(s, b, ot);

	}
	else {

		get_edge_intersection(a, b, l, r, ix, iy);

		Vertex * const new_vertex = new Vertex(ix, iy);

		constraints_new_vertices.push_back(new_vertex);

		ot = insert_vertex_into_edge(new_vertex, e);

		new_vertex->is_constrained = true;

		while (true) {

			if (!ot->contains(a))	ot = ot->neighbor_ccw(new_vertex);
			else					break;

		}

		ot->edge(new_vertex, a)->is_constrained = true;
		ot->edge(new_vertex, a)->marker = marker;

		// Denote Vertices of constrained edge as 'constrained' -> When smoothing, these vertices won't be displaced
		a->is_constrained = true;

		// This is typical situation, when the c.edge doesn't pass through any vertex but the 'new_vertex'
		insert_constraint_intersection<marker>(new_vertex, b, ot);



	}

};

template<GEOMETRIC_KERNEL GK, CHECK_KERNEL CK> template<E_MARKER marker>
void CDT<GK, CK>::insert_constraint_halving(Vertex * const a, Vertex * const b, Triangle * alpha) {


	// Serious problems with this method !!


	Triangle_handle const t = edge_location<marker>(a, b, a->adjacent_triangle);

	if (!t)
		return;


	double const x_mid = 0.5*(a->x + b->x);
	double const y_mid = 0.5*(a->y + b->y);

	Vertex_handle new_vertex = new Vertex(x_mid, y_mid);


	// If 'new_vertex' is duplicate, then 'insert_vertex(new_vertex)' returns the duplicate Vertex_handle, therefore the return value isn't the same as handle to 'new_vertex'
	Vertex_handle possible_duplicate = insert_vertex(new_vertex);

	// 'new_vertex' is discarded and it's coordinates (pretty much the same, but to ensure, that the edge is straight) are copied to the existing vertex
	if (new_vertex != possible_duplicate) {

		possible_duplicate->x = x_mid;
		possible_duplicate->y = y_mid;

		// free allocated memory
		delete new_vertex;

		new_vertex = possible_duplicate;
	
	}

	new_vertex->is_constrained = true;

	constraints_new_vertices.push_back(new_vertex);


	Triangle_handle const ot = edge_location<marker>(a, new_vertex, a->adjacent_triangle);
	Triangle_handle const ot2 = edge_location<marker>(new_vertex, b, new_vertex->adjacent_triangle);


	if (!ot && !ot2) {

		return;

	}

	//// If an edge through which the constrained edge have to go is also constrained, it can't be flipped and therefore intersection algorithm must be used. That one will insert vertex into blocking constrained edge
	//if (ot->edge(ot->vertex_index(a))->is_constrained) {
	//
	//	cout << "Halving algorithm encountered constrained edge. Intersection algorithm is used instead." << endl;
	//	insert_constraint_intersection<marker>(a, b, ot);
	//
	//	return;
	//
	//}

	if (ot) {

		 Vertex * const r = ot->vertex_ccw(a);
		 Vertex * const l = ot->vertex_cw(a);

		 if (on_edge_fast(a, new_vertex, l)) {

			 cout << "shit1" << endl;

			 ot->edge(a, l)->is_constrained = true;
			 ot->edge(a, l)->marker = marker;

			 l->is_constrained = true;

			 constraints_new_vertices.push_back(l);

			 insert_constraint_halving<marker>(l, b, ot);

			 return;

		 }
		 else if (on_edge_fast(a, new_vertex, r)) {

			 cout << "shit2" << endl;

			 ot->edge(a, r)->is_constrained = true;
			 ot->edge(a, r)->marker = marker;

			 r->is_constrained = true;

			 constraints_new_vertices.push_back(r);

			 insert_constraint_halving<marker>(r, b, ot);

			 return;

		 }

		insert_constraint_halving<marker>(a, new_vertex, ot);
		return;
	}
	if (ot2) {

		 Vertex * const r = ot2->vertex_ccw(new_vertex);
		 Vertex * const l = ot2->vertex_cw(new_vertex);

		 if (on_edge_fast(new_vertex, b, l)) {

			 cout << "shit1" << endl;

			 ot2->edge(new_vertex, l)->is_constrained = true;
			 ot2->edge(new_vertex, l)->marker = marker;

			 l->is_constrained = true;
		
			 constraints_new_vertices.push_back(l);

			 insert_constraint_halving<marker>(l, b, ot2);

			 return;

		 }
		 else if (on_edge_fast(new_vertex, b, r)) {

			 cout << "shit2" << endl;

			 ot2->edge(new_vertex, r)->is_constrained = true;
			 ot2->edge(new_vertex, r)->marker = marker;

			 r->is_constrained = true;
			
			 constraints_new_vertices.push_back(r);

			 insert_constraint_halving<marker>(r, b, ot2);

			 return;

		 }

		insert_constraint_halving<marker>(new_vertex, b, ot2);
		return;

	}

};
template<GEOMETRIC_KERNEL GK, CHECK_KERNEL CK> template<E_MARKER marker>
void CDT<GK, CK>::insert_constraint_flip(Vertex * const a, Vertex * const b, Triangle * alpha) {


	Triangle_handle const t = edge_location<marker>(a, b, a->adjacent_triangle);

	if (!t)
		return;



	// Have to implement a check whether two adjacent triangles form convex quadtrilateral or not ! -> possibly the easiest is 'Gift wrapping algorithm' also called 'Jarvis marrch'

};


/*****************************************************************************/
/*                                                                           */
/*    - Finds an starting vertex 'a' from which the constrained edge is		 */
/*		heading																 */
/*																			 */
/*****************************************************************************/
template<GEOMETRIC_KERNEL GK, CHECK_KERNEL CK> template<E_MARKER marker>
Triangle * CDT<GK, CK>::edge_location(Vertex * const a, Vertex * const b, Triangle * const alpha) {


	/*****************************************************************************/
	/*                                                                           */
	/*    - Locate triangle containg vertex 'a'                                  */
	/*                                                                           */
	/*****************************************************************************/
	//Triangle * t = alpha ? alpha : locate_triangle_straight_walk(triangles.back(), a);
	//Triangle * t = locate_triangle_lawson_walk(triangles.back(), a);
	Triangle * t = alpha == NULL ? locate_triangle_lawson_walk_remembering(triangles.back(), a) : alpha;


	/*****************************************************************************/
	/*                                                                           */
	/*    - For Debugging purpose :)                                             */
	/*                                                                           */
	/*****************************************************************************/
	if (!predicates.in_triangle_robust(t, a))
		std::cout << "'a' not in Triangle 't'" << std::endl;


	/*****************************************************************************/
	/*                                                                           */
	/*    - Circulate through triangles around the point 'a'					 */
	/*		so we can check if the edge exists, or to get the direction			 */
	/*		to the vertice 'b'													 */
	/*																			 */
	/*****************************************************************************/
	Vertex * const A = a;
	Vertex * L = t->vertex_cw(a);
	Vertex * R = t->vertex_ccw(a);

	// Do it with rotary traversal... it may be more efficient

	if (orientation_robust(R, A, b) <= 0.0) {

		while (orientation_robust(L, A, b) < 0.0) {

			R = L;
			t = t->neighbor(A, L);
			L = t->vertex_but(A, R);

		}
	}
	else {

		do {

			L = R;
			t = t->neighbor(A, R);
			R = t->vertex_but(A, L);

		} while (orientation_robust(R, A, b) > 0.0);
	}

	/*****************************************************************************/
	/*                                                                           */
	/*    - Now 't' is the starting triangle. Check its edges if there			 */
	/*		exists edge 'ab'													 */
	/*																			 */
	/*****************************************************************************/
	if (L == b || R == b) {

		std::cout << "Edge already exists. Skipping." << std::endl;

		t->edge(a, b)->is_constrained = true;
		t->edge(a, b)->marker = marker;

		return NULL;

	}


	return t;

};


/*****************************************************************************/
/*                                                                           */
/*    - 'make_individual_hole' : delete all triangles inside constrained	 */
/*								 polygon denoted by 'seed' vertex			 */
/*							   : triangle containing 'seed' is located		 */
/*								 by testing all triangles in triangulation   */
/*																			 */
/*																			 */
/*****************************************************************************/
template<GEOMETRIC_KERNEL GK, CHECK_KERNEL CK>
void CDT<GK, CK>::make_individual_hole(Vertex * const seed) {


	Triangle * t = NULL;

	for (size_t i = 0; i < triangles.size(); i++) {

		if (predicates.in_triangle(triangles[i], seed)) {

			t = triangles[i];
			break;

		}
	}


	if (!t)
		return;

	remove_triangles_in_hole(seed, t);

	num_vertices = vertices.size();
	num_edges = edges.size();
	num_triangles = triangles.size();

};
template<GEOMETRIC_KERNEL GK, CHECK_KERNEL CK>
void CDT<GK, CK>::make_individual_hole(Vertex seed) {


	Triangle * t = NULL;

	for (size_t i = 0; i < triangles.size(); i++) {

		if (predicates.in_triangle(triangles[i], &seed)) {

			t = triangles[i];
			break;

		}
	}

			
	if (!t)
		return;

	remove_triangles_in_hole(&seed, t);

	num_vertices = vertices.size();
	num_edges = edges.size();
	num_triangles = triangles.size();

};


/*****************************************************************************/
/*                                                                           */
/*    - 'make_hole' : delete all triangles inside constrained polygons		 */
/*                    which were created with 'insert_sequential_constraint' */
/*					: It is called after ! ALL ! closed						 */
/*					  polygons are inserted									 */
/*					: triangles containing 'seed's are located via		     */
/*					   walking algorithms									 */
/*																			 */
/*																			 */
/*****************************************************************************/
template<GEOMETRIC_KERNEL GK, CHECK_KERNEL CK> template<class input_iterator>
void CDT<GK, CK>::make_holes(input_iterator first_seed, input_iterator last_seed) {


	// Finds all initial triangles in all closed polygons.
	for (auto it = first_seed; it != last_seed; it++)
		make_individual_hole((*it));

		/*

		std::vector<Triangle *> seed_triangles;


	for (size_t i = 0; i < seeds.size(); i++)
		seed_triangles.push_back(locate_triangle_straight_walk(triangles.back(), seeds[i]));


	// Make holes in all closed polygons. (Spreading of a virus)
	for (size_t i = 0; i < seed_triangles.size(); i++)
		remove_triangles_in_hole(seeds[i], seed_triangles[i]);


	 'seed_triangles' are already destroyed in the proccess of digging a hole. Therefore do not 'delete' them !
	for (size_t i = 0; i < seed_triangles.size(); i++) {

		remove(seed_triangles, seed_triangles[i]);
		remove(seeds, seeds[i]);

		i--;

	}

	 */

};


/*****************************************************************************/
/*                                                                           */
/*    - 'remove_triangles_in_hole' : delete all triangles inside constrained */
/*									 polygon which are marked as 'OUTSIDE'	 */
/*								   : Any vertex situated inside a hole is	 */
/*									 deleted !!								 */
/*																			 */
/*																			 */
/*****************************************************************************/
template<GEOMETRIC_KERNEL GK, CHECK_KERNEL CK>
void CDT<GK, CK>::remove_triangles_in_hole(Vertex * const v, Triangle * const alpha) {

	
	// Mark triangles in hole as 'OUTSIDE'
	mark_triangle(alpha);


	std::vector<Vertex * > vertices_to_remove;

	// Delete Triangles with 'OUTSIDE' mark
	for (size_t i = 0; i < triangles.size(); i++) {


		Triangle * const t = triangles[i];

		if (t->marker == T_MARKER::OUTSIDE) {

			// Delete edges inside hole, which are not constrained
			for (unsigned j = 0; j < 3; j++) {


				Edge * const e = t->edges[j];

				// If edge is already deleted from other iteration 
				if (!e)
					continue;

				Triangle * const ot = t->neighbors[j];

				Vertex_handle const a = e->A;
				Vertex_handle const b = e->B;

				if (!e->is_constrained) {


					t->edges[j] = NULL;
					ot->edges[ot->edge_index(e)] = NULL;

					remove(edges, e);
					delete e;

					// Vertex inside a hole is pushed into 'vertices_to_remove' container, so that is uniquely present there
					if (!a->is_constrained) {

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
		
					if (!b->is_constrained) {

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
				else {

					// This is connected with 'super_triangle_remove'. Practically it's useless, because the mesh should never contain solo edges. (more like cosmetic fix, so it wont crash)
					// When there are cleared triangles (make_hole, etc.), whose neighbors are triangles connected to super triangle, 
					// there occurs a problem when we chcek for the NULL-ness of neighbors -> it is fixed by this

					if (ot)
						ot->neighbors[ot->edge_index(e)] = NULL;

					// Update new neighborhood adjacency of vertices (otherwise it may point to deleted triangle) !! Think about this little more
					a->adjacent_triangle = ot;
					b->adjacent_triangle = ot;

				}

			}

			remove(triangles, t);
			delete t;

			i--;

		}
	}


	//std::sort(vertices_to_remove.begin(), vertices_to_remove.end());
	//std::vector<Vertex *>::iterator it = std::unique(vertices_to_remove.begin(), vertices_to_remove.end());
	//
	//vertices_to_remove.resize(std::distance(vertices_to_remove.begin(), it));

	// Delete vertices inside hole
	for (size_t i = 0; i < vertices_to_remove.size(); i++) {

		Vertex * const v = vertices_to_remove[i];

		remove(vertices, v);
		delete v;

		num_deleted_vertices++;

	}

	// Re-index vertices
	if (!vertices_to_remove.empty()) {

		for (size_t i = 0; i < vertices.size(); i++)
			vertices[i]->index = i;

	}

	vertices_to_remove.clear();
		
};




/*****************************************************************************/
/*                                                                           */
/*    - 'mark_triangle' : Spreads a virus marking all triangles bounded by	 */
/*						  constrained polygon 'OUTSIDE'						 */
/*																			 */
/*****************************************************************************/
template<GEOMETRIC_KERNEL GK, CHECK_KERNEL CK>
void CDT<GK, CK>::mark_triangle(Triangle * const t) {


	t->marker = T_MARKER::OUTSIDE;

	for (unsigned i = 0; i < 3; i++) {


		// Same purpose as in 'remove_triangles_in_hole'.
		if (!t->neighbors[i])
			continue;

		if (t->neighbors[i]->marker == T_MARKER::OUTSIDE || t->edges[i]->is_constrained)
			continue;
		else
			mark_triangle(t->neighbors[i]);
			

	}

};






/*****************************************************************************/
/*                                                                           */
/*    - 'laplacian_smoothing' : smooth positions of vertices using			 */
/*								Laplace smoothing algorithm					 */
/*																			 */
/*****************************************************************************/
template<GEOMETRIC_KERNEL GK, CHECK_KERNEL CK>
void CDT<GK, CK>::laplacian_smoothing(unsigned const N) {


	// This have to be, otherwise problems with vertices of super_triangle (dependently on the function which lookup for neighbor vertices -> there those vertices may or may not exist)
	std::vector<bool> is_on_convex_hull;

	set_constrained_convex_hull(is_on_convex_hull);


	size_t nv = vertices.size();
	size_t nt = triangles.size();

	std::vector<Vertex> vertices_temp;

	vertices_temp.reserve(nv);

	// Copy initial coordinates of vertices. Allocaate once
	for (size_t i = 0; i < nv; i++) {

		Vertex * const v = vertices[i];

		double const x = v->x;
		double const y = v->y;

		vertices_temp.push_back(Vertex(x, y, i));

	}

	for (unsigned iteration = 0; iteration < N; iteration++) {


		// Smoothing each non-constrained vertex
		for (size_t i = 0; i < nv; i++) {


			Vertex * const v = vertices[i];

			if (v->is_constrained || is_on_convex_hull[i])
				continue;


			double const X = v->x;
			double const Y = v->y;


			unsigned const v_index = v->index;

			std::vector<Triangle *> adjacent_triangles;
			std::vector<Vertex *> adjacent_vertices;

			size_t const n = get_adjacent_triangles_and_vertices(v, adjacent_triangles, adjacent_vertices);
			size_t const nt = adjacent_triangles.size();

			double sx = 0.0;
			double sy = 0.0;


			for (size_t j = 0; j < adjacent_vertices.size(); j++) {

				unsigned const k = adjacent_vertices[j]->index;

				sx += vertices_temp[k].x;
				sy += vertices_temp[k].y;

			}


			double beta = 1.0;

			double x = (1.0 - beta)*X + beta * sx / n;
			double y = (1.0 - beta)*Y + beta * sy / n;

			// New position of vertex 'v' which we want to displace
			Vertex new_position(x, y);
			Vertex * const old_position = & vertices_temp[v_index];

			bool repeat;

			do {

				repeat = false;

				for (size_t k = 0; k < adjacent_vertices.size(); k++) {


					Vertex * const p = adjacent_vertices[k];


					for (size_t j = 0; j < nt; j++) {


						Triangle * const t = adjacent_triangles[j];

						//  Edge 'ab' can't be from initial triangulation (vertices_temp), because then we don't know, if the edge is NOW intersected
						Vertex * const a = t->vertex_cw(v);
						Vertex * const b = t->vertex_ccw(v);


						if (t->contains(p))
							continue;

						
						bool const self_intersection		= predicates.intersecting_edges(a, b, old_position, &new_position);
						bool const neighbor_intersection	= predicates.intersecting_edges(a, b, p, &new_position);

						if (self_intersection || neighbor_intersection) {

							repeat = true;

							beta = 0.5*beta;


							// There must be some problem with intersection / iteration through neighboring vertices and triangles (through k,j) / search for neighboring 'v' and 't'
							if (beta < 1.e-8) {

								cout << "beta under 1.e-8. Should not happen." << endl;
								
								x = X;
								y = Y;

								repeat = false;
								break;

							}
								

							x = (1.0 - beta)*X + beta * sx / n;
							y = (1.0 - beta)*Y + beta * sy / n;

							new_position.x = x;
							new_position.y = y;

							k--;

							break;

						}
						else
							repeat = false;
					}
				}


			} while (repeat);

			// Set new coordinates
			v->x = x;
			v->y = y;

			vertices_temp[v_index].x = x;
			vertices_temp[v_index].y = y;

			//for (size_t i = 0; i < adjacent_triangles.size(); i++)
			//	legalize(adjacent_triangles[i], v);

		}

		// Copy initial coordinates of vertices
		for (size_t i = 0; i < nv; i++) {

			Vertex * const v = vertices[i];

			double const x = v->x;
			double const y = v->y;

			Vertex * const p = &vertices_temp[i];

			p->x = x;
			p->y = y;

		}


		legalize_vertices();

	}


};


/*****************************************************************************/
/*                                                                           */
/*    -	Finds adjacent triangles and vertices neighboring vertex 'v'		 */
/*	  - Dependently on implementation, it either finds all triangles		 */
/*		and vertices, or, it finds all but super_triangle's triangles		 */
/*		and vertices														 */
/*																			 */
/*****************************************************************************/
template<GEOMETRIC_KERNEL GK, CHECK_KERNEL CK>
unsigned CDT<GK, CK>::get_adjacent_triangles_and_vertices(Vertex * const v, std::vector<Triangle *> & triangles, std::vector<Vertex *> & vertices) {


	Triangle * const stop_triangle = v->adjacent_triangle;
	Triangle * t = stop_triangle;

	unsigned nv = 0;


	bool direction_forward = true;

	do {

		triangles.push_back(t);

		//if ((!t->contains(st_v0) && !t->contains(st_v1) && !t->contains(st_v2)) || (t->vertex_cw(v) != st_v0 && t->vertex_cw(v) != st_v1 && t->vertex_cw(v) != st_v2)) {
		if (!t->contains(st_v0) && !t->contains(st_v1) && !t->contains(st_v2)) {

			vertices.push_back(t->vertex_cw(v));
			
			nv++;

		}

		t = t->neighbor_cw(v);

		if (!t)
			direction_forward = false;


	} while (direction_forward && t != stop_triangle);


	// If we have reached NULL triangle (thanks to some hole) we have to change the direction of rotation.
	if (!direction_forward) {

		t = stop_triangle->neighbor_ccw(v);

		while ( t && t != stop_triangle) {

			triangles.push_back(t);

			//if ((!t->contains(st_v0) && !t->contains(st_v1) && !t->contains(st_v2)) || (t->vertex_cw(v) != st_v0 && t->vertex_cw(v) != st_v1 && t->vertex_cw(v) != st_v2)) {
			if (!t->contains(st_v0) && !t->contains(st_v1) && !t->contains(st_v2)) {

				//triangles.push_back(t);
				vertices.push_back(t->vertex_cw(v));

				nv++;

			}

			if (!t->neighbor_ccw(v) && (t->vertex_ccw(v) != st_v0 && t->vertex_ccw(v) != st_v1 && t->vertex_ccw(v) != st_v2)) {

				vertices.push_back(t->vertex_ccw(v));
				nv++;

			}
				

			t = t->neighbor_ccw(v);

		}

	}



	return nv;

};


/*****************************************************************************/
/*                                                                           */
/*    - 'set_constrained_convex_hull' : denoted vertices on the convex hull	 */
/*								        as constrained when using laplace	 */
/*										smoothing							 */
/*																			 */
/*****************************************************************************/
template<GEOMETRIC_KERNEL GK, CHECK_KERNEL CK>
void CDT<GK, CK>::set_constrained_convex_hull(std::vector<bool> & is_on_convex_hull) {


	size_t nt = triangles.size();
	size_t nv = vertices.size();

	for (size_t i = 0; i < nv; i++)
		is_on_convex_hull.push_back(false);


	Triangle * t = st_v2->adjacent_triangle;

	while (!t->contains(st_v1))
		t = t->neighbor_cw(st_v2);


	do {

		t = t->neighbor_cw(st_v1);
		is_on_convex_hull[t->vertex_ccw(st_v1)->index] = true;

	} while (!t->contains(st_v0));

	do {

		t = t->neighbor_cw(st_v0);
		is_on_convex_hull[t->vertex_ccw(st_v0)->index] = true;

	} while (!t->contains(st_v2));

	do {

		t = t->neighbor_cw(st_v2);
		is_on_convex_hull[t->vertex_ccw(st_v2)->index] = true;

	} while (!t->contains(st_v1));


};




















template<GEOMETRIC_KERNEL GK, CHECK_KERNEL CK>
void CDT<GK, CK>::refinement_ruppert(double const Bound) {



	// Get all constrained segments up to now
	std::vector<Edge *> segments;

	for (size_t i = 0; i < edges.size(); i++) {

		Edge * const e = edges[i];

		if (e->is_constrained && e != st_e0 && e != st_e1 && e != st_e2)
			segments.push_back(e);

	}

	// While any segment 'e' is encroached : split_segment
	for (size_t i = 0; i < segments.size(); i++) {


		Edge * const e = segments[i];


		for (size_t j = 0; j < vertices.size(); j++) {

			Vertex * const v = vertices[j];

			if (e->contains(v))
				continue;

			if (is_encroeched(e, v)) {


				Vertex * const A = e->A;
				Vertex * const B = e->B;

				double const a_x = A->x;
				double const a_y = A->y;

				double const b_x = B->x;
				double const b_y = B->y;


				double const mid_x = 0.5*(a_x + b_x);
				double const mid_y = 0.5*(a_y + b_y);

				Vertex * const new_vertex = new Vertex(mid_x, mid_y);

				remove(segments, e);

				Triangle * k = insert_vertex_into_edge(new_vertex, e);

				new_vertex->is_constrained = true;


				while (true) {

					if (!k->contains(A))	k = k->neighbor_ccw(new_vertex);
					else					break;

				}

				segments.push_back(k->edge(new_vertex, A));

				while (true) {

					if (!k->contains(B))	k = k->neighbor_ccw(new_vertex);
					else					break;

				}

				segments.push_back(k->edge(new_vertex, B));


				// It is sufficient to start over with the next edge after this deleted edge 'e' -> because it is deleted, the iterator 'i' must be  lowered by 1
				i--;

				break;

			}
		}
	}


	double const L = shortest_edge_length();

	//double const BOUND = 1.0 / (2.0 * sin(Bound));
	//double const BOUND = 1.0;
	//cout << "Bound B : " << BOUND << endl;


	for (size_t i = 0; i < triangles.size(); i++) {



		Triangle * const t = triangles[i];

		if (t->contains(st_v0) || t->contains(st_v1) || t->contains(st_v2))
			continue;


		double x_c;
		double y_c;

		double r;
		double shortest_edge;

		find_circumcenter(t->vertices[0], t->vertices[1], t->vertices[2], x_c, y_c, r, shortest_edge);

		double const ratio = r / shortest_edge;		// Makes non-uniform size of triangles
		//double const ratio = r / L;					// Makes uniform size of triangles

		//double const B = sqr(BOUND);
		//double const B = sqr(sqrt(2.0));
		double const B = 0.69;// sqr(1.0);

		// radius 'r' and length of the shortest edge 'L' are squared -> postive numbers and exponent > 0, therefore inequalities are conserved
		if (ratio > B) {


			std::vector<Edge * > encroached_segments;
			Vertex circum_center(x_c, y_c);


			for (size_t j = 0; j < segments.size(); j++) {


				Edge * const e = segments[j];

				if (is_encroeched(e, &circum_center))
					encroached_segments.push_back(e);

			}

			if (!encroached_segments.empty()) {

				for (size_t j = 0; j < encroached_segments.size(); j++) {


					Edge * const e = encroached_segments[j];

					Vertex * const A = e->A;
					Vertex * const B = e->B;

					double const a_x = A->x;
					double const a_y = A->y;

					double const b_x = B->x;
					double const b_y = B->y;


					double const mid_x = 0.5*(a_x + b_x);
					double const mid_y = 0.5*(a_y + b_y);

					Vertex * const new_vertex = new Vertex(mid_x, mid_y);

					remove(segments, e);

					Triangle * k = insert_vertex_into_edge(new_vertex, e);

					new_vertex->is_constrained = true;


					while (true) {

						if (!k->contains(A))	k = k->neighbor_ccw(new_vertex);
						else					break;

					}

					segments.push_back(k->edge(new_vertex, A));

					while (true) {

						if (!k->contains(B))	k = k->neighbor_ccw(new_vertex);
						else					break;

					}

					segments.push_back(k->edge(new_vertex, B));

				}
				/*
				// =================================================================================
					// While any segment 'e' is encroached : split_segment
				for (size_t k = 0; k < segments.size(); k++) {


					Edge * const e = segments[k];


					for (size_t j = 0; j < vertices.size(); j++) {

						Vertex * const v = vertices[j];

						if (e->contains(v))
							continue;

						if (is_encroeched(e, v)) {


							Vertex * const A = e->A;
							Vertex * const B = e->B;

							double const a_x = A->x;
							double const a_y = A->y;

							double const b_x = B->x;
							double const b_y = B->y;


							double const mid_x = 0.5*(a_x + b_x);
							double const mid_y = 0.5*(a_y + b_y);

							Vertex * const new_vertex = new Vertex(mid_x, mid_y);

							remove(segments, e);

							Triangle * k = insert_vertex_into_edge(new_vertex, e);

							new_vertex->is_constrained = true;


							while (true) {

								if (!k->contains(A))	k = k->neighbor_ccw(new_vertex);
								else					break;

							}

							segments.push_back(k->edge(new_vertex, A));

							while (true) {

								if (!k->contains(B))	k = k->neighbor_ccw(new_vertex);
								else					break;

							}

							segments.push_back(k->edge(new_vertex, B));


							// It is sufficient to start over with the next edge after this deleted edge 'e' -> because it is deleted, the iterator 'i' must be  lowered by 1
							k--;

							break;

						}
					}
				}*/
				// =================================================================================

			}
			else {

				Vertex * const new_vertex = new Vertex(x_c, y_c);

				Vertex * const possible_duplicate = input_vertex(new_vertex, t); // fastest way is to begin from this triangle whose circumcenter i add ?

				//if (possible_duplicate != new_vertex)
				//	break;


			}
		}
	}




};