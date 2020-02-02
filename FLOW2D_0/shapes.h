#pragma once


#include "enumerators.h"

#include <vector>


class Vertex;
class Edge;
class Triangle;


class Vertex {


	friend class Triangle;

	friend class Edge;



	friend class Mesh;

	friend class PlanarStraightLineGraph;

	template<GEOMETRIC_KERNEL GK>
	friend class Triangulation;

	template <GEOMETRIC_KERNEL GK>
	friend class Predicates;

	/*****************************************/
	/*										 */
	/*				Friend functions		 */
	/*										 */
	/*****************************************/
	friend bool intersecting_edges(Edge * e, Edge * f);
	friend bool intersecting_edges(Vertex * p1, Vertex * p2, Vertex * q1, Vertex * q2);
	friend bool get_edge_intersection(Vertex * p1, Vertex * p2, Vertex * p3, Vertex * p4, double & x, double & y);
	friend bool in_circle_fast(Triangle * t, Vertex * v);
	friend bool in_circle_robust(Triangle * t, Vertex * v);
	friend bool on_edge_fast(Vertex * a, Vertex * b, Vertex * v);
	friend bool on_edge_fast(Edge * e, Vertex * v);
	friend bool on_edge_robust(Vertex * a, Vertex * b, Vertex * v);
	friend bool on_edge_robust(Edge * e, Vertex * v);
	friend double orientation_fast(Vertex * a, Vertex * b, Vertex * v);
	friend double orientation_robust(Vertex * a, Vertex * b, Vertex * v);
	friend double orientationAdapt(Vertex * a, Vertex * b, Vertex * v, const double detsum);
	friend bool in_circleAdapt(Triangle * t, Vertex * v, const double permanent);
	friend void get_mid_point(Edge * e, double &x_mid, double &y_mid);
	friend bool is_encroached(Vertex * const a, Vertex * const b, Vertex * const v);


	friend bool t_compare_x(Triangle * const t1, Triangle * const t2);
	friend bool t_compare_y(Triangle * const t1, Triangle * const t2);
	friend bool v_compare_x(Vertex * const v1, Vertex * const v2);
	friend bool v_compare_y(Vertex * const v1, Vertex * const v2);
	friend bool e_compare_x(Edge * const e1, Edge * const e2);
	friend bool e_compare_y(Edge * const e1, Edge * const e2);


public:
//private:


	/*****************************************/
	/*										 */
	/*				Data members			 */
	/*										 */
	/*****************************************/
	double x;
	double y;

	int index = -1;

	// Pointer to unspecified adjacent triangle
	Triangle * adjacent_triangle = NULL;

	// Marker denoting if the vertex is Constrained (can't be moved or deleted) or Free (can be moved e.g. by Laplacian smoothing or deleted e.g. by some refinement method)
	V_MARKER marker = V_MARKER::FREE;



	/*****************************************/
	/*										 */
	/*				Methods					 */
	/*										 */
	/*****************************************/
	void set(double X, double Y);
	bool is_almost_equal(Vertex * const v);
	//Vertex operator+(Vertex * const v);
	//Vertex & operator=(Vertex * const v);


public:


	Vertex(double X, double Y);
	~Vertex();

	
};

class Edge {


	friend class Vertex;

	friend class Triangle;



	friend class Mesh;

	template<GEOMETRIC_KERNEL GK>
	friend class Triangulation;

	template <GEOMETRIC_KERNEL GK>
	friend class Predicates;

	/*****************************************/
	/*										 */
	/*				Friend functions		 */
	/*										 */
	/*****************************************/
	friend bool intersecting_edges(Edge * e, Edge * f);
	friend bool on_edge_fast(Vertex * a, Vertex * b, Vertex * v);
	friend bool on_edge_fast(Edge * e, Vertex * v);
	friend bool on_edge_robust(Edge * e, Vertex * v);
	friend bool is_encroached(Vertex * const a, Vertex * const b, Vertex * const v);
	friend void get_mid_point(Edge * e, double &x_mid, double &y_mid);


	friend bool e_compare_x(Edge * const e1, Edge * const e2);
	friend bool e_compare_y(Edge * const e1, Edge * const e2);
	friend bool marker_compare_neumann(Edge * const e1, Edge * const e2);
	friend bool marker_compare_dirichlet(Edge * const e1, Edge * const e2);

public:
//private:


	/*****************************************/
	/*										 */
	/*				Data members			 */
	/*										 */
	/*****************************************/
	Vertex * a;
	Vertex * b;

	int index = -1;

	// Adjacent triangles on each of the edge
	Triangle * neighbors[2];

	// Marker denoting if there is some boundary condition
	E_MARKER marker = E_MARKER::NONE;

	// If the edge is constrained, then it can't be flipped
	bool is_constrained = false;
		

	Edge(Vertex * const A, Vertex * const B);
	~Edge();

	
	/*****************************************/
	/*										 */
	/*				Methods					 */
	/*										 */
	/*****************************************/
	// Set the edge's ending vertices
	void set(Vertex * const A, Vertex * const B);

	// Get adjacent triangle on the side 'i' (i = 0/1)
	Triangle * neighbor(unsigned i);

	// Set edge's adjacent triangle on the side 'i'
	Triangle *& set_neighbor(unsigned i);

	// Check if the edge contains vertex 'v'
	bool contains(Vertex * const v);

	// Check if the edge consists of vertices 'A' and 'B'
	bool contains(Vertex * const A, Vertex * const B);

	// Get length of this edge
	double length();

	Vertex * get_shared_vertex(Edge * e);


};

class Triangle {

	friend class Vertex;

	friend class Edge;

	friend class Mesh;
	
	template <GEOMETRIC_KERNEL GK>
	friend class Triangulation;

	template <GEOMETRIC_KERNEL GK>
	friend class Predicates;

	/*****************************************/
	/*										 */
	/*				Friend functions		 */
	/*										 */
	/*****************************************/
	friend bool in_triangle_fast(Triangle * t, Vertex * v);
	friend bool in_triangle_robust(Triangle * t, Vertex * v);
	friend bool in_circle_fast(Triangle * t, Vertex * v);
	friend bool in_circle_robust(Triangle * t, Vertex * v);
	friend bool in_circleAdapt(Triangle * t, Vertex * v, const double permanent);


	friend bool t_compare_x(Triangle * const t1, Triangle * const t2);
	friend bool t_compare_y(Triangle * const t1, Triangle * const t2);



public:
//private:


	/*****************************************/
	/*										 */
	/*				Data members			 */
	/*										 */
	/*****************************************/
	Vertex	 * vertices[3];
	Edge	 * edges[3];
	Triangle * neighbors[3];

	int index = -1;


	// Marker denoting if the triangle is Inside the triangulated domain or Outside of the domain and thus is to be deleted
	T_MARKER marker = T_MARKER::NONE;


	Triangle(Vertex * const a, Vertex * const b, Vertex * const c);
	~Triangle();



	/*****************************************/
	/*										 */
	/*				Methods					 */
	/*										 */
	/*****************************************/
	bool contains(Vertex * const v) const;
	bool contains(Vertex * const a, Vertex * const b) const;
	bool contains(Edge * const e) const;


	Vertex * get_vertex(unsigned i) const;
	Vertex * get_vertex(Edge *const  e) const;
	Vertex * get_vertex_cw(Vertex * const v)  const;
	Vertex * get_vertex_ccw(Vertex * const v) const;
	Vertex * get_vertex_but(Vertex * const a, Vertex * const b) const;

	Edge * get_edge(unsigned i) const;
	Edge * get_edge(Vertex * const a, Vertex * const b) const;
	Edge * get_edge_cw(Vertex * const v) const;
	Edge * get_edge_ccw(Vertex * const v) const;

	void set_edge(Edge * const e);
	Edge * & set_edge(int i);

	unsigned get_vertex_index(Vertex * const v) const;
	unsigned get_edge_index(Edge * const e) const;
	unsigned get_edge_index(Vertex * const a, Vertex * const b) const;
	

	Triangle * get_neighbor(unsigned i) const;
	Triangle * get_neighbor(Vertex * const a, Vertex * const b) const;
	Triangle * get_neighbor_cw(Vertex * const v) const;
	Triangle * get_neighbor_ccw(Vertex * const v) const;
	unsigned get_neighbor_index(Triangle * const t) const;

	void set_neighbor(Triangle * const t, unsigned i);
	void set_neighbor(Vertex * const a, Vertex * const b, Triangle * const t);
	void set_neighbor(Edge * const e, Triangle * const t);
	void set_neighbor(Triangle * const t);

	void null_neighbors();


	void rotate_triangle_cw(Vertex * const v, Vertex * const ov);

	Vertex   * get_opposite_vertex(Triangle * const t, Vertex * const v) const;
	Triangle * get_opposite_triangle(Vertex * const v) const;						// Probably useless - the same will do neighbor(Vertex * v) <- to do


	double area() const;
	void circum_center(double & x_c, double & y_c) const;
	double circum_radius_squared() const;
	double shortest_edge_squared() const;
	double ratio_squared() const;

	inline double orientation_determinant() const;

	bool is_bad(double const angle_bound, double const area_bound);

};
