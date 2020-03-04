#pragma once


#include "enumerators.h"

#include <vector>


class Vertex;
class Edge;
class Triangle;

typedef Vertex *		v_pointer;
typedef Edge *			e_pointer;
typedef Triangle *		t_pointer;


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
	friend bool intersecting_edges(e_pointer e, e_pointer f);
	friend bool intersecting_edges(v_pointer p1, v_pointer p2, v_pointer q1, v_pointer q2);
	friend bool get_edge_intersection(v_pointer p1, v_pointer p2, v_pointer p3, v_pointer p4, double & x, double & y);
	friend bool in_circle_fast(t_pointer t, v_pointer v);
	friend bool in_circle_robust(t_pointer t, v_pointer v);
	friend bool on_edge_fast(v_pointer a, v_pointer b, v_pointer v);
	friend bool on_edge_fast(e_pointer e, v_pointer v);
	friend bool on_edge_robust(v_pointer a, v_pointer b, v_pointer v);
	friend bool on_edge_robust(e_pointer e, v_pointer v);
	friend double orientation_fast(v_pointer a, v_pointer b, v_pointer v);
	friend double orientation_robust(v_pointer a, v_pointer b, v_pointer v);
	friend double orientationAdapt(v_pointer a, v_pointer b, v_pointer v, const double detsum);
	friend bool in_circleAdapt(t_pointer t, v_pointer v, const double permanent);
	friend void get_mid_point(e_pointer e, double &x_mid, double &y_mid);
	friend bool is_encroached(v_pointer const a, v_pointer const b, v_pointer const v);


	friend bool t_compare_x(t_pointer const t1, t_pointer const t2);
	friend bool t_compare_y(t_pointer const t1, t_pointer const t2);
	friend bool v_compare_x(v_pointer const v1, v_pointer const v2);
	friend bool v_compare_y(v_pointer const v1, v_pointer const v2);
	friend bool e_compare_x(e_pointer const e1, e_pointer const e2);
	friend bool e_compare_y(e_pointer const e1, e_pointer const e2);


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
	t_pointer adjacent_triangle = NULL;

	// Marker denoting if the vertex is Constrained (can't be moved or deleted) or Free (can be moved e.g. by Laplacian smoothing or deleted e.g. by some refinement method)
	V_MARKER marker = V_MARKER::FREE;



	/*****************************************/
	/*										 */
	/*				Methods					 */
	/*										 */
	/*****************************************/
	void set(double X, double Y);
	bool is_almost_equal(v_pointer const v);
	//Vertex operator+(v_pointer const v);
	//Vertex & operator=(v_pointer const v);


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
	friend bool intersecting_edges(e_pointer e, e_pointer f);
	friend bool on_edge_fast(v_pointer a, v_pointer b, v_pointer v);
	friend bool on_edge_fast(e_pointer e, v_pointer v);
	friend bool on_edge_robust(e_pointer e, v_pointer v);
	friend bool is_encroached(v_pointer const a, v_pointer const b, v_pointer const v);
	friend void get_mid_point(e_pointer e, double &x_mid, double &y_mid);


	friend bool e_compare_x(e_pointer const e1, e_pointer const e2);
	friend bool e_compare_y(e_pointer const e1, e_pointer const e2);
	friend bool marker_compare_neumann(e_pointer const e1, e_pointer const e2);
	friend bool marker_compare_dirichlet(e_pointer const e1, e_pointer const e2);

public:
//private:


	/*****************************************/
	/*										 */
	/*				Data members			 */
	/*										 */
	/*****************************************/
	v_pointer a;
	v_pointer b;

	int index = -1;

	// Adjacent triangles on each of the edge
	t_pointer neighbors[2];

	// Marker denoting if there is some boundary condition
	E_MARKER marker = E_MARKER::NONE;

	// If the edge is constrained, then it can't be flipped
	bool is_constrained = false;
		

	Edge(v_pointer const A, v_pointer const B);
	~Edge();

	
	/*****************************************/
	/*										 */
	/*				Methods					 */
	/*										 */
	/*****************************************/
	// Set the edge's ending vertices
	void set(v_pointer const A, v_pointer const B);

	// Get adjacent triangle on the side 'i' (i = 0/1)
	t_pointer neighbor(unsigned i);

	// Set edge's adjacent triangle on the side 'i'
	t_pointer& set_neighbor(unsigned i);

	// Check if the edge contains vertex 'v'
	bool contains(v_pointer const v);

	// Check if the edge consists of vertices 'A' and 'B'
	bool contains(v_pointer const A, v_pointer const B);

	// Get length of this edge
	double length();

	v_pointer get_shared_vertex(e_pointer e);


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
	friend bool in_triangle_fast(t_pointer t, v_pointer v);
	friend bool in_triangle_robust(t_pointer t, v_pointer v);
	friend bool in_circle_fast(t_pointer t, v_pointer v);
	friend bool in_circle_robust(t_pointer t, v_pointer v);
	friend bool in_circleAdapt(t_pointer t, v_pointer v, const double permanent);


	friend bool t_compare_x(t_pointer const t1, t_pointer const t2);
	friend bool t_compare_y(t_pointer const t1, t_pointer const t2);



public:
//private:


	/*****************************************/
	/*										 */
	/*				Data members			 */
	/*										 */
	/*****************************************/
	v_pointer vertices[3];
	e_pointer edges[3];
	t_pointer neighbors[3];

	int index = -1;


	// Marker denoting if the triangle is Inside the triangulated domain or Outside of the domain and thus is to be deleted
	T_MARKER marker = T_MARKER::NONE;


	Triangle(v_pointer const a, v_pointer const b, v_pointer const c);
	~Triangle();



	/*****************************************/
	/*										 */
	/*				Methods					 */
	/*										 */
	/*****************************************/
	bool contains(v_pointer const v) const;
	bool contains(v_pointer const a, v_pointer const b) const;
	bool contains(e_pointer const e) const;


	v_pointer get_vertex(unsigned i) const;
	v_pointer get_vertex(e_pointer const  e) const;
	v_pointer get_vertex_cw(v_pointer const v)  const;
	v_pointer get_vertex_ccw(v_pointer const v) const;
	v_pointer get_vertex_but(v_pointer const a, v_pointer const b) const;

	e_pointer get_edge(unsigned i) const;
	e_pointer get_edge(v_pointer const a, v_pointer const b) const;
	e_pointer get_edge_cw(v_pointer const v) const;
	e_pointer get_edge_ccw(v_pointer const v) const;

	void set_edge(e_pointer const e);
	e_pointer & set_edge(int i);

	unsigned get_vertex_index(v_pointer const v) const;
	unsigned get_edge_index(e_pointer const e) const;
	unsigned get_edge_index(v_pointer const a, v_pointer const b) const;
	

	t_pointer get_neighbor(unsigned i) const;
	t_pointer get_neighbor(v_pointer const a, v_pointer const b) const;
	t_pointer get_neighbor_cw(v_pointer const v) const;
	t_pointer get_neighbor_ccw(v_pointer const v) const;
	unsigned get_neighbor_index(t_pointer const t) const;

	void set_neighbor(t_pointer const t, unsigned i);
	void set_neighbor(v_pointer const a, v_pointer const b, t_pointer const t);
	void set_neighbor(e_pointer const e, t_pointer const t);
	void set_neighbor(t_pointer const t);

	void null_neighbors();


	void rotate_triangle_cw(v_pointer const v, v_pointer const ov);

	v_pointer get_opposite_vertex(t_pointer const t, v_pointer const v) const;
	t_pointer get_opposite_triangle(v_pointer const v) const;						// Probably useless - the same will do neighbor(v_pointer v) <- to do


	double area() const;
	void circum_center(double & x_c, double & y_c) const;
	double circum_radius_squared() const;
	double shortest_edge_squared() const;
	double ratio_squared() const;

	inline double orientation_determinant() const;

	bool is_bad(double const angle_bound, double const area_bound);

};
