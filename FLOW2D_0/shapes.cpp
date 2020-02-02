

#include "shapes.h"
#include "geometric.h"

#include <cassert>
#include <algorithm>
#include <iostream>


/*****************************************************************************/
/*																			 */
/*									VERTEX									 */
/*																			 */
/*****************************************************************************/
Vertex::Vertex(double X, double Y) {

	x = X;
	y = Y;

};
Vertex::~Vertex() {

};


void Vertex::set(double X, double Y) {

	x = X;
	y = Y;

}
bool Vertex::is_almost_equal(Vertex * const v) {

	//return (abs(x - v->x) < EPSILON_EQUALITY && abs(y - v->y) < EPSILON_EQUALITY);
	// return sqr(x - v->x) + sqr(y - v->y) < sqr(EPSILON_EQUALITY);	// if e.g.  eps < 1.e-8 then sqr(eps) < 1.e-16 -> overflow

	return sqrt(sqr(x - v->x) + sqr(y - v->y)) < EPSILON_EQUALITY;

};
/*
Vertex Vertex::operator+(Vertex * const v) {


	double const ax = x;
	double const ay = y;

	double const vx = v->x;
	double const vy = v->y;

	Vertex result(ax + vx, ay + vy);

	return result;

};
Vertex & Vertex::operator=(Vertex * const v) {

	return *v;

};
*/



/*****************************************************************************/
/*																			 */
/*									EDGE									 */
/*																			 */
/*****************************************************************************/
Edge::Edge(Vertex * const A, Vertex * const B) {

	a = A;
	b = B;

	neighbors[0] = NULL;
	neighbors[1] = NULL;

};
Edge::~Edge() {

	a = NULL;
	b = NULL;

	neighbors[0] = NULL;
	neighbors[1] = NULL;

};


void Edge::set(Vertex * const A, Vertex * const B) {

	a = A;
	b = B;

};
Triangle * Edge::neighbor(unsigned i) {

	return neighbors[i];

};
Triangle *& Edge::set_neighbor(unsigned i) {

	return neighbors[i];

};
bool Edge::contains(Vertex * const v) {

	return a == v || b == v;

};
bool Edge::contains(Vertex * const A, Vertex * const B) {

	return contains(A) && contains(B);

};
double Edge::length() {

	double const a_x = a->x;
	double const a_y = a->y;

	double const b_x = b->x;
	double const b_y = b->y;

	return sqrt(sqr(b_x - a_x) + sqr(b_y - a_y));

};
Vertex * Edge::get_shared_vertex(Edge * e) {

	if (a == e->a || a == e->b) return a;
	if (b == e->b || b == e->a) return b;
	
	assert(0 && "Edge::get_shared_vertex(Edge)");

	return NULL;

};


/*****************************************************************************/
/*																			 */
/*									TRIANGLE								 */
/*																			 */
/*****************************************************************************/
Triangle::Triangle(Vertex * const a, Vertex * const b, Vertex * const c) {


	neighbors[0] = NULL;
	neighbors[1] = NULL;
	neighbors[2] = NULL;

	edges[0] = NULL;
	edges[1] = NULL;
	edges[2] = NULL;

	if (orientation_robust(a, b, c) < 0.0) {

		vertices[0] = a;
		vertices[1] = c;
		vertices[2] = b;

	}
	else {

		vertices[0] = a;
		vertices[1] = b;
		vertices[2] = c;

	}

	if (orientation_robust(vertices[0], vertices[1], vertices[2]) <= 0.0) {

		// ORIENTATION ROBUST PROBABLY NOT WORKING ???????????

		//double oo = orientation_robust(vertices[0], vertices[1], vertices[2]);
		std::cout << "Big Phucking Problem" << std::endl;
		//oo = orientation_fast(vertices[0], vertices[1], vertices[2]);


	}
		

};
Triangle::~Triangle() {

	vertices[0] = vertices[1] = vertices[2] = NULL;
	edges[0] = edges[1] = edges[2] = NULL;
	neighbors[0] = neighbors[1] = neighbors[2] = NULL;

};


bool Triangle::contains(Vertex * const v) const {

	return vertices[0] == v || vertices[1] == v || vertices[2] == v;

};
bool Triangle::contains(Vertex * const a, Vertex * const b)	const {

	return contains(a) && contains(b);

};
bool Triangle::contains(Edge * const e) const {

	return contains(e->a, e->b);

};


Vertex * Triangle::get_vertex(unsigned i) const {

	return vertices[i];

};
Vertex * Triangle::get_vertex(Edge * const  e) const {

	if (e == edges[0]) return vertices[0];
	else if (e == edges[1]) return vertices[1];
	else if (e == edges[2]) return vertices[2];

	assert(0 && "Triangle::get_vertex(Edge)");

	return NULL;

};
Vertex * Triangle::get_vertex_cw(Vertex * const v) const {

	if		(vertices[0] == v) return vertices[2];
	else if (vertices[1] == v) return vertices[0];
	else if (vertices[2] == v) return vertices[1];
	
	assert(0 && "Triangle::vertex_cw(Vertex)");

	return NULL;

};
Vertex * Triangle::get_vertex_ccw(Vertex * const v) const {

	if		(vertices[0] == v) return vertices[1];
	else if (vertices[1] == v) return vertices[2];
	else if (vertices[2] == v) return vertices[0];

	assert(0 && "Triangle::vertex_ccw(Vertex)");

	return NULL;

};
Vertex * Triangle::get_vertex_but(Vertex * const a, Vertex * const b) const {

	if		((vertices[0] == a && vertices[1] == b) || (vertices[1] == a && vertices[0] == b)) return vertices[2];
	else if ((vertices[1] == a && vertices[2] == b) || (vertices[2] == a && vertices[1] == b)) return vertices[0];
	else if ((vertices[2] == a && vertices[0] == b) || (vertices[0] == a && vertices[2] == b)) return vertices[1];

	assert(0 && "Triangle::vertex_but(Vertex, Vertex)");

	return NULL;

};

Edge * Triangle::get_edge(unsigned i) const {

	return edges[i];

};
Edge * Triangle::get_edge(Vertex * const a, Vertex * const b) const {

	if		((a == vertices[0] && b == vertices[1]) || (a == vertices[1] && b == vertices[0])) return this->edges[2];
	else if ((a == vertices[1] && b == vertices[2]) || (a == vertices[2] && b == vertices[1])) return this->edges[0];
	else if ((a == vertices[2] && b == vertices[0]) || (a == vertices[0] && b == vertices[2])) return this->edges[1];

	assert(0 && "Triangle::edge(Vertex, Vertex)");

	return NULL;

};
Edge * Triangle::get_edge_cw(Vertex * const v) const {

	if		(vertices[0] == v) return edges[1];
	else if (vertices[1] == v) return edges[2];
	else if (vertices[2] == v) return edges[0];

	assert(0 && "Triangle::edge_cw(Vertex)");

	return NULL;

};
Edge * Triangle::get_edge_ccw(Vertex * const v) const {

	if		(vertices[0] == v) return edges[2];
	else if (vertices[1] == v) return edges[0];
	else if (vertices[2] == v) return edges[1];

	assert(0 && "Triangle::edge_ccw(Vertex)");

	return NULL;

};

void Triangle::set_edge(Edge * const e) {


	if		(e->contains(vertices[1], vertices[2])) edges[0] = e;
	else if (e->contains(vertices[0], vertices[2])) edges[1] = e;
	else if (e->contains(vertices[0], vertices[1])) edges[2] = e;

	else assert(0 && "Triangle::set_edge(Edge)");

	return;

};
Edge * & Triangle::set_edge(int i) {

	return edges[i];

};

unsigned Triangle::get_vertex_index(Vertex * const v) const {

	if		(v == vertices[0]) return 0;
	else if (v == vertices[1]) return 1;
	else if (v == vertices[2]) return 2;

	assert(0 && "Triangle::vertex_index(Vertex, Vertex)");

	return -1;

};
unsigned Triangle::get_edge_index(Edge * const e) const {

	if		(e == edges[0]) return 0;
	else if (e == edges[1]) return 1;
	else if (e == edges[2]) return 2;

	assert(0 && "Triangle::edge_index(Edge)");

	return -1;

};
unsigned Triangle::get_edge_index(Vertex * const a, Vertex * const b) const {

	if		((a == vertices[0] && b == vertices[1]) || (a == vertices[1] && b == vertices[0])) return 2;
	else if ((a == vertices[1] && b == vertices[2]) || (a == vertices[2] && b == vertices[1])) return 0;
	else if ((a == vertices[2] && b == vertices[0]) || (a == vertices[0] && b == vertices[2])) return 1;

	assert(0 && "Triangle::edge_index(Vertex, Vertex)");

	return -1;

};


Triangle * Triangle::get_neighbor(unsigned i) const {

	return neighbors[i];

};
Triangle * Triangle::get_neighbor(Vertex * const a, Vertex * const b) const {

	if		((a == vertices[2] && b == vertices[1]) || (a == vertices[1] && b == vertices[2]))	return neighbors[0];
	else if ((a == vertices[0] && b == vertices[2]) || (a == vertices[2] && b == vertices[0]))	return neighbors[1];
	else if ((a == vertices[0] && b == vertices[1]) || (a == vertices[1] && b == vertices[0]))	return neighbors[2];

	assert(0 && "Triangle::neighbor(Vertex, Vertex)");

	return NULL;

};
Triangle * Triangle::get_neighbor_cw(Vertex * const v) const {

	if		(vertices[0] == v) return neighbors[1];
	else if (vertices[1] == v) return neighbors[2];
	else if (vertices[2] == v) return neighbors[0];

	assert(0 && "Triangle::neighbor_cw(Vertex)");

	return NULL;

};
Triangle * Triangle::get_neighbor_ccw(Vertex * const v) const {

	if		(vertices[0] == v) return neighbors[2];
	else if (vertices[1] == v) return neighbors[0];
	else if (vertices[2] == v) return neighbors[1];

	assert(0 && "Triangle::neighbor_ccw(Vertex)");

	return NULL;

};
unsigned Triangle::get_neighbor_index(Triangle * const t) const {

	if		(neighbors[0] == t) return 0;
	else if (neighbors[1] == t) return 1;
	else if (neighbors[2] == t) return 2;
	
	assert(0 && "Triangle::neighbor_index(Triangle)");

	return -1;

};


void Triangle::set_neighbor(Triangle * const t, unsigned i) {

	neighbors[i] = t;

};
void Triangle::set_neighbor(Vertex * const a, Vertex * const b, Triangle * const t) {

	if		((a == vertices[2] && b == vertices[1]) || (a == vertices[1] && b == vertices[2]))	neighbors[0] = t;
	else if ((a == vertices[0] && b == vertices[2]) || (a == vertices[2] && b == vertices[0]))	neighbors[1] = t;
	else if ((a == vertices[0] && b == vertices[1]) || (a == vertices[1] && b == vertices[0]))	neighbors[2] = t;
	else assert(0 && "Triangle::set_neighbor(Vertex, Vertex, Triangle)");

	return;

};
void Triangle::set_neighbor(Edge * const e, Triangle * const t) {

	if		(e == edges[0]) neighbors[0] = t;
	else if (e == edges[1]) neighbors[1] = t;
	else if (e == edges[2]) neighbors[2] = t;
	else assert(0 && "Triangle::set_neighbor(Edge, Triangle)");

	return;

};
void Triangle::set_neighbor(Triangle * const t) {

	if (!t) 
		return;

	if (t->contains(vertices[1], vertices[2])) {

		neighbors[0] = t;
		t->set_neighbor(vertices[1], vertices[2], this);

	}
	else if (t->contains(vertices[0], vertices[2])) {

		neighbors[1] = t;
		t->set_neighbor(vertices[0], vertices[2], this);

	}
	else if (t->contains(vertices[0], vertices[1])) {

		neighbors[2] = t;
		t->set_neighbor(vertices[0], vertices[1], this);

	}
	else	
		assert(0 && "Triangle::set_neighbor(Triangle)");

	return;

};

void Triangle::null_neighbors() {

	neighbors[0] = NULL;
	neighbors[1] = NULL;
	neighbors[2] = NULL;

};



void Triangle::rotate_triangle_cw(Vertex * const v, Vertex * const ov) {

	if (v == vertices[0]) {

		vertices[1] = v; // vertices[0];
		vertices[0] = vertices[2];
		vertices[2] = ov;

	}
	else if (v == vertices[1]) {

		vertices[2] = v; // vertices[1];
		vertices[1] = vertices[0];
		vertices[0] = ov;

	}
	else if (v == vertices[2]) {

		vertices[0] = v; // vertices[2];
		vertices[2] = vertices[1];
		vertices[1] = ov;
	}
	else {

		assert(0 && "Triangle::rotate_triangle_cw(Vertex, Vertex)");

		return;

	}

	return;

};

Vertex * Triangle::get_opposite_vertex(Triangle * const t, Vertex * const v) const {

	Vertex * cw = t->get_vertex_cw(v);
	return get_vertex_cw(cw);

};
Triangle * Triangle::get_opposite_triangle(Vertex * const v) const {

	return get_neighbor(get_vertex_index(v));

};


double Triangle::area() const {

	return 0.5* abs(orientation_fast(vertices[0], vertices[1], vertices[2]));

};
void Triangle::circum_center(double & x_c, double & y_c) const {


	double const a_x = vertices[0]->x;
	double const a_y = vertices[0]->y;

	double const b_x = vertices[1]->x;
	double const b_y = vertices[1]->y;

	double const c_x = vertices[2]->x;
	double const c_y = vertices[2]->y;


	double const A = b_x - a_x;
	double const B = b_y - a_y;
	double const C = c_x - a_x;
	double const D = c_y - a_y;

	double const H = c_x - b_x;
	double const I = c_y - b_y;

	double const E = A * (a_x + b_x) + B * (a_y + b_y);
	double const F = C * (a_x + c_x) + D * (a_y + c_y);
	double const G = 2.0 * (A * I - B * H);

	// Circumcenter coordinates
	x_c = (D*E - B * F) / G;
	y_c = (A*F - C * E) / G;

};
double Triangle::circum_radius_squared() const {


	double const a_x = vertices[0]->x;
	double const a_y = vertices[0]->y;

	double const b_x = vertices[1]->x;
	double const b_y = vertices[1]->y;

	double const c_x = vertices[2]->x;
	double const c_y = vertices[2]->y;

	double const A = b_x - a_x;
	double const B = b_y - a_y;
	double const C = c_x - a_x;
	double const D = c_y - a_y;

	double const H = c_x - b_x;
	double const I = c_y - b_y;

	// Squared lenght of edges : |ab| , |bc| , |ca|
	double const AB = A * A + B * B;
	double const BC = H * H + I * I;
	double const CA = C * C + D * D;

	double const area_squared = sqr(C * I - D * H) / 4.0;

	// Formula for circum radius : sqrt( |ab| * |bc| * |ca| ) / ( 4 * Area )
	return AB * BC * CA / (16.0 * area_squared);

};
double Triangle::shortest_edge_squared() const {


	double const a_x = vertices[0]->x;
	double const a_y = vertices[0]->y;

	double const b_x = vertices[1]->x;
	double const b_y = vertices[1]->y;

	double const c_x = vertices[2]->x;
	double const c_y = vertices[2]->y;


	double const A = b_x - a_x;
	double const B = b_y - a_y;
	double const C = c_x - a_x;
	double const D = c_y - a_y;

	double const H = c_x - b_x;
	double const I = c_y - b_y;

	// Squared lenght of edges : |ab| , |bc| , |ca|
	double const AB = A * A + B * B;
	double const BC = H * H + I * I;
	double const CA = C * C + D * D;

	// Squared lenght of the shortest edge
	return std::fmin(AB, std::fmin(BC, CA));

};
double Triangle::ratio_squared() const {


	double const a_x = vertices[0]->x;
	double const a_y = vertices[0]->y;

	double const b_x = vertices[1]->x;
	double const b_y = vertices[1]->y;

	double const c_x = vertices[2]->x;
	double const c_y = vertices[2]->y;

	double const A = b_x - a_x;
	double const B = b_y - a_y;
	double const C = c_x - a_x;
	double const D = c_y - a_y;

	double const H = c_x - b_x;
	double const I = c_y - b_y;

	// Squared lenght of edges : |ab| , |bc| , |ca|
	double const AB = A * A + B * B;
	double const BC = H * H + I * I;
	double const CA = C * C + D * D;

	// Squared area
	double const area_squared = sqr(C * I - D * H) / 4.0;


	double const circum_radius_squared = AB * BC * CA / (16.0 * area_squared);
	double const shortest_edge_squared = std::fmin(AB, std::fmin(BC, CA));


	return circum_radius_squared / shortest_edge_squared;

};


inline double Triangle::orientation_determinant() const {

	double const x0 = vertices[0]->x;
	double const y0 = vertices[0]->y;

	double const x1 = vertices[1]->x;
	double const y1 = vertices[1]->y;

	double const x2 = vertices[2]->x;
	double const y2 = vertices[2]->y;

	return (x1 - x0)*(y2 - y0) - (x2 - x0)*(y1 - y0);

};


bool Triangle::is_bad(double const angle_bound, double const area_bound) {


	double const Bound = 0.5 / (sin(Pi * angle_bound / 180.0));
	
	if (this->ratio_squared() > sqr(Bound) || this->circum_radius_squared() > sqr(area_bound))
		return true;

	return false;

}