#include "mesh.h"




bool t_compare_x(Triangle * const t1, Triangle * const t2) {


	double const t1_min = std::fmin(std::fmin(t1->vertices[0]->x, t1->vertices[1]->x), t1->vertices[2]->x);
	double const t2_min = std::fmin(std::fmin(t2->vertices[0]->x, t2->vertices[1]->x), t2->vertices[2]->x);

	return t1_min < t2_min;

};
bool t_compare_y(Triangle * const t1, Triangle * const t2) {


	double const t1_min = std::fmin(std::fmin(t1->vertices[0]->y, t1->vertices[1]->y), t1->vertices[2]->y);
	double const t2_min = std::fmin(std::fmin(t2->vertices[0]->y, t2->vertices[1]->y), t2->vertices[2]->y);

	return t1_min < t2_min;

};
bool v_compare_x(Vertex * const v1, Vertex * const v2) {


	double const v1_min = v1->x;
	double const v2_min = v2->x;

	return v1_min < v2_min;

};
bool v_compare_y(Vertex * const v1, Vertex * const v2) {


	double const v1_min = v1->y;
	double const v2_min = v2->y;

	return v1_min < v2_min;

};
bool e_compare_x(Edge * const e1, Edge * const e2) {


	double const e1_min = std::fmin(e1->a->x, e1->b->x);
	double const e2_min = std::fmin(e2->a->x, e2->b->x);

	return e1_min < e2_min;

};
bool e_compare_y(Edge * const e1, Edge * const e2) {


	double const e1_min = std::fmin(e1->a->y, e1->b->y);
	double const e2_min = std::fmin(e2->a->y, e2->b->y);

	return e1_min < e2_min;

};


bool marker_compare_neumann(Edge * const e1, Edge * const e2) {


	E_MARKER const e1_m = e1->marker;
	E_MARKER const e2_m = e2->marker;

	if (e1_m == E_MARKER::NEUMANN && e2_m != E_MARKER::NEUMANN)
		return true;

	return false;

};
bool marker_compare_dirichlet(Edge * const e1, Edge * const e2) {


	E_MARKER const e1_m = e1->marker;
	E_MARKER const e2_m = e2->marker;

	if (e1_m == E_MARKER::DIRICHLET && e2_m != E_MARKER::DIRICHLET)
		return true;

	return false;

};





Mesh::~Mesh() {

	clear();

};



Vertex * Mesh::get_vertex(unsigned i) const {

	return vertices[i];

}
Edge * Mesh::get_edge(unsigned i) const {

	return edges[i];

}
Triangle * Mesh::get_triangle(unsigned i) const {

	return triangles[i];

}



unsigned Mesh::get_number_of_vertices() const {

	return num_vertices;

};
unsigned Mesh::get_number_of_edges() const {

	return num_edges;

};
unsigned Mesh::get_number_of_triangles() const {

	return num_triangles;

};


unsigned Mesh::get_number_of_dirichlet_edges() const {

	return num_dirichlet;

};
unsigned Mesh::get_number_of_neumann_edges() const {

	return num_neumann;

};


void Mesh::export_vertices(std::ofstream & stream) const {

	size_t const n = vertices.size();

	for (size_t i = 0; i < n; i++) {

		double const v_x = vertices[i]->x;
		double const v_y = vertices[i]->y;

		stream << v_x << "  " << v_y << std::endl;

	}

};
void Mesh::export_edges(std::ofstream & stream) const {

	size_t const n = edges.size();


	//// Gnuplot
	for (size_t i = 0; i < n; i++) {

		//if (edges[i]->is_constrained) {

		double const v0_x = edges[i]->a->x;
		double const v0_y = edges[i]->a->y;

		double const v1_x = edges[i]->b->x;
		double const v1_y = edges[i]->b->y;

		stream << v0_x << "  " << v0_y << std::endl;
		stream << v1_x << "  " << v1_y << std::endl << std::endl;

		//}

	}

	//// Matlab
	//for (size_t i = 0; i < n; i++) {

	//	//if (edges[i]->is_constrained) {

	//		double const v0_x = edges[i]->a->x;
	//		double const v0_y = edges[i]->a->y;

	//		double const v1_x = edges[i]->b->x;
	//		double const v1_y = edges[i]->b->y;

	//		stream << v0_x << "  " << v0_y << "  " << v1_x << "  " << v1_y << std::endl;

	//	//}

	//}

};
void Mesh::export_triangles(std::ofstream & stream) const {

	size_t const n = triangles.size();

	for (size_t i = 0; i < n; i++) {

		Triangle * t = triangles[i];

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


void Mesh::clear() {


	Vertex	 * v = NULL;
	Edge	 * e = NULL;
	Triangle * t = NULL;


	// Free allocated memory
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

	// Clear the vectors of the rubbish
	vertices.clear();
	edges.clear();
	triangles.clear();

};




template<NUMBERING algorithm>
void Mesh::apply_numbering() {


	switch (algorithm) {


	case NUMBERING::CM:

		cuthill_mckee();
		break;

	}


};

void Mesh::cuthill_mckee() {






};

