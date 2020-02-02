#pragma once



#include "shapes.h"
#include "enumerators.h"
//#include "triangulation.h"

#include <iostream>
#include <vector>

#include <list>


typedef Vertex *	v_pointer;
typedef Edge *		e_pointer;
typedef Triangle *	t_pointer;


class PlanarStraightLineGraph {


	template<GEOMETRIC_KERNEL GK>
	friend class Triangulation;


private:

	std::vector<v_pointer>		vertices_pslg;
	std::vector<unsigned *>		segments_pslg;		// 0. index of first vertex 1. index of second vertex 2. edge_marker 3. constraint_method

public:

	PlanarStraightLineGraph();
	~PlanarStraightLineGraph();

	template<V_MARKER marker>
	v_pointer insert_vertex(Vertex v);

	v_pointer get_vertex(unsigned i);

	// Inputed vertices are automatically marked as Constrained
	template<E_MARKER edge_marker, INSERT_CONSTRAINT constraint_method = INSERT_CONSTRAINT::INTERSECTION>
	void insert_constraint(v_pointer const a, v_pointer const b);

	// Inputed vertices are automatically marked as Constrained
	template<E_MARKER edge_marker, INSERT_CONSTRAINT constraint_method = INSERT_CONSTRAINT::INTERSECTION>
	void insert_constraint(unsigned i_a, unsigned i_b);


	unsigned get_number_of_vertices() const;
	unsigned get_number_of_segments() const;

	void clear();
	
};


template<V_MARKER marker>
v_pointer PlanarStraightLineGraph::insert_vertex(Vertex v) {


	size_t const nv = vertices_pslg.size();

	for (size_t i = 0; i < nv; i++) {

		if (v.is_almost_equal(vertices_pslg[i])) {

			//std::cout << "Input vertex is duplicate to some other. Returning its pointer." << std::endl;
			return vertices_pslg[i];

		}
	}

	v_pointer const new_vertex = new Vertex(v.x, v.y);

	new_vertex->index = vertices_pslg.size();
	new_vertex->marker = marker;

	vertices_pslg.push_back(new_vertex);

	return new_vertex;

};

template<E_MARKER edge_marker, INSERT_CONSTRAINT constraint_method>
void PlanarStraightLineGraph::insert_constraint(v_pointer const a, v_pointer const b) {


	unsigned const i_a = a->index;
	unsigned const i_b = b->index;

	if (i_a == i_b)
		std::cout << "Input vertices are duplicates." << std::endl;

	unsigned * constraint = new unsigned[4];

	constraint[0] = i_a;
	constraint[1] = i_b;
	constraint[2] = (unsigned)edge_marker;
	constraint[3] = (unsigned)constraint_method;

	segments_pslg.push_back(constraint);

};
template<E_MARKER edge_marker, INSERT_CONSTRAINT constraint_method>
void PlanarStraightLineGraph::insert_constraint(unsigned i_a, unsigned i_b) {


	size_t const nv = vertices_pslg.size();

	if (i_a == i_b || i_a >= nv || i_b >= nv)
		std::cout << "Invalid input." << std::endl;


	unsigned * constraint = new unsigned[4];

	constraint[0] = i_a;
	constraint[1] = i_b;
	constraint[2] = (unsigned)edge_marker;
	constraint[3] = (unsigned)constraint_method;

	segments_pslg.push_back(constraint);

};