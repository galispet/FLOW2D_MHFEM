
#include "PSLG.h"


PlanarStraightLineGraph::PlanarStraightLineGraph() {

};
PlanarStraightLineGraph::~PlanarStraightLineGraph() {

	clear();

};


v_pointer PlanarStraightLineGraph::get_vertex(unsigned i) {

	return vertices_pslg[i];

}




unsigned PlanarStraightLineGraph::get_number_of_vertices() const {

	return vertices_pslg.size();

};
unsigned PlanarStraightLineGraph::get_number_of_segments() const {

	return segments_pslg.size();

};

void PlanarStraightLineGraph::clear() {


	size_t const nv = vertices_pslg.size();
	size_t const ns = segments_pslg.size();

	for (size_t i = 0; i < nv; i++)
		delete vertices_pslg[i];

	for (size_t i = 0; i < ns; i++)
		delete segments_pslg[i];

	vertices_pslg.clear();
	segments_pslg.clear();

};