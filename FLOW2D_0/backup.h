#pragma once





/*
template <GEOMETRIC_KERNEL GK>
void Triangulation<GK>::mark_triangle(t_pointer const t) {


	t->marker = T_MARKER::OUTSIDE;


	for (unsigned j = 0; j < 3; j++) {


		e_pointer const e = t->edges[j];

		t_pointer const neighbor = t->neighbors[j];


		if (e) {


			if (!e->is_constrained) {

				remove(edges_tri, e);
				delete e;

			}
			else {

				if (neighbor) {

					unsigned const e_index = neighbor->get_edge_index(e);

					neighbor->edges[e_index] = NULL;
					neighbor->neighbors[e_index] = NULL;

				}

				continue;

			}



			if (!neighbor)
				continue;

			unsigned const e_index = neighbor->get_edge_index(e);

			neighbor->edges[e_index] = NULL;
			neighbor->neighbors[e_index] = NULL;

		}

		if (!neighbor || neighbor->marker == T_MARKER::OUTSIDE)
			continue;

		mark_triangle(neighbor);

	}

	remove(triangles_tri, t);
	delete t;

};
*/


/*

template<GEOMETRIC_KERNEL GK>
bool Triangulation<GK>::split_encroached_segments(v_pointer const p) {


	bool has_occured_splitting = false;


	for (size_t i = 0; i < segments_tri.size(); i++) {


		e_pointer const e = segments_tri[i];


		if (e->contains(p))
			continue;


		double x_mid = INFINITY;
		double y_mid = INFINITY;


		if (is_encroeched(e, p, x_mid, y_mid)) {


			v_pointer const new_vertex = new Vertex(x_mid, y_mid);

			new_vertex->marker = V_MARKER::CONSTRAINED;


			remove(segments_tri, e);

			e_pointer subsegment1 = NULL;
			e_pointer subsegment2 = NULL;

			insert_vertex_into_edge(e, new_vertex, subsegment1, subsegment2);


			segments_tri.push_back(subsegment1);
			segments_tri.push_back(subsegment2);


			has_occured_splitting = true;

			// It is sufficient to start over with the next edge after this deleted edge 'e' -> because it is deleted, the iterator 'i' must be  lowered by 1
			//i--;

			// Inserted vertex into subsegment may some other lead to be encroached ...
			i = -1;

			break;

		}
	}

	return has_occured_splitting;

};

template<GEOMETRIC_KERNEL GK>
void Triangulation<GK>::split_encroached_segments(std::vector<e_pointer> &encroached_segments, v_pointer const p) {


	size_t const ns = encroached_segments.size();

	std::vector<e_pointer> new_subsegments;
	new_subsegments.reserve(2 * ns);

	for (size_t i = 0; i < encroached_segments.size(); i++) {


		e_pointer const e = encroached_segments[i];

		double x_mid;
		double y_mid;

		get_mid_point(e, x_mid, y_mid);


		v_pointer const new_vertex = new Vertex(x_mid, y_mid);

		new_vertex->marker = V_MARKER::CONSTRAINED;
		vertices_tri.push_back(new_vertex);


		//remove(encroached_segments, e);
		remove(segments_tri, e);


		e_pointer subsegment1 = NULL;
		e_pointer subsegment2 = NULL;

		insert_vertex_into_edge(e, new_vertex, subsegment1, subsegment2);

		new_subsegments.push_back(subsegment1);
		new_subsegments.push_back(subsegment2);

		//i--;

		// Recursively split sugmenets if needed. Thats why there is vertex 'p' in the argumet list !!!!!!1

	}


	for (size_t i = 0; i < new_subsegments.size(); i++)
		segments_tri.push_back(new_subsegments[i]);


};
*/


/*

double find_circumcenter(Vertex * const a, Vertex * const b, Vertex * const c, double &x_c, double &y_c) {


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

	double const E = A * (a_x + b_x) + B * (a_y + b_y);
	double const F = C * (a_x + c_x) + D * (a_y + c_y);
	double const G = 2.0 * (A * I - B * H);


	// Circumcenter coordinates
	x_c = (D*E - B * F) / G;
	y_c = (A*F - C * E) / G;

	double const dx = x_c - a_x;
	double const dy = y_c - a_y;

	// Circumcenter radius squared
	return dx * dx + dy * dy;

};

double lenght_shortest_edge(Triangle * const t) {


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

	// Squared lenght of edges : |ab| , |bc| , |ca|
	double const AB = A * A + B * B;
	double const BC = H * H + I * I;
	double const CA = C * C + D * D;

	// Squared lenght of the shortest edge
	return std::fmin(AB, std::fmin(BC, CA));

};

*/


/*
template<GEOMETRIC_KERNEL GK>
t_pointer Triangulation<GK>::get_skinny_triangle_ruppert() {


	double const Area_Max = Area_Bound;
	double const Bound = 0.5 / (sin(Pi * Angle_Bound / 180.0));


	std::sort(triangles_tri.begin(), triangles_tri.end(), [](t_pointer const t, t_pointer const k) { return t->ratio_squared() < k->ratio_squared(); });



	for (size_t i = 0; i < triangles_tri.size(); i++) {


		t_pointer const t = triangles_tri[i];

		double const circum_radius_squared = t->circum_radius_squared();
		double const ratio_squared = t->ratio_squared();


		if (ratio_squared > sqr(Bound) || circum_radius_squared > sqr(Area_Max))
			return t;

	}

	return NULL;

};
void refinement_ruppert(double const angle_bound, double const area_max = INFINITY);
template<GEOMETRIC_KERNEL GK>
void Triangulation<GK>::refinement_ruppert(double const angle_bound, double const area_max) {


	Angle_Bound = angle_bound;
	Area_Bound = area_max;

	bool continue_refining = true;
	bool split_segments = true;

	while (true) {


		continue_refining = false;


		// Split every encroached segments / subsegments there are
		if (split_segments) {


			for (size_t i = 0; i < vertices_tri.size(); i++) {

				// !!!! Check if the vertex is visible. If not, do not add the segment
				v_pointer const v = vertices_tri[i];

				std::vector<e_pointer> encroached_segments;

				// Maybe we should chek encroached subsegments one by one, so there won't be visibility problem
				if (get_encroached_segments(v, encroached_segments)) {

					continue_refining = true;

					for (size_t j = 0; j < encroached_segments.size(); j++)
						split_encroached_segment(encroached_segments[j], v);

				}
			}
		}



		// Get any skinny triangle and its circumcenter
		t_pointer const skinny_triangle = get_skinny_triangle_ruppert();


		// If there are no skinny triangles left, we are done
		if (!skinny_triangle && !continue_refining)
			break;


		// Circum center coordinates
		double x_c;
		double y_c;

		skinny_triangle->circum_center(x_c, y_c);


		// Circumcenter of skinny triangle
		Vertex circum_center(x_c, y_c);

		std::vector<e_pointer> encroached_segments2;


		// !!!! Check if the circumcenter is visible. If not, do not add the segment
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
			insert_vertex(new Vertex(x_c, y_c), skinny_triangle);

		}

	}

};
*/


/*

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


		if (neighbor->marker == T_MARKER::OUTSIDE) {

			// If the neighbor is already processed
			continue;

		}
		else if (t->edges[i]->is_constrained) {

			// If the neighbor is part of triangulation which should not be deleted (domain boundaded by constrained segments)
			//e_pointer const e = t->edges[i];

			//v_pointer const a = e->a;
			//v_pointer const b = e->b;

			// Set the neighbor of the triangle on the boundary to NULL (it will be deleted)
			//neighbor->neighbors[neighbor->get_neighbor_index(t)] = NULL;


			// Set the adjacency of the constrained edge
			//e->neighbors[0] = neighbor;
			//e->neighbors[1] = NULL;

			// Update new neighborhood adjacency of vertices (otherwise it may point to marked triangle to be deleted)
			//a->adjacent_triangle = neighbor;
			//b->adjacent_triangle = neighbor;

			continue;

		}

		// Recurrent depth-search
		mark_triangle(neighbor, neighbor->get_neighbor_index(t));

	}

	// Insert this triangle into stack of triangles to be deleted (marked as outside)
	triangles_outside.push_back(t);

};

*/