#pragma once





/*
void solver::assemble_χ() {


	// Quadrature weights and points on reference segment [-1,1]
	gauss_quadrature_1D const quad(quadrature_order);
	unsigned const num_quad_points = quad.number_of_points;


	//Matrix<real> normals(2, 3);
	//Vector<real> parametrization(2);
	//Vector<real> parametrizationDerivative(2);
	//Matrix<real> basisRaviartThomas(2, 8);
	//Vector<real> basisEdgePolynomial(2);

	Eigen::MatrixXd normals(2, 3);
	Eigen::VectorXd parametrization(2);
	Eigen::VectorXd parametrizationDerivative(2);
	Eigen::MatrixXd basisRaviartThomas(2, 8);
	Eigen::VectorXd basisEdgePolynomial(2);



	χ.setZero();

	for (unsigned k = 0; k < nk; k++) {


		t_pointer const K = mesh->get_triangle(k);
		unsigned const k_index = K->index;


		//--------------------------------

		//Matrix<real> JF(2, 2);
		//
		//v_pointer const a = K->vertices[0];
		//v_pointer const b = K->vertices[1];
		//v_pointer const c = K->vertices[2];
		//
		//real const x0 = (real)a->x;
		//real const y0 = (real)a->y;
		//
		//real const x1 = (real)b->x;
		//real const y1 = (real)b->y;
		//
		//real const x2 = (real)c->x;
		//real const y2 = (real)c->y;
		//
		//real const detJF = (real) ((x1 - x0)*(y2 - y0) - (x2 - x0)*(y1 - y0));
		//
		//JF(0, 0) = x1 - x0;
		//JF(0, 1) = x2 - x0;
		//JF(1, 0) = y1 - y0;
		//JF(1, 1) = y2 - y0;

		//--------------------------------


		real orientations[3] = { 1.0 , 1.0 , 1.0 };
		real orientations3[3] = { 1.0 , 1.0 , 1.0 };
		//real orientations2[3] = { 1.0 , 1.0 , 1.0 };

		for (unsigned e = 0; e < 3; e++)
			orientations[e] = edgeOrientation(k_index, 0, e, 0);
		for (unsigned e = 0; e < 3; e++)
			orientations3[e] = edgeOrientation2(k_index, 0, e, 0);


		evaluate_edge_normal(normals, orientations3);
		//evaluate_edge_normal(normals, orientations3);


		for (unsigned n = 0; n < num_quad_points; n++) {

			for (unsigned El = 0; El < 3; El++) {


				real const orientation = orientations[El];
				real const orientation2 = orientations3[El];


				real const a = (real) 0.0;
				real const b = (real) El != 0 ? 1.0 : sqrt(2.0);

				//real a;
				//real b;
				//if (orientation == 1.0) {
				//	a = (real) 0.0;
				//	b = (real) El != 0 ? 1.0 : sqrt(2.0);
				//}
				//else {
				//	b = (real) 0.0;
				//	a = (real) El != 0 ? 1.0 : sqrt(2.0);
				//}


				real const c = (real)(b - a) / 2.0;
				real const d = (real)(b + a) / 2.0;

				real const x = (real)quad.points[n] * c + d;
				real const w = (real)quad.weigths[n] * c;

				//for (unsigned i = 0; i < num_quad_points; i++) {
				//
				//	real y = (real)quad.points[i] * c + d;
				//	evaluate_edge_parametrization(y, El, parametrization, orientation);
				//
				//	real const s = parametrization(0);
				//	real const t = parametrization(1);
				//
				//	std::cout << s << "    " << t << std::endl;
				//
				//}



				//real orientation2 = 1.0;
				evaluate_edge_parametrization(x, El, parametrization, orientation);
				evaluate_edge_parametrization_derivative(x, El, parametrizationDerivative, orientation);

				//Vector<real> const normal = normals.getColumn(El);
				//Vector<real> const r = parametrization;
				//Vector<real> const dr = parametrizationDerivative;

				Eigen::VectorXd const normal = normals.col(El);
				Eigen::VectorXd const r = parametrization;
				Eigen::VectorXd const dr = parametrizationDerivative;


				real const s = r(0);
				real const t = r(1);
				//real const drNorm = dr.norm();
				real const drNorm = 1.0;

				// Values of basis function on reference triangle' edge
				evaluate_raviartthomas_basis(s, t, basisRaviartThomas, orientations);
				evaluate_edge_polynomial_basis(x, El, basisEdgePolynomial, orientation);

				//Vector<real> const edgeBasis = basisEdgePolynomial;
				Eigen::VectorXd const edgeBasis = basisEdgePolynomial;


				for (unsigned m = 0; m < 8; m++) {


					//Vector<real> const Wm = basisRaviartThomas.getColumn(m);
					Eigen::VectorXd const Wm = basisRaviartThomas.col(m);

					//real const dotProduct = dot(Wm, normal);
					real const dotProduct = Wm.dot(normal);



					for (unsigned s = 0; s < 2; s++) {

						//real const varPhis = s == 1 ? orientation * edgeBasis(s) : edgeBasis(s);
						//real const varPhis = s == 1 ? orientation2 * edgeBasis(s) : edgeBasis(s);
						//real const varPhis = orientation * edgeBasis(s);
						real const varPhis = edgeBasis(s);

						χ.setCoeff(k_index, m, El, s) = χ(k_index, m, El, s) + w * dotProduct * varPhis * drNorm;

					}
				}
			}
		}

		for (unsigned El = 0; El < 3; El++)
			for (unsigned m = 0; m < 8; m++)
				for (unsigned j = 0; j < 2; j++)
					χ.setCoeff(k_index, m, El, j) = abs(χ(k_index, m, El, j)) < INTEGRAL_PRECISION ? 0.0 : χ(k_index, m, El, j);


		//for (unsigned m = 0; m < 8; m++) {
		////
		//	Matrix<real> integral(3, 2);
		//	integral.setZero();
		////
		//	for (unsigned El = 0; El < 3; El++)
		//		for (unsigned j = 0; j < 2; j++)
		//			integral(El, j) = χ(k_index, m, El, j);
		////
		//	std::cout << integral << std::endl;
		////
		//}

	}

};
*/

/*
void solver::assembleM() {


	M_j1_s1.setZero();
	M_j1_s2.setZero();

	M_j2_s1.setZero();
	M_j2_s2.setZero();


	for (unsigned e = 0; e < ne; e++) {


		e_pointer const E = mesh->get_edge(e);
		unsigned const e_index = E->index;
		E_MARKER const e_marker = E->marker;



		if (e_marker == E_MARKER::DIRICHLET) {

			M_j1_s1.coeffRef(e_index, e_index) = -1.0;
			M_j1_s2.coeffRef(e_index, e_index) =  0.0;

			M_j2_s1.coeffRef(e_index, e_index) =  0.0;
			M_j2_s2.coeffRef(e_index, e_index) = -1.0;

			continue;

		}

		//unsigned dof0;
		//unsigned dof1;
		//if (E->neighbors[0]) {
		//	dof0 = LI(E->neighbors[0], E, 0);
		//	dof1 = LI(E->neighbors[0], E, 1);
		//}
		//else {
		//	dof0 = LI(E->neighbors[1], E, 0);
		//	dof1 = LI(E->neighbors[1], E, 1);
		//}


		for (unsigned neighborElement = 0; neighborElement < 2; neighborElement++) {


			t_pointer const K = E->neighbors[neighborElement];

			if (!K)
				continue;


			unsigned const k_index = K->index;

			unsigned const dof0 = LI(K, E, 0);
			unsigned const dof1 = LI(K, E, 1);


			real const coeffChi0 = χ(k_index, LI(K, E, 0), K->get_edge_index(E), 0);
			real const coeffChi1 = χ(k_index, LI(K, E, 1), K->get_edge_index(E), 1);

			real const orient = edgeOrientation(K->index, 0, e_index, 0);

			// Loop over edges
			for (unsigned El = 0; El < 3; El++) {


				e_pointer const E_local = K->edges[El];

				unsigned const e_local_index_local = K->get_edge_index(E_local);	// Local index of local edge
				unsigned const e_local_index_global = E_local->index;				// Global index of local edge


				if (E_local == E) {

					for (unsigned i = 0; i < 8; i++) {

						M_j1_s1.coeffRef(e_index, e_index) += coeffChi0 * α(k_index, 0, i, dof0) * χ(k_index, i, e_local_index_local, 0) / viscosities[k_index];
						M_j1_s2.coeffRef(e_index, e_index) += coeffChi0 * α(k_index, 0, i, dof0) * χ(k_index, i, e_local_index_local, 1) / viscosities[k_index];

						M_j2_s1.coeffRef(e_index, e_index) += coeffChi1 * α(k_index, 0, i, dof1) * χ(k_index, i, e_local_index_local, 0) / viscosities[k_index];
						M_j2_s2.coeffRef(e_index, e_index) += coeffChi1 * α(k_index, 0, i, dof1) * χ(k_index, i, e_local_index_local, 1) / viscosities[k_index];

					}

					continue;

				}

				real ACHI_j1_s1 = 0.0;
				real ACHI_j1_s2 = 0.0;
				real ACHI_j2_s1 = 0.0;
				real ACHI_j2_s2 = 0.0;

				for (unsigned i = 0; i < 8; i++) {

					ACHI_j1_s1 += coeffChi0 * α(k_index, 0, i, dof0) * χ(k_index, i, e_local_index_local, 0);
					ACHI_j1_s2 += coeffChi0 * α(k_index, 0, i, dof0) * χ(k_index, i, e_local_index_local, 1);

					ACHI_j2_s1 += coeffChi1 * α(k_index, 0, i, dof1) * χ(k_index, i, e_local_index_local, 0);
					ACHI_j2_s2 += coeffChi1 * α(k_index, 0, i, dof1) * χ(k_index, i, e_local_index_local, 1);

				}


				//real const val1 = ACHI_j1_s1 / viscosities[k_index];
				//real const val2 = ACHI_j1_s2 / viscosities[k_index];
				//real const val3 = ACHI_j2_s1 / viscosities[k_index];
				//real const val4 = ACHI_j2_s2 / viscosities[k_index];
				//
				//M_j1_s1.coeffRef(e_index, e_local_index_global) = abs(val1) < INTEGRAL_PRECISION ? 0.0 : val1;
				//M_j1_s2.coeffRef(e_index, e_local_index_global) = abs(val2) < INTEGRAL_PRECISION ? 0.0 : val2;
				//
				//M_j2_s1.coeffRef(e_index, e_local_index_global) = abs(val3) < INTEGRAL_PRECISION ? 0.0 : val3;
				//M_j2_s2.coeffRef(e_index, e_local_index_global) = abs(val4) < INTEGRAL_PRECISION ? 0.0 : val4;


				M_j1_s1.coeffRef(e_index, e_local_index_global) = ACHI_j1_s1 / viscosities[k_index];
				M_j1_s2.coeffRef(e_index, e_local_index_global) = ACHI_j1_s2 / viscosities[k_index];

				M_j2_s1.coeffRef(e_index, e_local_index_global) = ACHI_j2_s1 / viscosities[k_index];
				M_j2_s2.coeffRef(e_index, e_local_index_global) = ACHI_j2_s2 / viscosities[k_index];

			}
		}
	}



	//for (unsigned e = 0; e < ne; e++) {
	//
	//	unsigned const e_index = mesh->get_edge(e)->index;
	//
	//	real const val1 = M_j1_s1.coeff(e_index, e_index);
	//	real const val2 = M_j1_s2.coeff(e_index, e_index);
	//	real const val3 = M_j2_s1.coeff(e_index, e_index);
	//	real const val4 = M_j2_s2.coeff(e_index, e_index);
	//
	//	M_j1_s1.coeffRef(e_index, e_index) = abs(val1) < INTEGRAL_PRECISION ? 0.0 : val1;
	//	M_j1_s2.coeffRef(e_index, e_index) = abs(val2) < INTEGRAL_PRECISION ? 0.0 : val2;
	//
	//	M_j2_s1.coeffRef(e_index, e_index) = abs(val3) < INTEGRAL_PRECISION ? 0.0 : val3;
	//	M_j2_s2.coeffRef(e_index, e_index) = abs(val4) < INTEGRAL_PRECISION ? 0.0 : val4;
	//
	//}
	//
	//
	//M_j1_s1 *= -1.0;
	//M_j1_s2 *= -1.0;
	//M_j2_s1 *= -1.0;
	//M_j2_s2 *= -1.0;
	//
	// j = 1
	//std::cout << M_j1_s1.toDense() << std::endl;
	//std::cout << M_j1_s2.toDense() << std::endl;
	//
	//// j = 2
	//std::cout << M_j2_s1.toDense() << std::endl;
	//std::cout << M_j2_s2.toDense() << std::endl;



};
*/

/*
void solver::upwindConcentration(t_pointer const K, double const t, CoeffMatrix<1, 3, 1> & ksi, double(&out)[3]) {


	Matrix <real> normal(2, 3);
	Matrix<real> parametrization(2, 3);
	Matrix<real> parametrizationDerivative(2, 8);
	Matrix<real> basisRaviartThomas(2, 8);
	Vector<real> basisPolynomial(3);
	Matrix<real> JF(2, 2);

	unsigned const k_index = K->index;

	v_pointer const a = K->vertices[0];
	v_pointer const b = K->vertices[1];
	v_pointer const c = K->vertices[2];

	real const x0 = (real)a->x;
	real const y0 = (real)a->y;

	real const x1 = (real)b->x;
	real const y1 = (real)b->y;

	real const x2 = (real)c->x;
	real const y2 = (real)c->y;

	JF(0, 0) = x1 - x0;
	JF(0, 1) = x2 - x0;
	JF(1, 0) = y1 - y0;
	JF(1, 1) = y2 - y0;

	evaluate_edge_normal(normal);


	for (unsigned Ei = 0; Ei < 3; Ei++) {


		e_pointer const E = K->get_edge(Ei);

		//unsigned const e_index_local = K->get_edge_index(E);
		E_MARKER const edge_marker = E->marker;

		t_pointer const Kn = K->neighbors[Ei];
		unsigned const kn_index = Kn->index;

		real const orientation[3] = { 1.0, 1.0, 1.0 };

		Vector<real> const v = normal.getColumn(Ei);




		real integral = 0.0;

		gauss_quadrature_1D quad(quadrature_order);
		unsigned const num_quad_points = quad.number_of_points;

		for (unsigned n = 0; n < num_quad_points; n++) {


			real const a = (real) 0.0;
			real const b = (real) Ei != 0 ? 1.0 : sqrt(2.0);

			real const c = (real) (b - a) / 2.0;
			real const d = (real) (b + a) / 2.0;

			real const x = (real) quad.points[n] * c + d;
			real const w = (real) quad.weigths[n] * c;

			evaluate_edge_parametrization(x, parametrization);
			evaluate_edge_parametrization_derivative(x, parametrizationDerivative);

			Vector<real> const r = parametrization.getColumn(Ei);
			Vector<real> const dr = parametrizationDerivative.getColumn(Ei);

			real const s = r(0);
			real const t = r(1);
			real const drNorm = dr.norm();

			// Values of basis function on reference triangle' edge at the gauss point (s,t)
			evaluate_polynomial_basis(s, t, basisPolynomial);
			evaluate_raviartthomas_basis(s, t, basisRaviartThomas, orientation);


			real velocitySign = 0.0;
			real concentration = 0.0;

			for (unsigned j = 0; j < 8; j++)
				velocitySign += velocities(K->index, 0, j, 0) * dot(v, JF *basisRaviartThomas.getColumn(j));

			// Upwinded concentration at the gauss point (s,t) of the edge E of the element K
			if (velocitySign >= 0.0)
				for (unsigned l = 0; l < 3; l++)
					concentration += ksi(k_index, 0, l, 0)*basisPolynomial(l);
			// Upwinded concentration at the gauss point (s,t) of the edge E of the element Kn
			else
				for (unsigned l = 0; l < 3; l++)
					concentration += ksi(kn_index, 0, l, 0)*basisPolynomial(l);


			integral += w * concentration * basisRaviartThomas.getColumn() * basisPolynomial(m);



			if (velocitySign >= 0.0)
				for (unsigned m = 0; m < 3; m++)
					integral += weight * ksi(K->index, 0, m, 0)*basisPolynomial(m);
			else
				for (unsigned m = 0; m < 3; m++)
					integral += weight * ksi(Kn->index, 0, m, 0)*basisPolynomial(m);



		}








		real integral = 0.0;



		out[Ei] = integral / E->length();

	}


	//unsigned const edge_index = K->get_edge_index(e);
	//E_MARKER const edge_marker = e->marker;
	//
	//t_pointer const Kn = K->neighbors[edge_index];
	//
	//v_pointer a = K->vertices[0];
	//v_pointer b = K->vertices[1];
	//v_pointer c = K->vertices[2];
	//
	//double const x0 = a->x;
	//double const y0 = a->y;
	//
	//double const x1 = b->x;
	//double const y1 = b->y;
	//
	//double const x2 = c->x;
	//double const y2 = c->y;
	//
	//double const iarea2x = 0.5 / K->area();
	//
	//// Centroid on element 'K'
	//double const centroid_x = (x0 + x1 + x2) / 3.0;
	//double const centroid_y = (y0 + y1 + y2) / 3.0;
	//
	//
	//// Transformed centroid onto the reference triangle
	//double const ref_centroid_x = iarea2x * (+(y2 - y0)*(centroid_x - x0) - (x2 - x0)*(centroid_y - y0));
	//double const ref_centroid_y = iarea2x * (-(y1 - y0)*(centroid_x - x0) + (x1 - x0)*(centroid_y - y0));
	//
	//
	//// 1. Physical velocity on edge 'e' with index in (0,1,2)  2.Global index of element K  3. row  4 column
	//double const velocity = velocities(K->index, edge_index, 0, 0);
	//
	//
	//double const CK = ksi(K->index, 0, 0, 0)*phi1(ref_centroid_x, ref_centroid_y);
	//
	//double aKn = 0.0;
	//double bKn = 0.0;
	//
	//double Knmid = 0.0;
	//double CKn = 0.0;
	//
	//double slope = 0.0;
	//double C = 0.0;
	//
	//
	//if (edge_marker == E_MARKER::NEUMANN) {
	//
	//	CKn = DIRICHLET_GAMMA_Q_concentration(x, t);
	//
	//	if (velocity < 0.0)
	//		return CKn;
	//
	//}
	//else if (edge_marker == E_MARKER::DIRICHLET) {
	//
	//	CKn = DIRICHLET_GAMMA_P_concentration(x, t);
	//
	//	if (velocity < 0.0)
	//		return CKn;
	//
	//}
	//else {
	//
	//	aKn = Kn->nodes[0]->x;
	//	bKn = Kn->nodes[1]->x;
	//
	//	Knmid = 0.5*(aKn + bKn);
	//
	//	CKn = ksi(Kn, 0, 0)*phi1(Knmid, aKn, bKn);
	//
	//}
	//
	//
	//double const Cmin = std::min(CK, CKn);
	//double const Cmax = std::max(CK, CKn);
	//
	//
	//if (velocity >= 0.0) {
	//
	//	slope = ksi(K, 0, 1)*phi2(x, aK, bK);
	//	C = CK + slope;
	//
	//}
	//else {
	//
	//	slope = ksi(Kn, 0, 1)*phi2(x, aKn, bKn);
	//	C = CKn + slope;
	//
	//}
	//
	////return C;
	//
	//if (Cmin <= C && C <= Cmax) {
	//
	//	//std::cout << "no limiter" << std::endl;
	//	return C;
	//
	//}
	////std::cout << "limiter" << std::endl;
	//
	//if (C < Cmin)
	//	return Cmin;
	//
	//return Cmax;
	//


};
*/

/*
void solver::assemble_χ() {


	unsigned const nk = _nk;


	//for (unsigned k = 0; k < nk; k++) {
	//
	//
	//	t_pointer const K = mesh->get_triangle(k);
	//	unsigned const k_index = K->index;
	//
	//	for (unsigned m = 0; m < 8; m++)
	//		for (unsigned e = 0; e < 3; e++)
	//			for (unsigned j = 0; j < 2; j++)
	//				χ(k_index, m, e, j) = (real) 0.0;
	//
	//	// We know the matrices, see definition.
	//	χ(k_index, 0, 0, 0) = (real) 1.0;	// w1 -> constant flux
	//	χ(k_index, 1, 1, 0) = (real) 1.0;	// w2 -> constant flux
	//	χ(k_index, 2, 2, 0) = (real) 1.0;	// w3 -> constant flux
	//
	//	χ(k_index, 3, 0, 1) = (real) 1.0;	// w4 -> linear flux
	//	χ(k_index, 4, 1, 1) = (real) 1.0;	// w5 -> linear flux
	//	χ(k_index, 5, 2, 1) = (real) 1.0;	// w6 -> linear flux
	//
	//	// w7 , w8 zero flux -> zeroed above
	//
	//
	//	//for (unsigned m = 0; m < 8; m++) {
	//
	//	//	Matrix<real> integral(3, 2);
	//	//	integral.setZero();
	//
	//	//	for (unsigned El = 0; El < 3; El++)
	//	//		for (unsigned j = 0; j < 2; j++)
	//	//			integral(El, j) = χ(k_index, m, El, j);
	//
	//	//	std::cout << integral << std::endl;
	//
	//	//}
	//}



	// Quadrature weights and points on reference segment [-1,1]
	gauss_quadrature_1D const quad(quadrature_order);
	unsigned const num_quad_points = quad.number_of_points;

	Matrix <real> normals(2, 3);
	Matrix<real> parametrization(2, 3);
	Matrix<real> parametrizationDerivative(2, 8);
	Matrix<real> basisRaviartThomas(2, 8);
	Matrix<real> basisEdgePolynomial(2, 3);


	evaluate_edge_normal(normals);


	for (unsigned k = 0; k < nk; k++) {


		t_pointer const K = mesh->get_triangle(k);
		unsigned const k_index = K->index;


		//--------------------------------
		//Matrix<real> JF(2, 2);

		//v_pointer const a = K->vertices[0];
		//v_pointer const b = K->vertices[1];
		//v_pointer const c = K->vertices[2];

		//real const x0 = (real)a->x;
		//real const y0 = (real)a->y;

		//real const x1 = (real)b->x;
		//real const y1 = (real)b->y;

		//real const x2 = (real)c->x;
		//real const y2 = (real)c->y;

		//JF(0, 0) = x1 - x0;
		//JF(0, 1) = x2 - x0;
		//JF(1, 0) = y1 - y0;
		//JF(1, 1) = y2 - y0;
		//--------------------------------

		real orientations[3];

		for (unsigned e = 0; e < 3; e++)
			orientations[e] = edgeOrientation(k_index, 0, e, 0);


		for (unsigned n = 0; n < num_quad_points; n++) {

			for (unsigned El = 0; El < 3; El++) {


				real const orientation = orientations[El];


				real const a = (real) 0.0;
				real const b = (real) El != 0 ? 1.0 : sqrt(2.0);

				real const c = (real) (b - a) / 2.0;
				real const d = (real) (b + a) / 2.0;

				real const x = (real) quad.points[n] * c + d;
				real const w = (real) quad.weigths[n] * c;

				//incorporate orientation into edge parametrization?

				evaluate_edge_parametrization(x, parametrization);
				evaluate_edge_parametrization_derivative(x, parametrizationDerivative);

				Vector<real> const normal = normals.getColumn(El);
				Vector<real> const r = parametrization.getColumn(El);
				Vector<real> const dr = parametrizationDerivative.getColumn(El);

				real const s = r(0);
				real const t = r(1);
				real const drNorm = dr.norm();

				// Values of basis function on reference triangle' edge
				evaluate_raviartthomas_basis(s, t, basisRaviartThomas, orientations);
				evaluate_edge_polynomial_basis(x, basisEdgePolynomial, orientation*orientation); // get rid of that orientation. Every time a compute it with the global edge orientation

				Vector<real> const edgeBasis = basisEdgePolynomial.getColumn(El);



				for (unsigned m = 0; m < 8; m++) {


					Vector<real> const Wm = basisRaviartThomas.getColumn(m);
					//Vector<real> const Wm = JF * basisRaviartThomas.getColumn(m);
					//real const dotProduct = dot(Wm, normal);
					real const dotProduct = orientation * dot(Wm, normal);

					for (unsigned j = 0; j < 2; j++) {

						real const varPhij = edgeBasis(j);

						χ(k_index, m, El, j) = χ(k_index, m, El, j) + w * dotProduct * varPhij * drNorm;

					}
				}
			}
		}

		for (unsigned El = 0; El < 3; El++)
			for (unsigned m = 0; m < 8; m++)
				for (unsigned j = 0; j < 2; j++)
					χ.setCoeff(k_index, m, El, j) = abs(χ(k_index, m, El, j)) < INTEGRAL_PRECISION ? 0.0 : χ(k_index, m, El, j);


		//for (unsigned m = 0; m < 8; m++) {

		//	Matrix<real> integral(3, 2);
		//	integral.setZero();

		//	for (unsigned El = 0; El < 3; El++)
		//		for (unsigned j = 0; j < 2; j++)
		//			integral(El, j) = χ(k_index, m, El, j);

		//	std::cout << integral << std::endl;

		//}

	}

};
*/


/*
void solver::assemble_δ() {


	unsigned const nk = _nk;

	// Quadrature weights and points on reference segment [-1,1]
	gauss_quadrature_1D const quad(quadrature_order);
	unsigned const num_quad_points = quad.number_of_points;

	Matrix <real> normal(2, 3);
	Matrix<real> integral(8, 3);
	Matrix<real> parametrization(2, 3);
	Matrix<real> parametrizationDerivative(2, 8);
	Matrix<real> basisRaviartThomas(2, 8);
	Vector<real> basisPolynomial(3);


	evaluate_edge_normal(normal);


	for (unsigned k = 0; k < nk; k++) {


		t_pointer const K = mesh->get_triangle(k);
		unsigned const k_index = K->index;

		real const orientation[3] = { 1.0, 1.0, 1.0 };



		for (unsigned n = 0; n < num_quad_points; n++) {

			for (unsigned El = 0; El < 3; El++) {


				real const a = (real) 0.0;
				real const b = (real) El != 0 ? 1.0 : sqrt(2.0);

				real const c = (real) (b - a) / 2.0;
				real const d = (real) (b + a) / 2.0;

				real const ksi = (real) quad.points[n] * c + d;
				real const weight = (real) quad.weigths[n] * c;

				evaluate_edge_parametrization(ksi, parametrization);
				evaluate_edge_parametrization_derivative(ksi, parametrizationDerivative);

				Vector<real> const v = normal.getColumn(El);
				Vector<real> const r = parametrization.getColumn(El);
				Vector<real> const dr = parametrizationDerivative.getColumn(El);

				real const s = r(0);
				real const t = r(1);
				real const drNorm = dr.norm();

				// Values of basis function on reference triangle' edge
				evaluate_raviartthomas_basis(s, t, basisRaviartThomas, orientation);
				evaluate_polynomial_basis(s, t, basisPolynomial);


				integral.setZero();

				for (unsigned j = 0; j < 8; j++) {


					Vector<real> const Wj = basisRaviartThomas.getColumn(j);
					real const dotProduct = dot(Wj, v);


					for (unsigned m = 0; m < 3; m++) {


						real const Phim = basisPolynomial(m);

						δ.setCoeff(k_index, j, El, m) = δ(k_index, j, El, m) + weight * dotProduct * Phim * drNorm;

					}
				}
			}
		}

		for (unsigned El = 0; El < 3; El++)
			for (unsigned j = 0; j < 8; j++)
				for (unsigned m = 0; m < 3; m++)
					δ.setCoeff(k_index, j, El, m) = abs(δ(k_index, j, El, m)) < INTEGRAL_PRECISION ? 0.0 : δ(k_index, j, El, m);


		//for (unsigned El = 0; El < 3; El++) {

		//	for (unsigned j = 0; j < 8; j++)
		//		for (unsigned m = 0; m < 3; m++)
		//			integral(j, m) = δ(k_index, j, El, m);

		//	std::cout << integral << std::endl;

		//}
	}



	//for (unsigned k = 0; k < nk; k++) {
	//
	//
	//	t_pointer const K = mesh->get_triangle(k);
	//	unsigned const k_index = K->index;
	//
	//	real const orient[3] = { 1.0, 1.0, 1.0 };
	//
	//
	//	for (unsigned j = 0; j < 8; j++) {
	//
	//		for (unsigned Ei = 0; Ei < 3; Ei++) {
	//
	//			for (unsigned m = 0; m < 3; m++) {
	//
	//
	//				real integral = 0.0;
	//
	//				for (unsigned n = 0; n < num_quad_points; n++) {
	//
	//
	//					// Quadrature point on reference segment [-1,1] and its weights
	//					real const ksi = quad.points[n];
	//					real const weight = quad.weigths[n];
	//
	//
	//					if (Ei == 0) {
	//
	//
	//						// E1 : ksi points need to be transformed to  ksi = [ 0 , sqrt(2) ]
	//						const double a = 0.0;
	//						const double b = sqrt(2.0);
	//
	//						const double c = (b - a) / 2.0;
	//						const double d = (b + a) / 2.0;
	//
	//						// Quadrature points and weight on [ 0 , sqrt(2) ]
	//						double const ksi = quad.points[n] * c + d;
	//						double const weight = quad.weigths[n] * c;
	//
	//
	//						// coordinates s,t on the edge E1 parametrized with ksi
	//						double const si = r11(ksi);
	//						double const ti = r12(ksi);
	//
	//						double const wj1 = evaluate_basis_reference2(j, 0, si, ti, orient);
	//						double const wj2 = evaluate_basis_reference2(j, 1, si, ti, orient);
	//
	//						double const d_r11 = dr11(ksi);
	//						double const d_r12 = dr12(ksi);
	//
	//						double const Phi_m = Phi(m, si, ti);
	//
	//
	//						integral += weight * (wj1 * d_r11 + wj2 * d_r12) * Phi_m;
	//
	//					}
	//					else if (Ei == 1) {
	//
	//
	//						// E2 : ksi points need to be transformed to  ksi = [ 0 , 1 ]
	//						const double a = 0.0;
	//						const double b = 1.0;
	//
	//						const double c = (b - a) / 2.0;
	//						const double d = (b + a) / 2.0;
	//
	//						// Quadrature points and weight on [ 0 , 1 ]
	//						double const ksi = quad.points[n] * c + d;
	//						double const weight = quad.weigths[n] * c;
	//
	//
	//						// coordinates s,t on the edge E2 parametrized with ksi
	//						double const si = r21(ksi);
	//						double const ti = r22(ksi);
	//
	//						double const wj1 = evaluate_basis_reference2(j, 0, si, ti, orient);
	//						double const wj2 = evaluate_basis_reference2(j, 1, si, ti, orient);
	//
	//						double const d_r21 = dr21(ksi);
	//						double const d_r22 = dr22(ksi);
	//
	//						double const Phi_m = Phi(m, si, ti);
	//
	//
	//						integral += weight * (wj1 * d_r21 + wj2 * d_r22) * Phi_m;
	//
	//					}
	//					else {
	//
	//
	//						// E3 : ksi points need to be transformed to  ksi = [ 0 , 1 ]
	//						const double a = 0.0;
	//						const double b = 1.0;
	//
	//						const double c = (b - a) / 2.0;
	//						const double d = (b + a) / 2.0;
	//
	//						// Quadrature points and weight on [ 0 , 1 ]
	//						double const ksi = quad.points[n] * c + d;
	//						double const weight = quad.weigths[n] * c;
	//
	//
	//						// coordinates s,t on the edge E3 parametrized with ksi
	//						double const si = r31(ksi);
	//						double const ti = r32(ksi);
	//
	//						double const wj1 = evaluate_basis_reference2(j, 0, si, ti, orient);
	//						double const wj2 = evaluate_basis_reference2(j, 1, si, ti, orient);
	//
	//						double const d_r31 = dr31(ksi);
	//						double const d_r32 = dr32(ksi);
	//
	//						double const Phi_m = Phi(m, si, ti);
	//
	//
	//						integral += weight * (wj1 * d_r31 + wj2 * d_r32) * Phi_m;
	//
	//					}
	//
	//				}
	//
	//				δ(k_index, j, Ei, m) = integral;
	//
	//			}
	//		}
	//	}
	//}


};*/


/*


		double basis[8][2];

		for (unsigned e = 0; e < 3; e++) {


			E = K->edges[e];

			unsigned const indexE = E->index;

			v_pointer a = E->a;
			v_pointer b = E->b;

			// Edge E is oriented ccw around elemnt K
			if (b != K->get_vertex_ccw(a)) {

				a = E->b;
				b = E->a;

			}



			double const x0 = a->x;
			double const y0 = a->y;

			double const x1 = b->x;
			double const y1 = b->y;

			double const dx = x1 - x0;
			double const dy = y1 - y0;

			// Normals
			//double const x_n = +dy;
			//double const y_n = -dx;

			//double const x_n = -dy;
			//double const y_n = +dx;

			double x_n;
			double y_n;

			// Normals on reference triangle
			if (indexE == 0) {

				x_n = 1.0 / sqrt(2.0);
				y_n = 1.0 / sqrt(2.0);

			}
			else if (indexE == 1) {

				x_n = -1.0;
				y_n = 0.0;;

			}
			else {

				x_n = 0.0;
				y_n = -1.0;;

			}

			for (unsigned n1 = 0; n1 < quadratureRule; n1++) {


				// Gauss Points, Weights on [-1,1]
				double const w = quad.weigths[n1];
				double const y = quad.points[n1];


				// Values of basis function from reference triangle to square.
				evaluate_basis2(s, t, K, basis, orient);


				for (unsigned m = 0; m < 8; m++) {

					// Phi1, Phi2, Phi3 are defined on reference triangle, therefore we need to input s,t and not x,y
					beta[m][0] += weight * divW[m] * phi1(s, t);
					beta[m][1] += weight * divW[m] * phi2(s, t);
					beta[m][2] += weight * divW[m] * phi3(s, t);

				}
			}


		}

		for (unsigned m = 0; m < 8; m++)
			for (unsigned e = 0; e < 3; e++)
				for (unsigned j = 0; j < 3; j++)
					χ(index, m, e, j);

					*/


/*
					unsigned const ne = _ne;


					for (unsigned e = 0; e < ne; e++) {


						e_pointer const E = mesh->get_edge(e);
						unsigned const e_index = E->index;

						t_pointer const K1 = E->neighbors[0];
						t_pointer const K2 = E->neighbors[1];


						E_MARKER const e_marker = E->marker;



						if (e_marker == E_MARKER::DIRICHLET) {

							M1.coeffRef(2 * e_index + 0, 2 * e_index + 0) = -1.0;
							M1.coeffRef(2 * e_index + 1, 2 * e_index + 1) = -1.0;

							M2.coeffRef(2 * e_index + 0, 2 * e_index + 0) = -1.0;
							M2.coeffRef(2 * e_index + 1, 2 * e_index + 1) = -1.0;

							continue;

						}

						real commonSum_j1_s1 = 0.0;
						real commonSum_j2_s1 = 0.0;
						real commonSum_j1_s2 = 0.0;
						real commonSum_j2_s2 = 0.0;

						if (K1) {


							unsigned const k_index = K1->index;

							// Loop over edges
							for (unsigned El = 0; El < 3; El++) {


								e_pointer const E_local = K1->edges[El];
								unsigned const e_index_local = K1->get_edge_index(E_local);
								unsigned const e_index_global = E_local->index;


								if (E_local == E) {

									// j = 1, s = 1
									for (unsigned i = 0; i < 8; i++)
										commonSum_j1_s1 += α(k_index, 0, i, LI(K1, E, 0))*χ(k_index, i, K1->get_edge_index(E), 0) / viscosities[k_index];
									// j = 2, s = 1
									for (unsigned i = 0; i < 8; i++)
										commonSum_j2_s1 += α(k_index, 0, i, LI(K1, E, 1))*χ(k_index, i, K1->get_edge_index(E), 0) / viscosities[k_index];
									// j = 1, s = 2
									for (unsigned i = 0; i < 8; i++)
										commonSum_j1_s2 += α(k_index, 0, i, LI(K1, E, 0))*χ(k_index, i, K1->get_edge_index(E), 1) / viscosities[k_index];
									// j = 2, s = 2
									for (unsigned i = 0; i < 8; i++)
										commonSum_j2_s2 += α(k_index, 0, i, LI(K1, E, 1))*χ(k_index, i, K1->get_edge_index(E), 1) / viscosities[k_index];

									continue;

								}

								real ACHI_j1_s1 = 0.0;
								real ACHI_j1_s2 = 0.0;
								real ACHI_j2_s1 = 0.0;
								real ACHI_j2_s2 = 0.0;

								for (unsigned i = 0; i < 8; i++) {

									// j = 1, s = 1
									ACHI_j1_s1 += α(k_index, 0, i, LI(K1, E_local, 0))*χ(k_index, i, e_index_local, 0);
									// j = 1, s = 2
									ACHI_j1_s2 += α(k_index, 0, i, LI(K1, E_local, 0))*χ(k_index, i, e_index_local, 1);

									// j = 2, s = 1
									ACHI_j2_s1 += α(k_index, 0, i, LI(K1, E_local, 1))*χ(k_index, i, e_index_local, 0);
									// j = 2, s = 2
									ACHI_j2_s2 += α(k_index, 0, i, LI(K1, E_local, 1))*χ(k_index, i, e_index_local, 1);

								}


								M1.coeffRef(2 * e_index + 0, 2 * e_index_global) = ACHI_j1_s1 / viscosities[k_index];
								M1.coeffRef(2 * e_index + 1, 2 * e_index_global) = ACHI_j2_s1 / viscosities[k_index];

								M2.coeffRef(2 * e_index + 0, 2 * e_index_global) = ACHI_j1_s2 / viscosities[k_index];
								M2.coeffRef(2 * e_index + 1, 2 * e_index_global) = ACHI_j2_s2 / viscosities[k_index];

							}

						}
						if (K2) {


							unsigned const k_index = K2->index;

							// Loop over edges
							for (unsigned El = 0; El < 3; El++) {


								e_pointer const E_local = K2->edges[El];
								unsigned const e_index_local = K2->get_edge_index(E_local);
								unsigned const e_index_global = E_local->index;


								if (E_local == E) {

									// j = 1, s = 1
									for (unsigned i = 0; i < 8; i++)
										commonSum_j1_s1 += α(k_index, 0, i, LI(K2, E, 0))*χ(k_index, i, K2->get_edge_index(E), 0) / viscosities[k_index];
									// j = 2, s = 1
									for (unsigned i = 0; i < 8; i++)
										commonSum_j2_s1 += α(k_index, 0, i, LI(K2, E, 1))*χ(k_index, i, K2->get_edge_index(E), 0) / viscosities[k_index];
									// j = 1, s = 2
									for (unsigned i = 0; i < 8; i++)
										commonSum_j1_s2 += α(k_index, 0, i, LI(K2, E, 0))*χ(k_index, i, K2->get_edge_index(E), 1) / viscosities[k_index];
									// j = 2, s = 2
									for (unsigned i = 0; i < 8; i++)
										commonSum_j2_s2 += α(k_index, 0, i, LI(K2, E, 1))*χ(k_index, i, K2->get_edge_index(E), 1) / viscosities[k_index];


									continue;

								}

								real ACHI_j1_s1 = 0.0;
								real ACHI_j1_s2 = 0.0;
								real ACHI_j2_s1 = 0.0;
								real ACHI_j2_s2 = 0.0;

								for (unsigned i = 0; i < 8; i++) {

									// j = 1, s = 1
									ACHI_j1_s1 += α(k_index, 0, i, LI(K2, E_local, 0))*χ(k_index, i, e_index_local, 0);
									// j = 1, s = 2
									ACHI_j1_s2 += α(k_index, 0, i, LI(K2, E_local, 0))*χ(k_index, i, e_index_local, 1);

									// j = 2, s = 1
									ACHI_j2_s1 += α(k_index, 0, i, LI(K2, E_local, 1))*χ(k_index, i, e_index_local, 0);
									// j = 2, s = 2
									ACHI_j2_s2 += α(k_index, 0, i, LI(K2, E_local, 1))*χ(k_index, i, e_index_local, 1);

								}


								M1.coeffRef(2 * e_index + 0, 2 * e_index_global) = ACHI_j1_s1 / viscosities[k_index];
								M1.coeffRef(2 * e_index + 1, 2 * e_index_global) = ACHI_j2_s1 / viscosities[k_index];

								M2.coeffRef(2 * e_index + 0, 2 * e_index_global) = ACHI_j1_s2 / viscosities[k_index];
								M2.coeffRef(2 * e_index + 1, 2 * e_index_global) = ACHI_j2_s2 / viscosities[k_index];

							}

						}



						M1.coeffRef(2 * e_index + 0, 2 * e_index) = commonSum_j1_s1;
						M1.coeffRef(2 * e_index + 1, 2 * e_index) = commonSum_j2_s1;

						M2.coeffRef(2 * e_index + 0, 2 * e_index) = commonSum_j1_s2;
						M2.coeffRef(2 * e_index + 1, 2 * e_index) = commonSum_j2_s2;


					}
					*/



/*
					for (unsigned e = 0; e < ne; e++) {


						e_pointer const E = mesh->get_edge(e);

						t_pointer const K1 = E->neighbors[0];
						t_pointer const K2 = E->neighbors[1];

						unsigned const e_index = E->index;
						E_MARKER const e_marker = E->marker;



						if (e_marker == E_MARKER::DIRICHLET) {

							M_j1_s1.coeffRef(e_index, 2 * e_index) = -1.0;
							M_j1_s2.coeffRef(e_index, 2 * e_index) = 0.0;

							M_j2_s1.coeffRef(e_index, 2 * e_index + 0) = 0.0;
							M_j2_s2.coeffRef(e_index, 2 * e_index + 1) = 0.0;

							continue;

						}


						real sum_j1_s1 = 0.0;
						real sum_j1_s2 = 0.0;
						real sum_j2_s1 = 0.0;
						real sum_j2_s2 = 0.0;

						for (unsigned k = 0; k < 1; k++) {


							t_pointer const K = E->neighbors[k];

							if (!K)
								continue;

							unsigned const k_index = K->index;


							for (unsigned l = 0; l < 3; l++) {


								e_pointer const El = K->edges[l];

								unsigned const e_index_local = K->get_edge_index(El);
								unsigned const e_index_global = El->index;


								real ACHI_j1_s1 = 0.0;
								real ACHI_j1_s2 = 0.0;
								real ACHI_j2_s1 = 0.0;
								real ACHI_j2_s2 = 0.0;

								for (unsigned i = 0; i < 8; i++) {

									// j = 1, s = 1
									ACHI_j1_s1 += α(k_index, 0, i, LI(K, El, 0))*χ(k_index, i, e_index_local, 0);
									// j = 1, s = 2
									ACHI_j1_s2 += α(k_index, 0, i, LI(K, El, 0))*χ(k_index, i, e_index_local, 1);

									// j = 2, s = 1
									ACHI_j2_s1 += α(k_index, 0, i, LI(K, El, 1))*χ(k_index, i, e_index_local, 0);
									// j = 2, s = 2
									ACHI_j2_s2 += α(k_index, 0, i, LI(K, El, 1))*χ(k_index, i, e_index_local, 1);

									//ACHI1 += α(t_index, 0, i, LI(K1, El, 0))*χ(t_index, i, l, 0);
									//ACHI2 += α(t_index, 0, i, LI(K1, El, 1))*χ(t_index, i, l, 1);

								}

								if (El == E) {

									sum_j1_s1 += ACHI_j1_s1 / viscosities[k_index];
									sum_j1_s2 += ACHI_j1_s2 / viscosities[k_index];

									sum_j2_s1 += ACHI_j2_s1 / viscosities[k_index];
									sum_j2_s2 += ACHI_j2_s2 / viscosities[k_index];

									continue;

								}


								M_j1_s1.coeffRef(e_index, e_index_global) = sum_j1_s1 / viscosities[k_index];
								M_j1_s2.coeffRef(e_index, e_index_global) = sum_j1_s2 / viscosities[k_index];

								M_j2_s1.coeffRef(e_index, e_index_global) = sum_j2_s1 / viscosities[k_index];
								M_j2_s2.coeffRef(e_index, e_index_global) = sum_j2_s2 / viscosities[k_index];


								M1.coeffRef(2 * e_index + 0, 2 * e_index_global) = sum_j1_s1 / viscosities[k_index];
								M1.coeffRef(2 * e_index + 1, 2 * e_index_global) = sum_j2_s1 / viscosities[k_index];

								M2.coeffRef(2 * e_index + 0, 2 * e_index_global) = sum_j1_s2 / viscosities[k_index];
								M2.coeffRef(2 * e_index + 1, 2 * e_index_global) = sum_j2_s2 / viscosities[k_index];

							}

						}

						M_j1_s1.coeffRef(e_index, e_index) = sum_j1_s1;
						M_j1_s2.coeffRef(e_index, e_index) = sum_j1_s2;

						M_j2_s1.coeffRef(e_index, e_index) = sum_j2_s1;
						M_j2_s2.coeffRef(e_index, e_index) = sum_j2_s2;


						M1.coeffRef(2 * e_index + 0, e_index) = sum_j1_s1;
						M1.coeffRef(2 * e_index + 1, e_index) = sum_j2_s1;

						M2.coeffRef(2 * e_index + 0, e_index) = sum_j1_s2;
						M2.coeffRef(2 * e_index + 1, e_index) = sum_j2_s2;

					}
					*/




/*
							if (K1) {


								unsigned const k_index = K1->index;

								for (unsigned m = 0; m < 3; m++) {


									real AB1 = 0.0;
									real AB2 = 0.0;

									for (unsigned i = 0; i < 8; i++) {

										AB1 += α(k_index, 0, i, LI(K1, E, 0))*β(k_index, 0, i, m);
										AB2 += α(k_index, 0, i, LI(K1, E, 1))*β(k_index, 0, i, m);

									}

									real const val1 = AB1 / viscosities[k_index];
									real const val2 = AB2 / viscosities[k_index];

									R1.coeffRef(e_index, 3 * k_index + m) = abs(val1) < INTEGRAL_PRECISION ? 0.0 : val1;
									R2.coeffRef(e_index, 3 * k_index + m) = abs(val2) < INTEGRAL_PRECISION ? 0.0 : val2;

								}
							}
							if (K2) {


								unsigned const k_index = K2->index;

								for (unsigned m = 0; m < 3; m++) {


									double AB1 = 0.0;
									double AB2 = 0.0;

									for (unsigned i = 0; i < 8; i++) {

										AB1 += α(k_index, 0, i, LI(K2, E, 0))*β(k_index, 0, i, m);
										AB2 += α(k_index, 0, i, LI(K2, E, 1))*β(k_index, 0, i, m);

									}


									real const val1 = AB1 / viscosities[k_index];
									real const val2 = AB2 / viscosities[k_index];

									R1.coeffRef(e_index, 3 * k_index + m) = abs(val1) < INTEGRAL_PRECISION ? 0.0 : val1;
									R2.coeffRef(e_index, 3 * k_index + m) = abs(val2) < INTEGRAL_PRECISION ? 0.0 : val2;

								}
}*/



/*
void solver::assemble_χ() {


	// Quadrature weights and points on reference segment [-1,1]
	gauss_quadrature_1D const quad(quadrature_order);
	unsigned const num_quad_points = quad.number_of_points;


	//Matrix<real> normals(2, 3);
	//Vector<real> parametrization(2);
	//Vector<real> parametrizationDerivative(2);
	//Matrix<real> basisRaviartThomas(2, 8);
	//Vector<real> basisEdgePolynomial(2);

	Eigen::MatrixXd normals(2, 3);
	Eigen::VectorXd parametrization(2);
	Eigen::VectorXd parametrizationDerivative(2);
	Eigen::MatrixXd basisRaviartThomas(2, 8);
	Eigen::VectorXd basisEdgePolynomial(2);



	χ.setZero();

	for (unsigned k = 0; k < nk; k++) {


		t_pointer const K = mesh->get_triangle(k);
		unsigned const k_index = K->index;


		//--------------------------------

		//Matrix<real> JF(2, 2);
		//
		//v_pointer const a = K->vertices[0];
		//v_pointer const b = K->vertices[1];
		//v_pointer const c = K->vertices[2];
		//
		//real const x0 = (real)a->x;
		//real const y0 = (real)a->y;
		//
		//real const x1 = (real)b->x;
		//real const y1 = (real)b->y;
		//
		//real const x2 = (real)c->x;
		//real const y2 = (real)c->y;
		//
		//real const detJF = (real) ((x1 - x0)*(y2 - y0) - (x2 - x0)*(y1 - y0));
		//
		//JF(0, 0) = x1 - x0;
		//JF(0, 1) = x2 - x0;
		//JF(1, 0) = y1 - y0;
		//JF(1, 1) = y2 - y0;

		//--------------------------------


		real orientations[3] = { 1.0 , 1.0 , 1.0 };
		real orientations3[3] = { 1.0 , 1.0 , 1.0 };
		//real orientations2[3] = { 1.0 , 1.0 , 1.0 };

		for (unsigned e = 0; e < 3; e++)
			orientations[e] = edgeOrientation(k_index, 0, e, 0);
		for (unsigned e = 0; e < 3; e++)
			orientations3[e] = edgeOrientation2(k_index, 0, e, 0);


		evaluate_edge_normal(normals, orientations3);
		//evaluate_edge_normal(normals, orientations3);


		for (unsigned n = 0; n < num_quad_points; n++) {

			for (unsigned El = 0; El < 3; El++) {


				real const orientation = orientations[El];
				real const orientation2 = orientations3[El];


				real const a = (real) 0.0;
				real const b = (real) El != 0 ? 1.0 : sqrt(2.0);

				//real a;
				//real b;
				//if (orientation == 1.0) {
				//	a = (real) 0.0;
				//	b = (real) El != 0 ? 1.0 : sqrt(2.0);
				//}
				//else {
				//	b = (real) 0.0;
				//	a = (real) El != 0 ? 1.0 : sqrt(2.0);
				//}


				real const c = (real)(b - a) / 2.0;
				real const d = (real)(b + a) / 2.0;

				real const x = (real)quad.points[n] * c + d;
				real const w = (real)quad.weigths[n] * c;

				//for (unsigned i = 0; i < num_quad_points; i++) {
				//
				//	real y = (real)quad.points[i] * c + d;
				//	evaluate_edge_parametrization(y, El, parametrization, orientation);
				//
				//	real const s = parametrization(0);
				//	real const t = parametrization(1);
				//
				//	std::cout << s << "    " << t << std::endl;
				//
				//}



				//real orientation2 = 1.0;
				evaluate_edge_parametrization(x, El, parametrization, orientation);
				evaluate_edge_parametrization_derivative(x, El, parametrizationDerivative, orientation);

				//Vector<real> const normal = normals.getColumn(El);
				//Vector<real> const r = parametrization;
				//Vector<real> const dr = parametrizationDerivative;

				Eigen::VectorXd const normal = normals.col(El);
				Eigen::VectorXd const r = parametrization;
				Eigen::VectorXd const dr = parametrizationDerivative;


				real const s = r(0);
				real const t = r(1);
				//real const drNorm = dr.norm();
				real const drNorm = 1.0;

				// Values of basis function on reference triangle' edge
				evaluate_raviartthomas_basis(s, t, basisRaviartThomas, orientations);
				evaluate_edge_polynomial_basis(x, El, basisEdgePolynomial, orientation);

				//Vector<real> const edgeBasis = basisEdgePolynomial;
				Eigen::VectorXd const edgeBasis = basisEdgePolynomial;


				for (unsigned m = 0; m < 8; m++) {


					//Vector<real> const Wm = basisRaviartThomas.getColumn(m);
					Eigen::VectorXd const Wm = basisRaviartThomas.col(m);

					//real const dotProduct = dot(Wm, normal);
					real const dotProduct = Wm.dot(normal);



					for (unsigned s = 0; s < 2; s++) {

						//real const varPhis = s == 1 ? orientation * edgeBasis(s) : edgeBasis(s);
						//real const varPhis = s == 1 ? orientation2 * edgeBasis(s) : edgeBasis(s);
						//real const varPhis = orientation * edgeBasis(s);
						real const varPhis = edgeBasis(s);

						χ.setCoeff(k_index, m, El, s) = χ(k_index, m, El, s) + w * dotProduct * varPhis * drNorm;

					}
				}
			}
		}

		for (unsigned El = 0; El < 3; El++)
			for (unsigned m = 0; m < 8; m++)
				for (unsigned j = 0; j < 2; j++)
					χ.setCoeff(k_index, m, El, j) = abs(χ(k_index, m, El, j)) < INTEGRAL_PRECISION ? 0.0 : χ(k_index, m, El, j);


		//for (unsigned m = 0; m < 8; m++) {
		////
		//	Matrix<real> integral(3, 2);
		//	integral.setZero();
		////
		//	for (unsigned El = 0; El < 3; El++)
		//		for (unsigned j = 0; j < 2; j++)
		//			integral(El, j) = χ(k_index, m, El, j);
		////
		//	std::cout << integral << std::endl;
		////
		//}

	}

};
*/


/*
	real const X = -0.5;

	std::cout << "----------------- Velocity on Neumann Edges -----------------" << std::endl;
	for (unsigned e = 0; e < ne; e++) {

		//e_pointer const E = mesh->get_edge(e);

		//t_pointer const K1 = E->neighbors[0];
		//t_pointer const K2 = E->neighbors[1];

		//if (K1 && E->marker == E_MARKER::NEUMANN) {

		//	std::cout << "K index: " << K1->index << std::endl;
		//	std::cout << "dof0: " << velocities(K1->index, 0, LI(K1, E, 0), 0) << std::endl;
		//	std::cout << "dof1: " << velocities(K1->index, 0, LI(K1, E, 1), 0) << std::endl;
		//	std::cout << "norm: " << velocity_in_normal_direction(K1, E, X) << std::endl << std::endl;

		//}
		//if (K2 && E->marker == E_MARKER::NEUMANN) {


		//	std::cout << "K index: " << K2->index << std::endl;
		//	std::cout << "dof0: " << velocities(K2->index, 0, LI(K2, E, 0), 0) << std::endl;
		//	std::cout << "dof1: " << velocities(K2->index, 0, LI(K2, E, 1), 0) << std::endl;
		//	std::cout << "norm: " << velocity_in_normal_direction(K2, E, X) << std::endl << std::endl;

		//}

	}
	std::cout << "----------------- Velocity on Dirichlet Edges -----------------" << std::endl;
	for (unsigned e = 0; e < ne; e++) {

		//e_pointer const E = mesh->get_edge(e);

		//t_pointer const K1 = E->neighbors[0];
		//t_pointer const K2 = E->neighbors[1];

		//if (K1 && E->marker == E_MARKER::DIRICHLET) {

		//	std::cout << "K index: " << K1->index << std::endl;
		//	std::cout << "dof0: " << velocities(K1->index, 0, LI(K1, E, 0), 0) << std::endl;
		//	std::cout << "dof1: " << velocities(K1->index, 0, LI(K1, E, 1), 0) << std::endl;
		//	std::cout << "norm: " << velocity_in_normal_direction(K1, E, X) << std::endl << std::endl;

		//}
		//if (K2 && E->marker == E_MARKER::DIRICHLET) {

		//	std::cout << "K index: " << K2->index << std::endl;
		//	std::cout << "dof0: " << velocities(K2->index, 0, LI(K2, E, 0), 0) << std::endl;
		//	std::cout << "dof1: " << velocities(K2->index, 0, LI(K2, E, 1), 0) << std::endl;
		//	std::cout << "norm: " << velocity_in_normal_direction(K2, E, X) << std::endl << std::endl;

		//}

	}
	std::cout << "----------------- Velocity on Inner Edges -----------------" << std::endl;
	for (unsigned e = 0; e < ne; e++) {

		e_pointer const E = mesh->get_edge(e);

		t_pointer const K1 = E->neighbors[0];
		t_pointer const K2 = E->neighbors[1];

		if (K1 && K2) {

			std::cout << K1->index << "  " << K2->index << std::endl;
			std::cout << "dof0: " << velocities(K1->index, 0, LI(K1, E, 0), 0) << "         " << velocities(K2->index, 0, LI(K2, E, 0), 0) << std::endl;
			std::cout << "dof1: " << velocities(K1->index, 0, LI(K1, E, 1), 0) << "         " << velocities(K2->index, 0, LI(K2, E, 1), 0) << std::endl;
			std::cout << "norm: " << velocity_in_normal_direction(K1, E, X) << "         " << velocity_in_normal_direction(K2, E, X) << std::endl << std::endl;

		}

	}
	std::cout << std::endl << std::endl;
	std::cout << "----------------- Value of continuity of normal component -----------------" << std::endl;
	for (unsigned e = 0; e < ne; e++) {

		//e_pointer const E = mesh->get_edge(e);

		//t_pointer const K1 = E->neighbors[0];
		//t_pointer const K2 = E->neighbors[1];

		//real const x = 0.7;

		//if (K1 && K2) {

		//	//std::cout << velocities(K1->index, 0, LI(K1, E, 0), 0) << "         " << velocities(K2->index, 0, LI(K2, E, 0), 0) << std::endl;
		//	std::cout << "K index: " << K1->index << "  " << K2->index << std::endl;
		//	std::cout << continuity_of_normal_component(K1, E, x) << "         " << continuity_of_normal_component(K2, E, x) << std::endl << std::endl;
		//	//std::cout << velocities(K1->index, 0, LI(K1, E, 1), 0) << "         " << velocities(K2->index, 0, LI(K2, E, 1), 0) << std::endl << std::endl;

		//}

	}

	std::cout << std::endl << std::endl;
	*/


/*
real solver::upwindConcentration(real const s, real const t, t_pointer const K, Vector<real> const normal, unsigned const El, CoeffMatrix<1, 3, 1> & ksi) {


	real const time = (nt + 1) * dt;
	unsigned const k_index = K->index;

	//Matrix<real> basisRaviartThomas(2, 8);
	//Vector<real> basisPolynomial(3);

	Eigen::MatrixXd basisRaviartThomas(2, 8);
	Eigen::VectorXd basisPolynomial(3);


	// Values of basis function on reference triangle' edge at the gauss point (s,t)
	evaluate_polynomial_basis(s, t, basisPolynomial);
	evaluate_raviartthomas_basis(s, t, basisRaviartThomas);



	real dotVelocityNormal = 0.0;
	real concentration = 0.0;

	if (K->edges[El]->marker == E_MARKER::NEUMANN) {

		//dotVelocityNormal = NEUMANN_GAMMA_Q_velocity(K->edges[El], time);
		dotVelocityNormal = NEUMANN_GAMMA_Q_velocity(K->edges[El], time);

	}
	else {


		//Matrix<real> JF(2, 2);
		Eigen::MatrixXd JF(2, 2);

		v_pointer const a = K->vertices[0];
		v_pointer const b = K->vertices[1];
		v_pointer const c = K->vertices[2];

		real const x0 = (real)a->x;
		real const y0 = (real)a->y;

		real const x1 = (real)b->x;
		real const y1 = (real)b->y;

		real const x2 = (real)c->x;
		real const y2 = (real)c->y;

		real const detJF = (real)abs((x1 - x0)*(y2 - y0) - (x2 - x0)*(y1 - y0));

		//JF(0, 0) = x1 - x0;
		//JF(0, 1) = x2 - x0;
		//JF(1, 0) = y1 - y0;
		//JF(1, 1) = y2 - y0;

		//Matrix<real> itJF = JF.inverse();
		//itJF.transposeInPlace();
		//Vector<real> origNormal = itJF * normal / (itJF * normal).norm();



		Eigen::Vector2d normal_eigen;
		normal_eigen(0) = normal(0);
		normal_eigen(1) = normal(1);

		JF(0, 0) = x1 - x0;
		JF(0, 1) = x2 - x0;
		JF(1, 0) = y1 - y0;
		JF(1, 1) = y2 - y0;

		Eigen::MatrixXd itJF = (JF.inverse()).transpose();
		Eigen::VectorXd origNormal = itJF * normal_eigen / (itJF * normal_eigen).norm();


		for (unsigned dof = 0; dof < 2; dof++) {

			//evaluate_raviartthomas_basis(s, t, basisRaviartThomas, orientations, dof);

			unsigned const j = LI(K, K->edges[El], dof);
			//real const val = velocities(k_index, 0, j, 0) * dot(normal, basisRaviartThomas.getColumn(j));

			//real const val = velocities(k_index, 0, dof, 0) * dot(normal, basisRaviartThomas.getColumn(dof));


			// Maybe try complete transformation both for w and n. But it is the same up to positive constant  idetJF / (itJF * normal).norm()
			//real const val = velocities(k_index, 0, j, 0) * dot(origNormal, JF * basisRaviartThomas.getColumn(j)) * idetJF;
			//real const val2 = velocities(k_index, 0, dof, 0) * dot(origNormal, JF * basisRaviartThomas.getColumn(dof)) * idetJF;

			real const val = velocities(k_index, 0, j, 0) * origNormal.dot(JF * basisRaviartThomas.col(j)) / detJF;
			//real const val2 = velocities(k_index, 0, dof, 0) * dot(origNormal, JF * basisRaviartThomas.getColumn(dof)) * idetJF;

			dotVelocityNormal += val;

		}
	}



	if (dotVelocityNormal >= 0.0) {

		
		////v_pointer const a = K->vertices[0];
		////v_pointer const b = K->vertices[1];
		////v_pointer const c = K->vertices[2];
		//
		////real const x0 = (real) a->x;
		////real const y0 = (real) a->y;
		//
		////real const x1 = (real) b->x;
		////real const y1 = (real) b->y;
		//
		////real const x2 = (real) c->x;
		////real const y2 = (real) c->y;
		//
		////real const idetJF = (real) 1.0 / abs((x1 - x0)*(y2 - y0) - (x2 - x0)*(y1 - y0));
		//
		////JF(0, 0) = x1 - x0;
		////JF(0, 1) = x2 - x0;
		////JF(1, 0) = y1 - y0;
		////JF(1, 1) = y2 - y0;
		//
		////Vector<real> const JPhim = JF * basisPolynomial;
		//
		////for (unsigned m = 0; m < 3; m++)
		////	std::cout << JPhim(m) << std::endl;
		////std::cout << std::endl;
		//
		////for (unsigned m = 0; m < 3; m++)
		////	concentration += idetJF * ksi(k_index, 0, m, 0) * JPhim(m);
		//
		

		//for (unsigned m = 0; m < 3; m++)
		//	concentration += ksi(k_index, 0, m, 0) * basisPolynomial(m);

concentration = ksi(k_index, 0, 0, 0) * basisPolynomial(0);

	}
	else {


	if (K->edges[El]->marker == E_MARKER::DIRICHLET)
		return barenblatt(s, t, time);

	//if (K->edges[El]->marker == E_MARKER::DIRICHLET)
	//	return DIRICHLET_GAMMA_P_concentration(K->edges[El], time);

	unsigned const kn_index = K->neighbors[El]->index;

	
	////v_pointer const a = Kn->vertices[0];
	////v_pointer const b = Kn->vertices[1];
	////v_pointer const c = Kn->vertices[2];
	//
	////real const x0 = (real) a->x;
	////real const y0 = (real) a->y;
	//
	////real const x1 = (real) b->x;
	////real const y1 = (real) b->y;
	//
	////real const x2 = (real) c->x;
	////real const y2 = (real) c->y;
	//
	////real const idetJF = (real) 1.0 / abs((x1 - x0)*(y2 - y0) - (x2 - x0)*(y1 - y0));
	//
	////JF(0, 0) = x1 - x0;
	////JF(0, 1) = x2 - x0;
	////JF(1, 0) = y1 - y0;
	////JF(1, 1) = y2 - y0;
	//
	////Vector<real> const JPhim = JF * basisPolynomial;
	//
	////for (unsigned m = 0; m < 3; m++)
	////	std::cout << JPhim(m) << std::endl;
	////std::cout << std::endl;
	//
	////for (unsigned m = 0; m < 3; m++)
	////	concentration += idetJF * ksi(kn_index, 0, m, 0) * JPhim(m);
	

	//for (unsigned m = 0; m < 3; m++)
	//	concentration += ksi(kn_index, 0, m, 0) * basisPolynomial(m);

	concentration = ksi(kn_index, 0, 0, 0) * basisPolynomial(0);

	}

	return concentration;

};
real solver::upwindConcentration(real const s, real const t, t_pointer const K, Eigen::VectorXd const normal, unsigned const El, CoeffMatrix<1, 3, 1> & ksi) {


	real const time = (nt + 1) * dt;
	unsigned const k_index = K->index;

	Eigen::VectorXd basisPolynomial(3);
	Eigen::MatrixXd basisRaviartThomas(2, 8);


	e_pointer const E = K->edges[El];


	evaluate_polynomial_basis(s, t, basisPolynomial);
	evaluate_raviartthomas_basis(s, t, basisRaviartThomas);


	real dotVelocityNormal = 0.0;
	real concentration = 0.0;

	if (E->marker == E_MARKER::NEUMANN) {

		dotVelocityNormal = NEUMANN_GAMMA_Q_velocity(E, time);

	}
	else {

		Eigen::MatrixXd JF(2, 2);

		v_pointer const a = K->vertices[0];
		v_pointer const b = K->vertices[1];
		v_pointer const c = K->vertices[2];

		real const x0 = (real)a->x;
		real const y0 = (real)a->y;

		real const x1 = (real)b->x;
		real const y1 = (real)b->y;

		real const x2 = (real)c->x;
		real const y2 = (real)c->y;

		real const detJF = (real)abs((x1 - x0)*(y2 - y0) - (x2 - x0)*(y1 - y0));

		JF(0, 0) = x1 - x0;
		JF(0, 1) = x2 - x0;
		JF(1, 0) = y1 - y0;
		JF(1, 1) = y2 - y0;

		Eigen::MatrixXd itJF = (JF.inverse()).transpose();
		Eigen::VectorXd origNormal = itJF * normal / (itJF * normal).norm();


		for (unsigned dof = 0; dof < 2; dof++) {


			unsigned const j = LI(K, E, dof);

			dotVelocityNormal += velocities(k_index, 0, j, 0) * origNormal.dot(JF * basisRaviartThomas.col(j)) / detJF;
			//dotVelocityNormal += velocities(k_index, 0, j, 0) * normal.dot(basisRaviartThomas.col(j));

		}
	}



	if (dotVelocityNormal >= 0.0) {

		//for (unsigned m = 0; m < 3; m++)
		//	concentration += ksi(k_index, 0, m, 0) * basisPolynomial(m);

		//// Total limiter
		concentration = ksi(k_index, 0, 0, 0) * basisPolynomial(0);

	}
	else {


		if (E->marker == E_MARKER::DIRICHLET)
			return barenblatt(s, t, time);

		//if (E->marker == E_MARKER::DIRICHLET)
		//	return DIRICHLET_GAMMA_P_concentration(E, time);

		unsigned const kn_index = K->neighbors[El]->index;


		//for (unsigned m = 0; m < 3; m++)
		//	concentration += ksi(kn_index, 0, m, 0) * basisPolynomial(m);

		//// Total limiter
		concentration = ksi(kn_index, 0, 0, 0) * basisPolynomial(0);

	}

	return concentration;

};



*/










//delta

/*

void solver::assemble_δ() {


	// Quadrature weights and points on reference segment [-1,1]
	gauss_quadrature_1D const quad(quadrature_order);
	unsigned const number_of_quadrature_points = quad.number_of_points;

	//Matrix<real> normals(2, 3);
	//Vector<real> parametrization(2);
	////Vector<real> parametrizationDerivative(2);
	//Matrix<real> basisRaviartThomas(2, 8);
	//Vector<real> basisPolynomial(3);


	Eigen::MatrixXd normals(2, 3);
	Eigen::VectorXd parametrization(2);
	//Eigen::VectorXd parametrizationDerivative(2);
	Eigen::MatrixXd basisRaviartThomas(2, 8);
	Eigen::VectorXd basisPolynomial(3);

	evaluate_edge_normal(normals);


	δ.setZero();

	for (unsigned k = 0; k < nk; k++) {


		t_pointer const K = mesh->get_triangle(k);
		unsigned const k_index = K->index;


		for (unsigned El = 0; El < 3; El++) {


			real const a = (real) 0.0;
			real const b = (real)(El != 0) ? 1.0 : sqrt(2.0);

			real const c = (real)(b - a) / 2.0;
			real const d = (real)(b + a) / 2.0;

			//Vector<real> const normal = normals.getColumn(El);
			Eigen::VectorXd const referenceNormal = normals.col(El);


			for (unsigned n = 0; n < number_of_quadrature_points; n++) {


				real const x = (real) quad.points[n] * c + d;
				real const w = (real) quad.weigths[n] * c;


				evaluate_edge_parametrization(x, El, parametrization);
				//evaluate_edge_parametrization_derivative(x, El, parametrizationDerivative);

				//Vector<real> const r = parametrization;
				//Vector<real> const dr = parametrizationDerivative;

				Eigen::VectorXd const r = parametrization;
				//Eigen::VectorXd const dr = parametrizationDerivative;

				real const s = r(0);
				real const t = r(1);
				//real const drNorm = dr.norm();
				real const drNorm = 1.0;


				evaluate_raviartthomas_basis(s, t, basisRaviartThomas);
				evaluate_polynomial_basis(s, t, basisPolynomial);


				real const tC = upwindConcentration(s, t, K, referenceNormal, El, n);


				for (unsigned m = 0; m < 3; m++) {


					real const Phim = basisPolynomial(m);

					for (unsigned j = 0; j < 8; j++) {

						//Vector<real> const Wj = basisRaviartThomas.getColumn(j);
						Eigen::VectorXd const Wj = basisRaviartThomas.col(j);

						//real const dotProduct = dot(Wj, normal);
						real const dotProduct = Wj.dot(referenceNormal);

						δ.setCoeff(k_index, m, El, j) = δ(k_index, m, El, j) + w * tC * dotProduct * Phim * drNorm;

					}
				}
			}
		}

		//for (unsigned El = 0; El < 3; El++)
		//	for (unsigned j = 0; j < 8; j++)
		//		for (unsigned m = 0; m < 3; m++)
		//			δ.setCoeff(k_index, m, El, j) = abs(δ(k_index, m, El, j)) < INTEGRAL_PRECISION ? 0.0 : δ(k_index, m, El, j);

	}

};
void solver::assemble_δ0() {


	gauss_quadrature_1D const quad(quadrature_order);
	unsigned const number_of_quadrature_points = quad.number_of_points;

	//Matrix<real> normals(2, 3);
	//Vector<real> parametrization(2);

	Eigen::MatrixXd normals(2, 3);
	Eigen::VectorXd parametrization(2);

	evaluate_edge_normal(normals);




	δ.setZero();

	for (unsigned k = 0; k < nk; k++) {


		t_pointer const K = mesh->get_triangle(k);
		unsigned const k_index = K->index;


		//CoeffMatrixOfElement<3, 3, 8> delta_temp(δ, k_index);


		for (unsigned El = 0; El < 3; El++) {


			real const a = (real) 0.0;
			real const b = (real)(El != 0) ? 1.0 : sqrt(2.0);

			real const c = (real)(b - a) / 2.0;
			real const d = (real)(b + a) / 2.0;

			//Vector<real> const normal = normals.getColumn(El);
			Eigen::VectorXd const referenceNormal = normals.col(El);


			for (unsigned n = 0; n < number_of_quadrature_points; n++) {


				real const x = (real)quad.points[n] * c + d;
				real const w = (real)quad.weigths[n] * c;


				evaluate_edge_parametrization(x, El, parametrization);
				//evaluate_edge_parametrization_derivative(x, El, parametrizationDerivative);


				real const s = parametrization(0);
				real const t = parametrization(1);
				//real const drNorm = parametrizationDerivative.norm();
				//real const drNorm = 1.0;


				//real const tC = upwindConcentration(s, t, K, basisRaviartThomas, basisPolynomial, referenceNormal, denominator, El, x);
				real const tC = upwindConcentration(s, t, K, referenceNormal, El, n);
				//real const tC = upwindConcentration2(K, referenceNormal, El, n);
				//real const tC = upwindConcentration3(K, referenceNormal, El, n);


				for (unsigned m = 0; m < 3; m++)
					for (unsigned j = 0; j < 8; j++)
						δ.setCoeff(k_index, m, El, j) = δ(k_index, m, El, j) + tC * raviartThomasBasisDotNormalTimesPolynomialBasis(n, j, El, m);
					  //delta_temp.setCoeff(m, El, j) = delta_temp(m, El, j) + tC * raviartThomasBasisDotNormalTimesPolynomialBasis(n, j, El, m);

			}
		}

		//for (unsigned El = 0; El < 3; El++)
		//	for (unsigned j = 0; j < 8; j++)
		//		for (unsigned m = 0; m < 3; m++)
		//			δ.setCoeff(k_index, m, El, j) = abs(δ(k_index, m, El, j)) < INTEGRAL_PRECISION ? 0.0 : δ(k_index, m, El, j);

	}

};
void solver::assemble_δ1() {

	//raviartThomasBasis_quadPoints.setNumberOfElements(quadrature_order)
	//polynomialBasis_quadPoints.setNumberOfElements(quadrature_order)

	// Quadrature weights and points on reference segment [-1,1]
	gauss_quadrature_1D const quad(quadrature_order);
	//unsigned const number_of_quadrature_points = quad.number_of_points;
	//double * const quad_points = quadrature_points_1D;
	//double * const quad_weights = quadrature_weights_1D;
	unsigned const number_of_quadrature_points = number_of_quad_points;


	//double * const quad_points_e0 = quadrature_points_1D_e0;
	//double * const quad_weights_e0 = quadrature_weights_1D_e0;

	//double * const quad_points_e12 = quadrature_points_1D_e12;
	//double * const quad_weights_e12 = quadrature_weights_1D_e12;


	//Matrix<real> normals(2, 3);
	//Vector<real> parametrization(2);

	Eigen::MatrixXd normals(2, 3);
	//Eigen::VectorXd parametrization(2);

	evaluate_edge_normal(normals);


	δ.setZero();

	for (unsigned k = 0; k < nk; k++) {


		t_pointer const K = mesh->get_triangle(k);
		unsigned const k_index = K->index;


		for (unsigned El = 0; El < 3; El++) {


			//real const a = (real) 0.0;
			//real const b = (real)(El != 0) ? 1.0 : sqrt(2.0);

			//real const c = (real)(b - a) / 2.0;
			//real const d = (real)(b + a) / 2.0;

			//Vector<real> const normal = normals.getColumn(El);
			Eigen::VectorXd const referenceNormal = normals.col(El);

			//double * quad_points;
			//double * quad_weights;

			//if (El == 0) {
			//	quad_points = quadrature_points_1D_e0;
			//	quad_weights = quadrature_weights_1D_e0;
			//}
			//else {
			//	quad_points = quadrature_points_1D_e12;
			//	quad_weights = quadrature_weights_1D_e12;
			//}


			for (unsigned n = 0; n < number_of_quadrature_points; n++) {


				//real const x = (real)quad.points[n] * c + d;
				//real const w = (real)quad.weigths[n] * c;

				//real const x = (real)quad_points[n] * c + d;
				//real const w = (real)quad_weights[n] * c;

				//real const x = quad_points[n];
				//real const w = quad_weights[n];

				//real x;
				//real w;

				//if (El == 0) {
				//	x = quadrature_points_1D_e0[n];
				//	w = quadrature_weights_1D_e0[n];
				//}
				//else {
				//	x = quadrature_points_1D_e12[n];
				//	w = quadrature_weights_1D_e12[n];
				//}


				//evaluate_edge_parametrization(x, El, parametrization);
				//evaluate_edge_parametrization_derivative(x, El, parametrizationDerivative);

				//Vector<real> const r = parametrization;
				//Vector<real> const dr = parametrizationDerivative;

				//Eigen::VectorXd const r = parametrization;
				//Eigen::VectorXd const dr = parametrizationDerivative;

				//real const s = r(0);
				//real const t = r(1);
				//real const drNorm = dr.norm();
				//real const drNorm = 1.0;

				//real const tC = upwindConcentration(s, t, K, basisRaviartThomas, basisPolynomial, referenceNormal, denominator, El, x);
				//real const tC = upwindConcentration(s, t, K, referenceNormal, El, n);
				//real const tC = upwindConcentration2(K, referenceNormal, El, n);
				real const tC = upwindConcentration3(K, referenceNormal, El, n);


				for (unsigned m = 0; m < 3; m++)
					for (unsigned j = 0; j < 8; j++)
						δ.setCoeff(k_index, m, El, j) = δ(k_index, m, El, j) + tC * raviartThomasBasisDotNormalTimesPolynomialBasis(n, j, El, m);

			}
		}

		//for (unsigned El = 0; El < 3; El++)
		//	for (unsigned j = 0; j < 8; j++)
		//		for (unsigned m = 0; m < 3; m++)
		//			δ.setCoeff(k_index, m, El, j) = abs(δ(k_index, m, El, j)) < INTEGRAL_PRECISION ? 0.0 : δ(k_index, m, El, j);

	}

};
void solver::assemble_δ2() {


	// Quadrature weights and points on reference segment [-1,1]
	gauss_quadrature_1D const quad(quadrature_order);
	unsigned const number_of_quadrature_points = quad.number_of_points;

	//Matrix<real> normals(2, 3);
	//Vector<real> parametrization(2);
	////Vector<real> parametrizationDerivative(2);
	//Matrix<real> basisRaviartThomas(2, 8);
	//Vector<real> basisPolynomial(3);


	Eigen::MatrixXd normals(2, 3);
	Eigen::VectorXd parametrization(2);
	//Eigen::VectorXd parametrizationDerivative(2);
	Eigen::MatrixXd basisRaviartThomas(2, 8);
	Eigen::VectorXd basisPolynomial(3);

	evaluate_edge_normal(normals);


	δ.setZero();

	for (unsigned k = 0; k < nk; k++) {


		t_pointer const K = mesh->get_triangle(k);
		unsigned const k_index = K->index;


		Eigen::MatrixXd JF(2, 2);

		v_pointer const a = K->vertices[0];
		v_pointer const b = K->vertices[1];
		v_pointer const c = K->vertices[2];

		real const x0 = (real)a->x;
		real const y0 = (real)a->y;

		real const x1 = (real)b->x;
		real const y1 = (real)b->y;

		real const x2 = (real)c->x;
		real const y2 = (real)c->y;

		real const detJF = (real) abs((x1 - x0)*(y2 - y0) - (x2 - x0)*(y1 - y0));

		JF(0, 0) = x1 - x0;
		JF(0, 1) = x2 - x0;
		JF(1, 0) = y1 - y0;
		JF(1, 1) = y2 - y0;

		Eigen::MatrixXd const itJF = (JF.inverse()).transpose();



		for (unsigned El = 0; El < 3; El++) {


			real const a = (real) 0.0;
			real const b = (real)(El != 0) ? 1.0 : sqrt(2.0);

			real const c = (real)(b - a) / 2.0;
			real const d = (real)(b + a) / 2.0;

			//Vector<real> const normal = normals.getColumn(El);
			Eigen::VectorXd const referenceNormal = normals.col(El);

			real const denominator = detJF * (itJF * referenceNormal).norm();

			e_pointer const E = K->edges[El];


			for (unsigned n = 0; n < number_of_quadrature_points; n++) {


				real const x = (real)quad.points[n] * c + d;
				real const w = (real)quad.weigths[n] * c;


				evaluate_edge_parametrization(x, El, parametrization);
				//evaluate_edge_parametrization_derivative(x, El, parametrizationDerivative);

				//Vector<real> const r = parametrization;
				//Vector<real> const dr = parametrizationDerivative;

				Eigen::VectorXd const r = parametrization;
				//Eigen::VectorXd const dr = parametrizationDerivative;

				real const s = r(0);
				real const t = r(1);
				//real const drNorm = dr.norm();
				real const drNorm = 1.0;


				evaluate_raviartthomas_basis(s, t, basisRaviartThomas);
				evaluate_polynomial_basis(s, t, basisPolynomial);

				real const tC = upwindConcentration(s, t, K, basisRaviartThomas, basisPolynomial, referenceNormal, denominator, El, quad.points[n]);


				for (unsigned m = 0; m < 3; m++) {


					real const Phim = basisPolynomial(m);

					for (unsigned j = 0; j < 8; j++) {

						//Vector<real> const Wj = basisRaviartThomas.getColumn(j);
						Eigen::VectorXd const Wj = basisRaviartThomas.col(j);

						//real const dotProduct = dot(Wj, normal);
						real const dotProduct = Wj.dot(referenceNormal);

						δ.setCoeff(k_index, m, El, j) = δ(k_index, m, El, j) + w * tC * dotProduct * Phim * drNorm;

					}
				}
			}
		}

		for (unsigned El = 0; El < 3; El++)
			for (unsigned j = 0; j < 8; j++)
				for (unsigned m = 0; m < 3; m++)
					δ.setCoeff(k_index, m, El, j) = abs(δ(k_index, m, El, j)) < INTEGRAL_PRECISION ? 0.0 : δ(k_index, m, El, j);

	}

};
void solver::assemble_δ3() {


	// Quadrature weights and points on reference segment [-1,1]
	gauss_quadrature_1D const quad(quadrature_order);
	unsigned const number_of_quadrature_points = quad.number_of_points;

	//Matrix<real> normals(2, 3);
	//Vector<real> parametrization(2);
	////Vector<real> parametrizationDerivative(2);
	//Matrix<real> basisRaviartThomas(2, 8);
	//Vector<real> basisPolynomial(3);


	Eigen::MatrixXd normals(2, 3);
	Eigen::VectorXd parametrization(2);
	//Eigen::VectorXd parametrizationDerivative(2);
	Eigen::MatrixXd basisRaviartThomas(2, 8);
	Eigen::VectorXd basisPolynomial(3);

	evaluate_edge_normal(normals);


	δ.setZero();

	for (unsigned k = 0; k < nk; k++) {


		t_pointer const K = mesh->get_triangle(k);
		unsigned const k_index = K->index;


		Eigen::MatrixXd JF(2, 2);

		v_pointer const a = K->vertices[0];
		v_pointer const b = K->vertices[1];
		v_pointer const c = K->vertices[2];

		real const x0 = (real)a->x;
		real const y0 = (real)a->y;

		real const x1 = (real)b->x;
		real const y1 = (real)b->y;

		real const x2 = (real)c->x;
		real const y2 = (real)c->y;

		real const detJF = (real)abs((x1 - x0)*(y2 - y0) - (x2 - x0)*(y1 - y0));

		JF(0, 0) = x1 - x0;
		JF(0, 1) = x2 - x0;
		JF(1, 0) = y1 - y0;
		JF(1, 1) = y2 - y0;

		Eigen::MatrixXd const itJF = (JF.inverse()).transpose();



		for (unsigned El = 0; El < 3; El++) {


			real const a = (real) 0.0;
			real const b = (real)(El != 0) ? 1.0 : sqrt(2.0);

			real const c = (real)(b - a) / 2.0;
			real const d = (real)(b + a) / 2.0;

			//Vector<real> const normal = normals.getColumn(El);
			Eigen::VectorXd const referenceNormal = normals.col(El);

			real const denominator = detJF * (itJF * referenceNormal).norm();

			e_pointer const E = K->edges[El];


			for (unsigned n = 0; n < number_of_quadrature_points; n++) {


				real const x = (real)quad.points[n] * c + d;
				real const w = (real)quad.weigths[n] * c;


				evaluate_edge_parametrization(x, El, parametrization);
				//evaluate_edge_parametrization_derivative(x, El, parametrizationDerivative);

				//Vector<real> const r = parametrization;
				//Vector<real> const dr = parametrizationDerivative;

				Eigen::VectorXd const r = parametrization;
				//Eigen::VectorXd const dr = parametrizationDerivative;

				real const s = r(0);
				real const t = r(1);
				//real const drNorm = dr.norm();
				real const drNorm = 1.0;


				evaluate_raviartthomas_basis(s, t, basisRaviartThomas);
				evaluate_polynomial_basis(s, t, basisPolynomial);

				real const tC = upwindConcentration(s, t, K, basisRaviartThomas, basisPolynomial, referenceNormal, denominator, El, quad.points[n]);


				for (unsigned m = 0; m < 3; m++)
					for (unsigned j = 0; j < 8; j++)
						δ.setCoeff(k_index, m, El, j) = δ(k_index, m, El, j) + tC * raviartThomasBasisDotNormalTimesPolynomialBasis(n, j, El, m);

			}
		}

		//for (unsigned El = 0; El < 3; El++)
		//	for (unsigned j = 0; j < 8; j++)
		//		for (unsigned m = 0; m < 3; m++)
		//			δ.setCoeff(k_index, m, El, j) = abs(δ(k_index, m, El, j)) < INTEGRAL_PRECISION ? 0.0 : δ(k_index, m, El, j);

	}

};
void solver::assemble_δ4() {


	// Quadrature weights and points on reference segment [-1,1]
	gauss_quadrature_1D const quad(quadrature_order);
	unsigned const number_of_quadrature_points = quad.number_of_points;

	//Matrix<real> normals(2, 3);
	//Vector<real> parametrization(2);
	////Vector<real> parametrizationDerivative(2);
	//Matrix<real> basisRaviartThomas(2, 8);
	//Vector<real> basisPolynomial(3);


	Eigen::MatrixXd normals(2, 3);
	Eigen::VectorXd parametrization(2);
	//Eigen::VectorXd parametrizationDerivative(2);
	//Eigen::MatrixXd basisRaviartThomas(2, 8);
	//Eigen::VectorXd basisPolynomial(3);

	evaluate_edge_normal(normals);


	δ.setZero();

	for (unsigned k = 0; k < nk; k++) {


		t_pointer const K = mesh->get_triangle(k);
		unsigned const k_index = K->index;


		Eigen::MatrixXd JF(2, 2);

		v_pointer const a = K->vertices[0];
		v_pointer const b = K->vertices[1];
		v_pointer const c = K->vertices[2];

		real const x0 = (real)a->x;
		real const y0 = (real)a->y;

		real const x1 = (real)b->x;
		real const y1 = (real)b->y;

		real const x2 = (real)c->x;
		real const y2 = (real)c->y;

		real const detJF = (real)abs((x1 - x0)*(y2 - y0) - (x2 - x0)*(y1 - y0));

		JF(0, 0) = x1 - x0;
		JF(0, 1) = x2 - x0;
		JF(1, 0) = y1 - y0;
		JF(1, 1) = y2 - y0;

		Eigen::MatrixXd const itJF = (JF.inverse()).transpose();



		for (unsigned El = 0; El < 3; El++) {


			real const a = (real) 0.0;
			real const b = (real)(El != 0) ? 1.0 : sqrt(2.0);

			real const c = (real)(b - a) / 2.0;
			real const d = (real)(b + a) / 2.0;

			//Vector<real> const normal = normals.getColumn(El);
			Eigen::VectorXd const referenceNormal = normals.col(El);

			real const denominator = detJF * (itJF * referenceNormal).norm();

			e_pointer const E = K->edges[El];


			for (unsigned n = 0; n < number_of_quadrature_points; n++) {


				real const x = (real)quad.points[n] * c + d;
				real const w = (real)quad.weigths[n] * c;


				evaluate_edge_parametrization(x, El, parametrization);
				//evaluate_edge_parametrization_derivative(x, El, parametrizationDerivative);

				//Vector<real> const r = parametrization;
				//Vector<real> const dr = parametrizationDerivative;

				Eigen::VectorXd const r = parametrization;
				//Eigen::VectorXd const dr = parametrizationDerivative;

				real const s = r(0);
				real const t = r(1);
				//real const drNorm = dr.norm();
				real const drNorm = 1.0;


				//evaluate_raviartthomas_basis(s, t, basisRaviartThomas);
				//evaluate_polynomial_basis(s, t, basisPolynomial);




				//////real normalVelocity = 0.0;
				//////if (E->marker == E_MARKER::NEUMANN)
				//////	normalVelocity = NEUMANN_GAMMA_Q_velocity(E, nt*dt);
				//////else
				//////	normalVelocity = velocityInNormalDirection(K, E, basisRaviartThomas, referenceNormal, denominator);


				//real const tC = upwindConcentration(s, t, K, basisRaviartThomas, basisPolynomial, referenceNormal, denominator, El, quad.points[n]);
				real const tC = upwindConcentration(s, t, K, referenceNormal, El, n);

				//Eigen::Matrix<real, 3, 8> integral;
				//integral.setZero();

				for (unsigned m = 0; m < 3; m++) {


					//real const Phim = basisPolynomial(m);

					for (unsigned j = 0; j < 8; j++) {

						//Vector<real> const Wj = basisRaviartThomas.getColumn(j);
						//Eigen::VectorXd const Wj = basisRaviartThomas.col(j);

						//real const dotProduct = dot(Wj, normal);
						//real const dotProduct = Wj.dot(referenceNormal);


						real const value = raviartThomasBasisDotNormalTimesPolynomialBasis(n, j, El, m);
						δ.setCoeff(k_index, m, El, j) = δ(k_index, m, El, j) + tC * value;


						//integral.coeffRef(m, j) += w * tC * dotProduct * Phim * drNorm;

						//δ.setCoeff(k_index, m, El, j) = δ(k_index, m, El, j) + w * tC * dotProduct * Phim * drNorm;

					}
				}
			}
		}

		for (unsigned El = 0; El < 3; El++)
			for (unsigned j = 0; j < 8; j++)
				for (unsigned m = 0; m < 3; m++)
					δ.setCoeff(k_index, m, El, j) = abs(δ(k_index, m, El, j)) < INTEGRAL_PRECISION ? 0.0 : δ(k_index, m, El, j);

	}

};
void solver::assemble_δ5() {

	//raviartThomasBasis_quadPoints.setNumberOfElements(quadrature_order)
	//polynomialBasis_quadPoints.setNumberOfElements(quadrature_order)

	// Quadrature weights and points on reference segment [-1,1]
	gauss_quadrature_1D const quad(quadrature_order);
	//unsigned const number_of_quadrature_points = quad.number_of_points;
	//double * const quad_points = quadrature_points_1D;
	//double * const quad_weights = quadrature_weights_1D;
	unsigned const number_of_quadrature_points = number_of_quad_points;


	//double * const quad_points_e0 = quadrature_points_1D_e0;
	//double * const quad_weights_e0 = quadrature_weights_1D_e0;

	//double * const quad_points_e12 = quadrature_points_1D_e12;
	//double * const quad_weights_e12 = quadrature_weights_1D_e12;


	//Matrix<real> normals(2, 3);
	//Vector<real> parametrization(2);

	Eigen::MatrixXd normals(2, 3);
	Eigen::VectorXd parametrization(2);

	evaluate_edge_normal(normals);


	δ.setZero();

	for (unsigned k = 0; k < nk; k++) {


		t_pointer const K = mesh->get_triangle(k);
		unsigned const k_index = K->index;


		for (unsigned El = 0; El < 3; El++) {


			real const a = (real) 0.0;
			real const b = (real)(El != 0) ? 1.0 : sqrt(2.0);

			real const c = (real)(b - a) / 2.0;
			real const d = (real)(b + a) / 2.0;

			//Vector<real> const normal = normals.getColumn(El);
			Eigen::VectorXd const referenceNormal = normals.col(El);

			//double * quad_points;
			//double * quad_weights;

			//if (El == 0) {
			//	quad_points = quadrature_points_1D_e0;
			//	quad_weights = quadrature_weights_1D_e0;
			//}
			//else {
			//	quad_points = quadrature_points_1D_e12;
			//	quad_weights = quadrature_weights_1D_e12;
			//}


			for (unsigned n = 0; n < number_of_quadrature_points; n++) {


				real const x = (real)quad.points[n] * c + d;
				real const w = (real)quad.weigths[n] * c;

				//real const x = (real)quad_points[n] * c + d;
				//real const w = (real)quad_weights[n] * c;

				//real const x = quad_points[n];
				//real const w = quad_weights[n];

				//real x;
				//real w;

				//if (El == 0) {
				//	x = quadrature_points_1D_e0[n];
				//	w = quadrature_weights_1D_e0[n];
				//}
				//else {
				//	x = quadrature_points_1D_e12[n];
				//	w = quadrature_weights_1D_e12[n];
				//}


				evaluate_edge_parametrization(x, El, parametrization);
				//evaluate_edge_parametrization_derivative(x, El, parametrizationDerivative);

				//Vector<real> const r = parametrization;
				//Vector<real> const dr = parametrizationDerivative;

				Eigen::VectorXd const r = parametrization;
				//Eigen::VectorXd const dr = parametrizationDerivative;

				real const s = r(0);
				real const t = r(1);
				//real const drNorm = dr.norm();
				//real const drNorm = 1.0;

				//real const tC = upwindConcentration(s, t, K, basisRaviartThomas, basisPolynomial, referenceNormal, denominator, El, x);
				real const tC = upwindConcentration(s, t, K, referenceNormal, El, n);
				//real const tC = upwindConcentration2(K, referenceNormal, El, n);
				//real const tC = upwindConcentration3(K, referenceNormal, El, n);


				for (unsigned m = 0; m < 3; m++)
					for (unsigned j = 0; j < 8; j++)
						δ.setCoeff(k_index, m, El, j) = δ(k_index, m, El, j) + tC * raviartThomasBasisDotNormalTimesPolynomialBasis(n, j, El, m);

			}
		}

		//for (unsigned El = 0; El < 3; El++)
		//	for (unsigned j = 0; j < 8; j++)
		//		for (unsigned m = 0; m < 3; m++)
		//			δ.setCoeff(k_index, m, El, j) = abs(δ(k_index, m, El, j)) < INTEGRAL_PRECISION ? 0.0 : δ(k_index, m, El, j);

	}

};
void solver::assemble_δ5() {

	//raviartThomasBasis_quadPoints.setNumberOfElements(quadrature_order)
	//polynomialBasis_quadPoints.setNumberOfElements(quadrature_order)

	// Quadrature weights and points on reference segment [-1,1]
	gauss_quadrature_1D const quad(quadrature_order);
	//unsigned const number_of_quadrature_points = quad.number_of_points;
	//double * const quad_points = quadrature_points_1D;
	//double * const quad_weights = quadrature_weights_1D;
	unsigned const number_of_quadrature_points = number_of_quad_points;


	//double * const quad_points_e0 = quadrature_points_1D_e0;
	//double * const quad_weights_e0 = quadrature_weights_1D_e0;

	//double * const quad_points_e12 = quadrature_points_1D_e12;
	//double * const quad_weights_e12 = quadrature_weights_1D_e12;


	//Matrix<real> normals(2, 3);
	//Vector<real> parametrization(2);

	Eigen::MatrixXd normals(2, 3);
	Eigen::VectorXd parametrization(2);

	evaluate_edge_normal(normals);


	δ.setZero();

	for (unsigned k = 0; k < nk; k++) {


		t_pointer const K = mesh->get_triangle(k);
		unsigned const k_index = K->index;


		for (unsigned El = 0; El < 3; El++) {


			real const a = (real) 0.0;
			real const b = (real)(El != 0) ? 1.0 : sqrt(2.0);

			real const c = (real)(b - a) / 2.0;
			real const d = (real)(b + a) / 2.0;

			//Vector<real> const normal = normals.getColumn(El);
			Eigen::VectorXd const referenceNormal = normals.col(El);

			//double * quad_points;
			//double * quad_weights;

			//if (El == 0) {
			//	quad_points = quadrature_points_1D_e0;
			//	quad_weights = quadrature_weights_1D_e0;
			//}
			//else {
			//	quad_points = quadrature_points_1D_e12;
			//	quad_weights = quadrature_weights_1D_e12;
			//}


			for (unsigned n = 0; n < number_of_quadrature_points; n++) {


				real const x = (real)quad.points[n] * c + d;
				real const w = (real)quad.weigths[n] * c;

				//real const x = (real)quad_points[n] * c + d;
				//real const w = (real)quad_weights[n] * c;

				//real const x = quad_points[n];
				//real const w = quad_weights[n];

				//real x;
				//real w;

				//if (El == 0) {
				//	x = quadrature_points_1D_e0[n];
				//	w = quadrature_weights_1D_e0[n];
				//}
				//else {
				//	x = quadrature_points_1D_e12[n];
				//	w = quadrature_weights_1D_e12[n];
				//}


				evaluate_edge_parametrization(x, El, parametrization);
				//evaluate_edge_parametrization_derivative(x, El, parametrizationDerivative);

				//Vector<real> const r = parametrization;
				//Vector<real> const dr = parametrizationDerivative;

				Eigen::VectorXd const r = parametrization;
				//Eigen::VectorXd const dr = parametrizationDerivative;

				real const s = r(0);
				real const t = r(1);
				//real const drNorm = dr.norm();
				//real const drNorm = 1.0;

				//real const tC = upwindConcentration(s, t, K, basisRaviartThomas, basisPolynomial, referenceNormal, denominator, El, x);
				real const tC = upwindConcentration(s, t, K, referenceNormal, El, n);
				//real const tC = upwindConcentration2(K, referenceNormal, El, n);
				//real const tC = upwindConcentration3(K, referenceNormal, El, n);


				for (unsigned m = 0; m < 3; m++)
					for (unsigned j = 0; j < 8; j++)
						δ.setCoeff(k_index, m, El, j) = δ(k_index, m, El, j) + tC * raviartThomasBasisDotNormalTimesPolynomialBasis(n, j, El, m);

			}
		}

		//for (unsigned El = 0; El < 3; El++)
		//	for (unsigned j = 0; j < 8; j++)
		//		for (unsigned m = 0; m < 3; m++)
		//			δ.setCoeff(k_index, m, El, j) = abs(δ(k_index, m, El, j)) < INTEGRAL_PRECISION ? 0.0 : δ(k_index, m, El, j);

	}

};


*/


// upwind
/*

real solver::upwindConcentration(real const & s, real const & t, t_pointer const & K, Eigen::MatrixXd const & basisRaviartThomas, Eigen::VectorXd const & basisPolynomial, Eigen::VectorXd const & normal, real const & denominator, unsigned const & El, real const & ksi) {


	real const time = (nt + 1) * dt;
	unsigned const k_index = K->index;


	e_pointer const E = K->edges[El];
	E_MARKER const e_marker = E->marker;


	//real const dotVelocityNormal = velocityInNormalDirection(K, E, basisRaviartThomas, normal, denominator);
	real dotVelocityNormal = 0.0;
	real concentration = 0.0;

	if (e_marker == E_MARKER::NEUMANN) {

		dotVelocityNormal = NEUMANN_GAMMA_Q_velocity(E, time);

		if (dotVelocityNormal < 0.0)
			return DIRICHLET_GAMMA_Q_concentration(s, t, time);

	}
	else {

		for (unsigned dof = 0; dof < 2; dof++) {

			unsigned const j = LI(K, E, dof);

			dotVelocityNormal += velocities(k_index, 0, j, 0) * normal.dot(basisRaviartThomas.col(j)) / denominator;
			//dotVelocityNormal += velocities(k_index, 0, j, 0) * Normal.dot(basisRaviartThomas.col(j));

		}
	}


	/////---------------------------------------------------------------------------
	//
	//real CK = 0.0;
	//real CKn = 0.0;
	//
	//
	//if (e_marker == E_MARKER::NEUMANN) {
	//
	//	CKn = DIRICHLET_GAMMA_Q_concentration(s, t, time);
	//
	//	if (dotVelocityNormal < 0.0)
	//		return CKn;
	//
	//}
	//else if (e_marker == E_MARKER::DIRICHLET) {
	//
	//	CKn = DIRICHLET_GAMMA_P_concentration(s, t, time);
	//
	//	if (dotVelocityNormal < 0.0)
	//		return CKn;
	//
	//}
	//else {
	//
	//	// Edge is not Dirichlet nor Neumann, therefore it must have TWO neighbors => There exists both CK, CKn
	//
	//	unsigned const kn_index = K->neighbors[El]->index;
	//
	//	for (unsigned m = 0; m < 3; m++)
	//		concentration += ksi(kn_index, 0, m, 0) * basisPolynomial(m);
	//
	//	//// Total limiter
	//	//concentration = ksi(kn_index, 0, 0, 0) * basisPolynomial(0);
	//}
	/////---------------------------------------------------------------------------



	if (dotVelocityNormal >= 0.0) {

		for (unsigned m = 0; m < 3; m++)
			concentration += ξ_prev(k_index, 0, m, 0) * basisPolynomial(m);

		//// Total limiter
		//concentration = ksi(k_index, 0, 0, 0) * basisPolynomial(0);

	}
	else {


		if (E->marker == E_MARKER::DIRICHLET)
			return DIRICHLET_GAMMA_P_concentration(s, t, time);


		t_pointer const Kn = K->neighbors[El];

		unsigned const kn_index = Kn->index;
		unsigned const e_index_Kn = Kn->get_edge_index(E);


		real const a = (real) 0.0;
		real const b = (real) e_index_Kn != 0 ? 1.0 : sqrt(2.0);

		real const c = (real) (b - a) / 2.0;
		real const d = (real) (b + a) / 2.0;

		real const x = (real) ksi * c + d;


		Eigen::VectorXd parametrization(2);
		Eigen::VectorXd basisPolynomial_local(3);


		evaluate_edge_parametrization(x, e_index_Kn, parametrization);

		real const ss = parametrization(0);
		real const tt = parametrization(1);

		evaluate_polynomial_basis(ss, tt, basisPolynomial_local);


		for (unsigned m = 0; m < 3; m++)
			concentration += ξ_prev(kn_index, 0, m, 0) * basisPolynomial_local(m);

		//// Total limiter
		//concentration = ksi(kn_index, 0, 0, 0) * basisPolynomial(0);

	}

	return concentration;

};

real solver::upwindConcentration(real const s, real const t, t_pointer const & K, Vector<real> const normal, unsigned const El, CoeffMatrix<1, 3, 1> & ksi) {

	real const time = (nt + 1) * dt;
	unsigned const k_index = K->index;

	Matrix<real> basisRaviartThomas(2, 8);
	Vector<real> basisPolynomial(3);


	e_pointer const E = K->edges[El];


	evaluate_polynomial_basis(s, t, basisPolynomial);
	evaluate_raviartthomas_basis(s, t, basisRaviartThomas);


	real dotVelocityNormal = 0.0;
	real concentration = 0.0;

	if (E->marker == E_MARKER::NEUMANN) {

		dotVelocityNormal = NEUMANN_GAMMA_Q_velocity(E, time);

	}
	else {


		Matrix<real> JF(2, 2);

		v_pointer const a = K->vertices[0];
		v_pointer const b = K->vertices[1];
		v_pointer const c = K->vertices[2];

		real const x0 = (real)a->x;
		real const y0 = (real)a->y;

		real const x1 = (real)b->x;
		real const y1 = (real)b->y;

		real const x2 = (real)c->x;
		real const y2 = (real)c->y;

		real const detJF = (real)abs((x1 - x0)*(y2 - y0) - (x2 - x0)*(y1 - y0));

		JF(0, 0) = x1 - x0;
		JF(0, 1) = x2 - x0;
		JF(1, 0) = y1 - y0;
		JF(1, 1) = y2 - y0;

		Matrix<real> itJF = JF.inverse();
		itJF.transposeInPlace();
		Vector<real> origNormal = itJF * normal / (itJF * normal).norm();


		for (unsigned dof = 0; dof < 2; dof++) {


			unsigned const j = LI(K, E, dof);

			dotVelocityNormal += velocities(k_index, 0, j, 0) * dot(JF * basisRaviartThomas.getColumn(j), origNormal) / detJF;
			//dotVelocityNormal += velocities(k_index, 0, j, 0) * normal.dot(basisRaviartThomas.col(j));

		}
	}



	if (dotVelocityNormal >= 0.0) {

		for (unsigned m = 0; m < 3; m++)
			concentration += ksi(k_index, 0, m, 0) * basisPolynomial(m);

		//// Total limiter
		//concentration = ksi(k_index, 0, 0, 0) * basisPolynomial(0);

	}
	else {


		if (E->marker == E_MARKER::DIRICHLET)
			return DIRICHLET_GAMMA_P_concentration(s, t, time);

		//if (E->marker == E_MARKER::DIRICHLET)
		//	return DIRICHLET_GAMMA_P_concentration(E, time);

		unsigned const kn_index = K->neighbors[El]->index;


		for (unsigned m = 0; m < 3; m++)
			concentration += ksi(kn_index, 0, m, 0) * basisPolynomial(m);

		//// Total limiter
		//concentration = ksi(kn_index, 0, 0, 0) * basisPolynomial(0);

	}

	return concentration;

};
real solver::upwindConcentration(real const s, real const t, t_pointer const & K, Eigen::VectorXd const normal, unsigned const El, unsigned const n) {


	real const time = (nt + 1) * dt;
	unsigned const k_index = K->index;

	Eigen::VectorXd basisPolynomial(3);
	Eigen::MatrixXd basisRaviartThomas(2, 8);


	e_pointer const E = K->edges[El];
	E_MARKER const e_marker = E->marker;


	evaluate_polynomial_basis(s, t, basisPolynomial);
	evaluate_raviartthomas_basis(s, t, basisRaviartThomas);


	real dotVelocityNormal = 0.0;
	real concentration = 0.0;

	if (e_marker == E_MARKER::NEUMANN) {

		dotVelocityNormal = NEUMANN_GAMMA_Q_velocity(E, time);

		if (dotVelocityNormal < 0.0) {

			gauss_quadrature_1D const quad(quadrature_order);

			real const a = (real) 0.0;
			real const b = (real) El != 0 ? 1.0 : sqrt(2.0);

			real const c = (real)(b - a) / 2.0;
			real const d = (real)(b + a) / 2.0;

			real const ksi = (real)quad.points[n] * c + d;

			Eigen::VectorXd parametrization(2);
			evaluate_edge_parametrization(ksi, El, parametrization);

			real const ss = parametrization(0);
			real const tt = parametrization(1);

			Eigen::Vector2d const refPoint(ss, tt);

			Eigen::MatrixXd JF(2, 2);

			JF(0, 0) = affineMappingMatrix[4 * k_index + 0];
			JF(0, 1) = affineMappingMatrix[4 * k_index + 1];
			JF(1, 0) = affineMappingMatrix[4 * k_index + 2];
			JF(1, 1) = affineMappingMatrix[4 * k_index + 3];

			real const x = K->vertices[0]->x + (JF * refPoint)(0);
			real const y = K->vertices[0]->y + (JF * refPoint)(1);

			return DIRICHLET_GAMMA_Q_concentration(x, y, time);

		}

	}
	else {

		Eigen::MatrixXd JF(2, 2);

		v_pointer const a = K->vertices[0];
		v_pointer const b = K->vertices[1];
		v_pointer const c = K->vertices[2];

		real const x0 = (real)a->x;
		real const y0 = (real)a->y;

		real const x1 = (real)b->x;
		real const y1 = (real)b->y;

		real const x2 = (real)c->x;
		real const y2 = (real)c->y;

		real const detJF = (real)abs((x1 - x0)*(y2 - y0) - (x2 - x0)*(y1 - y0));

		JF(0, 0) = x1 - x0;
		JF(0, 1) = x2 - x0;
		JF(1, 0) = y1 - y0;
		JF(1, 1) = y2 - y0;

		Eigen::MatrixXd itJF = (JF.inverse()).transpose();
		Eigen::VectorXd origNormal = itJF * normal / (itJF * normal).norm();


		for (unsigned dof = 0; dof < 2; dof++) {


			unsigned const j = LI(K, E, dof);

			dotVelocityNormal += velocities(k_index, 0, j, 0) * origNormal.dot(JF * basisRaviartThomas.col(j)) / detJF;
			//dotVelocityNormal += velocities(k_index, 0, j, 0) * normal.dot(basisRaviartThomas.col(j));

		}
	}



	if (dotVelocityNormal >= 0.0) {

		for (unsigned m = 0; m < 3; m++)
			concentration += ξ_prev(k_index, 0, m, 0) * basisPolynomial(m);

		//// Total limiter
		//concentration = ksi(k_index, 0, 0, 0) * basisPolynomial(0);

		//CK = concentration;

	}
	else {


		if (E->marker == E_MARKER::DIRICHLET) {

			gauss_quadrature_1D const quad(quadrature_order);

			real const a = (real) 0.0;
			real const b = (real) El != 0 ? 1.0 : sqrt(2.0);

			real const c = (real)(b - a) / 2.0;
			real const d = (real)(b + a) / 2.0;

			real const ksi = (real)quad.points[n] * c + d;

			Eigen::VectorXd parametrization(2);
			evaluate_edge_parametrization(ksi, El, parametrization);

			real const ss = parametrization(0);
			real const tt = parametrization(1);

			Eigen::Vector2d const refPoint(ss, tt);

			Eigen::MatrixXd JF(2, 2);

			JF(0, 0) = affineMappingMatrix[4 * k_index + 0];
			JF(0, 1) = affineMappingMatrix[4 * k_index + 1];
			JF(1, 0) = affineMappingMatrix[4 * k_index + 2];
			JF(1, 1) = affineMappingMatrix[4 * k_index + 3];

			real const x = K->vertices[0]->x + (JF * refPoint)(0);
			real const y = K->vertices[0]->y + (JF * refPoint)(1);

			return DIRICHLET_GAMMA_P_concentration(x, y, time);

		}


		unsigned const kn_index = K->neighbors[El]->index;
		unsigned const e_index_Kn = K->neighbors[El]->get_edge_index(E);

		real const a = (real) 0.0;
		real const b = (real) e_index_Kn != 0 ? 1.0 : sqrt(2.0);

		real const c = (real)(b - a) / 2.0;
		real const d = (real)(b + a) / 2.0;

		gauss_quadrature_1D const quad(quadrature_order);

		real const ksi = (real)quad.points[n] * c + d;


		Eigen::VectorXd parametrization(2);
		evaluate_edge_parametrization(ksi, e_index_Kn, parametrization);

		real const ss = parametrization(0);
		real const tt = parametrization(1);

		evaluate_polynomial_basis(ss, tt, basisPolynomial);


		for (unsigned m = 0; m < 3; m++)
			concentration += ξ_prev(kn_index, 0, m, 0) * basisPolynomial(m);

		//// Total limiter
		//concentration = ksi(kn_index, 0, 0, 0) * basisPolynomial(0);

		//CKn = concentration;

	}

	return concentration;

};
real solver::upwindConcentration2(t_pointer const & K, Eigen::VectorXd const normal, unsigned const El, unsigned const n) {


	real const time = (nt + 0) * dt;
	unsigned const k_index = K->index;

	//Eigen::VectorXd basisPolynomial(3);
	//Eigen::MatrixXd basisRaviartThomas(2, 8);


	e_pointer const E = K->edges[El];
	E_MARKER const e_marker = E->marker;


	//evaluate_polynomial_basis(s, t, basisPolynomial);
	//evaluate_raviartthomas_basis(s, t, basisRaviartThomas);


	real dotVelocityNormal = 0.0;
	real concentration = 0.0;

	if (e_marker == E_MARKER::NEUMANN) {


		dotVelocityNormal = NEUMANN_GAMMA_Q_velocity(E, time);

		if (dotVelocityNormal < 0.0) {


			gauss_quadrature_1D const quad(quadrature_order);

			real const a = (real) 0.0;
			real const b = (real) El != 0 ? 1.0 : sqrt(2.0);

			real const c = (real)(b - a) / 2.0;
			real const d = (real)(b + a) / 2.0;

			real const ksi = (real) quad.points[n] * c + d;

			Eigen::VectorXd parametrization(2);
			evaluate_edge_parametrization(ksi, El, parametrization);


			real const s = parametrization(0);
			real const t = parametrization(1);

			Eigen::Vector2d const refPoint(s, t);

			Eigen::MatrixXd JF(2, 2);

			JF(0, 0) = affineMappingMatrix[4 * k_index + 0];
			JF(0, 1) = affineMappingMatrix[4 * k_index + 1];
			JF(1, 0) = affineMappingMatrix[4 * k_index + 2];
			JF(1, 1) = affineMappingMatrix[4 * k_index + 3];

			real const x = K->vertices[0]->x + (JF * refPoint)(0);
			real const y = K->vertices[0]->y + (JF * refPoint)(1);

			return DIRICHLET_GAMMA_Q_concentration(x, y, time);

		}

	}
	else {

		//Eigen::MatrixXd JF(2, 2);
		//
		//v_pointer const a = K->vertices[0];
		//v_pointer const b = K->vertices[1];
		//v_pointer const c = K->vertices[2];
		//
		//real const x0 = (real)a->x;
		//real const y0 = (real)a->y;
		//
		//real const x1 = (real)b->x;
		//real const y1 = (real)b->y;
		//
		//real const x2 = (real)c->x;
		//real const y2 = (real)c->y;
		//
		//real const detJF = (real)abs((x1 - x0)*(y2 - y0) - (x2 - x0)*(y1 - y0));
		//
		//JF(0, 0) = x1 - x0;
		//JF(0, 1) = x2 - x0;
		//JF(1, 0) = y1 - y0;
		//JF(1, 1) = y2 - y0;
		//
		//Eigen::MatrixXd itJF = (JF.inverse()).transpose();
		//Eigen::VectorXd origNormal = itJF * normal / (itJF * normal).norm();
		//Eigen::VectorXd origNormal = normal / (itJF * normal).norm() / detJF;;

		Eigen::VectorXd const origNormal = normal;
		//Eigen::VectorXd const origNormal = normal / denominators[3 * k_index + El];

		//Eigen::VectorXd Wj(2);

		for (unsigned dof = 0; dof < 2; dof++) {


			unsigned const j = LI(K, E, dof);

			//Wj.coeffRef(0) = raviartThomasBasis_quadPoints(n, El, j, 0);
			//Wj.coeffRef(1) = raviartThomasBasis_quadPoints(n, El, j, 1);

			//dotVelocityNormal += velocities(k_index, 0, j, 0) * origNormal.dot(Wj);
			dotVelocityNormal += velocities(k_index, 0, j, 0) * (origNormal(0) * raviartThomasBasis_quadPoints(n, El, j, 0) + origNormal(1) * raviartThomasBasis_quadPoints(n, El, j, 1));
			//dotVelocityNormal += velocities(k_index, 0, j, 0) * origNormal.dot(JF * basisRaviartThomas.col(j)) / detJF;
			//dotVelocityNormal += velocities(k_index, 0, j, 0) * normal.dot(basisRaviartThomas.col(j));

		}
	}


	//real CK = 0.0;
	//real CKn = 0.0;
	//
	//
	//if (e_marker == E_MARKER::NEUMANN) {
	//
	//	CKn = DIRICHLET_GAMMA_Q_concentration(s, t, time);
	//
	//	if (dotVelocityNormal < 0.0)
	//		return CKn;
	//
	//}
	//else if (e_marker == E_MARKER::DIRICHLET) {
	//
	//	CKn = DIRICHLET_GAMMA_P_concentration(s, t, time);
	//
	//	if (dotVelocityNormal < 0.0)
	//		return CKn;
	//
	//}
	//else {
	//
	//	// Edge is not Dirichlet nor Neumann, therefore it must have TWO neighbors => There exists both CK, CKn
	//
	//	unsigned const kn_index = K->neighbors[El]->index;
	//
	//	for (unsigned m = 0; m < 3; m++)
	//		concentration += ksi(kn_index, 0, m, 0) * basisPolynomial(m);
	//
	//	//// Total limiter
	//	//concentration = ksi(kn_index, 0, 0, 0) * basisPolynomial(0);
	//}


	if (dotVelocityNormal >= 0.0) {

		//for (unsigned m = 0; m < 3; m++)
		//	concentration += ξ_prev(k_index, 0, m, 0) * basisPolynomial(m);

		for (unsigned m = 0; m < 3; m++)
			concentration += ξ_prev(k_index, 0, m, 0) * polynomialBasis_quadPoints(n, 0, El, m);

		//// Total limiter
		//concentration = ksi(k_index, 0, 0, 0) * basisPolynomial(0);

		//CK = concentration;

	}
	else {


		if (E->marker == E_MARKER::DIRICHLET) {


			gauss_quadrature_1D const quad(quadrature_order);

			real const a = (real) 0.0;
			real const b = (real) El != 0 ? 1.0 : sqrt(2.0);

			real const c = (real)(b - a) / 2.0;
			real const d = (real)(b + a) / 2.0;

			real const ksi = (real)quad.points[n] * c + d;


			Eigen::VectorXd parametrization(2);
			evaluate_edge_parametrization(ksi, El, parametrization);

			real const s = parametrization(0);
			real const t = parametrization(1);

			Eigen::Vector2d const refPoint(s, t);

			Eigen::MatrixXd JF(2, 2);

			JF(0, 0) = affineMappingMatrix[4 * k_index + 0];
			JF(0, 1) = affineMappingMatrix[4 * k_index + 1];
			JF(1, 0) = affineMappingMatrix[4 * k_index + 2];
			JF(1, 1) = affineMappingMatrix[4 * k_index + 3];

			real const x = K->vertices[0]->x + (JF * refPoint)(0);
			real const y = K->vertices[0]->y + (JF * refPoint)(1);

			return DIRICHLET_GAMMA_P_concentration(x, y, time);

		}

		unsigned const kn_index = K->neighbors[El]->index;
		unsigned const e_index_Kn = K->neighbors[El]->get_edge_index(E);

		//real const a = (real) 0.0;
		//real const b = (real)(e_index_Kn != 0) ? 1.0 : sqrt(2.0);
		//
		//real const c = (real)(b - a) / 2.0;
		//real const d = (real)(b + a) / 2.0;
		//
		//gauss_quadrature_1D const quad(quadrature_order);
		//
		//real const x = (real)quad.points[n] * c + d;
		//
		//
		//Eigen::VectorXd parametrization(2);
		//evaluate_edge_parametrization(x, e_index_Kn, parametrization);
		//
		//Eigen::VectorXd const r = parametrization;
		//
		//real const ss = r(0);
		//real const tt = r(1);
		//
		//evaluate_polynomial_basis(ss, tt, basisPolynomial);


		//for (unsigned m = 0; m < 3; m++)
		//	concentration += ξ_prev(kn_index, 0, m, 0) * basisPolynomial(m);


		for (unsigned m = 0; m < 3; m++)
			concentration += ξ_prev(kn_index, 0, m, 0) * polynomialBasis_quadPoints(n, 0, e_index_Kn, m);

		//// Total limiter
		//concentration = ksi(kn_index, 0, 0, 0) * basisPolynomial(0);

		//CKn = concentration;

	}

	return concentration;

};
real solver::upwindConcentration3(t_pointer const & K, Eigen::VectorXd const normal, unsigned const El, unsigned const n) {


	real const time = (nt + 0) * dt;
	unsigned const k_index = K->index;

	e_pointer const E = K->edges[El];
	E_MARKER const e_marker = E->marker;


	real dotVelocityNormal = 0.0;
	real concentration = 0.0;

	if (e_marker == E_MARKER::NEUMANN) {

		dotVelocityNormal = NEUMANN_GAMMA_Q_velocity(E, time);

		if (dotVelocityNormal < 0.0) {


			real const ksi = El != 0 ? quadrature_points_1D_e12[n] : quadrature_points_1D_e0[n];

			Eigen::VectorXd parametrization(2);
			evaluate_edge_parametrization(ksi, El, parametrization);

			real const s = parametrization(0);
			real const t = parametrization(1);

			Eigen::Vector2d const refPoint(s, t);

			Eigen::MatrixXd JF(2, 2);

			JF(0, 0) = affineMappingMatrix[4 * k_index + 0];
			JF(0, 1) = affineMappingMatrix[4 * k_index + 1];
			JF(1, 0) = affineMappingMatrix[4 * k_index + 2];
			JF(1, 1) = affineMappingMatrix[4 * k_index + 3];

			real const x = K->vertices[0]->x + (JF * refPoint)(0);
			real const y = K->vertices[0]->y + (JF * refPoint)(1);

			return DIRICHLET_GAMMA_Q_concentration(x, y, time);

		}


	}
	else {

		Eigen::Vector2d const physicalNormal = normal;
		//Eigen::Vector2d const physicalNormal = normal / denominators[3 * k_index + El];

		for (unsigned dof = 0; dof < 2; dof++) {

			unsigned const j = LI(K, E, dof);

			dotVelocityNormal += velocities(k_index, 0, j, 0) * (physicalNormal(0) * raviartThomasBasis_quadPoints(n, El, j, 0) + physicalNormal(1) * raviartThomasBasis_quadPoints(n, El, j, 1));

		}
	}

	if (dotVelocityNormal >= 0.0) {

		for (unsigned m = 0; m < 3; m++)
			concentration += ξ_prev(k_index, 0, m, 0) * polynomialBasis_quadPoints(n, 0, El, m);

	}
	else {


		if (E->marker == E_MARKER::DIRICHLET) {

			real const ksi = El != 0 ? quadrature_points_1D_e12[n] : quadrature_points_1D_e0[n];

			Eigen::VectorXd parametrization(2);
			evaluate_edge_parametrization(ksi, El, parametrization);

			real const s = parametrization(0);
			real const t = parametrization(1);

			Eigen::Vector2d const refPoint(s, t);

			Eigen::MatrixXd JF(2, 2);

			JF(0, 0) = affineMappingMatrix[4 * k_index + 0];
			JF(0, 1) = affineMappingMatrix[4 * k_index + 1];
			JF(1, 0) = affineMappingMatrix[4 * k_index + 2];
			JF(1, 1) = affineMappingMatrix[4 * k_index + 3];

			real const x = K->vertices[0]->x + (JF * refPoint)(0);
			real const y = K->vertices[0]->y + (JF * refPoint)(1);

			return DIRICHLET_GAMMA_P_concentration(x, y, time);

		}


		unsigned const kn_index = K->neighbors[El]->index;
		unsigned const e_index_Kn = K->neighbors[El]->get_edge_index(E);


		for (unsigned m = 0; m < 3; m++)
			concentration += ξ_prev(kn_index, 0, m, 0) * polynomialBasis_quadPoints(n, 0, e_index_Kn, m);

	}

	return concentration;

};
real solver::upwindConcentration4(t_pointer const & K, unsigned const El, unsigned const n) {


	real const time = (nt + 0) * dt;

	unsigned const k_index = K->index;
	e_pointer const E = K->edges[El];
	E_MARKER const e_marker = E->marker;


	real dotVelocityNormal = 0.0;
	real concentration = 0.0;

	if (e_marker == E_MARKER::NEUMANN) {


		dotVelocityNormal = NEUMANN_GAMMA_Q_velocity(E, time);


		if (dotVelocityNormal < 0.0) {

			real const ksi = El != 0 ? quadrature_points_1D_e12[n] : quadrature_points_1D_e0[n];

			Eigen::VectorXd parametrization(2);
			evaluate_edge_parametrization(ksi, El, parametrization);

			real const s = parametrization(0);
			real const t = parametrization(1);

			Eigen::Vector2d const refPoint(s, t);

			Eigen::Matrix2d JF(2, 2);

			JF(0, 0) = affineMappingMatrix[4 * k_index + 0];
			JF(0, 1) = affineMappingMatrix[4 * k_index + 1];
			JF(1, 0) = affineMappingMatrix[4 * k_index + 2];
			JF(1, 1) = affineMappingMatrix[4 * k_index + 3];

			real const x = K->vertices[0]->x + (JF * refPoint)(0);
			real const y = K->vertices[0]->y + (JF * refPoint)(1);

			return DIRICHLET_GAMMA_Q_concentration(x, y, time);

		}
	}
	else {

		Eigen::Vector2d const physicalNormal = evaluate_edge_normal(El);
		//Eigen::Vector2d const physicalNormal = evaluate_edge_normal(El) / denominators[3 * k_index + El];

		for (unsigned dof = 0; dof < 2; dof++) {

			unsigned const j = LI(K, E, dof);

			dotVelocityNormal += velocities(k_index, 0, j, 0) * (physicalNormal(0) * raviartThomasBasis_quadPoints(n, El, j, 0) + physicalNormal(1) * raviartThomasBasis_quadPoints(n, El, j, 1));

		}
	}

	if (dotVelocityNormal >= 0.0) {

		for (unsigned m = 0; m < 3; m++)
			concentration += ξ_prev(k_index, 0, m, 0) * polynomialBasis_quadPoints(n, 0, El, m);

	}
	else {


		if (E->marker == E_MARKER::DIRICHLET) {

			real const ksi = El != 0 ? quadrature_points_1D_e12[n] : quadrature_points_1D_e0[n];

			Eigen::VectorXd parametrization(2);
			evaluate_edge_parametrization(ksi, El, parametrization);

			real const s = parametrization(0);
			real const t = parametrization(1);

			Eigen::Vector2d const refPoint(s, t);

			Eigen::Matrix2d JF(2, 2);

			JF(0, 0) = affineMappingMatrix[4 * k_index + 0];
			JF(0, 1) = affineMappingMatrix[4 * k_index + 1];
			JF(1, 0) = affineMappingMatrix[4 * k_index + 2];
			JF(1, 1) = affineMappingMatrix[4 * k_index + 3];

			real const x = K->vertices[0]->x + (JF * refPoint)(0);
			real const y = K->vertices[0]->y + (JF * refPoint)(1);

			return DIRICHLET_GAMMA_P_concentration(x, y, time);

		}


		unsigned const kn_index = K->neighbors[El]->index;
		unsigned const e_index_Kn = K->neighbors[El]->get_edge_index(E);


		for (unsigned m = 0; m < 3; m++)
			concentration += ξ_prev(kn_index, 0, m, 0) * polynomialBasis_quadPoints(n, 0, e_index_Kn, m);

	}

	return concentration;

};
real solver::upwindConcentration5(t_pointer const & K, unsigned const El, unsigned const n) {


	real const time = (nt + 0) * dt;

	unsigned const k_index = K->index;
	e_pointer const E = K->edges[El];
	E_MARKER const e_marker = E->marker;


	real dotVelocityNormal = 0.0;
	real concentration = 0.0;

	if (e_marker == E_MARKER::NEUMANN) {


		dotVelocityNormal = NEUMANN_GAMMA_Q_velocity(E, time);


		if (dotVelocityNormal < 0.0) {

			real const ksi = El != 0 ? quadrature_points_1D_e12[n] : quadrature_points_1D_e0[n];

			Eigen::VectorXd parametrization(2);
			evaluate_edge_parametrization(ksi, El, parametrization);

			real const s = parametrization(0);
			real const t = parametrization(1);

			real const JF00 = affineMappingMatrix[4 * k_index + 0];
			real const JF01 = affineMappingMatrix[4 * k_index + 1];
			real const JF10 = affineMappingMatrix[4 * k_index + 2];
			real const JF11 = affineMappingMatrix[4 * k_index + 3];

			real const x = K->vertices[0]->x + JF00 * s + JF01 * t;
			real const y = K->vertices[0]->y + JF10 * s + JF11 * t;

			return DIRICHLET_GAMMA_Q_concentration(x, y, time);

		}
	}
	else {

		Eigen::Vector2d const physicalNormal = evaluate_edge_normal(El);
		//Eigen::Vector2d const physicalNormal = evaluate_edge_normal(El) / denominators[3 * k_index + El];

		for (unsigned dof = 0; dof < 2; dof++) {

			unsigned const j = LI(K, E, dof);

			dotVelocityNormal += velocities(k_index, 0, j, 0) * (physicalNormal(0) * raviartThomasBasis_quadPoints(n, El, j, 0) + physicalNormal(1) * raviartThomasBasis_quadPoints(n, El, j, 1));

		}
	}

	if (dotVelocityNormal >= 0.0) {

		for (unsigned m = 0; m < 3; m++)
			concentration += ξ_prev(k_index, 0, m, 0) * polynomialBasis_quadPoints(n, 0, El, m);

	}
	else {


		if (E->marker == E_MARKER::DIRICHLET) {

			real const ksi = El != 0 ? quadrature_points_1D_e12[n] : quadrature_points_1D_e0[n];

			Eigen::VectorXd parametrization(2);
			evaluate_edge_parametrization(ksi, El, parametrization);

			real const s = parametrization(0);
			real const t = parametrization(1);

			real const JF00 = affineMappingMatrix[4 * k_index + 0];
			real const JF01 = affineMappingMatrix[4 * k_index + 1];
			real const JF10 = affineMappingMatrix[4 * k_index + 2];
			real const JF11 = affineMappingMatrix[4 * k_index + 3];

			real const x = K->vertices[0]->x + JF00 * s + JF01 * t;
			real const y = K->vertices[0]->y + JF10 * s + JF11 * t;

			return DIRICHLET_GAMMA_P_concentration(x, y, time);

		}


		unsigned const kn_index = K->neighbors[El]->index;
		unsigned const e_index_Kn = K->neighbors[El]->get_edge_index(E);


		for (unsigned m = 0; m < 3; m++)
			concentration += ξ_prev(kn_index, 0, m, 0) * polynomialBasis_quadPoints(n, 0, e_index_Kn, m);

	}

	return concentration;

};
real solver::upwindConcentration6(t_pointer const & K, unsigned const El, unsigned const n) {


	real const time = (nt + 0) * dt;

	unsigned const k_index = K->index;
	e_pointer const E = K->edges[El];
	E_MARKER const e_marker = E->marker;


	real dotVelocityNormal = 0.0;
	real concentration = 0.0;

	if (e_marker == E_MARKER::NEUMANN) {


		dotVelocityNormal = NEUMANN_GAMMA_Q_velocity(E, time);


		if (dotVelocityNormal < 0.0) {

			real const ksi = El != 0 ? quadrature_points_1D_e12[n] : quadrature_points_1D_e0[n];

			Eigen::VectorXd parametrization(2);
			evaluate_edge_parametrization(ksi, El, parametrization);

			real const s = parametrization(0);
			real const t = parametrization(1);

			//real const s = edgeQuadraturePoints_s[number_of_quad_points * El + n];
			//real const t = edgeQuadraturePoints_t[number_of_quad_points * El + n];

			real const JF00 = affineMappingMatrix[4 * k_index + 0];
			real const JF01 = affineMappingMatrix[4 * k_index + 1];
			real const JF10 = affineMappingMatrix[4 * k_index + 2];
			real const JF11 = affineMappingMatrix[4 * k_index + 3];

			real const x = K->vertices[0]->x + JF00 * s + JF01 * t;
			real const y = K->vertices[0]->y + JF10 * s + JF11 * t;

			return DIRICHLET_GAMMA_Q_concentration(x, y, time);

		}
	}
	else {

		Eigen::Vector2d const physicalNormal = evaluate_edge_normal(El);
		//Eigen::Vector2d const physicalNormal = evaluate_edge_normal(El) / denominators[3 * k_index + El];

		for (unsigned dof = 0; dof < 2; dof++) {

			unsigned const j = LI(K, E, dof);

			dotVelocityNormal += velocities(k_index, 0, j, 0) * (physicalNormal(0) * raviartThomasBasis_quadPoints(n, El, j, 0) + physicalNormal(1) * raviartThomasBasis_quadPoints(n, El, j, 1));

		}
	}

	if (dotVelocityNormal >= 0.0) {

		for (unsigned m = 0; m < 3; m++)
			concentration += ξ_prev(k_index, 0, m, 0) * polynomialBasis_quadPoints(n, 0, El, m);

	}
	else {


		if (E->marker == E_MARKER::DIRICHLET) {


			real const ksi = El != 0 ? quadrature_points_1D_e12[n] : quadrature_points_1D_e0[n];

			Eigen::VectorXd parametrization(2);
			evaluate_edge_parametrization(ksi, El, parametrization);

			real const s = parametrization(0);
			real const t = parametrization(1);

			//real const s = edgeQuadraturePoints_s[number_of_quad_points * El + n];
			//real const t = edgeQuadraturePoints_t[number_of_quad_points * El + n];

			real const JF00 = affineMappingMatrix[4 * k_index + 0];
			real const JF01 = affineMappingMatrix[4 * k_index + 1];
			real const JF10 = affineMappingMatrix[4 * k_index + 2];
			real const JF11 = affineMappingMatrix[4 * k_index + 3];

			real const x = K->vertices[0]->x + JF00 * s + JF01 * t;
			real const y = K->vertices[0]->y + JF10 * s + JF11 * t;

			return DIRICHLET_GAMMA_P_concentration(x, y, time);

		}


		unsigned const kn_index = K->neighbors[El]->index;
		unsigned const e_index_Kn = K->neighbors[El]->get_edge_index(E);


		for (unsigned m = 0; m < 3; m++)
			concentration += ξ_prev(kn_index, 0, m, 0) * polynomialBasis_quadPoints(n, 0, e_index_Kn, m);

	}

	return concentration;

};
real solver::upwindConcentration7(t_pointer const & K, unsigned const El, unsigned const n) {


	real const time = (nt + 0) * dt;

	unsigned const number_of_quad_points_loc = number_of_quad_points;

	unsigned const k_index = K->index;
	e_pointer const E = K->edges[El];
	E_MARKER const e_marker = E->marker;


	real dotVelocityNormal = 0.0;
	real concentration = 0.0;

	if (e_marker == E_MARKER::NEUMANN) {


		dotVelocityNormal = NEUMANN_GAMMA_Q_velocity(E, time);


		if (dotVelocityNormal < 0.0) {

			unsigned const index = n + number_of_quad_points_loc * (El + 3 * k_index);

			real const x = edgeQuadraturePoints_x[index];
			real const y = edgeQuadraturePoints_y[index];

			return DIRICHLET_GAMMA_Q_concentration(x, y, time);

		}
	}
	else {

		Eigen::Vector2d const physicalNormal = evaluate_edge_normal(El);
		//Eigen::Vector2d const physicalNormal = evaluate_edge_normal(El) / denominators[3 * k_index + El];

		for (unsigned dof = 0; dof < 2; dof++) {

			unsigned const j = LI(K, E, dof);

			dotVelocityNormal += velocities(k_index, 0, j, 0) * (physicalNormal(0) * raviartThomasBasis_quadPoints(n, El, j, 0) + physicalNormal(1) * raviartThomasBasis_quadPoints(n, El, j, 1));

		}
	}

	if (dotVelocityNormal >= 0.0) {

		for (unsigned m = 0; m < 3; m++)
			concentration += ξ_prev(k_index, 0, m, 0) * polynomialBasis_quadPoints(n, 0, El, m);

	}
	else {


		if (E->marker == E_MARKER::DIRICHLET) {

			unsigned const index = n + number_of_quad_points_loc * (El + 3 * k_index);

			real const x = edgeQuadraturePoints_x[index];
			real const y = edgeQuadraturePoints_y[index];

			return DIRICHLET_GAMMA_P_concentration(x, y, time);

		}


		unsigned const kn_index = K->neighbors[El]->index;
		unsigned const e_index_Kn = K->neighbors[El]->get_edge_index(E);


		for (unsigned m = 0; m < 3; m++)
			concentration += ξ_prev(kn_index, 0, m, 0) * polynomialBasis_quadPoints(n, 0, e_index_Kn, m);

	}

	return concentration;

};
real solver::upwindConcentration8(t_pointer const & K, unsigned const El, unsigned const n) {


	real const time = (nt + 0) * dt;

	unsigned const number_of_quad_points_loc = number_of_quad_points;

	unsigned const k_index = K->index;
	e_pointer const E = K->edges[El];
	E_MARKER const e_marker = E->marker;


	real dotVelocityNormal = 0.0;
	real concentration = 0.0;

	if (e_marker == E_MARKER::NEUMANN) {


		dotVelocityNormal = NEUMANN_GAMMA_Q_velocity(E, time);


		if (dotVelocityNormal < 0.0) {

			unsigned const index = n + number_of_quad_points_loc * (El + 3 * k_index);

			real const x = edgeQuadraturePoints_x[index];
			real const y = edgeQuadraturePoints_y[index];

			return DIRICHLET_GAMMA_Q_concentration(x, y, time);

		}
	}
	else {

		for (unsigned dof = 0; dof < 2; dof++) {

			unsigned const j = LI(K, E, dof);

			dotVelocityNormal += velocities(k_index, 0, j, 0) * physicalNormalDotPhysicalRaviartThomasBasis_quadPoints[j + 8 * (n + number_of_quad_points_loc * (El + 3 * k_index))];

		}
	}

	if (dotVelocityNormal >= 0.0) {

		for (unsigned m = 0; m < 3; m++)
			concentration += ξ_prev(k_index, 0, m, 0) * polynomialBasis_quadPoints(n, 0, El, m);

	}
	else {


		if (E->marker == E_MARKER::DIRICHLET) {

			unsigned const index = n + number_of_quad_points_loc * (El + 3 * k_index);

			real const x = edgeQuadraturePoints_x[index];
			real const y = edgeQuadraturePoints_y[index];

			return DIRICHLET_GAMMA_P_concentration(x, y, time);

		}


		unsigned const kn_index = K->neighbors[El]->index;
		unsigned const e_index_Kn = K->neighbors[El]->get_edge_index(E);


		for (unsigned m = 0; m < 3; m++)
			concentration += ξ_prev(kn_index, 0, m, 0) * polynomialBasis_quadPoints(n, 0, e_index_Kn, m);

	}

	return concentration;

};

*/


// constructor / destructor
/*

solver::solver(Mesh & m, unsigned nt_0, double dt_0) : nk(m.get_number_of_triangles()), ne(m.get_number_of_edges()) {


	mesh = &m;

	nt = nt_0;
	dt = dt_0;




viscosities = new double[nk];
porosities = new double[nk];




betas = new double[nk];
betas_prev = new double[nk];




// Internal Pressures
π.setNumberOfElements(nk);

// Trace pressures
tπ.setNumberOfElements(nk);

// Concentrations
ξ.setNumberOfElements(nk);

// Velocities
velocities.setNumberOfElements(nk);

// Sink / sources
sources.setNumberOfElements(nk);


π.setZero();
tπ.setZero();
ξ.setZero();
velocities.setZero();
sources.setZero();

π_n.setZero();
π_prev.setZero();
ξ_n.setZero();
rkFp.setZero();
rkFp_n.setZero();
rkFc.setZero();
rkFc_n.setZero();




π_n.setNumberOfElements(nk);
π_prev.setNumberOfElements(nk);

ξ_n.setNumberOfElements(nk);
ξ_prev.setNumberOfElements(nk);

rkFp.setNumberOfElements(nk);
rkFp_n.setNumberOfElements(nk);

rkFc.setNumberOfElements(nk);
rkFc_n.setNumberOfElements(nk);




α.setNumberOfElements(nk);
β.setNumberOfElements(nk);
χ.setNumberOfElements(nk);
η.setNumberOfElements(nk);
τ.setNumberOfElements(nk);
γ.setNumberOfElements(nk);
δ.setNumberOfElements(nk);

σ.setNumberOfElements(nk);
λ.setNumberOfElements(nk);
φ.setNumberOfElements(nk);


// Precomputed values of basis function in quadrature points
raviartThomasBasis_quadPoints.setNumberOfElements(quadrature_order);
polynomialBasis_quadPoints.setNumberOfElements(quadrature_order);

raviartThomasBasisDotNormalTimesPolynomialBasis.setNumberOfElements(quadrature_order);


α.setZero();
β.setZero();
χ.setZero();
η.setZero();
τ.setZero();
γ.setZero();
δ.setZero();

σ.setZero();
λ.setZero();
φ.setZero();


// Precomputed values of basis function in quadrature points
raviartThomasBasis_quadPoints.setZero();
polynomialBasis_quadPoints.setZero();
raviartThomasBasisDotNormalTimesPolynomialBasis.setZero();


// Jacovian matrix of the affine mapping F : reference triangle -> phzsical triangle
affineMappingMatrixJF.setNumberOfElements(nk);

affineMappingMatrixDeterminant = new double[nk];



// ================================




affineMappingMatrixJF.setZero();
//edgeOrientation.setZero();
//edgeOrientation2.setZero();

System.resize(ne, ne);
Rhs.resize(ne);






internalPressureSystem.resize(2 * ne, 2 * ne);
pressureSystemRhs.resize(2 * ne);

// Inverse matrix of the matrix for Internal Pressures
iD.resize(3 * nk, 3 * nk);

// Matrices for Trace Pressures: H1 = mean pressure (dof 1), H2 = linear pressure (dof 2)
H1.resize(3 * nk, ne);
H2.resize(3 * nk, ne);

// Right-hand side of the Pressure equation
G.resize(3 * nk);

// Resulting Internal Pressures: 1. Mean Pressure, 2. Linear Pressure in x-direction, 3. Linear Pressure in y-direction
π_eigen.resize(3 * nk);




// 'tracePressureSystem_LU' is needed only locally. Its LU decomposition is stored in the 'sparseLUsolver_TracePressureSystem'
smat tracePressureSystem_LU;

traceSystem.resize(ne, ne);
traceRHS1.resize(ne, 3 * nk);
traceRHS2.resize(ne, ne);

tracePressureSystem_LU.resize(2 * ne, 2 * ne);
traceSystemRhs.resize(2 * ne);

// Matrices for Internal Pressures
R1.resize(ne, 3 * nk);
R2.resize(ne, 3 * nk);

// Matrices for Trace Pressures
M_j1_s1.resize(ne, ne);
M_j1_s2.resize(ne, ne);
M_j2_s1.resize(ne, ne);
M_j2_s2.resize(ne, ne);

// Right-hand side of the system
V1.resize(ne);
V2.resize(ne);

// Resulting Trace Pressures: Tp1 = mean pressure (dof 1), Tp2 = linear pressure (dof 2)
Tp1.resize(ne);
Tp2.resize(ne);



initializeValues();



assemble_α();
assemble_β();
assemble_χ();
assemble_η();
assemble_τ();



assembleR();
assembleM();


traceRHS2 = M_j2_s1 * M_j1_s1.cwiseInverse();
traceRHS1 = R2 - traceRHS2 * R1;
Eigen::MatrixXd const traceSystem = M_j2_s2 - traceRHS2 * M_j1_s2;

denseLUsolver_tracePressureSystem.compute(traceSystem);
sparseLUsolver_M11.compute(M_j1_s1);




std::vector<Eigen::Triplet<real>> triplet;

for (unsigned i = 0; i < ne; i++) {
	for (unsigned j = 0; j < ne; j++) {

		real const M11 = M_j1_s1.coeff(i, j);
		real const M12 = M_j1_s2.coeff(i, j);
		real const M21 = M_j2_s1.coeff(i, j);
		real const M22 = M_j2_s2.coeff(i, j);

		Eigen::Triplet<real> const T1(i, j, M11);
		Eigen::Triplet<real> const T2(i, j + ne, M12);
		Eigen::Triplet<real> const T3(i + ne, j, M21);
		Eigen::Triplet<real> const T4(i + ne, j + ne, M22);

		triplet.push_back(T1);
		triplet.push_back(T2);
		triplet.push_back(T3);
		triplet.push_back(T4);

	}
}

// Assembly sparse matrix
tracePressureSystem_LU.setFromTriplets(triplet.begin(), triplet.end());

// Compute LU decomposition and store it in the 'sparseLUsolver_TracePressureSystem'
sparseLUsolver_TracePressureSystem.compute(tracePressureSystem_LU);



//std::cout << M_j1_s1.toDense() << std::endl << std::endl;
//std::cout << M_j1_s2.toDense() << std::endl << std::endl;
//std::cout << M_j2_s1.toDense() << std::endl << std::endl;
//std::cout << M_j2_s2.toDense() << std::endl << std::endl;
//std::cout << R1.toDense() << std::endl;
//std::cout << R2.toDense() << std::endl;
//std::cout << inverseTraceSystem << std::endl << std::endl;





Eigen::MatrixXd normals(2, 3);
Eigen::VectorXd parametrization(2);
//Eigen::VectorXd parametrizationDerivative(2);
Eigen::MatrixXd basisRaviartThomas(2, 8);
Eigen::VectorXd basisPolynomial(3);


gauss_quadrature_1D const quadrature(quadrature_order);
unsigned const number_of_quadrature_points = quadrature.number_of_points;

number_of_quad_points = number_of_quadrature_points;
quadrature_points_1D_e0 = new double[number_of_quad_points];
quadrature_weights_1D_e0 = new double[number_of_quad_points];
quadrature_points_1D_e12 = new double[number_of_quad_points];
quadrature_weights_1D_e12 = new double[number_of_quad_points];


evaluate_edge_normal(normals);

//denominators = new double[3 * nk];
affineMappingMatrix = new double[4 * nk];

//edgeQuadraturePoints_s = new double[3 * number_of_quadrature_points];
//edgeQuadraturePoints_t = new double[3 * number_of_quadrature_points];

edgeQuadraturePoints_x = new double[nk * 3 * number_of_quadrature_points];
edgeQuadraturePoints_y = new double[nk * 3 * number_of_quadrature_points];

physicalNormalDotPhysicalRaviartThomasBasis_quadPoints = new double[8 * nk * 3 * number_of_quadrature_points];

//indecesOfDirichletEdges = new unsigned[mesh->get_number_of_dirichlet_edges()];
//edgeQuadraturePointsDirichlet_x = new double[mesh->get_number_of_dirichlet_edges()];
//edgeQuadraturePointsDirichlet_y = new double[mesh->get_number_of_dirichlet_edges()];

//unsigned count = 0;
//for (unsigned e = 0; e < ne; e++) {

//	if (mesh->get_edge(e)->marker == E_MARKER::DIRICHLET)
//		indecesOfDirichletEdges[count++] = mesh->get_edge(e)->index;

//}


//indeces_element = new unsigned[nk];
//indeces_edge = new unsigned[ne];

//for (unsigned k = 0; k < nk; k++)
//	indeces_element[k] = mesh->get_triangle(k)->index;

//for (unsigned e = 0; e < ne; e++)
//	indeces_edge[e] = mesh->get_edge(e)->index;


for (unsigned k = 0; k < nk; k++) {


	t_pointer const K = mesh->get_triangle(k);
	unsigned const k_index = K->index;

	Eigen::MatrixXd JF(2, 2);

	v_pointer const a = K->vertices[0];
	v_pointer const b = K->vertices[1];
	v_pointer const c = K->vertices[2];

	real const x0 = (real)a->x;
	real const y0 = (real)a->y;

	real const x1 = (real)b->x;
	real const y1 = (real)b->y;

	real const x2 = (real)c->x;
	real const y2 = (real)c->y;

	real const detJF = (real)abs((x1 - x0)*(y2 - y0) - (x2 - x0)*(y1 - y0));

	JF(0, 0) = x1 - x0;
	JF(0, 1) = x2 - x0;
	JF(1, 0) = y1 - y0;
	JF(1, 1) = y2 - y0;

	affineMappingMatrix[4 * k_index + 0] = JF(0, 0);
	affineMappingMatrix[4 * k_index + 1] = JF(0, 1);
	affineMappingMatrix[4 * k_index + 2] = JF(1, 0);
	affineMappingMatrix[4 * k_index + 3] = JF(1, 1);

	Eigen::MatrixXd itJF = (JF.inverse()).transpose();
	//Eigen::VectorXd origNormal = itJF * normal / (itJF * normal).norm();

	//for (unsigned El = 0; El < 3; El++)
	//	denominators[3 * k_index + El] = (itJF * normals.col(El)).norm() * detJF;

}

for (unsigned k = 0; k < nk; k++) {

	t_pointer const K = mesh->get_triangle(k);
	unsigned const k_index = K->index;


	Eigen::MatrixXd JF(2, 2);

	v_pointer const a = K->vertices[0];
	v_pointer const b = K->vertices[1];
	v_pointer const c = K->vertices[2];

	real const x0 = (real)a->x;
	real const y0 = (real)a->y;

	real const x1 = (real)b->x;
	real const y1 = (real)b->y;

	real const x2 = (real)c->x;
	real const y2 = (real)c->y;

	real const detJF = (real)abs((x1 - x0)*(y2 - y0) - (x2 - x0)*(y1 - y0));

	JF(0, 0) = x1 - x0;
	JF(0, 1) = x2 - x0;
	JF(1, 0) = y1 - y0;
	JF(1, 1) = y2 - y0;

	affineMappingMatrix[4 * k_index + 0] = JF(0, 0);
	affineMappingMatrix[4 * k_index + 1] = JF(0, 1);
	affineMappingMatrix[4 * k_index + 2] = JF(1, 0);
	affineMappingMatrix[4 * k_index + 3] = JF(1, 1);

	Eigen::MatrixXd itJF = (JF.inverse()).transpose();
	//Eigen::VectorXd origNormal = itJF * normal / (itJF * normal).norm();



	for (unsigned El = 0; El < 3; El++) {


		e_pointer const E = K->edges[El];

		real const a = (real) 0.0;
		real const b = (real)El != 0 ? 1.0 : sqrt(2.0);

		real const c = (real)(b - a) / 2.0;
		real const d = (real)(b + a) / 2.0;

		//Vector<real> const normal = normals.getColumn(El);
		Eigen::VectorXd const referenceNormal = normals.col(El);

		real const denominator = (itJF * normals.col(El)).norm() * detJF;



		for (unsigned n = 0; n < number_of_quadrature_points; n++) {


			real const x = (real)quadrature.points[n] * c + d;
			real const w = (real)quadrature.weigths[n] * c;

			evaluate_edge_parametrization(x, El, parametrization);
			//evaluate_edge_parametrization_derivative(x, El, parametrizationDerivative);

			//Vector<real> const r = parametrization;
			//Vector<real> const dr = parametrizationDerivative;

			Eigen::VectorXd const r = parametrization;
			//Eigen::VectorXd const dr = parametrizationDerivative;

			real const s = r(0);
			real const t = r(1);
			//real const drNorm = dr.norm();
			real const drNorm = 1.0;


			evaluate_raviartthomas_basis(s, t, basisRaviartThomas);


			
			Eigen::Vector2d const physicalNormal = referenceNormal / denominator;

			for (unsigned j = 0; j < 8; j++) {

				real const dotNormal = physicalNormal.dot(basisRaviartThomas.col(j));

				physicalNormalDotPhysicalRaviartThomasBasis_quadPoints[j + 8 * (n + number_of_quadrature_points * (El + 3 * k_index))] = dotNormal;

			}

			


		}
	}

}



for (unsigned El = 0; El < 3; El++) {


	real const a = (real) 0.0;
	real const b = (real)El != 0 ? 1.0 : sqrt(2.0);

	real const c = (real)(b - a) / 2.0;
	real const d = (real)(b + a) / 2.0;

	//Vector<real> const normal = normals.getColumn(El);
	Eigen::VectorXd const referenceNormal = normals.col(El);


	for (unsigned n = 0; n < number_of_quadrature_points; n++) {


		real const x = (real)quadrature.points[n] * c + d;
		real const w = (real)quadrature.weigths[n] * c;

		if (El == 0) {
			quadrature_points_1D_e0[n] = x;
			quadrature_weights_1D_e0[n] = w;
		}
		else {
			quadrature_points_1D_e12[n] = x;
			quadrature_weights_1D_e12[n] = w;
		}



		evaluate_edge_parametrization(x, El, parametrization);
		//evaluate_edge_parametrization_derivative(x, El, parametrizationDerivative);

		//Vector<real> const r = parametrization;
		//Vector<real> const dr = parametrizationDerivative;

		Eigen::VectorXd const r = parametrization;
		//Eigen::VectorXd const dr = parametrizationDerivative;

		real const s = r(0);
		real const t = r(1);
		//real const drNorm = dr.norm();
		real const drNorm = 1.0;


		evaluate_raviartthomas_basis(s, t, basisRaviartThomas);
		evaluate_polynomial_basis(s, t, basisPolynomial);


		for (unsigned m = 0; m < 3; m++) {


			real const Phim = basisPolynomial(m);

			polynomialBasis_quadPoints.setCoeff(n, 0, El, m) = Phim;

			for (unsigned j = 0; j < 8; j++) {

				//Vector<real> const Wj = basisRaviartThomas.getColumn(j);
				Eigen::VectorXd const Wj = basisRaviartThomas.col(j);

				//real const dotProduct = dot(Wj, normal);
				real const dotProduct = Wj.dot(referenceNormal);


				real const value = w * dotProduct * Phim * drNorm;

				raviartThomasBasis_quadPoints.setCoeff(n, El, j, 0) = Wj(0);
				raviartThomasBasis_quadPoints.setCoeff(n, El, j, 1) = Wj(1);

				raviartThomasBasisDotNormalTimesPolynomialBasis.setCoeff(n, j, El, m) = abs(value) < INTEGRAL_PRECISION ? 0.0 : value;

			}
		}
	}
}

//for (unsigned El = 0; El < 3; El++) {
//	for (unsigned n = 0; n < number_of_quadrature_points; n++) {

//		real const x = El != 0 ? quadrature_points_1D_e12[n] : quadrature_points_1D_e0[n];

//		Eigen::VectorXd parametrization(2);
//		evaluate_edge_parametrization(x, El, parametrization);

//		real const s = parametrization(0);
//		real const t = parametrization(1);

//		edgeQuadraturePoints_s[number_of_quadrature_points * El + n] = s;
//		edgeQuadraturePoints_t[number_of_quadrature_points * El + n] = t;

//	}
//}
for (unsigned k = 0; k < nk; k++) {

	t_pointer const K = mesh->get_triangle(k);
	unsigned const k_index = K->index;

	for (unsigned El = 0; El < 3; El++) {
		for (unsigned n = 0; n < number_of_quadrature_points; n++) {


			real const ksi = El != 0 ? quadrature_points_1D_e12[n] : quadrature_points_1D_e0[n];

			Eigen::VectorXd parametrization(2);
			evaluate_edge_parametrization(ksi, El, parametrization);

			real const s = parametrization(0);
			real const t = parametrization(1);

			//real const s = edgeQuadraturePoints_s[number_of_quadrature_points * El + n];
			//real const t = edgeQuadraturePoints_t[number_of_quadrature_points * El + n];

			real const JF00 = affineMappingMatrix[4 * k_index + 0];
			real const JF01 = affineMappingMatrix[4 * k_index + 1];
			real const JF10 = affineMappingMatrix[4 * k_index + 2];
			real const JF11 = affineMappingMatrix[4 * k_index + 3];

			real const x = K->vertices[0]->x + JF00 * s + JF01 * t;
			real const y = K->vertices[0]->y + JF10 * s + JF11 * t;

			edgeQuadraturePoints_x[n + number_of_quadrature_points * (El + 3 * k_index)] = x;
			edgeQuadraturePoints_y[n + number_of_quadrature_points * (El + 3 * k_index)] = y;

		}
	}
}



};
solver::~solver() {


	delete[] quadrature_points_1D_e0;
	delete[] quadrature_weights_1D_e0;
	delete[] quadrature_points_1D_e12;
	delete[] quadrature_weights_1D_e12;

	//delete[] denominators;
	delete[] affineMappingMatrix;

	//delete[] edgeQuadraturePoints_s;
	//delete[] edgeQuadraturePoints_t;

	delete[] edgeQuadraturePoints_x;
	delete[] edgeQuadraturePoints_y;

	//delete[] edgeQuadraturePointsDirichlet_x;
	//delete[] edgeQuadraturePointsDirichlet_y;

	//delete[] indecesOfDirichletEdges;

	//delete[] indeces_element;
	//delete[] indeces_edge;

	delete[] physicalNormalDotPhysicalRaviartThomasBasis_quadPoints;




	delete[] affineMappingMatrixDeterminant;

	delete[] viscosities;
	delete[] porosities;

	delete[] betas;
	delete[] betas_prev;

};


*/

/*


solver::solver(Mesh & m, unsigned nt_0, double dt_0) : nk(m.get_number_of_triangles()), ne(m.get_number_of_edges()) {


	mesh = &m;

	nt = nt_0;
	dt = dt_0;



viscosities = new double[nk];
porosities = new double[nk];





betas = new double[nk];
betas_prev = new double[nk];




// Internal Pressures
π.setNumberOfElements(nk);

// Trace pressures
tπ.setNumberOfElements(nk);

// Concentrations
ξ.setNumberOfElements(nk);

// Velocities
velocities.setNumberOfElements(nk);

// Sink / sources
sources.setNumberOfElements(nk);


π.setZero();
tπ.setZero();
ξ.setZero();
velocities.setZero();
sources.setZero();

π_n.setZero();
π_prev.setZero();
ξ_n.setZero();
rkFp.setZero();
rkFp_n.setZero();
rkFc.setZero();
rkFc_n.setZero();



π_n.setNumberOfElements(nk);
π_prev.setNumberOfElements(nk);

ξ_n.setNumberOfElements(nk);
ξ_prev.setNumberOfElements(nk);

rkFp.setNumberOfElements(nk);
rkFp_n.setNumberOfElements(nk);

rkFc.setNumberOfElements(nk);
rkFc_n.setNumberOfElements(nk);




α.setNumberOfElements(nk);
β.setNumberOfElements(nk);
χ.setNumberOfElements(nk);
η.setNumberOfElements(nk);
τ.setNumberOfElements(nk);
γ.setNumberOfElements(nk);
δ.setNumberOfElements(nk);

σ.setNumberOfElements(nk);
λ.setNumberOfElements(nk);
φ.setNumberOfElements(nk);


// Precomputed values of basis function in quadrature points
raviartThomasBasis_quadPoints.setNumberOfElements(quadrature_order);
polynomialBasis_quadPoints.setNumberOfElements(quadrature_order);

raviartThomasBasisDotNormalTimesPolynomialBasis.setNumberOfElements(quadrature_order);


α.setZero();
β.setZero();
χ.setZero();
η.setZero();
τ.setZero();
γ.setZero();
δ.setZero();

σ.setZero();
λ.setZero();
φ.setZero();


// Precomputed values of basis function in quadrature points
raviartThomasBasis_quadPoints.setZero();
polynomialBasis_quadPoints.setZero();
raviartThomasBasisDotNormalTimesPolynomialBasis.setZero();


// Jacovian matrix of the affine mapping F : reference triangle -> phzsical triangle
affineMappingMatrixJF.setNumberOfElements(nk);

affineMappingMatrixDeterminant = new double[nk];



// ================================




affineMappingMatrixJF.setZero();
//edgeOrientation.setZero();
//edgeOrientation2.setZero();

System.resize(ne, ne);
Rhs.resize(ne);




internalPressureSystem.resize(2 * ne, 2 * ne);
pressureSystemRhs.resize(2 * ne);

// Inverse matrix of the matrix for Internal Pressures
iD.resize(3 * nk, 3 * nk);

// Matrices for Trace Pressures: H1 = mean pressure (dof 1), H2 = linear pressure (dof 2)
H1.resize(3 * nk, ne);
H2.resize(3 * nk, ne);

// Right-hand side of the Pressure equation
G.resize(3 * nk);

// Resulting Internal Pressures: 1. Mean Pressure, 2. Linear Pressure in x-direction, 3. Linear Pressure in y-direction
π_eigen.resize(3 * nk);



// 'tracePressureSystem_LU' is needed only locally. Its LU decomposition is stored in the 'sparseLUsolver_TracePressureSystem'
smat tracePressureSystem_LU;

traceSystem.resize(ne, ne);
traceRHS1.resize(ne, 3 * nk);
traceRHS2.resize(ne, ne);

tracePressureSystem_LU.resize(2 * ne, 2 * ne);
traceSystemRhs.resize(2 * ne);

// Matrices for Internal Pressures
R1.resize(ne, 3 * nk);
R2.resize(ne, 3 * nk);

// Matrices for Trace Pressures
M_j1_s1.resize(ne, ne);
M_j1_s2.resize(ne, ne);
M_j2_s1.resize(ne, ne);
M_j2_s2.resize(ne, ne);

// Right-hand side of the system
V1.resize(ne);
V2.resize(ne);

// Resulting Trace Pressures: Tp1 = mean pressure (dof 1), Tp2 = linear pressure (dof 2)
Tp1.resize(ne);
Tp2.resize(ne);





initializeValues();




assemble_α();
assemble_β();
assemble_χ();
assemble_η();
assemble_τ();




assembleR();
assembleM();


traceRHS2 = M_j2_s1 * M_j1_s1.cwiseInverse();
traceRHS1 = R2 - traceRHS2 * R1;
Eigen::MatrixXd const traceSystem = M_j2_s2 - traceRHS2 * M_j1_s2;

denseLUsolver_tracePressureSystem.compute(traceSystem);
sparseLUsolver_M11.compute(M_j1_s1);


std::vector<Eigen::Triplet<real>> triplet;

for (unsigned i = 0; i < ne; i++) {
	for (unsigned j = 0; j < ne; j++) {

		real const M11 = M_j1_s1.coeff(i, j);
		real const M12 = M_j1_s2.coeff(i, j);
		real const M21 = M_j2_s1.coeff(i, j);
		real const M22 = M_j2_s2.coeff(i, j);

		Eigen::Triplet<real> const T1(i, j, M11);
		Eigen::Triplet<real> const T2(i, j + ne, M12);
		Eigen::Triplet<real> const T3(i + ne, j, M21);
		Eigen::Triplet<real> const T4(i + ne, j + ne, M22);

		triplet.push_back(T1);
		triplet.push_back(T2);
		triplet.push_back(T3);
		triplet.push_back(T4);

	}
}

// Assembly sparse matrix
tracePressureSystem_LU.setFromTriplets(triplet.begin(), triplet.end());

// Compute LU decomposition and store it in the 'sparseLUsolver_TracePressureSystem'
sparseLUsolver_TracePressureSystem.compute(tracePressureSystem_LU);



//std::cout << M_j1_s1.toDense() << std::endl << std::endl;
//std::cout << M_j1_s2.toDense() << std::endl << std::endl;
//std::cout << M_j2_s1.toDense() << std::endl << std::endl;
//std::cout << M_j2_s2.toDense() << std::endl << std::endl;
//std::cout << R1.toDense() << std::endl;
//std::cout << R2.toDense() << std::endl;
//std::cout << inverseTraceSystem << std::endl << std::endl;





Eigen::MatrixXd normals(2, 3);
Eigen::VectorXd parametrization(2);
//Eigen::VectorXd parametrizationDerivative(2);
Eigen::MatrixXd basisRaviartThomas(2, 8);
Eigen::VectorXd basisPolynomial(3);


gauss_quadrature_1D const quadrature(quadrature_order);
unsigned const number_of_quadrature_points = quadrature.number_of_points;

number_of_quad_points = number_of_quadrature_points;
quadrature_points_1D_e0 = new double[number_of_quad_points];
quadrature_weights_1D_e0 = new double[number_of_quad_points];
quadrature_points_1D_e12 = new double[number_of_quad_points];
quadrature_weights_1D_e12 = new double[number_of_quad_points];


evaluate_edge_normal(normals);

//denominators = new double[3 * nk];
affineMappingMatrix = new double[4 * nk];

//edgeQuadraturePoints_s = new double[3 * number_of_quadrature_points];
//edgeQuadraturePoints_t = new double[3 * number_of_quadrature_points];

edgeQuadraturePoints_x = new double[nk * 3 * number_of_quadrature_points];
edgeQuadraturePoints_y = new double[nk * 3 * number_of_quadrature_points];

physicalNormalDotPhysicalRaviartThomasBasis_quadPoints = new double[8 * nk * 3 * number_of_quadrature_points];

//indecesOfDirichletEdges = new unsigned[mesh->get_number_of_dirichlet_edges()];
//edgeQuadraturePointsDirichlet_x = new double[mesh->get_number_of_dirichlet_edges()];
//edgeQuadraturePointsDirichlet_y = new double[mesh->get_number_of_dirichlet_edges()];

//unsigned count = 0;
//for (unsigned e = 0; e < ne; e++) {

//	if (mesh->get_edge(e)->marker == E_MARKER::DIRICHLET)
//		indecesOfDirichletEdges[count++] = mesh->get_edge(e)->index;

//}


//indeces_element = new unsigned[nk];
//indeces_edge = new unsigned[ne];

//for (unsigned k = 0; k < nk; k++)
//	indeces_element[k] = mesh->get_triangle(k)->index;

//for (unsigned e = 0; e < ne; e++)
//	indeces_edge[e] = mesh->get_edge(e)->index;


for (unsigned k = 0; k < nk; k++) {


	t_pointer const K = mesh->get_triangle(k);
	unsigned const k_index = K->index;

	Eigen::MatrixXd JF(2, 2);

	v_pointer const a = K->vertices[0];
	v_pointer const b = K->vertices[1];
	v_pointer const c = K->vertices[2];

	real const x0 = (real)a->x;
	real const y0 = (real)a->y;

	real const x1 = (real)b->x;
	real const y1 = (real)b->y;

	real const x2 = (real)c->x;
	real const y2 = (real)c->y;

	real const detJF = (real)abs((x1 - x0)*(y2 - y0) - (x2 - x0)*(y1 - y0));

	JF(0, 0) = x1 - x0;
	JF(0, 1) = x2 - x0;
	JF(1, 0) = y1 - y0;
	JF(1, 1) = y2 - y0;

	affineMappingMatrix[4 * k_index + 0] = JF(0, 0);
	affineMappingMatrix[4 * k_index + 1] = JF(0, 1);
	affineMappingMatrix[4 * k_index + 2] = JF(1, 0);
	affineMappingMatrix[4 * k_index + 3] = JF(1, 1);

	Eigen::MatrixXd itJF = (JF.inverse()).transpose();
	//Eigen::VectorXd origNormal = itJF * normal / (itJF * normal).norm();

	//for (unsigned El = 0; El < 3; El++)
	//	denominators[3 * k_index + El] = (itJF * normals.col(El)).norm() * detJF;

}

for (unsigned k = 0; k < nk; k++) {

	t_pointer const K = mesh->get_triangle(k);
	unsigned const k_index = K->index;


	Eigen::MatrixXd JF(2, 2);

	v_pointer const a = K->vertices[0];
	v_pointer const b = K->vertices[1];
	v_pointer const c = K->vertices[2];

	real const x0 = (real)a->x;
	real const y0 = (real)a->y;

	real const x1 = (real)b->x;
	real const y1 = (real)b->y;

	real const x2 = (real)c->x;
	real const y2 = (real)c->y;

	real const detJF = (real)abs((x1 - x0)*(y2 - y0) - (x2 - x0)*(y1 - y0));

	JF(0, 0) = x1 - x0;
	JF(0, 1) = x2 - x0;
	JF(1, 0) = y1 - y0;
	JF(1, 1) = y2 - y0;

	affineMappingMatrix[4 * k_index + 0] = JF(0, 0);
	affineMappingMatrix[4 * k_index + 1] = JF(0, 1);
	affineMappingMatrix[4 * k_index + 2] = JF(1, 0);
	affineMappingMatrix[4 * k_index + 3] = JF(1, 1);

	Eigen::MatrixXd itJF = (JF.inverse()).transpose();
	//Eigen::VectorXd origNormal = itJF * normal / (itJF * normal).norm();



	for (unsigned El = 0; El < 3; El++) {


		e_pointer const E = K->edges[El];

		real const a = (real) 0.0;
		real const b = (real)El != 0 ? 1.0 : sqrt(2.0);

		real const c = (real)(b - a) / 2.0;
		real const d = (real)(b + a) / 2.0;

		//Vector<real> const normal = normals.getColumn(El);
		Eigen::VectorXd const referenceNormal = normals.col(El);

		real const denominator = (itJF * normals.col(El)).norm() * detJF;



		for (unsigned n = 0; n < number_of_quadrature_points; n++) {


			real const x = (real)quadrature.points[n] * c + d;
			real const w = (real)quadrature.weigths[n] * c;

			evaluate_edge_parametrization(x, El, parametrization);
			//evaluate_edge_parametrization_derivative(x, El, parametrizationDerivative);

			//Vector<real> const r = parametrization;
			//Vector<real> const dr = parametrizationDerivative;

			Eigen::VectorXd const r = parametrization;
			//Eigen::VectorXd const dr = parametrizationDerivative;

			real const s = r(0);
			real const t = r(1);
			//real const drNorm = dr.norm();
			real const drNorm = 1.0;


			evaluate_raviartthomas_basis(s, t, basisRaviartThomas);


		
			Eigen::Vector2d const physicalNormal = referenceNormal / denominator;

			for (unsigned j = 0; j < 8; j++) {

				real const dotNormal = physicalNormal.dot(basisRaviartThomas.col(j));

				physicalNormalDotPhysicalRaviartThomasBasis_quadPoints[j + 8 * (n + number_of_quadrature_points * (El + 3 * k_index))] = dotNormal;

			}

		


		}
	}

}



for (unsigned El = 0; El < 3; El++) {


	real const a = (real) 0.0;
	real const b = (real)El != 0 ? 1.0 : sqrt(2.0);

	real const c = (real)(b - a) / 2.0;
	real const d = (real)(b + a) / 2.0;

	//Vector<real> const normal = normals.getColumn(El);
	Eigen::VectorXd const referenceNormal = normals.col(El);


	for (unsigned n = 0; n < number_of_quadrature_points; n++) {


		real const x = (real)quadrature.points[n] * c + d;
		real const w = (real)quadrature.weigths[n] * c;

		if (El == 0) {
			quadrature_points_1D_e0[n] = x;
			quadrature_weights_1D_e0[n] = w;
		}
		else {
			quadrature_points_1D_e12[n] = x;
			quadrature_weights_1D_e12[n] = w;
		}



		evaluate_edge_parametrization(x, El, parametrization);
		//evaluate_edge_parametrization_derivative(x, El, parametrizationDerivative);

		//Vector<real> const r = parametrization;
		//Vector<real> const dr = parametrizationDerivative;

		Eigen::VectorXd const r = parametrization;
		//Eigen::VectorXd const dr = parametrizationDerivative;

		real const s = r(0);
		real const t = r(1);
		//real const drNorm = dr.norm();
		real const drNorm = 1.0;


		evaluate_raviartthomas_basis(s, t, basisRaviartThomas);
		evaluate_polynomial_basis(s, t, basisPolynomial);


		for (unsigned m = 0; m < 3; m++) {


			real const Phim = basisPolynomial(m);

			polynomialBasis_quadPoints.setCoeff(n, 0, El, m) = Phim;

			for (unsigned j = 0; j < 8; j++) {

				//Vector<real> const Wj = basisRaviartThomas.getColumn(j);
				Eigen::VectorXd const Wj = basisRaviartThomas.col(j);

				//real const dotProduct = dot(Wj, normal);
				real const dotProduct = Wj.dot(referenceNormal);


				real const value = w * dotProduct * Phim * drNorm;

				raviartThomasBasis_quadPoints.setCoeff(n, El, j, 0) = Wj(0);
				raviartThomasBasis_quadPoints.setCoeff(n, El, j, 1) = Wj(1);

				raviartThomasBasisDotNormalTimesPolynomialBasis.setCoeff(n, j, El, m) = abs(value) < INTEGRAL_PRECISION ? 0.0 : value;

			}
		}
	}
}

//for (unsigned El = 0; El < 3; El++) {
//	for (unsigned n = 0; n < number_of_quadrature_points; n++) {

//		real const x = El != 0 ? quadrature_points_1D_e12[n] : quadrature_points_1D_e0[n];

//		Eigen::VectorXd parametrization(2);
//		evaluate_edge_parametrization(x, El, parametrization);

//		real const s = parametrization(0);
//		real const t = parametrization(1);

//		edgeQuadraturePoints_s[number_of_quadrature_points * El + n] = s;
//		edgeQuadraturePoints_t[number_of_quadrature_points * El + n] = t;

//	}
//}
for (unsigned k = 0; k < nk; k++) {

	t_pointer const K = mesh->get_triangle(k);
	unsigned const k_index = K->index;

	for (unsigned El = 0; El < 3; El++) {
		for (unsigned n = 0; n < number_of_quadrature_points; n++) {


			real const ksi = El != 0 ? quadrature_points_1D_e12[n] : quadrature_points_1D_e0[n];

			Eigen::VectorXd parametrization(2);
			evaluate_edge_parametrization(ksi, El, parametrization);

			real const s = parametrization(0);
			real const t = parametrization(1);

			//real const s = edgeQuadraturePoints_s[number_of_quadrature_points * El + n];
			//real const t = edgeQuadraturePoints_t[number_of_quadrature_points * El + n];

			real const JF00 = affineMappingMatrix[4 * k_index + 0];
			real const JF01 = affineMappingMatrix[4 * k_index + 1];
			real const JF10 = affineMappingMatrix[4 * k_index + 2];
			real const JF11 = affineMappingMatrix[4 * k_index + 3];

			real const x = K->vertices[0]->x + JF00 * s + JF01 * t;
			real const y = K->vertices[0]->y + JF10 * s + JF11 * t;

			edgeQuadraturePoints_x[n + number_of_quadrature_points * (El + 3 * k_index)] = x;
			edgeQuadraturePoints_y[n + number_of_quadrature_points * (El + 3 * k_index)] = y;

		}
	}
}


};


*/

/*


	//indecesOfDirichletEdges = new unsigned[mesh->get_number_of_dirichlet_edges()];
	//edgeQuadraturePointsDirichlet_x = new double[mesh->get_number_of_dirichlet_edges()];
	//edgeQuadraturePointsDirichlet_y = new double[mesh->get_number_of_dirichlet_edges()];

	//unsigned count = 0;
	//for (unsigned e = 0; e < ne; e++) {

	//	if (mesh->get_edge(e)->marker == E_MARKER::DIRICHLET)
	//		indecesOfDirichletEdges[count++] = mesh->get_edge(e)->index;

	//}


	//indeces_element = new unsigned[nk];
	//indeces_edge = new unsigned[ne];

	//for (unsigned k = 0; k < nk; k++)
	//	indeces_element[k] = mesh->get_triangle(k)->index;

	//for (unsigned e = 0; e < ne; e++)
	//	indeces_edge[e] = mesh->get_edge(e)->index;


	for (unsigned k = 0; k < nk; k++) {

		t_pointer const K = mesh->get_triangle(k);
		unsigned const k_index = K->index;


		Eigen::MatrixXd JF(2, 2);

		v_pointer const a = K->vertices[0];
		v_pointer const b = K->vertices[1];
		v_pointer const c = K->vertices[2];

		real const x0 = (real)a->x;
		real const y0 = (real)a->y;

		real const x1 = (real)b->x;
		real const y1 = (real)b->y;

		real const x2 = (real)c->x;
		real const y2 = (real)c->y;

		real const detJF = (real)abs((x1 - x0)*(y2 - y0) - (x2 - x0)*(y1 - y0));

		JF(0, 0) = x1 - x0;
		JF(0, 1) = x2 - x0;
		JF(1, 0) = y1 - y0;
		JF(1, 1) = y2 - y0;

		affineMappingMatrix[4 * k_index + 0] = JF(0, 0);
		affineMappingMatrix[4 * k_index + 1] = JF(0, 1);
		affineMappingMatrix[4 * k_index + 2] = JF(1, 0);
		affineMappingMatrix[4 * k_index + 3] = JF(1, 1);

		Eigen::MatrixXd itJF = (JF.inverse()).transpose();
		//Eigen::VectorXd origNormal = itJF * normal / (itJF * normal).norm();



		for (unsigned El = 0; El < 3; El++) {


			e_pointer const E = K->edges[El];

			real const a = (real) 0.0;
			real const b = (real)El != 0 ? 1.0 : sqrt(2.0);

			real const c = (real)(b - a) / 2.0;
			real const d = (real)(b + a) / 2.0;

			//Vector<real> const normal = normals.getColumn(El);
			Eigen::VectorXd const referenceNormal = normals.col(El);

			real const denominator = (itJF * normals.col(El)).norm() * detJF;



			for (unsigned n = 0; n < number_of_quadrature_points; n++) {


				real const x = (real)quadrature.points[n] * c + d;
				real const w = (real)quadrature.weigths[n] * c;

				evaluate_edge_parametrization(x, El, parametrization);
				//evaluate_edge_parametrization_derivative(x, El, parametrizationDerivative);

				//Vector<real> const r = parametrization;
				//Vector<real> const dr = parametrizationDerivative;

				Eigen::VectorXd const r = parametrization;
				//Eigen::VectorXd const dr = parametrizationDerivative;

				real const s = r(0);
				real const t = r(1);
				//real const drNorm = dr.norm();
				real const drNorm = 1.0;


				evaluate_raviartthomas_basis(s, t, basisRaviartThomas);


			
Eigen::Vector2d const physicalNormal = referenceNormal / denominator;

for (unsigned j = 0; j < 8; j++) {

	real const dotNormal = physicalNormal.dot(basisRaviartThomas.col(j));

	physicalNormalDotPhysicalRaviartThomasBasis_quadPoints[j + 8 * (n + number_of_quadrature_points * (El + 3 * k_index))] = dotNormal;

}




			}
		}

	}

	for (unsigned El = 0; El < 3; El++) {


		real const a = (real) 0.0;
		real const b = (real)El != 0 ? 1.0 : sqrt(2.0);

		real const c = (real)(b - a) / 2.0;
		real const d = (real)(b + a) / 2.0;

		//Vector<real> const normal = normals.getColumn(El);
		Eigen::VectorXd const referenceNormal = normals.col(El);


		for (unsigned n = 0; n < number_of_quadrature_points; n++) {


			real const x = (real)quadrature.points[n] * c + d;
			real const w = (real)quadrature.weigths[n] * c;

			if (El == 0) {
				quadrature_points_1D_e0[n] = x;
				quadrature_weights_1D_e0[n] = w;
			}
			else {
				quadrature_points_1D_e12[n] = x;
				quadrature_weights_1D_e12[n] = w;
			}



			evaluate_edge_parametrization(x, El, parametrization);
			//evaluate_edge_parametrization_derivative(x, El, parametrizationDerivative);

			//Vector<real> const r = parametrization;
			//Vector<real> const dr = parametrizationDerivative;

			Eigen::VectorXd const r = parametrization;
			//Eigen::VectorXd const dr = parametrizationDerivative;

			real const s = r(0);
			real const t = r(1);
			//real const drNorm = dr.norm();
			real const drNorm = 1.0;


			evaluate_raviartthomas_basis(s, t, basisRaviartThomas);
			evaluate_polynomial_basis(s, t, basisPolynomial);


			for (unsigned m = 0; m < 3; m++) {


				real const Phim = basisPolynomial(m);

				polynomialBasis_quadPoints.setCoeff(n, 0, El, m) = Phim;

				for (unsigned j = 0; j < 8; j++) {

					//Vector<real> const Wj = basisRaviartThomas.getColumn(j);
					Eigen::VectorXd const Wj = basisRaviartThomas.col(j);

					//real const dotProduct = dot(Wj, normal);
					real const dotProduct = Wj.dot(referenceNormal);


					real const value = w * dotProduct * Phim * drNorm;

					raviartThomasBasis_quadPoints.setCoeff(n, El, j, 0) = Wj(0);
					raviartThomasBasis_quadPoints.setCoeff(n, El, j, 1) = Wj(1);

					raviartThomasBasisDotNormalTimesPolynomialBasis.setCoeff(n, j, El, m) = abs(value) < INTEGRAL_PRECISION ? 0.0 : value;

				}
			}
		}
	}

	//for (unsigned El = 0; El < 3; El++) {
	//	for (unsigned n = 0; n < number_of_quadrature_points; n++) {

	//		real const x = El != 0 ? quadrature_points_1D_e12[n] : quadrature_points_1D_e0[n];

	//		Eigen::VectorXd parametrization(2);
	//		evaluate_edge_parametrization(x, El, parametrization);

	//		real const s = parametrization(0);
	//		real const t = parametrization(1);

	//		edgeQuadraturePoints_s[number_of_quadrature_points * El + n] = s;
	//		edgeQuadraturePoints_t[number_of_quadrature_points * El + n] = t;

	//	}
	//}


	for (unsigned k = 0; k < nk; k++) {

		t_pointer const K = mesh->get_triangle(k);
		unsigned const k_index = K->index;

		for (unsigned El = 0; El < 3; El++) {

			for (unsigned n = 0; n < number_of_quadrature_points; n++) {


				real const ksi = El != 0 ? quadrature_points_1D_e12[n] : quadrature_points_1D_e0[n];

				Eigen::VectorXd parametrization(2);
				evaluate_edge_parametrization(ksi, El, parametrization);

				real const s = parametrization(0);
				real const t = parametrization(1);

				//real const s = edgeQuadraturePoints_s[number_of_quadrature_points * El + n];
				//real const t = edgeQuadraturePoints_t[number_of_quadrature_points * El + n];

				real const JF00 = affineMappingMatrix[4 * k_index + 0];
				real const JF01 = affineMappingMatrix[4 * k_index + 1];
				real const JF10 = affineMappingMatrix[4 * k_index + 2];
				real const JF11 = affineMappingMatrix[4 * k_index + 3];

				real const x = K->vertices[0]->x + JF00 * s + JF01 * t;
				real const y = K->vertices[0]->y + JF10 * s + JF11 * t;

				edgeQuadraturePoints_x[n + number_of_quadrature_points * (El + 3 * k_index)] = x;
				edgeQuadraturePoints_y[n + number_of_quadrature_points * (El + 3 * k_index)] = y;

			}
		}
	}

*/



// assemble R1 * inverseD * H1
/*

void solver::assembleRInverseDH() {



	//return;
	CoeffMatrix2D<3, 3> iDH1_matrix;
	CoeffMatrix2D<3, 3> iDH2_matrix;

	CoeffMatrix2D<3, 3> R1_matrix;
	CoeffMatrix2D<3, 3> R2_matrix;

	iDH1_matrix.setNumberOfElements(nk);
	iDH2_matrix.setNumberOfElements(nk);

	R1_matrix.setNumberOfElements(nk);
	R2_matrix.setNumberOfElements(nk);


	real const coefficient = θ * dt;


	Eigen::Matrix3d block;
	Eigen::Matrix3d block1;
	Eigen::Matrix3d block2;


	for (unsigned k = 0; k < nk; k++) {



		t_pointer const K = mesh->get_triangle(k);

		unsigned const k_index = K->index;



for (unsigned r = 0; r < 3; r++)
	for (unsigned s = 0; s < 3; s++)
		block(r, s) = δij(r, s) - coefficient * σ(k_index, 0, r, s);

Eigen::Matrix3d const inverse_block = block.inverse();



for (unsigned m = 0; m < 3; m++) {
	for (unsigned Ei = 0; Ei < 3; Ei++) {

		block1(m, Ei) = λ(k_index, 0, m, Ei);
		block2(m, Ei) = λ(k_index, 1, m, Ei);

	}
}

block1 *= -coefficient;
block2 *= -coefficient;




Eigen::Matrix3d const iDH1block = inverse_block * block1;
Eigen::Matrix3d const iDH2block = inverse_block * block2;




for (unsigned i = 0; i < 3; i++) {

	iDH1_matrix.setCoeff(k_index, i, 0) = iDH1block(i, 0);
	iDH1_matrix.setCoeff(k_index, i, 1) = iDH1block(i, 1);
	iDH1_matrix.setCoeff(k_index, i, 2) = iDH1block(i, 2);

	iDH2_matrix.setCoeff(k_index, i, 0) = iDH2block(i, 0);
	iDH2_matrix.setCoeff(k_index, i, 1) = iDH2block(i, 1);
	iDH2_matrix.setCoeff(k_index, i, 2) = iDH2block(i, 2);

}



for (unsigned Ei = 0; Ei < 3; Ei++) {

	e_pointer const E = K->edges[Ei];
	unsigned const e_index = E->index;

	unsigned const dof0 = LI(K, E, 0);
	unsigned const dof1 = LI(K, E, 1);


	for (unsigned m = 0; m < 3; m++) {


		real AB1 = 0.0;
		real AB2 = 0.0;

		for (unsigned i = 0; i < 8; i++) {

			AB1 += α(k_index, 0, i, dof0)*β(k_index, 0, i, m);
			AB2 += α(k_index, 0, i, dof1)*β(k_index, 0, i, m);

		}

		real const val1 = AB1 / viscosities[k_index];
		real const val2 = AB2 / viscosities[k_index];

		//delete R1, R2 and make it in this format only ?

		R1_matrix.setCoeff(k_index, Ei, m) = abs(val1) < INTEGRAL_PRECISION ? 0.0 : val1;
		R2_matrix.setCoeff(k_index, Ei, m) = abs(val2) < INTEGRAL_PRECISION ? 0.0 : val2;

	}
}
	}



	R1iDH1.setZero();
	R1iDH2.setZero();
	R2iDH1.setZero();
	R2iDH2.setZero();

	for (unsigned k = 0; k < nk; k++) {



		t_pointer const K = mesh->get_triangle(k);
		unsigned const k_index = K->index;


		for (unsigned ei = 0; ei < 3; ei++) {


			e_pointer const Ei = K->edges[ei];
			unsigned const e_index_i = Ei->index;


			if (Ei->marker == E_MARKER::DIRICHLET) {

				R1iDH1.coeffRef(e_index_i, e_index_i) = 0.0;
				R1iDH2.coeffRef(e_index_i, e_index_i) = 0.0;

				R2iDH1.coeffRef(e_index_i, e_index_i) = 0.0;
				R2iDH2.coeffRef(e_index_i, e_index_i) = 0.0;

				continue;

			}



			for (unsigned ej = 0; ej < 3; ej++) {


				e_pointer const Ej = K->edges[ej];
				unsigned const e_index_j = Ej->index;


				real sum11 = 0.0;
				real sum12 = 0.0;
				real sum21 = 0.0;
				real sum22 = 0.0;

				for (unsigned i = 0; i < 3; i++) {

					sum11 += R1_matrix(k_index, ei, i) * iDH1_matrix(k_index, i, ej);
					sum12 += R1_matrix(k_index, ei, i) * iDH2_matrix(k_index, i, ej);
					sum21 += R2_matrix(k_index, ei, i) * iDH1_matrix(k_index, i, ej);
					sum22 += R2_matrix(k_index, ei, i) * iDH2_matrix(k_index, i, ej);

				}

				R1iDH1.coeffRef(e_index_i, e_index_j) += sum11;
				R1iDH2.coeffRef(e_index_i, e_index_j) += sum12;
				R2iDH1.coeffRef(e_index_i, e_index_j) += sum21;
				R2iDH2.coeffRef(e_index_i, e_index_j) += sum22;

			}
		}
	}


	



	//for (unsigned e = 0; e < ne; e++) {


	//	e_pointer const E = mesh->get_edge(e);
	//	unsigned const e_index = E->index;
	//	E_MARKER const e_marker = E->marker;


	//	if (e_marker == E_MARKER::DIRICHLET) {


	//		for (unsigned neighborElement = 0; neighborElement < 2; neighborElement++) {

	//			t_pointer const K = E->neighbors[neighborElement];

	//			if (!K)
	//				continue;

	//			R1iDH1.coeffRef(e_index, e_index) = 0.0;
	//			R1iDH2.coeffRef(e_index, e_index) = 0.0;

	//			R2iDH1.coeffRef(e_index, e_index) = 0.0;
	//			R2iDH2.coeffRef(e_index, e_index) = 0.0;

	//			for (unsigned ee = 0; ee < 3; ee++) {

	//				e_pointer const EE = K->edges[ee];

	//				if (EE == E)
	//					continue;

	//				unsigned const e_index_loc = EE->index;


	//				R1iDH1.coeffRef(e_index, e_index_loc) = 0.0;
	//				R1iDH2.coeffRef(e_index, e_index_loc) = 0.0;

	//				R2iDH1.coeffRef(e_index, e_index_loc) = 0.0;
	//				R2iDH2.coeffRef(e_index, e_index_loc) = 0.0;

	//			}
	//		}

	//		continue;

	//	}

	//	real sum1_a_1 = 0.0;
	//	real sum2_a_1 = 0.0;
	//	real sum3_a_1 = 0.0;

	//	real sum1_b_1 = 0.0;
	//	real sum2_b_1 = 0.0;
	//	real sum3_b_1 = 0.0;

	//	real sum1_c_1 = 0.0;
	//	real sum2_c_1 = 0.0;
	//	real sum3_c_1 = 0.0;

	//	real sum1_a_1 = 0.0;
	//	real sum2_a_1 = 0.0;
	//	real sum3_a_1 = 0.0;

	//	real sum1_a_2 = 0.0;
	//	real sum2_a_2 = 0.0;
	//	real sum3_a_2 = 0.0;



	//	real sum1_1 = 0.0;
	//	real sum2_1 = 0.0;
	//	real sum3_1 = 0.0;

	//	for (unsigned neighborElement = 0; neighborElement < 2; neighborElement++) {

	//		t_pointer const K = E->neighbors[neighborElement];

	//		if (!K)
	//			continue;

	//		unsigned const k_index = K->index;

	//		for (unsigned ee = 0; ee < 3; ee++) {


	//			e_pointer const EE = K->edges[ee];

	//			if (EE == E)
	//				continue;

	//			unsigned const e_index_loc = EE->index;


	//			for (unsigned i = 0; i < 3; i++) {

	//				sum1_1 = R1_matrix(k_index, ee, i)*iDH1_matrix(k_index, i, 0);
	//				sum2_1 = R1_matrix(k_index, ee, i)*iDH1_matrix(k_index, i, 1);
	//				sum3_1 = R1_matrix(k_index, ee, i)*iDH1_matrix(k_index, i, 2);


	//				//sum1_a_1 = R1_matrix(k_index, 0, i)*iDH1_matrix(k_index, i, 0);
	//				//sum2_a_1 = R1_matrix(k_index, 0, i)*iDH1_matrix(k_index, i, 1);
	//				//sum3_a_1 = R1_matrix(k_index, 0, i)*iDH1_matrix(k_index, i, 2);

	//				//sum1_b_1 = R1_matrix(k_index, 1, i)*iDH1_matrix(k_index, i, 0);
	//				//sum2_b_1 = R1_matrix(k_index, 1, i)*iDH1_matrix(k_index, i, 1);
	//				//sum3_b_1 = R1_matrix(k_index, 1, i)*iDH1_matrix(k_index, i, 2);

	//				//sum1_c_1 = R1_matrix(k_index, 2, i)*iDH1_matrix(k_index, i, 0);
	//				//sum2_c_1 = R1_matrix(k_index, 2, i)*iDH1_matrix(k_index, i, 1);
	//				//sum3_c_1 = R1_matrix(k_index, 2, i)*iDH1_matrix(k_index, i, 2);

	//			}


	//			R1iDH1.coeffRef(e_index, e_index_loc) += sum1_1;
	//			//	R1iDH2;
	//			//R2iDH1;
	//			//R2iDH2;

	//		}
	//	}

	//	//R1iDH1.coeffRef(e_index,)
	//	//R1iDH2;
	//	//R2iDH1;
	//	//R2iDH2;


	//}
	//


	
	//for (unsigned e = 0; e < ne; e++) {


	//	e_pointer const E = mesh->get_edge(e);
	//	unsigned const e_index = E->index;
	//	E_MARKER const e_marker = E->marker;


	//	if (e_marker == E_MARKER::DIRICHLET) {


	//		for (unsigned neighborElement = 0; neighborElement < 2; neighborElement++) {

	//			t_pointer const K = E->neighbors[neighborElement];

	//			if (!K)
	//				continue;

	//			R1iDH1.coeffRef(e_index, e_index) = 0.0;
	//			R1iDH2.coeffRef(e_index, e_index) = 0.0;

	//			R2iDH1.coeffRef(e_index, e_index) = 0.0;
	//			R2iDH2.coeffRef(e_index, e_index) = 0.0;

	//			for (unsigned ee = 0; ee < 3; ee++) {

	//				e_pointer const EE = K->edges[ee];

	//				if (EE == E)
	//					continue;

	//				unsigned const e_index_loc = EE->index;


	//				R1iDH1.coeffRef(e_index, e_index_loc) = 0.0;
	//				R1iDH2.coeffRef(e_index, e_index_loc) = 0.0;

	//				R2iDH1.coeffRef(e_index, e_index_loc) = 0.0;
	//				R2iDH2.coeffRef(e_index, e_index_loc) = 0.0;

	//			}
	//		}

	//		continue;

	//	}

	//	real sum1_a_1 = 0.0;
	//	real sum2_a_1 = 0.0;
	//	real sum3_a_1 = 0.0;

	//	real sum1_b_1 = 0.0;
	//	real sum2_b_1 = 0.0;
	//	real sum3_b_1 = 0.0;

	//	real sum1_c_1 = 0.0;
	//	real sum2_c_1 = 0.0;
	//	real sum3_c_1 = 0.0;

	//	real sum1_a_1 = 0.0;
	//	real sum2_a_1 = 0.0;
	//	real sum3_a_1 = 0.0;

	//	real sum1_a_2 = 0.0;
	//	real sum2_a_2 = 0.0;
	//	real sum3_a_2 = 0.0;



	//	real sum1_1 = 0.0;
	//	real sum2_1 = 0.0;
	//	real sum3_1 = 0.0;

	//	for (unsigned neighborElement = 0; neighborElement < 2; neighborElement++) {

	//		t_pointer const K = E->neighbors[neighborElement];

	//		if (!K)
	//			continue;

	//		unsigned const k_index = K->index;

	//		for (unsigned ee = 0; ee < 3; ee++) {


	//			e_pointer const EE = K->edges[ee];

	//			if (EE == E)
	//				continue;

	//			unsigned const e_index_loc = EE->index;


	//			for (unsigned i = 0; i < 3; i++) {

	//				sum1_1 = R1_matrix(k_index, ee, i)*iDH1_matrix(k_index, i, 0);
	//				sum2_1 = R1_matrix(k_index, ee, i)*iDH1_matrix(k_index, i, 1);
	//				sum3_1 = R1_matrix(k_index, ee, i)*iDH1_matrix(k_index, i, 2);


	//				//sum1_a_1 = R1_matrix(k_index, 0, i)*iDH1_matrix(k_index, i, 0);
	//				//sum2_a_1 = R1_matrix(k_index, 0, i)*iDH1_matrix(k_index, i, 1);
	//				//sum3_a_1 = R1_matrix(k_index, 0, i)*iDH1_matrix(k_index, i, 2);

	//				//sum1_b_1 = R1_matrix(k_index, 1, i)*iDH1_matrix(k_index, i, 0);
	//				//sum2_b_1 = R1_matrix(k_index, 1, i)*iDH1_matrix(k_index, i, 1);
	//				//sum3_b_1 = R1_matrix(k_index, 1, i)*iDH1_matrix(k_index, i, 2);

	//				//sum1_c_1 = R1_matrix(k_index, 2, i)*iDH1_matrix(k_index, i, 0);
	//				//sum2_c_1 = R1_matrix(k_index, 2, i)*iDH1_matrix(k_index, i, 1);
	//				//sum3_c_1 = R1_matrix(k_index, 2, i)*iDH1_matrix(k_index, i, 2);

	//			}


	//			R1iDH1.coeffRef(e_index, e_index_loc) += sum1_1;
	//			//	R1iDH2;
	//			//R2iDH1;
	//			//R2iDH2;

	//		}
	//	}

	//	//R1iDH1.coeffRef(e_index,)
	//	//R1iDH2;
	//	//R2iDH1;
	//	//R2iDH2;


	//}
	//


	



	//for (unsigned e = 0; e < ne; e++) {


	//	e_pointer const E = mesh->get_edge(e);
	//	unsigned const e_index = E->index;

	//	for (unsigned neighborElement = 0; neighborElement < 2; neighborElement++) {

	//		t_pointer const K = E->neighbors[neighborElement];

	//		if (!K)
	//			continue;


	//		unsigned const k_index = K->index;

	//		unsigned const dof0 = LI(K, E, 0);
	//		unsigned const dof1 = LI(K, E, 1);


	//		for (unsigned m = 0; m < 3; m++) {


	//			real AB1 = 0.0;
	//			real AB2 = 0.0;

	//			for (unsigned i = 0; i < 8; i++) {

	//				AB1 += α(k_index, 0, i, dof0)*β(k_index, 0, i, m);
	//				AB2 += α(k_index, 0, i, dof1)*β(k_index, 0, i, m);

	//			}

	//			real const val1 = AB1 / viscosities[k_index];
	//			real const val2 = AB2 / viscosities[k_index];

	//			//delete R1, R2 and make it in this format only ?

	//			R1_matrix.setCoeff(e_index, neighborElement, m) = abs(val1) < INTEGRAL_PRECISION ? 0.0 : val1;
	//			R2_matrix.setCoeff(e_index, neighborElement, m) = abs(val2) < INTEGRAL_PRECISION ? 0.0 : val2;

	//		}
	//	}
	//}


	//for (unsigned e = 0; e < ne; e++) {


	//	e_pointer const E = mesh->get_edge(e);
	//	unsigned const e_index = E->index;
	//	E_MARKER const e_marker = E->marker;


	//	if (e_marker == E_MARKER::DIRICHLET) {


	//		for (unsigned neighborElement = 0; neighborElement < 2; neighborElement++) {

	//			t_pointer const K = E->neighbors[neighborElement];

	//			if (!K)
	//				continue;

	//			R1iDH1.coeffRef(e_index, e_index) = 0.0;
	//			R1iDH2.coeffRef(e_index, e_index) = 0.0;

	//			R2iDH1.coeffRef(e_index, e_index) = 0.0;
	//			R2iDH2.coeffRef(e_index, e_index) = 0.0;

	//			for (unsigned ee = 0; ee < 3; ee++) {

	//				e_pointer const EE = K->edges[ee];

	//				if (EE == E)
	//					continue;

	//				unsigned const e_index_loc = EE->index;


	//				R1iDH1.coeffRef(e_index, e_index_loc) = 0.0;
	//				R1iDH2.coeffRef(e_index, e_index_loc) = 0.0;

	//				R2iDH1.coeffRef(e_index, e_index_loc) = 0.0;
	//				R2iDH2.coeffRef(e_index, e_index_loc) = 0.0;

	//			}
	//		}

	//		continue;

	//	}


	//	for (unsigned neighborElement = 0; neighborElement < 2; neighborElement++) {

	//		t_pointer const K = E->neighbors[neighborElement];

	//		if (!K)
	//			continue;

	//		unsigned const k_index = K->index;

	//		for (unsigned ee = 0; ee < 3; ee++) {


	//			e_pointer const EE = K->edges[ee];

	//			if (EE == E)
	//				continue;

	//			unsigned const e_index_loc = EE->index;


	//			real sum1 = 0.0;

	//			for (unsigned i = 0; i < 3; i++)
	//				sum1 += R1_matrix(e_index, neighborElement, i)*iDH1_matrix(k_index, i, neighborElement);


	//		}
	//		R1iDH1;
	//		R1iDH2;
	//		R2iDH1;
	//		R2iDH2;



	//	}


	//}













	//for (unsigned e = 0; e < ne; e++) {


	//	e_pointer const E = mesh->get_edge(e);
	//	unsigned const e_index = E->index;
	//	E_MARKER const e_marker = E->marker;


	//	if (e_marker == E_MARKER::DIRICHLET) {


	//		for (unsigned neighborElement = 0; neighborElement < 2; neighborElement++) {

	//			t_pointer const K = E->neighbors[neighborElement];

	//			if (!K)
	//				continue;

	//			R1iDH1.coeffRef(e_index, e_index) = 0.0;
	//			R1iDH2.coeffRef(e_index, e_index) = 0.0;

	//			R2iDH1.coeffRef(e_index, e_index) = 0.0;
	//			R2iDH2.coeffRef(e_index, e_index) = 0.0;

	//			for (unsigned ee = 0; ee < 3; ee++) {

	//				e_pointer const EE = K->edges[ee];

	//				if (EE == E)
	//					continue;

	//				unsigned const e_index_loc = EE->index;


	//				R1iDH1.coeffRef(e_index, e_index_loc) = 0.0;
	//				R1iDH2.coeffRef(e_index, e_index_loc) = 0.0;

	//				R2iDH1.coeffRef(e_index, e_index_loc) = 0.0;
	//				R2iDH2.coeffRef(e_index, e_index_loc) = 0.0;

	//			}
	//		}

	//		continue;

	//	}





	//	real sum1_1 = 0.0;
	//	real sum2_1 = 0.0;
	//	real sum3_1 = 0.0;

	//	real sum1_2 = 0.0;
	//	real sum2_2 = 0.0;
	//	real sum3_2 = 0.0;

	//	for (unsigned neighborElement = 0; neighborElement < 2; neighborElement++) {

	//		t_pointer const K = E->neighbors[neighborElement];

	//		if (!K)
	//			continue;

	//		unsigned const k_index = K->index;
	//		unsigned const start_index = 3 * k_index;

	//		unsigned const dof0 = LI(K, E, 0);
	//		unsigned const dof1 = LI(K, E, 1);

	//		for (unsigned ee = 0; ee < 3; ee++) {


	//			e_pointer const EE = K->edges[ee];

	//			if (EE == E)
	//				continue;


	//			unsigned const e_index_loc = EE->index;

	//			for (unsigned r = 0; r < 3; r++)
	//				for (unsigned s = 0; s < 3; s++)
	//					block_σ(r, s) = δij(r, s) - coefficient * σ(k_index, 0, r, s);


	//			Eigen::Matrix3d inverse_block = block_σ.inverse();

	//			for (unsigned i = 0; i < 3; i++)
	//				for (unsigned j = 0; j < 3; j++)
	//					inverse_block.coeffRef(i, j) = abs(inverse_block(i, j)) < INTEGRAL_PRECISION ? 0.0 : inverse_block(i, j);



	//			for (unsigned m = 0; m < 3; m++) {
	//				for (unsigned Ei = 0; Ei < 3; Ei++) {

	//					block_λ1(m, Ei) = λ(k_index, 0, m, Ei);
	//					block_λ2(m, Ei) = λ(k_index, 1, m, Ei);

	//				}
	//			}

	//			block_λ1 *= -coefficient;
	//			block_λ2 *= -coefficient;



	//			Eigen::Matrix3d const iDH1block = inverse_block * block_λ1;
	//			Eigen::Matrix3d const iDH2block = inverse_block * block_λ2;

	//			//
	//			//for (unsigned i = 0; i < 3; i++) {
	//			//
	//			//
	//			//	_iDH1.coeffRef(3 * k_index + i, e_index0) = iDH1block.coeff(i, 0);
	//			//	_iDH1.coeffRef(3 * k_index + i, e_index1) = iDH1block.coeff(i, 1);
	//			//	_iDH1.coeffRef(3 * k_index + i, e_index2) = iDH1block.coeff(i, 2);
	//			//
	//			//	_iDH2.coeffRef(3 * k_index + i, e_index0) = iDH2block.coeff(i, 0);
	//			//	_iDH2.coeffRef(3 * k_index + i, e_index1) = iDH2block.coeff(i, 1);
	//			//	_iDH2.coeffRef(3 * k_index + i, e_index2) = iDH2block.coeff(i, 2);
	//			//
	//			//}



	//			Eigen::Matrix3d R1block;
	//			Eigen::Matrix3d R2block;


	//			Eigen::Vector3d vector_horizontal1_1;
	//			Eigen::Vector3d vector_horizontal2_1;
	//			Eigen::Vector3d vector_horizontal3_1;

	//			Eigen::Vector3d vector_horizontal1_2;
	//			Eigen::Vector3d vector_horizontal2_2;
	//			Eigen::Vector3d vector_horizontal3_2;



	//			for (unsigned m = 0; m < 3; m++) {


	//				real AB1 = 0.0;
	//				real AB2 = 0.0;

	//				for (unsigned i = 0; i < 8; i++) {

	//					AB1 += α(k_index, 0, i, dof0)*β(k_index, 0, i, m);
	//					AB2 += α(k_index, 0, i, dof1)*β(k_index, 0, i, m);

	//				}

	//				real const val1 = AB1 / viscosities[k_index];
	//				real const val2 = AB2 / viscosities[k_index];

	//				vector_horizontal1_1.coeffRef(m) = abs(val1) < INTEGRAL_PRECISION ? 0.0 : val1;
	//				vector_horizontal1_2.coeffRef(m) = abs(val2) < INTEGRAL_PRECISION ? 0.0 : val2;

	//			}

	//			sum1_1 += vector_horizontal1_1.transpose() * iDH1block.col(0);
	//			sum1_2 += vector_horizontal1_2.transpose() * iDH1block.col(0);




	//			R1iDH1.coeffRef(e_index, e_index_loc) = vec1.transpose() * iDH1block.col(0);



	//			sum1 += val1 * iDH1block(j, m);
	//			sum2 += val1 * iDH2block(j, m);

	//			sum3 += val2 * iDH1block(j, m);
	//			sum4 += val2 * iDH2block(j, m);

	//			R1iDH2;
	//			R2iDH1;
	//			R2iDH2;



	//			for (unsigned m = 0; m < 3; m++) {


	//				real AB1 = 0.0;
	//				real AB2 = 0.0;

	//				for (unsigned i = 0; i < 8; i++) {

	//					AB1 += α(k_index, 0, i, dof0)*β(k_index, 0, i, m);
	//					AB2 += α(k_index, 0, i, dof1)*β(k_index, 0, i, m);

	//				}

	//				real const val1 = AB1 / viscosities[k_index];
	//				real const val2 = AB2 / viscosities[k_index];

	//				R1.coeffRef(e_index, 3 * k_index + m) = abs(val1) < INTEGRAL_PRECISION ? 0.0 : val1;
	//				R2.coeffRef(e_index, 3 * k_index + m) = abs(val2) < INTEGRAL_PRECISION ? 0.0 : val2;

	//			}


	//		}
	//	}
	//}



































	//real const coefficient = θ * dt;

	//iD.setZero();
	//H1.setZero();
	//H2.setZero();

	//assemble_λ();
	//assemble_σ();

	//Eigen::Matrix3d block;
	//Eigen::Matrix3d block1;
	//Eigen::Matrix3d block2;

	//for (unsigned k = 0; k < nk; k++) {



	//	t_pointer const K = mesh->get_triangle(k);

	//	unsigned const k_index = K->index;
	//	unsigned const start_index = 3 * k_index;

	//	e_pointer const E0 = K->edges[0];
	//	e_pointer const E1 = K->edges[1];
	//	e_pointer const E2 = K->edges[2];

	//	unsigned const e_index0 = E0->index;
	//	unsigned const e_index1 = E1->index;
	//	unsigned const e_index2 = E2->index;



	//	for (unsigned r = 0; r < 3; r++)
	//		for (unsigned s = 0; s < 3; s++)
	//			block(r, s) = δij(r, s) - coefficient * σ(k_index, 0, r, s);


	//	Eigen::Matrix3d inverse_block = block.inverse();

	//	for (unsigned i = 0; i < 3; i++)
	//		for (unsigned j = 0; j < 3; j++)
	//			inverse_block.coeffRef(i, j) = abs(inverse_block(i, j)) < INTEGRAL_PRECISION ? 0.0 : inverse_block(i, j);




	//	for (unsigned m = 0; m < 3; m++) {
	//		for (unsigned Ei = 0; Ei < 3; Ei++) {

	//			block1(m, Ei) = λ(k_index, 0, m, Ei);
	//			block2(m, Ei) = λ(k_index, 1, m, Ei);

	//		}
	//	}

	//	block1 *= -coefficient;
	//	block2 *= -coefficient;




	//	Eigen::Matrix3d const iDH1block = inverse_block * block1;
	//	Eigen::Matrix3d const iDH2block = inverse_block * block2;



	//	for (unsigned i = 0; i < 3; i++) {


	//		_iDH1.coeffRef(3 * k_index + i, e_index0) = iDH1block.coeff(i, 0);
	//		_iDH1.coeffRef(3 * k_index + i, e_index1) = iDH1block.coeff(i, 1);
	//		_iDH1.coeffRef(3 * k_index + i, e_index2) = iDH1block.coeff(i, 2);

	//		_iDH2.coeffRef(3 * k_index + i, e_index0) = iDH2block.coeff(i, 0);
	//		_iDH2.coeffRef(3 * k_index + i, e_index1) = iDH2block.coeff(i, 1);
	//		_iDH2.coeffRef(3 * k_index + i, e_index2) = iDH2block.coeff(i, 2);

	//	}
	//}
















	//for (unsigned e = 0; e < ne; e++) {


	//	e_pointer const E = mesh->get_edge(e);
	//	unsigned const e_index = E->index;
	//	E_MARKER const e_marker = E->marker;


	//	if (e_marker == E_MARKER::DIRICHLET) {


	//		for (unsigned neighborElement = 0; neighborElement < 2; neighborElement++) {

	//			t_pointer const K = E->neighbors[neighborElement];

	//			if (!K)
	//				continue;

	//			unsigned const k_index = K->index;

	//			for (unsigned m = 0; m < 3; m++) {

	//				R1.coeffRef(e_index, 3 * k_index + m) = 0.0;
	//				R2.coeffRef(e_index, 3 * k_index + m) = 0.0;

	//			}
	//		}

	//		continue;

	//	}





	//	real sum1 = 0.0;
	//	real sum2 = 0.0;
	//	real sum3 = 0.0;
	//	real sum4 = 0.0;

	//	for (unsigned neighborElement = 0; neighborElement < 2; neighborElement++) {

	//		t_pointer const K = E->neighbors[neighborElement];

	//		if (!K)
	//			continue;

	//		unsigned const k_index = K->index;

	//		unsigned const dof0 = LI(K, E, 0);
	//		unsigned const dof1 = LI(K, E, 1);

	//		for (unsigned ee = 0; ee < 3; ee++) {


	//			e_pointer const EE = K->edges[ee];

	//			if (EE == E)
	//				continue;

	//			unsigned const e_index_loc = EE->index;



	//			unsigned const start_index = 3 * k_index;

	//			//e_pointer const E0 = K->edges[0];
	//			//e_pointer const E1 = K->edges[1];
	//			//e_pointer const E2 = K->edges[2];

	//			//unsigned const e_index0 = E0->index;
	//			//unsigned const e_index1 = E1->index;
	//			//unsigned const e_index2 = E2->index;


	//			for (unsigned r = 0; r < 3; r++)
	//				for (unsigned s = 0; s < 3; s++)
	//					block(r, s) = δij(r, s) - coefficient * σ(k_index, 0, r, s);


	//			Eigen::Matrix3d inverse_block = block.inverse();

	//			for (unsigned i = 0; i < 3; i++)
	//				for (unsigned j = 0; j < 3; j++)
	//					inverse_block.coeffRef(i, j) = abs(inverse_block(i, j)) < INTEGRAL_PRECISION ? 0.0 : inverse_block(i, j);




	//			for (unsigned m = 0; m < 3; m++) {
	//				for (unsigned Ei = 0; Ei < 3; Ei++) {

	//					block1(m, Ei) = λ(k_index, 0, m, Ei);
	//					block2(m, Ei) = λ(k_index, 1, m, Ei);

	//				}
	//			}

	//			block1 *= -coefficient;
	//			block2 *= -coefficient;



	//			Eigen::Matrix3d const iDH1block = inverse_block * block1;
	//			Eigen::Matrix3d const iDH2block = inverse_block * block2;


	//			//for (unsigned i = 0; i < 3; i++) {


	//			//	_iDH1.coeffRef(3 * k_index + i, e_index0) = iDH1block.coeff(i, 0);
	//			//	_iDH1.coeffRef(3 * k_index + i, e_index1) = iDH1block.coeff(i, 1);
	//			//	_iDH1.coeffRef(3 * k_index + i, e_index2) = iDH1block.coeff(i, 2);

	//			//	_iDH2.coeffRef(3 * k_index + i, e_index0) = iDH2block.coeff(i, 0);
	//			//	_iDH2.coeffRef(3 * k_index + i, e_index1) = iDH2block.coeff(i, 1);
	//			//	_iDH2.coeffRef(3 * k_index + i, e_index2) = iDH2block.coeff(i, 2);

	//			//}



	//			Eigen::Vector3d vec1;
	//			Eigen::Vector3d vec2;

	//			for (unsigned m = 0; m < 3; m++) {


	//				real AB1 = 0.0;
	//				real AB2 = 0.0;

	//				for (unsigned i = 0; i < 8; i++) {

	//					AB1 += α(k_index, 0, i, dof0)*β(k_index, 0, i, m);
	//					AB2 += α(k_index, 0, i, dof1)*β(k_index, 0, i, m);

	//				}

	//				real const val1 = AB1 / viscosities[k_index];
	//				real const val2 = AB2 / viscosities[k_index];

	//				vec1.coeffRef(m) = abs(val1) < INTEGRAL_PRECISION ? 0.0 : val1;
	//				vec2.coeffRef(m) = abs(val2) < INTEGRAL_PRECISION ? 0.0 : val2;

	//			}

	//			sum1 += vec1.transpose() * iDH1block.col(0);
	//			sum2 += vec1.transpose() * iDH1block.col(0);
	//			sum3 += 0.0;
	//			//sum4 += 0.0;

	//			R1iDH1.coeffRef(e_index, e_index_loc) = vec1.transpose() * iDH1block.col(0);



	//			sum1 += val1 * iDH1block(j, m);
	//			sum2 += val1 * iDH2block(j, m);

	//			sum3 += val2 * iDH1block(j, m);
	//			sum4 += val2 * iDH2block(j, m);

	//				R1iDH2;
	//			R2iDH1;
	//			R2iDH2;



	//			for (unsigned m = 0; m < 3; m++) {


	//				real AB1 = 0.0;
	//				real AB2 = 0.0;

	//				for (unsigned i = 0; i < 8; i++) {

	//					AB1 += α(k_index, 0, i, dof0)*β(k_index, 0, i, m);
	//					AB2 += α(k_index, 0, i, dof1)*β(k_index, 0, i, m);

	//				}

	//				real const val1 = AB1 / viscosities[k_index];
	//				real const val2 = AB2 / viscosities[k_index];

	//				R1.coeffRef(e_index, 3 * k_index + m) = abs(val1) < INTEGRAL_PRECISION ? 0.0 : val1;
	//				R2.coeffRef(e_index, 3 * k_index + m) = abs(val2) < INTEGRAL_PRECISION ? 0.0 : val2;

	//			}


	//		}
	//	}
	//}

	
//};


*/

/*

void solver::assembleRInverseDH_whole() {


	real const coefficient = θ * dt;

	Eigen::Matrix3d block;
	Eigen::Matrix3d block1;
	Eigen::Matrix3d block2;


	R1iDH1.setZero();
	R1iDH2.setZero();
	R2iDH1.setZero();
	R2iDH2.setZero();

	for (unsigned k = 0; k < nk; k++) {



		t_pointer const K = mesh->get_triangle(k);

		unsigned const k_index = K->index;




for (unsigned r = 0; r < 3; r++)
	for (unsigned s = 0; s < 3; s++)
		block(r, s) = δij(r, s) - coefficient * σ(k_index, 0, r, s);

Eigen::Matrix3d const inverse_block = block.inverse();




for (unsigned m = 0; m < 3; m++) {
	for (unsigned Ei = 0; Ei < 3; Ei++) {

		block1(m, Ei) = λ(k_index, 0, m, Ei);
		block2(m, Ei) = λ(k_index, 1, m, Ei);

	}
}

block1 *= -coefficient;
block2 *= -coefficient;




Eigen::Matrix3d const iDH1block = inverse_block * block1;
Eigen::Matrix3d const iDH2block = inverse_block * block2;





for (unsigned ei = 0; ei < 3; ei++) {


	e_pointer const Ei = K->edges[ei];
	unsigned const e_index_i = K->edges[ei]->index;


	if (Ei->marker == E_MARKER::DIRICHLET) {

		R1iDH1.coeffRef(e_index_i, e_index_i) = 0.0;
		R1iDH2.coeffRef(e_index_i, e_index_i) = 0.0;

		R2iDH1.coeffRef(e_index_i, e_index_i) = 0.0;
		R2iDH2.coeffRef(e_index_i, e_index_i) = 0.0;

		continue;

	}


	for (unsigned ej = 0; ej < 3; ej++) {


		unsigned const e_index_j = K->edges[ej]->index;


		real sum11 = 0.0;
		real sum12 = 0.0;
		real sum21 = 0.0;
		real sum22 = 0.0;

		// Number of degrees of freedom of internal pressure
		for (unsigned m = 0; m < 3; m++) {

			sum11 += R1_matrix(k_index, ei, m) * iDH1block(m, ej);
			sum12 += R1_matrix(k_index, ei, m) * iDH2block(m, ej);
			sum21 += R2_matrix(k_index, ei, m) * iDH1block(m, ej);
			sum22 += R2_matrix(k_index, ei, m) * iDH2block(m, ej);

		}

		R1iDH1.coeffRef(e_index_i, e_index_j) += sum11;
		R1iDH2.coeffRef(e_index_i, e_index_j) += sum12;
		R2iDH1.coeffRef(e_index_i, e_index_j) += sum21;
		R2iDH2.coeffRef(e_index_i, e_index_j) += sum22;

	}
}
	}

};

*/

/*

void solver::assembleRInverseDH_whole_immidiate() {


	real const coefficient = θ * dt;

	Eigen::Matrix3d block;
	Eigen::Matrix3d block1;
	Eigen::Matrix3d block2;


	//R1iDH1.setZero();
	//R1iDH2.setZero();
	//R2iDH1.setZero();
	//R2iDH2.setZero();

	internalPressureSystem_smat.setZero();

	for (unsigned k = 0; k < nk; k++) {



		t_pointer const K = mesh->get_triangle(k);

		unsigned const k_index = K->index;




for (unsigned r = 0; r < 3; r++)
	for (unsigned s = 0; s < 3; s++)
		block(r, s) = δij(r, s) - coefficient * σ(k_index, 0, r, s);

Eigen::Matrix3d const inverse_block = block.inverse();



for (unsigned m = 0; m < 3; m++) {
	for (unsigned Ei = 0; Ei < 3; Ei++) {

		block1(m, Ei) = λ(k_index, 0, m, Ei);
		block2(m, Ei) = λ(k_index, 1, m, Ei);

	}
}

block1 *= -coefficient;
block2 *= -coefficient;



Eigen::Matrix3d const iDH1block = inverse_block * block1;
Eigen::Matrix3d const iDH2block = inverse_block * block2;




for (unsigned ei = 0; ei < 3; ei++) {


	e_pointer const Ei = K->edges[ei];
	unsigned const e_index_i = K->edges[ei]->index;


	if (Ei->marker == E_MARKER::DIRICHLET) {


		//R1iDH1.coeffRef(e_index_i, e_index_i) = 0.0;
		//R1iDH2.coeffRef(e_index_i, e_index_i) = 0.0;

		//R2iDH1.coeffRef(e_index_i, e_index_i) = 0.0;
		//R2iDH2.coeffRef(e_index_i, e_index_i) = 0.0;


		//R1iDH1.coeffRef(e_index_i, e_index_i) = M_j1_s1.coeff(e_index_i, e_index_i);
		//R1iDH2.coeffRef(e_index_i, e_index_i) = M_j1_s2.coeff(e_index_i, e_index_i);

		//R2iDH1.coeffRef(e_index_i, e_index_i) = M_j2_s1.coeff(e_index_i, e_index_i);
		//R2iDH2.coeffRef(e_index_i, e_index_i) = M_j2_s2.coeff(e_index_i, e_index_i);



		internalPressureSystem_smat.coeffRef(e_index_i + 0, e_index_i + 0) = 0.0;
		internalPressureSystem_smat.coeffRef(e_index_i + 0, e_index_i + ne) = 0.0;
		internalPressureSystem_smat.coeffRef(e_index_i + ne, e_index_i + 0) = 0.0;
		internalPressureSystem_smat.coeffRef(e_index_i + ne, e_index_i + ne) = 0.0;

		//internalPressureSystem_smat.coeffRef(e_index_i + 0, e_index_i + 0) = M_j1_s1.coeff(e_index_i, e_index_i);
		//internalPressureSystem_smat.coeffRef(e_index_i + 0, e_index_i + ne) = M_j1_s2.coeff(e_index_i, e_index_i);
		//internalPressureSystem_smat.coeffRef(e_index_i + ne, e_index_i + 0) = M_j2_s1.coeff(e_index_i, e_index_i);
		//internalPressureSystem_smat.coeffRef(e_index_i + ne, e_index_i + ne) = M_j2_s2.coeff(e_index_i, e_index_i);



		continue;

	}


	for (unsigned ej = 0; ej < 3; ej++) {


		unsigned const e_index_j = K->edges[ej]->index;


		real sum11 = 0.0;
		real sum12 = 0.0;
		real sum21 = 0.0;
		real sum22 = 0.0;

		// Number of degrees of freedom of internal pressure
		for (unsigned m = 0; m < 3; m++) {

			sum11 += R1_matrix(k_index, ei, m) * iDH1block(m, ej);
			sum12 += R1_matrix(k_index, ei, m) * iDH2block(m, ej);
			sum21 += R2_matrix(k_index, ei, m) * iDH1block(m, ej);
			sum22 += R2_matrix(k_index, ei, m) * iDH2block(m, ej);

		}

		//R1iDH1.coeffRef(e_index_i, e_index_j) += sum11;
		//R1iDH2.coeffRef(e_index_i, e_index_j) += sum12;
		//R2iDH1.coeffRef(e_index_i, e_index_j) += sum21;
		//R2iDH2.coeffRef(e_index_i, e_index_j) += sum22;


		internalPressureSystem_smat.coeffRef(e_index_i + 0, e_index_j + 0) += sum11;
		internalPressureSystem_smat.coeffRef(e_index_i + 0, e_index_j + ne) += sum12;
		internalPressureSystem_smat.coeffRef(e_index_i + ne, e_index_j + 0) += sum21;
		internalPressureSystem_smat.coeffRef(e_index_i + ne, e_index_j + ne) += sum22;


		//R1iDH1.coeffRef(e_index_i, e_index_j) += sum11 + M_j1_s1.coeff(e_index_i, e_index_j);
		//R1iDH2.coeffRef(e_index_i, e_index_j) += sum12 + M_j1_s2.coeff(e_index_i, e_index_j);
		//R2iDH1.coeffRef(e_index_i, e_index_j) += sum21 + M_j2_s1.coeff(e_index_i, e_index_j);
		//R2iDH2.coeffRef(e_index_i, e_index_j) += sum22 + M_j2_s2.coeff(e_index_i, e_index_j);


	}
}
	}



	//for (unsigned e = 0; e < ne; e++) {


	//	e_pointer const E = mesh->get_edge(e);
	//	unsigned const e_index = E->index;
	//	E_MARKER const e_marker = E->marker;



	//	//if (e_marker == E_MARKER::DIRICHLET) {

	//	//	continue;

	//	//}


	//	for (unsigned neighborElement = 0; neighborElement < 2; neighborElement++) {


	//		t_pointer const K = E->neighbors[neighborElement];

	//		if (!K)
	//			continue;


	//		// Loop over edges
	//		for (unsigned El = 0; El < 3; El++) {

	//			unsigned const e_local_index_global = K->edges[El]->index;				// Global index of local edge

	//			//if (K->edges[El] == E)
	//			//	continue;


	//			internalPressureSystem_smat.coeffRef(e_index, e_local_index_global) += M_j1_s1.coeff(e_index, e_local_index_global);
	//			internalPressureSystem_smat.coeffRef(e_index, e_local_index_global + ne) += M_j1_s2.coeff(e_index, e_local_index_global);
	//			internalPressureSystem_smat.coeffRef(e_index + ne, e_local_index_global) += M_j2_s1.coeff(e_index, e_local_index_global);
	//			internalPressureSystem_smat.coeffRef(e_index + ne, e_local_index_global + ne) += M_j2_s2.coeff(e_index, e_local_index_global);

	//		}
	//	}
	//}




};


*/

/*


void solver::assembleInverseDH() {



	real const coefficient = θ * dt;

	Eigen::Matrix3d block;
	Eigen::Matrix3d block1;
	Eigen::Matrix3d block2;

	for (unsigned k = 0; k < nk; k++) {



		t_pointer const K = mesh->get_triangle(k);

		unsigned const k_index = K->index;
		unsigned const start_index = 3 * k_index;

		e_pointer const E0 = K->edges[0];
		e_pointer const E1 = K->edges[1];
		e_pointer const E2 = K->edges[2];

		unsigned const e_index0 = E0->index;
		unsigned const e_index1 = E1->index;
		unsigned const e_index2 = E2->index;


for (unsigned r = 0; r < 3; r++)
	for (unsigned s = 0; s < 3; s++)
		block(r, s) = δij(r, s) - coefficient * σ(k_index, 0, r, s);

Eigen::Matrix3d inverse_block = block.inverse();

for (unsigned i = 0; i < 3; i++)
	for (unsigned j = 0; j < 3; j++)
		inverse_block.coeffRef(i, j) = abs(inverse_block(i, j)) < INTEGRAL_PRECISION ? 0.0 : inverse_block(i, j);




for (unsigned m = 0; m < 3; m++) {
	for (unsigned Ei = 0; Ei < 3; Ei++) {

		block1(m, Ei) = λ(k_index, 0, m, Ei);
		block2(m, Ei) = λ(k_index, 1, m, Ei);

	}
}

block1 *= -coefficient;
block2 *= -coefficient;





Eigen::Matrix3d const iDH1block = inverse_block * block1;
Eigen::Matrix3d const iDH2block = inverse_block * block2;



for (unsigned i = 0; i < 3; i++) {


	_iDH1.coeffRef(3 * k_index + i, e_index0) = iDH1block.coeff(i, 0);
	_iDH1.coeffRef(3 * k_index + i, e_index1) = iDH1block.coeff(i, 1);
	_iDH1.coeffRef(3 * k_index + i, e_index2) = iDH1block.coeff(i, 2);

	_iDH2.coeffRef(3 * k_index + i, e_index0) = iDH2block.coeff(i, 0);
	_iDH2.coeffRef(3 * k_index + i, e_index1) = iDH2block.coeff(i, 1);
	_iDH2.coeffRef(3 * k_index + i, e_index2) = iDH2block.coeff(i, 2);

}
	}

};
void solver::assembleRInverseDH() {


	real const coefficient = θ * dt;

	Eigen::Matrix3d block;
	Eigen::Matrix3d block1;
	Eigen::Matrix3d block2;


	for (unsigned k = 0; k < nk; k++) {



		t_pointer const K = mesh->get_triangle(k);

		unsigned const k_index = K->index;



		for (unsigned r = 0; r < 3; r++)
			for (unsigned s = 0; s < 3; s++)
				block(r, s) = δij(r, s) - coefficient * σ(k_index, 0, r, s);

		Eigen::Matrix3d const inverse_block = block.inverse();


		for (unsigned m = 0; m < 3; m++) {
			for (unsigned Ei = 0; Ei < 3; Ei++) {

				block1(m, Ei) = λ(k_index, 0, m, Ei);
				block2(m, Ei) = λ(k_index, 1, m, Ei);

			}
		}

		block1 *= -coefficient;
		block2 *= -coefficient;



		Eigen::Matrix3d const iDH1block = inverse_block * block1;
		Eigen::Matrix3d const iDH2block = inverse_block * block2;




		for (unsigned i = 0; i < 3; i++) {

			iDH1_matrix.setCoeff(k_index, i, 0) = iDH1block(i, 0);
			iDH1_matrix.setCoeff(k_index, i, 1) = iDH1block(i, 1);
			iDH1_matrix.setCoeff(k_index, i, 2) = iDH1block(i, 2);

			iDH2_matrix.setCoeff(k_index, i, 0) = iDH2block(i, 0);
			iDH2_matrix.setCoeff(k_index, i, 1) = iDH2block(i, 1);
			iDH2_matrix.setCoeff(k_index, i, 2) = iDH2block(i, 2);

		}



		for (unsigned Ei = 0; Ei < 3; Ei++) {

			e_pointer const E = K->edges[Ei];
			unsigned const e_index = E->index;

			unsigned const dof0 = LI(K, E, 0);
			unsigned const dof1 = LI(K, E, 1);


			for (unsigned m = 0; m < 3; m++) {


				real AB1 = 0.0;
				real AB2 = 0.0;

				for (unsigned i = 0; i < 8; i++) {

					AB1 += α(k_index, 0, i, dof0)*β(k_index, 0, i, m);
					AB2 += α(k_index, 0, i, dof1)*β(k_index, 0, i, m);

				}

				real const val1 = AB1 / viscosities[k_index];
				real const val2 = AB2 / viscosities[k_index];

				//delete R1, R2 and make it in this format only ?

				R1_matrix.setCoeff(k_index, Ei, m) = abs(val1) < INTEGRAL_PRECISION ? 0.0 : val1;
				R2_matrix.setCoeff(k_index, Ei, m) = abs(val2) < INTEGRAL_PRECISION ? 0.0 : val2;

			}
		}
	}





	R1iDH1.setZero();
	R1iDH2.setZero();
	R2iDH1.setZero();
	R2iDH2.setZero();

	for (unsigned k = 0; k < nk; k++) {



		t_pointer const K = mesh->get_triangle(k);
		unsigned const k_index = K->index;


		for (unsigned ei = 0; ei < 3; ei++) {


			e_pointer const Ei = K->edges[ei];
			unsigned const e_index_i = Ei->index;


			if (Ei->marker == E_MARKER::DIRICHLET) {

				R1iDH1.coeffRef(e_index_i, e_index_i) = 0.0;
				R1iDH2.coeffRef(e_index_i, e_index_i) = 0.0;

				R2iDH1.coeffRef(e_index_i, e_index_i) = 0.0;
				R2iDH2.coeffRef(e_index_i, e_index_i) = 0.0;

				continue;

			}



			for (unsigned ej = 0; ej < 3; ej++) {


				e_pointer const Ej = K->edges[ej];
				unsigned const e_index_j = Ej->index;


				real sum11 = 0.0;
				real sum12 = 0.0;
				real sum21 = 0.0;
				real sum22 = 0.0;

				for (unsigned i = 0; i < 3; i++) {

					sum11 += R1_matrix(k_index, ei, i) * iDH1_matrix(k_index, i, ej);
					sum12 += R1_matrix(k_index, ei, i) * iDH2_matrix(k_index, i, ej);
					sum21 += R2_matrix(k_index, ei, i) * iDH1_matrix(k_index, i, ej);
					sum22 += R2_matrix(k_index, ei, i) * iDH2_matrix(k_index, i, ej);

				}

				R1iDH1.coeffRef(e_index_i, e_index_j) += sum11;
				R1iDH2.coeffRef(e_index_i, e_index_j) += sum12;
				R2iDH1.coeffRef(e_index_i, e_index_j) += sum21;
				R2iDH2.coeffRef(e_index_i, e_index_j) += sum22;

			}
		}
	}

};
void solver::assembleRInverseDH_whole() {


	real const coefficient = θ * dt;

	Eigen::Matrix3d block;
	Eigen::Matrix3d block1;
	Eigen::Matrix3d block2;

	std::vector<Eigen::Triplet<double>> tri;
	//tri.reserve(estimation_of_entries);


	// It is sufficient to zero only diagonal elements. Therefore, there is no need for += in the sequel
	for (unsigned e = 0; e < ne; e++) {

		//R1iDH1.coeffRef(e, e) = 0.0;
		//R1iDH2.coeffRef(e, e) = 0.0;
		//R2iDH1.coeffRef(e, e) = 0.0;
		//R2iDH2.coeffRef(e, e) = 0.0;

		//R1iDH1.coeffRef(e, e) = M_j1_s1.coeff(e, e);
		//R1iDH2.coeffRef(e, e) = M_j1_s2.coeff(e, e);
		//R2iDH1.coeffRef(e, e) = M_j2_s1.coeff(e, e);
		//R2iDH2.coeffRef(e, e) = M_j2_s2.coeff(e, e);


		//internalPressureSystem_smat.coeffRef(e, e) = M_j1_s1.coeff(e, e);
		//internalPressureSystem_smat.coeffRef(e, e + ne) = M_j1_s2.coeff(e, e);
		//internalPressureSystem_smat.coeffRef(e + ne, e) = M_j2_s1.coeff(e, e);
		//internalPressureSystem_smat.coeffRef(e + ne, e + ne) = M_j2_s2.coeff(e, e);


		//internalPressureSystem_smat.coeffRef(e, e) = 0.0;
		//internalPressureSystem_smat.coeffRef(e, e + ne) = 0.0;
		//internalPressureSystem_smat.coeffRef(e + ne, e) = 0.0;
		//internalPressureSystem_smat.coeffRef(e + ne, e + ne) = 0.0;

		real const M11 = M_j1_s1.coeff(e, e);
		real const M12 = M_j1_s2.coeff(e, e);
		real const M21 = M_j2_s1.coeff(e, e);
		real const M22 = M_j2_s2.coeff(e, e);

		Eigen::Triplet<real> const T1(e, e, M11);
		Eigen::Triplet<real> const T2(e, e + ne, M12);
		Eigen::Triplet<real> const T3(e + ne, e, M21);
		Eigen::Triplet<real> const T4(e + ne, e + ne, M22);

		tri.push_back(T1);
		tri.push_back(T2);
		tri.push_back(T3);
		tri.push_back(T4);

	}

	for (unsigned k = 0; k < nk; k++) {



		t_pointer const K = mesh->get_triangle(k);

		unsigned const k_index = K->index;


		for (unsigned r = 0; r < 3; r++)
			for (unsigned s = 0; s < 3; s++)
				block(r, s) = δij(r, s) - coefficient * σ(k_index, 0, r, s);

		Eigen::Matrix3d const inverse_block = block.inverse();


		for (unsigned m = 0; m < 3; m++) {
			for (unsigned Ei = 0; Ei < 3; Ei++) {

				block1(m, Ei) = λ(k_index, 0, m, Ei);
				block2(m, Ei) = λ(k_index, 1, m, Ei);

			}
		}

		block1 *= -coefficient;
		block2 *= -coefficient;


		Eigen::Matrix3d const iDH1block = inverse_block * block1;
		Eigen::Matrix3d const iDH2block = inverse_block * block2;



		for (unsigned ei = 0; ei < 3; ei++) {


			e_pointer const Ei = K->edges[ei];
			unsigned const e_index_i = K->edges[ei]->index;


			// Diagonal elements are already zeroed
			if (Ei->marker == E_MARKER::DIRICHLET)
				continue;


			for (unsigned ej = 0; ej < 3; ej++) {


				unsigned const e_index_j = K->edges[ej]->index;


				real sum11 = 0.0;
				real sum12 = 0.0;
				real sum21 = 0.0;
				real sum22 = 0.0;

				// Number of degrees of freedom of internal pressure
				for (unsigned m = 0; m < 3; m++) {

					sum11 += R1_matrix(k_index, ei, m) * iDH1block(m, ej);
					sum12 += R1_matrix(k_index, ei, m) * iDH2block(m, ej);
					sum21 += R2_matrix(k_index, ei, m) * iDH1block(m, ej);
					sum22 += R2_matrix(k_index, ei, m) * iDH2block(m, ej);

				}

				// Because diagonal elements were zeroed at the beginning, the += operator is needed only here
				if (e_index_i == e_index_j) {

					//R1iDH1.coeffRef(e_index_i, e_index_i) += sum11;
					//R1iDH2.coeffRef(e_index_i, e_index_i) += sum12;
					//R2iDH1.coeffRef(e_index_i, e_index_i) += sum21;
					//R2iDH2.coeffRef(e_index_i, e_index_i) += sum22;


					//internalPressureSystem_smat.coeffRef(e_index_i, e_index_i) += sum11;
					//internalPressureSystem_smat.coeffRef(e_index_i, e_index_i + ne) += sum12;
					//internalPressureSystem_smat.coeffRef(e_index_i + ne, e_index_i) += sum21;
					//internalPressureSystem_smat.coeffRef(e_index_i + ne, e_index_i + ne) += sum22;

					Eigen::Triplet<real> const T1(e_index_i, e_index_i, sum11);
					Eigen::Triplet<real> const T2(e_index_i, e_index_i + ne, sum12);
					Eigen::Triplet<real> const T3(e_index_i + ne, e_index_i, sum21);
					Eigen::Triplet<real> const T4(e_index_i + ne, e_index_i + ne, sum22);

					tri.push_back(T1);
					tri.push_back(T2);
					tri.push_back(T3);
					tri.push_back(T4);

					continue;

				}

				//R1iDH1.coeffRef(e_index_i, e_index_j) = sum11;
				//R1iDH2.coeffRef(e_index_i, e_index_j) = sum12;
				//R2iDH1.coeffRef(e_index_i, e_index_j) = sum21;
				//R2iDH2.coeffRef(e_index_i, e_index_j) = sum22;

				//R1iDH1.coeffRef(e_index_i, e_index_j) = sum11 + M_j1_s1.coeff(e_index_i, e_index_j);
				//R1iDH2.coeffRef(e_index_i, e_index_j) = sum12 + M_j1_s2.coeff(e_index_i, e_index_j);
				//R2iDH1.coeffRef(e_index_i, e_index_j) = sum21 + M_j2_s1.coeff(e_index_i, e_index_j);
				//R2iDH2.coeffRef(e_index_i, e_index_j) = sum22 + M_j2_s2.coeff(e_index_i, e_index_j);

				//internalPressureSystem_smat.coeffRef(e_index_i, e_index_j) = sum11 + M_j1_s1.coeff(e_index_i, e_index_j);
				//internalPressureSystem_smat.coeffRef(e_index_i, e_index_j + ne) = sum12 + M_j1_s2.coeff(e_index_i, e_index_j);
				//internalPressureSystem_smat.coeffRef(e_index_i + ne, e_index_j) = sum21 + M_j2_s1.coeff(e_index_i, e_index_j);
				//internalPressureSystem_smat.coeffRef(e_index_i + ne, e_index_j + ne) = sum22 + M_j2_s2.coeff(e_index_i, e_index_j);


				real const M11 = sum11 + M_j1_s1.coeff(e_index_i, e_index_j);
				real const M12 = sum12 + M_j1_s2.coeff(e_index_i, e_index_j);
				real const M21 = sum21 + M_j2_s1.coeff(e_index_i, e_index_j);
				real const M22 = sum22 + M_j2_s2.coeff(e_index_i, e_index_j);


				Eigen::Triplet<real> const T1(e_index_i, e_index_j, M11);
				Eigen::Triplet<real> const T2(e_index_i, e_index_j + ne, M12);
				Eigen::Triplet<real> const T3(e_index_i + ne, e_index_j, M21);
				Eigen::Triplet<real> const T4(e_index_i + ne, e_index_j + ne, M22);

				tri.push_back(T1);
				tri.push_back(T2);
				tri.push_back(T3);
				tri.push_back(T4);

			}
		}
	}

	internalPressureSystem_smat.setFromTriplets(tri.begin(), tri.end());

};

*/

/*

template<unsigned QuadraturePrecision>
void solver<QuadraturePrecision>::assemblePressureSystemMatrix3() {


	real const coefficient = θ * dt;

	Eigen::Matrix3d block;
	Eigen::Matrix3d block1;
	Eigen::Matrix3d block2;


	// It is sufficient to zero only diagonal elements. Therefore, there is no need for += in the sequel
	for (unsigned e = 0; e < ne; e++) {

		internalPressureSystem.coeffRef(e, e) = M_j1_s1.coeff(e, e);
		internalPressureSystem.coeffRef(e, e + ne) = M_j1_s2.coeff(e, e);
		internalPressureSystem.coeffRef(e + ne, e) = M_j2_s1.coeff(e, e);
		internalPressureSystem.coeffRef(e + ne, e + ne) = M_j2_s2.coeff(e, e);

	}

	for (unsigned k = 0; k < nk; k++) {



		t_pointer const K = mesh->get_triangle(k);

		unsigned const k_index = K->index;




for (unsigned r = 0; r < 3; r++)
	for (unsigned s = 0; s < 3; s++)
		block(r, s) = δij(r, s) - coefficient * σ(k_index, r, s);

Eigen::Matrix3d const inverse_block = block.inverse();




for (unsigned m = 0; m < 3; m++) {
	for (unsigned Ei = 0; Ei < 3; Ei++) {

		block1(m, Ei) = λ(k_index, 0, m, Ei);
		block2(m, Ei) = λ(k_index, 1, m, Ei);

	}
}

block1 *= -coefficient;
block2 *= -coefficient;




Eigen::Matrix3d const iDH1block = inverse_block * block1;
Eigen::Matrix3d const iDH2block = inverse_block * block2;





for (unsigned ei = 0; ei < 3; ei++) {


	e_pointer const Ei = K->edges[ei];
	unsigned const e_index_i = K->edges[ei]->index;


	// Diagonal elements are already zeroed
	if (Ei->marker == E_MARKER::DIRICHLET)
		continue;


	for (unsigned ej = 0; ej < 3; ej++) {


		unsigned const e_index_j = K->edges[ej]->index;


		real sum11 = 0.0;
		real sum12 = 0.0;
		real sum21 = 0.0;
		real sum22 = 0.0;

		// Number of degrees of freedom of internal pressure
		for (unsigned m = 0; m < 3; m++) {

			sum11 += R1_block(k_index, ei, m) * iDH1block(m, ej);
			sum12 += R1_block(k_index, ei, m) * iDH2block(m, ej);
			sum21 += R2_block(k_index, ei, m) * iDH1block(m, ej);
			sum22 += R2_block(k_index, ei, m) * iDH2block(m, ej);

		}

		// Because diagonal elements were zeroed at the beginning, the += operator is needed only here
		if (e_index_i == e_index_j) {

			internalPressureSystem.coeffRef(e_index_i, e_index_i) += sum11;
			internalPressureSystem.coeffRef(e_index_i, e_index_i + ne) += sum12;
			internalPressureSystem.coeffRef(e_index_i + ne, e_index_i) += sum21;
			internalPressureSystem.coeffRef(e_index_i + ne, e_index_i + ne) += sum22;

			continue;

		}

		internalPressureSystem.coeffRef(e_index_i, e_index_j) = sum11 + M_j1_s1.coeff(e_index_i, e_index_j);
		internalPressureSystem.coeffRef(e_index_i, e_index_j + ne) = sum12 + M_j1_s2.coeff(e_index_i, e_index_j);
		internalPressureSystem.coeffRef(e_index_i + ne, e_index_j) = sum21 + M_j2_s1.coeff(e_index_i, e_index_j);
		internalPressureSystem.coeffRef(e_index_i + ne, e_index_j + ne) = sum22 + M_j2_s2.coeff(e_index_i, e_index_j);

	}
}
	}

};
template<unsigned QuadraturePrecision>
void solver<QuadraturePrecision>::assemblePressureSystemMatrix4() {


	real const coefficient = θ * dt;

	Eigen::Matrix3d block;
	Eigen::Matrix3d block1;
	Eigen::Matrix3d block2;

	std::vector<Eigen::Triplet<real>> tri;


	// It is sufficient to zero only diagonal elements. Therefore, there is no need for += in the sequel
	for (unsigned e = 0; e < ne; e++) {

		real const M11 = M_j1_s1.coeff(e, e);
		real const M12 = M_j1_s2.coeff(e, e);
		real const M21 = M_j2_s1.coeff(e, e);
		real const M22 = M_j2_s2.coeff(e, e);

		Eigen::Triplet<real> const T1(e, e, M11);
		Eigen::Triplet<real> const T2(e, e + ne, M12);
		Eigen::Triplet<real> const T3(e + ne, e, M21);
		Eigen::Triplet<real> const T4(e + ne, e + ne, M22);

		tri.push_back(T1);
		tri.push_back(T2);
		tri.push_back(T3);
		tri.push_back(T4);

	}

	for (unsigned k = 0; k < nk; k++) {



		t_pointer const K = mesh->get_triangle(k);

		unsigned const k_index = K->index;



		for (unsigned r = 0; r < 3; r++)
			for (unsigned s = 0; s < 3; s++)
				block(r, s) = δij(r, s) - coefficient * σ(k_index, r, s);

		Eigen::Matrix3d const inverse_block = block.inverse();


		for (unsigned m = 0; m < 3; m++) {
			for (unsigned Ei = 0; Ei < 3; Ei++) {

				block1(m, Ei) = λ(k_index, 0, m, Ei);
				block2(m, Ei) = λ(k_index, 1, m, Ei);

			}
		}

		block1 *= -coefficient;
		block2 *= -coefficient;



		Eigen::Matrix3d const iDH1block = inverse_block * block1;
		Eigen::Matrix3d const iDH2block = inverse_block * block2;




		for (unsigned ei = 0; ei < 3; ei++) {


			e_pointer const Ei = K->edges[ei];
			unsigned const e_index_i = K->edges[ei]->index;


			// Diagonal elements are already zeroed/or there is Mij already
			if (Ei->marker == E_MARKER::DIRICHLET)
				continue;


			for (unsigned ej = 0; ej < 3; ej++) {


				unsigned const e_index_j = K->edges[ej]->index;


				real sum11 = 0.0;
				real sum12 = 0.0;
				real sum21 = 0.0;
				real sum22 = 0.0;

				// Number of degrees of freedom of internal pressure
				for (unsigned m = 0; m < 3; m++) {

					sum11 += R1_block(k_index, ei, m) * iDH1block(m, ej);
					sum12 += R1_block(k_index, ei, m) * iDH2block(m, ej);
					sum21 += R2_block(k_index, ei, m) * iDH1block(m, ej);
					sum22 += R2_block(k_index, ei, m) * iDH2block(m, ej);

				}

				// Because diagonal elements were zeroed at the beginning, the += operator is needed only here
				if (e_index_i == e_index_j) {

					Eigen::Triplet<real> const T1(e_index_i, e_index_i, sum11);
					Eigen::Triplet<real> const T2(e_index_i, e_index_i + ne, sum12);
					Eigen::Triplet<real> const T3(e_index_i + ne, e_index_i, sum21);
					Eigen::Triplet<real> const T4(e_index_i + ne, e_index_i + ne, sum22);

					tri.push_back(T1);
					tri.push_back(T2);
					tri.push_back(T3);
					tri.push_back(T4);

					continue;

				}

				real const M11 = sum11 + M_j1_s1.coeff(e_index_i, e_index_j);
				real const M12 = sum12 + M_j1_s2.coeff(e_index_i, e_index_j);
				real const M21 = sum21 + M_j2_s1.coeff(e_index_i, e_index_j);
				real const M22 = sum22 + M_j2_s2.coeff(e_index_i, e_index_j);

				Eigen::Triplet<real> const T1(e_index_i, e_index_j, M11);
				Eigen::Triplet<real> const T2(e_index_i, e_index_j + ne, M12);
				Eigen::Triplet<real> const T3(e_index_i + ne, e_index_j, M21);
				Eigen::Triplet<real> const T4(e_index_i + ne, e_index_j + ne, M22);

				tri.push_back(T1);
				tri.push_back(T2);
				tri.push_back(T3);
				tri.push_back(T4);

			}
		}
	}

	internalPressureSystem.setFromTriplets(tri.begin(), tri.end());

};


*/

//pressure
/*
void solver::computePressureEquation() {


	assembleInverseD();
	assembleH();

	assembleG();
	assembleV();


	//internalPressureSystem.block(0, 0, ne, ne) = R1 * iD * H1 + M_j1_s1;
	//internalPressureSystem.block(0, ne, ne, ne) = R1 * iD * H2 + M_j1_s2;
	//internalPressureSystem.block(ne, 0, ne, ne) = R2 * iD * H1 + M_j2_s1;
	//internalPressureSystem.block(ne, ne, ne, ne) = R2 * iD * H2 + M_j2_s2;

	//pressureSystemRhs.head(ne) = R1 * iD * G - V1;
	//pressureSystemRhs.tail(ne) = R2 * iD * G - V2;

	//denseLUsolver_InternalPressureSystem.compute(internalPressureSystem);
	//vec const solution = denseLUsolver_InternalPressureSystem.solve(pressureSystemRhs);

	//Tp1 = solution.head(ne);
	//Tp2 = solution.tail(ne);



	
	//Eigen::MatrixXd const R1iD = R1 * iD;
	//Eigen::MatrixXd const R2iD = R2 * iD;

	//Eigen::MatrixXd const R1iDH1 = R1iD * H1;
	//Eigen::MatrixXd const R1iDH2 = R1iD * H2;
	//Eigen::MatrixXd const R2iDH1 = R2iD * H1;

	//Eigen::MatrixXd R2iDH1M21 = R2iDH1 + M_j2_s1;
	//Eigen::MatrixXd R1iDH2M12 = R1iDH2 + M_j1_s2;

	//Eigen::MatrixXd const R1iDH1 = R1iD * H1;
	//Eigen::MatrixXd const R1iDH2 = R1iD * H2;
	//Eigen::MatrixXd const R2iDH1 = R2iD * H1;

	//Eigen::MatrixXd const R2iDH1M21 = R2iDH1 + M_j2_s1;
	//Eigen::MatrixXd const R1iDH2M12 = R1iDH2 + M_j1_s2;
	//Eigen::MatrixXd inverseR1iDH1M11 = (R1iDH1 + M_j1_s1);

	//inverseR1iDH1M11 = inverseR1iDH1M11.inverse();

	//Eigen::MatrixXd const System1 = (R2iD * H2 + M_j2_s2) - R2iDH1M21 * inverseR1iDH1M11 * R1iDH2M12;

	//Eigen::VectorXd const R1iDGV1 = R1iD * G - V1;
	//Eigen::VectorXd const RHS1 = (R2iD * G - V2) - R2iDH1M21 * inverseR1iDH1M11 * R1iDGV1;

	//denseLUsolver_InternalPressureSystem.compute(System1);
	//Tp2 = denseLUsolver_InternalPressureSystem.solve(RHS1);
	//Tp1 = inverseR1iDH1M11 * (R1iDGV1 - R1iDH2M12 * Tp2);

	//π_eigen = iD * (G - H1 * Tp1 - H2 * Tp2);
	


	//Eigen::VectorXd RHS1a;
	//Eigen::VectorXd R1iDGV1;

	//Eigen::MatrixXd R1iDH2;

	//Eigen::MatrixXd R2iDH1M21;
	//Eigen::MatrixXd R1iDH2M12;

	////smat System1a;
	////smat System1b;

	//Eigen::MatrixXd System1a;
	//Eigen::MatrixXd System1b;

	//Eigen::MatrixXd inverseR1iDH1M11;

	//Eigen::MatrixXd const R1iD = R1 * iD;
	//Eigen::MatrixXd const R2iD = R2 * iD;

	//Eigen::MatrixXd const R2iDH1 = R2iD * H1;
	//Eigen::MatrixXd const R1iDH1 = R1iD * H1;




	//omp_set_num_threads(4);
	//#pragma omp parallel
	//{
	//#pragma omp sections
	//	{
	//#pragma omp section
	//		{
	//			inverseR1iDH1M11 = (R1iDH1 + M_j1_s1);
	//			inverseR1iDH1M11 = inverseR1iDH1M11.inverse();
	//		}
	//#pragma omp section
	//		{
	//			R1iDH2 = R1iD * H2;
	//			R1iDH2M12 = R1iDH2 + M_j1_s2;
	//		}
	//#pragma omp section
	//		{
	//			R2iDH1M21 = R2iDH1 + M_j2_s1;
	//			R1iDGV1 = R1iD * G - V1;
	//		}
	//#pragma omp section
	//		{
	//			System1a = (R2iD * H2 + M_j2_s2);
	//			RHS1a = (R2iD * G - V2);
	//		}
	//	}
	//}


	//System1b = R2iDH1M21 * inverseR1iDH1M11 * R1iDH2M12;



	//Eigen::MatrixXd const System1 = System1a - System1b;
	//Eigen::VectorXd const RHS1 = RHS1a - R2iDH1M21 * inverseR1iDH1M11 * R1iDGV1;


	////smat const System1 = System1a - System1b;
	////Eigen::SparseLU<smat> solver;
	////solver.compute(System1);
	////Tp2 = solver.solve(RHS1);
	////Tp1 = inverseR1iDH1M11 * (R1iDGV1 - R1iDH2M12 * Tp2);

	//denseLUsolver_InternalPressureSystem.compute(System1);
	//Tp2 = denseLUsolver_InternalPressureSystem.solve(RHS1);
	//Tp1 = inverseR1iDH1M11 * (R1iDGV1 - R1iDH2M12 * Tp2);




	//smat R1iDH1M11;
	//smat R1iDH2M12;
	//smat R2iDH1M21;
	//smat R2iDH2M22;

smat const R1iD = R1 * iD;
smat const R2iD = R2 * iD;


//omp_set_num_threads(4);
//#pragma omp parallel 
//{
//	#pragma omp sections
//	{
//		#pragma omp section
//		{
//			R1iDH1M11 = R1iD * H1 + M_j1_s1;
//		}
//		#pragma omp section
//		{
//			R1iDH2M12 = R1iD * H2 + M_j1_s2;
//		}
//		#pragma omp section
//		{
//			R2iDH1M21 = R2iD * H2 + M_j2_s2;
//		}
//		#pragma omp section
//		{
//			R2iDH2M22 = R2iD * H2 + M_j2_s2;
//		}
//	}
//}

//omp_set_num_threads(2);
//#pragma omp parallel 
//{
//	#pragma omp sections
//		{
//		#pragma omp section
//		{
//			R1iDH1M11 = R1iD * H1 + M_j1_s1;
//		}
//		#pragma omp section
//		{
//			R1iDH2M12 = R1iD * H2 + M_j1_s2;
//		}
//	}
//}

//omp_set_num_threads(2);
//#pragma omp parallel 
//{
//	#pragma omp sections
//	{
//		#pragma omp section
//		{
//			R2iDH1M21 = R2iD * H2 + M_j2_s2;
//		}
//	#pragma omp section
//		{
//			R2iDH2M22 = R2iD * H2 + M_j2_s2;
//		}
//	}
//}

//R1iDH1M11 = R1iD * H1 + M_j1_s1;
//R1iDH2M12 = R1iD * H2 + M_j1_s2;
//R2iDH1M21 = R2iD * H1 + M_j2_s1;
//R2iDH2M22 = R2iD * H2 + M_j2_s2;




//smat _R1iDH1M11;
//smat _R1iDH2M12;
//smat _R2iDH1M21;
//smat _R2iDH2M22;

//assembleInverseDH();
assembleRInverseDH_whole();

//assembleRInverseDH();
//assembleRInverseDH_whole_immidiate();


//std::cout << (R1 * _iDH1).toDense() << std::endl << std::endl << std::endl;
//std::cout << (R1iDH1).toDense() << std::endl << std::endl << std::endl << std::endl << std::endl;

//std::cout << (R1 * _iDH2).toDense() << std::endl << std::endl << std::endl;
//std::cout << (R1iDH2).toDense() << std::endl << std::endl << std::endl << std::endl << std::endl;

//std::cout << (R2 * _iDH1).toDense() << std::endl << std::endl << std::endl;
//std::cout << (R2iDH1).toDense() << std::endl << std::endl << std::endl << std::endl << std::endl;

//std::cout << (R2 * _iDH2).toDense() << std::endl << std::endl << std::endl;
//std::cout << (R2iDH2).toDense() << std::endl << std::endl << std::endl << std::endl << std::endl;

//std::cout << (R1 * _iDH1).toDense() - (R1iDH1).toDense() << std::endl << std::endl << std::endl;
//std::cout << (R1 * _iDH2).toDense() - (R1iDH2).toDense() << std::endl << std::endl << std::endl;
//std::cout << (R2 * _iDH1).toDense() - (R2iDH1).toDense() << std::endl << std::endl << std::endl;
//std::cout << (R2 * _iDH2).toDense() - (R2iDH2).toDense() << std::endl << std::endl << std::endl;


//omp_set_num_threads(2);
//#pragma omp parallel 
//{
//	#pragma omp sections
//	{
//		#pragma omp section
//		{
//			_R1iDH1M11 = R1 * _iDH1 + M_j1_s1;
//			_R1iDH2M12 = R1 * _iDH2 + M_j1_s2;
//		}
//		#pragma omp section
//		{
//			_R2iDH1M21 = R2 * _iDH1 + M_j2_s1;
//			_R2iDH2M22 = R2 * _iDH2 + M_j2_s2;
//		}
//	}
//}

//smat const _R1iDH1M11 = R1 * _iDH1 + M_j1_s1;
//smat const _R1iDH2M12 = R1 * _iDH2 + M_j1_s2;
//smat const _R2iDH1M21 = R2 * _iDH1 + M_j2_s1;
//smat const _R2iDH2M22 = R2 * _iDH2 + M_j2_s2;

//smat const _R1iDH1M11 = R1iDH1 + M_j1_s1;
//smat const _R1iDH2M12 = R1iDH2 + M_j1_s2;
//smat const _R2iDH1M21 = R2iDH1 + M_j2_s1;
//smat const _R2iDH2M22 = R2iDH2 + M_j2_s2;

//smat const _R1iDH1M11 = R1iDH1;
//smat const _R1iDH2M12 = R1iDH2;
//smat const _R2iDH1M21 = R2iDH1;
//smat const _R2iDH2M22 = R2iDH2;



//smat _R1iDH1M11;
//smat _R1iDH2M12;
//smat _R2iDH1M21;
//smat _R2iDH2M22;

//	omp_set_num_threads(4);
//#pragma omp parallel 
//{
//	#pragma omp sections
//	{
//		#pragma omp section
//		{
//			_R1iDH1M11 = R1iDH1 + M_j1_s1;
//		}
//		#pragma omp section
//		{
//			_R1iDH2M12 = R1iDH2 + M_j1_s2;
//		}
//		#pragma omp section
//		{
//			_R2iDH1M21 = R2iDH1 + M_j2_s1;
//		}
//		#pragma omp section
//		{
//			_R2iDH2M22 = R2iDH2 + M_j2_s2;
//		}
//	}
//}


//Turn this into function assembleRInverseDH_whole -> it could be much faster. But do not delete previous versions -> will be needed for linux (lapack etc.)
//std::vector<Eigen::Triplet<double>> tri;

//for (int i = 0; i < _R1iDH1M11.outerSize(); i++)
//	for (Eigen::SparseMatrix<double>::InnerIterator it(_R1iDH1M11, i); it; ++it)
//		tri.emplace_back(it.row() + 0, it.col() + 0, it.value());

//for (int i = 0; i < _R1iDH2M12.outerSize(); i++)
//	for (Eigen::SparseMatrix<double>::InnerIterator it(_R1iDH2M12, i); it; ++it)
//		tri.emplace_back(it.row() + 0, it.col() + ne, it.value());

//for (int i = 0; i < _R2iDH1M21.outerSize(); i++)
//	for (Eigen::SparseMatrix<double>::InnerIterator it(_R2iDH1M21, i); it; ++it)
//		tri.emplace_back(it.row() + ne, it.col() + 0, it.value());

//for (int i = 0; i < _R2iDH2M22.outerSize(); i++)
//	for (Eigen::SparseMatrix<double>::InnerIterator it(_R2iDH2M22, i); it; ++it)
//		tri.emplace_back(it.row() + ne, it.col() + ne, it.value());


//internalPressureSystem_smat.setFromTriplets(tri.begin(), tri.end());



//	for (int i = 0; i < M_j1_s1.outerSize(); i++)
//	for (Eigen::SparseMatrix<double>::InnerIterator it(M_j1_s1, i); it; ++it)
//		internalPressureSystem_smat.coeffRef(it.row() + 0, it.col() + 0) += it.value();

//for (int i = 0; i < M_j1_s2.outerSize(); i++)
//	for (Eigen::SparseMatrix<double>::InnerIterator it(M_j1_s2, i); it; ++it)
//		internalPressureSystem_smat.coeffRef(it.row() + 0, it.col() + ne) += it.value();

//for (int i = 0; i < M_j2_s1.outerSize(); i++)
//	for (Eigen::SparseMatrix<double>::InnerIterator it(M_j2_s1, i); it; ++it)
//		internalPressureSystem_smat.coeffRef(it.row() + ne, it.col() + 0) += it.value();

//for (int i = 0; i < M_j2_s2.outerSize(); i++)
//	for (Eigen::SparseMatrix<double>::InnerIterator it(M_j2_s2, i); it; ++it)
//		internalPressureSystem_smat.coeffRef(it.row() + ne, it.col() + ne) += it.value();


Eigen::SparseLU<smat> solver;

//omp_set_num_threads(2);
//#pragma omp parallel 
//{
//	#pragma omp sections
//	{
//		#pragma omp section
//		{
//			pressureSystemRhs.head(ne) = R1iD * G - V1;
//			pressureSystemRhs.tail(ne) = R2iD * G - V2;
//		}
//		#pragma omp section
//		{
//			solver.compute(internalPressureSystem_smat);
//		}
//	}
//}

pressureSystemRhs.head(ne) = R1iD * G - V1;
pressureSystemRhs.tail(ne) = R2iD * G - V2;

solver.compute(internalPressureSystem_smat);

vec const solution = solver.solve(pressureSystemRhs);

Tp1 = solution.head(ne);
Tp2 = solution.tail(ne);






// Internal Pressures
π_eigen = iD * (G - H1 * Tp1 - H2 * Tp2);

// Copy Internal Pressure eigen-solution to each element
for (unsigned k = 0; k < nk; k++) {

	unsigned const k_index = mesh->get_triangle(k)->index;

	for (unsigned m = 0; m < 3; m++)
		π.setCoeff(k_index, 0, m, 0) = π_eigen[3 * k_index + m];

}

// Copy Trace Pressure solution to each element
for (unsigned e = 0; e < ne; e++) {


	e_pointer const E = mesh->get_edge(e);

	real const tpval1 = Tp1[E->index];
	real const tpval2 = Tp2[E->index];

	for (unsigned neighbor = 0; neighbor < 2; neighbor++) {

		t_pointer const K = E->neighbors[neighbor];

		if (!K)
			continue;

		unsigned const k_index = K->index;
		unsigned const e_index_local = K->get_edge_index(E);

		tπ.setCoeff(k_index, 0, e_index_local, 0) = tpval1;
		tπ.setCoeff(k_index, 0, e_index_local, 1) = tpval2;

	}
}

assemble_λ();
assemble_σ();

//// This is needed for the next time level
//for (unsigned k = 0; k < nk; k++) {

//	unsigned const k_index = mesh->get_triangle(k)->index;

//	for (unsigned m = 0; m < 3; m++) {

//		real val1 = 0.0;
//		real val2 = 0.0;

//		for (unsigned j = 0; j < 3; j++)
//			val1 += σ(k_index, 0, m, j) * π(k_index, 0, j, 0);

//		for (unsigned El = 0; El < 3; El++)
//			for (unsigned s = 0; s < 2; s++)
//				val2 += λ(k_index, s, m, El) * tπ(k_index, 0, El, s);

//		rkFp.setCoeff(k_index, 0, m, 0) = val1 + val2;

//	}
//}


};

*/



