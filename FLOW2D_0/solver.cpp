
#include "solver.h"
#include "integration.h"
#include "mesh.h"
#include "misc.h"
#include "matrix.h"

#include <iostream>
#include <fstream>
#include <iomanip>

#include <float.h>
#include <omp.h>

#include<thread>



typedef double real;
double const INTEGRAL_PRECISION = 1e-11;


std::vector<Eigen::Triplet<double>> to_triplets(Eigen::SparseMatrix<double> & M, unsigned const ii, unsigned const jj) {

	std::vector<Eigen::Triplet<double>> v;

	for (int i = 0; i < M.outerSize(); i++)
		for (Eigen::SparseMatrix<double>::InnerIterator it(M, i); it; ++it)
			v.emplace_back(it.row() + ii, it.col() + jj, it.value());

	return v;
};


solver::solver(Mesh & m, unsigned nt_0, double dt_0) : nk(m.get_number_of_triangles()), ne(m.get_number_of_edges()) {


	mesh = &m;

	nt = nt_0;
	dt = dt_0;


	/*****************************************************************************/
	/*                                                                           */
	/*    - Memory allocation for the physical quantities   			         */
	/*                                                                           */
	/*****************************************************************************/

	viscosities = new double[nk];
	porosities = new double[nk];


	/*****************************************************************************/
	/*                                                                           */
	/*    - Memory allocation for the beta coefficient. This coefficient         */
	/*      links model equations and EOS                                        */
	/*                                                                           */
	/*****************************************************************************/

	betas = new double[nk];
	betas_prev = new double[nk];


	/*****************************************************************************/
	/*                                                                           */
	/*    - Memory allocation for the unknowns:   							     */
	/*      Concentration, Pressure, Trace pressure, Velocity field              */
	/*                                                                           */
	/*****************************************************************************/

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


	/*****************************************************************************/
	/*                                                                           */
	/*    - Memory allocation for the auxillary variables used in			     */
	/*      time discretization                                                  */
	/*    - 'n'    = current time level                                          */
	/*      'prev' = previous iteration in the inner loop                        */
	/*                                                                           */
	/*****************************************************************************/

	π_n.setNumberOfElements(nk);
	π_prev.setNumberOfElements(nk);

	ξ_n.setNumberOfElements(nk);
	ξ_prev.setNumberOfElements(nk);

	rkFp.setNumberOfElements(nk);
	rkFp_n.setNumberOfElements(nk);

	rkFc.setNumberOfElements(nk);
	rkFc_n.setNumberOfElements(nk);


	/*****************************************************************************/
	/*                                                                           */
	/*    - Memory allocation for the integral coefficients				         */
	/*                                                                           */
	/*****************************************************************************/

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




	/*****************************************************************************/
	/*                                                                           */
	/*    - Memory allocation for the pressure system					         */
	/*                                                                           */
	/*****************************************************************************/

	internalPressureSystem.resize(2 * ne, 2 * ne);
	pressureSystemRhs.resize(2 * ne);

	// Inverse matrix of the matrix for Internal Pressures
	iD.resize(3 * nk, 3 * nk);

	// Matrices for Trace Pressures: H1 = mean pressure (dof 1), H2 = linear pressure (dof 2)
	H1.resize(3 * nk, ne);
	H2.resize(3 * nk, ne);


	R1iDH1.resize(ne, ne);
	R1iDH2.resize(ne, ne);
	R2iDH1.resize(ne, ne);
	R2iDH2.resize(ne, ne);

	R1iDH1M11.resize(ne, ne);
	R1iDH2M12.resize(ne, ne);
	R2iDH1M21.resize(ne, ne);
	R2iDH2M22.resize(ne, ne);

	iDH1_matrix.setNumberOfElements(nk);
	iDH2_matrix.setNumberOfElements(nk);

	R1_matrix.setNumberOfElements(nk);
	R2_matrix.setNumberOfElements(nk);

	// Right-hand side of the Pressure equation
	G.resize(3 * nk);

	// Resulting Internal Pressures: 1. Mean Pressure, 2. Linear Pressure in x-direction, 3. Linear Pressure in y-direction
	π_eigen.resize(3 * nk);


	/*****************************************************************************/
	/*                                                                           */
	/*    - Memory allocation for the Trace Pressure system				         */
	/*                                                                           */
	/*****************************************************************************/

	smat tracePressureSystem_LU(2 * ne, 2 * ne);
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

	   
	/*****************************************************************************/
	/*                                                                           */
	/*    - Initilize initial pressure/concentration condition				     */
	/*      and physical quantities: viscosity, porosity, etc.                   */
	/*                                                                           */
	/*****************************************************************************/

	initializeValues();


	/*****************************************************************************/
	/*                                                                           */
	/*    - Initilize integral coefficient matrices							     */
	/*                                                                           */
	/*****************************************************************************/

	assemble_α();
	assemble_β();
	assemble_χ();
	assemble_η();
	assemble_τ();


	/*****************************************************************************/
	/*                                                                           */
	/*    - Initilize constant matrices										     */
	/*                                                                           */
	/*****************************************************************************/

	assembleR();
	assembleM();


	/*****************************************************************************/
	/*                                                                           */
	/*    - Assembly of the trace pressure system and its LU decomposition       */
	/*                                                                           */
	/*****************************************************************************/

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








	for (unsigned k = 0; k < nk; k++) {


		t_pointer const K = mesh->get_triangle(k);
		unsigned const k_index = K->index;


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





	/****************************************************************/
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

		real const detJF = (real) abs((x1 - x0)*(y2 - y0) - (x2 - x0)*(y1 - y0));

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


				/**/
				Eigen::Vector2d const physicalNormal = referenceNormal / denominator;

				for (unsigned j = 0; j < 8; j++) {

					real const dotNormal = physicalNormal.dot(basisRaviartThomas.col(j));

					physicalNormalDotPhysicalRaviartThomasBasis_quadPoints[j + 8 * (n + number_of_quadrature_points * (El + 3 * k_index))] = dotNormal;

				}

				/**/


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



void solver::initializeValues() {


	double const time = nt * dt;


	for (unsigned k = 0; k < nk; k++) {


		t_pointer const K = mesh->get_triangle(k);
		unsigned const k_index = K->index;


		real const area = K->area();


		viscosities[k_index] = integrate_triangle(K, viscosity) / area;
		porosities[k_index] = integrate_triangle(K, porosity) / area;


		π.setCoeff(k_index, 0, 0, 0) = integrate_triangle(K, time, barenblatt) / area;
		π.setCoeff(k_index, 0, 1, 0) = 0.0;
		π.setCoeff(k_index, 0, 2, 0) = 0.0;


		ξ.setCoeff(k_index, 0, 0, 0) = integrate_triangle(K, time, barenblatt) / area;
		ξ.setCoeff(k_index, 0, 1, 0) = 0.0;
		ξ.setCoeff(k_index, 0, 2, 0) = 0.0;

		sources.setCoeff(k_index, 0, 0, 0) = F1(K, time);
		sources.setCoeff(k_index, 0, 1, 0) = F2(K, time);
		sources.setCoeff(k_index, 0, 2, 0) = F3(K, time);

		rkFc.setCoeff(k_index, 0, 0, 0) = 0.0;
		rkFc.setCoeff(k_index, 0, 1, 0) = 0.0;
		rkFc.setCoeff(k_index, 0, 2, 0) = 0.0;

		rkFc_n.setCoeff(k_index, 0, 0, 0) = rkFc(k_index, 0, 0, 0);
		rkFc_n.setCoeff(k_index, 0, 1, 0) = rkFc(k_index, 0, 1, 0);
		rkFc_n.setCoeff(k_index, 0, 2, 0) = rkFc(k_index, 0, 2, 0);


		rkFp.setCoeff(k_index, 0, 0, 0) = 0.0;
		rkFp.setCoeff(k_index, 0, 1, 0) = 0.0;
		rkFp.setCoeff(k_index, 0, 2, 0) = 0.0;

		rkFp_n.setCoeff(k_index, 0, 0, 0) = rkFp(k_index, 0, 0, 0);
		rkFp_n.setCoeff(k_index, 0, 1, 0) = rkFp(k_index, 0, 1, 0);
		rkFp_n.setCoeff(k_index, 0, 2, 0) = rkFp(k_index, 0, 2, 0);


		v_pointer const a = K->vertices[0];
		v_pointer const b = K->vertices[1];
		v_pointer const c = K->vertices[2];

		real const x0 = (real)a->x;
		real const y0 = (real)a->y;

		real const x1 = (real)b->x;
		real const y1 = (real)b->y;

		real const x2 = (real)c->x;
		real const y2 = (real)c->y;

		affineMappingMatrixJF.setCoeff(k_index, 0, 0, 0) = x1 - x0;
		affineMappingMatrixJF.setCoeff(k_index, 0, 0, 1) = x2 - x0;
		affineMappingMatrixJF.setCoeff(k_index, 0, 1, 0) = y1 - y0;
		affineMappingMatrixJF.setCoeff(k_index, 0, 1, 1) = y2 - y0;

		affineMappingMatrixDeterminant[k_index] = (x1 - x0)*(y2 - y0) - (x2 - x0)*(y1 - y0);

	}

};
void solver::computeBetas() {


	t_pointer K = NULL;

	for (unsigned k = 0; k < nk; k++) {

		//K = mesh->getElement(k);

		betas[k] = 1.0;// Deos(ξ(K, 0, 0));

	}
	//	betas[i] = (eos(concentrations[i]) - eos(concentrations_temporary[i])) / (concentrations[i] - concentrations_temporary[i]);

};

bool solver::stopCriterion() {


	double sP1 = 0.0;
	double sP2 = 0.0;

	double sC1 = 0.0;
	double sC2 = 0.0;

	double sB1 = 0.0;
	double sB2 = 0.0;

	t_pointer K = NULL;


	for (unsigned k = 0; k < nk; k++) {

		unsigned const k_index = mesh->get_triangle(k)->index;

		sP1 += sqr(π(k_index, 0, 0, 0) - π_prev(k_index, 0, 0, 0));
		sP2 += sqr(π(k_index, 0, 0, 0));

	}

	double const val_P = sP1 / sP2;

	if (val_P > TOL)
		return false;


	for (unsigned k = 0; k < nk; k++) {

		unsigned const k_index = mesh->get_triangle(k)->index;

		sC1 += sqr(ξ(k_index, 0, 0, 0) - ξ_prev(k_index, 0, 0, 0));
		sC2 += sqr(ξ(k_index, 0, 0, 0));

	}

	double const val_C = sC1 / sC2;

	if (val_C > TOL)
		return false;


	for (unsigned k = 0; k < nk; k++) {

		sB1 += sqr(betas[k] - betas_prev[k]);
		sB2 += sqr(betas[k]);

	}

	double const val_B = sB1 / sB2;

	if (val_B > TOL)
		return false;

	return true;


	/*
	for (unsigned k = 0; k < nk; k++) {


		K = mesh->getElement(k);

		double const a = K->nodes[0]->x;
		double const b = K->nodes[1]->x;


		quadrature const quad(error_quadRule, a, b, METHOD::GAUSS);

		double y = 0.0;
		double w = 0.0;

		double val = 0.0;


		for (unsigned j = 0; j < error_quadRule; j++) {


			y = quad.points[j];
			w = quad.weigths[j];

			val = ((π(K, 0, 0)*phi1(y, a, b) + π(K, 0, 1)*phi2(y, a, b)) - (π_prev(K, 0, 0)*phi1(y, a, b) + π_prev(K, 0, 1) * phi2(y, a, b)));

			sP1 += w * sqr(val);
			sP2 += w * sqr(π(K, 0, 0)*phi1(y, a, b) + π(K, 0, 1)*phi2(y, a, b));

		}

	}

	double const val_P = sP1 / sP2;

	if (val_P > TOL)
		return false;


	return true;
	*/

};
void solver::concentrationCorrection() {

	//unsigned const nk = nk;

	//double const eps = sqrt(DBL_EPSILON);

	//double s1 = 0.0;
	//double s2 = 0.0;

	//element * K = NULL;

	//for (unsigned k = 0; k < nk; k++) {

	//	K = mesh->getElement(k);

	//	s1 += sqr(ξ(K, 0, 0) - ξ_n(K, 0, 0));
	//	s2 += sqr(ξ(K, 0, 0));

	//}

	//if (sqrt(s1) / sqrt(s2) < eps) {

	//	for (unsigned k = 0; k < nk; k++) {

	//		std::cout << "Corrected" << std::endl;
	//		K = mesh->getElement(k);

	//		ξ.setCoeff(K, 0, 0) = ξ_n(K, 0, 0) + DBL_EPSILON * ξ(K, 0, 0);

	//	}

	//}

};

void solver::computeTracePressures() {


	// Assembly eigen-vector of internal pressures
	for (unsigned k = 0; k < nk; k++) {

		unsigned const k_index = mesh->get_triangle(k)->index;

		for (unsigned m = 0; m < 3; m++)
			π_eigen[3 * k_index + m] = π(k_index, 0, m, 0);

	}

	assembleV();


	traceSystemRhs.head(ne) = R1 * π_eigen - V1;
	traceSystemRhs.tail(ne) = R2 * π_eigen - V2;

	vec const solution = sparseLUsolver_TracePressureSystem.solve(traceSystemRhs);

	Tp1 = solution.head(ne);
	Tp2 = solution.tail(ne);


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

};
void solver::computePressureEquation() {


	assembleInverseD();
	assembleH();

	assembleG();
	assembleV();

	
	smat const R1iD = R1 * iD;
	smat const R2iD = R2 * iD;


	//assemblePressureSystemMatrix1();
	//assemblePressureSystemMatrix2();
	//assemblePressureSystemMatrix3();
	assemblePressureSystemMatrix4();

	//Eigen::SparseLU<smat> solver;
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

	Eigen::SparseLU<smat> const solver(internalPressureSystem);

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

void solver::computeVelocities() {


	real const time = (nt + 1)* dt;

	for (unsigned k = 0; k < nk; k++) {


		t_pointer const K = mesh->get_triangle(k);
		unsigned const k_index = K->index;


		//CoeffMatrixOfElement<1, 8, 1> velocity_temp(velocities, k_index);
		//CoeffMatrixOfElement<1, 8, 8> alpha_temp(α, k_index);
		//CoeffMatrixOfElement<1, 8, 3> beta_temp(β, k_index);
		//CoeffMatrixOfElement<1, 3, 1> pi_temp(π, k_index);
		//CoeffMatrixOfElement<8, 3, 2> chi_temp(χ, k_index);
		//CoeffMatrixOfElement<1, 3, 2> tpi_temp(tπ, k_index);


		// Loop over element's edges and their degrees of freedom
		for (unsigned El = 0; El < 3; El++) {


			e_pointer const E = K->edges[El];

			if (E->marker == E_MARKER::NEUMANN) {

				velocities.setCoeff(k_index, 0, LI(K, E, 0), 0) = NEUMANN_GAMMA_Q_velocity(E, time);
				velocities.setCoeff(k_index, 0, LI(K, E, 1), 0) = 0.0;

				//velocity_temp.setCoeff(0, LI(K, E, 0), 0) = NEUMANN_GAMMA_Q_velocity(E, time);
				//velocity_temp.setCoeff(0, LI(K, E, 1), 0) = 0.0;

				continue;

			}

			// Loop over two degrees of freedom on each edge
			for (unsigned dof = 0; dof < 2; dof++) {

				unsigned const j = LI(K, E, dof);

				real val1 = 0.0;
				real val2 = 0.0;

				for (unsigned m = 0; m < 8; m++) {

					real PB = 0.0;
					real ECHITP = 0.0;

					real const alpha = α(k_index, 0, m, j);

					for (unsigned i = 0; i < 3; i++)
						PB += β(k_index, 0, m, i) * π(k_index, 0, i, 0);

					val1 += alpha * PB;


					for (unsigned l = 0; l < 3; l++) {

						real CHITP = 0.0;

						for (unsigned s = 0; s < 2; s++)
							CHITP += χ(k_index, m, l, s) * tπ(k_index, 0, l, s);

						ECHITP += CHITP;

					}

					/*for (unsigned i = 0; i < 3; i++)
						PB += beta_temp(0, m, i) * pi_temp(0, i, 0);

					val1 += alpha * PB;


					for (unsigned l = 0; l < 3; l++) {

						real CHITP = 0.0;

						for (unsigned s = 0; s < 2; s++)
							CHITP += chi_temp(m, l, s) * tpi_temp(0, l, s);

						ECHITP += CHITP;

					}*/

					val2 += alpha * ECHITP;

				}

				velocities.setCoeff(k_index, 0, j, 0) = (val1 - val2) / viscosities[k_index];
				//velocity_temp.setCoeff(0, j, 0) = (val1 - val2) / viscosities[k_index];

			}
		}


		// compute Bubble velocities inside element
		for (unsigned dof = 6; dof < 8; dof++) {

			real val1 = 0.0;
			real val2 = 0.0;

			for (unsigned m = 0; m < 8; m++) {


				real PB = 0.0;
				real ECHITP = 0.0;

				real const alpha = α(k_index, 0, m, dof);


				for (unsigned i = 0; i < 3; i++)
					PB += π(k_index, 0, i, 0)*β(k_index, 0, m, i);

				val1 += alpha * PB;

				for (unsigned l = 0; l < 3; l++) {

					real CHITP = 0.0;
					for (unsigned s = 0; s < 2; s++)
						CHITP += χ(k_index, m, l, s) * tπ(k_index, 0, l, s);

					ECHITP += CHITP;

				}

				/*for (unsigned i = 0; i < 3; i++)
					PB += beta_temp(0, m, i) * pi_temp(0, i, 0);

				val1 += alpha * PB;


				for (unsigned l = 0; l < 3; l++) {

					real CHITP = 0.0;

					for (unsigned s = 0; s < 2; s++)
						CHITP += chi_temp(m, l, s) * tpi_temp(0, l, s);

					ECHITP += CHITP;

				}*/

				val2 += alpha * ECHITP;

			}

			velocities.setCoeff(k_index, 0, dof, 0) = (val1 - val2) / viscosities[k_index];
			//velocity_temp.setCoeff(0, dof, 0) = (val1 - val2) / viscosities[k_index];

		}
	}

};
void solver::updateConcentrations_explicit() {


	for (unsigned k = 0; k < nk; k++) {


		unsigned const k_index = mesh->get_triangle(k)->index;


		//CoeffMatrixOfElement<1, 3, 3> eta_temp(η, k_index);
		//CoeffMatrixOfElement<1, 3, 8> gamma_temp(γ, k_index);
		//CoeffMatrixOfElement<1, 8, 1> velocity_temp(velocities, k_index);
		//CoeffMatrixOfElement<1, 3, 1> rkFc_temp(rkFc, k_index);
		//CoeffMatrixOfElement<1, 3, 1> ksi_temp(ξ, k_index);
		//CoeffMatrixOfElement<1, 3, 1> ksi_n_temp(ξ_n, k_index);
		//CoeffMatrixOfElement<1, 3, 1> rkFc_n_temp(rkFc_n, k_index);


		real val0 = 0.0;
		real val1 = 0.0;
		real val2 = 0.0;

		for (unsigned j = 0; j < 8; j++) {


			real etaGamma0 = 0.0;
			real etaGamma1 = 0.0;
			real etaGamma2 = 0.0;

			for (unsigned l = 0; l < 3; l++) {

				etaGamma0 += η(k_index, 0, 0, l)*γ(k_index, 0, l, j);
				etaGamma1 += η(k_index, 0, 1, l)*γ(k_index, 0, l, j);
				etaGamma2 += η(k_index, 0, 2, l)*γ(k_index, 0, l, j);

				//etaGamma0 += eta_temp(0, 0, l)*gamma_temp(0, l, j);
				//etaGamma1 += eta_temp(0, 1, l)*gamma_temp(0, l, j);
				//etaGamma2 += eta_temp(0, 2, l)*gamma_temp(0, l, j);

			}

			val0 += velocities(k_index, 0, j, 0)*etaGamma0;
			val1 += velocities(k_index, 0, j, 0)*etaGamma1;
			val2 += velocities(k_index, 0, j, 0)*etaGamma2;

			//val0 += velocity_temp(0, j, 0)*etaGamma0;
			//val1 += velocity_temp(0, j, 0)*etaGamma1;
			//val2 += velocity_temp(0, j, 0)*etaGamma2;

		}

		rkFc.setCoeff(k_index, 0, 0, 0) = -val0 / porosities[k_index];
		rkFc.setCoeff(k_index, 0, 1, 0) = -val1 / porosities[k_index];
		rkFc.setCoeff(k_index, 0, 2, 0) = -val2 / porosities[k_index];

		//rkFc_temp.setCoeff(0, 0, 0) = -val0 / porosities[k_index];
		//rkFc_temp.setCoeff(0, 1, 0) = -val1 / porosities[k_index];
		//rkFc_temp.setCoeff(0, 2, 0) = -val2 / porosities[k_index];

		// Backward Euler
		//ξ.setCoeff(K, 0, 0) = ξ_n(K, 0, 0) + dt * rkF(K, 0, 0);
		//ξ.setCoeff(K, 0, 1) = ξ_n(K, 0, 1) + dt * rkF(K, 0, 1);

		// Crank-Nicolson
		//ξ.setCoeff(K, 0, 0) = ξ_n(K, 0, 0) + 0.5*dt * (rkF_n(K, 0, 0) + rkF(K, 0, 0));
		//ξ.setCoeff(K, 0, 1) = ξ_n(K, 0, 1) + 0.5*dt * (rkF_n(K, 0, 1) + rkF(K, 0, 1));


		// General θ scheme
		ξ.setCoeff(k_index, 0, 0, 0) = ξ_n(k_index, 0, 0, 0) + dt * (θ * rkFc(k_index, 0, 0, 0) + (1.0 - θ) * rkFc_n(k_index, 0, 0, 0));
		ξ.setCoeff(k_index, 0, 1, 0) = ξ_n(k_index, 0, 1, 0) + dt * (θ * rkFc(k_index, 0, 1, 0) + (1.0 - θ) * rkFc_n(k_index, 0, 1, 0));
		ξ.setCoeff(k_index, 0, 2, 0) = ξ_n(k_index, 0, 2, 0) + dt * (θ * rkFc(k_index, 0, 2, 0) + (1.0 - θ) * rkFc_n(k_index, 0, 2, 0));

		//ksi_temp.setCoeff(0, 0, 0) = ksi_n_temp(0, 0, 0) + dt * (θ * rkFc_temp(0, 0, 0) + (1.0 - θ) * rkFc_n_temp(0, 0, 0));
		//ksi_temp.setCoeff(0, 1, 0) = ksi_n_temp(0, 1, 0) + dt * (θ * rkFc_temp(0, 1, 0) + (1.0 - θ) * rkFc_n_temp(0, 1, 0));
		//ksi_temp.setCoeff(0, 2, 0) = ksi_n_temp(0, 2, 0) + dt * (θ * rkFc_temp(0, 2, 0) + (1.0 - θ) * rkFc_n_temp(0, 2, 0));

	}

};


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

real solver::upwindConcentration(real const s, real const t, t_pointer const & K, Vector<real> const normal, unsigned const El, CoeffMatrix3D<1, 3, 1> & ksi) {

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


	real const time = (nt + 1) * dt;

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

void solver::assembleR() {


	for (unsigned e = 0; e < ne; e++) {


		e_pointer const E = mesh->get_edge(e);
		unsigned const e_index = E->index;
		E_MARKER const e_marker = E->marker;


		if (e_marker == E_MARKER::DIRICHLET) {


			for (unsigned neighborElement = 0; neighborElement < 2; neighborElement++) {

				t_pointer const K = E->neighbors[neighborElement];

				if (!K)
					continue;

				unsigned const k_index = K->index;

				for (unsigned m = 0; m < 3; m++) {

					R1.coeffRef(e_index, 3 * k_index + m) = 0.0;
					R2.coeffRef(e_index, 3 * k_index + m) = 0.0;

				}
			}

			continue;

		}


		for (unsigned neighborElement = 0; neighborElement < 2; neighborElement++) {

			t_pointer const K = E->neighbors[neighborElement];

			if (!K)
				continue;

			unsigned const k_index = K->index;

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

				R1.coeffRef(e_index, 3 * k_index + m) = abs(val1) < INTEGRAL_PRECISION ? 0.0 : val1;
				R2.coeffRef(e_index, 3 * k_index + m) = abs(val2) < INTEGRAL_PRECISION ? 0.0 : val2;

			}
		}
	}

};
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
			M_j1_s2.coeffRef(e_index, e_index) = 0.0;

			M_j2_s1.coeffRef(e_index, e_index) = 0.0;
			M_j2_s2.coeffRef(e_index, e_index) = -1.0;

			//R1iDH1.insert(e_index, e_index) = 0.0;
			//R1iDH2.insert(e_index, e_index) = 0.0;
			//R2iDH1.insert(e_index, e_index) = 0.0;
			//R2iDH2.insert(e_index, e_index) = 0.0;

			continue;

		}


		for (unsigned neighborElement = 0; neighborElement < 2; neighborElement++) {


			t_pointer const K = E->neighbors[neighborElement];

			if (!K)
				continue;


			unsigned const k_index = K->index;

			unsigned const dof0 = LI(K, E, 0);
			unsigned const dof1 = LI(K, E, 1);


			// Loop over edges
			for (unsigned El = 0; El < 3; El++) {


				e_pointer const E_local = K->edges[El];

				unsigned const e_local_index_local = K->get_edge_index(E_local);	// Local index of local edge
				unsigned const e_local_index_global = E_local->index;				// Global index of local edge


				real ACHI_j1_s1 = 0.0;
				real ACHI_j1_s2 = 0.0;
				real ACHI_j2_s1 = 0.0;
				real ACHI_j2_s2 = 0.0;

				for (unsigned i = 0; i < 8; i++) {

					ACHI_j1_s1 += α(k_index, 0, i, dof0) * χ(k_index, i, e_local_index_local, 0);
					ACHI_j1_s2 += α(k_index, 0, i, dof0) * χ(k_index, i, e_local_index_local, 1);

					ACHI_j2_s1 += α(k_index, 0, i, dof1) * χ(k_index, i, e_local_index_local, 0);
					ACHI_j2_s2 += α(k_index, 0, i, dof1) * χ(k_index, i, e_local_index_local, 1);

				}

				real const val1 = ACHI_j1_s1 / viscosities[k_index];
				real const val2 = ACHI_j1_s2 / viscosities[k_index];

				real const val3 = ACHI_j2_s1 / viscosities[k_index];
				real const val4 = ACHI_j2_s2 / viscosities[k_index];

				M_j1_s1.coeffRef(e_index, e_local_index_global) += val1;
				M_j1_s2.coeffRef(e_index, e_local_index_global) += val2;

				M_j2_s1.coeffRef(e_index, e_local_index_global) += val3;
				M_j2_s2.coeffRef(e_index, e_local_index_global) += val4;

				//R1iDH1.insert(e_index, e_local_index_global) = 0.0;
				//R1iDH2.insert(e_index, e_local_index_global) = 0.0;
				//R2iDH1.insert(e_index, e_local_index_global) = 0.0;
				//R2iDH2.insert(e_index, e_local_index_global) = 0.0;

			}
		}
	}

};
void solver::assembleV() {


	real const time = (nt + 1)*dt;

	for (unsigned e = 0; e < ne; e++) {


		e_pointer const E = mesh->get_edge(e);
		unsigned const e_index = E->index;
		E_MARKER const marker = E->marker;


		V1[e_index] = 0.0;
		V2[e_index] = 0.0;


		if (marker == E_MARKER::NEUMANN) {

			real const val = NEUMANN_GAMMA_Q_velocity(E, time);
			V1[e_index] = NEUMANN_GAMMA_Q_velocity(E, time);

		}
		else if (marker == E_MARKER::DIRICHLET) {

			real const val = DIRICHLET_GAMMA_P_pressure(E, time);
			V1[e_index] = DIRICHLET_GAMMA_P_pressure(E, time);

			//real const val = DIRICHLET_GAMMA_P_pressure(E, time);
			//V1[e_index] = DIRICHLET_GAMMA_P_pressure(E, time);

		}
			

	}

};

void solver::assembleInverseD() {


	iD.setZero();

	assemble_σ();

	real const coeff = θ * dt;
	Eigen::Matrix3d block;


	for (unsigned k = 0; k < nk; k++) {


		t_pointer const K = mesh->get_triangle(k);

		unsigned const k_index = K->index;
		unsigned const start_index = 3 * k_index;


		// General θ scheme
		for (unsigned r = 0; r < 3; r++)
			for (unsigned s = 0; s < 3; s++)
				block(r, s) = δij(r, s) - coeff * σ(k_index, 0, r, s);


		Eigen::Matrix3d inverse_block = block.inverse();

		for (unsigned i = 0; i < 3; i++)
			for (unsigned j = 0; j < 3; j++)
				inverse_block.coeffRef(i, j) = abs(inverse_block(i, j)) < INTEGRAL_PRECISION ? 0.0 : inverse_block(i, j);


		for (unsigned i = 0; i < 3; i++)
			for (unsigned j = 0; j < 3; j++)
				iD.coeffRef(start_index + i, start_index + j) = inverse_block(i, j);

	}

	//std::cout << iD.toDense() << std::endl;

};
void solver::assembleH() {


	H1.setZero();
	H2.setZero();

	assemble_λ();

	Eigen::Matrix3d block1;
	Eigen::Matrix3d block2;

	for (unsigned k = 0; k < nk; k++) {


		t_pointer const K = mesh->get_triangle(k);
		unsigned const k_index = K->index;

		e_pointer const E0 = K->edges[0];
		e_pointer const E1 = K->edges[1];
		e_pointer const E2 = K->edges[2];

		unsigned const e_index0 = E0->index;
		unsigned const e_index1 = E1->index;
		unsigned const e_index2 = E2->index;


		for (unsigned m = 0; m < 3; m++) {
			for (unsigned Ei = 0; Ei < 3; Ei++) {

				block1(m, Ei) = λ(k_index, 0, m, Ei);
				block2(m, Ei) = λ(k_index, 1, m, Ei);

			}
		}

		real const coeff = -θ * dt;

		block1 *= coeff;
		block2 *= coeff;

		// s = 1
		H1.coeffRef(3 * k_index + 0, e_index0) = block1.coeff(0, 0);
		H1.coeffRef(3 * k_index + 0, e_index1) = block1.coeff(0, 1);
		H1.coeffRef(3 * k_index + 0, e_index2) = block1.coeff(0, 2);

		H1.coeffRef(3 * k_index + 1, e_index0) = block1.coeff(1, 0);
		H1.coeffRef(3 * k_index + 1, e_index1) = block1.coeff(1, 1);
		H1.coeffRef(3 * k_index + 1, e_index2) = block1.coeff(1, 2);

		H1.coeffRef(3 * k_index + 2, e_index0) = block1.coeff(2, 0);
		H1.coeffRef(3 * k_index + 2, e_index1) = block1.coeff(2, 1);
		H1.coeffRef(3 * k_index + 2, e_index2) = block1.coeff(2, 2);

		// s = 2
		H2.coeffRef(3 * k_index + 0, e_index0) = block2.coeff(0, 0);
		H2.coeffRef(3 * k_index + 0, e_index1) = block2.coeff(0, 1);
		H2.coeffRef(3 * k_index + 0, e_index2) = block2.coeff(0, 2);

		H2.coeffRef(3 * k_index + 1, e_index0) = block2.coeff(1, 0);
		H2.coeffRef(3 * k_index + 1, e_index1) = block2.coeff(1, 1);
		H2.coeffRef(3 * k_index + 1, e_index2) = block2.coeff(1, 2);

		H2.coeffRef(3 * k_index + 2, e_index0) = block2.coeff(2, 0);
		H2.coeffRef(3 * k_index + 2, e_index1) = block2.coeff(2, 1);
		H2.coeffRef(3 * k_index + 2, e_index2) = block2.coeff(2, 2);

	}

	//std::cout << H1.toDense() << std::endl << std::endl;
	//std::cout << H2.toDense() << std::endl;

};
void solver::assembleG() {


	real const coeff = dt * (1.0 - θ);

	for (unsigned k = 0; k < nk; k++) {

		unsigned const k_index = mesh->get_triangle(k)->index;

		// General θ scheme
		for (unsigned m = 0; m < 3; m++)
			G[3 * k_index + m] = π_n(k_index, 0, m, 0) + coeff * rkFp_n(k_index, 0, m, 0);

	}

};


void solver::assemblePressureSystemMatrix1() {


	real const coefficient = θ * dt;

	Eigen::Matrix3d block;
	Eigen::Matrix3d block1;
	Eigen::Matrix3d block2;


	// It is sufficient to zero only diagonal elements. Therefore, there is no need for += in the sequel
	for (unsigned e = 0; e < ne; e++) {

		R1iDH1.coeffRef(e, e) = 0.0;
		R1iDH2.coeffRef(e, e) = 0.0;
		R2iDH1.coeffRef(e, e) = 0.0;
		R2iDH2.coeffRef(e, e) = 0.0;

	}

	for (unsigned k = 0; k < nk; k++) {



		t_pointer const K = mesh->get_triangle(k);

		unsigned const k_index = K->index;


		/*****************************************************************************/
		/*                                                                           */
		/*    - Inverse of the matrix D										         */
		/*                                                                           */
		/*****************************************************************************/

		for (unsigned r = 0; r < 3; r++)
			for (unsigned s = 0; s < 3; s++)
				block(r, s) = δij(r, s) - coefficient * σ(k_index, 0, r, s);

		Eigen::Matrix3d const inverse_block = block.inverse();


		/*****************************************************************************/
		/*                                                                           */
		/*    - H1, H2														         */
		/*                                                                           */
		/*****************************************************************************/

		for (unsigned m = 0; m < 3; m++) {
			for (unsigned Ei = 0; Ei < 3; Ei++) {

				block1(m, Ei) = λ(k_index, 0, m, Ei);
				block2(m, Ei) = λ(k_index, 1, m, Ei);

			}
		}

		block1 *= -coefficient;
		block2 *= -coefficient;


		/*****************************************************************************/
		/*                                                                           */
		/*    - H1, H2 blocks multiplied by the blocks Inverse of D			         */
		/*                                                                           */
		/*****************************************************************************/

		Eigen::Matrix3d const iDH1block = inverse_block * block1;
		Eigen::Matrix3d const iDH2block = inverse_block * block2;



		/*****************************************************************************/
		/*                                                                           */
		/*    - Assemble of the resulting matrix R1 * D.inverse * H1		         */
		/*                                                                           */
		/*****************************************************************************/

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

					R1iDH1.coeffRef(e_index_i, e_index_i) += sum11;
					R1iDH2.coeffRef(e_index_i, e_index_i) += sum12;
					R2iDH1.coeffRef(e_index_i, e_index_i) += sum21;
					R2iDH2.coeffRef(e_index_i, e_index_i) += sum22;

					continue;

				}

				R1iDH1.coeffRef(e_index_i, e_index_j) = sum11;
				R1iDH2.coeffRef(e_index_i, e_index_j) = sum12;
				R2iDH1.coeffRef(e_index_i, e_index_j) = sum21;
				R2iDH2.coeffRef(e_index_i, e_index_j) = sum22;

			}
		}
	}

	R1iDH1M11 = R1iDH1 + M_j1_s1;
	R1iDH2M12 = R1iDH2 + M_j1_s2;
	R2iDH1M21 = R2iDH1 + M_j2_s1;
	R2iDH2M22 = R2iDH2 + M_j2_s2;

	/*
	smat _R1iDH1M11;
	smat _R1iDH2M12;
	smat _R2iDH1M21;
	smat _R2iDH2M22;

	omp_set_num_threads(4);
	#pragma omp parallel 
	{
		#pragma omp sections
		{
			#pragma omp section
			{
				_R1iDH1M11 = R1iDH1 + M_j1_s1;
			}
			#pragma omp section
			{
				_R1iDH2M12 = R1iDH2 + M_j1_s2;
			}
			#pragma omp section
			{
				_R2iDH1M21 = R2iDH1 + M_j2_s1;
			}
			#pragma omp section
			{
				_R2iDH2M22 = R2iDH2 + M_j2_s2;
			}
		}
	}
	*/


	std::vector<Eigen::Triplet<double>> tri;

	for (int i = 0; i < R1iDH1M11.outerSize(); i++)
		for (Eigen::SparseMatrix<double>::InnerIterator it(R1iDH1M11, i); it; ++it)
			tri.emplace_back(it.row() + 0, it.col() + 0, it.value());

	for (int i = 0; i < R1iDH2M12.outerSize(); i++)
		for (Eigen::SparseMatrix<double>::InnerIterator it(R1iDH2M12, i); it; ++it)
			tri.emplace_back(it.row() + 0, it.col() + ne, it.value());

	for (int i = 0; i < R2iDH1M21.outerSize(); i++)
		for (Eigen::SparseMatrix<double>::InnerIterator it(R2iDH1M21, i); it; ++it)
			tri.emplace_back(it.row() + ne, it.col() + 0, it.value());

	for (int i = 0; i < R2iDH2M22.outerSize(); i++)
		for (Eigen::SparseMatrix<double>::InnerIterator it(R2iDH2M22, i); it; ++it)
			tri.emplace_back(it.row() + ne, it.col() + ne, it.value());

	/*
	for (int i = 0; i < R1iDH1M11.outerSize(); i++) {

		for (Eigen::SparseMatrix<double>::InnerIterator it(R1iDH1M11, i); it; ++it)
			tri.emplace_back(it.row() + 0, it.col() + 0, it.value());
		for (Eigen::SparseMatrix<double>::InnerIterator it(R1iDH2M12, i); it; ++it)
			tri.emplace_back(it.row() + 0, it.col() + ne, it.value());
		for (Eigen::SparseMatrix<double>::InnerIterator it(R2iDH1M21, i); it; ++it)
			tri.emplace_back(it.row() + ne, it.col() + 0, it.value());
		for (Eigen::SparseMatrix<double>::InnerIterator it(R2iDH2M22, i); it; ++it)
			tri.emplace_back(it.row() + ne, it.col() + ne, it.value());

	}
	*/

	internalPressureSystem.setFromTriplets(tri.begin(), tri.end());

};
void solver::assemblePressureSystemMatrix2() {


	real const coefficient = θ * dt;

	Eigen::Matrix3d block;
	Eigen::Matrix3d block1;
	Eigen::Matrix3d block2;


	// It is sufficient to zero only diagonal elements. Therefore, there is no need for += in the sequel
	for (unsigned e = 0; e < ne; e++) {

		R1iDH1M11.coeffRef(e, e) = M_j1_s1.coeff(e, e);
		R1iDH2M12.coeffRef(e, e) = M_j1_s2.coeff(e, e);
		R2iDH1M21.coeffRef(e, e) = M_j2_s1.coeff(e, e);
		R2iDH2M22.coeffRef(e, e) = M_j2_s2.coeff(e, e);

	}

	for (unsigned k = 0; k < nk; k++) {



		t_pointer const K = mesh->get_triangle(k);

		unsigned const k_index = K->index;


		/*****************************************************************************/
		/*                                                                           */
		/*    - Inverse of the matrix D										         */
		/*                                                                           */
		/*****************************************************************************/

		for (unsigned r = 0; r < 3; r++)
			for (unsigned s = 0; s < 3; s++)
				block(r, s) = δij(r, s) - coefficient * σ(k_index, 0, r, s);

		Eigen::Matrix3d const inverse_block = block.inverse();


		/*****************************************************************************/
		/*                                                                           */
		/*    - H1, H2														         */
		/*                                                                           */
		/*****************************************************************************/

		for (unsigned m = 0; m < 3; m++) {
			for (unsigned Ei = 0; Ei < 3; Ei++) {

				block1(m, Ei) = λ(k_index, 0, m, Ei);
				block2(m, Ei) = λ(k_index, 1, m, Ei);

			}
		}

		block1 *= -coefficient;
		block2 *= -coefficient;


		/*****************************************************************************/
		/*                                                                           */
		/*    - H1, H2 blocks multiplied by the blocks Inverse of D			         */
		/*                                                                           */
		/*****************************************************************************/

		Eigen::Matrix3d const iDH1block = inverse_block * block1;
		Eigen::Matrix3d const iDH2block = inverse_block * block2;



		/*****************************************************************************/
		/*                                                                           */
		/*    - Assemble of the resulting matrix R1 * D.inverse * H1		         */
		/*                                                                           */
		/*****************************************************************************/

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

					R1iDH1M11.coeffRef(e_index_i, e_index_i) += sum11;
					R1iDH2M12.coeffRef(e_index_i, e_index_i) += sum12;
					R2iDH1M21.coeffRef(e_index_i, e_index_i) += sum21;
					R2iDH2M22.coeffRef(e_index_i, e_index_i) += sum22;

					continue;

				}

				R1iDH1M11.coeffRef(e_index_i, e_index_j) = sum11 + M_j1_s1.coeff(e_index_i, e_index_j);
				R1iDH2M12.coeffRef(e_index_i, e_index_j) = sum12 + M_j1_s2.coeff(e_index_i, e_index_j);
				R2iDH1M21.coeffRef(e_index_i, e_index_j) = sum21 + M_j2_s1.coeff(e_index_i, e_index_j);
				R2iDH2M22.coeffRef(e_index_i, e_index_j) = sum22 + M_j2_s2.coeff(e_index_i, e_index_j);

			}
		}
	}

	std::vector<Eigen::Triplet<double>> tri;

	for (int i = 0; i < R1iDH1M11.outerSize(); i++)
		for (Eigen::SparseMatrix<double>::InnerIterator it(R1iDH1M11, i); it; ++it)
			tri.emplace_back(it.row() + 0, it.col() + 0, it.value());

	for (int i = 0; i < R1iDH2M12.outerSize(); i++)
		for (Eigen::SparseMatrix<double>::InnerIterator it(R1iDH2M12, i); it; ++it)
			tri.emplace_back(it.row() + 0, it.col() + ne, it.value());

	for (int i = 0; i < R2iDH1M21.outerSize(); i++)
		for (Eigen::SparseMatrix<double>::InnerIterator it(R2iDH1M21, i); it; ++it)
			tri.emplace_back(it.row() + ne, it.col() + 0, it.value());

	for (int i = 0; i < R2iDH2M22.outerSize(); i++)
		for (Eigen::SparseMatrix<double>::InnerIterator it(R2iDH2M22, i); it; ++it)
			tri.emplace_back(it.row() + ne, it.col() + ne, it.value());

	/*
	for (int i = 0; i < R1iDH1M11.outerSize(); i++) {

		for (Eigen::SparseMatrix<double>::InnerIterator it(R1iDH1M11, i); it; ++it)
			tri.emplace_back(it.row() + 0, it.col() + 0, it.value());
		for (Eigen::SparseMatrix<double>::InnerIterator it(R1iDH2M12, i); it; ++it)
			tri.emplace_back(it.row() + 0, it.col() + ne, it.value());
		for (Eigen::SparseMatrix<double>::InnerIterator it(R2iDH1M21, i); it; ++it)
			tri.emplace_back(it.row() + ne, it.col() + 0, it.value());
		for (Eigen::SparseMatrix<double>::InnerIterator it(R2iDH2M22, i); it; ++it)
			tri.emplace_back(it.row() + ne, it.col() + ne, it.value());

	}
	*/

	internalPressureSystem.setFromTriplets(tri.begin(), tri.end());

};
void solver::assemblePressureSystemMatrix3() {


	real const coefficient = θ * dt;

	Eigen::Matrix3d block;
	Eigen::Matrix3d block1;
	Eigen::Matrix3d block2;


	// It is sufficient to zero only diagonal elements. Therefore, there is no need for += in the sequel
	for (unsigned e = 0; e < ne; e++) {

		internalPressureSystem.coeffRef(e,		e)		= M_j1_s1.coeff(e, e);
		internalPressureSystem.coeffRef(e,		e + ne)	= M_j1_s2.coeff(e, e);
		internalPressureSystem.coeffRef(e + ne, e)		= M_j2_s1.coeff(e, e);
		internalPressureSystem.coeffRef(e + ne, e + ne)	= M_j2_s2.coeff(e, e);

	}

	for (unsigned k = 0; k < nk; k++) {



		t_pointer const K = mesh->get_triangle(k);

		unsigned const k_index = K->index;


		/*****************************************************************************/
		/*                                                                           */
		/*    - Inverse of the matrix D										         */
		/*                                                                           */
		/*****************************************************************************/

		for (unsigned r = 0; r < 3; r++)
			for (unsigned s = 0; s < 3; s++)
				block(r, s) = δij(r, s) - coefficient * σ(k_index, 0, r, s);

		Eigen::Matrix3d const inverse_block = block.inverse();


		/*****************************************************************************/
		/*                                                                           */
		/*    - H1, H2														         */
		/*                                                                           */
		/*****************************************************************************/

		for (unsigned m = 0; m < 3; m++) {
			for (unsigned Ei = 0; Ei < 3; Ei++) {

				block1(m, Ei) = λ(k_index, 0, m, Ei);
				block2(m, Ei) = λ(k_index, 1, m, Ei);

			}
		}

		block1 *= -coefficient;
		block2 *= -coefficient;


		/*****************************************************************************/
		/*                                                                           */
		/*    - H1, H2 blocks multiplied by the blocks Inverse of D			         */
		/*                                                                           */
		/*****************************************************************************/

		Eigen::Matrix3d const iDH1block = inverse_block * block1;
		Eigen::Matrix3d const iDH2block = inverse_block * block2;



		/*****************************************************************************/
		/*                                                                           */
		/*    - Assemble of the resulting matrix R1 * D.inverse * H1		         */
		/*                                                                           */
		/*****************************************************************************/

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

					internalPressureSystem.coeffRef(e_index_i,		e_index_i)		+= sum11;
					internalPressureSystem.coeffRef(e_index_i,		e_index_i + ne)	+= sum12;
					internalPressureSystem.coeffRef(e_index_i + ne, e_index_i)		+= sum21;
					internalPressureSystem.coeffRef(e_index_i + ne, e_index_i + ne)	+= sum22;

					continue;

				}

				internalPressureSystem.coeffRef(e_index_i,		e_index_j)		= sum11 + M_j1_s1.coeff(e_index_i, e_index_j);
				internalPressureSystem.coeffRef(e_index_i,		e_index_j + ne)	= sum12 + M_j1_s2.coeff(e_index_i, e_index_j);
				internalPressureSystem.coeffRef(e_index_i + ne, e_index_j)		= sum21 + M_j2_s1.coeff(e_index_i, e_index_j);
				internalPressureSystem.coeffRef(e_index_i + ne, e_index_j + ne)	= sum22 + M_j2_s2.coeff(e_index_i, e_index_j);

			}
		}
	}

};
void solver::assemblePressureSystemMatrix4() {


	real const coefficient = θ * dt;

	Eigen::Matrix3d block;
	Eigen::Matrix3d block1;
	Eigen::Matrix3d block2;

	std::vector<Eigen::Triplet<double>> tri;
	//tri.reserve(estimation_of_entries);


	// It is sufficient to zero only diagonal elements. Therefore, there is no need for += in the sequel
	for (unsigned e = 0; e < ne; e++) {

		real const M11 = M_j1_s1.coeff(e, e);
		real const M12 = M_j1_s2.coeff(e, e);
		real const M21 = M_j2_s1.coeff(e, e);
		real const M22 = M_j2_s2.coeff(e, e);

		Eigen::Triplet<real> const T1(e,		e,		M11);
		Eigen::Triplet<real> const T2(e,		e + ne, M12);
		Eigen::Triplet<real> const T3(e + ne,	e,		M21);
		Eigen::Triplet<real> const T4(e + ne,	e + ne, M22);

		tri.push_back(T1);
		tri.push_back(T2);
		tri.push_back(T3);
		tri.push_back(T4);

	}

	for (unsigned k = 0; k < nk; k++) {



		t_pointer const K = mesh->get_triangle(k);

		unsigned const k_index = K->index;


		/*****************************************************************************/
		/*                                                                           */
		/*    - Inverse of the matrix D										         */
		/*                                                                           */
		/*****************************************************************************/

		for (unsigned r = 0; r < 3; r++)
			for (unsigned s = 0; s < 3; s++)
				block(r, s) = δij(r, s) - coefficient * σ(k_index, 0, r, s);

		Eigen::Matrix3d const inverse_block = block.inverse();


		/*****************************************************************************/
		/*                                                                           */
		/*    - H1, H2														         */
		/*                                                                           */
		/*****************************************************************************/

		for (unsigned m = 0; m < 3; m++) {
			for (unsigned Ei = 0; Ei < 3; Ei++) {

				block1(m, Ei) = λ(k_index, 0, m, Ei);
				block2(m, Ei) = λ(k_index, 1, m, Ei);

			}
		}

		block1 *= -coefficient;
		block2 *= -coefficient;


		/*****************************************************************************/
		/*                                                                           */
		/*    - H1, H2 blocks multiplied by the blocks Inverse of D			         */
		/*                                                                           */
		/*****************************************************************************/

		Eigen::Matrix3d const iDH1block = inverse_block * block1;
		Eigen::Matrix3d const iDH2block = inverse_block * block2;



		/*****************************************************************************/
		/*                                                                           */
		/*    - Assemble of the resulting matrix R1 * D.inverse * H1		         */
		/*                                                                           */
		/*****************************************************************************/

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

					Eigen::Triplet<real> const T1(e_index_i,		e_index_i,		sum11);
					Eigen::Triplet<real> const T2(e_index_i,		e_index_i + ne,	sum12);
					Eigen::Triplet<real> const T3(e_index_i + ne,	e_index_i,		sum21);
					Eigen::Triplet<real> const T4(e_index_i + ne,	e_index_i + ne,	sum22);

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

				Eigen::Triplet<real> const T1(e_index_i,		e_index_j,		M11);
				Eigen::Triplet<real> const T2(e_index_i,		e_index_j + ne, M12);
				Eigen::Triplet<real> const T3(e_index_i + ne,	e_index_j,		M21);
				Eigen::Triplet<real> const T4(e_index_i + ne,	e_index_j + ne, M22);

				tri.push_back(T1);
				tri.push_back(T2);
				tri.push_back(T3);
				tri.push_back(T4);

			}
		}
	}

	internalPressureSystem.setFromTriplets(tri.begin(), tri.end());

};



void solver::assemble_α() {


	// Quadrature weights and points on reference triangle
	quadrature_triangle const quad(quadrature_order);
	unsigned const number_of_quadrature_points = quad.number_of_points;


	//Matrix<real> integral(8, 8);
	//Matrix<real> basisRaviartThomas(2, 8);
	//Matrix<real> JF(2, 2);

	Eigen::MatrixXd integral(8, 8);
	Eigen::MatrixXd basisRaviartThomas(2, 8);
	Eigen::MatrixXd JF(2, 2);

	α.setZero();

	for (unsigned k = 0; k < nk; k++) {


		t_pointer const K = mesh->get_triangle(k);
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

		real const detJF = (real) abs((x1 - x0)*(y2 - y0) - (x2 - x0)*(y1 - y0));

		//JF(0, 0) = x1 - x0;
		//JF(0, 1) = x2 - x0;
		//JF(1, 0) = y1 - y0;
		//JF(1, 1) = y2 - y0;


		JF.coeffRef(0, 0) = x1 - x0;
		JF.coeffRef(0, 1) = x2 - x0;
		JF.coeffRef(1, 0) = y1 - y0;
		JF.coeffRef(1, 1) = y2 - y0;

		/*
		real const idetJF = (real) 1.0 / abs(affineMappingMatrixDeterminant[k_index]);

		JF(0, 0) = affineMappingMatrixJF(k_index, 0, 0, 0);
		JF(0, 1) = affineMappingMatrixJF(k_index, 0, 0, 1);
		JF(1, 0) = affineMappingMatrixJF(k_index, 0, 1, 0);
		JF(1, 1) = affineMappingMatrixJF(k_index, 0, 1, 1);
		*/

		integral.setZero();

		for (unsigned n = 0; n < number_of_quadrature_points; n++) {


			real const s = (real) quad.points_x[n];
			real const t = (real) quad.points_y[n];
			real const w = (real) 0.5 * quad.weigths[n];

			// Corresponding coordinates on the element K
			real const x = x0 + (x1 - x0)*s + (x2 - x0)*t;
			real const y = y0 + (y1 - y0)*s + (y2 - y0)*t;

			// Get inverse of permeability tensor

			/*
			//Matrix<real> iK(2, 2);

			//permeability(x, y, iK);
			//iK.inverseInPlace();
			*/

			Eigen::Matrix2d iK(2, 2);

			iK.coeffRef(0, 0) = 1.0;
			iK.coeffRef(0, 1) = 0.0;
			iK.coeffRef(1, 0) = 0.0;
			iK.coeffRef(1, 1) = 1.0;

			evaluate_raviartthomas_basis(s, t, basisRaviartThomas);


			/*
			for (unsigned i = 0; i < 8; i++) {


				Vector<real> const JWi = JF * basisRaviartThomas.getColumn(i);

				for (unsigned j = i; j < 8; j++) {


					Vector<real> const JWj = JF * basisRaviartThomas.getColumn(j);

					real const quadraticForm = dot(JWi, iK*JWj);

					integral(i, j) = integral(i, j) + idetJF * weight * quadraticForm;

					if (i != j)
						integral(j, i) = integral(i, j);

				}
			}
			*/


			for (unsigned i = 0; i < 8; i++) {


				Eigen::Vector2d const JWi = JF * basisRaviartThomas.col(i);

				//real const w0 = JWi(0);
				//real const w1 = JWi(1);

				for (unsigned j = 0; j < 8; j++) {


					Eigen::Vector2d const JWj = JF * basisRaviartThomas.col(j);

					//real const w00 = JWj(0);
					//real const w11 = JWj(1);

					real const quadraticForm = JWi.dot(iK*JWj);

					integral.coeffRef(i, j) = integral(i, j) + w * quadraticForm;

					//if (i != j)
					//	integral.coeffRef(j, i) = integral(i, j);

				}
			}


		}

		integral = integral / detJF;

		/*
		//std::cout << integral << std::endl;


		//for (unsigned i = 0; i < 8; i++)
		//	for (unsigned j = 0; j < 8; j++)
		//		if (abs(integral(i, j)) <= INTEGRAL_PRECISION)
		//			integral(i, j) = (real) 0.0;

		//
		//integral.inverseInPlace();

		//for (unsigned i = 0; i < 8; i++)
		//	for (unsigned j = 0; j < 8; j++)
		//		if (abs(integral(i, j)) <= INTEGRAL_PRECISION)
		//			integral(i, j) = (real) 0.0;


		//for (unsigned i = 0; i < 8; i++)
		//	for (unsigned j = 0; j < 8; j++)
		//		α.setCoeff(k_index, 0, i, j) = integral(i, j);
		*/

		for (unsigned i = 0; i < 8; i++)
			for (unsigned j = 0; j < 8; j++)
				integral.coeffRef(i, j) = abs(integral(i, j)) < INTEGRAL_PRECISION ? 0.0 : integral(i, j);

		integral = integral.inverse();

		for (unsigned i = 0; i < 8; i++)
			for (unsigned j = 0; j < 8; j++)
				integral.coeffRef(i, j) = abs(integral(i, j)) < INTEGRAL_PRECISION ? 0.0 : integral(i, j);


		for (unsigned i = 0; i < 8; i++)
			for (unsigned j = 0; j < 8; j++)
				α.setCoeff(k_index, 0, i, j) = integral(i, j);

	}

};
void solver::assemble_β() {


	// Quadrature weights and points on reference triangle
	quadrature_triangle const quad(quadrature_order);
	unsigned const number_of_quadrature_points = quad.number_of_points;

	//Matrix<real> integral(8, 3);
	//Vector<real> basisPolynomial(3);
	//Vector<real> basisRaviartThomasDivergence(8);


	Eigen::Matrix<real, 8, 3> integral;
	Eigen::VectorXd basisPolynomial(3);
	Eigen::VectorXd basisRaviartThomasDivergence(8);


	β.setZero();

	for (unsigned k = 0; k < nk; k++) {


		unsigned const k_index = mesh->get_triangle(k)->index;


		integral.setZero();

		for (unsigned n = 0; n < number_of_quadrature_points; n++) {


			real const s = (real) quad.points_x[n];
			real const t = (real) quad.points_y[n];
			real const w = (real) 0.5 * quad.weigths[n];


			evaluate_raviartthomas_basis_divergence(s, t, basisRaviartThomasDivergence);
			evaluate_polynomial_basis(s, t, basisPolynomial);


			for (unsigned i = 0; i < 8; i++) {


				real const dWi = basisRaviartThomasDivergence(i);

				for (unsigned j = 0; j < 3; j++) {


					real const Phij = basisPolynomial(j);

					//integral(i, j) = integral(i, j) + weight * dWi * Phij;
					integral.coeffRef(i, j) = integral(i, j) + w * dWi * Phij;

				}
			}
		}

		//for (unsigned i = 0; i < 8; i++)
		//	for (unsigned j = 0; j < 3; j++)
		//		integral(i, j) = (real)abs(integral(i, j)) < INTEGRAL_PRECISION ? 0.0 : integral(i, j);

		for (unsigned i = 0; i < 8; i++)
			for (unsigned j = 0; j < 3; j++)
				integral.coeffRef(i, j) = (real)abs(integral(i, j)) < INTEGRAL_PRECISION ? 0.0 : integral(i, j);

		for (unsigned m = 0; m < 8; m++)
			for (unsigned j = 0; j < 3; j++)
				β.setCoeff(k_index, 0, m, j) = integral(m, j);

	}

};
void solver::assemble_χ() {


	// Quadrature weights and points on reference segment [-1,1]
	gauss_quadrature_1D const quad(quadrature_order);
	unsigned const number_of_quadrature_points = quad.number_of_points;


	Eigen::MatrixXd normals(2, 3);
	Eigen::VectorXd parametrization(2);
	//Eigen::VectorXd parametrizationDerivative(2);
	Eigen::MatrixXd basisRaviartThomas(2, 8);
	Eigen::VectorXd basisEdgePolynomial(2);

	evaluate_edge_normal(normals);


	χ.setZero();

	for (unsigned k = 0; k < nk; k++) {


		unsigned const k_index = mesh->get_triangle(k)->index;


		for (unsigned El = 0; El < 3; El++) {


			real const a = (real) 0.0;
			real const b = (real)El != 0 ? 1.0 : sqrt(2.0);

			real const c = (real)(b - a) / 2.0;
			real const d = (real)(b + a) / 2.0;

			Eigen::VectorXd const normal = normals.col(El);

			for (unsigned n = 0; n < number_of_quadrature_points; n++) {


				real const x = (real) quad.points[n] * c + d;
				real const w = (real) quad.weigths[n] * c;


				evaluate_edge_parametrization(x, El, parametrization);
				//evaluate_edge_parametrization_derivative(x, El, parametrizationDerivative);

				Eigen::VectorXd const r = parametrization;
				//Eigen::VectorXd const dr = parametrizationDerivative;


				real const s = r(0);
				real const t = r(1);
				//real const drNorm = dr.norm();
				real const drNorm = 1.0;


				evaluate_raviartthomas_basis(s, t, basisRaviartThomas);
				evaluate_edge_polynomial_basis(x, El, basisEdgePolynomial);


				for (unsigned m = 0; m < 8; m++) {


					Eigen::VectorXd const Wm = basisRaviartThomas.col(m);

					real const dotProduct = Wm.dot(normal);


					for (unsigned s = 0; s < 2; s++) {

						real const varPhis = basisEdgePolynomial(s);

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
void solver::assemble_η() {


	// Quadrature weights and points on reference triangle
	quadrature_triangle const quad(quadrature_order);
	unsigned const number_of_quadrature_points = quad.number_of_points;

	//Matrix<real> integral(3, 3);
	//Vector<real> basisPolynomial(3);

	Eigen::MatrixXd integral(3, 3);
	Eigen::VectorXd basisPolynomial(3);


	η.setZero();

	for (unsigned k = 0; k < nk; k++) {


		t_pointer const K = mesh->get_triangle(k);
		unsigned const k_index = K->index;


		v_pointer const a = K->vertices[0];
		v_pointer const b = K->vertices[1];
		v_pointer const c = K->vertices[2];

		real const x0 = a->x;
		real const y0 = a->y;

		real const x1 = b->x;
		real const y1 = b->y;

		real const x2 = c->x;
		real const y2 = c->y;

		real const detJF = (real)abs((x1 - x0)*(y2 - y0) - (x2 - x0)*(y1 - y0));


		integral.setZero();

		for (unsigned n = 0; n < number_of_quadrature_points; n++) {


			real const s = (real) quad.points_x[n];
			real const t = (real) quad.points_y[n];
			real const w = (real) 0.5 * quad.weigths[n];

			evaluate_polynomial_basis(s, t, basisPolynomial);


			for (unsigned m = 0; m < 3; m++) {


				real const Phim = basisPolynomial(m);

				for (unsigned j = 0; j < 3; j++) {


					real const Phij = basisPolynomial(j);

					//integral(i, j) = integral(i, j) + weight * Phii * Phij;

					//if (i != j)
					//	integral(j, i) = integral(i, j);

					integral.coeffRef(m, j) = integral(m, j) + w * Phim * Phij;

					//if (m != j)
					//	integral.coeffRef(j, m) = integral(m, j);

				}
			}
		}

		integral *= detJF;

		//integral.inverseInPlace();

		for (unsigned i = 0; i < 3; i++)
			for (unsigned j = 0; j < 3; j++)
				integral.coeffRef(i, j) = abs(integral(i, j)) < INTEGRAL_PRECISION ? 0.0 : integral(i, j);

		integral = integral.inverse();

		for (unsigned i = 0; i < 3; i++)
			for (unsigned j = 0; j < 3; j++)
				integral.coeffRef(i, j) = abs(integral(i, j)) < INTEGRAL_PRECISION ? 0.0 : integral(i, j);

		for (unsigned i = 0; i < 3; i++)
			for (unsigned j = 0; j < 3; j++)
				η.setCoeff(k_index, 0, i, j) = integral(i, j);

	}

};
void solver::assemble_τ() {


	// Quadrature weights and points on reference triangle
	quadrature_triangle const quad(quadrature_order);
	unsigned const number_of_quadrature_points = quad.number_of_points;

	//Matrix<real> basisRaviartThomas(2, 8);
	//Vector<real> basisPolynomial(3);
	//Matrix<real> basisPolynomialGradient(2, 3);

	Eigen::MatrixXd basisRaviartThomas(2, 8);
	Eigen::VectorXd basisPolynomial(3);
	Eigen::MatrixXd basisPolynomialGradient(2, 3);


	τ.setZero();

	for (unsigned k = 0; k < nk; k++) {


		unsigned const k_index = mesh->get_triangle(k)->index;


		for (unsigned n = 0; n < number_of_quadrature_points; n++) {


			real const s = (real) quad.points_x[n];
			real const t = (real) quad.points_y[n];
			real const w = (real) 0.5 * quad.weigths[n];


			evaluate_raviartthomas_basis(s, t, basisRaviartThomas);
			evaluate_polynomial_basis(s, t, basisPolynomial);
			evaluate_polynomial_basis_gradient(s, t, basisPolynomialGradient);


			for (unsigned m = 0; m < 3; m++) {

				//Vector<real> const dPhim = basisPolynomialGradient.getColumn(m);
				Eigen::VectorXd const dPhim = basisPolynomialGradient.col(m);

				for (unsigned j = 0; j < 8; j++) {


					//Vector<real> const Wj = basisRaviartThomas.getColumn(j);
					//real const dotProduct = dot(Wj, dPhim);

					Eigen::VectorXd const Wj = basisRaviartThomas.col(j);
					real const dotProduct = Wj.dot(dPhim);

					for (unsigned l = 0; l < 3; l++) {

						real const Phil = basisPolynomial(l);

						τ.setCoeff(k_index, m, j, l) = τ(k_index, m, j, l) + w * dotProduct * Phil;

					}
				}
			}
		}

		for (unsigned m = 0; m < 3; m++)
			for (unsigned i = 0; i < 8; i++)
				for (unsigned j = 0; j < 3; j++)
					τ.setCoeff(k_index, m, i, j) = abs(τ(k_index, m, i, j)) < INTEGRAL_PRECISION ? 0.0 : τ(k_index, m, i, j);

	}

};
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
			Eigen::VectorXd const normal = normals.col(El);


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


				real const tC = upwindConcentration(s, t, K, normal, El, n);


				for (unsigned m = 0; m < 3; m++) {


					real const Phim = basisPolynomial(m);

					for (unsigned j = 0; j < 8; j++) {

						//Vector<real> const Wj = basisRaviartThomas.getColumn(j);
						Eigen::VectorXd const Wj = basisRaviartThomas.col(j);

						//real const dotProduct = dot(Wj, normal);
						real const dotProduct = Wj.dot(normal);

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
*/

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


	//gauss_quadrature_1D const quad(quadrature_order);
	//unsigned const number_of_quadrature_points = quad.number_of_points;
	unsigned const number_of_quadrature_points = number_of_quad_points;

	//Matrix<real> normals(2, 3);
	//Vector<real> parametrization(2);

	//Eigen::MatrixXd normals(2, 3);
	//Eigen::VectorXd parametrization(2);

	//evaluate_edge_normal(normals);


	δ.setZero();

	for (unsigned k = 0; k < nk; k++) {


		t_pointer const K = mesh->get_triangle(k);
		unsigned const k_index = K->index;


		//CoeffMatrixOfElement<3, 3, 8> delta_temp(δ, k_index);


		for (unsigned El = 0; El < 3; El++) {


			//real const a = (real) 0.0;
			//real const b = (real)(El != 0) ? 1.0 : sqrt(2.0);

			//real const c = (real)(b - a) / 2.0;
			//real const d = (real)(b + a) / 2.0;

			//Vector<real> const normal = normals.getColumn(El);
			//Eigen::VectorXd const referenceNormal = normals.col(El);


			for (unsigned n = 0; n < number_of_quadrature_points; n++) {


				//real const x = (real)quad.points[n] * c + d;
				//real const w = (real)quad.weigths[n] * c;


				//evaluate_edge_parametrization(x, El, parametrization);
				//evaluate_edge_parametrization_derivative(x, El, parametrizationDerivative);


				//real const s = parametrization(0);
				//real const t = parametrization(1);
				//real const drNorm = parametrizationDerivative.norm();
				//real const drNorm = 1.0;


				//real const tC = upwindConcentration(s, t, K, basisRaviartThomas, basisPolynomial, referenceNormal, denominator, El, x);
				//real const tC = upwindConcentration(s, t, K, referenceNormal, El, n);
				//real const tC = upwindConcentration2(K, referenceNormal, El, n);
				//real const tC = upwindConcentration3(K, referenceNormal, El, n);
				//real const tC = upwindConcentration4(K, El, n);
				//real const tC = upwindConcentration5(K, El, n);
				//real const tC = upwindConcentration6(K, El, n);
				//real const tC = upwindConcentration7(K, El, n);
				real const tC = upwindConcentration8(K, El, n);

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
void solver::assemble_δ00() {


	unsigned const number_of_quadrature_points = number_of_quad_points;


	δ.setZero();

	for (unsigned k = 0; k < nk; k++) {


		t_pointer const K = mesh->get_triangle(k);
		unsigned const k_index =  K->index;


		for (unsigned El = 0; El < 3; El++) {

			for (unsigned n = 0; n < number_of_quadrature_points; n++) {

				//real const tC = upwindConcentration4(K, El, n);
				//real const tC = upwindConcentration5(K, El, n);
				//real const tC = upwindConcentration6(K, El, n);
				//real const tC = upwindConcentration7(K, El, n);
				real const tC = upwindConcentration8(K, El, n);


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
void solver::assemble_δ11() {


	unsigned const number_of_quadrature_points = number_of_quad_points;

	//Matrix<real> normals(2, 3);
	Eigen::MatrixXd normals(2, 3);

	evaluate_edge_normal(normals);


	δ.setZero();

	for (unsigned k = 0; k < nk; k++) {


		t_pointer const K = mesh->get_triangle(k);
		unsigned const k_index = K->index;


		for (unsigned El = 0; El < 3; El++) {

			//Vector<real> const normal = normals.getColumn(El);
			Eigen::VectorXd const referenceNormal = normals.col(El);

			for (unsigned n = 0; n < number_of_quadrature_points; n++) {

				real const tC = upwindConcentration2(K, referenceNormal, El, n);
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
void solver::assemble_δ12() {


	unsigned const number_of_quadrature_points = number_of_quad_points;

	//Matrix<real> normals(2, 3);
	Eigen::MatrixXd normals(2, 3);

	evaluate_edge_normal(normals);


	δ.setZero();

	for (unsigned k = 0; k < nk; k++) {


		t_pointer const K = mesh->get_triangle(k);
		unsigned const k_index = K->index;


		for (unsigned El = 0; El < 3; El++) {

			//Vector<real> const normal = normals.getColumn(El);
			Eigen::VectorXd const referenceNormal = normals.col(El);

			for (unsigned n = 0; n < number_of_quadrature_points; n++) {

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
void solver::assemble_δ5() {


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
void solver::assemble_δ6() {

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
void solver::assemble_δ7() {

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


		Eigen::MatrixXd JF(2, 2);

		v_pointer const a = K->vertices[0];
		v_pointer const b = K->vertices[1];
		v_pointer const c = K->vertices[2];

		real const x0 = (real) a->x;
		real const y0 = (real) a->y;

		real const x1 = (real) b->x;
		real const y1 = (real) b->y;

		real const x2 = (real) c->x;
		real const y2 = (real) c->y;

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
				//real const drNorm = 1.0;


				evaluate_raviartthomas_basis(s, t, basisRaviartThomas);
				evaluate_polynomial_basis(s, t, basisPolynomial);



				
				//real normalVelocity = 0.0;  
				//if (E->marker == E_MARKER::NEUMANN)
				//	normalVelocity = NEUMANN_GAMMA_Q_velocity(E, nt*dt);
				//else
				//	normalVelocity = velocityInNormalDirection(K, E, basisRaviartThomas, referenceNormal, denominator);
				

				real const tC = upwindConcentration(s, t, K, basisRaviartThomas, basisPolynomial, referenceNormal, denominator, El, quad.points[n]);


				//Eigen::Matrix<real, 3, 8> integral;
				//integral.setZero();

				for (unsigned m = 0; m < 3; m++) {


					real const Phim = basisPolynomial(m);


					for (unsigned j = 0; j < 8; j++) {

						//Vector<real> const Wj = basisRaviartThomas.getColumn(j);
						//Eigen::VectorXd const Wj = basisRaviartThomas.col(j);

						//real const dotProduct = dot(Wj, normal);
						//real const dotProduct = Wj.dot(referenceNormal);

						real const value = raviartThomasBasisDotNormalTimesPolynomialBasis(n, j, El, m);
						δ.setCoeff(k_index, m, El, j) = δ(k_index, m, El, j) + tC * value;

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
*/
void solver::assemble_γ() {


	for (unsigned k = 0; k < nk; k++) {


		t_pointer const K = mesh->get_triangle(k);
		unsigned const k_index = K->index;


		//CoeffMatrixOfElement<3, 3, 8> delta_temp(δ, k_index);
		//CoeffMatrixOfElement<3, 8, 3> tau_temp(τ, k_index);
		//CoeffMatrixOfElement<1, 3, 8> gamma_temp(γ, k_index);


		for (unsigned m = 0; m < 3; m++) {

			for (unsigned j = 0; j < 8; j++) {


				real val1 = 0.0;
				real val2 = 0.0;

				for (unsigned l = 0; l < 3; l++)
					val1 += δ(k_index, m, l, j);

				for (unsigned l = 0; l < 3; l++)
					val2 += τ(k_index, m, j, l) * ξ_prev(k_index, 0, l, 0);

				γ.setCoeff(k_index, 0, m, j) = val1 - val2;

				//for (unsigned l = 0; l < 3; l++)
				//	val1 += delta_temp(m, l, j);

				//for (unsigned l = 0; l < 3; l++)
				//	val2 += tau_temp(m, j, l) * ξ_prev(k_index, 0, l, 0);


				//gamma_temp.setCoeff(0, m, j) = val1 - val2;

			}
		}
	}

};


void solver::assemble_σ() {


	σ.setZero();

	for (unsigned k = 0; k < nk; k++) {


		t_pointer const K = mesh->get_triangle(k);
		unsigned const k_index = K->index;

		real const coeff = betas_prev[k_index] / (viscosities[k_index] * porosities[k_index]);


		for (unsigned m = 0; m < 3; m++) {

			for (unsigned l = 0; l < 3; l++) {


				real val = 0.0;
				for (unsigned q = 0; q < 3; q++) {


					real GAB = 0.0;
					for (unsigned j = 0; j < 8; j++) {


						real AB = 0.0;
						for (unsigned i = 0; i < 8; i++)
							AB += α(k_index, 0, j, i) * β(k_index, 0, i, l);

						GAB += γ(k_index, 0, q, j) * AB;

					}

					val += η(k_index, 0, m, q) * GAB;

				}

				σ.setCoeff(k_index, 0, m, l) = -coeff * val;

			}
		}
	}

};
void solver::assemble_λ() {


	for (unsigned k = 0; k < nk; k++) {


		t_pointer const K = mesh->get_triangle(k);
		unsigned const k_index = K->index;

		real const coeff = betas_prev[k_index] / (viscosities[k_index] * porosities[k_index]);


		for (unsigned s = 0; s < 2; s++) {

			for (unsigned m = 0; m < 3; m++) {

				for (unsigned El = 0; El < 3; El++) {


					real val = 0.0;
					for (unsigned q = 0; q < 3; q++) {


						real YACHI = 0.0;
						for (unsigned j = 0; j < 8; j++) {


							real ACHI = 0.0;
							for (unsigned i = 0; i < 8; i++)
								ACHI += α(k_index, 0, i, j)*χ(k_index, i, El, s);

							YACHI += γ(k_index, 0, q, j)*ACHI;

						}

						val += η(k_index, 0, m, q)*YACHI;

					}

					λ.setCoeff(k_index, s, m, El) = coeff * val;

				}
			}
		}
	}

};
void solver::assemble_φ() {};
void solver::assemble_ψ() {};


void solver::getSolution() {


	// Compute initial values for iterations
	computeTracePressures();
	computeVelocities();
	computeBetas();

	// Set iteration level l := 0
	π_n = π;
	π_prev = π;

	ξ_n = ξ;
	ξ_prev = ξ;

	rkFc_n = rkFc;
	rkFp_n = rkFp;

	hard_copy(betas_prev, betas, nk);



	unsigned counter = 0;

	while (counter < MAX_IT) {


		// 1. assemble_δ0(); 
		// 2. assemble_δ11();
		// 3. assemble_δ12();

		assemble_δ00();

		//assemble_δ(); 
		//assemble_δ1(); 
		//assemble_δ2(); 
		//assemble_δ3(); 
		//assemble_δ5();

		assemble_γ();

		computePressureEquation();
		computeVelocities();

		updateConcentrations_explicit();

		//concentrationCorrection();

		computeBetas();


		if (stopCriterion()) {

			break;

		}
			


		// Set new iteration level	l := l+1
		π_prev = π;
		ξ_prev = ξ;

		//std::ofstream txtFile;
		//txtFile.open("C:\\Users\\pgali\\Desktop\\output_pressure.txt");
		//exportSolution(txtFile);
		//txtFile.close();


		//hard_copy(betas_prev, betas, nk);

		//for (unsigned k = 0; k < nk; k++)
		//	betas_prev[k] = betas[k];


		counter++;

	}

	assemble_λ();
	assemble_σ();

	// This is needed for the next time level
	for (unsigned k = 0; k < nk; k++) {

		unsigned const k_index = mesh->get_triangle(k)->index;

		for (unsigned m = 0; m < 3; m++) {

			real val1 = 0.0;
			real val2 = 0.0;

			for (unsigned j = 0; j < 3; j++)
				val1 += σ(k_index, 0, m, j) * π(k_index, 0, j, 0);

			for (unsigned El = 0; El < 3; El++)
				for (unsigned s = 0; s < 2; s++)
					val2 += λ(k_index, s, m, El) * tπ(k_index, 0, El, s);

			rkFp.setCoeff(k_index, 0, m, 0) = val1 + val2;

		}
	}


	if (nt % 1000 == 0)
		std::cout << nt << " - Iterations : " << counter << std::endl;

};




void solver::exportSolution(std::ofstream & txtFile) {


	double const t = nt * dt;


	for (unsigned k = 0; k < nk; k++) {


		t_pointer const	K = mesh->get_triangle(k);

		unsigned const index = K->index;

		v_pointer const a = K->vertices[0];
		v_pointer const b = K->vertices[1];
		v_pointer const c = K->vertices[2];

		double const x0 = a->x;
		double const y0 = a->y;

		double const x1 = b->x;
		double const y1 = b->y;

		double const x2 = c->x;
		double const y2 = c->y;

		double const x[3] = { x0, x1, x2 };
		double const y[3] = { y0, y1, y2 };


		Eigen::MatrixXd iJF(2, 2);

		iJF(0, 0) = x1 - x0;
		iJF(0, 1) = x2 - x0;
		iJF(1, 0) = y1 - y0;
		iJF(1, 1) = y2 - y0;

		iJF = iJF.inverse();

		Eigen::Vector2d P0;
		Eigen::Vector2d P1;
		Eigen::Vector2d P2;

		P0.coeffRef(0) = x0 - x0;
		P0.coeffRef(1) = y0 - y0;

		P1.coeffRef(0) = x1 - x0;
		P1.coeffRef(1) = y1 - y0;

		P2.coeffRef(0) = x2 - x0;
		P2.coeffRef(1) = y2 - y0;

		real const S0 = (iJF*P0)(0);
		real const T0 = (iJF*P0)(1);

		real const S1 = (iJF*P1)(0);
		real const T1 = (iJF*P1)(1);

		real const S2 = (iJF*P2)(0);
		real const T2 = (iJF*P2)(1);

		double const S[3] = { S0, S1, S2 };
		double const T[3] = { T0, T1, T2 };


		for (unsigned i = 0; i < 3; i++) {

			//// Pressures
			//real const value = π(index, 0, 0, 0) * phi1(S[i], T[i]) + π(index, 0, 1, 0) * phi2(S[i], T[i]) + π(index, 0, 2, 0) * phi3(S[i], T[i]);
			//// Concentrations
			real const value = ξ(index, 0, 0, 0) * phi1(S[i], T[i]) + ξ(index, 0, 1, 0) * phi2(S[i], T[i]) + ξ(index, 0, 2, 0) * phi3(S[i], T[i]);

			txtFile << std::setprecision(20) << x[i] << "	" << y[i] << "	" << value << "	" << index << std::endl;

		}

		txtFile << std::endl;

	}

};







void solver::compute_error(std::ofstream & txtFile) {


	//unsigned const nk = nk;
	//double const t = nt * dt;

	//element * K = NULL;

	//double errorL1 = 0.0;
	//double errorL2 = 0.0;
	//double errorMAX = 0.0;



	////	π	 ξ

	//for (unsigned k = 0; k < nk; k++) {


	//	K = mesh->getElement(k);

	//	double const a = K->nodes[0]->x;
	//	double const b = K->nodes[1]->x;


	//	quadrature quad(error_quadRule, a, b, METHOD::GAUSS);

	//	double y = 0.0;
	//	double w = 0.0;

	//	double s1 = 0.0;
	//	double s2 = 0.0;

	//	double pressure = 0.0;
	//	double analytic = 0.0;

	//	double val = 0.0;


	//	for (unsigned j = 0; j < error_quadRule; j++) {


	//		y = quad.points[j];
	//		w = quad.weigths[j];

	//		pressure = π(K, 0, 0)*phi1(y, a, b) + π(K, 0, 1)*phi2(y, a, b);
	//		analytic = barenblatt(y, t);

	//		val = abs(pressure - analytic);

	//		s1 += w * val;
	//		s2 += w * sqr(val);

	//		if (val > errorMAX)
	//			errorMAX = val;

	//	}

	//	errorL1 += s1;
	//	errorL2 += s2;

	//}

	//txtFile <<  "#L1 L2 MAX" << std::endl;
	//txtFile << std::setprecision(20) << errorL1 << std::endl;
	//txtFile << std::setprecision(20) << sqrt(errorL2) << std::endl;
	//txtFile << std::setprecision(20) << errorMAX << std::endl;
	//
	////std::cout << "Error L1 : " << errorL1 << std::endl;
	////std::cout << "Error L2 : " << sqrt(errorL2) << std::endl;
	////std::cout << "Error Max : " << errorMAX << std::endl;


};



























































































































/*






real solver::compute_α(t_pointer const K, unsigned const i, unsigned const j) {


	// Quadrature weights and points on reference triangle
	quadrature_triangle const quad(quadrature_order);
	unsigned const num_quad_points = quad.number_of_points;


	Matrix<real> integral(8, 8);
	Matrix<real> basisRaviartThomas(2, 8);
	Matrix<real> JF(2, 2);


	unsigned const k_index = K->index;

	real orientations[3];

	for (unsigned e = 0; e < 3; e++) {

		if (K->edges[e]->marker == E_MARKER::NEUMANN || K->edges[e]->marker == E_MARKER::DIRICHLET)
			orientations[e] = 1.0;
		else
			orientations[e] = edgeOrientation(k_index, 0, e, 0);

	}


	v_pointer const a = K->vertices[0];
	v_pointer const b = K->vertices[1];
	v_pointer const c = K->vertices[2];

	real const x0 = (real)a->x;
	real const y0 = (real)a->y;

	real const x1 = (real)b->x;
	real const y1 = (real)b->y;

	real const x2 = (real)c->x;
	real const y2 = (real)c->y;

	real const idetJF = (real) 1.0 / abs((x1 - x0)*(y2 - y0) - (x2 - x0)*(y1 - y0));

	JF(0, 0) = x1 - x0;
	JF(0, 1) = x2 - x0;
	JF(1, 0) = y1 - y0;
	JF(1, 1) = y2 - y0;


	integral.setZero();

	for (unsigned n = 0; n < num_quad_points; n++) {


		real const s = (real)quad.points_x[n];
		real const t = (real)quad.points_y[n];

		real const w = (real)quad.weigths[n];

		// There must be 1/2. See 'integrate_triangle' at "misc.h"
		real const weight = (real) 0.5 * w;

		// Corresponding coordinates on the element K
		real const x = x0 + (x1 - x0)*s + (x2 - x0)*t;
		real const y = y0 + (y1 - y0)*s + (y2 - y0)*t;

		// Get inverse of permeability tensor
		Matrix<real> iK(2, 2);

		permeability(x, y, iK);
		iK.inverseInPlace();

		// Evaluate Raviart-Thomas basis functions on the reference triangle at quadrature points s,t
		evaluate_raviartthomas_basis(s, t, basisRaviartThomas, orientations);


		for (unsigned i = 0; i < 8; i++) {


			Vector<real> const JWi = JF * basisRaviartThomas.getColumn(i);

			for (unsigned j = i; j < 8; j++) {


				Vector<real> const JWj = JF * basisRaviartThomas.getColumn(j);

				real const quadraticForm = dot(JWi, iK*JWj);

				integral(i, j) = integral(i, j) + idetJF * weight * quadraticForm;

				if (i != j)
					integral(j, i) = integral(i, j);

			}
		}


	}

	//std::cout << integral << std::endl;


	for (unsigned i = 0; i < 8; i++)
		for (unsigned j = 0; j < 8; j++)
			if (abs(integral(i, j)) <= INTEGRAL_PRECISION)
				integral(i, j) = (real) 0.0;


	integral.inverseInPlace();

	for (unsigned i = 0; i < 8; i++)
		for (unsigned j = 0; j < 8; j++)
			if (abs(integral(i, j)) <= INTEGRAL_PRECISION)
				integral(i, j) = (real) 0.0;


	for (unsigned i = 0; i < 8; i++)
		for (unsigned j = 0; j < 8; j++)
			α.setCoeff(k_index, 0, i, j) = integral(i, j);



};
real solver::compute_β(t_pointer const K, unsigned const m, unsigned const j) {


	// Quadrature weights and points on reference triangle
	quadrature_triangle const quad(quadrature_order);
	unsigned const num_quad_points = quad.number_of_points;

	Matrix<real> integral(8, 3);
	Vector<real> basisPolynomial(3);
	Vector<real> basisRaviartThomasDivergence(8);



	for (unsigned k = 0; k < nk; k++) {


		t_pointer const K = mesh->get_triangle(k);
		unsigned const k_index = K->index;


		real orientations[3];

		for (unsigned e = 0; e < 3; e++) {

			if (K->edges[e]->marker == E_MARKER::NEUMANN || K->edges[e]->marker == E_MARKER::DIRICHLET)
				orientations[e] = 1.0;
			else
				orientations[e] = edgeOrientation(k_index, 0, e, 0);

		}

		//for (unsigned e = 0; e < 3; e++)
		//	orientations[e] = edgeOrientation(k_index, 0, e, 0);


		v_pointer const a = K->vertices[0];
		v_pointer const b = K->vertices[1];
		v_pointer const c = K->vertices[2];

		real const x0 = (real)a->x;
		real const y0 = (real)a->y;

		real const x1 = (real)b->x;
		real const y1 = (real)b->y;

		real const x2 = (real)c->x;
		real const y2 = (real)c->y;

		//real const sig = (real) abs(((x1 - x0)*(y2 - y0) - (x2 - x0)*(y1 - y0))) / ((x1 - x0)*(y2 - y0) - (x2 - x0)*(y1 - y0));

		integral.setZero();

		for (unsigned n = 0; n < num_quad_points; n++) {


			real const s = (real)quad.points_x[n];
			real const t = (real)quad.points_y[n];
			real const w = (real)quad.weigths[n];

			real const weight = (real) 0.5 * w;

			// Values of divergence of basis function on reference triangle
			evaluate_raviartthomas_basis_divergence(s, t, basisRaviartThomasDivergence, orientations);
			// Phi1, Phi2, Phi3 are defined on reference triangle, therefore we need to input s,t and not x,y
			evaluate_polynomial_basis(s, t, basisPolynomial);


			for (unsigned i = 0; i < 8; i++) {


				real const dWi = basisRaviartThomasDivergence(i);

				for (unsigned j = 0; j < 3; j++) {


					real const Phij = basisPolynomial(j);

					integral(i, j) = integral(i, j) + weight * dWi * Phij;
					//β(k_index, 0, i, j) += weight * dWi * Phij;

				}
			}
		}

		for (unsigned i = 0; i < 8; i++)
			for (unsigned j = 0; j < 3; j++)
				if (abs(integral(i, j)) <= INTEGRAL_PRECISION)
					integral(i, j) = (real) 0.0;

		for (unsigned m = 0; m < 8; m++)
			for (unsigned j = 0; j < 3; j++)
				β.setCoeff(k_index, 0, m, j) = integral(m, j);

		//std::cout << integral << std::endl;

	}

};
real solver::compute_χ(t_pointer const K, unsigned const m, unsigned const El, unsigned const j) {


	// Quadrature weights and points on reference segment [-1,1]
	gauss_quadrature_1D const quad(quadrature_order);
	unsigned const num_quad_points = quad.number_of_points;

	Matrix <real> normals(2, 3);
	Matrix<real> parametrization(2, 3);
	Matrix<real> parametrizationDerivative(2, 8);
	Matrix<real> basisRaviartThomas(2, 8);
	Matrix<real> basisEdgePolynomial(2, 3);


	evaluate_edge_normal(normals);


	χ.setZero();

	for (unsigned k = 0; k < nk; k++) {


		t_pointer const K = mesh->get_triangle(k);
		unsigned const k_index = K->index;


		//--------------------------------
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

		real const detJF = (real)((x1 - x0)*(y2 - y0) - (x2 - x0)*(y1 - y0));

		JF(0, 0) = x1 - x0;
		JF(0, 1) = x2 - x0;
		JF(1, 0) = y1 - y0;
		JF(1, 1) = y2 - y0;
		//--------------------------------

		real orientations[3];

		//for (unsigned e = 0; e < 3; e++)
		//	orientations[e] = edgeOrientation(k_index, 0, e, 0);

		for (unsigned e = 0; e < 3; e++) {

			if (K->edges[e]->marker == E_MARKER::NEUMANN || K->edges[e]->marker == E_MARKER::DIRICHLET)
				orientations[e] = 1.0;
			else
				orientations[e] = edgeOrientation(k_index, 0, e, 0);

		}


		for (unsigned n = 0; n < num_quad_points; n++) {

			for (unsigned El = 0; El < 3; El++) {


				real const orientation = orientations[El];


				real const a = (real) 0.0;
				real const b = (real)El != 0 ? 1.0 : sqrt(2.0);

				real const c = (real)(b - a) / 2.0;
				real const d = (real)(b + a) / 2.0;

				real const x = (real)quad.points[n] * c + d;
				real const w = (real)quad.weigths[n] * c;

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
				evaluate_edge_polynomial_basis(x, basisEdgePolynomial, orientation); // get rid of that orientation. Every time a compute it with the global edge orientation

				Vector<real> const edgeBasis = basisEdgePolynomial.getColumn(El);


				//real const coeff = El == 0 ? sqrt(2.0) : 1.0;


				for (unsigned m = 0; m < 8; m++) {


					Vector<real> const Wm = basisRaviartThomas.getColumn(m);
					//real const dotProduct = orientation * dot(Wm, normal);
					real const dotProduct = dot(Wm, normal);

					for (unsigned j = 0; j < 2; j++) {

						real const varPhij = edgeBasis(j);

						χ.setCoeff(k_index, m, El, j) = χ(k_index, m, El, j) + w * dotProduct * varPhij * drNorm;
						//χ.setCoeff(k_index, m, El, j) = χ(k_index, m, El, j) + coeff * w * dotProduct * varPhij * drNorm;

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

		//	//std::cout << integral << std::endl;

		//}

	}

};
real solver::compute_η(t_pointer const K, unsigned const m, unsigned const j) {


	// Quadrature weights and points on reference triangle
	quadrature_triangle const quad(quadrature_order);
	unsigned const num_quad_points = quad.number_of_points;

	Matrix<real> integral(3, 3);
	Vector<real> basisPolynomial(3);


	for (unsigned k = 0; k < nk; k++) {


		t_pointer const K = mesh->get_triangle(k);
		unsigned const k_index = K->index;


		v_pointer const a = K->vertices[0];
		v_pointer const b = K->vertices[1];
		v_pointer const c = K->vertices[2];

		real const x0 = a->x;
		real const y0 = a->y;

		real const x1 = b->x;
		real const y1 = b->y;

		real const x2 = c->x;
		real const y2 = c->y;

		real const detJF = (real)abs((x1 - x0)*(y2 - y0) - (x2 - x0)*(y1 - y0));

		integral.setZero();


		for (unsigned n = 0; n < num_quad_points; n++) {


			real const s = (real)quad.points_x[n];
			real const t = (real)quad.points_y[n];
			real const w = (real)quad.weigths[n];

			real const weight = 0.5 * w;

			// Phi1, Phi2, Phi3 are defined on reference triangle, therefore we need to input s,t and not x,y
			evaluate_polynomial_basis(s, t, basisPolynomial);


			for (unsigned i = 0; i < 3; i++) {


				real const Phii = basisPolynomial(i);

				for (unsigned j = i; j < 3; j++) {


					real const Phij = basisPolynomial(j);

					integral(i, j) = integral(i, j) + weight * Phii * Phij;

					if (i != j)
						integral(j, i) = integral(i, j);

				}
			}
		}

		integral = detJF * integral;

		// Inverse of symmetric 3x3 matrix
		integral.inverseInPlace();

		for (unsigned i = 0; i < 3; i++)
			for (unsigned j = 0; j < 3; j++)
				η.setCoeff(k_index, 0, i, j) = integral(i, j);

		//std::cout << integral << std::endl;

	}

};
real solver::compute_τ(t_pointer const K, unsigned const m, unsigned const j, unsigned const l) {


	// Quadrature weights and points on reference triangle
	quadrature_triangle const quad(quadrature_order);
	unsigned const num_quad_points = quad.number_of_points;

	Matrix<real> basisRaviartThomas(2, 8);
	Vector<real> basisPolynomial(3);
	Matrix<real> basisPolynomialGradient(3, 2);


	τ.setZero();

	for (unsigned k = 0; k < nk; k++) {


		t_pointer const K = mesh->get_triangle(k);
		unsigned const k_index = K->index;


		real orientations[3];

		//for (unsigned e = 0; e < 3; e++)
		//	orientations[e] = edgeOrientation(k_index, 0, e, 0);

		for (unsigned e = 0; e < 3; e++) {

			if (K->edges[e]->marker == E_MARKER::NEUMANN || K->edges[e]->marker == E_MARKER::DIRICHLET)
				orientations[e] = 1.0;
			else
				orientations[e] = edgeOrientation(k_index, 0, e, 0);

		}


		for (unsigned n = 0; n < num_quad_points; n++) {


			real const s = (real)quad.points_x[n];
			real const t = (real)quad.points_y[n];
			real const w = (real)quad.weigths[n];

			real const weight = 0.5 * w;

			// Values of basis function on reference triangle
			evaluate_raviartthomas_basis(s, t, basisRaviartThomas, orientations);
			evaluate_polynomial_basis(s, t, basisPolynomial);
			evaluate_polynomial_basis_gradient(s, t, basisPolynomialGradient);


			for (unsigned m = 0; m < 3; m++) {


				Vector<real> const dPhim = basisPolynomialGradient.getColumn(m);

				for (unsigned i = 0; i < 8; i++) {


					Vector<real> const Wi = basisRaviartThomas.getColumn(i);
					real const dotProduct = dot(Wi, dPhim);

					for (unsigned j = 0; j < 3; j++) {

						real const Phij = basisPolynomial(j);

						τ.setCoeff(k_index, m, i, j) += weight * dotProduct * Phij;

					}
				}
			}
		}

		for (unsigned m = 0; m < 3; m++)
			for (unsigned i = 0; i < 8; i++)
				for (unsigned j = 0; j < 3; j++)
					τ.setCoeff(k_index, m, i, j) = abs(τ(k_index, m, i, j)) < INTEGRAL_PRECISION ? 0.0 : τ(k_index, m, i, j);


		//for (unsigned m = 0; m < 3; m++) {

		//	Matrix<real> integral(8, 3);

		//	for (unsigned i = 0; i < 8; i++)
		//		for (unsigned j = 0; j < 3; j++)
		//			integral(i, j) = τ(k_index, m, i, j);

		//	std::cout << integral << std::endl;


		//}

	}

};
real solver::compute_δ(t_pointer const K, unsigned const j, unsigned const El, unsigned const m) {


	// Quadrature weights and points on reference segment [-1,1]
	gauss_quadrature_1D const quad(quadrature_order);
	unsigned const num_quad_points = quad.number_of_points;

	Matrix <real> normals(2, 3);
	Matrix<real> parametrization(2, 3);
	Matrix<real> parametrizationDerivative(2, 8);
	Matrix<real> basisRaviartThomas(2, 8);
	Vector<real> basisPolynomial(3);


	evaluate_edge_normal(normals);





		t_pointer const K = mesh->get_triangle(k);
		unsigned const k_index = K->index;


		real orientations[3];

		//for (unsigned e = 0; e < 3; e++)
		//	orientations[e] = edgeOrientation(k_index, 0, e, 0);

		for (unsigned e = 0; e < 3; e++) {

			if (K->edges[e]->marker == E_MARKER::NEUMANN || K->edges[e]->marker == E_MARKER::DIRICHLET)
				orientations[e] = 1.0;
			else
				orientations[e] = edgeOrientation(k_index, 0, e, 0);

		}



		for (unsigned n = 0; n < num_quad_points; n++) {

			for (unsigned El = 0; El < 3; El++) {


				real const orientation = orientations[El];


				real const a = (real) 0.0;
				real const b = (real)El != 0 ? 1.0 : sqrt(2.0);

				real const c = (real)(b - a) / 2.0;
				real const d = (real)(b + a) / 2.0;

				real const x = (real)quad.points[n] * c + d;
				real const w = (real)quad.weigths[n] * c;

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
				evaluate_polynomial_basis(s, t, basisPolynomial);


				real const tC = upwindConcentration(s, t, K, normal, El, ξ_prev);
				//real const tC = upwindConcentration(s, t, K, normal, El, ξ_prev) / K->edges[El]->length();



				//real const coeff = El == 0 ? sqrt(2.0) : 1.0;


				for (unsigned m = 0; m < 3; m++) {


					real const Phim = basisPolynomial(m);

					for (unsigned j = 0; j < 8; j++) {

						Vector<real> const Wj = basisRaviartThomas.getColumn(j);
						real const dotProduct = dot(Wj, normal);
						//real const dotProduct = coeff * dot(Wj, normal);
						//real const dotProduct = orientation * dot(Wj, normal);

						δ.setCoeff(k_index, m, El, j) = δ(k_index, m, El, j) + w * tC * dotProduct * Phim * drNorm;

					}
				}
			}
		}

		for (unsigned El = 0; El < 3; El++)
			for (unsigned j = 0; j < 8; j++)
				for (unsigned m = 0; m < 3; m++)
					δ.setCoeff(k_index, m, El, j) = abs(δ(k_index, m, El, j)) < INTEGRAL_PRECISION ? 0.0 : δ(k_index, m, El, j);


};
void solver::compute_γ() {


	for (unsigned k = 0; k < nk; k++) {


		t_pointer const K = mesh->get_triangle(k);
		unsigned const k_index = K->index;


		for (unsigned m = 0; m < 3; m++) {

			for (unsigned j = 0; j < 8; j++) {


				real val1 = 0.0;
				real val2 = 0.0;

				for (unsigned l = 0; l < 3; l++)
					val1 += δ(k_index, m, l, j);

				for (unsigned l = 0; l < 3; l++)
					val2 += τ(k_index, m, j, l) * ξ_prev(k_index, 0, l, 0);


				γ.setCoeff(k_index, 0, m, j) = val1 - val2;

			}
		}
	}

};

*/







