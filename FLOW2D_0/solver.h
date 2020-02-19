#pragma once


#include "coefficient_matrix.h"
#include "integration.h"
#include "mesh.h"
#include "matrix.h"
#include <omp.h>

#include <Eigen/Sparse>
#include <Eigen/Dense>


typedef double real;

real const INTEGRAL_PRECISION = 1e-11;

typedef Eigen::SparseMatrix<real>			SparseMatrix;
typedef Eigen::VectorXd						DenseVector;


enum SCHEME { IMPLICIT, EXPLICIT };
enum PRINT { PRESSURE, CONCENTRATION, TRACE_PRESSURE };


template<unsigned i>
constexpr unsigned const get_number_of_quadrature_points_edge() {

	switch (i) {

	case 2:
		return 2;
	case 3:
		return 3;
	case 4:
		return 4;
	case 5:
		return 5;
	case 6:
		return 6;
	case 7:
		return 7;
	case 8:
		return 8;
	case 9:
		return 9;
	case 11:
		return 11;
	case 13:
		return 13;

	default:
		return 0;

	}

};
template<unsigned i>
constexpr unsigned const get_number_of_quadrature_points_triangle() {

	switch (i) {

	case 2:
		return 3;
	case 3:
		return 6;
	case 4:
		return 7;
	case 5:
		return 7;
	case 6:
		return 12;
	case 7:
		return 13;
	case 8:
		return 19;
	case 9:
		return 19;
	case 11:
		return 28;
	case 13:
		return 37;

	default:
		return 0;

	}

};

template <unsigned QuadraturePrecision>
class solver {


public:

	solver(Mesh & m, unsigned nt_0, double dt_0);
	~solver();

	void getSolution();

	void setTimeStep(double _dt) { this->dt = _dt; };
	void setTimeLevel(int _nt) { this->nt = _nt; };

	void exportSolution(std::ofstream & txtFile);

	void compute_error(std::ofstream & txtFile);



private:




	// Triangulation mesh
	Mesh * mesh;

	// nk = number of elements in the mesh, ne = number of edges in the mesh
	unsigned const nk;
	unsigned const ne;


	//Eigen::MatrixXd * π2 = new Eigen::MatrixXd[50];

	Eigen::VectorXd		p;
	Eigen::VectorXd		tp;
	Eigen::VectorXd		c;
	CoeffMatrix1D<3>	v;

	// Qunatities on time levels. Used in theta-scheme current quantites on time level 'n', last iteration in l 'prev'
	Eigen::VectorXd p_n;
	Eigen::VectorXd p_prev;
	Eigen::VectorXd c_n;
	Eigen::VectorXd c_prev;
	Eigen::VectorXd rkFp;
	Eigen::VectorXd rkFp_n;
	Eigen::VectorXd rkFc;
	Eigen::VectorXd rkFc_n;


	// Integral (geometric) coeffiecients
	CoeffMatrix2D<3, 3>		α;
	real *					σ;
	CoeffMatrix1D<3>		λ;

	static unsigned const NumberOfQuadraturePointsEdge		= get_number_of_quadrature_points_edge<QuadraturePrecision>();
	static unsigned const NumberOfQuadraturePointsTriangle	= get_number_of_quadrature_points_triangle<QuadraturePrecision>();


	unsigned * ElementIndeces = NULL;



	real * viscosities;
	real * porosities;

	real * betas;
	real * betas_prev;

	// Trace pressure system matrix
	SparseMatrix R;
	SparseMatrix M;
	DenseVector V;

	// Pressure system matrix
	SparseMatrix iD;
	SparseMatrix H;
	DenseVector G;


	// Solver for the computation of LU factorization
	Eigen::SparseLU<SparseMatrix>	sparseLUsolver_TracePressureSystem;

	DenseVector pressureSystemRhs;
	DenseVector traceSystemRhs;

	// Time level 'nt' 
	int nt;
	// Time increment 'dt'
	double dt;


	void initializeValues();
	void computeBetas();

	bool stopCriterion();
	void concentrationCorrection();


	void computePressureEquation();

	void computeTracePressures();
	void computeVelocities();
	void updateConcentrations_explicit();

	real upwindConcentration(t_pointer const & K, unsigned const El);


	void assembleR();
	void assembleM();
	void assembleV();

	void assembleInverseD();
	void assembleH();
	void assembleG();

	void assemble_α();
	void assemble_σ();
	void assemble_λ();

};











template<unsigned QuadraturePrecision>
solver<QuadraturePrecision>::solver(Mesh & m, unsigned nt_0, double dt_0)

	: nk(m.get_number_of_triangles()), ne(m.get_number_of_edges())

{


	mesh = &m;

	nt = nt_0;
	dt = dt_0;


	/*****************************************************************************/
	/*                                                                           */
	/*    - Memory allocation for the physical quantities   			         */
	/*                                                                           */
	/*****************************************************************************/

	viscosities = new double[nk];
	porosities	= new double[nk];


	/*****************************************************************************/
	/*                                                                           */
	/*    - Memory allocation for the beta coefficient. This coefficient         */
	/*      links model equations and EOS                                        */
	/*                                                                           */
	/*****************************************************************************/

	betas		= new double[nk];
	betas_prev	= new double[nk];


	/*****************************************************************************/
	/*                                                                           */
	/*    - Memory allocation for the unknowns : p - Internal Pressures		     */
	/*										   : tp - Trace pressures		     */
	/*										   : c - Concentrations     	     */
	/*										   : v - Velocities				     */
	/*                                                                           */
	/*****************************************************************************/

	p	.resize(nk);
	tp	.resize(ne);
	c	.resize(nk);
	v	.setNumberOfElements(nk);

	v.setZero();

	/*****************************************************************************/
	/*                                                                           */
	/*    - Memory allocation for the auxillary variables used in			     */
	/*      time discretization                                                  */
	/*    - 'n'    = current time level                                          */
	/*      'prev' = previous iteration in the inner loop                        */
	/*                                                                           */
	/*****************************************************************************/

	p_n		.resize(nk);
	p_prev	.resize(nk);
	c_n		.resize(nk);
	c_prev	.resize(nk);
	rkFp	.resize(nk);
	rkFp_n	.resize(nk);
	rkFc	.resize(nk);
	rkFc_n	.resize(nk);

	/*****************************************************************************/
	/*                                                                           */
	/*    - Memory allocation for the integral coefficients				         */
	/*                                                                           */
	/*****************************************************************************/

	α.setNumberOfElements(nk);
	λ.setNumberOfElements(nk);
	σ = new real[nk];
	

	α.setZero();
	λ.setZero();

	/*****************************************************************************/
	/*                                                                           */
	/*    - Memory allocation for the pressure system					         */
	/*                                                                           */
	/*****************************************************************************/

	pressureSystemRhs		.resize(ne);
	traceSystemRhs			.resize(ne);

	iD	.resize(nk, nk);
	H	.resize(nk, ne);
	G	.resize(nk);

	M	.resize(ne, ne);
	R	.resize(ne, nk);
	V	.resize(ne);


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

	/*****************************************************************************/
	/*                                                                           */
	/*    - Initilize constant matrices										     */
	/*    - Assembly of the trace pressure system and compute					 */
	/*		its LU decomposition												 */
	/*                                                                           */
	/*****************************************************************************/

	assembleR();
	assembleM();


	/*****************************************************************************/
	/*                                                                           */
	/*    - Memory allocation for precomputed values of basis function			 */
	/* 		in quadrature points												 */
	/*                                                                           */
	/*****************************************************************************/

	ElementIndeces = new unsigned[nk];

	/*****************************************************************************/
	/*                                                                           */
	/*    - Precomputation of values defined on the elements				     */
	/*                                                                           */
	/*****************************************************************************/

	for (unsigned k = 0; k < nk; k++) {


		t_pointer const K = mesh->get_triangle(k);
		unsigned const k_index = K->index;

		ElementIndeces[k] = k_index;

	}


};
template<unsigned QuadraturePrecision>
solver<QuadraturePrecision>::~solver() {


	delete[] ElementIndeces;

	delete[] σ;

	delete[] viscosities;
	delete[] porosities;

	delete[] betas;
	delete[] betas_prev;

};




template<unsigned QuadraturePrecision>
void solver<QuadraturePrecision>::initializeValues() {


	double const time = nt * dt;


	for (unsigned k = 0; k < nk; k++) {


		t_pointer const K = mesh->get_triangle(k);
		unsigned const k_index = K->index;

		real const area = K->area();


		p[k_index] = integrate_triangle(K, time, barenblatt) / area;
		c[k_index] = integrate_triangle(K, time, barenblatt) / area;

		/*****************************************************************************/
		/*                                                                           */
		/*    - Mean values of the viscosity, porosity on each element			     */
		/*                                                                           */
		/*****************************************************************************/

		viscosities[k_index]	= integrate_triangle(K, viscosity) / area;
		porosities[k_index]		= integrate_triangle(K, porosity) / area;



		/*****************************************************************************/
		/*                                                                           */
		/*    - Auxilary variebles												     */
		/*                                                                           */
		/*****************************************************************************/

		rkFc[k_index] = 0.0;
		rkFc_n[k_index] = rkFc[k_index];

		rkFp[k_index] = 0.0;
		rkFp_n[k_index] = rkFp[k_index];

	}

	//std::cout << p << std::endl;
	
};
template<unsigned QuadraturePrecision>
void solver<QuadraturePrecision>::computeBetas() {


	t_pointer K = NULL;

	for (unsigned k = 0; k < nk; k++) {

		//K = mesh->getElement(k);

		betas[k] = 1.0;// Deos(ξ(K, 0, 0));

	}
	//	betas[i] = (eos(concentrations[i]) - eos(concentrations_temporary[i])) / (concentrations[i] - concentrations_temporary[i]);

};


template<unsigned QuadraturePrecision>
bool solver<QuadraturePrecision>::stopCriterion() {


	real ErrorPressure = 0.0;
	real NormPressure = 0.0;

	real ErrorConcentration = 0.0;
	real NormConcentration = 0.0;

	real ErrorBeta = 0.0;
	real NormBeta = 0.0;

	for (unsigned k = 0; k < nk; k++) {


		t_pointer const	K = mesh->get_triangle(k);
		unsigned const k_index = K->index;

		real const IntegralError = sqr(p[k_index] - p_prev[k_index]);
		real const IntegralNorm = sqr(p[k_index]);

		ErrorPressure	+= IntegralError;
		NormPressure	+= IntegralNorm;

	}

	real const eP = ErrorPressure / NormPressure;

	if (eP > TOL)
		return false;


	for (unsigned k = 0; k < nk; k++) {


		t_pointer const	K = mesh->get_triangle(k);
		unsigned const k_index = K->index;

		real const IntegralError = sqr(c[k_index] - c_prev[k_index]);
		real const IntegralNorm = sqr(c[k_index]);

		ErrorConcentration	+= IntegralError;
		NormConcentration	+= IntegralNorm;

	}

	real const eC = ErrorConcentration / NormConcentration;

	if (eC > TOL)
		return false;

	
	/*for (unsigned k = 0; k < nk; k++) {

		sB1 += sqr(betas[k] - betas_prev[k]);
		sB2 += sqr(betas[k]);

	}

	double const val_B = sB1 / sB2;

	if (val_B > TOL)
		return false;*/

	return true;


};
template<unsigned QuadraturePrecision>
void solver<QuadraturePrecision>::concentrationCorrection() {

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



template<unsigned QuadraturePrecision>
void solver<QuadraturePrecision>::computeTracePressures() {


	assembleV();

	//traceSystemRhs	= R * p - V;
	tp				= sparseLUsolver_TracePressureSystem.solve(R * p - V);

	//std::cout << tp << std::endl;
	//
	//std::ofstream txtFile;
	//
	//txtFile.open("C:\\Users\\pgali\\Desktop\\flow2d\\tp.txt");
	//
	//for (unsigned e = 0; e < ne; e++) {
	//
	//
	//	e_pointer const E = mesh->get_edge(e);
	//
	//	v_pointer const va = E->a;
	//	v_pointer const vb = E->b;
	//
	//	real const x0 = (real)va->x;
	//	real const y0 = (real)va->y;
	//
	//	real const x1 = (real)vb->x;
	//	real const y1 = (real)vb->y;
	//
	//	txtFile << x0 << " " << y0 << " " << std::setprecision(20) << tp[e] << std::endl;
	//	txtFile << x1 << " " << y1 << " " << std::setprecision(20) << tp[e] << std::endl;
	//	txtFile << std::endl;
	//}
	//
	//txtFile.close();

};
template<unsigned QuadraturePrecision>
void solver<QuadraturePrecision>::computePressureEquation() {


	assembleInverseD();
	assembleH();

	assembleG();
	assembleV();

	
	Eigen::SparseLU<SparseMatrix> const solver(R * iD * H + M);
	//pressureSystemRhs = R * iD * G - V;
	tp = solver.solve(R * iD * G - V);

	p = iD * (G - H * tp);

};

template<unsigned QuadraturePrecision>
void solver<QuadraturePrecision>::computeVelocities() {


	real const time = (nt + 1)* dt;

	for (unsigned k = 0; k < nk; k++) {


		t_pointer const K = mesh->get_triangle(k);
		unsigned const k_index = K->index;


		for (unsigned El = 0; El < 3; El++) {


			e_pointer const E = K->edges[El];

			//unsigned const e_index_local = K->get_edge_index(E);

			if (E->marker == E_MARKER::NEUMANN) {

				//v.setCoeff(k_index, El) = NEUMANN_GAMMA_Q_velocity(E, time);
				v.setCoeff(k_index, El) = NEUMANN_GAMMA_Q_velocity(E, time);
				continue;

			}


			real AlphaP = 0.0;
			real AlphaTP = 0.0;

			for (unsigned j = 0; j < 3; j++)
				AlphaP += α(k_index, El, j);

			for (unsigned l = 0; l < 3; l++)
				AlphaTP += α(k_index, El, l) * tp(K->edges[l]->index);
				//AlphaTP += α(k_index, El, l) * tp(E->index);

			v.setCoeff(k_index, El) = (AlphaP * p[k_index] - AlphaTP) / viscosities[k_index];

		}
	}

};
template<unsigned QuadraturePrecision>
void solver<QuadraturePrecision>::updateConcentrations_explicit() {


	for (unsigned k = 0; k < nk; k++) {


		t_pointer const K = mesh->get_triangle(k);
		unsigned const k_index = ElementIndeces[k];


		real Value = 0.0;

		for (unsigned El = 0; El < 3; El++)
			Value += v(k_index, El) * upwindConcentration(K, El);


		rkFc[k_index] = -Value / (porosities[k_index] * K->area());

		c[k_index] = c_n[k_index] + dt * (1.0 - θ) * rkFc_n[k_index] + dt * θ * rkFc[k_index];

	}

};


template<unsigned QuadraturePrecision>
real solver<QuadraturePrecision>::upwindConcentration(t_pointer const & K, unsigned const El) {


	real const time = (nt + 1) * dt;

	unsigned const k_index = K->index;

	e_pointer const E = K->edges[El];
	E_MARKER const e_marker = E->marker;


	real VelocityDotNormal = 0.0;
	real Concentration = 0.0;

	if (e_marker == E_MARKER::NEUMANN) {


		VelocityDotNormal = NEUMANN_GAMMA_Q_velocity(E, time);

		if (VelocityDotNormal < 0.0)
			return DIRICHLET_GAMMA_Q_concentration(E, time);

	}
	else {
		VelocityDotNormal = v(k_index, El);
	}


	if (VelocityDotNormal >= 0.0) {

		Concentration = c_prev[k_index];

	}
	else {


		if (E->marker == E_MARKER::DIRICHLET)
			return DIRICHLET_GAMMA_P_concentration(E, time);


		unsigned const kn_index = K->neighbors[El]->index;

		Concentration = c_prev[kn_index];
		//VelocityDotNormal = v(kn_index, K->neighbors[El]->get_edge_index(E));

	}

	return Concentration;
	//return Concentration / E->length();;

};


template<unsigned QuadraturePrecision>
void solver<QuadraturePrecision>::assembleR() {


	for (unsigned e = 0; e < ne; e++) {


		e_pointer const E = mesh->get_edge(e);
		unsigned const e_index = E->index;
		E_MARKER const e_marker = E->marker;


		if (e_marker == E_MARKER::DIRICHLET) {

			//continue;

			for (unsigned neighborElement = 0; neighborElement < 2; neighborElement++) {

				t_pointer const K = E->neighbors[neighborElement];

				if (!K)
					continue;

				unsigned const k_index = K->index;

				R.coeffRef(e_index, k_index) = 0.0;

			}

			continue;

		}


		for (unsigned neighborElement = 0; neighborElement < 2; neighborElement++) {

			t_pointer const K = E->neighbors[neighborElement];

			if (!K)
				continue;

			unsigned const dof = LI(K, E);
			unsigned const k_index = K->index;

			real Alpha = 0.0;

			for (unsigned j = 0; j < 3; j++)
				Alpha += α(k_index, dof, j);

			real const Value = Alpha / viscosities[k_index];

			R.coeffRef(e_index, k_index) = abs(Value) < INTEGRAL_PRECISION ? 0.0 : Value;

		}
	}

	//std::cout << R.toDense() << std::endl << std::endl;

};
template<unsigned QuadraturePrecision>
void solver<QuadraturePrecision>::assembleM() {


	std::vector<Eigen::Triplet<real>> triplet;

	M.setZero();

	for (unsigned e = 0; e < ne; e++) {


		e_pointer const E = mesh->get_edge(e);
		unsigned const e_index = E->index;
		E_MARKER const e_marker = E->marker;


		if (e_marker == E_MARKER::DIRICHLET) {

			Eigen::Triplet<real> const T(e_index, e_index, -1.0);
			triplet.push_back(T);

			M.coeffRef(e_index, e_index) = -1.0;

			continue;

		}

		for (unsigned neighborElement = 0; neighborElement < 2; neighborElement++) {


			t_pointer const K = E->neighbors[neighborElement];

			if (!K)
				continue;


			unsigned const dof = LI(K, E);
			unsigned const k_index = K->index;


			for (unsigned El = 0; El < 3; El++) {


				e_pointer const E_local = K->edges[El];

				unsigned const e_local_index_local = K->get_edge_index(E_local);	// Local index of local edge ==== El
				unsigned const e_local_index_global = E_local->index;				// Global index of local edge


				real const Value = α(k_index, dof, e_local_index_local) / viscosities[k_index];

				Eigen::Triplet<real> const T(e_index, e_local_index_global, Value);
				triplet.push_back(T);

				M.coeffRef(e_index, e_local_index_global) += Value;
				

			}
		}
	}

	//std::cout << M.toDense() << std::endl << std::endl << std::endl;

	SparseMatrix tracePressureSystem_LU(ne, ne);

	tracePressureSystem_LU.setFromTriplets(triplet.begin(), triplet.end());

	//std::cout << (M- tracePressureSystem_LU).toDense() << std::endl << std::endl << std::endl << std::endl;
	//std::cout << tracePressureSystem_LU.toDense() << std::endl << std::endl << std::endl;

	sparseLUsolver_TracePressureSystem.compute(tracePressureSystem_LU);


};
template<unsigned QuadraturePrecision>
void solver<QuadraturePrecision>::assembleV() {


	real const time = (nt + 1) * dt;

	for (unsigned e = 0; e < ne; e++) {


		e_pointer const E = mesh->get_edge(e);
		unsigned const e_index = E->index;
		E_MARKER const marker = E->marker;


		V[e_index] = 0.0;

		if (marker == E_MARKER::NEUMANN)
			V[e_index] = NEUMANN_GAMMA_Q_velocity(E, time);
		else if (marker == E_MARKER::DIRICHLET)
			V[e_index] = DIRICHLET_GAMMA_P_pressure(E, time);

	}

	//std::cout << V << std::endl;
};


template<unsigned QuadraturePrecision>
void solver<QuadraturePrecision>::assembleInverseD() {


	assemble_σ();

	real const TimeCoeff = θ * dt;


	for (unsigned k = 0; k < nk; k++) {


		t_pointer const K = mesh->get_triangle(k);
		unsigned const k_index = K->index;


		real const DiagValue = 1.0 - TimeCoeff * σ[k_index];

		iD.coeffRef(k_index, k_index) = 1.0 / DiagValue;

	}

};
template<unsigned QuadraturePrecision>
void solver<QuadraturePrecision>::assembleH() {


	assemble_λ();

	real const TimeCoeff = -θ * dt;


	for (unsigned k = 0; k < nk; k++) {


		t_pointer const K = mesh->get_triangle(k);
		unsigned const k_index = K->index;


		for (unsigned El = 0; El < 3; El++) {

			e_pointer const E = K->edges[El];
			unsigned const e_index = E->index;

			H.coeffRef(k_index, e_index) = TimeCoeff * λ(k_index, El);
			//H.coeffRef(k_index, e_index) = TimeCoeff * λ(k_index, K->get_edge_index(E));

		}
	}

};
template<unsigned QuadraturePrecision>
void solver<QuadraturePrecision>::assembleG() {


	real const TimeCoeff = dt * (1.0 - θ);

	for (unsigned k = 0; k < nk; k++) {

		unsigned const k_index = ElementIndeces[k];

		G[k_index] = p_n[k_index] + TimeCoeff * rkFp_n[k_index];

	}

};


template<unsigned QuadraturePrecision>
void solver<QuadraturePrecision>::assemble_α() {


	quadrature_triangle const QuadratureOnTriangle(QuadraturePrecision);

	Eigen::MatrixXd Integral(3, 3);
	Eigen::MatrixXd BasisRaviartThomas(2, 3);
	Eigen::MatrixXd JF(2, 2);

	for (unsigned k = 0; k < nk; k++) {


		t_pointer const K = mesh->get_triangle(k);
		unsigned const k_index = K->index;

		v_pointer const va = K->vertices[0];
		v_pointer const vb = K->vertices[1];
		v_pointer const vc = K->vertices[2];

		real const x0 = (real)va->x;
		real const y0 = (real)va->y;

		real const x1 = (real)vb->x;
		real const y1 = (real)vb->y;

		real const x2 = (real)vc->x;
		real const y2 = (real)vc->y;

		real const detJF = (real)abs((x1 - x0)*(y2 - y0) - (x2 - x0)*(y1 - y0));

		JF.coeffRef(0, 0) = x1 - x0;
		JF.coeffRef(0, 1) = x2 - x0;
		JF.coeffRef(1, 0) = y1 - y0;
		JF.coeffRef(1, 1) = y2 - y0;


		Integral.setZero();

		real const X[3] = { x0,x1,x2 };
		real const Y[3] = { y0,y1,y2 };

		/*
		for (unsigned n = 0; n < 3; n++) {


			real const s = (real) X[n];
			real const t = (real) Y[n];
			real const w = (real) 1.0 / 3.0;

			// Corresponding coordinates on the element K
			real const x = x0 + JF(0, 0) * s + JF(0, 1) * t;
			real const y = y0 + JF(1, 0) * s + JF(1, 1) * t;

			// Get the inverse of the permeability tensor
			Eigen::MatrixXd iK(2, 2);

			//K << 1.0, 0.0,
			//	 0.0, 1.0;

			iK.coeffRef(0, 0) = 1.0;
			iK.coeffRef(0, 1) = 0.0;
			iK.coeffRef(1, 0) = 0.0;
			iK.coeffRef(1, 1) = 1.0;

			iK = iK.inverse();

			evaluate_raviartthomas_basis(s, t, BasisRaviartThomas);


			for (unsigned i = 0; i < 3; i++) {


				Eigen::VectorXd const JFWi = JF * (BasisRaviartThomas.col(i));

				for (unsigned j = 0; j < 3; j++) {


					Eigen::VectorXd const JFWj = JF * (BasisRaviartThomas.col(j));

					Integral.coeffRef(i, j) += w * JFWi.dot(iK * JFWj);

				}
			}
		}
		*/

		
		for (unsigned n = 0; n < NumberOfQuadraturePointsTriangle; n++) {


			real const s = (real)QuadratureOnTriangle.points_x[n];
			real const t = (real)QuadratureOnTriangle.points_y[n];
			real const w = (real)0.5 * QuadratureOnTriangle.weights[n];

			// Corresponding coordinates on the element K
			real const x = x0 + JF(0, 0) * s + JF(0, 1) * t;
			real const y = y0 + JF(1, 0) * s + JF(1, 1) * t;

			// Get the inverse of the permeability tensor
			Eigen::MatrixXd iK(2, 2);

			//K << 1.0, 0.0,
			//	 0.0, 1.0;

			iK.coeffRef(0, 0) = 1.0;
			iK.coeffRef(0, 1) = 0.0;
			iK.coeffRef(1, 0) = 0.0;
			iK.coeffRef(1, 1) = 1.0;

			iK = iK.inverse();

			evaluate_raviartthomas_basis(s, t, BasisRaviartThomas);


			for (unsigned i = 0; i < 3; i++) {


				Eigen::VectorXd const JFWi = JF * (BasisRaviartThomas.col(i));

				for (unsigned j = 0; j < 3; j++) {


					Eigen::VectorXd const JFWj = JF * (BasisRaviartThomas.col(j));

					Integral.coeffRef(i, j) += w * JFWi.dot(iK * JFWj);

				}
			}
		}
		

		//std::cout << Integral << std::endl << std::endl;

		Integral = Integral / detJF;

		//Integral = K->area() * Integral / detJF;

		//std::cout << Integral << std::endl << std::endl;

		for (unsigned i = 0; i < 3; i++)
			for (unsigned j = 0; j < 3; j++)
				Integral.coeffRef(i, j) = abs(Integral(i, j)) < INTEGRAL_PRECISION ? 0.0 : Integral(i, j);

		//real const r0 = Integral.row(0).sum();
		//real const r1 = Integral.row(1).sum();
		//real const r2 = Integral.row(2).sum();
		//
		//std::cout << Integral << std::endl << std::endl;
		//
		//Integral.setZero();
		//
		//Integral(0, 0) = r0;
		//Integral(1, 1) = r1;
		//Integral(2, 2) = r2;


		//std::cout << Integral << std::endl << std::endl;
		//Eigen::MatrixXd iii = Integral;

		Integral = Integral.inverse();

		//std::cout << Integral* iii << std::endl << std::endl;

		for (unsigned i = 0; i < 3; i++)
			for (unsigned j = 0; j < 3; j++)
				α.setCoeff(k_index, i, j) = abs(Integral(i, j)) < INTEGRAL_PRECISION ? 0.0 : Integral(i, j);

	}

};
template<unsigned QuadraturePrecision>
void solver<QuadraturePrecision>::assemble_σ() {



	for (unsigned k = 0; k < nk; k++) {

		t_pointer const K = mesh->get_triangle(k);
		unsigned const k_index = ElementIndeces[k];
		real const Area = K->area();

		real const ElementCoeff = betas[k_index] / (porosities[k_index] * viscosities[k_index] * Area);


		real Val = 0.0;

		for (unsigned El = 0; El < 3; El++) {


			real Alpha = 0.0;

			for (unsigned j = 0; j < 3; j++)
				Alpha += α(k_index, El, j);

			real const tC = upwindConcentration(K, El);

			Val += tC * Alpha;

		}

		σ[k_index] = -ElementCoeff * Val;

	}

};
template<unsigned QuadraturePrecision>
void solver<QuadraturePrecision>::assemble_λ() {


	for (unsigned k = 0; k < nk; k++) {

		t_pointer const K = mesh->get_triangle(k);
		unsigned const k_index = ElementIndeces[k];
		real const Area = K->area();

		real const ElementCoeff = betas[k_index] / (porosities[k_index] * viscosities[k_index] * Area);


		for (unsigned l = 0; l < 3; l++) {

			real Value = 0.0;

			for (unsigned El = 0; El < 3; El++)
				Value += upwindConcentration(K, El) * α(k_index, El, l);

			λ.setCoeff(k_index, l) = ElementCoeff * Value;

		}
	}

};



template<unsigned QuadraturePrecision>
void solver<QuadraturePrecision>::getSolution() {


	// Compute initial values for iterations
	computeTracePressures();
	computeVelocities();
	computeBetas();

	std::ofstream txtFile;
	txtFile.open("C:\\Users\\pgali\\Desktop\\flow2d\\v.txt");
	for (unsigned k = 0; k < nk; k++) {

		for (unsigned j = 0; j < 3; j++)
			txtFile << v(mesh->get_triangle(k)->index, j) << " ";

		txtFile << std::endl;
	}
	txtFile.close();

	// Set iteration level l := 0
	p_n		= p;
	p_prev	= p;

	c_n		= c;
	c_prev	= c;

	rkFc_n = rkFc;
	rkFp_n = rkFp;

	//hard_copy(betas_prev, betas, nk);
	for (unsigned k = 0; k < nk; k++)
		betas_prev[k] = betas[k];



	unsigned counter = 0;

	while (counter < MAX_IT) {


		computePressureEquation();
		computeVelocities();

		updateConcentrations_explicit();

		//concentrationCorrection();

		computeBetas();


		if (stopCriterion())
			break;



		// Set new iteration level	l := l+1
		p_prev = p;
		c_prev = c;

		//std::ofstream txtFile;
		//txtFile.open("C:\\Users\\pgali\\Desktop\\output_pressure.txt");
		//exportSolution(txtFile);
		//txtFile.close();


		//hard_copy(betas_prev, betas, nk);
		for (unsigned k = 0; k < nk; k++)
			betas_prev[k] = betas[k];


		counter++;

	}


	//std::ofstream txtFile;
	//txtFile.open("C:\\Users\\pgali\\Desktop\\flow2d\\v.txt");
	//for (unsigned k = 0; k < nk; k++) {

	//	for (unsigned j = 0; j < 3; j++)
	//		txtFile << v(mesh->get_triangle(k)->index, j) << " ";

	//	txtFile << std::endl;
	//}
	//txtFile.close();


	std::cout << counter << std::endl;

	assemble_λ();
	assemble_σ();

	// This is needed for the next time level
	for (unsigned k = 0; k < nk; k++) {


		t_pointer const K = mesh->get_triangle(k);
		unsigned const k_index = ElementIndeces[k];

		
		real const val1 = σ[k_index] * p[k_index];
		real val2 = 0.0;

		for (unsigned El = 0; El < 3; El++)
			val2 += λ(k_index, El) * tp(K->edges[El]->index);

		rkFp[k_index] = val1 + val2;

	}


	if (nt % 1000 == 0)
		std::cout << nt << " - Iterations : " << counter << std::endl;

};
template<unsigned QuadraturePrecision>
void solver<QuadraturePrecision>::exportSolution(std::ofstream & txtFile) {

	//Eigen::MatrixXd JF(2, 2);
	//real const time = nt * dt;


	for (unsigned k = 0; k < nk; k++) {


		t_pointer const	K = mesh->get_triangle(k);

		unsigned const k_index = K->index;

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


		//JF.coeffRef(0, 0) = x1 - x0;
		//JF.coeffRef(0, 1) = x2 - x0;
		//JF.coeffRef(1, 0) = y1 - y0;
		//JF.coeffRef(1, 1) = y2 - y0;

		//real const X0 = x0 + JF(0, 0) * 0.0 + JF(0, 1) * 0.0;
		//real const Y0 = y0 + JF(1, 0) * 0.0 + JF(1, 1) * 0.0;

		//real const X1 = x0 + JF(0, 0) * 1.0 + JF(0, 1) * 0.0;
		//real const Y1 = y0 + JF(1, 0) * 1.0 + JF(1, 1) * 0.0;

		//real const X2 = x0 + JF(0, 0) * 0.0 + JF(0, 1) * 1.0;
		//real const Y2 = y0 + JF(1, 0) * 0.0 + JF(1, 1) * 1.0;

		//real const XX[3] = { X0,X1,X2 };
		//real const YY[3] = { Y0,Y1,Y2 };

		for (unsigned i = 0; i < 3; i++) {

			//// Pressures
			real const Value = p[k_index];
			//// Concentrations
			//real const Value = c[k_index];


			//real const Value = barenblatt(XX[i], YY[i], time);


			txtFile << std::setprecision(20) << x[i] << "\t" << y[i] << "\t" << Value << std::endl;
			//txtFile << std::setprecision(20) << x[i] << "	" << y[i] << "	" << value << "	" << index << std::endl;

		}


		txtFile << std::endl;

	}

};



template<unsigned QuadraturePrecision>
void solver<QuadraturePrecision>::compute_error(std::ofstream & txtFile) {


	real const time = nt * dt;

	quadrature_triangle const QuadratureOnTriangle(9);
	unsigned const NumberOfQuadraturePoints = QuadratureOnTriangle.NumberOfPoints;

	Eigen::MatrixXd JF(2, 2);

	real errorL1 = 0.0;
	real errorL2 = 0.0;
	real errorMAX = 0.0;

	for (unsigned k = 0; k < nk; k++) {


		t_pointer const	K = mesh->get_triangle(k);
		unsigned const k_index = K->index;

		v_pointer const va = K->vertices[0];
		v_pointer const vb = K->vertices[1];
		v_pointer const vc = K->vertices[2];

		real const x0 = (real)va->x;
		real const y0 = (real)va->y;

		real const x1 = (real)vb->x;
		real const y1 = (real)vb->y;

		real const x2 = (real)vc->x;
		real const y2 = (real)vc->y;

		real const detJF = (real)abs((x1 - x0)*(y2 - y0) - (x2 - x0)*(y1 - y0));

		JF.coeffRef(0, 0) = x1 - x0;
		JF.coeffRef(0, 1) = x2 - x0;
		JF.coeffRef(1, 0) = y1 - y0;
		JF.coeffRef(1, 1) = y2 - y0;


		real Integral = 0.0;

		for (unsigned n = 0; n < NumberOfQuadraturePoints; n++) {


			real const s = (real) QuadratureOnTriangle.points_x[n];
			real const t = (real) QuadratureOnTriangle.points_y[n];
			real const w = (real) 0.5 * QuadratureOnTriangle.weights[n];

			real const X = x0 + JF(0, 0) * s + JF(0, 1) * t;
			real const Y = y0 + JF(1, 0) * s + JF(1, 1) * t;

			real const PressureK = p[k_index];


			Integral += w * sqr(PressureK - barenblatt(X, Y, time));

		}


		errorL2 += detJF * Integral;

	}

	//txtFile << "#L1 L2 MAX" << std::endl;
	//txtFile << std::setprecision(20) << errorL1 << std::endl;
	//txtFile << std::setprecision(20) << sqrt(errorL2) << std::endl;
	//txtFile << std::setprecision(20) << errorMAX << std::endl;

	//std::cout << "Error L1 : " << errorL1 << std::endl;
	std::cout << "Error L2 : " << sqrt(errorL2) << std::endl;
	//std::cout << "Error Max : " << errorMAX << std::endl;


};