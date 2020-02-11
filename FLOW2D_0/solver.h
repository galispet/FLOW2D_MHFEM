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

	Eigen::VectorXd p;
	Eigen::MatrixXd tp;
	Eigen::VectorXd c;
	Eigen::MatrixXd v;

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


	// Quadrature points on the edges of the physical triangle
	CoeffMatrix2D<3, NumberOfQuadraturePointsEdge>	 QuadraturePoints_Edge_x;
	CoeffMatrix2D<3, NumberOfQuadraturePointsEdge>	 QuadraturePoints_Edge_y;


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
	DenseVector iD;
	SparseMatrix H;
	DenseVector G;


	// Solver for the computation of LU factorization
	Eigen::SparseLU<SparseMatrix>	sparseLUsolver_TracePressureSystem;
	SparseMatrix					internalPressureSystem;

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

	real upwindConcentration(t_pointer const & K, unsigned const El, unsigned const n);
	real upwindConcentration_limiter(t_pointer const & K, unsigned const El, unsigned const n);


	void assembleR();
	void assembleM();
	void assembleV();

	void assembleInverseD();
	void assembleH();
	void assembleG();


	void assemblePressureSystem();

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
	tp	.resize(nk, 3);
	c	.resize(nk);
	v	.resize(nk);

	σ = new real[nk];
	λ.setNumberOfElements(nk);

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
	α.setZero();

	/*****************************************************************************/
	/*                                                                           */
	/*    - Memory allocation for the pressure system					         */
	/*                                                                           */
	/*****************************************************************************/

	internalPressureSystem	.resize(ne, ne);
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

	gauss_quadrature_1D const GaussQuadratureOnEdge(QuadraturePrecision);

	QuadraturePoints_RaviartThomasBasis									.setNumberOfElements(NumberOfQuadraturePointsEdge);
	QuadraturePoints_PhysicalNormalDotPhysicalRaviartThomasBasis		.setNumberOfElements(nk);

	QuadraturePoints_Edge_x												.setNumberOfElements(nk);
	QuadraturePoints_Edge_y												.setNumberOfElements(nk);

	AffineMappingMatrixDeterminant = new real[nk];

	ElementIndeces = new unsigned[nk];



	QuadraturePoints_RaviartThomasBasis									.setZero();
	QuadraturePoints_PhysicalNormalDotPhysicalRaviartThomasBasis		.setZero();

	QuadraturePoints_Edge_x												.setZero();
	QuadraturePoints_Edge_y												.setZero();



	Eigen::MatrixXd JF(2, 2);
	Eigen::MatrixXd ReferenceNormals(2, 3);
	Eigen::VectorXd Parametrization(2);
	Eigen::MatrixXd BasisRaviartThomas(2, 3);

	evaluate_edge_normal(ReferenceNormals);



	/*****************************************************************************/
	/*                                                                           */
	/*    - Precomputation of values defined on the elements				     */
	/*                                                                           */
	/*****************************************************************************/

	for (unsigned k = 0; k < nk; k++) {


		t_pointer const K = mesh->get_triangle(k);
		unsigned const k_index = K->index;

		ElementIndeces[k] = k_index;


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


		/*****************************************************************************/
		/*                                                                           */
		/*    - Precomputed values of : Determinant of the affine mapping matrix JF  */
		/*							  : Affine mapping matrix JF				     */
		/*                                                                           */
		/*****************************************************************************/

		AffineMappingMatrixDeterminant[k_index] = detJF;

		JF(0, 0) = x1 - x0;
		JF(0, 1) = x2 - x0;
		JF(1, 0) = y1 - y0;
		JF(1, 1) = y2 - y0;

		Eigen::MatrixXd const itJF = (JF.inverse()).transpose();


		for (unsigned El = 0; El < 3; El++) {


			e_pointer const E = K->edges[El];
			unsigned const e_index = E->index;
			E_MARKER const e_marker = E->marker;

			unsigned const dof0 = LI(K, E, 0);
			unsigned const dof1 = LI(K, E, 1);


			real const a = (real) 0.0;
			real const b = (real)El != 0 ? 1.0 : sqrt(2.0);

			real const c = (real)(b - a) / 2.0;
			real const d = (real)(b + a) / 2.0;


			//Vector<real> const normal = normals.getColumn(El);
			Eigen::Vector2d const PhysicalNormal = (itJF * ReferenceNormals.col(El)) / (itJF * ReferenceNormals.col(El)).norm();



			/*****************************************************************************/
			/*                                                                           */
			/*    - Precomputed values of : Physical normal * Physical RT basis		     */
			/*							  : Quadrature points on each edge of the		 */
			/*							    reference triangle						     */
			/*                                                                           */
			/*    - Used in : Upwind concentration                                       */
			/*                                                                           */
			/*****************************************************************************/

			for (unsigned n = 0; n < NumberOfQuadraturePointsEdge; n++) {


				real const x = (real)GaussQuadratureOnEdge.points[n] * c + d;
				real const w = (real)GaussQuadratureOnEdge.weights[n] * c;


				evaluate_edge_parametrization(x, El, Parametrization);

				real const s = Parametrization(0);
				real const t = Parametrization(1);

				real const X = x0 + JF(0, 0) * s + JF(0, 1) * t;
				real const Y = y0 + JF(1, 0) * s + JF(1, 1) * t;

				evaluate_raviartthomas_basis(s, t, BasisRaviartThomas);


				QuadraturePoints_Edge_x.setCoeff(k_index, El, n) = abs(X) < INTEGRAL_PRECISION ? 0.0 : X;
				QuadraturePoints_Edge_y.setCoeff(k_index, El, n) = abs(Y) < INTEGRAL_PRECISION ? 0.0 : Y;


				for (unsigned j = 0; j < 3; j++)
					QuadraturePoints_PhysicalNormalDotPhysicalRaviartThomasBasis.setCoeff(k_index, El, n, j) = PhysicalNormal.dot(JF * BasisRaviartThomas.col(j)) / detJF;

			}
		}
	}



	/*****************************************************************************/
	/*                                                                           */
	/*    - Precomputation of values defined on the edges					     */
	/*                                                                           */
	/*****************************************************************************/

	for (unsigned El = 0; El < 3; El++) {


		real const a = (real) 0.0;
		real const b = (real)El != 0 ? 1.0 : sqrt(2.0);

		real const c = (real)(b - a) / 2.0;
		real const d = (real)(b + a) / 2.0;


		for (unsigned n = 0; n < NumberOfQuadraturePointsEdge; n++) {


			real const x = (real)GaussQuadratureOnEdge.points[n] * c + d;
			real const w = (real)GaussQuadratureOnEdge.weights[n] * c;


			/*****************************************************************************/
			/*                                                                           */
			/*    - Precomputed values of : Quadrature points on each edge of the	     */
			/*								reference triangle						     */
			/*                                                                           */
			/*****************************************************************************/

			evaluate_edge_parametrization(x, El, Parametrization);

			real const s = Parametrization(0);
			real const t = Parametrization(1);

			evaluate_raviartthomas_basis(s, t, BasisRaviartThomas);

			/*****************************************************************************/
			/*                                                                           */
			/*    - Precomputed values of : RT basis on the reference triangle		     */
			/*                                                                           */
			/*****************************************************************************/

			for (unsigned j = 0; j < 8; j++) {

				Eigen::VectorXd const Wj = BasisRaviartThomas.col(j);

				QuadraturePoints_RaviartThomasBasis.setCoeff(n, El, j, 0) = Wj(0);
				QuadraturePoints_RaviartThomasBasis.setCoeff(n, El, j, 1) = Wj(1);

			}
		}
	}


};
template<unsigned QuadraturePrecision>
solver<QuadraturePrecision>::~solver() {

	delete[] AffineMappingMatrixDeterminant;

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


		v_pointer const va = K->vertices[0];
		v_pointer const vb = K->vertices[1];
		v_pointer const vc = K->vertices[2];

		real const x0 = (real)va->x;
		real const y0 = (real)va->y;

		real const x1 = (real)vb->x;
		real const y1 = (real)vb->y;

		real const x2 = (real)vc->x;
		real const y2 = (real)vc->y;


		p.setCoeff[k_index] = integrate_triangle(K, time, barenblatt) / area;
		c.setCoeff[k_index] = integrate_triangle(K, time, barenblatt) / area;

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


/*
template<unsigned QuadraturePrecision>
bool solver<QuadraturePrecision>::stopCriterion() {


	double sP1 = 0.0;
	double sP2 = 0.0;

	double sC1 = 0.0;
	double sC2 = 0.0;

	double sB1 = 0.0;
	double sB2 = 0.0;

	t_pointer K = NULL;


	for (unsigned k = 0; k < nk; k++) {

		//unsigned const k_index = mesh->get_triangle(k)->index;
		unsigned const k_index = ElementIndeces[k];

		sP1 += sqr(π(k_index, 0) - π_prev(k_index, 0));
		sP2 += sqr(π(k_index, 0));

	}

	double const val_P = sP1 / sP2;

	if (val_P > TOL)
		return false;


	for (unsigned k = 0; k < nk; k++) {

		//unsigned const k_index = mesh->get_triangle(k)->index;
		unsigned const k_index = ElementIndeces[k];

		sC1 += sqr(ξ(k_index, 0) - ξ_prev(k_index, 0));
		sC2 += sqr(ξ(k_index, 0));

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


	
	//for (unsigned k = 0; k < nk; k++) {
	//	K = mesh->getElement(k);
	//	double const a = K->nodes[0]->x;
	//	double const b = K->nodes[1]->x;
	//	quadrature const quad(error_quadRule, a, b, METHOD::GAUSS);
	//	double y = 0.0;
	//	double w = 0.0;
	//	double val = 0.0;
	//	for (unsigned j = 0; j < error_quadRule; j++) {
	//		y = quad.points[j];
	//		w = quad.weights[j];
	//		val = ((π(K, 0, 0)*phi1(y, a, b) + π(K, 0, 1)*phi2(y, a, b)) - (π_prev(K, 0, 0)*phi1(y, a, b) + π_prev(K, 0, 1) * phi2(y, a, b)));
	//		sP1 += w * sqr(val);
	//		sP2 += w * sqr(π(K, 0, 0)*phi1(y, a, b) + π(K, 0, 1)*phi2(y, a, b));
	//	}
	//}
	//double const val_P = sP1 / sP2;
	//if (val_P > TOL)
	//	return false;
	//return true;
	

};*/
template<unsigned QuadraturePrecision>
bool solver<QuadraturePrecision>::stopCriterion() {


	real const time = nt * dt;

	quadrature_triangle const QuadratureOnTriangle(9);
	unsigned const NumberOfQuadraturePoints = QuadratureOnTriangle.NumberOfPoints;

	Eigen::VectorXd BasisPolynomial(3);
	Eigen::MatrixXd JF(2, 2);

	real ErrorPressure = 0.0;
	real NormPressure = 0.0;

	real ErrorConcentration = 0.0;
	real NormConcentration = 0.0;

	real ErrorBeta = 0.0;
	real NormBeta = 0.0;

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


		real IntegralError = 0.0;
		real IntegralNorm = 0.0;

		for (unsigned n = 0; n < NumberOfQuadraturePoints; n++) {


			real const s = (real) QuadratureOnTriangle.points_x[n];
			real const t = (real) QuadratureOnTriangle.points_y[n];
			real const w = (real) 0.5 * QuadratureOnTriangle.weights[n];

			evaluate_polynomial_basis(s, t, BasisPolynomial);


			real Difference = 0.0;
			real Norm = 0.0;

			for (unsigned j = 0; j < 3; j++) {

				Difference	+= (π(k_index, j) - π_prev(k_index, j)) * BasisPolynomial(j);
				Norm		+= π(k_index, j) * BasisPolynomial(j);

			}
				

			IntegralError	+= w * sqr(Difference);
			IntegralNorm	+= w * sqr(Norm);

		}

		ErrorPressure	+= detJF * IntegralError;
		NormPressure	+= detJF * IntegralNorm;

	}

	real const eP = ErrorPressure / NormPressure;

	if (eP > TOL)
		return false;


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


		real IntegralError = 0.0;
		real IntegralNorm = 0.0;

		for (unsigned n = 0; n < NumberOfQuadraturePoints; n++) {


			real const s = (real)QuadratureOnTriangle.points_x[n];
			real const t = (real)QuadratureOnTriangle.points_y[n];
			real const w = (real) 0.5 * QuadratureOnTriangle.weights[n];

			evaluate_polynomial_basis(s, t, BasisPolynomial);


			real Difference = 0.0;
			real Norm = 0.0;

			for (unsigned j = 0; j < 3; j++) {

				Difference	+= (ξ(k_index, j) - ξ_prev(k_index, j)) * BasisPolynomial(j);
				Norm		+= ξ(k_index, j) * BasisPolynomial(j);

			}


			IntegralError += w * sqr(Difference);
			IntegralNorm += w * sqr(Norm);

		}

		ErrorConcentration	+= detJF * IntegralError;
		NormConcentration	+= detJF * IntegralNorm;

	}

	real const eC = ErrorPressure / NormPressure;

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

	traceSystemRhs	= R * p - V;
	tp				= sparseLUsolver_TracePressureSystem.solve(traceSystemRhs);


};
template<unsigned QuadraturePrecision>
void solver<QuadraturePrecision>::computePressureEquation() {


	assembleInverseD();
	assembleH();

	assembleG();
	assembleV();

	
	Eigen::SparseLU<SparseMatrix> const solver(R * iD.asDiagonal() * H + M);
	pressureSystemRhs = R * iD.asDiagonal() * G - V;
	tp = solver.solve(pressureSystemRhs);

	p = iD.asDiagonal() * (G - H * tp);

};


template<unsigned QuadraturePrecision>
void solver<QuadraturePrecision>::computeVelocities() {


	real const time = (nt + 1)* dt;

	v.setZero();

	for (unsigned k = 0; k < nk; k++) {


		t_pointer const K = mesh->get_triangle(k);
		unsigned const k_index = K->index;


		for (unsigned El = 0; El < 3; El++) {


			e_pointer const E = K->edges[El];

			if (E->marker == E_MARKER::NEUMANN) {

				v.setCoeff(k_index, El) = NEUMANN_GAMMA_Q_velocity(E, time);
				continue;

			}

			for (unsigned m = 0; m < 3; m++) {


				real AlphaP		= 0.0;
				real AlphaTP	= 0.0;

				for (unsigned j = 0; j < 3; j++)
					AlphaP += α(k_index, m, j);

				for (unsigned l = 0; l < 3; l++)
					AlphaTP += α(k_index, m, l) * tp(k_index, l);

				v.setCoeff(k_index, m) = (AlphaP * p[k_index] - AlphaTP) / viscosities[k_index];

			}
		}
	}

};

template<unsigned QuadraturePrecision>
void solver<QuadraturePrecision>::updateConcentrations_explicit() {


	for (unsigned k = 0; k < nk; k++) {


		t_pointer const K = mesh->get_triangle(k);
		unsigned const k_index = ElementIndeces[k];


		real const ElementCoeff = 1.0 / (porosities[k_index] * K->area());



		real Val = 0.0;

		for (unsigned El = 0; El < 3; El++) {


			real const tC = upwindConcentration(K, El);
			Val += v(k_index, El) * tC;

		}

		rkFc.setCoeff[k_index] = -Val * ElementCoeff;

		c.setCoeff[k_index] = c_n[k_index] + dt * (θ * rkFc[k_index] + (1.0 - θ) * rkFc_n[k_index]);

	}

};


template<unsigned QuadraturePrecision>
real solver<QuadraturePrecision>::upwindConcentration(t_pointer const & K, unsigned const El, unsigned const n) {


	real const time = (nt + 1) * dt;

	unsigned const k_index = K->index;

	e_pointer const E = K->edges[El];
	E_MARKER const e_marker = E->marker;


	real VelocityDotNormal = 0.0;
	real Concentration = 0.0;

	if (e_marker == E_MARKER::NEUMANN) {


		VelocityDotNormal = NEUMANN_GAMMA_Q_velocity(E, time);

		if (VelocityDotNormal < 0.0) {

			real const X = QuadraturePoints_Edge_x(k_index, El, n);
			real const Y = QuadraturePoints_Edge_y(k_index, El, n);

			return DIRICHLET_GAMMA_Q_concentration(X, Y, time);

		}
	}
	else {

		for (unsigned dof = 0; dof < 2; dof++) {

			unsigned const j = LI(K, E, dof);

			VelocityDotNormal += υ(k_index, j) * QuadraturePoints_PhysicalNormalDotPhysicalRaviartThomasBasis(k_index, El, n, j);

		}
	}

	
	if (VelocityDotNormal >= 0.0) {

		for (unsigned m = 0; m < 3; m++)
			Concentration += ξ_prev(k_index, m) * QuadraturePoints_PolynomialBasis(n, El, m);

		//Concentration = ξ_prev(k_index, 0) * QuadraturePoints_PolynomialBasis(n, El, 0);

	}
	else {


		if (E->marker == E_MARKER::DIRICHLET) {

			real const X = QuadraturePoints_Edge_x(k_index, El, n);
			real const Y = QuadraturePoints_Edge_y(k_index, El, n);

			return DIRICHLET_GAMMA_P_concentration(X, Y, time);

		}


		unsigned const kn_index = K->neighbors[El]->index;
		unsigned const e_index_Kn_loc = K->neighbors[El]->get_edge_index(E);


		for (unsigned m = 0; m < 3; m++)
			Concentration += ξ_prev(kn_index, m) * QuadraturePoints_PolynomialBasis(n, e_index_Kn_loc, m);

		//Concentration = ξ_prev(kn_index, 0) * QuadraturePoints_PolynomialBasis(n, e_index_Kn_loc, 0);

	}
	
	return Concentration;

};
template<unsigned QuadraturePrecision>
real solver<QuadraturePrecision>::upwindConcentration_limiter(t_pointer const & K, unsigned const El, unsigned const n) {


	real const time = (nt + 1) * dt;

	unsigned const k_index = K->index;

	e_pointer const E = K->edges[El];
	E_MARKER const e_marker = E->marker;


	real VelocityDotNormal = 0.0;
	real Concentration = 0.0;

	if (e_marker == E_MARKER::NEUMANN) {


		VelocityDotNormal = NEUMANN_GAMMA_Q_velocity(E, time);

		if (VelocityDotNormal < 0.0) {

			real const X = QuadraturePoints_Edge_x(k_index, El, n);
			real const Y = QuadraturePoints_Edge_y(k_index, El, n);

			return DIRICHLET_GAMMA_Q_concentration(X, Y, time);

		}
	}
	else {

		for (unsigned dof = 0; dof < 2; dof++) {

			unsigned const j = LI(K, E, dof);

			VelocityDotNormal += υ(k_index, j) * QuadraturePoints_PhysicalNormalDotPhysicalRaviartThomasBasis(k_index, El, n, j);

		}
	}

	if (VelocityDotNormal >= 0.0) {

		real const CKmean = ξ_prev(k_index, 0) * QuadraturePoints_PolynomialBasis(n, El, 0);

		for (unsigned m = 0; m < 3; m++)
			Concentration += ξ_prev(k_index, m) * QuadraturePoints_PolynomialBasis(n, El, m);

		if (K->neighbors[El]) {

			unsigned const kn_index = K->neighbors[El]->index;
			unsigned const e_index_Kn_loc = K->neighbors[El]->get_edge_index(E);

			real const CKNmean = ξ_prev(kn_index, 0) * QuadraturePoints_PolynomialBasis(n, e_index_Kn_loc, 0);


			real const Cmin = std::min(CKmean, CKNmean);
			real const Cmax = std::max(CKmean, CKNmean);


			if (Concentration < Cmin)
				return Cmin;
			else if (Concentration > Cmax)
				return Cmax;

		}

		return Concentration;

	}
	else {


		if (E->marker == E_MARKER::DIRICHLET) {

			real const X = QuadraturePoints_Edge_x(k_index, El, n);
			real const Y = QuadraturePoints_Edge_y(k_index, El, n);

			return DIRICHLET_GAMMA_P_concentration(X, Y, time);

		}

		unsigned const kn_index = K->neighbors[El]->index;
		unsigned const e_index_Kn_loc = K->neighbors[El]->get_edge_index(E);


		real const CKmean = ξ_prev(k_index, 0) * QuadraturePoints_PolynomialBasis(n, El, 0);
		real const CKNmean = ξ_prev(kn_index, 0) * QuadraturePoints_PolynomialBasis(n, e_index_Kn_loc, 0);

		real const Cmin = std::min(CKmean, CKNmean);
		real const Cmax = std::max(CKmean, CKNmean);


		Concentration = 0.0;
		for (unsigned m = 0; m < 3; m++)
			Concentration += ξ_prev(kn_index, m) * QuadraturePoints_PolynomialBasis(n, e_index_Kn_loc, m);

		if (Concentration < Cmin)
			return Cmin;
		else if (Concentration > Cmax)
			return Cmax;


		return Concentration;


	}

	/*
	if (VelocityDotNormal >= 0.0) {

		for (unsigned m = 0; m < 3; m++)
			Concentration += ξ_prev(k_index, m) * QuadraturePoints_PolynomialBasis(n, El, m);

		//Concentration = ξ_prev(k_index, 0) * QuadraturePoints_PolynomialBasis(n, El, 0);

	}
	else {


		if (E->marker == E_MARKER::DIRICHLET) {

			real const X = QuadraturePoints_Edge_x(k_index, El, n);
			real const Y = QuadraturePoints_Edge_y(k_index, El, n);

			return DIRICHLET_GAMMA_P_concentration(X, Y, time);

		}


		unsigned const kn_index = K->neighbors[El]->index;
		unsigned const e_index_Kn_loc = K->neighbors[El]->get_edge_index(E);


		for (unsigned m = 0; m < 3; m++)
			Concentration += ξ_prev(kn_index, m) * QuadraturePoints_PolynomialBasis(n, e_index_Kn_loc, m);

		//Concentration = ξ_prev(kn_index, 0) * QuadraturePoints_PolynomialBasis(n, e_index_Kn_loc, 0);

	}
	*/

	
	//return Concentration;

};



template<unsigned QuadraturePrecision>
void solver<QuadraturePrecision>::assembleR() {


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

				R.coeffRef(e_index, k_index) = 0.0;

			}

			continue;

		}


		for (unsigned neighborElement = 0; neighborElement < 2; neighborElement++) {

			t_pointer const K = E->neighbors[neighborElement];

			if (!K)
				continue;

			unsigned const k_index = K->index;

			real Alpha = 0.0;

			for (unsigned i = 0; i < 3; i++)
				Alpha += α(k_index, LI(K, E), i);

			real const val = Alpha / viscosities[k_index];

			R.coeffRef(e_index, k_index) = abs(val) < INTEGRAL_PRECISION ? 0.0 : val;

		}
	}

};
template<unsigned QuadraturePrecision>
void solver<QuadraturePrecision>::assembleM() {


	std::vector<Eigen::Triplet<real>> triplet;

	for (unsigned e = 0; e < ne; e++) {


		e_pointer const E = mesh->get_edge(e);
		unsigned const e_index = E->index;
		E_MARKER const e_marker = E->marker;


		if (e_marker == E_MARKER::DIRICHLET) {

			Eigen::Triplet<real> const T(e_index, e_index, -1.0);
			tripletM.push_back(T);

			continue;

		}

		for (unsigned neighborElement = 0; neighborElement < 2; neighborElement++) {


			t_pointer const K = E->neighbors[neighborElement];

			if (!K)
				continue;


			unsigned const k_index = K->index;


			for (unsigned El = 0; El < 3; El++) {


				e_pointer const E_local = K->edges[El];

				unsigned const e_local_index_local = K->get_edge_index(E_local);	// Local index of local edge
				unsigned const e_local_index_global = E_local->index;				// Global index of local edge

				real Alpha = 0.0;

				for (unsigned i = 0; i < 3; i++)
					Alpha += α(k_index, LI(K, E), i);

				Eigen::Triplet<real> const T(e_index, e_local_index_global, Alpha / viscosities[k_index]);
				tripletM.push_back(T);

			}
		}
	}

	SparseMatrix tracePressureSystem_LU(ne, ne);

	tracePressureSystem_LU.setFromTriplets(triplet.begin(), triplet.end());

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

};


template<unsigned QuadraturePrecision>
void solver<QuadraturePrecision>::assembleInverseD() {


	assemble_σ();

	real const TimeCoeff = θ * dt;


	for (unsigned k = 0; k < nk; k++) {


		t_pointer const K = mesh->get_triangle(k);
		unsigned const k_index = K->index;
		real const Area = K->area();

		iD[k_index] = 1.0 / (1.0 + TimeCoeff * σ[k_index]);

	}

};
template<unsigned QuadraturePrecision>
void solver<QuadraturePrecision>::assembleH() {


	assemble_λ();

	real const TimeCoeff = -θ * dt;


	for (unsigned k = 0; k < nk; k++) {


		t_pointer const K = mesh->get_triangle(k);
		unsigned const k_index = K->index;

		e_pointer const E0 = K->edges[0];
		e_pointer const E1 = K->edges[1];
		e_pointer const E2 = K->edges[2];

		unsigned const e_index0 = E0->index;
		unsigned const e_index1 = E1->index;
		unsigned const e_index2 = E2->index;


		H.coeffRef(k_index, e_index0) = TimeCoeff * λ(k_index, 0);
		H.coeffRef(k_index, e_index1) = TimeCoeff * λ(k_index, 1);
		H.coeffRef(k_index, e_index2) = TimeCoeff * λ(k_index, 2);

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

	Eigen::Matrix3d Integral;
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

		for (unsigned n = 0; n < NumberOfQuadraturePointsTriangle; n++) {


			real const s = (real)QuadratureOnTriangle.points_x[n];
			real const t = (real)QuadratureOnTriangle.points_y[n];
			real const w = (real) 0.5 * QuadratureOnTriangle.weights[n];

			// Corresponding coordinates on the element K
			real const x = x0 + JF(0, 0) * s + JF(0, 1) * t;
			real const y = y0 + JF(1, 0) * s + JF(1, 1) * t;

			// Get the inverse of the permeability tensor
			Eigen::Matrix2d K(2, 2);

			K << 1.0, 0.0,
				 0.0, 1.0;

			//K.coeffRef(0, 0) = 1.0;
			//K.coeffRef(0, 1) = 0.0;
			//K.coeffRef(1, 0) = 0.0;
			//K.coeffRef(1, 1) = 1.0;

			Eigen::Matrix2d const iK = K.inverse();

			evaluate_raviartthomas_basis(s, t, BasisRaviartThomas);


			for (unsigned i = 0; i < 3; i++) {


				Eigen::Vector2d const JFWi = JF * BasisRaviartThomas.col(i);

				for (unsigned j = 0; j < 3; j++) {


					Eigen::Vector2d const JFWj = JF * BasisRaviartThomas.col(j);

					Integral.coeffRef(i, j) += w * JFWi.dot(iK * JFWj);

				}
			}
		}

		Integral = Integral / detJF;

		for (unsigned i = 0; i < 3; i++)
			for (unsigned j = 0; j < 3; j++)
				Integral.coeffRef(i, j) = abs(Integral(i, j)) < INTEGRAL_PRECISION ? 0.0 : Integral(i, j);

		Integral = Integral.inverse();

		for (unsigned i = 0; i < 3; i++)
			for (unsigned j = 0; j < 3; j++)
				α.setCoeff(k_index, i, j) = abs(Integral(i, j)) < INTEGRAL_PRECISION ? 0.0 : Integral(i, j);

	}

};
template<unsigned QuadraturePrecision>
void solver<QuadraturePrecision>::assemble_σ() {



	for (unsigned k = 0; k < nk; k++) {

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

		σ.setCoeff(k_index, m, l) = ElementCoeff * Val;

	}

};
template<unsigned QuadraturePrecision>
void solver<QuadraturePrecision>::assemble_λ() {


	for (unsigned k = 0; k < nk; k++) {

		unsigned const k_index = ElementIndeces[k];
		real const Area = K->area();

		real const ElementCoeff = betas[k_index] / (porosities[k_index] * viscosities[k_index] * Area);


		for (unsigned j = 0; j < 3; j++) {

			real Val = 0.0;

			for (unsigned El = 0; El < 3; El++) {


				real const Alpha = α(k_index, El, j);

				real const tC = upwindConcentration(K, El);

				Val += tC * Alpha;

			}

			λ.setCoeff(k_index, j) = ElementCoeff * Val;

		}
	}

};



template<unsigned QuadraturePrecision>
void solver<QuadraturePrecision>::getSolution() {


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

	//hard_copy(betas_prev, betas, nk);
	for (unsigned k = 0; k < nk; k++)
		betas_prev[k] = betas[k];



	unsigned counter = 0;

	while (counter < MAX_IT) {


		assemble_δ();
		assemble_γ();

		computePressureEquation();
		computeVelocities();

		updateConcentrations_explicit();

		//concentrationCorrection();

		computeBetas();


		if (stopCriterion())
			break;



		// Set new iteration level	l := l+1
		π_prev = π;
		ξ_prev = ξ;

		//std::ofstream txtFile;
		//txtFile.open("C:\\Users\\pgali\\Desktop\\output_pressure.txt");
		//exportSolution(txtFile);
		//txtFile.close();


		//hard_copy(betas_prev, betas, nk);
		for (unsigned k = 0; k < nk; k++)
			betas_prev[k] = betas[k];


		counter++;

	}

	std::cout << counter << std::endl;

	assemble_λ();
	assemble_σ();

	// This is needed for the next time level
	for (unsigned k = 0; k < nk; k++) {

		//unsigned const k_index = mesh->get_triangle(k)->index;
		unsigned const k_index = ElementIndeces[k];

		for (unsigned m = 0; m < 3; m++) {

			real val1 = 0.0;
			real val2 = 0.0;

			for (unsigned j = 0; j < 3; j++)
				val1 += σ(k_index, m, j) * π(k_index, j);

			for (unsigned El = 0; El < 3; El++)
				for (unsigned s = 0; s < 2; s++)
					val2 += λ(k_index, s, m, El) * tπ(k_index, El, s);

			rkFp.setCoeff(k_index, m) = val1 + val2;

		}
	}


	if (nt % 1000 == 0)
		std::cout << nt << " - Iterations : " << counter << std::endl;

};
template<unsigned QuadraturePrecision>
void solver<QuadraturePrecision>::exportSolution(std::ofstream & txtFile) {


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

		Eigen::Vector2d const P0(x0 - x0, y0 - y0);
		Eigen::Vector2d const P1(x1 - x0, y1 - y0);
		Eigen::Vector2d const P2(x2 - x0, y2 - y0);

		real const S0 = (iJF*P0)(0);
		real const T0 = (iJF*P0)(1);

		real const S1 = (iJF*P1)(0);
		real const T1 = (iJF*P1)(1);

		real const S2 = (iJF*P2)(0);
		real const T2 = (iJF*P2)(1);

		double const S[3] = { S0, S1, S2 };
		double const T[3] = { T0, T1, T2 };

		//double const S[4] = { S0, S1, S2, S0 };
		//double const T[4] = { T0, T1, T2, T0 };


		for (unsigned i = 0; i < 3; i++) {

			//// Pressures
			real const value = π(index, 0) * phi1(S[i], T[i]) + π(index, 1) * phi2(S[i], T[i]) + π(index, 2) * phi3(S[i], T[i]);
			//// Concentrations
			//real const value = ξ(index, 0) * phi1(S[i], T[i]) + ξ(index, 1) * phi2(S[i], T[i]) + ξ(index, 2) * phi3(S[i], T[i]);

			txtFile << std::setprecision(20) << x[i] << "\t" << y[i] << "\t" << value << std::endl;
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

	Eigen::VectorXd BasisPolynomial(3);
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

			evaluate_polynomial_basis(s, t, BasisPolynomial);

			real PressureK = 0.0;

			for (unsigned j = 0; j < 3; j++)
				PressureK += π(k_index, j) * BasisPolynomial(j);

			Integral += w * sqr(PressureK - barenblatt(X, Y, time));

		}


		errorL2 += detJF * Integral;

		/*
		for (unsigned j = 0; j < error_quadRule; j++) {


			y = quad.points[j];
			w = quad.weigths[j];

			pressure = π(K, 0, 0)*phi1(y, a, b) + π(K, 0, 1)*phi2(y, a, b);
			analytic = barenblatt(y, t);

			val = abs(pressure - analytic);

			s1 += w * val;
			s2 += w * sqr(val);

			if (val > errorMAX)
				errorMAX = val;

		}

		errorL1 += s1;
		errorL2 += s2;
		*/

	}

	//txtFile << "#L1 L2 MAX" << std::endl;
	//txtFile << std::setprecision(20) << errorL1 << std::endl;
	//txtFile << std::setprecision(20) << sqrt(errorL2) << std::endl;
	//txtFile << std::setprecision(20) << errorMAX << std::endl;

	//std::cout << "Error L1 : " << errorL1 << std::endl;
	std::cout << "Error L2 : " << sqrt(errorL2) << std::endl;
	//std::cout << "Error Max : " << errorMAX << std::endl;


};