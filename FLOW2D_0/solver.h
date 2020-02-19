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

	// Unknows vectors of quantities
	CoeffMatrix1D<3>	π;
	CoeffMatrix1D<3>	ξ;
	CoeffMatrix2D<3, 2> tπ;
	CoeffMatrix1D<8>	υ;
	CoeffMatrix1D<3>	sources;

	// Qunatities on time levels. Used in theta-scheme current quantites on time level 'n', last iteration in l 'prev'
	CoeffMatrix1D<3> π_n;
	CoeffMatrix1D<3> π_prev;

	CoeffMatrix1D<3> ξ_n;
	CoeffMatrix1D<3> ξ_prev;

	CoeffMatrix1D<3> rkFc;
	CoeffMatrix1D<3> rkFc_n;

	CoeffMatrix1D<3> rkFp;
	CoeffMatrix1D<3> rkFp_n;


	// Integral (geometric) coeffiecients
	CoeffMatrix2D<8, 8>				α;
	CoeffMatrix3D<3, 3, 8>			δ;
	CoeffMatrix2D<3, 8>				γ;
	SingleCoeffMatrix2D<8, 3>		β;
	SingleCoeffMatrix3D<8, 3, 2>	χ;
	SingleCoeffMatrix2D<3, 3>		η;
	SingleCoeffMatrix3D<3, 8, 3>	τ;


	CoeffMatrix2D<3, 3>		σ;
	CoeffMatrix3D<2, 3, 3>	λ;
	CoeffMatrix1D<3>		φ;
	CoeffMatrix1D<3>		ψ;
	CoeffMatrix1D<8>		Ω;



	static unsigned const NumberOfQuadraturePointsEdge = get_number_of_quadrature_points_edge<QuadraturePrecision>();
	static unsigned const NumberOfQuadraturePointsTriangle = get_number_of_quadrature_points_triangle<QuadraturePrecision>();


	// Precomputed values of basis function in quadrature points
	CoeffMatrix3D<3, 8, 2>								QuadraturePoints_RaviartThomasBasis;
	CoeffMatrix2D<3, 3>									QuadraturePoints_PolynomialBasis;
	CoeffMatrix3D<8, 3, 3>								QuadraturePoints_RaviartThomasBasisDotNormalTimesPolynomialBasis;
	CoeffMatrix3D<3, NumberOfQuadraturePointsEdge, 8>	QuadraturePoints_PhysicalNormalDotPhysicalRaviartThomasBasis;
	CoeffMatrix3D<8, 3, 2>								AlphaTimesChi;
	CoeffMatrix2D<8, 3>									AlphaTimesBeta;

	// Quadrature points on the edges of the physical triangle
	CoeffMatrix2D<3, NumberOfQuadraturePointsEdge>	 QuadraturePoints_Edge_x;
	CoeffMatrix2D<3, NumberOfQuadraturePointsEdge>	 QuadraturePoints_Edge_y;


	// Determinant of the Affine mapping matrix JF
	real * AffineMappingMatrixDeterminant = NULL;

	unsigned * ElementIndeces = NULL;

	CoeffMatrix1D<3> edgeOrientation;


	real * viscosities;
	real * porosities;

	real * betas;
	real * betas_prev;

	// Trace pressure system matrix
	SparseMatrix R1;
	SparseMatrix R2;
	SparseMatrix M_j1_s1;
	SparseMatrix M_j1_s2;
	SparseMatrix M_j2_s1;
	SparseMatrix M_j2_s2;
	DenseVector V1;
	DenseVector V2;
	DenseVector Tp1;
	DenseVector Tp2;
	DenseVector π_eigen;

	SparseMatrix R1iD;
	SparseMatrix R2iD;

	// Pressure system matrix
	SparseMatrix iD;
	SparseMatrix H1;
	SparseMatrix H2;
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


	void assemblePressureSystemMatrix1();
	void assemblePressureSystemMatrix2();
	void assemblePressureSystemMatrix3();
	void assemblePressureSystemMatrix4();
	void assemblePressureSystemMatrixAndInverseDAndH_1();
	void assemblePressureSystemMatrixAndInverseDAndH_2();


	CoeffMatrix2D<3, 3> R1_matrix;
	CoeffMatrix2D<3, 3> R2_matrix;


	void assemble_α();
	void assemble_β();
	void assemble_χ();
	void assemble_η();
	void assemble_τ();
	void assemble_δ();
	void assemble_γ();

	void assemble_σ();
	void assemble_λ();
	void assemble_φ();
	void assemble_ψ();



	real velocity_in_normal_direction(t_pointer const K, e_pointer const E, real const x);
	real continuity_of_normal_component(t_pointer const K, e_pointer const E, real const x);



};











template<unsigned QuadraturePrecision>
solver<QuadraturePrecision>::solver(Mesh & m, unsigned nt_0, double dt_0)

	: nk(m.get_number_of_triangles()), ne(m.get_number_of_edges())
	//,
	//  NumberOfQuadraturePointsEdge(get_number_of_quadrature_points_edge()),
	//  NumberOfQuadraturePointsTriangle(get_number_of_quadrature_points_triangle())

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
	/*    - Memory allocation for the unknowns : π - Internal Pressures		     */
	/*										   : tπ - Trace pressures		     */
	/*										   : ξ - Concentrations     	     */
	/*										   : υ - Velocities				     */
	/*										   : π - Internal Pressures		     */
	/*                                                                           */
	/*****************************************************************************/

	π		.setNumberOfElements(nk);
	tπ		.setNumberOfElements(nk);
	ξ		.setNumberOfElements(nk);
	υ		.setNumberOfElements(nk);
	sources	.setNumberOfElements(nk);	// Sink / sources

	π		.setZero();
	tπ		.setZero();
	ξ		.setZero();
	υ		.setZero();
	sources	.setZero();


	/*****************************************************************************/
	/*                                                                           */
	/*    - Memory allocation for the auxillary variables used in			     */
	/*      time discretization                                                  */
	/*    - 'n'    = current time level                                          */
	/*      'prev' = previous iteration in the inner loop                        */
	/*                                                                           */
	/*****************************************************************************/

	π_n		.setNumberOfElements(nk);
	π_prev	.setNumberOfElements(nk);
	ξ_n		.setNumberOfElements(nk);
	ξ_prev	.setNumberOfElements(nk);
	rkFp	.setNumberOfElements(nk);
	rkFp_n	.setNumberOfElements(nk);
	rkFc	.setNumberOfElements(nk);
	rkFc_n	.setNumberOfElements(nk);

	π_n		.setZero();
	π_prev	.setZero();
	ξ_n		.setZero();
	ξ_prev	.setZero();
	rkFp	.setZero();
	rkFp_n	.setZero();
	rkFc	.setZero();
	rkFc_n	.setZero();


	/*****************************************************************************/
	/*                                                                           */
	/*    - Memory allocation for the integral coefficients				         */
	/*                                                                           */
	/*****************************************************************************/

	α.setNumberOfElements(nk);
	γ.setNumberOfElements(nk);
	δ.setNumberOfElements(nk);

	σ.setNumberOfElements(nk);
	λ.setNumberOfElements(nk);
	φ.setNumberOfElements(nk);


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


	//For All sparse matrix call matrix.data().squeeze()
	//Maybe up to now, there is still much unused memory?????


	/*****************************************************************************/
	/*                                                                           */
	/*    - Memory allocation for the pressure system					         */
	/*                                                                           */
	/*****************************************************************************/

	internalPressureSystem	.resize(2 * ne, 2 * ne);
	pressureSystemRhs		.resize(2 * ne);

	//internalPressureSystem.reserve(Eigen::VectorXi::Constant(2 * ne, 10));


	// Inverse matrix of the matrix for Internal Pressures
	iD.resize(3 * nk, 3 * nk);

	H1.resize(3 * nk, ne);
	H2.resize(3 * nk, ne);


	R1_matrix.setNumberOfElements(nk);
	R2_matrix.setNumberOfElements(nk);

	R1iD.resize(ne, 3 * nk);
	R2iD.resize(ne, 3 * nk);


	// Right-hand side of the Pressure equation
	G.resize(3 * nk);

	// Resulting Internal Pressures: 1. Mean Pressure, 2. Linear Pressure in x-direction, 3. Linear Pressure in y-direction
	π_eigen.resize(3 * nk);


	/*****************************************************************************/
	/*                                                                           */
	/*    - Memory allocation for the Trace Pressure system				         */
	/*                                                                           */
	/*****************************************************************************/

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
	edgeOrientation.setNumberOfElements(nk);

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
	QuadraturePoints_PolynomialBasis									.setNumberOfElements(NumberOfQuadraturePointsEdge);
	QuadraturePoints_RaviartThomasBasisDotNormalTimesPolynomialBasis	.setNumberOfElements(NumberOfQuadraturePointsEdge);
	QuadraturePoints_PhysicalNormalDotPhysicalRaviartThomasBasis		.setNumberOfElements(nk);
	AlphaTimesChi														.setNumberOfElements(nk);
	AlphaTimesBeta														.setNumberOfElements(nk);

	QuadraturePoints_Edge_x												.setNumberOfElements(nk);
	QuadraturePoints_Edge_y												.setNumberOfElements(nk);

	AffineMappingMatrixDeterminant = new real[nk];

	ElementIndeces = new unsigned[nk];






	QuadraturePoints_RaviartThomasBasis									.setZero();
	QuadraturePoints_PolynomialBasis									.setZero();
	QuadraturePoints_RaviartThomasBasisDotNormalTimesPolynomialBasis	.setZero();
	QuadraturePoints_PhysicalNormalDotPhysicalRaviartThomasBasis		.setZero();
	AlphaTimesChi														.setZero();
	AlphaTimesBeta														.setZero();

	QuadraturePoints_Edge_x												.setZero();
	QuadraturePoints_Edge_y												.setZero();



	Eigen::MatrixXd JF(2, 2);
	Eigen::MatrixXd ReferenceNormals(2, 3);
	Eigen::VectorXd Parametrization(2);
	Eigen::MatrixXd BasisRaviartThomas(2, 8);
	Eigen::VectorXd BasisPolynomial(3);

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




		/*****************************************************************************/
		/*                                                                           */
		/*    - Precomputed values of : Alpha * Beta							     */
		/*                                                                           */
		/*    - Used in : Sigma computation											 */
		/*                                                                           */
		/*****************************************************************************/


		for (unsigned l = 0; l < 3; l++) {

			for (unsigned j = 0; j < 8; j++) {

				real AB = 0.0;

				for (unsigned i = 0; i < 8; i++)
					AB += α(k_index, j, i) * β(i, l);

				AlphaTimesBeta.setCoeff(k_index, j, l) = AB;

			}
		}


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


				for (unsigned j = 0; j < 8; j++)
					QuadraturePoints_PhysicalNormalDotPhysicalRaviartThomasBasis.setCoeff(k_index, El, n, j) = PhysicalNormal.dot(JF * BasisRaviartThomas.col(j)) / detJF;

			}




			/*****************************************************************************/
			/*                                                                           */
			/*    - Precomputed values of : Alpha * Chi								     */
			/*                                                                           */
			/*    - Used in : Lambda computation				                           */
			/*                                                                           */
			/*****************************************************************************/

			for (unsigned s = 0; s < 2; s++) {

				for (unsigned j = 0; j < 8; j++) {

					real ACHI = 0.0;

					for (unsigned i = 0; i < 8; i++)
						ACHI += α(k_index, j, i) * χ(i, El, s);

					AlphaTimesChi.setCoeff(k_index, j, El, s) = ACHI;

				}
			}

			/*****************************************************************************/
			/*                                                                           */
			/*    - Precomputed values of : Matrices R1, R2 in the block form		     */
			/*                                                                           */
			/*    - Used in : Trace / Internal pressure system                           */
			/*                                                                           */
			/*****************************************************************************/

			if (e_marker == E_MARKER::DIRICHLET) {

				R1_matrix.setCoeff(k_index, El, 0) = 0.0;
				R1_matrix.setCoeff(k_index, El, 1) = 0.0;
				R1_matrix.setCoeff(k_index, El, 2) = 0.0;

				R2_matrix.setCoeff(k_index, El, 0) = 0.0;
				R2_matrix.setCoeff(k_index, El, 1) = 0.0;
				R2_matrix.setCoeff(k_index, El, 2) = 0.0;

				continue;

			}

			for (unsigned m = 0; m < 3; m++) {


				real AB1 = 0.0;
				real AB2 = 0.0;

				for (unsigned i = 0; i < 8; i++) {

					AB1 += α(k_index, dof0, i)*β(i, m);
					AB2 += α(k_index, dof1, i)*β(i, m);

				}

				real const val1 = AB1 / viscosities[k_index];
				real const val2 = AB2 / viscosities[k_index];

				R1_matrix.setCoeff(k_index, El, m) = abs(val1) < INTEGRAL_PRECISION ? 0.0 : val1;
				R2_matrix.setCoeff(k_index, El, m) = abs(val2) < INTEGRAL_PRECISION ? 0.0 : val2;

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

		//Vector<real> const normal = ReferenceNormals.getColumn(El);
		Eigen::VectorXd const ReferenceNormal = ReferenceNormals.col(El);


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
			evaluate_polynomial_basis(s, t, BasisPolynomial);



			/*****************************************************************************/
			/*                                                                           */
			/*    - Precomputed values of : RT basis on the reference triangle		     */
			/*							  : Polynomial basis on the reference triangle	 */
			/*							  : RT basis * Normal * Polynomial basis		 */
			/*                                                                           */
			/*****************************************************************************/

			for (unsigned j = 0; j < 8; j++) {


				//Vector<real> const Wj = BasisRaviartThomas.getColumn(j);
				Eigen::VectorXd const Wj = BasisRaviartThomas.col(j);

				//real const dotProduct = dot(Wj, normal);
				real const DotProduct = Wj.dot(ReferenceNormal);


				QuadraturePoints_RaviartThomasBasis.setCoeff(n, El, j, 0) = Wj(0);
				QuadraturePoints_RaviartThomasBasis.setCoeff(n, El, j, 1) = Wj(1);

				for (unsigned m = 0; m < 3; m++) {


					real const Phim		= BasisPolynomial(m);
					real const Value	= w * DotProduct * Phim;

					QuadraturePoints_PolynomialBasis									.setCoeff(n, El, m)		= Phim;
					QuadraturePoints_RaviartThomasBasisDotNormalTimesPolynomialBasis	.setCoeff(n, j, El, m)	= abs(Value) < INTEGRAL_PRECISION ? 0.0 : Value;


				}
			}
		}
	}


};
template<unsigned QuadraturePrecision>
solver<QuadraturePrecision>::~solver() {

	delete[] AffineMappingMatrixDeterminant;

	delete[] ElementIndeces;

	//delete[] QuadraturePoints_Edge0;
	//delete[] QuadraturePoints_Edge12;
	//delete[] QuadratureWeights_Edge0;
	//delete[] QuadratureWeights_Edge12;

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


		/*****************************************************************************/
		/*                                                                           */
		/*    - Interpolant of the Barenblatt solution for the initial condition     */
		/*                                                                           */
		/*****************************************************************************/

		real const B0 = barenblatt(x0, y0, time);
		real const B1 = barenblatt(x1, y1, time);
		real const B2 = barenblatt(x2, y2, time);

		Eigen::Vector3d const B(B0, B1, B2);
		Eigen::Matrix3d M;
		M << 1.0, -1.0, -1.0,
			 1.0, +1.0, -1.0,
			 1.0, -1.0, +1.0;
		Eigen::Vector3d const Solution = M.colPivHouseholderQr().solve(B);

		π.setCoeff(k_index, 0) = Solution(0);
		π.setCoeff(k_index, 1) = Solution(1);
		π.setCoeff(k_index, 2) = Solution(2);

		ξ.setCoeff(k_index, 0) = Solution(0);
		ξ.setCoeff(k_index, 1) = Solution(1);
		ξ.setCoeff(k_index, 2) = Solution(2);

		sources.setCoeff(k_index, 0) = F1(K, time);
		sources.setCoeff(k_index, 1) = F2(K, time);
		sources.setCoeff(k_index, 2) = F3(K, time);


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

		rkFc.setCoeff(k_index, 0) = 0.0;
		rkFc.setCoeff(k_index, 1) = 0.0;
		rkFc.setCoeff(k_index, 2) = 0.0;

		rkFc_n.setCoeff(k_index, 0) = rkFc(k_index, 0);
		rkFc_n.setCoeff(k_index, 1) = rkFc(k_index, 1);
		rkFc_n.setCoeff(k_index, 2) = rkFc(k_index, 2);


		rkFp.setCoeff(k_index, 0) = 0.0;
		rkFp.setCoeff(k_index, 1) = 0.0;
		rkFp.setCoeff(k_index, 2) = 0.0;

		rkFp_n.setCoeff(k_index, 0) = rkFp(k_index, 0);
		rkFp_n.setCoeff(k_index, 1) = rkFp(k_index, 1);
		rkFp_n.setCoeff(k_index, 2) = rkFp(k_index, 2);


		// Set edge orientations of edges on each element K with respect to GLOBAL orientation of respective edge
		for (unsigned e = 0; e < 3; e++) {

			e_pointer const E = K->edges[e];

			v_pointer const a = E->a;
			v_pointer const b = E->b;

			v_pointer const v = K->get_vertex_cw(a);
			v_pointer const p = K->get_vertex(e);

			bool const isBoundaryEdge = E->marker == E_MARKER::NEUMANN || E->marker == E_MARKER::DIRICHLET;


			if (v != p && !isBoundaryEdge)
				edgeOrientation.setCoeff(k_index, e) = -1.0;
			else
				edgeOrientation.setCoeff(k_index, e) = +1.0;

		}

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


	// Assembly eigen-vector of internal pressures
	for (unsigned k = 0; k < nk; k++) {

		//unsigned const k_index = mesh->get_triangle(k)->index;
		unsigned const k_index = ElementIndeces[k];

		for (unsigned m = 0; m < 3; m++)
			π_eigen[3 * k_index + m] = π(k_index, m);

	}

	assembleV();


	traceSystemRhs.head(ne) = R1 * π_eigen - V1;
	traceSystemRhs.tail(ne) = R2 * π_eigen - V2;

	DenseVector const solution = sparseLUsolver_TracePressureSystem.solve(traceSystemRhs);

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

			tπ.setCoeff(k_index, e_index_local, 0) = tpval1;
			tπ.setCoeff(k_index, e_index_local, 1) = tpval2;

		}
	}

};
template<unsigned QuadraturePrecision>
void solver<QuadraturePrecision>::computePressureEquation() {



	assembleG();
	assembleV();

	//assemblePressureSystemMatrixAndInverseDAndH_1();
	assemblePressureSystemMatrixAndInverseDAndH_2();

	pressureSystemRhs.head(ne) = R1iD * G - V1;
	pressureSystemRhs.tail(ne) = R2iD * G - V2;


	Eigen::SparseLU<SparseMatrix> const solver(internalPressureSystem);

	DenseVector const solution = solver.solve(pressureSystemRhs);

	Tp1 = solution.head(ne);
	Tp2 = solution.tail(ne);



	// Internal Pressures
	π_eigen = iD * (G - H1 * Tp1 - H2 * Tp2);

	// Copy Internal Pressure eigen-solution to each element
	for (unsigned k = 0; k < nk; k++) {

		//unsigned const k_index = mesh->get_triangle(k)->index;
		unsigned const k_index = ElementIndeces[k];

		for (unsigned m = 0; m < 3; m++)
			π.setCoeff(k_index, m) = π_eigen[3 * k_index + m];

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

			tπ.setCoeff(k_index, e_index_local, 0) = tpval1;
			tπ.setCoeff(k_index, e_index_local, 1) = tpval2;

		}
	}

	/*
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
*/

};



template<unsigned QuadraturePrecision>
void solver<QuadraturePrecision>::computeVelocities() {


	real const time = (nt + 1)* dt;

	for (unsigned k = 0; k < nk; k++) {


		t_pointer const K = mesh->get_triangle(k);
		unsigned const k_index = K->index;
		//unsigned const k_index = ElementIndeces[k];


		// Loop over element's edges and their degrees of freedom
		for (unsigned El = 0; El < 3; El++) {


			e_pointer const E = K->edges[El];

			if (E->marker == E_MARKER::NEUMANN) {

				υ.setCoeff(k_index, LI(K, E, 0)) = NEUMANN_GAMMA_Q_velocity(E, time);
				υ.setCoeff(k_index, LI(K, E, 1)) = 0.0;

				continue;

			}

			// Loop over the two degrees of freedom on each edge
			for (unsigned dof = 0; dof < 2; dof++) {


				unsigned const j = LI(K, E, dof);

				real val1 = 0.0;
				real val2 = 0.0;

				for (unsigned m = 0; m < 8; m++) {


					real const alpha = α(k_index, m, j);

					real PB = 0.0;
					real ECHITP = 0.0;

					for (unsigned i = 0; i < 3; i++)
						PB += β(m, i) * π(k_index, i);

					val1 += alpha * PB;


					for (unsigned l = 0; l < 3; l++) {

						real CHITP = 0.0;

						for (unsigned s = 0; s < 2; s++)
							CHITP += χ(m, l, s) * tπ(k_index, l, s);

						ECHITP += CHITP;

					}

					val2 += alpha * ECHITP;

				}

				υ.setCoeff(k_index, j) = (val1 - val2) / viscosities[k_index];

			}
		}


		// compute Bubble velocities inside element
		for (unsigned dof = 6; dof < 8; dof++) {


			real val1 = 0.0;
			real val2 = 0.0;

			for (unsigned m = 0; m < 8; m++) {


				real const alpha = α(k_index, dof, m);

				real PB = 0.0;
				real ECHITP = 0.0;

				for (unsigned i = 0; i < 3; i++)
					PB += β(m, i) * π(k_index, i);

				val1 += alpha * PB;


				for (unsigned l = 0; l < 3; l++) {

					real CHITP = 0.0;

					for (unsigned s = 0; s < 2; s++)
						CHITP += χ(m, l, s) * tπ(k_index, l, s);

					ECHITP += CHITP;

				}

				val2 += alpha * ECHITP;

			}

			υ.setCoeff(k_index, dof) = (val1 - val2) / viscosities[k_index];

		}
	}

};
template<unsigned QuadraturePrecision>
void solver<QuadraturePrecision>::updateConcentrations_explicit() {


	for (unsigned k = 0; k < nk; k++) {


		unsigned const k_index = ElementIndeces[k];

		real val0 = 0.0;
		real val1 = 0.0;
		real val2 = 0.0;

		for (unsigned j = 0; j < 8; j++) {


			real etaGamma0 = 0.0;
			real etaGamma1 = 0.0;
			real etaGamma2 = 0.0;

			for (unsigned l = 0; l < 3; l++) {

				etaGamma0 += η(0, l) * γ(k_index, l, j);
				etaGamma1 += η(1, l) * γ(k_index, l, j);
				etaGamma2 += η(2, l) * γ(k_index, l, j);

			}

			real const detJF = AffineMappingMatrixDeterminant[k_index];

			// In the computation of Eta, there is coefficient detJF. When inverting, the coefficient is inverted
			val0 += υ(k_index, j) * etaGamma0 / detJF;
			val1 += υ(k_index, j) * etaGamma1 / detJF;
			val2 += υ(k_index, j) * etaGamma2 / detJF;

		}

		real const Porosity = porosities[k_index];

		rkFc.setCoeff(k_index, 0) = -val0 / Porosity;
		rkFc.setCoeff(k_index, 1) = -val1 / Porosity;
		rkFc.setCoeff(k_index, 2) = -val2 / Porosity;

		// General θ scheme
		ξ.setCoeff(k_index, 0) = ξ_n(k_index, 0) + dt * (θ * rkFc(k_index, 0) + (1.0 - θ) * rkFc_n(k_index, 0));
		ξ.setCoeff(k_index, 1) = ξ_n(k_index, 1) + dt * (θ * rkFc(k_index, 1) + (1.0 - θ) * rkFc_n(k_index, 1));
		ξ.setCoeff(k_index, 2) = ξ_n(k_index, 2) + dt * (θ * rkFc(k_index, 2) + (1.0 - θ) * rkFc_n(k_index, 2));

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

					AB1 += α(k_index, dof0, i)*β(i, m);
					AB2 += α(k_index, dof1, i)*β(i, m);

				}

				real const val1 = AB1 / viscosities[k_index];
				real const val2 = AB2 / viscosities[k_index];

				R1.coeffRef(e_index, 3 * k_index + m) = abs(val1) < INTEGRAL_PRECISION ? 0.0 : val1;
				R2.coeffRef(e_index, 3 * k_index + m) = abs(val2) < INTEGRAL_PRECISION ? 0.0 : val2;

			}
		}
	}

};
template<unsigned QuadraturePrecision>
void solver<QuadraturePrecision>::assembleM() {


	//M_j1_s1.setZero();
	//M_j1_s2.setZero();

	//M_j2_s1.setZero();
	//M_j2_s2.setZero();

	std::vector<Eigen::Triplet<real>> triplet;


	for (unsigned e = 0; e < ne; e++) {


		e_pointer const E = mesh->get_edge(e);
		unsigned const e_index = E->index;
		E_MARKER const e_marker = E->marker;



		if (e_marker == E_MARKER::DIRICHLET) {

			M_j1_s1.coeffRef(e_index, e_index) = -1.0;
			M_j1_s2.coeffRef(e_index, e_index) = 0.0;

			M_j2_s1.coeffRef(e_index, e_index) = 0.0;
			M_j2_s2.coeffRef(e_index, e_index) = -1.0;

			Eigen::Triplet<real> const T1(e_index,		e_index,		-1.0);
			Eigen::Triplet<real> const T2(e_index,		e_index + ne,	0.0);
			Eigen::Triplet<real> const T3(e_index + ne, e_index,		0.0);
			Eigen::Triplet<real> const T4(e_index + ne, e_index + ne,	-1.0);

			triplet.push_back(T1);
			triplet.push_back(T2);
			triplet.push_back(T3);
			triplet.push_back(T4);

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

					ACHI_j1_s1 += α(k_index, dof0, i) * χ(i, e_local_index_local, 0);
					ACHI_j1_s2 += α(k_index, dof0, i) * χ(i, e_local_index_local, 1);

					ACHI_j2_s1 += α(k_index, dof1, i) * χ(i, e_local_index_local, 0);
					ACHI_j2_s2 += α(k_index, dof1, i) * χ(i, e_local_index_local, 1);

				}

				real const val1 = ACHI_j1_s1 / viscosities[k_index];
				real const val2 = ACHI_j1_s2 / viscosities[k_index];

				real const val3 = ACHI_j2_s1 / viscosities[k_index];
				real const val4 = ACHI_j2_s2 / viscosities[k_index];

				
				/*if (e_local_index_global == e_index) {
					M_j1_s1.coeffRef(e_index, e_local_index_global) += val1;
					M_j1_s2.coeffRef(e_index, e_local_index_global) += val2;
					M_j2_s1.coeffRef(e_index, e_local_index_global) += val3;
					M_j2_s2.coeffRef(e_index, e_local_index_global) += val4;
					continue;
				}
				M_j1_s1.coeffRef(e_index, e_local_index_global) = val1;
				M_j1_s2.coeffRef(e_index, e_local_index_global) = val2;
				M_j2_s1.coeffRef(e_index, e_local_index_global) = val3;
				M_j2_s2.coeffRef(e_index, e_local_index_global) = val4;*/
				

				M_j1_s1.coeffRef(e_index, e_local_index_global) += val1;
				M_j1_s2.coeffRef(e_index, e_local_index_global) += val2;

				M_j2_s1.coeffRef(e_index, e_local_index_global) += val3;
				M_j2_s2.coeffRef(e_index, e_local_index_global) += val4;

				Eigen::Triplet<real> const T1(e_index,		e_local_index_global,		val1);
				Eigen::Triplet<real> const T2(e_index,		e_local_index_global + ne,	val2);
				Eigen::Triplet<real> const T3(e_index + ne, e_local_index_global,		val3);
				Eigen::Triplet<real> const T4(e_index + ne, e_local_index_global + ne,	val4);

				triplet.push_back(T1);
				triplet.push_back(T2);
				triplet.push_back(T3);
				triplet.push_back(T4);

			}
		}
	}

	SparseMatrix tracePressureSystem_LU(2 * ne, 2 * ne);

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


		V1[e_index] = 0.0;
		V2[e_index] = 0.0;

		if (marker == E_MARKER::NEUMANN)
			V1[e_index] = NEUMANN_GAMMA_Q_velocity(E, time);
		//else if (marker == E_MARKER::DIRICHLET)
		//	V1[e_index] = DIRICHLET_GAMMA_P_pressure(E, time);
		else if (marker == E_MARKER::DIRICHLET) {


			real const x0 = (real) E->a->x;
			real const y0 = (real) E->a->y;

			real const x1 = (real) E->b->x;
			real const y1 = (real) E->b->y;


			/*****************************************************************************/
			/*                                                                           */
			/*    - Interpolant of the Barenblatt solution for the Dirichlet Edge	     */
			/*                                                                           */
			/*****************************************************************************/

			real const B0 = barenblatt(x0, y0, time);
			real const B1 = barenblatt(x1, y1, time);

			Eigen::Vector2d const B(B0, B1);
			Eigen::Matrix2d M;
			M << +1.0, -1.0,
				 +1.0, +1.0;
			Eigen::Vector2d const Solution = M.colPivHouseholderQr().solve(B);

			V1[e_index] = Solution(0);
			V2[e_index] = Solution(1);

		}
	}

};




template<unsigned QuadraturePrecision>
void solver<QuadraturePrecision>::assembleInverseD() {


	assemble_σ();

	real const coeff = θ * dt;
	Eigen::Matrix3d block;


	for (unsigned k = 0; k < nk; k++) {


		t_pointer const K = mesh->get_triangle(k);

		unsigned const k_index = K->index;
		unsigned const start_index = 3 * k_index;


		for (unsigned r = 0; r < 3; r++)
			for (unsigned s = 0; s < 3; s++)
				block(r, s) = δij(r, s) - coeff * σ(k_index, r, s);


		Eigen::Matrix3d const inverse_block = block.inverse();


		for (unsigned i = 0; i < 3; i++)
			for (unsigned j = 0; j < 3; j++)
				iD.coeffRef(start_index + i, start_index + j) = abs(inverse_block(i, j)) < INTEGRAL_PRECISION ? 0.0 : inverse_block(i, j);

	}

};
template<unsigned QuadraturePrecision>
void solver<QuadraturePrecision>::assembleH() {


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


		/*
		for (unsigned m = 0; m < 3; m++) {
			H1.coeffRef(3 * k_index + m, e_index0) = block1.coeff(m, 0);
			H1.coeffRef(3 * k_index + m, e_index1) = block1.coeff(m, 1);
			H1.coeffRef(3 * k_index + m, e_index2) = block1.coeff(m, 2);
			H2.coeffRef(3 * k_index + m, e_index0) = block2.coeff(m, 0);
			H2.coeffRef(3 * k_index + m, e_index1) = block2.coeff(m, 1);
			H2.coeffRef(3 * k_index + m, e_index2) = block2.coeff(m, 2);
		}
		*/

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
template<unsigned QuadraturePrecision>
void solver<QuadraturePrecision>::assembleG() {


	real const coeff = dt * (1.0 - θ);

	for (unsigned k = 0; k < nk; k++) {

		unsigned const k_index = ElementIndeces[k];

		for (unsigned m = 0; m < 3; m++)
			G[3 * k_index + m] = π_n(k_index, m) + coeff * rkFp_n(k_index, m);

	}

};





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


		/*****************************************************************************/
		/*                                                                           */
		/*    - Inverse of the matrix D										         */
		/*                                                                           */
		/*****************************************************************************/

		for (unsigned r = 0; r < 3; r++)
			for (unsigned s = 0; s < 3; s++)
				block(r, s) = δij(r, s) - coefficient * σ(k_index, r, s);

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
	//tri.reserve(estimation_of_entries);


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


		/*****************************************************************************/
		/*                                                                           */
		/*    - Inverse of the matrix D										         */
		/*                                                                           */
		/*****************************************************************************/

		for (unsigned r = 0; r < 3; r++)
			for (unsigned s = 0; s < 3; s++)
				block(r, s) = δij(r, s) - coefficient * σ(k_index, r, s);

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

					sum11 += R1_matrix(k_index, ei, m) * iDH1block(m, ej);
					sum12 += R1_matrix(k_index, ei, m) * iDH2block(m, ej);
					sum21 += R2_matrix(k_index, ei, m) * iDH1block(m, ej);
					sum22 += R2_matrix(k_index, ei, m) * iDH2block(m, ej);

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


template<unsigned QuadraturePrecision>
void solver<QuadraturePrecision>::assemblePressureSystemMatrixAndInverseDAndH_1() {


	real const coefficient = θ * dt;

	Eigen::Matrix3d block;
	Eigen::Matrix3d block1;
	Eigen::Matrix3d block2;

	assemble_σ();
	assemble_λ();

	//omp_set_num_threads(2);
	//#pragma omp parallel 
	//{
	//	#pragma omp sections
	//	{
	//		#pragma omp section
	//		{
	//			assemble_σ();
	//		}
	//		#pragma omp section
	//		{
	//			assemble_λ();
	//		}
	//	}
	//}


	// It is sufficient to zero only diagonal elements. Therefore, there is no need for += in the sequel
	for (unsigned e = 0; e < ne; e++) {

		internalPressureSystem.coeffRef(e,		e)			= M_j1_s1.coeff(e, e);
		internalPressureSystem.coeffRef(e,		e + ne)		= M_j1_s2.coeff(e, e);
		internalPressureSystem.coeffRef(e + ne, e)			= M_j2_s1.coeff(e, e);
		internalPressureSystem.coeffRef(e + ne, e + ne)		= M_j2_s2.coeff(e, e);

	}

	for (int k = 0; k < nk; k++) {



		t_pointer const K = mesh->get_triangle(k);

		unsigned const k_index = K->index;
		unsigned const start_index = 3 * k_index;


		/*****************************************************************************/
		/*                                                                           */
		/*    - Inverse of the matrix D										         */
		/*                                                                           */
		/*****************************************************************************/

		for (unsigned r = 0; r < 3; r++)
			for (unsigned s = 0; s < 3; s++)
				block(r, s) = δij(r, s) - coefficient * σ(k_index, r, s);

		Eigen::Matrix3d const inverse_block = block.inverse();

		for (unsigned i = 0; i < 3; i++)
			for (unsigned j = 0; j < 3; j++)
				iD.coeffRef(start_index + i, start_index + j) = inverse_block(i, j);



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



			for (unsigned j = 0; j < 3; j++) {


				/*****************************************************************************/
				/*                                                                           */
				/*    - Assemble Matrices H1, H2									         */
				/*                                                                           */
				/*****************************************************************************/

				H1.coeffRef(start_index + j, e_index_i) = block1.coeff(j, ei);
				H2.coeffRef(start_index + j, e_index_i) = block2.coeff(j, ei);


				/*****************************************************************************/
				/*                                                                           */
				/*    - Assemble Matrices R1 * D.inverse, R2 * D.inverse			         */
				/*                                                                           */
				/*****************************************************************************/

				real sumA = 0.0;
				real sumB = 0.0;

				for (unsigned m = 0; m < 3; m++) {

					sumA += R1_matrix(k_index, ei, m) * inverse_block(m, j);
					sumB += R2_matrix(k_index, ei, m) * inverse_block(m, j);

				}

				R1iD.coeffRef(e_index_i, start_index + j) = sumA;
				R2iD.coeffRef(e_index_i, start_index + j) = sumB;

			}


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
void solver<QuadraturePrecision>::assemblePressureSystemMatrixAndInverseDAndH_2() {


	real const coefficient = θ * dt;

	Eigen::Matrix3d block;
	Eigen::Matrix3d block1;
	Eigen::Matrix3d block2;

	std::vector<Eigen::Triplet<real>> tri;
	std::vector<Eigen::Triplet<real>> triR1;
	std::vector<Eigen::Triplet<real>> triR2;

	assemble_σ();
	assemble_λ();

	//omp_set_num_threads(2);
	//#pragma omp parallel 
	//{
	//	#pragma omp sections
	//	{
	//		#pragma omp section
	//		{
	//			assemble_σ();
	//		}
	//		#pragma omp section
	//		{
	//			assemble_λ();
	//		}
	//	}
	//}

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

	//#pragma omp parallel for private(block,block1,block2) shared(tri,triR1,triR2)
	for (unsigned k = 0; k < nk; k++) {



		t_pointer const K = mesh->get_triangle(k);

		unsigned const k_index = K->index;
		unsigned const start_index = 3 * k_index;


		/*****************************************************************************/
		/*                                                                           */
		/*    - Inverse of the matrix D										         */
		/*                                                                           */
		/*****************************************************************************/

		for (unsigned r = 0; r < 3; r++)
			for (unsigned s = 0; s < 3; s++)
				block(r, s) = δij(r, s) - coefficient * σ(k_index, r, s);

		Eigen::Matrix3d const inverse_block = block.inverse();


		// Assemble iD
		for (unsigned i = 0; i < 3; i++)
			for (unsigned j = 0; j < 3; j++)
				iD.coeffRef(start_index + i, start_index + j) = inverse_block(i, j);

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
		/*    - Assembly of the resulting matrix R1 * D.inverse * H1		         */
		/*                                                                           */
		/*****************************************************************************/

		for (unsigned ei = 0; ei < 3; ei++) {


			e_pointer const Ei = K->edges[ei];
			unsigned const e_index_i = K->edges[ei]->index;


			for (unsigned j = 0; j < 3; j++) {


				/*****************************************************************************/
				/*                                                                           */
				/*    - Assembly of the Matrices H1, H2								         */
				/*                                                                           */
				/*****************************************************************************/

				H1.coeffRef(start_index + j, e_index_i) = block1.coeff(j, ei);
				H2.coeffRef(start_index + j, e_index_i) = block2.coeff(j, ei);


				/*****************************************************************************/
				/*                                                                           */
				/*    - Assemble Matrices R1 * D.inverse, R2 * D.inverse			         */
				/*                                                                           */
				/*****************************************************************************/

				real sumA = 0.0;
				real sumB = 0.0;

				for (unsigned m = 0; m < 3; m++) {

					sumA += R1_matrix(k_index, ei, m) * inverse_block(m, j);
					sumB += R2_matrix(k_index, ei, m) * inverse_block(m, j);

				}

				Eigen::Triplet<real> const T1(e_index_i, start_index + j, sumA);
				Eigen::Triplet<real> const T2(e_index_i, start_index + j, sumB);

				triR1.push_back(T1);
				triR2.push_back(T2);

			}


			// Diagonal elements are already zeroed/or there is Mij already
			if (Ei->marker == E_MARKER::DIRICHLET)
				continue;

//#pragma omp parallel for
			for (int ej = 0; ej < 3; ej++) {


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

	R1iD.setFromTriplets(triR1.begin(), triR1.end());
	R2iD.setFromTriplets(triR2.begin(), triR2.end());

};






template<unsigned QuadraturePrecision>
void solver<QuadraturePrecision>::assemble_α() {


	quadrature_triangle const QuadratureOnTriangle(QuadraturePrecision);

	Eigen::MatrixXd Integral(8, 8);
	Eigen::MatrixXd BasisRaviartThomas(2, 8);
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

			K.coeffRef(0, 0) = 1.0;
			K.coeffRef(0, 1) = 0.0;
			K.coeffRef(1, 0) = 0.0;
			K.coeffRef(1, 1) = 1.0;

			Eigen::Matrix2d const iK = K.inverse();

			evaluate_raviartthomas_basis(s, t, BasisRaviartThomas);


			for (unsigned i = 0; i < 8; i++) {


				Eigen::Vector2d const JFWi = JF * BasisRaviartThomas.col(i);

				for (unsigned j = 0; j < 8; j++) {


					Eigen::Vector2d const JFWj = JF * BasisRaviartThomas.col(j);

					Integral.coeffRef(i, j) += w * JFWi.dot(iK * JFWj);

				}
			}
		}

		Integral = Integral / detJF;

		for (unsigned i = 0; i < 8; i++)
			for (unsigned j = 0; j < 8; j++)
				Integral.coeffRef(i, j) = abs(Integral(i, j)) < INTEGRAL_PRECISION ? 0.0 : Integral(i, j);

		Integral = Integral.inverse();

		for (unsigned i = 0; i < 8; i++)
			for (unsigned j = 0; j < 8; j++)
				α.setCoeff(k_index, i, j) = abs(Integral(i, j)) < INTEGRAL_PRECISION ? 0.0 : Integral(i, j);

	}

};
template<unsigned QuadraturePrecision>
void solver<QuadraturePrecision>::assemble_β() {


	quadrature_triangle const QuadratureOnTriangle(QuadraturePrecision);

	Eigen::MatrixXd Integral(8, 3);
	Eigen::VectorXd BasisPolynomial(3);
	Eigen::VectorXd BasisRaviartThomasDivergence(8);

	Integral.setZero();

	for (unsigned n = 0; n < NumberOfQuadraturePointsTriangle; n++) {


		real const s = (real) QuadratureOnTriangle.points_x[n];
		real const t = (real) QuadratureOnTriangle.points_y[n];
		real const w = (real) 0.5 * QuadratureOnTriangle.weights[n];


		evaluate_raviartthomas_basis_divergence(s, t, BasisRaviartThomasDivergence);
		evaluate_polynomial_basis(s, t, BasisPolynomial);


		for (unsigned i = 0; i < 8; i++) {

			real const dWi = BasisRaviartThomasDivergence(i);

			for (unsigned j = 0; j < 3; j++) {

				real const Phij = BasisPolynomial(j);

				Integral.coeffRef(i, j) += w * dWi * Phij;

			}
		}
	}

	for (unsigned m = 0; m < 8; m++)
		for (unsigned j = 0; j < 3; j++)
			β.setCoeff(m, j) = abs(Integral(m, j)) < INTEGRAL_PRECISION ? 0.0 : Integral(m, j);


};
template<unsigned QuadraturePrecision>
void solver<QuadraturePrecision>::assemble_χ() {


	gauss_quadrature_1D const QuadratureOnEdge(QuadraturePrecision);


	Eigen::MatrixXd ReferenceNormals(2, 3);
	Eigen::VectorXd Parametrization(2);
	Eigen::MatrixXd BasisRaviartThomas(2, 8);
	Eigen::VectorXd BasisEdgePolynomial(2);

	evaluate_edge_normal(ReferenceNormals);


	for (unsigned El = 0; El < 3; El++) {


		real const a = (real) 0.0;
		real const b = (real)El != 0 ? 1.0 : sqrt(2.0);

		real const c = (real)(b - a) / 2.0;
		real const d = (real)(b + a) / 2.0;

		Eigen::VectorXd const ReferenceNormal = ReferenceNormals.col(El);


		for (unsigned n = 0; n < NumberOfQuadraturePointsEdge; n++) {


			real const x = (real)QuadratureOnEdge.points[n] * c + d;
			real const w = (real)QuadratureOnEdge.weights[n] * c;


			evaluate_edge_parametrization(x, El, Parametrization);

			real const s = Parametrization(0);
			real const t = Parametrization(1);
			real const drNorm = 1.0;


			evaluate_raviartthomas_basis(s, t, BasisRaviartThomas);
			evaluate_edge_polynomial_basis(x, El, BasisEdgePolynomial);


			for (unsigned m = 0; m < 8; m++) {


				Eigen::VectorXd const Wm = BasisRaviartThomas.col(m);

				real const dotProduct = Wm.dot(ReferenceNormal);


				for (unsigned s = 0; s < 2; s++) {

					real const varPhis = BasisEdgePolynomial(s);

					χ.setCoeff(m, El, s) = χ(m, El, s) + w * dotProduct * varPhis * drNorm;

				}
			}
		}
	}

	for (unsigned m = 0; m < 8; m++)
		for (unsigned El = 0; El < 3; El++)
			for (unsigned j = 0; j < 2; j++)
				χ.setCoeff(m, El, j) = abs(χ(m, El, j)) < INTEGRAL_PRECISION ? 0.0 : χ(m, El, j);

	//for (unsigned m = 0; m < 8; m++) {
	//	//
	//	Matrix<real> integral(3, 2);
	//	integral.setZero();
	//	//
	//	for (unsigned El = 0; El < 3; El++)
	//		for (unsigned j = 0; j < 2; j++)
	//			integral(El, j) = χ(m, El, j);
	//	//
	//	std::cout << integral << std::endl;
	//	//
	//}


};
template<unsigned QuadraturePrecision>
void solver<QuadraturePrecision>::assemble_η() {


	quadrature_triangle const QuadratureOnTriangle(QuadraturePrecision);

	Eigen::MatrixXd Integral(3, 3);
	Eigen::VectorXd BasisPolynomial(3);


	Integral.setZero();

	for (unsigned n = 0; n < NumberOfQuadraturePointsTriangle; n++) {


		real const s = (real)QuadratureOnTriangle.points_x[n];
		real const t = (real)QuadratureOnTriangle.points_y[n];
		real const w = (real) 0.5 * QuadratureOnTriangle.weights[n];

		evaluate_polynomial_basis(s, t, BasisPolynomial);


		for (unsigned m = 0; m < 3; m++) {


			real const Phim = BasisPolynomial(m);

			for (unsigned j = 0; j < 3; j++) {


				real const Phij = BasisPolynomial(j);

				Integral.coeffRef(m, j) += w * Phim * Phij;

			}
		}
	}

	for (unsigned i = 0; i < 3; i++)
		for (unsigned j = 0; j < 3; j++)
			Integral.coeffRef(i, j) = abs(Integral(i, j)) < INTEGRAL_PRECISION ? 0.0 : Integral(i, j);

	Integral = Integral.inverse();

	for (unsigned i = 0; i < 3; i++)
		for (unsigned j = 0; j < 3; j++)
			η.setCoeff(i, j) = abs(Integral(i, j)) < INTEGRAL_PRECISION ? 0.0 : Integral(i, j);

};
template<unsigned QuadraturePrecision>
void solver<QuadraturePrecision>::assemble_τ() {


	quadrature_triangle const QuadratureOnTriangle(QuadraturePrecision);

	Eigen::MatrixXd BasisRaviartThomas(2, 8);
	Eigen::VectorXd BasisPolynomial(3);
	Eigen::MatrixXd BasisPolynomialGradient(2, 3);


	for (unsigned n = 0; n < NumberOfQuadraturePointsTriangle; n++) {


		real const s = (real)QuadratureOnTriangle.points_x[n];
		real const t = (real)QuadratureOnTriangle.points_y[n];
		real const w = (real) 0.5 * QuadratureOnTriangle.weights[n];


		evaluate_raviartthomas_basis(s, t, BasisRaviartThomas);
		evaluate_polynomial_basis(s, t, BasisPolynomial);
		evaluate_polynomial_basis_gradient(s, t, BasisPolynomialGradient);


		for (unsigned m = 0; m < 3; m++) {


			Eigen::VectorXd const dPhim = BasisPolynomialGradient.col(m);

			for (unsigned j = 0; j < 8; j++) {


				Eigen::VectorXd const Wj = BasisRaviartThomas.col(j);
				real const dotProduct = Wj.dot(dPhim);

				for (unsigned l = 0; l < 3; l++) {

					real const Phil = BasisPolynomial(l);

					τ.setCoeff(m, j, l) = τ(m, j, l) + w * dotProduct * Phil;

				}
			}
		}
	}

	for (unsigned m = 0; m < 3; m++)
		for (unsigned i = 0; i < 8; i++)
			for (unsigned j = 0; j < 3; j++)
				τ.setCoeff(m, i, j) = abs(τ(m, i, j)) < INTEGRAL_PRECISION ? 0.0 : τ(m, i, j);

};
template<unsigned QuadraturePrecision>
void solver<QuadraturePrecision>::assemble_γ() {


	for (unsigned k = 0; k < nk; k++) {


		unsigned const k_index = ElementIndeces[k];

		for (unsigned m = 0; m < 3; m++) {

			for (unsigned j = 0; j < 8; j++) {


				real val1 = 0.0;
				real val2 = 0.0;

				for (unsigned l = 0; l < 3; l++)
					val1 += δ(k_index, m, l, j);

				for (unsigned l = 0; l < 3; l++)
					val2 += τ(m, j, l) * ξ_prev(k_index, l);

				γ.setCoeff(k_index, m, j) = val1 - val2;

			}
		}
	}

};
template<unsigned QuadraturePrecision>
void solver<QuadraturePrecision>::assemble_δ() {



	δ.setZero();

	for (unsigned k = 0; k < nk; k++) {


		t_pointer const K = mesh->get_triangle(k);
		unsigned const k_index = K->index;


		for (unsigned El = 0; El < 3; El++) {

			for (unsigned n = 0; n < NumberOfQuadraturePointsEdge; n++) {


				//real const tC = upwindConcentration_limiter(K, El, n);
				real const tC = upwindConcentration(K, El, n);

				for (unsigned m = 0; m < 3; m++)
					for (unsigned j = 0; j < 8; j++)
						δ.setCoeff(k_index, m, El, j) = δ(k_index, m, El, j) + tC * QuadraturePoints_RaviartThomasBasisDotNormalTimesPolynomialBasis(n, j, El, m);

			}
		}

		/*for (unsigned El = 0; El < 3; El++) {
			for (unsigned m = 0; m < 3; m++) {
				for (unsigned j = 0; j < 8; j++) {
					real integral = 0.0;
					for (unsigned n = 0; n < NumberOfQuadraturePointsEdge; n++) {
						real const tC = upwindConcentration8(K, El, n);
						integral += tC * QuadraturePoints_RaviartThomasBasisDotNormalTimesPolynomialBasis(n, j, El, m);
					}
					δ.setCoeff(k_index, m, El, j) = abs(integral) < INTEGRAL_PRECISION ? 0.0 : integral;
				}
			}
		}*/

		for (unsigned El = 0; El < 3; El++)
			for (unsigned m = 0; m < 3; m++)
				for (unsigned j = 0; j < 8; j++)
					δ.setCoeff(k_index, m, El, j) = abs(δ(k_index, m, El, j)) < INTEGRAL_PRECISION ? 0.0 : δ(k_index, m, El, j);

	}

};

template<unsigned QuadraturePrecision>
void solver<QuadraturePrecision>::assemble_σ() {



	for (unsigned k = 0; k < nk; k++) {

		unsigned const k_index = ElementIndeces[k];

		real const coeff = betas_prev[k_index] / (viscosities[k_index] * porosities[k_index]);


		for (unsigned m = 0; m < 3; m++) {

			for (unsigned l = 0; l < 3; l++) {


				real val = 0.0;
				for (unsigned q = 0; q < 3; q++) {


					real GAB = 0.0;
					for (unsigned j = 0; j < 8; j++) {

						real const AB = AlphaTimesBeta.setCoeff(k_index, j, l);

						GAB += γ(k_index, q, j) * AB;

					}

					val += η(m, q) * GAB;

				}

				// In the computation of Eta, there is coefficient detJF. When inverting, the coefficient is inverted
				σ.setCoeff(k_index, m, l) = -coeff * val / AffineMappingMatrixDeterminant[k_index];

			}
		}
	}

};
template<unsigned QuadraturePrecision>
void solver<QuadraturePrecision>::assemble_λ() {


	for (unsigned k = 0; k < nk; k++) {

		unsigned const k_index = ElementIndeces[k];

		real const coeff = betas_prev[k_index] / (viscosities[k_index] * porosities[k_index]);


		for (unsigned s = 0; s < 2; s++) {

			for (unsigned m = 0; m < 3; m++) {

				for (unsigned El = 0; El < 3; El++) {


					real val = 0.0;
					for (unsigned q = 0; q < 3; q++) {


						real YACHI = 0.0;
						for (unsigned j = 0; j < 8; j++) {

							real const ACHI = AlphaTimesChi(k_index, j, El, s);

							YACHI += γ(k_index, q, j) * ACHI;

						}

						val += η(m, q) * YACHI;

					}

					// In the computation of Eta, there is coefficient detJF. When inverting, the coefficient is inverted
					λ.setCoeff(k_index, s, m, El) = coeff * val / AffineMappingMatrixDeterminant[k_index];

				}
			}
		}
	}

};
template<unsigned QuadraturePrecision>
void solver<QuadraturePrecision>::assemble_φ() {};
template<unsigned QuadraturePrecision>
void solver<QuadraturePrecision>::assemble_ψ() {};




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


	std::ofstream txtFile;
	txtFile.open("C:\\Users\\pgali\\Desktop\\flow2d\\v.txt");
	for (unsigned k = 0; k < nk; k++) {
	
		for (unsigned j = 0; j < 8; j++)
			txtFile << υ(mesh->get_triangle(k)->index, j) << " ";

		txtFile << std::endl;
	}
	txtFile.close();

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