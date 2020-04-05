
#pragma once



/*****************************************************************************/
/*                                                                           */
/*    - Preprocessor flags for eigen										 */
/*                                                                           */
/*****************************************************************************/
#define EIGEN_NO_DEBUG
#define NDEBUG
#define EIGEN_UNROLLING_LIMIT 500




#include "miscellaneous.h"
#include "coefficient_matrix.h"
#include "matrix.h"
#include "mesh.h"

//#include <Eigen/UmfPackSupport>
#include <Eigen/Sparse>
#include <omp.h>
#include <immintrin.h>


enum scheme { CRANK_NICOLSON, EULER_BACKWARD };

	typedef Eigen::SparseMatrix<Real>				SparseMatrix;
	typedef Eigen::Matrix<Real, Eigen::Dynamic, 1>	DenseVector;


/*****************************************************************************/
/*                                                                           */
/*    - Get number of quadrature points in compilation						 */
/*                                                                           */
/*****************************************************************************/
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
	case 10:
		return 10;	// order 10 not implemented -> using order 9
	case 11:
		return 11;
	case 12:
		return 12;	// order 12 not implemented -> using order 11
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
	case 10:
		return 19;	// order 10 not implemented -> using order 9
	case 11:
		return 28;
	case 12:
		return 28;	// order 12 not implemented -> using order 11
	case 13:
		return 37;

	default:
		return 0;

	}

};



template <unsigned QuadraturePrecision = 6, scheme TimeScheme = CRANK_NICOLSON>
class solver {






	/*****************************************************************************/
	/*                                                                           */
	/*    - Numerical model parameters									         */
	/*                                                                           */
	/*		- TimeSchemeParameter : 0.5 = Crank-Nicolson, 1.0 = Backward Euler	 */
	/*		- INTEGRAL_PRECISION  : Values lower than this are treated as zero	 */
	/*							  : Used in numerical integration etc.			 */
	/*		- TOLERANCE			  : Tolerance for the inner loop				 */
	/*							  : Succesive updates of P,C,Theta differ		 */
	/*							    less than this bound						 */
	/*		- MAX_ITERATIONS	  : Maximumnumber of iterations of the inner	 */
	/*							    loop, which updates P,C,Theta				 */
	/*                                                                           */
	/*****************************************************************************/
	Real const TimeSchemeParameter		= TimeScheme == CRANK_NICOLSON ? 0.5 : 1.0;
	Real const INTEGRAL_PRECISION		= 1e-11;
	Real const TOLERANCE				= DBL_EPSILON;
	unsigned const MAX_ITERATIONS		= 100;



public:


	void getSolution();

	void setTimeStep(Real const _dt) { this->dt = _dt; };
	void setTimeLevel(int const _nt) { this->nt = _nt; };


	void exportPressures(std::string const & fileName);
	void exportConcentrations(std::string const & fileName);

	void exportTracePressures(std::string const & fileName);

	void exportVelocityField(std::string const & fileName);
	void exportVelocities(std::string const & fileName);

	void computeError(std::string const & fileName);


	solver(MESH & mesh, int const nt0, Real const dt0);
	~solver();


private:


	/*****************************************************************************/
	/*                                                                           */
	/*    - nk : number of elements in the Mesh							         */
	/*      ne : number of edges in the Mesh							         */
	/*                                                                           */
	/*    - nt : Current time level of the solution						         */
	/*      dt : Time step												         */
	/*                                                                           */
	/*    - NumberOfQuadraturePointsEdge : number of quad points on edge		 */
	/*      NumberOfQuadraturePointsTriangle : number of quad points on element  */
	/*                                                                           */
	/*****************************************************************************/
	unsigned const nk;
	unsigned const ne;

	int		nt;
	Real	dt;

	static unsigned const NumberOfQuadraturePointsEdge		= get_number_of_quadrature_points_edge<QuadraturePrecision>();
	static unsigned const NumberOfQuadraturePointsTriangle	= get_number_of_quadrature_points_triangle<QuadraturePrecision>();



	/*****************************************************************************/
	/*                                                                           */
	/*    - Unknowns quantities of the model   : Pi			- Internal Pressures */
	/*										   : TPi		- Trace pressures	 */
	/*										   : Xi			- Concentrations     */
	/*										   : Velocity	- Velocities		 */
	/*                                                                           */
	/*										   : Sources	- Sources / Sinks	 */
	/*                                                                           */
	/*****************************************************************************/
	CoeffMatrix1D<3>	Pi;
	CoeffMatrix1D<3>	Xi;
	CoeffMatrix2D<3, 2> TPi;
	CoeffMatrix1D<8>	Velocity;
	CoeffMatrix1D<3>	Sources;


	/*****************************************************************************/
	/*                                                                           */
	/*    - Unknowns quantities on different time levels					     */
	/*			: 'n'		- current time level                                 */
	/*			: 'prev'	- previous iteration in the inner loop               */
	/*                                                                           */
	/*****************************************************************************/
	CoeffMatrix1D<3> Pi_n;
	CoeffMatrix1D<3> Pi_prev;

	CoeffMatrix1D<3> Xi_n;
	CoeffMatrix1D<3> Xi_prev;

	CoeffMatrix1D<3> rkFc;
	CoeffMatrix1D<3> rkFc_n;

	CoeffMatrix1D<3> rkFp_n;


	/*****************************************************************************/
	/*                                                                           */
	/*    - Integral geometric coefficients								         */
	/*                                                                           */
	/*****************************************************************************/
	CoeffMatrix2D<8, 8>				Alpha;
	CoeffMatrix3D<3, 3, 8>			Delta;
	CoeffMatrix2D<3, 8>				Gamma;
	SingleCoeffMatrix2D<8, 3>		Beta;
	CoeffMatrix3D<8, 3, 2>			Chi;
	SingleCoeffMatrix2D<3, 3>		Eta;
	SingleCoeffMatrix3D<3, 8, 3>	Tau;


	CoeffMatrix2D<3, 3>		Sigma;
	CoeffMatrix3D<2, 3, 3>	Lambda;
	//CoeffMatrix1D<3>		BigPhi;
	//CoeffMatrix1D<3>		BigPsi;
	//CoeffMatrix1D<8>		BigOmega;


	/*****************************************************************************/
	/*                                                                           */
	/*    - Precomputed values of basis function in quadrature points	         */
	/*                                                                           */
	/*    - Quadrature points on the edges of the element				         */
	/*       : QuadraturePoints_RaviartThomasBasis								 */
	/*       : QuadraturePoints_PolynomialBasis									 */
	/*       : QuadraturePoints_PolynomialBasis									 */
	/*       : QuadraturePoints_PolynomialBasisOpposite							 */
	/*                * Used in Upwinding from the neighboring element		     */
	/*                * n-th quad point is on the opposite sides of an edge      */
	/*       : QuadraturePoints_RaviartThomasBasisDotNormalTimesPolynomialBasis  */
	/*                                                                           */
	/*    - Quadrature points on the element							         */
	/*       : QuadraturePoints_PhysicalNormalDotPhysicalRaviartThomasBasis	     */
	/*       : AlphaTimesChi													 */
	/*       : AlphaTimesBeta													 */
	/*                                                                           */
	/*****************************************************************************/
	CoeffMatrix1D<3>									QuadraturePointsAndWeightsOnReferenceTriangle;
	CoeffMatrix1D<3>									QuadraturePoints_PolynomialBasisOnReferenceTriangle;

	//CoeffMatrix1D<3, NumberOfQuadraturePointsTriangle>	QuadraturePointsAndWeightsOnReferenceTriangle;
	//CoeffMatrix1D<3, NumberOfQuadraturePointsTriangle>	QuadraturePoints_PolynomialBasisOnReferenceTriangle;


	CoeffMatrix3D<3, 8, 2>								QuadraturePoints_RaviartThomasBasis;
	CoeffMatrix2D<3, 3>									QuadraturePoints_PolynomialBasis;
	CoeffMatrix2D<3, 3>									QuadraturePoints_PolynomialBasisOpposite;
	CoeffMatrix3D<8, 3, 3>								QuadraturePoints_RaviartThomasBasisDotNormalTimesPolynomialBasis;
	CoeffMatrix3D<3, NumberOfQuadraturePointsEdge, 8>	QuadraturePoints_PhysicalNormalDotPhysicalRaviartThomasBasis;
	CoeffMatrix3D<8, 3, 2>								AlphaTimesChi;
	CoeffMatrix2D<8, 3>									AlphaTimesBeta;


	/*****************************************************************************/
	/*                                                                           */
	/*    - Quadrature points on the edges of the physical triangle		         */
	/*                                                                           */
	/*****************************************************************************/
	CoeffMatrix2D<3, NumberOfQuadraturePointsEdge>	 QuadraturePoints_Edge_x;
	CoeffMatrix2D<3, NumberOfQuadraturePointsEdge>	 QuadraturePoints_Edge_y;


	/*****************************************************************************/
	/*                                                                           */
	/*    - Some auxilary things					   							 */
	/*                                                                           */
	/*****************************************************************************/
	CoeffMatrix1D<3> edgeOrientation;

	tm_pointer	* MeshElements				= NULL;
	em_pointer	* MeshEdges					= NULL;
	unsigned	* MeshElementIndeces		= NULL;
	unsigned	* MeshEdgesIndeces			= NULL;

	Real * AffineMappingMatrixDeterminant	= NULL;
	Real * PorosityViscosityDeterminant		= NULL;

	/*
	em_pointer * EdgesDirichlet = NULL;
	em_pointer * EdgesNeumann	= NULL;

	unsigned * EdgesDirichlet_Indeces	= NULL;
	unsigned * EdgesNeumann_Indeces		= NULL;

	unsigned DirichletCount = 0;
	unsigned NeumannCount	= 0;

	for (unsigned e = 0; e < ne; e++){

		em_pointer const E = Mesh->get_edge(e);

		E_MARKER const e_marker = E->marker;
		unsigned const e_index = E->index;

		if (e_marker == E_MARKER::DIRICHLET){

			EdgesDirichlet[DirichletCount]			= E;
			EdgesDirichlet_Indeces[DirichletCount]	= e_index;

			DirichletCount++;

		}
		else if (e_marker == E_MARKER::NEUMANN){

			EdgesNeumann[NeumannCount]			= E;
			EdgesNeumann_Indeces[NeumannCount]	= e_index;

			NeumannCount++;

		}
	}
	*/


	/*****************************************************************************/
	/*                                                                           */
	/*    - Porous medium physical quantities and coefficients Thetas, which	 */
	/*		are connected to Equation Of State (EOS)							 */
	/*                                                                           */
	/*****************************************************************************/
	Real * Viscosities;
	Real * Porosities;

	Real * Thetas;
	Real * Thetas_prev;


	/*****************************************************************************/
	/*                                                                           */
	/*    - Pressure and Trace Pressure system matrices							 */
	/*                                                                           */
	/*****************************************************************************/
	SparseMatrix	R1;
	SparseMatrix	R2;
	SparseMatrix	M_j1_s1;
	SparseMatrix	M_j1_s2;
	SparseMatrix	M_j2_s1;
	SparseMatrix	M_j2_s2;

	DenseVector		V1;
	DenseVector		V2;
	DenseVector		G;

	DenseVector		Tp;
	DenseVector		PiTemp;

	SparseMatrix	iD;
	SparseMatrix	H1;
	SparseMatrix	H2;

	SparseMatrix	R1iD;
	SparseMatrix	R2iD;
	SparseMatrix	iDH1;
	SparseMatrix	iDH2;

//DenseVector		iDG;
//DenseVector		R1iDG;
//DenseVector		R2iDG;



	CoeffMatrix2D<3, 3> R1_block;
	CoeffMatrix2D<3, 3> R2_block;


	/*****************************************************************************/
	/*                                                                           */
	/*    - Eigen solver for the computation of LU factorization				 */
	/*                                                                           */
	/*****************************************************************************/
	//Eigen::UmfPackLU<SparseMatrix>	LUFacorizationUMFPACK;

	Eigen::BiCGSTAB<SparseMatrix, Eigen::DiagonalPreconditioner<Real>>	BiConjugateGradientSolver;
	Eigen::SparseLU<SparseMatrix>										sparseLUsolver_TracePressureSystem;
	Eigen::SparseLU<SparseMatrix>										sparseLUsolver_PressureSystem;
	SparseMatrix														PressureSystem;

	std::vector<Eigen::Triplet<Real>>	TripletVector;
	std::vector<Eigen::Triplet<Real>>	TripletVectorInverseD;

	DenseVector traceSystemRhs;
	DenseVector pressureSystemRhs;

	//Eigen::COLAMDOrdering<int>										SparseOrdering;
	//Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic, int>	PermutationMatrix;


	/*****************************************************************************/
	/*                                                                           */
	/*    - Class methods														 */
	/*                                                                           */
	/*****************************************************************************/
	void initializeValues();
	void computeThetas();

	bool stopCriterion();
	void concentrationCorrection();

	void computeTracePressures();
	void computePressureEquation();

	void computeVelocities();
	void updateConcentrations();

	Real upwindConcentration(tm_pointer const & K, unsigned const El, unsigned const n);
	Real upwindConcentration_limiter(tm_pointer const & K, unsigned const El, unsigned const n);

	void assembleR();
	void assembleM();
	void assembleV();

	void assembleInverseD();
	void assembleH();
	void assembleG();

	void assemblePressureSystem();
	void getSparsityPatternOfThePressureSystem();

	void assemble_Alpha();
	void assemble_Beta();
	void assemble_Chi();
	void assemble_Eta();
	void assemble_Tau();
	void assemble_Delta();
	void assemble_Gamma();

	void assemble_Sigma();
	void assemble_Lambda();
	void assemble_BigPhi();
	void assemble_BigPsi();
	void assemble_BigOmega();


	/*****************************************************************************/
	/*                                                                           */
	/*    - Evaluation of Raviart-Thomas, P1 basis function, Edge basis P1(E),	 */
	/*      edges parametrization, normal vectors							     */
	/*                                                                           */
	/*****************************************************************************/
	inline void evaluate_raviartthomas_basis(Real const s, Real const t, Eigen::Matrix<Real, 2, 8> & out);
	inline void evaluate_raviartthomas_basis_divergence(Real const s, Real const t, Eigen::Matrix<Real, 8, 1> & out);

	inline void evaluate_polynomial_basis(Real const s, Real const t, Eigen::Matrix<Real, 3, 1> & out);
	inline void evaluate_polynomial_basis_gradient(Real const s, Real const t, Eigen::Matrix<Real, 2, 3> & out);

	inline void evaluate_edge_polynomial_basis(Real const ksi, unsigned const El, Eigen::Matrix<Real, 2, 1> & out, Real const & orientation);

	inline void evaluate_edge_normal(Eigen::Matrix<Real, 2, 3> & out);
	inline void evaluate_edge_parametrization(Real const ksi, unsigned const El, Eigen::Matrix<Real, 2, 1> & out);
	inline void evaluate_edge_parametrization_derivative(Real const ksi, unsigned const El, Eigen::Matrix<Real, 2, 1> & out);
	inline void evaluate_edge_parametrization_opposite(Real const ksi, unsigned const El, Eigen::Matrix<Real, 2, 1> & out);
	inline void evaluate_edge_parametrization_opposite_derivative(Real const ksi, unsigned const El, Eigen::Matrix<Real, 2, 1> & out);


	inline Real phi0(Real const s, Real const t);
	inline Real phi1(Real const s, Real const t);
	inline Real phi2(Real const s, Real const t);


	/*****************************************************************************/
	/*                                                                           */
	/*    - Integral of Source * Polynomial basis function						 */
	/*                                                                           */
	/*****************************************************************************/
	Real F0(tm_pointer const K, Real const time);
	Real F1(tm_pointer const K, Real const time);
	Real F2(tm_pointer const K, Real const time);

	void evaluate_source_integrals();



	/*****************************************************************************/
	/*                                                                           */
	/*    - Return the number of the basis function on the E of the K			 */
	/*			: (0 if DOF == 0), (1 if DOF == 1)								 */
	/*			: E0 -> 0,3														 */
	/*			: E1 -> 1,4														 */
	/*			: E2 -> 2,5														 */
	/*	  - 6,7 are defined inside element and are not connected to any edge	 */
	/*                                                                           */
	/*****************************************************************************/
	inline unsigned LI(tm_pointer const & K, em_pointer const & E, unsigned const & DOF) {

		return K->get_edge_index(E) + 3 * DOF;

	};
	inline void LI(tm_pointer const & K, em_pointer const & E, unsigned(&out)[2]) {

		unsigned const e_index_local = K->get_edge_index(E);

		out[0] = e_index_local;
		out[1] = e_index_local + 3;

	};


	SparseMatrix P00;
	SparseMatrix P01;
	SparseMatrix P10;
	SparseMatrix P11;

	std::vector<Eigen::Triplet<Real>> SystemTriplet;
};



template<unsigned QuadraturePrecision, scheme TimeScheme>
solver<QuadraturePrecision, TimeScheme>::solver(MESH & mesh, int const nt0, Real const dt0) : nk(mesh.get_number_of_triangles()), ne(mesh.get_number_of_edges()), nt(nt0), dt(dt0) {


	/*****************************************************************************/
	/*                                                                           */
	/*    - Memory allocation for the physical quantities   			         */
	/*                                                                           */
	/*****************************************************************************/
	Viscosities = new Real[nk];
	Porosities	= new Real[nk];


	/*****************************************************************************/
	/*                                                                           */
	/*    - Memory allocation for the beta coefficient. This coefficient         */
	/*      links model equations and EOS                                        */
	/*                                                                           */
	/*****************************************************************************/
	Thetas		= new Real[nk];
	Thetas_prev = new Real[nk];


	/*****************************************************************************/
	/*                                                                           */
	/*    - Memory allocation for the unknowns : Pi			- Internal Pressures */
	/*										   : TPi		- Trace pressures	 */
	/*										   : Xi			- Concentrations     */
	/*										   : Velocity	- Velocities		 */
	/*                                                                           */
	/*										   : Sources	- Sources / Sinks	 */
	/*                                                                           */
	/*****************************************************************************/
	Pi		.setNumberOfElements(nk);
	TPi		.setNumberOfElements(nk);
	Xi		.setNumberOfElements(nk);
	Velocity.setNumberOfElements(nk);
	Sources	.setNumberOfElements(nk);

	Pi		.setZero();
	TPi		.setZero();
	Xi		.setZero();
	Velocity.setZero();
	Sources	.setZero();


	/*****************************************************************************/
	/*                                                                           */
	/*    - Memory allocation for the quantities on different time levels	     */
	/*			: 'n'		- current time level                                 */
	/*			: 'prev'	- previous iteration in the inner loop               */
	/*                                                                           */
	/*****************************************************************************/
	Pi_n	.setNumberOfElements(nk);
	Pi_prev	.setNumberOfElements(nk);
	Xi_n	.setNumberOfElements(nk);
	Xi_prev	.setNumberOfElements(nk);
	rkFc	.setNumberOfElements(nk);
	rkFc_n	.setNumberOfElements(nk);
	rkFp_n	.setNumberOfElements(nk);

	Pi_n	.setZero();
	Pi_prev	.setZero();
	Xi_n	.setZero();
	Xi_prev	.setZero();
	rkFc	.setZero();
	rkFc_n	.setZero();
	rkFp_n	.setZero();


	/*****************************************************************************/
	/*                                                                           */
	/*    - Memory allocation for the integral geometric coefficients	         */
	/*                                                                           */
	/*****************************************************************************/
	Alpha	.setNumberOfElements(nk);
	Gamma	.setNumberOfElements(nk);
	Delta	.setNumberOfElements(nk);
	Chi		.setNumberOfElements(nk);

	Sigma	.setNumberOfElements(nk);
	Lambda	.setNumberOfElements(nk);

	//BigPhi	.setNumberOfElements(nk);
	//BigPsi	.setNumberOfElements(nk);
	//BigOmega.setNumberOfElements(nk);

	Alpha	.setZero();
	Beta	.setZero();
	Chi		.setZero();
	Eta		.setZero();
	Tau		.setZero();
	Gamma	.setZero();
	Delta	.setZero();

	Sigma	.setZero();
	Lambda	.setZero();

	//BigPhi	.setZero();
	//BigPsi	.setZero();
	//BigOmega.setZero();


	/*****************************************************************************/
	/*                                                                           */
	/*    - Memory allocation for the linear system matrices			         */
	/*                                                                           */
	/*****************************************************************************/
	PressureSystem		.resize(2 * ne, 2 * ne);
	pressureSystemRhs	.resize(2 * ne);
	traceSystemRhs		.resize(2 * ne);

	R1					.resize(ne, 3 * nk);
	R2					.resize(ne, 3 * nk);
	M_j1_s1				.resize(ne, ne);
	M_j1_s2				.resize(ne, ne);
	M_j2_s1				.resize(ne, ne);
	M_j2_s2				.resize(ne, ne);

	V1					.resize(ne);
	V2					.resize(ne);
	G					.resize(3 * nk);

	Tp					.resize(2 * ne);
	PiTemp			.resize(3 * nk);

	iD					.resize(3 * nk, 3 * nk);
	H1		.resize(3 * nk, ne);
	H2		.resize(3 * nk, ne);
	R1iD				.resize(ne, 3 * nk);
	R2iD				.resize(ne, 3 * nk);

	iDH1				.resize(3 * nk, ne);
	iDH2				.resize(3 * nk, ne);

	R1_block			.setNumberOfElements(nk);
	R2_block			.setNumberOfElements(nk);

	//iDG.resize(3 * nk);
	//R1iDG.resize(ne);
	//R2iDG.resize(ne);



	P00.resize(ne, ne);
	P01.resize(ne, ne);
	P10.resize(ne, ne);
	P11.resize(ne, ne);


	/*****************************************************************************/
	/*                                                                           */
	/*    - Some auxilary and test things					   			         */
	/*                                                                           */
	/*****************************************************************************/
	MeshElements					= new tm_pointer[nk];
	MeshEdges						= new em_pointer[ne];
	MeshElementIndeces				= new unsigned[nk];
	MeshEdgesIndeces				= new unsigned[ne];
	AffineMappingMatrixDeterminant	= new Real[nk];
	PorosityViscosityDeterminant	= new Real[nk];

	edgeOrientation.setNumberOfElements(nk);


	/*****************************************************************************/
	/*                                                                           */
	/*    - This auxilary variable helps compute Edge integral coefficient Chi.	 */
	/*      Edge basis function are defined globaly on each edge. From each      */
	/*      element these basis functions must be the same                       */
	/*                                                                           */
	/*    - Make local copies of pointers to the mesh primitives 			     */
	/*                                                                           */
	/*****************************************************************************/
	for (unsigned k = 0; k < nk; k++) {


		tm_pointer const K			= mesh.get_triangle(k);
		unsigned const	 k_index	= K->index;

		MeshElements[k]			= K;
		MeshElementIndeces[k]	= k_index;

		for (unsigned El = 0; El < 3; El++) {


			em_pointer const E		 = K->edges[El];
			unsigned const	 e_index = E->index;

			vm_pointer const v = K->get_vertex_cw(E->a);
			vm_pointer const p = K->get_vertex(El);

			if (v != p) edgeOrientation.setCoeff(k_index, El) = -1.0;
			else		edgeOrientation.setCoeff(k_index, El) = +1.0;

		}
	}


	/*****************************************************************************/
	/*                                                                           */
	/*    - Make local copies of pointers to the mesh primitives 			     */
	/*                                                                           */
	/*****************************************************************************/
	for (unsigned e = 0; e < ne; e++) {

		em_pointer const E = mesh.get_edge(e);

		MeshEdges[e]		= E;
		MeshEdgesIndeces[e] = E->index;

	}


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
	assemble_Alpha();
	assemble_Beta();
	assemble_Chi();
	assemble_Eta();
	assemble_Tau();


	/*****************************************************************************/
	/*                                                                           */
	/*    - Memory allocation for precomputed values of basis function			 */
	/* 		in quadrature points												 */
	/*                                                                           */
	/*****************************************************************************/
	QuadraturePoints_RaviartThomasBasis									.setNumberOfElements(NumberOfQuadraturePointsEdge);
	QuadraturePoints_PolynomialBasis									.setNumberOfElements(NumberOfQuadraturePointsEdge);
	QuadraturePoints_PolynomialBasisOpposite							.setNumberOfElements(NumberOfQuadraturePointsEdge);
	QuadraturePoints_RaviartThomasBasisDotNormalTimesPolynomialBasis	.setNumberOfElements(NumberOfQuadraturePointsEdge);
	QuadraturePoints_PhysicalNormalDotPhysicalRaviartThomasBasis		.setNumberOfElements(nk);
	AlphaTimesChi														.setNumberOfElements(nk);
	AlphaTimesBeta														.setNumberOfElements(nk);

	QuadraturePoints_Edge_x												.setNumberOfElements(nk);
	QuadraturePoints_Edge_y												.setNumberOfElements(nk);


	QuadraturePoints_RaviartThomasBasis									.setZero();
	QuadraturePoints_PolynomialBasis									.setZero();
	QuadraturePoints_PolynomialBasisOpposite							.setZero();
	QuadraturePoints_RaviartThomasBasisDotNormalTimesPolynomialBasis	.setZero();
	QuadraturePoints_PhysicalNormalDotPhysicalRaviartThomasBasis		.setZero();
	AlphaTimesChi														.setZero();
	AlphaTimesBeta														.setZero();

	QuadraturePoints_Edge_x												.setZero();
	QuadraturePoints_Edge_y												.setZero();


	/*****************************************************************************/
	/*                                                                           */
	/*    - Some auxilary variables							   			         */
	/*                                                                           */
	/*****************************************************************************/
	quadrature_triangle<Real> const GaussQuadratureOnTriangle(QuadraturePrecision);
	gauss_quadrature_1D<Real> const GaussQuadratureOnEdge(QuadraturePrecision);

	Eigen::Matrix<Real, 2, 2> JF;
	Eigen::Matrix<Real, 2, 3> ReferenceNormals;
	Eigen::Matrix<Real, 2, 8> BasisRaviartThomas;

	Eigen::Matrix<Real, 2, 1> Parametrization;
	Eigen::Matrix<Real, 2, 1> ParametrizationOpposite;
	Eigen::Matrix<Real, 3, 1> BasisPolynomial;
	Eigen::Matrix<Real, 3, 1> BasisPolynomialOpposite;

	evaluate_edge_normal(ReferenceNormals);


	/*****************************************************************************/
	/*                                                                           */
	/*    - Precomputation of the quadrature points on the reference triangle	 */
	/*                                                                           */
	/*****************************************************************************/
	QuadraturePointsAndWeightsOnReferenceTriangle			.setNumberOfElements(NumberOfQuadraturePointsTriangle);
	QuadraturePoints_PolynomialBasisOnReferenceTriangle		.setNumberOfElements(NumberOfQuadraturePointsTriangle);

	QuadraturePointsAndWeightsOnReferenceTriangle			.setZero();
	QuadraturePoints_PolynomialBasisOnReferenceTriangle		.setZero();

	for (unsigned n = 0; n < NumberOfQuadraturePointsTriangle; n++) {

		Real const s = GaussQuadratureOnTriangle.points_x[n];
		Real const t = GaussQuadratureOnTriangle.points_y[n];
		Real const w = GaussQuadratureOnTriangle.weights[n];

		QuadraturePointsAndWeightsOnReferenceTriangle.setCoeff(n, 0) = s;
		QuadraturePointsAndWeightsOnReferenceTriangle.setCoeff(n, 1) = t;
		QuadraturePointsAndWeightsOnReferenceTriangle.setCoeff(n, 2) = w;

		evaluate_polynomial_basis(s, t, BasisPolynomial);

		QuadraturePoints_PolynomialBasisOnReferenceTriangle.setCoeff(n, 0) = BasisPolynomial(0);
		QuadraturePoints_PolynomialBasisOnReferenceTriangle.setCoeff(n, 1) = BasisPolynomial(1);
		QuadraturePoints_PolynomialBasisOnReferenceTriangle.setCoeff(n, 2) = BasisPolynomial(2);

	}


	/*****************************************************************************/
	/*                                                                           */
	/*    - Precomputation of values defined on the Elements				     */
	/*                                                                           */
	/*****************************************************************************/
	for (unsigned k = 0; k < nk; k++) {


		tm_pointer const K			= MeshElements[k];
		unsigned const	 k_index	= K->index;


		vm_pointer const va = K->vertices[0];
		vm_pointer const vb = K->vertices[1];
		vm_pointer const vc = K->vertices[2];

		Real const x0 = va->x;
		Real const y0 = va->y;

		Real const x1 = vb->x;
		Real const y1 = vb->y;

		Real const x2 = vc->x;
		Real const y2 = vc->y;

		Real const detJF = abs((x1 - x0)*(y2 - y0) - (x2 - x0)*(y1 - y0));


		/*****************************************************************************/
		/*                                                                           */
		/*    - Precomputed values of : Determinant of the affine mapping matrix JF  */
		/*							  : Affine mapping matrix JF				     */
		/*                                                                           */
		/*****************************************************************************/
		AffineMappingMatrixDeterminant[k_index] = detJF;
		PorosityViscosityDeterminant[k_index]	= 1.0 / (Viscosities[k_index] * Porosities[k_index] * detJF);

		JF(0, 0) = x1 - x0;
		JF(0, 1) = x2 - x0;
		JF(1, 0) = y1 - y0;
		JF(1, 1) = y2 - y0;

		Eigen::Matrix<Real, 2, 2> const itJF = (JF.inverse()).transpose();


		/*****************************************************************************/
		/*                                                                           */
		/*    - Precomputed values of : Alpha * Beta							     */
		/*                                                                           */
		/*    - Used in : Sigma computation											 */
		/*                                                                           */
		/*****************************************************************************/
		for (unsigned j = 0; j < 8; j++) {
			for (unsigned l = 0; l < 3; l++) {

				Real AlphaBeta = 0.0;

				for (unsigned i = 0; i < 8; i++)
					AlphaBeta += Alpha(k_index, j, i) * Beta(i, l);

				AlphaTimesBeta.setCoeff(k_index, j, l) = AlphaBeta;

			}
		}


		for (unsigned El = 0; El < 3; El++) {


			em_pointer const E			= K->edges[El];
			unsigned const	 e_index	= E->index;
			E_MARKER const	 e_marker	= E->marker;


			Real const a = 0.0;
			Real const b = El != 0 ? 1.0 : sqrt(2.0);

			Real const c = 0.5 * (b - a);
			Real const d = 0.5 * (b + a);


			//Vector<real> const normal = normals.getColumn(El);
			Eigen::Matrix<Real, 2, 1> const PhysicalNormal = (itJF * ReferenceNormals.col(El)) / (itJF * ReferenceNormals.col(El)).norm();



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


				Real const x = GaussQuadratureOnEdge.points[n] * c + d;
				Real const w = GaussQuadratureOnEdge.weights[n] * c;


				evaluate_edge_parametrization(x, El, Parametrization);


				Real const s = Parametrization(0);
				Real const t = Parametrization(1);

				Real const X = x0 + JF(0, 0) * s + JF(0, 1) * t;
				Real const Y = y0 + JF(1, 0) * s + JF(1, 1) * t;

				evaluate_raviartthomas_basis(s, t, BasisRaviartThomas);


				QuadraturePoints_Edge_x.setCoeff(k_index, El, n) = abs(X) < INTEGRAL_PRECISION ? 0.0 : X;
				QuadraturePoints_Edge_y.setCoeff(k_index, El, n) = abs(Y) < INTEGRAL_PRECISION ? 0.0 : Y;


				for (unsigned j = 0; j < 8; j++) {

					Real const Value = PhysicalNormal.dot(JF * BasisRaviartThomas.col(j)) / detJF;
					QuadraturePoints_PhysicalNormalDotPhysicalRaviartThomasBasis.setCoeff(k_index, El, n, j) = abs(Value) < INTEGRAL_PRECISION ? 0.0 : Value;

				}
			}



			/*****************************************************************************/
			/*                                                                           */
			/*    - Precomputed values of : Alpha * Chi								     */
			/*                                                                           */
			/*    - Used in : computation of Lambda				                         */
			/*                                                                           */
			/*****************************************************************************/
			for (unsigned j = 0; j < 8; j++) {
				for (unsigned s = 0; s < 2; s++) {

					Real AlphaChi = 0.0;

					for (unsigned i = 0; i < 8; i++)
						AlphaChi += Alpha(k_index, j, i) * Chi(k_index, i, El, s);

					AlphaTimesChi.setCoeff(k_index, j, El, s) = AlphaChi;

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

				R1_block.setCoeff(k_index, El, 0) = 0.0;
				R1_block.setCoeff(k_index, El, 1) = 0.0;
				R1_block.setCoeff(k_index, El, 2) = 0.0;

				R2_block.setCoeff(k_index, El, 0) = 0.0;
				R2_block.setCoeff(k_index, El, 1) = 0.0;
				R2_block.setCoeff(k_index, El, 2) = 0.0;

				continue;

			}


			unsigned const dof0 = LI(K, E, 0);
			unsigned const dof1 = LI(K, E, 1);

			Real const ChiCoeff0 = Chi(k_index, dof0, El, 0);
			Real const ChiCoeff1 = Chi(k_index, dof1, El, 1);


			for (unsigned m = 0; m < 3; m++) {


				Real AlphaBeta1 = 0.0;
				Real AlphaBeta2 = 0.0;

				for (unsigned i = 0; i < 8; i++) {

					AlphaBeta1 += Alpha(k_index, dof0, i) * Beta(i, m);
					AlphaBeta2 += Alpha(k_index, dof1, i) * Beta(i, m);

				}

				Real const Value1 = ChiCoeff0 * AlphaBeta1 / Viscosities[k_index];
				Real const Value2 = ChiCoeff1 * AlphaBeta2 / Viscosities[k_index];

				R1_block.setCoeff(k_index, El, m) = abs(Value1) < INTEGRAL_PRECISION ? 0.0 : Value1;
				R2_block.setCoeff(k_index, El, m) = abs(Value2) < INTEGRAL_PRECISION ? 0.0 : Value2;

			}
		}
	}



	/*****************************************************************************/
	/*                                                                           */
	/*    - Precomputation of values defined on the edges					     */
	/*                                                                           */
	/*****************************************************************************/
	for (unsigned El = 0; El < 3; El++) {


		Real const a = 0.0;
		Real const b = El != 0 ? 1.0 : sqrt(2.0);

		Real const c = 0.5 * (b - a);
		Real const d = 0.5 * (b + a);

		//Vector<real> const normal = ReferenceNormals.getColumn(El);
		Eigen::Matrix<Real, 2, 1> const ReferenceNormal = ReferenceNormals.col(El);


		for (unsigned n = 0; n < NumberOfQuadraturePointsEdge; n++) {


			Real const x = GaussQuadratureOnEdge.points[n] * c + d;
			Real const w = GaussQuadratureOnEdge.weights[n] * c;


			/*****************************************************************************/
			/*                                                                           */
			/*    - Precomputed values of : Quadrature points on each edge of the	     */
			/*								reference triangle						     */
			/*                                                                           */
			/*****************************************************************************/
			evaluate_edge_parametrization(x, El, Parametrization);
			evaluate_edge_parametrization_opposite(x, El, ParametrizationOpposite);

			Real const s = Parametrization(0);
			Real const t = Parametrization(1);

			Real const s_opposite = ParametrizationOpposite(0);
			Real const t_opposite = ParametrizationOpposite(1);

			evaluate_raviartthomas_basis(s, t, BasisRaviartThomas);
			evaluate_polynomial_basis(s, t, BasisPolynomial);
			evaluate_polynomial_basis(s_opposite, t_opposite, BasisPolynomialOpposite);


			/*****************************************************************************/
			/*                                                                           */
			/*    - Precomputed values of : RT basis on the reference triangle		     */
			/*							  : Polynomial basis on the reference triangle	 */
			/*							  : RT basis * Normal * Polynomial basis		 */
			/*                                                                           */
			/*****************************************************************************/
			for (unsigned j = 0; j < 8; j++) {


				//Vector<real> const Wj = BasisRaviartThomas.getColumn(j);
				Eigen::Matrix<Real, 2, 1> const Wj = BasisRaviartThomas.col(j);

				//real const dotProduct = dot(Wj, normal);
				Real const DotProduct = Wj.dot(ReferenceNormal);


				QuadraturePoints_RaviartThomasBasis.setCoeff(n, El, j, 0) = Wj(0);
				QuadraturePoints_RaviartThomasBasis.setCoeff(n, El, j, 1) = Wj(1);

				for (unsigned m = 0; m < 3; m++) {


					Real const Phim			 = BasisPolynomial(m);
					Real const Phim_opposite = BasisPolynomialOpposite(m);
					Real const Value		 = w * DotProduct * Phim;

					QuadraturePoints_PolynomialBasis								.setCoeff(n, El, m) = Phim;
					QuadraturePoints_PolynomialBasisOpposite						.setCoeff(n, El, m) = Phim_opposite;
					QuadraturePoints_RaviartThomasBasisDotNormalTimesPolynomialBasis.setCoeff(n, j, El, m) = abs(Value) < INTEGRAL_PRECISION ? 0.0 : Value;

				}
			}
		}
	}



	/*****************************************************************************/
	/*                                                                           */
	/*    - Evaluate integrals of Source * P1 basis function phi_m			     */
	/*                                                                           */
	/*****************************************************************************/
	evaluate_source_integrals();
	

	/*****************************************************************************/
	/*                                                                           */
	/*    - Initilize constant matrices										     */
	/*    - Assembly of the trace pressure system and compute					 */
	/*		its LU decomposition												 */
	/*                                                                           */
	/*    - Get sparsity pattern of the Pressure system	(to optimize Eigen)		 */
	/*			: This will compute permutation only once,                       */
	/*            then we can just keep calling factorize(A)                     */
	/*            and the pattern will be already known		                     */
	/*                                                                           */
	/*****************************************************************************/
	std::cout << std::endl;
	std::cout << "/*****************/" << std::endl;
	std::cout << "/*               */" << std::endl;
	std::cout << "/*  Assembling   */" << std::endl;
	std::cout << "/*  -R           */" << std::endl;
	std::cout << "/*               */" << std::endl;
	std::cout << "/*****************/" << std::endl;
	std::cout << std::endl;

	assembleR();

	std::cout << std::endl;
	std::cout << "/*****************/" << std::endl;
	std::cout << "/*               */" << std::endl;
	std::cout << "/*  Assembling & */" << std::endl;
	std::cout << "/*  Factorizing  */" << std::endl;
	std::cout << "/*  -M           */" << std::endl;
	std::cout << "/*               */" << std::endl;
	std::cout << "/*****************/" << std::endl;
	std::cout << std::endl;

	assembleM();

	std::cout << std::endl;
	std::cout << "/*****************/" << std::endl;
	std::cout << "/*               */" << std::endl;
	std::cout << "/*  Sparsity     */" << std::endl;
	std::cout << "/*  Pattern      */" << std::endl;
	std::cout << "/*               */" << std::endl;
	std::cout << "/*****************/" << std::endl;
	std::cout << std::endl;

	getSparsityPatternOfThePressureSystem();


};
template<unsigned QuadraturePrecision, scheme TimeScheme>
solver<QuadraturePrecision, TimeScheme>::~solver() {

	delete[] Viscosities;
	delete[] Porosities;

	delete[] Thetas;
	delete[] Thetas_prev;

	delete[] MeshElements;
	delete[] MeshEdges;
	delete[] MeshElementIndeces;
	delete[] MeshEdgesIndeces;

	delete[] AffineMappingMatrixDeterminant;
	delete[] PorosityViscosityDeterminant;

};




template<unsigned QuadraturePrecision, scheme TimeScheme>
void solver<QuadraturePrecision, TimeScheme>::initializeValues() {


	Real const time = nt * dt;


	for (unsigned k = 0; k < nk; k++) {


		//tm_pointer const K = Mesh->get_triangle(k);
		tm_pointer const K		 = MeshElements[k];
		unsigned const	 k_index = K->index;
		Real const		 area	 = K->area();


		vm_pointer const va = K->vertices[0];
		vm_pointer const vb = K->vertices[1];
		vm_pointer const vc = K->vertices[2];

		Real const x0 = va->x;
		Real const y0 = va->y;

		Real const x1 = vb->x;
		Real const y1 = vb->y;

		Real const x2 = vc->x;
		Real const y2 = vc->y;


		/*****************************************************************************/
		/*                                                                           */
		/*    - Interpolant of the Barenblatt solution for the initial condition     */
		/*                                                                           */
		/*****************************************************************************/
		Real const B0 = barenblatt(x0, y0, time);
		Real const B1 = barenblatt(x1, y1, time);
		Real const B2 = barenblatt(x2, y2, time);

		Eigen::Vector3d const B(B0, B1, B2);

		Eigen::Matrix3d iM;
		iM << 0.0, 0.5, 0.5,
			 -0.5, 0.5, 0.0,
			 -0.5, 0.0, 0.5;
		Eigen::Vector3d const Solution = iM * B;

		Xi.setCoeff(k_index, 0) = Solution(0);
		Xi.setCoeff(k_index, 1) = Solution(1);
		Xi.setCoeff(k_index, 2) = Solution(2);

		Pi.setCoeff(k_index, 0) = equationOfState(Xi(k_index, 0));
		Pi.setCoeff(k_index, 1) = equationOfState(Xi(k_index, 1));
		Pi.setCoeff(k_index, 2) = equationOfState(Xi(k_index, 2));

		//Pi.setCoeff(k_index, 0) = integrate_triangle<Real>(K, time, barenblatt) / area;
		//Pi.setCoeff(k_index, 1) = 0.0;
		//Pi.setCoeff(k_index, 2) = 0.0;

		//Xi.setCoeff(k_index, 0) = integrate_triangle<Real>(K, time, barenblatt) / area;
		//Xi.setCoeff(k_index, 1) = 0.0;
		//Xi.setCoeff(k_index, 2) = 0.0;



		/*****************************************************************************/
		/*                                                                           */
		/*    - Mean values of the viscosity, porosity on each element			     */
		/*                                                                           */
		/*****************************************************************************/
		Viscosities[k_index] = integrateTriangle<QuadraturePrecision>(K, viscosity) / area;
		Porosities[k_index]  = integrateTriangle<QuadraturePrecision>(K, porosity) / area;


		/*****************************************************************************/
		/*                                                                           */
		/*    - Auxilary variables												     */
		/*                                                                           */
		/*****************************************************************************/
		rkFc.setCoeff(k_index, 0) = 0.0;
		rkFc.setCoeff(k_index, 1) = 0.0;
		rkFc.setCoeff(k_index, 2) = 0.0;

		rkFc_n.setCoeff(k_index, 0) = rkFc(k_index, 0);
		rkFc_n.setCoeff(k_index, 1) = rkFc(k_index, 1);
		rkFc_n.setCoeff(k_index, 2) = rkFc(k_index, 2);

		rkFp_n.setCoeff(k_index, 0) = 0.0;
		rkFp_n.setCoeff(k_index, 1) = 0.0;
		rkFp_n.setCoeff(k_index, 2) = 0.0;

	}

};
template<unsigned QuadraturePrecision, scheme TimeScheme>
void solver<QuadraturePrecision, TimeScheme>::computeThetas() {

	for (unsigned k = 0; k < nk; k++)
		Thetas[k] = 1.0;

};

/*
template<typename Real, unsigned QuadraturePrecision, scheme TimeScheme>
bool solver<Real, QuadraturePrecision, TimeScheme>::stopCriterion() {


	double sP1 = 0.0;
	double sP2 = 0.0;

	double sC1 = 0.0;
	double sC2 = 0.0;

	double sB1 = 0.0;
	double sB2 = 0.0;

	t_pointer K = NULL;


	for (unsigned k = 0; k < nk; k++) {

		//unsigned const k_index = mesh->get_triangle(k)->index;
		unsigned const k_index = MeshElementIndeces[k];

		sP1 += square(Pi(k_index, 0) - Pi_prev(k_index, 0));
		sP2 += square(Pi(k_index, 0));

	}

	double const val_P = sP1 / sP2;

	if (val_P > TOL)
		return false;


	for (unsigned k = 0; k < nk; k++) {

		//unsigned const k_index = mesh->get_triangle(k)->index;
		unsigned const k_index = MeshElementIndeces[k];

		sC1 += square(Xi(k_index, 0) - Xi_prev(k_index, 0));
		sC2 += square(Xi(k_index, 0));

	}

	double const val_C = sC1 / sC2;

	if (val_C > TOL)
		return false;


	for (unsigned k = 0; k < nk; k++) {

		sB1 += square(Thetas[k] - Thetas_prev[k]);
		sB2 += square(Thetas[k]);

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
	//		val = ((Pi(K, 0, 0)*phi1(y, a, b) + Pi(K, 0, 1)*phi2(y, a, b)) - (Pi_prev(K, 0, 0)*phi1(y, a, b) + Pi_prev(K, 0, 1) * phi2(y, a, b)));
	//		sP1 += w * square(val);
	//		sP2 += w * square(Pi(K, 0, 0)*phi1(y, a, b) + Pi(K, 0, 1)*phi2(y, a, b));
	//	}
	//}
	//double const val_P = sP1 / sP2;
	//if (val_P > TOL)
	//	return false;
	//return true;


};
*/
/*

template<unsigned QuadraturePrecision>
bool solver<QuadraturePrecision>::stopCriterion() {


	real const time = nt * dt;

	quadrature_triangle const QuadratureOnTriangle(quadrature_order);
	unsigned const			  NumberOfQuadraturePoints = QuadratureOnTriangle.NumberOfPoints;

	Eigen::VectorXd BasisPolynomial(3);

	real ErrorPressure = 0.0;
	real NormPressure = 0.0;

	real ErrorConcentration = 0.0;
	real NormConcentration = 0.0;

	real ErrorBeta = 0.0;
	real NormBeta = 0.0;

	for (unsigned k = 0; k < nk; k++) {


		unsigned const	k_index = MeshElementIndeces[k];
		real const		detJF = AffineMappingMatrixDeterminant[k_index];


		real IntegralErrorPressure	= 0.0;
		real IntegralNormPressure	= 0.0;

		real IntegralErrorConcentration = 0.0;
		real IntegralNormConcentration	= 0.0;

		for (unsigned n = 0; n < NumberOfQuadraturePoints; n++) {


			real const s = (real)QuadratureOnTriangle.points_x[n];
			real const t = (real)QuadratureOnTriangle.points_y[n];
			real const w = (real) 0.5 * QuadratureOnTriangle.weights[n];

			evaluate_polynomial_basis(s, t, BasisPolynomial);


			real DifferencePressure = 0.0;
			real NormPressure		= 0.0;

			real DifferenceConcentration	= 0.0;
			real NormConcentration			= 0.0;

			for (unsigned j = 0; j < 3; j++) {

				DifferencePressure	+= (Pi(k_index, j) - Pi_prev(k_index, j)) * BasisPolynomial(j);
				NormPressure		+= Pi(k_index, j) * BasisPolynomial(j);

				DifferenceConcentration += (Xi(k_index, j) - Xi_prev(k_index, j)) * BasisPolynomial(j);
				NormConcentration		+= Xi(k_index, j) * BasisPolynomial(j);

			}

			IntegralErrorPressure	+= w * square(DifferencePressure);
			IntegralNormPressure	+= w * square(NormPressure);

			IntegralErrorConcentration	+= w * square(DifferenceConcentration);
			IntegralNormConcentration	+= w * square(NormConcentration);

		}

		ErrorPressure	+= detJF * IntegralErrorPressure;
		NormPressure	+= detJF * IntegralNormPressure;

		ErrorConcentration	+= detJF * IntegralErrorConcentration;
		NormConcentration	+= detJF * IntegralNormConcentration;

	}

	real const eP = ErrorPressure / NormPressure;
	real const eC = ErrorConcentration / NormConcentration;

	if (eP > TOL || eC > TOL)
		return false;

	return true;


};
*/

template<unsigned QuadraturePrecision, scheme TimeScheme>
bool solver<QuadraturePrecision, TimeScheme>::stopCriterion() {


	Eigen::Vector3d BasisPolynomial(3);

	Real ErrorPressure = 0.0;
	Real NormPressure = 0.0;

	Real ErrorConcentration = 0.0;
	Real NormConcentration = 0.0;

	Real ErrorThetas = 0.0;
	Real NormThetas = 0.0;

	for (unsigned k = 0; k < nk; k++) {


		unsigned const k_index	 = MeshElementIndeces[k];
		Real const	   HalfDetJF = 0.5 * AffineMappingMatrixDeterminant[k_index];


		Real IntegralError = 0.0;
		Real IntegralNorm  = 0.0;


		for (unsigned n = 0; n < NumberOfQuadraturePointsTriangle; n++) {


			//Real const s = QuadraturePointsAndWeightsOnReferenceTriangle(n, 0);
			//Real const t = QuadraturePointsAndWeightsOnReferenceTriangle(n, 1);
			Real const w = QuadraturePointsAndWeightsOnReferenceTriangle(n, 2);

			//evaluate_polynomial_basis(s, t, BasisPolynomial);

			Real Difference = 0.0;
			Real Norm = 0.0;

			for (unsigned j = 0; j < 3; j++) {

				Real const Pj = QuadraturePoints_PolynomialBasisOnReferenceTriangle(n, j);

				Difference	+= Pj * (Pi(k_index, j) - Pi_prev(k_index, j));
				Norm		+= Pj * Pi(k_index, j);
			
				//Difference  += BasisPolynomial(j) * (Pi(k_index, j) - Pi_prev(k_index, j));
				//Norm		+= BasisPolynomial(j) * Pi(k_index, j);
			
			}
			
			IntegralError	+= w * square(Difference);
			IntegralNorm	+= w * square(Norm);

		}

		ErrorPressure	+= HalfDetJF * IntegralError;
		NormPressure	+= HalfDetJF * IntegralNorm;

	}

	Real const eP = ErrorPressure / NormPressure;

	if (eP > TOLERANCE)
		return false;

	/*
	for (unsigned k = 0; k < nk; k++) {


		unsigned const	k_index = MeshElementIndeces[k];
		Real const		detJF = AffineMappingMatrixDeterminant[k_index];


		Real IntegralError = 0.0;
		Real IntegralNorm = 0.0;

		for (unsigned n = 0; n < NumberOfQuadraturePoints; n++) {


			Real const s = QuadratureOnTriangle.points_x[n];
			Real const t = QuadratureOnTriangle.points_y[n];
			Real const w = 0.5 * QuadratureOnTriangle.weights[n];

			evaluate_polynomial_basis(s, t, BasisPolynomial);


			Real Difference = 0.0;
			Real Norm = 0.0;

			for (unsigned j = 0; j < 3; j++) {

				Difference += (Xi(k_index, j) - Xi_prev(k_index, j)) * BasisPolynomial(j);
				Norm += Xi(k_index, j) * BasisPolynomial(j);

			}

			IntegralError += w * square(Difference);
			IntegralNorm += w * square(Norm);

		}

		ErrorConcentration += detJF * IntegralError;
		NormConcentration += detJF * IntegralNorm;

	}

	Real const eC = ErrorConcentration / NormConcentration;

	if (eC > TOLERANCE)
		return false;
		*/

		/*for (unsigned k = 0; k < nk; k++) {

			sB1 += square(betas[k] - betas_prev[k]);
			sB2 += square(betas[k]);

		}

		double const val_B = sB1 / sB2;

		if (val_B > TOLERANCE)
			return false;*/

	return true;


};

template<unsigned QuadraturePrecision, scheme TimeScheme>
void solver<QuadraturePrecision, TimeScheme>::concentrationCorrection() {

	//unsigned const nk = nk;

	//double const eps = sqrt(DBL_EPSILON);

	//double s1 = 0.0;
	//double s2 = 0.0;

	//element * K = NULL;

	//for (unsigned k = 0; k < nk; k++) {

	//	K = mesh->getElement(k);

	//	s1 += square(Xi(K, 0, 0) - Xi_n(K, 0, 0));
	//	s2 += square(Xi(K, 0, 0));

	//}

	//if (sqrt(s1) / sqrt(s2) < eps) {

	//	for (unsigned k = 0; k < nk; k++) {

	//		std::cout << "Corrected" << std::endl;
	//		K = mesh->getElement(k);

	//		Xi.setCoeff(K, 0, 0) = Xi_n(K, 0, 0) + DBL_EPSILON * Xi(K, 0, 0);

	//	}

	//}

};



template<unsigned QuadraturePrecision, scheme TimeScheme>
void solver<QuadraturePrecision, TimeScheme>::computeTracePressures() {



	/*****************************************************************************/
	/*                                                                           */
	/*    - Copy Internal Pressures into Eigen container						 */
	/*                                                                           */
	/*****************************************************************************/
	for (unsigned k = 0; k < nk; k++) {
	
		unsigned const k_index = MeshElementIndeces[k];
	
		for (unsigned m = 0; m < 3; m++)
			PiTemp[3 * k_index + m] = Pi(k_index, m);
	
	}

	assembleV();


	traceSystemRhs.head(ne) = R1 * PiTemp - V1;
	traceSystemRhs.tail(ne) = R2 * PiTemp - V2;

	Tp = sparseLUsolver_TracePressureSystem.solve(traceSystemRhs);


	/*****************************************************************************/
	/*                                                                           */
	/*    - Copy Trace Pressure solution to each elements's edges from			 */
	/*      the Eigen container                                                  */
	/*                                                                           */
	/*****************************************************************************/
	for (unsigned e = 0; e < ne; e++) {


		em_pointer const E = MeshEdges[e];

		Real const TpValue1 = Tp[E->index];
		Real const TpValue2 = Tp[E->index + ne];


		for (unsigned neighbor = 0; neighbor < 2; neighbor++) {

			tm_pointer const K = E->neighbors[neighbor];

			if (!K)
				continue;

			unsigned const k_index		 = K->index;
			unsigned const e_index_local = K->get_edge_index(E);

			TPi.setCoeff(k_index, e_index_local, 0) = TpValue1;
			TPi.setCoeff(k_index, e_index_local, 1) = TpValue2;

		}
	}

};
template<unsigned QuadraturePrecision, scheme TimeScheme>
void solver<QuadraturePrecision, TimeScheme>::computePressureEquation() {



	assembleG();
	assembleV();

	assemblePressureSystem();

	pressureSystemRhs.head(ne) = R1iD * G - V1;
	pressureSystemRhs.tail(ne) = R2iD * G - V2;


	/*****************************************************************************/
	/*                                                                           */
	/*    - Sparsity pattern is already known. Method factorize() is enough		 */
	/*                                                                           */
	/*****************************************************************************/
	//sparseLUsolver_PressureSystem.factorize(PressureSystem);
	BiConjugateGradientSolver.factorize(PressureSystem);

	Tp		= BiConjugateGradientSolver.solveWithGuess(pressureSystemRhs, Tp);
	PiTemp	= iD * G - (iDH1 * Tp.head(ne) + iDH2 * Tp.tail(ne));

	
	/*****************************************************************************/
	/*                                                                           */
	/*    - Copy Internal Pressures from Eigen container						 */
	/*                                                                           */
	/*****************************************************************************/
	for (unsigned k = 0; k < nk; k++) {

		unsigned const k_index = MeshElementIndeces[k];

		for (unsigned m = 0; m < 3; m++)
			Pi.setCoeff(k_index, m) = PiTemp[3 * k_index + m];

	}

	/*****************************************************************************/
	/*                                                                           */
	/*    - Copy Trace Pressure solution to each elements's edges from			 */
	/*      the Eigen container                                                  */
	/*                                                                           */
	/*****************************************************************************/
	for (unsigned e = 0; e < ne; e++) {


		em_pointer const E = MeshEdges[e];

		Real const TpValue1 = Tp[E->index];
		Real const TpValue2 = Tp[E->index + ne];


		for (unsigned neighbor = 0; neighbor < 2; neighbor++) {

			tm_pointer const K = E->neighbors[neighbor];

			if (!K)
				continue;

			unsigned const k_index		 = K->index;
			unsigned const e_index_local = K->get_edge_index(E);

			TPi.setCoeff(k_index, e_index_local, 0) = TpValue1;
			TPi.setCoeff(k_index, e_index_local, 1) = TpValue2;

		}
	}

};



template<unsigned QuadraturePrecision, scheme TimeScheme>
void solver<QuadraturePrecision, TimeScheme>::computeVelocities() {


	Real const time = (nt + 1) * dt;

//	omp_set_num_threads(2);
//#pragma omp parallel for 
	for (int k = 0; k < nk; k++) {


		tm_pointer const K		 = MeshElements[k];
		unsigned const	 k_index = MeshElementIndeces[k];


		/*****************************************************************************/
		/*                                                                           */
		/*    - Loop over edges of the element										 */
		/*                                                                           */
		/*****************************************************************************/
		for (unsigned El = 0; El < 3; El++) {


			em_pointer const E = K->edges[El];


			if (E->marker == E_MARKER::NEUMANN) {

				Velocity.setCoeff(k_index, LI(K, E, 0)) = NEUMANN_GAMMA_Q_velocity(E, time);
				Velocity.setCoeff(k_index, LI(K, E, 1)) = 0.0;

				continue;

			}

			/*****************************************************************************/
			/*                                                                           */
			/*    - Loop over the two degrees of freedom on each edge					 */
			/*                                                                           */
			/*****************************************************************************/
			for (unsigned dof = 0; dof < 2; dof++) {


				unsigned const j = LI(K, E, dof);

				Real Value1 = 0.0;
				Real Value2 = 0.0;

				for (unsigned m = 0; m < 8; m++) {


					Real const AlphaTemp = Alpha(k_index, m, j);

					Real BetaPi		= 0.0;
					Real ChiTracePi = 0.0;

					for (unsigned i = 0; i < 3; i++)
						BetaPi += Beta(m, i) * Pi(k_index, i);
										
					Value1 += AlphaTemp * BetaPi;


					for (unsigned l = 0; l < 3; l++)
						for (unsigned s = 0; s < 2; s++)
							ChiTracePi += Chi(k_index, m, l, s) * TPi(k_index, l, s);

					Value2 += AlphaTemp * ChiTracePi;

				}

				Velocity.setCoeff(k_index, j) = (Value1 - Value2) / Viscosities[k_index];

			}
		}


		/*****************************************************************************/
		/*                                                                           */
		/*    - Compute Bubble velocities inside element							 */
		/*                                                                           */
		/*****************************************************************************/
		for (unsigned dof = 6; dof < 8; dof++) {


			Real Value1 = 0.0;
			Real Value2 = 0.0;

			for (unsigned m = 0; m < 8; m++) {


				Real const AlphaTemp = Alpha(k_index, dof, m);

				Real BetaPi		= 0.0;
				Real ChiTracePi = 0.0;


				for (unsigned i = 0; i < 3; i++)
					BetaPi += Beta(m, i) * Pi(k_index, i);
				
				Value1 += AlphaTemp * BetaPi;


				for (unsigned l = 0; l < 3; l++)
					for (unsigned s = 0; s < 2; s++)
						ChiTracePi += Chi(k_index, m, l, s) * TPi(k_index, l, s);

				Value2 += AlphaTemp * ChiTracePi;

			}

			Velocity.setCoeff(k_index, dof) = (Value1 - Value2) / Viscosities[k_index];

		}

	}

};
template<unsigned QuadraturePrecision, scheme TimeScheme>
void solver<QuadraturePrecision, TimeScheme>::updateConcentrations() {


	Real const TimeCoefficient1 = dt * TimeSchemeParameter;
	Real const TimeCoefficient2 = dt * (1.0 - TimeSchemeParameter);

//	omp_set_num_threads(2);
//#pragma omp parallel for 
	for (int k = 0; k < nk; k++) {


		unsigned const	k_index		= MeshElementIndeces[k];
		Real const		Coefficient = 1.0 / (Porosities[k_index] * AffineMappingMatrixDeterminant[k_index]);

		Real Value0 = 0.0;
		Real Value1 = 0.0;
		Real Value2 = 0.0;

		for (unsigned j = 0; j < 8; j++) {


			Real EtaGamma0 = 0.0;
			Real EtaGamma1 = 0.0;
			Real EtaGamma2 = 0.0;

			for (unsigned l = 0; l < 3; l++) {

				EtaGamma0 += Eta(0, l) * Gamma(k_index, l, j);
				EtaGamma1 += Eta(1, l) * Gamma(k_index, l, j);
				EtaGamma2 += Eta(2, l) * Gamma(k_index, l, j);

			}

			Value0 += Velocity(k_index, j) * EtaGamma0;
			Value1 += Velocity(k_index, j) * EtaGamma1;
			Value2 += Velocity(k_index, j) * EtaGamma2;

		}

		// In the computation of Eta, there is coefficient detJF. When inverting, the coefficient is inverted
		rkFc.setCoeff(k_index, 0) = -Value0 * Coefficient;
		rkFc.setCoeff(k_index, 1) = -Value1 * Coefficient;
		rkFc.setCoeff(k_index, 2) = -Value2 * Coefficient;

		Xi.setCoeff(k_index, 0) = Xi_n(k_index, 0) + (TimeCoefficient1 * rkFc(k_index, 0) + TimeCoefficient2 * rkFc_n(k_index, 0));
		Xi.setCoeff(k_index, 1) = Xi_n(k_index, 1) + (TimeCoefficient1 * rkFc(k_index, 1) + TimeCoefficient2 * rkFc_n(k_index, 1));
		Xi.setCoeff(k_index, 2) = Xi_n(k_index, 2) + (TimeCoefficient1 * rkFc(k_index, 2) + TimeCoefficient2 * rkFc_n(k_index, 2));

	}

};



template<unsigned QuadraturePrecision, scheme TimeScheme>
Real solver<QuadraturePrecision, TimeScheme>::upwindConcentration(tm_pointer const & K, unsigned const El, unsigned const n) {


	Real const time = (nt + 1) * dt;

	unsigned const	 k_index	= K->index;
	em_pointer const E			= K->edges[El];
	E_MARKER const	 e_marker	= E->marker;


	Real VelocityDotNormal = 0.0;
	Real Concentration = 0.0;

	if (e_marker == E_MARKER::NEUMANN) {


		VelocityDotNormal = NEUMANN_GAMMA_Q_velocity(E, time);

		//Real const X = QuadraturePoints_Edge_x(k_index, El, n);
		//Real const Y = QuadraturePoints_Edge_y(k_index, El, n);
		//
		//VelocityDotNormal = NEUMANN_GAMMA_Q_velocity(X, Y, time);

		if (VelocityDotNormal < 0.0) {

			Real const X = QuadraturePoints_Edge_x(k_index, El, n);
			Real const Y = QuadraturePoints_Edge_y(k_index, El, n);

			return DIRICHLET_GAMMA_Q_concentration(X, Y, time);

		}
	}
	else {

		for (unsigned dof = 0; dof < 2; dof++) {

			unsigned const j = LI(K, E, dof);

			VelocityDotNormal += Velocity(k_index, j) * QuadraturePoints_PhysicalNormalDotPhysicalRaviartThomasBasis(k_index, El, n, j);

		}
	}


	if (VelocityDotNormal >= 0.0) {

		for (unsigned m = 0; m < 3; m++)
			Concentration += Xi_prev(k_index, m) * QuadraturePoints_PolynomialBasis(n, El, m);

	}
	else {


		if (E->marker == E_MARKER::DIRICHLET) {

			Real const X = QuadraturePoints_Edge_x(k_index, El, n);
			Real const Y = QuadraturePoints_Edge_y(k_index, El, n);

			return DIRICHLET_GAMMA_P_concentration(X, Y, time);

		}

		unsigned const kn_index			= K->neighbors[El]->index;
		unsigned const e_index_Kn_loc	= K->neighbors[El]->get_edge_index(E);

		// The linear combination of the polynomial basis must be computed backwards
		for (unsigned m = 0; m < 3; m++)
			Concentration += Xi_prev(kn_index, m) * QuadraturePoints_PolynomialBasisOpposite(n, e_index_Kn_loc, m);

	}

	return Concentration;

};
template<unsigned QuadraturePrecision, scheme TimeScheme>
Real solver<QuadraturePrecision, TimeScheme>::upwindConcentration_limiter(tm_pointer const & K, unsigned const El, unsigned const n) {


	Real const time = (nt + 1) * dt;

	unsigned const		k_index = K->index;
	em_pointer const	E = K->edges[El];
	E_MARKER const		e_marker = E->marker;


	Real VelocityDotNormal = 0.0;
	Real Concentration = 0.0;

	if (e_marker == E_MARKER::NEUMANN) {


		VelocityDotNormal = NEUMANN_GAMMA_Q_velocity(E, time);

		if (VelocityDotNormal < 0.0) {

			Real const X = QuadraturePoints_Edge_x(k_index, El, n);
			Real const Y = QuadraturePoints_Edge_y(k_index, El, n);

			return DIRICHLET_GAMMA_Q_concentration(X, Y, time);

		}
	}
	else {

		for (unsigned dof = 0; dof < 2; dof++) {

			unsigned const j = LI(K, E, dof);

			VelocityDotNormal += Velocity(k_index, j) * QuadraturePoints_PhysicalNormalDotPhysicalRaviartThomasBasis(k_index, El, n, j);

		}
	}

	if (VelocityDotNormal >= 0.0) {

		Real const CKmean = Xi_prev(k_index, 0) * QuadraturePoints_PolynomialBasis(n, El, 0);

		for (unsigned m = 0; m < 3; m++)
			Concentration += Xi_prev(k_index, m) * QuadraturePoints_PolynomialBasis(n, El, m);

		if (K->neighbors[El]) {

			unsigned const kn_index = K->neighbors[El]->index;
			unsigned const e_index_Kn_loc = K->neighbors[El]->get_edge_index(E);

			Real const CKNmean = Xi_prev(kn_index, 0) * QuadraturePoints_PolynomialBasis(n, e_index_Kn_loc, 0);


			Real const Cmin = std::min(CKmean, CKNmean);
			Real const Cmax = std::max(CKmean, CKNmean);


			if (Concentration < Cmin) return Cmin;
			else if (Concentration > Cmax) return Cmax;

		}

		return Concentration;

	}
	else {


		if (E->marker == E_MARKER::DIRICHLET) {

			Real const X = QuadraturePoints_Edge_x(k_index, El, n);
			Real const Y = QuadraturePoints_Edge_y(k_index, El, n);

			return DIRICHLET_GAMMA_P_concentration(X, Y, time);

		}

		unsigned const kn_index = K->neighbors[El]->index;
		unsigned const e_index_Kn_loc = K->neighbors[El]->get_edge_index(E);


		Real const CKmean = Xi_prev(k_index, 0) * QuadraturePoints_PolynomialBasis(n, El, 0);
		Real const CKNmean = Xi_prev(kn_index, 0) * QuadraturePoints_PolynomialBasis(n, e_index_Kn_loc, 0);

		Real const Cmin = std::min(CKmean, CKNmean);
		Real const Cmax = std::max(CKmean, CKNmean);


		Concentration = 0.0;

		for (unsigned m = 0; m < 3; m++)
			Concentration += Xi_prev(kn_index, m) * QuadraturePoints_PolynomialBasis(n, e_index_Kn_loc, m);

		if (Concentration < Cmin) return Cmin;
		else if (Concentration > Cmax) return Cmax;

		return Concentration;

	}

	/*
	if (VelocityDotNormal >= 0.0) {

		for (unsigned m = 0; m < 3; m++)
			Concentration += Xi_prev(k_index, m) * QuadraturePoints_PolynomialBasis(n, El, m);

		//Concentration = Xi_prev(k_index, 0) * QuadraturePoints_PolynomialBasis(n, El, 0);

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
			Concentration += Xi_prev(kn_index, m) * QuadraturePoints_PolynomialBasis(n, e_index_Kn_loc, m);

		//Concentration = Xi_prev(kn_index, 0) * QuadraturePoints_PolynomialBasis(n, e_index_Kn_loc, 0);

	}
	*/


	//return Concentration;

};



template<unsigned QuadraturePrecision, scheme TimeScheme>
void solver<QuadraturePrecision, TimeScheme>::assembleR() {



	for (unsigned e = 0; e < ne; e++) {


		//em_pointer const E			= Mesh->get_edge(e);
		em_pointer const E			= MeshEdges[e];
		unsigned const	 e_index	= E->index;
		E_MARKER const	 e_marker	= E->marker;


		if (e_marker == E_MARKER::DIRICHLET)
			continue;


		for (unsigned neighborElement = 0; neighborElement < 2; neighborElement++) {

			tm_pointer const K = E->neighbors[neighborElement];

			if (!K)
				continue;


			unsigned const k_index = K->index;

			unsigned const dof0 = LI(K, E, 0);
			unsigned const dof1 = LI(K, E, 1);

			Real const ChiCoeff0 = Chi(k_index, dof0, K->get_edge_index(E), 0);
			Real const ChiCoeff1 = Chi(k_index, dof1, K->get_edge_index(E), 1);


			for (unsigned l = 0; l < 3; l++) {

				Real Value1 = 0.0;
				Real Value2 = 0.0;

				for (unsigned j = 0; j < 8; j++) {

					Value1 += Alpha(k_index, dof0, j) * Beta(j, l);
					Value2 += Alpha(k_index, dof1, j) * Beta(j, l);

				}

				Value1 *= ChiCoeff0 / Viscosities[k_index];
				Value2 *= ChiCoeff1 / Viscosities[k_index];

				R1.coeffRef(e_index, 3 * k_index + l) = abs(Value1) < INTEGRAL_PRECISION ? 0.0 : Value1;
				R2.coeffRef(e_index, 3 * k_index + l) = abs(Value2) < INTEGRAL_PRECISION ? 0.0 : Value2;

			}
		}

		/*
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

						AB1 += Alpha(k_index, dof0, i)*Beta(i, m);
						AB2 += Alpha(k_index, dof1, i)*Beta(i, m);

					}

					real const val1 = AB1 / viscosities[k_index];
					real const val2 = AB2 / viscosities[k_index];

					R1.coeffRef(e_index, 3 * k_index + m) = abs(val1) < INTEGRAL_PRECISION ? 0.0 : val1;
					R2.coeffRef(e_index, 3 * k_index + m) = abs(val2) < INTEGRAL_PRECISION ? 0.0 : val2;

				}
			}
			*/

	}

};
template<unsigned QuadraturePrecision, scheme TimeScheme>
void solver<QuadraturePrecision, TimeScheme>::assembleM() {


	//unsigned const NumberOfDirichletEdges	= Mesh->get_number_of_dirichlet_edges();
	//unsigned const NumberOfNeumannEdges		= Mesh->get_number_of_neumann_edges();
	//unsigned const NumberOfBoundaryEdges	= NumberOfDirichletEdges + NumberOfNeumannEdges;

	//unsigned const NumberOfElements = 4 * (NumberOfDirichletEdges + (NumberOfBoundaryEdges - NumberOfDirichletEdges) * 3 + (ne - NumberOfBoundaryEdges) * 5 + ne - NumberOfBoundaryEdges);

	std::vector<Eigen::Triplet<Real>> triplet;
	//triplet.reserve(NumberOfElements);

	for (unsigned e = 0; e < ne; e++) {


		//em_pointer const E = Mesh->get_edge(e);
		em_pointer const E			= MeshEdges[e];
		unsigned const	 e_index	= E->index;
		E_MARKER const	 e_marker	= E->marker;



		if (e_marker == E_MARKER::DIRICHLET) {

			M_j1_s1.coeffRef(e_index, e_index) = -1.0;
			M_j1_s2.coeffRef(e_index, e_index) = +0.0;

			M_j2_s1.coeffRef(e_index, e_index) = +0.0;
			M_j2_s2.coeffRef(e_index, e_index) = -1.0;

			Eigen::Triplet<Real> const T_j1_s1(e_index,		 e_index,		-1.0);
			Eigen::Triplet<Real> const T_j1_s2(e_index,		 e_index + ne,	+0.0);
			Eigen::Triplet<Real> const T_j2_s1(e_index + ne, e_index,		+0.0);
			Eigen::Triplet<Real> const T_j2_s2(e_index + ne, e_index + ne,  -1.0);

			triplet.push_back(T_j1_s1);
			triplet.push_back(T_j1_s2);
			triplet.push_back(T_j2_s1);
			triplet.push_back(T_j2_s2);

			continue;

		}


		for (unsigned neighborElement = 0; neighborElement < 2; neighborElement++) {


			tm_pointer const K = E->neighbors[neighborElement];

			if (!K)
				continue;


			unsigned const k_index = K->index;

			unsigned const dof0 = LI(K, E, 0);
			unsigned const dof1 = LI(K, E, 1);

			Real const ChiCoeff0 = Chi(k_index, dof0, K->get_edge_index(E), 0);
			Real const ChiCoeff1 = Chi(k_index, dof1, K->get_edge_index(E), 1);


			for (unsigned El = 0; El < 3; El++) {


				em_pointer const E_local			  = K->edges[El];
				unsigned const	 e_local_index_global = E_local->index;

				Real ACHI_j1_s1 = 0.0;
				Real ACHI_j1_s2 = 0.0;
				Real ACHI_j2_s1 = 0.0;
				Real ACHI_j2_s2 = 0.0;

				for (unsigned j = 0; j < 8; j++) {

					ACHI_j1_s1 += Alpha(k_index, dof0, j) * Chi(k_index, j, El, 0);
					ACHI_j1_s2 += Alpha(k_index, dof0, j) * Chi(k_index, j, El, 1);

					ACHI_j2_s1 += Alpha(k_index, dof1, j) * Chi(k_index, j, El, 0);
					ACHI_j2_s2 += Alpha(k_index, dof1, j) * Chi(k_index, j, El, 1);

				}

				Real const Value11 = ChiCoeff0 * ACHI_j1_s1 / Viscosities[k_index];
				Real const Value12 = ChiCoeff0 * ACHI_j1_s2 / Viscosities[k_index];
				Real const Value21 = ChiCoeff1 * ACHI_j2_s1 / Viscosities[k_index];
				Real const Value22 = ChiCoeff1 * ACHI_j2_s2 / Viscosities[k_index];


				M_j1_s1.coeffRef(e_index, e_local_index_global) += Value11;
				M_j1_s2.coeffRef(e_index, e_local_index_global) += Value12;
				M_j2_s1.coeffRef(e_index, e_local_index_global) += Value21;
				M_j2_s2.coeffRef(e_index, e_local_index_global) += Value22;


				Eigen::Triplet<Real> const T_j1_s1(e_index,		 e_local_index_global,		Value11);
				Eigen::Triplet<Real> const T_j1_s2(e_index,		 e_local_index_global + ne, Value12);
				Eigen::Triplet<Real> const T_j2_s1(e_index + ne, e_local_index_global,		Value21);
				Eigen::Triplet<Real> const T_j2_s2(e_index + ne, e_local_index_global + ne, Value22);

				triplet.push_back(T_j1_s1);
				triplet.push_back(T_j1_s2);
				triplet.push_back(T_j2_s1);
				triplet.push_back(T_j2_s2);

			}
		}

		/*
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

					ACHI_j1_s1 += Alpha(k_index, dof0, i) * Chi(i, e_local_index_local, 0);
					ACHI_j1_s2 += Alpha(k_index, dof0, i) * Chi(i, e_local_index_local, 1);

					ACHI_j2_s1 += Alpha(k_index, dof1, i) * Chi(i, e_local_index_local, 0);
					ACHI_j2_s2 += Alpha(k_index, dof1, i) * Chi(i, e_local_index_local, 1);

				}

				real const val1 = ACHI_j1_s1 / viscosities[k_index];
				real const val2 = ACHI_j1_s2 / viscosities[k_index];

				real const val3 = ACHI_j2_s1 / viscosities[k_index];
				real const val4 = ACHI_j2_s2 / viscosities[k_index];


				//if (e_local_index_global == e_index) {
				//	M_j1_s1.coeffRef(e_index, e_local_index_global) += val1;
				//	M_j1_s2.coeffRef(e_index, e_local_index_global) += val2;
				//	M_j2_s1.coeffRef(e_index, e_local_index_global) += val3;
				//	M_j2_s2.coeffRef(e_index, e_local_index_global) += val4;
				//	continue;
				//}
				//M_j1_s1.coeffRef(e_index, e_local_index_global) = val1;
				//M_j1_s2.coeffRef(e_index, e_local_index_global) = val2;
				//M_j2_s1.coeffRef(e_index, e_local_index_global) = val3;
				//M_j2_s2.coeffRef(e_index, e_local_index_global) = val4;


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
		*/

	}
	

	SparseMatrix tracePressureSystem_LU(2 * ne, 2 * ne);

	/*****************************************************************************/
	/*                                                                           */
	/*    - Assembly of the matrix [ M11,M12 ; M21,M22 ] using Eigen Triplets    */
	/*                                                                           */
	/*****************************************************************************/
	tracePressureSystem_LU.setFromTriplets(triplet.begin(), triplet.end());


	/*****************************************************************************/
	/*                                                                           */
	/*    - Compute the LU decomposition of the [ M11,M12 ; M21,M22 ]			 */
	/*                                                                           */
	/*****************************************************************************/
	sparseLUsolver_TracePressureSystem.compute(tracePressureSystem_LU);


};
template<unsigned QuadraturePrecision, scheme TimeScheme>
void solver<QuadraturePrecision, TimeScheme>::assembleV() {


	Real const time = (nt + 1) * dt;

//	omp_set_num_threads(2);
//#pragma omp parallel for 
	for (int e = 0; e < ne; e++) {


		//em_pointer const E = Mesh->get_edge(e);
		em_pointer const E			= MeshEdges[e];
		unsigned const	 e_index	= E->index;
		E_MARKER const	 e_marker	= E->marker;


		V1[e_index] = 0.0;
		V2[e_index] = 0.0;


		if (e_marker == E_MARKER::NEUMANN) {


			Real ChiCoeff0 = 1.0;
			Real ChiCoeff1 = 1.0;

			for (unsigned neighborElement = 0; neighborElement < 2; neighborElement++) {

				tm_pointer const K = E->neighbors[neighborElement];

				if (!K)
					continue;

				unsigned const k_index = K->index;

				unsigned const dof0 = LI(K, E, 0);
				unsigned const dof1 = LI(K, E, 1);

				ChiCoeff0 = Chi(k_index, dof0, K->get_edge_index(E), 0);
				ChiCoeff1 = Chi(k_index, dof1, K->get_edge_index(E), 1);

				break;

			}

			V1[e_index] = ChiCoeff0 * NEUMANN_GAMMA_Q_velocity(E, time);
			V2[e_index] = ChiCoeff1 * 0.0;

		}

		else if (e_marker == E_MARKER::DIRICHLET) {

			Real const x0 = E->a->x;
			Real const y0 = E->a->y;

			Real const x1 = E->b->x;
			Real const y1 = E->b->y;


			/*****************************************************************************/
			/*                                                                           */
			/*    - Interpolant of the Barenblatt solution for the Dirichlet Edge	     */
			/*    - Optimized : System matrix is known, so we can write the solution     */
			/*					immediately using inverse matrix (which is also known)   */
			/*                                                                           */
			/*                : M = [1, -1 ; 1, 1]    M^(-1) = 0.5*[1, 1; -1, 1]         */
			/*                                                                           */
			/*****************************************************************************/
			Real const B0 = DIRICHLET_GAMMA_P_pressure(x0, y0, time);
			Real const B1 = DIRICHLET_GAMMA_P_pressure(x1, y1, time);

			V1[e_index] = 0.5 * (+B0 + B1);
			V2[e_index] = 0.5 * (-B0 + B1);

			//Eigen::Vector2d const Solution = [] { Eigen::Matrix2d Temp; Temp << 0.5, 0.5, -0.5, 0.5; return Temp; }() * Eigen::Vector2d (B0, B1);
			//
			//Eigen::Vector2d const B(B0, B1);
			//Eigen::Matrix2d M;
			//
			//M << +1.0, -1.0,
			//	 +1.0, +1.0;
			//
			//Eigen::Vector2d const Solution = M.colPivHouseholderQr().solve(B);
			//
			//V1[e_index] = Solution(0);
			//V2[e_index] = Solution(1);

		}
	}

};




template<unsigned QuadraturePrecision, scheme TimeScheme>
void solver<QuadraturePrecision, TimeScheme>::assembleInverseD() {


	//assemble_Sigma();

	Real const TimeCoefficient = TimeSchemeParameter * dt;
	Eigen::Matrix3d Block;


	for (unsigned k = 0; k < nk; k++) {


		unsigned const	 k_index	 = MeshElementIndeces[k];
		unsigned const	 start_index = 3 * k_index;


		for (unsigned r = 0; r < 3; r++)
			for (unsigned s = 0; s < 3; s++)
				Block(r, s) = kroneckerDelta(r, s) - TimeCoefficient * Sigma(k_index, r, s);

		Eigen::Matrix3d const InverseBlock = Block.inverse();

		for (unsigned i = 0; i < 3; i++)
			for (unsigned j = 0; j < 3; j++)
				iD.coeffRef(start_index + i, start_index + j) = InverseBlock(i, j);





	}

};
template<unsigned QuadraturePrecision, scheme TimeScheme>
void solver<QuadraturePrecision, TimeScheme>::assembleH() {


	//assemble_Lambda();

	Eigen::Matrix3d Block1;
	Eigen::Matrix3d Block2;

	Real const TimeCoefficient = -TimeSchemeParameter * dt;


	for (unsigned k = 0; k < nk; k++) {


		tm_pointer const K		 = MeshElements[k];
		unsigned const	 k_index = MeshElementIndeces[k];

		em_pointer const E0 = K->edges[0];
		em_pointer const E1 = K->edges[1];
		em_pointer const E2 = K->edges[2];

		unsigned const e_index0 = E0->index;
		unsigned const e_index1 = E1->index;
		unsigned const e_index2 = E2->index;


		for (unsigned m = 0; m < 3; m++) {
			for (unsigned Ei = 0; Ei < 3; Ei++) {

				Block1(m, Ei) = TimeCoefficient * Lambda(k_index, 0, m, Ei);
				Block2(m, Ei) = TimeCoefficient * Lambda(k_index, 1, m, Ei);

			}
		}


		for (unsigned m = 0; m < 3; m++) {

			H1.coeffRef(3 * k_index + m, e_index0) = Block1.coeff(m, 0);
			H1.coeffRef(3 * k_index + m, e_index1) = Block1.coeff(m, 1);
			H1.coeffRef(3 * k_index + m, e_index2) = Block1.coeff(m, 2);

			H2.coeffRef(3 * k_index + m, e_index0) = Block2.coeff(m, 0);
			H2.coeffRef(3 * k_index + m, e_index1) = Block2.coeff(m, 1);
			H2.coeffRef(3 * k_index + m, e_index2) = Block2.coeff(m, 2);

		}
	}

};
template<unsigned QuadraturePrecision, scheme TimeScheme>
void solver<QuadraturePrecision, TimeScheme>::assembleG() {


	Real const TimeCoefficient = dt * (1.0 - TimeSchemeParameter);

//	omp_set_num_threads(2);
//#pragma omp parallel for 
	for (int k = 0; k < nk; k++) {

		unsigned const k_index = MeshElementIndeces[k];

		for (unsigned m = 0; m < 3; m++)
			G[3 * k_index + m] = Pi_n(k_index, m) + TimeCoefficient * rkFp_n(k_index, m);

	}

};



template<unsigned QuadraturePrecision, scheme TimeScheme>
void solver<QuadraturePrecision, TimeScheme>::assemblePressureSystem() {


	Real const TimeCoefficient = TimeSchemeParameter * dt;

	Eigen::Matrix3d Block;
	Eigen::Matrix3d Block1;
	Eigen::Matrix3d Block2;

	unsigned count = 0;

	/*****************************************************************************/
	/*                                                                           */
	/*    - Assemble of the diagonal. Therefore, there is no need 				 */
	/*		for += operator in the sequel									     */
	/*                                                                           */
	/*****************************************************************************/
	for (unsigned e = 0; e < ne; e++) {

		Real const M11 = M_j1_s1.coeff(e, e);
		Real const M12 = M_j1_s2.coeff(e, e);
		Real const M21 = M_j2_s1.coeff(e, e);
		Real const M22 = M_j2_s2.coeff(e, e);

		Eigen::Triplet<Real> const T11(e,		e,		M11);
		Eigen::Triplet<Real> const T12(e,		e + ne, M12);
		Eigen::Triplet<Real> const T21(e + ne,	e,		M21);
		Eigen::Triplet<Real> const T22(e + ne,	e + ne, M22);

		TripletVector[count++] = T11;
		TripletVector[count++] = T12;
		TripletVector[count++] = T21;
		TripletVector[count++] = T22;

		//PressureSystem.coeffRef(e,		e)		= M11;
		//PressureSystem.coeffRef(e,		e + ne) = M12;
		//PressureSystem.coeffRef(e + ne, e)		= M21;
		//PressureSystem.coeffRef(e + ne, e + ne) = M22;

	}

	//#pragma omp parallel for
	for (unsigned k = 0; k < nk; k++) {

		tm_pointer const K			 = MeshElements[k];
		unsigned const	 k_index	 = MeshElementIndeces[k];
		unsigned const	 start_index = 3 * k_index;

		/*
		Real * const R1BlockArray0 = R1_block.getRowPointer(k_index, 0, 0);
		Real * const R1BlockArray1 = R1_block.getRowPointer(k_index, 1, 0);
		Real * const R1BlockArray2 = R1_block.getRowPointer(k_index, 2, 0);
		*/

		/*****************************************************************************/
		/*                                                                           */
		/*    - Inverse of the matrix D										         */
		/*                                                                           */
		/*****************************************************************************/
		for (unsigned r = 0; r < 3; r++)
			for (unsigned s = 0; s < 3; s++)
				Block(r, s) = kroneckerDelta(r, s) - TimeCoefficient * Sigma(k_index, r, s);
			
		Eigen::Matrix3d const InverseBlock = Block.inverse();

		for (unsigned i = 0; i < 3; i++)
			for (unsigned j = 0; j < 3; j++)
				iD.coeffRef(start_index + i, start_index + j) = InverseBlock(i, j);


		/*****************************************************************************/
		/*                                                                           */
		/*    - H1, H2 blocks multiplied by the blocks Inverse of D					 */
		/*                                                                           */
		/*****************************************************************************/
		Eigen::Matrix3d iDH1Block;
		Eigen::Matrix3d iDH2Block;

		for (unsigned m = 0; m < 3; m++) {


			Real const InverseBlockArray[3] = { InverseBlock.coeff(m, 0), InverseBlock.coeff(m, 1) , InverseBlock.coeff(m, 2)};

			for (unsigned Ei = 0; Ei < 3; Ei++) {

				Real Sum0 = 0.0;
				Real Sum1 = 0.0;

				Real const LambdaArray0[3] = { Lambda(k_index, 0, 0, Ei), Lambda(k_index, 0, 1, Ei) , Lambda(k_index, 0, 2, Ei)};
				Real const LambdaArray1[3] = { Lambda(k_index, 1, 0, Ei), Lambda(k_index, 1, 1, Ei) , Lambda(k_index, 1, 2, Ei)};

				for (unsigned s = 0; s < 3; s++) {

					//Sum0 += InverseBlockArray[s] * Lambda(k_index, 0, s, Ei);
					//Sum1 += InverseBlockArray[s] * Lambda(k_index, 1, s, Ei);

					//Sum0 += InverseBlock.coeff(m, s) * LambdaArray0[s];
					//Sum1 += InverseBlock.coeff(m, s) * LambdaArray1[s];

					//Sum0 += InverseBlock.coeff(m, s) * Lambda(k_index, 0, s, Ei);
					//Sum1 += InverseBlock.coeff(m, s) * Lambda(k_index, 1, s, Ei);

					Sum0 += InverseBlockArray[s] * LambdaArray0[s];
					Sum1 += InverseBlockArray[s] * LambdaArray1[s];

				}
				
				iDH1Block.coeffRef(m, Ei) = -TimeCoefficient * Sum0;
				iDH2Block.coeffRef(m, Ei) = -TimeCoefficient * Sum1;

			}
		}


		/*****************************************************************************/
		/*                                                                           */
		/*    - Assembly of the resulting matrix R1 * iD * H1						 */
		/*                                                                           */
		/*****************************************************************************/

		/*

		//Real * const R1BlockPointer = R1_block.getRowPointer(k_index, 0, 0);
		//Real * const R2BlockPointer = R2_block.getRowPointer(k_index, 0, 0);

		//Real R1BlockTimesInverseBlock[3][3];
		//Real R2BlockTimesInverseBlock[3][3];

		//for (unsigned ei = 0; ei < 3; ei++) {
		//	for (unsigned j = 0; j < 3; j++) {

		//		Real Sum0 = 0.0;
		//		Real Sum1 = 0.0;

		//		for (unsigned k = 0; k < 3; k++) {

		//			Sum0 += R1BlockPointer[3 * ei + k] * InverseBlock.coeff(k, j);
		//			Sum1 += R2BlockPointer[3 * ei + k] * InverseBlock.coeff(k, j);

		//		}

		//		R1BlockTimesInverseBlock[ei][j] = Sum0;
		//		R2BlockTimesInverseBlock[ei][j] = Sum1;
		//		
		//	}
		//}

		*/

		/*
		//__declspec(align(64)) Real const R1Row0[4] = { R1_block(k_index, 0, 0), R1_block(k_index, 0, 1), R1_block(k_index, 0, 2), 0.0 };
		//__declspec(align(64)) Real const R1Row1[4] = { R1_block(k_index, 1, 0), R1_block(k_index, 1, 1), R1_block(k_index, 1, 2), 0.0 };
		//__declspec(align(64)) Real const R1Row2[4] = { R1_block(k_index, 2, 0), R1_block(k_index, 2, 1), R1_block(k_index, 2, 2), 0.0 };

		//__declspec(align(64)) Real const R2Row0[4] = { R2_block(k_index, 0, 0), R2_block(k_index, 0, 1), R2_block(k_index, 0, 2), 0.0 };
		//__declspec(align(64)) Real const R2Row1[4] = { R2_block(k_index, 1, 0), R2_block(k_index, 1, 1), R2_block(k_index, 1, 2), 0.0 };
		//__declspec(align(64)) Real const R2Row2[4] = { R2_block(k_index, 2, 0), R2_block(k_index, 2, 1), R2_block(k_index, 2, 2), 0.0 };

		//__declspec(align(64)) Real const InverseBlockTransRow0[4] = { InverseBlock(0, 0), InverseBlock(1, 0), InverseBlock(2, 0), 0.0 };
		//__declspec(align(64)) Real const InverseBlockTransRow1[4] = { InverseBlock(0, 1), InverseBlock(1, 1), InverseBlock(2, 1), 0.0 };
		//__declspec(align(64)) Real const InverseBlockTransRow2[4] = { InverseBlock(0, 2), InverseBlock(1, 2), InverseBlock(2, 2), 0.0 };

		//Real R1TimesInverseBlock[9];
		//Real R2TimesInverseBlock[9];




		////for (unsigned ei = 0; ei < 3; ei++) {
		////	for (unsigned j = 0; j < 3; j++) {

		////		Real Sum1 = 0.0;
		////		Real Sum2 = 0.0;

		////		for (unsigned m = 0; m < 4; m++) {

		////			Sum1 += R1_block(k_index, ei, m) * InverseBlock(m, j);
		////			Sum2 += R2_block(k_index, ei, m) * InverseBlock(m, j);

		////		}

		////		R1TimesInverseBlock[ei][j] = Sum1;
		////		R2TimesInverseBlock[ei][j] = Sum2;

		////	}
		////}


		//__declspec(align(64)) Real S1um00[4];
		//__declspec(align(64)) Real S1um01[4];
		//__declspec(align(64)) Real S1um02[4];

		//__declspec(align(64)) Real S1um10[4];
		//__declspec(align(64)) Real S1um11[4];
		//__declspec(align(64)) Real S1um12[4];

		//__declspec(align(64)) Real S1um20[4];
		//__declspec(align(64)) Real S1um21[4];
		//__declspec(align(64)) Real S1um22[4];


		//__declspec(align(64)) Real S2um00[4];
		//__declspec(align(64)) Real S2um01[4];
		//__declspec(align(64)) Real S2um02[4];

		//__declspec(align(64)) Real S2um10[4];
		//__declspec(align(64)) Real S2um11[4];
		//__declspec(align(64)) Real S2um12[4];

		//__declspec(align(64)) Real S2um20[4];
		//__declspec(align(64)) Real S2um21[4];
		//__declspec(align(64)) Real S2um22[4];

		//for (unsigned m = 0; m < 4; m++) {

		//	S1um00[m] = R1Row0[m] * InverseBlockTransRow0[m];
		//	S1um01[m] = R1Row0[m] * InverseBlockTransRow1[m];
		//	S1um02[m] = R1Row0[m] * InverseBlockTransRow2[m];

		//	S1um10[m] = R1Row1[m] * InverseBlockTransRow0[m];
		//	S1um11[m] = R1Row1[m] * InverseBlockTransRow1[m];
		//	S1um12[m] = R1Row1[m] * InverseBlockTransRow2[m];

		//	S1um20[m] = R1Row2[m] * InverseBlockTransRow0[m];
		//	S1um21[m] = R1Row2[m] * InverseBlockTransRow1[m];
		//	S1um22[m] = R1Row2[m] * InverseBlockTransRow2[m];


		//	S2um00[m] = R2Row0[m] * InverseBlockTransRow0[m];
		//	S2um01[m] = R2Row0[m] * InverseBlockTransRow1[m];
		//	S2um02[m] = R2Row0[m] * InverseBlockTransRow2[m];

		//	S2um10[m] = R2Row1[m] * InverseBlockTransRow0[m];
		//	S2um11[m] = R2Row1[m] * InverseBlockTransRow1[m];
		//	S2um12[m] = R2Row1[m] * InverseBlockTransRow2[m];

		//	S2um20[m] = R2Row2[m] * InverseBlockTransRow0[m];
		//	S2um21[m] = R2Row2[m] * InverseBlockTransRow1[m];
		//	S2um22[m] = R2Row2[m] * InverseBlockTransRow2[m];

		//}

		//R1TimesInverseBlock[0 + 0 * 3] = S1um00[0] + S1um00[1] + S1um00[2];
		//R1TimesInverseBlock[1 + 0 * 3] = S1um01[0] + S1um01[1] + S1um01[2];
		//R1TimesInverseBlock[2 + 0 * 3] = S1um02[0] + S1um02[1] + S1um02[2];

		//R1TimesInverseBlock[0 + 1 * 3] = S1um10[0] + S1um10[1] + S1um10[2];
		//R1TimesInverseBlock[1 + 1 * 3] = S1um11[0] + S1um11[1] + S1um11[2];
		//R1TimesInverseBlock[2 + 1 * 3] = S1um12[0] + S1um12[1] + S1um12[2];

		//R1TimesInverseBlock[0 + 2 * 3] = S1um20[0] + S1um20[1] + S1um20[2];
		//R1TimesInverseBlock[1 + 2 * 3] = S1um21[0] + S1um21[1] + S1um21[2];
		//R1TimesInverseBlock[2 + 2 * 3] = S1um22[0] + S1um22[1] + S1um22[2];


		//R2TimesInverseBlock[0 + 0 * 3] = S2um00[0] + S2um00[1] + S2um00[2];
		//R2TimesInverseBlock[1 + 0 * 3] = S2um01[0] + S2um01[1] + S2um01[2];
		//R2TimesInverseBlock[2 + 0 * 3] = S2um02[0] + S2um02[1] + S2um02[2];

		//R2TimesInverseBlock[0 + 1 * 3] = S2um10[0] + S2um10[1] + S2um10[2];
		//R2TimesInverseBlock[1 + 1 * 3] = S2um11[0] + S2um11[1] + S2um11[2];
		//R2TimesInverseBlock[2 + 1 * 3] = S2um12[0] + S2um12[1] + S2um12[2];

		//R2TimesInverseBlock[0 + 2 * 3] = S2um20[0] + S2um20[1] + S2um20[2];
		//R2TimesInverseBlock[1 + 2 * 3] = S2um21[0] + S2um21[1] + S2um21[2];
		//R2TimesInverseBlock[2 + 2 * 3] = S2um22[0] + S2um22[1] + S2um22[2];

		*/

		for (unsigned ei = 0; ei < 3; ei++) {


			em_pointer const Ei = K->edges[ei];
			unsigned const	 e_index_i = Ei->index;


			for (unsigned j = 0; j < 3; j++) {


				unsigned const start_index_j = start_index + j;

				/*****************************************************************************/
				/*                                                                           */
				/*    - Assemble Matrices H1, H2									         */
				/*                                                                           */
				/*****************************************************************************/
				//H1.coeffRef(start_index + j, e_index_i) = Block1.coeff(j, ei);
				//H2.coeffRef(start_index + j, e_index_i) = Block2.coeff(j, ei);


				/*****************************************************************************/
				/*                                                                           */
				/*    - Assembly of the Matrices iDH1, iDH2	: iD * H1, iD * H2			     */
				/*                                                                           */
				/*****************************************************************************/
				iDH1.coeffRef(start_index_j, e_index_i) = iDH1Block.coeff(j, ei);
				iDH2.coeffRef(start_index_j, e_index_i) = iDH2Block.coeff(j, ei);


				/*****************************************************************************/
				/*                                                                           */
				/*    - Assemble Matrices R1 * D.inverse, R2 * D.inverse			         */
				/*                                                                           */
				/*****************************************************************************/
				Real Sum1 = 0.0;
				Real Sum2 = 0.0;

				Real const InverseBlockArray[3] = { InverseBlock.coeff(0, j), InverseBlock.coeff(1, j) , InverseBlock.coeff(2, j) };

				for (unsigned m = 0; m < 3; m++) {

					/*
					//Sum1 += R1BlockPointer[3 * ei + m] * InverseBlockArray[m];
					//Sum2 += R2BlockPointer[3 * ei + m] * InverseBlockArray[m];
					*/

					//Sum1 += R1_block(k_index, ei, m) * InverseBlock(m, j);
					//Sum2 += R2_block(k_index, ei, m) * InverseBlock(m, j);

					Sum1 += R1_block(k_index, ei, m) * InverseBlockArray[m];
					Sum2 += R2_block(k_index, ei, m) * InverseBlockArray[m];

				}

				//R1iD.coeffRef(e_index_i, start_index_j) = R1TimesInverseBlock[j + 3 * ei];
				//R2iD.coeffRef(e_index_i, start_index_j) = R2TimesInverseBlock[j + 3 * ei];

				R1iD.coeffRef(e_index_i, start_index_j) = Sum1;
				R2iD.coeffRef(e_index_i, start_index_j) = Sum2;


				/*
				//Real const Temp0 = R1BlockTimesInverseBlock[ei][j];
				//Real const Temp1 = R2BlockTimesInverseBlock[ei][j];

				//R1iD.coeffRef(e_index_i, start_index_j) = Temp0;
				//R2iD.coeffRef(e_index_i, start_index_j) = Temp1;
				*/

			}


			// Diagonal elements are already zeroed/or there is Mij already
			if (Ei->marker == E_MARKER::DIRICHLET)
				continue;


			//#pragma omp parallel for
			for (unsigned ej = 0; ej < 3; ej++) {


				unsigned const e_index_j = K->edges[ej]->index;


				Real sum11 = 0.0;
				Real sum12 = 0.0;
				Real sum21 = 0.0;
				Real sum22 = 0.0;

				// Number of degrees of freedom of internal pressure
				for (unsigned m = 0; m < 3; m++) {

					/*
					//sum11 += R1BlockPointer[3 * ei + m] * iDH1Block(m, ej);
					//sum12 += R1BlockPointer[3 * ei + m] * iDH2Block(m, ej);

					//sum21 += R2BlockPointer[3 * ei + m] * iDH1Block(m, ej);
					//sum22 += R2BlockPointer[3 * ei + m] * iDH2Block(m, ej);

					*/

					sum11 += R1_block(k_index, ei, m) * iDH1Block(m, ej);
					sum12 += R1_block(k_index, ei, m) * iDH2Block(m, ej);

					sum21 += R2_block(k_index, ei, m) * iDH1Block(m, ej);
					sum22 += R2_block(k_index, ei, m) * iDH2Block(m, ej);

				}

				// Because diagonal elements were zeroed at the beginning, the += operator is needed only here (if there was any. Will be when migrate to linux and No Eigen will be used)
				if (e_index_i == e_index_j) {

					Eigen::Triplet<Real> const T11(e_index_i,		e_index_i,		sum11);
					Eigen::Triplet<Real> const T12(e_index_i,		e_index_i + ne, sum12);
					Eigen::Triplet<Real> const T21(e_index_i + ne,	e_index_i,		sum21);
					Eigen::Triplet<Real> const T22(e_index_i + ne,	e_index_i + ne, sum22);

					TripletVector[count++] = T11;
					TripletVector[count++] = T12;
					TripletVector[count++] = T21;
					TripletVector[count++] = T22;

					//PressureSystem.coeffRef(e_index_i,		e_index_i)		+= sum11;
					//PressureSystem.coeffRef(e_index_i,		e_index_i + ne) += sum12;
					//PressureSystem.coeffRef(e_index_i + ne, e_index_i)		+= sum21;
					//PressureSystem.coeffRef(e_index_i + ne, e_index_i + ne) += sum22;

					continue;

				}


				Real const M11 = sum11 + M_j1_s1.coeff(e_index_i, e_index_j);
				Real const M12 = sum12 + M_j1_s2.coeff(e_index_i, e_index_j);
				Real const M21 = sum21 + M_j2_s1.coeff(e_index_i, e_index_j);
				Real const M22 = sum22 + M_j2_s2.coeff(e_index_i, e_index_j);


				Eigen::Triplet<Real> const T11(e_index_i,		e_index_j,		M11);
				Eigen::Triplet<Real> const T12(e_index_i,		e_index_j + ne, M12);
				Eigen::Triplet<Real> const T21(e_index_i + ne,	e_index_j,		M21);
				Eigen::Triplet<Real> const T22(e_index_i + ne,	e_index_j + ne, M22);

				TripletVector[count++] = T11;
				TripletVector[count++] = T12;
				TripletVector[count++] = T21;
				TripletVector[count++] = T22;

				//PressureSystem.coeffRef(e_index_i,		e_index_j)		= M11;
				//PressureSystem.coeffRef(e_index_i,		e_index_j + ne) = M12;
				//PressureSystem.coeffRef(e_index_i + ne, e_index_j)		= M21;
				//PressureSystem.coeffRef(e_index_i + ne, e_index_j + ne) = M22;

			}
		}
	}

	PressureSystem.setFromTriplets(TripletVector.cbegin(), TripletVector.cend());

};
template<unsigned QuadraturePrecision, scheme TimeScheme>
void solver<QuadraturePrecision, TimeScheme>::getSparsityPatternOfThePressureSystem() {


	Real const DummyFillIn = 1.0;

	std::vector<Eigen::Triplet<Real>> triR1iD;
	std::vector<Eigen::Triplet<Real>> triR2iD;
	std::vector<Eigen::Triplet<Real>> triiDH1;
	std::vector<Eigen::Triplet<Real>> triiDH2;


	for (unsigned e = 0; e < ne; e++) {

		Eigen::Triplet<Real> const T11(e,		e,		DummyFillIn);
		Eigen::Triplet<Real> const T12(e,		e + ne, DummyFillIn);
		Eigen::Triplet<Real> const T21(e + ne,	e,		DummyFillIn);
		Eigen::Triplet<Real> const T22(e + ne,	e + ne, DummyFillIn);

		TripletVector.push_back(T11);
		TripletVector.push_back(T12);
		TripletVector.push_back(T21);
		TripletVector.push_back(T22);

	}

	for (unsigned k = 0; k < nk; k++) {


		tm_pointer		 K			 = MeshElements[k];
		unsigned const	 k_index	 = MeshElementIndeces[k];
		unsigned const	 start_index = 3 * k_index;


		for (unsigned ei = 0; ei < 3; ei++) {


			em_pointer const Ei			= K->edges[ei];
			unsigned const	 e_index_i  = Ei->index;


			for (unsigned j = 0; j < 3; j++) {


				/*****************************************************************************/
				/*                                                                           */
				/*    - Sparsity pattern of the Matrices iD * H1, iD * H2				     */
				/*                                                                           */
				/*****************************************************************************/
				Eigen::Triplet<Real> const TH1(start_index + j, e_index_i, DummyFillIn);
				Eigen::Triplet<Real> const TH2(start_index + j, e_index_i, DummyFillIn);

				triiDH1.push_back(TH1);
				triiDH2.push_back(TH2);


				/*****************************************************************************/
				/*                                                                           */
				/*    - Sparsity pattern of the Matrices R1 * iD, R2 * iD			         */
				/*                                                                           */
				/*****************************************************************************/
				Eigen::Triplet<Real> const TR1(e_index_i, start_index + j, DummyFillIn);
				Eigen::Triplet<Real> const TR2(e_index_i, start_index + j, DummyFillIn);

				triR1iD.push_back(TR1);
				triR2iD.push_back(TR2);


			}



			// Diagonal elements are already zeroed/or there is Mij already
			if (Ei->marker == E_MARKER::DIRICHLET)
				continue;


			/*****************************************************************************/
			/*                                                                           */
			/*    - Sparsity patter of the Matrix PressureSystem						 */
			/*                                                                           */
			/*****************************************************************************/
			for (unsigned ej = 0; ej < 3; ej++) {


				unsigned const e_index_j = K->edges[ej]->index;

				// Because diagonal elements were zeroed at the beginning, the += operator is needed only here (if there was any. Will be when migrate to linux and No Eigen will be used)
				if (e_index_i == e_index_j) {

					Eigen::Triplet<Real> const T11(e_index_i,		e_index_i,		DummyFillIn);
					Eigen::Triplet<Real> const T12(e_index_i,		e_index_i + ne, DummyFillIn);
					Eigen::Triplet<Real> const T21(e_index_i + ne,	e_index_i,		DummyFillIn);
					Eigen::Triplet<Real> const T22(e_index_i + ne,	e_index_i + ne, DummyFillIn);

					TripletVector.push_back(T11);
					TripletVector.push_back(T12);
					TripletVector.push_back(T21);
					TripletVector.push_back(T22);

					continue;

				}


				Eigen::Triplet<Real> const T11(e_index_i,		e_index_j,		DummyFillIn);
				Eigen::Triplet<Real> const T12(e_index_i,		e_index_j + ne, DummyFillIn);
				Eigen::Triplet<Real> const T21(e_index_i + ne,	e_index_j,		DummyFillIn);
				Eigen::Triplet<Real> const T22(e_index_i + ne,	e_index_j + ne, DummyFillIn);

				TripletVector.push_back(T11);
				TripletVector.push_back(T12);
				TripletVector.push_back(T21);
				TripletVector.push_back(T22);

			}
		}
	}

	R1iD			.setFromTriplets(triR1iD.cbegin(), triR1iD.cend());
	R2iD			.setFromTriplets(triR2iD.cbegin(), triR2iD.cend());
	iDH1			.setFromTriplets(triiDH1.cbegin(), triiDH1.cend());
	iDH2			.setFromTriplets(triiDH2.cbegin(), triiDH2.cend());

	PressureSystem.	setFromTriplets(TripletVector.cbegin(), TripletVector.cend());
	
	BiConjugateGradientSolver		.analyzePattern(PressureSystem);
	//sparseLUsolver_PressureSystem	.analyzePattern(PressureSystem);

	TripletVector.shrink_to_fit();
	
	//BiConjugateGradientSolver.setTolerance(1-10);

	R1iD.setZero();
	R2iD.setZero();
	iDH1.setZero();
	iDH2.setZero();

	PressureSystem.setZero();

};



template<unsigned QuadraturePrecision, scheme TimeScheme>
void solver<QuadraturePrecision, TimeScheme>::assemble_Alpha() {


	quadrature_triangle<Real> const QuadratureOnTriangle(QuadraturePrecision);

	Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> Integral(8, 8);
	Eigen::Matrix<Real, 2, 8> BasisRaviartThomas;
	Eigen::Matrix<Real, 2, 2> JF;

	for (unsigned k = 0; k < nk; k++) {


		//tm_pointer const K = Mesh->get_triangle(k);
		tm_pointer const K		 = MeshElements[k];
		unsigned const   k_index = K->index;

		vm_pointer const va = K->vertices[0];
		vm_pointer const vb = K->vertices[1];
		vm_pointer const vc = K->vertices[2];

		Real const x0 = va->x;
		Real const y0 = va->y;

		Real const x1 = vb->x;
		Real const y1 = vb->y;

		Real const x2 = vc->x;
		Real const y2 = vc->y;

		Real const detJF = abs((x1 - x0)*(y2 - y0) - (x2 - x0)*(y1 - y0));

		JF.coeffRef(0, 0) = x1 - x0;
		JF.coeffRef(0, 1) = x2 - x0;
		JF.coeffRef(1, 0) = y1 - y0;
		JF.coeffRef(1, 1) = y2 - y0;


		Integral.setZero();

		for (unsigned n = 0; n < NumberOfQuadraturePointsTriangle; n++) {


			Real const s = QuadratureOnTriangle.points_x[n];
			Real const t = QuadratureOnTriangle.points_y[n];
			Real const w = 0.5 * QuadratureOnTriangle.weights[n];

			Real const x = x0 + JF(0, 0) * s + JF(0, 1) * t;
			Real const y = y0 + JF(1, 0) * s + JF(1, 1) * t;


			evaluate_raviartthomas_basis(s, t, BasisRaviartThomas);


			// Get the inverse of the permeability tensor K
			Eigen::Matrix<Real, 2, 2> K;

			permeability(x, y, K);

			Eigen::Matrix<Real, 2, 2> const iK = K.inverse();


			for (unsigned i = 0; i < 8; i++) {


				Eigen::Matrix<Real, 2, 1> const JFWi = JF * BasisRaviartThomas.col(i);

				for (unsigned j = 0; j < 8; j++) {


					Eigen::Matrix<Real, 2, 1> const JFWj = JF * BasisRaviartThomas.col(j);

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
				Alpha.setCoeff(k_index, i, j) = abs(Integral(i, j)) < INTEGRAL_PRECISION ? 0.0 : Integral(i, j);

	}

};
template<unsigned QuadraturePrecision, scheme TimeScheme>
void solver<QuadraturePrecision, TimeScheme>::assemble_Beta() {


	quadrature_triangle<Real> const QuadratureOnTriangle(QuadraturePrecision);

	Eigen::Matrix<Real, 8, 3> Integral;
	Eigen::Matrix<Real, 3, 1> BasisPolynomial;
	Eigen::Matrix<Real, 8, 1> BasisRaviartThomasDivergence;


	Integral.setZero();

	for (unsigned n = 0; n < NumberOfQuadraturePointsTriangle; n++) {


		Real const s = QuadratureOnTriangle.points_x[n];
		Real const t = QuadratureOnTriangle.points_y[n];
		Real const w = 0.5 * QuadratureOnTriangle.weights[n];


		evaluate_raviartthomas_basis_divergence(s, t, BasisRaviartThomasDivergence);
		evaluate_polynomial_basis(s, t, BasisPolynomial);


		for (unsigned i = 0; i < 8; i++) {

			Real const dWi = BasisRaviartThomasDivergence(i);

			for (unsigned j = 0; j < 3; j++) {

				Real const Phij = BasisPolynomial(j);

				Integral.coeffRef(i, j) += w * dWi * Phij;

			}
		}
	}

	for (unsigned m = 0; m < 8; m++)
		for (unsigned j = 0; j < 3; j++)
			Beta.setCoeff(m, j) = abs(Integral(m, j)) < INTEGRAL_PRECISION ? 0.0 : Integral(m, j);


};
template<unsigned QuadraturePrecision, scheme TimeScheme>
void solver<QuadraturePrecision, TimeScheme>::assemble_Chi() {


	gauss_quadrature_1D<Real> const QuadratureOnEdge(QuadraturePrecision);


	Eigen::Matrix<Real, 2, 3> ReferenceNormals;
	Eigen::Matrix<Real, 2, 1> Parametrization;
	Eigen::Matrix<Real, 2, 8> BasisRaviartThomas;
	Eigen::Matrix<Real, 2, 1> BasisEdgePolynomial;

	evaluate_edge_normal(ReferenceNormals);

	Chi.setZero();

	for (unsigned e = 0; e < ne; e++) {


		em_pointer const E		 = MeshEdges[e];
		unsigned const	 e_index = E->index;


		for (unsigned neighborElement = 0; neighborElement < 2; neighborElement++) {


			tm_pointer const K = E->neighbors[neighborElement];

			if (!K)
				continue;


			unsigned const k_index = K->index;
			unsigned const e_index_local = K->get_edge_index(E);

			Real const orientation = edgeOrientation(k_index, e_index_local);

			Real const a = 0.0;
			Real const b = e_index_local != 0 ? 1.0 : sqrt(2.0);

			Real const c = 0.5 * (b - a);
			Real const d = 0.5 * (b + a);


			Eigen::Matrix<Real, 2, 1> const ReferenceNormal = ReferenceNormals.col(e_index_local);


			for (unsigned n = 0; n < NumberOfQuadraturePointsEdge; n++) {


				Real const x = QuadratureOnEdge.points[n] * c + d;
				Real const w = QuadratureOnEdge.weights[n] * c;


				evaluate_edge_parametrization(x, e_index_local, Parametrization);

				Real const s = Parametrization(0);
				Real const t = Parametrization(1);
				Real const drNorm = 1.0;


				evaluate_raviartthomas_basis(s, t, BasisRaviartThomas);
				evaluate_edge_polynomial_basis(x, e_index_local, BasisEdgePolynomial, orientation);


				for (unsigned m = 0; m < 8; m++) {


					Eigen::Matrix<Real, 2, 1> const	Wm = BasisRaviartThomas.col(m);
					Real const						dotProduct = Wm.dot(ReferenceNormal);


					for (unsigned s = 0; s < 2; s++) {

						Real const varPhis = BasisEdgePolynomial(s);

						Chi.setCoeff(k_index, m, e_index_local, s) = Chi(k_index, m, e_index_local, s) + w * dotProduct * varPhis * drNorm;

					}
				}
			}


			for (unsigned m = 0; m < 8; m++)
				for (unsigned s = 0; s < 2; s++)
					Chi.setCoeff(k_index, m, e_index_local, s) = abs(Chi(k_index, m, e_index_local, s)) < INTEGRAL_PRECISION ? 0.0 : Chi(k_index, m, e_index_local, s);

		}
	}

};
template<unsigned QuadraturePrecision, scheme TimeScheme>
void solver<QuadraturePrecision, TimeScheme>::assemble_Eta() {


	quadrature_triangle<Real> const QuadratureOnTriangle(QuadraturePrecision);

	Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> Integral(3, 3);
	Eigen::Matrix<Real, 3, 1>							BasisPolynomial;


	Integral.setZero();

	for (unsigned n = 0; n < NumberOfQuadraturePointsTriangle; n++) {


		Real const s = QuadratureOnTriangle.points_x[n];
		Real const t = QuadratureOnTriangle.points_y[n];
		Real const w = 0.5 * QuadratureOnTriangle.weights[n];

		evaluate_polynomial_basis(s, t, BasisPolynomial);


		for (unsigned m = 0; m < 3; m++) {


			Real const Phim = BasisPolynomial(m);

			for (unsigned j = 0; j < 3; j++) {


				Real const Phij = BasisPolynomial(j);

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
			Eta.setCoeff(i, j) = abs(Integral(i, j)) < INTEGRAL_PRECISION ? 0.0 : Integral(i, j);

};
template<unsigned QuadraturePrecision, scheme TimeScheme>
void solver<QuadraturePrecision, TimeScheme>::assemble_Tau() {


	quadrature_triangle<Real> const QuadratureOnTriangle(QuadraturePrecision);

	Eigen::Matrix<Real, 2, 8> BasisRaviartThomas;
	Eigen::Matrix<Real, 3, 1> BasisPolynomial;
	Eigen::Matrix<Real, 2, 3> BasisPolynomialGradient;


	for (unsigned n = 0; n < NumberOfQuadraturePointsTriangle; n++) {


		Real const s = QuadratureOnTriangle.points_x[n];
		Real const t = QuadratureOnTriangle.points_y[n];
		Real const w = 0.5 * QuadratureOnTriangle.weights[n];


		evaluate_raviartthomas_basis(s, t, BasisRaviartThomas);
		evaluate_polynomial_basis(s, t, BasisPolynomial);
		evaluate_polynomial_basis_gradient(s, t, BasisPolynomialGradient);


		for (unsigned m = 0; m < 3; m++) {


			Eigen::Matrix<Real, 2, 1> const dPhim = BasisPolynomialGradient.col(m);

			for (unsigned j = 0; j < 8; j++) {


				Eigen::Matrix<Real, 2, 1> const	Wj = BasisRaviartThomas.col(j);
				Real const						dotProduct = Wj.dot(dPhim);

				for (unsigned l = 0; l < 3; l++) {

					Real const Phil = BasisPolynomial(l);

					Tau.setCoeff(m, j, l) = Tau(m, j, l) + w * dotProduct * Phil;

				}
			}
		}
	}

	for (unsigned m = 0; m < 3; m++)
		for (unsigned i = 0; i < 8; i++)
			for (unsigned j = 0; j < 3; j++)
				Tau.setCoeff(m, i, j) = abs(Tau(m, i, j)) < INTEGRAL_PRECISION ? 0.0 : Tau(m, i, j);

};
template<unsigned QuadraturePrecision, scheme TimeScheme>
void solver<QuadraturePrecision, TimeScheme>::assemble_Delta() {


	Delta.setZero();

	Eigen::Matrix<Real, 3, 8> Integral;

	for (unsigned k = 0; k < nk; k++) {


		tm_pointer const K		 = MeshElements[k];
		unsigned const	 k_index = K->index;


		for (unsigned El = 0; El < 3; El++) {


			Integral.setZero();

			for (unsigned n = 0; n < NumberOfQuadraturePointsEdge; n++) {


				Real const tC = upwindConcentration(K, El, n);
				
				for (unsigned m = 0; m < 3; m++) {

					Integral.coeffRef(m, 0) += tC * QuadraturePoints_RaviartThomasBasisDotNormalTimesPolynomialBasis(n, 0, El, m);
					Integral.coeffRef(m, 1) += tC * QuadraturePoints_RaviartThomasBasisDotNormalTimesPolynomialBasis(n, 1, El, m);
					Integral.coeffRef(m, 2) += tC * QuadraturePoints_RaviartThomasBasisDotNormalTimesPolynomialBasis(n, 2, El, m);
					Integral.coeffRef(m, 3) += tC * QuadraturePoints_RaviartThomasBasisDotNormalTimesPolynomialBasis(n, 3, El, m);
					Integral.coeffRef(m, 4) += tC * QuadraturePoints_RaviartThomasBasisDotNormalTimesPolynomialBasis(n, 4, El, m);
					Integral.coeffRef(m, 5) += tC * QuadraturePoints_RaviartThomasBasisDotNormalTimesPolynomialBasis(n, 5, El, m);
					Integral.coeffRef(m, 6) += tC * QuadraturePoints_RaviartThomasBasisDotNormalTimesPolynomialBasis(n, 6, El, m);
					Integral.coeffRef(m, 7) += tC * QuadraturePoints_RaviartThomasBasisDotNormalTimesPolynomialBasis(n, 7, El, m);

					//for (unsigned j = 0; j < 8; j++)
					//	Integral.coeffRef(m, j) += tC * QuadraturePoints_RaviartThomasBasisDotNormalTimesPolynomialBasis(n, j, El, m);

				}
				

			}

			for (unsigned m = 0; m < 3; m++)
				for (unsigned j = 0; j < 8; j++)
					Delta.setCoeff(k_index, m, El, j) = Integral(m, j);

		}
	}

};
template<unsigned QuadraturePrecision, scheme TimeScheme>
void solver<QuadraturePrecision, TimeScheme>::assemble_Gamma() {

//	omp_set_num_threads(2);
//#pragma omp parallel for 
	for (int k = 0; k < nk; k++) {


		unsigned const k_index = MeshElementIndeces[k];

		for (unsigned m = 0; m < 3; m++) {

			for (unsigned j = 0; j < 8; j++) {


				/*
				Real Value1 = 0.0;
				Real Value2 = 0.0;

				for (unsigned El = 0; El < 3; El++)
					Value1 += Delta(k_index, m, El, j);

				for (unsigned l = 0; l < 3; l++)
					Value2 += Tau(m, j, l) * Xi_prev(k_index, l);
				*/

				/*****************************************************************************/
				/*                                                                           */
				/*    - Sum over element's edges El = 0,1,2									 */
				/*                                                                           */
				/*****************************************************************************/
				Real const Value1 = Delta(k_index, m, 0, j) + Delta(k_index, m, 1, j) + Delta(k_index, m, 2, j);

				/*****************************************************************************/
				/*                                                                           */
				/*    - Sum over DG degrees of freedom l = 0,1,2							 */
				/*                                                                           */
				/*****************************************************************************/
				Real const Value2 = Tau(m, j, 0) * Xi_prev(k_index, 0) + Tau(m, j, 1) * Xi_prev(k_index, 1) + Tau(m, j, 2) * Xi_prev(k_index, 2);

				Gamma.setCoeff(k_index, m, j) = Value1 - Value2;

			}
		}
	}

};



template<unsigned QuadraturePrecision, scheme TimeScheme>
void solver<QuadraturePrecision, TimeScheme>::assemble_Sigma() {

//	omp_set_num_threads(2);
//#pragma omp parallel for 
	for (int k = 0; k < nk; k++) {


		unsigned const	k_index		= MeshElementIndeces[k];
		Real const		Coefficient = -Thetas_prev[k_index] * PorosityViscosityDeterminant[k_index];


		for (unsigned m = 0; m < 3; m++) {
			for (unsigned l = 0; l < 3; l++) {


				/*
				Real Value = 0.0;

				for (unsigned q = 0; q < 3; q++) {

					Real GammaAlphaBeta = 0.0;

					for (unsigned j = 0; j < 8; j++)
						GammaAlphaBeta += Gamma(k_index, q, j) * AlphaTimesBeta(k_index, j, l);

					Value += Eta(m, q) * GammaAlphaBeta;

				}
				*/

				/*****************************************************************************/
				/*                                                                           */
				/*    - Sum over DG degrees of freedom q = 0,1,2							 */
				/*                                                                           */
				/*****************************************************************************/
				Real GammaAlphaBeta0 = 0.0;
				Real GammaAlphaBeta1 = 0.0;
				Real GammaAlphaBeta2 = 0.0;

				for (unsigned j = 0; j < 8; j++) {

					Real const AlphaBeta = AlphaTimesBeta(k_index, j, l);

					GammaAlphaBeta0 += Gamma(k_index, 0, j) * AlphaBeta;
					GammaAlphaBeta1 += Gamma(k_index, 1, j) * AlphaBeta;
					GammaAlphaBeta2 += Gamma(k_index, 2, j) * AlphaBeta;

				}

				Real const Value = Eta(m, 0) * GammaAlphaBeta0 + Eta(m, 1) * GammaAlphaBeta1 + Eta(m, 2) * GammaAlphaBeta2;

				// In the computation of Eta, there is coefficient detJF. When inverting, the coefficient is inverted
				Sigma.setCoeff(k_index, m, l) = Coefficient * Value;

			}
		}
	}

};
template<unsigned QuadraturePrecision, scheme TimeScheme>
void solver<QuadraturePrecision, TimeScheme>::assemble_Lambda() {

//	omp_set_num_threads(2);
//#pragma omp parallel for 
	for (int k = 0; k < nk; k++) {

		unsigned const	k_index		= MeshElementIndeces[k];
		Real const		Coefficient = Thetas_prev[k_index] * PorosityViscosityDeterminant[k_index];


		for (unsigned s = 0; s < 2; s++) {
			for (unsigned m = 0; m < 3; m++) {
				for (unsigned El = 0; El < 3; El++) {


					/*
					Real Value = 0.0;

					for (unsigned q = 0; q < 3; q++) {

						Real GammaAlphaChi = 0.0;

						for (unsigned j = 0; j < 8; j++)
							GammaAlphaChi += Gamma(k_index, q, j) * AlphaTimesChi(k_index, j, El, s);

						Value += Eta(m, q) * GammaAlphaChi;

					}
					*/

					/*****************************************************************************/
					/*                                                                           */
					/*    - Sum over DG degrees of freedom q = 0,1,2							 */
					/*                                                                           */
					/*****************************************************************************/
					Real GammaAlphaChi0 = 0.0;
					Real GammaAlphaChi1 = 0.0;
					Real GammaAlphaChi2 = 0.0;

					for (unsigned j = 0; j < 8; j++) {

						Real const AlphaChi = AlphaTimesChi(k_index, j, El, s);

						GammaAlphaChi0 += Gamma(k_index, 0, j) * AlphaChi;
						GammaAlphaChi1 += Gamma(k_index, 1, j) * AlphaChi;
						GammaAlphaChi2 += Gamma(k_index, 2, j) * AlphaChi;

					}

					Real const Value = Eta(m, 0) * GammaAlphaChi0 + Eta(m, 1) * GammaAlphaChi1 + Eta(m, 2) * GammaAlphaChi2;

					// In the computation of Eta, there is coefficient detJF. When inverting, the coefficient is inverted
					Lambda.setCoeff(k_index, s, m, El) = Coefficient * Value;

				}
			}
		}
	}

};
template<unsigned QuadraturePrecision, scheme TimeScheme>
void solver<QuadraturePrecision, TimeScheme>::assemble_BigPhi() {};
template<unsigned QuadraturePrecision, scheme TimeScheme>
void solver<QuadraturePrecision, TimeScheme>::assemble_BigPsi() {};
template<unsigned QuadraturePrecision, scheme TimeScheme>
void solver<QuadraturePrecision, TimeScheme>::assemble_BigOmega() {};


template<unsigned QuadraturePrecision, scheme TimeScheme>
void solver<QuadraturePrecision, TimeScheme>::getSolution() {



	/*****************************************************************************/
	/*                                                                           */
	/*    - Compute initial trace pressures from known inner pressures and	     */
	/*      velocity field and coefficients Thetas							     */
	/*                                                                           */
	/*    - This is initializing step. This is first guess for the               */
	/*      trace pressures and velocities on the time level nt = (n + 1),       */
	/*      therefore quantities are computed on the (n + 1)-th time level       */
	/*                                                                           */
	/*****************************************************************************/
	computeTracePressures();
	computeVelocities();
	computeThetas();

	//std::string fileName_velocity = "C:\\Users\\pgali\\Desktop\\eoc\\velocity_";
	//exportVelocities(fileName_velocity + std::to_string(nt) + ".txt");


	/*****************************************************************************/
	/*                                                                           */
	/*    - Set unknowns for inner iterations. Set iteration number: l = 0	     */
	/*			: Pi_n	     - pressures on the n-th time level				 	 */
	/*			: Pi_prev     - pressures on the (n+1),(l-1)-th time level		 */
	/*			: Xi_n	     - concentrations on the n-th time level		 	 */
	/*			: Xi_prev     - concentrations on the (n+1),(l-1)-th time level	 */
	/*			: Thetas      - Thetas on the n-th time level		 			 */
	/*			: Thetas_prev - Thetas on the (n+1),(l-1)-th time level			 */
	/*                                                                           */
	/*****************************************************************************/
	Pi_n = Pi;
	Pi_prev = Pi;

	Xi_n = Xi;
	Xi_prev = Xi;

	hardCopy(Thetas_prev, Thetas, nk);


	/*****************************************************************************/
	/*                                                                           */
	/*    - Right-hand sides of the equation y'(t) = F(x,t,y(t))			     */
	/*      arrising from Runge-Kutta schemes									 */
	/*			: rkFc_n - RHS for concentrations on the n-th time level	     */
	/*			: rkFc   - RHS for concentrations on (n+1),(l-1)-th time level	 */
	/*			: rkFp_n - RHS for pressures on the n-th time level				 */
	/*                                                                           */
	/*****************************************************************************/
	rkFc_n = rkFc;



	unsigned counter = 0;

	while (counter < MAX_ITERATIONS) {


		assemble_Delta();
		assemble_Gamma();

		assemble_Sigma();
		assemble_Lambda();


		//omp_set_num_threads(2);
		//#pragma omp parallel
		//{
		//	#pragma omp sections
		//	{
		//		#pragma omp section
		//		{
		//			assemble_Sigma();
		//		}
		//		#pragma omp section
		//		{
		//			assemble_Lambda();
		//		}
		//	}
		//}


		computePressureEquation();
		computeVelocities();

		updateConcentrations();

		//concentrationCorrection();
		//computeThetas();

		if (stopCriterion())
			break;

		/*****************************************************************************/
		/*                                                                           */
		/*    - Set next iteration number: l = l + 1								 */
		/*                                                                           */
		/*****************************************************************************/
		Pi_prev = Pi;
		Xi_prev = Xi;

		//HardCopy(Thetas_prev, Thetas, nk);

		counter++;

	}


	/*****************************************************************************/
	/*                                                                           */
	/*    - This is needed for the next time level as the part of the			 */
	/*      Right-Hand Side of the pressure system - Crank-Nicolson              */
	/*                                                                           */
	/*****************************************************************************/
	Xi_prev = Xi;

	assemble_Delta();
	assemble_Gamma();

	assemble_Lambda();
	assemble_Sigma();

	for (unsigned k = 0; k < nk; k++) {


		unsigned const k_index = MeshElementIndeces[k];

		for (unsigned m = 0; m < 3; m++) {

			Real Value1 = 0.0;
			Real Value2 = 0.0;
			
			for (unsigned j = 0; j < 3; j++)
				Value1 += Sigma(k_index, m, j) * Pi(k_index, j);

			for (unsigned El = 0; El < 3; El++)
				for (unsigned s = 0; s < 2; s++)
					Value2 += Lambda(k_index, s, m, El) * TPi(k_index, El, s);

			rkFp_n.setCoeff(k_index, m) = Value1 + Value2;

		}
	}


	//std::cout << counter << std::endl;

	//if (nt % 1000 == 0)
	//std::cout << nt << " - Iterations : " << counter << std::endl;

};







/*****************************************************************************/
/*                                                                           */
/*    - Export computed quantites into .txt file							 */
/*                                                                           */
/*****************************************************************************/
template<unsigned QuadraturePrecision, scheme TimeScheme>
void solver<QuadraturePrecision, TimeScheme>::exportPressures(std::string const & fileName) {


	std::ofstream OFSTxtFile(fileName);

	for (unsigned k = 0; k < nk; k++) {


		tm_pointer const K = MeshElements[k];
		unsigned const	 k_index = MeshElementIndeces[k];

		vm_pointer const a = K->vertices[0];
		vm_pointer const b = K->vertices[1];
		vm_pointer const c = K->vertices[2];

		Real const x0 = a->x;
		Real const y0 = a->y;

		Real const x1 = b->x;
		Real const y1 = b->y;

		Real const x2 = c->x;
		Real const y2 = c->y;

		Real const x[3] = { x0, x1, x2 };
		Real const y[3] = { y0, y1, y2 };

		Real const S[3] = { 0.0, 1.0, 0.0 };
		Real const T[3] = { 0.0, 0.0, 1.0 };

		//Real const S[4] = { 0.0, 1.0, 0.0, 0.0 };
		//Real const T[4] = { 0.0, 0.0, 1.0, 0.0 };
		
		for (unsigned i = 0; i < 3; i++)
			OFSTxtFile << std::setprecision(20) << x[i] << "\t" << y[i] << "\t" << Pi(k_index, 0) * phi0(S[i], T[i]) + Pi(k_index, 1) * phi1(S[i], T[i]) + Pi(k_index, 2) * phi2(S[i], T[i]) << std::endl;

		OFSTxtFile << std::endl;

	}

	OFSTxtFile.close();

};
template<unsigned QuadraturePrecision, scheme TimeScheme>
void solver<QuadraturePrecision, TimeScheme>::exportConcentrations(std::string const & fileName) {


	std::ofstream OFSTxtFile(fileName);

	for (unsigned k = 0; k < nk; k++) {


		tm_pointer const K = MeshElements[k];
		unsigned const	 k_index = MeshElementIndeces[k];

		vm_pointer const a = K->vertices[0];
		vm_pointer const b = K->vertices[1];
		vm_pointer const c = K->vertices[2];

		Real const x0 = a->x;
		Real const y0 = a->y;

		Real const x1 = b->x;
		Real const y1 = b->y;

		Real const x2 = c->x;
		Real const y2 = c->y;

		Real const x[3] = { x0, x1, x2 };
		Real const y[3] = { y0, y1, y2 };

		Real const S[3] = { 0.0, 1.0, 0.0 };
		Real const T[3] = { 0.0, 0.0, 1.0 };

		//Real const S[4] = { 0.0, 1.0, 0.0, 0.0 };
		//Real const T[4] = { 0.0, 0.0, 1.0, 0.0 };

		for (unsigned i = 0; i < 3; i++)
			OFSTxtFile << std::setprecision(20) << x[i] << "\t" << y[i] << "\t" << Xi(k_index, 0) * phi0(S[i], T[i]) + Xi(k_index, 1) * phi1(S[i], T[i]) + Xi(k_index, 2) * phi2(S[i], T[i]) << std::endl;

		OFSTxtFile << std::endl;

	}

	OFSTxtFile.close();

};

template<unsigned QuadraturePrecision, scheme TimeScheme>
void solver<QuadraturePrecision, TimeScheme>::exportTracePressures(std::string const & fileName) {


	//std::ofstream OFSTxtFile(fileName);
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
	//	OFSTxtFile << x0 << " " << y0 << " " << std::setprecision(20) << tp[e] << std::endl;
	//	OFSTxtFile << x1 << " " << y1 << " " << std::setprecision(20) << tp[e] << std::endl;
	//	OFSTxtFile << std::endl;
	//}
	//
	//OFSTxtFile.close();

};

template<unsigned QuadraturePrecision, scheme TimeScheme>
void solver<QuadraturePrecision, TimeScheme>::exportVelocityField(std::string const & fileName) {



	quadrature_triangle<Real> const QuadratureOnTriangle(4);
	unsigned const					NumberOfQuadraturePoints = QuadratureOnTriangle.NumberOfPoints;

	Eigen::MatrixXd BasisRaviartThomas(2, 8);
	Eigen::Matrix2d JF;

	std::ofstream OFSTxtFile(fileName);

	for (unsigned k = 0; k < nk; k++) {


		tm_pointer const K = MeshElements[k];
		unsigned const	 k_index = MeshElementIndeces[k];

		vm_pointer const a = K->vertices[0];
		vm_pointer const b = K->vertices[1];
		vm_pointer const c = K->vertices[2];

		Real const x0 = a->x;
		Real const y0 = a->y;

		Real const x1 = b->x;
		Real const y1 = b->y;

		Real const x2 = c->x;
		Real const y2 = c->y;

		Real const detJF = abs((x1 - x0)*(y2 - y0) - (x2 - x0)*(y1 - y0));

		JF(0, 0) = x1 - x0;
		JF(0, 1) = x2 - x0;
		JF(1, 0) = y1 - y0;
		JF(1, 1) = y2 - y0;


		for (unsigned n = 0; n < NumberOfQuadraturePoints; n++) {


			Real const s = QuadratureOnTriangle.points_x[n];
			Real const t = QuadratureOnTriangle.points_y[n];

			Real const x = x0 + JF(0, 0) * s + JF(0, 1) * t;
			Real const y = y0 + JF(1, 0) * s + JF(1, 1) * t;

			evaluate_raviartthomas_basis(s, t, BasisRaviartThomas);


			Eigen::Vector2d Velocity(0.0, 0.0);

			for (unsigned i = 0; i < 8; i++)
				Velocity += Velocity(k_index, i) * JF * BasisRaviartThomas.col(i) / detJF;


			OFSTxtFile << std::setprecision(20) << x << "\t" << y << "\t" << Velocity(0) << "\t" << Velocity(1) << std::endl;

		}
	}

	OFSTxtFile.close();

};
template<unsigned QuadraturePrecision, scheme TimeScheme>
void solver<QuadraturePrecision, TimeScheme>::exportVelocities(std::string const & fileName) {


	std::ofstream OFSTxtFile(fileName);

	for (unsigned k = 0; k < nk; k++) {

		for (unsigned j = 0; j < 8; j++)
			OFSTxtFile << Velocity(MeshElementIndeces[k], j) << " ";

		OFSTxtFile << std::endl;
	}

	OFSTxtFile.close();

};

template<unsigned QuadraturePrecision, scheme TimeScheme>
void solver<QuadraturePrecision, TimeScheme>::computeError(std::string const & fileName) {


	/*****************************************************************************/
	/*                                                                           */
	/*    - When leaving time for-cycle, time is still on the NT-th level		 */
	/*      but the solution is already on (NT+1)-th level						 */
	/*                                                                           */
	/*****************************************************************************/
	Real const time = (nt + 1) * dt;


	quadrature_triangle<Real> const	QuadratureOnTriangle(QuadraturePrecision);
	unsigned const					NumberOfQuadraturePoints = QuadratureOnTriangle.NumberOfPoints;

	Eigen::Vector3d BasisPolynomial(3);
	Eigen::Matrix2d JF;

	Real ErrorL1 = 0.0;
	Real ErrorL2 = 0.0;
	Real ErrorMax = 0.0;

	for (unsigned k = 0; k < nk; k++) {


		tm_pointer const K = MeshElements[k];
		unsigned const	 k_index = MeshElementIndeces[k];

		vm_pointer const va = K->vertices[0];
		vm_pointer const vb = K->vertices[1];
		vm_pointer const vc = K->vertices[2];

		Real const x0 = va->x;
		Real const y0 = va->y;

		Real const x1 = vb->x;
		Real const y1 = vb->y;

		Real const x2 = vc->x;
		Real const y2 = vc->y;

		Real const HalfDetJF = 0.5 * abs((x1 - x0)*(y2 - y0) - (x2 - x0)*(y1 - y0));

		JF.coeffRef(0, 0) = x1 - x0;
		JF.coeffRef(0, 1) = x2 - x0;
		JF.coeffRef(1, 0) = y1 - y0;
		JF.coeffRef(1, 1) = y2 - y0;


		Real const B0 = barenblatt(x0, y0, time);
		Real const B1 = barenblatt(x1, y1, time);
		Real const B2 = barenblatt(x2, y2, time);

		Eigen::Vector3d const B(B0, B1, B2);
		Eigen::Matrix3d M;
		M << 1.0, -1.0, -1.0,
			1.0, +1.0, -1.0,
			1.0, -1.0, +1.0;
		Eigen::Vector3d const Solution = M.inverse() * B;


		Real L1NormOnElement = 0.0;
		Real L2NormOnElement = 0.0;
		Real MaxNormOnElement = 0.0;

		for (unsigned n = 0; n < NumberOfQuadraturePoints; n++) {


			//for (unsigned n = 0; n < NumberOfQuadraturePointsTriangle; n++) {

			//Real const s = QuadraturePointsAndWeightsOnReferenceTriangle(n, 0);
			//Real const t = QuadraturePointsAndWeightsOnReferenceTriangle(n, 1);
			//Real const w = QuadraturePointsAndWeightsOnReferenceTriangle(n, 2);


			Real const s = QuadratureOnTriangle.points_x[n];
			Real const t = QuadratureOnTriangle.points_y[n];
			Real const w = QuadratureOnTriangle.weights[n];

			Real const X = x0 + JF(0, 0) * s + JF(0, 1) * t;
			Real const Y = y0 + JF(1, 0) * s + JF(1, 1) * t;

			evaluate_polynomial_basis(s, t, BasisPolynomial);

			Real PressureK = 0.0;

			for (unsigned j = 0; j < 3; j++)
				PressureK += Pi(k_index, j) * BasisPolynomial(j);

			Real const Difference = abs(PressureK - barenblatt(X, Y, time));


			L1NormOnElement += w * Difference;
			L2NormOnElement += w * square(Difference);
			MaxNormOnElement = Difference > MaxNormOnElement ? Difference : MaxNormOnElement;


			//real BarenblattK = 0.0;
			//for (unsigned j = 0; j < 3; j++)
			//	BarenblattK += Solution(j) * BasisPolynomial(j);
			//
			//Integral += w * square(PressureK - BarenblattK);

		}

		ErrorL1 += HalfDetJF * L1NormOnElement;
		ErrorL2 += HalfDetJF * L2NormOnElement;
		ErrorMax = MaxNormOnElement;

	}

	std::ofstream OFSTxtFile(fileName);

	OFSTxtFile << "#L1 L2 MAX" << std::endl;

	OFSTxtFile << std::setprecision(20) << ErrorL1 << std::endl;
	OFSTxtFile << std::setprecision(20) << sqrt(ErrorL2) << std::endl;
	OFSTxtFile << std::setprecision(20) << ErrorMax << std::endl;

	OFSTxtFile.close();

	std::cout << "Error L1	: " << ErrorL1 << std::endl;
	std::cout << "Error L2	: " << sqrt(ErrorL2) << std::endl;
	std::cout << "Error Max	: " << ErrorMax << std::endl;

};





/*****************************************************************************/
/*                                                                           */
/*    - Evaluation of Raviart-Thomas and P1 basis function, normal vectors,  */
/*      edges parametrization, Edge basis P1(E)							     */
/*                                                                           */
/*****************************************************************************/
template<unsigned QuadraturePrecision, scheme TimeScheme>
inline void solver<QuadraturePrecision, TimeScheme>::evaluate_raviartthomas_basis(Real const s, Real const t, Eigen::Matrix<Real, 2, 8> & out) {


	out.coeffRef(0, 0) = -3.0*s + 4.0*s*t + 4.0*s * s;
	out.coeffRef(1, 0) = -3.0*t + 4.0*s*t + 4.0*t * t;

	out.coeffRef(0, 1) = -1.0 + 5.0*s - 4.0*s * s;
	out.coeffRef(1, 1) = t - 4.0*s * t;

	out.coeffRef(0, 2) = s - 4.0*s * t;
	out.coeffRef(1, 2) = -1.0 + 5.0*t - 4.0*t * t;

	out.coeffRef(0, 3) = s + 4.0*t * s - 4.0*s * s;
	out.coeffRef(1, 3) = -t - 4.0*s * t + 4.0*t * t;

	out.coeffRef(0, 4) = -3.0 + 7.0*s + 6.0*t - 8.0*t * s - 4.0*s * s;
	out.coeffRef(1, 4) = 5.0*t - 4.0*s * t - 8.0*t * t;

	out.coeffRef(0, 5) = -5.0*s + 4.0*t * s + 8.0*s * s;
	out.coeffRef(1, 5) = 3.0 - 6.0*s - 7.0*t + 8.0*s * t + 4.0*t * t;

	out.coeffRef(0, 6) = 16.0*s - 8.0*s * t - 16.0*s * s;
	out.coeffRef(1, 6) = 8.0*t - 16.0*s * t - 8.0*t * t;

	out.coeffRef(0, 7) = 8.0*s - 16.0*s * t - 8.0*s * s;
	out.coeffRef(1, 7) = 16.0*t - 8.0*s * t - 16.0*t * t;

};
template<unsigned QuadraturePrecision, scheme TimeScheme>
inline void solver<QuadraturePrecision, TimeScheme>::evaluate_raviartthomas_basis_divergence(Real const s, Real const t, Eigen::Matrix<Real, 8, 1> & out) {

	out.coeffRef(0) = -3.0 + 4.0*t + 8.0*s - 3.0 + 4.0*s + 8.0*t;
	out.coeffRef(1) = 5.0 - 8.0*s + 1.0 - 4.0*s;
	out.coeffRef(2) = 1.0 - 4.0*t + 5.0 - 8.0*t;
	out.coeffRef(3) = 1.0 + 4.0 * t - 8.0*s - 1.0 - 4.0*s + 8.0*t;
	out.coeffRef(4) = 7.0 - 8.0*t - 8.0*s + 5.0 - 4.0*s - 16.0*t;
	out.coeffRef(5) = -5.0 + 16.0*s + 4.0*t - 7.0 + 8.0*s + 8.0*t;
	out.coeffRef(6) = 16.0 - 8.0*t - 32.0*s + 8.0 - 16.0*s - 16.0*t;
	out.coeffRef(7) = 8.0 - 16.0*s - 16.0*t + 16.0 - 8.0*s - 32.0*t;

};


template<unsigned QuadraturePrecision, scheme TimeScheme>
inline void solver<QuadraturePrecision, TimeScheme>::evaluate_polynomial_basis(Real const s, Real const t, Eigen::Matrix<Real, 3, 1> & out) {

	out.coeffRef(0) = +1.0;
	out.coeffRef(1) = -1.0 + 2.0*s;
	out.coeffRef(2) = -1.0 + 2.0*t;

};
template<unsigned QuadraturePrecision, scheme TimeScheme>
inline void solver<QuadraturePrecision, TimeScheme>::evaluate_polynomial_basis_gradient(Real const s, Real const t, Eigen::Matrix<Real, 2, 3> & out) {

	out.coeffRef(0, 0) = 0.0;
	out.coeffRef(1, 0) = 0.0;

	out.coeffRef(0, 1) = 2.0;
	out.coeffRef(1, 1) = 0.0;

	out.coeffRef(0, 2) = 0.0;
	out.coeffRef(1, 2) = 2.0;

};


template<unsigned QuadraturePrecision, scheme TimeScheme>
inline void solver<QuadraturePrecision, TimeScheme>::evaluate_edge_polynomial_basis(Real const ksi, unsigned const El, Eigen::Matrix<Real, 2, 1> & out, Real const & orientation) {


	switch (El) {

	case 0:
		out.coeffRef(0) = 1.0;
		out.coeffRef(1) = orientation * (2.0 / sqrt(2.0)) * (ksi - sqrt(2.0) / 2.0);
		return;
	case 1:
		out.coeffRef(0) = 1.0;
		out.coeffRef(1) = orientation * 2.0 * (ksi - 0.5);
		return;
	case 2:
		out.coeffRef(0) = 1.0;
		out.coeffRef(1) = orientation * 2.0 * (ksi - 0.5);
		return;

	}

};


template<unsigned QuadraturePrecision, scheme TimeScheme>
inline void solver<QuadraturePrecision, TimeScheme>::evaluate_edge_parametrization(Real const ksi, unsigned const El, Eigen::Matrix<Real, 2, 1> & out) {


	switch (El) {

	case 0:

		out.coeffRef(0) = 1.0 - ksi / sqrt(2.0);
		out.coeffRef(1) = ksi / sqrt(2.0);
		return;

	case 1:

		out.coeffRef(0) = 0.0;
		out.coeffRef(1) = 1.0 - ksi;
		return;

	case 2:

		out.coeffRef(0) = ksi;
		out.coeffRef(1) = 0.0;
		return;


	}

};
template<unsigned QuadraturePrecision, scheme TimeScheme>
inline void solver<QuadraturePrecision, TimeScheme>::evaluate_edge_parametrization_derivative(Real const ksi, unsigned const El, Eigen::Matrix<Real, 2, 1> & out) {


	switch (El) {

	case 0:

		out.coeffRef(0) = -1.0 / sqrt(2.0);
		out.coeffRef(1) = 1.0 / sqrt(2.0);
		return;

	case 1:

		out.coeffRef(0) = 0.0;
		out.coeffRef(1) = -1.0;
		return;

	case 2:

		out.coeffRef(0) = 1.0;
		out.coeffRef(1) = 0.0;
		return;

	}

};
template<unsigned QuadraturePrecision, scheme TimeScheme>
inline void solver<QuadraturePrecision, TimeScheme>::evaluate_edge_parametrization_opposite(Real const ksi, unsigned const El, Eigen::Matrix<Real, 2, 1> & out) {


	switch (El) {

	case 0:

		out.coeffRef(0) = ksi / sqrt(2.0);
		out.coeffRef(1) = 1.0 - ksi / sqrt(2.0);
		return;

	case 1:

		out.coeffRef(0) = 0.0;
		out.coeffRef(1) = ksi;
		return;

	case 2:

		out.coeffRef(0) = 1.0 - ksi;
		out.coeffRef(1) = 0.0;
		return;

	}

};
template<unsigned QuadraturePrecision, scheme TimeScheme>
inline void solver<QuadraturePrecision, TimeScheme>::evaluate_edge_parametrization_opposite_derivative(Real const ksi, unsigned const El, Eigen::Matrix<Real, 2, 1> & out) {


	switch (El) {

	case 0:

		out.coeffRef(0) = 1.0 / sqrt(2.0);
		out.coeffRef(1) = -1.0 / sqrt(2.0);
		return;

	case 1:

		out.coeffRef(0) = 0.0;
		out.coeffRef(1) = 1.0;
		return;

	case 2:

		out.coeffRef(0) = -1.0;
		out.coeffRef(1) = 0.0;
		return;

	}

};


template<unsigned QuadraturePrecision, scheme TimeScheme>
inline void solver<QuadraturePrecision, TimeScheme>::evaluate_edge_normal(Eigen::Matrix<Real, 2, 3> & out) {

	out.coeffRef(0, 0) = 1.0 / sqrt(2.0);
	out.coeffRef(1, 0) = 1.0 / sqrt(2.0);

	out.coeffRef(0, 1) = -1.0;
	out.coeffRef(1, 1) = 0.0;

	out.coeffRef(0, 2) = 0.0;
	out.coeffRef(1, 2) = -1.0;

};


template<unsigned QuadraturePrecision, scheme TimeScheme>
inline Real solver<QuadraturePrecision, TimeScheme>::phi0(Real const s, Real const t) {

	return 1.0;

};
template<unsigned QuadraturePrecision, scheme TimeScheme>
inline Real solver<QuadraturePrecision, TimeScheme>::phi1(Real const s, Real const t) {

	return -1.0 + 2.0*s;

};
template<unsigned QuadraturePrecision, scheme TimeScheme>
inline Real solver<QuadraturePrecision, TimeScheme>::phi2(Real const s, Real const t) {

	return -1.0 + 2.0*t;

};



/*****************************************************************************/
/*                                                                           */
/*    - Integral of Source * P1 basis function phi_m						 */
/*                                                                           */
/*****************************************************************************/
template<unsigned QuadraturePrecision, scheme TimeScheme>
Real solver<QuadraturePrecision, TimeScheme>::F0(tm_pointer K, Real const time) {


	vm_pointer const a = K->vertices[0];
	vm_pointer const b = K->vertices[1];
	vm_pointer const c = K->vertices[2];

	Real const x0 = a->x;
	Real const y0 = a->y;

	Real const x1 = b->x;
	Real const y1 = b->y;

	Real const x2 = c->x;
	Real const y2 = c->y;

	Real const detJF = abs((x1 - x0)*(y2 - y0) - (x2 - x0)*(y1 - y0));

	//quadrature_triangle<Real> quad(quadrature_order);
	//unsigned const num_quad_points = quad.NumberOfPoints;

	Real integral = 0.0;

	for (unsigned i = 0; i < NumberOfQuadraturePointsTriangle; i++) {

		Real const s = QuadraturePointsAndWeightsOnReferenceTriangle(i, 0);
		Real const t = QuadraturePointsAndWeightsOnReferenceTriangle(i, 1);
		Real const w = QuadraturePointsAndWeightsOnReferenceTriangle(i, 2);

		//Real const s = quad.points_x[i];
		//Real const t = quad.points_y[i];
		//Real const w = 0.5 * quad.weights[i];

		Real const x = x0 + (x1 - x0)*s + (x2 - x0)*t;
		Real const y = y0 + (y1 - y0)*s + (y2 - y0)*t;

		integral += w * source(x, y, time) * phi0(s, t);

	}

	return 0.5 * detJF * integral;

};
template<unsigned QuadraturePrecision, scheme TimeScheme>
Real solver<QuadraturePrecision, TimeScheme>::F1(tm_pointer K, Real const time) {


	vm_pointer const a = K->vertices[0];
	vm_pointer const b = K->vertices[1];
	vm_pointer const c = K->vertices[2];

	Real const x0 = a->x;
	Real const y0 = a->y;

	Real const x1 = b->x;
	Real const y1 = b->y;

	Real const x2 = c->x;
	Real const y2 = c->y;

	Real const detJF = abs((x1 - x0)*(y2 - y0) - (x2 - x0)*(y1 - y0));

	//quadrature_triangle<Real> quad(quadrature_order);
	//unsigned const num_quad_points = quad.NumberOfPoints;

	Real integral = 0.0;

	for (unsigned i = 0; i < NumberOfQuadraturePointsTriangle; i++) {

		Real const s = QuadraturePointsAndWeightsOnReferenceTriangle(i, 0);
		Real const t = QuadraturePointsAndWeightsOnReferenceTriangle(i, 1);
		Real const w = QuadraturePointsAndWeightsOnReferenceTriangle(i, 2);

		//Real const s = quad.points_x[i];
		//Real const t = quad.points_y[i];
		//Real const w = 0.5 * quad.weights[i];

		Real const x = x0 + (x1 - x0)*s + (x2 - x0)*t;
		Real const y = y0 + (y1 - y0)*s + (y2 - y0)*t;

		integral += w * source(x, y, time) * phi1(s, t);

	}

	return 0.5 * detJF * integral;

};
template<unsigned QuadraturePrecision, scheme TimeScheme>
Real solver<QuadraturePrecision, TimeScheme>::F2(tm_pointer K, Real const time) {


	vm_pointer const a = K->vertices[0];
	vm_pointer const b = K->vertices[1];
	vm_pointer const c = K->vertices[2];

	Real const x0 = a->x;
	Real const y0 = a->y;

	Real const x1 = b->x;
	Real const y1 = b->y;

	Real const x2 = c->x;
	Real const y2 = c->y;

	Real const detJF = abs((x1 - x0)*(y2 - y0) - (x2 - x0)*(y1 - y0));

	//quadrature_triangle<Real> quad(quadrature_order);
	//unsigned const num_quad_points = quad.NumberOfPoints;

	Real integral = 0.0;

	for (unsigned i = 0; i < NumberOfQuadraturePointsTriangle; i++) {

		Real const s = QuadraturePointsAndWeightsOnReferenceTriangle(i, 0);
		Real const t = QuadraturePointsAndWeightsOnReferenceTriangle(i, 1);
		Real const w = QuadraturePointsAndWeightsOnReferenceTriangle(i, 2);

		//Real const s = quad.points_x[i];
		//Real const t = quad.points_y[i];
		//Real const w = 0.5 * quad.weights[i];

		Real const x = x0 + (x1 - x0)*s + (x2 - x0)*t;
		Real const y = y0 + (y1 - y0)*s + (y2 - y0)*t;

		integral += w * source(x, y, time) * phi2(s, t);

	}

	return 0.5 * detJF * integral;

};

template<unsigned QuadraturePrecision, scheme TimeScheme>
void solver<QuadraturePrecision, TimeScheme>::evaluate_source_integrals() {


	Real const time = nt * dt;

	for (unsigned k = 0; k < nk; k++) {


		//tm_pointer const K = Mesh->get_triangle(k);
		tm_pointer const K		 = MeshElements[k];
		unsigned const   k_index = K->index;

		vm_pointer const va = K->vertices[0];
		vm_pointer const vb = K->vertices[1];
		vm_pointer const vc = K->vertices[2];

		Real const x0 = va->x;
		Real const y0 = va->y;

		Real const x1 = vb->x;
		Real const y1 = vb->y;

		Real const x2 = vc->x;
		Real const y2 = vc->y;

		Real const x10 = x1 - x0;
		Real const x20 = x2 - x0;
		Real const y10 = y1 - y0;
		Real const y20 = y2 - y0;

		Real const HalfDetJF = 0.5 * abs(x10 * y20 - x20 * y10);


		Real Integral0 = 0.0;
		Real Integral1 = 0.0;
		Real Integral2 = 0.0;

		for (unsigned n = 0; n < NumberOfQuadraturePointsTriangle; n++) {


			Real const s = QuadraturePointsAndWeightsOnReferenceTriangle(n, 0);
			Real const t = QuadraturePointsAndWeightsOnReferenceTriangle(n, 1);
			Real const w = QuadraturePointsAndWeightsOnReferenceTriangle(n, 2);

			Real const x = x0 + (x1 - x0)*s + (x2 - x0)*t;
			Real const y = y0 + (y1 - y0)*s + (y2 - y0)*t;


			Real const SourceValueTimesWeight = w * source(x, y, time);

			Integral0 += SourceValueTimesWeight * phi0(s, t);
			Integral1 += SourceValueTimesWeight * phi1(s, t);
			Integral2 += SourceValueTimesWeight * phi2(s, t);

		}

		Sources.setCoeff(k_index, 0) = HalfDetJF * Integral0;
		Sources.setCoeff(k_index, 1) = HalfDetJF * Integral1;
		Sources.setCoeff(k_index, 2) = HalfDetJF * Integral2;

	}

};

