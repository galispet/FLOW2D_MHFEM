#pragma once



#include "coefficient_matrix.h"
#include "integration.h"
#include "mesh.h"

#include <omp.h>
#include <Eigen/Sparse>
#include <Eigen/Dense>




typedef Eigen::SparseMatrix<real>			SparseMatrix;
typedef Eigen::VectorXd						DenseVector;

enum TimeDiscretization { CrankNicolson, EulerBackward };



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



template <TimeDiscretization TimeScheme, unsigned QuadraturePrecision>
class Solver {



	double const TimeSchemeParameter			= (TimeScheme == TimeDiscretization::CrankNicolson) ? 0.5 : 1.0;
	double const NumericalZeroBound				= 1e-11;
	double const NonlinearSolverPrecision		= DBL_EPSILON;
	unsigned const NonlinearSolverMaxIterations = 100;



public:

	Solver(Mesh & m, unsigned nt_0, double dt_0);
	~Solver();

	void getSolution();

	void setTimeStep(double _dt) {	this->dt = _dt; };
	void setTimeLevel(int _nt) {	this->nt = _nt; };

	void exportPressures(std::string const & fileName);
	void exportConcentrations(std::string const & fileName);
	void exportVelocityField(std::string const & fileName);
	void computeError(std::string const & fileName);



private:




	// Triangulation mesh
	Mesh * mesh;

	// nk = number of elements in the mesh, ne = number of edges in the mesh
	unsigned nk;
	unsigned ne;


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

	CoeffMatrix2D<3, 3>		Alpha;
	CoeffMatrix1D<3>		Lambda;
	double *				Sigma;


	static unsigned const NumberOfQuadraturePointsEdge		= get_number_of_quadrature_points_edge<QuadraturePrecision>();
	static unsigned const NumberOfQuadraturePointsTriangle	= get_number_of_quadrature_points_triangle<QuadraturePrecision>();


	double * viscosities;
	double * porosities;

	double * thetas;
	double * thetas_prev;

	// Trace pressure system matrix
	SparseMatrix	R;
	SparseMatrix	M;
	DenseVector		V;

	// Pressure system matrix
	SparseMatrix	iD;
	SparseMatrix	H;
	DenseVector		G;



	Eigen::BiCGSTAB<SparseMatrix, Eigen::DiagonalPreconditioner<double>>	BiConjugateGradientSolver;
	Eigen::SparseLU<SparseMatrix, Eigen::COLAMDOrdering<int>>				sparseLUsolver_TracePressureSystem;



	int		nt;
	double	dt;


	void initializeValues();
	void computeThetas();

	bool stopCriterion();
	void concentrationCorrection();


	void computePressureEquation();

	void computeTracePressures();
	void computeVelocities();
	void updateConcentrations();

	double upwindConcentration(t_pointer const & K, unsigned const El);


	void assembleR();
	void assembleM();
	void assembleV();

	void assembleInverseD();
	void assembleH();
	void assembleG();

	void assemble_Alpha();
	void assemble_Sigma();
	void assemble_Lambda();

	


	/*****************************************************************************/
	/*                                                                           */
	/*    - Miscellaneous									   			         */
	/*                                                                           */
	/*****************************************************************************/
	unsigned LocalIndex(t_pointer const & K, e_pointer const & E) {

		return K->get_edge_index(E);

	};

	void evaluateBasisRaviartThomas(double const s, double const t, Eigen::Matrix<double, 2, 3> & out);

	static double Barenblatt(double x, double y, double time) {

		double const norm_squared = x * x + y * y;

		//return (1.0 / pow(time, (1.0 / 3.0)))*fmax(1.0 - norm_squared / (12.0 * pow(time, (2.0 / 3.0))), 0.0);

		return (1.0 / sqrt(time))*fmax(1.0 - norm_squared / (16.0 * sqrt(time)), 0.0);

	};

	double integrate_edge(e_pointer const e, double  time, double (*fun)(double , double , double )) {


		double const x0 = e->a->x;
		double const y0 = e->a->y;

		double const x1 = e->b->x;
		double const y1 = e->b->y;

		double const length = e->length();


		// Quadrature weights and points on reference segment [-1,1]
		gauss_quadrature_1D quad(QuadraturePrecision);
		unsigned const numberOfQuadraturePoints = quad.NumberOfPoints;


		double integral = 0.0;

		for (unsigned n = 0; n < numberOfQuadraturePoints; n++) {


			// Quadrature point on reference segment [-1,1] and its weights
			double const s = quad.points[n];
			double const w = quad.weights[n];

			double const x = x0 + 0.5*(1.0 + s)*(x1 - x0);
			double const y = y0 + 0.5*(1.0 + s)*(y1 - y0);

			integral += w * fun(x, y, time);

		}

		return 0.5 * length * integral;

	};
	double integrate_edge(e_pointer const e, double(*fun)(double , double  )) {


		double const x0 = e->a->x;
		double const y0 = e->a->y;
		double const x1 = e->b->x;
		double const y1 = e->b->y;

		double const length = e->length();



		// Quadrature weights and points on reference segment [-1,1]
		gauss_quadrature_1D quad(QuadraturePrecision);
		unsigned const numberOfQuadraturePoints = quad.NumberOfPoints;


		double integral = 0.0;

		for (unsigned n = 0; n < numberOfQuadraturePoints; n++) {


			// Quadrature point on reference segment [-1,1] and its weights
			double const s = quad.points[n];
			double const w = quad.weights[n];

			double const x = x0 + 0.5*(1.0 + s)*(x1 - x0);
			double const y = y0 + 0.5*(1.0 + s)*(y1 - y0);

			integral += w * fun(x, y);

		}


		return 0.5 * length * integral;

	};

	double integrate_triangle(t_pointer const K, double const time, double(*fun)(double, double, double)) {


		v_pointer const v0 = K->vertices[0];
		v_pointer const v1 = K->vertices[1];
		v_pointer const v2 = K->vertices[2];

		double const x0 = v0->x;
		double const y0 = v0->y;
		double const x1 = v1->x;
		double const y1 = v1->y;
		double const x2 = v2->x;
		double const y2 = v2->y;

		double const detJF = abs((x1 - x0)*(y2 - y0) - (x2 - x0)*(y1 - y0));


		// Quadrature weights and points on [-1,1]
		quadrature_triangle quad(QuadraturePrecision);
		unsigned const numberOfQuadraturePoints = quad.NumberOfPoints;


		double integral = 0.0;

		for (unsigned n = 0; n < numberOfQuadraturePoints; n++) {


			double const s = quad.points_x[n];
			double const t = quad.points_y[n];
			double const w = quad.weights[n];

			// Gauss points on given triangle. This is the transformation : [-1,1] x [-1,1] square --> reference triangle ([0,0],[1,0],[0,1]) --> given triangle ([x0y,0],[x1,y1],[x2,y2]) 
			double const x = x0 + (x1 - x0)*s + (x2 - x0)*t;
			double const y = y0 + (y1 - y0)*s + (y2 - y0)*t;

			integral += w * fun(x, y, time);

		}

		// Not quite sure, why there is the 1/2. See https://people.sc.fsu.edu/~jburkardt/datasets/quadrature_rules_tri/quadrature_rules_tri.html
		return 0.5 * detJF * integral;

	};
	double integrate_triangle(t_pointer const K, double (*fun)(double , double )) {


		v_pointer const v0 = K->vertices[0];
		v_pointer const v1 = K->vertices[1];
		v_pointer const v2 = K->vertices[2];

		double const x0 = v0->x;
		double const y0 = v0->y;
		double const x1 = v1->x;
		double const y1 = v1->y;
		double const x2 = v2->x;
		double const y2 = v2->y;

		double const detJF = abs((x1 - x0)*(y2 - y0) - (x2 - x0)*(y1 - y0));


		quadrature_triangle quad(QuadraturePrecision);
		unsigned const numberOfQuadraturePoints = quad.NumberOfPoints;


		double integral = 0.0;

		for (unsigned n = 0; n < numberOfQuadraturePoints; n++) {


			double const s = quad.points_x[n];
			double const t = quad.points_y[n];
			double const w = quad.weights[n];

			// Gauss points on given triangle. This is the transformation : [-1,1] x [-1,1] square --> reference triangle ([0,0],[1,0],[0,1]) --> given triangle ([x0y,0],[x1,y1],[x2,y2]) 
			double const x = x0 + (x1 - x0)*s + (x2 - x0)*t;
			double const y = y0 + (y1 - y0)*s + (y2 - y0)*t;

			integral += w * fun(x, y);

		}

		// Not quite sure, why there is the 1/2. See https://people.sc.fsu.edu/~jburkardt/datasets/quadrature_rules_tri/quadrature_rules_tri.html
		return 0.5 * detJF * integral;

	};

	
	static void Permeability(double const x, double const y, Eigen::Matrix<double, 2, 2> & out) {

		out.coeffRef(0, 0) = 1.0;
		out.coeffRef(0, 1) = 0.0;
		out.coeffRef(1, 0) = 0.0;
		out.coeffRef(1, 1) = 1.0;

	};
	static double Porosity(double const x, double const y) {

		return 1.0;

	};
	static double Viscosity(double const x, double const y) {

		return 0.5;

	};
	static double Source(double const x, double const y, double const t) {

		return 0.0;

	};


	double NEUMANN_GAMMA_Q_velocity(e_pointer const E, double const time) {

		return 0.0;

	};
	double DIRICHLET_GAMMA_Q_concentration(e_pointer const E, double const time) {

		return 0.0;

	};
	double DIRICHLET_GAMMA_P_concentration(e_pointer const E, double const time) {

		return integrate_edge(E, time, barenblatt) / E->length();

	};
	double DIRICHLET_GAMMA_P_pressure(e_pointer const E, double const time) {

		return integrate_edge(E, time, barenblatt) / E->length();

	};
	
	double DIRICHLET_GAMMA_Q_concentration(double const s, double const t, double const time) {

		return 0.0 * barenblatt(s, t, time);

	};
	double DIRICHLET_GAMMA_P_concentration(double const s, double const t, double const time) {

		return barenblatt(s, t, time);

	};
	double DIRICHLET_GAMMA_P_pressure(double const s, double const t, double const time) {

		return barenblatt(s, t, time);

	};

};




template <TimeDiscretization TimeScheme, unsigned QuadraturePrecision>
Solver<TimeScheme, QuadraturePrecision>::Solver(Mesh & m, unsigned nt_0, double dt_0)

	: nk(m.get_number_of_triangles()), ne(m.get_number_of_edges())

{


	mesh = & m;

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
	thetas		= new double[nk];
	thetas_prev	= new double[nk];


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
	v	.setZero();

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
	Alpha	.setNumberOfElements(nk);
	Lambda	.setNumberOfElements(nk);
	Sigma = new double[nk];

	Alpha.setZero();
	Lambda.setZero();

	/*****************************************************************************/
	/*                                                                           */
	/*    - Memory allocation for the pressure system					         */
	/*                                                                           */
	/*****************************************************************************/
	iD	.resize(nk, nk);
	H	.resize(nk, ne);
	G	.resize(nk);

	M.resize(ne, ne);
	R.resize(ne, nk);
	V.resize(ne);


	/*****************************************************************************/
	/*                                                                           */
	/*    - Initilize initial pressure,concentration condition				     */
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
	/*    - Get sparsity pattern for the pressure sytem. BiCG is faster then     */
	/*      LU decomposition                                                     */
	/*                                                                           */
	/*****************************************************************************/
	assemble_Sigma();
	assemble_Lambda();

	assembleInverseD();
	assembleH();

	BiConjugateGradientSolver.analyzePattern(R * iD * H + M);


};
template <TimeDiscretization TimeScheme, unsigned QuadraturePrecision>
Solver<TimeScheme, QuadraturePrecision>::~Solver() {


	delete[] Sigma;

	delete[] viscosities;
	delete[] porosities;

	delete[] thetas;
	delete[] thetas_prev;

};


template <TimeDiscretization TimeScheme, unsigned QuadraturePrecision>
void Solver<TimeScheme, QuadraturePrecision>::initializeValues() {


	double const time = nt * dt;


	for (unsigned k = 0; k < nk; k++) {


		t_pointer const K		= mesh->get_triangle(k);
		unsigned const k_index	= K->index;
		double const Area		= K->area();


		p[k_index] = integrate_triangle(K, time, &Solver::Barenblatt) / Area;
		c[k_index] = integrate_triangle(K, time, &Solver::Barenblatt) / Area;

		/*****************************************************************************/
		/*                                                                           */
		/*    - Mean values of the viscosity and porosity on each element		     */
		/*                                                                           */
		/*****************************************************************************/
		viscosities[k_index]	= integrate_triangle(K, &Solver::Viscosity) / Area;
		porosities[k_index]		= integrate_triangle(K, &Solver::Porosity) / Area;


		/*****************************************************************************/
		/*                                                                           */
		/*    - Auxilary variables used in temporal discretization				     */
		/*                                                                           */
		/*****************************************************************************/
		rkFc[k_index]	= 0.0;
		rkFc_n[k_index] = rkFc[k_index];

		rkFp[k_index]	= 0.0;
		rkFp_n[k_index] = rkFp[k_index];

	}


};
template <TimeDiscretization TimeScheme, unsigned QuadraturePrecision>
void Solver<TimeScheme, QuadraturePrecision>::computeThetas() {

	for (unsigned k = 0; k < nk; k++)
		thetas[k] = 1.0;

};



template <TimeDiscretization TimeScheme, unsigned QuadraturePrecision>
bool Solver<TimeScheme, QuadraturePrecision>::stopCriterion() {


	double ErrorPressure	= 0.0;
	double NormPressure		= 0.0;

	double ErrorConcentration	= 0.0;
	double NormConcentration	= 0.0;

	double ErrorBeta	= 0.0;
	double NormBeta		= 0.0;

	for (unsigned k = 0; k < nk; k++) {


		t_pointer const	K		= mesh->get_triangle(k);
		unsigned const k_index	= K->index;
		double const Area		= K->area();


		double const IntegralError	= sqr(p[k_index] - p_prev[k_index]);
		double const IntegralNorm	= sqr(p[k_index]);

		ErrorPressure	+= Area * IntegralError;
		NormPressure	+= Area * IntegralNorm;

	}

	double const eP = ErrorPressure / NormPressure;

	if (eP > NonlinearSolverPrecision)
		return false;


	for (unsigned k = 0; k < nk; k++) {


		t_pointer const	K		= mesh->get_triangle(k);
		unsigned const k_index	= K->index;
		double const Area		= K->area();

		double const IntegralError	= sqr(c[k_index] - c_prev[k_index]);
		double const IntegralNorm	= sqr(c[k_index]);

		ErrorConcentration	+= Area * IntegralError;
		NormConcentration	+= Area * IntegralNorm;

	}

	double const eC = ErrorConcentration / NormConcentration;

	if (eC > NonlinearSolverPrecision)
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
template<TimeDiscretization TimeScheme, unsigned QuadraturePrecision>
void Solver<TimeScheme, QuadraturePrecision>::concentrationCorrection() {


};



template <TimeDiscretization TimeScheme, unsigned QuadraturePrecision>
void Solver<TimeScheme, QuadraturePrecision>::computeTracePressures() {

	tp = sparseLUsolver_TracePressureSystem.solve(R * p - V);

};
template <TimeDiscretization TimeScheme, unsigned QuadraturePrecision>
void Solver<TimeScheme, QuadraturePrecision>::computePressureEquation() {


	BiConjugateGradientSolver.factorize(R * iD * H + M);
	tp = BiConjugateGradientSolver.solveWithGuess(R * iD * G - V, tp);

	/*Eigen::SparseLU<SparseMatrix> const solver(R * iD * H + M);
	tp = solver.solve(R * iD * G - V);*/

	p = iD * (G - H * tp);

};



template <TimeDiscretization TimeScheme, unsigned QuadraturePrecision>
void Solver<TimeScheme, QuadraturePrecision>::computeVelocities() {


	double const time = (nt + 1)* dt;

	for (unsigned k = 0; k < nk; k++) {


		t_pointer const K		= mesh->get_triangle(k);
		unsigned const k_index	= K->index;


		for (unsigned El = 0; El < 3; El++) {


			e_pointer const E = K->edges[El];

			if (E->marker == E_MARKER::NEUMANN) {

				v.setCoeff(k_index, El) = NEUMANN_GAMMA_Q_velocity(E, time);

				continue;

			}


			double AlphaSum	= 0.0;
			double AlphaTP	= 0.0;

			for (unsigned j = 0; j < 3; j++)
				AlphaSum += Alpha(k_index, El, j);

			for (unsigned l = 0; l < 3; l++)
				AlphaTP += Alpha(k_index, El, l) * tp(K->edges[l]->index);


			v.setCoeff(k_index, El) = (AlphaSum * p[k_index] - AlphaTP) / viscosities[k_index];

		}
	}

};
template <TimeDiscretization TimeScheme, unsigned QuadraturePrecision>
void Solver<TimeScheme, QuadraturePrecision>::updateConcentrations() {



	double const TimeCoefficient1 = dt * TimeSchemeParameter;
	double const TimeCoefficient2 = dt * (1.0 - TimeSchemeParameter);

	for (unsigned k = 0; k < nk; k++) {


		t_pointer const K		= mesh->get_triangle(k);
		unsigned const k_index	= K->index;
		double const Area		= K->area();


		double Value = 0.0;

		for (unsigned El = 0; El < 3; El++)
			Value += v(k_index, El) * upwindConcentration(K, El);


		rkFc[k_index] = -Value / (porosities[k_index] * Area);

		c[k_index] = c_n[k_index] + TimeCoefficient1 * rkFc[k_index] + TimeCoefficient2 * rkFc_n[k_index];

	}


};
template <TimeDiscretization TimeScheme, unsigned QuadraturePrecision>
double Solver<TimeScheme, QuadraturePrecision>::upwindConcentration(t_pointer const & K, unsigned const El) {



	double const time = (nt + 1) * dt;


	unsigned const k_index	= K->index;
	e_pointer const E		= K->edges[El];
	E_MARKER const e_marker = E->marker;


	double VelocityDotNormal	= 0.0;
	double Concentration		= 0.0;

	if (e_marker == E_MARKER::NEUMANN) {


		VelocityDotNormal = NEUMANN_GAMMA_Q_velocity(E, time);

		if (VelocityDotNormal < 0.0)
			return DIRICHLET_GAMMA_Q_concentration(E, time);


	}
	else
		VelocityDotNormal = v(k_index, El);
	


	if (VelocityDotNormal >= 0.0)
		Concentration = c_prev[k_index];
	else {


		if (E->marker == E_MARKER::DIRICHLET)
			return DIRICHLET_GAMMA_P_concentration(E, time);


		unsigned const kn_index = K->neighbors[El]->index;

		Concentration = c_prev[kn_index];
		
	}

	return Concentration;

};



template <TimeDiscretization TimeScheme, unsigned QuadraturePrecision>
void Solver<TimeScheme, QuadraturePrecision>::assembleR() {


	std::vector<Eigen::Triplet<double>> TripletsR;

	TripletsR.reserve(3 * nk);


	for (unsigned e = 0; e < ne; e++) {


		e_pointer const E		= mesh->get_edge(e);
		unsigned const e_index	= E->index;
		E_MARKER const e_marker = E->marker;


		if (e_marker == E_MARKER::DIRICHLET)
			continue;


		for (unsigned neighbor = 0; neighbor < 2; neighbor++) {


			t_pointer const K = E->neighbors[neighbor];

			if (!K)
				continue;


			unsigned const dof		= LocalIndex(K, E);
			unsigned const k_index	= K->index;

			double AlphaSum = 0.0;

			for (unsigned j = 0; j < 3; j++)
				AlphaSum += Alpha(k_index, dof, j);


			double const Value			= AlphaSum / viscosities[k_index];
			double const dumpedValue	= fabs(Value) < NumericalZeroBound ? 0.0 : Value;

			Eigen::Triplet<double> const T(e_index, k_index, dumpedValue);

			TripletsR.push_back(T);

		}
	}

	R.setFromTriplets(TripletsR.begin(), TripletsR.end());

};
template <TimeDiscretization TimeScheme, unsigned QuadraturePrecision>
void Solver<TimeScheme, QuadraturePrecision>::assembleM() {


	std::vector<Eigen::Triplet<double>> TripletsM;

	TripletsM.reserve(5 * ne);


	for (unsigned e = 0; e < ne; e++) {


		e_pointer const E		= mesh->get_edge(e);
		unsigned const e_index	= E->index;
		E_MARKER const e_marker = E->marker;


		if (e_marker == E_MARKER::DIRICHLET) {


			Eigen::Triplet<double> const T(e_index, e_index, -1.0);
			TripletsM.push_back(T);

			continue;

		}

		for (unsigned neighbor = 0; neighbor < 2; neighbor++) {


			t_pointer const K = E->neighbors[neighbor];

			if (!K)
				continue;


			unsigned const dof		= LocalIndex(K, E);
			unsigned const k_index	= K->index;


			for (unsigned El = 0; El < 3; El++) {


				e_pointer const E_local				= K->edges[El];
				unsigned const e_local_index_global = E_local->index;


				double const Value			= Alpha(k_index, dof, El) / viscosities[k_index];
				double const dumpedValue	= std::fabs(Value) < NumericalZeroBound ? 0.0 : Value;


				Eigen::Triplet<double> const T(e_index, e_local_index_global, dumpedValue);
				TripletsM.push_back(T);

			}
		}
	}


	TripletsM.shrink_to_fit();

	M.setFromTriplets(TripletsM.begin(), TripletsM.end());


	/*****************************************************************************/
	/*                                                                           */
	/*    - Compute the LU decomposition of the matrix M						 */
	/*                                                                           */
	/*****************************************************************************/	
	sparseLUsolver_TracePressureSystem.compute(M);

};
template <TimeDiscretization TimeScheme, unsigned QuadraturePrecision>
void Solver<TimeScheme, QuadraturePrecision>::assembleV() {


	double const time = (nt + 1) * dt;

	for (unsigned e = 0; e < ne; e++) {


		e_pointer const E		= mesh->get_edge(e);
		unsigned const e_index	= E->index;
		E_MARKER const marker	= E->marker;


		V[e_index] = 0.0;

		if (marker == E_MARKER::NEUMANN)
			V[e_index] = NEUMANN_GAMMA_Q_velocity(E, time);
		else if (marker == E_MARKER::DIRICHLET)
			V[e_index] = DIRICHLET_GAMMA_P_pressure(E, time);

	}

};


template <TimeDiscretization TimeScheme, unsigned QuadraturePrecision>
void Solver<TimeScheme, QuadraturePrecision>::assembleInverseD() {


	double const TimeCoefficient = TimeSchemeParameter * dt;

	std::vector<Eigen::Triplet<double>> TripletsInverseD;
	TripletsInverseD.reserve(nk);


	for (unsigned k = 0; k < nk; k++) {


		t_pointer const K		= mesh->get_triangle(k);
		unsigned const k_index	= K->index;


		double const Value = 1.0 / (1.0 - TimeCoefficient * Sigma[k_index]);

		Eigen::Triplet<double> const T(k_index, k_index, Value);

		TripletsInverseD.push_back(T);
	}

	TripletsInverseD.shrink_to_fit();

	iD.setFromTriplets(TripletsInverseD.begin(), TripletsInverseD.end());

};
template <TimeDiscretization TimeScheme, unsigned QuadraturePrecision>
void Solver<TimeScheme, QuadraturePrecision>::assembleH() {


	double const TimeCoefficient = -TimeSchemeParameter * dt;

	std::vector<Eigen::Triplet<double>> TripletsH;	
	TripletsH.reserve(3 * nk);


	for (unsigned k = 0; k < nk; k++) {


		t_pointer const K		= mesh->get_triangle(k);
		unsigned const k_index	= K->index;


		for (unsigned El = 0; El < 3; El++) {

			e_pointer const E		= K->edges[El];
			unsigned const e_index	= E->index;


			double const Value = TimeCoefficient * Lambda(k_index, El);

			Eigen::Triplet<double> const T(k_index, e_index, Value);

			TripletsH.push_back(T);

		}
	}

	TripletsH.shrink_to_fit();

	H.setFromTriplets(TripletsH.begin(), TripletsH.end());

};
template <TimeDiscretization TimeScheme, unsigned QuadraturePrecision>
void Solver<TimeScheme, QuadraturePrecision>::assembleG() {


	double const TimeCoefficient = dt * (1.0 - TimeSchemeParameter);

	for (unsigned k = 0; k < nk; k++) {

		unsigned const k_index = mesh->get_triangle(k)->index;

		G[k_index] = p_n[k_index] + TimeCoefficient * rkFp_n[k_index];

	}


};


template <TimeDiscretization TimeScheme, unsigned QuadraturePrecision>
void Solver<TimeScheme, QuadraturePrecision>::assemble_Alpha() {


	quadrature_triangle const QuadratureOnTriangle(QuadraturePrecision);

	Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>	Integral(3, 3);
	Eigen::Matrix<double, 2, 3>								BasisRaviartThomas;


	for (int k = 0; k < nk; k++) {



		t_pointer const K		= mesh->get_triangle(k);
		unsigned const k_index	= K->index;

		v_pointer const v0 = K->vertices[0];
		v_pointer const v1 = K->vertices[1];
		v_pointer const v2 = K->vertices[2];

		double const x0		= v0->x;
		double const y0		= v0->y;
		double const x10	= v1->x - x0;
		double const y10	= v1->y - y0;
		double const x20	= v2->x - x0;
		double const y20	= v2->y - y0;

		double const detJF = fabs(x10 * y20 - x20 * y10);


		Eigen::Matrix<double, 2, 2> JF;

		JF.coeffRef(0, 0) = x10;
		JF.coeffRef(0, 1) = x20;
		JF.coeffRef(1, 0) = y10;
		JF.coeffRef(1, 1) = y20;

		Integral.setZero();


		for (unsigned n = 0; n < NumberOfQuadraturePointsTriangle; n++) {


			double const s =	   QuadratureOnTriangle.points_x[n];
			double const t =	   QuadratureOnTriangle.points_y[n];
			double const w = 0.5 * QuadratureOnTriangle.weights[n];

			double const x = x0 + (x10)*s + (x20)*t;
			double const y = y0 + (y10)*s + (y20)*t;


			Eigen::Matrix2d iK(2, 2);
			Permeability(x, y, iK);
			iK = iK.inverse();

			evaluateBasisRaviartThomas(s, t, BasisRaviartThomas);


			for (unsigned i = 0; i < 3; i++) {
				
				Eigen::Matrix<double, 2, 1> const JFWi = JF * (BasisRaviartThomas.col(i));

				for (unsigned j = 0; j < 3; j++) {

					Eigen::Matrix<double, 2, 1> const JFWj = JF * (BasisRaviartThomas.col(j));

					Integral.coeffRef(i, j) += w * JFWi.dot(iK * JFWj);

				}
			}

		}

		Integral = Integral / detJF;


		for (unsigned i = 0; i < 3; i++)
			for (unsigned j = 0; j < 3; j++)
				Integral.coeffRef(i, j) = abs(Integral(i, j)) < NumericalZeroBound ? 0.0 : Integral(i, j);

		Integral = Integral.inverse();

		for (unsigned i = 0; i < 3; i++)
			for (unsigned j = 0; j < 3; j++)
				Alpha.setCoeff(k_index, i, j) = abs(Integral(i, j)) < NumericalZeroBound ? 0.0 : Integral(i, j);
		
	}

};
template <TimeDiscretization TimeScheme, unsigned QuadraturePrecision>
void Solver<TimeScheme, QuadraturePrecision>::assemble_Sigma() {


	for (int k = 0; k < nk; k++) {


		t_pointer const K		= mesh->get_triangle(k);
		unsigned const k_index	= K->index;
		double const Area		= K->area();

		double const Coefficient = -thetas[k_index] / (porosities[k_index] * viscosities[k_index] * Area);
		

		double Value = 0.0;

		for (unsigned El = 0; El < 3; El++) {

			// Loop unroll for the degrees of freedom (in the RT0: degree of freedom = edge)
			double const AlphaSum = Alpha(k_index, El, 0) +
									Alpha(k_index, El, 1) +
									Alpha(k_index, El, 2);

			double const tC	= upwindConcentration(K, El);

			Value += tC * AlphaSum;

		}

		Sigma[k_index] = Coefficient * Value;

	}

};
template <TimeDiscretization TimeScheme, unsigned QuadraturePrecision>
void Solver<TimeScheme, QuadraturePrecision>::assemble_Lambda() {


	for (unsigned k = 0; k < nk; k++) {

		t_pointer const K		= mesh->get_triangle(k);
		unsigned const k_index	= K->index;
		double const Area		= K->area();

		double const ElementCoeff = thetas[k_index] / (porosities[k_index] * viscosities[k_index] * Area);


		for (unsigned l = 0; l < 3; l++) {


			// Loop unroll for the edges
			double const Value = upwindConcentration(K, 0) * Alpha(k_index, 0, l) +
								 upwindConcentration(K, 1) * Alpha(k_index, 1, l) +
								 upwindConcentration(K, 2) * Alpha(k_index, 2, l);

			Lambda.setCoeff(k_index, l) = ElementCoeff * Value;

		}
	}

};



template <TimeDiscretization TimeScheme, unsigned QuadraturePrecision>
void Solver<TimeScheme, QuadraturePrecision>::getSolution() {


	/*****************************************************************************/
	/*                                                                           */
	/*    - Compute initial trace pressures from known inner pressures and	     */
	/*      velocity field and coefficients Thetas							     */
	/*                                                                           */
	/*    - This is initializing step. This is first guess for the               */
	/*      trace pressures and velocities on the time level nt = (n + 1)        */
	/*                                                                           */
	/*****************************************************************************/
	assembleV();

	computeTracePressures();
	computeVelocities();
	computeThetas();


	/*****************************************************************************/
	/*                                                                           */
	/*    - Set unknowns for inner iterations. Set iteration number: l = 0	     */
	/*                                                                           */
	/*          : _n	= previous n-th time level                               */
	/*          : _prev = this (n+1)-th time level, but previous iteration (l-1) */
	/*                                                                           */
	/*			: p_n	     - pressures on the n-th time level				 	 */
	/*			: p_prev     - pressures on the (n+1),(l-1)-th time level		 */
	/*			: c_n	     - concentrations on the n-th time level		 	 */
	/*			: c_prev     - concentrations on the (n+1),(l-1)-th time level	 */
	/*			: thetas      - Thetas on the n-th time level		 			 */
	/*			: thetas_prev - Thetas on the (n+1),(l-1)-th time level			 */
	/*                                                                           */
	/*****************************************************************************/
	p_n		= p;
	p_prev	= p;

	c_n		= c;
	c_prev	= c;

	for (unsigned i = 0; i < nk; i++)
		thetas_prev[i] = thetas[i];


	/*****************************************************************************/
	/*                                                                           */
	/*    - Right-hand sides of the equation y'(t) = F(x,t,y(t))			     */
	/*      arrising from Runge-Kutta schemes									 */
	/*                                                                           */
	/*			: rkFc_n - RHS for concentrations on the n-th time level	     */
	/*			: rkFc   - RHS for concentrations on (n+1),(l-1)-th time level	 */
	/*			: rkFp_n - RHS for pressures on the n-th time level				 */
	/*                                                                           */
	/*****************************************************************************/
	rkFc_n = rkFc;


	unsigned counter = 0;

	while (counter < NonlinearSolverMaxIterations) {

		
		/*****************************************************************************/
		/*                                                                           */
		/*    - Update axuilary variables which show up in matrices					 */
		/*                                                                           */
		/*****************************************************************************/
		assemble_Sigma();
		assemble_Lambda();


		/*****************************************************************************/
		/*                                                                           */
		/*    - Assembly of the matrices											 */
		/*                                                                           */
		/*****************************************************************************/
		assembleInverseD();
		assembleH();


		/*****************************************************************************/
		/*                                                                           */
		/*    - Assembly of the right-hand sides									 */
		/*                                                                           */
		/*****************************************************************************/
		assembleG();
		assembleV();


		computePressureEquation();
		computeVelocities();
		updateConcentrations();

		//concentrationCorrection();

		computeThetas();


		if (stopCriterion())
			break;


		/*****************************************************************************/
		/*                                                                           */
		/*    - Next iteration number: l = l + 1									 */
		/*                                                                           */
		/*****************************************************************************/
		p_prev = p;
		c_prev = c;


		for (unsigned i = 0; i < nk; i++)
			thetas_prev[i] = thetas[i];


		counter++;

	}


	/*****************************************************************************/
	/*                                                                           */
	/*    - The rest is needed for the next time level as the part of the		 */
	/*      Right-Hand Side of the pressure system - Euler / Crank-Nicolson      */
	/*                                                                           */
	/*****************************************************************************/
	c_prev = c;


	/*****************************************************************************/
	/*                                                                           */
	/*    - These are dependent on the concentration, therefore they must be	 */
	/*      updated																 */
	/*                                                                           */
	/*****************************************************************************/
	assemble_Lambda();
	assemble_Sigma();

	
	for (unsigned k = 0; k < nk; k++) {


		t_pointer const K		= mesh->get_triangle(k);
		unsigned const k_index	= K->index;


		double const val1	= Sigma[k_index] * p[k_index];
		double val2			= 0.0;

		for (unsigned El = 0; El < 3; El++)
			val2 += Lambda(k_index, El) * tp(K->edges[El]->index);

		/*****************************************************************************/
		/*                                                                           */
		/*    - Here we update the part on the n-th time level of the right-hand side*/
		/*      of the discretized pressure equation. Not 'rkFp' which is on		 */
		/*      the (n+1)-th time level !                                            */
		/*                                                                           */
		/*****************************************************************************/
		rkFp_n[k_index] = val1 + val2;

	}

	std::cout << nt << " - Iterations : " << counter << std::endl;

};




template <TimeDiscretization TimeScheme, unsigned QuadraturePrecision>
void Solver<TimeScheme, QuadraturePrecision>::exportPressures(std::string const & fileName) {


	std::ofstream OFSTxtFile(fileName);

	for (unsigned k = 0; k < nk; k++) {


		t_pointer const	K		= mesh->get_triangle(k);
		unsigned const k_index	= K->index;

		v_pointer const v0 = K->vertices[0];
		v_pointer const v1 = K->vertices[1];
		v_pointer const v2 = K->vertices[2];

		double const x0 = v0->x;
		double const y0 = v0->y;
		double const x1 = v1->x;
		double const y1 = v1->y;
		double const x2 = v2->x;
		double const y2 = v2->y;

		double const x[3] = { x0, x1, x2 };
		double const y[3] = { y0, y1, y2 };


		double const Value = p[k_index];

		for (unsigned i = 0; i < 3; i++)
			OFSTxtFile << std::setprecision(20) << x[i] << "\t" << y[i] << "\t" << Value << std::endl;


		OFSTxtFile << std::endl;

	}


	OFSTxtFile.close();

};
template <TimeDiscretization TimeScheme, unsigned QuadraturePrecision>
void Solver<TimeScheme, QuadraturePrecision>::exportConcentrations(std::string const & fileName) {


	std::ofstream OFSTxtFile(fileName);

	for (unsigned k = 0; k < nk; k++) {


		t_pointer const	K		= mesh->get_triangle(k);
		unsigned const k_index	= K->index;

		v_pointer const v0 = K->vertices[0];
		v_pointer const v1 = K->vertices[1];
		v_pointer const v2 = K->vertices[2];

		double const x0 = v0->x;
		double const y0 = v0->y;

		double const x1 = v1->x;
		double const y1 = v1->y;

		double const x2 = v2->x;
		double const y2 = v2->y;

		double const x[3] = { x0, x1, x2 };
		double const y[3] = { y0, y1, y2 };


		double const Value = c[k_index];

		for (unsigned i = 0; i < 3; i++)
			OFSTxtFile << std::setprecision(20) << x[i] << "\t" << y[i] << "\t" << Value << std::endl;


		OFSTxtFile << std::endl;

	}


	OFSTxtFile.close();

};
template <TimeDiscretization TimeScheme, unsigned QuadraturePrecision>
void Solver<TimeScheme, QuadraturePrecision>::exportVelocityField(std::string const & fileName) {



	quadrature_triangle const QuadratureOnTriangle(3);
	unsigned const NumberOfQuadraturePoints = QuadratureOnTriangle.NumberOfPoints;

	Eigen::MatrixXd BasisRaviartThomas(2, 8);
	Eigen::MatrixXd JF(2, 2);


	std::ofstream OFSTxtFile(fileName);


	for (unsigned k = 0; k < nk; k++) {


		t_pointer const	K		= mesh->get_triangle(k);
		unsigned const k_index	= K->index;

		v_pointer const v0 = K->vertices[0];
		v_pointer const v1 = K->vertices[1];
		v_pointer const v2 = K->vertices[2];

		double const x0 = v0->x;
		double const y0 = v0->y;
		double const x1 = v1->x;
		double const y1 = v1->y;
		double const x2 = v2->x;
		double const y2 = v2->y;

		double const detJF = abs((x1 - x0)*(y2 - y0) - (x2 - x0)*(y1 - y0));

		JF(0, 0) = x1 - x0;
		JF(0, 1) = x2 - x0;
		JF(1, 0) = y1 - y0;
		JF(1, 1) = y2 - y0;


		for (unsigned n = 0; n < NumberOfQuadraturePoints; n++) {


			double const s = QuadratureOnTriangle.points_x[n];
			double const t = QuadratureOnTriangle.points_y[n];

			double const x = x0 + JF(0, 0) * s + JF(0, 1) * t;
			double const y = y0 + JF(1, 0) * s + JF(1, 1) * t;

			evaluateBasisRaviartThomas(s, t, BasisRaviartThomas);


			Eigen::Vector2d Velocity(0.0, 0.0);

			for (unsigned i = 0; i < 3; i++)
				Velocity += v(k_index, i) * JF * BasisRaviartThomas.col(i) / detJF;


			OFSTxtFile << std::setprecision(20) << x << "\t" << y << "\t" << Velocity(0) << "\t" << Velocity(1) << std::endl;

		}

	}

	OFSTxtFile.close();

};
template <TimeDiscretization TimeScheme, unsigned QuadraturePrecision>
void Solver<TimeScheme, QuadraturePrecision>::computeError(std::string const & fileName) {


	/*****************************************************************************/
	/*                                                                           */
	/*    - When leaving time for-cycle, time is still on the NT-th level		 */
	/*      but the solution is already on (NT+1)-th level						 */
	/*                                                                           */
	/*****************************************************************************/
	double const time = (nt + 1) * dt;


	quadrature_triangle const QuadratureOnTriangle(QuadraturePrecision);
	unsigned const NumberOfQuadraturePoints = QuadratureOnTriangle.NumberOfPoints;

	Eigen::Matrix2d JF;


	double ErrorL1	= 0.0;
	double ErrorL2	= 0.0;
	double ErrorMax = 0.0;

	for (int k = 0; k < nk; k++) {


		t_pointer const	K		= mesh->get_triangle(k);
		unsigned const k_index	= K->index;

		v_pointer const v0 = K->vertices[0];
		v_pointer const v1 = K->vertices[1];
		v_pointer const v2 = K->vertices[2];

		double const x0 = v0->x;
		double const y0 = v0->y;
		double const x1 = v1->x;
		double const y1 = v1->y;
		double const x2 = v2->x;
		double const y2 = v2->y;

		double const detJF = fabs((x1 - x0)*(y2 - y0) - (x2 - x0)*(y1 - y0));

		JF.coeffRef(0, 0) = x1 - x0;
		JF.coeffRef(0, 1) = x2 - x0;
		JF.coeffRef(1, 0) = y1 - y0;
		JF.coeffRef(1, 1) = y2 - y0;


		double L1NormOnElement	= 0.0;
		double L2NormOnElement	= 0.0;
		double MaxNormOnElement = 0.0;


		//double const PressureK = p[k_index];
		double const PressureK = c[k_index];


		for (unsigned n = 0; n < NumberOfQuadraturePoints; n++) {


			double const s =	   QuadratureOnTriangle.points_x[n];
			double const t =	   QuadratureOnTriangle.points_y[n];
			double const w = 0.5 * QuadratureOnTriangle.weights[n];

			double const x = x0 + JF(0, 0) * s + JF(0, 1) * t;
			double const y = y0 + JF(1, 0) * s + JF(1, 1) * t;

			double const Difference = fabs(PressureK - barenblatt(x, y, time));

			
			L1NormOnElement		+= w * Difference;
			L2NormOnElement		+= w * square(Difference);
			MaxNormOnElement	 = Difference > MaxNormOnElement ? Difference : MaxNormOnElement;

		}

		ErrorL1		+= detJF * L1NormOnElement;
		ErrorL2		+= detJF * L2NormOnElement;
		ErrorMax	 = MaxNormOnElement;

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




template <TimeDiscretization TimeScheme, unsigned QuadraturePrecision>
void Solver<TimeScheme, QuadraturePrecision>::evaluateBasisRaviartThomas(double const s, double const t, Eigen::Matrix<double, 2, 3> & out) {

	out(0, 0) = s;
	out(1, 0) = t;

	out(0, 1) = s - 1.0;
	out(1, 1) = t;

	out(0, 2) = s;
	out(1, 2) = t - 1.0;

};