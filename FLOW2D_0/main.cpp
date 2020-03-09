/*****************************************************************************/
/*                                                                           */
/*    - To do											 			         */
/*                                                                           */
/*****************************************************************************/
// Neumann boundary : Add second Degree Of Freedom to the velocity prescribed on the Neumann edges
//					  Maybe do the same as for prescribed pressures on Dirichlet edges (Interpolant)
// Check the values at initializeValues() for rkFp, rkFc ... maybe I should compute it also from the initial condition? -> CN scheme will have 2nd order?


/*****************************************************************************/
/*                                                                           */
/*    - Speed up recommendations						 			         */
/*                                                                           */
/*****************************************************************************/
// quadrature_points_x/y calculate only for dirihlet edges/neuman somehow (hash table?)
// Allocate memory in the beginning of the assemble function (typically those Eigen::triplet make only once, not in the loop as const)->will it be faster?
// Try to reserve memory for matrices, see
		/*
			unsigned const NumberOfDirichletEdges	= Mesh->get_number_of_dirichlet_edges();
			unsigned const NumberOfNeumannEdges		= Mesh->get_number_of_neumann_edges();
			unsigned const NumberOfBoundaryEdges	= NumberOfDirichletEdges + NumberOfNeumannEdges;

			unsigned const NumberOfElements			= 4 * (NumberOfDirichletEdges + (NumberOfBoundaryEdges - NumberOfDirichletEdges) * 3 + (ne - NumberOfBoundaryEdges) * 5 + ne - NumberOfBoundaryEdges);

			std::vector<Eigen::Triplet<Real>> triplet;
			triplet.reserve(NumberOfElements);
		*/
		// Use another type of LI(K,e,out) which gives both DOF numbers in one step
		// Interchange the for loops (j-l) in updateConcentrations
		// Make special array of pointers to edges on the boundary (Dirichlet, Neumann) and for Inner edges (No need to check if the edge is neumann/dirichlet)
		// Fast Multiplication Of Sparse matrix and vector (iDH1 * Tp1, iDH2 * Tp2 in computePressure)
		// Precompute values of P(1) on the reference triangle. Now only values on edges are precomputed


#include "mesh.h"

#include "mesh2.h"
#include "triangulation.h"
#include "PSLG.h"
//#include "misc.h"

//#include "solver.h"
//#include "solver2.h"
#include "solver3.h"

#include <iostream>
#include <fstream>
#include <ctime>
#include <vector>
#include <iomanip>

void generate_vertices();
void exportMesh(std::ofstream & txtFile, Mesh const & m);



//const double a_x = 0.0;
//const double b_x = 25.0;
//
//const double a_y = 0.0;
//const double b_y = 25.0;

//const double a_x = -25.0;
//const double b_x = 25.0;
//
//const double a_y = -25.0;
//const double b_y = 25.0;

const double a_x = 0.0;
const double b_x = 40.0;

const double a_y = 0.0;
const double b_y = 40.0;

unsigned const refinement = 2 * 2;

const int N_x = 4 * refinement;
const int N_y = N_x;

//unsigned const nt0 = 25 * (refinement * refinement);
//unsigned const NT = 150 * (refinement * refinement);
//double const dt = 300.0 / (refinement * refinement);

unsigned const nt0 = 25 * (refinement);
unsigned const NT = 150 * (refinement);
double const dt = 300.0 / (refinement);




// const reference to an object T (T const & ) -> const is redunant, reerence is always cosnt (it will never be reference to anz=ythig else)




std::string fileName_mesh		= "C:\\Users\\pgali\\Desktop\\eoc\\mesh.txt";
std::string fileName_velocity	= "C:\\Users\\pgali\\Desktop\\eoc\\velocity";

std::string fileName_pressure		= "C:\\Users\\pgali\\Desktop\\eoc\\pressure";
std::string fileName_concentration	= "C:\\Users\\pgali\\Desktop\\eoc\\concentration";
std::string fileName_error			 = "C:\\Users\\pgali\\Desktop\\eoc\\error";



std::vector<Vertex>		vertices;

int main() {


	PlanarStraightLineGraph pslg;
	std::vector<Vertex>		seeds;
	GEOMETRIC_KERNEL const	GK = GEOMETRIC_KERNEL::INEXACT_PREDICATES_INEXACT_CONSTRUCT;


	/*****************************************************************************/
	/*                                                                           */
	/*    - Generate vertices for the mesh triangulation						 */
	/*                                                                           */
	/*****************************************************************************/
	generate_vertices();



	/*****************************************************************************/
	/*                                                                           */
	/*    - Insert vertices in to Planar Straight Line Graph and define			 */
	/*      constraints (boundary edges)										 */
	/*                                                                           */
	/*****************************************************************************/
	for (size_t i = 0; i < vertices.size(); i++)
		pslg.insert_vertex<V_MARKER::FREE>(vertices[i]);

	v_pointer const va = pslg.get_vertex(0);
	v_pointer const vb = pslg.get_vertex(N_y - 1);
	v_pointer const vc = pslg.get_vertex(N_x*N_y - 1);
	v_pointer const vd = pslg.get_vertex(N_x*N_y - 1 - (N_y - 1));


	pslg.insert_constraint<E_MARKER::NEUMANN>(va, vb);
	pslg.insert_constraint<E_MARKER::DIRICHLET>(vb, vc);
	pslg.insert_constraint<E_MARKER::DIRICHLET>(vc, vd);
	pslg.insert_constraint<E_MARKER::NEUMANN>(vd, va);

	//pslg.insert_constraint<E_MARKER::DIRICHLET>(va, vb);
	//pslg.insert_constraint<E_MARKER::DIRICHLET>(vb, vc);
	//pslg.insert_constraint<E_MARKER::DIRICHLET>(vc, vd);
	//pslg.insert_constraint<E_MARKER::DIRICHLET>(vd, va);



	/*****************************************************************************/
	/*                                                                           */
	/*    - Construct triangulation from the input PSLG and defined holes		 */
	/*			: Problem with Mesh (less probably triangulation)				 */
	/*            when inserting seed                                            */
	/*                                                                           */
	/*****************************************************************************/
	//seeds.push_back(Vertex(8.0, 8.0));
	Triangulation<GK>	triangulation(pslg, seeds);


	/*****************************************************************************/
	/*                                                                           */
	/*    - Print Triangulation	information										 */
	/*                                                                           */
	/*****************************************************************************/
	std::cout << "/*****************/" << std::endl;
	std::cout << "/*               */" << std::endl;
	std::cout << "/* Triangulation */" << std::endl;
	std::cout << "/*               */" << std::endl;
	std::cout << "/*****************/" << std::endl;
	std::cout << std::endl;
	cout << "No. Vertices  : " << triangulation.get_number_of_vertices() << endl;
	cout << "No. Edges     : " << triangulation.get_number_of_edges() << endl;
	cout << "No. Triangles : " << triangulation.get_number_of_triangles() << endl;
	cout << endl;



	/*****************************************************************************/
	/*                                                                           */
	/*    - Construct computation mesh from the triangulation					 */
	/*                                                                           */
	/*****************************************************************************/
	//Mesh				mesh(triangulation);
	Mesh2				mesh(triangulation);


	/*****************************************************************************/
	/*                                                                           */
	/*    - Create text file of the mesh: coordinates							 */
	/*                                                                           */
	/*****************************************************************************/
	//OFSTxtFile.open(directory_mesh);
	//exportMesh(OFSTxtFile, mesh);
	//OFSTxtFile.close();



	/*****************************************************************************/
	/*                                                                           */
	/*    - Create instance of the solver										 */
	/*                                                                           */
	/*****************************************************************************/
	//solver<quadrature_order> solution(mesh, nt0, dt);
	//solver2<double, 7, scheme::CRANK_NICOLSON> solution(mesh, nt0, dt);
	solver3<7, scheme::CRANK_NICOLSON> solution(mesh, nt0, dt);



	/*****************************************************************************/
	/*                                                                           */
	/*    - Create text file of the initial condition							 */
	/*                                                                           */
	/*****************************************************************************/
	solution.exportPressures(fileName_pressure + "_" + std::to_string(nt0) + ".txt");
	//solution.exportConcentrations(fileName_concentration + "_" + std::to_string(nt0) + ".txt");



	std::cout << std::endl;
	std::cout << "/*****************/" << std::endl;
	std::cout << "/*               */" << std::endl;
	std::cout << "/*  Computing... */" << std::endl;
	std::cout << "/*               */" << std::endl;
	std::cout << "/*****************/" << std::endl;
	std::cout << std::endl;




	clock_t const begin = clock();

	for (int nt = nt0 + 1; nt < NT + 1; nt++) {


		/*****************************************************************************/
		/*                                                                           */
		/*    - Solver is on the n-th time level									 */
		/*    - Compute solution on the (n+1)-th time level							 */
		/*    - Set the solver to the new (n+1)-th time level						 */
		/*                                                                           */
		/*****************************************************************************/
		solution.getSolution();
		solution.setTimeLevel(nt);


		/*****************************************************************************/
		/*                                                                           */
		/*    - Create text file of the solution on the (n+1)-th time level			 */
		/*                                                                           */
		/*****************************************************************************/
		solution.exportPressures(fileName_pressure + "_" + std::to_string(nt) + ".txt");
		//solution.exportConcentrations(fileName_concentration + "_" + std::to_string(nt) + ".txt");


		/*****************************************************************************/
		/*                                                                           */
		/*    - Create text file of the velocities on the (n+1)-th time level		 */
		/*                                                                           */
		/*****************************************************************************/
		//solution.exportVelocities(fileName_velocity + "_" + std::to_string(nt) + ".txt");


	}

	clock_t const end = clock();




	//solution.computeError(fileName_error + std::to_string(refinement));



	std::cout << std::endl;
	std::cout << "CPU clocks : " << end - begin << std::endl << std::endl;

	system("pause");

	return 0;

}



void generate_vertices() {

	for (int i = 0; i < N_x; i++) {

		double const x = a_x + i * (b_x - a_x) / (N_x - 1);

		for (int j = 0; j < N_y; j++) {

			double const y = a_y + j * (b_y - a_y) / (N_y - 1);

			vertices.push_back(Vertex(x, y));

		}
	}

};

void exportMesh(std::string const & fileName, Mesh2 const & m) {


	std::ofstream OFSTxtFile(fileName);

	for (unsigned k = 0; k < m.get_number_of_triangles(); k++) {


		tm_pointer const	K = m.get_triangle(k);

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


		for (unsigned i = 0; i < 3; i++)
			OFSTxtFile << std::setprecision(20) << x[i] << "\t" << y[i] << std::endl;

		OFSTxtFile << std::endl;

	}

	OFSTxtFile.close();

};