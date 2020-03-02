

#include "mesh.h"
#include "triangulation.h"
#include "PSLG.h"
#include "misc.h"

#include "solver.h"

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

unsigned const refinement = 2*2*2*2*2;

const int N_x = 4 * refinement;
const int N_y = N_x;

unsigned const nt0 = 125 * (refinement * refinement);
unsigned const NT = 150 * (refinement * refinement);
double const dt = 300.0 / (refinement * refinement);






std::string directory_mesh			= "C:\\Users\\pgali\\Desktop\\eoc\\mesh.txt";
std::string directory_velocities	= "C:\\Users\\pgali\\Desktop\\eoc\\velocities_";

std::string directory_solution		= "C:\\Users\\pgali\\Desktop\\eoc\\output_";
std::string directory_error			= "C:\\Users\\pgali\\Desktop\\eoc\\error_";

std::ofstream OFSTxtFile;


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
	Mesh				mesh(triangulation);



	
	/*****************************************************************************/
	/*                                                                           */
	/*    - Print Triangulation / Mesh information								 */
	/*                                                                           */
	/*****************************************************************************/
	cout << "*************** Triangulation ***************" << endl;
	cout << "No. Vertices  : " << triangulation.get_number_of_vertices() << endl;
	cout << "No. Edges     : " << triangulation.get_number_of_edges() << endl;
	cout << "No. Triangles : " << triangulation.get_number_of_triangles() << endl;
	cout << "*********************************************" << endl << endl;





	/*****************************************************************************/
	/*                                                                           */
	/*    - Create text file of the mesh: coordinates							 */
	/*                                                                           */
	/*****************************************************************************/
	OFSTxtFile.open(directory_mesh);
	exportMesh(OFSTxtFile, mesh);
	OFSTxtFile.close();



	/*****************************************************************************/
	/*                                                                           */
	/*    - Create instance of the solver										 */
	/*                                                                           */
	/*****************************************************************************/
	solver<quadrature_order> solution(mesh, nt0, dt);




	/*****************************************************************************/
	/*                                                                           */
	/*    - Create text file of the initial condition							 */
	/*                                                                           */
	/*****************************************************************************/
	OFSTxtFile.open(directory_solution + std::to_string(nt0) + ".txt");
	solution.exportSolution(OFSTxtFile);
	OFSTxtFile.close();

	

	clock_t begin = clock();

	for (unsigned nt = nt0 + 1; nt < NT + 1; nt++) {

		

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
		//OFSTxtFile.open(directory_solution + std::to_string(nt) + ".txt");
		//solution.exportSolution(OFSTxtFile);
		//OFSTxtFile.close();


		/*****************************************************************************/
		/*                                                                           */
		/*    - Create text file of the velocities on the (n+1)-th time level		 */
		/*                                                                           */
		/*****************************************************************************/
		//OFSTxtFile.open(directory_velocities + std::to_string(nt) + ".txt");
		//solution.exportSolution(OFSTxtFile);
		//OFSTxtFile.close();


	}

	clock_t end = clock();





	OFSTxtFile.open(directory_error + std::to_string(N_x) + ".txt");
	solution.compute_error(OFSTxtFile);
	OFSTxtFile.close();





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

void exportMesh(std::ofstream & txtFile, Mesh const & m) {


	for (unsigned k = 0; k < m.get_number_of_triangles(); k++) {


		t_pointer const	K = m.get_triangle(k);

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


		for (unsigned i = 0; i < 3; i++)
			txtFile << std::setprecision(20) << x[i] << "\t" << y[i] << std::endl;

		txtFile << std::endl;

	}


};