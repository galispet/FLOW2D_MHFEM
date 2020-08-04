

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



double const ax = 0.0;
double const bx = 40.0;

double const ay = 0.0;
double const by = 40.0;

unsigned const refinement = 2*2*2*2*2;

const int Nx = 4 * refinement;
const int Ny = Nx;

unsigned const NT0	= 25 * (refinement);
unsigned const NT	= 150 * (refinement);
double const dt		= 300.0 / (refinement);


std::string fileName_mesh			= "D:\\simulations\\multicomponentflow\\mesh.txt";
std::string fileName_velocity		= "D:\\simulations\\multicomponentflow\\velocity_";

std::string fileName_pressure		= "D:\\simulations\\multicomponentflow\\pressure_";
std::string fileName_concentration	= "D:\\simulations\\multicomponentflow\\concentration_";
std::string fileName_error			= "D:\\simulations\\multicomponentflow\\error_";

std::ofstream txtFile;





PlanarStraightLineGraph pslg;
std::vector<Vertex> seeds;
GEOMETRIC_KERNEL const GK = GEOMETRIC_KERNEL::INEXACT_PREDICATES_INEXACT_CONSTRUCT;

TimeDiscretization const TimeScheme			= TimeDiscretization::CrankNicolson;
unsigned const IntegrationPrecisionOrder	= 6;



void exportMesh(std::ofstream & txtFile, Mesh const & m);


int main() {


	std::vector<Vertex> vertices;

	for (int i = 0; i < Nx; i++) {
		for (int j = 0; j < Ny; j++) {

			double const x = ax + i * (bx - ax) / (Nx - 1);
			double const y = ay + j * (by - ay) / (Ny - 1);

			vertices.push_back(Vertex(x, y));

		}
	}

	for (size_t i = 0; i < vertices.size(); i++)
		pslg.insert_vertex<V_MARKER::FREE>(vertices[i]);

	v_pointer const va = pslg.get_vertex(0);
	v_pointer const vb = pslg.get_vertex(Ny - 1);
	v_pointer const vc = pslg.get_vertex(Nx*Ny - 1);
	v_pointer const vd = pslg.get_vertex(Nx*Ny - 1 - (Ny - 1));


	pslg.insert_constraint<E_MARKER::NEUMANN>(va, vb);
	pslg.insert_constraint<E_MARKER::DIRICHLET>(vb, vc);
	pslg.insert_constraint<E_MARKER::DIRICHLET>(vc, vd);
	pslg.insert_constraint<E_MARKER::NEUMANN>(vd, va);

	//pslg.insert_constraint<E_MARKER::DIRICHLET>(va, vb);
	//pslg.insert_constraint<E_MARKER::DIRICHLET>(vb, vc);
	//pslg.insert_constraint<E_MARKER::DIRICHLET>(vc, vd);
	//pslg.insert_constraint<E_MARKER::DIRICHLET>(vd, va);

	// Problem with Mesh when inserting seed
	//seeds.push_back(Vertex(8.0, 8.0));



	Triangulation<GK> triangulation(pslg, seeds);

	Mesh mesh(triangulation);

	

	cout << "*************** Triangulation ***************" << endl;
	cout << "No. Vertices  : " << triangulation.get_number_of_vertices() << endl;
	cout << "No. Edges     : " << triangulation.get_number_of_edges() << endl;
	cout << "No. Triangles : " << triangulation.get_number_of_triangles() << endl;
	cout << "*********************************************" << endl << endl;


	txtFile.open(fileName_mesh);
	exportMesh(txtFile, mesh);
	txtFile.close();


	Solver<TimeScheme, IntegrationPrecisionOrder> solution(mesh, NT0, dt);

	//solution.exportPressures(fileName_pressure + std::to_string(NT0) + ".txt");
	//solution.exportConcentrations(fileName_concentration + "_" + std::to_string(nt0) + ".txt");


	
	clock_t begin = clock();

	for (unsigned nt = NT0 + 1; nt < NT + 1; nt++) {


		solution.getSolution();
		solution.setTimeLevel(nt);

		//solution.exportPressures(fileName_pressure + std::to_string(nt) + ".txt");
		//solution.exportConcentrations(fileName_concentration + std::to_string(nt) + ".txt");

	}

	clock_t end = clock();


	//solution.exportPressures(fileName_pressure + std::to_string(NT) + ".txt");
	//solution.exportConcentrations(fileName_concentration + std::to_string(NT) + ".txt");


	solution.computeError(fileName_error + std::to_string(Nx) + ".txt");


	std::cout << std::endl;
	std::cout << "CPU clocks : " << end - begin << std::endl << std::endl;

	system("pause");


	return 0;

}




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