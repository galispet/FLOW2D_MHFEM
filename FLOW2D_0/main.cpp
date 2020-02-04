

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



//const double a_x = 0.0;
//const double b_x = 60.0;
//
//const double a_y = 0.0;
//const double b_y = 60.0;

const double a_x = -30.0;
const double b_x = 30.0;

const double a_y = -30.0;
const double b_y = 30.0;

const int N_x = 40;
const int N_y = 40;


unsigned const nt0 = 20;
unsigned const NT = 50;
double const dt = 50;



std::string directory_solution = "C:\\Users\\pgali\\Desktop\\flow2d\\output_c_";
std::string directory_error = "C:\\Users\\pgali\\Desktop\\error_";

std::ofstream txtFile;
std::ofstream txtFile_error;

std::vector<Vertex> vertices;


PlanarStraightLineGraph pslg;
std::vector<Vertex> seeds;
GEOMETRIC_KERNEL const GK = GEOMETRIC_KERNEL::INEXACT_PREDICATES_INEXACT_CONSTRUCT;


void generate_vertices();



int main() {




	// *************************** Generating mesh ***************************
	// ***********************************************************************

	// TODO
	// Mesh - sort by the wish of the user ->template parameter sorting (1.inner, 2. neumann, 3. dirichlet and permutations ....) 

	generate_vertices();

	for (size_t i = 0; i < vertices.size(); i++)
		pslg.insert_vertex<V_MARKER::FREE>(vertices[i]);

	v_pointer const va = pslg.get_vertex(0);
	v_pointer const vb = pslg.get_vertex(N_y - 1);
	v_pointer const vc = pslg.get_vertex(N_x*N_y - 1);
	v_pointer const vd = pslg.get_vertex(N_x*N_y - 1 - (N_y - 1));


	//pslg.insert_constraint<E_MARKER::NEUMANN>(va, vb);
	//pslg.insert_constraint<E_MARKER::DIRICHLET>(vb, vc);
	//pslg.insert_constraint<E_MARKER::DIRICHLET>(vc, vd);
	//pslg.insert_constraint<E_MARKER::NEUMANN>(vd, va);

	pslg.insert_constraint<E_MARKER::DIRICHLET>(va, vb);
	pslg.insert_constraint<E_MARKER::DIRICHLET>(vb, vc);
	pslg.insert_constraint<E_MARKER::DIRICHLET>(vc, vd);
	pslg.insert_constraint<E_MARKER::DIRICHLET>(vd, va);


	//Vertex  vaa(11.0, 5.5);
	//Vertex  vbb(11.0, 11.0);
	//Vertex  vcc(5.5, 11.0);
	//Vertex  vdd(5.5, 5.5);

	//pslg.insert_vertex<V_MARKER::FREE>(vaa);
	//pslg.insert_vertex<V_MARKER::FREE>(vbb);
	//pslg.insert_vertex<V_MARKER::FREE>(vcc);
	//pslg.insert_vertex<V_MARKER::FREE>(vdd);

	//v_pointer const _vaa = pslg.get_vertex(vertices.size()-4);
	//v_pointer const _vbb = pslg.get_vertex(vertices.size()-3);
	//v_pointer const _vcc = pslg.get_vertex(vertices.size()-2);
	//v_pointer const _vdd = pslg.get_vertex(vertices.size()-1);

	//pslg.insert_constraint<E_MARKER::DIRICHLET>(_vaa, _vbb);
	//pslg.insert_constraint<E_MARKER::DIRICHLET>(_vbb, _vcc);
	//pslg.insert_constraint<E_MARKER::DIRICHLET>(_vcc, _vdd);
	//pslg.insert_constraint<E_MARKER::DIRICHLET>(_vdd, _vaa);






	//pslg.insert_vertex<V_MARKER::FREE>(Vertex(0.0, 0.0));
	//pslg.insert_vertex<V_MARKER::FREE>(Vertex(50.0, 0.0));
	//pslg.insert_vertex<V_MARKER::FREE>(Vertex(100.0, 0.0));
	//pslg.insert_vertex<V_MARKER::FREE>(Vertex(0.0, 50.0));
	//pslg.insert_vertex<V_MARKER::FREE>(Vertex(50.0, 50.0));

	//v_pointer const va = pslg.get_vertex(0);
	//v_pointer const vb = pslg.get_vertex(1);
	//v_pointer const vc = pslg.get_vertex(2);
	//v_pointer const vd = pslg.get_vertex(3);
	//v_pointer const ve = pslg.get_vertex(4);

	//pslg.insert_constraint<E_MARKER::NEUMANN>(va, vb);
	//pslg.insert_constraint<E_MARKER::NEUMANN>(vb, vc);
	//pslg.insert_constraint<E_MARKER::NEUMANN>(vc, ve);
	//pslg.insert_constraint<E_MARKER::NEUMANN>(ve, vd);
	//pslg.insert_constraint<E_MARKER::NEUMANN>(vd, va);
	//pslg.insert_constraint<E_MARKER::NEUMANN>(vd, vb);
	//pslg.insert_constraint<E_MARKER::NEUMANN>(vb, ve);

	//seeds.push_back(Vertex(8.0, 8.0));

	// Problem with Mesh when inserting seed
	Triangulation<GK> triangulation(pslg, seeds);

	// Problem with Mesh when inserting seed
	Mesh mesh(triangulation);

	
	cout << "*************** Triangulation ***************" << endl;
	cout << "No. Vertices  : " << triangulation.get_number_of_vertices() << endl;
	cout << "No. Edges     : " << triangulation.get_number_of_edges() << endl;
	cout << "No. Triangles : " << triangulation.get_number_of_triangles() << endl;
	cout << "*********************************************" << endl << endl;
	// ***********************************************************************
	// ***********************************************************************


	solver<quadrature_order> solution(mesh, nt0, dt);





	txtFile.open(directory_solution + std::to_string(nt0) + ".txt");
	//txtFile.open("C:\\Users\\pgali\\Desktop\\output_pressure.txt");
	//txtFile.open("C:\\Users\\pgali\\Desktop\\flow2d\\output_c.txt");
	solution.exportSolution(txtFile);
	txtFile.close();

	
	clock_t begin = clock();

	for (unsigned nt = nt0 + 1; nt < NT + 1; nt++) {


		solution.getSolution();

		solution.setTimeLevel(nt);


		txtFile.open(directory_solution + std::to_string(nt) + ".txt");
		solution.exportSolution(txtFile);
		txtFile.close();

	}

	clock_t end = clock();

	//txtFile.open(directory_solution + std::to_string(NT) + ".txt");
	//solution.exportSolution(txtFile);
	//txtFile.close();


	//txtFile_error.open(directory_error + std::to_string(N_x*N_y) + ".txt");
	solution.compute_error(txtFile_error);
	//txtFile_error.close();





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