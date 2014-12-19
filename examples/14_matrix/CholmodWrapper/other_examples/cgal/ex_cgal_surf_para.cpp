// =====================================================================	
// An example to show how to use CholmodWrapper as a solver trati
// for Surface_mesh_parameterization in CGAL. 
// =====================================================================	

// =====================================================================	
#include <iostream>
#include <fstream>
#include <CGAL/Cartesian.h>
#include <CGAL/Timer.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/parameterize.h>
#include <CGAL/Parameterization_polyhedron_adaptor_3.h>
// border parameterizer headers
#include <CGAL/Two_vertices_parameterizer_3.h>
#include <CGAL/Circular_border_parameterizer_3.h>
#include <CGAL/Square_border_parameterizer_3.h>
// essential parameterizer headers
#include <CGAL/Barycentric_mapping_parameterizer_3.h>
#include <CGAL/Discrete_conformal_map_parameterizer_3.h>
#include <CGAL/Discrete_authalic_parameterizer_3.h>
#include <CGAL/Mean_value_coordinates_parameterizer_3.h>
#include <CGAL/LSCM_parameterizer_3.h>
// polyhedron io
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/IO/print_OFF.h>
// solver headers
#include <CholmodWrapper/Cholmod_solver_traits.h>
// =====================================================================	

// =====================================================================	

typedef CGAL::Cartesian<double> Kernel;
typedef CGAL::Polyhedron_3<Kernel> Mesh;
typedef CGAL::Parameterization_polyhedron_adaptor_3<Mesh> Para_mesh_adaptor;
typedef Cholmod_solver_traits<double> Solver;

// =====================================================================	
#if 1	// free para boundary
typedef CGAL::LSCM_parameterizer_3
	<Para_mesh_adaptor,
	CGAL::Two_vertices_parameterizer_3<Para_mesh_adaptor>, 
	Solver> Parameterizer;
#else	// fixed para boundary
// ---------------------------------------------  surface parameterization type
//typedef CGAL::Discrete_authalic_parameterizer_3
//typedef CGAL::Discrete_conformal_map_parameterizer_3
typedef CGAL::Barycentric_mapping_parameterizer_3
//typedef CGAL::Mean_value_coordinates_parameterizer_3
// ---------------------------------------------  surface parameterization type
	<Para_mesh_adaptor,
	// --------------------------------------------- boundry paramterization type
	//CGAL::Circular_border_arc_length_parameterizer_3<Para_mesh_adaptor>, 
	CGAL::Square_border_arc_length_parameterizer_3<Para_mesh_adaptor>, 			  	
	// --------------------------------------------- boundry paramterization type
	Solver> Parameterizer;
#endif

// =====================================================================	

bool mesh_parameterization(Mesh& input_mesh, Mesh& para_mesh)
{
	Para_mesh_adaptor mesh_adaptor(input_mesh);
	Parameterizer::Error_code err = 
		CGAL::parameterize(mesh_adaptor, Parameterizer());	

	if (err != Parameterizer::OK)
	{
		std::cerr << "FATAL ERROR: " << 
			Parameterizer::get_error_message(err) << std::endl;		
		return false;
	}

	para_mesh = input_mesh;
	for (Mesh::Vertex_iterator i_it = input_mesh.vertices_begin(),
		o_it = para_mesh.vertices_begin();
		i_it != input_mesh.vertices_end(); ++i_it, ++o_it)
	{
		o_it->point() = Mesh::Traits::Point_3(
				mesh_adaptor.info(i_it->halfedge())->uv().x(),
				mesh_adaptor.info(i_it->halfedge())->uv().y(),
				0);		
	}
	return true;
}

int main(int argc, char* argv[])
{
	const char* input_file	= "holes.off";
	const char* output_file = "holes_para.off";

	// ---------------------------------------------
	// input mesh
	std::ifstream ifs(input_file);
	if (ifs.bad()) 
	{
		std::cerr << "Error: cannot open the input file.\n";
		return -1;
	}
	Mesh input_mesh, para_mesh;
	ifs >> input_mesh;
	ifs.close();

	// ---------------------------------------------
	// surface parameterization
	CGAL::Timer task_timer;
	std::cout << "begin parameterization ......\n";
	task_timer.start();	
	if (!mesh_parameterization(input_mesh, para_mesh))
	{
		std::cerr << "Error: canot parameterize the mesh.\n";
		return -1;
	}
	std::cout << "parameterization end: " << task_timer.time() << "s\n";

	// ---------------------------------------------
	// output parameterized mesh
	std::ofstream ofs(output_file);
	if (ofs.bad())
	{
		std::cerr << "Error: cannot open the output file.\n";
		return -1;
	}
	CGAL::print_OFF (ofs, para_mesh);
	ofs.close();

	return 0;
}

