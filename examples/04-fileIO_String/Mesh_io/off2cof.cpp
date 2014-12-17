#pragma warning(disable : 4503 )

#pragma warning (push)
	#pragma warning(disable : 4267 4996 )
	#include <CGAL/Cartesian.h>
	#include <CGAL/Polyhedron_3.h>
	#include <CGAL/IO/Polyhedron_iostream.h>
	#include <DGAL/Polyhedron_3.h>
	#include <DGAL/IO/Polyhedron_scan_cof.h>
	#include <DGAL/IO/Polyhderon_outputer_cof.h>
#pragma warning (pop)

#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <fstream>

// ----------------------------------------------------------------------------
// Private types
// ----------------------------------------------------------------------------

typedef CGAL::Cartesian<double>             Kernel;
typedef DGAL::Polyhedron_3<Kernel>       Polyhedron;

// ----------------------------------------------------------------------------
// main()
// ----------------------------------------------------------------------------

int main(int argc, char * argv[])
{
    //***************************************
    // decode parameters
    //***************************************

    if (argc != 3)
    {
        std::cerr << "Usage: " << argv[0] << " input_file.off" 
			 << " output_file.cof" << std::endl;
		int i;
		std::cin >> i;
        return(EXIT_FAILURE);
    }

    const char* input_filename   = argv[1];
	const char* output_filename  = argv[2];

    //***************************************
    // Read the mesh in OFF format
    //***************************************

	std::ifstream istream(input_filename);
    if(!istream)
    {
        std::cerr << "FATAL ERROR: cannot open file " << input_filename << std::endl;
        return EXIT_FAILURE;
    }
    Polyhedron mesh;
	istream >> mesh;
	istream.close();

	//***************************************
    // Output the mesh in COF format without constraint
    //***************************************
	std::string ofname(output_filename);
	std::ofstream ostream(ofname.c_str());
	DGAL::Polyhderon_outputer_cof<Polyhedron> outputer(mesh, ostream);
	outputer.save();
	ostream.close();
}