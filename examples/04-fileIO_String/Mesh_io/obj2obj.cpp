// obj2obj.cpp

#pragma warning(disable : 4503 )

#pragma warning (push)
	#pragma warning(disable : 4267 4996 4819)
	#include <CGAL/Cartesian.h>
	#include <DGAL/Polyhedron_3.h>
	#include <DGAL/IO/Polyhedron_scan_obj.h>
	#include <DGAL/IO/Polyhedron_outputer_obj.h>

#pragma warning (pop)

#include <iostream>
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
    std::cout << "Read obj to an polyhedron and output the polyhedron to an obj file!" << std::endl;
   
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

    // File name is:
    const char* input_filename  = argv[1];
	const char* output_filename  = argv[2];

    //***************************************
    // Read the mesh
    //***************************************
	std::ifstream stream(input_filename);
    if(!stream)
    {
        std::cerr << "FATAL ERROR: cannot open file " << input_filename << std::endl;
        return EXIT_FAILURE;
    }
    Polyhedron mesh;
	DGAL::Polyhedron_scan_obj<Polyhedron::HDS> bd(stream);
	mesh.delegate(bd);
	stream.close();

    //***************************************
    // output statistics info
    //***************************************
	//todo 


    //***************************************
    // output the mesh
    //***************************************
	std::string ofname(output_filename);
	std::ofstream ostream(ofname.c_str());
	DGAL::Polyhedron_outputer_obj<Polyhedron> outputer(mesh, ostream);
	outputer.save();
	ostream.close();


	//***************************************
    // compare the input file and the output file
    //***************************************
	//todo    

    return 0;
}


