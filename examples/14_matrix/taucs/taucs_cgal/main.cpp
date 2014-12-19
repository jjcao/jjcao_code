#include <fstream>
#include "calculator.h"

#include <CGAL/IO/Polyhedron_iostream.h>
#include <vector>

void lu_symbolic_factor_solve(Matrix& A, Matrix& MA, Vector& B, Vector& X, std::ofstream& logger)
{
	logger << std::endl << "## Taucs in CGAL do not support lu_symbolic_factor!" << std::endl;
}	
void ooc_lu_factor_solve(Matrix& A,Vector& B, Vector& X, std::ofstream& logger)
{
	int status(0);
#ifdef USE_OPENNL
	std::string funcName = "BICGSTAB in OpenNL version";
#else
	std::string funcName = "ooc_lu_factor_solve in Taucs version";	
#endif
	logger << std::endl << "#### Begin " << funcName << std::endl;

	CGAL::Timer timer;
	timer.start();
	// Solve "A*X = B". On success, solution is (1/D) * X.
	Solver solver = Solver();
	double D;
	if ( !solver.linear_solver(A, B, X, D))
	{
		status = -1;//Base::ERROR_CANNOT_SOLVE_LINEAR_SYSTEM;
		logger << "  Error in solving a linear systems!" << std::endl;	}
	else{
		logger << "  solving a linear systems: " << timer.time() << " seconds." << std::endl;
		timer.reset();
	}

	///////////////////////////////////////
	///////////////////////////////////////
	/*if ( status==0)
	{
		int dim (A.row_dimension());
		for( int i = 0; i < dim; ++i)
		{
			logger << X[i] << ", ";
		}
	}
	logger << std::endl;*/
	logger << "End " << funcName << std::endl;
}
void ll_symbolic_factor_solve(SymmetricMatrix &SA, SymmetricMatrix &SMA, Vector& B, Vector& X, std::ofstream& logger)
{
	int status(0);
	logger << std::endl << "#### Begin ll_symbolic_factor_solve" << std::endl;

#ifdef USE_OPENNL
	logger << "Please not define USE_OPENNL" << std::endl;
#else
	CGAL::Timer timer;
	timer.start();

	int i = 0;
	void* L = taucs_ccs_factor_llt_symbolic((taucs_ccs_matrix*)SA.get_taucs_matrix());
	logger << "....taucs_ccs_factor_llt_symbolic: " << timer.time() << " seconds." << std::endl;
	timer.reset();

	i   = taucs_ccs_factor_llt_numeric((taucs_ccs_matrix*)SMA.get_taucs_matrix(),L);
	if (i)
	{
		logger << "Failed to numerical factorization!" << std::endl;
		status = -1;
	}
	else
	{
		logger << "....taucs_ccs_factor_llt_numeric: " << timer.time() << " seconds." << std::endl;
		timer.reset();
		status = taucs_supernodal_solve_llt(L,X.get_taucs_vector(),B.get_taucs_vector());
	}

	/*char* options[] = { "taucs.factor.LLT=true", NULL };
	void* F = NULL;
	int i = taucs_linsolve((taucs_ccs_matrix*) A.get_taucs_matrix(), 
						&F,0,
						NULL,
						NULL,
						options,NULL);  
	i = taucs_linsolve((taucs_ccs_matrix*) A.get_taucs_matrix(), 
						&F,1,
						X.get_taucs_vector(),
						B.get_taucs_vector(),
						NULL,NULL);  
	taucs_linsolve(NULL,&F,0,NULL,NULL,NULL,NULL);
	if (i != TAUCS_SUCCESS)
	{
		status = -1;
		printf ("Solution error.\n");
		if (i==TAUCS_ERROR)
			printf ("Generic error.");
		if (i==TAUCS_ERROR_NOMEM)
			printf ("NOMEM error.");    
		if (i==TAUCS_ERROR_BADARGS)
			printf ("BADARGS error.");   
		if (i==TAUCS_ERROR_MAXDEPTH)
			printf ("MAXDEPTH error.");    
		if (i==TAUCS_ERROR_INDEFINITE)
			printf ("NOT POSITIVE DEFINITE error."); 
	}*/
	
	///////////////////////////////////////
	///////////////////////////////////////
	/*if ( status==0)
	{
		logger << "  solving a linear systems: " << timer.time() << " seconds." << std::endl;
		for( int i = 0; i < dim; ++i)
		{
			logger << X[i] << ", ";
		}
	}
	logger << std::endl;*/
	logger << "....taucs_supernodal_solve_llt: " << timer.time() << " seconds." << std::endl;
	timer.reset();
	taucs_supernodal_factor_free_numeric(L);	
#endif
	logger << "End ll_symbolic_factor_solve!" << std::endl << std::endl;
}
void fillMatrices(Matrix& A, Matrix& MA, SymmetricMatrix &SA, SymmetricMatrix &SMA, Vector& B)
{
	//A = taucs_ccs_read_ijv("A.txt",TAUCS_DOUBLE);

	int dim (A.row_dimension());
	for( int i = 0; i < dim; ++i)
	{
		A.set_coef(i, i, 1); MA.set_coef(i, i,1);
		SA.set_coef(i, i, 1); SMA.set_coef(i, i,1);
		B[i] = 2;
	}
	MA.set_coef(1, 3, -0.25); MA.set_coef(3, 1, -0.25);
	MA.set_coef(1, 4, -0.25); MA.set_coef(4, 1, -0.25);
	MA.set_coef(1, 8, -0.25); MA.set_coef(8, 1, -0.25);
	MA.set_coef(1, 99, -0.25);MA.set_coef(99, 1, -0.25);
	MA.set_coef(81, 23, -0.2);
	MA.set_coef(81, 4, -0.2); 
	MA.set_coef(81, 38, -0.2); 
	MA.set_coef(81, 299, -0.2);
	MA.set_coef(81, 55, -0.2);
	A.set_coef(1, 3, -0.25); A.set_coef(3, 1, -0.25);
	A.set_coef(1, 4, -0.25); A.set_coef(4, 1, -0.25);
	A.set_coef(1, 8, -0.25); A.set_coef(8, 1, -0.25);
	A.set_coef(1, 99, -0.25);A.set_coef(99, 1, -0.25);
	A.set_coef(4, 13, -0.5); A.set_coef(13, 4, -0.5);
	A.set_coef(4, 199, -0.5);A.set_coef(199, 4, -0.5);
	A.set_coef(81, 23, -0.2); A.set_coef(23, 81, -0.2);
	A.set_coef(81, 4, -0.2); A.set_coef(4, 81, -0.2);
	A.set_coef(81, 38, -0.2); A.set_coef(38, 81, -0.2);
	A.set_coef(81, 299, -0.2);A.set_coef(299, 81, -0.2);
	A.set_coef(81, 55, -0.2);A.set_coef(55, 81, -0.2);

	SMA.set_coef(1, 3, -0.25); SMA.set_coef(3, 1, -0.25);
	SMA.set_coef(1, 4, -0.25); SMA.set_coef(4, 1, -0.25);
	SMA.set_coef(1, 8, -0.25); SMA.set_coef(8, 1, -0.25);
	SMA.set_coef(1, 99, -0.25);SMA.set_coef(99, 1, -0.25);
	SA.set_coef(1, 3, -0.25); SA.set_coef(3, 1, -0.25);
	SA.set_coef(1, 4, -0.25); SA.set_coef(4, 1, -0.25);
	SA.set_coef(1, 8, -0.25); SA.set_coef(8, 1, -0.25);
	SA.set_coef(1, 99, -0.25);SA.set_coef(99, 1, -0.25);
	SA.set_coef(4, 13, -0.5); SA.set_coef(13, 4, -0.5);
	SA.set_coef(4, 199, -0.5);SA.set_coef(199, 4, -0.5);
	SA.set_coef(81, 23, -0.2); SA.set_coef(23, 81, -0.2);
	SA.set_coef(81, 4, -0.2); SA.set_coef(4, 81, -0.2);
	SA.set_coef(81, 38, -0.2); SA.set_coef(38, 81, -0.2);
	SA.set_coef(81, 299, -0.2);SA.set_coef(299, 81, -0.2);
	SA.set_coef(81, 55, -0.2);SA.set_coef(55, 81, -0.2);

	//i=taucs_ccs_write_ijv((taucs_ccs_matrix*)A.get_taucs_matrix(), "A.txt" );
}
void testSimpleMatrix(int dim, std::ofstream &logger)
{
	CGAL::Timer timer;
	timer.start();
	logger << std::endl << "## Begin testSimpleMatrix: " << std::endl;

	if ( dim < 1) return;
	
	Matrix A(dim, dim);          	  Matrix MA(dim, dim);
	SymmetricMatrix SA(dim, dim);SymmetricMatrix SMA(dim, dim);
	Vector X(dim), B(dim);
	fillMatrices(A, MA, SA, SMA, B);	
	logger << ".... fill 4 matries (" << dim << " x " << dim << "): " << timer.time() << " seconds." << std::endl;
	timer.reset();
	
	ll_symbolic_factor_solve(SA, SMA, B, X, logger);	
	lu_symbolic_factor_solve(A, MA, B, X, logger);//u_symbolic_factor_solve(A, MA, B, X, logger);
#ifdef USE_OPENNL
	ooc_lu_factor_solve(MA, B, X, logger);//??, when use ooc_lu in Taucs of mt and mtd version.
#endif

	//create_write_matrix();//checked
	//test2();//checked
	//taucs_pseudo();// cause:A'A is always not symmetric positive definite
}
int testMeshMatrix(char* filename, char* solverType, std::ofstream &logger)
{
	CGAL::Timer timer;
	timer.start();

	logger << std::endl << "## Begin testMeshMatrix: " << std::endl;
	// Read the mesh
    std::ifstream stream(filename);
    Polyhedron mesh;
    stream >> mesh;
    if(!stream || !mesh.is_valid() || mesh.empty())
    {
        logger << "FATAL ERROR: cannot read OFF file " << filename << std::endl;
        return -1;
    }
	logger << "Read mesh: " << filename << " with " << mesh.size_of_vertices() << " vertices in " 
		   << timer.time() << " seconds" << std::endl;
	
	std::list<ConstrianedVertex*> vcMap;
	vcMap.push_back(new ConstrianedVertex(mesh.vertices_begin(),0));

	std::vector <double> rValues;
	rValues.reserve(mesh.size_of_vertices());
	for ( int i = 0; i<mesh.size_of_vertices(); ++i)
	{
		rValues.push_back(0.9);
	}

	Calculator calc(0);//0 MVC, 1 DCP, 5 Tutte
	calc.compute(&mesh, vcMap,rValues,solverType,logger);

	logger << "End testMeshMatrix: " << std::endl;
	return mesh.size_of_vertices();
}
int main(int argc, char* argv[])
{
	if (argc<3)
    {
        std::cout << "Usage: " << argv[0] << " <mesh file name (*.off)> " << " <SolverType>"<< std::endl;
		std::cout << "SolverType: " << std::endl;
		std::cout << "  Taucs::lu_symbolic"<< std::endl;
		std::cout << "  Taucs::ooc_lu"<< std::endl;
		std::cout << "  OpenNL::BICGSTAB"<< std::endl;
		std::cout << "  Taucs::ll_symbolic"<< std::endl;			
        return(EXIT_FAILURE);
    }
	std::ofstream logger("log.txt");	
	CGAL::Timer timer;
	timer.start();

#ifdef _DEBUG
	logger << "Begin main() in Debug mode.The solver is : " << argv[2] << ". The cpu time is:" << timer.time() << std::endl;
#else
	logger << "Begin main() in Release mode.The solver is : " << argv[2] << ". The cpu time is:" << timer.time() << std::endl;
#endif
	timer.reset();

	int dim = testMeshMatrix(argv[1], argv[2], logger);
	//testSimpleMatrix(dim, logger);
	logger << std::endl << "End main(), time: " << timer.time() << std::endl;
	logger.close();

	return EXIT_SUCCESS;
}