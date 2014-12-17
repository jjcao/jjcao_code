#include <fstream>
//#include <malloc.h>
//#include <taucsaddon.h>
#include "calculator.h"

#include <CGAL/IO/Polyhedron_iostream.h>
#include <vector>
void ooc_lu_factor_solve(Matrix& A,Vector& B, Vector& X, std::ofstream& logger)
{
	int status(0);	
#ifdef USE_OPENNL
	std::string funcName = "BICGSTAB in OpenNL version";
#else
	std::string funcName = "ooc_lu_factor_solve in Taucs version";	
#endif
	logger << std::endl << "#### Begin " << funcName << std::endl;

	double oldTime = taucs_ctime();
	logger << "start time: " << oldTime << std::endl;

	// Solve "A*X = B". On success, solution is (1/D) * X.
	Solver solver = Solver();
	double D;
	if ( !solver.linear_solver(A, B, X, D))
	{
		status = -1;//Base::ERROR_CANNOT_SOLVE_LINEAR_SYSTEM;
		logger << "  Error in solving a linear systems!" << std::endl;
	}
	else{
		logger << "  solving a linear systems: " << (taucs_ctime() - oldTime) << " seconds." << std::endl;
	}

	///////////////////////////////////////
	///////////////////////////////////////
	//if ( status==0)
	//{
	//	int dim (A.row_dimension());
	//	for( int i = 0; i < dim; ++i)
	//	{
	//		logger << X[i] << ", ";
	//	}
	//}
	//logger << std::endl;
	logger << "End " << funcName << std::endl;
}
void lu_symbolic_factor_solve(Matrix& A, Matrix& MA, Vector& B, Vector& X, std::ofstream& logger)
{
	int status(0);
	logger << std::endl << "#### Begin lu_symbolic_factor_solve" << std::endl;
#ifndef USE_OPENNL
	double oldTime = taucs_ctime() ;
	logger << "start time: " << oldTime << std::endl;

	///////////////////////////////////////
	///////////////////////////////////////
	int* rowperm = NULL;
    int* colperm = NULL;
	char msg[] = "colamd";
    taucs_ccs_order((taucs_ccs_matrix*)A.get_taucs_matrix(),&rowperm,&colperm, msg);
	logger << "  taucs_ccs_order: " << taucs_ctime() - oldTime << " seconds." << std::endl;
	oldTime = taucs_ctime() ;
	
	taucs_multilu_symbolic* LS = taucs_ccs_factor_lu_symbolic((taucs_ccs_matrix*)A.get_taucs_matrix(),colperm);
	logger << "  taucs_ccs_factor_lu_symbolic: " << taucs_ctime() - oldTime << " seconds." << std::endl;
	oldTime = taucs_ctime() ;

	taucs_multilu_factor* LV = taucs_ccs_factor_lu_numeric((taucs_ccs_matrix*)MA.get_taucs_matrix(),LS,100);
	logger << "  taucs_ccs_factor_lu_numeric: " << taucs_ctime() - oldTime << " seconds." << std::endl;
	oldTime = taucs_ctime() ;

	if (!LV)
	{
		logger << "Failed to numerical factorization!" << std::endl;
		status = -1;
	}
	else
	{
		status = taucs_multilu_solve(LV,X.get_taucs_vector(),B.get_taucs_vector());
		if ( status==0)
		{
			logger << "  solving a linear systems: " << taucs_ctime() - oldTime << " seconds." << std::endl;	
			/*int dim (A.row_dimension());
			for( int i = 0; i < dim; ++i)
			{
				logger << X[i] << ", ";
			}*/
		}
		logger << std::endl;
		taucs_multilu_symbolic_free(LS);
		taucs_multilu_factor_free(LV);		
	}
#endif
	logger << "End lu_symbolic_factor_solve" << std::endl << std::endl;
}
void ll_symbolic_factor_solve(SymmetricMatrix &SA, SymmetricMatrix &SMA, Vector& B, Vector& X, std::ofstream& logger)
{
	int status(0);
	logger << std::endl << "#### Begin ll_symbolic_factor_solve" << std::endl;
#ifdef USE_OPENNL
	logger << "Please not define USE_OPENNL" << std::endl;
#else
	double oldTime = taucs_ctime() ;

	int i = 0;
	void* L = taucs_ccs_factor_llt_symbolic((taucs_ccs_matrix*)SA.get_taucs_matrix());
	logger << "....taucs_ccs_factor_llt_symbolic: " << taucs_ctime()-oldTime << " seconds." << std::endl;
	oldTime = taucs_ctime();

	i   = taucs_ccs_factor_llt_numeric((taucs_ccs_matrix*)SMA.get_taucs_matrix(),L);
	if (i)
	{
		logger << "Failed to numerical factorization!" << std::endl;
		status = -1;
	}
	else
	{
		logger << "....taucs_ccs_factor_llt_numeric: " << taucs_ctime()-oldTime << " seconds." << std::endl;
		oldTime = taucs_ctime();
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
		logger << "  solving a linear systems: " << taucs_ctime()-oldTime << " seconds." << std::endl;
		oldTime = taucs_ctime();
		for( int i = 0; i < dim; ++i)
		{
			logger << X[i] << ", ";
		}
	}
	logger << std::endl;*/
	
	logger << "....taucs_supernodal_solve_llt: " << taucs_ctime()-oldTime << " seconds." << std::endl;
	oldTime = taucs_ctime();
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
	double oldTime = taucs_ctime();
	logger << std::endl << "## Begin testSimpleMatrix: " << oldTime << std::endl;

	if ( dim < 1) return;
	
	Matrix A(dim, dim);          	  Matrix MA(dim, dim);
	SymmetricMatrix SA(dim, dim);SymmetricMatrix SMA(dim, dim);
	Vector X(dim), B(dim);
	fillMatrices(A, MA, SA, SMA, B);	
	logger << ".... fill 4 matries (" << dim << " x " << dim << "): " << (taucs_ctime() - oldTime) << " seconds." << std::endl;

	
	ll_symbolic_factor_solve(SA, SMA, B, X, logger);	
	lu_symbolic_factor_solve(A, MA, B, X, logger);//u_symbolic_factor_solve(A, MA, B, X, logger);
#ifdef USE_OPENNL
	ooc_lu_factor_solve(MA, B, X, logger);//Ã»·´Ó¦, when use ooc_lu in Taucs of mt and mtd version.
#endif

	//create_write_matrix();//checked
	//test2();//checked
	//taucs_pseudo();// cause:A'A is always not symmetric positive definite
}
int testMeshMatrix(char* filename, char* solverType, std::ofstream &logger)
{
	double oldTime = taucs_ctime();

	logger << std::endl << "## Begin testMeshMatrix: " << oldTime << std::endl;
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
		   << taucs_ctime() - oldTime << " seconds" << std::endl;
	
	std::list<ConstrianedVertex*> vcMap;
	vcMap.push_back(new ConstrianedVertex(mesh.vertices_begin(),0));

	std::vector <double> rValues;
	rValues.reserve(mesh.size_of_vertices());
	for ( int i = 0; i<mesh.size_of_vertices(); ++i)
	{
		rValues.push_back(0.9);
	}

	Calculator calc(5);//0 MVC, 1 DCP, 5 Tutte
	calc.compute(&mesh, vcMap,rValues,solverType,logger);

	logger << "End testMeshMatrix: " << taucs_ctime() - oldTime << std::endl;
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
#ifdef _DEBUG
	logger << "Begin main() in Debug mode.The solver is : " << argv[2] << ". The cpu time is:" << taucs_ctime() << std::endl;
#else
	logger << "Begin main() in Release mode.The solver is : " << argv[2] << ". The cpu time is:" << taucs_ctime() << std::endl;
#endif

	int dim = testMeshMatrix(argv[1], argv[2], logger);
	//testSimpleMatrix(dim, logger);
	logger << std::endl << "End main(), time: " << taucs_ctime() << std::endl;
	logger.close();

	return EXIT_SUCCESS;
}

//void taucs_pseudo()
//{ 
//	taucsType* A_tB, *B, *X;
//
//	int i, j, m, n, val;
//
//	n = 4; m = 4; val  = 0;
//
//	taucs_ccs_matrix *A, *A_t, *A_tA;
//	std::vector<std::map<int, taucsType> > rowsC(n);
//
//
//	//fill a (sparse) column matrix rowsC(nxm) with values
//
//	for (i=0;i<m; i++)	
//		for( j=0;j<n; j++ )	  
//			rowsC[j][i] = i;			
//
//
//	//call taucs matrix creation and create the symmetric squared matrix A'A
//	A = CreateTaucsMatrixFromColumns(rowsC,m,TAUCS_DOUBLE);
//	i=taucs_ccs_write_ijv(A, "A.txt" );
//
//	A_t = MatrixTranspose(A);
//	i=taucs_ccs_write_ijv(A_t, "A_t.txt" );
//	A_tA = Mul2NonSymmMatSymmResult(A_t, A);
//	i=taucs_ccs_write_ijv(A_tA, "A_tA.txt" );
//
//	taucs_ccs_free(A);
//
//	B = new taucsType[m];
//	X = new taucsType[n];
//	A_tB = new taucsType[n];
//
//	//fill the right part of the equation AX = B;
//	for (i = 0; i<n; i++)		
//		B[i] = 0;//constrain;	
//
//
//	//create A'B
//	MulNonSymmMatrixVector(A_t, B, A_tB);		
//
//	char* options[] = { "taucs.factor.LU=true", NULL };
//
//	//solve A'Ax = A'B;
//	i = taucs_linsolve(A_tA, NULL, 1, X, A_tB, options, NULL);
//	if (i != TAUCS_SUCCESS)
//	{
//		printf ("Solution error.\n");
//		if (i==TAUCS_ERROR)
//			printf ("Generic error.");
//		if (i==TAUCS_ERROR_NOMEM)
//			printf ("NOMEM error.");    
//		if (i==TAUCS_ERROR_BADARGS)
//			printf ("BADARGS error.");   
//		if (i==TAUCS_ERROR_MAXDEPTH)
//			printf ("MAXDEPTH error.");    
//		if (i==TAUCS_ERROR_INDEFINITE)
//			printf ("NOT POSITIVE DEFINITE error.");         
//	}
//	else
//	{
//		printf ("Solution success.\n");	 
//		for (i = 0; i < n; i++)
//			printf ("%f\n",X[i]);      
//	}
//	taucs_ccs_free(A_t);
//	taucs_ccs_free(A_tA);
//
//	delete [] B;
//	delete [] X;
//	delete [] A_tB;
//}
//void create_write_matrix()
//{
//	int m = 4,n = 4,nnz = 8, i;
//
//	taucs_ccs_matrix * pMatrix = taucs_ccs_create( m, n, nnz, TAUCS_DOUBLE );
//	typedef double T;
//	//T* taucs_values = (T*)pMatrix->values.v;
//
//	pMatrix->colptr[0] = 0;
//	pMatrix->colptr[1] = 4;
//	pMatrix->colptr[2] = 8;
//	for (int j=0; j<n; j++) 
//	{
//		for (int i = pMatrix->colptr[j]; i < pMatrix->colptr[j+1]; i++) 
//		{
//			pMatrix->rowind[i] = j;
//			pMatrix->taucs_values[i] = i;//Aij
//		}		 
//	}
//
//	//for ( i = 0; i < nnz; i++ )
//	//{
//	//	pMatrix->rowind[i] = i;
//	//	pMatrix->taucs_values[i] = i;//taucs_values
//	//}
//	i = taucs_ccs_write_ijv(pMatrix, "test.txt" );
//	taucs_dccs_free(pMatrix);
//}
//void test2()
//{
//	int m = 4,n = 4,nnz = 7, i, used;
//	taucs_double x[4];  
//	taucs_double b[4];  
//	taucs_ccs_matrix * pMatrix;  
//	taucs_logfile ("stdout");
//	pMatrix = taucs_ccs_create( m, n, nnz, TAUCS_DOUBLE );
//	pMatrix->colptr[0] = 0;
//	pMatrix->colptr[1] = 2;  
//	pMatrix->colptr[2] = 4;  
//	pMatrix->colptr[3] = 6;  
//	pMatrix->colptr[4] = 7;    
//	used = 0;
//	for (i = 0; i < 4; i++ )    
//	{
//		pMatrix->rowind[used] = i;
//		pMatrix->taucs_values[used] = 1.0;
//		used++;
//		if (i+1 <= 3)
//		{
//			pMatrix->rowind[used] = i+1;
//			pMatrix->taucs_values[used] = 0.5;
//			used++;
//		}  
//	}
//	pMatrix->flags += TAUCS_SYMMETRIC;  
//	pMatrix->flags += TAUCS_LOWER;
//	i = taucs_ccs_write_ijv(pMatrix, "test.txt" );        
//
//	b[0] = 1.0; b[1] = 2.0; b[2] = 3.0; b[3] = 4.0;
//	char* options[] = { "taucs.factor.LLT=true","taucs.factor.ordering=metis",NULL };    
//	i = taucs_linsolve (pMatrix,NULL,1,x,b,options,NULL);  
//	if (i != TAUCS_SUCCESS)
//	{
//		printf ("Solution error.\n");
//		if (i==TAUCS_ERROR)
//			printf ("Generic error.");
//		if (i==TAUCS_ERROR_NOMEM)
//			printf ("NOMEM error.");    
//		if (i==TAUCS_ERROR_BADARGS)
//			printf ("BADARGS error.");   
//		if (i==TAUCS_ERROR_MAXDEPTH)
//			printf ("MAXDEPTH error.");    
//		if (i==TAUCS_ERROR_INDEFINITE)
//			printf ("NOT POSITIVE DEFINITE error.");         
//	}
//	else
//	{
//		printf ("Solution success.\n");	 
//		for (i = 0; i < 4; i++)
//			printf ("%f\n",x[i]);      
//	}
//	taucs_dccs_free(pMatrix);
//}