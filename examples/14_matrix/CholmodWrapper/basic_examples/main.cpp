#include <iostream>

#include <CholmodWrapper/Cholmod_solver_traits.h>
#include <CholmodWrapper/Sparse_coordinate_matrix.h>


typedef Cholmod_solver_traits<double>	Solver;
typedef Solver::Sparse_matrix			Sparse_matrix;
typedef Solver::Dense_matrix			Dense_matrix;
bool test2()
{
	int m = 3, n = 2, c = 1;
	Sparse_matrix A(m, n);

	// initialize the system matrix
	//		  [ 1	1 ]
	// A =    [ 1	-1]
	//		  [ 1	2 ]
	// 		  
	A.add_coef(0, 0, 1);
	A.add_coef(0, 1, 1);
	A.add_coef(1, 0, 1);
	A.add_coef(1, 1, -1);
	A.add_coef(2, 0, 1);
	A.add_coef(2, 1, 2);

	// initialize the right hand side of the sytem
	Dense_matrix  B(m, c);
	//		  [2]
	// b =    [0]
	//		  [4]
	B.set(0, 0, 2);
	B.set(1, 0, 0);
	B.set(2, 0, 4);

	Solver solver;
	Dense_matrix X(n, c);

	if (!solver.precompute(A))	// perform Cholesky factorization
	{
		std::cerr << "Error: cannot precompute the system matrix.\n";		
	}
	if (!solver.solve(B, X))// solve by back-substitution
	{
		std::cerr << "Error: the system matrix has not been precomputed yet.\n";	
	}

	//output
	std::cout << "solve AX=B" <<std::endl;
	std::cout << "A is :" << std::endl << A;
	std::cout << "B is :" << std::endl << B;
	std::cout << "The solution X is :" << std::endl << X;
	//		  [ 1.14286]
	//x =     [ 1.28571]
	std::cin.get();

	return true;
}
bool test1()
{
	int m = 3, n = 2, c = 1;
	Sparse_matrix A(m, n);

	// initialize the system matrix
	//		  [ 1	1 ]
	// A =    [ 1	-1]
	//		  [ 1	2 ]
	// 		  
	A.add_coef(0, 0, 1);
	A.add_coef(0, 1, 1);
	A.add_coef(1, 0, 1);
	A.add_coef(1, 1, -1);
	A.add_coef(2, 0, 1);
	A.add_coef(2, 1, 2);

	// initialize the right hand side of the sytem
	Dense_matrix  B(m, c);
	//		  [2]
	// b =    [0]
	//		  [4]
	B.set(0, 0, 2);
	B.set(1, 0, 0);
	B.set(2, 0, 4);


	Solver solver;
	Dense_matrix X(n, c);

	if (!solver.linear_solver(A, B, X))
	{
		std::cerr << "Error: cannot solve the linear system.\n";
		return false;
	}

	//if (!solver.precompute(A))	// perform Cholesky factorization
	//{
	//	std::cerr << "Error: cannot precompute the system matrix.\n";		
	//}
	//if (!solver.solve(B, X))// solve by back-substitution
	//{
	//	std::cerr << "Error: the system matrix has not been precomputed yet.\n";	
	//}

	//output
	std::cout << "solve AX=B" <<std::endl;
	std::cout << "A is :" << std::endl << A;
	std::cout << "B is :" << std::endl << B;
	std::cout << "The solution X is :" << std::endl << X;
	//		  [ 1.14286]
	//x =     [ 1.28571]
   std::cin.get();


	///*****************cout****************************************/
	//ofstream output("result.txt");
	//for (int i=0;i<n;++i)
		//output <<X.get(i,0)<<'\n';
	//output.close();

   return true;
}
int main(int argc, char* argv[])
{
	test1();
	test2();
	return 0;
}
