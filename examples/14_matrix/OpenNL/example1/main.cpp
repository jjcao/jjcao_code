#pragma warning(disable:4244 4267 4819 4996)

#include <CGAL/Cartesian.h>
#include <CGAL/Timer.h>
#include <sstream>
#include <fstream>
#include <string>

#include <OpenNL/linear_solver.h>
//typedef OpenNL::DefaultLinearSolverTraits<double> Solver;
typedef OpenNL::SymmetricLinearSolverTraits<double> Solver;
typedef Solver::Vector      Vector;
typedef Solver::Matrix      SparseMatrix;
typedef OpenNL::FullVector<double> FullVector;
typedef OpenNL::DefaultLinearSolverTraits< SparseMatrix, FullVector> SolverTraits ;
typedef OpenNL::LinearSolver<SolverTraits> Ls_solver ;

void solve(std::ofstream& logger)
{
	int status(0);
	CGAL::Timer timer;
	timer.start();

	Solver solver = Solver();
	int dim(999);
	SparseMatrix A(dim, dim);
	Vector X(dim), B(dim);

	///////////////////////////////////////
	///////////////////////////////////////
	logger << std::endl << "## Begin solve" << std::endl;
	for( int i = 0; i < dim; ++i)
	{
		A.set_coef(i, i, 1);
		B[i] = 2;
	}
	A.set_coef(1, 3, -0.25); A.set_coef(3, 1, -0.25);
	A.set_coef(1, 4, -0.25); A.set_coef(4, 1, -0.25);
	A.set_coef(1, 8, -0.25); A.set_coef(8, 1, -0.25);
	A.set_coef(1, 99, -0.25);A.set_coef(99, 1, -0.25);
	logger << "  matrix filling (" << dim << " x " << dim << "): " << timer.time() << " seconds." << std::endl;
	timer.reset();

	///////////////////////////////////////
	///////////////////////////////////////
	// Solve "A*X = B". On success, solution is (1/D) * X.
	double D;
	if ( !solver.linear_solver(A, B, X, D))
	{
		status = -1;//Base::ERROR_CANNOT_SOLVE_LINEAR_SYSTEM;
		logger << "  Error in solving a linear systems!" << std::endl;
	}
	else{
		logger << "  solving a linear systems: " << timer.time() << " seconds." << std::endl;
		timer.reset();
	}

	///////////////////////////////////////
	///////////////////////////////////////
	if ( status==0)
	{
		for( int i = 0; i < dim; ++i)
		{
			logger << X[i] << ", ";
		}
	}
	logger << std::endl;
	logger << "End solve!" << std::endl << std::endl;
}
void least_square()
{
	//Ls_solver solver(6);
	//solver.set_least_squares(true) ;
	//solver.begin_system() ;

	//solver.variable(0).set_value(1);
	//solver.variable(0).lock();
	//solver.variable(1).set_value(0);
	//solver.variable(2).set_value(0);
	//solver.variable(3).set_value(2);
	//solver.variable(3).lock();
	//solver.variable(4).set_value(0);
	//solver.variable(5).set_value(0);

	//solver.begin_row() ;
 //   solver.add_coefficient(0, 1)  ;
 //   solver_.add_coefficient(1,  b-d)  ;
 //   solver_.add_coefficient(u1_id,   -c)  ;
 //   solver_.add_coefficient(v1_id,    d)  ;
 //   solver_.add_coefficient(u2_id,    a) ;
 //   solver.end_row() ;

 //   setup_lscm() ;
 //   solver.end_system() ;
 //   std::cout << "Solving ..." << std::endl ;
 //   solver.solve() ;
}
int main(int argc, char* argv[])
{	
	std::ofstream logger("log.txt");
	std::cout.rdbuf(logger.rdbuf());

	solve(logger);	

	logger.close();
	return 0;
}