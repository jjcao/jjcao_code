#include "qttaucs.h"
#pragma warning(disable:4244 4267 4819 4996)

#include <CGAL/Cartesian.h>

#include <CGAL/Taucs_solver_traits.h>
typedef CGAL::Taucs_solver_traits<double>     Solver;
typedef CGAL::Taucs_symmetric_solver_traits<double>     SymSolver;

#include <sstream>
#include <fstream>
#include <string>

typedef Solver::Vector      Vector;
typedef Solver::Matrix      Matrix;
typedef SymSolver::Vector      SymVector;
typedef SymSolver::Matrix      SymMatrix;

void solve_SymbolicFactor(std::ofstream& logger)
{
	int status(0);
	logger << std::endl << "## Begin solve_SymbolicFactor" << std::endl;
	CGAL::Timer timer;
	timer.start();

	///////////////////////////////////////
	///////////////////////////////////////
	int dim(4);
	SymMatrix A(dim, dim);
	Matrix MA(dim, dim);
	SymVector X(dim), B(dim);
	for( int i = 0; i < dim; ++i)
	{
		A.set_coef(i, i, 1);
		MA.set_coef(i, i,1);
		B[i] = 2;
	}
	MA.set_coef(1, 3, 0); MA.set_coef(3, 1, -0.5);//-0.25
	//MA.set_coef(1, 4, -0.25); MA.set_coef(4, 1, -0.25);
	//MA.set_coef(1, 8, -0.25); MA.set_coef(8, 1, -0.25);
	//MA.set_coef(1, 99, -0.25);MA.set_coef(99, 1, -0.25);
	A.set_coef(1, 3, -0.25); A.set_coef(3, 1, -0.25);
	//A.set_coef(1, 4, -0.25); A.set_coef(4, 1, -0.25);
	//A.set_coef(1, 8, -0.25); A.set_coef(8, 1, -0.25);
	//A.set_coef(1, 99, -0.25);A.set_coef(99, 1, -0.25);
	//A.set_coef(4, 13, -0.5); A.set_coef(13, 4, -0.5);
	//A.set_coef(4, 199, -0.5);A.set_coef(199, 4, -0.5);
	//A.set_coef(81, 23, -0.2); A.set_coef(23, 81, -0.2);
	//A.set_coef(81, 4, -0.2); A.set_coef(4, 81, -0.2);
	//A.set_coef(81, 38, -0.2); A.set_coef(38, 81, -0.2);
	//A.set_coef(81, 299, -0.2);A.set_coef(299, 81, -0.2);
	//A.set_coef(81, 55, -0.2);A.set_coef(55, 81, -0.2);
	logger << "  matrix filling (" << dim << " x " << dim << "): " << timer.time() << " seconds." << std::endl;
	timer.reset();
	int i=taucs_ccs_write_ijv((taucs_ccs_matrix*)A.get_taucs_matrix(), "A.txt" );
	i=taucs_ccs_write_ijv((taucs_ccs_matrix*)MA.get_taucs_matrix(), "B.txt" );

	///////////////////////////////////////
	///////////////////////////////////////
	void* L = taucs_ccs_factor_llt_symbolic((taucs_ccs_matrix*)A.get_taucs_matrix());
	i   = taucs_ccs_factor_llt_numeric((taucs_ccs_matrix*)MA.get_taucs_matrix(),L);
	if (i)
	{
		logger << "Failed to numerical factorization!" << std::endl;
		status = -1;
	}
	else
	{
		status = taucs_supernodal_solve_llt(L,X.get_taucs_vector(),B.get_taucs_vector());
		taucs_supernodal_factor_free_numeric(L);
	}	

	///////////////////////////////////////
	///////////////////////////////////////
	if ( status==0)
	{
		logger << "  solving a linear systems: " << timer.time() << " seconds." << std::endl;
		timer.reset();
		for( int i = 0; i < dim; ++i)
		{
			logger << X[i] << ", ";
		}
	}
	logger << std::endl;
	logger << "End solve_SymbolicFactor!" << std::endl << std::endl;
}
qtTaucs::qtTaucs(QWidget *parent, Qt::WFlags flags)
	: QMainWindow(parent, flags)
{
	ui.setupUi(this);
}

qtTaucs::~qtTaucs()
{

}


void qtTaucs::on_actionTest_triggered()
{
	std::ofstream logger("log.txt");	

	solve_SymbolicFactor(logger);	

	logger.close();
}