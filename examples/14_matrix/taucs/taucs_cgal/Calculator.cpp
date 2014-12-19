#include "Calculator.h"

Calculator::Calculator(int weightOption)
{
	switch (weightOption){
	case 0: m_weightOption = new MvcOption;	break;
	case 1: m_weightOption = new DcpOption;	break;
	case 3: m_weightOption = new SpringOption1;	break;
	case 4: m_weightOption = new SpringOption2;	break;
	default: m_weightOption = new TutteOption; 
	}
}
Calculator::~Calculator()
{
	if (m_weightOption)
	{
		delete m_weightOption;
		m_weightOption = 0;
	}
}
int Calculator::index_mesh_vertices(Polyhedron* mesh)
{
	int i(0);
	for (Vertex_iterator itor = mesh->vertices_begin();	itor != mesh->vertices_end();++itor)
	{
		itor->is_parameterized(false);
		itor->s(-1.0);
		itor->index(i);
		++i;
	}
	return i;
}
void Calculator::setScalarForVc(Polyhedron* mesh, std::list<ConstrianedVertex*>& vcMap)
{
	for (std::list<ConstrianedVertex*>::iterator it = vcMap.begin(); it != vcMap.end(); ++it)
	{
		ConstrianedVertex* vc = *it;
		Vertex_handle vh = vc->m_vh;
		vh->is_parameterized(true);
		vh->s(vc->m_scalar);
	}
}
void Calculator::initialize_system_from_mesh_vc(Matrix& A, Vector& Bu, std::list<ConstrianedVertex*>& vcMap)
{
	for (std::list<ConstrianedVertex*>::iterator it = vcMap.begin(); it != vcMap.end(); ++it)
	{
		ConstrianedVertex* vc = *it;
		int index = vc->m_vh->index();
		double scalar = vc->m_scalar;

		// Write a as diagonal coefficient of A
		A.set_coef(index, index, 1);

		Bu[index] = scalar;
	}
}
// if cValue, Poisson; else Laplacian
int Calculator::setup_inner_vertex_relations(Matrix& A, Vector& Bu, 
											Polyhedron* mesh, Vertex_handle vh,
											std::vector<double> &cValues)
{
	int i = vh->index();

	// circulate over vertices around 'vertex' to compute w_ii and w_ijs
	double w_ii = 0;
	int vertexIndex = 0;
	Halfedge_vertex_circulator v_j = vh->vertex_begin();
	Halfedge_vertex_circulator end = v_j;
	CGAL_For_all(v_j, end)
	{
		// Call to virtual method to do the actual coefficient computation
		double w_ij = -1.0 * m_weightOption->operator()(mesh, vh, v_j);

		// w_ii = - sum of w_ijs
		w_ii -= w_ij;

		// Get j index
		int j = v_j->opposite()->vertex()->index();

		// Set w_ij in matrix
		A.set_coef(i,j, w_ij);

		vertexIndex++;
	}
	if (vertexIndex < 2)
		return -1;//Base::ERROR_NON_TRIANGULAR_MESH;

	// Set w_ii in matrix
	A.set_coef(i,i, w_ii);
	if ( !cValues.empty())// for Poisson
	{
		Bu[i] = std::abs(cValues[i]);
	}
	return i;//for OK
}
int Calculator::setup_inner_vertex_relations(SymmetricMatrix& A, Vector& Bu, 
											Polyhedron* mesh, Vertex_handle vh,
											std::vector<double> &cValues)
{
	int i = vh->index();

	// circulate over vertices around 'vertex' to compute w_ii and w_ijs
	double w_ii = 0;
	int vertexIndex = 0;
	Halfedge_vertex_circulator v_j = vh->vertex_begin();
	Halfedge_vertex_circulator end = v_j;
	CGAL_For_all(v_j, end)
	{
		// Call to virtual method to do the actual coefficient computation
		double w_ij = -1.0 * m_weightOption->operator()(mesh, vh, v_j);

		// w_ii = - sum of w_ijs
		w_ii -= w_ij;

		// Get j index
		int j = v_j->opposite()->vertex()->index();

		// Set w_ij in matrix
		A.set_coef(i,j, w_ij);

		vertexIndex++;
	}
	if (vertexIndex < 2)
		return -1;//Base::ERROR_NON_TRIANGULAR_MESH;

	// Set w_ii in matrix
	A.set_coef(i,i, w_ii);
	if ( !cValues.empty())// for Poisson
	{
		Bu[i] = std::abs(cValues[i]);
	}
	return i;//for OK
}
void Calculator::set_mesh_scalar_from_system(Polyhedron* mesh, Vector& Xu)
{
	Vertex_iterator vertexIt;
	for (vertexIt = mesh->vertices_begin();
		vertexIt != mesh->vertices_end();
		++vertexIt)
	{
		int index = vertexIt->index();

		double s = Xu[index];

		vertexIt->s(s);
		vertexIt->is_parameterized(true);
	}
}
int Calculator::lu_symbolic_factor_solve(Matrix& A, Matrix& MA, Vector& B, Vector& X, std::ofstream& logger)
{;
	logger << std::endl << "#### Begin Calculator::lu_symbolic_factor_solve" << std::endl;
	logger << std::endl << "## Taucs in CGAL do not support lu_symbolic_factor!" << std::endl;
	logger << "End lu_symbolic_factor_solve" << std::endl;
	return -1;
}
int Calculator::ooc_lu_factor_solve(Matrix& A,Vector& B, Vector& X, std::ofstream& logger)
{
	int status(0);	
#ifdef USE_OPENNL
	std::string funcName = "Calculator::OpenNL::BICGSTAB";
#else
	std::string funcName = "Calculator::Taucs::ooc_lu_factor_solveTaucs";	
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
		logger << "  Error in solving a linear systems!" << std::endl;
	}
	else{
		logger << "  solving a linear systems: " << timer.time() << " seconds." << std::endl;
		timer.reset();
	}

	logger << "End " << funcName << std::endl;

	return status;
}
int Calculator::ll_symbolic_factor_solve(SymmetricMatrix &SA, SymmetricMatrix &SMA, Vector& B, Vector& X, std::ofstream& logger)
{
	int status(0);
	logger << std::endl << "#### Begin Calculator::ll_symbolic_factor_solve" << std::endl;

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
	
	logger << "....taucs_supernodal_solve_llt: " << timer.time() << " seconds." << std::endl;
	timer.reset();
	taucs_supernodal_factor_free_numeric(L);

#endif
	logger << "End ll_symbolic_factor_solve!" << std::endl << std::endl;
	return status;
}
void Calculator::setupSystem(Polyhedron* mesh, Matrix&A, Vector&Bu, std::list<ConstrianedVertex*>& vcMap, std::vector<double> &cValues)
{
	int status(-1);
	// Initialize A, Bu 
	// Fill the constrained vertices' lines in the linear systems:
	// "s = constant"
	initialize_system_from_mesh_vc(A, Bu, vcMap);

	// Fill the matrix for the inner vertices v_i: compute A's coefficient
	// w_ij for each neighbor j; then w_ii = - sum of w_ijs
	for (Vertex_iterator vertexIt = mesh->vertices_begin(); vertexIt != mesh->vertices_end(); ++vertexIt)
	{
		// inner vertices only
		if ( vertexIt->is_parameterized())
			continue;

		// Compute the line i of matrix A for i inner vertex
		status = setup_inner_vertex_relations(A, Bu, mesh,vertexIt,cValues);
	}	
}
void Calculator::setupSymmetricSystem(Polyhedron* mesh, SymmetricMatrix&A, Vector&Bu, std::list<ConstrianedVertex*>& vcMap, std::vector<double> &cValues)
{
	int status(-1);
	// Fill the matrix for the inner vertices v_i: compute A's coefficient
	// w_ij for each neighbor j; then w_ii = - sum of w_ijs
	for (Vertex_iterator vertexIt = mesh->vertices_begin(); vertexIt != mesh->vertices_end(); ++vertexIt)
	{
		// Compute the line i of matrix A for i inner vertex
		status = setup_inner_vertex_relations(A, Bu, mesh,vertexIt,cValues);
	}	
}
int Calculator::compute(Polyhedron* mesh,std::list<ConstrianedVertex*>& vcMap,
						std::vector<double> &cValues, char* solverType, std::ofstream& logger)
{
	int status(0);//0 for success

#ifdef DEBUG_TRACE
	logger << "  Calculator::compute() begin!" <<std::endl;
	// Create timer for traces
	CGAL::Timer timer;
	timer.start();
#endif
	// Count and index vertices
	int nbVertices = index_mesh_vertices(mesh);

	// Mark all vertices as NOT "parameterized"
    Vertex_iterator vertexIt;
    for (vertexIt = mesh->vertices_begin(); vertexIt != mesh->vertices_end();
        vertexIt++)
    {
		vertexIt->is_parameterized(false);
    }

	// set scalar for constrained vertices and mark them as "parameterized"
	setScalarForVc(mesh, vcMap);
#ifdef DEBUG_TRACE
	logger << "....set scalar for constrained vertices: " << timer.time() << " seconds." << std::endl;
	timer.reset();
#endif

	// Create a sparse linear systems "A*Xu = Bu" (one line/column per vertex)
	Matrix A(nbVertices, nbVertices);
	Vector Xu(nbVertices), Bu(nbVertices);
	setupSystem(mesh, A, Bu, vcMap, cValues);

#ifdef DEBUG_TRACE
	logger << "....matrix filling (" << nbVertices << " x " << nbVertices << "): "
		<< timer.time() << " seconds." << std::endl;
	timer.reset();
#endif

	SymmetricMatrix SA(nbVertices, nbVertices);
	Vector Xv(nbVertices), Bv(nbVertices);
	setupSymmetricSystem(mesh,SA, Bv, vcMap, cValues);	

#ifdef DEBUG_TRACE
	logger << "....symmetric matrix filling (" << nbVertices << " x " << nbVertices << "): "
		<< timer.time() << " seconds." << std::endl;
	timer.reset();
#endif

	std::string st(solverType);
	if( st.compare("Taucs::lu_symbolic")==0)
	{
		status = lu_symbolic_factor_solve(A, A, Bu, Xu, logger);
	}
	else if (st.compare("Taucs::ooc_lu")==0 || st.compare("OpenNL::BICGSTAB")==0)
	{
		status = ooc_lu_factor_solve(A, Bu, Xu, logger);
	}
	else if (st.compare("Taucs::ll_symbolic")==0)
	{
		status = ll_symbolic_factor_solve(SA, SA, Bv, Xv, logger);
	}
	else
	{
		logger << "!! not support!" << std::endl;
		status = -1;
	}

	if (status < 0)
		return status;


#ifdef DEBUG_TRACE
	timer.reset();
#endif	
	// Copy scalar into s of each vertex
	set_mesh_scalar_from_system (mesh, Xu);
#ifdef DEBUG_TRACE
	logger << "....copy computed UVs to mesh :" << timer.time() << " seconds." << std::endl;
	timer.reset();
#endif

	return status;
}