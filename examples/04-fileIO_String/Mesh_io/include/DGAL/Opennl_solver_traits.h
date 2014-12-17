#ifndef DGAL_OPENNL_SOLVER_TRAITS_H
#define DGAL_OPENNL_SOLVER_TRAITS_H

#include <DGAL/config.h>
#include <CGAL/OpenNL/linear_solver.h>
#include <CholmodWrapper/Cholmod_dense_vector.h>
#include <CholmodWrapper/Cholmod_dense_matrix.h>
#include <CholmodWrapper/Sparse_coordinate_matrix.h>
#include <CholmodWrapper/Cholmod_conversion.h>

DGAL_BEGIN_NAMESPACE

template<class TRAITS> 
class Opennl_solver_traits : public OpenNL::DefaultLinearSolverTraits<TRAITS>
{
	// Public types
public:
	typedef typename TRAITS                       NT;
	typedef typename Cholmod_dense_matrix<NT>     Dense_matrix;
	// Public operations
public:
	void factor(cholmod_factor* factor){}
	void matrix(cholmod_sparse* matrix){}
	cholmod_factor* factor(){return 0;}
	cholmod_sparse* matrix(){return 0;}
	
	template <class Dense_matrix_T>
	bool linear_solver(const Matrix& A, const Dense_matrix_T& B, Dense_matrix_T& X, cholmod_factor* factor=0)
	{
		return false;
	}
	template <class Dense_matrix_T>
	bool solve(const Dense_matrix_T& B, Dense_matrix_T& X, cholmod_factor* factor=0){return false;}

};

DGAL_END_NAMESPACE

#endif // DGAL_OPENNL_SOLVER_TRAITS_H
