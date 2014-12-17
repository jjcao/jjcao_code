#ifndef DGAL_TAUCS_SOLVER_TRAITS_H
#define DGAL_TAUCS_SOLVER_TRAITS_H

#include <DGAL/config.h>
#include <CGAL/Taucs_solver_traits.h>
#include <CholmodWrapper/Cholmod_dense_vector.h>
#include <CholmodWrapper/Cholmod_dense_matrix.h>
#include <CholmodWrapper/Sparse_coordinate_matrix.h>
#include <CholmodWrapper/Cholmod_conversion.h>

DGAL_BEGIN_NAMESPACE
template<class T> 
class Taucs_solver_traits : public CGAL::Taucs_solver_traits<T>
{
// Public types
public:
	typedef typename CGAL::Taucs_solver_traits<T> Base;
    typedef typename CGAL::Taucs_matrix<T>             Matrix;
    typedef typename CGAL::Taucs_vector<T>             Vector;
	typedef Matrix            Sparse_matrix;
	typedef Vector             Vector;
    typedef T                           NT;

// Public operations
public:

	void factor(cholmod_factor* factor){}
	void matrix(cholmod_sparse* matrix){}
	cholmod_factor* factor(){return 0;}
	cholmod_sparse* matrix(){return 0;}
	template <class Dense_matrix_T>
	bool linear_solver(const Sparse_matrix& A, const Dense_matrix_T& B, Dense_matrix_T& X, cholmod_factor* factor=0)
	{
		return false;	
	}
	template <class Dense_matrix_T>
	bool solve(const Dense_matrix_T& B, Dense_matrix_T& X, cholmod_factor* factor=0){return false;}
};


DGAL_END_NAMESPACE

#endif // DGAL_TAUCS_SOLVER_TRAITS_H
