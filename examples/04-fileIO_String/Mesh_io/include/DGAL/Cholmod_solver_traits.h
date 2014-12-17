#ifndef DGAL_CHOLMOD_SOLVER_TRAITS_H
#define DGAL_CHOLMOD_SOLVER_TRAITS_H

#include <DGAL/config.h>
#include <cassert>
#include <CholmodWrapper/Cholmod_dense_vector.h>
#include <CholmodWrapper/Cholmod_dense_matrix.h>
#include <CholmodWrapper/Sparse_coordinate_matrix.h>
#include <CholmodWrapper/Cholmod_conversion.h>

DGAL_BEGIN_NAMESPACE

/// The class Cholmod_symmetric_solver_traits
/// is a traits class for solving GENERAL (aka unsymmetric) sparse linear systems
/// using SuiteSparse cholmod solvers family.
/// 
/// Cholmod_solver_traits: Model of the SparseLinearAlgebraTraits_d concept.
///						   Model of the Direct_linear_solver concept.
/// Sparse_matrix: Model of the SparseLinearAlgebraTraits_d::Matrix concept.
/// Dense_vector: Model of the SparseLinearAlgebraTraits_d::Vector concept.
template < 
	class T, 					///< Currently only double is tested.	
	class Dense_vector2   =		Cholmod_dense_vector<T>,
	class Dense_matrix2   =		Cholmod_dense_matrix<T>,
	class Sparse_matrix2  =		Sparse_coordinate_matrix<T>
>
class Cholmod_solver_traits
{
	// -- public types -----------------------------
public:
	typedef T								NT;
	// for SparseLinearAlgebraTraits_d concept
	typedef Sparse_matrix2					Matrix;
	typedef Dense_vector2					Vector;

	typedef Sparse_matrix2					Sparse_matrix;	
	typedef Dense_vector2					Dense_vector;
	typedef Dense_matrix2					Dense_matrix;

	// -- public operations ------------------------
public:
	Cholmod_solver_traits(bool verbose = false)
		: m_verbose (verbose)
		, m_symmetric (false)
		, m_num_of_unknows (0)
		, m_L (0)
		, m_At (0)
	{
		cholmod_start(&m_cholmod_common);
	}

	~Cholmod_solver_traits()
	{
		release_precomputation();
		cholmod_finish(&m_cholmod_common);
	}

	/// Solve the sparse linear system "A * X = B" in the least-squares sense.
	/// @param A the system matrix
	/// @param B the right hand side of the system. It can be either Dense_vector or Dense_matrix.
	/// @param X the final solution.
	/// @param D last entry of a homogeneous coordinate. Useless for Cholmod_solver_traits.
	/// @return true on success. 
	template <class Dense_matrix_T>
	bool linear_solver(
		const Sparse_matrix& A, const Dense_matrix_T& B, Dense_matrix_T& X, NT& D, cholmod_factor* factor=0)
	{
		D = 1; // CHOLMOD does not support homogeneous coordinates

		if (!has_precomputed())
		{
			if (!precompute(A, factor) && m_verbose)  // perform Cholesky factorization 
			{
				std::cout << "FATAL ERROR: cannot PRECOMPUTE the factorization.\n";
				return false;
			}
		}
 
		if (!solve(B, X, factor) && m_verbose)	 // solve using back-substituion
		{
			std::cout << "FATAL ERROR: cannot SOLVE the system after factorization.\n";
			return false;
		}
		return true;
	}

	template <class Dense_matrix_T>
	bool linear_solver(
		const Sparse_matrix& A, const Dense_matrix_T& B, Dense_matrix_T& X, cholmod_factor* factor=0)
	{
		NT D;
		return linear_solver(A, B, X, D, factor);
	}

	/// Peform precomputation (Cholesky decomposition) given a system matrix A.
	/// If A.is_symmetric(), A = L*L' is performed.
	/// Otherwise, A'*A = L*L' is performed. 
	/// @return true on success. 
	bool precompute(const Sparse_matrix& A, cholmod_factor* factor=0)
	{
		release_precomputation(); // release old stuffs

		m_symmetric = A.is_symmetric();


		// If A is symmetric, 
		// else A*A' is analyzed.
		if (m_symmetric) // A is analyzed.
		{
			// for symmetric matrices, At = A
			m_At = create_cholmod_sparse(A, &m_cholmod_common);			
		}
		else			// A*A' is analyzed: we need A'*A however.
		{
#ifdef OUTPUT_INFO
			CGAL::Timer timer;	timer.start();
#endif		
			cholmod_sparse* tmp_A = create_cholmod_sparse(A, &m_cholmod_common);
			assert (tmp_A != 0);
#ifdef OUTPUT_INFO
			std::cout << "matrix create: " << timer.time() << " seconds." << std::endl;
			timer.reset();
#endif
			m_At = cholmod_transpose(tmp_A, 1 /* array transpose */, &m_cholmod_common);
#ifdef OUTPUT_INFO
			std::cout << "matrix multiply: " << timer.time() << " seconds." << std::endl;
#endif			
			cholmod_free_sparse(&tmp_A, &m_cholmod_common);
		}
		
		assert (m_At != 0);
		if (!factor)
		{
#ifdef OUTPUT_INFO
			CGAL::Timer timer;	timer.start();
#endif	
			m_L = cholmod_analyze(m_At, &m_cholmod_common);
			if (m_L == 0) return false;

			cholmod_factorize(m_At, m_L, &m_cholmod_common);
			
#ifdef OUTPUT_INFO
			std::cout << "cholmod_factorize: " << timer.time() << " seconds." << std::endl;
#endif			
			/*cholmod_sparse * sm = cholmod_factor_to_sparse(m_L, &(m_cholmod_common));
			cholmod_triplet* triplet = cholmod_sparse_to_triplet(sm, &(m_cholmod_common));
			int* Ti		= static_cast<int*> (triplet->i);
			int* Tj		= static_cast<int*> (triplet->j);
			double* Tx = static_cast<double*>(triplet->x); 
			for(int i = 0; i < triplet->nzmax; ++i)
			{
				std::cout << Ti[i] <<", "<< Tj[i] <<", "<< Tx[i]<< std::endl;
			}*/
		}
		// number of unknowns
		m_num_of_unknows = static_cast<unsigned int>(m_At->nrow); // <==> m_A->mcol;

		return true;
	}

	/// Solve the system using back-substitution A X = B. 	
	/// X and B can be either Dense_vector or Dense_matrix. 
	/// When X and B are in Dense_matrix, it means solve 
	/// A [x1 x2 ... xn] = [b1 b2 ... bn] simultaneously. 
	/// Precondition: the system has been precomputed.
	template <class Dense_matrix_T>
	bool solve(const Dense_matrix_T& B, Dense_matrix_T& X, cholmod_factor* factor=0)
	{		
		if (!has_precomputed() && !factor) return false;

		cholmod_dense* b = B.create_cholmod_dense(&m_cholmod_common);
		cholmod_dense* x = 0;

		if (factor)
			m_L = factor;

		if (m_symmetric)	// A x = b
			x = cholmod_solve(CHOLMOD_A, m_L, b, &m_cholmod_common);
		else				// A'A x = A'b
		{
			cholmod_dense* Atb = cholmod_zeros(
				m_At->nrow, b->ncol, CHOLMOD_REAL, &m_cholmod_common);	
			double alpha[2] = {1, 0}, beta[2] = {0, 0};
			cholmod_sdmult(m_At, 0, alpha, beta, b, Atb, &m_cholmod_common);
//////////////////////////////////////////////////////////////////////////
//{			
//	FILE* file;
//	fopen_s(&file, "m_At.txt","w");
//	if(!file){
//		std::cout << "save matrix failed: m_At.txt" << std::endl;
//		return false;
//	}
//	cholmod_write_sparse(file, m_At, 0, 0, &m_cholmod_common);			
//	fclose(file);
//}
//{			
//	FILE* file;
//	fopen_s(&file, "Atb.txt","w");
//	if(!file){
//		std::cout << "save matrix failed: Atb.txt" << std::endl;		
//		return false;
//	}
//	cholmod_write_dense(file, Atb, 0, &m_cholmod_common);			
//	fclose(file);
//}
//{
//	std::cout << "matrix factor: " << std::endl;
//	cholmod_sparse * sm = cholmod_factor_to_sparse(m_L, &(m_cholmod_common));
//	cholmod_triplet* triplet = cholmod_sparse_to_triplet(sm, &(m_cholmod_common));
//	int* Ti		= static_cast<int*> (triplet->i);
//	int* Tj		= static_cast<int*> (triplet->j);
//	double* Tx = static_cast<double*>(triplet->x); 
//	for(int i = 0; i < triplet->nzmax; ++i)
//	{
//		std::cout << Ti[i] <<", "<< Tj[i] <<", "<< Tx[i]<< std::endl;
//	}
//}
////////////////////////////////////////////////////////////////////
#ifdef OUTPUT_INFO
CGAL::Timer timer;	timer.start();
#endif
			x = cholmod_solve(CHOLMOD_A, m_L, Atb, &m_cholmod_common);
#ifdef OUTPUT_INFO
std::cout << "cholmod_solve: " << timer.time() << " seconds." << std::endl;
#endif
		}	

		assert (x != 0);
		X.assign(x);

		if (factor)
			m_L = 0;

		return true;
	}

	/// Has the system precomputed? 
	bool has_precomputed() const { return m_L != 0; }
	void factor(cholmod_factor* factor){m_L = factor;}
	void matrix(cholmod_sparse* matrix){m_At = matrix;}
	cholmod_factor* factor(){return m_L;}
	cholmod_sparse* matrix(){return m_At;}

	// -- private operations -----------------------
private:
	void release_precomputation()
	{
		// it is safe to release NULL pointers in cholmod.
		cholmod_free_factor(&m_L,	&m_cholmod_common);
		cholmod_free_sparse(&m_At,	&m_cholmod_common);
	}

	// -- private variables ------------------------
private:	
	cholmod_factor* m_L;
	cholmod_sparse* m_At;

	cholmod_common m_cholmod_common;
	bool m_symmetric;
	bool m_verbose;
	unsigned int m_num_of_unknows;
};

typedef Sparse_coordinate_matrix<Num_type_default> Sparse_matrix_default;
typedef Cholmod_solver_traits<Num_type_default> Solver_defualt;

DGAL_END_NAMESPACE

#endif // DGAL_CHOLMOD_SOLVER_TRAITS_H
// =====================================================================	