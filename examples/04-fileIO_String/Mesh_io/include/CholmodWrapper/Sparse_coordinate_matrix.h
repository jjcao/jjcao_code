// Copyright (c) 2006  Hong Kong University of Science & Technology
// All rights reserved.
//
// This file is part of CholmodWrapper 
// (http://ihome.ust.hk/~fuhb/software.htm#cholmod_wrapper); 
// you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CholmodWrapper.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// Author(s) : Hongbo Fu, fuplus@gmail.com, 2006

// =====================================================================	
#ifndef SPARSE_COORDINATE_MATRIX_H
#define SPARSE_COORDINATE_MATRIX_H

#include <list>
#include <algorithm>

/// The class Sparse_coordinate_matrix is a sparse matrix representation 
/// in coordinate storage format (i, j, value). 
/// It supports duplicate entries, which will be summed finally.
///
/// Concept: Model of the SparseLinearAlgebraTraits_d::Matrix concept.
template <class T>
class Sparse_coordinate_matrix
{
	// -- private types ----------------------------
public:
	struct Triplet_entry ///< triplet: (i, j, value)
	{
		int i, j; // row and column index
		T value;
		Triplet_entry() : i(0), j(0), value(T()) {}
		Triplet_entry(int r, int c, T val)
			: i(r), j(c), value(val) {}
		void negate(){value = -value;}
	};
	/// Predicate: determine if the coordinates of two triplets are the same
	struct Triplet_entry_zero_element_predicate
	{		
		bool operator () (const Triplet_entry& right)
		{
			if (right.value == 0) return true;
			else return false;
		}
	};	
	typedef std::list<Triplet_entry>						Triplets_container;

	// -- public types -----------------------------
public:
	typedef T												NT;		
	typedef typename Triplets_container::iterator			Coordinate_iterator;
	typedef typename Triplets_container::const_iterator		Coordinate_const_iterator;

	// -- public operations ------------------------
public:
	/// Create a square matrix initialized with zeros.
	Sparse_coordinate_matrix(int dim, bool symmetric = false)		
		: m_symmetric (symmetric)
	{
		assert(dim > 0);
		m_row_dimension     = dim;
		m_column_dimension  = dim;
	}

	/// Create a rectangular matrix initialized with zeros.
	Sparse_coordinate_matrix(
		int  rows,                 ///< Matrix dimensions.
		int  columns, bool symmetric = false)		
		: m_symmetric (symmetric)
	{
		assert ( rows > 0 && columns > 0 );		
		if(symmetric) assert (rows == columns);	// must be a square matrix

		m_row_dimension     = rows;
		m_column_dimension  = columns;		
	}

	void setDimensions(int rows, int cols){m_row_dimension = rows; m_column_dimension = cols;}
	/// Return the matrix number of rows
	int row_dimension() const    { return m_row_dimension; }
	/// Return the matrix number of columns
	int column_dimension() const { return m_column_dimension; }

	/// Is a symmetric matrix?
	bool is_symmetric() const { return m_symmetric; }

	/// Write access to a matrix coefficient: a_ij <- val.
	/// 
	/// @note The efficiency of \c set_coef is much slower than that of \c add_coef
	/// So I strongly recommend using add_coef instead of set_coef.
	///
	/// Preconditions:
	/// - 0 <= i < row_dimension.
	/// - 0 <= j < column_dimension.
	void set_coef(int i, int j, T val) 
	{		
		reset_coef_to_zero(i, j);
		add_coef(i, j, val);		
	}

	/// Write access to a matrix coefficient: 
	/// if(exist(a_ij)) 
	///		a_ij <- a_ij + val;		
	/// else 
	///		a_ij <- val; 
	///
	/// Optimization for symmetric matrices:
	/// Sparse_coordinate_matrix stores only the lower triangle
	/// add_coef() does nothing if (i, j) belongs to the upper triangle.
	/// 
	/// Preconditions:
	/// - 0 <= i < row_dimension.
	/// - 0 <= j < column_dimension.	
	void add_coef(int i, int j, T val) 
	{
		assert ( i < m_row_dimension && j < m_column_dimension );
		
		if (val == 0) return; // do nothing
		if (m_symmetric && j > i) return; // do nothing		
		m_coordinates.push_back( Triplet_entry(i, j, val) );
	}

	/// Return an upper bound of number of non-zero elements.
	/// Since we allow duplicate entries, it is not efficient to directly count nnz.
	size_t upper_bound_of_nnz() const { return m_coordinates.size(); }

	/// Iterators for traverse the list of triplets
	Coordinate_iterator coords_begin() { return m_coordinates.begin(); }
	Coordinate_const_iterator coords_begin() const { return m_coordinates.begin(); }
	Coordinate_iterator coords_end() { return m_coordinates.end(); }
	Coordinate_const_iterator coords_end() const { return m_coordinates.end(); }

	void remove_bogus_nonzero_entries()
	{
 		m_coordinates.erase(
			// iterator of begin
			std::remove_if(
			m_coordinates.begin(), m_coordinates.end(),
			Triplet_entry_zero_element_predicate()),
			// iterator of end
 			m_coordinates.end());
	}


	void copy(Sparse_coordinate_matrix* in)
	{
		m_coordinates.clear();
		for (Coordinate_iterator it = in->coords_begin(); it != in->coords_end(); ++it)
			m_coordinates.push_back( Triplet_entry(it->i, it->j, it->value) );
	}

	// -- private operations -----------------------
private:
	/// set val(i, j) = 0;
	void reset_coef_to_zero(int i, int j)
	{
		assert ( i < m_row_dimension && j < m_column_dimension );
		// simulating removing effect by setting their values to zero
		for (Coordinate_iterator it = coords_begin(); it != coords_end(); ++it)
			if (it->i == i && it->j == j)
				it->value = 0;
	}	

	// -- private variables ------------------------
private:
	// Matrix dimensions
	int						m_row_dimension;
	int						m_column_dimension;
	Triplets_container		m_coordinates;
	bool					m_symmetric;
};

template <class T>
std::ostream& operator << (std::ostream& out, const Sparse_coordinate_matrix<T>& matrix)
{
	for (Sparse_coordinate_matrix<T>::Coordinate_const_iterator it = matrix.coords_begin(); it != matrix.coords_end(); ++it)
	{
		out << it->i << ", " << it->j << ", " << it->value << std::endl;
	}
	return out;
}

#endif // SPARSE_COORDINATE_MATRIX_H
// =====================================================================	