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
#ifndef CHOLMOD_DENSE_MATRIX_H	
#define CHOLMOD_DENSE_MATRIX_H

#include <vector>
#include <iostream>
#include <cholmod.h>

// =====================================================================	
#pragma warning (push)
#pragma warning (disable: 4018)

/// class Cholmod_dense_matrix is a simple array of array.
/// Recommend using set/get functions to access the elements.
template <class T>
class Cholmod_dense_matrix : public std::vector< std::vector<T> >
{	
	// -- public types -----------------------------
public:
	typedef T NT;

	// -- public operations ------------------------
public:
	// Create a vector initialized with zeros
	Cholmod_dense_matrix(int row, int col)
		: m_row_dimension(row), m_column_dimension(col)
	{
		assert (row > 0 && col > 0);
		resize(col, std::vector<T>(row, 0.0) );	// memory allocation
	}

 	/// Convert to a cholmod_dense matrix	
 	cholmod_dense* create_cholmod_dense(cholmod_common* c) const
 	{
 		cholmod_dense* output = 
			cholmod_zeros(m_row_dimension, m_column_dimension, CHOLMOD_REAL, c);
		for(int i = 0; i < output->nrow; ++i)
		{
			for(int j = 0; j < output->ncol; ++j)
				((T*)output->x)[i + j * output->d] = (*this)[j][i];
		}

 		return output;
 	}
 
 	/// Copy values of a single column or single row cholmod_dense 
	/// matrix to Cholmod_dense_matrix
 	void assign (const cholmod_dense* cd_matrix) 
 	{
 		assert (cd_matrix != 0);
		assert (cd_matrix->nrow == m_row_dimension && 
			cd_matrix->ncol == m_column_dimension);

		for(int i = 0; i < cd_matrix->nrow; ++i)
		{
			for(int j = 0; j < cd_matrix->ncol; ++j)
				(*this)[j][i] = ((T*)cd_matrix->x)[i + j * cd_matrix->d];
		}
 	}

	inline void set(int row, int col, T val) { (*this)[col][row] = val; }
	inline T& get(int row, int col) { return (*this)[col][row]; }
	inline const T& get(int row, int col) const { return (*this)[col][row]; }

	/// Return the matrix number of rows
	int row_dimension() const    { return m_row_dimension; }
	/// Return the matrix number of columns
	int column_dimension() const { return m_column_dimension; }

	// -- private variables ------------------------
private:
	int	m_row_dimension;
	int m_column_dimension;
};

template <class T>
std::ostream& operator << (std::ostream& out, const Cholmod_dense_matrix<T>& matrix)
{
	for (int i = 0; i < matrix.row_dimension(); ++i)
	{
		out << "Row " << i << ":";
		for (int j = 0; j < matrix.column_dimension(); ++j)
			out << " " << matrix.get(i, j);
		out << "\n";
	}
	return out;
}


#pragma warning(pop) 
// =====================================================================	

#endif // CHOLMOD_DENSE_MATRIX_H
// =====================================================================	