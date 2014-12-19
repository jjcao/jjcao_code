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
#ifndef CHOLMOD_DENSE_VECTOR_H
#define CHOLMOD_DENSE_VECTOR_H

#include <vector>
#include <cholmod.h>

/// class Cholmod_dense_vector is a simple array, derived from std::vector.
///
/// Concept: Model of the SparseLinearAlgebraTraits_d::Vector concept.
template <class T>
class Cholmod_dense_vector : public std::vector<T>
{	
private: 
	typedef std::vector<T> Super_class;
	
	// -- public types -----------------------------
public:
	typedef T NT;

	// -- public operations ------------------------
public:
	// Create a vector initialized with zeros
	Cholmod_dense_vector(int dimension)
		: std::vector<T>(dimension, 0.0)
	{
	}

	Cholmod_dense_vector(const std::vector<T>& vec)
	{
		resize(vec.size());
		std::copy(vec.begin(), vec.end(), this->begin());
	}
	Cholmod_dense_vector(const Cholmod_dense_vector<T>& vec)
	{
		resize(vec.size());
		std::copy(vec.begin(), vec.end(), this->begin());
	}

	void setDimension(int in){resize(in);}
	int dimension() const { return Super_class::size(); }
	inline void set(int i, T val) { (*this)[i] = val; }
	inline T get(int i) const { return (*this)[i]; }

	/// Convert to a cholmod_dense matrix	
	cholmod_dense* create_cholmod_dense(cholmod_common* c) const
	{
		cholmod_dense* output = cholmod_zeros(Super_class::size(), 1, CHOLMOD_REAL, c);
		std::copy(Super_class::begin(), Super_class::end(), (T*)output->x);

		return output;
	}

	/// Copy values of a single column or single row cholmod_dense matrix 
	/// to Cholmod_dense_vector
	void assign (const cholmod_dense* cd_matrix) 
	{
		assert (cd_matrix != 0);
		assert (
			(Super_class::size() == cd_matrix->nrow && cd_matrix->ncol == 1)	|| 
			(Super_class::size() == cd_matrix->ncol && cd_matrix->nrow == 1) );

		T* ptr = (T*)cd_matrix->x;
		std::copy(ptr, ptr + size(), Super_class::begin());
	}
};
template <class T>
std::ostream& operator << (std::ostream& out, const Cholmod_dense_vector<T>& vector)
{
	for (Cholmod_dense_vector<T>::const_iterator it = vector.begin(); it != vector.end(); ++it)
	{
		out << *it << std::endl;
	}
	return out;
}
#endif // CHOLMOD_DENSE_VECTOR_H
// =====================================================================	