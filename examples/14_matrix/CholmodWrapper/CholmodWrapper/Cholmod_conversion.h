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
#ifndef CHOLMOD_CONVERSION_H
#define CHOLMOD_CONVERSION_H

#include <cholmod.h>
#include "Sparse_coordinate_matrix.h"

/// Convert sparse matrix from coordinate format to cholmod_triplet.
/// Mainly for internal use.
template <class T>
cholmod_triplet* create_cholmod_triplet(
	const Sparse_coordinate_matrix<T>& coor_matrix, cholmod_common* c)
{
	Sparse_coordinate_matrix<T> local_copy(coor_matrix);
	local_copy.remove_bogus_nonzero_entries();
	size_t upper_bound_of_nnz = local_copy.upper_bound_of_nnz();

	// stype = 0: cholmod_triplet is "unsymmetric": 
	//			  use both upper and lower triangular parts
	// stype > 0: matrix is square and symmetric. 
	//			  use the lower triangular part to store non-zero entries
	int stype = local_copy.is_symmetric() ? 1 : 0;

	cholmod_triplet* triplet_matrix = cholmod_allocate_triplet(
		local_copy.row_dimension(), local_copy.column_dimension(), 
		upper_bound_of_nnz, stype, CHOLMOD_REAL, c);
	triplet_matrix->nnz = upper_bound_of_nnz;

	if (upper_bound_of_nnz == 0) return triplet_matrix; // empty

	int* Ti		= static_cast<int*>	(triplet_matrix->i);
	int* Tj		= static_cast<int*> (triplet_matrix->j);
	T* Tx		= static_cast<T*>	(triplet_matrix->x);
	int n = 0;
	for (typename Sparse_coordinate_matrix<T>::Coordinate_const_iterator it = 
		local_copy.coords_begin(); it != local_copy.coords_end(); ++it, ++n)
	{
		Ti[n] = it->i;	Tj[n] = it->j;	Tx[n] = it->value;
	}
	return triplet_matrix;
}

/// Convert sparse matrix from coordinate format to compressed-column form.
/// This conversion function is not implemented as a member function of 
/// Sparse_coordinate_matrix, as I want Sparse_coordinate_matrix more general
/// instead of depending on cholmod.
/// Mainly for internal use.
template <class T>
cholmod_sparse* create_cholmod_sparse(
	const Sparse_coordinate_matrix<T>& coor_matrix, cholmod_common* c)
{
	cholmod_triplet* triplet_matrix = create_cholmod_triplet(coor_matrix, c);
	if (triplet_matrix == 0) return 0; // empty matrix

	cholmod_sparse* output = cholmod_triplet_to_sparse(
		triplet_matrix, static_cast<int>(triplet_matrix->nnz), c);

	cholmod_free_triplet(&triplet_matrix, c);

	return output;
}

#endif // CHOLMOD_CONVERSION_H
// =====================================================================	