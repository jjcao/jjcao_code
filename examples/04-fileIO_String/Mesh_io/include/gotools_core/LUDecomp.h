//===========================================================================
// GoTools - SINTEF Geometry Tools version 1.0.1
//
// GoTools module: CORE
//
// Copyright (C) 2000-2005 SINTEF ICT, Applied Mathematics, Norway.
//
// This program is free software; you can redistribute it and/or          
// modify it under the terms of the GNU General Public License            
// as published by the Free Software Foundation version 2 of the License. 
//
// This program is distributed in the hope that it will be useful,        
// but WITHOUT ANY WARRANTY; without even the implied warranty of         
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          
// GNU General Public License for more details.                           
//
// You should have received a copy of the GNU General Public License      
// along with this program; if not, write to the Free Software            
// Foundation, Inc.,                                                      
// 59 Temple Place - Suite 330,                                           
// Boston, MA  02111-1307, USA.                                           
//
// Contact information: e-mail: tor.dokken@sintef.no                      
// SINTEF ICT, Department of Applied Mathematics,                         
// P.O. Box 124 Blindern,                                                 
// 0314 Oslo, Norway.                                                     
//
// Other licenses are also available for this software, notably licenses
// for:
// - Building commercial software.                                        
// - Building software whose source code you wish to keep private.        
//===========================================================================
#ifndef _LUDECOMP_NEW_H
#define _LUDECOMP_NEW_H

#include <vector>
#include <stdexcept>

namespace Go 
{
///\addtogroup utils
///\{


/** LU decomposition algorithm, based on Crout's algorithm
 * \param mat The matrix to be decomposed.  The actually used class must support the 
 *            operation [][] and return 'double'.
 * \param num_rows Number of rows (which is also equal to the number of columns).
 * \param perm Should point to an n-sized array where the permutation is written.
 * \param parity Value upon function completion is 'true' if the number of row 
 *               interchanges was pair.  'false' otherwise.
 */
//===========================================================================
template<typename SquareMatrix> 
void LUDecomp(SquareMatrix& mat, int num_rows, int* perm, bool& parity);
//===========================================================================

/** Solve the system Ax = b for x, using LU decomposition of the matrix A.
 * Upon successful completion of the function, the matrix A will be LU decomposed,
 * and the solution x will be computed and stored at the memory area pointed to by the
 * last argument.
 * \param A   The system matrix A.  The actually used class must support the operation
 *            [][] and return 'double'.
 * \param num_unknowns Number of unknowns in the equation system.
 * \param vec At function invocation, 'vec' should point to an array of T, containing 
 *            the vector 'b'.  On successful completion of the function, this array 
 *            will contain the solution for x.
 */
//===========================================================================
template<typename SquareMatrix, typename T>
void LUsolveSystem(SquareMatrix& A, int num_unknowns, T* vec);
//===========================================================================

/** Using forward substitution to calculate x on the system Lx = b, where L is a
 * lower triangular matrix with unitary diagonal.
 * \param L The lower triangular matrix.  The actually used class must support the
 *          operation [][] and return 'double'.
 * \param x At function invocation, x should point to an array of T, containing the
 *          vector 'b'.  On successful completion of the function, this array will 
 *          contain the solution for x.
 * \param num_unknowns The system size (number of unknowns).
 */
//===========================================================================
template<typename SquareMatrix, typename T>
void forwardSubstitution(const SquareMatrix& L, T* x, int num_unknowns);
//===========================================================================

/** Using backward substitution to calculate x on the system Ux = b, where U is an
 * upper triangular matrix with unitary diagonal.
 * \param U The upper triangular matrix.  The actually used class must support the
 *          operation [][] and return 'double'.
 * \param x At function invocation, x should point to an array of T, containing the 
 *          vector 'b'.  On successful completion of the function, this array will 
 *          contain the solution for x.
 * \param num_unknowns The system size (number of unknowns).
 */
//===========================================================================
template<typename SquareMatrix, typename T>
void backwardSubstitution(const SquareMatrix& U, T* x, int num_unknowns);
//===========================================================================


///\}
}; // end namespace Go

#include "LUDecomp_implementation.h"

#endif // _LUDECOMP_NEW_H


