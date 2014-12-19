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
#ifndef _FACTORIAL_H
#define _FACTORIAL_H


namespace Go {
///\addtogroup utils
///\{



/// Compile-time factorial calculations.


/// Calculate the factorial of a given integer N.
template<int N>
struct Factorial {
    enum { value = N * Factorial<N-1>::value };
};


/// \internal
template <>
struct Factorial<1> {
    enum { value = 1 };
};


/// Compute the inverse of the factorial of a given integer N, where
/// T is typically 'float' or 'double'.
template<typename T, int N>
struct InverseFactorial {
    static inline T val()
    {
	return T(1) / T(Factorial<N>::value);
    }
};
 

///\}
};
#endif // _FACTORIAL_H


