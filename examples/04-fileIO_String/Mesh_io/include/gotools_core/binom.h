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
#ifndef _BINOM_H
#define _BINOM_H

#include <stdexcept>
#include <vector>

namespace Go {
///\addtogroup utils
///\{



/// Computes the binomial coefficient: n! / (i! (n-i)!)
inline double binom(int n, int i)
{
    if (i < 0 || i > n)
	return 0.0;

    static std::vector<std::vector<double> > pascals_triangle(10);
    static bool first_time = true;
    if (first_time) {
	first_time = false;
	pascals_triangle.resize(1);
	pascals_triangle[0].resize(1);
	pascals_triangle[0][0] = 1;
    }
    int old_size = pascals_triangle.size();
    if (old_size < n+1) {
	// We must expand the triangle
	pascals_triangle.resize(n+1);
	// Compute the terms of the new rows
	for (int nr = old_size; nr < n+1; ++nr) {
	    pascals_triangle[nr].resize(nr+1);
	    pascals_triangle[nr][0] = 1.0;
	    for (int j = 1; j < nr; ++j) {
		pascals_triangle[nr][j]
		    = pascals_triangle[nr-1][j-1] + pascals_triangle[nr-1][j];
	    }
	    pascals_triangle[nr][nr] = 1.0;
	}
    }
    return pascals_triangle[n][i];
}

/// computes n! (n factorial)
inline double factorial(int n)
{
    double res = 1;
    for (int i = 2; i <= n; ++i) {
 	res *= i;
    }
    return res;
}

/// computes the trinomial coefficient: n! / (i! j! (n-i-j)!)
inline double trinomial(int n, int i, int j)
{
    if (i < 0 || i > n || j < 0 || j > n)
	return 0;

    return binom(n, i) * binom(n-i, j);
}


/// computes the quadrinomial coefficient: n! / (i! j! k! (n-i-j-k)!)
inline double quadrinomial(int n, int i, int j, int k)
{
    if (i < 0 || i > n || j < 0 || j > n || k < 0 || k > n)
	return 0;

    return binom(n, i) * binom(n-i, j) * binom(n-i-j, k);
}

///\}
}; // end Go

#endif // _BINOM_H


