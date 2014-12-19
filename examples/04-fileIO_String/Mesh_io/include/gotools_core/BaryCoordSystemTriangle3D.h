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
#ifndef _BARYCOORDSYSTEMTRIANGLE3D_H
#define _BARYCOORDSYSTEMTRIANGLE3D_H


#include "Array.h"
#include "Volumes.h"
#include <iostream>
#include <algorithm>


namespace Go {
///\addtogroup utils
///\{


    /** A barycentric coordinate system for a triangle (2-manifold) embedded in 3D.
     *  Note that this differs from what can be expressed using the 
     *  BaryCoordSystem template in that for the latter, the dimension of the simplex
     *  and the dimension of the space is equal.
     */

class BaryCoordSystemTriangle3D {
public:
    /// Empty default constructor
    BaryCoordSystemTriangle3D() { }

    /// Constructor. Takes an array of points that will become
    /// the corners of the coordinate simplex.
    BaryCoordSystemTriangle3D(const Array<double, 3>* corners)
    {
	std::copy(corners, corners+3, corners_);
	total_area_ = area(corners);
	normal_ = (corners[1]-corners[0]) % (corners[2]-corners[0]);
	normal_.normalize();
    }

    /// Input is a barycentric point, output is the corresponding
    /// cartesian point.
    template <typename T>
    Array<T, 3> baryToCart(const Array<T, 3>& bary_pt) const
    {
        Array<T, 3> cart_pt(T(0.0), T(0.0), T(0.0));
	for (int i = 0; i < 3; ++i) {
	    cart_pt[0] += corners_[i][0] * bary_pt[i];
	    cart_pt[1] += corners_[i][1] * bary_pt[i];
	    cart_pt[2] += corners_[i][2] * bary_pt[i];
	}
	return cart_pt;
    }

    /// Input is a cartesian point, output is the corresponding
    /// barycentric point.
    template <typename T>
    Array<T, 3> cartToBary(const Array<T, 3>& cart_pt) const
    {
        static Array<T, 3> subtriangle[3];
	int i;
	for (i = 1; i < 3; ++i) {
	    subtriangle[i] = corners_[i];
	}

	Array<T, 3> bary_pt;
	for (i = 0; i < 3; ++i) {
	    subtriangle[i] = cart_pt;
	    bary_pt[i] = signed_area(subtriangle, normal_);
	    subtriangle[i] = corners_[i];
	}
	return bary_pt / total_area_;
    }

private:
    Array<double, 3> corners_[3];
    Array<double, 3> normal_;
    double total_area_;
};


///\}
} // namespace Go

#endif // _BARYCOORDSYSTEMTRIANGLE3D_H


