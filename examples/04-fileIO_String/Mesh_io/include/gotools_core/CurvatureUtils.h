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
#ifndef _CURVATUREUTILS_
#define _CURVATUREUTILS_

#include "Point.h"
#include <vector>

namespace Go
{
///\addtogroup utils
///\{

    /** Help functions in related to curvature.
     */

    /// Given position, first and second derivative
    /// of a curve passing through a point, compute
    /// the unit tangent, curvature vector and curvature 
    /// radius of this curve.
    double curvatureRadius(const std::vector<Point>& der,
			   std::vector<Point>& unitder);

    /// Computes the step length along a curve based on radius of curvature
    /// at a point on the curve, and an absolute tolerance.
    double stepLenFromRadius(double radius, double aepsge);

    /// To create the tangent length for interpolating a
    /// circular arc with an almost equi-oscillating Hermit qubic.
    double tanLenFromRadius(double radius, double angle);

    /// Given position, first and second derivative in both ends of
    /// an Hermite segment, compute parameter interval and tangent lengths
    /// in order to stay close to a circular segment.
    void getHermiteData(const std::vector<Point>& der1,
			const std::vector<Point>& der2, 
			double& parint, double& len1, double& len2);

///\}
} // End of namespace Go


#endif


