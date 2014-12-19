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
#ifndef _GOHERMITEINTERPOLATOR_H
#define _GOHERMITEINTERPOLATOR_H


#include "Interpolator.h"
#include "BsplineBasis.h"
#include "Point.h"

namespace Go
{
///\addtogroup geometry
///\{


/** An Interpolator that generates a hermite spline
 *  curve through given points and tangents.
 */

class HermiteInterpolator : public Interpolator
{
 public:
    /// Constructor takes no arguments
    HermiteInterpolator() {}
 
    /// Virtual destructor enables safe inheritance
    virtual ~HermiteInterpolator();

    // inherited from Interpolator
    virtual const BsplineBasis& basis();

    /// Hermite interpolation of a sequence of points with
    /// associated tangents and parameter values
    /// \param num_points number of points to interpolate
    /// \param dimension dimension of points to interpolate (2D, 3D, etc..)
    /// \param param_start pointer to the start of the array where
    ///                    the parameter values of the points are stored.
    ///                    This should be a strictly increasing sequence
    ///                    of 'num_points' values.
    /// \param data_start pointer to the start of the array where the
    ///                   points and tangents to be interpolated are stored.
    ///                   Each point and tangent consist of 'dimension'
    ///                   coordinates, and each tangent is stored immediately
    ///                   after its corresponding point.
    /// \retval coefs The control points of the computed hermite interpolation
    ///               curve will be returned in this vector. (Use the basis() 
    ///               function to get the associated b-spline basis).
    virtual void interpolate(int num_points,
			     int dimension,
			     const double* param_start,
			     const double* data_start,
			     std::vector<double>& coefs);

    /// Hermite interpolation of a sequence of points with 
    /// associated tangents and parameter values.
    /// \param data This vector contains the set of points and 
    ///             tangents to be interpolated.  The size of the
    ///             vector is twice the total number of points.
    ///             Each entry on the form [2i] represents a point,
    ///             and the entry [2i+1] represents the associated
    ///             tangent.
    /// \param param This vector represent the parameterization of the
    ///              points, and should have one entry per point.
    ///              (Making it half the size of the 'data' vector).
    ///              The parameter sequence must be strictly increaing.
    /// \retval coefs The control points of the computed hermite 
    ///               interpolation curve will be returned in this 
    ///               vector.  (Use the basis() function to get the 
    ///               associated b-spline basis.
    void interpolate(const std::vector<Point>& data,
		     const std::vector<double>& param,
		     std::vector<double>& coefs);

 private:
    BsplineBasis basis_;
}; 

///\}
} // namespace Go

#endif // _GOHERMITEINTERPOLATOR_H


 

