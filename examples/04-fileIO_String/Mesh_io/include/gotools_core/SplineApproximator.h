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
#ifndef _GOSPLINEAPPROXIMATOR_H
#define _GOSPLINEAPPROXIMATOR_H


#include "Interpolator.h"
#include "BsplineBasis.h"

namespace Go
{
///\addtogroup geometry
///\{



    /** An Interpolator that generates a spline curve approximating
     *  the given dataset in the least squares sense.
     */
class SplineApproximator : public Interpolator
{
public:
    /// Constructor takes no arguments
    SplineApproximator()
	: num_coefs_(0)
    {
    }

    /// Virtual destructor ensures safe inheritance
    virtual ~SplineApproximator();

    // inherited from Interpolator
    virtual const BsplineBasis& basis();

    /// The interpolating function, as inherited by \ref Interpolator.  
    /// Prior to calling this function, the user must have specified:
    /// - \em either the number of control points to use in the 
    ///   approximating curve by \ref setNumCoefs().
    /// - \em or/and directly specified the BsplineBasis to use
    ///   by \ref setSplineSpace().
    /// A default BsplineBasis will be generated if the user has only set
    /// the number of control points prior to calling interpolate().
    /// For parameter list, see \ref Interpolator.
    virtual void interpolate(int num_points,
			     int dimension,
			     const double* param_start,
			     const double* data_start,
			     std::vector<double>& coefs);
    
    /// Specify the number of basis functions / control points to use
    /// in the approximating curve.
    void setNumCoefs(int num) {
	num_coefs_ = num;
    }

    /// Directly specify the spline space in which to search for the 
    /// approximating function.
    void setSplineSpace(const BsplineBasis& basis) {
	basis_ = basis;
    }


private:
    int num_coefs_;
    BsplineBasis basis_;
};


///\}
} // namespace Go


#endif // _GOSPLINEAPPROXIMATOR_H


