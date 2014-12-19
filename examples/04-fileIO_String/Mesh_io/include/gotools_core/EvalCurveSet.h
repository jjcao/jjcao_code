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
#ifndef _GOEVALSURFACE_H_
#define _GOEVALSURFACE_H_


#include "Point.h"

#include <vector>

namespace Go
{
///\addtogroup geometry
///\{


/// This abstract class provides an interface to a set of curves that can
/// be evaluated.
/// Representing the actual geometry, typically used when iteratively
/// approximating the set of curves on the same basis.
class EvalCurveSet
{
public:
    /// virtual destructor ensures safe inheritance
    virtual ~EvalCurveSet()
    {}
  
    /// Evaluate the curves.
    /// \param t parameter in which to evaluate.
    /// \return the evaluated points for the curve set.
    virtual std::vector<Point> eval(double t)=0;

    /// Evaluate the curve derivatives.
    /// \param t parameter in which to evaluate.
    /// \param n number of derivatives to compute.
    /// \param der the evaluated points up to the n'th derivative for the curve set.
    virtual void eval(double t, int n, std::vector<std::vector<Point> >& der)=0; // n = order of diff

    /// Start parameter of domain.
    /// \return start parameter of the spline space.
    virtual double start()=0;

    /// End parameter of domain.
    /// \return end parameter of the spline space.
    virtual double end()=0;

    /// The geometric dimension of the spline curves.
    /// \return geometric dimension of the space.
    virtual int dim() = 0;

    /// Whether the approximation is within tolerances in input parameter.
    /// \param par parameter in which to evaluate.
    /// \param approxpos whether the input points are within tolerance from the
    ///                  evaluated points (as given by eval()).
    /// \param tol1 tolerance used to decide approximation accuracy.
    /// \param tol2 tolerance used to decide approximation accuracy.
    /// \return whether the approximation is within tolerances in input parameter.
    virtual bool approximationOK(double par, const std::vector<Point>& approxpos,
				 double tol1, double tol2)=0;

    /// The number of curves in the curve set.
    /// \return the number of curves in the curve set.
    virtual int nmbCvs() = 0;
};

///\}
} // namespace Go

#endif // _GOEVALSURFACE_H_

