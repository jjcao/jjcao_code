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
#ifndef _GODOMAIN_H
#define _GODOMAIN_H


#include "Array.h"


namespace Go
{
///\addtogroup geometry
///\{

    /** Abstract base class representing a 2D parameter domain.
     *
     */

class Domain
{
 public:
    /// virtual destructor ensures safe inheritance
    virtual ~Domain() {}

    /// check whether a given parameter pair is located inside the domain
    /// \param point the (u,v)-pair that we want to test.
    /// \param tolerance the tolerance used (ruling what to do when 'point' is
    ///        located very near the edge of the domain).
    /// \return 'true' if 'point' is inside the domain, or within 'tolerance' 
    ///         from being inside the domain, 'false' otherwise'.
    virtual bool isInDomain(const Array<double, 2>& point, 
			    double tolerance) const = 0;

    /// check whether a given parameter pair is located on the domain boundary
    /// \param point the (u,v)-pair that we want to test
    /// \param tolerance the tolerance used (how 'far' from the boundary our
    ///        (u,v) pair can be and still be considered 'on' the boundary.
    /// \return 'true' if the point is considered to be on the boundary (within
    ///         'tolerance', 'false' otherwise.
    virtual bool isOnBoundary(const Array<double, 2>& point, 
			      double tolerance) const = 0;

    /// Find the (u, v) point in the Domain that is closest (using Euclidean distance
    /// in R^2) to a given (u, v) point.  If the given point is in the domain, then 
    /// the answer is obviously the same point.
    /// \param point the (u,v) parameter pair that we want to find the closest 
    ///              parameter pair to \em inside Domain.
    /// \param clo_pt the resulting closest parameter point.
    /// \param tolerance the tolerance used in defining whether the given point is 
    ///        already inside the domain.
    virtual void closestInDomain(const Array<double, 2>& point,
				 Array<double, 2>& clo_pt,
				 double tolerance) const = 0;
    
    /// Find the (u, v) point on the boundary of the Domain that is closest
    /// (using Euclidean distance in R^2) to a given (u, v) point.  If the
    /// point is already considered \em on the boundary, then the answer is obviously
    /// the same point.
    /// \param point the (u,v) parameter pair that we want to find to closest parameter
    ///              pair to \em on the Domain border.
    /// \param clo_bd_pt the resulting closest border point.
    /// \param tolerance the tolerance used in defining whether the given point is
    ////                 already inside the domain.
    virtual void closestOnBoundary(const Array<double, 2>& point,
				   Array<double, 2>& clo_bd_pt,
				   double tolerance) const = 0;

};


///\}
} // namespace Go


#endif // _GODOMAIN_H


