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
#ifndef _GOINTERSECTIONS_H
#define _GOINTERSECTIONS_H

/// \file GoIntersections.h
/// Declaration file for a set of free intersection functions operating on
/// object belonging to GoTools but using the functionality of SISL.

#include "SplineCurve.h"
#include <utility>

namespace Go
{
///\addtogroup geometry
///\{


/// Enumeration of various pretopologies that can be associated with
/// intersections.
enum {  pretop_UNDEF,
	pretop_IN,  
	pretop_OUT, 
	pretop_ON,  
	pretop_AT };


/// Intersect two 2D spline curves. Collect intersection parameters
/// and pretopology information.
/// \param cv1 pointer to the first 2D curve 
/// \param cv2 pointer to the second 2D curve
/// \param epsge geometrical tolerance
/// \retval intersections this vector will contain the parameter pairs
///         of the found intersections (one vector entry per intersection.
///         The two parameters in the pair<> correspond to the parameter value
///         in 'cv1' and 'cv2' for a particular intersection.
/// \retval pretopology vector containing a pretopology indicator for each
///                     detected intersection point.  There is one entry per
///                     intersection point.
void intersect2Dcurves(const ParamCurve* cv1, 
		       const ParamCurve* cv2, 
		       double epsge,
		       std::vector<std::pair<double,double> >& intersections,
		       std::vector<int>& pretopology);

///Intersect two spline curves. Collect intersection parameters.
/// \param cv1 pointer to the first spline curve
/// \param cv2 pointer to the second spline curve
/// \param epsge geometrical tolerance
/// \retval intersections this vector will contain the parameter pairs
///         of the found intersections (one vector entry per intersection.
///         The two parameters in the pair<> correspond to the parameter value
///         in 'cv1' and 'cv2' for a particular intersection.
void intersectcurves(SplineCurve* cv1, SplineCurve* cv2, double epsge,
		     std::vector<std::pair<double,double> >& intersections);

/// Compute the closest point between two curves
/// \param cv1 pointer to the first curve
/// \param cv2 pointer to the second curve
/// \param epsge geometrical tolerance
/// \retval par1 parameter of the closest point in the first curve
/// \retval par2 parameter of the closest point in the second curve
/// \retval dist distance between the curves at the closest point
void closestPtCurves(SplineCurve* cv1, SplineCurve* cv2, double epsge,
		     double& par1, double& par2, double& dist);

///\}
} // namespace Go


#endif // _GOINTERSECTIONS_H

