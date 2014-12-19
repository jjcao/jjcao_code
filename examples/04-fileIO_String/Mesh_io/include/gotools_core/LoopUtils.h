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
#ifndef _LOOPUTILS_H
#define _LOOPUTILS_H


#include "SplineCurve.h"
#include "CurveOnSurface.h"

#include <vector>
#include <boost/smart_ptr.hpp>

namespace Go{
///\addtogroup geometry
///\{


/// Functions for checking the orientation of loops (closed curves), and 
/// whether one loop on a surface encloses another.
///\}
namespace LoopUtils {
///\addtogroup geometry
///\{


    /// Check if a closed 2D-loop is oriented counterclockwise or not.
    /// \param simple_par_loop a sequence of 2D curves that are joined start-to-end and that form
    ///                        a closed loop in the plane.
    /// \param int_tol (geometric) tolerance used for internal computations (intersection detections)
    /// \return 'true' if the loop was found to be oriented CCW, otherwise 'false'.
    bool loopIsCCW(const std::vector<boost::shared_ptr<Go::SplineCurve> >& simple_par_loop, 
		   double int_tol);

    /// Check if a loop defined by CurveOnSurface s is oriented counterclockwise in the surface's 
    /// parametric domain.
    /// \param loop a sequence of CurveOnSurface s that are jointed start-to-end and that form
    ///             a closed loop on the surface.  
    /// \param int_tol (geometric) tolerance used for internal computations (intersection detections)
    /// \return 'true' if the loop was found to be oriented CCW, otherwise 'false'.
    bool paramIsCCW(const std::vector< boost::shared_ptr<Go::CurveOnSurface> >& loop, double int_tol);

    /// Loops expected to be disjoint, except possibly share part of boundary.

    /// Test whether one loop lies entirely within another.  This function does not work in the 
    /// general case; it makes the assumption that the loops do NOT intersect each other 
    /// transversally (their boundaries are allowed to tangentially touch though).  The algorithm
    /// works by testing a single point on the first loop for being inside the second loop, so 
    /// if the first loop lay partially inside, partially outside the second, the answer would be
    /// arbitrary.
    /// \param first_loop the first loop
    /// \param second_loop the second loop
    /// \param loop_tol the tolerance for defining coincidence between start/endpoints on the 
    ///                 consecutive curve segments consituting a loop.
    /// \param int_tol tolerance used for intersection calculations
    /// \return 'true' if 'first_loop' was found to be located inside 'second_loop' (given the 
    ///         assumptions above).  'false' otherwise.
    bool firstLoopInsideSecond(const std::vector<boost::shared_ptr<Go::CurveOnSurface> >& first_loop,
			       const std::vector<boost::shared_ptr<Go::CurveOnSurface> >& second_loop,
			       double loop_tol, double int_tol);


///\}
} // end namespace Go
///\addtogroup geometry
///\{

///\}
} // end namespace LoopUtils

#endif // _LOOPUTILS_H


