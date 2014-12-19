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
#ifndef _GOSURFACEGEN_H
#define _GOSURFACEGEN_H

#include "SplineSurface.h"
#include "SplineCurve.h"


/// This namespace contains functions used to create a Hahns Surface based on
/// input of a number of boundary curves (3<=nmb_cvs<=6).
/// Typically used for patching a model in which one or more surfaces are missing.

namespace HahnsSurfaceGen{
///\addtogroup geometry
///\{



    /// The routine which actually creates vector of blending surfaces, with given
    /// curves as their total bnd curve. # of return surfaces equals # of curves.
    /// All bnd curves share orientation (clockwise), and all cross curves point
    /// inwards. A missing cross curve is indicated by a 0 pointer.
    /// bnd_curves must form a loop. As we approximate modified cross curves,
    /// we include parameters to denote the exactness of the approximation.
    /// All curves expected to live in 3-dimensional space.
    /// \param bnd_curves edge curves for the Hahns Surface.
    /// \param cross_curves corresponding cross tangent curves.
    /// \param neighbour_tol allowed distance between corresponding end points.
    /// \param kink_tol allowed angle between corresponding tangents.
    /// \param knot_diff_tol parametric tolerance for equality of knots.
    /// \return vector containing the Hahns Surface.
    std::vector<boost::shared_ptr<Go::ParamSurface> >
    constructPolygonialSurface(std::vector<boost::shared_ptr<Go::ParamCurve> >& bnd_curves,
			       std::vector<boost::shared_ptr<Go::ParamCurve> >& cross_curves,
			       double neighbour_tol,
			       double kink_tol,
			       double knot_diff_tol);

    /// Given input of bnd_curves and cross_tangent curves, create the corresponding
    /// Hahns Surface.  Assuming input curves fulfill corner conditions.
    /// All curves expected to live in 3-dimensional space.
    /// \param bnd_curves edge curves for the Hahns Surface.
    /// \param mod_cross_curves corresponding cross tangent curves.
    /// \param neighbour_tol allowed distance between corresponding end points.
    /// \param kink_tol allowed angle between corresponding tangents.
    /// \param knot_diff_tol parametric tolerance for equality of knots.
    /// \return vector containing the Hahns Surface.
    std::vector<boost::shared_ptr<Go::ParamSurface> >
    constructHahnsSurface(std::vector<boost::shared_ptr<Go::SplineCurve> >& bnd_curves,
			  std::vector<boost::shared_ptr<Go::SplineCurve> >& mod_cross_curves,
			  double neighbour_tol,
			  double kink_tol,
			  double knot_diff_tol);

///\}
} // end of namespace



#endif // _GOSURFACEGEN_H

