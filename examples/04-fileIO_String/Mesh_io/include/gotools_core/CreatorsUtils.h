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
#ifndef _CREATORSUTILS_H
#define _CREATORSUTILS_H


#include "SplineCurve.h"
#include "SplineSurface.h"
#include "CurveOnSurface.h"
#include "BsplineBasis.h"
#include "boost/smart_ptr.hpp"

namespace Go {
///\addtogroup geometry
///\{


/// Related to the generation of cross tangent curves.
///\}
namespace CreatorsUtils
{
///\addtogroup geometry
///\{


    /// Helper function. Returns the parametric representation of the input
    /// curve.
    /// The returned curve is created (i.e. it should be be handled by a
    /// smart ptr)!
    /// This function takes as argument a vector of (shared pointers to)
    /// CurveOnSurfaces, which
    /// are assumed to lie on the same surface and to represent consecutive
    /// parts of a larger
    /// curve, and return the parametric representation of this larger curve.
    /// \param cv a vector of (shared pointers to) CurveOnSurfaces.
    ///           We suppose that these curves
    ///           are all lying on the same surface and that they constitute
    ///           consecutive parts
    ///           of a global curve.  
    /// \return if successful, a pointer to the newly generated SplineCurve,
    ///         which represents
    ///         the global curve that is the union of the input curves.
    ///         This SplineCurve is 
    ///         represented in the parameter plane of the surface, and is
    ///         therefore a 2D curve.
    ///         If it could not be constructed (due a failure of any of
    ///         the provided CurveOnSurfaces
    ///         to provide their internal parameter curve), a null pointer
    ///         is returned instead.
    ///         The user assumes ownership of the SplineCurve.
    SplineCurve* getParametricCurve
    (const std::vector<boost::shared_ptr<const CurveOnSurface> >& cv);

    /// Generate the inwards pointing cross-tangent curve along the input trim
    /// curve.
    /// \param cv the trim curve.
    /// \return the inwards cross tangent curve.
    boost::shared_ptr<Go::SplineCurve>
    createCrossTangent(const Go::CurveOnSurface& cv);

    /// The cross tangent cv along input cv is created. Length given by length
    /// of derivs in
    /// sf along cv or by input_cross_cv if it exists. Direction of created
    /// cross tangent cv
    /// is to the left of cv (i.e. it should point into the surface).
    /// If input_cross_cv exists so does basis_space_cv and they must share
    /// parametrization.
    /// If input_cross_cv exists, but not basis_space_cv, it shares
    /// parametrization with cv.
    /// \param cv the trim curve along which we are to compute the cross
    ///        tangent cv.
    /// \param basis_space_cv if != NULL we use the basis and parametrization
    ///                       of basis_space_cv.
    /// \param cross_cv_ref if != NULL the curve defines the angle between
    ///                     the boundary curve and the new cross tangent
    ///                     curve.
    ///                     Additionally it defines the length of the new
    ///                     cross tangent curve.
    /// \param appr_offset_cv whether the method should approximate the offset
    ///                       curve
    ///                       (as opposed to the cross tangent curve).
    /// \return the cross tangent (or offset) curve.
    boost::shared_ptr<Go::SplineCurve>
    createCrossTangent(const Go::CurveOnSurface& cv,
		       boost::shared_ptr<Go::SplineCurve> basis_space_cv,
		       const Go::SplineCurve* cross_cv_ref,
		       bool appr_offset_cv = true);

///\}
} // of namespace CreatorsUtils.
///\addtogroup geometry
///\{


///\}
}; // end namespace Go



#endif // _CREATORSUTILS_H

