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
#ifndef _GOSURFACECREATORS_H
#define _GOSURFACECREATORS_H


#include "SplineSurface.h"
#include "SplineCurve.h"
#include "CurveOnSurface.h"


namespace Go
{
///\addtogroup geometry
///\{


/// Various functions for generating SplineSurface s by approximation, blending,
/// etc.
///\}
namespace SurfaceCreators
{
///\addtogroup geometry
///\{


    /// Given input of two parametric surfaces, which intersect along an edge,
    /// create a smooth surface by moving boundary of surfs to inner part of
    /// current surfs. Normals of input surfs are assumed to be consistent.
    /// \param surfs the input surfaces, size of vector is 2.
    /// \param int_cvs intersection curves between the two surfaces.  Vector has size 2.
    ///                Both parameter curves should exist, the space curve should be the same.
    /// \param dist_0 the offset space distance in surfs[0].
    /// \param dist_1 the offset space distance in surfs[1].
    /// \param epsge the geometrical tolerance when offsetting.
    /// \param trim_crvs the offset trim curves: par_cv0, space_cv0, par_cv1, space_cv1.
    boost::shared_ptr<SplineSurface>
    createSmoothTransition(const std::vector<boost::shared_ptr<const ParamSurface> >& surfs,
			   const std::vector<boost::shared_ptr<const CurveOnSurface> >& int_cvs,
			   double dist_0, double dist_1, double epsge,
			   std::vector<boost::shared_ptr<SplineCurve> >& trim_crvs);

    /// Return the product of the two spline surfaces. Expecting them to be
    /// non-rational.
    /// \param sf1 the first surface.
    /// \param sf2 the second surface.
    /// \return the surface product.
    boost::shared_ptr<SplineSurface> mult1DSurfaces(const SplineSurface& sf1,
						    const SplineSurface& sf2);

    /// Return the product of the two Bezier patches. Expecting them to be
    /// non-rational.
    /// \param patch1 SplineSurface of Bezier type (i.e. no inner knots).
    /// \param patch2 SplineSurface of Bezier type (i.e. no inner knots).
    /// \return the surface product.
    boost::shared_ptr<SplineSurface> mult1DBezierPatches(const SplineSurface& patch1,
							 const SplineSurface& patch2);

    /// Return the rational surface 'nom_sf/den_sf'.
    /// \param nom_sf the nominator surface.
    /// \param den_sf the denominator surface.
    /// \param weights_in_first true if the coefficients of the nom_sf have been multiplied
    ///                         by the corresonding rational coefficients in den_sf.
    /// \return the rational surface.
    boost::shared_ptr<Go::SplineSurface> mergeRationalParts(const Go::SplineSurface& nom_sf,
							    const Go::SplineSurface& den_sf,
							    bool weights_in_first = false);

    /// Given input of 1d-sf, return the 3d visualization (u, v, f(u, v)).
    /// \param sf_1d 1-dimensional surface.
    /// \return the 3-dimensional surface.
    boost::shared_ptr<Go::SplineSurface> insertParamDomain(const Go::SplineSurface& sf_1d);

    /// Assuming the sf is rational, separate the geometric space from the homogenuous.
    /// \param sf the input spline surface.
    /// \return the non-rational parts of sf.
    std::vector<boost::shared_ptr<Go::SplineSurface> >
    separateRationalParts(const Go::SplineSurface& sf);

///\}
} // end of namespace SurfaceCreators
///\addtogroup geometry
///\{


///\}
} // end of namespace Go


#endif // _GOSURFACECREATORS_H


