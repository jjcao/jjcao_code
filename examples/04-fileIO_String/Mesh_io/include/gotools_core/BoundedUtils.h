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
#ifndef _GOBOUNDEDUTILS_H
#define _GOBOUNDEDUTILS_H


#include "BoundedSurface.h"
#include "SplineSurface.h"
#include "CurveOnSurface.h"
#include "LoopUtils.h"

#include <boost/smart_ptr.hpp>

using std::vector;
using boost::shared_ptr;
using Go::SplineSurface;
using Go::BoundedSurface;
using Go::CurveOnSurface;
using Go::Point;
using Go::ParamCurve;

namespace Go {
///\addtogroup geometry
///\{


/// Functions related to the trimming of surfaces, etc.  Also contains functions
/// for spatial transformations of such surfaces (rotation, translation, etc.)
///\}
namespace BoundedUtils {
///\addtogroup geometry
///\{


    /// Extract those parts of a given CurveOnSurface where its parameter curve
    /// lies inside the parameter domain of a BoundedSurface.
    /// \param curve the CurveOnSurface that we want to check
    /// \param bounded_surf the BoundedSurface whose parameter domain we will
    ///                     check against
    /// \param epsge geometric tolerance used in calculations
    /// \return a vector containing those segments of 'curve' that have parameter
    ///         descriptions inside the parameter domain of the 'bounded_surf'.
    std::vector<shared_ptr<CurveOnSurface> >
      intersectWithSurface(CurveOnSurface& curve,
			   BoundedSurface& bounded_surf, double epsge);

    /// We have two set of CurveOnSurface s, 'curves1' and 'curves2', and two
    /// BoundedSurface s, 'bd_sf1' and 'bd_sf2'.  We extract those segments of
    /// curves in 'curves1' and 'curves2' that have parameter curves in the 
    /// parameter domains of respective 'bd_sf1' and 'bd_sf2'.  Then we compare
    /// the resulting segments from the two sets against each others, and keep 
    /// those that overlap spatially.  Segments are split in two if necessary
    /// in order to make start and end points from one set coincide with those from
    /// the other set.  'curves1' and 'curves2' are then cleared and filled with the
    /// resulting curve segments.
    /// \param curves1 see above, we suppose that the curves are NOT self-intersecting
    /// \param bd_sf1 see above, we suppose that the underlying surface is the same
    ///               as the one refered to by the curves in 'curves1'.
    /// \param curves2 see above, we suppose that the curves are NOT self-intersecting
    /// \param bd_sf2 see above, we suppose that the underlying surface is the same
    ///               as the one refered to by the curves in 'curves2'.
    /// \param epsge geometric epsilon used in closest point calculations, checks for
    ///              coincidence, etc.
    void intersectWithSurfaces(vector<shared_ptr<CurveOnSurface> >& curves1,
			       shared_ptr<BoundedSurface>& bd_sf1,
			       vector<shared_ptr<CurveOnSurface> >& curves2,
			       shared_ptr<BoundedSurface>& bd_sf2,
			       double epsge);

    /// Given input of partial boundary curves, extract parts of boundary making it a
    /// boundary loop (or more). Input segments expected to be ordered, going in the same direction.
    /// part_bd_cvs should lie on sf (as oppsed to only parts of cvs).

    /// This function tries to complete "partial" boundary loops by filling out
    /// the missing parts using fragments from the domain boundary of a BoundedSurface.
    /// \param sf the surface whose domain boundary will be used
    /// \param part_bnd_cvs a vector of (shared pointers to) curve segments that represent
    ///                     incomplete loops.  Upon function return, this vector will be emptied.
    /// \return a vector contained the loops that the function was able to completely
    ///         close using curve segments from 'part_bnd_cvs' and the domain boundaries of 'sf'.
    vector<vector<shared_ptr<CurveOnSurface> > >
      getBoundaryLoops(const BoundedSurface& sf, 
		       vector<shared_ptr<CurveOnSurface> >& part_bnd_cvs);

    /// We intersect a parametric surface with a plane, and return the surface(s)
    /// consisting only of the part(s) of the surface that were located on the 
    /// positive side of the intersection.  If there was no intersection, an empty
    /// stl-vector is returned.  The plane is defined by its normal and a point 
    /// located on it.
    /// \param surf the parametric surface.  It must be either a BoundedSurface or 
    ///             a SplineSurface.
    /// \param point a point on the plane to intersect against
    /// \param normal normal of the plane to intersect against
    /// \param epsge geometric tolerance
    /// \return the surface(s) consisting of the part(s) of 'surf' that were located
    ///         on the positive side of the intersection.
    vector<shared_ptr<BoundedSurface> >
      trimWithPlane(const shared_ptr<Go::ParamSurface>& surf,
		    Point point, Point normal, double epsge);

    /// must be BoundedSurfaces or SplineSurface
    /// underlying surfaces must be of type SplineSurface
    
    /// If the argument surfaces intersect, and if the intersection curves result in
    /// new boundary loops being defined, then the new surface parts defined within these
    /// domains will be returned.  Otherwise, the return vector will be empty.
    /// \param sf1 the first surface to participate in the intersection
    /// \param sf2 the second surface to participate in the intersection
    /// \param epsge geometrical tolerance to be used in computations
    /// \return a vector with the surfaces representing the parts of the original surfaces
    ///         enclosed by new parametrical loops arising when combining existing loops
    ///         with the curves defined by the intersection.
    vector<shared_ptr<BoundedSurface> >
    trimSurfWithSurf(const shared_ptr<Go::ParamSurface>& sf1,
		     const shared_ptr<Go::ParamSurface>& sf2, double epsge);
    

    /// If surf already is a BoundedSurface, return clone. If SplineSurface,
    /// convert. Otherwise, Error. Return surface is created inside function.

    /// Convert a SplineSurface to a BoundedSurface.  All information is copied, nothing
    /// is shared. 
    /// \param surf the SplineSurface to convert
    /// \param space_epsilon the tolerance assigned to the newly created BoundedSurface
    /// \return (pointer to) a BoundedSurface that represent the same surface as 'surf'.  
    ///         The user assumes ownership of the object.
    BoundedSurface* convertToBoundedSurface(const SplineSurface& surf,
					      double space_epsilon);


   /// All input loops are expected to be simple, lying on surface. They are sorted based
   /// on orientation. No pair of curves with the same orientation may lie inside/outside eachother.
   /// It they do the outer/inner (ccw/cw) loop(s) will be erased.
   /// Furthermore assuming no pair of loops intersect (may touch tangentially).


    /// Define the surfaces that result from trimming a given SplineSurface with a set of 
    /// boundary loops.  Counterclockwise loops define the interior of a parameter domain, while
    /// clockwise loops define holes in the domain.
    /// \param loops each entry in the outermost vector represent a vector of curves that together
    ///              specify a loop in 2D parametrical space of the surface 'under_sf'.  These
    ///              are the trim curves.  The curves are expected to be simple, lying on the 
    ///              surface 'under_sf'.  No pair of curve loops with the same orientation may lie
    ///              inside/outside of each other; in that case, the irrelevant loop will be reased.
    ///              Furthermore, we assume that no pair of loops intersect transversally (they are
    ///              still allowed to touch tangentially).
    /// \param under_sf the underlying SplineSurface that we are going to trim with the curves
    ///                 given in 'loops'.
    /// \param epsgeo geometrical tolerance used in computations
    /// \return a vector containing BoundedSurface s that each represent a trimmed part of the
    ///         'under_sf' surface.
    vector<shared_ptr<BoundedSurface> >
     createTrimmedSurfs(vector<vector<boost::shared_ptr<CurveOnSurface> > >& loops,
			shared_ptr<SplineSurface>& under_sf, 
			double epsgeo);

    /// Find the intersection curve(s) between a SplineSurface and a given plane.
    /// The plane is defined by its normal and a point on the plane.
    /// \param surf the SplineSurface to intersect with the plane
    /// \param pnt a point lying on the plane
    /// \param normal the normal of the plane
    /// \param geom_tol geometrical tolerance to be used for intersection computations
    /// \return a vector with (shared pointers to) CurveOnSurface s, which represent
    ///         the intersection curves found.
    std::vector<shared_ptr<CurveOnSurface> >
      intersectWithPlane(shared_ptr<SplineSurface>& surf,
			 Point pnt, Point normal, double geom_tol);

    /// Find the intersction curve(s) between two spline surfaces
    /// \param sf1 the first surface to intersect
    /// \param sf2 the second surface to intersect
    /// \retval int_segments1 a vector of CurveOnSurface s, representing the intersection 
    ///                       curves as lying on 'sf1'.
    /// \retval int_segments2 a vector of CurveOnSurface s, representing the intersection
    ///                       curves as lying on 'sf2'.
    /// \param epsge geometrical tolerance to be used for intersection computations
    void getIntersectionCurve(shared_ptr<SplineSurface>& sf1,
			      shared_ptr<SplineSurface>& sf2,
			      vector<shared_ptr<CurveOnSurface> >& int_segments1,
			      vector<shared_ptr<CurveOnSurface> >& int_segments2,
			      double epsge);


    /// Translate a given BoundedSurface 
    /// \param trans_vec the vector specifying the translation to apply to the surface
    /// \param bd_sf the surface to translate
    /// \param deg_eps an epsilon value used when determining degenerate boundary loops
    void translateBoundedSurf(Point trans_vec, BoundedSurface& bd_sf,
			      double deg_eps);

    // Rotate a given BoundedSurface
    /// \param rot_axis a vector specifying the axis of rotation
    /// \param alpha the angle of rotation (in radians)
    /// \param bf_sf the surface to rotate
    /// \param deg_eps an epsilon value used when determining degenerate boundary loops
    void rotateBoundedSurf(Point rot_axis, double alpha,
			   BoundedSurface& bf_sf, double deg_eps);

///\}
} // namespace Go
///\}
} // namespace BoundedUtils

#endif // _GOBOUNDEDUTILS_H

