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
#ifndef _GOPARAMSURFACE_H
#define _GOPARAMSURFACE_H

#include "Domain.h"
#include "GeomObject.h"
#include "RectDomain.h"
#include "CurveLoop.h"
#include <vector>


namespace Go
{
///\addtogroup geometry
///\{


class CurveLoop; // forward declaration

  /** Base class for parametric surfaces in Go
   *
   */
class ParamSurface : public GeomObject
{
public:
    /// Virtual destructor, enables safe inheritance.
    virtual ~ParamSurface();

    /// make a clone of this surface and return a pointer to it (user is 
    /// responsible for clearing up memory afterwards).
    /// \return pointer to cloned object
    virtual ParamSurface* clone() const = 0;
    /// Return the parameter domain of the surface.  This may be a simple
    /// rectangular domain (\ref RectDomain) or any other subclass of
    /// \ref Domain (such as GoCurveBoundedDomain, found in the
    /// \c sisl_dependent module).
    /// \return a Domain object describing the parametric domain of the surface
    virtual const Domain& parameterDomain() const = 0;

    /// Get a rectangular parameter domain that is guaranteed to contain the
    /// surface's \ref parameterDomain().  It may be the same.  There is no
    /// guarantee that this is the smallest domain containing the actual domain.
    /// \return a RectDomain that is guaranteed to include the surface's total
    ///         parameter domain.
    virtual RectDomain containingDomain() const = 0;


    /// Returns the anticlockwise, outer boundary loop of the surface.
    /// \param degenerate_epsilon edges whose length is smaller than this value
    ///        are ignored.
    /// \return a CurveLoop describing the anticlockwise, outer boundary loop of 
    ///         the surface.
    virtual CurveLoop outerBoundaryLoop(double degenerate_epsilon
					  = DEFAULT_SPACE_EPSILON) const = 0;

    /// Returns the anticlockwise outer boundary loop of the surface, together with 
    /// clockwise loops of any interior boundaries, such that the surface always is
    /// 'to the left of' the loops.
    /// \param degenerate_epsilon edges whose length is smaller than this value are
    ///                           ignored.
    /// \return a vector containing CurveLoops.  The first of these describe the 
    ///         outer boundary of the surface (clockwise), whereas the others describe
    ///         boundaries of interior holes (clockwise).
    virtual std::vector<CurveLoop> allBoundaryLoops(double degenerate_epsilon
						      = DEFAULT_SPACE_EPSILON) const = 0;

    /// Creates a DirectionCone covering all normals to this surface.
    /// \return a DirectionCone (not necessarily the smallest) containing all normals 
    ///         to this surface.
    virtual DirectionCone normalCone() const = 0;
    
    /// Creates a DirectionCone covering all tangents to 
    /// this surface along a given parameter direction.
    /// \param pardir_is_u if 'true', then the DirectionCone will be defined on basis 
    ///        of the surface's tangents along the first parameter direction.  Otherwise
    ///        the second parameter direction will be used.
    /// \return a DirectionCone (not necessarily the smallest) containing all tangents
    ///         to this surface along the specified parameter direction.
    virtual DirectionCone tangentCone(bool pardir_is_u) const = 0;

    /// Creates a composite box enclosing the surface. The composite box
    /// consists of an inner and an edge box. The inner box is
    /// supposed to be made from the interior of the surface, while the
    /// edge box is made from the boundary curves. The default
    /// implementation simply makes both boxes identical to the
    /// regular bounding box.
    /// \return the CompositeBox of the surface, as specified above
    virtual CompositeBox compositeBox() const;
   
    /// Evaluates the surface's position for a given parameter pair.
    /// \param pt the result of the evaluation is written here 
    /// \param upar the first parameter
    /// \param vpar the second parameter
    virtual void point(Point& pt, double upar, double vpar) const = 0;

    /// Evaluates the surface's position and a certain number of derivatives
    /// for a given parameter pair.
    /// \param pts the vector containing the evaluated values.  Its size must be 
    ///            set by the user prior to calling this function, and should be
    ///            equal to (derivs+1) * (derivs+2) / 2.  Upon completion of the
    ///            function, its first entry is the surface's position at the given
    ///            parameter pair.  Then, if 'derivs' > 0, the two next entries will
    ///            be the surface tangents along the first and second parameter 
    ///            direction.  The next three entries are the second- and cross 
    ///            derivatives, in the order (du2, dudv, dv2), and similar for 
    ///            even higher derivatives.
    /// \param upar the first parameter 
    /// \param vpar the second parameter
    /// \param derivs number of requested derivatives
    /// \param u_from_right specify whether derivatives along the first parameter are
    ///                     to be calculated from the right ('true', default) or from
    ///                     the left ('false')
    /// \param v_from_right specify whether derivatives along the second parameter are
    ///                     to be calculated from the right ('true', default) or from
    ///                     the left ('false')
    /// \param resolution tolerance used when determining whether parameters are located 
    ///                   at special values of the parameter domain (in particualar; knot
    ///                   values in case of spline objects.
    virtual void point(std::vector<Point>& pts, 
		       double upar, double vpar,
		       int derivs,
		       bool u_from_right = true,
		       bool v_from_right = true,
		       double resolution = 1.0e-12) const = 0;

    /// Evaluate the surface's position at a certain parameter pair
    /// \param upar the first parameter
    /// \param vpar the second parameter
    /// \return the surface's position for this parameter pair.
    /// NB: This function is implemented in terms of the ParamSurface's virtual
    /// 'point(...)' function, but is itself not virtual.  If you make a concrete 
    /// subclass and wish to make this function visible to the user, you must
    /// put a 'using ParamSurface::point' statement in the class definition.
    Point point(double upar, double vpar) const;

    /// Evaluate the surface's position and a certain number of derivatives at
    /// a given parameter.
    /// \param upar the first parameter 
    /// \param vpar the second parameter
    /// \param derivs number of requested derivatives
    /// \return the vector containing the evaluated values.  Its size will be
    ///         equal to (derivs+1) * (derivs+2) / 2.  Its first entry is the surface's 
    ///         position at the given parameter pair.  Then, if 'derivs' > 0, the two next
    ////        entries will be the surface tangents along the first and second parameter 
    ///         direction.  The next three entries are the second- and cross derivatives, 
    ///         in the order (du2, dudv, dv2), and similar for even higher derivatives.
    /// NB: This function is implemented in terms of the ParamSurface's virtual
    /// 'point(...)' functions, but is itself not virtual.  If you make a concrete 
    /// subclass and wish to make this function visible to the user, you must
    /// put a 'using ParamCurve::point' in the class definition.
    std::vector<Point> point(double upar, double vpar,
			       int derivs) const;

    /// Evaluates the surface normal for a given parameter pair
    /// \param n the computed normal will be written to this variable
    /// \param upar the first parameter
    /// \param vpar the second parameter
    virtual void normal(Point& n, double upar, double vpar) const = 0;

    /// Get the curve(s) obtained by intersecting the surface with one of its constant
    /// parameter curves.  For surfaces without holes, this will be the parameter curve
    /// itself; for surfaces with interior holes this may be a collection of several, 
    /// disjoint curves.  
    /// \param parameter parameter value for the constant parameter (either u or v)
    /// \param pardir_is_u specify whether the \em moving parameter (as opposed to the 
    ///                    \em constant parameter) is the first ('true') or the second
    ///                    ('false') one.
    /// \return a vector containing shared pointers to the obtained, newly constructed
    ///          constant-parameter curves.
    virtual std::vector<boost::shared_ptr<ParamCurve> >
    constParamCurves(double parameter, bool pardir_is_u) const = 0;

    /// Get the surface(s) obtained by cropping the parameter domain of this surface
    /// between given values for the first and second parameter.  In general, for 
    /// surfaces with no interior holes, the result will be \em one surface; however,
    /// for surfaces with interior holes, the result might be \em several \em disjoint
    /// surfaces.
    /// \param from_upar lower value for the first parameter in the subdomain
    /// \param from_vpar lower value for the second parameter in the subdomain
    /// \param to_upar upper value for the first parameter in the subdomain
    /// \param to_vpar upper value for the second parameter in the subdomain
    /// \param fuzzy tolerance used when determining intersection with interior 
    ///        boundaries
    /// \return a vector contained shared pointers to the obtained, newly constructed
    ///         sub-surfaces.
    virtual std::vector<boost::shared_ptr<ParamSurface> >
    subSurfaces(double from_upar, double from_vpar,
		double to_upar, double to_vpar,
		double fuzzy = DEFAULT_PARAMETER_EPSILON) const = 0;

    /// Determine the parameter value of the start of the 'next
    /// segment' from a parameter value, along a given parameter
    /// direction.  A 'segment' is here defined as a parameter
    /// interval in which there will be no discontinuities in
    /// derivatives or other artifacts.  For spline objects, a segment
    /// will typically be the interval between two consecutive,
    /// non-coincident knots.
    /// \param dir the parameter direction in which we search for the
    /// next segment (0 or 1)
    /// \param par the parameter value starting from which we search
    /// for the start value of the next segment
    /// \param forward define whether we shall move forward ('true')
    /// or backwards when searching along this parameter
    /// \param tol tolerance used for determining whether the 'par' is
    /// already located \em on the next segment value
    /// \return the value of the start value of the next segment (or
    /// the end of the previous segment, if we are moving
    /// backwards...)
    virtual double nextSegmentVal(int dir, double par, bool forward, double tol) const = 0;

    /// Iterates to the closest point to pt on the surface. 
    /// \param pt the point to find the closest point to
    /// \param clo_u u parameter of the closest point
    /// \param clo_v v parameter of the closest point
    /// \param clo_pt the geometric position of the closest point
    /// \param clo_dist the distance between pt and clo_pt
    /// \param epsilon parameter tolerance (will in any case not be higher than
    ///                sqrt(machine_precision) x magnitude of solution
    /// \param domain_of_interest pointer to parameter domain in which to search for 
    ///                           closest point. If a NULL pointer is used, the entire
    ///                           surface is searched.
    /// \param seed pointer to parameter values where iteration starts.
    virtual void closestPoint(const Point& pt,
			      double&        clo_u,
			      double&        clo_v, 
			      Point&       clo_pt,
			      double&        clo_dist,
			      double         epsilon,
			      const RectDomain* domain_of_interest = NULL,
			      double   *seed = 0) const = 0;

    void singularity(double& sing_u,
		     double& sing_v, 
		     Point& sing_pt,
		     double& sing_dist,
		     double epsilon,
		     const RectDomain* rd = NULL,
		     double *seed = 0) const;

    /// Iterates to the closest point to pt on the boundary of the surface.
    /// \see closestPoint()
    virtual void closestBoundaryPoint(const Point& pt,
				      double&        clo_u,
				      double&        clo_v, 
				      Point&       clo_pt,
				      double&        clo_dist,
				      double epsilon,
				      const RectDomain* rd = NULL,
				      double *seed = 0) const = 0;

    /// Get the boundary curve segment between two points on the boundary, as 
    /// well as the cross-tangent curve.  If the given points are not positioned 
    /// on the same boundary (within a certain tolerance), no curves will be created.
    /// \param pt1 the first point on the boundary, given by the user
    /// \param pt2 the second point on the boundary, given by the user
    /// \param epsilon the tolerance used when determining whether the given points
    ///                are lying on a boundary, and if they do, whether they both lie
    ///                on the \em same boundary.
    /// \param cv upon return, this will point to a newly created curve representing
    ///           the boundary curve between 'pt1' and 'pt2'.  The user assumes ownership
    ///           of this object and is responsible for its deletion.  No curve is created
    ///           if the given points are not found to lie on the same boundary.
    /// \param crosscv upon return, this will point to a newly created curve representing
    ///                the cross-boundary curve between 'pt1' and 'pt2'  The user assumes
    ///                ownership of this object and is responsible for its deletion.
    ///                The direction is outwards from the surface.
    ///                No curve is created if the given points are not found to lie on the
    ///                same boundary.
    /// \param knot_tol tolerance used when working with the knot-vector, to specify how
    ///                 close a parameter value must be to a knot in order to be considered
    ///                 'on' the knot.
    virtual void getBoundaryInfo(Point& pt1, Point& pt2,
				 double epsilon, SplineCurve*& cv,
				 SplineCurve*& crosscv, double knot_tol = 1e-05) const = 0;

    /// Turns the direction of the normal of the surface.
    virtual void turnOrientation() = 0;

    /// Reverses the direction of the basis in input direction.
    /// \param direction_is_u if 'true', the first parameter direction will be reversed,
    ///                       otherwise, the second parameter direction will be reversed
    virtual void reverseParameterDirection(bool direction_is_u) = 0;

    /// Swaps the two parameter directions
    virtual void swapParameterDirection() = 0;

    /// The order of the edge indicators (bottom, right, top, left)
    /// matches the edge_number of edgeCurve().

    /// Query whether any of the four boundary curves are degenerate (zero length) within 
    /// a certain tolerance.  In the below, we refer to 'u' as the first parameter and 'v'
    /// as the second.
    /// \param b 'true' upon return of function if the boundary (v = v_min) is degenerate
    /// \param r 'true' upon return of function if the boundary (v = v_max) is degenerate
    /// \param t 'true' upon return of function if the boundary (u = u_min) is degenerate
    /// \param l 'true' upon return of function if the boundary (u = u_max) is degenerate
    /// \param tolerance boundaries are considered degenerate if their length is shorter
    ///        than this value, given by the user
    /// \return 'true' if at least one boundary curve was found to be degenerate, 'false'
    ///         otherwise.
    virtual bool isDegenerate(bool& b, bool& r,
			      bool& t, bool& l, double tolerance) const = 0;
};


///\}
} // namespace Go


#endif // _GOPARAMSURFACE_H


