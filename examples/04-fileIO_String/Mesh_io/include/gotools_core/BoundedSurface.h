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
#ifndef _GOBOUNDEDSURFACE_H
#define _GOBOUNDEDSURFACE_H

#include "ParamSurface.h"
#include "CurveOnSurface.h"
#include "CurveBoundedDomain.h"

using boost::shared_ptr;
using std::vector;

namespace Go
{
///\addtogroup geometry
///\{


/// The class representing trimmed surfaces in Go
class BoundedSurface : public ParamSurface
{
public:

    /// Create an empty BoundedSurface that can be assigned or read() into
    BoundedSurface();

    /// Create a BoundedSurface by specifying the underlying surface and a loop
    /// of curves that specifies the trimming of the surface.
    /// \param surf the created BoundedSurface will represent a trimmed version of 
    ///             this surface
    /// \param loop a vector of CurveOnSurface s that together describe the boundary
    ///             that defines the trimming of the  surface.  The curves in this
    ///             vector should all lie on the surface in question, and when placed
    ///             head-to-tail they should form a closed loop with counterclockwise
    ///             orientation.
    /// \param space_epsilon geometrical tolerance used when treating the loops.
    BoundedSurface(shared_ptr<ParamSurface> surf,
                     vector<shared_ptr<CurveOnSurface> > loop,
		     double space_epsilon);


    /// Create a BoundedSurface by specifying the underlying surface and a number of 
    /// loops of curves that specify the trimming of the surface.
    /// \param surf the created BoundedSurface will represent a trimmed version of 
    ///             this surface
    /// \param loops each entry in 'loop is a vector of CurveOnSurface s that describe 
    ///             a closed loop forming a part of the trimmed surface's boundary.  (Since
    ///             the surface may have internal holes, more than one loop might be 
    ///             required to describe its boundary).  The curve loops should all
    ///             lie on the surface in question.  The first entry in this vector
    ///             describes the outermost boundary, which should be oriented 
    ///             counterclockwise.  The other entries represent holes, and should
    ///             be oriented clockwise.
    /// \param space_epsilon geometrical tolerance used when treating the loops.
    BoundedSurface(shared_ptr<ParamSurface> surf,
		     vector<vector<shared_ptr<CurveOnSurface> > > loops,
		     double space_epsilon);

    /// Virtual destructor ensures safe inheritance
    virtual ~BoundedSurface();


    // From Streamable

    // @afr: These should not be called!
    /// read this BoundedSurface from a stream
    virtual void read (std::istream& is);
    /// write this BoundedSurface to a stream
    virtual void write (std::ostream& os) const;

    // From GeomObject

    /// Return the object's bounding box
    virtual BoundingBox boundingBox() const;

    /// Return the dimension of the space in which the object lies (usually 2 or 3)    
    virtual int dimension() const;

    /// Return the class type identifier of type BoundedSurface
    virtual ClassType instanceType() const;

    /// Return the class type identifier of type BoundedSurface
    static ClassType classType()
    { return Class_BoundedSurface; }

    /// clone this BoundedSurface and return a pointer to the clone
    virtual BoundedSurface* clone() const
    { return new BoundedSurface(*this); }

    // From ParamSurface

    /// Creates a DirectionCone covering all normals to this surface.
    /// \return a DirectionCone (not necessarily the smallest) containing all normals 
    ///         to this surface.
    virtual DirectionCone normalCone() const;

    /// Creates a DirectionCone covering all tangents to 
    /// this surface along a given parameter direction.
    /// \param pardir_is_u if 'true', then the DirectionCone will be defined on basis 
    ///        of the surface's tangents along the first parameter direction.  Otherwise
    ///        the second parameter direction will be used.
    /// \return a DirectionCone (not necessarily the smallest) containing all tangents
    ///         to this surface along the specified parameter direction.
    virtual DirectionCone tangentCone(bool pardir_is_u) const;

    /// Return the parameter domain of the surface.  This may be a simple
    /// rectangular domain (\ref RectDomain) or any other subclass of
    /// \ref Domain (such as CurveBoundedDomain, found in the
    /// \c sisl_dependent module).
    /// \return a Domain object describing the parametric domain of the surface
    virtual const CurveBoundedDomain& parameterDomain() const;

    /// Get a rectangular parameter domain that is guaranteed to contain the
    /// surface's \ref parameterDomain().  It may be the same.  There is no
    /// guarantee that this is the smallest domain containing the actual domain.
    /// \return a RectDomain that is guaranteed to include the surface's total
    ///         parameter domain.
    virtual RectDomain containingDomain() const;


    /// Returns the anticlockwise, outer boundary loop of the surface.
    /// \param degenerate_epsilon edges whose length is smaller than this value
    ///        are ignored.
    /// \return a CurveLoop describing the anticlockwise, outer boundary loop of 
    ///         the surface.
    virtual CurveLoop outerBoundaryLoop(double degenerate_epsilon
					  = DEFAULT_SPACE_EPSILON) const;

    /// Returns the anticlockwise outer boundary loop of the surface, together with 
    /// clockwise loops of any interior boundaries, such that the surface always is
    /// 'to the left of' the loops.
    /// \param degenerate_epsilon edges whose length is smaller than this value are
    ///                           ignored.
    /// \return a vector containing CurveLoops.  The first of these describe the 
    ///         outer boundary of the surface (counterclockwise), whereas the others 
    ///         describe boundaries of interior holes (clockwise).
    virtual std::vector<CurveLoop> allBoundaryLoops(double degenerate_epsilon
						      = DEFAULT_SPACE_EPSILON) const;

    /// Returns the anticlockwise outer boundary loop of the surface, together with 
    /// clockwise loops of any interior boundaries, such that the surface always is
    /// 'to the left of' the loops.  This function works like \ref allBoundaryLoops(),
    /// except that it includes degenerate edges.
    /// \return vector containing CurveLoops.  The first of these describe the
    ///         outer boundary of the surface (counterclockwise) whereas the others
    ///         describe boundaries of interior holes (clockwise).
    std::vector<CurveLoop> absolutelyAllBoundaryLoops() const;


    /// Evaluates the surface's position for a given parameter pair.
    /// \param pt the result of the evaluation is written here 
    /// \param upar the first parameter
    /// \param vpar the second parameter
    virtual void point(Point& pt, double upar, double vpar) const;


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
		       double resolution = 1.0e-12) const;

    /// Evaluates the surface's position and a certain number of derivatives for 
    /// a given parameter pair.
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
    virtual void point(std::vector<Point>& pts, 
		       double upar, double vpar,
		       int derivs) const;

    //    using ParamSurface::point;
    /// Evaluates the surface normal for a given parameter pair
    /// \param n the computed normal will be written to this variable
    /// \param upar the first parameter
    /// \param vpar the second parameter
    virtual void normal(Point& n, double upar, double vpar) const;

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
    constParamCurves(double parameter, bool pardir_is_u) const;

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
		double fuzzy = DEFAULT_PARAMETER_EPSILON) const;

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
    virtual double nextSegmentVal(int dir, double par, bool forward, double tol) const;

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
			      double   *seed = 0) const;


    /// Iterates to the closest point to pt on the boundary of the surface.
    /// \see closestPoint()
    virtual void closestBoundaryPoint(const Point& pt,
				      double&        clo_u,
				      double&        clo_v, 
				      Point&       clo_pt,
				      double&        clo_dist,
				      double epsilon,
				      const RectDomain* domain_of_interest = NULL,
				      double *seed = 0) const;

    /// Get the boundary curve segment between two points on the boundary, as 
    /// well as the cross-tangent curve.  If the given points are not positioned 
    /// on the same boundary (within a certain tolerance), no curves will be created.
    /// \b NB: This function has not yet been implemented!
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
				 SplineCurve*& crosscv, double knot_tol = 1e-05) const;


    /// Get the boundary curve segment between two points on the same boundary loop.
    /// If the given points are not positioned on the same boundary loop (within a certain
    /// tolerance), no curves will be retuned.
    /// \param pt1 the first point on the boundary, given by the user
    /// \param pt2 the second point on the boundary, given by the user
    /// \retval bd_cvs upon return, this will contain shared pointers to curves that, 
    ///                taken consecutively, describe the requested boundary segment in 
    ///                its entirety.
    void getBoundaryInfo(Point& pt1, Point& pt2,
			 vector<shared_ptr<CurveOnSurface> >& bd_cvs) const;

    /// Turns the direction of the normal of the surface.
    virtual void turnOrientation();

    /// Reverses the direction of the basis in input direction.
    /// \b NB: This function has not yet been implemented!
    /// \param direction_is_u if 'true', the first parameter direction will be reversed,
    ///                       otherwise, the second parameter direction will be reversed
    virtual void reverseParameterDirection(bool direction_is_u);

    // If a segment in a loop in boundary_loops_ is not G1, curve is split.
    // Function to be called prior to a topology builder relying on smooth segments.

    /// This function processes all the curves that participate in defining the surface's
    /// (trimmed) boundary.  Those curves that are not G1 within a certain tolerance
    /// are split into severalcurves, so that all G1-discontinuities will end up \em between
    /// consecutive curve segments.
    /// \param kink the tolerance to use for checking G1 continuity
    void makeBoundaryCurvesG1(double kink);

    /// Swaps the two parameter directions
    virtual void swapParameterDirection();

    /// Change the parameter domain of the underlying surface, and modify the
    /// boundary loops with respect to this change
    /// \param u1 new start value of first parameter
    /// \param u2 new end value of first parameter
    /// \param v1 new start value of second parameter
    /// \param v2 new end value of second parameter
    void setParameterDomain(double u1, double u2, double v1, double v2);

    // If a boundary loop is represented as a single curve, it is split into 3 parts.
    // Handy tue to current limitations in topology analysator.

    /// Split all boundary loops defined by only \em curve into three parts.
    /// (This somewhat exotic member function is included due to its handiness with the
    /// GoTools topology analysator).
    void splitSingleLoops();

    // Access functions

    /// Get a pointer to the underlying surface
    /// \return shared pointer to the underlying surface
    shared_ptr<ParamSurface> underlyingSurface()
    { return surface_; }

    /// Get a pointer to the underlying surface
    /// \return shared pointer to the underlying surface
    shared_ptr<const ParamSurface> underlyingSurface() const
    { return surface_; }

    /// Get the number of boundary loops that describe the trimmed surface.
    int numberOfLoops()
      { return boundary_loops_.size(); }

    /// Get a shared pointer to a specific boundary loop
    shared_ptr<CurveLoop> loop(int idx)
      { return boundary_loops_[idx]; }

    /// Get the space-curve resulting from fixing one of the surface's parameters and 
    /// moving the other along its allowed range (inside the trimmed domain).  If this
    /// results in several disjoint curves, an exception is thrown.
    /// \param parameter the parameter value of the fixed parameter
    /// \param direction_is_u if 'true' then the "free" parameter will be the first one,
    ///                       and the second parameter will be fixed.  If 'false', it is
    ///                       the other way around.
    /// \return a newly created SplineCurve representing the requested space-curve.
    ///         The ownership is assumed by the user.
    SplineCurve* constParamCurve(double parameter, bool direction_is_u) const;

    /// Query whether any of the four boundaries of the \em underlying \em surface
    /// are degenerate (zero length) within a certain tolerance.  In the below, we refer
    /// to 'u' as the first parameter and 'v' as the second.
    /// \param b 'true' upon return of function if the boundary (v = v_min) is degenerate
    /// \param r 'true' upon return of function if the boundary (v = v_max) is degenerate
    /// \param t 'true' upon return of function if the boundary (u = u_min) is degenerate
    /// \param l 'true' upon return of function if the boundary (u = u_max) is degenerate
    /// \param tolerance boundaries are considered degenerate if their length is shorter
    ///        than this value, given by the user
    /// \return 'true' if at least one boundary curve was found to be degenerate, 'false'
    ///         otherwise.
    virtual bool isDegenerate(bool& b, bool& r,
			      bool& t, bool& l, double tolerance) const;


private:
    /// The underlying surface
    shared_ptr<ParamSurface> surface_;

    /// The curves describing the boundaries of the surface.
    /// First element is the outer boundary loop (ordering is done by constructor).
    vector<shared_ptr<CurveLoop> > boundary_loops_;

    mutable CurveBoundedDomain domain_;

    /// Helper functions
    /// Order boundary_loops_ w/outer boundary loop first. Called from constructor.
    void orderBoundaryLoops(double degenerate_epsilon = DEFAULT_SPACE_EPSILON);

    // We want the boundary curve to be at least c1. To be called from public function.
    vector<shared_ptr<CurveOnSurface> >
    splitIntoC1Curves(shared_ptr<CurveOnSurface>& curve, double space_epsilon, double kink);

    // Run through the bounday loops, returning the smallest epsgeo.
    double getEpsGeo() const;

};


///\}
} // namespace Go


#endif // _GOBOUNDEDSURFACE_H


