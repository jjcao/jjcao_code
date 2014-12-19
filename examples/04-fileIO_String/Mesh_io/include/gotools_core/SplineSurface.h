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
#ifndef _GOSPLINESURFACE_H
#define _GOSPLINESURFACE_H


#include "ParamSurface.h"
#include "BsplineBasis.h"
#include "RectDomain.h"
#include "ScratchVect.h"


namespace Go
{
///\addtogroup geometry
///\{


class Interpolator;
class SplineCurve;
class DirectionCone;

/// \brief SplineSurface provides methodes for storing,
/// reading and manipulating rational and non-rational
/// B-spline surfaces.
///
/// Non-rational B-spline surfaces represented on the form
/// \f[ \sum_{i=1}^{n_1}\sum_{j=1}^{n_2} P_{i,j}
/// B_{i,o_1}(u) B_{j,o_2}(v), \f]
/// where the B-spline coefficients are stored in a vector
/// called coefs. The coefs is stored as
/// \f[ P_{0,0}, P_{0,1},  \ldots
/// P_{n_1, n_2}, \f]
/// where \f$ P_{i,j} \f$ is represented as dim doubles.
///
/// NURBS surfaces are represented on the form
/// \f[ \frac{\sum_{i=1}^{n_1}\sum_{j=1}^{n_2} w_{i,j}P_{i,j}
/// B_{i,o_1}(u) B_{j,o_2}(v)}
/// {\sum_{i=1}^{n_1}\sum_{j=1}^{n_2} w_{i,j} B_{i,o_1}(u) B_{j,o_2}(v)} \f]
/// where \f$ B_{i,o} \f$ is the <em> i</em>'th non-rational B-spline of
/// order \e o. The Projected coefficients are stored in the
/// coefs vector, i.e
/// \f[ w_{0,0}*P_{0,0}, w_{0,1}*P_{0,1}, \ldots
/// w_{n_1, n_2}*P_{n_1, n_2}. \f]
/// In addition the cefficients for the surface in projective
/// space are kept,  i.e
/// \f[ P_{0,0}, w_{0,0},  P_{0,1}, w_{0,1}, \ldots
/// P_{n_1, n_2}, w_{n_1, n_2}. \f]
/// With this representation surface is in the convex hull
/// of the coefficients stored in coefs both for rational
/// and non rational surfaces.
class SplineSurface : public ParamSurface
{
public:
    /// Creates an uninitialized SplineSurface, which can only be assigned to 
    /// or read(...) into.
    SplineSurface()
	: dim_(-1), rational_(false)
    {
    }

    /// Create a SplineSurface by explicitly providing all spline-related 
    /// information.
    /// \param number1 number of control points along the u-parameter
    /// \param number2 number of control points along the v-parameter
    /// \param order1 BsplineBasis order along the u-parameter
    /// \param order2 BsplineBasis order along the v-parameter
    /// \param knot1start pointer to the array describing the knotvector 
    ///                   for the u-parameter
    /// \param knot2start pointer to the array describing the knotvector
    ///                   for the v-parameter
    /// \param coefsstart pointer to the array where the control points
    ///                   are consecutively stored.  The storage order is
    ///                   such that control points along the u-parameter have
    ///                   the shortest stride (stored right after each other).
    ///                   If the surface is rational, pay attention to the
    ///                   comments below.
    /// \param dim        dimension of the space in which the surface lies 
    ///                   (usually 3).  
    /// \param rational Specify whether the surface is rational or not.
    ///                 If the surface is rational, coefficients must be in
    ///                 the following format: 
    ///                 wP1 wP2 .... wPdim w.   Ie. a (dim+1)-dimensional form.
    ///                 (This is the same form that is used within SISL).
    template <typename RandomIterator1,
	      typename RandomIterator2,
	      typename RandomIterator3>
    SplineSurface(int number1,
		    int number2,
		    int order1,
		    int order2,
		    RandomIterator1 knot1start,
		    RandomIterator2 knot2start,
		    RandomIterator3 coefsstart,
		    int dim,
		    bool rational = false)
	: dim_(dim), rational_(rational),
	  basis_u_(number1, order1, knot1start),
	  basis_v_(number2, order2, knot2start)
    {
	if (rational) {
	    int n = (dim+1)*number1*number2;
	    rcoefs_.resize(n);
	    std::copy(coefsstart, coefsstart + n, rcoefs_.begin());
	    coefs_.resize(dim*number1*number2);
	    updateCoefsFromRcoefs();
	} else {
	    int n = dim*number1*number2;
	    coefs_.resize(n);
	    std::copy(coefsstart, coefsstart + n, coefs_.begin());
	}
    }

    /// Create a SplineSurface by explicitly providing all spline-related 
    /// information.
    /// \param basis_u the BsplineBasis to be used along the u-direction
    /// \param basis_v the BsplineBasis to be used along the v-direction
    /// \param coefsstart pointer to the array where the control points
    ///                   are consecutively stored.  The storage order is
    ///                   such that control points along the u-parameter have
    ///                   the shortest stride (stored right after each other).
    ///                   If the surface is rational, pay attention to the
    ///                   comments below.
    /// \param dim dimension of the space in which the surface lies (usually 3).
    /// \param rational Specify whether the surface is rational or not.
    ///                 If the surface is rational, coefficients must be in the
    ///                 following format:
    ///                 wP1 wP2 .... wPdim w.   Ie. a (dim+1)-dimensional form.
    ///                 (This is the same form that is used within SISL).
    template <typename RandomIterator>
    SplineSurface(const BsplineBasis& basis_u,
		    const BsplineBasis& basis_v,
		    RandomIterator coefsstart,
		    int dim,
		    bool rational = false)
	: dim_(dim), rational_(rational),
	  basis_u_(basis_u),
	  basis_v_(basis_v)
    {
	int number1 = basis_u.numCoefs();
	int number2 = basis_v.numCoefs();
	if (rational) {
	    int n = (dim+1)*number1*number2;
	    rcoefs_.resize(n);
	    std::copy(coefsstart, coefsstart + n, rcoefs_.begin());
	    coefs_.resize(dim*number1*number2);
	    updateCoefsFromRcoefs();
	} else {
	    int n = dim*number1*number2;
	    coefs_.resize(n);
	    std::copy(coefsstart, coefsstart + n, coefs_.begin());
	}
    }

    /// Virtual destructor, enables safe inheritance.
    virtual ~SplineSurface();

    // inherited from Streamable
    virtual void read (std::istream& is);

    // inherited from Streamable
    virtual void write (std::ostream& os) const;

    // inherited from GeomObject
    virtual BoundingBox boundingBox() const;

    // inherited from GeomObject
    virtual int dimension() const;
    
    /// quick swap of two SplineSurface objects with each other
    void swap(SplineSurface& other);

    // inherited from GeomObject
    virtual ClassType instanceType() const;

    // inherited from GeomObject
    static ClassType classType()
    { return Class_SplineSurface; }
// #if ((_MSC_VER > 0) && (_MSC_VER < 1300))
//     virtual GeomObject* clone() const
//     { return new SplineSurface(*this); }
// #else
    virtual SplineSurface* clone() const
    { return new SplineSurface(*this); }
// #endif

    // inherited from ParamSurface
// #if ((_MSC_VER > 0) && (_MSC_VER < 1300))
//     virtual const Domain& parameterDomain() const;
// #else
    virtual const RectDomain& parameterDomain() const;
// #endif

    // inherited from ParamSurface
    virtual RectDomain containingDomain() const;

    // inherited from ParamSurface
    virtual CurveLoop outerBoundaryLoop(double degenerate_epsilon
					  = DEFAULT_SPACE_EPSILON) const;
    // inherited from ParamSurface
    virtual std::vector<CurveLoop> allBoundaryLoops(double degenerate_epsilon
						      = DEFAULT_SPACE_EPSILON) const;

    // inherited from ParamSurface
    virtual void point(Point& pt, double upar, double vpar) const;

    // Output: Partial derivatives up to order derivs (pts[0]=S(u,v),
    // pts[1]=dS/du=S_u, pts[2]=S_v, pts[3]=S_uu, pts[4]=S_uv, pts[5]=S_vv, ...)
    // inherited from ParamSurface
    virtual void point(std::vector<Point>& pts, 
		       double upar, double vpar,
		       int derivs,
		       bool u_from_right = true,
		       bool v_from_right = true,
		       double resolution = 1.0e-12) const;

    /// Get the start value for the u-parameter
    /// \return the start value for the u-parameter
    virtual double startparam_u() const;

    /// Get the start value for the v-parameter
    /// \return the start value for the v-parameter
    virtual double startparam_v() const;

    /// Get the end value for the u-parameter
    /// \return the end value for the u-parameter
    virtual double endparam_u() const;

    /// Get the end value for the v-parameter
    /// \return the end value for the v-parameter
    virtual double endparam_v() const;

    // inherited from ParamSurface
    virtual void normal(Point& n, double upar, double vpar) const;

    /// Enumerates the method for computing the normal cone
    enum NormalConeMethod { 
	SederbergMeyers = 0,
	SMCornersFirst = 1
    };
    /// Returns a cone that contains the convex hull of all normalized
    /// tangents of the surface. Note: It is an overestimate.
    /// \return the normal cone of the surface
    DirectionCone normalCone(NormalConeMethod method) const;

    /// Function that calls normalCone(NormalConeMethod) with method =
    /// SederbergMeyers. Needed because normalCone() is virtual! 
    /// (Inherited from ParamSurface).
    /// \return a DirectionCone (not necessarily the smallest) containing all normals 
    ///         to this surface.
    virtual DirectionCone normalCone() const;

    // inherited from ParamSurface
    virtual DirectionCone tangentCone(bool pardir_is_u) const;

    // Creates a composite box enclosing the surface. The composite box
    // consists of an inner and an edge box. The inner box is
    // made from the interior control points of the surface, while the
    // edge box is made from the boundary control points.
    // Inherited from ParamSurface
    virtual CompositeBox compositeBox() const;

    /// Not yet implemented
    SplineSurface* normal() const;

    /// Returns the normal surface corresponding to this surface, as 
    /// described in: \n
    /// <tt> Computing normal vector Bezier patches, Przemyslaw Kiciak,</tt>
    /// <tt> Computer Aided Design 18 (2001), 699-710</tt>.
    /// \return pointer to a newly created SplineSurface which is the 
    ///         normal surface.  User assumes ownership of this object,
    ///         and is responsible for destroying it.
    SplineSurface* normalSurface() const;

    /// Get the derivative surface.
    /// Expresses the i,j-th derivative of a spline surface as
    /// a spline surface. Ported from sisl routine s1386.
    /// \param ider1 what derivative to use along the u parameter
    /// \param ider2 what derivative to use along the v parameter
    /// \return pointer to a newly constructed SplineSurface which is 
    ///         the derivative surface.  User assumes ownership of this
    ///         object, and is responsible for destroying it.
    SplineSurface* derivSurface(int ider1, int ider2) const;

    /// Get a SplineSurface which represent a part of 'this' SplineSurface
    /// \param from_upar start value for u-parameter in the sub-surface to be 
    ///                  generated 
    /// \param from_vpar start value for v-parameter in the sub-surface to be
    ///                  generated
    /// \param to_upar end value for u-parameter in the sub-surface to be generated
    /// \param to_vpar end value for v-parameter in the sub-surface to be generated
    /// \param fuzzy tolerance used to determine whether given parameter values
    ///              are located directly \em knot values.
    SplineSurface* subSurface(double from_upar,
			      double from_vpar,
			      double to_upar,
			      double to_vpar,
			      double fuzzy =
			      DEFAULT_PARAMETER_EPSILON) const;

    // inherited from ParamSurface
    virtual std::vector<boost::shared_ptr<ParamSurface> >
    subSurfaces(double from_upar, double from_vpar,
		double to_upar, double to_vpar,
		double fuzzy = DEFAULT_PARAMETER_EPSILON) const;

    // inherited from ParamSurface
    virtual void closestPoint(const Point& pt,
			      double&        clo_u,
			      double&        clo_v, 
			      Point&         clo_pt,
			      double&        clo_dist,
			      double         epsilon,
			      const RectDomain* domain_of_interest = NULL,
			      double   *seed = 0) const;

    // inherited from ParamSurface
    virtual void closestBoundaryPoint(const Point& pt,
				      double&        clo_u,
				      double&        clo_v, 
				      Point&         clo_pt,
				      double&        clo_dist,
				      double         epsilon,
				      const RectDomain* rd = NULL,
				      double *seed = 0) const;

    /// Appends a surface to the end of 'this' surface along a specified parameter
    /// direction.  The knotvector and control point grid of 'this' surface will be
    /// extended, and the degree will be raised if necessary.  Note that there \em
    /// will be side-effects on the other surface; its order might be raised, its
    /// knotvector will become k-regular and reparametrized, one of its edges will be
    /// moved to coincide with an edge of 'this' surface, etc.
    /// \param sf the surface to append to 'this' surface
    /// \param join_dir the parameter direction along which to extend 'this' surface
    ///                 by appending the 'sf' surface.  If it is equal to 1, we will
    ///                 extend the surface along the u-parameter; if it equals 2, we
    ///                 will extend the surface along the v-parameter.
    /// \param continuity the desired level of continuity at the transition between
    ///                   the two surfaces (can be from -1 to order()).  The higher 
    ///                   the value, the more the surfaces will have to be locally 
    ///                   modified around the seam.
    /// \param dist upon function return, this variable will hold the estimated
    ///             maximum distorsion after the 'smoothing' of the seam in order 
    ///             to achieve the desired continuity.   
    /// \param repar The parametrization of the concerned parameter in the other 
    ///              surface will \em always be shifted so that it starts where the
    ///              parametrization of 'this' surface ends.  However, if 'repar' is 
    ///              set to 'true', it will also be \em scaled as a function of position
    ///              of control points close to the transition.
    void appendSurface(ParamSurface* sf, int join_dir,
		       int continuity, double& dist, bool repar=true);

    void appendSurface(SplineSurface* sf, int join_dir, int continuity, double& dist);

    /// Short hand function to call \ref appendSurface with C^1 continuity.
    /// \param sf the surface to append to 'this' surface
    /// \param join_dir the parameter direction along which to extend 'this' surface
    ///                 by appending the 'sf' surface.  If it is equal to 1, we will
    ///                 extend the surface along the u-parameter; if it equals 2, we
    ///                 will extend the surface along the v-parameter.
    /// \param repar The parametrization of the concerned parameter in the other 
    ///              surface will \em always be shifted so that it starts where the
    ///              parametrization of 'this' surface ends.  However, if 'repar' is 
    ///              set to 'true', it will also be \em scaled as a function of position
    ///              of control points close to the transition.
    void appendSurface(ParamSurface* sf, int join_dir, bool repar=true);

    // inherited from ParamSurface
    virtual void getBoundaryInfo(Point& pt1, Point& pt2, 
				 double epsilon, SplineCurve*& cv,
				 SplineCurve*& crosscv, double knot_tol = 1e-05) const;

    /// Get the boundary curve and the outward cross tangent curve between
    /// two parameter values on a given boundary.
    /// \param par1 first parameter along the boundary
    /// \param par2 second parameter along the boundary
    /// \param bdindex index of the boundary in question.  The boundaries are
    ///                indexed as in \ref getBoundaryIdx().
    /// \param cv returns a pointer to a generated SplineCurve representing the
    ///           boundary between the given parameters.
    ///           The user assumes ownership and is responsible for deletion.
    /// \param crosscv returns a pointer to a generated SplineCurve representing
    ///                the cross tangent curve between the given parameters on the
    ///                boundary. The direction is outwards from the surface.
    ///                The user assumes ownership and is responsible for
    ///                deletion.
    /// \param knot_tol tolerance used when working with the knot-vector, to specify how
    ///                 close a parameter value must be to a knot in order to be considered
    ///                 'on' the knot.
    void getBoundaryInfo(double par1, double par2,
			 int bdindex, SplineCurve*& cv,
			 SplineCurve*& crosscv, double knot_tol = 1e-05) const;

    /// Given two points on the surface boundary, find the number of the
    /// corresponding boundary and the curve parameter of the closest points
    /// on this surface boundary.
    /// 
    /// Ordering of boundaries: 
    /// \verbatim
    ///                      1 
    ///           ---------------------- 
    ///           |                    |
    ///         2 |                    | 3
    ///      v    |                    |
    ///      ^    ----------------------
    ///      |-> u            0
    /// \endverbatim
    /// \param pt1 point in geometry space
    /// \param pt2 point in geometry space
    /// \param epsilon geometric tolerance
    /// \param bdindex if the two points are on a common boundary, 
    ///                the index of the boundary, otherwise -1.
    /// \param par1 if bdindex is 0 or 1, the u parameter of pt1, 
    ///             otherwise the v parameter.
    /// \param par2 if bdindex is 0 or 1, the u parameter of pt2, 
    ///             otherwise the v parameter.
    /// \param knot_tol tolerance used when working with the knot-vector, to specify how
    ///                 close a parameter value must be to a knot in order to be considered
    ///                 'on' the knot.
    void getBoundaryIdx(Point& pt1, Point& pt2, 
			double epsilon, int& bdindex,
			double& par1, double& par2, double knot_tol = 1e-05) const;

    /// Given one point on the surface boundary, find the index number of the
    // corresponding boundary.  The Boundaries are indexed as in \getBoundaryIdx().
    /// \param pt1 the point we want to find the boundary for
    /// \param epsilon the geometric distance the point can have from a boundary and 
    ///                still be considered as laying on the boundary.
    /// \param bdindex the index of the boundary on which the point lies.  If the 
    ///                point was found not to lie on any boundary, the value will 
    ///                be -1.
    /// \param knot_tol tolerance used when working with the knot-vector, to specify how
    ///                 close a parameter value must be to a knot in order to be considered
    ///                 'on' the knot.
    void getBoundaryIdx(Point& pt1, double epsilon, int& bdindex,
			double knot_tol = 1e-05) const;

    // inherited from ParamSurface
    virtual void turnOrientation();

    // inherited from ParamSurface
    virtual void swapParameterDirection();

    // inherited from ParamSurface
    virtual void reverseParameterDirection(bool direction_is_u);

    /** If the two points pt1 and pt2 are on the same edge the
     edge index is returned, otherwise -1 is returned.
     Edges are counted anticlockwise, and starts with the
     edge (upar, 0). */

    /// find whether two points are lying on the same surface boundary; if
    /// this is the case, return the index of the boundary (as specified
    /// in \ref getBoundaryIdx()), else return 1
    /// \return The boundary index on which the points lie, if they are found
    ///         to both lie on the same boundary.  In the opposite case, -1 
    ///         is returned.
    int boundaryIndex(Point& param_pt1, Point& param_pt2) const;

    /// get a const reference to the BsplineBasis for the first parameter
    /// \return const reference to the BsplineBasis for the first parameter
    const BsplineBasis& basis_u() const
    { return basis_u_; }

    /// get a const reference to the BsplineBasis for the second parameter
    /// \return const reference to the BsplineBasis for the second parameter
    const BsplineBasis& basis_v() const
    { return basis_v_; }

    /// get one of the BsplineBasises of the surface
    /// \param i specify whether to return the BsplineBasis for the first 
    ///          parameter (0), or for the second parameter (1).
    /// \return const reference to the requested BsplineBasis.
    const BsplineBasis& basis(int i) const
    { return (i==0) ? basis_u_ : basis_v_; }

    /// Query the number of control points along the first parameter direction
    /// \return the number of control points along the first parameter direction
    int numCoefs_u() const
    { return basis_u_.numCoefs(); }

    /// Query the number of control points along the second parameter direction
    /// \return the number of control points along the second parameter direction
    int numCoefs_v() const
    { return basis_v_.numCoefs(); }

    /// Query the order of the BsplineBasis for the first parameter
    /// \return  the order of the BsplineBasis for the first parameter
    int order_u() const
    { return basis_u_.order(); }

    /// Query the order of the BsplineBasis for the second parameter
    /// \return the order of the BsplineBasis for the second parameter
    int order_v() const
    { return basis_v_.order(); }

    /// Query whether the surface is rational
    /// \return 'true' if the surface is rational, 'false' otherwise
    bool rational() const
    { return rational_; }
    
    /// Get an iterator to the start of the internal array of non-rational control 
    /// points.
    /// \return an (nonconst) iterator to the start of the internal array of non-
    ///         rational control points
    std::vector<double>::iterator coefs_begin()
    { return coefs_.begin(); }

    /// Get an iterator to the one-past-end position of the internal array of non-
    /// rational control points
    /// \return an (nonconst) iterator to the one-past-end position of the internal
    ///         array of non-rational control points
    std::vector<double>::iterator coefs_end()
    { return coefs_.end(); }

    /// Get a const iterator to the start of the internal array of non-rational
    /// control points.
    /// \return a const iterator to the start of the internal array of non-rational
    ///         control points.
    std::vector<double>::const_iterator coefs_begin() const
    { return coefs_.begin(); }

    /// Get a const iterator to the one-past-end position of the internal array of
    /// non-rational control points.
    /// \return a const iterator to the one-past-end position of the internal array of
    ///         non-rational control points.
    std::vector<double>::const_iterator coefs_end() const
    { return coefs_.end(); }

    /// Get an iterator to the start of the internal array of \em rational control 
    /// points.
    /// \return an (nonconst) iterator ro the start of the internal array of rational
    ///         control points.
    std::vector<double>::iterator rcoefs_begin()
    { return rcoefs_.begin(); }

    /// Get an iterator to the one-past-end position of the internal array of 
    /// \em rational control points.
    /// \return an (nonconst) iterator to the start of the internal array of rational
    ///         control points.
    std::vector<double>::iterator rcoefs_end()
    { return rcoefs_.end(); }

    /// Get a const iterator to the start of the internal array of \em rational
    /// control points.
    /// \return a const iterator to the start of the internal array of rational
    ///         control points
    std::vector<double>::const_iterator rcoefs_begin() const
    { return rcoefs_.begin(); }

    /// Get a const iterator to the one-past-end position of the internal array
    /// of \em rational control points.
    /// \return a const iterator to the one-past-end position of the internal array 
    ///         of rational control points.
    std::vector<double>::const_iterator rcoefs_end() const
    { return rcoefs_.end(); }

    // inherited from ParamSurface
    virtual bool isDegenerate(bool& b, bool& r,
		      bool& t, bool& l, double tolerance) const;

    /// set the parameter domain to a given rectangle
    /// \param u1 new min. value of first parameter span
    /// \param u2 new max. value of first parameter span
    /// \param v1 new min. value of second parameter span
    /// \param v2 new max. value of second parameter span
    void setParameterDomain(double u1, double u2, double v1, double v2);

    /// Insert a new knot in the knotvector of the first parameter
    /// \param apar the parameter value at which a new knot will be inserted
    void insertKnot_u(double apar);
    
    /// Insert new knots in the knotvector of the first parameter
    /// \param new_knots a vector containing the parameter values of the
    ///                  new knots to insert.
    void insertKnot_u(const std::vector<double>& new_knots);

    /// Insert a new knot in the knotvector of the second parameter
    /// \param apar the parameter value at which a new knot will be inserted
    void insertKnot_v(double apar);
    
    /// Insert new knots in the knotvector of the second parameter
    /// \param new_knots a vector containing the parameter values of the 
    ///                  new knots to insert.
    void insertKnot_v(const std::vector<double>& new_knots);

    /// Remove a knot from the knotvector of the first parameter.
    /// \param tpar the parameter value of the knot to be removed
    void removeKnot_u(double tpar);

    /// Remove a knot from the knotvector of the second parameter.
    /// \param tpar the parameter value of the knot to be removed.
    void removeKnot_v(double tpar);

    /// Inserts knots in the u knot vector, such that all knots
    /// have multiplicity order
    void makeBernsteinKnotsU();

    /// Inserts knots in the v knot vector, such that all knots
    /// have multiplicity order
    void makeBernsteinKnotsV();

    /// Returns the number of knot intervals in u knot vector.
    /// \return the number of knot intervals in the knotvector for the first 
    ///         parameter
    int numberOfPatches_u() const;

    /// Returns the number of knot intervals in v knot vector.
    /// \return the number of knot intervals in the knotvector for the second
    ///         parameter
    int numberOfPatches_v() const;

    /// Raise the order of the spline surface as indicated by parameters.
    /// \param raise_u the order of the BsplineBasis associated with the first
    ///                parameter will be raised this many times.
    /// \param raise_v the order of the BsplineBasis associated with the second
    ///                parameter will be raised this many times.
    void raiseOrder(int raise_u, int raise_v);

    /// Generate and return a SplineCurve that represents a constant parameter 
    /// curve on the surface
    /// \param parameter value of the fixed parameter
    /// \param pardir_is_u 'true' if the first parameter is the running parameter,
    ///                    'false' otherwise.
    /// \return pointer to a newly constructed SplineCurve representing the 
    ///         specified constant parameter curve.  It is the user's reponsibility
    ///         to delete it when it is no longer needed.
    SplineCurve* constParamCurve(double parameter,
				 bool pardir_is_u) const;

    /// Generate and return two SplineCurves, representing a constant parameter 
    /// curve on the surface as well as it cross-tangent curve.
    /// \param parameter value of the fixed parameter
    /// \param pardir_is_u 'true' if the first parameter is the running parameter,
    ///                    'false' otherwise.
    /// \param cv upon function completion, 'cv' will point to a newly constructed
    ///           SplineCurve representing the specified constant parameter curve.
    ///           It is the user's responsibility to delete it when it is no longer
    ///           needed.
    /// \param crosscv upon function completion, 'crosscv' will point to a newly
    ///                constructed SplineCurve representing the cross-tangent curve
    ///                of the specified constant parameter curve.  It is the user's
    ///                reponsibility to delete it when it is no longer needed.
    void constParamCurve(double parameter, 
			 bool pardir_is_u, 
			 SplineCurve*& cv, 
			 SplineCurve*& crosscv) const;

    // inherited from ParamSurface
    virtual std::vector< boost::shared_ptr<ParamCurve> >
    constParamCurves(double parameter, bool pardir_is_u) const;

    /// Here edge_number means:
    ///  0 -> bottom edge,
    ///  1 -> right edge,
    ///  2 -> top edge and
    ///  3 -> left edge.
    /// The returned pointer is new'ed and the responsibility
    /// of the caller. You should either delete it, or
    /// manage it with a smart pointer.


    /// Generate and return a SplineCurve representing one of the surface's four
    /// edges.
    /// \param ccw_edge_number indicates which edge the user requests. 
    /// \verbatim 
    /// 0 -> bottom edge
    /// 1 -> right edge
    /// 2 -> top edge
    /// 3 -> left edge \endverbatim
    /// \return A pointer to a newly constructed SplineCurve representing the
    ///         requested edge.  It is the user's responsibility to delete it when
    ///         it is no longer needed.
    SplineCurve* edgeCurve(int ccw_edge_number) const;

    /// Remake 'this'surface to interpolate (or approximate)  a given point grid.
    /// \param interpolator1 Interpolator to apply to the first parameter direction
    /// \param interpolator2 Interpolator to apply to the second parameter direction
    /// \param num_points1 number of points to interpolate along the first parameter 
    ///                    direction
    /// \param num_points2 number of points to interpolate along the second parameter
    ///                    direction
    /// \param dim dimension of the points to interpolate (usually 3)
    /// \param param1_start pointer to the start of the input grid's parametrization
    ///                     along the first parameter
    /// \param param2_start pointer to the start of the input grid's parametrization
    ///                     along the second parameter
    /// \param data_start pointer to the start of the storage array for the data point
    ///                   grid to interpolate
    void interpolate(Interpolator& interpolator1,
		     Interpolator& interpolator2,
		     int num_points1,
		     int num_points2,
		     int dim,
		     const double* param1_start,
		     const double* param2_start,
		     const double* data_start);

    /// Evaluate points and normals on an entire grid, taking computational advantage
    /// over calculating all these values simultaneously rather than one-by-one.
    /// \param num_u number of values to evaluate along first parameter direction
    /// \param num_v number of values to evaluate along second parameter direction
    /// \param points upon function return, this vector holds all the evaluated points
    /// \param normals upon function return, this vector holds all the evaluated normals
    /// \param param_u upon function return, this vector holds all the numerical values
    ///                for the first parameter where evaluation has taken place
    /// \param param_v upon function return, this vector holds all the numerical values
    ///                for the second parameter where evaluation has taken place.
    void gridEvaluator(int num_u, int num_v,
		       std::vector<double>& points,
		       std::vector<double>& normals,
		       std::vector<double>& param_u,
		       std::vector<double>& param_v) const;

    /// Evaluate points on an entire grid, taking computational advantage
    /// over calculating all these values simultaneously rather than one-by-one.
    /// \param num_u number of values to evaluate along first parameter direction
    /// \param num_v number of values to evaluate along second parameter direction
    /// \param points upon function return, this vector holds all the evaluated points
    /// \param param_u upon function return, this vector holds all the numerical values
    ///                for the first parameter where evaluation has taken place
    /// \param param_v upon function return, this vector holds all the numerical values
    ///                for the second parameter where evaluation has taken place.
    void gridEvaluator(int num_u, int num_v,
		       std::vector<double>& points,
		       std::vector<double>& param_u,
		       std::vector<double>& param_v) const;

    // inherited from ParamSurface
    virtual double nextSegmentVal(int dir, double par, bool forward, double tol) const;


private:

    struct degenerate_info
    {
	bool is_set_;
	bool b_, r_, t_, l_;
	
	degenerate_info()
	{ is_set_ = false; }
    };

    // Canonical data
    int dim_;
    bool rational_;
    BsplineBasis basis_u_;
    BsplineBasis basis_v_;
    std::vector<double> coefs_;   /// Like ecoef in SISL
    std::vector<double> rcoefs_;  /// Like rcoef in SISL, only used if rational

    // Generated data
    mutable RectDomain domain_;
    mutable CurveLoop spatial_boundary_;
    mutable degenerate_info degen_;


    // Helper functions
    void updateCoefsFromRcoefs();
    std::vector<double>& activeCoefs() { return rational_ ? rcoefs_ : coefs_; }
    bool normal_not_failsafe(Point& n, double upar, double vpar) const;
    bool search_for_normal(bool interval_in_u,
			   double fixed_parameter,
			   double interval_start, // normal is not defined here
			   double interval_end,
			   Point& normal) const;

    // Members new to this class
    void pointsGrid(int m1, int m2, int derivs,
		    const double* basisvals1,
		    const double* basisvals2,
		    const int* knotint1,
		    const int* knotint2,
		    double* result,
		    double* normals = 0) const;

    // Rewritten pointsGrid, to avoid reformatting results.
    // Does not return the derivatives, works only for a 3D non-rational
    // spline.
    void pointsGridNoDerivs(int m1, int m2,
			    const double* basisvals1,
			    const double* basisvals2,
			    const int* knotint1,
			    const int* knotint2,
			    double* result,
			    double* normals = 0) const;

    // Actually computes the closest point. Only difference is the explicit
    // robust_seedfind-parameter, which is always true in the virtual
    // closestPoint() function, when that function calls closestPointImpl().
    void closestPointImpl(const Point& pt,
			  double&        clo_u,
			  double&        clo_v, 
			  Point&         clo_pt,
			  double&        clo_dist,
			  double         epsilon,
			  const RectDomain* domain_of_interest = NULL,
			  bool robust_seedfind = true,
			  double   *seed = 0) const;


};


///\}
} // namespace Go




#endif // _GOSPLINESURFACE_H


