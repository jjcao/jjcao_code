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
#ifndef _GOSPLINECURVE_H
#define _GOSPLINECURVE_H


#include "DirectionCone.h"
#include "ParamCurve.h"
#include "BsplineBasis.h"


namespace Go
{
///\addtogroup geometry
///\{


class Interpolator;

/// \brief SplineCurve provides methodes for storing, reading and
/// manipulating rational and non-rational B-spline curves.
class SplineCurve : public ParamCurve
{
public:

    /// Creates an uninitialized SplineCurve, which can only be
    /// assigned to or read() into.
    SplineCurve() : dim_(-1), rational_(false) { }

    /// Create a SplineCurve by explicitly providing all
    /// spline-related information.
    /// \param number number of control points
    /// \param order order of b-spline basis
    /// \param knotstart pointer to the array describing the
    /// knotvector
    /// \param coefsstart pointer to the array where the control
    /// points are consecutively stored.
    /// \param dim dimension of the space in which the curve lies
    /// (usually 2 or 3).
    /// \param rational Specify whether the curve is rational or not.
    /// If the curve is rational, coefficients must be in the
    /// following format: wP1 wP2 .... wPdim w.  I.e., a
    /// (dim+1)-dimensional form.  (This is the same form that is used
    /// within SISL).
    template <typename RandomIterator1, typename RandomIterator2>
    SplineCurve(int number,
		  int order,
		  RandomIterator1 knotstart,
		  RandomIterator2 coefsstart,
		  int dim,
		  bool rational = false)
	: dim_(dim), rational_(rational),
	  basis_(number, order, knotstart)
    {
	if (rational) {
	    rcoefs_.resize((dim+1)*number);
	    std::copy(coefsstart, coefsstart+(dim+1)*number, 
		      rcoefs_.begin());
	    coefs_.resize(dim*number);
	    updateCoefsFromRcoefs();
	} else {
	    coefs_.resize(dim*number);
	    std::copy(coefsstart, coefsstart + dim*number, 
		      coefs_.begin());
	}
    }

    /// Construct a straight curve between two points (order 2 -
    /// linear).  The parameter span will be from 0 to T, where T is
    /// the euclidean distance between the two points.
    /// \param pnt1 first of the two points defining the straight
    /// curve
    /// \param pnt2 second of the two points defining the straight
    /// curve
    SplineCurve(const Point& pnt1, const Point& pnt2);

    /// Make a straight curve between two points (order 2 - linear).
    /// The parameter span is also given by the user.
    /// \param pnt1 first of the two points defining the straight
    /// curve
    /// \param startpar start parameter of curve
    /// \param pnt2 second of the two points defining the straight
    /// curve
    /// \param endpar end parameter of curve
    SplineCurve(const Point& pnt1, double startpar, 
		const Point& pnt2, double endpar);

    /// Virtual destructor, enables safe inheritance.
    virtual ~SplineCurve();

    // Inherited from Streamable
    virtual void read (std::istream& is);

    // Inherited from Streamable
    virtual void write (std::ostream& os) const;

    // Inherited from GeomObject
    virtual BoundingBox boundingBox() const;

    // Inherited from GeomObject
    virtual int dimension() const;

    // Inherited from GeomObject
    virtual ClassType instanceType() const;

    // Inherited from GeomObject
    static ClassType classType()
    { return Class_SplineCurve; }
    // Inherited from GeomObject
// #if ((_MSC_VER > 0) && (_MSC_VER < 1300))
//     virtual GeomObject* clone() const
//     { return new SplineCurve(*this); }
// #else
    virtual SplineCurve* clone() const
    { return new SplineCurve(*this); }
// #endif

    // Inherited from ParamCurve
    virtual void point(Point& pt, double tpar) const;

    // Inherited from ParamCurve
    virtual void point(std::vector<Point>& pts, 
		       double tpar,
		       int derivs,
		       bool from_right = true) const;

    // Inherited from ParamCurve
    virtual double startparam() const;

    // Inherited from ParamCurve
    virtual double endparam() const;

    // Inherited from ParamCurve
    virtual void reverseParameterDirection(bool switchparam = false);

    // Inherited from ParamCurve
    virtual SplineCurve* geometryCurve();

    // Inherited from ParamCurve
    virtual DirectionCone directionCone() const;

    // Inherited from ParamCurve
    virtual CompositeBox compositeBox() const;

    // Inherited from ParamCurve
    virtual bool isDegenerate(double degenerate_epsilon);

    /// Make a curve expressing the i'th derivative of 'this' curve,
    /// and return a pointer to it.
    /// \param ider the number of the derivative of which we want to
    /// make a curve
    /// \return a pointer to the newly generated curve.  User assumes
    /// ownership and is responsible for deletion.

    SplineCurve* derivCurve(int ider) const;

    // Inherited from ParamCurve
#if ((_MSC_VER > 0) && (_MSC_VER < 1300))
    virtual ParamCurve* subCurve(double from_par, double to_par,
				 double fuzzy 
				 = DEFAULT_PARAMETER_EPSILON) const;
#else
    virtual SplineCurve* subCurve(double from_par, double to_par,
				  double fuzzy 
				  = DEFAULT_PARAMETER_EPSILON) const;
#endif

    // Inherited from ParamCurve
    virtual void closestPoint(const Point& pt,
			      double         tmin,
			      double         tmax,
			      double&        clo_t,
			      Point&         clo_pt,
			      double&        clo_dist,
			      double const   *seed = 0) const;

    /// Inherited from ParamCurve.
    /// Appends a curve to end of this curve.  The knotvector and
    /// control point vector of 'this' curve will be extended, and the
    /// degree will be raised if necessary.  Note that there \em will
    /// be side-effects on the 'other_curve' (its order might be
    /// raised, its knotvector will become k-regular and evt.
    /// reparametrized, its start point will be moved to coincide with
    /// the end point of 'this' curve, etc.)
    ///
    /// \param other_curve the curve to append to this curve
    /// \param continuity which level of continuity we demand at the
    /// transition between the two curves (can be from -1 to order(),
    /// but the higher the value the more the curves will have to be
    /// locally modified.
    /// \param dist upon function return, this variable will hold the
    /// estimated maximum distorsion after 'smoothing' of joined curve
    /// to achieve the desired continuity.
    /// \param repar The parametrization of the 'other_curve' will \em
    /// always be shifted so that it starts where the parametrization
    /// of 'this' curve ends.  However, if 'repar' is set to 'true',
    /// it will also be \em scaled as a function of position of
    /// control points close to the transition.
    virtual void appendCurve(ParamCurve* other_curve,
			     int continuity, 
			     double& dist, 
			     bool repar=true);
    
    /// Inherited from ParamCurve.
    /// Short hand function to call \ref appendCurve with C^1
    /// continuity.
    /// \param cv the curve to append to this curve
    /// \param repar The parametrization of the 'other_curve' will \em
    /// always be shifted so that it starts where the parametrization
    /// of 'this' curve ends.  However, if 'repar' is set to 'true',
    /// it will also be \em scaled as a function of position of
    /// control points close to the transition.
    virtual void appendCurve(ParamCurve* cv, bool repar=true);

    /// Get a const reference to the BsplineBasis of the curve
    /// \return const reference to the curve's BsplineBasis.
    const BsplineBasis& basis() const
    { return basis_; }

    /// Get a reference to the BsplineBasis of the curve
    /// \return reference to the curve's BsplineBasis.
    BsplineBasis& basis()
    { return basis_; }

    /// Query the number of control points of the curve
    /// \return the number of control points of the curve.
    int numCoefs() const
    { return basis_.numCoefs(); }

    /// Query the order of the spline space in which the curve lies
    /// \return the order of the curve's spline space
    int order() const
    { return basis_.order(); }

    /// Query whether or not the curve is rational.
    /// \return 'true' if the curve is rational, 'false' otherwise.
    bool rational() const
    { return rational_; }

    /// Get an iterator to the beginning of the knot vector
    /// \return an iterator to the beginning of the knot vector
    std::vector<double>::iterator knotsBegin()
    { return basis_.begin(); }
    /// Get a one-past-end iterator to the knot vector
    /// \return an iterator to one-past-end of the knot vector
    std::vector<double>::iterator knotsEnd()
    { return basis_.end(); }
    /// Get a const iterator to the beginning of the knot vector
    /// \return a const iterator to the beginning of the knot vector
    std::vector<double>::const_iterator knotsBegin() const
    { return basis_.begin(); }
    /// Get a one-past-end const iterator to the knot vector
    /// \return a const iterator to one-past-end of the knot vector
    std::vector<double>::const_iterator knotsEnd() const
    { return basis_.end(); }


    /// Get an iterator to the start of the curve's internal,
    /// non-rational control point array
    /// \return an iterator to the start of the curves non-rational
    /// control point array
    std::vector<double>::iterator coefs_begin() 
    { return coefs_.begin(); }
    /// Get a one-past-end iterator to the curve's non-rational,
    /// internal control point array
    /// \return an iterator to one-past-end of the curve's
    /// non-rational, internal control point array
    std::vector<double>::iterator coefs_end() 
    { return coefs_.end(); }
    /// Get a const iterator to the start of the curve's non-rational,
    /// internal control point array
    /// \return a const iterator to the start of the curve's
    /// non-rational control point array.
    std::vector<double>::const_iterator coefs_begin() const 
    { return coefs_.begin(); }
    /// Get a one-past-end const iterator to the curve's non-rational,
    /// internal control point array
    /// \return a const iterator to one-past-end of the curve's
    /// non-rational, internal control point array
    std::vector<double>::const_iterator coefs_end() const 
    { return coefs_.end(); }
    /// Get an iterator to the start of the curve's internal, rational
    /// control point array
    /// \return an iterator to the start of the curves rational
    /// control point array
    std::vector<double>::iterator rcoefs_begin() 
    { return rcoefs_.begin(); }
    /// Get a one-past-end iterator to the curve's rational, internal
    /// control point array
    /// \return an iterator to one-past-end of the curve's rational,
    /// internal control point array
    std::vector<double>::iterator rcoefs_end() 
    { return rcoefs_.end(); }
    /// Get a const iterator to the start of the curve's rational,
    /// internal control point array
    /// \return a const iterator to the start of the curve's rational
    /// control point array.
    std::vector<double>::const_iterator rcoefs_begin() const 
    { return rcoefs_.begin(); }
    /// Get a one-past-end const iterator to the curve's rational,
    /// internal control point array
    /// \return a const iterator to one-past-end of the curve's
    /// rational, internal control point array
    std::vector<double>::const_iterator rcoefs_end() const 
    { return rcoefs_.end(); }

    /// Remake 'this' SplineCurve to interpolate (or approximate) a
    /// given sequence of data points.
    ///
    /// \param interpolator reference to the Interpolator object
    /// specifying the interpolation method to use
    /// \param num_points number of points to interpolate
    /// \param dim dimension of the Euclidean space in which the
    /// points lie (usually 2 or 3)
    /// \param param_start pointer to the start of an array expressing
    /// the parameter values of the given data points.
    /// \param data_start pointer to the array where the coordinates
    /// of the data points are stored.
    void interpolate(Interpolator& interpolator,
		     int num_points,
		     int dim,
		     const double* param_start,
		     const double* data_start);

    /// Insert a new knot into the curve's knotvector
    /// \param apar parameter value of the new knot
    void insertKnot(double apar);

    /// Insert several knots into the curve's knotvector
    /// \param new_knots vector containing the parameter values of the
    /// new knots to be inserted into the knotvector.
    void insertKnot(const std::vector<double>& new_knots);

    /// Insert knots into the knotvector such that all knots get
    /// multiplicity equal to the order of the b-spline basis.
    void makeBernsteinKnots();

    /// Rescale the knotvector so that the total parameter span now
    /// ranges from 't1' to 't2'.
    /// \param t1 new start value of the spline's parameter span
    /// \param t2 new end value of the spline's parameter span
    void setParameterInterval(double t1, double t2);

    /// Remove a knot from the knotvector
    /// \param tpar the parameter value of the knot to be removed.
    void removeKnot(double tpar);

    /// Raise the order of the curve's b-spline basis without changing
    /// the shape of the curve
    /// \param i specifies how many times the order will be raised.
    void raiseOrder(int i = 1);

    /// Make the knotvector k-regular at start.  Useful when
    /// k-regularity of a knotvector is to be assumed.
    void makeKnotStartRegular();

    /// Make the knotvector k-regular at end.  Useful when
    /// k-regularity of a knotvector is to be assumed.
    void makeKnotEndRegular();

    /// quick swap of 'this' SplineCurve with the 'other' one.
    /// \param other the SplineCurve to swap with 'this' one.
    void swap(SplineCurve& other);

    /// Inherited from ParamCurve.  Returns the value of the next
    /// knot.
    /// \param par the parameter from which we will look for the start
    /// of the next interval.
    /// \param forward specify whether we will look forwards or
    /// backwards.
    /// \param tol a tolerance specifying how close 'par' has to be to
    /// a knot to be considered 'on' the knot.
    virtual double nextSegmentVal(double par, bool forward, 
				  double tol) const;


private:
    // Canonical data
    int dim_;
    bool rational_;
    BsplineBasis basis_;
    std::vector<double> coefs_;   /// Like ecoef in SISL
    std::vector<double> rcoefs_;  /// Like rcoef in SISL, only used if
				  /// rational

    // Helper functions
    void updateCoefsFromRcoefs();
    /// Appends this curve to itself in a periodic fashion - that is,
    /// assuming the curve has a periodic structure wrt knots and
    /// coefs, the curve will be wound twice.
    ///
    /// This is a helper for subCurve().
    ///
    /// The caller is assumed to know that the curve is periodic in
    /// the sense tested by analyzePeriodicity() from GeometryTools.h,
    /// but for now we also test inside this function. This test may
    /// go away later.
    void appendSelfPeriodic();

};


///\}
} // namespace Go



#endif // _GOSPLINECURVE_H


