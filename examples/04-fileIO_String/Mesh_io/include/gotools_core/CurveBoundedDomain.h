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
#ifndef _GOCURVEBOUNDEDDOMAIN_H
#define _GOCURVEBOUNDEDDOMAIN_H

#include "Domain.h"
#include "RectDomain.h"
#include "CurveLoop.h"
#include "SplineSurface.h"
#include "CurveOnSurface.h"

using std::vector;
using std::pair;
using boost::shared_ptr;

namespace Go
{
///\addtogroup geometry
///\{


/// A 2D parameter domain represented by its boundaries.  The domain is
/// not necessarily connected.
class CurveBoundedDomain : public Domain
{
public:
    /// Constructor generating an empty domain
    CurveBoundedDomain()
    {}

    /// The curve loop must contain either 2D ParamCurve objects or
    /// ?D CurveOnSurface objects.
    /// One of the loops (the first, preferrably) should be ccw, the rest cw.


    /// Constructor generating a domain specified by one or more 2D 
    /// CurveLoops.  
    /// \param loops The loops defining the domain to be created.  The domain
    ///              is defined as the union of the interior of these loops.
    ///              The CurveLoop s must contain either 2D ParamCurve objects
    ///              XD CurveOnSurface objects.
    CurveBoundedDomain(vector<shared_ptr<CurveLoop> > loops);

    /// Constructor generating a (connected) domain specified by a CurveLoop
    /// \param ccw_loop the curve loop specifying the domain.  The domain will
    ///                 represent the interior of this (supposedly 
    ///                 counterclockwise loop).
    CurveBoundedDomain(shared_ptr<CurveLoop> ccw_loop);

    /// Virtual destructor, enables safe inheritance.
    virtual ~CurveBoundedDomain();

    // virtual DomainType domainType() const;
    // @@ Seems that direction of bd cv is not used. Fine by me.

    /// Query whether a given parameter pair is inside the domain or not.
    /// \param point array containing the parameter pair
    /// \param tolerance the tolerance to be used.  In order to be considered
    ///                  'inside', the point must be located inside one of the
    ///                  defining CurveLoop s, as well as being at a distance
    ///                  more than 'tolerance' from any point on that CurveLoop.
    /// \return 'true' if the point is found to be inside the domain, 'false' 
    ///         otherwise.
    virtual bool isInDomain(const Array<double, 2>& point, 
			    double tolerance) const;

    /// check whether a given parameter pair is located \em on the domain boundary
    /// \param point array containing the parameter pair 
    /// \param tolerance the tolerance used.  (how 'far' from the boundary our
    ///                  parameter pair can be and still be considered 'on' the 
    ///                  boundary.
    /// \return 'true' if the parameter pair is considered to be on the boundary,
    ///         'false' otherwise.
    virtual bool isOnBoundary(const Array<double, 2>& point, 
			      double tolerance) const;

    /// Find the parameter pair contained in the domain that is closest (using
    /// Euclidean distance in R^2) to a given parameter pair.  If the given parameter
    /// pair is already in the domain, then  the answer is obviously that parameter pair
    /// \param point the (u,v) parameter pair that we want to find the closest 
    ///              parameter pair to \em inside the CurveBoundedDomain.
    /// \retval clo_pt the resulting closest parameter point.
    /// \param tolerance the tolerance used in defining whether the given point is 
    ///        already inside the domain.
    virtual void closestInDomain(const Array<double, 2>& point,
				 Array<double, 2>& clo_pt,
				 double tolerance) const;

    /// Find the parameter pair on the boundary of the domain that is closest
    /// (using Euclidean distance in R^2) to a given parameter pair.  If the
    /// point is already considered \em on the boundary, then the answer is obviously
    /// the same point.
    /// \param point the (u,v) parameter pair that we want to find to closest parameter
    ///              pair to \em on the CurveBoundedDomain border.
    /// \retval clo_bd_pt the closest point on the border.
    /// \param tolerance the tolerance used in defining whether the given point is
    ////                 on the border.
    virtual void closestOnBoundary(const Array<double, 2>& point,
				   Array<double, 2>& clo_bd_pt,
				   double tolerance) const;

    /// Get a rectangular domain that is guaranteed to contain the CurveBoundedDomain.
    /// \return a RectDomain that contains this CurveBoundedDomain.  If the
    ///         loops defining the CurveBoundedDomain are specified by 2D parameter 
    ///         curves, the RectDomain will be the smallest rectangular domain containing
    ///         all the control points of the boundary curves.  In the opposite case,
    ///         the loops are assumed to be CurveOnSurface s, and the returned RectDomain
    ///         will be that of the underlying surface.
    RectDomain containingDomain() const;

    /// On the 2D parameter plane, consider the (iso)curve defined by fixing one of the
    /// parameters at a certain value.  Now, pick those parts of this curve that are
    /// covered by 'this' CurveBoundedDomain.  Then, lift these (parameter) curves
    /// into 3D space by means of a SplineSurface.  This function will return the
    /// CurveOnSurface s defined by this procedure.
    /// \param pardir this parameter specifies which of the two parameter directions
    ///               that should be 'free' (the other will be fixed).  'pardir' should
    ///               take the value of '1' or '2'.  A value of '1' means that the 
    ///               \em first parameter direction will be free, while '2' means that 
    ///               the \em second parameter direction will be free.
    /// \param parval the parameter value of the fixed parameter
    /// \param tolerance the tolerance when determining which parts of the isoparameter
    ///                  curve are located inside 'this' CurveBoundedDomain.
    /// \param srf The surface used to 'lift' the resulting isoparameter curve intervals
    /// \retval trim_pieces vector containing the resulting CurveOnSurface s.
    void clipWithDomain(int pardir, double parval, 
			double tolerance, shared_ptr<SplineSurface> srf,
			vector<shared_ptr<CurveOnSurface> >& trim_pieces) const;

    // Determine which intervals of a 2D spline curve lies inside the bounded
    // domain.  The start and end parameter values for the inside intervals are
    // found in the vector 'params_start_end_interval'.  Even entries in this 
    // vector marks the start of an inside interval, odd entries marks the end.


    /// Given a curve in the 2D parameter plane, determine those parts of the curve
    /// that are contained inside 'this' CurveBoundedDomain.
    /// \param curve this is the 2D curve that we want to examine
    /// \param tolerance this is the tolerance used when determining which parts of 
    ///                  of 'curve' are inside 'this' CurveBoundedDomain.
    /// \retval params_start_end_interval a vector containing the start- and end 
    ///                                   parameters of the curve segments that were
    ///                                   found to be inside this domain.  
    ///                                   An even indexed entry marks the start parameter
    ///                                   of a curve segment, while the following, odd
    ///                                   indexed entry marks the end parameter of 
    ///                                   the same curve segment.
    void findPcurveInsideSegments(const SplineCurve& curve,
				  double tolerance,
				  vector<double>& params_start_end_interval) const;

private:
    typedef struct intersection_point {
	double par1, par2;
	int pretop;
	
	intersection_point(double p1, double p2, int top)
	{ par1 = p1; par2 = p2; pretop = top;}
	
    } intersection_point;
    
    // Comparisement function to use in std::sort
    static bool par1_compare(const intersection_point& el1,
			     const intersection_point& el2)
    {
	if (el1.par1 < el2.par1)
	    return true;
	else
	    return false;
    }
    

    // We store a set of curve loops
    vector<shared_ptr<CurveLoop> > loops_;

    // Fetch all intervals in one parameter direction
    // going through a specific point lying inside the 
    // bounded domain.
    void getInsideIntervals(int pardir, double parval, double tolerance,
			    vector<pair<double, double> >& insideInts) const;

    // We return a pointer to a parameter curve defining boundary. If loops_
    // consists of CoCurveOnSurface's, the parameter domain curve is returned.
    // Otherwise we make sure that dimension really is 2.
    shared_ptr<ParamCurve> getParameterCurve(int loop_nmb, int curve_nmb) const;

    // Determine which intervals of a 2D spline curve lies inside the bounded 
    // domain.  The start and end parameter values for the inside intervals are 
    // found in the intersection_points in the 'intpt' vector.  Those
    // intersection_points with even indexes marks the start of such an interval,
    // each following (odd-indexed) intersection_point mark the end of that interval.
    void findPcurveInsideSegments(const SplineCurve& curve,
				  double tolerance, 
				  vector<intersection_point>& intpt) const;

};


///\}
} // namespace Go

#endif // _GOCURVEBOUNDEDDOMAIN_H


