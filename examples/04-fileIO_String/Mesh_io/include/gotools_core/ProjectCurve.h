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
#ifndef _GOPROJECTCURVE_
#define _GOPROJECTCURVE_

#include <boost/smart_ptr.hpp>

#include "Point.h"
#include "EvalCurve.h"
#include "SplineCurve.h"
#include "SplineSurface.h"


namespace Go
{
///\addtogroup geometry
///\{


/// This class represents the curve obtained by projecting a 
/// given 3D curve onto a given part of a given 3D surface.

class ProjectCurve : public EvalCurve
{
public:

    /// Constructor, taking one 3D curve and one 3D surface.  The user 
    /// may additionally provide explicit values for the start and end points
    /// of the curve (supposedly on the surface) as well as specifying
    /// the parameter domain of interest of the surface.
    /// \param space_crv the 3D space curve that will be projected onto a surface
    /// \param surf the 3D surface that we will project the space curve onto
    /// \param start_par_pt explicit position of start point (can be a zero pointer,
    ///                     in which case the start point will be evaluated by projection,
    ///                     just like any other point).
    /// \param end_par_pt explicit position of end point (can be a zero pointer,
    ///                   in which case the start point will be evaluated by projection,
    ///                   just like any other point).
    /// \param epsgeo geometric tolerance to use when projecting curve onto surface, and 
    ///               when using the approximationOK() function.
    /// \param domain_of_interest if the user wants to limit the surface to a certain 
    ///                           parametric domain, it can be specified here.
    ProjectCurve(boost::shared_ptr<Go::SplineCurve>& space_crv, 
		   boost::shared_ptr<Go::SplineSurface>& surf,
		   boost::shared_ptr<Go::Point>& start_par_pt, 
		   boost::shared_ptr<Go::Point>& end_par_pt,
		   double epsgeo,
		   const RectDomain* domain_of_interest = NULL);

    /// virtual destructor ensures safe inheritance
    virtual ~ProjectCurve();
    
    // Inherited from EvalCurve
    virtual Go::Point eval( double t) const;

    // Inherited from EvalCurve
    virtual void eval(double t, int n, Go::Point der[]) const;

    // Inherited from EvalCurve
    virtual double start() const;

    // Inherited from EvalCurve
    virtual double end() const;

    /// Inherited from EvalCurve::dim().  For this class, the returned dimension will be that
    /// of the surface parameter domain, ie. 2, NOT that of the space curve.
    virtual int dim() const;

    /// Inherited from EvalCurve::approximationOK().  For this class, the specified tolerances
    /// are not used; the internally stored 'epsgeo' value is used as tolerance (this value was
    /// specified in the constructor).
    /// \param par the parameter at which to check the curve
    /// \param approxpos the position we want to check whether or not the curve
    ///                  approximates for parameter 'par'.
    /// \param tol1 unused
    /// \param tol2 unused
    /// \return 'true' if the curve approximates the point at the parameter, 'false'
    ///         otherwise.
    virtual bool approximationOK(double par, Go::Point approxpos,
				 double tol1, double tol2) const; 

private:
    const boost::shared_ptr<Go::SplineCurve> space_crv_;
    const boost::shared_ptr<Go::SplineSurface> surf_;
    const boost::shared_ptr<Go::Point> start_par_pt_; // When projecting end pts may be of special interest.
    const boost::shared_ptr<Go::Point> end_par_pt_;
    const double epsgeo_;
    const RectDomain* domain_of_interest_;

    // Simple function, interpolates and pts.
    std::vector<double> createSeed(double tpar) const;

};


///\}
}

#endif //_GOPROJECTCURVE_

