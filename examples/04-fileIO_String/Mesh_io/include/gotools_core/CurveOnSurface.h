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
#ifndef _GOCURVEONSURFACE_H
#define _GOCURVEONSURFACE_H

#include "ParamSurface.h"
#include "ParamCurve.h"
#include <boost/smart_ptr.hpp>


namespace Go
{
///\addtogroup geometry
///\{


class SplineCurve;

    /** Missing doxygen documentation class description
     *
     */

class CurveOnSurface : public ParamCurve
{
public:
    /// Define an empty CurveOnSurface that can be assigned or \ref read()
    /// into
    CurveOnSurface();

    /// Construct a CurveOnSurface by specifying the surface and either a space curve
    /// or a curve in the parameter plane.  Does not clone any of the input, just 
    /// sets (smart) pointers.
    /// \param surf pointer to the underlying surface
    /// \param curve pointer to the curve specifying the CurveOnSurface.  This
    ///              curve may either be in the parametric domain of the surface
    ///              (2D curve), or a space curve that the user has assured 
    ///              to be coincident with the surface.  'preferparameter' specifies
    ///              which kind of curve this is.
    /// \param preferparameter if this is set to 'true', then 'curve' is assumed
    ///                        to be a curve in the parametric domain of the surface.
    ///                        Otherwise, it is assumed to be a space (3D) curve.
    CurveOnSurface(boost::shared_ptr<ParamSurface> surf,
		     boost::shared_ptr<ParamCurve> curve,
                     bool preferparameter);
    
    /// Construct a CurveOnSurface by specifying the surface, the space curve and the
    /// curve in the parameter plane.  The arguments are checked for consistency, and if
    /// they appear incoherent, an exception will be thrown.
    /// Does not clone any of the input, just sets (smart) pointers.
    /// \param surf pointer to the underlying surface
    /// \param pcurve pointer to the curve representing the CurveOnSurface in the
    ///               parametric domain of the surface.
    /// \param spacecurve pointer to the curve that is the spatial (3D) representation 
    ///                   of the CurveOnSurface.
    /// \param preferparameter specify whether the parametric curve or the space
    ///                        curve are preferred for internal computations.
    CurveOnSurface(boost::shared_ptr<ParamSurface> surf,
		     boost::shared_ptr<ParamCurve> pcurve,
                     boost::shared_ptr<ParamCurve> spacecurve,
                     bool preferparameter);

    /// Copy constructor.  The copy constructor will not clone() the underlying
    /// surface, but it will clone() both the parametric and the spatial curve.
    /// \param surface_curve the CurveOnSurface to copy into 'this' CurveOnSurface.
    CurveOnSurface(const CurveOnSurface& surface_curve);

    /// Assignment operator.  Like the copy constructor, the assignment operator
    /// clone()s the curves, and not the surface.
    /// \param other the CurveOnSurface to copy into 'this' CurveOnSurface.
    CurveOnSurface& operator= (const CurveOnSurface& other);
    
    /// Destructor.
    /// Trivial because memory is managed by boost::shared_ptr.
    virtual ~CurveOnSurface();


    // inherited from Streamable
    virtual void read (std::istream& is);
    virtual void write (std::ostream& os) const;

    // inherited from GeomObject
    // @afr: This one only returns the bounding box of the underlying surface
    virtual BoundingBox boundingBox() const;
    virtual DirectionCone directionCone() const;
    virtual int dimension() const;
    virtual ClassType instanceType() const;
    static ClassType classType()
    { return Class_CurveOnSurface; }
    
    // The function clone() calls the copy constructor,
    // so clone() also makes deep copies of the curves,
    // but not the underlying surface.
// #ifdef _MSC_VER
// #if _MSC_VER < 1300
//     virtual GeomObject* clone() const
//       { return new CurveOnSurface(*this); }
// #else
//     virtual CurveOnSurface* clone() const
//       { return new CurveOnSurface(*this); }
// #endif // _MSC_VER < 1300
// #else
    virtual CurveOnSurface* clone() const
    { return new CurveOnSurface(*this); }
// #endif

    // inherited from ParamCurve
    virtual void point(Point& pt, double tpar) const;

    /// Inherited from \ref ParamCurve. Only works for 'derivs' = 0 or 1.
    /// \see ParamCurve::point()
    virtual void point(std::vector<Point>& pts, 
		       double tpar,
		       int derivs, bool from_right = true) const;
    
    // inherited from ParamCurve
    virtual double startparam() const;

    // inherited from ParamCurve
    virtual double endparam() const;

    // inherited from ParamCurve
    virtual void reverseParameterDirection(bool switchparam = false);

    // inherited from ParamCurve
    virtual SplineCurve* geometryCurve();

    // inherited from ParamCurve
    virtual bool isDegenerate(double degenerate_epsilon);


// #if ((_MSC_VER > 0) && (_MSC_VER < 1300))
//     // inherited from ParamCurve
//     virtual ParamCurve* subCurve(double from_par, double to_par,
// 				   double fuzzy = DEFAULT_PARAMETER_EPSILON) const;
// #else
    // inherited from ParamCurve
    virtual CurveOnSurface* subCurve(double from_par, double to_par,
				       double fuzzy =
				       DEFAULT_PARAMETER_EPSILON) const;
// #endif

    // inherited from ParamCurve
    virtual void closestPoint(const Point& pt,
			      double         tmin,
			      double         tmax,
			      double&        clo_t,
			      Point&       clo_pt,
			      double&        clo_dist,
			      double const   *seed = 0) const;

    // inherited from ParamCurve.  NB: does not check whether the resulting ParamCurve
    // stays inside parameter domain (or that the space curve stays on surface).
    virtual void appendCurve(ParamCurve* cv, bool reparam=true);

    // inherited from ParamCurve.  NB: does not check whether the resulting ParamCurve
    // stays inside parameter domain (or that the space curve stays on surface).
    virtual void appendCurve(ParamCurve* cv,
			     int continuity, double& dist, bool reparam=true);

    /// Set the underlying surface to the one pointed to by the argument
    /// \param surface the pointer to the surface we will set as underlying for this
    ///                CurveOnSurface.
    void setUnderlyingSurface(boost::shared_ptr<ParamSurface> surface)
    {surface_ = surface;}

    /// Inherited from \ref ParamCurve.  If the parametric curve is set to be the
    /// 'prefered' one, this function will return the next segment value for the 
    /// parametric curve; otherwise it will return the next segment value for the 
    /// spatial 3D curve.
    /// See also \ref ParamCurve::nextSegmentVal()
    virtual double nextSegmentVal(double par, bool forward, double tol) const;


    /// Get a shared pointer to the underlying surface
    /// \return a shared pointer to the underlying surface
    boost::shared_ptr<ParamSurface> underlyingSurface()
    { return surface_; }

    /// Get a shared pointer to the curve in the parameter domain.
    /// \return a shared pointer to the curve in the parameter domain
    boost::shared_ptr<ParamCurve> parameterCurve()
    { return pcurve_; }

    /// Get a shared pointer to the space curve
    /// \return a shared pointer to the space curve.
    boost::shared_ptr<ParamCurve> spaceCurve()
    { return spacecurve_; }

    /// Get a constant, shared pointer to the underlying surface
    /// \return a const-pointer to the underlying surface
    boost::shared_ptr<const ParamSurface> underlyingSurface() const
    { return surface_; }

    /// Get a constant, shared pointer to the curve in the parameter domain.
    /// \return a const-pointer to the curve in the parameter domain.
    boost::shared_ptr<const ParamCurve> parameterCurve() const
    { return pcurve_; }

    /// Get a constant, shared pointer to the space curve
    /// \return a const-shared pointer to the space curve.
    boost::shared_ptr<const ParamCurve> spaceCurve() const
    { return spacecurve_; }
    
    /// Query whether the parameter curve or the space curve is prefered for computation
    /// in this object.
    bool parPref() const
    { return prefer_parameter_; }

    /// Get the rectangle enclosing the underlying surface's parametric domain.
    /// \return the RectDomain for the underlying surface.
    RectDomain containingDomain() const;

private:
    /// The underlying surface
    boost::shared_ptr<ParamSurface> surface_;
    /// The 2D curve in the parameter domain of the surface.
          // May point to null.
    boost::shared_ptr<ParamCurve> pcurve_;
    // An instance of the curve in the surface. May point to null.
    boost::shared_ptr<ParamCurve> spacecurve_;
    // Which representation to prefer if both exist
    bool prefer_parameter_;
};


///\}
} // namespace Go

#endif // _GOCURVEONSURFACE_H


