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
#ifndef _GOLIFTCURVE_
#define _GOLIFTCURVE_



#include "Point.h"
#include "EvalCurve.h"
#include "SplineCurve.h"
#include "SplineSurface.h"

#include <boost/smart_ptr.hpp>

namespace Go
{
///\addtogroup geometry
///\{


/// This class represents a "lift curve", generated from taking
/// a surface and a 2D curve and then evaluation the surface 
/// at the parameter values obtained by evaluating the 2D curve.

class LiftCurve : public EvalCurve
{
public:

  /// Constructor, taking a 2D parameter curve and a surface.
  /// \param parameter_crv the 2D parameter curve that will be 'lifted'
  /// \param surf the surface on which the resulting, 'lifted' curve will
  ///             lie
  /// \param epsgeo geometrical tolerance used when running the 'approximationOK'
  ///               function.
  LiftCurve(boost::shared_ptr<Go::SplineCurve>& parameter_crv,
	      boost::shared_ptr<Go::SplineSurface>& surf,
	      double epsgeo);

  /// virtual destructor enables safe inheritance
  virtual ~LiftCurve();

  // Inherited from EvalCurve
  virtual Point eval( double t) const;

  // Inherited from EvalCurve
  virtual void eval(double t, int n, Point der[]) const;

  // Inherited from EvalCurve
  virtual double start() const;

  // Inherited from EvalCurve
  virtual double end() const;

  /// Dimension of the lifted curve (i.e. 3).
  virtual int dim() const;

  /// Inherited from EvalCurve::approximationOK().  
  /// \param par the parameter at which to check the curve
  /// \param approxpos the position we want to check whether or not the curve
  ///                  approximates for parameter 'par'.
  /// \param tol1 unused
  /// \param tol2 unused
  /// \return 'true' if the curve approximates the point at the parameter
  ///         (within the tolerance given in the constructor, 'epsgeo'). 'false'
  ///         otherwise.
  virtual bool approximationOK(double par, Point approxpos,
			       double tol1, double tol2) const;

 private:
  const boost::shared_ptr<Go::SplineCurve> parameter_crv_;
  const boost::shared_ptr<Go::SplineSurface> surf_;
  const double epsgeo_;

};


///\}
} // namespace Go

#endif //_GOLIFTCURVE_

