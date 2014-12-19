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
#ifndef _GOHERMITEAPPC_H_
#define _GOHERMITEAPPC_H_

#include "EvalCurve.h"
#include "SplineCurve.h"
#include "HermiteGrid1D.h"


namespace Go
{
///\addtogroup geometry
///\{


/// This class is used to generate a SplineCurve from a EvalCurve
/// using Hermite interpolation.  The generated curve will approximate the
/// EvalCurve within specified tolerances.

class HermiteAppC
{
public:

    /// Constructor where the tolerances and the curve to approximate are specified.
    /// \param crv the curve that we want to generate a Hermite approximation of.
    ///            The curve is \em not copied, only pointed to by the HermiteAppC.
    /// \param tolerance1 the required geometrical accuracy of approximation
    /// \param tolerance2 another tolerance, used for some kinds of EvalCurves
    HermiteAppC(EvalCurve* crv, double tolerance1, double tolerance2);

    /// Constructor where the tolerances and the curve to approximate are specified,
    /// as well as the parameters for which we will sample the input curve before
    /// starting the approximating process.
    /// \param crv the curve that we want to generate a Hermite approximation of.
    ///            The curve is \em not copied, only pointed to by the HermiteAppC.
    /// \param initpars pointer to the array of parameter values for which we will
    ///        sample the input curve.
    /// \param n number of parameter values in the array 'initpars'.
    /// \param tolerance1 the required geometrical accuracy of approximation
    /// \param tolerance2 another tolerance, used for some kinds of EvalCurves
    HermiteAppC(EvalCurve* crv, double initpars[],int n, 
		  double tolerance1, double tolerance2);

    /// Empty destructor
    ~HermiteAppC(){}

    /// Refine the internal sampling of the curve to approximate such that
    /// the Hermite interpolated curve of this sampling approximates parametrically
    /// the original curve within the specified tolerance.
    void refineApproximation();	// Refine Hermite Grid.

    /// Return the cubic spline curve Hermite interpolating the grid.
    /// \return the cubic spline curve that Hermite interpolates the sampled 
    ///         points of the EvalCurve specified in the constructor.
    boost::shared_ptr<SplineCurve> getCurve();	

 private:
    EvalCurve* curve_;	// Pointer to original curve existing outside *this.
    const double tol1_;
    const double tol2_;
    const double min_interval_;	// Smaller intervals are not refined
    HermiteGrid1D grid_;
    boost::shared_ptr<SplineCurve> curve_approx_; // Spline representation of approximation

    bool testSegment(int j, double& new_knot);	// Distance to _original
    int bisectSegment(int);



};

///\}
} // namespace Go;

#endif

