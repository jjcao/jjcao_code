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
#ifndef _GOHERMITEAPPS_H_
#define _GOHERMITEAPPS_H_

#include "EvalCurve.h"
#include "SplineCurve.h"
#include "SplineSurface.h"
#include "HermiteGrid1DMulti.h"

namespace Go
{
///\addtogroup geometry
///\{


/// This class is used to generate a set of SplineCurves from a EvalCurveSet (which
/// itself represents a set of related curves) using Hermiet interpolation.  The generated
/// curves will approximate those defined by the EvalCurveSet within specified 
/// tolerances.  This class is really a generalization of HermiteAppC.

class HermiteAppS
{
public:
    /// Constructor where the tolerances and the curves to approximate are 
    /// specified.
    /// \param surface the curve set that we want to approximate by Hermite 
    ///                interpolation of sampled values.
    /// \param tolerance1 the required geometrical accuracy of approximation
    /// \param tolerance2 another tolerance, used for some kinds of EvalCurveSets.
    /// \param dims vector that specifies the spatial dimensions of the curves 
    ///             contained in the EvalCurveSet.  The size of the vector should
    ///             be equal to the total number of curves in 'surf', ie. the return
    ///             value of its EvalCurveSet::nmbCvs() function.
    HermiteAppS(EvalCurveSet* surface, 
		  double tolerance1, 
		  double tolerance2, 
		  std::vector<int> dims);

    /// Constructor where the tolerances and the curve set to approximate are 
    /// specified, as well as the parameters where they should be sampled prior to
    /// the Hermite interpolation.
    /// \param surface the curve set that we want to approximate by Hermite
    ///                interpolation of sampled values.
    /// \param initpars pointer to the array of parameter values for which we will
    ///                 sample the input curves.
    /// \param n number of parameter values in the array 'initpars[]'.
    /// \param tolerance1 the required geometrical accuracy of approximation
    /// \param tolerance2 another tolerance, used for some kinds of EvalCurveSets.
    /// \param dims vector that specifies the spatial dimensions of the curves 
    ///             contained in the EvalCurveSet.  The size of the vector should
    ///             be equal to the total number of curves in 'surf', ie. the return
    ///             value of its EvalCurveSet::nmbCvs() function.
    HermiteAppS(EvalCurveSet* surface, double initpars[], int n,
		  double tolerance1, double tolerance2, std::vector<int> dims);

    /// Empty destructor
    ~HermiteAppS(){}

    /// Refine the internal sampling of the curves such that the Hermite 
    /// interpolated curves parametrically approximates the original curve within
    /// a specified tolerance.
    void refineApproximation();	// Refine Hermite Grid.

    /// Return the cubic spline curves intepolating the grid (ie. approximating
    /// the original curve set). 
    /// \return a vector containing shared pointers to the newly created 
    ///         spline curves that Hermite interpolate the sampled points
    ///         of the EvalCurveSet (curve set) specified in the constructor.
    std::vector<boost::shared_ptr<SplineCurve> > getCurves();

private:
    EvalCurveSet* surface_;     // Pointer to original surface existing outside *this.
    HermiteGrid1DMulti grid_; // We define a grid for all curves, separating point values.
    const double tol1_;
    const double tol2_;
    const double min_interval_;	// Smaller intervals are not refined
    std::vector<boost::shared_ptr<SplineCurve> > curve_approx_; // Spline representations of approximation.
    /*   shared_ptr<SplineSurface> surface_approx_; // Spline representation of approximation */

    bool testSegment(int j, double& new_knot);	// Distance to _original
    int bisectSegment(int);


};


///\}
} // namespace Go

#endif // _GOHERMITEAPPS_H_

