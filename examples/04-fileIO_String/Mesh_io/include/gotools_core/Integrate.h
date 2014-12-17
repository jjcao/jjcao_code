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
#ifndef _GOINTEGRATE_H_
#define _GOINTEGRATE_H_

namespace Go
{
///\addtogroup geometry
///\{


class BsplineBasis;

  /// Functions used to compute integrals of inner products of B-splines

  /** Compute all definite integrals of inner products of
   *  derivatives of B-splines up to a given order where the
   *  differentiation is of the same order for both B-splines.
   *  The interval of integration are equal to the parameter
   *  intervals of the surface in the current par. dir.
   * @param basis B-spline basis.
   * @param ider Number of derivatives to compute.
   * @param lim1 Start of parameter interval.
   * @param lim2 End of parameter interval.
   * @param integral Computed integrals.
   */
    void GaussQuadInner(const BsplineBasis& basis, int ider, double lim1,
			double lim2, double*** integral);

  /** Compute all definite integrals of inner products of
   *  derivatives of B-splines up to a given order where the
   *  differentiation is of the same order for both B-splines.
   *  The interval of integration are equal to the parameter
   *  intervals of the surface in the current par. dir.
   * @param basis B-spline basis.
   * @param ider Number of derivatives to compute.
   * @param lim1 Start of parameter interval.
   * @param lim2 End of parameter interval.
   * @param integral Computed integrals.
   */
    void GaussQuadInner2(const BsplineBasis& basis, int ider, double lim1,
			 double lim2, double** integral);


///\}
};
#endif

