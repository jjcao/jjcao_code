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
#ifndef _CLOSESTPTCURVES_H
#define _CLOSESTPTCURVES_H

/** @file closestPtCurves.h
 * Declaration file for an intersection function operating on 
 * objects of class ParamCurve belonging to geometry.
 */

#include "ParamCurve.h"

namespace Go {
///\addtogroup geometry
///\{


  /** Compute a closest point or an intersection point between two curves.
   * Ported from the sisl-functions s1770, s1770_s9corr, s1770_s9dir
   * \param cv1   Curve number one.
   * \param cv2   Curve number two.
   * \retval par1 Parameter value of the first curve's point.
   * \retval par2 Parameter value of the second curve's point.
   * \retval dist Distance between the points.
   * \retval ptc1 Point on curve number one.
   * \retval ptc2 Point on curve number two.
   */
  void closestPtCurves(const ParamCurve* cv1, const ParamCurve* cv2,
		       double& par1, double& par2, double& dist,
		       Point& ptc1, Point& ptc2);

  /*
  void closestPtCurves(ParamCurve* cv1, ParamCurve*cv2, double tmin1,
		       double tmax1, double tmin2, double tmax2,
		       double& par1, double& par2, double& dist,
		       Point& ptc1, Point& ptc2)
  */

  /** Newton iteration on the distance function between two curves to
   * find a closest point or an intersection point.
   * The solution is restricted to the intervals [tmin1,tmax1] and 
   * [tmin2,tmax2] in the first and second curves's parameter domain.
   * Ported from the sisl function s1770.
   * \param cv1 Curve number one.
   * \param cv2 Curve number two.
   * \param tmin1 Min. parameter value. First curve  
   * \param tmax1 Max. parameter value. First curve
   * \param tmin2 Min. parameter value. Second curve
   * \param tmax2 Max. parameter value. Second curve
   * \param seed1 Start point for iteration along curve number one.
   * \param seed2 Start point for iteration along curve number two.
   * \retval par1 Parameter value of the first curve's point.
   * \retval par2 Parameter value of the second curve's point.
   * \retval dist Distance between the points.
   * \retval ptc1 Point on curve number one.
   * \retval ptc2 Point on curve number two.
   */
  void closestPtCurves(const ParamCurve* cv1, const ParamCurve* cv2, double tmin1,
		       double tmax1, double tmin2, double tmax2,
		       double seed1, double seed2, double& par1, double& par2,
		       double& dist, Point& ptc1, Point& ptc2);

 /** Computes initial start points for iteration along the curves.
   * \param pc1 Curve number one.
   * \param pc2 Curve number two.
   * \retval seed1 Start point for iteration along curve number one.
   * \retval seed2 Start point for iteration along curve number two.
   */
  void computeSeedCvCv(const SplineCurve* pc1, const SplineCurve* pc2,
		       double& seed1, double& seed2);


  /** Adjust delta to satisfy: \f[ astart \leq acoef+delta \leq aend \f]
   * Ported from the sisl function s1770_s9corr.
   */
  void insideParamDomain(double& delta, double acoef, double astart,
			 double aend);


  /** Computes the distance vector and value beetween
   * a point on the first curve and a point on the second
   * curve. And computes a next step on both curves.
   * This is equivalent to the nearest way to the
   * parameter plane in the tangent plane from a point in the
   * distance surface between two curves.
   * Ported from the sisl function s1770_s9dir.
   * METHOD : The method is to compute the parameter distance to the points
   * on both tangents which is closest to each other.
   * The difference vector beetween these points are orthogonal
   * to both tangents. If the distance vector beetween the two
   * points on the curve is "diff" and the two derivative vectors
   * are "der1" and "der2", and the two wanted parameter distances
   * are "dt1" and "dt2", then we get the following system of 
   * equations:
   * @verbatim
      <dt1*der1+dist-dt2*der2,der2> = 0
      <dt1*der1+dist-dt2*der2,der1> = 0
   	       This is further:
    
    | -<der1,der2>   <der2,der2> |  | dt1 |   | <diff,der2> |
    |                            |  |     | = |             |
    | -<der1,der1>   <der1,der2> |  | dt2 |   | <diff,der1> |
   @endverbatim
   * 
   * The solution of this matrix equation dt1,dt2 are returned in the
   * parameters cdiff1,cdiff2.
   *
   */
  void nextStep(double& cdist, double& cdiff1, double& cdiff2,
		std::vector<Point>& eval1, std::vector<Point>& eval2);

///\}
} // namespace Go

#endif // _CLOSESTPTCURVES_H

