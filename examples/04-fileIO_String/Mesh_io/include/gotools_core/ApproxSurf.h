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
#ifndef _GOAPPROXSURF_H_
#define _GOAPPROXSURF_H_

//   -----------------------------------------------------------------------
//      Interface file for class ApproxSurf
//   -----------------------------------------------------------------------
//
//       Approximate a set of points by a B-spline curve to
//       satisfy a given accuracy
//
//       Implementation of the member functions are given in the
//       following files:
//
//          1. ApproxSurf.C
//
//   -----------------------------------------------------------------------
//    Written by: Vibeke Skytt                           04-00
//   -----------------------------------------------------------------------

#include "SplineSurface.h"
#include "SplineCurve.h"
#include <vector>



class SmoothSurf;

namespace Go
{
///\addtogroup geometry
///\{

/// This class can generate a B-spline surface that approximates
/// a set of points for a given accuracy.
class ApproxSurf
{
 public:
    /// Constructor where the user specifies the boundary curves of the surface
    /// to generate, a parameter domain for the surface, the points to approximate
    /// and their parameter values, as well as the geometric tolerance.  The two 
    /// spline basises of the generated tensor product spline surface will be
    /// determined by unifying basises of opposing boundary curves.
    /// \param crvs the boundary curves of the surface to be generated.  This vector
    ///             should contain exactly \em four curves, whose endpoints are 
    ///             connected so that they form a loop.
    /// \param points vector containing the coordinates of the points that this 
    ///               surface should interpolate.  They are stored in
    ///               "xyzxyz...-fashion".
    /// \param parvals vector containing the parameter values of the points given in
    ///                the 'points' vector.  They are stored in "uvuv...-fashion".
    /// \param domain pointer to an array of four doubles specifying the parametric
    ///               domain for the surface to be generated.  They should be stored
    ///               as "u_min, u_max, v_min, v_max".
    /// \param dim spatial dimension of the points (usually 3).
    /// \param aepsge geometric tolerance to use internally
    /// \param constdir The points \em will be reparameterized internally according
    ///                 to their spatial position with respect to the surface that 
    ///                 shall be generated.  However, they might be reparameterized
    ///                 in both their u and v parameters, only in their u parameters
    ///                 or only in their v parameters.  The user can specify this 
    ///                 with 'constdir'.  If 'constdir' is set to 0, the points will
    ///                 be reparameterized in  both u and v.  If 'constdir' is set to
    ///                 1, they will only be reparameterized in the v parameter. 
    ///                 If 'constdir' is set to 2, they will only be reparameterized in
    ///                 the u parameter.
    ApproxSurf(std::vector<boost::shared_ptr<SplineCurve> > & crvs,
		 const std::vector<double>& points, 
		 const std::vector<double>& parvals,
		 double domain[],
		 int dim, 
		 double aepsge,
		 int constdir = 0);

    /// Constructor where the user specifies a spline surface that should be
    /// modified, the points to approximate and their parameter values, as well
    /// as the geometric tolerance.  The surface that is given as argument is not 
    /// copied internally, only pointed to, so it \em will be modified.
    /// \param srf the surface that will be modified to approximate the points.
    ///            Assumed to contain k-regular knots.
    /// \param points vector containing the coordinates of the points that this 
    ///               surface should interpolate.  They are stored in
    ///               "xyzxyz...-fashion".
    /// \param parvals vector containing the parameter values of the points given in
    ///                the 'points' vector.  They are stored in "uvuv...-fashion".
    /// \param dim spatial dimension of the points (usually 3).
    /// \param aepsge geometric tolerance to use internally
    /// \param constdir The points \em will be reparameterized internally according
    ///                 to their spatial position with respect to the surface that 
    ///                 shall be generated.  However, they might be reparameterized
    ///                 in both their u and v parameters, only in their u parameters
    ///                 or only in their v parameters.  The user can specify this 
    ///                 with 'constdir'.  If 'constdir' is set to 0, the points will
    ///                 be reparameterized in  both u and v.  If 'constdir' is set to
    ///                 1, they will only be reparameterized in the v parameter. 
    ///                 If 'constdir' is set to 2, they will only be reparameterized in
    ///                 the u parameter.
    ApproxSurf(boost::shared_ptr<SplineSurface>& srf,
		 const std::vector<double>& points, 
		 const std::vector<double>& parvals,
		 int dim, double aepsge, int constdir = 0);


    /// Destructor
    ~ApproxSurf();

    /// Sets the smoothing weight to something other than the default (1e-9).
    /// The value should lie in the unit interval, typically close to 0.
    /// \param smooth the new smoothing weight.
    void setSmoothingWeight(double smooth)
	{
	    ASSERT(smoothweight_ >= 0.0 && smoothweight_ <= 1.0);
	    smoothweight_ = smooth;
	}

    /// Decide whether or not the total bondary of the surface should be kept fixed
    /// (i.e. unchanged by approximation process).   Default is true.  Cross derivatives
    /// will not be kept fixed.  (If you want to keep cross derivatives fixed, use
    /// the edgeFix() member function instead).
    /// \param fix_boundary if 'true' the boundary of the surface will not be modified,
    ///                     if 'false' it will be open to modification.
    void setFixBoundary(bool fix_boundary)
	{
	    int fix = (fix_boundary) ? 1 : 0;
	    edge_derivs_[0] = edge_derivs_[1] = edge_derivs_[2] = edge_derivs_[3] = fix;
	}

    /// Decide whether specific edges of the surface's boundary should be kept fixed
    /// (i.e. unchanged by approximation process), as well as a certain number of cross-
    /// derivatives across these curves.
    /// \param edge_fix pointer to an array of four integers specifying to which extent
    ///                 each surface edge should be kept fixed during the approximation 
    ///                 process.
    ///                 A value of 0 means that it will not be kept fixed, 1 means that
    ///                 its position will be kept fixed, 2 means that its position and 
    ///                 cross-tangent will be kept fixed, etc. 
    ///                 The integers are associated with the surface edges starting 
    ///                 with the edge 'v=vmin' and moving counterclockwise.
    void edgeFix(int edge_fix[])  // CCV
	{
	    for (int ki=0; ki<4; ki++)
		edge_derivs_[ki] = edge_fix[ki];
	}

    /// Forces the surface to approximate certain normals at
    /// certain parameter values.
    /// \param points  this vector contains the normals that should be approximated.
    ///                They are stored in 'xyzxyz... fashion'.
    /// \param parvals this vector contains the parameter values of the normals that
    ///                should be approximated.  They are stored in 'uvuv... fashion'.
    void setNormalConditions(const std::vector<double>& points, 
			     const std::vector<double>& parvals)
	{
	    use_normals_ = true;
	    norm_points_ = points;
	    norm_parvals_ = parvals;
	}

    /// Fetch the approximating surface

    /// When everything else is set, this function can be used to run the 
    /// approximation process and fetch the approximated surface.
    /// \retval maxdist report the maximum distance between the approximated 
    ///                 surface and the data points
    /// \retval avdist report the average distance between the approximated 
    ///                surface and the datapoints
    /// \param nmb_out_eps report the number of data points that were found to
    ///                    lie outside the geometric tolerance (as specified by 
    ///                    the 'aepsge' argument to the ApproxSurf constructor.
    /// \param max_iter maximum number of iterations
    /// \return a shared pointer to the generated SplineCurve, approximating 
    ///         the points as specified.
    boost::shared_ptr<SplineSurface> getApproxSurf(double& maxdist, 
						   double& avdist,
						   int& nmb_out_eps, 
						   int max_iter=4);

    /// Reparameterize the data points by a closest point
    /// match against the current surface.
    int reParam();


 protected:
    /// Default constructor
    ApproxSurf();

 private:
    boost::shared_ptr<SplineSurface> curr_srf_;
    double maxdist_;
    double avdist_;
    int outsideeps_;
    double aepsge_;
    double smoothweight_;
    double smoothfac_;
    bool use_normals_;
    int edge_derivs_[4];

    int dim_;
    std::vector<double> points_;
    std::vector<double> parvals_;
    std::vector<double> norm_points_;
    std::vector<double> norm_parvals_;
    int constdir_;

    /// Generate an initial curve representing the spline space
    int makeInitSurf(std::vector<boost::shared_ptr<SplineCurve> > &crvs, 
		     double domain[]);

    /// Generate a smoothing surface
    int makeSmoothSurf();

    /// Check the accuracy of the current surface
    int checkAccuracy();

    /// Generate the approximating surface
    int doApprox(int max_iter);


};

///\}
}  // namespace Go

#endif







