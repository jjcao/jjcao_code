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
#ifndef _GOCROSSTANOFFDIST_
#define _GOCROSSTANOFFDIST_

#include "EvalCurve.h"
#include "SplineCurve.h"
#include <vector>
#include <boost/smart_ptr.hpp>


namespace Go
{
///\addtogroup geometry
///\{


/// This class defines an evaluator-based offset-curve.
/// We're blending between a set of curves seen as an evaluator
/// based curve.
class CrossTanOffDist : public EvalCurve
{
public:

    /// Constructor
    /// \param poscurve the curve to offset from.
    /// \param tangcv1 tangent curve along poscurve.
    /// \param tangcv2 cross tangent curve along poscurve.
    /// \param blend1 1-dimensional blending function for tangcv1.
    /// \param blend2 1-dimensional blending function for tangcv2.
    /// \param opposite1 space curve corresponding to poscurve,
    ///                  parametrized in the opposite direction.
    /// \param opposite2 space curve corresponding to poscurve,
    ///                  parametrized in the opposite direction.
    ///                  May be equal to opposite1.
    /// \param factor if != 1.0 the poscurve, opposite1 & opposite2
    ///               will be used to define the opposite curve, so
    ///               as to minimize differences in parametrization
    ///               along the cv.
    CrossTanOffDist(boost::shared_ptr<SplineCurve>& poscurve,
		    boost::shared_ptr<SplineCurve>& tangcv1,
		    boost::shared_ptr<SplineCurve>& tangcv2,
		    boost::shared_ptr<SplineCurve>& blend1,
		    boost::shared_ptr<SplineCurve>& blend2,
		    boost::shared_ptr<SplineCurve>& opposite1,
		    boost::shared_ptr<SplineCurve>& opposite2,
		    double factor);


    /// Destructor.
    virtual ~CrossTanOffDist();

    /// Evaluate a point on the curve for a given parameter
    /// \param t the parameter for which to evaluate the curve.
    /// \return the evaluated point
    virtual Point eval( double t) const;

    /// Evaluate a point and a certain number of derivatives 
    /// on the curve for a given parameter.
    /// \param t the parameter for which to evaluate the curve.
    /// \param n the number of derivatives (0 or more)
    /// \retval der pointer to an array of Points where the 
    ///         result will be written.  The position will be stored
    ///         first, then the first derivative (tangent), then the
    ///         second, etc..
    ///         \b NB: For most (all) derived classes of 'EvalCurve', 
    ///         the implementation actually only supports the computation of 
    ///         one derivative, i.e. if n > 1, only one derivative will be 
    ///         computed anyway.
    virtual void eval( double t, int n, Point der[]) const; // n = order of diff

    /// Get the start parameter of the curve.
    /// \return the start parameter of the curve.
    virtual double start() const;

    /// Get the end parameter of the curve.
    /// \return  the end parameter of the curve.
    virtual double end() const;

    /// Get the dimension of the space in which the curve lies.
    /// \return the space dimension of the curve.
    virtual int dim() const;

    /// Check if the curve, evaluated at a given parameter, approximates
    /// a given position within a given tolerance.
    /// \param par the parameter at which to check the curve
    /// \param approxpos the position we want to check whether or not the curve
    ///                  approximates for parameter 'par'.
    /// \param tol1 approximation tolerance,
    /// \param tol2 another approximation tolerance (its use is defined by some of
    ///             the derived classes.
    /// \return 'true' if the curve approximates the point at the parameter, 'false'
    ///         otherwise.
    virtual bool approximationOK(double par, Point approxpos,
				 double tol1, double tol2) const;

private:
    const boost::shared_ptr<SplineCurve> poscurve_;
    std::vector<boost::shared_ptr<SplineCurve> > tangcurves_;
    std::vector<boost::shared_ptr<SplineCurve> > blends_;
    boost::shared_ptr<SplineCurve> oppositepos_;
    boost::shared_ptr<SplineCurve> lengthfac_;
    boost::shared_ptr<SplineCurve> avcross_;
    const double epstol_;

    /// Evaluate up to n derivatives of th offset tangent in the input parameter.
    /// \param t parameter in which to evaluate.
    /// \param derblend the derivatives in the offset direction.
    /// \param derproj the projection of the difference vector between the
    ///                derivatives in the offset direction when using the blending
    ///                functions and when using a linear blend between end offset
    ///                tangents.
    void evalcrtan(double t, Point& blend, Point& projdiff) const ;

    /// Evaluate up to n derivatives of th offset tangent in the input parameter.
    /// \param t parameter in which to evaluate.
    /// \param n the number of derivatives to compute.
    /// \param derblend array of size 'n + 1' containing the derivatives in the
    ///                 offset direction.
    /// \param derproj array of the same size, containing the projection of the
    ///                difference vector between the derivatives in the offset
    ///                direction when using the blending functions and when using
    ///                a linear blend between end offset tangents.
    void evalcrtan(double t, int n, Point derblend[], Point derproj[]) const;

    /// Evaluate the offset point using the blending functions to create a linear
    /// combination of the input tangent curves.
    /// \param t the parameter in which to evaluate.
    /// \return the blended offset point.
    Point evalblend(double t) const;

    /// Evaluate the difference between poscurve_ & oppositepos_ in input parameter.
    /// \param t parameter in which to evaluate.
    /// \return the difference vector between poscurve_ & oppositepos_ in t.
    Point evaldiff(double t) const ;

    /// Evaluate the difference between poscurve_ & oppositepos_ in input parameter,
    /// with the given number of derivatives.
    /// \param t parameter in which to evaluate.
    /// \param n the number of derivatives to compute.
    /// \param der the difference vector between poscurve_ & oppositepos_ in t, up
    ///            to the n'th derivative. Size of vector is 'n + 1'.
    void evaldiff(double t, int n, Point der[]) const ;
};


///\}
} // namespace Go

#endif

