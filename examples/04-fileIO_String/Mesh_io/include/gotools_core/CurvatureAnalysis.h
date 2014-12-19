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
#ifndef _CURVATUREANALYSIS_H
#define _CURVATUREANALYSIS_H


#include "ParamSurface.h"

/// \file CurvatureAnalysis.h
///Functions computing fundamental forms and curvature.

namespace Go
{
///\addtogroup geometry
///\{

    /// Computes the coefficients of the first fundamental form.
    /// The returned values are stored 
    /// { E F G [Eu Fu Gu Ev Fv Gv
    ///   [Euu Fuu Guu Euv Fuv Guv Evv Fvv Gvv [...] ] ] }.
    /// Only one derivative is implemented so far.
    /// \param sf reference to the concerned surface
    /// \param u first parameter of point where we want to compute the first 
    ///          fundamental form
    /// \param v second parameter of point where we want to compute the first
    ///          fundamental form
    /// \param derivs number of (partial) derivatives of the first fundamental
    ///               form that we want to include in the result.  So far, 
    ///               only the computation of the first derivative is actually 
    ///               implemented.
    /// \param form the computed values, on the format described above
    void computeFirstFundamentalForm(const ParamSurface& sf,
				     double u, double v, int derivs,
				     std::vector<double>& form);

    /// Computes the coefficients of the first and second fundamental
    /// forms.
    /// The returned values are stored 
    /// { E F G } in form1 and { e f g } in form2.
    /// This function cannot compute derivatives.
    /// \param sf reference to the concerned surface
    /// \param u first parameter of point where we want to compute the fundamental
    ///          forms.
    /// \param v second parameter of point where we want to compute the fundamental
    ///          forms
    /// \param form1 the values associated with the first fundamental form will be 
    ///              returned here.
    /// \param form2 the values associated with the second fundamental form will be
    ///              returned here.
    void computeSecondFundamentalForm(const ParamSurface& sf,
				      double u, double v,
				      double form1[3],
				      double form2[3]);

    /// Computes the Gaussian (K) and mean (H) curvatures.
    /// \param sf reference to the concerned surface
    /// \param u first parameter of point where we want to carry out computation
    /// \param v second parameter of point where weh want to carry out computation
    /// \param K value of Gaussian curvature returned here
    /// \param H value of mean curvature returned here
    void curvatures(const ParamSurface& sf,
		    double u, double v,
		    double& K, double& H);
///\}
} // namespace Go

#endif // _CURVATUREANALYSIS_H


