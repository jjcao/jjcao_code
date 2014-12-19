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
#ifndef _CONSTRAINTDEFINITIONS_H
#define _CONSTRAINTDEFINITIONS_H

#include <vector>


namespace Go
{
///\addtogroup geometry
///\{


    /// Struct defining linear side constraints between control points in
    /// a surface.
    /// The elements of factor_ correspond to a linear combination of
    /// control points on the left side of the equation, whilst
    /// constant_term_ denotes the right side of the equation.
    typedef struct sideConstraint
    {
	int dim_;  // Dimension of coefficient (max 3).
	// For each coefficient involved in the constraint, the 
	// index of the coefficient is given and the factor corresponding
	// to the coefficient in the equation.
	std::vector<std::pair<int, double> > factor_;
	double constant_term_[3];  // The constant term in the current equation.
    } sideConstraint;


    /// Struct defining linear side constraints between control points in
    /// a set of surfaces.
    /// The elements of factor_ correspond to a linear combination of
    /// control points on the left side of the equation, whilst
    /// constant_term_ denotes the right side of the equation.
    typedef struct sideConstraintSet
    {
	int dim_;  // Dimension of coefficient (max 3).
	// For each coefficient involved in the constraint, the index
	// of the surface, the index of the coefficient and the 
	// factor corresponding to the coefficient in the equation is given.
	std::vector<std::pair<std::pair<int,int>, double> > factor_; 
	// The constant term in the current equation, on the right side of the
	//  equation. For dim < 3 not all elements are used.
	double constant_term_[3];
	/// Default constructor.
	sideConstraintSet()
	{
	    dim_ = 3;
	    constant_term_[0] = constant_term_[1] = constant_term_[2] = 0.0;
	}

	/// Constructor.
	/// \param dim dimension of geometric space.
	sideConstraintSet(int dim)
	{
	    dim_ = dim;
	    constant_term_[0] = constant_term_[1] = constant_term_[2] = 0.0;
	}
    } sideConstraintSet;

///\}
}

#endif // _CONSTRAINTDEFINITIONS_H


