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
#ifndef _COORDINATESYSTEM_H
#define _COORDINATESYSTEM_H

#include "Array.h"
#include "MatrixXD.h"

namespace Go
{
///\addtogroup utils
///\{


    /// Defines a Cartesian coordinate system

template <int Dim>
class CoordinateSystem
{
public:
    typedef MatrixXD<double, Dim> Matrix;
    typedef Array<double, Dim> Vector;

    /// The default constructor creates an object with 
    /// an identity rotation matrix and zero translation vector.
    CoordinateSystem()
    {
	rotation_.identity();
	translation_.zero();
    }

    /// Creates a CoordinateSystem with the specified 
    /// rotation matrix and translation vector.
    CoordinateSystem(const Matrix& rot, const Vector& trans)
	: rotation_(rot), translation_(trans)
    {}

    /// Get the 'global' coordinates of a vector 'v' expressed
    /// in the CoordinateSystem.
    Vector operator* (const Vector& v)
    {
	return rotation_*v + translation_;
    }

    /// Create a CoordinateSystem that is a composition of 'this' 
    /// CoordinateSystem and 'c'.
    CoordinateSystem operator* (const CoordinateSystem& c)
    {
	Matrix newrot = rotation_*c.rotation_;
	Vector newtrans = translation_ + rotation_*c.translation_;
	return CoordinateSystem(newrot, newtrans);
    }

    /// Get the rotation matrix.
    Matrix& rot() { return rotation_; }

    /// Get the rotation matrix. Const version.
    const Matrix& rot() const { return rotation_; }

    /// Get the translation vector.
    Vector& tr() { return translation_; }

    /// Get the translation vector. Const version.
    const Vector& tr() const { return translation_; }

private:
    Matrix rotation_;
    Vector translation_;
};


///\}
} // namespace Go



#endif // _COORDINATESYSTEM_H


