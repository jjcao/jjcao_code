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
#ifndef _SPLINEDEBUGUTILS_H
#define _SPLINEDEBUGUTILS_H


#include "SplineSurface.h"
#include "SplineCurve.h"
#include "boost/smart_ptr.hpp"


namespace Go
{
///\addtogroup geometry
///\{


    /// For debugging. Writes a parameter curve in the xy-plane. Remove when
    /// GoViewer handles 2D curves.

    /// For debugging.  Writes a parameter curve (2D) in the xy-plane, for a given z-value.
    /// It will be written (with header) to the specified stream as a 3D curve.
    /// \param pcurve the parameter curve we want to write as a 3D curve
    /// \param os the stream to which we want to write the 3D curve
    /// \param z the constant z-value for the generated curve
    void writeSpaceParamCurve(const SplineCurve& pcurve, std::ostream& os, double z = 0.0);

    /// writes the geometric object (with header) to the specified file name.
    /// \param geom_obj the object to write to file.
    /// \param to_file the file name to which the object will be written.
    void objToFile(GeomObject* geom_obj, char *to_file);

    /// writes the geometric objects (with header) to the specified file name.
    /// \param geom_objs the objects to write to file.
    /// \param to_file the file name to which the objects will be written.
    void objsToFile(std::vector<boost::shared_ptr<GeomObject> >& geom_objs,
		    char *to_file);

    /// Write a SplineCurve to a stream using the SISL file format (not the Go format).
    /// \param spline_cv the curve to write to a stream
    /// \param os the stream to which the curve will be written (in SISL format)x
    void writeSISLFormat(const SplineCurve& spline_cv, std::ostream& os);

///\}
} // End of namespace Go


#endif // _SPLINEDEBUGUTILS_H


