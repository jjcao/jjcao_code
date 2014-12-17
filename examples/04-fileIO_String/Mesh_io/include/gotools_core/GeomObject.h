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
#ifndef _GOGEOMOBJECT_H
#define _GOGEOMOBJECT_H

#include "BoundingBox.h"
#include "Streamable.h"
#include "ClassType.h"

namespace Go
{
///\addtogroup geometry
///\{


const int MAJOR_VERSION = 1;
const int MINOR_VERSION = 0;

    /** 
     *  Base class for geometrical objects (curves, surfaces, etc.) regrouping
     *  all the properties that they have in common.
     */

class GeomObject : public Streamable
{
public:
    virtual ~GeomObject();
    
    /// Return the object's bounding box
    virtual BoundingBox boundingBox() const = 0;
    
    /// Return the dimension of the space in which the object lies (usually 2 or 3)
    virtual int dimension() const = 0;

    /// Return the class type identifier of a given, derived instance of GeomObject
    virtual ClassType instanceType() const = 0;

    /// Return the class type identifier of a given class derived from GeomObject
    static ClassType classType();

    /// Clone the GeomObject and return a pointer to the clone.
    virtual GeomObject* clone() const = 0;

    /// Write header information of the GeomObject to stream.  This typically precedes
    /// the act of writing the object itself to a stream, to signal to the receiver
    /// what object is streamed.
    void writeStandardHeader(std::ostream& os) const;
};

///\}
} // namespace Go

#endif // _GOGEOMOBJECT_H





