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
#ifndef _GOSTREAMABLE_H
#define _GOSTREAMABLE_H

#include <iostream>
#include <iomanip>
#include "errormacros.h"

namespace Go
{
///\addtogroup geometry
///\{


    /** 
     *  Base class for streamable objects, ie. objects which can be read from 
     *  and written to a stream.
     */

class Streamable
{
public:
    virtual ~Streamable();

    /// read object from stream
    /// \param is stream from which object is read
    virtual void read (std::istream& is) = 0;
    /// write object to stream
    /// \param os stream to which object is written
    virtual void write (std::ostream& os) const = 0;

    // Exception class
    class EofException{};
};

inline std::istream& operator >> (std::istream& is, Go::Streamable& obj)
{
    ALWAYS_ERROR_IF(is.eof(), "End of file reached. Cannot read.");
    obj.read(is);
    return is;
}

inline std::ostream& operator << (std::ostream& os,
				  const Go::Streamable& obj)
{
    obj.write(os);
    return os;
}

///\}
} // namespace Go



#endif // _GOSTREAMABLE_H


