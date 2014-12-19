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
#ifndef _GOOBJECTHEADER_H
#define _GOOBJECTHEADER_H


#include <vector>
#include "ClassType.h"
#include "Streamable.h"

namespace Go
{
///\addtogroup geometry
///\{


    /** An object representing the "header" usually preceeding a GeomObject 
     *  in a stream.  This header contains information about the GeomObject,
     *  and this information can be read by ObjectHeader and accessed by its 
     *  member functions.
     */ 


class ObjectHeader : public Streamable
{
public:
    /// Default constructor (uninitialized header)
    ObjectHeader()
	: class_type_(Class_Unknown),
	  major_version_(0),
	  minor_version_(0)
    {}

    /// Constructor creating an ObjectHeader that is initialized
    /// with a given ClassType, major and minor version.
    /// \param t the ClassType of the GeomObject that this ObjectHeader
    ///          shall represent.
    /// \param major major version number
    /// \param minor minor version number
    ObjectHeader(ClassType t, int major, int minor)
	: class_type_(t),
	  major_version_(major),
	  minor_version_(minor)
    {}

    /// Constructor creating an ObjectHeader that is initialized with
    /// a given ClassType, major and minor version and auxiliary, class-
    /// specific data.
    /// \param t the ClassType of the GeomObject that this ObjectHeader
    ///          shall represent.
    /// \param major major version number
    /// \param minor minor version number
    /// \param auxdata auxiliary data, specific to the ClassType.
    ObjectHeader(ClassType t, int major, int minor,
		   const std::vector<int>& auxdata)
	: class_type_(t),
	  major_version_(major),
	  minor_version_(minor),
	  auxillary_data_(auxdata)
    {}

    /// Virtual destructor, allowing safe destruction of derived objects.
    virtual ~ObjectHeader();

    /// Read the ObjectHeader from an input stream
    /// \param is the input stream from which the ObjectHeader is read
    virtual void read (std::istream& is);

    /// Write the ObjectHeader to an output stream
    /// \param os the output stream to which the ObjectHeader is written
    virtual void write (std::ostream& os) const;

    /// Get the ClassType stored in this ObjectHeader
    ClassType classType() { return class_type_; }

    /// Get the major version number stored in this ObjectHeader
    int majorVersion() { return major_version_; }

    /// Get the minor version number stored in this ObjectHeader
    int minorVersion() { return minor_version_; }

    /// Get the size of the auxiliary data stored in this ObjectHeader 
    /// (size measured in number of ints).
    int auxdataSize() { return auxillary_data_.size(); }

    /// Get a certain piece of auxiliary data (an integer).
    /// \param i the requested integer's position in the auxiliary data vector
    int auxdata(int i) { return auxillary_data_[i]; }

private:
    ClassType class_type_;
    int major_version_;
    int minor_version_;
    std::vector<int> auxillary_data_;
};


///\}
} // namespace Go


#endif // _GOOBJECTHEADER_H


