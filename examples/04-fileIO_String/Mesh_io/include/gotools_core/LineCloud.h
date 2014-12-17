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
#ifndef _GOLINECLOUD_H
#define _GOLINECLOUD_H


#include "GeomObject.h"
#include "Array.h"
#include <vector>


namespace Go
{
///\addtogroup geometry
///\{


    /** GeomObject representing a collection of line segments in space.
     *
     */

class LineCloud : public GeomObject
{
public:

    /// Makes an unitialized LineCloud that can later be assigned or read() into.
    LineCloud()
    {};

    /** start is supposed to point to the start of the array to be copied,
	its valuetype should be convertible to double */

    /// Generate a LineCloud based on values stored in memory.  The lines should
    /// be stored as a sequence of point pairs indicating the start and end position
    /// of each line in the line cloud.  Each point is stored as (x_coord, y_coord, z_coord).
    /// The coordinates should be convertible to 'double'.  Ex: p0_start_x, p0_start_y, 
    /// p0_start_z, p0_end_x, p0_end_y, p0_end_z, p1_start_x, p1_start_y, p1_start_z,....
    /// \param start pointer to the start of the memory area from which point information
    ///              is to be copied
    /// \param numlines number of lines in the LineCloud
    template <typename ForwardIterator>
    LineCloud(ForwardIterator start, int numlines)
	: points_(numlines*2)
    {
	std::copy(start, start + 6*numlines, points_[0].begin());
    }

    /// Virtual destructor, enables safe inheritance.
    virtual ~LineCloud();

    /// Read a LineCloud from an input stream
    /// \param is the stream from which we will read the LineCloud
    virtual void read (std::istream& is);

    /// Write the LineCloud to stream
    /// \param os the stream that we will write the LineCloud to
    virtual void write (std::ostream& os) const;

    /// Get the BoundingBox of the LineCloud
    /// \return the BoundingBox enclosing the LineCloud
    virtual BoundingBox boundingBox() const;

    /// Query the dimension of the space in which the LineCloud is embedded 
    /// (currently, only dimension=3 is allowed...)
    /// \return the dimension of the space (always 3 for now)
    virtual int dimension() const
    { return 3; }

    /// Get the ClassType of this GeomObject (which is of course the ClassType identified
    /// LineCloud).
    /// \return the ClassType identifier of LineCloud
    virtual ClassType instanceType() const
    { return LineCloud::classType(); }

    /// Get the ClassType identifier of LineCloud.
    /// \return the ClassType identifier of LineCloud
    static ClassType classType()
    { return Class_LineCloud; }

    /// Clone this object
    /// \return a pointer to a cloned LineCloud.
// #ifdef _MSC_VER
// #if _MSC_VER < 1300
//     virtual GeomObject* clone() const
//     { return new LineCloud(*this); }
// #else
//     virtual LineCloud* clone() const
//     { return new LineCloud(*this); }
// #endif // _MSC_VER < 1300
// #else
    virtual LineCloud* clone() const
    { return new LineCloud(*this); }
// #endif

    /// Fill a line cloud with information read from memory.  The layout of the
    /// read information should be as for the LineCloud constructor: LineCloud(ForwardIterator 
    /// start, int numlines)
    /// \param points pointer to the memory area where the coordinates of the elements
    ///               in the LineCloud can be found. (This information will be copied).
    /// \param numlines number of lines in the LineCloud.
    void setCloud(const double* points, int numlines);

    /// Query the number of lines in the LineCloud
    /// \return the LineCloud's number of lines.
    int numLines() const { return points_.size()/2; }

    /// Get a start or end point from a line in the LineCloud (non-const version)
    /// \param i the index of the start/end point.  If 'i' is pair, then the returned
    ///          point is a start point, else it is an end point.
    /// \return a reference to the requested start/end point
    Vector3D& point(int i) { return points_[i]; }

    /// Get a start or end point from a line in the LineCloud (const version)
    /// \param i the index of the start/end point.  If 'i' is pair, then the returned
    ///          point is a start point, else it is an end point.
    /// \return a const-reference to the requested start/end point
    const Vector3D& point(int i) const { return points_[i]; }

    /// Get a pointer to the start of the internal memory area where line information
    /// is stored.
    /// \return a pointer to the beginning of the array where line information is stored.
    double* rawData() { return points_[0].begin(); }

private:
    std::vector<Vector3D> points_;
};


///\}
} // namespace Go



#endif // _GOLINECLOUD_H


