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
#ifndef _ROTATEDBOX_H
#define _ROTATEDBOX_H

#include "CompositeBox.h"
#include "MatrixXD.h"
#include <boost/shared_ptr.hpp>

namespace Go
{
///\addtogroup utils
///\{



    /** A rotated version of CompositeBox.
     *  It works in the same way, except that the boxes are
     *  aligned with an arbitrary (given) coordinate system.
     */

class RotatedBox
{
public:
    /// Given an array of dim-dimensional points and a
    /// coordinate system given by (axis[0], axis[1], axis[0] x axis[1]),
    /// construct a rotated box containing all points.
    /// If in 2D, coordinate system is (axis[0], Rot(Pi/2)*axis[0])
    template <typename RandomAccessIterator>
    RotatedBox(RandomAccessIterator start,
	       int dim,
	       int num_u,
	       int num_v,
	       const Point* axis)
    {
	// Make coordinate system.
	setCs(axis, dim);
	// Make box.
	setFromArray(start, dim, num_u, num_v);
    }

    /// Creates a RotatedBox with the Point low specifying
    /// the lower bound in all dimensions and high specifying
    /// the upper bound. The inner and edge boxes are equal.
    RotatedBox(const Point& low, const Point& high,
	       const Point* axis)
    {
	// Make coordinate system.
	setCs(axis, low.size());
	// Make box.
	setFromPoints(low, high);
    }

    /// Do not inherit from this class -- nonvirtual destructor.
    ~RotatedBox()
    {
    }

    /// Given an array of dim-dimensional points stored as doubles
    /// or floats, makes the smallest composite box containing all
    /// points in the array. The array must be like a control point
    /// grid, with num_u points in the fastest running direction,
    /// and num_v points in the other direction. For curves, use
    /// num_v = 1.
    template <typename RandomAccessIterator>
    void setFromArray(RandomAccessIterator start,
		      int dim,
		      int num_u,
		      int num_v)
    {
	// Make a temporary array of points
	std::vector<double> pts(dim*num_u*num_v);
	std::copy(start, start + dim*num_u*num_v, &pts[0]);
	// Transform the points.
	if (dim == 2) {
	    for (int i = 0; i < num_u*num_v; ++i) {
		Point p(&pts[0] + i*2,
			&pts[0] + (i+1)*2, false);
		p = cs2_*p;
		pts[i*2] = p[0];
		pts[i*2+1] = p[1];
	    }
	} else if (dim == 3) {
	    for (int i = 0; i < num_u*num_v; ++i) {
		Point p(&pts[0] + i*3,
			&pts[0] + (i+1)*3, false);
		p = cs3_*p;
		pts[i*3] = p[0];
		pts[i*3+1] = p[1];
		pts[i*3+2] = p[2];
	    }
	} else {
	    THROW("Only supports 2 and 3 dimensions.");
	}
	// Make the composite box.
	box_.reset(new CompositeBox(&pts[0], dim, num_u, num_v));
    }

    /// Makes the bounding box have lower bounds as specified in low
    /// and upper bounds as specified in high.
    void setFromPoints(const Point& low, const Point& high)
    {
	int dim = low.size();
	// Transform the points.
	if (dim == 2) {
	    Point p1 = cs2_*low;
	    Point p2 = cs2_*high;
	    // Make the composite box.
	    box_.reset(new CompositeBox(p1, p2));
	} else if (dim == 3) {
	    Point p1 = cs3_*low;
	    Point p2 = cs3_*high;
	    // Make the composite box.
	    box_.reset(new CompositeBox(p1, p2));
	} else {
	    THROW("Only supports 2 and 3 dimensions.");
	}
    }

    /// The dimension of the rotated box.
    int  dimension()  const
    {
	return box_->dimension();
    }

    /// The composite box. WARNING: The coordinates of this box must
    /// be interpreted in the coordinate system given by coordsystem().
    const CompositeBox& box() const
    {
	return *box_;
    }

    /// Returns true if the point pt is inside the
    /// box, up to tolerances. Tolerances may be specified
    /// separately for inner and edge boxes.
    bool containsPoint(const Point& pt,
		       double toli = 0.0,
		       double tole = 0.0) const
    {
	int dim = dimension();
	Point rotp;
	if (dim == 2) {
	    rotp = cs2_*pt;
	} else if (dim == 3) {
	    rotp = cs2_*pt;
	} else {
	    THROW("Only supports 2 and 3 dimensions.");
	}
	return box_->containsPoint(rotp, toli, tole);
    }

    /// Returns true if the two boxes overlap, up to 
    /// tolerances. Tolerances may be specified
    /// separately for inner and edge boxes.
    bool overlaps(const RotatedBox& box,
		  double toli = 0.0,
		  double tole = 0.0) const
    {
	MatrixXD<double, 2> m2 = cs2_;
	m2 += -box.cs2_;
	if (m2.frobeniusNorm() > 1e-12) {
	    THROW("The two RotatedBox objects have different coordinate systems.");
	}
	MatrixXD<double, 3> m3 = cs3_;
	m3 += -box.cs3_;
	if (m3.frobeniusNorm() > 1e-12) {
	    THROW("The two RotatedBox objects have different coordinate systems.");
	}
	return box_->overlaps(box.box(), toli, tole);
    }

    /// Returns true if this box contain the box passed 
    /// as a parameter,up to tolerances. Tolerances may be
    /// specified separately for inner and edge boxes.
    bool containsBox(const RotatedBox& box,
		     double toli = 0.0,
		     double tole = 0.0) const
    {
	MatrixXD<double, 2> m2 = cs2_;
	m2 += -box.cs2_;
	if (m2.frobeniusNorm() > 1e-12) {
	    THROW("The two RotatedBox objects have different coordinate systems.");
	}
	MatrixXD<double, 3> m3 = cs3_;
	m3 += -box.cs3_;
	if (m3.frobeniusNorm() > 1e-12) {
	    THROW("The two RotatedBox objects have different coordinate systems.");
	}
	return box_->containsBox(box.box(), toli, tole);
    }


private:
    void setCs(const Point* axis, int dim)
    {
	// What we actually compute is the inverse coordinate system.
	// This way when we multiply a vector with cs we get its
	// coordinates in the system given by (axis[0], Rot(Pi/2)*axis[0]).
	// Initiate the other array to zero to avoid an exception
	if (dim == 2) {
	    cs2_(0,0) = axis[0][0];
	    cs2_(0,1) = axis[0][1];
	    cs2_(1,0) = -axis[0][1];
	    cs2_(1,1) = axis[0][0];

	    cs3_(0,0) = cs3_(0,1) = cs3_(0,2) = 0.0;
	    cs3_(1,0) = cs3_(1,1) = cs3_(1,2) = 0.0;
	    cs3_(2,0) = cs3_(2,1) = cs3_(2,2) = 0.0;
	} else if (dim == 3) {
	    Point zaxis = axis[0] % axis[1];
	    zaxis.normalize();
	    cs3_(0,0) = axis[0][0];
	    cs3_(0,1) = axis[0][1];
	    cs3_(0,2) = axis[0][2];
	    cs3_(1,0) = axis[1][0];
	    cs3_(1,1) = axis[1][1];
	    cs3_(1,2) = axis[1][2];
	    cs3_(2,0) = zaxis[0];
	    cs3_(2,1) = zaxis[1];
	    cs3_(2,2) = zaxis[2];

	    cs2_(0,0) = cs2_(0,1) = 0.0;
	    cs2_(1,0) = cs2_(1,1) = 0.0;
	} else {
	    THROW("Only supports 2 and 3 dimensions.");
	}
    }

    // Data members
    boost::shared_ptr<CompositeBox> box_;
    MatrixXD<double, 2> cs2_;
    MatrixXD<double, 3> cs3_;
};

///\}
} // namespace Go



#endif // _ROTATEDBOX_H


