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
#ifndef _RECTGRID_H
#define _RECTGRID_H

#include "GeomObject.h"

namespace Go
{
///\addtogroup geometry
///\{


/// This class represent the simplest of all quadrangulations: a rectangular
/// grid.
class RectGrid : public GeomObject
{
public:

    /// Constructor that defines an empty grid, which can only be assigned or 
    /// \ref read() into.
    RectGrid()
	: numu_(0), numv_(0), dim_(-1)
    {
    }
    
    /// Constructor explicitly specifying a RectGrid
    /// \param numu number of points along the first grid direction
    /// \param numv number of points along the second grid direction
    /// \param dim dimension of the points (usually 2 or 3)
    /// \param pts pointer to the array where the points are stored.
    ///            The points should be stored so that the first grid 
    ///            has lowest stride (ie., coefficients are stored 
    ///            (u1v1, u2v1,... unv1, u1v2, u2v2, ...unv2,.... unvm)
    RectGrid(int numu, int numv, int dim, double* pts)
	: numu_(numu), numv_(numv), dim_(dim),
	  points_(pts, pts + numu*numv*dim)
    {
    }
    /// virtual destructor assures safe inheritacne
    virtual ~RectGrid();

    // Inherited from GeomObject
    virtual BoundingBox boundingBox() const;

    // Inherited from GeomObject
    virtual int dimension() const;

    // Inherited from GeomObject
    virtual ClassType instanceType() const;

    // Inherited from GeomObject
    static ClassType classType()
    {
	return Class_RectGrid;
    }

    // Inherited from GeomObject
    virtual GeomObject* clone() const
    {
	return new RectGrid(*this);
    }

    // Inherited from Streamable
    void read(std::istream& is);

    // Inherited from Streamable
    void write(std::ostream& os) const;

    /// Set grid according to the user-supplied data
    /// \param numu number of points along the first grid direction
    /// \param numv number of points along the second grid direction
    /// \param dim dimension of the points (usually 2 or 3)
    /// \param pts pointer to the array where the points are stored.
    ///            The points should be stored so that the first grid 
    ///            has lowest stride (ie., coefficients are stored 
    ///            (u1v1, u2v1,... unv1, u1v2, u2v2, ...unv2,.... unvm)
    void setGrid(int numu, int numv, int dim, const double* pts)
    {
	std::vector<double> tmp(pts, pts + numu*numv*dim);
	points_.swap(tmp);
	numu_ = numu;
	numv_ = numv;
	dim_ = dim;
    }

    /// Get number of points along the first grid direction
    /// \return the number of points along the first grid direction
    int numCoefs_u() const
    {
	return numu_;
    }

    /// Get number of points along the second grid direction
    /// \return the number of points along the second grid direction
    int numCoefs_v() const
    {
	return numv_;
    }

    /// Get a pointer to the start of the array where grid point coordinates
    /// are stored.
    /// \return a pointer to the storage array for grid points.
    double* rawData()
    {
	return &points_[0];
    }

    /// Get a const pointer to the start of the array where grid point coordinates
    /// are stored.
    /// \return a const pointer to the storage array for grid points.
    const double* rawData() const
    {
	return &points_[0];
    }

private:
    int numu_;
    int numv_;
    int dim_;
    std::vector<double> points_;
};

///\}
} // namespace Go


#endif // _RECTGRID_H


