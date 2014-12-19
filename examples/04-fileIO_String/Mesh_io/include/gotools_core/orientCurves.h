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
#ifndef _ORIENTCURVES_H
#define _ORIENTCURVES_H

namespace Go
{
///\addtogroup geometry
///\{


//===========================================================================
/// This function sorts and orients a set of curves so that curves whose endpoints 
/// coincide will be ordered consecutively, and eventually 'reversed' so that 
/// startpoints meet endpoints.  It is assumed that the set of curves given 
/// constitute one or several manifolds, ie., that there will not be cases where
/// three or more endpoints meet at a common point (this would make the above
/// mentionned ordering impossible).
/// \param curves a vector containing the set of (pointers to) curves to be 
///        analysed.  It is expected that they constitute one or more 1-manifolds.
/// \param permutation upon return, this vector will contain a permutation of the indices
///        to the input curves such that after this permutation, curves whose endpoints 
///        are connected will become 'consecutive'.
/// \param reversed upon return, this vector, whose length will be equal to that of
///        'curves' and 'permutation', will contain bool values.  If curve at position
///        'i' in the 'curves' vector needs to be reversed in order to connect 
///        startpoint-to-endpoint with its neighbours, then the corresponding value in 
///        'reversed' will be 'true', and 'false' if no reversal is necessary.
/// \param neighbour_tol the tolerance used when checking for coincident points.
/// \param assume_manifold if the user specifies 'true', then the manifold property
///                        will be assumed, but will not be explicitly checked.  If 
///                        'false' is specified, then the input \em will be checked
///                        for consistency with respect to this, and an exception will
///                        be cast if the manifold condition is violated.
template <class PtrToCurveType>
inline void orientCurves(const std::vector<PtrToCurveType>& curves,
			 std::vector<int>& permutation,
			 std::vector<bool>& reversed,
			 double neighbour_tol,
			 bool assume_manifold = true)
//===========================================================================
{
    // registering all points
    int i, num_seg = curves.size();
    std::vector<Point> endpoints(2 * num_seg);
    for (i = 0; i < num_seg; ++i) {
	curves[i]->point(endpoints[2 * i + 0], curves[i]->startparam());
	curves[i]->point(endpoints[2 * i + 1], curves[i]->endparam());
    }

    // detecting connections
    std::vector<int> connected_to(2 * num_seg, -1);
    for (i = 0; i < 2 * num_seg; ++i) {
	if (connected_to[i] == -1 || !assume_manifold) {
	    Point& p1 = endpoints[i];
	    bool found = (connected_to[i] != -1);
	    // this point has not been checked for connections yet
	    for (int j = i+1; j < 2 * num_seg; ++j) {
		if (connected_to[j] == -1) {
		    if (p1.dist(endpoints[j]) < neighbour_tol) {
			// we consider these points to be connected
			connected_to[i] = j;
			connected_to[j] = i;
			if (assume_manifold) {
			    break;
			} else if (found == true) {
			    THROW("Multiple connected points detected in "
				  " orientCurves().  Not a 1-manifold!" );
			} else { 
			    found = true;
			}
		    }
		}
	    }
	}
    }

    // Making the permutation std::vector
    permutation.clear();
    permutation.reserve(num_seg);
    reversed.clear();
    reversed.reserve(num_seg);
    std::vector<bool> visited(2*num_seg, false);
    // First handle all open chains
    for (i = 0; i < 2*num_seg; ++i) {
	if (!visited[i] && connected_to[i] == -1) {
	    int current = i;
	    do {
		permutation.push_back(current/2);
		bool current_parity = (current%2 == 0);
		int partner = current_parity ? current + 1 : current - 1;
		reversed.push_back(!current_parity);
		visited[current] = true;
		visited[partner] = true;
		current = connected_to[partner];
	    } while (current != -1);
	}
    }
    // Then all closed chains (loops)
    for (i = 0; i < 2*num_seg; ++i) {
	if (!visited[i]) {
	    int current = i;
	    do {
		permutation.push_back(current/2);
		bool current_parity = (current%2 == 0);
		int partner = current_parity ? current + 1 : current - 1;
		reversed.push_back(!current_parity);
		visited[current] = true;
		visited[partner] = true;
		current = connected_to[partner];
	    } while (current != i);
	}
    }
}


///\}
} // namespace Go

#endif // _ORIENTCURVES_H


