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
#ifndef _GOSISLCONVERSION_H
#define _GOSISLCONVERSION_H

//#include "sisl.h"
#include "SplineCurve.h"
#include "SplineSurface.h"

/// \file SISLconversion.h
/// Declaration file for a set of free conversion functions
/// between SISL and Spline curves and surfaces.

struct SISLCurve;
struct SISLSurf;

namespace Go
{
///\addtogroup geometry
///\{

/// Convert a SplineCurve to a SISLCurve
/// \param cv the SplineCurve to convert
/// \param copy if 'true', then the generated SISLCurve will have its own
///             copy of the coefficient information contained in 'cv'.
///             Otherwise, it will share this information with 'cv' (ie. 
///             it will only contain a pointer into the corresponding
///             storage array in 'cv'.
/// \return A newly generated SISLCurve that describes the same curve as 
///         'cv'.  The user assumes ownership and is responsible for 
///         cleaning up (which means calling the SISL function freeCurve(...)
///         on the pointer when it should be destroyed).
SISLCurve* Curve2SISL( const SplineCurve& cv, bool copy = true);

/// Convert a SISLCurve to a SplineCurve
/// \param cv the SISLcurve to convert
/// \return A newly generated SplineCurve that describes the same curve
///         as 'cv'.  The user assumes ownership and is responsible for 
///         cleaning up by calling the \c delete function.
SplineCurve* SISLCurve2Go( const SISLCurve* const cv);

/// Convert a SplineSurface to a SISLSurface
/// \param sf the SplineSurface to convert
/// \param copy if 'true', then the generated SISLSurf will have its own
///             copy of the coefficient information contained in 'sf'. 
///             Otherwise, it will share this information with 'sf' (ie.
///             it will only contain a pointer into the corresponding 
///             storage array in 'sf'.  
/// \return a newly generated SISLSurf that describes the same surface
///         as 'sf'.  The user assumes ownership and is responsible for 
///         cleaning up (which means calling the SISL function freeSurf(...)
///         on the pointer when it should be destroyed).    
SISLSurf* GoSurf2SISL( const SplineSurface& sf, bool copy = true);

/// Convert a SISLSurface to a SplineSurface
/// \param sf the SISLSurf to convert
/// \return A newly generated SplineSurface that describes the same surface
///         as 'sf'.  The user assumes ownership and is responsible for 
///         cleaning up by calling the \c delete function.
SplineSurface* SISLSurf2Go( SISLSurf* sf);

///\}
} // namespace Go

#endif // _GOSISLCONVERSION_H


