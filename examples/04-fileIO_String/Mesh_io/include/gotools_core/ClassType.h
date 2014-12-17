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
#ifndef _GOCLASSTYPE_H
#define _GOCLASSTYPE_H

namespace Go
{
///\addtogroup geometry
///\{


    /** All concrete classes that inherit GeomObject should
     *  have a matching enum in ClassType. It is necessary to
     *  have a central location for the type identification numbers,
     *  so that two classes in different modules (for example) does
     *  not try to use the same number. At the same time, it is
     *  something we'd like to avoid, because it forces you to
     *  update this file every time you add a new GeomObject-
     *  inheriting class.
     */
    enum ClassType
    {
	Class_Unknown = 0,

	Class_SplineCurve = 100,
	Class_CurveOnSurface = 110,

	Class_SplineSurface = 200,
	Class_BoundedSurface = 210,
	Class_GoBaryPolSurface = 220,
	Class_GoHBSplineParamSurface = 230,

	Class_Go3dsObject = 300,
	Class_GoHeTriang = 310,
	Class_GoSdTriang = 320,
	Class_GoQuadMesh = 330,
	Class_GoHybridMesh = 340,
	Class_ParamTriang = 350,
	Class_GoVrmlGeometry = 360,
	Class_PointCloud = 400,
	Class_LineCloud = 410,
	Class_GoTriangleSets = 500,
	Class_RectGrid = 510
    };

///\}
}


#endif // _GOCLASSTYPE_H


