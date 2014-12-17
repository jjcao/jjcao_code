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
#ifndef _GOVALUES_H_
#define _GOVALUES_H_

#ifndef MICROSOFT
#include <values.h>
#endif

#ifndef MAXINT
//const int MAXINT = 2147483647;
#define MAXINT 2147483647;
#endif

#ifndef MAXDOUBLE
#define MAXDOUBLE   1.79769313486231570e+308
#endif

#ifndef M_PI
const double M_PI = 3.14159265358979323846;
#endif

#ifndef DEFAULT_SPACE_EPSILON
#define DEFAULT_SPACE_EPSILON 1e-10
#endif 

#ifndef DEFAULT_PARAMETER_EPSILON
#define DEFAULT_PARAMETER_EPSILON 1e-12
#endif


#endif

/** @file Values.h
 * Defines the constants  MAXINT, MAXDOUBLE and M_PI if they are not
 * defined by system.
 */

