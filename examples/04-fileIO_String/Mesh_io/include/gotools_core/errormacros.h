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
#ifndef _ERRORMACROS_H
#define _ERRORMACROS_H

#include <iostream>

/// Usage: REPORT;
/// Usage: MESSAGE("Message string.");
#ifdef NVERBOSE // Not verbose mode
#  ifndef REPORT
#    define REPORT
#  endif
#  ifndef MESSAGE
#    define MESSAGE(x)
#  endif
#  ifndef MESSAGE_IF
#    define MESSAGE_IF(cond, m)
#  endif
#else // Verbose mode
#  ifndef REPORT
#    define REPORT std::cerr << "\nIn file " << __FILE__ << ", line " << __LINE__ << std::endl
#  endif
#  ifndef MESSAGE
#    define MESSAGE(x) std::cerr << "\nIn file " << __FILE__ << ", line " << __LINE__ << ": " << x << std::endl
#  endif
#  ifndef MESSAGE_IF
#    define MESSAGE_IF(cond, m) do {if(cond) MESSAGE(m);} while(0)
#  endif
#endif

/// Usage: THROW("Error message string.");
#ifndef THROW
#  define THROW(x) MESSAGE(x), throw std::exception()
#endif

#define ALWAYS_ERROR_IF(condition, message) do {if(condition){ THROW(message);}} while(0)

/// Usage: ASSERT(condition)
/// Usage: ASSERT2(condition, "Error message string.")
/// Usage: DEBUG_ERROR_IF(condition, "Error message string.");
#ifdef NDEBUG // Not in debug mode
#  ifndef ASSERT
#    define ASSERT(x)
#  endif
#  ifndef ASSERT2
#    define ASSERT2(cond, x)
#  endif
#  ifndef DEBUG_ERROR_IF
#    define DEBUG_ERROR_IF(cond, x)
#  endif
#else // Debug mode
#  ifndef ASSERT
#    define ASSERT(cond) if (!(cond)) THROW("Assertion \'" #cond "\' failed.")
#  endif
#  ifndef ASSERT2
#    define ASSERT2(cond, x) do { if (!(cond)) THROW(x);} while(0)
#  endif
#  ifndef DEBUG_ERROR_IF
//#    define DEBUG_ERROR_IF(cond, x) if (cond) THROW(x) 
#    define DEBUG_ERROR_IF(cond, x) do { if (cond) THROW(x); } while(0)
#  endif
#endif


#endif // _ERRORMACROS_H





