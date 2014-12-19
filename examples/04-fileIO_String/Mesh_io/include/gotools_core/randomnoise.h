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
#ifndef _RANDOMNOISE_H
#define _RANDOMNOISE_H

namespace Go {
///\addtogroup utils
///\{


//===========================================================================
//                    FUNCTIONS FOR RANDOM DATA
//===========================================================================
//! Gives a certain number of random samples drawn from the normal distribution.
//! \param res a pointer to the array where the resulting samples should be written.
//! \param mean_err the sigma parameter to the normal distribution.
//! \param num_samples the desired number of samples (should also be the size of the
//! array pointed to by \em res.
void normalNoise(double* res, double mean_err, int num_samples);

//! Gives a certain number of random samples drawn from the uniform distribution.
//! \param res a pointer to the array where the resulting samples should be written.
//! \param lval lower bound of the range from which the samples can take their values.
//! \param uval upper bound of the range from which the samples can take their values.
//! \param num_samples the desired number of samples (should also be the size of the
//! array pointed to by \em res.
void uniformNoise(double* res, double lval, double uval, int num_samples);

///\}
}; // namespace Go


#endif // _RANDOMNOISE_H


