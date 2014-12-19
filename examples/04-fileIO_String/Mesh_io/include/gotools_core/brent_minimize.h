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
#ifndef _BRENT_MINIMIZE_H
#define _BRENT_MINIMIZE_H

#include "GeneralFunctionMinimizer.h"

namespace Go
{
///\addtogroup utils
///\{

    /** Brief description. 
     *  Detailed description.
     */

template<class Functor>
class Fun2Fun
{
public:
    Fun2Fun(const Functor& f, double a, double b) : f_(f), a_(a), b_(b) {}
    double operator()(const double* arg) const
    {
	return f_(*arg);
    }
/// void grad(const double* arg, double* grad) const;
    double minPar(int n) const { return a_; }
    double maxPar(int n) const { return b_; }
private:
    Functor f_;
    double a_;
    double b_;
};

//===========================================================================
template<class Functor>
inline double brent_minimize(const Functor& f,
			     double a, double b, double c,
			     double& parmin,
			     const double rel_tolerance = std::sqrt(std::numeric_limits<double>::epsilon()))
//===========================================================================
{
    Fun2Fun<Functor> f2(f, a, c);
    Go::FunctionMinimizer<Fun2Fun<Functor> > fmin(1, f2, &a, rel_tolerance);
    Go::Point dir(1);
    dir[0] = 1.0;
    double bracket[3];
    double fval_brak[3];
    bracket[0] = 0.0;
    bracket[1] = b-a;
    bracket[2] = c-a;
    fval_brak[0] = f(a);
    fval_brak[1] = f(b);
    fval_brak[2] = f(c);
    double minimum = fmin.linminBrent(dir, bracket, fval_brak);
    parmin = fmin.getPar(0);
    return minimum;
}


///\}
} // namespace Go

#endif // _BRENT_MINIMIZE_H


