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
#ifndef _RATIONAL_H
#define _RATIONAL_H

#include <iostream>

namespace Go {
///\addtogroup utils
///\{

    
    /** Class representing rational numbers
     *
     */


class Rational
{
public:
    
    /// Construct a rational number whose value is 0.
    Rational() : p_(0), q_(1) {}

    /// Construct a rational number whose value is 'p' (integer)
    Rational(int p) : p_(p), q_(1) {}

    /// Construct a rational number whose value is 'p'/'q' (fraction)
    Rational(int p, int q) : p_(p), q_(q) {}

    /// Add a different rational number to this rational number
    Rational& operator += (const Rational& other)
    {
	p_ = p_*other.q_ + q_*other.p_;
	q_ = q_*other.q_;
	simplify();
	return *this;
    }

    /// Subtract a different rational number from this rational number
    Rational& operator -= (const Rational& other)
    {
	Rational tmp = -other;
	(*this) += tmp;
	return *this;
    }

    /// Multiply this rational number with a different rational number
    Rational& operator *= (const Rational& other)
    {
	p_ = p_*other.p_;
	q_ = q_*other.q_;
	simplify();
	return *this;
    }

    /// Divide this rational number by a different rational number
    Rational& operator /= (const Rational& other)
    {
	p_ = p_*other.q_;
	q_ = q_*other.p_;
	simplify();
	return *this;
    }

    /// Return the additive inverse of this rational number
    Rational operator- () const
    {
	return Rational(-p_, q_);
    }

    /// Test this rational number for equality with another rational number
    bool operator == (const Rational r)
    {
	return (p_ == r.p_ && q_ == r.q_);
    }

    /// Test this rational number for difference with another rational number
    bool operator != (const Rational r)
    {
	return (p_ != r.p_ || q_ != r.q_);
    }

    /// Write this rational number to a stream
    void write(std::ostream& os) const
    {
	os << p_ << '/' << q_;
    }

    /// Simplify the internal fractional expression of this rational number
    void simplify()
    {
	int n = std::min(abs(p_), abs(q_));
	for (int i = 2; i <= n; ++i) {
	    while (p_%i==0 && q_%i==0) {
		p_ /= i;
		q_ /= i;
		n /= i;
	    }
	}
    }

private:
    int p_;
    int q_;
};

Rational operator + (const Rational& r1, const Rational r2)
{
    Rational res = r1;
    res += r2;
    return res;
}

Rational operator - (const Rational& r1, const Rational r2)
{
    Rational res = r1;
    res -= r2;
    return res;
}

Rational operator * (const Rational& r1, const Rational r2)
{
    Rational res = r1;
    res *= r2;
    return res;
}

Rational operator / (const Rational& r1, const Rational r2)
{
    Rational res = r1;
    res /= r2;
    return res;
}

std::ostream& operator << (std::ostream& os, const Rational& p)
{
    p.write(os);
    return os;
}

///\}
}; // end namespace Go
#endif // _RATIONAL_H


