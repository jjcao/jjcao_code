
#ifndef DGAL_Matrix_color_h_
#define DGAL_Matrix_color_h_

#include <DGAL/config.h>
#include <iostream>
#include <algorithm>
using std::istream;
using std::ostream;
using std::max;
using std::min;
#include <math.h>

DGAL_BEGIN_NAMESPACE

template<class Number_type>
class Color
{
public:
	void set_color(Number_type currValue, Number_type totalValue)
	{
		Number_type ncv = currValue/totalValue;
		
		if ( ncv < 0.5)//blue->green
		{
			Number_type t = 2 * ncv;
			r = 0;
			g = t;
			b = 1-t;
		}
		else//green->red
		{
			Number_type t = 2 * (ncv-0.5);
			r = t;
			g = 1-t;
			b = 0;
		}		
	}
	/*void fromXYZ(Number_type x, Number_type y, Number_type z);
    void toXYZ  (Number_type& x, Number_type& y, Number_type& z);

    void fromYIQ(Number_type q, Number_type i, Number_type y);
    void toYIQ  (Number_type& q, Number_type& i, Number_type& y);

    void fromHSV(Number_type h, Number_type s, Number_type v);
    void toHSV  (Number_type& h, Number_type& s, Number_type& v);*/
public:
	void set(Number_type R=0.0, Number_type G=0.0, Number_type B=0.0, Number_type A=1.0){
		r = R; g = G; b = B; a = A;
	}
	template<class InType>
	void to_array(InType* in){
		in[0] = InType(r);
		in[1] = InType(g);
		in[2] = InType(b);
		in[3] = InType(a);
	}
	Color(Number_type R=0, Number_type G=0, Number_type B=0, Number_type A=1) : r(R),g(G),b(B),a(A){} 
	Color& operator=(const Color& in) { r=in.r ; g=in.g; b=in.b ; a=in.a; return *this ; }
public:
	Number_type r,g,b, a;	
};
  
typedef Color<double> ColorF;
 
 // /*!
 //   \brief transforms from the XYZ color space to the RGB colorspace

 //   \author Philippe Lavoie
 //   \date 17 March 1999
 // */
 // inline void Color::fromXYZ(double x, double y, double z){
 //   r = (unsigned char)(255.0*(3.240479*x-1.537510*y-0.498535*z));
 //   g = (unsigned char)(255.0*(-0.969256*x+1.875992*y+0.041556*z));
 //   b = (unsigned char)(255.0*(0.055648*x-0.204043*y+1.057311*z));
 // }

 // inline void Color::toXYZ(double &x, double& y, double& z){
 //   
 // }

 // /*!
 //   \brief transforms from the YIQ color space to the RGB colorspace

 //   This is the same color space used by NTSC

 //   \param y the luminicance
 //   \param i the chromacity
 //   \param q the chromacity

 //   \author Philippe Lavoie
 //   \date 14 May 1999
 // */
 // inline void Color::fromYIQ(double y, double i, double q){
 //   r = (unsigned char)(255.0*(1.0030893*y+0.954849*i+0.6178597*q));
 //   g = (unsigned char)(255.0*(0.9967760*y-0.27070623*i-0.64478833*q));
 //   b = (unsigned char)(255.0*(1.0084978*y-1.11048518*i+1.69956753125));
 // }

 // /*!
 //   \brief transforms to the YIQ color space to the RGB colorspace

 //   This is the same color space used by NTSC

 //   \param y the luminicance
 //   \param i the chromacity
 //   \param q the chromacity

 //   \author Philippe Lavoie
 //   \date 14 May 1999
 // */
 // inline void Color::toYIQ(double &y, double &i, double &q){
 //   double R= double(R)/255.0 ;
 //   double G= double(R)/255.0 ;
 //   double B= double(R)/255.0 ;
 //   y = 0.299*R + 0.587*G + 0.114*B ;
 //   i = 0.596*R - 0.275*G - 0.321*B ;
 //   q = 0.212*R - 0.528*G + 0.311*B ;
 // }

 // /*!
 //   \brief from the HSV color space
 //   
 //   \param h hue valid inside [0,360]
 //   \param s saturation valid inside [0,1]
 //   \param v value valid inside [0,1]

 //   \author Philippe Lavoie
 //   \date 14 May 1999
 //  */
 // inline void Color::fromHSV(double h, double s, double v){
 //   if(s==0.0){
 //     r=g=b=0;
 //     return;
 //   }
 //   if(h>=360)
 //     h = 0.0 ;
 //   if(h<=0.0)
 //     h = 0.0 ;
 //   h /= 60.0 ;
 //   int i = int(floor(h)) ;
 //   double f = h-double(i);
 //   double p = v*(1.0-s);
 //   double q = v*(1-(s*f));
 //   double t = v*(1-(s*(1-f)));
 //   switch(i){
 //   case 0: 
 //     r = (unsigned char)(255.0*v) ; 
 //     g = (unsigned char)(255.0*t) ; 
 //     b = (unsigned char)(255.0*p) ; break ; 
 //   case 1: 
 //     r = (unsigned char)(255.0*q) ; 
 //     g = (unsigned char)(255.0*v) ; 
 //     b = (unsigned char)(255.0*p) ; break ; 
 //   case 2: 
 //     r = (unsigned char)(255.0*p) ; 
 //     g = (unsigned char)(255.0*v) ; 
 //     b = (unsigned char)(255.0*t) ; break ; 
 //   case 3: 
 //     r = (unsigned char)(255.0*p) ; 
 //     g = (unsigned char)(255.0*q) ; 
 //     b = (unsigned char)(255.0*v) ; break ; 
 //   case 4: 
 //     r = (unsigned char)(255.0*t) ; 
 //     g = (unsigned char)(255.0*p) ; 
 //     b = (unsigned char)(255.0*v) ; break ; 
 //   case 5: 
 //   default:
 //     r = (unsigned char)(255.0*v) ; 
 //     g = (unsigned char)(255.0*p) ; 
 //     b = (unsigned char)(255.0*q) ; break ; 
 //   }
 // }

 // /*!
 //   \brief to the HSV color space

 //   \param h hue valid inside [0,360]
 //   \param s saturation valid inside [0,1]
 //   \param v value valid inside [0,1]

 //   \author Philippe Lavoie
 //   \date 14 May 1999
 //  */
 // inline void Color::toHSV(double &h, double &s, double &v){
 //   double R = double(r)/255.0;
 //   double G = double(g)/255.0;
 //   double B = double(b)/255.0;
 //   
	//double max1 = max(R,max(G,B));
 //   double min1 = min(R,min(G,B));
 //   
 //   int maxI = max(r,max(g,b));

 //   v = max1 ; 
 //   s = (max1>0) ? (max1-min1)/max1 : 0 ; 
 //   h = 0 ;
 //   if(s>0){
 //     double delta = max1-min1 ;
 //     if(r==maxI){
	//h = (G-B)/delta ;
 //     }
 //     else{
	//if(g==maxI){
	//  h = 2.0 + (B-R)/delta ;
	//}
	//else{
	//  h = 4.0 + (R-G)/delta ;
	//}
 //     }
 //     h *= 60.0 ;
 //     if(h<0)
	//h += 360.0 ; 
 //   }
 // }

DGAL_END_NAMESPACE

#endif//DGAL_Matrix_color_h_
