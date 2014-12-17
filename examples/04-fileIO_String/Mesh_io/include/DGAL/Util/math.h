#ifndef DGAL_MATH__H__
#define DGAL_MATH__H__

#include <DGAL/config.h>
#include <numeric>
DGAL_BEGIN_NAMESPACE

//typedef CGAL::Cartesian<Num_type_default> 
/// Return cotangent of (P,Q,R) corner (ie cotan of QP,QR angle).
template<
	class Kernel = Kernel_default
>
class Math{
public:
	typedef typename Kernel::Point_2 Point_2;
	typedef typename Kernel::Point_3 Point_3;
	typedef typename Kernel::Vector_3 Vector_3;

	static double cotangent(const Point_3& P,const Point_3& Q,const Point_3& R)
	{
		Vector_3 u = P - Q;
		Vector_3 v = R - Q;
		// (u . v)/((u x v).len)
		double dot = (u*v);
		Vector_3 cross_vector = CGAL::cross_product(u,v);
		double cross_norm = std::sqrt(cross_vector*cross_vector);
		if(cross_norm != 0.0)
			return (dot/cross_norm);
		else
			return 0.0; // undefined
	}
	static double cotangent(const Point_2& P,const Point_2& Q, const Point_2& R)
	{
		Point_3 p(P.x(), P.y(), 0.0);
		Point_3 q(Q.x(), Q.y(), 0.0);
		Point_3 r(R.x(), R.y(), 0.0);
	
		return cotangent(p,q,r);
	}
	
	//                                                  -> ->
	/// Return tangent of (P,Q,R) corner (ie tangent of QP,QR angle).
	static double tangent(const Point_3& P,const Point_3& Q,const Point_3& R)
	{
		Vector_3 u = P - Q;
		Vector_3 v = R - Q;
		// (u . v)/((u x v).len)
		double dot = (u*v);
		CGAL_surface_mesh_parameterization_assertion(dot != 0.0);
		Vector_3 cross_vector = CGAL::cross_product(u,v);
		double cross_norm = std::sqrt(cross_vector*cross_vector);
		if(dot != 0.0)
			return (cross_norm/dot);
		else
			return 0.0; // undefined
	}
	
	//                                                     -> ->
	/// Return angle (in radians) of of (P,Q,R) corner (ie QP,QR angle).
	static double compute_angle_rad(const Point_3& P,
									const Point_3& Q,
									const Point_3& R)
	{
		static const double PI = 3.14159265359;
	
		Vector_3 u = P - Q;
		Vector_3 v = R - Q;
	
		// check
		double product = std::sqrt(u*u) * std::sqrt(v*v);
		if(product == 0)
			return 0.0;
	
		// cosine
		double dot = (u*v);
		double cosine = dot / product;
	
		// sine
		Vector_3 w = CGAL::cross_product(u,v);
		double AbsSine = std::sqrt(w*w) / product;
	
		if(cosine >= 0)
			return std::asin(fix_sine(AbsSine));
		else
			return PI-std::asin(fix_sine(AbsSine));
	}
	/// Fix sine.
	static double fix_sine(double sine)
	{
		if(sine >= 1)
			return 1;
		else if(sine <= -1)
			return -1;
		else
			return sine;
	}
};


DGAL_END_NAMESPACE
#endif//DGAL_MATH__H__