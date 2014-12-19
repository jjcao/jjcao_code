#ifndef DGAL_BORDER_PARIMETERIZE__H__
#define DGAL_BORDER_PARIMETERIZE__H__

#include <DGAL/config.h>
#include <DGAL/Util/math.h>

enum BORDER_STRATEGY {
	CIRCULAR_ARC  = 0,
	SQUARE_ARC,
	MEASURED
};

DGAL_BEGIN_NAMESPACE


template<
	class Polyhedron_3
>
class Border_parameterizer
{
public:
	typedef typename Polyhedron_3::NT            NT;
	typedef typename Polyhedron_3::Traits Kernel;
	typedef typename DGAL::Math<Kernel> Math;
		
	typedef typename Polyhedron_3::Point_3 Point_3;
	typedef typename Polyhedron_3::Vector_3 Vector_3;	
	typedef typename Polyhedron_3::Point_2       Point_2;
	typedef typename Polyhedron_3::Vector_2      Vector_2;
	
	typedef typename Polyhedron_3::Facet         Facet;
	typedef typename Polyhedron_3::Facet_handle  Facet_handle;
	typedef typename Polyhedron_3::Facet_const_handle	Facet_const_handle;
	typedef typename Polyhedron_3::Facet_iterator Facet_iterator;
	typedef typename Polyhedron_3::Facet_const_iterator	Facet_const_iterator;
	
	typedef typename Polyhedron_3::Vertex        Vertex;
	typedef typename Polyhedron_3::Vertex_handle Vertex_handle;
	typedef typename Polyhedron_3::Vertex_const_handle	Vertex_const_handle;
	typedef typename Polyhedron_3::Vertex_iterator Vertex_iterator;
	typedef typename Polyhedron_3::Vertex_const_iterator Vertex_const_iterator;
	
	typedef typename Polyhedron_3::Halfedge_around_facet_circulator
		Vertex_around_facet_circulator;
	typedef typename Polyhedron_3::Halfedge_around_facet_const_circulator
		Vertex_around_facet_const_circulator;
	typedef typename Polyhedron_3::Halfedge_around_vertex_circulator
		Vertex_around_vertex_circulator;
	typedef typename Polyhedron_3::Halfedge_around_vertex_const_circulator
		Vertex_around_vertex_const_circulator;		

	typedef typename std::list<Vertex_handle>::iterator Border_vertex_iterator;
		
	virtual ~Border_parameterizer(){}
	virtual bool parameterize_border(Border_vertex_iterator first, Border_vertex_iterator last, 
			Vertex_handle center_vertex)=0;
};

template<
	class Polyhedron_3
>
class Border_polar_map : public Border_parameterizer<Polyhedron_3>
{
public:
	bool parameterize_border(Border_vertex_iterator first, Border_vertex_iterator last, Vertex_handle center_vertex)
	{
		// Nothing to do if no border
		if (first == last)
			return false;

		Point_3 p0 = center_vertex->point();

		std::list<double> angles;
		Border_vertex_iterator it1 = first;
		Border_vertex_iterator it2 = it1;
		++it2;

		for (; it2 != last;++it1,++it2)
		{       
			Point_3 p1 = (*it1)->point();
			Point_3 p2 = (*it2)->point();

			angles.push_back(std::abs(Math::compute_angle_rad(p1,p0,p2)));		
		}
		it2 = first;
		angles.push_back(std::abs(Math::compute_angle_rad((*it1)->point(),p0,(*it2)->point())));	

		double sumSpaceAngle = std::accumulate(angles.begin(), angles.end(), 0.0);
		static const double PI = 3.14159265359;
		double tmp = 2*PI/sumSpaceAngle;
		double currAngle(0.0);
		std::list<double>::iterator it = angles.begin();
		for (; it!=angles.end(); ++it)
		{
			*it = (*it) * tmp;
			currAngle += *it;
			*it = currAngle;
		}

		it1 = first;	
		Vertex_handle vh = *it1;
		double radius1 = (*it1)->s();
		vh->uv(1.0,0.5);
		//mesh.set_vertex_parameterized(vh, true);
		++it1;

		for (it = angles.begin();it1 != last;++it1,++it)
		{       
			vh = *it1;
			double r = vh->s()/(radius1*2);
			double u = 0.5 + std::cos(*it) * r;
			double v = 0.5 + std::sin(*it) * r;
			vh->uv(u,v);
			//mesh.set_vertex_parameterized(vh, true);
		}	
		return true;
	}
};

//
// Class Circular_border_parameterizer_3
//
template<class Polyhedron_3>           //< 3D surface
class Circular_border_parameterizer_3 : public Border_parameterizer<Polyhedron_3>
{
public:
	/// Destructor of base class should be virtual.
	virtual ~Circular_border_parameterizer_3() {}

	// Default constructor, copy constructor and operator =() are fine

	/// Assign to mesh's border vertices a 2D position (i.e. a (u,v) pair)
	/// on border's shape. Mark them as "parameterized".
	bool parameterize_border(Border_vertex_iterator first, Border_vertex_iterator last, Vertex_handle center_vertex=NULL)
	{
		// Nothing to do if no border
		if (first == last)
			return false;

		// Compute the total border length
		double total_len = compute_border_length(first, last);
		if (total_len == 0)
			return false;

		const double PI = 3.14159265359;
		const double tmp = 2*PI/total_len;
		double len = 0.0;           // current position on circle in [0, total_len]
		for(Border_vertex_iterator it = first;it != last;++it)
		{
			Vertex_handle vh = *it;
			double angle = len*tmp; // current position on the circle in radians

			// map vertex on unit circle
			Point_2 uv;
			uv = Point_2(0.5+0.5*std::cos(-angle),0.5+0.5*std::sin(-angle));
			
			vh->uv(uv.x(), uv.y());
			//mesh.set_vertex_parameterized(it, true);

			// Get next iterator (looping)
			Border_vertex_iterator next = it;
			next++;
			if(next == last)
				next = first;

			// Add 'length' of it -> next vector to 'len'
			len += compute_edge_length(*it, *next);
		}
		
		return true;
	}

	// Protected operations
protected:
	/// Compute the length of an edge.
	virtual double compute_edge_length(Vertex_const_handle source,Vertex_const_handle target) = 0;

	// Private operations
private:
	/// Compute the total length of the border
	double compute_border_length(Border_vertex_iterator first, Border_vertex_iterator last)
	{
		double len = 0.0;
		for(Border_vertex_iterator it =first;it != last;++it)
		{
			// Get next iterator (looping)
			Border_vertex_iterator next = it;
			next++;
			if(next == last)
				next = first;

			// Add 'length' of it -> next vector to 'len'
			len += compute_edge_length(*it, *next);
		}
		return len;
	}
};

template<class Polyhedron_3>      //< 3D surface
class Circular_border_uniform_parameterizer_3
	: public Circular_border_parameterizer_3<Polyhedron_3>
{
protected:
	/// Compute the length of an edge.
	virtual double compute_edge_length(Vertex_const_handle source,Vertex_const_handle target)
	{
		/// Uniform border parameterization: points are equally spaced.
		return 1;
	}
};
template<class Polyhedron_3>           //< 3D surface
class Circular_border_arc_length_parameterizer_3
	: public Circular_border_parameterizer_3<Polyhedron_3>
{
protected:
	/// Compute the length of an edge.
	virtual double compute_edge_length(Vertex_const_handle source,Vertex_const_handle target)
	{
		/// Arc-length border parameterization: (u,v) values are
		/// proportional to the length of border edges.
		Vector_3 v = target->point() - source->point();
		return std::sqrt(v*v);
	}
};

DGAL_END_NAMESPACE
#endif//DGAL_BORDER_PARIMETERIZE__H__