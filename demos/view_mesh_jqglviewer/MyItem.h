#ifndef MY_ITEM_H__
#define MY_ITEM_H__

#include <CGAL/Cartesian.h>
#include <CGAL/IO/Color.h>
#include <CGAL/Polyhedron_3.h>

// CGAL kernel
typedef double number_type;
typedef CGAL::Cartesian<number_type>       Kernel;

// functor: compute face normal 
struct Face_normal
{
	template <class Face>
	void operator()(Face& f)
	{
		Kernel::Vector_3 sum = CGAL::NULL_VECTOR;
		Face::Halfedge_around_facet_circulator h = f.facet_begin();
		do
		{
			Kernel::Vector_3 normal = CGAL::cross_product(
				h->next()->vertex()->point() - h->vertex()->point(),
				h->next()->next()->vertex()->point() - h->next()->vertex()->point());
			double sqnorm = normal * normal;
			if(sqnorm != 0)
				normal = normal / (float)std::sqrt(sqnorm);
			sum = sum + normal;
		}
		while(++h != f.facet_begin());
		float sqnorm = (float) (sum * sum);
		if(sqnorm != 0.0)
			f.normal_ = sum / std::sqrt(sqnorm);
		else
		{
			f.normal_ = CGAL::NULL_VECTOR;
			//TRACE("degenerate face\n");
		}
	}
};


// functor: compute vertex normal 
struct Vertex_normal
{
	template <class Vertex>
	void operator()(Vertex& v)
	{
		Kernel::Vector_3 normal = CGAL::NULL_VECTOR;
		Vertex::Halfedge_around_vertex_const_circulator pHalfedge = v.vertex_begin();
		Vertex::Halfedge_around_vertex_const_circulator begin = pHalfedge;
		CGAL_For_all(pHalfedge,begin) 
			if(!pHalfedge->is_border())
				normal = normal + pHalfedge->facet()->normal_;
		float sqnorm = (float) (normal * normal);
		if(sqnorm != 0.0f)
			v.normal_ = normal / (float)std::sqrt(sqnorm);
		else
			v.normal_ = CGAL::NULL_VECTOR;
	}
};

// a refined facet with a normal, a tag
template <class Refs>
class My_face : public CGAL::HalfedgeDS_face_base<Refs>
{
public:
	// tag
	int tag_; //-1 for // uninitialized
	// normal
	Kernel::Vector_3 normal_;
	My_face():tag_(-1){}
};

// a refined vertex with a normal, a tag, an index, a scalar 
template <class Refs, class T, class P>
class My_vertex : public CGAL::HalfedgeDS_vertex_base<Refs, T, P>
{
public:
	// index
	int index_;
	number_type scalar_;
	// misc
	int tag_;
	// normal
	Kernel::Vector_3 normal_;

public:
	// life cycle
	My_vertex(){init();}
	// repeat mandatory constructors
	My_vertex(const P& pt)	: CGAL::HalfedgeDS_vertex_base<Refs, T, P>(pt){
		init();
	}

	void init(){
		index_ = -1;           // uninitialized
		tag_ = -1;             // uninitialized		
	}
};

struct My_items : public CGAL::Polyhedron_items_3
{
	// wrap vertex
	template <class Refs, class Traits>
	struct Vertex_wrapper
	{
		typedef typename Traits::Point_3  Point;
		typedef My_vertex<Refs, CGAL::Tag_true, Point> Vertex;		
	};

	// wrap face
	template <class Refs, class Traits>
	struct Face_wrapper
	{
		typedef My_face<Refs> Face;
	};
};

typedef CGAL::Polyhedron_3<Kernel, My_items>          Polyhedron;
//typedef CGAL::Polyhedron_3<Kernel>          Polyhedron;

#endif
