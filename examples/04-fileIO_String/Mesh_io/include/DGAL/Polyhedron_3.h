#ifndef DGAL_POLYHEDRON_H
#define DGAL_POLYHEDRON_H

#pragma warning (push)
#pragma warning(disable : 4244 4267 4503 4819 4996 )
#include <CGAL/Cartesian.h>	
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <DGAL/config.h>
#include <list>

DGAL_BEGIN_NAMESPACE

// CGAL kernel
typedef double number_type;
typedef CGAL::Cartesian<number_type> Kernel;

// a refined facet with a normal and a tag
template <class Refs, class T, class P, class Norm>
class Facet : public CGAL::HalfedgeDS_face_base<Refs, T>
{
	typedef Kernel::Vector_3 Vector_3;
	typedef Kernel::Point_3  Point_3;
private:
	// tag
	int tag_; //-1 for // uninitialized

	// normal
	Norm normal_;

	Point_3 center_;

public:

	// life cycle
	// no constructors to repeat, since only
	// default constructor mandatory

	Facet():tag_(-1)
	{
	}

	// tag
	const int&  tag() { return tag_; }
	void tag(const int& t)  { tag_ = t; }

	// center
	Point_3& center() { return center_; }
	const Point_3& center() const { return center_; }

	// normal
	typedef Norm Normal_3;
	Normal_3& normal() { return normal_; }
	const Normal_3& normal() const { return normal_; }

	// distance
	double distance(Point_3& point) const
	{
		Vector_3 vec = (point-center_);
		return std::sqrt(vec*vec);
	}
};
template <class Refs, class T, class P, class Norm>
class Vertex : public CGAL::HalfedgeDS_vertex_base<Refs, T, P>
{
private:
	bool is_constrained_;
	bool is_parameterized_;
	double u_;   // parameter
	double v_;	
	double s_;
	int index_;
	// tag
	int tag_;
	// normal
	Norm normal_;
public:
	// life cycle
	// -1 means uninitialized
	Vertex():
	  is_constrained_(false),is_parameterized_(false),
		  u_(-1.0), v_(-1.0),s_(-1.0),
		  index_(-1),tag_(-1)
	  {
	  }
	  // repeat mandatory constructors
	  Vertex(const P& pt):
	  CGAL::HalfedgeDS_vertex_base<Refs, T, P>(pt),
		  is_constrained_(false),is_parameterized_(false),
		  u_(-1.0), v_(-1.0),s_(-1.0),
		  index_(-1),tag_(-1)
	  {		
	  }

	  bool is_parameterized() const {	return is_parameterized_;}
	  void is_parameterized(bool in){is_parameterized_ = in;}

	  // normal
	  typedef Norm Normal_3;
	  Normal_3& normal() { return normal_; }
	  const Normal_3& normal() const { return normal_; }

	  // u,v texture
	  const double u() const {return u_;}
	  const double v() const {  return v_; }
	  const double s() const {  return s_; }
	  void uv(double u, double v)  { u_=u; v_=v; }
	  void u(double in){u_=in;}
	  void v(double in){v_=in;}
	  void s(double in){s_=in;}

	  bool is_constrained() const {
		  return is_constrained_; 
	  }
	  void is_constrained(bool in) { 
		  is_constrained_ = in; 
	  }

	  int index() const {return index_;}
	  void index(int in){ index_ = in;}
	  int tag() const {return tag_;}
	  void tag(int in){ tag_ = in;}
};

// a refined halfedge with a general tag, error, uv, is_parameterized, index
template <class Refs, class Tprev, class Tvertex, class Tface, class Norm>
class Halfedge : public CGAL::HalfedgeDS_halfedge_base<Refs,Tprev,Tvertex,Tface>
{
public:
	double err_;
private:

	// tag
	int tag_; 
	double length_;

	// parameterization
	bool is_parameterized_;
	double u_;                 // texture coordinates
	double v_;
	int index_;                // for parameterization

public:

	// life cycle
	Halfedge()
	{
		tag_ = -1;             // uninitialized

		u_ = 0.0;
		v_ = 0.0;
		index_ = -1;           // uninitialized
		is_parameterized_ = false;
	}

	// tag
	const int& tag() const { return tag_;  }
	int& tag() { return tag_;  }
	void tag(const int& t)  { tag_ = t; }

	// precomputed distance
	double length() const { return length_; }
	void length(double length) { length_ = length; }

	// texture coordinates
	double u() const { return u_; }
	double v() const { return v_; }
	void uv(double u, double v) { u_ = u; v_ = v; }
	void u(double in){u_=in;}
	void v(double in){v_=in;}

	// param.
	bool is_parameterized() const { return is_parameterized_; }
	void is_parameterized(bool is)  { is_parameterized_ = is; }

	// index
	int index() const { return index_; }
	void index(int i) { index_ = i; }
};

// An items type using my face.
struct Polyhedron_items_3 : public CGAL::Polyhedron_items_3 {
	// wrap vertex
	template <class Refs, class Traits>
	struct Vertex_wrapper
	{
		typedef typename Traits::Point_3  Point;
		typedef typename Traits::Vector_3 Normal;
		typedef Vertex<Refs, CGAL::Tag_true,Point, Normal> Vertex;
	};
	// wrap halfedge
	template <class Refs, class Traits>
	struct Halfedge_wrapper
	{
		typedef typename Traits::Vector_3 Normal;
		typedef Halfedge<Refs,
			CGAL::Tag_true,
			CGAL::Tag_true,
			CGAL::Tag_true,
			Normal> Halfedge;
	};
	// wrap face
	template <class Refs, class Traits>
	struct Face_wrapper
	{
		typedef typename Traits::Point_3  Point;
		typedef typename Traits::Vector_3 Normal;
		typedef Facet<Refs,
			CGAL::Tag_true,
			Point,
			Normal> Face;
	};
};



template
< 
	class PolyhedronTraits_3,
    class PolyhedronItems_3 = Polyhedron_items_3,
#ifndef CGAL_CFG_NO_TMPL_IN_TMPL_PARAM
    template < class T, class I, class A>
#endif
	class T_HDS = CGAL::HalfedgeDS_default, 
	class Alloc = CGAL_ALLOCATOR(int)
>
class Polyhedron_3 
	: public CGAL::Polyhedron_3<PolyhedronTraits_3, PolyhedronItems_3, T_HDS, Alloc >
{
public:
	typedef typename PolyhedronTraits_3                    Traits;
	typedef typename PolyhedronItems_3                     Items;
	typedef typename Polyhedron_3<Traits, Items, T_HDS, Alloc>
		                                                   Self;
	typedef typename CGAL::Polyhedron_3<Traits, Items, T_HDS, Alloc >
		                                                   Base;

	typedef CGAL::I_Polyhedron_derived_items_3<Items>      Derived_items;
#ifndef CGAL_CFG_NO_TMPL_IN_TMPL_PARAM
    typedef T_HDS< Traits, Derived_items, Alloc>           HDS;
#else
    typedef typename T_HDS::template HDS< Traits, Derived_items, Alloc>  
		                                                   HDS;
#endif

	typedef typename HDS::Vertex                           Vertex;
	typedef typename HDS::Halfedge                         Halfedge;
	typedef typename HDS::Face                             Face;

	// Number type
	typedef typename Traits::FT                           NT;
	typedef typename Traits::Point_2                      Point_2;
	typedef typename Traits::Point_3                      Point_3;
	typedef typename Traits::Vector_2                     Vector_2;
	typedef typename Traits::Vector_3                     Vector_3;
	typedef typename Traits::Line_2                        Line_2;
	typedef typename Traits::Line_3                        Line_3;
	typedef typename Traits::Plane_3                       Plane_3;
	typedef typename Traits::Direction_3                   Direction_3;
	typedef typename Traits::Aff_transformation_3          Aff_transformation_3;
	typedef typename Traits::Iso_cuboid_3                  Iso_cuboid_3;
		
	// compute facet center
	struct Facet_center  // (functor)
	{
		template<class Facet>
		void operator()(Facet& f)
		{
			Vector_3 vec(0.0,0.0,0.0);
			int degree = 0;
			typedef typename Facet::Halfedge_around_facet_circulator Circ;
			Circ h = f.facet_begin(); Circ end = h;
			do
			{
				vec = vec + (h->vertex()->point()-CGAL::ORIGIN);
				degree++;
			}
			while (++h != end);
			f.center() = CGAL::ORIGIN + (vec/degree);
		}
	};
	// compute facet normal 
	struct Facet_normal // (functor)
	{
		template <class Facet>
		void operator()(Facet& f)
		{
			typename Facet::Normal_3 sum = CGAL::NULL_VECTOR;
			typename Facet::Halfedge_around_facet_circulator h = f.facet_begin();
			do
			{
				typename Facet::Normal_3 normal = CGAL::cross_product(
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
				f.normal() = sum / std::sqrt(sqnorm);
			else
			{
				f.normal() = CGAL::NULL_VECTOR;
				//TRACE("degenerate face\n");
			}
		}
	};

	// compute vertex normal 
	struct Vertex_normal // (functor)
	{
		template <class Vertex>
		void operator()(Vertex& v)
		{
			typename Vertex::Normal_3 normal = CGAL::NULL_VECTOR;
			Vertex::Halfedge_around_vertex_circulator pHalfedge = v.vertex_begin();
			Vertex::Halfedge_around_vertex_circulator begin = pHalfedge;
			CGAL_For_all(pHalfedge,begin) 
				if(!pHalfedge->is_border())
					normal = normal + pHalfedge->facet()->normal();
			float sqnorm = (float) (normal * normal);
			if(sqnorm != 0.0f)
				v.normal() = normal / (float)std::sqrt(sqnorm);
			else
				v.normal() = CGAL::NULL_VECTOR;
		}
	};

	struct Index_predicate{
		Index_predicate(int ind):m_ind(ind){}
		template<class Element>
		bool operator()(Element& elem){
			return (m_ind == elem.index());
		}
		int m_ind;
	};
	struct Generate_index{
		Generate_index():m_ind(-1){}
		template<class Element>
		void operator()(Element& elem){
			elem.index(++m_ind);
		}
		int m_ind;
	};
	struct Set_tag{
		Set_tag(int v):m_tag(v){}
		template<class Element>
		void operator()(Element& elem){
			elem.tag(m_tag);
		}
		void operator()(Vertex_handle& elem){
			elem->tag(m_tag);
		}
		int m_tag;
	};
	struct Edge_length
	{
		template <class Halfedge_T>
		void operator()(Halfedge_T& he)
		{
			Vector_3 v = he.vertex()->point() - he.opposite()->vertex()->point();
			double len = std::sqrt(v*v);
			he.length(len); he.opposite()->length(len);
		}
	};
	template<class Vector_T> 
	struct Edge_mean_length
	{
		Edge_mean_length(Vector_T& result):m_result(result){}

		template <class Vertex_T>
		void operator()(Vertex_T& vertex)
		{
			Vertex_T::Halfedge_around_vertex_circulator hc, end;
			hc = vertex.vertex_begin();
			end  =  hc;	

			double s(0),degree(0);
			CGAL_For_all(hc,end)
			{
				s += hc->length();
				++degree;
			}
			s = s/degree;
			m_result[vertex.index()] = s;
		}
		Vector_T& m_result;
	};
public:
	void index_vertices(){
		std::for_each(vertices_begin(), vertices_end(), Generate_index());
	}
	void index_halfedges(){
		std::for_each(halfedges_begin(), halfedges_end(), Generate_index());
	}
	//void index_halfedges(){std::for_each(halfedges_begin(), halfedges_end(), Generate_index());}
	//void index_edges(){std::generate(edges_begin(), edges_end(), Generate_index(-1));}
	//void index_facets(){std::generate(facets_begin(), facets_end(), Generate_index(-1));}	

	/// Preconditions:
	/// index vertices first
	Vertex_handle get_vertex(int index)
	{
		Vertex_iterator result = std::find_if(vertices_begin(), vertices_end(), Index_predicate(index));
		return result;
	}

	// tag (vertices, halfedges, facets)
	void tag_vertices(const int tag){
		std::for_each(vertices_begin(), vertices_end(), Set_tag(tag));
	}
	void tag_halfedges(const int tag){
		std::for_each(halfedges_begin(), halfedges_end(), Set_tag(tag));
	}
	void tag_facets(const int tag){
		std::for_each(facets_begin(), facets_end(), Set_tag(tag));
	}
	// normals (per facet, then per vertex)
	void compute_normals_per_facet()
	{
		std::for_each(facets_begin(),facets_end(),Facet_normal());
	}
	void compute_normals_per_vertex()
	{
		std::for_each(vertices_begin(),vertices_end(),Vertex_normal());
	}
	void compute_normals()
	{
		compute_normals_per_facet();
		compute_normals_per_vertex();
	}

	bool is_parameterized(){return vertices_begin()->is_parameterized();}
	bool is_border(Vertex_handle vh){
		Halfedge_around_vertex_circulator he = vh->vertex_begin();
		Halfedge_around_vertex_circulator end = he;
		if(he == NULL) // isolated vertex
			return true;
		CGAL_For_all(he,end)
			if(he->is_border())
				return true;
		return false;
	}
	bool is_inner_facet(Facet_handle fh)
	{
		typedef Halfedge_around_facet_circulator circ;
		circ h = fh->facet_begin();
		circ end = h;
		do
		{
			if(h->opposite()->is_border())
				return false;
		}
		while(++h != end);
		return true;
	}


	///?要不要constrained vertex呢？
	void add_constrained_vertex(Vertex_handle vh)
	{
		std::list<Vertex_handle>::iterator it = std::find(m_constrained_vertices.begin(),
			m_constrained_vertices.end(), vh);
		if ( it == m_constrained_vertices.end()) 
		{
			m_constrained_vertices.push_back(vh);
			vh->is_constrained(true);
		}
		else 
		{
			m_constrained_vertices.erase(it);
			vh->is_constrained(false);
		}
	}
	void set_constrained_vertices(std::list<Vertex_handle>& cvs)
	{
		m_constrained_vertices.clear();
		for(std::list<Vertex_handle>::iterator it = cvs.begin(); it!=cvs.end(); ++it)
		{
			Vertex_handle cv = *it;
			m_constrained_vertices.push_back(cv);
			cv->is_constrained(true);
		}
	}
	void clear_constrained_vertices()
	{	
		for(std::list<Vertex_handle>::iterator it = m_constrained_vertices.begin(); it!=m_constrained_vertices.end(); ++it)
		{
			Vertex_handle vh = *it;
			vh->is_constrained(false);
		}
		m_constrained_vertices.clear();
	}

	std::list<Vertex_handle>& constrained_vertices(){return m_constrained_vertices;}

	Vertex_handle vertex_min(int coord, double &min){
		struct Min_vertex{
			Min_vertex(int coord, double min):m_coord(coord),m_min(min){}
			void operator()(Vertex& v){
				double value = v.point()[m_coord];
				if(value < m_min)
				{
					m_min = value;
					m_vh = &v;
				}
			}
			int m_coord;
			double m_min;
			Vertex_handle m_vh;
		};

		//////////////////////////////////////////////////////////////////////////
		CGAL_assertion(size_of_vertices() > 0);
		min = vertices_begin()->point()[coord];
		Min_vertex mv = std::for_each(vertices_begin(), vertices_end(), Min_vertex(coord,min));
		min = mv.m_min;
		return mv.m_vh;
	}
	Vertex_handle vertex_max(int coord, double &max)
	{
		struct Max_vertex{
			Max_vertex(int coord, double max):m_coord(coord),m_max(max){}
			void operator()(Vertex& v){
				double value = v.point()[m_coord];
				if(value > m_max)
				{
					m_max = value;
					m_vh = &v;
				}
			}
			int m_coord;
			double m_max;
			Vertex_handle m_vh;
		};

		//////////////////////////////////////////////////////////////////////////
		CGAL_assertion(size_of_vertices() > 0);
		max = vertices_begin()->point()[coord];
		Max_vertex mv = std::for_each(vertices_begin(), vertices_end(), Max_vertex(coord,max));
		max = mv.m_max;
		return mv.m_vh;
	}


	void farthest_point_aligned(Vertex_handle &pVertexMin, Vertex_handle &pVertexMax)
	{
		double xmin,xmax,ymin,ymax,zmin,zmax;
		Vertex_handle pVertex_xMin = vertex_min(0,xmin);
		Vertex_handle pVertex_xMax = vertex_max(0,xmax);
		Vertex_handle pVertex_yMin = vertex_min(1,ymin);
		Vertex_handle pVertex_yMax = vertex_max(1,ymax);
		Vertex_handle pVertex_zMin = vertex_min(2,zmin);
		Vertex_handle pVertex_zMax = vertex_max(2,zmax);
		double xdiff = xmax-xmin;
		double ydiff = ymax-ymin;
		double zdiff = zmax-zmin;
		if (xdiff >= max(ydiff,zdiff))
		{
			pVertexMin = pVertex_xMin;
			pVertexMax = pVertex_xMax;
		}
		else if (ydiff >= max(xdiff,zdiff))
		{
			pVertexMin = pVertex_yMin;
			pVertexMax = pVertex_yMax;
		}
		else
		{
			pVertexMin = pVertex_zMin;
			pVertexMax = pVertex_zMax;
		}
	}
	Halfedge_handle get_halfedge(Vertex_handle source, Vertex_handle target)
	{
		assert(source != NULL);
		assert(target != NULL);

		Halfedge_around_vertex_circulator cir = target->vertex_begin(),cir_end = cir;
		CGAL_For_all(cir, cir_end)
			if (cir->opposite()->vertex() == source)
				return cir;

		assert(false);              // error if we reach this point
		return NULL;
	}
	// compute distance from facet center to halfedge center
	double distance(Facet_handle pFacet, Halfedge_handle pHalfedge)
	{
		// we assume
		Point_3 center_facet = pFacet->center();

		Vector_3 v = (pHalfedge->opposite()->vertex()->point()
			- pHalfedge->vertex()->point());
		Point_3 center_halfedge = pHalfedge->vertex()->point() + (v/2);
		Vector_3 d = center_facet-center_halfedge;
		return std::sqrt(d*d);
	}
	void transform(Aff_transformation_3 t){
		for(Vertex_iterator vi = vertices_begin(); vi !=  vertices_end(); ++vi){
			Point_3& p = vi->point();
			pVertex->point() = p.transform(t);
		}
	}

	// get closest inner facet
	Facet_handle get_closest_inner_facet(Point_3& point)
	{
		Facet_iterator pFace = facets_begin();
		Facet_handle pClosest = pFace;
		double min = pFace->distance(point);
		for(;pFace != facets_end();++pFace){
			if(is_inner(pFace)){
				double distance = pFace->distance(point);
				if(distance < min){
					pClosest = pFace;
					min = distance;
				}
			}
		}
		return pClosest;
	}

	template<class Vector_T> void compute_mean_spokeLen(Vector_T& result)
	{
		std::for_each(edges_begin(), edges_end(), Edge_length());
		std::for_each(vertices_begin(), vertices_end(), Edge_mean_length<Vector_T>(result));
	}
	void compute_spoke_len(std::vector<double> & result, double threshold, int option)
	{
		std::for_each(edges_begin(), edges_end(), Edge_length());

		result.reserve(size_of_vertices());
		double s, mn(0.0), mx(0.0);
		for(Vertex_iterator vi = vertices_begin(); vi!=vertices_end(); ++vi)
		{
			Halfedge_vertex_circulator  hc = vi->vertex_begin();
			Halfedge_vertex_circulator end  =  hc;	
			mn = threshold; mx=0;
			CGAL_For_all(hc,end)
			{
				if ( hc->u() < mn)
				{
					mn = hc->u();
				}
				if (hc->u() > mx)
				{
					mx = hc->u();
				}
			}
			if (option==0)
			{
				s = mn;
			}
			else
			{
				s = (mn+mx)/2.0;
			}
			result.push(s);
		}
	}

	
	//////////////////////////////////////////////////////////////////////////
	/// is valid in geometry for 2D case.
	bool is_domain_valid()
	{
		Halfedge_around_facet_circulator h;
		Vertex_handle vhs[3];
		
		/////////////////////////////////////
		Facet_iterator iter = facets_begin();
		h = iter->facet_begin();
		int i(0);
		do
		{
			vhs[i] = h->vertex();
			++i;
		}
		while(++h != iter->facet_begin());
		Vector_3 v1 = vhs[1]->point()-vhs[0]->point();
		Vector_3 v2 = vhs[2]->point()-vhs[1]->point();
		Vector_3 v3 = CGAL::cross_product(v1,v2);
		bool isPositive(v3.z()>0);

		++iter;
        for(; iter != facets_end(); ++iter)
		{
			h = iter->facet_begin();
			int i(0);
			do
			{
				vhs[i] = h->vertex();
				++i;
			}
			while(++h != iter->facet_begin());
			v1 = vhs[1]->point()-vhs[0]->point();
		    v2 = vhs[2]->point()-vhs[1]->point();
			v3 = CGAL::cross_product(v1,v2);

			bool tmp(v3.z()>0);
			if (isPositive != tmp)
			{
				return false;
			}
		}
        
        return true;
	}
	//precondition: each vertex has uv
	void set_halfedges_uv()
	{
		for(Vertex_iterator pVertex = vertices_begin();pVertex != vertices_end(); ++pVertex)
		{
			// Loop over all incident halfedges
			Halfedge_around_vertex_circulator cir     = pVertex->vertex_begin(),
											  cir_end = cir;
			CGAL_For_all(cir, cir_end)
			{
				cir->uv(pVertex->u(), pVertex->v());
				cir->is_parameterized(true);
			}
		}
	}
	void set_vertices_uv()
	{
		for(Vertex_iterator pVertex = vertices_begin();pVertex != vertices_end(); ++pVertex)
		{
			Halfedge_around_vertex_circulator cir     = pVertex->vertex_begin();
			pVertex->uv(cir->u(), cir->v());
			pVertex->is_parameterized(true);
		}
	}
	void setUVByMesh(Polyhedron_3* domain)
	{
		Vertex_iterator it1 = vertices_begin();
		Vertex_iterator it2 = domain->vertices_begin();	
		for ( ; it1 != vertices_end(); ++it1,++it2)
		{
			Point_3 p3 = it2->point();
			it1->uv(p3.x(), p3.y());
			it1->is_parameterized(true);
		}
	}
	void setUVBySelf()
	{
		Vertex_iterator it1 = vertices_begin();
		for ( ; it1 != vertices_end(); ++it1)
		{
			Point_3 p3 = it1->point();
			it1->uv(p3.x(), p3.y());
			it1->is_parameterized(true);
		}
	}
	
	//typename std::list<Vertex_handle>::const_iterator main_border_begin() const { return m_main_border.begin();}
	//typename std::list<Vertex_handle>::const_iterator main_border_end()   const { return m_main_border.end();}
	typename std::list<Vertex_handle>::iterator main_border_begin() { return m_main_border.begin();}
	typename std::list<Vertex_handle>::iterator main_border_end()   { return m_main_border.end();}
	typename std::list<Vertex_handle> main_border()   { return m_main_border;}

	std::list<Vertex_handle> extract_longest_border(){
		std::list<Vertex_handle> longest_border;    // returned list
		double                   max_len = 0;       // length of longest_border

		// Tag all vertices as unprocessed
		const int tag_free = 0;
		const int tag_done = 1;
		std::for_each(vertices_begin(), vertices_end(), Set_tag(tag_free));

		// find all closed borders and keep longest one
		int nb = 0;
		while (1)
		{
			// Find a border tagged as "free" and tag it as "processed"
			std::list<Vertex_handle> border = find_free_border(tag_free, tag_done);
			if(border.empty())
				break;

			// compute  total len of 'border'
			double len = 0.0;
			std::list<Vertex_handle>::iterator it;
			for(it = border.begin(); it != border.end(); it++)
			{
				// Get next iterator (looping)
				std::list<Vertex_handle>::iterator next = it;
				next++;
				if (next == border.end())
					next = border.begin();

				Vector_3 vect = (*next)->point() - (*it)->point();
				len += std::sqrt(vect*vect);
			}

			// Keep 'border' if longer
			if (len > max_len)
			{
				longest_border = border;
				max_len = len;
			}

			++nb;
		}

		return longest_border;
	}
	std::list<Vertex_handle>& extract_main_border_by_length()
	{
		m_main_border = extract_longest_border();
		return m_main_border;
	}
	Iso_cuboid_3 compute_bounding_box()
	{
		if(size_of_vertices() == 0){
			return m_bbox;
		}

		NT xmin,xmax,ymin,ymax,zmin,zmax;
		Vertex_iterator pVertex = vertices_begin();
		xmin = xmax = pVertex->point().x();
		ymin = ymax = pVertex->point().y();
		zmin = zmax = pVertex->point().z();
		for(;pVertex !=  vertices_end();++pVertex)
		{
			const Point_3& p = pVertex->point();

			xmin =  min(xmin,p.x());
			ymin =  min(ymin,p.y());
			zmin =  min(zmin,p.z());

			xmax =  max(xmax,p.x());
			ymax =  max(ymax,p.y());
			zmax =  max(zmax,p.z());
		}
		m_bbox = Iso_cuboid_3(xmin,ymin,zmin,
			xmax,ymax,zmax);
		return m_bbox;
	}
	Iso_cuboid_3& bbox(){return m_bbox;}

	void compare_deformation(Polyhedron_3* initial,int from)
	{
		if (1==from)
		{
			//edge difference
			Halfedge_iterator h1=halfedges_begin();
			Halfedge_iterator h2=initial->halfedges_begin();
			for(;h1!= halfedges_end();++h1,++h2)
			{
				//edge
				Vector_3 vec=h1->vertex()->point()-h1->opposite()->vertex()->point();
				double l1=CGAL::sqrt(vec*vec);
				vec=h2->vertex()->point()-h2->opposite()->vertex()->point();
				double l2=CGAL::sqrt(vec*vec);
				h1->u(abs(l1-l2)/l1);
			}
			//save to vertex
			Vertex_iterator v1 = vertices_begin();	
			for ( ; v1 != vertices_end(); ++v1)
			{
				v1->s(0);
				Halfedge_vertex_circulator pHalfEdge = v1->vertex_begin();
				Halfedge_vertex_circulator end=pHalfEdge;
				int num=0;
				CGAL_For_all(pHalfEdge,end)
				{
					v1->s(v1->s()+pHalfEdge->u());
					++num;
				}
				v1->s(v1->s()/num/2);
			}
		}
		else if (2==from)
		{
			//angle difference
			Halfedge_iterator h1=halfedges_begin();
			Halfedge_iterator h2=initial->halfedges_begin();
			compute_normals_per_facet();
			initial->compute_normals_per_facet();
			for(;h1!= halfedges_end();++h1,++h2)
			{
				//angle
				if (!h1->is_border()&&!h1->opposite()->is_border())
				{
					double l1=acos(h1->facet()->normal()*h1->opposite()->facet()->normal());
					double l2=acos(h2->facet()->normal()*h2->opposite()->facet()->normal());
					h1->u(abs(l1-l2));
				}		
			}
			//save to vertex
			Vertex_iterator v1 = vertices_begin();	
			for ( ; v1 != vertices_end(); ++v1)
			{
				v1->s(0);
				Halfedge_vertex_circulator pHalfEdge = v1->vertex_begin();
				Halfedge_vertex_circulator end=pHalfEdge;
				int num=0;
				CGAL_For_all(pHalfEdge,end)
				{
					v1->s(v1->s()+pHalfEdge->u());
					++num;
				}
				v1->s(v1->s()/num);
			}
		} 
		else
		{
			//vertex difference
			Vertex_iterator v1 = vertices_begin();
			Vertex_iterator v2 = initial->vertices_begin();	
			for ( ; v1 != vertices_end(); ++v1,++v2)
			{
				Vector_3 vec =v1->point()-v2->point();
				v1->s(CGAL::sqrt(vec*vec));
			}		
		}	
	}

	
	Polyhedron_3* deep_copy(bool b_create_domain=false)
	{
		Polyhedron_3* result = new Polyhedron_3();
		Copy_Construct_Mesh<Polyhedron_3> bd(this, b_create_domain);
		result->delegate(bd);	
		return result;
	}
	Polyhedron_3* build_domain(){
		return deep_copy(true);    
	}

private:
	/// Find a border tagged as "free" and tag it as "processed".
    /// Return an empty list if not found.
	std::list<Vertex_handle> find_free_border(int tag_free, int tag_done)
	{
		std::list<Vertex_handle> border;    // returned list

		// get any border vertex with "free" tag
		Vertex_handle seed_vertex = NULL;
		// 为了和DGAL::Parameterization_polyhedron_adaptor_3和CGAL::Parameterization_polyhedron_adaptor_3的兼容
		// 放弃下面直接从边界取点的方法，采用比较复杂的
		//for (Edge_iterator ei=border_edges_begin(); ei!=edges_end(); ++ei)
		//{
		//	Vertex_handle vh = ei->vertex();
		//	if( vh->tag() ==tag_free)
		//	{
		//		seed_vertex = vh;
		//           break;
		//	}
		//}
		for (Vertex_iterator vi=vertices_begin(); vi!=vertices_end(); ++vi)
		{
			Vertex_handle vh = vi;
			if( is_border(vh) && vh->tag() ==tag_free)
			{
				seed_vertex = vh;
				break;
			}
		}

		if (seed_vertex == NULL)
			return border;                  // return empty list

		// Get the border containing seed_vertex
		border = get_border(seed_vertex);

		// Tag border vertices as "processed"
		std::list<Vertex_handle>::iterator it;
		for(it = border.begin(); it != border.end(); it++)
		{
			Vertex_handle vh = *it;
			vh->tag(tag_done);
		}

		return border;
	}
	// Return the border containing seed_vertex.
    /// Return an empty list if not found.
    std::list<Vertex_handle> get_border(Vertex_handle seed_vertex)
	{
		std::list<Vertex_handle> border;    // returned list

		Halfedge_around_vertex_circulator pHalfedge = seed_vertex->vertex_begin();
		Halfedge_around_vertex_circulator end       = pHalfedge;

		// if isolated vertex
		if (pHalfedge == NULL) {
			border.push_back(seed_vertex);
			return border;
		}

		// Get seed_vertex' border halfedge
		Halfedge_handle  seed_halfedge = NULL;
		CGAL_For_all(pHalfedge,end) {
			if(pHalfedge->is_border()) {
				seed_halfedge = pHalfedge;
				break;
			}
		}

		// if inner vertex
		if (seed_halfedge == NULL)
			return border;                  // return empty list

		// Add seed vertex
		border.push_back(seed_vertex);

		// fill border
		int size = 1;
		Halfedge_handle current_halfedge = seed_halfedge;
		do
		{
			// Stop if end of loop
			Halfedge_handle next_halfedge = current_halfedge->next();
			Vertex_handle next_vertex = next_halfedge->vertex();
			if(next_vertex == seed_vertex)
				break;

			// Add vertex
			border.push_back(next_vertex);

			current_halfedge = next_halfedge;
			size++;
		}
		while(1);

		return border;
	}

private:
	// bounding box
	Iso_cuboid_3 m_bbox;
	std::list<Vertex_handle> m_constrained_vertices;
	std::list<Vertex_handle> m_main_border;
};

template <class Polyhedron>
class Copy_Construct_Mesh : public CGAL::Modifier_base<typename Polyhedron::HDS> {
public:
	typedef typename Polyhedron::HDS               HDS;
	typedef typename Polyhedron::Vertex            Vertex;	
	typedef typename Vertex::Point                 Point;
	typedef typename Polyhedron::Vertex_handle     Vertex_handle;	
	typedef typename Polyhedron::Vertex_iterator   Vertex_iterator;	

	typedef typename Polyhedron::Facet_handle      Facet_handle;
	typedef typename Polyhedron::Facet_iterator    Facet_iterator;	

	typedef typename Polyhedron::Halfedge_around_vertex_circulator   Halfedge_vertex_circulator;
	typedef typename Polyhedron::Halfedge_around_facet_circulator    Halfedge_facet_circulator;

	typedef typename Polyhedron::Halfedge_handle                     Halfedge_handle;

	typedef typename CGAL::Polyhedron_incremental_builder_3<HDS>              builder;

public:
	bool is_build_domain_;
	Polyhedron* mesh_;//original mesh.	
public:
	Copy_Construct_Mesh(Polyhedron *mesh, bool is_build_domain=false)
		:mesh_(mesh),is_build_domain_(is_build_domain){}

	void operator()( HDS& hds){
		// Postcondition: `hds' is a valid polyhedral surface.
		builder B( hds, true);

		int vs(mesh_->size_of_vertices());
		int fs(mesh_->size_of_facets());
		int n_h(mesh_->size_of_halfedges());
		B.begin_surface( vs, fs, n_h);

		add_vertices(B);
		add_facets(B);

		B.end_surface();
	}
private:	
	void add_vertices(builder &B){
		int i(0);
		for(Vertex_iterator vit = mesh_->vertices_begin();vit != mesh_->vertices_end(); ++vit,++i)
		{
			Vertex_handle vh;
			if ( is_build_domain_)
				vh = B.add_vertex( Point(vit->u(), vit->v(), 0.0) );
			else
				vh = B.add_vertex(vit->point());

			vh->uv( vit->u(), vit->v());
			vh->s( vit->s());

			vit->index(i);
		}
	}
	void add_facets(builder &B){
		for(Facet_iterator fit = mesh_->facets_begin();fit != mesh_->facets_end(); ++fit)
		{		
			B.begin_facet();

			Halfedge_facet_circulator j = fit ->facet_begin();
			do {//访问pos的3个顶点
				B.add_vertex_to_facet(j->vertex()->index());
			} while ( ++j != fit->facet_begin());			

			Halfedge_handle hh = B.end_facet();
		}
	}

};

typedef Polyhedron_3<Kernel_default> Polyhedron_default;

DGAL_END_NAMESPACE
#pragma warning (pop)

#endif//DGAL_POLYHEDRON_H