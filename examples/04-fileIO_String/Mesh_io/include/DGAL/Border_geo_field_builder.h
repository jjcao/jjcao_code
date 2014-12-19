#ifndef BOEDER_GEO_FIELD_BUILDER__H
#define BOEDER_GEO_FIELD_BUILDER__H

#include <DGAL/config.h>
DGAL_BEGIN_NAMESPACE

/// Compute a geodesic distance field from specified vertices. 
/// If these vertices are all border vertices, it builds a boundary geodesic distance field.
/// Reffer to Wang C, Wang Y, Tang K, Yuen MMF. Reduce the stretch in surface flattening
/// by finding cutting paths to the surface boundary. Computer-Aided Design 2004;
/// 36(8):665-677.
template<
class Polyhedron_3
>
class Border_geo_field_builder
{
public:
	typedef typename Polyhedron_3 Polyhedron;
	typedef typename Polyhedron::Vertex_handle Vertex_handle;
	typedef typename Polyhedron::Vertex_handle Vertex_iterator;	
	typedef typename Polyhedron::Halfedge_around_vertex_circulator Halfedge_vertex_circulator;
	typedef typename Polyhedron::Vector_3 Vector_3;

	enum {TAG_FREE=0,TAG_DONE};
	struct Set_scalar{
		Set_scalar(double v):s(v){}
		template<class Element>
		void operator()(Element& elem){
			elem.s(s);
		}
		void operator()(Vertex_handle& elem){
			elem->s(s);
		}
		double s;
	};
	struct Edge_length_limit
	{
		Edge_length_limit():m_min_edge_len(std::numeric_limits<int>::max()),m_max_edge_len(0){}
		template <class Halfedge_T>
		void operator()(Halfedge_T& he)
		{
			Vector_3 v = he.vertex()->point() - he.opposite()->vertex()->point();
			double len = std::sqrt(v*v);
			he.length(len); he.opposite()->length(len);

			if (len>m_max_edge_len)
			{
				m_max_edge_len = len;
			}
			else if (len<m_min_edge_len)
			{
				m_min_edge_len = len;
			}
		}
		double m_max_edge_len;
		double m_min_edge_len;
	};

	Border_geo_field_builder(double threshold=0.1):m_minEdgelen(-1),m_maxEdgelen(-1),
							m_threshold(threshold){}
public:
	bool compute(Polyhedron* mesh){
		return compute(mesh,mesh->main_border());
	}

	bool compute(Polyhedron* mesh, std::list<Vertex_handle>& cvs)
	{
		bool result(false);
#ifdef OUTPUT_INFO
		CGAL::Timer timer;	timer.start();
#endif	
		double tmp(getEdgeRatio(mesh));
		if ( tmp> m_threshold)
			result = computeUniform(mesh, cvs);
		else
			result = computeNonuniform(mesh, cvs);
#ifdef OUTPUT_INFO
		std::cout << "border geo field build: " << timer.time() << " seconds." << std::endl;
		timer.reset();
#endif
		return result;
	}
	double getEdgeRatio(Polyhedron* mesh){
		if(m_maxEdgelen<0){//calculate min and max length of all edges			
			Edge_length_limit ell = std::for_each(mesh->edges_begin(), mesh->edges_end(), Edge_length_limit());
			m_maxEdgelen = ell.m_max_edge_len;
			m_minEdgelen = ell.m_min_edge_len;
		}

		return m_minEdgelen/m_maxEdgelen;
	}
	double threshold(){return m_threshold;}
	void threshold(double in){m_threshold=in;}
	void* matrix(){return 0;}
	void* factor(){return 0;}
protected:
	bool computeUniform(Polyhedron* mesh, std::list<Vertex_handle>& cvs){
		std::list<Vertex_handle> lv;
		init(mesh, cvs, lv);		

		float dist(0.0);
		float s(0.0);
		do{
			std::list<Vertex_handle> lv1;
			for(std::list<Vertex_handle>::iterator itor=lv.begin(); itor!=lv.end(); ++itor)
			{
				Vertex_handle vh = *itor;
				s = vh->s();
				Vertex_handle ovh;

				Halfedge_vertex_circulator  he  = vh->vertex_begin();
				Halfedge_vertex_circulator  edgeE  = he;
				CGAL_For_all(he,edgeE)
				{
					ovh = he->opposite()->vertex();
					dist = s + he->length();
					if ( dist < ovh->s())
					{
						ovh->s(dist);
					}				
					if ( ovh->halfedge()->tag()==TAG_FREE)
					{
						lv1.push_back(ovh);
						ovh->halfedge()->tag(TAG_DONE);
					}				
				}
			}        
			lv.clear();
			lv=lv1;
		}while(!lv.empty());

		return true;
	}
	bool computeNonuniform(Polyhedron* mesh, std::list<Vertex_handle>& cvs){
		std::list<Vertex_handle> lv;
		init(mesh, cvs,lv);

		float lamda(m_minEdgelen); 
		float dist(0.0);
		float s(0.0);
		do{
			std::list<Vertex_handle> lv1;

			for(std::list<Vertex_handle>::iterator itor=lv.begin(); itor!=lv.end(); ++itor)
			{
				Vertex_handle vh = *itor;
				s = vh->s();
				Vertex_handle ovh;

				Halfedge_vertex_circulator  he  = vh->vertex_begin();
				Halfedge_vertex_circulator  edgeE  = he;
				CGAL_For_all(he,edgeE)
				{
					ovh = he->opposite()->vertex();
					dist = s + he->length();
					if ( dist < ovh->s())
					{
						ovh->s(dist);
					}				
					if (( ovh->halfedge()->tag()==TAG_FREE)&&(ovh->s()<lamda))
					{
						lv1.push_back(ovh);
						ovh->halfedge()->tag(TAG_DONE);
					}

				}
			}

			for(std::list<Vertex_handle>::iterator itor=lv.begin(); itor!=lv.end(); ++itor)
			{

				Vertex_handle vh = *itor;
				Halfedge_vertex_circulator  he  = vh->vertex_begin();
				Halfedge_vertex_circulator  edgeE  = he;
				CGAL_For_all(he,edgeE)
				{
					if (he->opposite()->vertex()->halfedge()->tag()==TAG_FREE)
					{
						lv1.push_back (he->vertex());
						break;
					}
				}
			}

			lv.clear();
			lv=lv1;
			lamda+=m_minEdgelen;
		}while(!lv.empty());

		return true;
	}

	void init(Polyhedron* mesh,std::list<Vertex_handle>& cvs, std::list<Vertex_handle>& lv){
		double maxF(std::numeric_limits<double>::max());

		for(Vertex_iterator vh = mesh->vertices_begin(); vh != mesh->vertices_end(); ++vh)
		{
			vh->halfedge()->tag(TAG_FREE);
			vh->s(maxF);
		}

		for (std::list<Vertex_handle>::iterator it = cvs.begin(); it!=cvs.end();++it)
		{
			Vertex_handle vh = *it;
			lv.push_back(vh);
			vh->halfedge()->tag(TAG_DONE);
			vh->s(0.0);
		}
	}

private:
	double m_threshold;
	double m_minEdgelen;
	double m_maxEdgelen;
};


DGAL_END_NAMESPACE
#endif//BOEDER_GEO_FIELD_BUILDER__H