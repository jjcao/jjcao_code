#ifndef DGAL_WEIGHT_STRATEGY__H__
#define DGAL_WEIGHT_STRATEGY__H__

#include <DGAL/config.h>
#include <DGAL/Cholmod_solver_traits.h>
#include <DGAL/Util/math.h>

enum WEIGHT_STRATEGY {
	TUTTE  = 0,
	SPRING,
	DCP,
	MVC
};

DGAL_BEGIN_NAMESPACE

//#pragma warning(disable : 4244 4267 4503 4819 4996 )


/// Postconditions:
/// L is normalized;
template<
	class Polyhedron_3,
	class Sparse_matrix
>
class Weight_abstract
{
public:
	typedef typename Polyhedron_3::Point_3 Point_3;
	typedef typename Polyhedron_3::Vector_3 Vector_3;
	typedef typename Polyhedron_3::Vertex_handle Vertex_handle;
	typedef typename Polyhedron_3::Facet_handle Facet_handle;	
	typedef typename Polyhedron_3::Halfedge_around_vertex_circulator Halfedge_vertex_circulator;
	typedef typename Polyhedron_3::Halfedge_around_facet_circulator Halfedge_facet_circulator;
	
	typedef typename Sparse_matrix::NT NT;

	struct Triplet_entry ///< triplet: (i, j, value)
	{
		int i, j; // row and column index
		NT value;
		Triplet_entry() : i(0), j(0), value(NT()) {}
		Triplet_entry(int r, int c, NT val)
			: i(r), j(c), value(val) {}
		void negate(){value = -value;}
	};
	
	struct Fill_matrix{
		Fill_matrix(Sparse_matrix* L, double sum, bool normalize):L(L),normalize(normalize),sum(sum){}
		void operator()(Triplet_entry& entry){
			if(normalize)
				L->add_coef(entry.i,entry.j, (entry.value)/sum);
			else
				L->add_coef(entry.i,entry.j, entry.value);
		}
		Sparse_matrix* L;
		double sum;
		bool normalize;
	};
	//////////////////////////////////////////////////////////////////////////
public:		
	Weight_abstract(WEIGHT_STRATEGY type):m_type(type),m_bNormalize(true){}
	/// Preconditions:
	/// row i has never been set before.
	void set_row(Vertex_handle& vh, Sparse_matrix* L)
	{
		std::list<Triplet_entry> row;//current row
		double sum = compute(row, vh);
		std::for_each(row.begin(), row.end(), Fill_matrix(L,sum,m_bNormalize));			
	}
	/// Preconditions:
	/// row i has never been set before.
	void set_constrained_row(Vertex_handle& vh, Sparse_matrix* L)
	{
		L->add_coef(vh->index(), vh->index(), 1);	
	}
	void normailze(bool in){m_bNormalize=in;}
	bool is_normailized(){return m_bNormalize;} 
	WEIGHT_STRATEGY type(){return m_type;}
protected:
	bool is_border(Vertex_handle vh)
	{
		Halfedge_vertex_circulator pHalfedge = vh->vertex_begin();
		Halfedge_vertex_circulator end = pHalfedge;
		if(pHalfedge == NULL) // isolated vertex
			return true;
		CGAL_For_all(pHalfedge,end)
			if(pHalfedge->is_border())
				return true;
		return false;
	}
	void get_other_vertex(Vertex_handle vi, Facet_handle fh, Vertex_handle& vj, Vertex_handle& vk)
	{
		std::list<Vertex_handle> vhs;
		Halfedge_facet_circulator he = fh->facet_begin();
		Halfedge_facet_circulator he_end = he;
		CGAL_For_all(he, he_end)
		{
			vhs.push_back(he->vertex());
		}
		std::remove(vhs.begin(), vhs.end(), vi);

		std::list<Vertex_handle>::iterator it=vhs.begin();
		vj = *it;
		++it;
		vk = *it;
	}
	virtual double compute(std::list<Triplet_entry>& row, Vertex_handle& vh)=0;
	WEIGHT_STRATEGY m_type;
	bool m_bNormalize;
};
template<
	class Polyhedron_3,
	class Sparse_matrix
>
class Weight_Tutte : public Weight_abstract<Polyhedron_3,Sparse_matrix>
{
public:
	Weight_Tutte():Weight_abstract(TUTTE){}
protected:
	double compute(std::list<Triplet_entry>& row, Vertex_handle& vh)
	{
		double result(0.0);
		int i = vh->index();
		Facet_handle dfh;
		Halfedge_vertex_circulator he  =  vh->vertex_begin();
		Halfedge_vertex_circulator end  =  he;

		if(is_border(vh))
		{
			CGAL_For_all(he,end) 
			{
				Facet_handle fh = he->facet();
				if(fh!=dfh)
				{
					Vertex_handle vj,vk;
					get_other_vertex(vh, fh, vj, vk);
					int j = vj->index(), k=vk->index();

					if(is_border(vj))
					{
						row.push_back(Triplet_entry(i,j,-1));
						result+=1;
					}else{
						row.push_back(Triplet_entry(i,j,-0.5));
						result+=0.5;
					}
					if(is_border(vk))
					{
						row.push_back(Triplet_entry(i,k,-1));
						result+=1;
					}else{
						row.push_back(Triplet_entry(i,k,-0.5));
						result+=0.5;
					}
				}				
			}
		}else{
			CGAL_For_all(he,end) 
			{
				Facet_handle fh = he->facet();
				Vertex_handle vj,vk;
				get_other_vertex(vh, fh, vj, vk);
				int j = vj->index(), k=vk->index();

				row.push_back(Triplet_entry(i,j,-0.5));
				row.push_back(Triplet_entry(i,k,-0.5));
				result+=1;			
			}
		}

		row.push_back(Triplet_entry(i,i,result));

		return result;
	}
};
template<
	class Polyhedron_3,
	class Sparse_matrix
>
class Weight_spring : public Weight_abstract<Polyhedron_3,Sparse_matrix>
{
public:
	Weight_spring():Weight_abstract(SPRING){}
protected:
	double length(Vertex_handle& vi, Vertex_handle& vj)
	{
		Vector_3 v = vi->point() - vj->point();
		return std::sqrt(v*v);
	}
	double compute(std::list<Triplet_entry>& row, Vertex_handle& vh)
	{
		double result(0.0);
		int i = vh->index();
		Facet_handle dfh;
		Halfedge_vertex_circulator he  =  vh->vertex_begin();
		Halfedge_vertex_circulator end  =  he;

		if(is_border(vh))
		{
			double len(0.0);
			CGAL_For_all(he,end) 
			{
				Facet_handle fh = he->facet();
				if(fh!=dfh)
				{
					Vertex_handle vj,vk;
					get_other_vertex(vh, fh, vj, vk);
					int j = vj->index(), k=vk->index();

					len = 1.0/length(vh, vj);
					if(is_border(vj))
					{						
						row.push_back(Triplet_entry(i,j,-len));						
						result+=len;
					}else{
						row.push_back(Triplet_entry(i,j,-0.5*len));
						result+=0.5*len;
					}
					
					len = length(vh, vk);
					if(is_border(vk))
					{
						row.push_back(Triplet_entry(i,k,-len));	
						result+=len;
					}else{
						row.push_back(Triplet_entry(i,k,-0.5*len));
						result+=0.5*len;
					}
				}				
			}
		}else{
			double lenj(0.0), lenk(0);
			CGAL_For_all(he,end) 
			{
				Facet_handle fh = he->facet();
				Vertex_handle vj,vk;
				get_other_vertex(vh, fh, vj, vk);
				int j = vj->index(), k=vk->index();

				lenj = 0.5/length(vh, vj); lenk = 0.5/length(vh, vk);
				row.push_back(Triplet_entry(i,j,-lenj));
				row.push_back(Triplet_entry(i,k,-lenk));
				result+=lenj;result+=lenk;
			}
		}

		row.push_back(Triplet_entry(i,i,result));

		return result;
	}
};

/// todo
template<
	class Polyhedron_3,
	class Sparse_matrix
>
class Weight_DCP : public Weight_abstract<Polyhedron_3,Sparse_matrix>
{
public:
	typedef typename Polyhedron_3::Traits Kernel;
	typedef typename DGAL::Math<Kernel> Math;
	Weight_DCP():Weight_abstract(DCP){}
protected:
	double length(Vertex_handle& vi, Vertex_handle& vj)
	{
		Vector_3 v = vi->point() - vj->point();
		return std::sqrt(v*v);
	}
	double compute(std::list<Triplet_entry>& row, Vertex_handle& vh)
	{
		double result(0.0);
		int i = vh->index();
		Facet_handle dfh;
		Halfedge_vertex_circulator he  =  vh->vertex_begin();
		Halfedge_vertex_circulator end  =  he;

		CGAL_For_all(he,end) 
		{
			Facet_handle fh = he->facet();
			if(fh!=NULL)
			{
				Vertex_handle vj,vk;
				get_other_vertex(vh, fh, vj, vk);				
				
				double cotg_alpha_ij = Math::cotangent(vj->point(),vk->point(),vh->point());//0.5*
				double cotg_alpha_ik = Math::cotangent(vh->point(),vj->point(),vk->point());//0.5*
				
				int j = vj->index(), k=vk->index();
				row.push_back(Triplet_entry(i,j,-cotg_alpha_ij));
				row.push_back(Triplet_entry(i,k,-cotg_alpha_ik));
				result+=cotg_alpha_ij;result+=cotg_alpha_ik;
			}
		}

		row.push_back(Triplet_entry(i,i,result));

		return result;
	}
};

DGAL_END_NAMESPACE
#endif//DGAL_WEIGHT_STRATEGY__H__