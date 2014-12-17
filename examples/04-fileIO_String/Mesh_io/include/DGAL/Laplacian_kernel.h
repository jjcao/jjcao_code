#ifndef DGAL_LAPLACIAN_KERNEL__H__
#define DGAL_LAPLACIAN_KERNEL__H__

#include <DGAL/config.h>
#include <DGAL/Polyhedron_3.h>
#include <DGAL/Cholmod_solver_traits.h>
#include <DGAL/Weight_strategy.h>
#include <DGAL/Border_parameterizer.h>
DGAL_BEGIN_NAMESPACE

#pragma warning(disable : 4244 4267 4503 4819 4996 )

template<
	class Polyhedron_3,	
	class Weight_strategy,
	class Solver
>
class Laplacian_kernel
{
public:
	typedef typename Polyhedron_3                    Polyhedron;
	typedef typename Polyhedron::Vertex_handle       Vertex_handle;
	typedef typename Solver::Matrix    Matrix;
	typedef typename Solver::Vector     Vector;	
	typedef typename Polyhedron::Vertex              Vertex;
	typedef typename Polyhedron::Vertex_handle       Vertex_handle;
	struct Set_row
	{
		Set_row(Weight_strategy& ws,Matrix* matrix):matrix(matrix),ws(ws){}
		void operator() (Vertex& vb){
			Vertex_handle vh = &vb;
			if (vh->tag())
			{
				ws.set_constrained_row(vh, matrix);
			}
			else
			{
				ws.set_row(vh, matrix);
			}		
		}
		Weight_strategy& ws;
		Matrix* matrix;
	};
public:
	Laplacian_kernel():m_matrix(0){}
	virtual ~Laplacian_kernel(){
		if(m_matrix) delete m_matrix;
	}
	template<class Const_iterator>
	void compute_laplacian_matrix(Polyhedron* mesh, Const_iterator first, Const_iterator last)
	{
		const int tag_free = 0;
		const int tag_done = 1;

		//////////////////////////////////////////////////////////////////////////
		if (m_matrix)
			delete m_matrix;

#ifdef OUTPUT_INFO
		CGAL::Timer timer;	timer.start();
#endif		
		m_matrix = new Matrix(mesh->size_of_vertices());
		m_constraints.clear();
		std::copy(first,last, std::back_inserter(m_constraints));

		// Tag all vertices as unprocessed
		mesh->tag_vertices(tag_free);
		// tag constrained vertices as processed
		std::for_each(m_constraints.begin(),m_constraints.end(),Polyhedron::Set_tag(tag_done));

		//setup weight both for constrained and free vertices
		std::for_each(mesh->vertices_begin(),mesh->vertices_end(),Set_row(m_weight_strategy, m_matrix));
#ifdef OUTPUT_INFO
		std::cout << "matrix setup: " << timer.time() << " seconds." << std::endl;
		timer.reset();
#endif
	}
	//void set_constraints(){}
	void factor(cholmod_factor* factor){m_solver.factor(factor);}
	void matrix(cholmod_sparse* matrix){m_solver.matrix(matrix);}
	cholmod_factor* factor(){return m_solver.factor();}
	cholmod_sparse* matrix(){return m_solver.matrix();}
	Weight_strategy& weight_strategy(){return m_weight_strategy;}

protected:
	Weight_strategy m_weight_strategy;
	Matrix* m_matrix;
	std::list<Vertex_handle> m_constraints;
	Solver m_solver;
};

template<
	class Polyhedron_3,
	class Weight_strategy,
	class Solver
>
class Laplacian_field_builder : public Laplacian_kernel<Polyhedron_3, Weight_strategy, Solver>
{
public:
	struct Set_value
	{
		Set_value(Vector& b, double v):b(b),v(v){}
		void operator()(Vertex_handle vh)
		{
			b[vh->index()] = v;
		}
		Vector& b;
		double v;
	};
	struct Set_scalar
	{
		Set_scalar(Vector& s):scalars(s){}
		void operator() (Vertex& vh){
			double s = scalars[vh.index()];
			if (abs(s)<1.0e-015)	
				s = 0;
			vh.s(s);
		}
		void operator()(Vertex_handle& elem){
			this->operator ()(*elem);
		}
		Vector& scalars;
	};

public:	
	template <class Dense_matrix_T>
	bool solve(Dense_matrix_T& B, Dense_matrix_T& X, cholmod_factor* factor=0)
	{
		if(factor) 
			return m_solver.solve(B, X, factor);
		else
			return m_solver.linear_solver(*m_matrix, B, X, factor);
	}
};

template<
	class Polyhedron_3,
	class Border_strategy,
	class Weight_strategy,
	class Solver
>
class Planar_parameterizer : public Laplacian_kernel<Polyhedron_3, Weight_strategy,Solver>
{
	typedef typename Laplacian_kernel<Polyhedron_3, Weight_strategy,Solver> Base;
	typedef typename Polyhedron_3::NT NT;
	typedef typename Cholmod_dense_matrix<NT> Dense_matrix;
public:
	struct Set_uv2matrix
	{
		Set_uv2matrix(Dense_matrix& b):b(b){}
		void operator()(Vertex_handle vh)
		{
			b.set(vh->index(),0,vh->u());
			b.set(vh->index(),1,vh->v());
		}
		Dense_matrix& b;
	};
	struct Set_uv2mesh
	{
		Set_uv2mesh(Dense_matrix& x):x(x){}
		void operator()(Vertex& vh)
		{
			vh.u(x.get(vh.index(),0));
			vh.v(x.get(vh.index(),1));
		}
		Dense_matrix& x;
	};

	bool parameterize(Polyhedron* mesh, Vertex_handle v1 = NULL){
		compute_laplacian_matrix(mesh, mesh->main_border_begin(), mesh->main_border_end());

		Border_strategy bs;
		bs.parameterize_border(mesh->main_border_begin(), mesh->main_border_end(),v1);

		// compute right hand
		Dense_matrix b(mesh->size_of_vertices(),2);
		Dense_matrix x(mesh->size_of_vertices(),2);
		std::for_each(mesh->main_border_begin(), mesh->main_border_end(),Set_uv2matrix(b));

		// solve
		if(!m_solver.linear_solver(*m_matrix,b,x))
		{
			std::cout << "solved failed!" << std::endl;
			return false;
		}

		std::for_each(mesh->vertices_begin(), mesh->vertices_end(),Set_uv2mesh(x));
		return true;
	}
	void parameterize_border(){}
};

DGAL_END_NAMESPACE
#endif //DGAL_LAPLACIAN_KERNEL__H__