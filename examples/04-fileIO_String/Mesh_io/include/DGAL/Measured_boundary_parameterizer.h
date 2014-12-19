#ifndef MEASURED_BOUNDARY_PARAMETERIZER__H
#define MEASURED_BOUNDARY_PARAMETERIZER__H

#pragma warning(disable :  4996 )

#include <DGAL/config.h>

#include <DGAL/Laplacian_kernel.h>
#include <DGAL/Weight_strategy.h>
#include <CGAL/Timer.h>
DGAL_BEGIN_NAMESPACE

/// Poisson distance field builder solves a Poisson system with:
/// constraints: constraints
/// weight: spring
/// right hand: mean spoke length, set constrained line as 0
template<
	class Polyhedron, 
	class Weight_strategy = Weight_spring<Polyhedron>,
	class Solver = Cholmod_solver_traits<Kernel_default::FT>
>
class Poisson_distance_field_builder{
public:
	typedef typename Solver::Matrix    Matrix;
	typedef typename Solver::Vector     Vector;
	typedef Laplacian_field_builder<Polyhedron,Weight_strategy,Solver> PDF_builder;
public:
	Poisson_distance_field_builder(){m_builder.weight_strategy().normailze(false);}
	bool compute(Polyhedron* mesh,std::list<Vertex_handle>& constraints)
	{
		m_builder.compute_laplacian_matrix(mesh, constraints.begin(), constraints.end());

		// compute right hand
		int n(mesh->size_of_vertices());
		Vector b(n),x(n);
		mesh->compute_mean_spokeLen(b);	
		std::for_each(constraints.begin(),constraints.end(),PDF_builder::Set_value(b,0));

		// solve
		if(!m_builder.solve(b,x)){
			std::cout << "solved failed!" << std::endl;
			return false;
		}

		std::for_each(mesh->vertices_begin(),mesh->vertices_end(),PDF_builder::Set_scalar(x));
		return true;
	}

	cholmod_sparse* matrix(){return m_builder.matrix();}
	PDF_builder& builder(){return m_builder;}
protected:
	PDF_builder m_builder;
};

/// Poisson border distance field builder solves a Poisson system with:
/// constraints: boundary vertices;
/// weight: dcp;
/// right hand: mean spoke length, set constrained line as 0
template<
class Polyhedron, 
class Weight_strategy = Weight_DCP<Polyhedron>,//should be DCP
class Solver = Cholmod_solver_traits<Kernel_default::FT>
>
class Poisson_border_field_builder : public Poisson_distance_field_builder<Polyhedron, Weight_strategy, Solver>
{
public:	
	bool compute(Polyhedron* mesh){
		return compute(mesh,mesh->main_border());
	}
	bool compute(Polyhedron* mesh, std::list<Vertex_handle>& cvs)
	{
		m_builder.compute_laplacian_matrix(mesh, cvs.begin(), cvs.end());

		// compute right hand
		int n(mesh->size_of_vertices());
		Vector b(n),x(n);
		mesh->compute_mean_spokeLen(b);		
		std::for_each(mesh->main_border_begin(), mesh->main_border_end(),PDF_builder::Set_value(b,0));

		// solve
		if(!m_builder.solve(b,x)){
			std::cout << "solved failed!" << std::endl;
			return false;
		}

		std::for_each(mesh->vertices_begin(),mesh->vertices_end(),PDF_builder::Set_scalar(x));

		return true;
	}
	void* matrix(){return m_builder.matrix();}
	void* factor(){return m_builder.factor();}
};

template<
class Border_field_builder,
class Distance_field_builder,
class Border_parameterizer,
class Solver
>
class Measured_boundary_parameterizer
{
public:

	typedef typename Solver::Matrix    Matrix;
	typedef typename Solver::Vector     Vector;
	
	//typedef typename Solver::Dense_matrix     Dense_matrix;
	typedef typename Weight_DCP<Polyhedron,Matrix>   Weight_DCP;

	typedef DGAL::Planar_parameterizer<Polyhedron, Border_parameterizer,Weight_DCP, Solver> MB_parameterizer;	
public:

	bool parameterize(Polyhedron* mesh)
	{
		// 1. compute disk-like mesh's center vertex using Border_field_builder
		//    If Border_field_builder is Poisson_border_field_builder, then its matrix and factor would be
		//    reused in step 3 (parameterization) to speed up.
		Border_field_builder bfd;
		if(!bfd.compute(mesh,mesh->main_border())) return false;
		if(center()==0)
		{			
			m_center = compute_center(mesh);		
		}		
		
		std::list<Vertex_handle> constraints(1,m_center);

		// 2. compute distance map from specified center
		Distance_field_builder dfd;
		if(!dfd.compute(mesh, constraints)) return false;		

		// 3. measured boundary parameterization
		// constraints: boundary vertices
		// weight: dcp
		// right hand: zero for right hand, set constrained line as specified
		MB_parameterizer mbp;
		mbp.matrix(static_cast<cholmod_sparse*>(bfd.matrix()));
		mbp.factor(static_cast<cholmod_factor*>(bfd.factor()));
		bool result = mbp.parameterize(mesh,m_center);
		mbp.matrix(0);
		mbp.factor(0);
		return result;
	}
	struct Vertex_with_max_scalar
	{
		typedef Polyhedron::Vertex Vertex;
		Vertex_with_max_scalar(Vertex_handle vh):m_vh(vh){}
		void operator() (Vertex& vh){
			if (vh.s() > m_vh->s())	
			{
				m_vh = &vh;
			}
		}
		operator Vertex_handle(){return m_vh;}
		Vertex_handle m_vh;
	};
	/// compute disk-like mesh's center vertex using distance field saved in vertex.s();	
	Vertex_handle compute_center(Polyhedron* mesh)
	{	
		Vertex_handle vcenter = mesh->vertices_begin();
		vcenter = std::for_each(mesh->vertices_begin(),mesh->vertices_end(),Vertex_with_max_scalar(vcenter));
		return vcenter;		
	}	

	void center(Vertex_handle in){m_center = in;}
	Vertex_handle center(){return m_center;}
	
private:
	Vertex_handle m_center;//computed center vertex
};
DGAL_END_NAMESPACE
#endif //MEASURED_BOUNDARY_PARAMETERIZER__H