#ifndef DGAL_MEASURED_BOUNDARY__H
#define DGAL_MEASURED_BOUNDARY__H

#include <CGAL/circulator.h>
#include <CGAL/Timer.h>
#include <CGAL/OpenNL/linear_solver.h>

#include <CGAL/Parameterizer_traits_3.h>
#include <CGAL/Circular_border_parameterizer_3.h>
#include <CGAL/Parameterization_mesh_feature_extractor.h>
#include <CGAL/surface_mesh_parameterization_assertions.h>

#include <DGAL/config.h>
#include <DGAL/Parameterization_policies_3.h>
#include <DGAL/Parameterization_constrain_policies_3.h>
#include <iostream>
#include <numeric>
DGAL_BEGIN_NAMESPACE

template<
	class ParameterizationMesh_3,
	class WeightPolicy  	        = DCP_policy<ParameterizationMesh_3>,
	class SolverPolicy              = OpenNL::DefaultLinearSolverTraits<typename ParameterizationMesh_3::NT>
>
class MB : public CGAL::Parameterizer_traits_3<ParameterizationMesh_3>
{
	// Private types
private:
	// Superclass
	typedef Parameterizer_traits_3<ParameterizationMesh_3>
		Base;

	// Public types
public:
	// We have to repeat the types exported by superclass
	typedef typename Base::Error_code       Error_code;
	typedef ParameterizationMesh_3          Adaptor;
	// Private types
private:
	typedef typename SolverPolicy           Solver;
	// Mesh_Adaptor_3 subtypes:
	typedef typename Adaptor::NT            NT;
	typedef typename Adaptor::Point_2       Point_2;
	typedef typename Adaptor::Point_3       Point_3;
	typedef typename Adaptor::Vector_2      Vector_2;
	typedef typename Adaptor::Vector_3      Vector_3;
	typedef typename Adaptor::Facet         Facet;
	typedef typename Adaptor::Facet_handle  Facet_handle;
	typedef typename Adaptor::Facet_const_handle
		Facet_const_handle;
	typedef typename Adaptor::Facet_iterator Facet_iterator;
	typedef typename Adaptor::Facet_const_iterator
		Facet_const_iterator;
	typedef typename Adaptor::Vertex        Vertex;
	typedef typename Adaptor::Vertex_handle Vertex_handle;
	typedef typename Adaptor::Vertex_const_handle
		Vertex_const_handle;
	typedef typename Adaptor::Vertex_iterator Vertex_iterator;
	typedef typename Adaptor::Vertex_const_iterator
		Vertex_const_iterator;
	typedef typename Adaptor::Border_vertex_iterator
		Border_vertex_iterator;
	typedef typename Adaptor::Border_vertex_const_iterator
		Border_vertex_const_iterator;
	typedef typename Adaptor::Vertex_around_facet_circulator
		Vertex_around_facet_circulator;
	typedef typename Adaptor::Vertex_around_facet_const_circulator
		Vertex_around_facet_const_circulator;
	typedef typename Adaptor::Vertex_around_vertex_circulator
		Vertex_around_vertex_circulator;
	typedef typename Adaptor::Vertex_around_vertex_const_circulator
		Vertex_around_vertex_const_circulator;

	// SolverPolicy subtypes:
	typedef typename SolverPolicy::Vector      Vector;
	typedef typename SolverPolicy::Matrix      Matrix;

	// Private types
private:
	// Mesh_Adaptor_3 subtypes:
	typedef typename Adaptor::NT            NT;

public:
	MB(){}
	Error_code  parameterize(Adaptor& mesh);
private:
	WeightPolicy weightPolicy;
	// Protected operations
protected:
	/// Check parameterize() preconditions:
	/// - 'mesh' must be a surface with one connected component.
	/// - 'mesh' must be a triangular mesh.
	/// - the mesh border must be mapped onto a convex polygon.
	Error_code  check_parameterize_preconditions(Adaptor& mesh);

	/// set u=v=0 for the start boundary vertex and mark them as "parameterized"
	Error_code parameterize_border(Adaptor& mesh,Vertex_handle startVertex);

	/// Initialize A, Bu and Bv after border parameterization.
	/// Fill the first border vertex's lines in both linear systems:
	/// "u = constant" and "v = constant".
	///
	/// Preconditions:
	/// - vertices must be indexed.
	/// - A, Bu and Bv must be allocated.
	/// - the first border vertex must be parameterized.
	void  initialize_system_from_mesh_border (Matrix& A, Vector& Bu, Vector& Bv,
		const Adaptor& mesh);

	/// Compute the line i of matrix A for i inner vertex:
	/// - call compute_w_ij() to compute the A coefficient w_ij for each neighbor v_j.
	/// - compute w_ii = - sum of w_ijs.
	///
	/// Preconditions:
	/// - vertices must be indexed.
	/// - vertex i musn't be already parameterized.
	/// - line i of A must contain only zeros.
	Error_code setup_inner_vertex_relations(Matrix& A,
		Vector& Bu,
		Vector& Bv,
		Adaptor& mesh,
		Vertex_handle vertex);

	/// Copy Xu and Xv coordinates into the (u,v) pair of each surface vertex.
	void  set_mesh_uv_from_system (Adaptor& mesh,
		const Vector& Xu, const Vector& Xv);

	/// Check parameterize() postconditions:
	/// - 3D -> 2D mapping is one-to-one.
	Error_code check_parameterize_postconditions(const Adaptor& mesh,
		const Matrix& A,
		const Vector& Bu,
		const Vector& Bv);

	/// Check if 3D -> 2D mapping is one-to-one.
	/// The default implementation checks each normal.
	/// todo It should be optimized for different weight policy and border policy
	bool  is_one_to_one_mapping(const Adaptor& mesh,
		const Matrix& A,
		const Vector& Bu,
		const Vector& Bv);
};

template<class Adaptor, class WeightPolicy,class SolverPolicy>
inline
typename CGAL::Parameterizer_traits_3<Adaptor>::Error_code
MB<Adaptor, WeightPolicy, SolverPolicy>::
check_parameterize_preconditions(Adaptor& mesh)
{
	Error_code status = Base::OK;	    // returned value

	// Helper class to compute genus or count borders, vertices, ...
	typedef CGAL::Parameterization_mesh_feature_extractor<Adaptor>
		Mesh_feature_extractor;
	Mesh_feature_extractor feature_extractor(mesh);

	// Allways check that mesh is not empty
	if (mesh.mesh_vertices_begin() == mesh.mesh_vertices_end())
		status = Base::ERROR_EMPTY_MESH;
	if (status != Base::OK)
		return status;

	// The whole surface parameterization package is restricted to triangular meshes
	CGAL_surface_mesh_parameterization_expensive_precondition_code(            \
		status = mesh.is_mesh_triangular() ? Base::OK                         \
		: Base::ERROR_NON_TRIANGULAR_MESH; \
		);
	if (status != Base::OK)
		return status;

	// The whole package is restricted to surfaces: genus = 0,
	// one connected component and at least one border
	CGAL_surface_mesh_parameterization_expensive_precondition_code(         \
		int genus = feature_extractor.get_genus();                          \
		int nb_borders = feature_extractor.get_nb_borders();                \
		int nb_components = feature_extractor.get_nb_connex_components();   \
		status = (genus == 0 && nb_borders >= 1 && nb_components == 1)      \
		? Base::OK                                                   \
		: Base::ERROR_NO_SURFACE_MESH;                               \
		);
	if (status != Base::OK)
		return status;

	return status;
}

template<class Adaptor, class WeightPolicy,class SolverPolicy>
inline
void MB<Adaptor, WeightPolicy, SolverPolicy>::
initialize_system_from_mesh_border (Matrix& A, Vector& Bu, Vector& Bv,
									const Adaptor& mesh)
{
	for (Border_vertex_const_iterator it = mesh.mesh_main_border_vertices_begin();
		it != mesh.mesh_main_border_vertices_end();
		it++)
	{
		CGAL_surface_mesh_parameterization_assertion(mesh.is_vertex_parameterized(it));

		// Get vertex index in sparse linear system
		int index = mesh.get_vertex_index(it);

		// Write a as diagonal coefficient of A
		A.set_coef(index, index, 1);

		// Write constant in Bu and Bv
		Point_2 uv = mesh.get_vertex_uv(it);
		Bu[index] = uv.x();
		Bv[index] = uv.y();
	}
}

template<class Adaptor, class WeightPolicy,class SolverPolicy>
inline
typename CGAL::Parameterizer_traits_3<Adaptor>::Error_code
MB<Adaptor, WeightPolicy, SolverPolicy>::
setup_inner_vertex_relations(Matrix& A,
							 Vector& Bu,
							 Vector& Bv,
							 Adaptor& mesh,
							 Vertex_handle vertex)
{
	CGAL_surface_mesh_parameterization_assertion( ! mesh.is_vertex_on_main_border(vertex) );
	CGAL_surface_mesh_parameterization_assertion( ! mesh.is_vertex_parameterized(vertex) );

	int i = mesh.get_vertex_index(vertex);

	// circulate over vertices around 'vertex' to compute w_ii and w_ijs
	NT w_ii = 0;
	int vertexIndex = 0;
	Vertex_around_vertex_circulator v_j = mesh.vertices_around_vertex_begin(vertex);
	Vertex_around_vertex_circulator end = v_j;
	CGAL_For_all(v_j, end)
	{
		// Call to virtual method to do the actual coefficient computation
		NT w_ij = -1.0 * weightPolicy.compute_w_ij(mesh, vertex, v_j);

		// w_ii = - sum of w_ijs
		w_ii -= w_ij;

		// Get j index
		int j = mesh.get_vertex_index(v_j);

		// Set w_ij in matrix
		A.set_coef(i,j, w_ij);

		vertexIndex++;
	}
	if (vertexIndex < 2)
		return Base::ERROR_NON_TRIANGULAR_MESH;

	// Set w_ii in matrix
	A.set_coef(i,i, w_ii);

	return Base::OK;
}

template<class Adaptor, class WeightPolicy,class SolverPolicy>
inline
void MB<Adaptor, WeightPolicy, SolverPolicy>::
set_mesh_uv_from_system(Adaptor& mesh,
						const Vector& Xu, const Vector& Xv)
{
	Vertex_iterator vertexIt;
	for (vertexIt = mesh.mesh_vertices_begin();
		vertexIt != mesh.mesh_vertices_end();
		vertexIt++)
	{
		int index = mesh.get_vertex_index(vertexIt);

		NT u = Xu[index];
		NT v = Xv[index];

		// Fill vertex (u,v) and mark it as "parameterized"
		mesh.set_vertex_uv(vertexIt, Point_2(u,v));
		mesh.set_vertex_parameterized(vertexIt, true);
	}
}
template<class Adaptor, class WeightPolicy,class SolverPolicy>
inline
typename CGAL::Parameterizer_traits_3<Adaptor>::Error_code
MB<Adaptor, WeightPolicy, SolverPolicy>::
check_parameterize_postconditions(const Adaptor& mesh,
								  const Matrix& A,
								  const Vector& Bu,
								  const Vector& Bv)
{
	Error_code status = Base::OK;

	// Check if 3D -> 2D mapping is one-to-one
	CGAL_surface_mesh_parameterization_postcondition_code(  \
		status = is_one_to_one_mapping(mesh, A, Bu, Bv)     \
		? Base::OK                                   \
		: Base::ERROR_NO_1_TO_1_MAPPING;             \
		);
	if (status != Base::OK)
		return status;

	return status;
}

template<class Adaptor, class WeightPolicy,class SolverPolicy>
inline
bool MB<Adaptor, WeightPolicy, SolverPolicy>::
is_one_to_one_mapping(const Adaptor& mesh,
					  const Matrix& A,
					  const Vector& Bu,
					  const Vector& Bv)
{
	Vector_3    first_triangle_normal;

	for (Facet_const_iterator facetIt = mesh.mesh_facets_begin();
		facetIt != mesh.mesh_facets_end();
		facetIt++)
	{
		// Get 3 vertices of the facet
		Vertex_const_handle v0, v1, v2;
		int vertexIndex = 0;
		Vertex_around_facet_const_circulator cir = mesh.facet_vertices_begin(facetIt),
			end = cir;
		CGAL_For_all(cir, end)
		{
			if (vertexIndex == 0)
				v0 = cir;
			else if (vertexIndex == 1)
				v1 = cir;
			else if (vertexIndex == 2)
				v2 = cir;

			vertexIndex++;
		}
		CGAL_surface_mesh_parameterization_assertion(vertexIndex >= 3);

		// Get the 3 vertices position IN 2D
		Point_2 p0 = mesh.get_vertex_uv(v0) ;
		Point_2 p1 = mesh.get_vertex_uv(v1) ;
		Point_2 p2 = mesh.get_vertex_uv(v2) ;

		// Compute the facet normal
		Point_3 p0_3D(p0.x(), p0.y(), 0);
		Point_3 p1_3D(p1.x(), p1.y(), 0);
		Point_3 p2_3D(p2.x(), p2.y(), 0);
		Vector_3 v01_3D = p1_3D - p0_3D;
		Vector_3 v02_3D = p2_3D - p0_3D;
		Vector_3 normal = CGAL::cross_product(v01_3D, v02_3D);

		// Check that all normals are oriented the same way
		// => no 2D triangle is flipped
		if (cir == mesh.facet_vertices_begin(facetIt))
		{
			first_triangle_normal = normal;
		}
		else
		{
			if (first_triangle_normal * normal < 0)
				return false;
		}
	}

	return true;            // OK if we reach this point
}
template<class Adaptor, class WeightPolicy,class SolverPolicy>
inline
typename CGAL::Parameterizer_traits_3<Adaptor>::Error_code
MB<Adaptor, WeightPolicy, SolverPolicy>::
parameterize_border(Adaptor& mesh,Vertex_handle startVertex)
{
	// Nothing to do if no border
	if (mesh.mesh_main_border_vertices_begin() == mesh.mesh_main_border_vertices_end())
		return Parameterizer_traits_3<Adaptor>::ERROR_BORDER_TOO_SHORT;

	Point_3 p0 = startVertex->point();

	std::list<double> angles;
	Border_vertex_iterator it1 = mesh.mesh_main_border_vertices_begin();
	Border_vertex_iterator it2 = mesh.mesh_main_border_vertices_begin();
	++it2;

	for (; it2 != mesh.mesh_main_border_vertices_end();++it1,++it2)
	{       
		Point_3 p1 = it1->point();
		Point_3 p2 = it2->point();

		angles.push_back(weightPolicy.compute_angle_rad(p1,p0,p2));		
	}
	it2 = mesh.mesh_main_border_vertices_begin();
	angles.push_back(weightPolicy.compute_angle_rad(it1->point(),p0,it2->point()));	

	double sumSpaceAngle = std::accumulate(angles.begin(), angles.end(), 0.0);
	static const double PI = 3.14159265359;
	double tmp = 2*PI/sumSpaceAngle;
	double currAngle(0.0);
	std::list<double>::iterator it = angles.begin();
	for (; it!=angles.end(); ++it)
	{
		*it = std::abs(*it) * tmp;
		currAngle += *it;
		*it = currAngle;
	}

	it1 = mesh.mesh_main_border_vertices_begin();	
	Vertex_handle vh = it1;
	double radius1 = it1->s();
	mesh.set_vertex_uv(vh, Point_2(1.0,0.5)); 
	mesh.set_vertex_parameterized(vh, true);
	++it1;

	for (it = angles.begin();it1 != mesh.mesh_main_border_vertices_end();++it1,++it)
	{       
		vh = it1;
		double r = vh->s()/(radius1*2);
		double u = 0.5 + std::cos(*it) * r;
		double v = 0.5 + std::sin(*it) * r;
		mesh.set_vertex_uv(vh, Point_2(u,v)); 
		mesh.set_vertex_parameterized(vh, true);
	}	

#ifdef DEBUG_TRACE
	std::cerr << "  map one vertex..." << std::endl;
	// std::cerr << "    #" << mesh.get_vertex_index(vh) << "(" << vh->vertex()->index() << ") parameterized " << std::endl;
#endif

	return Parameterizer_traits_3<Adaptor>::OK;
}
template<class Adaptor, class WeightPolicy,class SolverPolicy>
inline
typename CGAL::Parameterizer_traits_3<Adaptor>::Error_code
MB<Adaptor, WeightPolicy, SolverPolicy>::
parameterize(Adaptor& mesh)
{
	Error_code status = CGAL::Parameterizer_traits_3<Adaptor>::OK;
#ifdef DEBUG_TRACE
	// Create timer for traces
	CGAL::Timer timer;
	timer.start();
	// Check preconditions
	status = check_parameterize_preconditions(mesh);
	std::cerr << "  parameterization preconditions: " << timer.time() << " seconds." << std::endl;
	timer.reset();
	if (status != Base::OK)
		return status;
#endif

	// Count vertices
	int nbVertices = mesh.count_mesh_vertices();

	// Index vertices from 0 to nbVertices-1
	mesh.index_mesh_vertices();

	// Mark all vertices as NOT "parameterized"
	Vertex_iterator vertexIt;
	for (vertexIt = mesh.mesh_vertices_begin();
		vertexIt != mesh.mesh_vertices_end();
		vertexIt++)
	{
		mesh.set_vertex_parameterized(vertexIt, false);
	}

	//
	std::list<Vertex_handle>::iterator it = mesh.constrained_vertices().begin();
	// set u=v=0 for the start boundary vertex and mark them as "parameterized"
	status = parameterize_border(mesh, *it);
#ifdef DEBUG_TRACE
	std::cerr << "  border vertices parameterization: " << timer.time() << " seconds." << std::endl;
	timer.reset();
#endif
	if (status != Base::OK)
		return status;

	// Create two sparse linear systems "A*Xu = Bu" and "A*Xv = Bv" (one line/column per vertex)
	Matrix A(nbVertices, nbVertices);
	Vector Xu(nbVertices), Xv(nbVertices), Bu(nbVertices), Bv(nbVertices);

	// Initialize A, Xu, Xv, Bu and Bv after border parameterization
	// Fill the border vertices' lines in both linear systems:
	// "u = constant" and "v = constant"
	initialize_system_from_mesh_border (A, Bu, Bv, mesh);

	// Fill the matrix for the inner vertices v_i: compute A's coefficient
	// w_ij for each neighbor j; then w_ii = - sum of w_ijs
	for (vertexIt = mesh.mesh_vertices_begin();
		vertexIt != mesh.mesh_vertices_end();
		vertexIt++)
	{		
		CGAL_surface_mesh_parameterization_assertion(mesh.is_vertex_on_main_border(vertexIt)
			== mesh.is_vertex_parameterized(vertexIt));

		// inner vertices only
		if ( mesh.is_vertex_on_main_border(vertexIt))
		{
			continue;
		}		

		// Compute the line i of matrix A for i inner vertex
		status = setup_inner_vertex_relations(A, Bu, Bv,
			mesh,
			vertexIt);
		if (status != Base::OK)
			return status;
	}
#ifdef DEBUG_TRACE
	std::cerr << "  matrix filling (" << nbVertices << " x " << nbVertices << "): "
		<< timer.time() << " seconds." << std::endl;
	timer.reset();
#endif

	// Solve "A*Xu = Bu". On success, solution is (1/Du) * Xu.
	// Solve "A*Xv = Bv". On success, solution is (1/Dv) * Xv.
	NT Du, Dv;
	Solver solver = Solver();
	if ( !solver.linear_solver(A, Bu, Xu, Du) ||
		!solver.linear_solver(A, Bv, Xv, Dv) )
	{
		status = Base::ERROR_CANNOT_SOLVE_LINEAR_SYSTEM;
	}
#ifdef DEBUG_TRACE
	std::cerr << "  solving two linear systems: "
		<< timer.time() << " seconds." << std::endl;
	timer.reset();
#endif
	if (status != Base::OK)
		return status;

	// Copy Xu and Xv coordinates into the (u,v) pair of each vertex
	set_mesh_uv_from_system (mesh, Xu, Xv);
#ifdef DEBUG_TRACE
	std::cerr << "  copy computed UVs to mesh :"
		<< timer.time() << " seconds." << std::endl;
	timer.reset();
	// Check postconditions
	status = check_parameterize_postconditions(mesh, A, Bu, Bv);
	std::cerr << "  parameterization postconditions: " << timer.time() << " seconds." << std::endl;
	if (status != Base::OK)
		return status;
#endif
	return status;
}
DGAL_END_NAMESPACE
#endif//DGAL_MEASURED_BOUNDARY__H