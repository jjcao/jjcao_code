// Copyright (c) 2007-2010  Dalian University of Technology (China). 
// All rights reserved.
//
// This file is part of DGAL; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL: http://jjcao1231.googlepages.com $
// $Id: planar_parameterizer_3.h 2007-04-04 $
//
//
// Author(s)     : JJCAO <jjcao@hotmail.com>

#ifndef DGAL_PLANAR_PARAMETERIZER_3_H
#define DGAL_PLANAR_PARAMETERIZER_3_H

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

DGAL_BEGIN_NAMESPACE

// ------------------------------------------------------------------------------------
// Declaration
// ------------------------------------------------------------------------------------

/// The class Planar_parameterizer_3
/// is a class of planar parameterization methods supporting different border condition and parameterization,
/// different constrains (including hard constrains and no constrains), 
/// different parameterization methods belong to the family of affine combination, and
/// different linear system solver.
/// 
/// Anyway, you can modify the behavior of this class by construcing it with different policies.
/// - Borderpolicy is ...
/// - WeightPolicy is responsible for computing w_ij = (i, j) coefficient of matrix A for j neighbor vertex of i.
/// + WeightPolicy and BorderPolicy together determine the final case of is_one_to_one_mapping().
/// - ConstrainPolicy is responsible for ...
/// + If the hard constrain policy is applied, the linear system will be constructed with hard constrians. 
/// - SolverPolicy is responsible for ...
///
/// @todo Planar_parameterizer_3 should have the ability to add Steiner points for
/// always getting a solution.
///
/// Concept:
/// Model of the CGAL::ParameterizerTraits_3 concept.
///
/// Design Pattern:
/// Planar_parameterizer_3<ParameterizationMesh_3, ...> class is a 
/// copyist of CGAL::Fixed_border_parameterizer_3, but change the original design to 
/// the Policy-Based Design, and this class is a host of four policies.
/// Planar_parameterizer_3<ParameterizationMesh_3, ...> class is a
/// Strategy [GHJV95]: it implements (part of) a strategy of surface parameterization
/// for models of ParameterizationMesh_3.

template
<
    class ParameterizationMesh_3,       ///< 3D surface mesh
    class BorderPolicy         ///< Strategy to parameterize the surface border
		        = CGAL::Circular_border_arc_length_parameterizer_3<ParameterizationMesh_3>,
	class WeightPolicy
		        = DCP_policy<ParameterizationMesh_3>,
	class SolverPolicy   ///< Traits class to solve a sparse linear system
                = OpenNL::DefaultLinearSolverTraits<typename ParameterizationMesh_3::NT>,
	class ConstrainPolicy
		        = No_constrian_policy<ParameterizationMesh_3,
				                      SolverPolicy::Vector,SolverPolicy::Matrix>    
>
class Planar_parameterizer_3
	: public CGAL::Parameterizer_traits_3<ParameterizationMesh_3>,
	         BorderPolicy,
			 WeightPolicy,
			 SolverPolicy,
			 ConstrainPolicy			 
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

// Public operations
public:
	/// Constructor
    Planar_parameterizer_3(){}
	Error_code  parameterize(Adaptor& mesh);

// Protected operations
protected:
    /// Check parameterize() preconditions:
    /// - 'mesh' must be a surface with one connected component.
    /// - 'mesh' must be a triangular mesh.
    /// - the mesh border must be mapped onto a convex polygon.
    Error_code  check_parameterize_preconditions(Adaptor& mesh);

    /// Initialize A, Bu and Bv after border parameterization.
    /// Fill the border vertices' lines in both linear systems:
    /// "u = constant" and "v = constant".
    ///
    /// Preconditions:
    /// - vertices must be indexed.
    /// - A, Bu and Bv must be allocated.
    /// - border vertices must be parameterized.
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

template<class Adaptor, class BorderPolicy, class WeightPolicy,class ConstrainPolicy,class SolverPolicy>
inline
typename CGAL::Parameterizer_traits_3<Adaptor>::Error_code
Planar_parameterizer_3<Adaptor, BorderPolicy, WeightPolicy, ConstrainPolicy, SolverPolicy>::
parameterize(Adaptor& mesh)
{
#ifdef DEBUG_TRACE
    // Create timer for traces
    CGAL::Timer timer;
    timer.start();
#endif

    // Check preconditions
    Error_code status = check_parameterize_preconditions(mesh);
#ifdef DEBUG_TRACE
    std::cerr << "  parameterization preconditions: " << timer.time() << " seconds." << std::endl;
    timer.reset();
#endif
    if (status != Base::OK)
        return status;

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

    // Compute (u,v) for border vertices
    // and mark them as "parameterized"
    status = parameterize_border(mesh);
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
    //
    // So the linear system is not a symmetric definite positive matrix formally 
    // even for Tutte Barycentric Mapping and Discrete Conformal Map algorithms.
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
		
		// Compute the line i of matrix A for i inner vertex with constriant
		if (mesh.is_constrained_vertex(vertexIt))
		{
			setup_constrained_vertex_relations(A, Bu, Bv,
                                              mesh,
                                              vertexIt);
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
    if ( !linear_solver(A, Bu, Xu, Du) ||
         !linear_solver(A, Bv, Xv, Dv) )
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

    // WARNING: this package does not support homogeneous coordinates!
    CGAL_surface_mesh_parameterization_assertion(Du == 1.0);
    CGAL_surface_mesh_parameterization_assertion(Dv == 1.0);

    // Copy Xu and Xv coordinates into the (u,v) pair of each vertex
    set_mesh_uv_from_system (mesh, Xu, Xv);
#ifdef DEBUG_TRACE
    std::cerr << "  copy computed UVs to mesh :"
              << timer.time() << " seconds." << std::endl;
    timer.reset();
#endif


    // Check postconditions
    status = check_parameterize_postconditions(mesh, A, Bu, Bv);
#ifdef DEBUG_TRACE
    std::cerr << "  parameterization postconditions: " << timer.time() << " seconds." << std::endl;
#endif
    if (status != Base::OK)
        return status;

    return status;
}

template<class Adaptor, class BorderPolicy, class WeightPolicy,class ConstrainPolicy,class SolverPolicy>
inline
typename CGAL::Parameterizer_traits_3<Adaptor>::Error_code
Planar_parameterizer_3<Adaptor, BorderPolicy, WeightPolicy, ConstrainPolicy, SolverPolicy>::
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

    // One-to-one mapping is guaranteed if all w_ij coefficients are > 0 (for j vertex neighbor of i)
    // and if the surface border is mapped onto a 2D convex polygon
    CGAL_surface_mesh_parameterization_precondition_code(       \
        status = is_border_convex()  \
               ? Base::OK                                       \
               : Base::ERROR_NON_CONVEX_BORDER;                    \
    );
    if (status != Base::OK)
        return status;

    return status;
}

template<class Adaptor, class BorderPolicy, class WeightPolicy,class ConstrainPolicy,class SolverPolicy>
inline
void Planar_parameterizer_3<Adaptor, BorderPolicy, WeightPolicy, ConstrainPolicy, SolverPolicy>::
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

template<class Adaptor, class BorderPolicy, class WeightPolicy,class ConstrainPolicy,class SolverPolicy>
inline
typename CGAL::Parameterizer_traits_3<Adaptor>::Error_code
Planar_parameterizer_3<Adaptor, BorderPolicy, WeightPolicy, ConstrainPolicy, SolverPolicy>::
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
        NT w_ij = -1.0 * compute_w_ij(mesh, vertex, v_j);

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

template<class Adaptor, class BorderPolicy, class WeightPolicy,class ConstrainPolicy,class SolverPolicy>
inline
void Planar_parameterizer_3<Adaptor, BorderPolicy, WeightPolicy, ConstrainPolicy, SolverPolicy>::
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

template<class Adaptor, class BorderPolicy, class WeightPolicy,class ConstrainPolicy,class SolverPolicy>
inline
typename CGAL::Parameterizer_traits_3<Adaptor>::Error_code
Planar_parameterizer_3<Adaptor, BorderPolicy, WeightPolicy, ConstrainPolicy, SolverPolicy>::
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

template<class Adaptor, class BorderPolicy, class WeightPolicy,class ConstrainPolicy,class SolverPolicy>
inline
bool Planar_parameterizer_3<Adaptor, BorderPolicy, WeightPolicy, ConstrainPolicy, SolverPolicy>::
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

DGAL_END_NAMESPACE

#endif//DGAL_PLANAR_PARAMETERIZER_3_H