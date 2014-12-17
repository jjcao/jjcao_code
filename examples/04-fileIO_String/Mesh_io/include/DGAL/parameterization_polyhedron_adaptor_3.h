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
// $Id: Parameterization_polyhedron_adaptor_3.h 2007-04-04 $
//
//
// Author(s)     : JJCAO <jjcao@hotmail.com>

#ifndef DGAL_PARAMETERIZATION_POLYHEDRON_ADAPTOR_H
#define DGAL_PARAMETERIZATION_POLYHEDRON_ADAPTOR_H

#include <DGAL/config.h>
#include <CGAL/Parameterization_polyhedron_adaptor_3.h>
#include <DGAL/Parameterization_constrain_policies_3.h>

DGAL_BEGIN_NAMESPACE

template
<
	class Polyhedron_3_
>
class Parameterization_polyhedron_adaptor_3
{

// Forward references
protected:
    struct                                           Project_halfedge_vertex;
    struct                                           Project_vertex_handle_vertex;
    struct                                           Project_opposite_halfedge_vertex;

// Private types
public:
	typedef typename Polyhedron_3_                   Polyhedron;
    // Halfedge
    typedef typename Polyhedron::Halfedge            Halfedge;
    typedef typename Polyhedron::Halfedge_handle     Halfedge_handle;
    typedef typename Polyhedron::Halfedge_const_handle
                                                     Halfedge_const_handle;
    typedef typename Polyhedron::Halfedge_iterator   Halfedge_iterator;
    typedef typename Polyhedron::Halfedge_const_iterator
                                                     Halfedge_const_iterator;
    typedef typename Polyhedron::Halfedge_around_vertex_circulator
                                                     Halfedge_around_vertex_circulator;
    typedef typename Polyhedron::Halfedge_around_vertex_const_circulator
                                                     Halfedge_around_vertex_const_circulator;
    typedef typename Polyhedron::Halfedge_around_facet_circulator
                                                     Halfedge_around_facet_circulator;
    typedef typename Polyhedron::Halfedge_around_facet_const_circulator
                                                     Halfedge_around_facet_const_circulator;

// Public types
public:

    //******************************************************************
    // LEVEL 1 INTERFACE:
    // for "normal" classes that do not deal with virtual seams
    // Example: all parameterization methods
    //******************************************************************

    // Number type
    typedef typename Polyhedron::Traits::FT          NT;
    // Points and vectors
    typedef typename Polyhedron::Traits::Point_2     Point_2;
    typedef typename Polyhedron::Traits::Point_3     Point_3;
    typedef typename Polyhedron::Traits::Vector_2    Vector_2;
    typedef typename Polyhedron::Traits::Vector_3    Vector_3;

    // Facet
    typedef typename Polyhedron::Facet               Facet;
    typedef typename Polyhedron::Facet_handle        Facet_handle;
    typedef typename Polyhedron::Facet_const_handle  Facet_const_handle;
    // Iterator over all mesh facets
    typedef typename Polyhedron::Facet_iterator      Facet_iterator;
    typedef typename Polyhedron::Facet_const_iterator
                                                     Facet_const_iterator;

    // Vertex
    typedef typename Polyhedron::Vertex              Vertex;
    typedef typename Polyhedron::Vertex_handle       Vertex_handle;
    typedef typename Polyhedron::Vertex_const_handle Vertex_const_handle;
    // Iterator over all mesh vertices
    typedef typename Polyhedron::Vertex_iterator     Vertex_iterator;
    typedef typename Polyhedron::Vertex_const_iterator
                                                     Vertex_const_iterator;

    // Iterator over mesh border vertices
    typedef CGAL::Convertible_iterator_project<typename std::list<Vertex_handle>::iterator,
                                               Project_vertex_handle_vertex,
                                               Vertex_const_handle,
                                               Vertex_handle>
                                                     Border_vertex_iterator;
    typedef CGAL::Convertible_iterator_project<typename std::list<Vertex_handle>::const_iterator,
                                               Project_vertex_handle_vertex,
                                               Vertex_const_handle>
                                                     Border_vertex_const_iterator;
    // Counter-clockwise circulator over a facet's vertices
    typedef CGAL::Convertible_circulator_project<typename Halfedge_around_facet_circulator,
                                                 Project_halfedge_vertex,
                                                 Vertex&,
                                                 Vertex*,
                                                 Vertex_const_handle,
                                                 Vertex_handle>
                                                     Vertex_around_facet_circulator;
    typedef CGAL::Convertible_circulator_project<typename Halfedge_around_facet_const_circulator,
                                                 Project_halfedge_vertex,
                                                 const Vertex&,
                                                 const Vertex*,
                                                 Vertex_const_handle>
                                                     Vertex_around_facet_const_circulator;
    // Clockwise circulator over the vertices incident to a vertex
    typedef CGAL::Convertible_circulator_project<typename Halfedge_around_vertex_circulator,
                                                 Project_opposite_halfedge_vertex,
                                                 Vertex&,
                                                 Vertex*,
                                                 Vertex_const_handle,
                                                 Vertex_handle>
                                                     Vertex_around_vertex_circulator;
    typedef CGAL::Convertible_circulator_project<typename Halfedge_around_vertex_const_circulator,
                                                 Project_opposite_halfedge_vertex,
                                                 const Vertex&,
                                                 const Vertex*,
                                                 Vertex_const_handle>
                                                     Vertex_around_vertex_const_circulator;

//Public operations
public:
	Parameterization_polyhedron_adaptor_3(Polyhedron& mesh)
	: m_polyhedron(mesh)
	{
		// Extract mesh's longest border
        m_main_border = extract_longest_border(mesh);

#ifndef NDEBUG
        // Index vertices right away to ease debugging
        index_mesh_vertices();
#endif
	}

	// Get the adapted mesh
    Polyhedron&       get_adapted_mesh()       { return m_polyhedron; }
    const Polyhedron& get_adapted_mesh() const { return m_polyhedron; }

	virtual bool  is_constrained_vertex(Vertex_const_handle vertex) const 
	{
		return false;
    }

    // Get halfedge from source and target vertices
    // Will assert if such an halfedge doesn't exist
    Halfedge_const_handle get_halfedge(
        Vertex_const_handle source, Vertex_const_handle target) const
    {
        assert(source != NULL);
        assert(target != NULL);

        Halfedge_around_vertex_const_circulator cir     = target->vertex_begin(),
                                                cir_end = cir;
        CGAL_For_all(cir, cir_end)
            if (cir->opposite()->vertex() == source)
                return cir;

        assert(false);              // error if we reach this point
        return NULL;
    }
    Halfedge_handle get_halfedge(Vertex_handle source, Vertex_handle target)
    {
        Halfedge_const_handle halfedge = get_halfedge((Vertex_const_handle)source,
                                                      (Vertex_const_handle)target);
        return const_cast<Halfedge*>(&*halfedge);
    }

    //******************************************************************
    // LEVEL 1 INTERFACE:
    // for "normal" classes that do not deal with virtual seams
    // Example: all parameterization methods
    //******************************************************************

    // MESH INTERFACE

    // Get iterator over first vertex of mesh
    Vertex_iterator  mesh_vertices_begin() {
        return m_polyhedron.vertices_begin();
    }
    Vertex_const_iterator  mesh_vertices_begin() const {
        return m_polyhedron.vertices_begin();
    }

    // Get iterator over past-the-end vertex of mesh
    Vertex_iterator  mesh_vertices_end() {
        return m_polyhedron.vertices_end();
    }
    Vertex_const_iterator  mesh_vertices_end() const {
        return m_polyhedron.vertices_end();
    }

    // Count the number of vertices of the mesh
    int  count_mesh_vertices() const {
        int index = 0;
        for (Vertex_const_iterator it=mesh_vertices_begin(); it!=mesh_vertices_end(); it++)
            index++;
        return index;
    }

    // Index vertices of the mesh from 0 to count_mesh_vertices()-1
    void  index_mesh_vertices ()
    {
        //fprintf(stderr,"  index Parameterization_polyhedron_adaptor_3 vertices:\n");
        int index = 0;
        for (Vertex_iterator it=mesh_vertices_begin(); it!=mesh_vertices_end(); it++)
        {
            Point_3 position = get_vertex_position(it);
            //fprintf(stderr, "    %d=(%f,%f,%f)\n",
            //                index,
            //                (float)position.x(),
            //                (float)position.y(),
            //                (float)position.z());
            set_vertex_index(it, index++);
        }
        //fprintf(stderr,"    ok\n");
    }

    // Get iterator over first vertex of mesh's main border
    Border_vertex_iterator  mesh_main_border_vertices_begin() {
        return Border_vertex_iterator(m_main_border.begin());
    }
    Border_vertex_const_iterator  mesh_main_border_vertices_begin() const {
        return Border_vertex_const_iterator(m_main_border.begin());
    }

    // Get iterator over past-the-end vertex of mesh's main border
    Border_vertex_iterator  mesh_main_border_vertices_end() {
        return Border_vertex_iterator(m_main_border.end());
    }
    Border_vertex_const_iterator  mesh_main_border_vertices_end() const {
        return Border_vertex_const_iterator(m_main_border.end());
    }

    Halfedge_iterator  mesh_halfedges_begin() {
        return m_polyhedron.halfedges_begin();
    }
    Halfedge_const_iterator  mesh_halfedges_begin() const {
        return m_polyhedron.halfedges_begin();
    }

    Halfedge_iterator  mesh_halfedges_end() {
        return m_polyhedron.halfedges_end();
    }
    Halfedge_const_iterator  mesh_halfedges_end() const {
        return m_polyhedron.halfedges_end();
    }

    // Return the border containing seed_vertex
    // Return an empty list if not found
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

    // Get iterator over first facet of mesh
    Facet_iterator  mesh_facets_begin() {
        return m_polyhedron.facets_begin();
    }
    Facet_const_iterator  mesh_facets_begin() const {
        return m_polyhedron.facets_begin();
    }

    // Get iterator over past-the-end facet of mesh
    Facet_iterator  mesh_facets_end() {
        return m_polyhedron.facets_end();
    }
    Facet_const_iterator  mesh_facets_end() const {
        return m_polyhedron.facets_end();
    }

    // Count the number of facets of the mesh
    int  count_mesh_facets() const {
        int index = 0;
        for (Facet_const_iterator it=mesh_facets_begin(); it!=mesh_facets_end(); it++)
            index++;
        return index;
    }

    // Return true of all mesh's facets are triangles
    bool  is_mesh_triangular() const {
        for (Facet_const_iterator it = mesh_facets_begin(); it != mesh_facets_end(); it++)
            if (count_facet_vertices(it) != 3)
                return false;
        return true;            // mesh is triangular if we reach this point
    }

    // Count the number of halfedges of the mesh
    int  count_mesh_halfedges() const {
        int index = 0;
        for (Halfedge_iterator pHalfedge = m_polyhedron.halfedges_begin();
             pHalfedge != m_polyhedron.halfedges_end();
             pHalfedge++)
        {
            index++;
        }
        return index;
    }

    // FACET INTERFACE

    // Get circulator over facet's vertices
    Vertex_around_facet_circulator  facet_vertices_begin(Facet_handle facet) {
        return Vertex_around_facet_circulator(facet->facet_begin());
    }
    Vertex_around_facet_const_circulator  facet_vertices_begin(Facet_const_handle facet) const {
        return Vertex_around_facet_const_circulator(facet->facet_begin());
    }

    // Count the number of vertices of a facet
    int  count_facet_vertices(Facet_const_handle facet) const {
        int index = 0;
        Vertex_around_facet_const_circulator cir     = facet_vertices_begin(facet),
                                             cir_end = cir;
        CGAL_For_all(cir, cir_end)
            index++;
        return index;
    }

    // VERTEX INTERFACE

    // Get the 3D position of a vertex
    Point_3 get_vertex_position(Vertex_const_handle vertex) const {
        return vertex->point();
    }
	Point_3 get_vertex_position(Vertex_handle vertex) const {
		//vertex->
        return vertex->point();
    }

    // Get/set the 2D position (u/v pair) of a vertex. Default value is undefined.
    // (stored in halfedges sharing the same vertex)
    Point_2  get_vertex_uv(Vertex_const_handle vertex) const {
        return get_corners_uv(vertex, NULL, NULL);
    }
    void  set_vertex_uv(Vertex_handle vertex, const Point_2& uv) {
        set_corners_uv(vertex, NULL, NULL, uv);
    }

    // Get/set "is parameterized" field of vertex. Default value is undefined.
    // (stored in halfedges sharing the same vertex)
    bool  is_vertex_parameterized(Vertex_const_handle vertex) const {
        return are_corners_parameterized(vertex, NULL, NULL);
    }
    void  set_vertex_parameterized(Vertex_handle vertex, bool parameterized) {
        set_corners_parameterized(vertex, NULL, NULL, parameterized);
    }

    // Get/set vertex index. Default value is undefined.
    // (stored in Polyhedron vertex for debugging purpose)
    int  get_vertex_index(Vertex_const_handle vertex) const {
        //return get_corners_index(vertex, NULL, NULL);
        return vertex->index();
    }
    void  set_vertex_index(Vertex_handle vertex, int index)
    {
        //set_corners_index(vertex, NULL, NULL, index);
        vertex->index(index);
    }

    // Get/set vertex' all purpose tag. Default value is undefined.
    // (stored in halfedges sharing the same vertex)
    int  get_vertex_tag(Vertex_const_handle vertex) const {
        return get_corners_tag(vertex, NULL, NULL);
    }
    void set_vertex_tag(Vertex_handle vertex, int tag) {
        set_corners_tag(vertex, NULL, NULL, tag);
    }

    // Return true if a vertex belongs to ANY mesh's border
    bool  is_vertex_on_border(Vertex_const_handle vertex) const {
        return m_polyhedron.is_border(vertex);
    }

    // Return true if a vertex belongs to the UNIQUE mesh's main border,
    // ie the mesh's LONGEST border
    bool  is_vertex_on_main_border(Vertex_const_handle vertex) const {
        return std::find(m_main_border.begin(),
                         m_main_border.end(),
                         (Vertex*)&*vertex) != m_main_border.end();
    }

    // Get circulator over the vertices incident to 'vertex'
    // 'start_position' defines the optional initial position of the circulator
    Vertex_around_vertex_circulator vertices_around_vertex_begin(
                            Vertex_handle vertex,
                            Vertex_handle start_position = Vertex_handle())
    {
        if (start_position == NULL)
            return Vertex_around_vertex_circulator(vertex->vertex_begin());
        else
            return Vertex_around_vertex_circulator(
            		Halfedge_around_vertex_circulator(
                        	get_halfedge(start_position, vertex)));
    }
    Vertex_around_vertex_const_circulator vertices_around_vertex_begin(
                            Vertex_const_handle vertex,
                            Vertex_const_handle start_position = Vertex_const_handle()) const
    {
        if (start_position == NULL)
            return Vertex_around_vertex_const_circulator(vertex->vertex_begin());
        else
            return Vertex_around_vertex_const_circulator(
            		Halfedge_around_vertex_const_circulator(
                        	get_halfedge(start_position, vertex)));
    }

    //******************************************************************
    // LEVEL 2 INTERFACE:
    // for classes that deal with virtual seams
    // Example: Parameterization_mesh_patch_3
    //******************************************************************

    // VERTEX INTERFACE

    // Get/set vertex seaming flag. Default value is undefined.
    int  get_vertex_seaming(Vertex_const_handle vertex) const {
        return vertex->seaming();
    }
    void set_vertex_seaming(Vertex_handle vertex, int seaming) {
        vertex->seaming(seaming);
    }

    // EDGE INTERFACE

    // Get/set oriented edge's seaming flag, ie position of the oriented edge
    // wrt to the UNIQUE main border
    int  get_halfedge_seaming(Vertex_const_handle source, Vertex_const_handle target) const {
        return get_halfedge(source, target)->seaming();
    }
    void set_halfedge_seaming(Vertex_handle source, Vertex_handle target, int seaming) {
        get_halfedge(source, target)->seaming(seaming);
    }

    // CORNER INTERFACE

    // Get/set the 2D position (= (u,v) pair) of corners at the "right"
    // of the prev_vertex -> vertex -> next_vertex line.
    // Default value is undefined.
    // (stored in incident halfedges)
    Point_2 get_corners_uv(Vertex_const_handle vertex,
                           Vertex_const_handle prev_vertex,
                           Vertex_const_handle next_vertex) const
    {
        // if inner vertex
        if (prev_vertex == NULL && next_vertex == NULL)
        {
            // get (u,v) pair from any incident halfedge
            return Point_2(vertex->halfedge()->u(), vertex->halfedge()->v());
        }
        else // if seam vertex
        {
            assert(prev_vertex != NULL);
            assert(next_vertex != NULL);

            // get (u,v) pair from first inner halfedge (clockwise)
            Halfedge_around_vertex_const_circulator cir(
                                get_halfedge(next_vertex, vertex) );
            return Point_2(cir->u(), cir->v());
        }
    }
    void set_corners_uv(Vertex_handle vertex,
                        Vertex_const_handle prev_vertex,
                        Vertex_const_handle next_vertex,
                        const Point_2& uv)
    {
        // if inner vertex
        if (prev_vertex == NULL && next_vertex == NULL)
        {
            // Loop over all incident halfedges
            Halfedge_around_vertex_circulator cir     = vertex->vertex_begin(),
                                              cir_end = cir;
            CGAL_For_all(cir, cir_end)
                cir->uv(uv.x(), uv.y());
        }
        else // if seam vertex
        {
            assert(prev_vertex != NULL);
            assert(next_vertex != NULL);

            // first inner halfedge (for a clockwise rotation)
            Halfedge_around_vertex_circulator cir(
                                get_halfedge((Vertex*)&*next_vertex, vertex) );

            // past-the-end inner halfedge (for a clockwise rotation)
            Halfedge_around_vertex_circulator cir_end(
                                get_halfedge((Vertex*)&*prev_vertex, vertex) );

            // Loop over incident halfedges at the "right"
            // of the prev_vertex -> vertex -> next_vertex line
            CGAL_For_all(cir, cir_end)
                cir->uv(uv.x(), uv.y());
        }
    }

    // Get/set "is parameterized" field of corners at the "right"
    // of the prev_vertex -> vertex -> next_vertex line.
    // Default value is undefined.
    // (stored in incident halfedges)
    bool are_corners_parameterized(Vertex_const_handle vertex,
                                   Vertex_const_handle prev_vertex,
                                   Vertex_const_handle next_vertex) const
    {
        // if inner vertex
        if (prev_vertex == NULL && next_vertex == NULL)
        {
            // get "is parameterized" field from any incident halfedge
            return vertex->halfedge()->is_parameterized();
        }
        else // if seam vertex
        {
            assert(prev_vertex != NULL);
            assert(next_vertex != NULL);

            // get "is parameterized" field from first inner halfedge (clockwise)
            Halfedge_around_vertex_const_circulator cir(
                                get_halfedge(next_vertex, vertex) );
            return cir->is_parameterized();
        }
    }
    void set_corners_parameterized(Vertex_handle vertex,
                                   Vertex_const_handle prev_vertex,
                                   Vertex_const_handle next_vertex,
                                   bool parameterized)
    {
        // if inner vertex
        if (prev_vertex == NULL && next_vertex == NULL)
        {
            // Loop over all incident halfedges
            Halfedge_around_vertex_circulator cir     = vertex->vertex_begin(),
                                              cir_end = cir;
            CGAL_For_all(cir, cir_end)
                cir->is_parameterized(parameterized);
        }
        else // if seam vertex
        {
            assert(prev_vertex != NULL);
            assert(next_vertex != NULL);

            // first inner halfedge (for a clockwise rotation)
            Halfedge_around_vertex_circulator cir(
                                get_halfedge((Vertex*)&*next_vertex, vertex) );

            // past-the-end inner halfedge (for a clockwise rotation)
            Halfedge_around_vertex_circulator cir_end(
                                get_halfedge((Vertex*)&*prev_vertex, vertex) );

            // Loop over incident halfedges at the "right"
            // of the prev_vertex -> vertex -> next_vertex line
            CGAL_For_all(cir, cir_end)
                cir->is_parameterized(parameterized);
        }
    }

    // Get/set index of corners at the "right"
    // of the prev_vertex -> vertex -> next_vertex line.
    // Default value is undefined.
    // (stored in incident halfedges)
    int get_corners_index(Vertex_const_handle vertex,
                          Vertex_const_handle prev_vertex,
                          Vertex_const_handle next_vertex) const
    {
        // if inner vertex
        if (prev_vertex == NULL && next_vertex == NULL)
        {
            // get index from any incident halfedge
            return vertex->halfedge()->index();
        }
        else // if seam vertex
        {
            assert(prev_vertex != NULL);
            assert(next_vertex != NULL);

            // get index from first inner halfedge (clockwise)
            Halfedge_around_vertex_const_circulator cir(
                                get_halfedge(next_vertex, vertex) );
            return cir->index();
        }
    }
    void set_corners_index(Vertex_handle vertex,
                           Vertex_const_handle prev_vertex,
                           Vertex_const_handle next_vertex,
                           int index)
    {
        // if inner vertex
        if (prev_vertex == NULL && next_vertex == NULL)
        {
            // Loop over all incident halfedges
            Halfedge_around_vertex_circulator cir     = vertex->vertex_begin(),
                                              cir_end = cir;
            CGAL_For_all(cir, cir_end)
                cir->index(index);
        }
        else // if seam vertex
        {
            assert(prev_vertex != NULL);
            assert(next_vertex != NULL);

            // first inner halfedge (for a clockwise rotation)
            Halfedge_around_vertex_circulator cir(
                                get_halfedge((Vertex*)&*next_vertex, vertex) );

            // past-the-end inner halfedge (for a clockwise rotation)
            Halfedge_around_vertex_circulator cir_end(
                                get_halfedge((Vertex*)&*prev_vertex, vertex) );

            // Loop over incident halfedges at the "right"
            // of the prev_vertex -> vertex -> next_vertex line
            CGAL_For_all(cir, cir_end)
                cir->index(index);
        }
    }

    // Get/set all purpose tag of corners at the "right"
    // of the prev_vertex -> vertex -> next_vertex line.
    // Default value is undefined.
    // (stored in incident halfedges)
    int get_corners_tag(Vertex_const_handle vertex,
                        Vertex_const_handle prev_vertex,
                        Vertex_const_handle next_vertex) const
    {
        // if inner vertex
        if (prev_vertex == NULL && next_vertex == NULL)
        {
            // get tag from any incident halfedge
            return vertex->halfedge()->tag();
        }
        else // if seam vertex
        {
            assert(prev_vertex != NULL);
            assert(next_vertex != NULL);

            // get tag from first inner halfedge (clockwise)
            Halfedge_around_vertex_const_circulator cir(
                                get_halfedge(next_vertex, vertex) );
            return cir->tag();
        }
    }
    void set_corners_tag(Vertex_handle vertex,
                         Vertex_const_handle prev_vertex,
                         Vertex_const_handle next_vertex,
                         int tag)
    {
        // if inner vertex
        if (prev_vertex == NULL && next_vertex == NULL)
        {
            // Loop over all incident halfedges
            Halfedge_around_vertex_circulator cir     = vertex->vertex_begin(),
                                              cir_end = cir;
            CGAL_For_all(cir, cir_end)
                cir->tag(tag);
        }
        else // if seam vertex
        {
            assert(prev_vertex != NULL);
            assert(next_vertex != NULL);

            // first inner halfedge (for a clockwise rotation)
            Halfedge_around_vertex_circulator cir(
                                get_halfedge((Vertex*)&*next_vertex, vertex) );

            // past-the-end inner halfedge (for a clockwise rotation)
            Halfedge_around_vertex_circulator cir_end(
                                get_halfedge((Vertex*)&*prev_vertex, vertex) );

            // Loop over incident halfedges at the "right"
            // of the prev_vertex -> vertex -> next_vertex line
            CGAL_For_all(cir, cir_end)
                cir->tag(tag);
        }
    }


// Private operations
protected:

    // Extract mesh's longest border
    std::list<Vertex_handle> extract_longest_border(Polyhedron& mesh)
    {
        std::list<Vertex_handle> longest_border;    // returned list
        double                   max_len = 0;       // length of longest_border

        // Tag all vertices as unprocessed
        const int tag_free = 0;
        const int tag_done = 1;
        for (Vertex_iterator it=mesh_vertices_begin(); it!=mesh_vertices_end(); it++)
             set_vertex_tag(it, tag_free);

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
            std::list<Vertex_handle>::const_iterator it;
            for(it = border.begin(); it != border.end(); it++)
            {
                // Get next iterator (looping)
                std::list<Vertex_handle>::const_iterator next = it;
                next++;
                if (next == border.end())
                    next = border.begin();

                Vector_3 vect = get_vertex_position(*next) - get_vertex_position(*it);
                len += std::sqrt(vect*vect);
            }

            // Keep 'border' if longer
            if (len > max_len)
            {
                longest_border = border;
                max_len = len;
            }

            nb++;
        }

        return longest_border;
    }

    // Find a border tagged as "free" and tag it as "processed"
    // Return an empty list if not found
    std::list<Vertex_handle> find_free_border(int tag_free, int tag_done)
    {
        std::list<Vertex_handle> border;    // returned list

        // get any border vertex with "free" tag
        Vertex_handle seed_vertex = NULL;
        for (Vertex_iterator pVertex = mesh_vertices_begin();
             pVertex != mesh_vertices_end();
             pVertex++)
        {
            if (is_vertex_on_border(pVertex) && get_vertex_tag(pVertex) == tag_free) {
                seed_vertex = pVertex;
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
            set_vertex_tag(*it, tag_done);

        return border;
    }

// Fields
protected:

    // The adapted mesh (cannot be NULL)
    Polyhedron&                      m_polyhedron;

    // Main border of a topological disc inside m_polyhedron (may be empty)
    std::list<Vertex_handle>         m_main_border;

	// Private types
private:

    // Utility class to generate the Vertex_around_facet_circulator type
    struct Project_halfedge_vertex {
		typedef Halfedge                            argument_type;
		typedef typename Parameterization_polyhedron_adaptor_3::Vertex
			Vertex;
		typedef Vertex                              result_type;

		/// Get the target vertex of a halfedge
		Vertex&       operator()(Halfedge& h)       const {
			return *(h.vertex());
		}
		const Vertex& operator()(const Halfedge& h) const {
			return *(h.vertex());
		}		
    };
    friend struct Project_halfedge_vertex; // SUN's CC 5.50 requires that

    // Utility class to generate the Border_vertex_iterator type
    struct Project_vertex_handle_vertex {
        typedef Vertex_handle                       argument_type;
        typedef typename Parameterization_polyhedron_adaptor_3::Vertex
                                                    Vertex;
        typedef Vertex                              result_type;

        // Convert Vertex_handle to Vertex
        Vertex&       operator()(Vertex_handle& vh)       const { return *vh; }
        const Vertex& operator()(const Vertex_handle& vh) const { return *vh; }
    };
    friend struct Project_vertex_handle_vertex; // SUN's CC 5.50 requires that

    // This class is used to generate the Vertex_around_vertex_circulator type
    struct Project_opposite_halfedge_vertex {
        typedef Halfedge                            argument_type;
        typedef typename Parameterization_polyhedron_adaptor_3::Vertex
                                                    Vertex;
        typedef Vertex                              result_type;

        // Get the source vertex of a halfedge
        Vertex&       operator()(Halfedge& h)       const {
            return *(h.opposite()->vertex());
        }
        const Vertex& operator()(const Halfedge& h) const {
            return *(h.opposite()->vertex());
        }
    };
    friend struct Project_opposite_halfedge_vertex; // SUN's CC 5.50 requires that
};

DGAL_END_NAMESPACE
#endif//DGAL_PARAMETERIZATION_POLYHEDRON_ADAPTOR_H