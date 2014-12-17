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
// This file is adapted from mesh_cutter of Pierre Alliez (pierre.alliez@sophia.inria.fr).
//
// $URL: http://jjcao1231.googlepages.com $
// $Id: Mesh_cutter_simple.h 2007-04-08 $
//
//
// Author(s)     : JJCAO <jjcao@hotmail.com>

#ifndef DGAL_MESH_CUTTER_SIMPLE_H
#define DGAL_MESH_CUTTER_SIMPLE_H

#include <list>


DGAL_BEGIN_NAMESPACE

template
<
	class Polyhedron_3_
>
class Mesh_cutter_simple
{
// Public types
public:
    typedef Polyhedron_3_                               Polyhedron;
	typedef double                                      Number_type;
	typedef typename Polyhedron::Halfedge_handle        Halfedge_handle;
    typedef std::list<Halfedge_handle>                  Backbone;
	typedef typename Polyhedron::Vertex_handle          Vertex_handle;
	typedef typename Polyhedron::Facet_handle           Facet_handle;

// Private types
private:
    typedef typename Polyhedron::Vector_3                        Vector_3;
    typedef typename Polyhedron::Point_3                         Point_3;
    enum {FREE,DONE,FIXED};

// Public operations
public:
    // life cycle
	Mesh_cutter_simple(Polyhedron& polyhedron)
		:m_pPolyhedron(&polyhedron),m_pBackbone(NULL)
    {
    }
    ~Mesh_cutter_simple() {}

    void cut(Backbone& backbone);
    void cut_genus(Backbone& backbone);

// Private operations
private:

    // genus > 0
    bool init();
    bool simplify();
    bool extend();
    void precompute_distances();
	void recursive_tag(Facet_handle pFacet,int index);

	//***************************************************
	// pick_best_halfedge
	//***************************************************
    Halfedge_handle 
		pick_best_halfedge( typename std::list< Halfedge_handle>::iterator &pos)
	{
		Halfedge_handle pBest = NULL;
		double min_distance = 1e308; //

		// cleanup
		std::list<Halfedge_handle>::iterator iter;
		for(iter  = m_pBackbone->begin();
			iter != m_pBackbone->end();
			iter++)
		{
			Halfedge_handle pHalfedge = (*iter);
			Halfedge_handle opposite = pHalfedge->opposite();
			Facet_handle pFacet = opposite->facet();

			// check
			if(pHalfedge->is_border() ||
				pFacet == NULL)
				continue;

			if(pFacet->tag() == DONE)
				continue;

			// no border vertex
			Vertex_handle pVertex = opposite->next()->vertex();
			if(m_pPolyhedron->is_border(pVertex))
				continue;

			// precomputed distance
			double distance = pHalfedge->distance();
			if(distance < min_distance)
			{
				pos = iter;
				pBest = pHalfedge;
				min_distance = distance;
			}
		}
		return pBest;
	}
    

// Fields
private:

    Polyhedron*                 m_pPolyhedron;  // the model to cut
    Backbone*                   m_pBackbone;    // the backbone to fill
    Facet_handle                m_pSeedFacet;
    Vertex_handle               m_pSeedVertex;
};

//***************************************************
// simple cut for genus 0 mesh
//***************************************************
template<class Polyhedron>
inline
void Mesh_cutter_simple<Polyhedron>::cut(Backbone& backbone)
{
    m_pBackbone = &backbone;

    // special init -> tag all vertices, but two
    m_pPolyhedron->tag_vertices(FREE);
    Polyhedron::Vertex_handle pVertexMin,pVertexMax;
    m_pPolyhedron->farthest_point_aligned(pVertexMin,pVertexMax);
    pVertexMin->tag(FIXED);
    pVertexMax->tag(FIXED);
    init();

    // cutting
    while(extend()) {}
}

/////////////////////////////////////////////////////
// GENUS > 0
/////////////////////////////////////////////////////

//***************************************************
// cut for genus>0 mesh
//***************************************************
template<class Polyhedron>
inline
void Mesh_cutter_simple<Polyhedron>::cut_genus(Backbone& backbone)
{
    m_pBackbone = &backbone;

    // init
    m_pPolyhedron->tag_vertices(FREE); // all free
    init();

    // cutting
    while(extend()) {}
}

//***************************************************
// init
//***************************************************
template<class Polyhedron>
inline
bool Mesh_cutter_simple<Polyhedron>::init()
{
    // tag facets
    m_pPolyhedron->tag_facets(FREE);

    // compute bounding box and center
    double xmin = m_pPolyhedron->xmin();
    double ymin = m_pPolyhedron->ymin();
    double zmin = m_pPolyhedron->zmin();
    double xmax = m_pPolyhedron->xmax();
    double ymax = m_pPolyhedron->ymax();
    double zmax = m_pPolyhedron->zmax();
    double xcenter = 0.5*(xmin+xmax);
    double ycenter = 0.5*(ymin+ymax);
    double zcenter = 0.5*(zmin+zmax);
    Point_3 center(xcenter,ycenter,zcenter);

    // get closest facet
    m_pSeedFacet = m_pPolyhedron->get_closest_inner_facet(center);
    CGAL_assertion(m_pSeedFacet != NULL);

    Polyhedron::Halfedge_handle he = m_pSeedFacet->halfedge();
    CGAL_assertion(he != NULL);
    CGAL_assertion(m_pBackbone != NULL);
    m_pBackbone->push_back(he);
    m_pBackbone->push_back(he->next());
    m_pBackbone->push_back(he->next()->next());

    precompute_distances();
    m_pSeedFacet->tag(DONE);

    return true;
}

//***************************************************
// extend
//***************************************************
template<class Polyhedron>
inline
bool Mesh_cutter_simple<Polyhedron>::extend()
{
    std::list<Polyhedron::Halfedge_handle>::iterator pos;
    Polyhedron::Halfedge_handle pHalfedge = pick_best_halfedge(pos);
    if(pHalfedge == NULL)
    return false;

    // flag facet
    pHalfedge->opposite()->facet()->tag(DONE);

    // insert halfedge
    std::list<Polyhedron::Halfedge_handle>::iterator tmp =
    m_pBackbone->insert(pos,pHalfedge->opposite()->next()->next());
    m_pBackbone->insert(tmp,pHalfedge->opposite()->next());

    // remove this one
    m_pBackbone->remove(pHalfedge);

    // simplify current backbone
    while(simplify());
    return true;
}

//***************************************************
// simplify
//***************************************************
template<class Polyhedron>
inline
bool Mesh_cutter_simple<Polyhedron>::simplify()
{
    // cleanup
    std::list<Polyhedron::Halfedge_handle>::iterator iter;
    for(iter  = m_pBackbone->begin();
        iter != m_pBackbone->end();
        iter++)
    {
        Polyhedron::Halfedge_handle pHalfedge = (*iter);
        Polyhedron::Halfedge_handle opposite = pHalfedge->opposite();

        // get next halfedge in the list
        iter++;
        Polyhedron::Halfedge_handle pNext = NULL;
        if(iter == m_pBackbone->end()) // loop
            pNext = (*m_pBackbone->begin());
        else
            pNext = (*iter);

        if(pNext == opposite &&
            pHalfedge->vertex()->tag() == FREE)
        {
            m_pBackbone->remove(pHalfedge);
            m_pBackbone->remove(opposite);
            return true;
        }

        iter--; // restore
    }
    return false;
}

//***************************************************
// precompute_distances
//***************************************************
template<class Polyhedron>
inline
void Mesh_cutter_simple<Polyhedron>::precompute_distances()
{
    Polyhedron::Halfedge_iterator pHalfedge;
    for(pHalfedge = m_pPolyhedron->halfedges_begin();
        pHalfedge != m_pPolyhedron->halfedges_end();
        pHalfedge++)
    pHalfedge->distance(m_pPolyhedron->distance(m_pSeedFacet,pHalfedge));
}

// Cut the mesh to make it homeomorphic to a disk
// or extract a region homeomorphic to a disc.
// Return the border of this region (empty on error)
//
// CAUTION:
// This method is provided "as is". It is very buggy and simply part of this example.
// Developers using this package should implement a more robust cut algorithm!
template <class Parameterization_polyhedron_adaptor>
std::list<typename Parameterization_polyhedron_adaptor::Vertex_handle>
cut_mesh(Parameterization_polyhedron_adaptor& mesh_adaptor)
{
	typedef Parameterization_polyhedron_adaptor       Polyhedron_adaptor;
	typedef typename std::list<Polyhedron_adaptor::Vertex_handle> 
		                                              Seam;
    // Helper class to compute genus or extract borders
    typedef CGAL::Parameterization_mesh_feature_extractor<Polyhedron_adaptor>
                                                      Mesh_feature_extractor;
    typedef Mesh_feature_extractor::Border            Border;
	typedef typename Polyhedron_adaptor::Polyhedron   WrappedMesh;
	typedef Mesh_cutter_simple<WrappedMesh>           Mesh_cutter;
    typedef Mesh_cutter::Backbone                     Backbone;

    Seam seam;              // returned list

    // Get refererence to Polyhedron_3 mesh
    Polyhedron& mesh = mesh_adaptor.get_adapted_mesh();

    // Extract mesh borders and compute genus
    Mesh_feature_extractor feature_extractor(mesh_adaptor);
    int nb_borders = feature_extractor.get_nb_borders();
    int genus = feature_extractor.get_genus();

    // If mesh is a topological disk
    if (genus == 0 && nb_borders > 0)
    {
        // Pick the longest border
        seam = feature_extractor.get_longest_border();
    }
    else // if mesh is NOT a topological disk, create a virtual cut
    {
        Backbone seamingBackbone;           // result of cutting
        Backbone::iterator he;

        // Compute a cutting path that makes the mesh a "virtual" topological disk
        mesh.compute_facet_centers();
        Mesh_cutter cutter(mesh);
        if (genus == 0)
        {
            // no border, we need to cut the mesh
            assert (nb_borders == 0);
            cutter.cut(seamingBackbone);    // simple cut
        }
        else // genus > 0 -> cut the mesh
        {
            cutter.cut_genus(seamingBackbone);
        }

        // The Mesh_cutter_simple class is quite buggy
        // => we check that seamingBackbone is valid
        //
        // 1) Check that seamingBackbone is not empty
        if (seamingBackbone.begin() == seamingBackbone.end())
            return seam;                    // return empty list
        //
        // 2) Check that seamingBackbone is a loop and
        //    count occurences of seam halfedges
        mesh.tag_halfedges(0);              // Reset counters
        for (he = seamingBackbone.begin(); he != seamingBackbone.end(); he++)
        {
            // Get next halfedge iterator (looping)
            Backbone::iterator next_he = he;
            next_he++;
            if (next_he == seamingBackbone.end())
                next_he = seamingBackbone.begin();

            // Check that seamingBackbone is a loop: check that
            // end of current HE == start of next one
            if ((*he)->vertex() != (*next_he)->opposite()->vertex())
                return seam;                // return empty list

            // Increment counter (in "tag" field) of seam halfedges
            (*he)->tag( (*he)->tag()+1 );
        }
        //
        // 3) check that the seamingBackbone is a two-way list
        for (he = seamingBackbone.begin(); he != seamingBackbone.end(); he++)
        {
            // Counter of halfedge and opposite halfedge must be 1
            if ((*he)->tag() != 1 || (*he)->opposite()->tag() != 1)
                return seam;                // return empty list
        }

        // Convert list of halfedges to a list of vertices
        for (he = seamingBackbone.begin(); he != seamingBackbone.end(); he++)
            seam.push_back((*he)->vertex());
    }

    return seam;
}

DGAL_END_NAMESPACE

#endif // DGAL_MESH_CUTTER_SIMPLE_H
