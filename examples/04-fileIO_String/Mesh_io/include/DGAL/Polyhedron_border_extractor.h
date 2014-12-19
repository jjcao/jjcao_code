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
// $Id: Polyhedron_border_extractor.h 2007-04-08 $
//
//
// Author(s)     : JJCAO <jjcao@hotmail.com>


#ifndef DGAL_BORDER_EXTRACTOR_H
#define DGAL_BORDER_EXTRACTOR_H

#include <DGAL/config.h>
#include <list>
#include <vector>

DGAL_BEGIN_NAMESPACE

template <class Polyhedron_3_>
class Polyhedron_border_extractor
{
private:
	typedef Polyhedron_3_                         Polyhedron;
	typedef typename Polyhedron::HDS              HDS;
	typedef typename HDS::Vertex_handle           Vertex_handle;
	typedef typename HDS::Halfedge_handle         Halfedge_handle;
	typedef typename HDS::Vertex_iterator         Vertex_iterator;
	typedef typename HDS::Halfedge_iterator       Halfedge_iterator;
	typedef typename Polyhedron::Point_3          Point_3;
	typedef typename Polyhedron::Vector_3         Vector_3;
	typedef typename Polyhedron::Halfedge_around_vertex_circulator
		                                          Halfedge_around_vertex_circulator;
public:
	typedef typename std::list< Vertex_handle >   Border_list;
	typedef typename std::vector< Border_list* >  Border_lists;
public:
	~Polyhedron_border_extractor()
	{
		clear();
	}
	void clear()
	{
		for (Border_lists::iterator itor = m_border_lists.begin();
			 itor != m_border_lists.end();
			 ++itor)
		{
			delete *itor;
		}

		m_border_lists.erase(m_border_lists.begin(), m_border_lists.end());
	}
	// Build numBorder border list from the border of the mesh 
	void extract(Polyhedron &mesh, int numBorder=4)
	{		
		clear();
		double border_len;

		Border_list wbl = extract_longest_border(mesh, &border_len);
		Border_list::const_iterator itor = wbl.begin();

		//create numBorder border list and push them all into m_border_lists
		for ( int i = 0; i < numBorder; ++i)
		{
			//by zhao
			if( i < numBorder - 1)
			{
				double pre_len = 0.0;
				double len = 0.0;
				Border_list::const_iterator it;
				for(it = itor; it != wbl.end(); ++it)
				{
					Border_list::const_iterator next = it;
					next++;
					if (next == wbl.end())
						next = wbl.begin();
					
					Vector_3 vect = (*next)->point() - (*it)->point();
					pre_len = len;
					len = pre_len + std::sqrt(vect*vect);
					
					if( pre_len < (border_len)/numBorder && len >= (border_len)/numBorder )
					{
						Border_list* bl = new Border_list();
						bl->insert(bl->end(), itor, it);
						m_border_lists.push_back(bl);
						itor = it;
						break;
					}
				}
			}
			else
			{
				Border_list* bl = new Border_list();
				Border_list::const_iterator itt = wbl.end();
				bl->insert(bl->end(), itor, itt);
				m_border_lists.push_back(bl);
			}
		}
	}

	Border_list* border_list(int indexBorder)
	{
		Border_list* result(0);

		if (indexBorder < 0) return result;

		if (indexBorder < (int)m_border_lists.size())
		{
			result = m_border_lists[indexBorder];
		}

		return result;
	}

private:
	 // Extract mesh's longest border
    Border_list extract_longest_border(Polyhedron& mesh, double *max_len)
    {
        Border_list              longest_border;    // returned list
								 *max_len = 0;       // length of longest_border

        // Tag all vertices as unprocessed
        const int tag_free = 0;
        const int tag_done = 1;
        for (Vertex_iterator it=mesh.vertices_begin(); it!=mesh.vertices_end(); it++)
             set_vertex_tag(it, tag_free);

        // find all closed borders and keep longest one
        int nb = 0;
        while (1)
        {
            // Find a border tagged as "free" and tag it as "processed"
            Border_list border = find_free_border(mesh, tag_free, tag_done);
            if(border.empty())
                break;

            // compute  total len of 'border'
            double len = 0.0;
            Border_list::const_iterator it;
            for(it = border.begin(); it != border.end(); it++)
            {
                // Get next iterator (looping)
                Border_list::const_iterator next = it;
                next++;
                if (next == border.end())
                    next = border.begin();

                Vector_3 vect = (*next)->point() - (*it)->point();
                len += std::sqrt(vect*vect);
            }

            // Keep 'border' if longer
            if (len > *max_len)
            {
                longest_border = border;
                *max_len = len;
            }

            nb++;
        }

        return longest_border;
    }
	void set_vertex_tag(Vertex_handle vertex, int tag)
    {
        // Loop over all incident halfedges
        Halfedge_around_vertex_circulator cir     = vertex->vertex_begin(),
                                          cir_end = cir;
        CGAL_For_all(cir, cir_end)
            cir->tag(tag);       
    }
	// Find a border tagged as "free" and tag it as "processed"
    // Return an empty list if not found
    Border_list find_free_border(Polyhedron& mesh, int tag_free, int tag_done)
    {
        Border_list border;    // returned list

        // get any border vertex with "free" tag
        Vertex_handle seed_vertex = NULL;
        for (Vertex_iterator pVertex = mesh.vertices_begin();
             pVertex != mesh.vertices_end();
             pVertex++)
        {
            if (mesh.is_border(pVertex) && pVertex->halfedge()->tag() == tag_free) {
                seed_vertex = pVertex;
                break;
            }
        }
        if (seed_vertex == NULL)
            return border;                  // return empty list

        // Get the border containing seed_vertex
        border = get_border(mesh, seed_vertex);

        // Tag border vertices as "processed"
        Border_list::iterator it;
        for(it = border.begin(); it != border.end(); it++)
            set_vertex_tag(*it, tag_done);

        return border;
    }
	// Return the border containing seed_vertex
    // Return an empty list if not found
    Border_list get_border(Polyhedron& mesh, Vertex_handle seed_vertex)
    {
        Border_list border;    // returned list

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
	Border_lists m_border_lists;
};

DGAL_END_NAMESPACE

#endif//DGAL_BORDER_EXTRACTOR_H