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
// $Id: Polyhderon_outputer_cof.h 2007-04-07 $
//
//
// Author(s)     : JJCAO <jjcao@hotmail.com>


#ifndef DGAL_POLYHEDRON_OUTPUTER_COF_H
#define DGAL_POLYHEDRON_OUTPUTER_COF_H

#include <DGAL/config.h>
#include <sstream>
#include <fstream>

DGAL_BEGIN_NAMESPACE

template <class Polyhedron_3>
class Polyhderon_outputer_cof
{
public:
    typedef Polyhedron_3                                 Polyhedron;
private:
	typedef typename Polyhedron::Vertex_handle           Vertex_handle;
	typedef typename Polyhedron::Vertex_const_handle     Vertex_const_handle;
    typedef typename Polyhedron::Vertex_const_iterator   Vertex_const_iterator;
	typedef typename Polyhedron::Vertex_handle_itor      Vertex_handle_itor;//
	typedef typename Polyhedron::Halfedge_const_iterator Halfedge_const_iterator;
	typedef typename Polyhedron::Facet_const_iterator    Facet_const_iterator;
	typedef typename Polyhedron::Halfedge_around_facet_const_circulator    
	                                                     Halfedge_around_facet_const_circulator;
private:
    Polyhedron&    m_polyhedron;
    std::ofstream& m_stream;//alloacated outside
public:
	Polyhderon_outputer_cof(Polyhedron &mesh, std::ofstream &stream):m_polyhedron(mesh), m_stream(stream)
	{
	}
	~Polyhderon_outputer_cof() {}

	void save()
	{
		CGAL::set_ascii_mode(m_stream);

		// Index all mesh vertices following the order of vertices_begin() iterator
		m_polyhedron.precompute_vertex_indices();
		// Index all mesh half edges following the order of halfedges_begin() iterator
		m_polyhedron.precompute_halfedge_indices();

		// output coordinates
		Vertex_const_iterator pVertex;
		for(pVertex = m_polyhedron.vertices_begin(); 
			pVertex != m_polyhedron.vertices_end(); ++pVertex)
		{
			m_stream << "v " << pVertex->point().x() << " "
						<< pVertex->point().y() << " "
						<< pVertex->point().z() << std::endl;
		}

		// Write UVs (1 UV / vertex)
		for(pVertex = m_polyhedron.vertices_begin(); 
			pVertex != m_polyhedron.vertices_end(); ++pVertex)
		{
			Vertex_const_handle vh(pVertex);
			if (vh->is_parameterized())
				m_stream << "vt " << vh->u() << " " << vh->v() << std::endl;
			//else
			//	m_stream << "vt " << 0.0 << " " << 0.0 << std::endl;
		}

		// Write facets
		Facet_const_iterator pFacet;
		for(pFacet = m_polyhedron.facets_begin(); pFacet != m_polyhedron.facets_end(); ++pFacet)
		{
			Halfedge_around_facet_const_circulator h = pFacet->facet_begin();
			m_stream << "f";
			do {
				m_stream << " " << h->vertex()->index();
				if (h->is_parameterized())
					m_stream <<  "/" << h->index();
			}
			while(++h != pFacet->facet_begin());
			m_stream << std::endl;
		}

		// Write constraint
		for(Vertex_handle_itor pVh = m_polyhedron.constrained_vertices_begin();
			pVh != m_polyhedron.constrained_vertices_end(); ++pVh)
		{
			Vertex_handle vh = *pVh;
			m_stream << "cv " << vh->index();
			if (vh->is_parameterized())
			{
				m_stream <<  " " << vh->u() << " " << vh->v();
			}
			m_stream << std::endl;
		}

	}

};

DGAL_END_NAMESPACE

#endif//DGAL_POLYHEDRON_OUTPUTER_COF_H