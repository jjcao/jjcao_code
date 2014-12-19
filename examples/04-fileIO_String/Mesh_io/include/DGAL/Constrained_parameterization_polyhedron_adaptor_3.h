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
// $Id: Constrained_parameterization_polyhedron_adaptor_3.h 2007-04-04 $
//
//
// Author(s)     : JJCAO <jjcao@hotmail.com>

#ifndef DGAL_CONSTRAINED_PARAMETERIZATION_POLYHEDRON_ADAPTOR_H
#define DGAL_CONSTRAINED_PARAMETERIZATION_POLYHEDRON_ADAPTOR_H

#include <DGAL/config.h>
#include <DGAL/Parameterization_constrain_policies_3.h>
#include <DGAL/Parameterization_polyhedron_adaptor_3.h>

DGAL_BEGIN_NAMESPACE

template
<
	class Polyhedron_3_
>
class Constrained_parameterization_polyhedron_adaptor_3
	:public Parameterization_polyhedron_adaptor_3<Polyhedron_3_>
{
// Private types
private:
	// Superclass
	typedef Parameterization_polyhedron_adaptor_3<Polyhedron_3_>  
		                                             Base;
public:
	typedef typename std::list<Vertex_handle>::iterator Constrain_vertex_iterator;
	typedef typename std::list<Vertex_handle>::const_iterator 
		                                                Constrain_vertex_const_iterator;

//Public operations
public:
	Constrained_parameterization_polyhedron_adaptor_3(Polyhedron& mesh)
	: Parameterization_polyhedron_adaptor_3(mesh)
	{
		//set uv for constrained vertices
		for (Constrain_vertex_iterator it = constrained_vertices().begin();
             it != constrained_vertices().end();    ++it)
		{		
			Vertex_handle vh = *it;
			set_vertex_uv(vh, Point_2(vh->u(), vh->v()));
		}
	}
	virtual bool  is_constrained_vertex(Vertex_const_handle vertex) const 
	{
		return vertex->is_constrained();
    }
	std::list<Vertex_handle>& constrained_vertices(){return m_polyhedron.constrained_vertices();}
};

DGAL_END_NAMESPACE
#endif//DGAL_CONSTRAINED_PARAMETERIZATION_POLYHEDRON_ADAPTOR_H