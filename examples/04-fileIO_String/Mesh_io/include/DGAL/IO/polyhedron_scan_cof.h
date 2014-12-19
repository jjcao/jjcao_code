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
// $Id: Polyhedron_scan_cof.h 2007-04-03 $
//
//
// Author(s)     : JJCAO <jjcao@hotmail.com>


#ifndef DGAL_BUILDER_COF_H
#define DGAL_BUILDER_COF_H

#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <DGAL/config.h>
#include <DGAL/UTIL/String_util.h>
#include <fstream>

DGAL_BEGIN_NAMESPACE

template <class Polyhedron_3_>
class Polyhedron_scan_cof : public CGAL::Modifier_base< typename Polyhedron_3_::HDS>
{
private:
	typedef Polyhedron_3_               Polyhedron;
	typedef typename Polyhedron::HDS    HDS;
	typedef typename HDS::Vertex_handle Vertex_handle;
	typedef typename HDS::Vertex::Point Point;
	typedef typename CGAL::Polyhedron_incremental_builder_3<HDS> 
		                                Builder;
private:
	std::ifstream&                      m_stream;//alloacated outside
	Polyhedron&                         m_polyhedron;
	std::vector<Vertex_handle>          m_vh; 
  
public:
	Polyhedron_scan_cof(Polyhedron &mesh, std::ifstream &stream):m_polyhedron(mesh),m_stream(stream)
	{
	}
	~Polyhedron_scan_cof() {}

	void operator()(HDS& hds)
	{
		Builder builder(hds,true);
		builder.begin_surface(3,1,6);
			ingore_addtional_info();
			read_vertices(builder);
			read_para_vertices(builder);
			read_facets(builder);
			read_constrained_vertices(builder);
		builder.end_surface();
	}
private:
	void ingore_addtional_info()
	{
		std::string cur_line;
		while(std::getline(m_stream,cur_line))
		{    
			int line_length = (int) cur_line.size();
			if (line_length == 0) continue;

			if( get_substr(cur_line) == "v")
			{
				roll_back_line(m_stream, line_length);
				break;
			}
		}
	}
	void read_vertices(Builder &builder)
	{
	  std::string cur_line;
	  while(std::getline(m_stream,cur_line))
	  {
		  int line_length = (int) cur_line.size();
		  if (line_length == 0) continue;

		  if( get_substr(cur_line) == "v")
		  {//v -0.00829827 0.583053 0.070454
			  double x, y, z;
			  x = str2double( get_substr(cur_line) );			 
			  y = str2double( get_substr(cur_line) );
			  z = str2double( get_substr(cur_line) );

			  m_vh.push_back( builder.add_vertex( Point(x, y, z) ) );
		  }
		  else
		  {
			  roll_back_line(m_stream, line_length);
			  break;
		  }	  
	  }
	}

// read parametered vertex coordinates
void read_para_vertices(Builder &builder) 
{
	int i(0);
	std::string cur_line;
	while(std::getline(m_stream,cur_line))
	{
	  int line_length = (int) cur_line.size();
	  if( get_substr(cur_line) == "vt")
	  {//vt 1.1 2.2
		  double u,v;
		  u = str2double( get_substr(cur_line) );			 
		  v = str2double( get_substr(cur_line) );
		  m_vh[i]->uv(u,v);
		  m_vh[i]->is_parameterized(true);
		  ++i;
	  }
	  else
	  {
		  roll_back_line(m_stream, line_length);
		  break;
	  }	  
	}
}

	// read facets and uv coordinates per halfedge
	void read_facets(Builder &builder)
	{
	  std::string cur_line;
	  while(std::getline(m_stream,cur_line))
	  {
		  int line_length = (int) cur_line.size();
		  if( get_substr(cur_line) == "f")
		  {//f 3 4 0
			  int v1,v2,v3;
			  v1 = str2int( get_substr(cur_line) );			 
			  v2 = str2int( get_substr(cur_line) );	
			  v3 = str2int( get_substr(cur_line) );	

			  builder.begin_facet();
				  builder.add_vertex_to_facet( v1);
				  builder.add_vertex_to_facet( v2);
				  builder.add_vertex_to_facet( v3);
			  builder.end_facet();
		  }
		  else
		  {
			  roll_back_line(m_stream, line_length);
			  break;
		  }	 
	  }
	}

	// read constrained vertex
	void read_constrained_vertices(Builder &builder)
	{
	  std::string cur_line;
	  while(std::getline(m_stream,cur_line))
	  {
		  int line_length = (int) cur_line.size();
		  if ( line_length == 0) continue;
		  if( get_substr(cur_line) == "cv")
		  {//"cv 0 0.5 0.5" or "cv 2"
			  int v1;
			  v1 = str2int( get_substr(cur_line) );	

			  Vertex_handle vh = m_vh[v1];
			  if ( !cur_line.empty())
			  {
				  double u, v;
				  u = str2double( get_substr(cur_line) );	
				  v = str2double( get_substr(cur_line) );
				  vh->uv(u,v);
			  }
			  vh->is_constrained(true);
			  m_polyhedron.add_constrained_vertex(vh);
		  }	  
		  else
		  {
			  std::cout<<"unkown row in cof!"<<std::endl;
		  }
	  }
	}
};

DGAL_END_NAMESPACE

#endif//DGAL_BUILDER_COF_H