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
// $Id:  Polyhedron_scan_obj.h 2007-04-03 $
//
//
// Author(s)     : JJCAO <jjcao@hotmail.com>


#ifndef DGAL_POLYHEDRON_SCAN_OBJ_H
#define DGAL_POLYHEDRON_SCAN_OBJ_H

#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <DGAL/config.h>
#include <DGAL/UTIL/String_util.h>
#include <fstream>
#include <locale>

DGAL_BEGIN_NAMESPACE

// The file format
// # some text
// Line is a comment until the end of the line
// v float float float
// A single vertex's geometric position in space. The first vertex listed in the file has index 1, and subsequent vertices are numbered sequentially.
// vn float float float
// A normal. The first normal in the file is index 1, and subsequent normals are numbered sequentially.
// vt float float
// A texture coordinate. The first texture coordinate in the file is index 1, and subsequent textures are numbered sequentially.
// f int int int ...
// or
// f int/int int/int int/int . . .
// or
// f int/int/int int/int/int int/int/int ...

template <class HDS>
class Polyhedron_scan_obj : public CGAL::Modifier_base<HDS>
{
private:
  typedef typename HDS::Vertex::Point Point;
  typedef typename HDS::Vertex_handle Vertex_handle;
  typedef typename CGAL::Polyhedron_incremental_builder_3<HDS> Builder;

private:
	std::ifstream&                      m_stream;//alloacated outside
	std::vector<Vertex_handle>          m_vh;
	std::vector<double> m_u;
	std::vector<double> m_v;

public:
	Polyhedron_scan_obj( std::ifstream &stream):m_stream(stream)
	{
	}
	~Polyhedron_scan_obj() {}

  void operator()(HDS& hds)
  { 
	  Builder builder(hds,true);
	  builder.begin_surface(3,1,6);
		read_vertices(builder);		
		read_para_vertices(builder);
		read_facets(builder);		
	  builder.end_surface();
  }

private:
// read vertex coordinates
  void read_vertices(Builder &builder)
  {
	  std::string cur_line;
	  //ingore some comment lines before "v *"
	  while(std::getline(m_stream,cur_line))
	  {
		  int line_length = (int) cur_line.size();
		  if (line_length > 0 && get_substr(cur_line) == "v")
		  {
			  roll_back_line(m_stream, line_length);
			  break;
		  }	 
	  }

	  //begin read vertices
	  while(std::getline(m_stream,cur_line))
	  {
		  int line_length = (int) cur_line.size();
		  if( line_length > 0 && get_substr(cur_line) == "v")
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

	//ingore some comment lines before "vt *"
	while(std::getline(m_stream,cur_line))
	{
		int line_length = (int) cur_line.size();
		if (line_length > 0)
		{
			std::string tmp = get_substr(cur_line);
			if ( tmp == "vt")
			{
				roll_back_line(m_stream, line_length);
				break;
			}
			else if (tmp == "f")
			{
				roll_back_line(m_stream, line_length);
				return;
			}
		}
	}

	//begin read parameter vertices
	while(std::getline(m_stream,cur_line))
	{
	  int line_length = (int) cur_line.size();
	  if( line_length > 0 && get_substr(cur_line) == "vt")
	  {//vt 1.1 2.2
		  double u,v;
		  u = str2double( get_substr(cur_line) );			 
		  v = str2double( get_substr(cur_line) );
		  m_u.push_back(u);
		  m_v.push_back(v);		  
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

	//ingore some comment lines before "f *"
	while(std::getline(m_stream,cur_line))
	{
		int line_length = (int) cur_line.size();
		if (line_length > 0 && get_substr(cur_line) == "f")
		{
			roll_back_line(m_stream, line_length);
			break;
		}
	}

	//begin read faces
	while(std::getline(m_stream,cur_line))
	{
		int line_length = (int) cur_line.size();
		if( line_length > 0 && cur_line.at(0) == 'f')
		{//f 3/3 2/2 1/1 or f 3 2 1

			int p[3], uv[3];
			parse_face_line(cur_line, p, uv);
			
			for (int i = 0; i < 3; ++i)
			{
				m_vh[p[i]]->uv(m_u[uv[i]],m_v[uv[i]]);	
				m_vh[p[i]]->is_parameterized(true);
			}
			
			builder.begin_facet();
			builder.add_vertex_to_facet( p[0]);
			builder.add_vertex_to_facet( p[1]);
			builder.add_vertex_to_facet( p[2]);
			builder.end_facet();
		}
		else
		{
			break;
		}	 
	}

}
private:
	void parse_face_line(std::string& in, int* p, int* uv)
	{
		//remove "f "
		in.erase(0, 2);

		//to 3 groups by spaces
		std::string::size_type idx(0);
		std::string str[3];
		for ( int i = 0; i < 3; ++i)
		{
			idx = in.find_first_of(" ");
			if ( idx == std::string::npos)
			{
				str[i] = in;
			}
			else
			{
				str[i] = in.substr(0, idx);
				in.erase(0, idx+1);
			}			
		} 
		
		//is contain texture
		std::locale loc1;
		for ( int i = 0; i < 3; ++i)
		{
			const char *string = str[i].c_str();

			const char* tmp = std::use_facet<std::ctype<char> > ( loc1 ).scan_not
			( std::ctype_base::digit, &string[0], &string[strlen(&string[0])] );

			idx = tmp - string;
			if ( idx)
			{
				std::string tmp = str[i].substr(0, idx);
				p[i] = str2int(tmp)-1;
				
				tmp = str[i].substr(idx+1);
				uv[i] = str2int(tmp)-1;
			}
			else
			{
				p[i] = str2int(str[i])-1;
				uv[i] = p[i];
			}			
		}
	}
};

DGAL_END_NAMESPACE

#endif//DGAL_POLYHEDRON_SCAN_OBJ_H