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
// $Id: Polyhedron_io.h 2007-04-09 $
//
//
// Author(s)     : JJCAO <jjcao@hotmail.com>


#ifndef DGAL_POLYHEDRON_IO_H
#define DGAL_POLYHEDRON_IO_H

#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/IO/Polyhedron_inventor_ostream.h> 
#include <CGAL/IO/Polyhedron_VRML_1_ostream.h> 
#include <CGAL/IO/Polyhedron_VRML_2_ostream.h> 
#include <CGAL/IO/Polyhedron_geomview_ostream.h>


#include <DGAL/config.h>
#include <DGAL/IO/Polyhedron_scan_cof.h>
#include <DGAL/IO/Polyhedron_scan_obj.h>
#include <DGAL/IO/Polyhderon_outputer_cof.h>
#include <DGAL/IO/Polyhedron_outputer_obj.h>
#include <DGAL/UTIL/String_util.h>

#include <sstream>
#include <fstream>
#include <string>

DGAL_BEGIN_NAMESPACE

class Polyhedron_io
{
public:
	 enum Error_code {
		OK,                             ///< Success
		UNSUPPORTED_FORMAT,
		FAILED_TO_OPEN_FILE,
		UNKNOWN_ERROR
		};
public:	 
	//
	
	template<class Polyhedron>
	static Error_code scan(Polyhedron& mesh, const char* filename)
	{
		std::string basename;		std::string extname; 

		if ( !parse_filename(filename, basename, extname))
		{
			return UNSUPPORTED_FORMAT;
		}

		//////////////////////////////////
		std::ifstream stream(filename);
		if(!stream)
		{
			return FAILED_TO_OPEN_FILE;
		}

		Error_code result = OK;
		std::transform(extname.begin(), extname.end(), extname.begin(), tolower);
		if ( extname == "cof")
		{
			Polyhedron_scan_cof<Polyhedron> bd(mesh, stream);
			mesh.delegate(bd);
		}
		else if ( extname == "obj")
		{
			Polyhedron_scan_obj<Polyhedron::HDS> bd(stream);
			mesh.delegate(bd);
		}
		else if ( extname == "off")
		{
			stream >> mesh;
		}
		else
		{
			result = UNSUPPORTED_FORMAT;
		}

		stream.close();
		return result;
	}

	template<class Polyhedron>
	static Error_code output(Polyhedron& mesh, const char* filename)
	{
		std::string basename;		std::string extname; 

		if ( !parse_filename(filename, basename, extname))
		{
			return UNSUPPORTED_FORMAT;
		}

		//////////////////////////////////
		std::ofstream ostream(filename);
		if(!ostream) return FAILED_TO_OPEN_FILE;

		Error_code result = OK;
		std::transform(extname.begin(), extname.end(), extname.begin(), tolower);
		if ( extname == "cof")
		{
			/*Polyhderon_outputer_cof<Polyhedron> outputer(mesh, ostream);
			outputer.save();*/
		}
		else if ( extname == "obj")
		{
			Polyhedron_outputer_obj<Polyhedron> outputer(mesh, ostream);
			outputer.save();
		}
		else if ( extname == "off")
		{	
			ostream << mesh;
		}
		else if ( extname == "iv")
		{	
			CGAL::Inventor_ostream out( ostream);
			out << mesh;
		}
		else if ( extname == "wrl")
		{	
			CGAL::VRML_2_ostream out( ostream);
			out << mesh;
		}
		//else if ( extname == "gv")
		//{	
		//	CGAL::Geomview_stream out;
		//	out << mesh;
		//}
		else
		{
			result = UNSUPPORTED_FORMAT;
		}

		ostream.close();
		return result;
	}

};


DGAL_END_NAMESPACE

#endif//DGAL_POLYHEDRON_IO_H