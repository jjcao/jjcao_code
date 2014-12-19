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
// $Id: String_util.h 2007-04-07 $
//
//
// Author(s)     : JJCAO <jjcao@hotmail.com>


#ifndef DGAL_STRING_UTIL_H
#define DGAL_STRING_UTIL_H

#include <fstream>
#include <sstream>
#include <string>

namespace jj
{
bool parse_filename(const char* filename, std::string &basename, std::string &extname)
{
	basename = filename;
	extname.clear();

	std::string::size_type idx = basename.rfind('.');
	if ( idx == std::string::npos)//filename does not contain any peroid	
	{	
		return false;
	}
	
	//split filename to base name and extension; 
	//base name contains call characters before the period
	//extension contains all characters after the period
	extname = basename.substr(idx+1); 
	basename = basename.substr(0, idx);	
	
	if ( extname.empty())
	{
		false;
	}

	return true;
}

void roll_back_line(std::ifstream& stream, int len)
{
	for(int i = 0;i <= len; ++i)
	{
		stream.unget();
	}
}


//in will be truncated
std::string get_substr(std::string& in, const std::string& lable=" ")//, bool bRemoveSpace = true
{
	std::string::size_type idx(0), idx2(0);
	idx = in.find_first_not_of(lable, idx);
	idx2 = in.find_first_of(lable, idx);

	std::string result;
	if ( idx2 == std::string::npos)
	{
		if ( idx == std::string::npos)
			result = in;
		else
			result = in.substr(idx);
		in = "";
	}
	else
	{
		result = in.substr(idx, idx2);
		in = in.substr(idx2);
	}
	
	return result;
}
double str2double(std::string in)
{
	std::stringstream ss( in);
	double result;
	ss >> result;
	return result;
}
int str2int(std::string in)
{
	std::stringstream ss( in);
	int result;
	ss >> result;
	return result;
}
}

#endif//DGAL_STRING_UTIL_H