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
// $Id: config.h 2007-03-28 $
// 
//
// Author(s)     : JJCAO <jjcao@hotmail.com>

#ifndef DGAL_CONFIG_H
#define DGAL_CONFIG_H

#define DGAL_VERSION 0.0.1

#define DGAL_BEGIN_NAMESPACE namespace DGAL {
#define DGAL_END_NAMESPACE }

#include <CGAL/Cartesian.h>	
#include <CGAL/Timer.h>

DGAL_BEGIN_NAMESPACE
enum Vertex_type {
	FIRST_BORDER_VERTEX  = -12,
	SECOND_BORDER_VERTEX  = -11,	
	CONTROL_VERTEX  = -4,
	MARKED_VERTEX  = 1000
};

typedef double Num_type_default;
typedef CGAL::Cartesian<Num_type_default> Kernel_default;
DGAL_END_NAMESPACE

#endif // DGAL_CONFIG_H
