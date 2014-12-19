#pragma once

// ***************************************************************
//  cgalPointToVtk   version:  1.0   ¡¤  date: 09/30/2008
//  -------------------------------------------------------------
//                           HuiWang
//  -------------------------------------------------------------
//  Copyright (C) 2008 - All Rights Reserved
// ***************************************************************
//                Display the CGAL point in VTK
// ***************************************************************

#include "stdafx.h"

#include <vtkSphereSource.h>
#include <vtkPolyDataMapper.h>
#include <vtkLODActor.h>
#include <vtkRenderer.h>
#include <vtkProperty.h>



class cgalPointToVtk
{
public:
	vtkSphereSource* m_sphereSource;
	vtkPolyDataMapper* m_polyDataMapper;
	vtkLODActor* m_actor;
	vtkRenderer* m_renderer;

public:
	cgalPointToVtk() {}
	cgalPointToVtk(vtkRenderer* renderer,double p[3],double radius,double color[3]);
	~cgalPointToVtk();

private:
	void clear();
};