#include "stdafx.h"
#include "cgalPointToVtk.h"

cgalPointToVtk::cgalPointToVtk(vtkRenderer* renderer,Point_3 p,double radius,double color[3])
{
	m_renderer = renderer;
	double center[3] = {p.x(),p.y(),p.z()};
	m_sphereSource = vtkSphereSource::New();
	m_sphereSource->SetCenter(center);
	m_sphereSource->SetRadius(radius);
	m_polyDataMapper = vtkPolyDataMapper::New();
	m_polyDataMapper->SetInputConnection(m_sphereSource->GetOutputPort() );
	m_actor = vtkLODActor::New();
	m_actor->SetMapper(m_polyDataMapper);
	//m_actor->PickableOff();
	m_actor->GetProperty()->SetColor(color[0],color[2],color[3]);
	m_renderer->AddActor(m_actor);

}

cgalPointToVtk::~cgalPointToVtk()
{
	clear();
}

void cgalPointToVtk::clear()
{
	m_sphereSource->Delete();
	m_polyDataMapper->Delete();
	m_actor->Delete();
	/*if(m_sphereSource)
	{
		m_sphereSource->Delete();
		m_sphereSource = 0;
	}
	if(m_polyDataMapper)
	{
		m_polyDataMapper->Delete();
		m_polyDataMapper =0;
	}
	*/
}