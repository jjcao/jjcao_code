#ifndef QT_VTK_H
#define QT_VTK_H

#include "qt_vtk_global.h"

#include <vtkAxesActor.h>
#include <vtkOrientationMarkerWidget.h>
#include <vtkRenderWindowInteractor.h>

template<class Polyhedron>	
vtkPolyData* cgal2vtk(Polyhedron* mesh)
{
	vtkPoints *points = vtkPoints::New();
	int pointNum(0);
	for ( Vertex_iterator itor = mesh->vertices_begin(); itor != mesh->vertices_end(); ++itor)
	{
		Point_3 p3 = itor->point();
		points->InsertNextPoint( p3.x(), p3.y(), p3.z());
		itor->index(pointNum);
		++pointNum;
	}
	vtkCellArray *polys = vtkCellArray::New();
	for ( Face_iterator itor = mesh->facets_begin(); itor != mesh->facets_end(); ++itor)
	{
		polys->InsertNextCell(itor->facet_degree());
		Halfedge_facet_circulator j = itor->facet_begin();        
		do {            
			polys->InsertCellPoint(j->vertex()->index());
		} while ( ++j != itor->facet_begin());		
	}

	vtkPolyData* result = vtkPolyData::New();
	result->SetPoints(points);
	result->SetPolys(polys);
	points->Delete();
	polys->Delete();

	return result;
}
class QT_VTK_EXPORT qt_vtk
{
public:
	static vtkOrientationMarkerWidget* createOrientationMarker(vtkRenderWindowInteractor *iren, bool bOn=true);
	static vtkAxesActor* createAxis();
};

#endif // QT_VTK_H
