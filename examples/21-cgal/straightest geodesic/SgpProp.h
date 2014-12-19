#ifndef SGPPROP_H
#define SGPPROP_H

#include "proplib_global.h"
#pragma warning (push)
#pragma warning(disable : 4244 4267 4503 4819 4996 )
#include <CGAL/Cartesian.h>	
#include <CGAL/Plane_3.h>
//#include <CGAL/Vector_3.h>
#pragma warning (pop)

#include "MyPolyhedron.h"
#include "c2vPolyLines.h"
#include <vtkActor.h>
#include <vtkRenderer.h>
#include <QObject>
#include <list>

typedef CGAL::Cartesian<double>  K;
typedef CGAL::Vector_3<K>  Vector_3;
typedef CGAL::Plane_3<K> Plane_3;

//Prop of straightest geodesic path
//refer to: "Straightest Paths on Meshes by Cutting Planes"
class PROPLIB_EXPORT SgpProp : public QObject
{
	Q_OBJECT
public:
	SgpProp(vtkRenderer* render, QObject *parent);
	~SgpProp();
signals:
	void updateVTK();
	void setStatusBarInfo(QString);
public slots:
	// methods
	void createPlane(Vertex_handle &center, Polyhedron* mesh);
	void createBorderSgp(Vertex_handle &center, Polyhedron* mesh);
	void createSgp(Vertex_handle &center, Polyhedron* mesh, int idx=-1);
	void createBorderPgp(Vertex_handle &center, Polyhedron* mesh);
	void createPgp(Vertex_handle &center, Polyhedron* mesh, int idx=-1);
	void createBorderRK(Vertex_handle &center, Polyhedron* mesh, vtkPolyData* meshData);
	void createRK(Vertex_handle &center, Polyhedron* mesh, vtkPolyData* meshData, int idx=-1);
	void clearActors();
	// properties
	void setVisibility(bool in);
	void setBasePlaneVisibility(bool in);
public:
	c2vPolyLines* getSgpLines(){return m_lines;}
private:
	void createSgp(Vertex_handle &center, Polyhedron* mesh, std::list<Vertex_handle>& vhs);
	void createPgp(Vertex_handle &center, Polyhedron* mesh, std::list<Vertex_handle>& vhs);
	void createRK(Vertex_handle &center, Polyhedron* mesh, vtkPolyData* meshData, std::list<Vertex_handle>& vhs);
	Vector_3 m_basePlaneNormal;
	c2vPolyLines* m_lines;// all straight lines from border to center
	vtkActor* m_basePlaneActor;
	vtkRenderer* m_render;
	std::list<std::list<Point_3> > m_sgps;

	std::list<vtkActor*> m_actors;
};

#endif // SGPPROP_H
