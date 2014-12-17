#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>

typedef CGAL::Simple_cartesian<double>     Kernel;
typedef CGAL::Polyhedron_3<Kernel>         Polyhedron;
typedef Polyhedron::Halfedge_handle        Halfedge_handle;
typedef Polyhedron::Vertex_iterator        Vertex_iterator;
typedef Polyhedron::Point_3                Point_3;

#include <iostream>
#include <sstream>
using namespace std;

void modify_vertex_position()
{
	Polyhedron P;
    Halfedge_handle h = P.make_tetrahedron();
    if ( P.is_tetrahedron(h))
	{
		int i(0);
		for(Vertex_iterator vi = P.vertices_begin(); vi != P.vertices_end(); ++vi,++i)
		{
		  std::cout << "before changing vertex " << i << ": " << vi->point().x() << vi->point().y() << vi->point().z() << endl;
		  Point_3 pt(1, 0, 0);
		  vi->point() = pt;
		  std::cout << "after changing vertex " << i << ": " << vi->point().x() << vi->point().y() << vi->point().z() << endl;
		}
	}
}

int main()
{
	modify_vertex_position();
	return 0;
}