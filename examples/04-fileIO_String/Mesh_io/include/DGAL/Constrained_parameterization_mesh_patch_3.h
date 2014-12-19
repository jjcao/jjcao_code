//Constrained_parameterization_mesh_patch_3.h
#ifndef DGAL_CONSTRAINED_PARAMETERIZATION_MESH_PATCH_H
#define DGAL_CONSTRAINED_PARAMETERIZATION_MESH_PATCH_H

#include <DGAL/config.h>
//#include <DGAL/Parameterization_constrain_policies_3.h>
#include <CGAL/Parameterization_mesh_patch_3.h>

DGAL_BEGIN_NAMESPACE

template
<
	class ParameterizationPatchableMesh_3
>
class Constrained_parameterization_mesh_patch_3
	:public CGAL::Parameterization_mesh_patch_3<ParameterizationPatchableMesh_3>
{
// Private types
private:
	// Superclass
	typedef CGAL::Parameterization_mesh_patch_3<ParameterizationPatchableMesh_3>
		                                             Base;
public:
	typedef ParameterizationPatchableMesh_3 Adaptor;
	typedef typename Adaptor::Constrain_vertex_iterator 
		                                            Constrain_vertex_iterator;
	typedef typename Adaptor::Constrain_vertex_const_iterator
		                                            Constrain_vertex_const_iterator;

//Public operations
public:
	template<class InputIterator>
	Constrained_parameterization_mesh_patch_3(Adaptor& mesh,
                                  InputIterator first_seam_vertex,
                                  InputIterator end_seam_vertex)
		: Base(mesh, first_seam_vertex, end_seam_vertex)
	{
	}
	bool  is_constrained_vertex(Vertex_const_handle vertex) const 
	{
		return vertex->vertex()->is_constrained();
		//get_decorated_mesh().
    }

};

DGAL_END_NAMESPACE
#endif //DGAL_CONSTRAINED_PARAMETERIZATION_MESH_PATCH_H