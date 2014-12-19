#ifndef DGAL_PARAMETERIZATION_CONSTRAIN_POLICIES_3_H
#define DGAL_PARAMETERIZATION_CONSTRAIN_POLICIES_3_H

#include <DGAL/config.h>

DGAL_BEGIN_NAMESPACE

enum ParamConstrainType//ParameterizationConstrainType
{
	NO_CONSTRAIN,
	HARD_CONSTRAIN,
	SOFT_CONSTRAIN
};

//template
//<
//    ParamConstrainType constrainType 
//>
//class Constrian_policy//general
//{
//};

template
<
    class ParameterizationMesh_3,
	class Vector,
	class Matrix
>
class Hard_constrian_policy
{
public:
	typedef ParameterizationMesh_3 Adaptor;
	typedef typename Adaptor::Vertex        Vertex;
    typedef typename Adaptor::Vertex_handle Vertex_handle;
    typedef typename Adaptor::Vertex_const_handle
                                            Vertex_const_handle;
    typedef typename Adaptor::Vertex_iterator Vertex_iterator;
    typedef typename Adaptor::Vertex_const_iterator
                                            Vertex_const_iterator;
	typedef typename Adaptor::Constrain_vertex_iterator
		                                    Constrain_vertex_iterator;
	typedef typename Adaptor::Point_2       Point_2;
protected:
	void  setup_constrained_vertex_relations (Matrix& A, Vector& Bu, Vector& Bv,
											  const Adaptor& mesh,
		                                      Vertex_handle vh)
	{
		bool tmp = mesh.is_vertex_parameterized(vh);
		CGAL_surface_mesh_parameterization_assertion(!tmp);

		// Get vertex index in sparse linear system
		int index = mesh.get_vertex_index(vh);

		// Write a as diagonal coefficient of A
		A.set_coef(index, index, 1);

		// Write constant in Bu and Bv
		Point_2 uv = mesh.get_vertex_uv(vh);
		Bu[index] = uv.x();
		Bv[index] = uv.y();
	}

};

template
<
    class ParameterizationMesh_3,
	class Vector,
	class Matrix
>
class No_constrian_policy
{
public:
	typedef ParameterizationMesh_3 Adaptor;
	typedef typename Adaptor::Vertex        Vertex;
    typedef typename Adaptor::Vertex_handle Vertex_handle;
    typedef typename Adaptor::Vertex_const_handle
                                            Vertex_const_handle;
    typedef typename Adaptor::Vertex_iterator Vertex_iterator;
    typedef typename Adaptor::Vertex_const_iterator
                                            Vertex_const_iterator;
protected:
	void  setup_constrained_vertex_relations (Matrix& A, Vector& Bu, Vector& Bv,
											  const Adaptor& mesh,
		                                      Vertex_handle vh)
	{
	}
};

DGAL_END_NAMESPACE

#endif //DGAL_PARAMETERIZATION_CONSTRAIN_POLICIES_3_H
