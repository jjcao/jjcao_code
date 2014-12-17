#ifndef DGAL_PARAMETERIZATION_POLICIES_3_H
#define DGAL_PARAMETERIZATION_POLICIES_3_H

#include <DGAL/config.h>

DGAL_BEGIN_NAMESPACE

template<class ParameterizationMesh_3>       //< 3D surface
class Parameterization_policy_3
{
// Public types
public:
    typedef ParameterizationMesh_3          Adaptor;

// Private types
protected:
    // Mesh_Adaptor_3 subtypes:
    typedef typename Adaptor::NT            NT;
    typedef typename Adaptor::Point_2       Point_2;
    typedef typename Adaptor::Point_3       Point_3;
    typedef typename Adaptor::Vector_2      Vector_2;
    typedef typename Adaptor::Vector_3      Vector_3;

// Public operations
public:
	Parameterization_policy_3(int type):m_type(type){}
    /// Destructor of base class should be virtual.
    virtual ~Parameterization_policy_3() {}

	int getType(){return m_type;}
    //                                                  -> ->
    /// Return cotangent of (P,Q,R) corner (ie cotan of QP,QR angle).
    double cotangent(const Point_3& P,
                     const Point_3& Q,
                     const Point_3& R)
    {
        Vector_3 u = P - Q;
        Vector_3 v = R - Q;
        // (u . v)/((u x v).len)
        double dot = (u*v);
        Vector_3 cross_vector = CGAL::cross_product(u,v);
        double cross_norm = std::sqrt(cross_vector*cross_vector);
        //CGAL_surface_mesh_parameterization_assertion(cross_norm != 0.0);
        if(cross_norm != 0.0)
            return (dot/cross_norm);
        else
            return 0.0; // undefined
    }
	double cotangent(const Point_2& P,
                     const Point_2& Q,
                     const Point_2& R)
    {
		Point_3 p(P.x(), P.y(), 0.0);
		Point_3 q(Q.x(), Q.y(), 0.0);
		Point_3 r(R.x(), R.y(), 0.0);

		return cotangent(p,q,r);
    }

    //                                                  -> ->
    /// Return tangent of (P,Q,R) corner (ie tangent of QP,QR angle).
    double tangent(const Point_3& P,
                   const Point_3& Q,
                   const Point_3& R)
    {
        Vector_3 u = P - Q;
        Vector_3 v = R - Q;
        // (u . v)/((u x v).len)
        double dot = (u*v);
        CGAL_surface_mesh_parameterization_assertion(dot != 0.0);
        Vector_3 cross_vector = CGAL::cross_product(u,v);
        double cross_norm = std::sqrt(cross_vector*cross_vector);
        if(dot != 0.0)
            return (cross_norm/dot);
        else
            return 0.0; // undefined
    }

    //                                                     -> ->
    /// Return angle (in radians) of of (P,Q,R) corner (ie QP,QR angle).
    static double compute_angle_rad(const Point_3& P,
                                    const Point_3& Q,
                                    const Point_3& R)
    {
        static const double PI = 3.14159265359;

        Vector_3 u = P - Q;
        Vector_3 v = R - Q;

        // check
        double product = std::sqrt(u*u) * std::sqrt(v*v);
        if(product == 0)
            return 0.0;

        // cosine
        double dot = (u*v);
        double cosine = dot / product;

        // sine
        Vector_3 w = CGAL::cross_product(u,v);
        double AbsSine = std::sqrt(w*w) / product;

        if(cosine >= 0)
            return std::asin(fix_sine(AbsSine));
        else
            return PI-std::asin(fix_sine(AbsSine));
    }
/// @endcond

// Private operations
private:
    /// Fix sine.
    static double fix_sine(double sine)
    {
        if(sine >= 1)
            return 1;
        else if(sine <= -1)
            return -1;
        else
            return sine;
    }
	int m_type;
};

//Discrete conformal mapping
template
<
    class ParameterizationMesh_3
>
class DCP_policy
	: public Parameterization_policy_3<ParameterizationMesh_3>
{
public:
	DCP_policy():Parameterization_policy_3(1){}
    typedef ParameterizationMesh_3          Adaptor;
// Private types
private:
    // Mesh_Adaptor_3 subtypes:
    typedef typename Adaptor::NT            NT;
    typedef typename Adaptor::Point_3       Point_3;
    typedef typename Adaptor::Vertex        Vertex;
    typedef typename Adaptor::Vertex_handle Vertex_handle;
    typedef typename Adaptor::Vertex_const_handle
                                            Vertex_const_handle;
    typedef typename Adaptor::Vertex_around_vertex_circulator
                                            Vertex_around_vertex_circulator;
    typedef typename Adaptor::Vertex_around_vertex_const_circulator
                                            Vertex_around_vertex_const_circulator;
public:
	/// Compute w_ij = (i, j) coefficient of matrix A for j neighbor vertex of i.
    /// Implementation note: Subclasses must at least implement compute_w_ij().
    NT compute_w_ij(Adaptor& mesh,
                            Vertex_handle main_vertex_v_i,
							Vertex_around_vertex_circulator neighbor_vertex_v_j)
	{		
		Point_3 position_v_i = mesh.get_vertex_position(main_vertex_v_i);
        Vertex_handle vh = neighbor_vertex_v_j;
        Point_3 position_v_j = mesh.get_vertex_position(vh);

        // Compute cotangent of (v_i,v_k,v_j) corner (ie cotan of v_k corner)
        // if v_k is the vertex before v_j when circulating around v_i
        Vertex_around_vertex_circulator previous_vertex_v_k = neighbor_vertex_v_j;
        previous_vertex_v_k --;
        Point_3 position_v_k = mesh.get_vertex_position(Vertex_handle(previous_vertex_v_k));
        double cotg_beta_ij  = cotangent(position_v_i, position_v_k, position_v_j);

        // Compute cotangent of (v_j,v_l,v_i) corner (ie cotan of v_l corner)
        // if v_l is the vertex after v_j when circulating around v_i
        Vertex_around_vertex_circulator next_vertex_v_l = neighbor_vertex_v_j;
        next_vertex_v_l ++;
        Point_3 position_v_l = mesh.get_vertex_position(Vertex_handle(next_vertex_v_l));
        double cotg_alpha_ij = cotangent(position_v_j, position_v_l, position_v_i);

        double weight = cotg_beta_ij+cotg_alpha_ij;
        return weight;
	}
	NT compute_w_ij_uv(Adaptor& mesh,
                            Vertex_handle main_vertex_v_i,
							Vertex_around_vertex_circulator neighbor_vertex_v_j)
	{		
		Point_2 position_v_i = mesh.get_vertex_uv(main_vertex_v_i);
        Vertex_handle vh = neighbor_vertex_v_j;
        Point_2 position_v_j = mesh.get_vertex_uv(vh);

        // Compute cotangent of (v_i,v_k,v_j) corner (ie cotan of v_k corner)
        // if v_k is the vertex before v_j when circulating around v_i
        Vertex_around_vertex_circulator previous_vertex_v_k = neighbor_vertex_v_j;
        previous_vertex_v_k --;
        Point_2 position_v_k = mesh.get_vertex_uv(Vertex_handle(previous_vertex_v_k));
        double cotg_beta_ij  = cotangent(position_v_i, position_v_k, position_v_j);

        // Compute cotangent of (v_j,v_l,v_i) corner (ie cotan of v_l corner)
        // if v_l is the vertex after v_j when circulating around v_i
        Vertex_around_vertex_circulator next_vertex_v_l = neighbor_vertex_v_j;
        next_vertex_v_l ++;
        Point_2 position_v_l = mesh.get_vertex_uv(Vertex_handle(next_vertex_v_l));
        double cotg_alpha_ij = cotangent(position_v_j, position_v_l, position_v_i);

        double weight = cotg_beta_ij+cotg_alpha_ij;
        return weight;
	}
};

//Discrete authalic mapping
template
<
    class ParameterizationMesh_3
>
class DAP_policy
	: public Parameterization_policy_3<ParameterizationMesh_3>
{
public:
    typedef ParameterizationMesh_3          Adaptor;
// Private types
private:
    // Mesh_Adaptor_3 subtypes:
    typedef typename Adaptor::NT            NT;
    typedef typename Adaptor::Point_3       Point_3;
    typedef typename Adaptor::Vertex        Vertex;
    typedef typename Adaptor::Vertex_handle Vertex_handle;
    typedef typename Adaptor::Vertex_const_handle
                                            Vertex_const_handle;
    typedef typename Adaptor::Vertex_around_vertex_circulator
                                            Vertex_around_vertex_circulator;
    typedef typename Adaptor::Vertex_around_vertex_const_circulator
                                            Vertex_around_vertex_const_circulator;
public:
	/// Compute w_ij = (i, j) coefficient of matrix A for j neighbor vertex of i.
    /// Implementation note: Subclasses must at least implement compute_w_ij().
    NT compute_w_ij(const Adaptor& mesh,
                            Vertex_handle main_vertex_v_i,
							Vertex_around_vertex_circulator neighbor_vertex_v_j)
	{
		Point_3 position_v_i = mesh.get_vertex_position(main_vertex_v_i);
        Vertex_handle vh = neighbor_vertex_v_j;
        Point_3 position_v_j = mesh.get_vertex_position(vh);

        // Compute the square norm of v_j -> v_i vector
        Vector_3 edge = position_v_i - position_v_j;
        double square_len = edge*edge;

        // Compute cotangent of (v_k,v_j,v_i) corner (ie cotan of v_j corner)
        // if v_k is the vertex before v_j when circulating around v_i
        Vertex_around_vertex_circulator previous_vertex_v_k = neighbor_vertex_v_j;
        previous_vertex_v_k --;
        Point_3 position_v_k = mesh.get_vertex_position(Vertex_handle(previous_vertex_v_k));
        double cotg_psi_ij  = cotangent(position_v_k, position_v_j, position_v_i);

        // Compute cotangent of (v_i,v_j,v_l) corner (ie cotan of v_j corner)
        // if v_l is the vertex after v_j when circulating around v_i
        Vertex_around_vertex_circulator next_vertex_v_l = neighbor_vertex_v_j;
        next_vertex_v_l ++;
        Point_3 position_v_l = mesh.get_vertex_position(Vertex_handle(next_vertex_v_l));
        double cotg_theta_ij = cotangent(position_v_i, position_v_j, position_v_l);

        double weight = 0.0;
        CGAL_surface_mesh_parameterization_assertion(square_len != 0.0);    // two points are identical!
        if(square_len != 0.0)
            weight = (cotg_psi_ij+cotg_theta_ij)/square_len;

        return weight;
	}
};

//Mean value coordinate mapping of Floater
template
<
    class ParameterizationMesh_3
>
class MVC_policy
	: public Parameterization_policy_3<ParameterizationMesh_3>
{
public:
	MVC_policy():Parameterization_policy_3(0){}
    typedef ParameterizationMesh_3          Adaptor;
// Private types
private:
    // Mesh_Adaptor_3 subtypes:
    typedef typename Adaptor::NT            NT;
    typedef typename Adaptor::Point_3       Point_3;
    typedef typename Adaptor::Vertex        Vertex;
    typedef typename Adaptor::Vertex_handle Vertex_handle;
    typedef typename Adaptor::Vertex_const_handle
                                            Vertex_const_handle;
    typedef typename Adaptor::Vertex_around_vertex_circulator
                                            Vertex_around_vertex_circulator;
    typedef typename Adaptor::Vertex_around_vertex_const_circulator
                                            Vertex_around_vertex_const_circulator;
// Protected operations
public:
    /// Compute w_ij = (i, j) coefficient of matrix A for j neighbor vertex of i.
    virtual NT compute_w_ij(Adaptor& mesh,
                            Vertex_handle main_vertex_v_i,
                            Vertex_around_vertex_circulator neighbor_vertex_v_j)
    {
        Point_3 position_v_i = mesh.get_vertex_position(main_vertex_v_i);
		Vertex_handle vh = neighbor_vertex_v_j;
        Point_3 position_v_j = mesh.get_vertex_position(vh);

        // Compute the norm of v_j -> v_i vector
        Vector_3 edge = position_v_i - position_v_j;
        double len = std::sqrt(edge*edge);

        // Compute angle of (v_j,v_i,v_k) corner (ie angle of v_i corner)
        // if v_k is the vertex before v_j when circulating around v_i
        Vertex_around_vertex_circulator previous_vertex_v_k = neighbor_vertex_v_j;
        previous_vertex_v_k --;
        Point_3 position_v_k = mesh.get_vertex_position(Vertex_handle(previous_vertex_v_k));
        double gamma_ij  = compute_angle_rad(position_v_j, position_v_i, position_v_k);

        // Compute angle of (v_l,v_i,v_j) corner (ie angle of v_i corner)
        // if v_l is the vertex after v_j when circulating around v_i
        Vertex_around_vertex_circulator next_vertex_v_l = neighbor_vertex_v_j;
        next_vertex_v_l ++;
        Point_3 position_v_l = mesh.get_vertex_position(Vertex_handle(next_vertex_v_l));
        double delta_ij = compute_angle_rad(position_v_l, position_v_i, position_v_j);

        double weight = 0.0;
        CGAL_surface_mesh_parameterization_assertion(len != 0.0);    // two points are identical!
        if(len != 0.0)
            weight = (std::tan(0.5*gamma_ij) + std::tan(0.5*delta_ij)) / len;
        CGAL_surface_mesh_parameterization_assertion(weight > 0);

        return weight;
    }
};

//Barycentric mapping of Floater
template
<
    class ParameterizationMesh_3
>
class Tutte_policy
	: public Parameterization_policy_3<ParameterizationMesh_3>
{
public:
	Tutte_policy():Parameterization_policy_3(2){}
    typedef ParameterizationMesh_3          Adaptor;
// Private types
private:
    // Mesh_Adaptor_3 subtypes:
    typedef typename Adaptor::NT            NT;
    typedef typename Adaptor::Point_3       Point_3;
    typedef typename Adaptor::Vertex        Vertex;
    typedef typename Adaptor::Vertex_handle Vertex_handle;
    typedef typename Adaptor::Vertex_const_handle
                                            Vertex_const_handle;
    typedef typename Adaptor::Vertex_around_vertex_circulator
                                            Vertex_around_vertex_circulator;
    typedef typename Adaptor::Vertex_around_vertex_const_circulator
                                            Vertex_around_vertex_const_circulator;
// Protected operations
public:
    /// Compute w_ij = (i, j) coefficient of matrix A for j neighbor vertex of i.
    virtual NT compute_w_ij(Adaptor& mesh,
                            Vertex_handle main_vertex_v_i,
                            Vertex_around_vertex_circulator neighbor_vertex_v_j)
    {
        Point_3 position_v_i = mesh.get_vertex_position(main_vertex_v_i);
        
        Vertex_around_vertex_circulator vvc = neighbor_vertex_v_j;
        int i (0);

		do {//访问一环顶点
			++i;
		} while ( ++vvc != neighbor_vertex_v_j);

        double weight = 1.0/i;

        return weight;
    }
};

template
<
    class ParameterizationMesh_3
>
class Spring_policy
	: public Parameterization_policy_3<ParameterizationMesh_3>
{
public:
	Spring_policy():Parameterization_policy_3(3){}
    typedef ParameterizationMesh_3          Adaptor;
// Private types
private:
    // Mesh_Adaptor_3 subtypes:
    typedef typename Adaptor::NT            NT;
    typedef typename Adaptor::Point_3       Point_3;
    typedef typename Adaptor::Vertex        Vertex;
    typedef typename Adaptor::Vertex_handle Vertex_handle;
    typedef typename Adaptor::Vertex_const_handle
                                            Vertex_const_handle;
	typedef typename Adaptor::Halfedge_handle
                                            Halfedge_handle;
	typedef typename Adaptor::Halfedge_const_handle
                                            Halfedge_const_handle;	
    typedef typename Adaptor::Vertex_around_vertex_circulator
                                            Vertex_around_vertex_circulator;
    typedef typename Adaptor::Vertex_around_vertex_const_circulator
                                            Vertex_around_vertex_const_circulator;
	typedef typename Adaptor::Halfedge_around_vertex_circulator
												Halfedge_around_vertex_circulator;
// Protected operations
public:
    /// Compute w_ij = (i, j) coefficient of matrix A for j neighbor vertex of i.
    virtual NT compute_w_ij(Adaptor& mesh,
                            Vertex_handle main_vertex_v_i,
                            Vertex_around_vertex_circulator neighbor_vertex_v_j)
    {
        Point_3 position_v_i = mesh.get_vertex_position(main_vertex_v_i);
		Vertex_handle vh = neighbor_vertex_v_j;
        Point_3 position_v_j = mesh.get_vertex_position(vh);
		 
		double tmp = main_vertex_v_i->r_;
		Halfedge_handle heh = mesh.get_halfedge(main_vertex_v_i,neighbor_vertex_v_j);
		//Halfedge_handle hh = static_cast<Halfedge_handle>(heh);
		tmp = heh->distance();
		/*Halfedge_around_vertex_circulator cir = main_vertex_v_i->vertex_begin(),
										cir_end = cir;
		CGAL_For_all(cir, cir_end)
			cir->distance();*/

        Vector_3 edge = position_v_i - position_v_j;
        double len = std::sqrt(edge*edge);//edge*edge;//

        double weight = 0.0;
        CGAL_surface_mesh_parameterization_assertion(len != 0.0);    // two points are identical!
        if(len != 0.0)
            weight = 1.0/len;
        CGAL_surface_mesh_parameterization_assertion(weight > 0);

        return weight;
    }
};


DGAL_END_NAMESPACE

#endif //DGAL_PARAMETERIZATION_POLICIES_3_H
