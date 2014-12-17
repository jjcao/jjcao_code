#ifndef DGAL_PARAMETERIZATION_MEASURER_3_H
#define DGAL_PARAMETERIZATION_MEASURER_3_H

#include <DGAL/config.h>

DGAL_BEGIN_NAMESPACE

#include <CGAL/Kernel/global_functions.h>

static const double PAI = 3.14159265259;
static const double EPSILON = 0.0000000001;

//stretch distortion
//P. V. Sander, J. Snyder, S. J. Gortler, and H. Hoppe. Texture mapping progressive meshes. 
//In Proceedings of ACMSIGGRAPH 2001, pages 409–416, 2001

//angle distortion
//Equation 1 of SURFACE PARAMETERIZATION FOR MESHING BY TRIANGULATION FLATTENING_00
template<class ParameterizationMesh_3>       //< 3D surface
class Parameterization_measurer_3
{
public:
	typedef ParameterizationMesh_3                             Mesh;
    typedef typename Mesh::Point_2                             Point_2;
    typedef typename Mesh::Point_3                             Point_3;
    typedef typename Mesh::Vector_2                            Vector_2;
    typedef typename Mesh::Vector_3                            Vector_3;
	typedef typename Mesh::Facet_iterator                      Facet_iterator;
	typedef typename Mesh::Facet_handle                        Facet_handle;
	typedef typename Mesh::Vertex_iterator                     Vertex_iterator;
	typedef typename Mesh::Halfedge_iterator                   Halfedge_iterator;
	typedef typename Mesh::Vertex_handle                       Vertex_handle;
	typedef typename Mesh::Halfedge_handle                     Halfedge_handle;
	typedef typename Mesh::Halfedge_around_facet_circulator    Halfedge_facet_circulator;
	typedef typename Mesh::Halfedge_around_vertex_circulator   Halfedge_vertex_circulator;
public:
	Parameterization_measurer_3():m_numFlipFacets(0)
	{}
private:
	int m_numFlipFacets;//has meaning only after calling l_2_distortion or l_infinite_distortion
    double m_angleDistortion;//for the whole mesh
	double m_l2Distortion;//for the whole mesh
	double m_l_infiniteDistortion;//for the whole mesh
	std::list<double> m_as;//angle_distortion per vertex
	std::list<double> m_l2s;//l_2_distortion per face
public:
	std::list<double>& getAngleDistortionPerV(){return m_as;}
	std::list<double>& getL2DistortionPerF(){return m_l2s;}
	double getAngleDistortion(){return m_angleDistortion;}
	double getL2Distortion(){return m_l2Distortion;}
	double getL_infiniteDistortion(){return m_l_infiniteDistortion;}
public:
	void computeStretchDistortion(Mesh* mesh)
	{
		m_l2Distortion = 0;
	    m_l_infiniteDistortion = 0;
		m_l2s.clear();

		m_numFlipFacets = 0;
		double uv_area(0.0);
		double mesh_area(0.0);

		for(Facet_iterator fi = mesh->facets_begin();fi != mesh->facets_end(); ++fi)
		{	
			double areaD; double areaM;
			double a_c; double a_b_c;
			prepare_l(fi->facet_begin(), a_c, a_b_c, areaD, areaM);

			double l2 = std::sqrt( a_c);
			double li = std::sqrt( a_c + a_b_c);
			mesh_area += areaM;
			uv_area += areaD;

			double tmp = l2*l2*areaM;
			m_l2s.push_back(tmp);
			m_l2Distortion += tmp;
			if ( li > m_l_infiniteDistortion) m_l_infiniteDistortion = li;
		}

		int tmp = mesh->size_of_facets() - m_numFlipFacets;
		if ( tmp>0)
		{		
			if (tmp<m_numFlipFacets)
				m_numFlipFacets = tmp;
		}
		else
		{
			m_numFlipFacets = 0;
		}

		m_l2Distortion = std::sqrt(m_l2Distortion/mesh_area);
		m_l2Distortion *= std::sqrt(uv_area/mesh_area);

		m_l_infiniteDistortion = std::sqrt(m_l_infiniteDistortion/mesh_area);
		m_l_infiniteDistortion *= std::sqrt(uv_area/mesh_area);
		
		if (m_numFlipFacets)
		{
			m_l2Distortion = -m_l2Distortion;
			m_l_infiniteDistortion = -m_l_infiniteDistortion;
		}
	}
	double computeAngleDistortion(Mesh* mesh)
	{
		m_angleDistortion = 0;
		m_as.clear();

		for (Vertex_iterator vi = mesh->vertices_begin(); vi != mesh->vertices_end(); ++vi)
		{
			std::list<double> anglesM; std::list<double> anglesD;
			double angle = computeAllAngles(vi, anglesM,anglesD);

			std::list<double>::iterator itm = anglesM.begin();
			std::list<double>::iterator itd = anglesD.begin();
			double alpha, phi,omega;

			int j = vi->index();
			double tmp = 0;
			if ( mesh->is_border(vi))
			{		
				for(; itm!=anglesM.end(); ++itm,++itd)
				{
					alpha = *itd;
					phi = *itm;// the only diff line between border and inner vertex
					if (phi!=0)
					{
						omega = 1/(phi*phi);
						tmp += (alpha-phi)*(alpha-phi)* omega;
					}
				}
				m_as.push_back(tmp);
				m_angleDistortion += tmp;
			}
			else
			{
				for(; itm!=anglesM.end(); ++itm,++itd)
				{
					alpha = *itd;
					phi = (*itm)*2*PAI/angle;// the only diff line between border and inner vertex
					omega = 1/(phi*phi);
					tmp += (alpha-phi)*(alpha-phi)* omega;					
				}
				m_as.push_back(tmp);
				m_angleDistortion += tmp;
			}
		}

		m_angleDistortion = m_angleDistortion/mesh->size_of_vertices();
		return m_angleDistortion;
	}
	double angular_distortion(Mesh* mesh)//from Graphite 2.1
	{
		double distortion(0);
		long count(0);
		double a3d(0), a2d(0);
		for (Halfedge_iterator hi = mesh->halfedges_begin(); hi != mesh->halfedges_end(); ++hi)
		{
			if(!hi->is_border())
			{				
				angle(hi, hi->next(), a2d, a3d);				
				distortion += (a3d - a2d)*(a3d - a2d) ;
				++count;		
			}			
		}
		
		distortion /= double(count);
		return distortion;			
	}
	//outdated call, use computeStretchDistoriton instead of this function
	double l_2_distortion(Mesh* mesh, bool is_normalized = true)
	{
		m_numFlipFacets = 0;
		double result(0.0);	
		double uv_area(0.0);
		double mesh_area(0.0);

		for(Facet_iterator pFacet = mesh->facets_begin();pFacet != mesh->facets_end(); ++pFacet)
		{	
			double area(0.0); double area_1(0.0);
			double l = l_2(pFacet->facet_begin(), area, area_1);
			
			mesh_area += area_1;
			uv_area += area;
			result += l*l*area_1;
		}

		int tmp = mesh->size_of_facets() - m_numFlipFacets;
		if ( tmp>0)
		{		
			if (tmp<m_numFlipFacets)
				m_numFlipFacets = tmp;
		}
		else
		{
			m_numFlipFacets = 0;
		}

		if (m_numFlipFacets)
		{
			result = 99999.9;
		}
		else
		{
			result = std::sqrt(result/mesh_area);
			if ( is_normalized)	result *= std::sqrt(uv_area/mesh_area);
		}

		return result;
	}
	//outdated call, use computeStretchDistoriton instead of this function
	double l_infinite_distortion(Mesh* mesh, bool is_normalized = true)
	{
		m_numFlipFacets = 0;
		double result(0.0);	
		double uv_area(0.0);
		double mesh_area(0.0);

		for(Facet_iterator pFacet = mesh->facets_begin();pFacet != mesh->facets_end(); ++pFacet)
		{	
			double area(0.0); double area_1(0.0);
			double l = l_infinite(pFacet ->facet_begin(), area, area_1);
			
			mesh_area += area_1;
			uv_area += area;

			if ( l > result) result = l;
		}

		/*if ( (long(mesh_area + 0.5) - mesh_area) < EPSILON)
			mesh_area = double( long(mesh_area + 0.5));
		if ( (long(uv_area + 0.5) - uv_area) < EPSILON)
			uv_area = double( long(uv_area + 0.5));*/

		int tmp = mesh->size_of_facets() - m_numFlipFacets;
		if ( tmp>0)
		{		
			if (tmp<m_numFlipFacets)
				m_numFlipFacets = tmp;
		}
		else
		{
			m_numFlipFacets = 0;
		}
		if (m_numFlipFacets)
		{
			result = 99999.9;
		}
		else
		{
			if ( is_normalized)	result *= std::sqrt(uv_area/mesh_area);
		}

		return result;
	}
	//outdated call, use computeAngleDistoriton instead of this function
	double angle_distortion(Mesh* mesh, bool is_normalized = true)
	{
		double result(0.0);
		for(Facet_iterator pFacet = mesh->facets_begin();pFacet != mesh->facets_end(); ++pFacet)
		{
			Halfedge_facet_circulator hfc = pFacet ->facet_begin();
			do {//访问pos的3个顶点
				Vertex_handle vh = hfc->vertex();
				//compute alpha in parameter domain of mesh
				Point_3 p(vh->u(),vh->v(),0.0);
				Halfedge_facet_circulator tc(hfc);
				++tc;
				Point_3 pl( tc->vertex()->u(),tc->vertex()->v(),0.0);
				++tc;
				Point_3 pr( tc->vertex()->u(),tc->vertex()->v(),0.0);
				double alpha = compute_angle_rad(pl, p, pr);

				//compute phi in mesh
				double phi = compute_phi(vh, pFacet, mesh->is_border(vh) );
				double omega = 1/(phi*phi);

				result += (alpha-phi)*(alpha-phi)* omega;

			} while ( ++hfc != pFacet->facet_begin());		
		}

		if ( is_normalized)	result = result/mesh->size_of_facets()/3.0;

		return result;
	}
	int getNumFlipFacets(){return m_numFlipFacets;}
protected:
	//
	////p for uv; q for original mesh
	void prepare_l(Halfedge_facet_circulator hfc, double& a_c, double& a_b_c, double& uv_area, double& mesh_area)
	{
		Vertex_handle vh = hfc->vertex();
		Vector_3 q1( vh->point().x(), vh->point().y(), vh->point().z()); 
		Point_2 p1( vh->u(), vh->v());
		++hfc; vh = hfc->vertex();
		Vector_3 q2( vh->point().x(), vh->point().y(), vh->point().z()); 
		Point_2 p2( vh->u(), vh->v());
		++hfc; vh = hfc->vertex();
		Vector_3 q3( vh->point().x(), vh->point().y(), vh->point().z()); 
		Point_2 p3( vh->u(), vh->v());

		uv_area = CGAL::area(p1, p2, p3);
		if ( uv_area<0)
			++m_numFlipFacets;

		uv_area = std::abs( uv_area);
		double two_a = uv_area * 2.0;
		Vector_3 s_partial_s = q1*(p2.y()-p3.y()) + q2*(p3.y()-p1.y()) + q3*(p1.y()-p2.y());
		s_partial_s = s_partial_s/two_a;
		Vector_3 s_partial_t = q1*(p3.x()-p2.x()) + q2*(p1.x()-p3.x()) + q3*(p2.x()-p1.x());
		s_partial_t = s_partial_t/two_a;

		double a = s_partial_s * s_partial_s;
		double b = s_partial_s * s_partial_t;
		double c = s_partial_t * s_partial_t;
		a_c = 0.5*(a + c);
		a_b_c = 0.5*std::sqrt( (a-c)*(a-c) + 4*b*b );
		//double tau = std::sqrt( a_c + a_b_c); // max singular value
		//double gamma = std::sqrt( a_c - a_b_c); // min singular value

		Vector_3 v1(q2 - q1); Vector_3 v2(q3 - q2);
		Vector_3 tmp = CGAL::cross_product(v1, v2);
		mesh_area = 0.5 * std::sqrt(tmp * tmp);
	}	
	double l_2(Halfedge_facet_circulator hfc, double& uv_area, double& mesh_area){
		double a_c; double a_b_c;
		prepare_l(hfc, a_c, a_b_c, uv_area, mesh_area);
		return std::sqrt( a_c); 
	}	
	double l_infinite(Halfedge_facet_circulator hfc, double& uv_area, double& mesh_area){
		double a_c; double a_b_c;
		prepare_l(hfc, a_c, a_b_c, uv_area, mesh_area);
		return std::sqrt( a_c + a_b_c);
	}

	/// Preconditions:
	///		h1->next() == h2
	void angle(Halfedge_handle h1, Halfedge_handle h2, double &a2d, double &a3d)
	{
		Vertex_handle vh = h1->opposite()->vertex();
		Point_3 p3l = vh->point();			
		Point_3 p2l(vh->u(),vh->v(),0.0);
		
		vh = h1->vertex();
		Point_3 p3 = vh->point();			
		Point_3 p2(vh->u(),vh->v(),0.0);		

		vh = h2->vertex();
		Point_3 p3r = vh->point();			
		Point_3 p2r(vh->u(),vh->v(),0.0);		

		a3d = compute_angle_rad(p3l, p3, p3r);
		a2d = compute_angle_rad(p2l, p2, p2r);	
	}
	
	//anglesM angles around vh on the mesh
	//anglesD angles around vh on the domain
	//return total angles aournd the vh on the mesh
	double computeAllAngles(Vertex_handle vh, std::list<double> &anglesM,std::list<double> &anglesD)
	{
		double result(0.0);
		
		Point_3 p = vh->point();
		Point_3 pd(vh->u(),vh->v(),0.0);
		Halfedge_vertex_circulator hvc = vh->vertex_begin();
		Halfedge_vertex_circulator end = hvc;
		CGAL_For_all(hvc, end)
		{
			Vertex_handle vht = hvc->opposite()->vertex();
			Point_3 pl = vht->point();			
			Point_3 pdl(vht->u(),vht->v(),0.0);

			vht = hvc->next()->vertex();			
			Point_3 pr = vht->point();
			Point_3 pdr(vht->u(),vht->v(),0.0);				

			double alphaM = compute_angle_rad(pl, p, pr);
			double alphaD = compute_angle_rad(pdl, pd, pdr);
			anglesM.push_back(alphaM);
			anglesD.push_back(alphaD);

			result += alphaM;
		}

		return result;
	}
	//
	//vh is a vertex of the face fh.
	double compute_phi(Vertex_handle vh, Facet_handle fh, bool is_border){
		double result(0.0);
		double beta(0.0);

		Halfedge_vertex_circulator hvc = vh->vertex_begin();
		do {
			Facet_handle cfh = hvc->facet();

			//(pl,p,pr)
			Halfedge_vertex_circulator hvc1 = hvc;
			++hvc1;
			Point_3 pl = hvc1->opposite()->vertex()->point();
			Point_3 p = vh->point();
			Point_3 pr = hvc->opposite()->vertex()->point();
			double alpha = compute_angle_rad(pl, p, pr);

			if ( cfh == fh){
				beta = alpha;
				if ( is_border) return beta;
			}
			result += alpha;

		} while ( ++hvc != vh->vertex_begin());

		result = (beta * 2 * PAI)/result;

		return result;
	}
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
        if(cross_norm != 0.0)
            return (dot/cross_norm);
		else{
			std::cerr << "cross norm = 0.0; lead to undefined behavior!" ;
            return 0.0; // undefined
		}
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
        Vector_3 cross_vector = CGAL::cross_product(u,v);
        double cross_norm = std::sqrt(cross_vector*cross_vector);
        if(dot != 0.0)
            return (cross_norm/dot);
		else{
			std::cerr << "denominator dot = 0.0; lead to undefined behavior!" ;
            return 0.0; // undefined
		}
    }

	//Triangle_3
    //                                                     -> ->
    /// Return angle (in radians) of of (P,Q,R) corner (ie QP,QR angle).
    double compute_angle_rad(const Point_3& P,
                                    const Point_3& Q,
                                    const Point_3& R)
    {
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
            return PAI-std::asin(fix_sine(AbsSine));
    }

// Private operations
private:
    /// Fix sine.
    double fix_sine(double sine)
    {
        if(sine >= 1)
            return 1;
        else if(sine <= -1)
            return -1;
        else
            return sine;
    }
};

DGAL_END_NAMESPACE

#endif//DGAL_PARAMETERIZATION_MEASURER_3_H