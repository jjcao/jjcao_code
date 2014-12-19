// ======================================================================
//
// Copyright (c) 2002 SurfLab of CISE of University of Florida
//
// File          : libs/src/cgalExt/Polyhedron_subdivision_rules.h
// Description   : Provides the subdivision rules used to template the 
//                 subdivision functions.
// Creation_date : 29 Jan 2002
// Author(s)     : Le-Jeng Shiue <sle-jeng@cise.ufl.edu>
//
// ======================================================================

// $Id: Polyhedron_subdivision_rules.h,v 1.13 2004/05/25 21:23:48 sle-jeng Exp $

#ifndef _POLYHEDRON_SUBDIVISION_RULES_H_01292002
#define _POLYHEDRON_SUBDIVISION_RULES_H_01292002

#include <CGAL/circulator.h>
#include <CGAL/Vector_3.h>
#include <DGAL/config.h>

DGAL_BEGIN_NAMESPACE

//#include <point_2.h>
// ======================================================================
template <class _Poly>
class quadralize_rule {
public:
	typedef _Poly                                         Polyhedron;

	typedef typename Polyhedron::Vertex_handle           Vertex_handle;
	typedef typename Polyhedron::Halfedge_handle         Halfedge_handle;
	typedef typename Polyhedron::Facet_handle            Facet_handle;
	typedef typename Polyhedron::Point_3                 Point;
	//typedef typename Polyhedron::Point_2                 Point_2;


	/**@name Class Methods */
	//@{
public:
	void face_point_rule(Facet_handle, Point&) {};
	void edge_point_rule(Halfedge_handle, Point&) {};
	void vertex_point_rule(Vertex_handle, Point&) {};

	void border_point_rule(Halfedge_handle, Point&, Point&) {};
	//@}
};


// ======================================================================
///
template <class _Poly>
class average_rule : public quadralize_rule<_Poly> {
public:
	typedef _Poly                                         Polyhedron;

	typedef typename Polyhedron::Traits                  Traits;
	typedef typename Traits::Kernel                      Kernel;

	typedef typename Polyhedron::Vertex_handle           Vertex_handle;
	typedef typename Polyhedron::Halfedge_handle         Halfedge_handle;
	typedef typename Polyhedron::Facet_handle            Facet_handle;

	typedef typename Polyhedron::Halfedge_around_facet_circulator  
		Halfedge_around_facet_circulator;
	typedef typename Polyhedron::Halfedge_around_vertex_circulator 
		Halfedge_around_vertex_circulator;

	typedef typename Polyhedron::Point_3                 Point;
	//typedef typename Polyhedron::Point_2                 Point_UV;

	typedef typename Kernel::FT                          FT;

	/**@name Class Methods */
	//@{
public:
	void face_point_rule(Facet_handle facet, Point& pt) {
#ifdef USE_CGAL_POINT
		Halfedge_around_facet_circulator hcir = facet->facet_begin();
		Vector vec = hcir->vertex()->point() - CGAL::ORIGIN;
		++hcir;
		do {
			vec = vec + hcir->vertex()->point();
		} while (++hcir != facet->facet_begin());
		pt = CGAL::ORIGIN + vec/circulator_size(hcir);
#else
		Halfedge_around_facet_circulator hcir = facet->facet_begin();
		int n = 0;
		FT p[] = {0,0,0};
		do {
			Point t = hcir->vertex()->point();
			p[0] += t[0], p[1] += t[1], p[2] += t[2]; 
			++n;
		} while (++hcir != facet->facet_begin());
		pt = Point(p[0]/n, p[1]/n, p[2]/n);
#endif
	}

	void edge_point_rule(Halfedge_handle edge, Point& pt) {
		Point p1 = edge->vertex()->point();
		Point p2 = edge->opposite()->vertex()->point();
		pt = Point((p1[0]+p2[0])/2, (p1[1]+p2[1])/2, (p1[2]+p2[2])/2);
	}

	void vertex_point_rule(Vertex_handle vertex, Point& pt) {
		pt = vertex->point();
	}

	void border_point_rule(Halfedge_handle edge, Point& ept, Point& vpt) {
		edge_point_rule(edge, ept);
	}
	//@}
};

// ======================================================================
///
template <class _Poly>
class CatmullClark_rule : public average_rule<_Poly> {
public:
	typedef _Poly                                        Polyhedron;

	typedef typename Polyhedron::Vertex_handle           Vertex_handle;
	typedef typename Polyhedron::Halfedge_handle         Halfedge_handle;
	typedef typename Polyhedron::Facet_handle            Facet_handle;

	typedef typename Polyhedron::Halfedge_around_facet_circulator  
		Halfedge_around_facet_circulator;
	typedef typename Polyhedron::Halfedge_around_vertex_circulator 
		Halfedge_around_vertex_circulator;

	typedef typename Polyhedron::Traits                  Traits;
	typedef typename Traits::Kernel                      Kernel;

	typedef typename Polyhedron::Point_3                 Point;

	/**@name Class Methods */
	//@{
public:
	void edge_point_rule(Halfedge_handle edge, Point& pt) {
		Point p1 = edge->vertex()->point();
		Point p2 = edge->opposite()->vertex()->point();
		Point f1, f2;
		face_point_rule(edge->facet(), f1);
		face_point_rule(edge->opposite()->facet(), f2);
		pt = Point((p1[0]+p2[0]+f1[0]+f2[0])/4,
			(p1[1]+p2[1]+f1[1]+f2[1])/4,
			(p1[2]+p2[2]+f1[2]+f2[2])/4 );
	}

	void vertex_point_rule(Vertex_handle vertex, Point& pt) {
		Halfedge_around_vertex_circulator vcir = vertex->vertex_begin();
		int n = circulator_size(vcir);    

		float Q[] = {0.0, 0.0, 0.0}, R[] = {0.0, 0.0, 0.0};
		Point& S = vertex->point();

		Point q;
		for (int i = 0; i < n; i++, ++vcir) {
			Point& p2 = vcir->opposite()->vertex()->point();
			R[0] += (S[0]+p2[0])/2;
			R[1] += (S[1]+p2[1])/2;
			R[2] += (S[2]+p2[2])/2;
			face_point_rule(vcir->facet(), q);
			Q[0] += q[0];      
			Q[1] += q[1];      
			Q[2] += q[2];
		}
		R[0] /= n;    R[1] /= n;    R[2] /= n;
		Q[0] /= n;    Q[1] /= n;    Q[2] /= n;

		pt = Point((Q[0] + 2*R[0] + S[0]*(n-3))/n,
			(Q[1] + 2*R[1] + S[1]*(n-3))/n,
			(Q[2] + 2*R[2] + S[2]*(n-3))/n );
	}

	void border_point_rule(Halfedge_handle edge, Point& ept, Point& vpt) {
		Point& ep1 = edge->vertex()->point();
		Point& ep2 = edge->opposite()->vertex()->point();
		ept = Point((ep1[0]+ep2[0])/2, (ep1[1]+ep2[1])/2, (ep1[2]+ep2[2])/2);

		Halfedge_around_vertex_circulator vcir = edge->vertex_begin();
		Point& vp1  = vcir->opposite()->vertex()->point();
		Point& vp0  = vcir->vertex()->point();
		Point& vp_1 = (--vcir)->opposite()->vertex()->point();
		vpt = Point((vp_1[0] + 6*vp0[0] + vp1[0])/8,
			(vp_1[1] + 6*vp0[1] + vp1[1])/8,
			(vp_1[2] + 6*vp0[2] + vp1[2])/8 );
	}
	//@}
};


// ======================================================================
///Loop subdivision rules
template <class _Poly>
class Loop_rule : public quadralize_rule<_Poly> {
public:
	typedef _Poly                                        Polyhedron;

	typedef typename Polyhedron::Vertex_handle           Vertex_handle;
	typedef typename Polyhedron::Halfedge_handle         Halfedge_handle;
	typedef typename Polyhedron::Facet_handle            Facet_handle;

	typedef typename Polyhedron::Halfedge_around_facet_circulator  
		Halfedge_around_facet_circulator;
	typedef typename Polyhedron::Halfedge_around_vertex_circulator 
		Halfedge_around_vertex_circulator;

	typedef typename Polyhedron::Traits                  Traits;
	typedef typename Traits::Kernel                      Kernel;

	typedef typename Polyhedron::Point_3                 Point;

	/**@name Class Methods */
	//@{
public:
	void edge_point_rule(Halfedge_handle edge, Point& pt) {
		Point& p1 = edge->vertex()->point();
		Point& p2 = edge->opposite()->vertex()->point();
		Point& f1 = edge->next()->vertex()->point();
		Point& f2 = edge->opposite()->next()->vertex()->point();

		pt = Point((3*(p1[0]+p2[0])+f1[0]+f2[0])/8,
			(3*(p1[1]+p2[1])+f1[1]+f2[1])/8,
			(3*(p1[2]+p2[2])+f1[2]+f2[2])/8 );
	}

	void vertex_point_rule(Vertex_handle vertex, Point& pt) {
		Halfedge_around_vertex_circulator vcir = vertex->vertex_begin();
		int n = circulator_size(vcir);    

		float Q[] = {0.0, 0.0, 0.0}, R[] = {0.0, 0.0, 0.0};
		Point& S = vertex->point();

		for (int i = 0; i < n; i++, ++vcir) {
			Point& p = vcir->opposite()->vertex()->point();
			R[0] += p[0]; 	R[1] += p[1]; 	R[2] += p[2];
		}
		if (n == 6) {
			pt = Point((10*S[0]+R[0])/16, (10*S[1]+R[1])/16, (10*S[2]+R[2])/16);
		} else {
			double Cn = 5.0/8.0 - std::sqrt(3+2*std::cos(6.283/n))/64.0;
			double Sw = n*(1-Cn)/Cn;
			double W = n/Cn;
			pt = Point((Sw*S[0]+R[0])/W, (Sw*S[1]+R[1])/W, (Sw*S[2]+R[2])/W);
		}
	}

	void border_point_rule(Halfedge_handle edge, Point& ept, Point& vpt) {
		Point& ep1 = edge->vertex()->point();
		Point& ep2 = edge->opposite()->vertex()->point();
		ept = Point((ep1[0]+ep2[0])/2, (ep1[1]+ep2[1])/2, (ep1[2]+ep2[2])/2);

		Halfedge_around_vertex_circulator vcir = edge->vertex_begin();
		Point& vp1  = vcir->opposite()->vertex()->point();
		Point& vp0  = vcir->vertex()->point();
		Point& vp_1 = (--vcir)->opposite()->vertex()->point();
		vpt = Point((vp_1[0] + 6*vp0[0] + vp1[0])/8,
			(vp_1[1] + 6*vp0[1] + vp1[1])/8,
			(vp_1[2] + 6*vp0[2] + vp1[2])/8 );
	}
	//@}
};

///
// ======================================================================
///Loop subdivision rules, the parameters will be subdivision together
template <class _Poly>
class Loop_Param_rule : public quadralize_rule<_Poly> {
public:
	typedef _Poly                                        Polyhedron;

	typedef typename Polyhedron::Vertex_handle           Vertex_handle;
	typedef typename Polyhedron::Halfedge_handle         Halfedge_handle;
	typedef typename Polyhedron::Facet_handle            Facet_handle;

	typedef typename Polyhedron::Halfedge_around_facet_circulator  
		Halfedge_around_facet_circulator;
	typedef typename Polyhedron::Halfedge_around_vertex_circulator 
		Halfedge_around_vertex_circulator;

	typedef typename Polyhedron::Traits                  Traits;
	typedef typename Traits::Kernel                      Kernel;

	typedef typename Polyhedron::Point_3                 Point;
	//typedef typename Polyhedron::Point_2                 Point_2;

	/**@name Class Methods */
	//@{
public:
	void edge_point_rule(Halfedge_handle edge, Point& pt, Point& uv) {
		Point& p1 = edge->vertex()->point();
		Point& p2 = edge->opposite()->vertex()->point();
		Point& f1 = edge->next()->vertex()->point();
		Point& f2 = edge->opposite()->next()->vertex()->point();


		//get uv 
		uv = Point((edge->vertex()->u()+edge->opposite()->vertex()->u())/2,
			(edge->vertex()->v()+edge->opposite()->vertex()->v())/2,0);


		pt = Point((3*(p1[0]+p2[0])+f1[0]+f2[0])/8,
			(3*(p1[1]+p2[1])+f1[1]+f2[1])/8,
			(3*(p1[2]+p2[2])+f1[2]+f2[2])/8 );
	}

	void vertex_point_rule(Vertex_handle vertex, Point& pt, Point& uv) 
	{
		Halfedge_around_vertex_circulator vcir = vertex->vertex_begin();
		int n = circulator_size(vcir);    

		//parameterization
		//get u, v
		uv = Point(vertex->u(), vertex->v(),0);

		float Q[] = {0.0, 0.0, 0.0}, R[] = {0.0, 0.0, 0.0};
		Point& S = vertex->point();

		for (int i = 0; i < n; i++, ++vcir) {
			Point& p = vcir->opposite()->vertex()->point();
			R[0] += p[0]; 	R[1] += p[1]; 	R[2] += p[2];
		}
		if (n == 6) {
			pt = Point((10*S[0]+R[0])/16, (10*S[1]+R[1])/16, (10*S[2]+R[2])/16);
		} else {
			double Cn = 5.0/8.0 - std::sqrt(3+2*std::cos(6.283/n))/64.0;
			double Sw = n*(1-Cn)/Cn;
			double W = n/Cn;
			pt = Point((Sw*S[0]+R[0])/W, (Sw*S[1]+R[1])/W, (Sw*S[2]+R[2])/W);
		}
	}

	void border_point_rule(Halfedge_handle edge, Point& ept, Point& vpt, Point& vuv, Point& euv)
	{

		Point& ep1 = edge->vertex()->point();
		Point& ep2 = edge->opposite()->vertex()->point();
		ept = Point((ep1[0]+ep2[0])/2, (ep1[1]+ep2[1])/2, (ep1[2]+ep2[2])/2);

		//parameterization
		//get u, v		
		euv = Point((edge->vertex()->u()+edge->opposite()->vertex()->u())/2,
			(edge->vertex()->v()+edge->opposite()->vertex()->v())/2,0);
		vuv = Point(edge->vertex()->u(), edge->vertex()->v(),0);

		Halfedge_around_vertex_circulator vcir = edge->vertex_begin();
		Point& vp1  = vcir->opposite()->vertex()->point();
		Point& vp0  = vcir->vertex()->point();
		Point& vp_1 = (--vcir)->opposite()->vertex()->point();
		vpt = Point((vp_1[0] + 6*vp0[0] + vp1[0])/8,
			(vp_1[1] + 6*vp0[1] + vp1[1])/8,
			(vp_1[2] + 6*vp0[2] + vp1[2])/8 );
	}
	//@}
};

///
template <class _Poly>
class Modify_Butterfly_rule : public quadralize_rule<_Poly> {
public:
	typedef _Poly                                        Polyhedron;

	typedef typename Polyhedron::Vertex_handle           Vertex_handle;
	typedef typename Polyhedron::Halfedge_handle         Halfedge_handle;
	typedef typename Polyhedron::Facet_handle            Facet_handle;

	typedef typename Polyhedron::Halfedge_around_facet_circulator  
		Halfedge_around_facet_circulator;
	typedef typename Polyhedron::Halfedge_around_vertex_circulator 
		Halfedge_around_vertex_circulator;

	typedef typename Polyhedron::Traits                  Traits;
	typedef typename Traits::Kernel                      Kernel;

	typedef typename Polyhedron::Point_3                 Point;

	/**@name Class Methods */
	//@{
public:
	void edge_point_rule(Halfedge_handle edge, Point& pt) {
		Point& p1 = edge->vertex()->point();
		Point& p2 = edge->opposite()->vertex()->point();

		int result = check_edge_type(edge);

		//border case and irrgular case
		//	TRACE("\n%d",result);
		if(result == 1 || result == 2)
		{
			pt = Point((p1[0]+p2[0])/2,(p1[1]+p2[1])/2,(p1[2]+p2[2])/2);
			return;
		}

		//regular case
		if(result == 3)
		{
			Point& f1 = edge->next()->vertex()->point();  
			Point& f2 = edge->opposite()->next()->vertex()->point();

			Point& vp1 = edge->next()->opposite()->next()->vertex()->point();
			Point& vp2 = edge->prev()->opposite()->next()->vertex()->point();
			Point& vp3 = edge->opposite()->next()->opposite()->next()->vertex()->point();
			Point& vp4 = edge->opposite()->prev()->opposite()->next()->vertex()->point();

			pt = Point(((p1[0]+p2[0])/2 + (f1[0]+f2[0])/8 - (vp1[0]+vp2[0]+vp3[0]+vp4[0])/16),
				((p1[1]+p2[1])/2 + (f1[1]+f2[1])/8 - (vp1[1]+vp2[1]+vp3[1]+vp4[1])/16),
				((p1[2]+p2[2])/2 + (f1[2]+f2[2])/8 - (vp1[2]+vp2[2]+vp3[2]+vp4[2])/16)); 

			return;		
		}

		Vertex_handle V1 = edge->vertex();
		Vertex_handle V2 = edge->opposite()->vertex(); 

		//irrgular case 4
		if(result == 4)
		{
			float coord[3];
			//ZeroMemory(coord, sizeof(float)*3);

			for(int j=0;j<3;j++)
				coord[j] = ButterflyIrregu(edge, j);

			pt = Point(coord[0],coord[1],coord[2]);
			return;
		}
		//irrgular case 5

		if(result == 5 )
		{
			float coord[3];
			//ZeroMemory(coord, sizeof(float)*3);

			for(int j=0;j<3;j++)
				coord[j] = ButterflyIrregu( edge->opposite(), j);

			pt = Point(coord[0],coord[1],coord[2]);
			return;

		}
		//irrgular case 6	
		if(result == 6)
		{
			float coord1[3],coord2[3];
			//ZeroMemory(coord1, sizeof(float)*3);
			//ZeroMemory(coord2, sizeof(float)*3);

			for(int j=0;j<3;j++){
				coord1[j] = ButterflyIrregu(edge, j);
				coord2[j] = ButterflyIrregu(edge->opposite(), j);
			}

			pt = Point((coord1[0]+coord2[0])/2,(coord1[1]+coord2[1])/2, (coord1[2]+coord2[2])/2);
			return;

		}


	}

	//irregular case4
	float ButterflyIrregu(Halfedge_handle edge, int i)
	{
		float coord = 0.0;

		Vertex_handle V = edge->vertex();
		Point& p1 = V->point();

		std::size_t nb = V->vertex_degree();

		Halfedge_around_vertex_circulator vcir = edge->vertex_begin();

		for(int j=0; j<nb; j++,vcir++){
			Point& temp = vcir->opposite()->vertex()->point();
			coord += temp[i]*ButterflyParam(nb,j);
		}

		coord += p1[i]*0.75;
		return coord;
	}

	//get parameter of modify butterfly subdivision scheme
	// n: degree of the current vertex
	// i: index
	float ButterflyParam(int n, int i)
	{
		float para;
		if( n == 3)
		{
			if(i == 0)
				para = 5.0f/12.0f;
			else
				para = -1.0f/12.0f;
		}
		else if( n == 4 )
		{
			if(i == 0)
				para = 3.0f/8.0f;
			else if(i == 2)
				para = -1.0f/8.0f;
			else
				para = 0.0f;
		}
		else
		{
			para = 0.25f + (float)cos(2.0f * 3.14159265359f*(float)i/(float)n) +
				0.5f * (float)cos(4.0f * 3.14159265359f*(float)i/(float)n);
			para = para/(float)n;

		}


		return para;
	}

	// check the type of edge
	int check_edge_type(Halfedge_handle edge)
	{

		//case 2 other edge of the border face 
		if(edge->next()->is_border_edge() || edge->prev()->is_border_edge())
			return 1;

		//vertices of edge
		Vertex_handle V1 = edge->vertex();
		Vertex_handle V2 = edge->opposite()->vertex();

		std::size_t nbv1 = V1->vertex_degree();
		std::size_t nbv2 = V2->vertex_degree();

		//case 3, mask is not complete
		int result = 0;
		Halfedge_around_vertex_circulator vcir = V1->vertex_begin();
		for(int i=0; i< nbv1; i++,vcir++)
			result += int(vcir->is_border_edge());
		vcir = V2->vertex_begin();
		for(int i=0; i< nbv2; i++,vcir++)
			result += int(vcir->is_border_edge());
		if(result > 0)
			return 2;

		//case 4; regulare case
		if(nbv1==6 && nbv2 ==6)
			return 3;

		//case5. current vertex is not irregular, but next is.
		else if(nbv1!=6 && nbv2 == 6)
			return 4;

		//case6  next vertex is not irrgular, but this is
		else if(nbv1==6 && nbv2 != 6)
			return 5;

		//default
		else
			return 6;

	}


	void vertex_point_rule(Vertex_handle vertex, Point& pt) {

		Halfedge_around_vertex_circulator vcir = vertex->vertex_begin();
		int n = circulator_size(vcir);    

		float Q[] = {0.0, 0.0, 0.0}, R[] = {0.0, 0.0, 0.0};
		Point& S = vertex->point();


		pt = Point(S[0], S[1], S[2]);
		//  }

	}

	void border_point_rule(Halfedge_handle edge, Point& ept, Point& vpt) {
		Point& ep1 = edge->vertex()->point();
		Point& ep2 = edge->opposite()->vertex()->point();
		vpt = Point(ep1[0],ep1[1],ep1[2]);

		ept = Point((ep1[0]+ep2[0])/2,(ep1[1]+ep2[1])/2,(ep1[2]+ep2[2])/2);


		/*	Halfedge_around_vertex_circulator vcir = edge->vertex_begin();
		Point& vp0  = (--vcir)->opposite()->vertex()->point();

		Halfedge_around_vertex_circulator vcir1 = edge->opposite()->vertex_begin();
		Point& vp1  = (--vcir1)->opposite()->vertex()->point();


		ept = Point(((ep1[0]+ep2[0])*9/16-(vp0[0]+vp1[0])/16), ((ep1[1]+ep2[1])*9/16-(vp0[1]+vp1[1])/16), ((ep1[2]+ep2[2])*9/16-(vp0[2]+vp1[2])/16));
		*/
	}
	//@}
};




/***********************************************************************************
*Modify_Butterfly subdivision rules, the parameters will be subdivision together
*
***********************************************************************************/
template <class _Poly>
class Modify_Butterfly_Param_rule : public quadralize_rule<_Poly> {
public:
	typedef _Poly                                        Polyhedron;

	typedef typename Polyhedron::Vertex_handle           Vertex_handle;
	typedef typename Polyhedron::Halfedge_handle         Halfedge_handle;
	typedef typename Polyhedron::Facet_handle            Facet_handle;

	typedef typename Polyhedron::Halfedge_around_facet_circulator  
		Halfedge_around_facet_circulator;
	typedef typename Polyhedron::Halfedge_around_vertex_circulator 
		Halfedge_around_vertex_circulator;

	typedef typename Polyhedron::Traits                  Traits;
	typedef typename Traits::Kernel                      Kernel;

	typedef typename Polyhedron::Point_3                 Point;
	//typedef typename Polyhedron::Point_3                 Point_3;

	/**@name Class Methods */
	//@{
public:
	void edge_point_rule(Halfedge_handle edge, Point& pt, Point& uv) {
		Point& p1 = edge->vertex()->point();
		Point& p2 = edge->opposite()->vertex()->point();

		uv = Point((edge->vertex()->u()+edge->opposite()->vertex()->u())/2,
			(edge->vertex()->v()+edge->opposite()->vertex()->v())/2, 0);

		//others case
		int result = check_edge_type(edge);

		// irrgular case
		if(result == 1 || result == 2)
		{
			pt = Point((p1[0]+p2[0])/2,(p1[1]+p2[1])/2,(p1[2]+p2[2])/2);
			return;
		}

		//regular case
		if(result == 3)
		{
			Point& f1 = edge->next()->vertex()->point();  
			Point& f2 = edge->opposite()->next()->vertex()->point();

			Point& vp1 = edge->next()->opposite()->next()->vertex()->point();
			Point& vp2 = edge->prev()->opposite()->next()->vertex()->point();
			Point& vp3 = edge->opposite()->next()->opposite()->next()->vertex()->point();
			Point& vp4 = edge->opposite()->prev()->opposite()->next()->vertex()->point();

			pt = Point(((p1[0]+p2[0])/2 + (f1[0]+f2[0])/8 - (vp1[0]+vp2[0]+vp3[0]+vp4[0])/16),
				((p1[1]+p2[1])/2 + (f1[1]+f2[1])/8 - (vp1[1]+vp2[1]+vp3[1]+vp4[1])/16),
				((p1[2]+p2[2])/2 + (f1[2]+f2[2])/8 - (vp1[2]+vp2[2]+vp3[2]+vp4[2])/16)); 

			return;		
		}

		Vertex_handle V1 = edge->vertex();
		Vertex_handle V2 = edge->opposite()->vertex(); 

		//irrgular case 4
		if(result == 4)
		{
			float coord[3];
			//ZeroMemory(coord, sizeof(float)*3);

			for(int j=0;j<3;j++)
				coord[j] = ButterflyIrregu(edge, j);

			pt = Point(coord[0],coord[1],coord[2]);
			return;
		}
		//irrgular case 5

		if(result == 5 )
		{
			float coord[3];
			//ZeroMemory(coord, sizeof(float)*3);

			for(int j=0;j<3;j++)
				coord[j] = ButterflyIrregu( edge->opposite(), j);

			pt = Point(coord[0],coord[1],coord[2]);
			return;

		}
		//irrgular case 6	
		if(result == 6)
		{
			float coord1[3],coord2[3];
			//ZeroMemory(coord1, sizeof(float)*3);
			//ZeroMemory(coord2, sizeof(float)*3);

			for(int j=0;j<3;j++){
				coord1[j] = ButterflyIrregu(edge, j);
				coord2[j] = ButterflyIrregu(edge->opposite(), j);
			}

			pt = Point((coord1[0]+coord2[0])/2,(coord1[1]+coord2[1])/2, (coord1[2]+coord2[2])/2);
			return;

		}


	}

	//irregular case4
	float ButterflyIrregu(Halfedge_handle edge, int i)
	{
	
		static float Param[28][30] = {{0.416667,-0.083333,-0.083333,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},
		{0.375000,0.000000,-0.125000,0.000000,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},
		{0.350000,0.030902,-0.080902,-0.080902,0.030902,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},
		{0.291667,0.083333,-0.083333,-0.041667,-0.083333,0.083333,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},
		{0.250000,0.108890,-0.060429,-0.048461,-0.048461,-0.060429,0.108890,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},
		{0.218750,0.119638,-0.031250,-0.057138,-0.031250,-0.057138,-0.031250,0.119638,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},
		{0.194444,0.122541,-0.005133,-0.055556,-0.034074,-0.034074,-0.055556,-0.005133,0.122541,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},
		{0.175000,0.121353,0.015451,-0.046353,-0.040451,-0.025000,-0.040451,-0.046353,0.015451,0.121353,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},
		{0.159091,0.118087,0.030726,-0.033824,-0.043274,-0.026261,-0.026261,-0.043274,-0.033824,0.030726,0.118087,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},
		{0.145833,0.113835,0.041667,-0.020833,-0.041667,-0.030502,-0.020833,-0.030502,-0.041667,-0.020833,0.041667,0.113835,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},
		{0.134615,0.109191,0.049289,-0.008841,-0.036835,-0.033711,-0.021401,-0.021401,-0.033711,-0.036835,-0.008841,0.049289,0.109191,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},
		{0.125000,0.104480,0.054445,0.001574,-0.030215,-0.034625,-0.024230,-0.017857,-0.024230,-0.034625,-0.030215,0.001574,0.054445,0.104480,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},
		{0.116667,0.099874,0.057791,0.010301,-0.022907,-0.033333,-0.026967,-0.018092,-0.018092,-0.026967,-0.033333,-0.022907,0.010301,0.057791,0.099874,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},
		{0.109375,0.095465,0.059819,0.017446,-0.015625,-0.030390,-0.028569,-0.020020,-0.015625,-0.020020,-0.028569,-0.030390,-0.015625,0.017446,0.059819,0.095465,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},
		{0.102941,0.091293,0.060891,0.023201,-0.008778,-0.026398,-0.028792,-0.022197,-0.015690,-0.015690,-0.022197,-0.028792,-0.026398,-0.008778,0.023201,0.060891,0.091293,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},
		{0.097222,0.087373,0.061270,0.027778,-0.002567,-0.021861,-0.027778,-0.023846,-0.017037,-0.013889,-0.017037,-0.023846,-0.027778,-0.021861,-0.002567,0.027778,0.061270,0.087373,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},
		{0.092105,0.083705,0.061152,0.031374,0.002934,-0.017145,-0.025807,-0.024662,-0.018737,-0.013866,-0.013866,-0.018737,-0.024662,-0.025807,-0.017145,0.002934,0.031374,0.061152,0.083705,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},
		{0.087500,0.080278,0.060676,0.034164,0.007725,-0.012500,-0.023176,-0.024615,-0.020225,-0.014827,-0.012500,-0.014827,-0.020225,-0.024615,-0.023176,-0.012500,0.007725,0.034164,0.060676,0.080278,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},
		{0.083333,0.077081,0.059948,0.036297,0.011848,-0.008080,-0.020143,-0.023810,-0.021223,-0.016154,-0.012431,-0.012431,-0.016154,-0.021223,-0.023810,-0.020143,-0.008080,0.011848,0.036297,0.059948,0.077081,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},
		{0.079545,0.074096,0.059044,0.037896,0.015363,-0.003974,-0.016912,-0.022402,-0.021637,-0.017434,-0.013130,-0.011364,-0.013130,-0.017434,-0.021637,-0.022402,-0.016912,-0.003974,0.015363,0.037896,0.059044,0.074096,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},
		{0.076087,0.071310,0.058020,0.039062,0.018336,-0.000224,-0.013634,-0.020554,-0.021483,-0.018434,-0.014171,-0.011271,-0.011271,-0.014171,-0.018434,-0.021483,-0.020554,-0.013634,-0.000224,0.018336,0.039062,0.058020,0.071310,0.0,0.0,0.0,0.0,0.0,0.0,0.0},
		{0.072917,0.068706,0.056918,0.039879,0.020833,0.003159,-0.010417,-0.018410,-0.020833,-0.019046,-0.015251,-0.011788,-0.010417,-0.011788,-0.015251,-0.019046,-0.020833,-0.018410,-0.010417,0.003159,0.020833,0.039879,0.056918,0.068706,0.0,0.0,0.0,0.0,0.0,0.0},
		{0.070000,0.066269,0.055769,0.040415,0.022917,0.006180,-0.007331,-0.016091,-0.019780,-0.019245,-0.016180,-0.012612,-0.010313,-0.010313,-0.012612,-0.016180,-0.019245,-0.019780,-0.016091,-0.007331,0.006180,0.022917,0.040415,0.055769,0.066269,0.0,0.0,0.0,0.0,0.0},
		{0.067308,0.063987,0.054596,0.040722,0.024645,0.008860,-0.004421,-0.013693,-0.018418,-0.019053,-0.016855,-0.013516,-0.010701,-0.009615,-0.010701,-0.013516,-0.016855,-0.019053,-0.018418,-0.013693,-0.004421,0.008860,0.024645,0.040722,0.054596,0.063987,0.0,0.0,0.0,0.0},
		{0.064815,0.061847,0.053415,0.040847,0.026065,0.011221,-0.001711,-0.011288,-0.016835,-0.018519,-0.017234,-0.014350,-0.011358,-0.009508,-0.009508,-0.011358,-0.014350,-0.017234,-0.018519,-0.016835,-0.011288,-0.001711,0.011221,0.026065,0.040847,0.053415,0.061847,0.0,0.0,0.0},
		{0.062500,0.059836,0.052240,0.040825,0.027222,0.013291,0.000787,-0.008929,-0.015107,-0.017701,-0.017313,-0.015020,-0.012115,-0.009802,-0.008929,-0.009802,-0.012115,-0.015020,-0.017313,-0.017701,-0.015107,-0.008929,0.000787,0.013291,0.027222,0.040825,0.052240,0.059836,0.0,0.0},
		{0.060345,0.057945,0.051078,0.040685,0.028155,0.015097,0.003072,-0.006653,-0.013297,-0.016660,-0.017112,-0.015480,-0.012850,-0.010331,-0.008822,-0.008822,-0.010331,-0.012850,-0.015480,-0.017112,-0.016660,-0.013297,-0.006653,0.003072,0.015097,0.028155,0.040685,0.051078,0.057945,0.0},
		{0.058333,0.056164,0.049937,0.040451,0.028896,0.016667,0.005150,-0.004485,-0.011453,-0.015451,-0.016667,-0.015713,-0.013484,-0.010966,-0.009046,-0.008333,-0.009046,-0.010966,-0.013484,-0.015713,-0.016667,-0.015451,-0.011453,-0.004485,0.005150,0.016667,0.028896,0.040451,0.049937,0.056164}};

		float coord = 0.0;

		Vertex_handle V = edge->vertex();
		Point& p1 = V->point();

		std::size_t nb = V->vertex_degree();

		Halfedge_around_vertex_circulator vcir = edge->vertex_begin();

		for(int j=0; j<nb; j++,vcir++){
			Point& temp = vcir->opposite()->vertex()->point();
			coord += temp[i]*Param[nb-3][j];
		}

		coord += p1[i]*0.75;
		return coord;
	}

	//get parameter of modify butterfly subdivision scheme
	// n: degree of the current vertex
	// i: index
	float ButterflyParam(int n, int i)
	{
		float para;
		if( n == 3)
		{
			if(i == 0)
				para = 0.41666666666667f;
			else
				para = -0.083333333333333f;
		}
		else if( n == 4 )
		{
			if(i == 0)
				para = 0.375f;
			else if(i == 2)
				para = -0.125f;
			else
				para = 0.0f;
		}
		else
		{
			para = 0.25f + (float)cos(2.0f * 3.14159265359f*(float)i/(float)n) +
				0.5f * (float)cos(4.0f * 3.14159265359f*(float)i/(float)n);
			para = para/(float)n;

		}


		return para;
	}

	// check the type of edge
	int check_edge_type(Halfedge_handle edge)
	{

		//case 2 other edge of the border face 
		if(edge->next()->is_border_edge() || edge->prev()->is_border_edge())
			return 1;

		//vertices of edge
		Vertex_handle V1 = edge->vertex();
		Vertex_handle V2 = edge->opposite()->vertex();

		std::size_t nbv1 = V1->vertex_degree();
		std::size_t nbv2 = V2->vertex_degree();

		//case 3, mask is not complete
		int result = 0;
		Halfedge_around_vertex_circulator vcir = V1->vertex_begin();
		for(int i=0; i< nbv1; i++,vcir++)
			result += int(vcir->is_border_edge());
		vcir = V2->vertex_begin();
		for(int i=0; i< nbv2; i++,vcir++)
			result += int(vcir->is_border_edge());
		if(result > 0)
			return 2;

		//case 4; regulare case
		if(nbv1==6 && nbv2 ==6)
			return 3;

		//case5. current vertex is not irregular, but next is.
		else if(nbv1!=6 && nbv2 == 6)
			return 4;

		//case6  next vertex is not irrgular, but this is
		else if(nbv1==6 && nbv2 != 6)
			return 5;

		//default
		else
			return 6;

	}


	void vertex_point_rule(Vertex_handle vertex, Point& pt, Point& uv) {

		Halfedge_around_vertex_circulator vcir = vertex->vertex_begin();
		int n = circulator_size(vcir);    

		float Q[] = {0.0, 0.0, 0.0}, R[] = {0.0, 0.0, 0.0};
		Point& S = vertex->point();


		//coord
		pt = Point(S[0], S[1], S[2]);

		//parameterization
		//get u, v
		uv = Point(vertex->u(), vertex->v(),0);
		//  }

	}

	void border_point_rule(Halfedge_handle edge, Point& ept, Point& vpt, Point& vuv, Point& euv) {
		Point& ep1 = edge->vertex()->point();
		Point& ep2 = edge->opposite()->vertex()->point();
		vpt = Point(ep1[0],ep1[1],ep1[2]);

		ept = Point((ep1[0]+ep2[0])/2,(ep1[1]+ep2[1])/2,(ep1[2]+ep2[2])/2);


		//get u,v

		euv = Point((edge->vertex()->u()+edge->opposite()->vertex()->u())/2,
			(edge->vertex()->v()+edge->opposite()->vertex()->v())/2, 0);
		vuv = Point(edge->vertex()->u(),edge->vertex()->v(),0);
		/*	Halfedge_around_vertex_circulator vcir = edge->vertex_begin();
		Point& vp0  = (--vcir)->opposite()->vertex()->point();

		Halfedge_around_vertex_circulator vcir1 = edge->opposite()->vertex_begin();
		Point& vp1  = (--vcir1)->opposite()->vertex()->point();


		ept = Point(((ep1[0]+ep2[0])*9/16-(vp0[0]+vp1[0])/16), ((ep1[1]+ep2[1])*9/16-(vp0[1]+vp1[1])/16), ((ep1[2]+ep2[2])*9/16-(vp0[2]+vp1[2])/16));
		*/
	}
	//@}
};




//==========================================================================
template <class _Poly>
class dualize_rule {
public:
	typedef _Poly                                        Polyhedron;
	typedef typename Polyhedron::Halfedge_handle         Halfedge_handle;

	typedef typename Polyhedron::Traits                  Traits;
	typedef typename Traits::Kernel                      Kernel;
	typedef typename Polyhedron::Point_3                 Point;

public:
	void point_rule(Halfedge_handle edge, Point& pt) {};
};


// ======================================================================
///
template <class _Poly>
class DooSabin_rule : public dualize_rule<_Poly> {
public:
	typedef _Poly                                        Polyhedron;
	typedef typename Polyhedron::Halfedge_handle         Halfedge_handle;
	typedef typename Polyhedron::Halfedge_around_facet_circulator  
		Halfedge_around_facet_circulator;

	typedef typename Polyhedron::Traits                  Traits;
	typedef typename Traits::Kernel                      Kernel;
	typedef typename Polyhedron::Point_3                 Point;
	typedef CGAL::Vector_3<Kernel>                       Vector;

public:
	void point_rule(Halfedge_handle he, Point& pt) {
		int n =  CGAL::circulator_size(he->facet()->facet_begin()); 

		Vector cv(0,0,0);
		if (n == 4) {
			cv = cv + (he->vertex()->point()-CGAL::ORIGIN)*9;
			cv = cv + (he->next()->vertex()->point()-CGAL::ORIGIN)*3;
			cv = cv + (he->next()->next()->vertex()->point()-CGAL::ORIGIN);
			cv = cv + (he->prev()->vertex()->point()-CGAL::ORIGIN)*3;
			cv = cv/16;
		} else {
			double a;
			for (int k = 0; k < n; ++k, he = he->next()) {
				if (k == 0) a = ((double)5/n) + 1;
				else a = (3+2*std::cos(2*k*3.141593/n))/n;
				cv = cv + (he->vertex()->point()-CGAL::ORIGIN)*a;
			}
			cv = cv/4;
		}
		pt = CGAL::ORIGIN + cv;
	}
};

DGAL_END_NAMESPACE

#endif //_POLYHEDRON_SUBDIVISION_RULES_H_01292002
