// Copyright (c) 2007-2010  Dalian University of Technology (China). 
// All rights reserved.
//
// This file is part of DGAL; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL: http://jjcao1231.googlepages.com $
// $Id: Gl_polyhedron_3_policy.h 2007-04-06 $
//
//
// Author(s)     : JJCAO <jjcao@hotmail.com>

#ifndef DGAL_GL_POLYHEDRON_3_POLICY_H
#define DGAL_GL_POLYHEDRON_3_POLICY_H

#pragma warning (push)
	#pragma warning(disable : 4267 4503 4996)
	#include <CGAL/Cartesian.h>
	#include <CGAL/Polyhedron_3.h>
	#include <DGAL/config.h>
	#include <Fl/Gl.h>
	#include <Gl/Glu.h>
	#include <algorithm>
	#include <list>
#pragma warning (pop)

DGAL_BEGIN_NAMESPACE

template 
<
	class PolyhedronTraits_3, 
	class PolyhedronItems_3,
#ifndef CGAL_CFG_NO_TMPL_IN_TMPL_PARAM
       template < class T, class I, class A>
#endif
	class T_HDS,
	class Alloc
>
class Gl_polyhedron_3_policy
{
public:
	typedef typename PolyhedronTraits_3    Traits;
	typedef typename PolyhedronItems_3     Items;
	typedef typename Gl_polyhedron_3_policy<Traits, Items, T_HDS, Alloc>
		                                   Self;

	typedef CGAL::I_Polyhedron_derived_items_3<Items>   Derived_items;
#ifndef CGAL_CFG_NO_TMPL_IN_TMPL_PARAM
    typedef T_HDS< Traits, Derived_items, Alloc>  HDS;
#else
    typedef typename T_HDS::template HDS< Traits, Derived_items, Alloc>  HDS;
#endif

    // Container stuff.
    typedef typename HDS::size_type               size_type;
    typedef typename HDS::difference_type         difference_type;
    typedef typename HDS::iterator_category       iterator_category;
    typedef typename HDS::Supports_removal        Supports_removal;

    // Geometry
	typedef typename Traits::FT                   FT;
    typedef typename Traits::Point_3              Point_3;
	typedef typename Traits::Point_2              Point_2;
	typedef typename Traits::Vector_3             Vector_3;
    typedef typename Traits::Plane_3              Plane_3;
	typedef typename Traits::Iso_cuboid_3         Iso_cuboid;

    // Items
    typedef typename HDS::Vertex                  Vertex;
    typedef typename HDS::Halfedge                Halfedge;
    typedef typename HDS::Face                    Face;

    typedef typename Vertex::Base                 VBase;
    typedef typename Halfedge::Base               HBase;
    typedef typename Face::Base                   FBase;

	typedef typename HBase::Supports_halfedge_prev  Supports_prev;

    // Handles and Iterators
    typedef typename HDS::Vertex_handle           Vertex_handle;
    typedef typename HDS::Halfedge_handle         Halfedge_handle;
    typedef typename HDS::Face_handle             Face_handle;
    typedef typename HDS::Vertex_iterator         Vertex_iterator;
    typedef typename HDS::Halfedge_iterator       Halfedge_iterator;
    typedef typename HDS::Face_iterator           Face_iterator;

    typedef typename HDS::Vertex_const_handle     Vertex_const_handle;
    typedef typename HDS::Halfedge_const_handle   Halfedge_const_handle;
    typedef typename HDS::Face_const_handle       Face_const_handle;
    typedef typename HDS::Vertex_const_iterator   Vertex_const_iterator;
    typedef typename HDS::Halfedge_const_iterator Halfedge_const_iterator;
    typedef typename HDS::Face_const_iterator     Face_const_iterator;

	/////////////
	// Circulator category.
	typedef CGAL::HalfedgeDS_circulator_traits<Supports_prev> 
		                                          Circulator_traits;
    typedef typename Circulator_traits::iterator_category 
		                                          Circulator_category;

    // Circulators around a vertex and around a facet.
    typedef CGAL::I_HalfedgeDS_facet_circ< Halfedge_handle, Circulator_category>
                                                  Halfedge_around_face_cior;

    typedef CGAL::I_HalfedgeDS_vertex_circ< Halfedge_handle, Circulator_category>
                                                  Halfedge_around_vertex_cior;

    typedef CGAL::I_HalfedgeDS_facet_circ<Halfedge_const_handle, Circulator_category>
		                                          Halfedge_around_face_const_cior;

    typedef CGAL::I_HalfedgeDS_vertex_circ<Halfedge_const_handle, Circulator_category>
		                                          Halfedge_around_vertex_const_cior;


	/////////////
	typedef typename std::list<Face_handle>::iterator         Face_itor;
	typedef typename std::list<Vertex_handle>::iterator       Vertex_itor;
	typedef typename std::list<Halfedge_handle>::iterator     Halfedge_itor;
	typedef CGAL::N_step_adaptor_derived<Halfedge_iterator, 2>
                                                  Edge_itor;
    typedef CGAL::N_step_adaptor_derived<Halfedge_const_iterator, 2>
                                                  Edge_const_itor;
public:
	
	void draw_facets(Face_itor begin, Face_itor end,
		             bool smooth_shading,bool use_normals, 
					 bool use_texture=false)
	{
		for(Face_itor pFacet = begin;	pFacet != end; ++pFacet)
		{
			Face_handle fh = *pFacet;
			// begin polygon assembly
			::glBegin(GL_POLYGON);
			draw_facet(fh,smooth_shading,use_normals, use_texture);
			::glEnd(); // end polygon assembly
		}
		glFlush();
	}
	void draw_facets(Face_const_iterator begin, Face_const_iterator end,
		             bool smooth_shading,bool use_normals, 
					 bool use_texture=false, int pick_mode=0)
	{
		int i(0);
		for(Face_const_iterator pFacet = begin; pFacet != end; ++pFacet, ++i)
		{
			if(pick_mode==3)
				glLoadName(i);

			// begin polygon assembly
			::glBegin(GL_POLYGON);
			draw_facet(pFacet,smooth_shading,use_normals, use_texture);
			::glEnd(); // end polygon assembly
		}
		glFlush();
	}
	virtual void draw_facet(Face_const_handle pFace, bool smooth_shading,
		            bool use_normals, bool use_texture=false)
	{
		// one normal per face
		if(use_normals && !smooth_shading)
		{
			const Face::Normal_3& normal = pFace->normal();
			::glNormal3f( (float)normal[0], (float)normal[1], (float)normal[2]);
		}

		// revolve around current face to get vertices
		Halfedge_around_face_const_cior pHalfedge = pFace->facet_begin();
		do
		{
			Vertex_const_handle vh = pHalfedge->vertex();

			// one normal per vertex
			if(use_normals && smooth_shading)
			{
				const Face::Normal_3& normal = vh->normal();
				::glNormal3f( (float)normal[0], (float)normal[1], (float)normal[2]);
			}

			// polygon assembly is performed per vertex
			const Point_3& point  = vh->point();		  
			if (use_texture)
			{				
				::glTexCoord2f((float)vh->u(), (float)vh->v());	
				::glVertex3d(point[0],point[1],point[2]);
			}
			else ::glVertex3d(point[0],point[1],point[2]);
		}
		while(++pHalfedge != pFace->facet_begin());
	}

	template<class Iterator>
	void draw_edges(Iterator begin, Iterator end)
	{
		::glBegin(GL_LINES);
		for(Iterator h = begin; h != end; ++h)
		{
			//draw line segment
			Halfedge_handle hh = *h;
			const Point_3& p1 = hh->prev()->vertex()->point();
			const Point_3& p2 = hh->vertex()->point();
			::glVertex3f( (float) p1[0],(float) p1[1],(float) p1[2]);
			::glVertex3f( (float) p2[0],(float) p2[1],(float) p2[2]);
		}		
		::glEnd();
	}
	void draw_edges(Edge_const_itor begin, Edge_const_itor end,
		            bool voronoi_edge = false, int pick_mode=0)
	{
		if(voronoi_edge)
		{
			::glBegin(GL_LINES);
			for(Edge_const_itor h = begin; h != end; ++h)
			{
				Halfedge_const_handle hh = h;
				Face_const_handle pFace1 = hh->facet();
				Face_const_handle pFace2 = hh->opposite()->facet();
				if(pFace1 == NULL || pFace2 == NULL)
				  continue;

				const Point_3 &p1 = hh->vertex()->point();
				const Point_3 &p2 = hh->next()->vertex()->point();
				const Point_3 &p3 = hh->next()->next()->vertex()->point();

				Kernel k;
				Point_3 d1 = k.construct_circumcenter_3_object()(p1,p2,p3);
				::glVertex3f(d1[0],d1[1],d1[2]);

				const Point_3 &pp1 = hh->opposite()->vertex()->point();
				const Point_3 &pp2 = hh->opposite()->next()->vertex()->point();
				const Point_3 &pp3 = hh->opposite()->next()->next()->vertex()->point();
				Point_3 d2 = k.construct_circumcenter_3_object()(pp1,pp2,pp3);
				::glVertex3f(d2[0],d2[1],d2[2]);
			}
			::glEnd();
			return;
		}

		int i(0);
		for(Edge_const_itor h = begin; h != end; ++h, ++i)
		{
			if(pick_mode==2)
				glLoadName(i);
			//draw line segment
			Halfedge_const_handle hh = h;
			const Point_3& p1 = hh->prev()->vertex()->point();
			const Point_3& p2 = hh->vertex()->point();
			::glBegin(GL_LINES);
			::glVertex3f( (float) p1[0],(float) p1[1],(float) p1[2]);
			::glVertex3f( (float) p2[0],(float) p2[1],(float) p2[2]);
			::glEnd();
		}
	}
	
	template<class Iterator>
	void draw_vertices(double scale,Iterator begin, Iterator end, int pick_mode)
	{
		GLUquadricObj* pQuadric = gluNewQuadric();
		::glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);

		int i(0); 
		for(Iterator it = begin;it !=  end; ++it, ++i)
		{
			if(pick_mode==1) glLoadName(i);

			Vertex_handle vh = *it;
			double radius = average_edge_length_around(vh);
			::glPushMatrix();
				::glTranslated(vh->point().x(),vh->point().y(), vh->point().z());
				::gluSphere(pQuadric,scale * radius,24,24); 
			::glPopMatrix();
		}

		gluDeleteQuadric(pQuadric);	
	}

	template<> 
	void draw_vertices <Vertex_const_iterator> (double scale,
		               Vertex_const_iterator begin, Vertex_const_iterator end,
					   int pick_mode)
	{
		GLUquadricObj* pQuadric = gluNewQuadric();
		::glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);

		int i(0); 
		for(Vertex_const_iterator it = begin;it !=  end; ++it, ++i)
		{
			if(pick_mode==1) glLoadName(i);

			Vertex_handle vh = it;
			double radius = average_edge_length_around(vh);
			::glPushMatrix();
				::glTranslated(vh->point().x(),vh->point().y(), vh->point().z());
				::gluSphere(pQuadric,scale * radius,24,24); 
			::glPopMatrix();
		}

		gluDeleteQuadric(pQuadric);
	}	

	template<> 
	void draw_vertices <Vertex_iterator> (double scale,
		               Vertex_iterator begin, Vertex_iterator end,
					   int pick_mode)
	{
		GLUquadricObj* pQuadric = gluNewQuadric();
		::glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);

		int i(0); 
		for(Vertex_iterator it = begin;it !=  end; ++it,++i)
		{
			if(pick_mode==1) glLoadName(i);

			Vertex_handle vh = it;
			double radius = average_edge_length_around(vh);
			::glPushMatrix();
				::glTranslated(vh->point().x(),vh->point().y(), vh->point().z());
				::gluSphere(pQuadric,scale * radius,24,24); 
			::glPopMatrix();
		}

		gluDeleteQuadric(pQuadric);
	}	

	// draw bounding box
	void draw_bounding_box(Iso_cuboid & bbox)
	{
	::glBegin(GL_LINES);

		// along x axis
		::glVertex3f(bbox.xmin(),bbox.ymin(),bbox.zmin());
		::glVertex3f(bbox.xmax(),bbox.ymin(),bbox.zmin());
		::glVertex3f(bbox.xmin(),bbox.ymin(),bbox.zmax());
		::glVertex3f(bbox.xmax(),bbox.ymin(),bbox.zmax());
		::glVertex3f(bbox.xmin(),bbox.ymax(),bbox.zmin());
		::glVertex3f(bbox.xmax(),bbox.ymax(),bbox.zmin());
		::glVertex3f(bbox.xmin(),bbox.ymax(),bbox.zmax());
		::glVertex3f(bbox.xmax(),bbox.ymax(),bbox.zmax());

		// along y axis
		::glVertex3f(bbox.xmin(),bbox.ymin(),bbox.zmin());
		::glVertex3f(bbox.xmin(),bbox.ymax(),bbox.zmin());
		::glVertex3f(bbox.xmin(),bbox.ymin(),bbox.zmax());
		::glVertex3f(bbox.xmin(),bbox.ymax(),bbox.zmax());
		::glVertex3f(bbox.xmax(),bbox.ymin(),bbox.zmin());
		::glVertex3f(bbox.xmax(),bbox.ymax(),bbox.zmin());
		::glVertex3f(bbox.xmax(),bbox.ymin(),bbox.zmax());
		::glVertex3f(bbox.xmax(),bbox.ymax(),bbox.zmax());

		// along z axis
		::glVertex3f(bbox.xmin(),bbox.ymin(),bbox.zmin());
		::glVertex3f(bbox.xmin(),bbox.ymin(),bbox.zmax());
		::glVertex3f(bbox.xmin(),bbox.ymax(),bbox.zmin());
		::glVertex3f(bbox.xmin(),bbox.ymax(),bbox.zmax());
		::glVertex3f(bbox.xmax(),bbox.ymin(),bbox.zmin());
		::glVertex3f(bbox.xmax(),bbox.ymin(),bbox.zmax());
		::glVertex3f(bbox.xmax(),bbox.ymax(),bbox.zmin());
		::glVertex3f(bbox.xmax(),bbox.ymax(),bbox.zmax());

	::glEnd();
	//draw an axis	
	float olw[] = {1.0};
	::glGetFloatv(GL_LINE_WIDTH, olw);
	glLineWidth(2.0);
	::glBegin(GL_LINES);
		glColor3f(1.0f,0.0f,0.0f);//red: x
		::glVertex3f(bbox.xmin(),bbox.ymin(),bbox.zmin());
		::glVertex3f(bbox.xmin()+1,bbox.ymin(),bbox.zmin());
		::glVertex3f(bbox.xmin()+0.9,bbox.ymin()+0.1,bbox.zmin());
		::glVertex3f(bbox.xmin()+1,bbox.ymin(),bbox.zmin());
		::glVertex3f(bbox.xmin()+0.9,bbox.ymin()-0.1,bbox.zmin());
		::glVertex3f(bbox.xmin()+1,bbox.ymin(),bbox.zmin());

		glColor3f(0.0f,1.0f,0.0f);//green: y
		::glVertex3f(bbox.xmin(),bbox.ymin(),bbox.zmin());
		::glVertex3f(bbox.xmin(),bbox.ymin()+1,bbox.zmin());
		::glVertex3f(bbox.xmin(),bbox.ymin()+0.9,bbox.zmin()+0.1);
		::glVertex3f(bbox.xmin(),bbox.ymin()+1,bbox.zmin());
		::glVertex3f(bbox.xmin(),bbox.ymin()+0.9,bbox.zmin()-0.1);
		::glVertex3f(bbox.xmin(),bbox.ymin()+1,bbox.zmin());

		glColor3f(0.0f,0.0f,1.0f);//blue: z
		::glVertex3f(bbox.xmin(),bbox.ymin(),bbox.zmin());
		::glVertex3f(bbox.xmin(),bbox.ymin(),bbox.zmin()+1);
		::glVertex3f(bbox.xmin()+0.1,bbox.ymin(),bbox.zmin()+0.9);
		::glVertex3f(bbox.xmin(),bbox.ymin(),bbox.zmin()+1);
		::glVertex3f(bbox.xmin()-0.1,bbox.ymin(),bbox.zmin()+0.9);
		::glVertex3f(bbox.xmin(),bbox.ymin(),bbox.zmin()+1);
		
	::glEnd();
	glLineWidth(olw[0]);
	}

	
private:
	// compute average edge length around a vertex
	FT average_edge_length_around(Vertex_handle pVertex)
	{
		FT sum = 0.0;
		Halfedge_around_vertex_cior pHalfEdge = pVertex->vertex_begin();
		Halfedge_around_vertex_cior end = pHalfEdge;
		Vector_3 vec(0.0,0.0,0.0);
		int degree = 0;
		CGAL_For_all(pHalfEdge,end)
		{
			Vector_3 vec = pHalfEdge->vertex()->point()-
				pHalfEdge->opposite()->vertex()->point();
			sum += std::sqrt(vec*vec);
			degree++;
		}
		return sum / (FT) degree;
	}
	
};

DGAL_END_NAMESPACE

#endif //DGAL_GL_POLYHEDRON_3_POLICY_H