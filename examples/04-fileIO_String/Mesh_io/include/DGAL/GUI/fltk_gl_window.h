#ifndef DGAL_FLTK_GL_WINDOW_H
#define DGAL_FLTK_GL_WINDOW_H

#pragma warning (push)
	#pragma warning(disable : 4083 4204 4244 4311 4312)
	#include <CGAL/Cartesian.h>
	#include <DGAL/Parameterization_factory.h>

	//for fltk
	#include <FL/Fl.H>
	#include <Fl/Fl_Gl_Window.h>
	#include <FL/Fl_Help_Dialog.H>

	#include <FL/gl.h>//lead to fatal error to use someboay.max() or std::max()

	#include <arcball/camera.h>
	#include <arcball/viewport.h>
	#include <arcball/arcball.h>
	#include <DGAL/UTIL/color.h>
	#include <Simple_gl_texture/Texture.H>
#pragma warning (pop)

DGAL_BEGIN_NAMESPACE

template
<
	class Polyhedron_3_
>
class Fltk_gl_window 
	: public Fl_Gl_Window
{
public:
	typedef typename Polyhedron_3_                 Polyhedron;
	typedef typename Polyhedron::Facet_iterator    Facet_iterator;
	typedef typename Polyhedron::Facet_handle      Facet_handle;
	typedef typename Polyhedron::Vertex_iterator   Vertex_iterator;
	typedef typename Polyhedron::Vertex_handle     Vertex_handle;
	typedef typename Polyhedron::Vertex_const_handle
		                                           Vertex_const_handle;
	typedef typename Polyhedron::Halfedge_iterator Halfedge_iterator;
	typedef typename Polyhedron::Halfedge_handle   Halfedge_handle;
	typedef typename Polyhedron::Edge_iterator     Edge_iterator;
	typedef typename Polyhedron::Aff_transformation_3     
												   Aff_transformation_3;
	typedef typename Polyhedron::FT                FT;
	typedef typename Polyhedron::Point_3           Point_3;
	typedef typename Polyhedron::Vector_3          Vector_3;
	typedef typename Polyhedron::Iso_cuboid_3      Iso_cuboid_3;

public:
	enum Pick_mode {
		PICK_VERTICES  = 1,
		PICK_EDGES     = 2,
		PICK_FACES     = 3
	};
public:
	Fltk_gl_window(int x, int y, int width, int height, const char* title=0)
		: Fl_Gl_Window(x, y, width, height, title),mesh_(0),domain_(0),tmp_mesh_(0)
	{
		init();
	}
	virtual ~Fltk_gl_window(){clear_mesh();}
protected:
	virtual void init();
	virtual void draw();
	virtual int handle(int event1);
	
public:	
	virtual void switch_mesh_domain()
	{
		if ( mesh_ == domain_)//mesh_ saves domain
		{
			mesh_ = tmp_mesh_;
		}
		else
		{
			mesh_ = domain_;
		}

		view_all(false);
	}

protected:
	void init_gl();
	void init_camera();
	void handle_mouse(int x, int y);
	void pick(int x, int y, int key_state);

public:	
	float view_all(bool check_first = true, bool trans=false);
	Fl_Widget* root(Fl_Widget* o) 
	{
		Fl_Widget* t=o->parent();
		while(t) {o=t;t=t->parent();}
		return o;
	}
	void clear_picked()
	{
		mesh_->picked_vertices().clear();
		mesh_->picked_edges().clear();
		mesh_->picked_facets().clear();
	}
		
	void clear_mesh()
	{
		if(mesh_ == domain_)
		{
			if(mesh_)    {delete mesh_;     mesh_=0;}
			if(tmp_mesh_){delete tmp_mesh_;   tmp_mesh_=0;}
			domain_=0;
		}
		else
		{
			if(mesh_)    {delete mesh_;     mesh_=0;}
			if(domain_)  {delete domain_;   domain_=0;}
			tmp_mesh_=0;
		}
	}

public:	
	virtual Polyhedron* mesh(){return mesh_;}
	void virtual mesh(Polyhedron* in)
	{
		clear_mesh();
		mesh_ = in;
		tmp_mesh_ = in;

		Vertex_const_handle vh( mesh_->vertices_begin( ));
		if ( vh->is_parameterized())
		{
			is_parameterized(true);
			domain_ = mesh_->build_domain();
		}
		else
			is_parameterized(false);
	}
	void set_mesh(Polyhedron* in){mesh_ = in;}
	Polyhedron* domain(){return domain_;}	

	Parameterize_option & parametrize_option(){return parametrize_option_;}
	bool is_parameterized(){
		return parametrize_option_.is_parameterized();
	}
	void is_parameterized(bool in){
		parametrize_option_.is_parameterized(in);
	}
	void set_texture(char* name, bool repeat_u=true, bool repeat_v=true){fetchTexture(name, repeat_u, repeat_v);}
	
	#define v3f(x) glVertex3fv(x)

	/* Draw a colored cube */
	void draw_cube(int wire= GL_POLYGON) 
	{
		float v0[3] = {0.0, 0.0, 0.0};
		float v1[3] = {1.0, 0.0, 0.0};
		float v2[3] = {1.0, 1.0, 0.0};
		float v3[3] = {0.0, 1.0, 0.0};
		float v4[3] = {0.0, 0.0, 1.0};
		float v5[3] = {1.0, 0.0, 1.0};
		float v6[3] = {1.0, 1.0, 1.0};
		float v7[3] = {0.0, 1.0, 1.0};

		glBegin(wire ? GL_LINE_LOOP : GL_POLYGON);
		glColor3ub(0,0,255);
		v3f(v0); v3f(v1); v3f(v2); v3f(v3);
		glEnd();
		glBegin(wire ? GL_LINE_LOOP : GL_POLYGON);
		glColor3ub(0,255,255); v3f(v4); v3f(v5); v3f(v6); v3f(v7);
		glEnd();
		glBegin(wire ? GL_LINE_LOOP : GL_POLYGON);
		glColor3ub(255,0,255); v3f(v0); v3f(v1); v3f(v5); v3f(v4);
		glEnd();
		glBegin(wire ? GL_LINE_LOOP : GL_POLYGON);
		glColor3ub(255,255,0); v3f(v2); v3f(v3); v3f(v7); v3f(v6);
		glEnd();
		glBegin(wire ? GL_LINE_LOOP : GL_POLYGON);
		glColor3ub(0,255,0); v3f(v0); v3f(v4); v3f(v7); v3f(v3);
		glEnd();
		glBegin(wire ? GL_LINE_LOOP : GL_POLYGON);
		glColor3ub(255,0,0); v3f(v1); v3f(v2); v3f(v6); v3f(v5);
		glEnd();
	}
	Pick_mode pick_mode(){return pick_mode_;}
	void pick_mode(Pick_mode in)
	{
		if ( pick_mode_ == PICK_VERTICES) mesh_->picked_vertices().clear();
		else if ( pick_mode_ == PICK_EDGES) mesh_->picked_edges().clear();
		else mesh_->picked_facets().clear();

		pick_mode_ = in;
	}

	bool smooth_shading(){return smooth_shading_;}		
	void smooth_shading(bool in)
	{
		smooth_shading_ = in;	
		if(smooth_shading_)	glShadeModel(GL_SMOOTH);
		else     			glShadeModel(GL_FLAT);
	}
	bool culling(){return culling_;}	
	void culling(bool in)
	{
		culling_ = in;	
		if(culling_)   glEnable(GL_CULL_FACE);
		else           glDisable(GL_CULL_FACE);
	}
	bool antialiasing(){return antialiasing_;}
	void antialiasing(bool in)
	{
		antialiasing_ = in;
		if(antialiasing_)
		{
			glEnable(GL_LINE_SMOOTH);
			glEnable(GL_BLEND);
			glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
			glHint(GL_LINE_SMOOTH_HINT,GL_NICEST);
			glLineWidth(antialiasing_edge_thickness_);//1.5f
		}
		else
		{
			glDisable(GL_LINE_SMOOTH);
			glDisable(GL_BLEND);
			glLineWidth(1.0f);
		}
	}
	
	bool show_texture(){ return show_texture_;}
	void show_texture(bool in)
	{
		show_texture_ = in; 		
	}
	float m_shininess(){ return m_shininess_;}
	void m_shininess(float in)
	{
		m_shininess_ = in; 
		float shininess[] = {m_shininess_};
		glMaterialfv( GL_FRONT, GL_SHININESS, shininess);
	}
	bool front_face(){ return front_face_;}
	void front_face(bool in)
	{
		front_face_ = in;
		if ( front_face_) glFrontFace(GL_CCW);
		else              glFrontFace(GL_CW);
	}
	ColorF l_ambient_color(){return l_ambient_color_;}
	ColorF l_diffuse_color(){return l_diffuse_color_;}
	ColorF l_specular_color(){return l_specular_color_;}
	ColorF m_ambient_color(){return m_ambient_color_;}
	ColorF m_diffuse_color(){return m_diffuse_color_;}
	ColorF m_specular_color(){return m_specular_color_;}
	ColorF m_emission_color(){return m_emission_color_;}
	void l_ambient_color(ColorF& in)
	{
		l_ambient_color_ = in;
		float lap[]  = { float(l_ambient_color_.r), float(l_ambient_color_.g), float(l_ambient_color_.b), float(l_ambient_color_.a)};
		glLightfv( GL_LIGHT0, GL_AMBIENT,  lap);
	}
	void l_diffuse_color(ColorF& in)
	{
		l_diffuse_color_ = in;
		float ldp[]  = { float(l_diffuse_color_.r), float(l_diffuse_color_.g), float(l_diffuse_color_.b), float(l_diffuse_color_.a)};
		glLightfv( GL_LIGHT0, GL_DIFFUSE,  ldp);
	}
	void l_specular_color(ColorF& in)
	{
		l_specular_color_ = in;
		float lsp[] = { float(l_specular_color_.r), float(l_specular_color_.g), float(l_specular_color_.b), float(l_specular_color_.a)};
		glLightfv( GL_LIGHT0, GL_SPECULAR, lsp);
	}
	void m_ambient_color(ColorF& in)
	{
		m_ambient_color_ = in;
		float mp[]  = { float(m_ambient_color_.r), float(m_ambient_color_.g), float(m_ambient_color_.b), float(m_ambient_color_.a)};
		glMaterialfv( GL_FRONT, GL_AMBIENT,   mp);
	}
	void m_diffuse_color(ColorF& in)
	{
		m_diffuse_color_ = in;
		float mdp[]  = { float(m_diffuse_color_.r), float(m_diffuse_color_.g), float(m_diffuse_color_.b), float(m_diffuse_color_.a)};
		glMaterialfv( GL_FRONT, GL_DIFFUSE,   mdp);
	}
	void m_specular_color(ColorF& in)
	{
		m_specular_color_ = in;
		float msp[]  = { float(m_specular_color_.r), float(m_specular_color_.g), float(m_specular_color_.b), float(m_specular_color_.a)};
		glMaterialfv( GL_FRONT, GL_SPECULAR,  msp);	
	}
	void m_emission_color(ColorF& in)
	{
		m_emission_color_ = in;
		float mep[]  = { float(m_emission_color_.r), float(m_emission_color_.g), float(m_emission_color_.b), float(m_emission_color_.a)};
		glMaterialfv( GL_FRONT, GL_EMISSION,  mep);
	}
	void background_color(){
		glClearColor(background_color_.r,background_color_.g,background_color_.b,1);
	}

	float light_x(){return l_pos_[0];}
	float light_y(){return l_pos_[1];}
	float light_z(){return l_pos_[2];}
	float light_w(){return l_pos_[3];}

	void light_x(float in){l_pos_[0] = in;}
	void light_y(float in){l_pos_[1] = in;}
	void light_z(float in){l_pos_[2] = in;}
	void light_w(float in){l_pos_[3] = in;}

protected:
	Polyhedron* mesh_;//represent a mesh	
	Polyhedron* tmp_mesh_;//another pointer point to the (*mesh_)
	Polyhedron* domain_; //represent the uv domain of the mesh
	Parameterize_option parametrize_option_;

	Pick_mode   pick_mode_;
	bool        smooth_shading_;
	bool        culling_;
	bool        antialiasing_;
	bool        show_texture_;
	bool        front_face_;

	CCamera     camera_;
	CViewport   viewport_;
	CArcball    arcball_;	

	float       l_pos_[4];//light position
	float       m_shininess_;
	ColorF      l_ambient_color_;//light ambient color
	ColorF      l_diffuse_color_;
	ColorF      l_specular_color_;
	ColorF      m_ambient_color_;//material ambient color
	ColorF      m_diffuse_color_;
	ColorF      m_specular_color_;
	ColorF      m_emission_color_;
	ColorF      background_color_;
	
public:	
	bool        first_view_;//for open mesh
	int         polygon_mode_;	
	bool        superimpose_edges_;
	bool        draw_voronoi_edges_;
	bool        superimpose_vertices_;	
	bool        draw_bounding_box_when_moving_;
	bool        draw_bounding_box_;
	float       super_point_radius_;
	float       super_edge_thickness_;
	float       antialiasing_edge_thickness_;//todo set by GUI
	bool        lighting_;

	// mouse
	bool        moving_;
	bool        l_button_down_;
	bool        r_button_down_;
	//color
	ColorF      edge_color_;//for super edge	
	ColorF      vertex_color_;//for super vertex
	ColorF      mesh_color_;//for mesh, include point, line and fill
	ColorF      picked_color_;
	ColorF      vertex_constrained_color_;
	ColorF      seam_color_;

	//ColorF env_color_;//environment color
};

// Setting up the camera to view the whole mesh.
//If trans == ture then the mesh is "unitized" first by translating it to the origin and scaling it to fit in a unit cube around the origin.
// Returns the scalefactor used.
template<class Polyhedron> inline float 
Fltk_gl_window<Polyhedron>::
view_all(bool check_first, bool trans)
{
	float result(1);

	if(mesh_ == 0)
		return result;

	if(!first_view_ && check_first) return result;
	first_view_ = false;

	// set up the camera to visualize the whole object
	mesh_->compute_bounding_box();
	if (trans)
	{
		Iso_cuboid_3 bbox = mesh_->bbox();

		/* calculate center of the model */
		FT cx = (bbox.xmin() + bbox.xmax()) / 2.0;
		FT cy = (bbox.ymin() + bbox.ymax()) / 2.0;
		FT cz = (bbox.zmin() + bbox.zmax()) / 2.0;

		Aff_transformation_3 tran(CGAL::TRANSLATION, Vector_3(-cx,-cy,-cz,1));

		/* calculate unitizing scale factor */
		FT x = bbox.xmax() - bbox.xmin();FT y = bbox.ymax() - bbox.ymin();FT z = bbox.zmax() - bbox.zmin();		
		FT scale = 2.0/std::sqrt(x*x+y*y+z*z);	
		result = (float)scale;
		Aff_transformation_3 sca(CGAL::SCALING, scale, 1);

		Aff_transformation_3 aff( sca * tran);
		mesh_->transform(aff);
		mesh_->compute_bounding_box();
	}

	CMatrix44 ArcballMatrix = arcball_.GetMatrix();

	CVector3d minBound, maxBound;
	minBound.Set(mesh_->xmin(),mesh_->ymin(),mesh_->zmin());
	maxBound.Set(mesh_->xmax(),mesh_->ymax(),mesh_->zmax());
	minBound = ArcballMatrix * minBound;
	maxBound = ArcballMatrix * maxBound;
	camera_.ViewAll(minBound[0],maxBound[0],minBound[1],
		maxBound[1],minBound[2],maxBound[2],viewport_);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	camera_.glDraw(viewport_);
	glMatrixMode(GL_MODELVIEW);

	return result;
}

template<class Polyhedron> inline void 
Fltk_gl_window<Polyhedron>::
init()         
{	
	pick_mode_ = PICK_VERTICES;
	first_view_ = false;
	front_face_ = true;

	l_pos_[0]=-5.0; l_pos_[1]=5.0; l_pos_[2]=5.0; l_pos_[3]=1.0;

	mode(FL_RGB | FL_ALPHA | FL_DEPTH | FL_DOUBLE);

	polygon_mode_ = GL_LINE;
	smooth_shading_ = true;
	culling_ = false;
	antialiasing_ = false;
	superimpose_edges_ = false;
	draw_voronoi_edges_ = false;
	superimpose_vertices_ = false;
	draw_bounding_box_when_moving_ = false;
	draw_bounding_box_ = false;

	moving_ = false;
	l_button_down_ = false;
	r_button_down_ = false;

	super_point_radius_ = 0.1f;
	super_edge_thickness_ = 1.5f;
	antialiasing_edge_thickness_ = 1.5f;

	background_color_.set        (0.0f, 0.0f, 0.0f);
	mesh_color_.set              (0.5f, 0.5f, 0.5f);
	edge_color_.set              (0.0f, 1.0f ,1.0f);
	vertex_color_.set            (1.0f, 1.0f, 0.0f);
	vertex_constrained_color_.set(0.0f, 0.0f ,1.0f);
	seam_color_.set              (0.0f, 1.0f ,0.0f);
	picked_color_.set            (1.0f, 0.0f, 0.0f);

	lighting_ = false;
	l_ambient_color_.set (0.0f, 0.0f, 0.0f);
	l_diffuse_color_.set (1.0f, 1.0f, 1.0f);
	l_specular_color_.set(1.0f, 1.0f, 1.0f);

	m_ambient_color_.set (0.11f, 0.06f, 0.11f);
	m_diffuse_color_.set (0.43f, 0.47f, 0.54f);
	m_specular_color_.set(0.33f, 0.33f, 0.52f);
	m_emission_color_.set(0.00f, 0.00f, 0.00f, 0.00f);
	m_shininess_ = 10;

	show_texture_ = false;
}

template<class Polyhedron> inline void 
Fltk_gl_window<Polyhedron>::
init_gl()
{
	background_color();//sets the clearing color
	glEnable(GL_DEPTH_TEST);

	glClearDepth(1.0f);									// Depth Buffer Setup
	glEnable(GL_DEPTH_TEST);							// Enables Depth Testing
	glDepthFunc(GL_LEQUAL);								// The Type Of Depth Testing To Do
	glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);	// Really Nice Perspective Calculations
	
	// lighting
	glEnable(GL_LIGHT0);
	glEnable(GL_LIGHTING);	
	glLightfv(GL_LIGHT0,GL_POSITION,l_pos_);
	glLightModelf(GL_LIGHT_MODEL_TWO_SIDE, 1.0f);//0.0-1.0

	// Lights properties
	float lap[]  = { float(l_ambient_color_.r), float(l_ambient_color_.g), float(l_ambient_color_.b), float(l_ambient_color_.a)};
	float ldp[]  = { float(l_diffuse_color_.r), float(l_diffuse_color_.g), float(l_diffuse_color_.b), float(l_diffuse_color_.a)};
	float lsp[] = { float(l_specular_color_.r), float(l_specular_color_.g), float(l_specular_color_.b), float(l_specular_color_.a)};

	glLightfv( GL_LIGHT0, GL_AMBIENT,  lap);
	glLightfv( GL_LIGHT0, GL_DIFFUSE,  ldp);
	glLightfv( GL_LIGHT0, GL_SPECULAR, lsp);

	//  material 
	float mp[]  = { float(m_ambient_color_.r), float(m_ambient_color_.g), float(m_ambient_color_.b), float(m_ambient_color_.a)};
	float mdp[]  = { float(m_diffuse_color_.r), float(m_diffuse_color_.g), float(m_diffuse_color_.b), float(m_diffuse_color_.a)};
	float msp[]  = { float(m_specular_color_.r), float(m_specular_color_.g), float(m_specular_color_.b), float(m_specular_color_.a)};
	float mep[]  = { float(m_emission_color_.r), float(m_emission_color_.g), float(m_emission_color_.b), float(m_emission_color_.a)};
	float shininess[] = {m_shininess_};
    glMaterialfv( GL_FRONT, GL_AMBIENT,   mp);
	glMaterialfv( GL_FRONT, GL_DIFFUSE,   mdp);
	glMaterialfv( GL_FRONT, GL_SPECULAR,  msp);	
    glMaterialfv( GL_FRONT, GL_EMISSION,  mep);
	glMaterialfv( GL_FRONT, GL_SHININESS, shininess);

	// Enable to set up color material
	glEnable(GL_COLOR_MATERIAL);
	// Set up the color material to be two sided and only use diffuse component
	glColorMaterial(GL_FRONT_AND_BACK, GL_DIFFUSE);

	// Set up the texture environment
	glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
	// Initially load the image in
	set_texture("textures", (bool)GL_REPEAT, (bool)GL_REPEAT);

	init_camera();
}

template<class Polyhedron> inline void 
Fltk_gl_window<Polyhedron>::
init_camera()
{
	// set viewport and camera
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	viewport_.SetOrigin(0,0);
	viewport_.SetSize(w(),h());
	camera_.SetHeightAngle(45.);
	camera_.SetPosition(0.,0.,5.);
	camera_.SetOrientation(0.,1.,0.,0);
	camera_.SetNearDistance(.1);
	camera_.SetFarDistance(1000.);
	viewport_.glDraw();
	camera_.glDraw(viewport_);
	glDrawBuffer(GL_BACK);
}

template<class Polyhedron> inline void 
Fltk_gl_window<Polyhedron>::
draw()
{
	if (!valid())//for window creat and resize 
	{
		init_gl();
	}// if 	

	if ( show_texture_) 
		glEnable(GL_TEXTURE_2D);
	else                
		glDisable(GL_TEXTURE_2D);
	
	// clears the color buffer and depth buffer using the current clearing color 
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT); 
	if (mesh_ == 0){	draw_cube();	return;	}

	view_all();// setup camera (only once)

	//begin of drawing
	glPushMatrix();//to avoid cumulate effect

	/////////////////////////////////////////////////////////////////////////
	//setup viewpoint from current arcball
	arcball_.glDraw();
	
	if(superimpose_edges_ || superimpose_vertices_)
	{
		// enable polygon offset
		glEnable(GL_POLYGON_OFFSET_FILL);
		glPolygonOffset(3.0f,1.0f);
	}

	// draw the bounding box 
	if((moving_ && draw_bounding_box_when_moving_) || 
		 draw_bounding_box_)
	{
		glColor3f(1.0f,0.0f,0.0f);
		glDisable(GL_LIGHTING);
		mesh_->draw_bounding_box(mesh_->bbox());
	}

	bool use_normals;
	// lighting option
	if(lighting_)
	{
		use_normals = true;
		glEnable(GL_LIGHTING);
		glLightfv(GL_LIGHT0,GL_POSITION,l_pos_);
	}
	else
	{
		use_normals = false;
		glDisable(GL_LIGHTING);
	}

	/////////////////////////////////////////////////////////////////////////
	//draw faces
	if(!moving_ || !draw_bounding_box_when_moving_)
	{
		glPolygonMode(GL_FRONT_AND_BACK, polygon_mode_);// polygon mode (GL_POINT, GL_LINE or GL_FILL)	

		glColor3d(mesh_color_.r,mesh_color_.g,mesh_color_.b);//set mesh color
		if ( is_parameterized() && show_texture_)
			mesh_->draw(smooth_shading_,use_normals, true);	
		else
			mesh_->draw(smooth_shading_,use_normals, false);
	}

	// disable lighting
	if(superimpose_vertices_ || superimpose_edges_ || pick_mode_>0)
		glDisable(GL_LIGHTING);

	//draw selected faces
	if(!moving_ || !draw_bounding_box_when_moving_)
	{
		glColor3d(picked_color_.r,picked_color_.g,picked_color_.b);
		mesh_->draw_picked_facets(smooth_shading_, use_normals, false);
	}

	// draw the mesh once again with a few options desactivated
	/////////////////////////////////////////////////////////////////////////
	//draw edges 	
	if(superimpose_edges_ && 
		 !(moving_ && draw_bounding_box_when_moving_))
	{
		// set line mode
		glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);

		glLineWidth(super_edge_thickness_);
		glColor3d(edge_color_.r,edge_color_.g,edge_color_.b);
		mesh_->draw_edges(draw_voronoi_edges_);
	} // end superimpose edges
	glLineWidth(float(super_edge_thickness_*1.5));
	glColor3d(picked_color_.r,picked_color_.g,picked_color_.b);
	mesh_->draw_picked_edges();

	// superimpose vertices
	if(superimpose_vertices_ &&
		 !(moving_ && draw_bounding_box_when_moving_))
	{
		//draw vertices
		glPolygonMode(GL_FRONT_AND_BACK,GL_POINT);// set point mode

		glColor3d(vertex_color_.r,vertex_color_.g,vertex_color_.b);
		mesh_->draw_vertices(this->super_point_radius_);
	} // end superimpose vertices
	glColor3d(picked_color_.r,picked_color_.g,picked_color_.b);
	mesh_->draw_picked_vertices(this->super_point_radius_*1.5);

	// disable polygon offset
	if(superimpose_vertices_ || superimpose_edges_)
		glDisable(GL_POLYGON_OFFSET_FILL);

	glPopMatrix();
	//end of drawing
}

template<class Polyhedron> inline int 
Fltk_gl_window<Polyhedron>::
handle(int et)
{
	static int x(this->x()), y(this->y()), w(this->w()), h(this->h());
	static int fx(root(this)->x()), fy(root(this)->y()), fw(root(this)->w()), fh(root(this)->h());
	static bool is_full(false);
	switch(et) 
	{
	case FL_FOCUS:
		Fl::focus(this);
		return 1;
	case FL_UNFOCUS :
		return 1;
	case FL_KEYBOARD:
	{
		switch (Fl::event_key()) 
		{
		case FL_Control_L & 'a'://Left Ctrl + a
			first_view_=true; redraw();// view all; 
			std::cout << "\"Left Ctrl + a\" pressed: view all." << std::endl;
			return 1;
		case FL_Escape:
			exit(0);
		case FL_F+1: //F1
			{
			std::cout << "\"F1\" pressed: show help." << std::endl;
			Fl_Help_Dialog* help = new Fl_Help_Dialog;
			help->load("help/help.html");
			help->show();
			}
			return 1;	
		case 'v':
			{
			std::cout << "\"v\" pressed: switch to pick vertice mode." << std::endl;
			pick_mode(PICK_VERTICES);
			redraw();
			}
			return 1;
		case 'e':
			{
			std::cout << "\"e\" pressed: switch to pick edge mode." << std::endl;
			pick_mode(PICK_EDGES);
			redraw();
			}
			return 1;
		case 'f':
			{
			std::cout << "\"f\" pressed: switch to pick face mode." << std::endl;
			pick_mode(PICK_FACES);
			redraw();
			}
			return 1;
		}
	}
		
	case FL_PUSH:
	{
	  int bt=Fl::event_button();
	  int x= Fl::event_x();
	  int y= Fl::event_y();
	  if(bt==FL_LEFT_MOUSE)
	  {
		  l_button_down_ = TRUE;
		  handle_mouse(x,y);
	  }
	  else if(bt==FL_RIGHT_MOUSE)
	  {
		  r_button_down_ = TRUE;
		  handle_mouse(x,y);
	  }
	  Fl::focus(this);
	  return 1;
	}	
	case FL_DRAG:
	{			 
	  int bt=Fl::event_button();
	  int x= Fl::event_x();
	  int y= Fl::event_y();
	  if(l_button_down_ || r_button_down_)
	  {
		  CVector3d vec = arcball_.Intersect(x,viewport_.yRes()-y,camera_,viewport_);
		  arcball_.Motion(vec);
		  moving_ = true;
		  redraw();
	  }
	  else
		  moving_ = false;

	  return 1;
	}
	case FL_RELEASE:   
	{
	  int bt=Fl::event_button();
	  int x= Fl::event_x(); int y= Fl::event_y();
	  l_button_down_ = FALSE;	
	  r_button_down_ = FALSE;	
	  moving_ = false;

	  if (bt==FL_LEFT_MOUSE)
		  handle_mouse(x,y);
	  else if (bt==FL_RIGHT_MOUSE)	
	  {
		  handle_mouse(x,y);
		  if ( Fl::event_state(FL_CTRL) ) 
		  {
			  this->pick(x,y, FL_CTRL);
		  }
		  else if ( Fl::event_state(FL_SHIFT) ) 
		  {
			  this->pick(x,y, FL_SHIFT);
		  }
	  }
	  
	  redraw();
	  return 1;
	}
	default:    
		return Fl_Gl_Window::handle(et);
	}
}

template<class Polyhedron> inline void 
Fltk_gl_window<Polyhedron>::
pick(int x, int y, int key_state)
{
	if(mesh_ ==0) return;

	GLuint	buffer[512];
	GLint	hits;

	// The Size Of The Viewport. [0] Is <x>, [1] Is <y>, [2] Is <length>, [3] Is <width>
	GLint	viewport[4];
	GLdouble modelview[16];
	GLdouble projection[16];

	// This Sets The Array <viewport> To The Size And Location Of The Screen Relative To The Window
    glGetIntegerv(GL_VIEWPORT, viewport);
	glGetDoublev(GL_MODELVIEW_MATRIX, modelview);
	glGetDoublev(GL_PROJECTION_MATRIX, projection);

	glSelectBuffer(512, buffer);	// Tell OpenGL To Use Our Array For Selection
	// Puts OpenGL In Selection Mode. Nothing Will Be Drawn.  Object ID's and Extents Are Stored In The Buffer.
	(void) glRenderMode(GL_SELECT);

	glInitNames();                  // Initializes The Name Stack
	glPushName((unsigned)-1);	    // Push 0 (At Least One Entry) Onto The Stack

	glMatrixMode(GL_PROJECTION);	// Selects The Projection Matrix					
	glPushMatrix();				    // Push The Projection Matrix				
		glLoadIdentity();		    // Resets The Matrix
		// This Creates A Matrix That Will Zoom Up To A Small Portion Of The Screen, Where The Mouse Is.
		gluPickMatrix((GLdouble) x, (GLdouble) (viewport[3]-y), 1.0f, 1.0f, viewport);
		glMultMatrixd(projection);  // Inherit projection from current projection matrix

		glMatrixMode(GL_MODELVIEW);	// Select The Modelview Matrix
		glLoadIdentity();           // Resets The Matrix
		glMultMatrixd(modelview);   // Inherit view from current view matrix

		glPushMatrix();			
			arcball_.glDraw();  //setup viewpoint from current arcball
			// Render The Targets To The Selection Buffer
			if(pick_mode_ == PICK_VERTICES)
			{
				mesh_->draw_vertices(super_point_radius_, pick_mode_);
			}
			else if (pick_mode_ == PICK_EDGES)
			{
				glLineWidth(super_edge_thickness_);
				mesh_->draw_edges(false, pick_mode_);
			}
			else if (pick_mode_ == PICK_FACES)
			{
				mesh_->draw(false, false, false, pick_mode_);
			}
			//else 
		glPopMatrix();			
		
		glMatrixMode(GL_PROJECTION);
	glPopMatrix();
	glMatrixMode(GL_MODELVIEW);	

	hits=glRenderMode(GL_RENDER);
	std::cout << "pick mode: " << pick_mode_ << "; and hits = " << hits << std::endl;

	if (hits > 0)												// If There Were More Than 0 Hits
	{
		int	choose = buffer[3];									// Make Our Selection The First Object
		int depth = buffer[1];									// Store How Far Away It Is 

		for (int loop = 1; loop < hits; loop++)					// Loop Through All The Detected Hits
		{
			// If This Object Is Closer To Us Than The One We Have Selected
			if (buffer[loop*4+1] < GLuint(depth))
			{
				choose = buffer[loop*4+3];						// Select The Closer Object
				depth = buffer[loop*4+1];						// Store How Far Away It Is
			}       
		}

		if(pick_mode_ == PICK_VERTICES)
		{
			std::cout << "Picked vertex name is " << choose << std::endl;
			int numVertices = mesh_->size_of_vertices();
			if (choose >= numVertices ) return;

			int i(0); 
			for(Vertex_iterator vi = mesh_->vertices_begin();vi!=mesh_->vertices_end(); ++vi,++i)
			{
				if ( i == choose)
				{				
					Vertex_handle vh = vi;
					if ( key_state == FL_CTRL)
						mesh_->picked_vertices().clear();

					mesh_->add_picked_vertex(vh);
					break;
				}
			 }
		}
		else if (pick_mode_ == PICK_EDGES)
		{
			std::cout << "Picked edge name is " << choose << std::endl;

			int numEdges = mesh_->size_of_halfedges()/2;
			if (choose >= numEdges ) return;
			
			int i(0);
			for(Edge_iterator h = mesh_->edges_begin(); h != mesh_->edges_end(); ++h,++i)
			{
				if ( i == choose)
				{	
					Halfedge_handle hh = h;
					if ( key_state == FL_CTRL)
						mesh_->picked_edges().clear();

					mesh_->add_picked_edge(hh);
					break;
				}
			}
		}
		else if (pick_mode_ == PICK_FACES)
		{
			std::cout << "Picked face name is " << choose << std::endl;
			int numFacets = mesh_->size_of_facets();
			if (choose >= numFacets ) return;

			int i(0);
			for(Facet_iterator pFacet = mesh_->facets_begin();pFacet != mesh_->facets_end(); ++pFacet, ++i)
			{
				if ( i == choose)
				{
					Facet_handle fh = pFacet;
					if ( key_state == FL_CTRL)
						mesh_->picked_facets().clear();
					mesh_->add_picked_facet(fh);
					break;
				}
			}
		}
	}
}

//flag: 1 for left button; 2 for right button; 3 for both button
template<class Polyhedron> inline void 
Fltk_gl_window<Polyhedron>::
handle_mouse(int x, int y)
{ 	
	CVector3d vec = arcball_.Intersect(x,viewport_.yRes()-y,camera_,viewport_);
	arcball_.EndDrag(vec);
	arcball_.SetMode(l_button_down_+2*r_button_down_);
	vec = arcball_.Intersect(x,viewport_.yRes()-y,camera_,viewport_);
	arcball_.BeginDrag(vec);	
}

DGAL_END_NAMESPACE

#endif//DGAL_FLTK_GL_WINDOW_H
