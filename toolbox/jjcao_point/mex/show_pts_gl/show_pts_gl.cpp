#pragma warning(disable : 4244 4290 4305 4800 4996)

#include <yvals.h>
#if (_MSC_VER >= 1600)
#define __STDC_UTF_16__
#endif

#include <mex.h>
#pragma comment(lib, "libmex.lib")
#pragma comment(lib, "libmx.lib")

#include <iostream>
#include <sstream>
#include <numeric>
#include <GL/glut.h>
#include "myutil.h"
#include "MyHeaps.h"

using namespace std;

void print(const char * msg){
	mexWarnMsgTxt(msg);
}

// for picking and draging, added by jjcao
#define BUFSIZE 512
GLuint selectBuf[BUFSIZE];
double snapDist(0xffffffff);

// GUI Variables
float clear_color[3] = { 1, 1, 1 }; //white
int window_height = 800; //1100
int window_width  = 800; //1100
double eyepos[3] = {0,0,1}; // from what position I look at the object

// Lights
float light_ambient[4] = { 0.3f, 0.3f, 0.3f, 1.0f };
float light_diffuse[4] = { 0.6f, 0.6f, 0.6f, 1.0f };
float light0_position[4] = { 0.0f, 0.0f, 5.0f, 0.0f };
float light1_position[4] = { 1.0f, 1.0f, 1.0f, 0.0f };

// Colors
//float surfelcolor[4] = {1, .73, .0, 0.5f}; // gold
//float surfelcolor[4] = {0.275, .337, .60, 0.5f}; // blue
float surfelcolor[4] = {1.00, .65, .35, 0.5f}; // orange
float blue[3] = {0, 0, 1};
float red[3] = {1, 0, 0};
float vertColor[3] = {1, .8, .8};//added by jjcao
float edgeColor[3] = {1, 0, 0};//added by jjcao
float edgePickColor[3] = {0, 0, 0.5}; // added by jjcao
float vertPickColor[3] = {0, 0, 1}; // added by jjcao
// Modelview
//double modelview_rot[16] = { 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1 }; //standard
double modelview_rot[16] = { -0.998294, 0.0230628, -0.0535224, 0,
-0.0174278, 0.758218, 0.651759, 0,
0.0556148, 0.651586, -0.756529, 0, 0, 0, 0, 1 }; //for Pecaso
//double modelview_rot[16] = {-0.886703, -0.156514, -0.435038, 0,
//										-0.26335,   0.944216,  0.197471,  0,
//										0.379864,  0.289753, -0.878491,  0, 0, 0, 0, 1 };//for fertility_1n
//double modelview_rot[16] = {0.589297, -0.325771, 0.739324, 0,
//										-0.72799, 0.182738, 0.660784, 0,
//										-0.350367, -0.927621, -0.129471,0,0, 0, 0, 1};//for runman1s1 reconstruction
//double modelview_rot[16] = {0.751328, -0.352342, 0.557993, 0,
//										-0.597814, -0.00526998, 0.801616, 0,
//										-0.279503, -0.9355857, -0.214595, 0, 0, 0, 0, 1};//for runman1s1 skeleton edit
// modelview data
enum ModelViewStatus {NORMAL, ROTATING, TRANSLATING, ZOOMING, PICKING} modelview_status;
int pick_vid = -1;// added by jjcao
double modelview_sx = 0, modelview_sy = 0; // modelview I/O starting position
double modelview_sz = 0; // added by jjcao
double modelview_tx=0, modelview_ty=0; // x/y axis translation
double modelview_zoom = 1; // zoom

// Flow control variables
bool showsplats			= false;
bool tranparencyEnabled = false;
int  show_backpointing = 0;
int  colorwhat         = 0;
double markersize      = 0.033241;

class Pcloud;
// Data models
Pcloud* pcloud(0); // the point cloud
double *pickedPtsIdx; // picked points' id, results of the program

// Local classes
class Pcloud{
public:
	vector<double> points;	// point  coordinates
	vector<double> normals; // normals
	vector<double> colors;   // color of cloud
	int 	npoints; 		// number of points
	int 	ndims;   		// dimensionality of the points

	// Constructor
	// not handle radius by now.
	Pcloud(int npts, const double *pts, const double *color=0, const double *normal=0, const double *radius=0, int radius_dim=0):npoints(npts),ndims(3){ 
		this->points.resize(npoints*ndims);

		for (int i = 0; i < npoints; ++i){
			setPoint(i, pts[i], pts[npoints + i], pts[2 * npoints + i]);
		}

		if (color){
			this->colors.resize(npoints*ndims);
			for (int i = 0; i < npoints; ++i){
				setColor(i, color[i], color[npoints + i], color[2 * npoints + i]);
			}
		}
		if (normal){
			this->normals.resize(npoints*ndims);
			for (int i = 0; i < npoints; ++i){
				setNormal(i, normal[i], normal[npoints + i], normal[2 * npoints + i]);
			}
		}

		//normalize();
	}

	int getNPoints(){return this->npoints;}
	/// Set coordinates of point "i" with x,y,z
	void setPoint( unsigned int i, double x, double y, double z ){
		assert( i<points.size() );
		points[i*ndims + 0] = x;
		points[i*ndims + 1] = y;
		points[i*ndims + 2] = z;
	}
	// added by jjcao
	vector<double> getPoint(int i){
		assert( i<npoints );
		vector<double> p(3);
		unsigned int j(i*ndims);
		p[0] = points[j];
		p[1] = points[j+1];
		p[2] = points[j+2];
		return p;
	}
	/// Set normal for point "i" with vx,vy,vz
	void setNormal( unsigned int i, double vx, double vy, double vz ){
		assert( i<normals.size() );
		normals[i*ndims + 0] = vx;
		normals[i*ndims + 1] = vy;
		normals[i*ndims + 2] = vz;
	}
	// added by jjcao
	void normalize(){// scale to unitBox and move to origin		
		double minx(0xffffffff), maxx(-minx);
		double miny(0xffffffff), maxy(-miny);
		double minz(0xffffffff), maxz(-minz);
		double x,y,z;

		for (int i = 0; i < this->npoints; ++i){
			x = this->points[i*ndims];
			y = this->points[i*ndims+1];
			z = this->points[i*ndims+2];
			if (x<minx) minx = x;
			if (x>maxx) maxx = x;
			if (y<miny) miny = y;
			if (y>maxy) maxy = y;
			if (z<minz) minz = z;
			if (z>maxz) maxz = z;
		}

		vector<double> m1(3), m2(3), c(3), d(3);
		m1[0] = minx; m1[1] = miny; m1[2] = minz;
		m2[0] = maxx; m2[1] = maxy; m2[2] = maxz;

		sum(m1,m2,c);
		scale(c, 0.5);
		diff(m2,m1,d);
		double s = max(d[0], d[1]);
		s = max(s, d[2]);
		s = 1.6 / s;

		int tmp;
		for (int i = 0; i < this->npoints; ++i){
			tmp = i*ndims;
			this->points[tmp] = (this->points[tmp]-c[0])*s; ++tmp;
			this->points[tmp] = (this->points[tmp]-c[1])*s; ++tmp;
			this->points[tmp] = (this->points[tmp]-c[2])*s; 
		}
	}
	// added by jjcao
	bool hasNormal(){
		//std::vector<double>::iterator it = this->normals.begin();
		//++it; ++it; ++it;
		//double res = std::accumulate(this->normals.begin(), it, 0.0);
		double res = std::accumulate(this->normals.begin(), this->normals.end(), 0.0);
		return bool(res);
	}
	void setColor(unsigned int i, double cx, double cy, double cz ){
		assert( i<colors.size() );
		colors[i*ndims + 0] = cx;
		colors[i*ndims + 1] = cy;
		colors[i*ndims + 2] = cz;
	}
	bool hasNormals() const{
		return normals.size();
	}
	inline double operator()(int i, int dim){
		assert(i<npoints );
		return points[ i*ndims + dim ];
	}
	void draw_transp(double disksize, int drawbackpoint=0, int colorwhat=0){
		glEnable(GL_BLEND);

		/// create a sorting field for the points to be drawn
		MinHeap<double> heap( npoints );
		double mv[16], pj[16]; int mw[16];
		glGetDoublev(GL_MODELVIEW_MATRIX, mv);
		glGetDoublev(GL_PROJECTION_MATRIX, pj);
		glGetIntegerv(GL_VIEWPORT, mw);
		vector<double> depth( npoints, 0 );
		vector<double> point(3,0);
		vector<double> wpoint(3,0);
		for (int i = 0; i < npoints; i++) {
			gluProject( points[ i*ndims + 0 ],
				points[ i*ndims + 1 ],
				points[ i*ndims + 2 ],
				mv, pj, mw,
				&wpoint[0], &wpoint[1], &wpoint[2]);
			depth[i] = -wpoint[2];
			heap.push( depth[i], i );
		}
		vector<int> srtidxs(npoints, 0);
		heap.heapsort( srtidxs );

		// extract camera vector
		float model[16];
		glGetFloatv (GL_MODELVIEW_MATRIX, model);
		vector<double> camera(3,0);
		camera[0] = model[2];
		camera[1] = model[6];
		camera[2] = model[10];

		// COLOR/QUADRIC SETUP
		//float diskfront[4] = {.81, .42, .11, 0.5f};
		float diskback[4] = {0,0,0,1};
		GLUquadricObj *q = gluNewQuadric();
		gluQuadricNormals (q,GLU_TRUE);
		glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, surfelcolor);
		glMaterialfv(GL_BACK, GL_AMBIENT_AND_DIFFUSE, diskback);

		/// DRAW BACKFACING
		vector<double> zaxis (3,0); zaxis[2] = 1;
		vector<double> normal(3,0);
		vector<double> crossv(3,0);
		double theta = 0;
		for (int i=0,j; i < npoints; i++){
			j = srtidxs[i];
			normal[ 0 ] = normals[ j*ndims + 0 ];
			normal[ 1 ] = normals[ j*ndims + 1 ];
			normal[ 2 ] = normals[ j*ndims + 2 ];
			if( dot(camera, normal)>=0 )
				continue;

			// estimate rotation parameters
			theta = acos( dot( zaxis, normal ) ) * 180 / PI;
			cross( zaxis, normal, crossv );
			glPushMatrix();
			glTranslatef( (*this)(j,0), (*this)(j,1), (*this)(j,2) );
			glRotatef( theta, crossv[0], crossv[1], crossv[2] );
			gluDisk( q, 0, disksize, 20, 1 );
			glPopMatrix();
		}

		// setup stencil
		glEnable(GL_STENCIL_TEST);
		glStencilFunc( GL_EQUAL, 0, 1 );
		glStencilOp( GL_KEEP, GL_KEEP, GL_INCR );

		/// DRAW FRONTFACING
		glEnable(GL_DEPTH_TEST);
		glDepthMask(GL_FALSE);


		for (int i=0,j; i < npoints; i++){
			j = srtidxs[npoints-i-1];
			normal[ 0 ] = normals[ j*ndims + 0 ];
			normal[ 1 ] = normals[ j*ndims + 1 ];
			normal[ 2 ] = normals[ j*ndims + 2 ];
			if( dot(camera, normal)<0 )
				continue;

			// estimate rotation parameters
			theta = acos( dot( zaxis, normal ) ) * 180 / PI;
			cross( zaxis, normal, crossv );
			glPushMatrix();
			glTranslatef( (*this)(j,0), (*this)(j,1), (*this)(j,2) );
			glRotatef( theta, crossv[0], crossv[1], crossv[2] );
			gluDisk( q, 0, disksize/2, 20, 1 );
			glPopMatrix();
		}
		glEnable(GL_DEPTH_TEST);
		glDepthMask(GL_TRUE);

		// turn off stencil
		glDisable(GL_STENCIL_TEST);

		// turn off transparent display
		glDisable(GL_BLEND);
	}

	void drawPickedPoints(float disksize, int drawbackpoint=0, int colorwhat=0, bool showsplats=0){
		if (pick_vid < 0) return;

		GLUquadricObj *q = gluNewQuadric();
		gluQuadricNormals (q,GLU_FALSE);

		vector<double> v = getPoint(pick_vid);

		// Are we using splats?
		if( hasNormals() && showsplats ){
			// retrieve modelview matrix
			float model[16];
			glGetFloatv (GL_MODELVIEW_MATRIX, model);
			vector<double> camera(3,0);
			camera[0] = model[2];
			camera[1] = model[6];
			camera[2] = model[10];

			vector<double> zaxis (3,0); zaxis[2] = 1;
			vector<double> normal(3,0);
			vector<double> crossv(3,0);
			double theta = 0;

			normal[ 0 ] = normals[ pick_vid*ndims + 0 ];
			normal[ 1 ] = normals[ pick_vid*ndims + 1 ];
			normal[ 2 ] = normals[ pick_vid*ndims + 2 ];
			theta = acos( dot( zaxis, normal ) ) * 180 / PI;
			cross( zaxis, normal, crossv );
			if( drawbackpoint==2 && dot(camera, normal)<0 )
				return;

			glPushMatrix();
			glTranslatef( v[0], v[1], v[2] );	
			glRotatef( theta, crossv[0], crossv[1], crossv[2] );
			glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, vertPickColor);
			gluDisk( q, 0, disksize*0.7, 20, 1 );
			glPopMatrix();
		}
		// Are we using spheres?
		if( showsplats == false ){
			glPushMatrix();
			glTranslatef( v[0], v[1], v[2] );	
			glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, vertPickColor);
			gluSphere(q, disksize*0.7, 10, 10);
			glPopMatrix();
		}
		
	}	
	// added by jjcao, for picking points
	void drawPoints(float cylindersize=.01){
		GLUquadricObj *q = gluNewQuadric();
		gluQuadricNormals (q,GLU_FALSE);

		// For every vertices in the graph
		for (int i = 0; i < this->npoints; ++i) {
			vector<double> v = getPoint(i);

			glPushMatrix();
			glPushName(i);
			glTranslatef( v[0], v[1], v[2] );	
			gluSphere(q, 0.5*cylindersize, 10, 10);
			glPopName();
			glPopMatrix();
		}
	}
	// draw a point cloud using OpenGLz
	void draw(double disksize, int drawbackpoint=0, int colorwhat=0, bool showsplats=0 ){
		// DEBUG: draw a simple sphere
		//		GLUquadricObj *q = gluNewQuadric();
		//		gluQuadricNormals (q,GLU_TRUE);
		//		glColor3f(1,1,1);

		// retrieve modelview matrix
		float model[16];
		glGetFloatv (GL_MODELVIEW_MATRIX, model);
		vector<double> camera(3,0);
		camera[0] = model[2];
		camera[1] = model[6];
		camera[2] = model[10];

		// draw a disk
		// float diskfront[3] = {.81, .42, .11};
		float diskback[3] = {0,0,0};
		GLUquadricObj *q = gluNewQuadric();
		gluQuadricNormals (q,GLU_TRUE);
		glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, surfelcolor);
		if( drawbackpoint == 0 )
			glMaterialfv(GL_BACK, GL_AMBIENT_AND_DIFFUSE, diskback);
		else if( drawbackpoint == 1 )
			glMaterialfv(GL_BACK, GL_AMBIENT_AND_DIFFUSE, surfelcolor);

		vector<double> zaxis (3,0); zaxis[2] = 1;
		vector<double> normal(3,0);
		vector<double> crossv(3,0);
		double theta = 0;
		for (int i = 0; i < npoints; i++){
			// Need to apply per splat color?
			if( colorwhat==1 ){
				float currcolor[3];
				currcolor[0] = colors[ i*ndims + 0 ];
				currcolor[1] = colors[ i*ndims + 1 ];
				currcolor[2] = colors[ i*ndims + 2 ];
				glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, currcolor);
			}

			// Are we using splats?
			if( hasNormals() && showsplats ){
				normal[ 0 ] = normals[ i*ndims + 0 ];
				normal[ 1 ] = normals[ i*ndims + 1 ];
				normal[ 2 ] = normals[ i*ndims + 2 ];
				theta = acos( dot( zaxis, normal ) ) * 180 / PI;
				cross( zaxis, normal, crossv );
				if( drawbackpoint==2 && dot(camera, normal)<0 )
					continue;

				glPushMatrix();
				glTranslatef( (*this)(i,0), (*this)(i,1), (*this)(i,2) );
				glRotatef( theta, crossv[0], crossv[1], crossv[2] );
				gluDisk( q, 0, disksize/2, 20, 1 );
				glPopMatrix();
			}
			// Are we using spheres?
			if( showsplats == false ){
				glPushMatrix();
				glTranslatef( (*this)(i,0), (*this)(i,1), (*this)(i,2) );
				gluSphere (q, disksize/2, 10, 10);
				glPopMatrix();
			}
		}
	}
};

void displayCallback() {
	// clear background
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);

	// ANIMATION CODE: Here is best place to change view if you want to spin around
	//glMatrixMode(GL_PROJECTION);
	//glLoadIdentity();
	//glOrtho(-1, 1, -1, 1, -10, 10);
	//gluLookAt( eyepos[0],eyepos[1],eyepos[2],0,0,0,0,1,0 );
	// END OF ANIMATION CODE

	// model drawing
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glTranslatef( modelview_tx, modelview_ty, 0 );
	glScalef( modelview_zoom, modelview_zoom, modelview_zoom );

	// ANIMATION CODE: Here is best place to change view if you want to spin the object around
	gluLookAt( eyepos[0],eyepos[1],eyepos[2],0,0,0,0,1,0 );
	// END OF ANIMATION CODE

	glMultMatrixd( modelview_rot );

	// this is for the example of rotational symmetry video
	//gluLookAt( eyepos[0],eyepos[1],eyepos[2],0,0,0,0,1,0 );

	//draw_axis_simple(); // added by jjcao

	// Draw the point cloud
	if( tranparencyEnabled )
		pcloud->draw_transp( markersize, show_backpointing, colorwhat );
	else
		pcloud->draw( markersize, show_backpointing, colorwhat, showsplats );
	pcloud->drawPickedPoints(markersize, show_backpointing, colorwhat, showsplats );

	// Draw the axis
	//draw_axis();
	glutSwapBuffers();
}

//added by jjcao
void startPicking(int cursorX, int cursorY) {
	// The Size Of <viewport>. [0] Is <x>, [1] Is <y>, [2] Is <length>, [3] Is <width>
	GLint	viewport[4];

	glSelectBuffer(BUFSIZE,selectBuf);
	glRenderMode(GL_SELECT);

	// This sets  <viewport> to the Size and Location Of The Screen Relative To The Window
	glGetIntegerv(GL_VIEWPORT, viewport);

	glMatrixMode(GL_PROJECTION);
	glPushMatrix();

	glGetIntegerv(GL_VIEWPORT,viewport);
	gluPickMatrix(cursorX,viewport[3]-cursorY,2,2,viewport);

	glMatrixMode(GL_MODELVIEW);	// Select The Modelview Matrix

	glInitNames();
}
// added by jjcao
void processHits (GLint hits, GLuint buffer[])
{
	GLuint names, *ptr, minZ,*ptrNames, numberOfNames;

	//cout << "hits =" << hits << endl;
	ptr = (GLuint *) buffer;
	minZ = 0xffffffff;
	for (int i = 0; i < hits; i++) {	
		names = *ptr;
		ptr++;
		if (*ptr < minZ) {
			numberOfNames = names;
			minZ = *ptr;
			ptrNames = ptr+2;
		}

		ptr += names+2;
	}
	//cout << "The closest hit names are ";
	ptr = ptrNames;
	pick_vid = *ptr;
	stringstream ss;	ss << "point picked: " << pick_vid << " points size: " << pcloud->getNPoints();
	print(ss.str().c_str());

	modelview_sz = minZ;

	if (numberOfNames<1){
		pick_vid = -1;
	}
	//for (unsigned int j = 0; j < numberOfNames; j++,ptr++) {
	//	cout << *ptr << ", ";
	//}
	//cout << endl;

	glutPostRedisplay();
}
// added by jjcao
void stopPicking() {
	int hits;

	// restoring the original projection matrix
	glMatrixMode(GL_PROJECTION);
	glPopMatrix();
	glMatrixMode(GL_MODELVIEW);
	glFlush();

	// returning to normal rendering mode
	hits = glRenderMode(GL_RENDER);

	// if there are hits process them
	if (hits != 0)
		processHits(hits,selectBuf);
}

// added by jjcao
void picking()
{
	// Render The Targets To The Selection Buffer
	pcloud -> drawPoints();
}
void mouseCallback(int button, int state, int x, int y) {
	//cout << "mousestat is " << state <<  " down: " << GLUT_DOWN << endl;
	//cout << x << ", " << y << endl; // added by jjcao

	// change modelview status
	if (state == GLUT_UP){
		// begin of added by jjcao
		if (modelview_status == PICKING && pick_vid>-1){
			vector<double> v = pcloud->getPoint(pick_vid);
			// drag, merge, or separate?
			int mod = glutGetModifiers();			
			if (mod == GLUT_ACTIVE_SHIFT){// merge								
			}else if (mod== GLUT_ACTIVE_CTRL){// separate, add a vertex				
			}
		}
		// end of added by jjcao

		modelview_status = NORMAL;
	}
	else if (state == GLUT_DOWN) {
		//cout << "button down: " << button << endl;

		if (button == GLUT_RIGHT_BUTTON){
			//cout << "left button" << endl;
			modelview_status = TRANSLATING;
		}
		else if (button == GLUT_LEFT_BUTTON){
			//cout << "right button" << endl;
			modelview_status = ROTATING;
		}
		else{ // added by jjcao
			//cout << "middle button" << endl;
			modelview_status = PICKING;

			startPicking(x,y);
			picking();
			stopPicking();
		}
		// set starting point
		modelview_sx = x;
		modelview_sy = y;
	}
}

void motionCallback(int x, int y) {
	double dx, dy, ang, len, v[3], rot[3];
	double z_axis[3] = { 0, 0, 1 };

	if (modelview_status != NORMAL) {

		// Get movement differentials
		dx = x - modelview_sx;
		dy = y - modelview_sy;

		//cout << "motion with status: " << modelview_status << endl;

		// If difference is significant
		if ( (dx != 0.0) || (dy != 0.0)) {
			if (modelview_status == ROTATING) {
				v[0] = dx;
				v[1] = -dy;
				v[2] = 0.0;

				len = length(v);
				ang = -((len * PI *10 )/180.0);

				normalize(v);
				cross(v, z_axis, rot);
				modelview_sx = x;
				modelview_sy = y;

				glLoadIdentity();
				glRotated(ang, rot[0], rot[1], rot[2]);
				glMultMatrixd(modelview_rot);
				glGetDoublev(GL_MODELVIEW_MATRIX, modelview_rot);
			}
			if (modelview_status == TRANSLATING) {
				modelview_sx = x;
				modelview_sy = y;

				modelview_tx += dx / window_width;
				modelview_ty -= dy / window_height;
			}
			if (modelview_status == PICKING && pick_vid>-1) {
				GLdouble modelMatrix[16];
				GLdouble projMatrix[16];
				int viewport[4];		

				glGetDoublev(GL_MODELVIEW_MATRIX,modelMatrix);				
				glGetDoublev(GL_PROJECTION_MATRIX,projMatrix);				
				glGetIntegerv(GL_VIEWPORT,viewport);

				//////////////////////////////////////////////////////////////////////////
				// modelview_sx,modelview_sy are the previous mouse position. modelview_sz is the
				//   Z value of the selected object. Together, they form the window
				//   coordinate position of the mouse in 3D window coordinate space.
				//   Vack-transform them into object coordinates.
				vector<double>	sp(3), p(3);
				GLint projected;
				projected = gluUnProject( modelview_sx, viewport[3]-modelview_sy, modelview_sz/(double)0xffffffff,
					modelMatrix,projMatrix,viewport,&sp[0], &sp[1],&sp[2]);
				assert( projected == GL_TRUE );

				// So the same with the current window coordnate mouse position.
				projected = gluUnProject( x, viewport[3]-y, modelview_sz/(double)0xffffffff,
					modelMatrix,projMatrix,viewport,&p[0], &p[1],&p[2]);
				assert( projected == GL_TRUE );

				// (p-sp) is the object coordinate delta motion of the mouse.
				//   Set this in the scene to reposition the selected object.
				vector<double> v = pcloud->getPoint(pick_vid);
				sum(v, p);
				diff(v,sp);
				pcloud->setPoint(pick_vid, v[0], v[1], v[2]);

				//////////////////////////////////////////////////////////////////////////
				modelview_sx = x;
				modelview_sy = y;
			}
			glutPostRedisplay();
		}//if ( (dx != 0.0) || (dy != 0.0)) 
	}
}

void reshapeCallback(int x_new, int y_new) {
	// update window dimension (1/1 ratio)
	window_width  = x_new;
	window_height = y_new; //changed by jjcao from x_new to y_new

	// re-set viewport
	glViewport(0, 0, window_width, window_height);

	// refresh
	glutPostRedisplay();
}
// keyboard
void keyboardCallback(unsigned char key, int x, int y) {
	stringstream ss;
	switch (key) {
	case 'q':
		exit(0);
		print("viewer quit.\n");
		break;
	case 's':
		showsplats = !showsplats;
		break;
	case '=':
		modelview_zoom *= 1.1;
		break;
	case 'm':
		print("modelview");

		for (int i = 0; i < 16; i++)
			ss << modelview_rot[i] << " ";

		print(ss.str().c_str());
		break;
	case '-':
		modelview_zoom *= 0.9;
		break;
	case 't':
		tranparencyEnabled = !tranparencyEnabled;
		break;
	case 'c':
		colorwhat = (colorwhat+1)%3;
		break;
	case 'h':
		print("EXAMPLE SYNTAX: pickedPtsIdx = show_pts_gl(verts) \n");
		print("EXAMPLE SYNTAX: pickedPtsIdx = show_pts_gl(verts, colors, normals, radius) \n");
		print("h: show this help\n");
		print("= / -: zoom/unzoom the model\n");
		print("[ / ]: decrease/increase size of splats or spherelets\n");
		print("t: toggle transparent/solid mode (requires a files with splats)\n");
		print("b: switch between color specification (available only in spherelets mode)\n");
		print("s: toggle splats/spheres\n");
		print("m: output the rotation matrix\n");
		print("l: clear picked points. \n"); // added by jjcao
		print("middle button down: picking vertex or edge of the skeleton\n");// added by jjcao
		print("middle button down and move: drag picked vertex of the skeleton\n");// added by jjcao		
		print("Up function key + picked vertex: increase radius of the vertex\n");// added by jjcao
		print("Down function key + picked vertex: decrease radius of the vertex\n");// added by jjcao
		break;
	case 'b':
		show_backpointing = (show_backpointing+1)%3;
		break;
	case '[':
		markersize*=.95;
		ss << "markersize " << markersize;
		print(ss.str().c_str());
		break;
	case ']':
		markersize/=.95;
		ss << "markersize " << markersize;
		break;
	case 'l': //clear picked vertex and edge
		pick_vid = -1;
		break;
	}
	glutPostRedisplay();
}

// Matlab exits when the mex exits. I do not know why? Maybe it is the problem of glutMainLoop
void mexFunction(int nlhs, mxArray*plhs[], 
				 int nrhs, const mxArray*prhs[])
{
	/* Check for proper number of arguments. */
	if (nrhs < 1) {
		mexErrMsgTxt("At least one parameter required.");
	}
	mexWarnMsgTxt("begin \n");
	mexPrintf("begin \n");

	/* Handle parameters and outputs. */
	int npts, radius_dim(0);
	double *pts(0), *colors(0), *normals(0), *radius(0);

	npts = mxGetM(prhs[0]);
	pts = mxGetPr(prhs[0]);
	if (nrhs > 1) {
		colors = mxGetPr(prhs[1]);
	}
	if (nrhs > 2) {
		normals = mxGetPr(prhs[2]);
	}
	if (nrhs > 3) {
		radius = mxGetPr(prhs[3]);
		radius_dim = mxGetN(prhs[3]);
	}

	plhs[0] = mxCreateDoubleMatrix(npts, 1, mxREAL);
	pickedPtsIdx = mxGetPr(plhs[0]);

	/* build data structure. */
	pcloud = new Pcloud(npts, pts, colors, normals, radius, radius_dim); // added by jjcao

	/* initialize glut stuff and create window. */
	int argc(1); char* argv[1]; argv[0]="show_pts_gl.exe";
	glutInit( &argc, argv);
	glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH | GLUT_STENCIL );
	// viewport is set automatically
	glutInitWindowSize(window_width, window_height);
	int winid = glutCreateWindow("show_pts_gl");	

	// register callbacks
	glutDisplayFunc(displayCallback);
	glutKeyboardFunc(keyboardCallback);
	//glutSpecialFunc(specialKeyboardCallback);
	glutMotionFunc(motionCallback);
	glutMouseFunc(mouseCallback);
	glutReshapeFunc(reshapeCallback);

	// default drawing settings
	glClearColor(clear_color[0], clear_color[1], clear_color[2], 0.0);
	glClearStencil( 0x0 );

	glLineWidth(1.0);
	glEnable(GL_DEPTH_TEST); // delete covered objects

	// FIXED lights positions/setup
	glEnable(GL_LIGHTING);
	glLightModeli( GL_LIGHT_MODEL_TWO_SIDE, true );
	glEnable(GL_LIGHT0);
	glLightfv(GL_LIGHT0, GL_AMBIENT, light_ambient);
	glLightfv(GL_LIGHT0, GL_DIFFUSE, light_diffuse);
	glLightf(GL_LIGHT0, GL_POSITION, light0_position[0]);

	// setup blending options
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	// Corrects changes in shading during zoom
	glEnable(GL_NORMALIZE);

	// set the projection mode
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glOrtho(-1, 1, -1, 1, -10, 10);
	// glFrustum(-1, 1, -1, 1, -10, 10);
	gluLookAt( eyepos[0],eyepos[1],eyepos[2],0,0,0,0,1,0 );

	// Shape material
	// glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, red);

	glutMainLoop();
	glutDestroyWindow(winid);
	mexPrintf("pointlab terminated correctly.\n");
}
