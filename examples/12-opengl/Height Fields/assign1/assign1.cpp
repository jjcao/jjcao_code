// assign1.cpp : Defines the entry point for the console application.
//

/*
  CSCI 480 Computer Graphics
  Assignment 1: Height Fields
  C++ starter code
*/

/*
 modified to CImg 2, libjpeg 9 and vc2015 by jjcao.github.io @ 2017
*/

#include "stdafx.h"
#define cimg_use_jpeg /* see <your jpeg-9 directory>/jpeglib.h for using libjpeg 9 */
#include "CImg.h" 
using namespace cimg_library;

#include <windows.h>
#include <stdlib.h>
#include <GL/glu.h>
#include <GL/glut.h>

CImg<unsigned char> heightData;
int winWidth(640), winHeight(480);

int g_iMenuId;

int g_vMousePos[2] = {0, 0};
int g_iLeftMouseButton = 0;    /* 1 if pressed, 0 if not */
int g_iMiddleMouseButton = 0;
int g_iRightMouseButton = 0;

typedef enum { ROTATE, TRANSLATE, SCALE } CONTROLSTATE;

CONTROLSTATE g_ControlState = ROTATE;

/* state of the world */
float g_vLandRotate[3] = {0.0, 0.0, 0.0};
float g_vLandTranslate[3] = {0.0, 0.0, 0.0};
float g_vLandScale[3] = {1.0, 1.0, 1.0};

void convertToNonInterleaved(int w, int h, unsigned char* tangled, unsigned char* untangled) {
	//Take string in format R1 G1 B1 R2 G2 B2... and re-write it 
	//in the format format R1 R2 G1 G2 B1 B2... 
	//Assume 8 bit values for red, green and blue color channels.
	//Assume there are no other channels
	//tangled is a pointer to the input string and untangled 
	//is a pointer to the output string. This method assumes that 
	//memory has already been allocated for the output string.

	int numPixels = w*h;
	int numColors = 3;
	for (int i = 0; i<numPixels; ++i) {
		int indexIntoInterleavedTuple = numColors*i;
		//Red
		untangled[i] = tangled[indexIntoInterleavedTuple];
		//Green
		untangled[numPixels + i] = tangled[indexIntoInterleavedTuple + 1];
		//Blue
		untangled[2 * numPixels + i] = tangled[indexIntoInterleavedTuple + 2];
	}
}

/* Write a screenshot to the specified filename */
void saveScreenshot (char *filename)
{
  int bytes = winWidth * winHeight * 3; //Color space is RGB
  unsigned char* buffer;
  buffer = (GLubyte *)malloc(bytes);
  glReadPixels(0, 0, winWidth, winHeight, GL_RGB, GL_UNSIGNED_BYTE, buffer);

  unsigned char* p = (unsigned char *)malloc(bytes);
  convertToNonInterleaved(winWidth, winHeight, buffer, p);

  CImg<unsigned char> img(p, winWidth, winHeight, 1, 3);
  free(p);
  printf("File to save to: %s\n", filename);
  img.mirror('y');
  //CImgDisplay main_disp(img, "Snapshot");
  //while (!main_disp.is_closed()) {
	 // main_disp.wait();
  //}

  img.save(filename);
}

void myinit()
{
  /* setup gl view here */
  glClearColor(0.0, 0.0, 0.0, 0.0);
  glEnable(GL_DEPTH_TEST);            // enable depth buffering
  glShadeModel(GL_SMOOTH);            // interpolate colors during rasterization
  
}


///////////////////////////////////////////////////////////////////////////////
// set the projection matrix as perspective
///////////////////////////////////////////////////////////////////////////////
void toPerspective()
{
	// set viewport to be the entire window
	glViewport(0, 0, (GLsizei)winWidth, (GLsizei)winHeight);

	//// set perspective viewing frustum
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();	
	float aspect = (float)(winWidth) / winHeight;
	float fovy(60.0f);
	gluPerspective(fovy, aspect, 1.0f, 1000.0f); // FOV, AspectRatio, NearClip, FarClip

	// switch to modelview matrix in order to set scene
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	
}
void reshapeCB(int w, int h)
{
	winWidth = w;
	winHeight = h;
	toPerspective();
	
}

void display()
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	// Reset transformations
	glLoadIdentity();

	gluLookAt(0.0, 0.0, 125.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0);
	glTranslatef(g_vLandTranslate[0], g_vLandTranslate[1], g_vLandTranslate[2]);
	glRotatef(g_vLandRotate[0], 1.0f, 0.0f, 0.0f);
	glRotatef(g_vLandRotate[1], 0.0f, 1.0f, 0.0f);
	glRotatef(g_vLandRotate[2], 0.0f, 0.0f, 1.0f);
	glRotatef(1.0f, g_vLandRotate[0], g_vLandRotate[1], g_vLandRotate[2]);
	glScalef(g_vLandScale[0], g_vLandScale[1], g_vLandScale[2]);

  /* draw 1x1 cube about origin */
  /* replace this code with your height field implementation */
  /* you may also want to precede it with your rotation/translation/scaling */
	float hh = heightData.height()*0.5;
	float hw = heightData.width()*0.5;
  glBegin(GL_POLYGON);
  glColor3f(1.0, 0.0, 0.0);
  glVertex3f(-hw, -hh, 0.0);
  glColor3f(1.0, 0.0, 0.0);
  glVertex3f(-hw, hh, 0.0);

  glColor3f(0.0, 1.0, 0.0);
  glVertex3f(hw, hh, 0.0);
  glColor3f(1.0, 0.0, 0.0);
  glVertex3f(hw, -hh, 0.0);
  glEnd();


	
	CImg<float> z0, z1;
	float x, y0, y1;
	int step(1);
	float hScale(25.0f);
	//glFrontFace(GL_CW);
  for (int i = 0; i < heightData.height() - step; i = i+step) {
	  y0 = i - hh;
	  y1 = i + step - hh;
	  glBegin(GL_TRIANGLE_STRIP);
	  for (int j = 0; j < heightData.width(); j = j + step) {
		  x = j - hw;
		  z0 = heightData.get_vector_at(j, i)/255.0 ; // 'top' vertex
		  z1 = heightData.get_vector_at(j, i+step) / 255.0;// 'bottom' vertex

		  // sequential top,bottom vert pairs generates a tri-strip
		  //glColor3f(z0(0), z0(1), z0(2));
		  glColor3f(0.0f, 0.0f, z0(0));
		  glVertex3f(x, y0, z0(0)* hScale);
		  //glVertex3f(x, y0, 1.0);
		  glColor3f(0.0f, 0.0f, z1(0));
		  glVertex3f(x, y1, z1(0)* hScale);
		  //glVertex3f(x, y1, 1.0);
		}// next pixel in current row
	  glEnd();
}// next row

  glutSwapBuffers();
}

void menufunc(int value)
{
  switch (value)
  {
    case 0:
      exit(0);
      break;
	case 1:
		// render vertices
		glPolygonMode(GL_FRONT_AND_BACK, GL_POINT);
		break;
	case 2:
		// render wireframe
		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
		break;
	case 3:
		// render solid triangles
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
		break;
	case 4:
		reshapeCB(winWidth, winHeight);

		g_vLandTranslate[0] = 0.0f, g_vLandTranslate[1] = 0.0f, g_vLandTranslate[2] = 0.0f;
		g_vLandRotate[0] = 0.0f, g_vLandRotate[1] = 0.0f, g_vLandRotate[2] = 0.0f;
		g_vLandScale[0] = 1.0f, g_vLandScale[1] = 1.0f, g_vLandScale[2] = 1.0f;
		break;
  }
}

void doIdle()
{
  /* do some stuff... */

  /* make the screen update */
  glutPostRedisplay();
}

/* converts mouse drags into information about 
rotation/translation/scaling */
void mousedrag(int x, int y)
{
  int vMouseDelta[2] = {x-g_vMousePos[0], y-g_vMousePos[1]};
  
  switch (g_ControlState)
  {
    case TRANSLATE:  
      if (g_iLeftMouseButton)
      {
        g_vLandTranslate[0] += vMouseDelta[0]*0.1;
        g_vLandTranslate[1] += vMouseDelta[1]*0.1;
      }
      if (g_iMiddleMouseButton)
      {
        g_vLandTranslate[2] += vMouseDelta[1]*0.1;
      }
	  
      break;
    case ROTATE:
      if (g_iLeftMouseButton)
      {
        g_vLandRotate[0] += vMouseDelta[1];
        g_vLandRotate[1] += vMouseDelta[0];
      }
      if (g_iMiddleMouseButton)
      {
        g_vLandRotate[2] += vMouseDelta[1];
      }
	  
      break;
    case SCALE:
      if (g_iLeftMouseButton)
      {
        g_vLandScale[0] *= 1.0+vMouseDelta[0]*0.001;
        g_vLandScale[1] *= 1.0+vMouseDelta[1]*0.001;
      }
      if (g_iMiddleMouseButton)
      {
        g_vLandScale[2] *= 1.0+vMouseDelta[1]*0.01;
      }
	  
      break;
  }
  g_vMousePos[0] = x;
  g_vMousePos[1] = y;
}

void mouseidle(int x, int y)
{
  g_vMousePos[0] = x;
  g_vMousePos[1] = y;
}

void mousebutton(int button, int state, int x, int y)
{

  switch (button)
  {
    case GLUT_LEFT_BUTTON:
      g_iLeftMouseButton = (state==GLUT_DOWN);
      break;
    case GLUT_MIDDLE_BUTTON:
      g_iMiddleMouseButton = (state==GLUT_DOWN);
      break;
    case GLUT_RIGHT_BUTTON:
      g_iRightMouseButton = (state==GLUT_DOWN);
      break;
  }
 
  switch(glutGetModifiers())
  {
    case GLUT_ACTIVE_CTRL:
      g_ControlState = TRANSLATE;
      break;
    case GLUT_ACTIVE_SHIFT:
      g_ControlState = SCALE;
      break;
    default:
      g_ControlState = ROTATE;
      break;
  }

  g_vMousePos[0] = x;
  g_vMousePos[1] = y;
}

int main(int argc, char* argv[])
{
	// I've set the argv[1] to spiral.jpg.
	// To change it, on the "Solution Explorer",
	// right click "assign1", choose "Properties",
	// go to "Configuration Properties", click "Debugging",
	// then type your texture name for the "Command Arguments"
	if (argc<2)
	{  
		printf ("usage: %s heightfield.jpg\n", argv[0]);
		exit(1);
	}

	heightData.load((char*)argv[1]);
	//heightData.display();




	/*
	create a window here..should be double buffered and use depth testing
	*/
	glutInit(&argc,(char**)argv);	
	glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGBA);
	glutInitWindowPosition(100, 100);
	glutInitWindowSize(winWidth, winHeight);	
	glutCreateWindow("Height Field");
	glutReshapeFunc(reshapeCB);

	/* tells glut to use a particular display function to redraw */
	glutDisplayFunc(display);
  
	// jjcao
	/* allow the user to quit using the right mouse button menu */
	g_iMenuId = glutCreateMenu(menufunc);
	glutSetMenu(g_iMenuId);
	glutAddMenuEntry("Quit",0);
	glutAddMenuEntry("Point", 1);
	glutAddMenuEntry("Wireframe", 2);
	glutAddMenuEntry("Solid", 3);
	glutAddMenuEntry("Restore", 4);
	glutAttachMenu(GLUT_RIGHT_BUTTON);
  
	/* replace with any animate code */
	glutIdleFunc(doIdle);

	/* callback for mouse drags */
	glutMotionFunc(mousedrag);
	/* callback for idle mouse movement */
	glutPassiveMotionFunc(mouseidle);
	/* callback for mouse button changes */
	glutMouseFunc(mousebutton);

	/* do initialization */
	myinit();

	glutMainLoop();
	return 0;
}