// assign2.cpp : Defines the entry point for the console application.
//

/*
	CSCI 480 Computer Graphics
	Assignment 2: Simulating a Roller Coaster
	C++ starter code
*/

#include "stdafx.h"
#define cimg_use_jpeg /* see <your jpeg-9 directory>/jpeglib.h for using libjpeg 9 */
#include "CImg.h" 
using namespace cimg_library;
#include "Utils.h"
#include <windows.h>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <GL/glu.h>
#include <GL/glut.h>



GLuint texture[1];

int winWidth(640), winHeight(480);

/* represents one control point along the spline */
struct Point {
	double x;
	double y;
	double z;
};

/* spline struct which contains how many control points, and an array of control points */
struct Spline {
	int numControlPoints;
	struct Point *points;
};

/* the spline array */
struct Spline *g_Splines;

/* total number of splines */
int g_iNumOfSplines;


int loadSplines(char *argv) {
	char *cName = (char *)malloc(128 * sizeof(char));
	FILE *fileList;
	FILE *fileSpline;
	int iType, i = 0, j, iLength;

	/* load the track file */
	fileList = fopen(argv, "r");
	if (fileList == NULL) {
		printf ("can't open file\n");
		exit(1);
	}
  
	/* stores the number of splines in a global variable */
	fscanf(fileList, "%d", &g_iNumOfSplines);

	g_Splines = (struct Spline *)malloc(g_iNumOfSplines * sizeof(struct Spline));

	/* reads through the spline files */
	for (j = 0; j < g_iNumOfSplines; j++) {
		i = 0;
		fscanf(fileList, "%s", cName);
		fileSpline = fopen(cName, "r");

		if (fileSpline == NULL) {
			printf ("can't open file\n");
			exit(1);
		}

		/* gets length for spline file */
		fscanf(fileSpline, "%d %d", &iLength, &iType);

		/* allocate memory for all the points */
		g_Splines[j].points = (struct Point *)malloc(iLength * sizeof(struct Point));
		g_Splines[j].numControlPoints = iLength;

		/* saves the data to the struct */
		while (fscanf(fileSpline, "%lf %lf %lf", 
			&g_Splines[j].points[i].x, 
			&g_Splines[j].points[i].y, 
			&g_Splines[j].points[i].z) != EOF) {
			i++;
		}
	}

	free(cName);

	return 0;
}
void texload(int i, char *filename)
{
	CImg<unsigned char> img(filename);
	//img.permute_axes("cxyz"); // no time to find why it does not work.

	unsigned char *tex = new unsigned char[img.height()*img.width() * 4];
	int p = 0;
	for (int i = 0; i < img.height(); ++i) {
		for (int j = 0; j < img.width(); ++j) {
			tex[p++] = img(i, j, 0, 0);
			tex[p++] = img(i, j, 0, 1);
			tex[p++] = img(i, j, 0, 2);
			tex[p++] = 255;
		}
	}


	glBindTexture(GL_TEXTURE_2D, texture[i]);
	//glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, img.width(), img.height(), 0, GL_RGBA, GL_UNSIGNED_BYTE, img);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, img.width(), img.height(), 0, GL_RGBA, GL_UNSIGNED_BYTE, tex);
	delete[] tex;
}
void myinit()
{
	/* setup gl view here */
	glClearColor(0.0, 0.0, 0.0, 0.0);
	glEnable(GL_DEPTH_TEST);            // enable depth buffering
	glShadeModel(GL_SMOOTH);            // interpolate colors during rasterization
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

	glGenTextures(1, texture);
	texload(0, "test.jpg");//OhioPyle-128.jpg

}

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

// Catmall-Rom Spline
// points: 4 points from p[sid] to p[sid+3]
// sid: id of the start control points in points
// return: a vector of 1 * 3
Vector crSpline(float u, Point* points, int sid = 0, float s=0.5)
{
	CImg<float> U(4, 1, 1, 1); //Vector U(u*u*u, u*u, u, 1);
	int j = 0;
	U(0, j) = u*u*u; U(1,j) = u*u; U(2,j) = u; U(3,j) = 1;
	
	CImg<float> C(3, 4, 1, 1); //Matrix C; // control matrix;
	for (int i = 0; i < 4; ++i) 
	{
		C(0,i) = points[sid+i].x; C(1, i) = points[sid+i].y; C(2, i) = points[sid+i].z;
	}
	
	CImg<float> B(4, 4, 1, 1); //Matrix B; // basis matrix;
	j = 0;
	B(0,j) = -s; B(1, j) = 2-s; B(2, j) = s-2; B(3, j) = s;
	j = 1;
	B(0, j) = 2*s; B(1, j) = s-3; B(2, j) = 3 - 2*s; B(3, j) = -s;
	j = 2;
	B(0, j) = -s; B(1, j) = 0; B(2, j) = s; B(3, j) = 0;
	j = 3;
	B(0, j) = 0; B(1, j) = 1; B(2, j) = 0; B(3, j) = 0;
	
	CImg<float> val = U * B * C;
	//CImg<float> val1 = U * B;
	//CImg<float> val2 = B * C;
	Vector tmp(val(0, 0), val(1, 0), val(2, 0));
	return tmp;
}
// sid: id of the start control points in points
void subdivide(float u0, float u1, Point* points, int sid=0, float maxLen=0.01)
{
	float um = 0.5*(u0 + u1);
	Vector x0 = crSpline(u0, points, sid);
	Vector x1 = crSpline(u1, points, sid);
	if ((x0 - x1).Magnitude() > maxLen)
	{
		subdivide(u0, um, points, sid, maxLen);
		subdivide(um, u1, points, sid, maxLen);
	}
	else
	{
		glVertex3f(x0.x, x0.y, x0.z);
		glVertex3f(x1.x, x1.y, x1.z);//glVertex3f(x1(0, 0), x1(0, 1), x1(0, 2));
	}
}
void display()
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	// Reset transformations
	glLoadIdentity();

	gluLookAt(0.0, 0.0, 5.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0);
	

	///////////////////////////////////////////////////
	// draw some texture mapped objects
	///////////////////////////////////////////////////
	// turn on texture mapping (this disables standard OpenGL lighting, unless in GL_MODULATE mode)
	glEnable(GL_TEXTURE_2D);
	glBindTexture(GL_TEXTURE_2D, texture[0]);
	//glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
	glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	
	// back
	glBegin(GL_POLYGON);
	glTexCoord2f(1.0, 0.0);
	glVertex3f(1.0, -1.0, -1.0);
	glTexCoord2f(0.0, 0.0);
	glVertex3f(-1.0, -1.0, -1.0);
	glTexCoord2f(0.0, 1.0);
	glVertex3f(-1.0, 1.0, -1.0);
	glTexCoord2f(1.0, 1.0);
	glVertex3f(1.0, 1.0, -1.0);
	glEnd();

	// bottom
	glBegin(GL_POLYGON);
	glTexCoord2f(1.0, 0.0);
	glVertex3f(1.0, -1.0, 1.0);
	glTexCoord2f(0.0, 0.0);
	glVertex3f(-1.0, -1.0, 1.0);
	glTexCoord2f(0.0, 1.0);
	glVertex3f(-1.0, -1.0, -1.0);
	glTexCoord2f(1.0, 1.0);
	glVertex3f(1.0, -1.0, -1.0);
	glEnd();
	glDisable(GL_TEXTURE_2D);

	// turn off texture mapping
	glDisable(GL_TEXTURE_2D);


	///////////////////////////////////////////////////
	// draw some non-texture mapped objects
	///////////////////////////////////////////////////
	for (int i = 0; i < g_iNumOfSplines; ++i)
	{
		Spline sp = g_Splines[i];
		// draw control points
		glColor3f(0.0, 1.0, 0.0);
		//glEnable(GL_PROGRAM_POINT_SIZE);
		glPointSize(3.0);
		glBegin(GL_POINTS);	glBegin(GL_LINE_STRIP);
		for (int j = 0; j < sp.numControlPoints; ++j)
			glVertex3f(sp.points[j].x, sp.points[j].y, sp.points[j].z);
		glEnd();

		///////////////////
		glColor3f(1.0, 0.0, 0.0);
		glLineWidth(1.0);
		glBegin(GL_LINES);		
		for (int j = 0; j < sp.numControlPoints-3; ++j)
		{	
			subdivide(0, 1, sp.points, j);
		}
		glEnd();
	}

	glutSwapBuffers();
}
int _tmain(int argc, _TCHAR* argv[])
{
	// I've set the argv[1] to track.txt.
	// To change it, on the "Solution Explorer",
	// right click "assign1", choose "Properties",
	// go to "Configuration Properties", click "Debugging",
	// then type your track file name for the "Command Arguments"
	if (argc<2)
	{  
		printf ("usage: %s <trackfile>\n", argv[0]);
		exit(0);
	}

	loadSplines(argv[1]);

	/*
	create a window here..should be double buffered and use depth testing
	*/
	glutInit(&argc, (char**)argv);
	glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGBA);
	glutInitWindowPosition(100, 100);
	glutInitWindowSize(winWidth, winHeight);
	glutCreateWindow("Roller Coaster");
	glutReshapeFunc(reshapeCB);

	/* tells glut to use a particular display function to redraw */
	glutDisplayFunc(display);

	// jjcao
	/* replace with any animate code */
	//glutIdleFunc(doIdle);

	/* callback for mouse drags */
	//glutMotionFunc(mousedrag);
	/* callback for idle mouse movement */
	//glutPassiveMotionFunc(mouseidle);
	/* callback for mouse button changes */
	//glutMouseFunc(mousebutton);

	/* do initialization */
	myinit();

	glutMainLoop();
	return 0;
}