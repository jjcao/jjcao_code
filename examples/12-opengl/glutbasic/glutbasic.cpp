// the headers
#include <iostream>

#include <GL/glut.h>
#include <GL/gl.h>
#include <GL/glu.h>

#pragma comment(lib, "glut32.lib")
void testModelViewMatrix()
{
	
	std::cout << "test modelview matrix begin: " << std::endl;
	glMatrixMode(GL_MODELVIEW);
	int depth;
	glGetIntegerv(GL_MODELVIEW_STACK_DEPTH, &depth); 
	std::cout << "matrix stack depth: " << depth << std::endl;

	float mat[16];
	glGetFloatv(GL_MODELVIEW_MATRIX, mat);
	std::cout << "matrix elements: " << std::endl;
	for (int i = 0; i < 4; ++i)
		std::cout << mat[4 * i] << ", " << mat[4 * i + 1] << ", " << mat[4 * i + 2] << ", " << mat[4 * i + 3] << "; " << std::endl;

	std::cout << "test modelview matrix end: " << std::endl << std::endl;
}

// called before main loop
void init()
{
	glClearColor(0.0, 0.0, 0.0, 0.0);   // set background color
	glEnable(GL_DEPTH_TEST);            // enable depth buffering
	glShadeModel(GL_SMOOTH);            // interpolate colors during rasterization
}


// display a frame
void display()
{
	testModelViewMatrix();
	//gluLookAt(0.0, 0.0, 125.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0);
	//testModelViewMatrix();

	// clear buffers
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	// draw a triangle
	glBegin(GL_TRIANGLES);
	//glColor4f(1.0,0.0,0.0,1.0);
	glVertex3f(0.0, 0.0, -0.0);
	//glColor4f(0.0,1.0,0.0,1.0);
	glVertex3f(2.0, 0.0, -0.0);
	//glColor4f(0.0,0.0,1.0,1.0);
	glVertex3f(0.0, 2.0, -0.0);
	glEnd();

	glutSwapBuffers(); // double buffer flush
}

// called every time window is resized to update projection matrix
void reshape(int w, int h)
{
	// setup image size
	//glViewport(0, 0.5*h, 0.5*w, 0.5*h);
	glViewport(0, 0, w, h);
	glMatrixMode(GL_PROJECTION);

	// setup camera
	//glFrustum(-0.1, 0.1, -float(h) / (10.0*float(w)), float(h) / (10.0*float(w)), 0.5, 1000.0);
	gluOrtho2D(0.0, 3.0, 0.0, 3.0 * (GLfloat)h / (GLfloat)w);
	//gluLookAt(0, 0, 0, 0, 0, -12, 0, 1, 0);

	glMatrixMode(GL_MODELVIEW);
}

// entry point
int main(int argc, char **argv)
{

	// initialize GLUT
	glutInit(&argc, argv);

	// request double buffer
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_DEPTH | GLUT_RGBA);

	// set window size
	glutInitWindowSize(500, 500);

	// set window position
	glutInitWindowPosition(0, 0);

	// creates a window
	glutCreateWindow("Ahahaha!");

	// initialize states
	init();	

	// GLUT callbacks
	glutReshapeFunc(reshape);
	glutDisplayFunc(display);


	// start GLUT program
	glutMainLoop();
	return 0;
}


