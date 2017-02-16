// the headers

#include <GL/glut.h>
#include <GL/gl.h>
#include <GL/glu.h>

#pragma comment(lib, "glut32.lib")

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
	// clear buffers
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glLoadIdentity(); // reset transformation

					  // draw a triangle
	glBegin(GL_TRIANGLES);
	//glColor4f(1.0,0.0,0.0,1.0);
	glVertex3f(0.0, 0.0, -10.0);
	//glColor4f(0.0,1.0,0.0,1.0);
	glVertex3f(1.0, 0.0, -10.0);
	//glColor4f(0.0,0.0,1.0,1.0);
	glVertex3f(0.0, 1.0, -10.0);
	glEnd();

	glutSwapBuffers(); // double buffer flush
}

// called every time window is resized to update projection matrix
void reshape(int w, int h)
{
	// setup image size
	glViewport(0, 0, w, h);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();

	// setup camera
	glFrustum(-0.1, 0.1, -float(h) / (10.0*float(w)), float(h) / (10.0*float(w)), 0.5, 1000.0);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
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
