// the headers

#include <GL/glut.h>
#include <GL/gl.h>
#include <GL/glu.h>

#pragma comment(lib, "glut32.lib")

GLfloat mv_matrix[16]; 

// called before main loop
void init()
{
	glClearColor(0.0, 0.0, 0.0, 0.0);   // set background color
	glEnable(GL_DEPTH_TEST);            // enable depth buffering
	glShadeModel(GL_SMOOTH);            // interpolate colors during rasterization
}

void drawAxis()
{
	glPointSize(5.0);
	glBegin(GL_LINES);
	glColor4f(1.0, 0.0, 0.0, 1.0);
	glVertex3f(0.0, 0.0, 0.0);
	glVertex3f(2.0, 0.0, 0.0);

	glColor4f(0.0, 1.0, 0.0, 1.0);
	glVertex3f(0.0, 0.0, 0.0);
	glVertex3f(0.0, 2.0, 0.0);

	glColor4f(0.0, 0.0, 1.0, 1.0);
	glVertex3f(0.0, 0.0, 0.0);
	glVertex3f(0.0, 0.0, 2.0);
	glEnd();
}
// display a frame
void display()
{
	// clear buffers
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	// move camera
	//gluLookAt(0, 0, 2, 0, 0, 0, 0, 1, 0);
	gluLookAt(0, 0, 2, 0, 0, 0, 1, 0, 0);
	glGetFloatv(GL_MODELVIEW_MATRIX, mv_matrix);

	//glTranslatef(0, 0, -2);
	drawAxis();
	glColor4f(1.0, 1.0, 1.0, 0.2);
	glutSolidSphere(0.5f, 20, 20);

	glutSwapBuffers(); // double buffer flush
}

// called every time window is resized to update projection matrix
void reshape(int w, int h)
{
	// setup image size
	glViewport(0, 0, w, h);

	// setup camera
	glMatrixMode(GL_PROJECTION);	
	glFrustum(-1, 1, -1, 1, 0.5, 100.0);
	glMatrixMode(GL_MODELVIEW);
}

// entry point
int main(int argc, char **argv)
{
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
