#include <GL/glut.h>
#include <GL/gl.h>
#include <GL/glu.h>

#pragma comment(lib, "glut32.lib")
static int width;
static int height;

static void displayTri()
{
	// draw a triangle
	glBegin(GL_TRIANGLES);
	//glColor4f(1.0,0.0,0.0,1.0);
	glVertex3f(0.0, 0.0, -10.0);
	//glColor4f(0.0,1.0,0.0,1.0);
	glVertex3f(2.0, 0.0, -10.0);
	//glColor4f(0.0,0.0,1.0,1.0);
	glVertex3f(0.0, 2.0, -10.0);
	glEnd();
}
static void display(void) {
	glClear(GL_COLOR_BUFFER_BIT);
	glColor3f(1.0f, 0.0f, 0.0f);

	glViewport(0, 0, width / 2, height / 2);
	glLoadIdentity();
	//gluLookAt(0.0, 0.0, -3.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0);
	displayTri();

	glViewport(width / 2, 0, width / 2, height / 2);
	glLoadIdentity();
	//gluLookAt(0.0, 0.0, 3.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0);
	displayTri();

	glViewport(0, height / 2, width / 2, height / 2);
	glLoadIdentity();
	gluLookAt(0.0, 3.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0);
	glutWireTeapot(1);

	glViewport(width / 2, height / 2, width / 2, height / 2);
	glLoadIdentity();
	gluLookAt(0.0, -3.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0);
	glutWireTeapot(1);

	glFlush();
}

static void reshape(int w, int h) {
	width = w;
	height = h;
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glFrustum(-1.0, 1.0, -1.0, 1.0, 1.5, 20.0);
	glMatrixMode(GL_MODELVIEW);
}

int main(int argc, char** argv) {
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB);
	glutInitWindowSize(500, 500);
	glutInitWindowPosition(100, 100);
	glutCreateWindow(argv[0]);
	glClearColor(0.0, 0.0, 0.0, 0.0);
	glShadeModel(GL_FLAT);
	glutDisplayFunc(display);
	glutReshapeFunc(reshape);
	glutMainLoop();
	return 0;
}