#include "GL/glut.h"
#include <math.h>
#include <stdio.h>
#include <string.h>

#define MAX_VERTEX 50000
#define MAX_SURFACE 50000

float gWidth;   // width of the render window
float gHeight;  // height of the render window

GLfloat gLightDiffuse[] = {0.0, 0.0, 1.0, 1.0};     // light diffuse
GLfloat gLightPosition[] = {0.0, 0.0, 1.0, 1.0};    // light position

GLfloat gDataVertex[MAX_VERTEX][3];     // vertex data: x, y, z
int gDataSurface[MAX_SURFACE][5];   // surface data: 3 or 4, first vertex, second vertex, third vertex, forth vertex
GLfloat gMinVertex[3];
GLfloat gMaxVertex[3];
GLfloat gCenterVertex[3];
char gTitle[8];
int gCountVertex;
int gCountSurface;
int gCount;

void init()
{
    char kind;
    int i, j, k;
    FILE* file = fopen("data.obj", "r");
    gCountVertex = 1;
    gCountSurface = 0;
    gCount = 0;
    gMinVertex[0] = gMinVertex[1] = gMinVertex[2] = 1.0;
    gMaxVertex[0] = gMaxVertex[1] = gMaxVertex[2] = 0.0;
    while(!feof(file))
    {
        fscanf(file, "%c", &kind);
        if(kind == 'v')
        {
            for(j = 0; j < 3; j++)
            {
                fscanf(file, "%f", &gDataVertex[gCountVertex][j]);
                if(gDataVertex[gCountVertex][j] < gMinVertex[j])
                {
                    gMinVertex[j] = gDataVertex[gCountVertex][j];
                }
                if(gDataVertex[gCountVertex][j] > gMaxVertex[j])
                {
                    gMaxVertex[j] = gDataVertex[gCountVertex][j];
                }
            }
            gCountVertex++;
        }
        else if(kind == 'f')
        {
            gDataSurface[gCountSurface][0] = 3;
            for(j = 0; j < 3; j++)
            {
                fscanf(file, "%d", &gDataSurface[gCountSurface][j + 1]);
            }
            gCountSurface++;
        }
    }
    gCenterVertex[0] = (gMinVertex[0] + gMaxVertex[0]) / 2;
    gCenterVertex[1] = (gMinVertex[1] + gMaxVertex[1]) / 2;
    gCenterVertex[2] = (gMinVertex[2] + gMaxVertex[2]) / 2;
    fclose(file);
}

void shape(int aFill)
{
    int i;
    GLfloat* v[3];
    GLfloat edge[3][3];
    GLfloat normal[3];
    GLfloat length;
    if(aFill)
    {
        glBegin(GL_TRIANGLES);
        for(i = 0; i < gCountSurface; i++)
        {
            v[0] = gDataVertex[gDataSurface[i][1]];
            v[1] = gDataVertex[gDataSurface[i][2]];
            v[2] = gDataVertex[gDataSurface[i][3]];
            // 3 edges
            edge[0][0] = v[1][0] - v[0][0];
            edge[0][1] = v[1][1] - v[0][1];
            edge[0][2] = v[1][2] - v[0][2];

            edge[1][0] = v[2][0] - v[1][0];
            edge[1][1] = v[2][1] - v[1][1];
            edge[1][2] = v[2][2] - v[1][2];
            
            edge[2][0] = v[0][0] - v[2][0];
            edge[2][1] = v[0][1] - v[2][1];
            edge[2][2] = v[0][2] - v[2][2];
            // normal
            normal[0] = edge[0][1] * edge[1][2] - edge[0][2] * edge[1][1];
            normal[1] = edge[0][2] * edge[1][0] - edge[0][0] * edge[1][2];
            normal[2] = edge[0][0] * edge[1][1] - edge[0][1] * edge[1][0];
            length = sqrt(normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2]);
            normal[0] /= length;
            normal[1] /= length;
            normal[2] /= length;
            glNormal3f(normal[0], normal[1], normal[2]);
            glVertex3f(v[0][0], v[0][1], v[0][2]);
            glVertex3f(v[1][0], v[1][1], v[1][2]);
            glVertex3f(v[2][0], v[2][1], v[2][2]);
        }
        glEnd();
    }
    else
    {
        for(i = 0; i < gCountSurface; i++)
        {
            v[0] = gDataVertex[gDataSurface[i][1]];
            v[1] = gDataVertex[gDataSurface[i][2]];
            v[2] = gDataVertex[gDataSurface[i][3]];
            glBegin(GL_LINE_LOOP);
            glVertex3f(v[0][0], v[0][1], v[0][2]);
            glVertex3f(v[1][0], v[1][1], v[1][2]);
            glVertex3f(v[2][0], v[2][1], v[2][2]);
            glEnd();
        }
    }
}

void reshape(int aWidth, int aHeight)
{
    // record the width and height of the current render window
    gWidth = aWidth;
    gHeight = aHeight;
}

void render(int aShape)
{
    switch(aShape)
    {
    case 0:
        {
            glTranslatef(-gCenterVertex[0], -gCenterVertex[1], -gCenterVertex[2]);
			 glRotatef(90, 1.0, 0.0, 0.0);
            glEnable(GL_LIGHTING);
            shape(0); // wire shape
            glDisable(GL_LIGHTING);
            glTranslatef(gCenterVertex[0], gCenterVertex[1], gCenterVertex[2]);
            break;
        }
    case 1:
        {
            glTranslatef(-gCenterVertex[0], -gCenterVertex[1], -gCenterVertex[2]);
            glEnable(GL_LIGHTING);
            shape(1); // solid shape
            glDisable(GL_LIGHTING);
            glTranslatef(gCenterVertex[0], gCenterVertex[1], gCenterVertex[2]);
            break;
        }
    }
}

void display(void)
{
    int row;
    int col;
    // clear the render window
    glViewport(0, 0, gWidth, gHeight);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    for(row = 0; row < 1; row++)
    {
        for(col = 0; col < 2; col++)
        {
            // render each geometric shape
            glViewport(gWidth / 2 * col, gHeight * row, gWidth / 2, gHeight);
            render(col);
        }
    }
    glFlush();
}

int main(int argc, char** argv)
{
    // initialize the render window
    glutInitWindowSize(800, 400);
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DEPTH | GLUT_SINGLE | GLUT_RGB);
    glutCreateWindow("GLUT NURBS");
    glutDisplayFunc(display);
    glutReshapeFunc(reshape);

    init();

    // set color
    glClearColor(0.5, 0.5, 0.5, 1.0);
    glColor3f(0.0, 0.0, 1.0);
    glLightfv(GL_LIGHT0, GL_DIFFUSE, gLightDiffuse);
    glLightfv(GL_LIGHT0, GL_POSITION, gLightPosition);
    glEnable(GL_LIGHT0);
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_AUTO_NORMAL);
    glEnable(GL_NORMALIZE);
    glShadeModel(GL_SMOOTH);

    // set OpenGL view model
    glMatrixMode(GL_PROJECTION);
    gluPerspective(
        22.0,   // field of view in degree
        1.0,    // aspect ratio
        1.0,    // Z near
        50.0);  // Z far
    glMatrixMode(GL_MODELVIEW);
    gluLookAt(
        0.0, 0.0, 35.0,  // eye is at (0,0,1.5)
        0.0, 0.0, 0.0,  // center is at (0,0,0)
        0.0, 1.0, 0.0); // up is in positive Y direction

    glutMainLoop();
    return 0;
}
