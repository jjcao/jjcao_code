/************************************************************************
     File:        MyWindow.cpp

     Author:     
                  Steven J Chenney, schenney@cs.wisc.edu
     Modifier
                  Yu-Chi Lai, yu-chi@cs.wisc.edu
     
     Comment:     For window for displaying the OpenGL operation
	                (c) 2004 Steven J. Chenney
   
     Platform:    Visio Studio.Net 2003
*************************************************************************/
#include <iostream>
#include <math.h>

#include "MyWindow.h"

// Include the texture manager
#include "Texture.H"

//**************************************************************************
//
// * When window has computing resource, it will do this callback function
//   Normally use for animation and it is similar to the time_out callback function
//==========================================================================
void IdleCallback(void* pData)
//======================================================================
{
  if (pData != NULL)  {
		MyWindow* pWindow = reinterpret_cast<MyWindow*>(pData);
		if (pWindow->animating)		{
				pWindow->rotation += pWindow->rotationIncrement / 100;
				pWindow->redraw();
		}
  }
}




//************************************************************************
//
// * Constructor
//========================================================================
MyWindow::MyWindow(int width, int height, char* title) 
				: Fl_Gl_Window(width, height, title)
//========================================================================
{
	// Set up the GL window's properties
  mode(FL_RGB | FL_ALPHA | FL_DEPTH | FL_DOUBLE);

	// Set up the initial condition for rotation angle, increment and animation flag
  rotation          = 0.f;
  rotationIncrement = 10.f;
  animating         = false;

	// Add the IDle callback function
  Fl::add_idle(IdleCallback, this);
}


//************************************************************************
//
// * Destructor
//========================================================================
MyWindow::~MyWindow(void)
//========================================================================
{}


//************************************************************************
//
// * Initialize necessary GL operation
//========================================================================
void MyWindow::InitializeGL(void)
//========================================================================
{
	// Which color we should fill in the color buffer when glClear is called
  glClearColor(.1f, .1f, .1f, 1);
	// Enable the depth test of GL operation
  glEnable(GL_DEPTH_TEST);
	// Enable the lighting functionality
  glEnable(GL_LIGHTING);
	// Enable the light source 0
  glEnable(GL_LIGHT0);

	// Set up the diffuse light color
  float lightColor[4] = {1, 1, 1, 1};
  glLightfv(GL_LIGHT0, GL_DIFFUSE, lightColor);

	// Enable to set up color material
  glEnable(GL_COLOR_MATERIAL);
	// Set up the color material to be two sided and only use diffuse component
  glColorMaterial(GL_FRONT_AND_BACK, GL_DIFFUSE);

	// Enable the 2D textures
  glEnable(GL_TEXTURE_2D);
	// Set up the texture environment
  glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
}


//************************************************************************
//
// * Overwrite the window draw function
//========================================================================
void MyWindow::draw()
//========================================================================
{
 
  static bool firstTime = true; // Flag to control whether this is first time
  if (firstTime)   { // The first time, we need to set up the necessary GL parameters
    InitializeGL();
    firstTime = false;
  }// if

	// clear the color and depth buffer
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);      

	//***********************************************************************************
	//
	// * Set up the GL camera
	//
	//***********************************************************************************
  // Set up the projection transformation
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  glFrustum(-1, 1, -1, 1, 1, 100);

	// Set up the Model view transformation
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  gluLookAt(0, 0, 3, 0, 0, 0, 0, 1, 0);

	// Set up the light position
	float lightPosition[4] = {5, 5, 5, 1};
  glLightfv(GL_LIGHT0, GL_POSITION, lightPosition);

	// Do the rotation
	glRotatef(rotation, 0, 1, 0);

	// Enable the 2D textures
  glEnable(GL_TEXTURE_2D);
	// Set up the texture environment
  glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);

	// Initially load the image in
  fetchTexture("cs559", GL_REPEAT, GL_REPEAT);
  // draw something
  DrawCube();
}

//************************************************************************
//
// * Draw a cube which is center at (0, 0, 0)
//========================================================================
void MyWindow::DrawCube()
//========================================================================
{
  glBegin(GL_QUADS);
    // front
    glNormal3f(0, 0, 1);
    glColor3f(1, 0, 0);
    glTexCoord2f(0, 1);      glVertex3f(-1, 1, 1);
    glTexCoord2f(0, 0);      glVertex3f(-1, -1, 1);
    glTexCoord2f(1, 0);      glVertex3f(1, -1, 1);
    glTexCoord2f(1, 1);      glVertex3f(1, 1, 1);

    // back
    glNormal3f(0, 0, -1);
    glColor3f(0, 1, 0);
    glTexCoord2f(1, 1);      glVertex3f(-1, 1, -1);
    glTexCoord2f(0, 1);      glVertex3f(1, 1, -1);
    glTexCoord2f(0, 0);      glVertex3f(1, -1, -1);
    glTexCoord2f(1, 0);      glVertex3f(-1, -1, -1);

    // top
    glNormal3f(0, 1, 0);
    glColor3f(0, 0, 1);
    glTexCoord2f(0, 1);      glVertex3f(-1, 1, -1);
    glTexCoord2f(0, 0);      glVertex3f(-1, 1, 1);
    glTexCoord2f(1, 0);      glVertex3f(1, 1, 1);
    glTexCoord2f(1, 1);      glVertex3f(1, 1, -1);

    // bottom
    glNormal3f(0, -1, 0);
    glColor3f(1, 1, 0);
    glTexCoord2f(0, 0);      glVertex3f(-1, -1, -1);
    glTexCoord2f(1, 0);      glVertex3f(1, -1, -1);
    glTexCoord2f(1, 1);      glVertex3f(1, -1, 1);
    glTexCoord2f(0, 1);      glVertex3f(-1, -1, 1);

    // left
    glNormal3f(-1, 0, 0);
    glColor3f(0, 1, 1);
    glTexCoord2f(0, 1);      glVertex3f(-1, 1, -1);
    glTexCoord2f(0, 0);      glVertex3f(-1, -1, -1);
    glTexCoord2f(1, 0);      glVertex3f(-1, -1, 1);
    glTexCoord2f(1, 1);      glVertex3f(-1, 1, 1);

    // right
    glNormal3f(1, 0, 0);
    glColor3f(1, 0, 1);
    glTexCoord2f(0, 1);      glVertex3f(1, 1, 1);
    glTexCoord2f(0, 0);      glVertex3f(1, -1, 1);
    glTexCoord2f(1, 0);      glVertex3f(1, -1, -1);
    glTexCoord2f(1, 1);      glVertex3f(1, 1, -1);
  glEnd();
}


//**********************************************************************
//
// * Handle the event
//======================================================================
int MyWindow::handle(int event)
//======================================================================
{
  switch (event)  {
		case FL_FOCUS:  // When the window become selected
		case FL_UNFOCUS: // When the window become unselected
			return 1;

		case FL_KEYBOARD:  // When key has been pressed
			int key = Fl::event_key();
			switch (key)	{
				case FL_Left: // The left key
					rotation -= rotationIncrement;
					redraw();
					return 1;

				case FL_Right: // The right key
					rotation += rotationIncrement;
					redraw();
					return 1;

				case ' ': // The space key
					animating = !animating;
					return 1;
		  }
	}

  return Fl_Gl_Window::handle(event);
}

