/************************************************************************
     File:        MyWindow.h

     Author:     
                  Steven J Chenney, schenney@cs.wisc.edu
     Modifier
                  Yu-Chi Lai, yu-chi@cs.wisc.edu
     
     Comment:     For window for displaying the OpenGL operation
	                (c) 2004 Steven J. Chenney
   
     Platform:    Visio Studio.Net 2003
*************************************************************************/
#ifndef MY_WINDOW_H
#define MY_WINDOW_H

#include <Fl/Fl_Gl_Window.h>

class MyWindow : public Fl_Gl_Window
{
	public:
		// Constructor and destructor
    MyWindow(int width, int height, char* title);
    virtual ~MyWindow(void);
    void InitializeGL(void);
    virtual void draw(void);
    void DrawCube(void);
    virtual int handle(int event);

    float          rotation;
    float          rotationIncrement;
    bool           animating;
    unsigned int   textureId;
};

#endif
