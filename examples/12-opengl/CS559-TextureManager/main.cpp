/************************************************************************
     File:        main.cpp

     Author:     
                  Steven J Chenney, schenney@cs.wisc.edu
     Modifier
                  Yu-Chi Lai, yu-chi@cs.wisc.edu
     
     Comment:     For demoing the texture in OpenGL
	                (c) 2004 Steven J. Chenney
   
     Platform:    Visio Studio.Net 2003
*************************************************************************/
#include <Fl/Fl.h>
#include "MyWindow.h"

int main(int argc, char** args)
{
	// Create the window
   MyWindow myWindow(400, 400, "CS559 Tutorial");
	 // Show the window
   myWindow.show();
   // Run it
   Fl::run();

   return 0;
}
