/************************************************************************
     File:        Texture.H

     Author:     
                  Michael L. Gleicher, gleicher@cs.wisc.edu
     Modifier
                  Yu-Chi Lai, yu-chi@cs.wisc.edu
     
     Comment:     Simple texture manager
	              (c) 2000 Michael L. Gleicher
                  This code may be modified and re-used providing that proper credit
                  is given to the author. 

				  Management of textures
                  probably more complicated than it needs to be, but it makes some
                  other things simple

                  The trick is to only have to load each texture once, no matter how
                  many things use it

     Compiler:    g++

     Platform:    Visio Studio.Net 2003
*************************************************************************/
#ifndef TEXTURE_H
#define TEXTURE_H

#include <Fl/Fl.h>
#include <Fl/Gl.h>
#include <Gl/Glu.h>
#include <stdio.h>
class Texture {
	public:
		// What we give to LibTarga
		char* name;                   //
		char* fname;					// the file name

		// What libtarga gives back...
		void* bits;
		int   width;
		int   height;

		GLuint texName;				// GL texture name

		// keep everything on a linked list
		Texture* next;
};

// The head of the texture link
extern Texture* theTextures;

// Get a texture- will hunt it down and load it if its not already loaded
// this even does the bind texture for you
// Note: the name should not include the ".tga"
// Note: the search path is (presently) hard coded into the .cpp
// file. this should be fixed someday.
GLuint fetchTexture(char* name, bool wrapS=false, bool wrapT=false);
#endif
