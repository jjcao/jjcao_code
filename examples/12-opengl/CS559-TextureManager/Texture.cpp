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
#include "Texture.H"
#include "libtarga.h"

// the list of all textures loaded in already
Texture* theTextures = 0;

//**********************************************************************************
// The places to look for textures...
//    this should be modified to look where you might find your textures.
//    ultimately, there should be a better way than just sticking this in
//    the code.
//**********************************************************************************
char* texturePaths[] = {
  ".",
	"..",
	"Textures/signs",
  "../Textures/signs",
  "Textures/Textures",
  "../Textures/Textures",
  "Textures/Objects",
  "../Textures/Objects",
  "c:/src/GraphicsTown/Textures/signs",
  "c:/src/GraphicsTown/Textures/Textures",
  0 
};

//*************************************************************************************
//
// Utility routine
//
//*************************************************************************************

//*************************************************************************************
//
// * Copy and create a new string from the s
//=====================================================================================
static char* newString(char* s)
//=====================================================================================
{
  char* n = new char [strlen(s)+1];
  strcpy(n,s);
  return n;
}

// 
//*************************************************************************************
//
// * The core guts of the thing...
//=====================================================================================
GLuint fetchTexture(char* name, bool wrapS, bool wrapT)
//=====================================================================================
{
  char buf[200]; // String buffer for all possible filenames

  int   w, h;    // Width and height of the loading in image
  void* b=0;     // The buffer data for the image

  int i;

  //************************************************************
  // Check to see if we've loaded the texture already
  //************************************************************
  for(Texture* t=theTextures; t; t=t->next) {
		if (!strcmp(t->name,name)) {
			glBindTexture(GL_TEXTURE_2D,t->texName);
			return t->texName;
		}	  
  }

  //************************************************************
  // We need to load the file
  // Search everywhere on the path until we find it
  //************************************************************
  for(i=0; texturePaths[i] && !b; i++) {
		_snprintf(buf,200,"%s/%s.tga", texturePaths[i],name);
		printf("Trying: `%s'\n", buf);
		b = tga_load(buf, &w, &h, TGA_TRUECOLOR_32);
	}

  // If we found an image, load it into Texture data structure and put it into the
  // link list
  if (b) {
    // Create the structure
		Texture* t  = new Texture;
		// Set up the field in the data structure
		t->name     = newString(name);
		t->fname    = newString(buf);
		t->bits     = b;
		t->width    = w;
		t->height   = h;
		t->next     = theTextures;
		theTextures = t;
	
		// Generate a unique ID inside opengl
		glGenTextures(1,&t->texName);
		
		// Print out the creating message
		printf("Created texture %u for %s\n",t->texName,t->name);

		// Why bind twice, this need to be corrected
		glBindTexture(GL_TEXTURE_2D, t->texName);
	
		// Set up the texture clamping information in s direction
		glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_WRAP_S, 
						wrapS ? GL_REPEAT : GL_CLAMP);
		// Set up the texture clamping information in t direction
		glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_WRAP_T, 
						wrapT ? GL_REPEAT : GL_CLAMP);
		// Set up the alignment inside the texture
		glPixelStorei(GL_UNPACK_ALIGNMENT,1);

		// Generate the mip maps
		int e = gluBuild2DMipmaps(GL_TEXTURE_2D,GL_RGBA,w,h,GL_RGBA,GL_UNSIGNED_BYTE,b);

		printf("Build MipMap(%d %d) returns %d = %s (%d)\n",w,h,e,gluErrorString(e));
		return t->texName;
  } else {
		printf("!! Couldn't load texture `%s'~\n",name);
		return 0;
  }
}
