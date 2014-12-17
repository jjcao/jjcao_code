#ifndef _libtarga_h_
#define _libtarga_h_


/* uncomment this line if you're compiling on a big-endian machine */
/* #define WORDS_BIGENDIAN */


/* make sure these types reflect your system's type sizes. */
#define byte    char
#define int32   int
#define int16   short

#define ubyte   unsigned byte
#define uint32  unsigned int32
#define uint16  unsigned int16



/*  
    Truecolor images supported:

    bits            breakdown   components
    --------------------------------------
    32              8-8-8-8     RGBA
    24              8-8-8       RGB
    16              5-6-5       RGB
    15              5-5-5-1     RGB (ignore extra bit)


    Paletted images supported:
    
    index size      palette entry   breakdown   components
    ------------------------------------------------------
    8               <any of above>  <same as above> ..
    16              <any of above>  <same as above> ..
    24              <any of above>  <same as above> ..

*/



/*

   Targa files are read in and converted to
   any of these three for you -- you choose which you want.

   This is the 'format' argument to tga_create/load/write.
   
   For create and load, format is what you want the data
   converted to.

   For write, format is what format the data you're writing
   is in. (NOT the format you want written)

   Only TGA_TRUECOLOR_32 supports an alpha channel.

*/

#define TGA_TRUECOLOR_32      (4)
#define TGA_TRUECOLOR_24      (3)


/*
   Image data will start in the low-left corner
   of the image.
*/


#ifdef __cplusplus
extern "C" {
#endif


/* Error handling routines */
int             tga_get_last_error();
const char *    tga_error_string( int error_code );


/* Creating/Loading images  --  a return of NULL indicates a fatal error */
void * tga_create( int width, int height, unsigned int format );
void * tga_load( const char * file, int * width, int * height, unsigned int format );


/* Writing images to file  --  a return of 1 indicates success, 0 indicates error*/
int tga_write_raw( const char * file, int width, int height, unsigned char * dat, unsigned int format );
int tga_write_rle( const char * file, int width, int height, unsigned char * dat, unsigned int format );



#ifdef __cplusplus
}
#endif


#endif /* _libtarga_h_ */