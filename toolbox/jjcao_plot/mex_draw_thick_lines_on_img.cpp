/*
//
// Draw multiple thick lines on the image
// 
// Inputs: 
//      InputImg: Input image (Grayscale or Color)
//      CoordPnt: Coordinates of end points of lines [r1 r2 c1 c2] 
//          (n x 4 double array, n: # of lines)
//      Thickness: Thickness of lines (The value is integer, but the type is double.)
//          (n x 1 double array, n: # of lines)
//      LineColor: Line colors (The values should be integers from 0 to 255)
//          (n x 4 double array, n: # of lines)
//          (Data type: double, ex. [255 255 255]: white )
//          (The channels of 'InputImg' and 'LineColor' should be consistent)
// 
//
// Outputs: 
//      OutputImg: Output image (The same format with InputImg) 
// 
//
// Copyright, Gunhee Kim (gunhee@cs.cmu.edu)
// Computer Science Department, Carnegie Mellon University, 
// October 19 2009
*/

#include <mex.h>
#include <math.h>
//#include <string.h>

#define UINT8 unsigned char 
#define UINT32 unsigned int 


#define 	max(a, b)   ((a) > (b) ? (a) : (b))
#define 	min(a, b)   ((a) < (b) ? (a) : (b))
#define 	abs(x)      ((x) >= 0 ? (x) : -(x))
#define 	pow2(x)     ((x)*(x))


void Bresenham(int x1, int x2, int y1, int yt)  ;
void plot(int x1, int y1)  ;
void draw_thick_line(int x1, int x2, int y1, int y2) ;

double *LineColor, *LineThickness ;
UINT8 *OutputImg ;
int nPnt, nChannel, indLine ;
int dimR, dimC, nDim ;

void mexFunction(
				 int nlhs,              // Number of left hand side (output) arguments
				 mxArray *plhs[],       // Array of left hand side arguments
				 int nrhs,              // Number of right hand side (input) arguments
				 const mxArray *prhs[]  // Array of right hand side arguments
				 ) {

	UINT8 *InputImg ;
    double *CoordPnt, *InCoord ;
    int i, dimPnt, add_r1, add_r2, add_c1, add_c2 ;
    
    /* Check for proper number of arguments. */
    if (nrhs <4)
    { 
        mexErrMsgTxt("Four inputs are required.");
    } 
    else if (nlhs > 1) 
    {
        mexErrMsgTxt("Too many output assigned."); 
    } 
            
    InputImg = (UINT8 *)mxGetPr(prhs[0]);       // Input image (Data type: uint8*)            
    
    InCoord = (double *) mxGetPr(prhs[1]);     // [r1 r2 c1 c2] (Data type: double*)
    nPnt = (int) mxGetM(prhs[1]);
    dimPnt = (int) mxGetN(prhs[1]);
    
    LineThickness = (double *)mxGetPr(prhs[2]);     // Line Thickness (Data type: double*)
    LineColor = (double *)mxGetPr(prhs[3]);         // Line color (Data type: double*)
    nChannel = (int)mxGetN(prhs[3]);
    
    // Dimensions of input image 
    dimR = (int)mxGetM(prhs[0]);
    dimC = (int)mxGetN(prhs[0]);
    nDim = (int)mxGetNumberOfDimensions(prhs[0]);

	/*    
    // If the channels of Input image and line color are different, terminate it. 
    if (((nDim==2) && (nChannel==3)) || ((nDim==3) && (nChannel==1)) ) {
        mexErrMsgTxt("The channels of Input image and line color are different !");
    }
    */
    // If the input image is color image
    if(nChannel==3)
        dimC = dimC/nChannel;
    
    //mexPrintf("%d %d %d %d %d %d\n", nChannel, dimR, dimC, nDim, dimPnt, nPnt) ;
    
    // Set output
    plhs[0] = mxCreateNumericArray(mxGetNumberOfDimensions(prhs[0]), 
            mxGetDimensions(prhs[0]),mxUINT8_CLASS,mxREAL);
    OutputImg = (UINT8*)mxGetPr(plhs[0]);

    // copy InputImg to OutputImg
    for (i=0; i<dimR*dimC*nChannel; i++) 
        OutputImg[i] = InputImg[i] -1 ;

	
    // decrease coordinates by 1
    CoordPnt = (double *) mxCalloc(nPnt*dimPnt, sizeof(double));
	// If the coordinate< 0 or > image dimension, reduce it.
    for (i=0; i<nPnt*dimPnt/2; i++)
	{
		CoordPnt[i] = InCoord[i] -1 ;
		if (CoordPnt[i]<0) CoordPnt[i] = 0 ; 
		if (CoordPnt[i]>dimR) CoordPnt[i] = dimR ; 
	}
    for (i=nPnt*dimPnt/2; i<nPnt*dimPnt; i++)
	{
		CoordPnt[i] = InCoord[i] -1 ;
		if (CoordPnt[i]<0) CoordPnt[i] = 0 ; 
		if (CoordPnt[i]>dimC) CoordPnt[i] = dimC ; 
	}
    
    
    // adding indices for each dim. 
    add_r1 = 0 ; 
    add_r2 = nPnt ;
    add_c1 = 2*nPnt ;
    add_c2 = 3*nPnt ;

    // Main loop
    for (i=0; i<nPnt; i++) {
        // lineIdx
        indLine = i ;
        
        draw_thick_line((int)CoordPnt[i+add_c1],    // c1 = x1
            (int)CoordPnt[i+add_c2],                // c2 = x2
            (int)CoordPnt[i+add_r1],                // r1 = y1
            (int)CoordPnt[i+add_r2]                 // r2 = y2
        ) ;
    }

    mxFree(CoordPnt) ;
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////
// 	End of the function "Main" function
//////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Set the values
void plot(int x1, int y1) 
{
    int idx = y1 + dimR*x1 ;
    for (int i=0; i<nChannel; i++)
        OutputImg[dimR*dimC*i+idx] = (UINT8)LineColor[indLine+i*nPnt] ;
}

// draw thick lines
// (reference) Computer Graphics by A.P. Godse
void draw_thick_line(int x1, int x2, int y1, int y2) 
{
    int i;
    float wy, wx ;
    int thickness = (int)LineThickness[indLine] ;
    Bresenham(x1, x2, y1, y2) ;
    
    if (abs((float)(y2-y1)/(float)(x2-x1)) < 1) {
        wy = (float)(thickness-1)*sqrt((float)(pow2(x2-x1)+pow2(y2-y1))) / (2.*(float)abs(x2-x1)) ;
        //mexPrintf("wy: %f %d %f \n", 2.*(float)abs(x2-x1), thickness, wy) ;
        for (i=0; i<wy; i++) {
            if ( ((y1-i)>-1) && ((y2-i)>-1) )
                Bresenham(x1, x2, y1-i, y2-i) ;
            if ( ((y1+i)<dimR) && ((y2+i)<dimR) )
                Bresenham(x1, x2, y1+i, y2+i) ;
        }
    }
    else {
        wx = (float)(thickness-1)*sqrt((float)(pow2(x2-x1)+pow2(y2-y1))) / (2.*(float)abs(y2-y1)) ;
        //mexPrintf("wx: %f %d %f \n", 2.*(float)abs(y2-y1), thickness, wx) ;
        for (i=0; i<wx; i++) {
            if ( ((x1-i)>-1) && ((x2-i)>-1) )
                Bresenham(x1-i, x2-i, y1, y2) ;
            if ( ((x1+i)<dimC) && ((x2+i)<dimC) )
                Bresenham(x1+i, x2+i, y1, y2) ;
        }
    }
    
}

// The following code is Bresenham's Line Algorithm.
// (Reference) http://roguebasin.roguelikedevelopment.org/index.php?title=Bresenham's_Line_Algorithm
////////////////////////////////////////////////////////////////////////////////
void Bresenham(int x1,
    int x2,
    int y1,
    int y2)
{
    int delta_x = abs(x2 - x1) << 1;
    int delta_y = abs(y2 - y1) << 1;

    // if x1 == x2 or y1 == y2, then it does not matter what we set here
    signed char ix = x2 > x1?1:-1;
    signed char iy = y2 > y1?1:-1;

    plot(x1, y1);

    if (delta_x >= delta_y)
    {
        // error may go below zero
        int error = delta_y - (delta_x >> 1);

        while (x1 != x2)
        {
            if (error >= 0)
            {
                if (error || (ix > 0))
                {
                    y1 += iy;
                    error -= delta_x;
                }
                // else do nothing
            }
            // else do nothing

            x1 += ix;
            error += delta_y;

            plot(x1, y1);
        }
    }
    else
    {
        // error may go below zero
        int error = delta_x - (delta_y >> 1);

        while (y1 != y2)
        {
            if (error >= 0)
            {
                if (error || (iy > 0))
                {
                    x1 += ix;
                    error -= delta_y;
                }
                // else do nothing
            }
            // else do nothing

            y1 += iy;
            error += delta_x;

            plot(x1, y1);
        }
    }
}
