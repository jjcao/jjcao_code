//mex mex_get_depth.cpp -lOpenGL32

#include <mex.h>
#include "windows.h"
//#include "glut.h"
#include <gl/gl.h>
#include <gl/glu.h>

#define GL_VIEWPORT                       0x0BA2
#define GL_DEPTH_COMPONENT                0x1902
#define GL_FLOAT                          0x1406

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    int viewport[4], i, x, y;
    int colLen;
    float *data;
    double *matrix;

    glGetIntegerv(GL_VIEWPORT, viewport);
    data = (float*)malloc(viewport[2]*viewport[3]*sizeof(float));
    glReadPixels(0, 0, viewport[2], viewport[3], GL_DEPTH_COMPONENT, GL_FLOAT, data);

    plhs[0] = mxCreateNumericMatrix(viewport[3], viewport[2], mxDOUBLE_CLASS, mxREAL);
    matrix = mxGetPr(plhs[0]);
    colLen = mxGetM(plhs[0]);

    for(x = 0; x < viewport[2]; ++ x) {
        for(y = 0; y < viewport[3]; ++ y)
           matrix[x * colLen + y] = data[x + (viewport[3] - 1 - y) * viewport[2]];
    }

    free(data);
    return;
}