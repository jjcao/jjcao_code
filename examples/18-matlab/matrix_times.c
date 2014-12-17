#include <mex.h>
void mexFunction(int nlhs, mxArray*plhs[], 
				 int nrhs, const mxArray*prhs[])
{
	double *A, *B;
	int m, n;
	/* check for arguments here, omitted */
	m= mxGetM(prhs[0]);
	n= mxGetN(prhs[0]);
    A = mxGetPr(prhs[0]);
	plhs[0] = mxCreateDoubleMatrix(m, n, mxREAL);	
	B = mxGetPr(plhs[0]);
    for(int i=0; i<m*n; ++i)
    {
        *(B+i)=*(A+i)*2.0;
    }
}
