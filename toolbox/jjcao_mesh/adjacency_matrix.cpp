/*=================================================================
*
*
* usage: 
          [I J] = adjacency_matrix(faces);
          A = sparse(I,J,ones(size(I)));
		  A = double(A>0);
		  A = max(A, A');
* actually use the following matlab code directly, maybe faster!
		  A = sparse([f(:,1); f(:,1); f(:,2); f(:,2); f(:,3); f(:,3)], ...
                   [f(:,2); f(:,3); f(:,1); f(:,3); f(:,1); f(:,2)], ...
                   1.0);
        % avoid double links
        A = double(A>0);
*
* Adapted by JJCAO, 2012
* Hui Wang, wanghui19841109@gmail.com, Dec. 11, 2011
*
*=================================================================*/

#include <mex.h>

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray*prhs[])
{
	// input
	int rowf = mxGetM(prhs[0]);
	int colf = mxGetN(prhs[0]);
	if(colf != 3)
		mexErrMsgTxt("The mesh must be triangle mesh!");

	double *faces = mxGetPr(prhs[0]);

	// output
	plhs[0] = mxCreateDoubleMatrix(3 * rowf,1,mxREAL);// the size is guessed
	plhs[1] = mxCreateDoubleMatrix(3 * rowf,1,mxREAL);
	double *Ir = mxGetPr(plhs[0]);
	double *Jr = mxGetPr(plhs[1]);

	// compute
	int k = 0;
	int v1,v2,v3;
	for(int i = 0; i < rowf; ++i)
	{
		v1 = (int)faces[i];
		v2 = (int)faces[rowf + i];
		v3 = (int)faces[2 * rowf + i];
        
        Ir[k] = v1;
		Jr[k] = v2;
		++k;

		Ir[k] = v2;
		Jr[k] = v3;
		++k;

		Ir[k] = v3;
		Jr[k] = v1;
		++k;
	}	
}