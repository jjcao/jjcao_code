/*=================================================================
 * The index of the sparse weight matrix with cot
 *
 * Usage: L is the Laplacian Matrix
   [I,J,Vr] = cotangentLaplacianMatrix(V,F);
	Q = sparse(I,J,Vr);
	W = 0.5*(Q + Q');
	diagValue = sum(W,2);
	D = sparse(1:length(diagValue),1:length(diagValue),diagValue);
	L = D - W;
 *
 *
 * =================================================================*/

#include <mex.h>
#include <math.h>

double cotangent(double P[], double Q[], double R[])
{
    double dot = (P[0] - Q[0]) * (R[0] - Q[0]) + (P[1] - Q[1]) * (R[1] - Q[1])
    + (P[2] - Q[2]) * (R[2] - Q[2]);
    
    double d[3];
    double cross_norm;
    
    d[0] = (Q[1] - P[1]) * (R[2] - P[2]) - (Q[2] - P[2]) * (R[1] - P[1]);
    d[1] = (Q[2] - P[2]) * (R[0] - P[0]) - (Q[0] - P[0]) * (R[2] - P[2]);
    d[2] = (Q[0] - P[0]) * (R[1] - P[1]) - (Q[1] - P[1]) * (R[0] - P[0]);
    cross_norm = sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
    
    if(cross_norm != 0.0)
        return (dot/cross_norm);
    else
        return 0.0;
}


void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray*prhs[])
{
    int i,rowv,rowf,colf,n1,n2,n3,k;
    double *V,*F,*Ir,*Jr,*Vr;
    double v1[3], v2[3], v3[3];
    
    rowv = mxGetM(prhs[0]);
    rowf = mxGetM(prhs[1]);
    colf = mxGetN(prhs[1]);
    if(colf != 3)
        mexErrMsgTxt("The mesh must be triangle mesh!");
    
    V = mxGetPr(prhs[0]);
    F = mxGetPr(prhs[1]);
    
    plhs[0] = mxCreateDoubleMatrix(3 * rowf,1,mxREAL);
    plhs[1] = mxCreateDoubleMatrix(3 * rowf,1,mxREAL);
    plhs[2] = mxCreateDoubleMatrix(3 * rowf,1,mxREAL);
    Ir = mxGetPr(plhs[0]);
    Jr = mxGetPr(plhs[1]);
    Vr = mxGetPr(plhs[2]);
    
    k = 0;
    for(i = 0;i < rowf;i++)
    {
        n1 = (int)F[i]-1;
        n2 = (int)F[rowf + i]-1;
        n3 = (int)F[2 * rowf + i]-1;
        v1[0] = V[n1];
        v1[1] = V[rowv + n1];
        v1[2] = V[2*rowv + n1];
        
        v2[0] = V[n2];
        v2[1] = V[rowv + n2];
        v2[2] = V[2*rowv + n2];
        
        v3[0] = V[n3];
        v3[1] = V[rowv + n3];
        v3[2] = V[2*rowv + n3];
        
        
        Ir[k] = n1+1;
        Jr[k] = n2+1;
        Vr[k] = cotangent(v1,v3,v2);
        k++;
        Ir[k] = n3+1;
        Jr[k] = n1+1;
        Vr[k] = cotangent(v3,v2,v1);
        k++;
        Ir[k] = n2+1;
        Jr[k] = n3+1;
        Vr[k] = cotangent(v2,v1,v3);
        k++;
    }
}

