#include <mex.h>


void nearcolor(double *sp,const double *R,const double *G,const double *B, const int colorNum,int row,int col ,int depth)
{
	double *dis = new double [colorNum]; 
	for(int i = 0; i<row*col; ++i)
	{
		int min = 255*255*255;
		int pos = 0;
		for(int j = 0; j<colorNum; ++j)
		{
			dis[j] = (sp[i]-R[j])*(sp[i]-R[j])+(sp[i+row*col]-G[j])*(sp[i+row*col]-G[j])+(sp[i+row*col*2]-B[j])*(sp[i+row*col*2]-B[j]);
			if(min>dis[j])
			{
				min = dis[j];
				pos = j;
			}

		}
		sp[i] = R[pos];
		sp[i+row*col] = G[pos];
		sp[i+row*col*2] = B[pos];
		
	}
}
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
	double *sp = mxGetPr(prhs[0]);
	double *R = mxGetPr(prhs[1]);
	double *G = mxGetPr(prhs[2]);
	double *B = mxGetPr(prhs[3]);
	double *colorNum = mxGetPr(prhs[4]);
	double *row = mxGetPr(prhs[5]);
	double *col = mxGetPr(prhs[6]);
	double *depth = mxGetPr(prhs[7]);
	 nearcolor( sp,R,G,B, (int) colorNum[0],(int) row[0],(int) col[0] ,(int) depth[0]);
	double sum = row[0]*col[0]*depth[0];
	 double *img;
	 plhs[0] = mxCreateDoubleMatrix(1,(int)sum,mxREAL);
	 img = mxGetPr(plhs[0]);
	 for (int i = 0;i<(int)sum;++i)
	 {
		 img[i] = sp[i];
	 }
}