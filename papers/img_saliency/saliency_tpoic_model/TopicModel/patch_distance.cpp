#include <yvals.h>
#if (_MSC_VER >= 1600)
#define __STDC_UTF_16__
#endif
#include<mex.h>
#include<iostream>
#include<mex.h>
#include<iostream>
#include<list>
#include<vector>
#include<cmath>
using namespace std;
#define dim 3

typedef struct N
{
	double  num;
	double  colorL;
	double  colorA;
	double  colorB;

bool operator < (N &N)
{
	return ( num >= N.num );
}
bool operator > (N &N)
{
	return ( num <= N.num );
}
} N; 

double distanceIJ(vector<N> colorArrayI,vector<N> colorArrayJ)
{
	float distance = 0;
	for(int i=0;i<colorArrayI.size();++i)
	{
		for(int j=0;j<colorArrayJ.size();++j)
		{
			double deltaL = colorArrayI[i].colorL-colorArrayJ[j].colorL;
			double deltaA = colorArrayI[i].colorA-colorArrayJ[j].colorA;
			double deltaB = colorArrayI[i].colorB-colorArrayJ[j].colorB;
			double tmpdistance = deltaL*deltaL+deltaA*deltaA+deltaB*deltaB;
			distance += sqrt(tmpdistance)*colorArrayI[i].num*colorArrayJ[j].num;
		}
	}
	return distance;
}

void patch_distance( const double *pixelNumLab,const double *colorL,const double *colorA, const double *colorB,const double *colorNum,const double *count,const int maxlabel ,vector< vector<double> > &distance)
{
	vector< vector<N> > colorArray(maxlabel,0);
	int tmpcount = 0;
	for(int i=0;i<maxlabel;++i)
	{
		colorArray[i].resize((int)colorNum[i]);
		
		for(int j=0;j<(int)colorNum[i];++j)
		{
			colorArray[i][j].num = pixelNumLab[j+tmpcount]/count[i];
			colorArray[i][j].colorL = colorL[j+tmpcount];
			colorArray[i][j].colorA = colorA[j+tmpcount];
			colorArray[i][j].colorB = colorB[j+tmpcount];
		}
		tmpcount +=(int)colorNum[i];
	}
	// vector装每个patch的颜色Lab和比例num
	

	for(int i=0;i<maxlabel;++i)
	{
		for(int j=0;j<maxlabel;++j)
		{
			distance[i][j] = distanceIJ(colorArray[i],colorArray[j]);
		}
	}

}


void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
	double *pixelNumLab = mxGetPr(prhs[0]);
	double *L = mxGetPr(prhs[1]);
	double *a = mxGetPr(prhs[2]);
	double *b = mxGetPr(prhs[3]);
	double *colorNumEachPatch = mxGetPr(prhs[4]);
	double *pixelNumEachPatch = mxGetPr(prhs[5]);
	double *numPatches = mxGetPr(prhs[6]);
	
	vector< vector<double> > distance( (int)numPatches[0],vector<double>((int)numPatches[0],0) );

	patch_distance(pixelNumLab, L, a, b, colorNumEachPatch, pixelNumEachPatch,(int)numPatches[0] , distance);
	
	plhs[0] = mxCreateDoubleMatrix((int)numPatches[0],(int)numPatches[0],mxREAL);
	double *dis= mxGetPr(plhs[0]);
	for(int i=0;i<(int)numPatches[0];++i)
		for(int j=0;j<(int)numPatches[0];++j)
		{
			dis[(int)numPatches[0]*i+j] = distance[i][j];
		}
	
}

	



