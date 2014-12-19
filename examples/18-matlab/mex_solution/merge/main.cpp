#include<iostream>
#include<stdlib.h>
#include<mex.h>
#pragma comment(lib,"libmx.lib")
#pragma comment(lib,"libmex.lib")
#pragma comment(lib,"libmat.lib")

using namespace std;

void Merge(double Array[],double TempArray[],int left,int right,int middle)// jjcao, function name, lowercase
{
	for(int j = left;j <= right;j++)
		TempArray[j] = Array[j];
	int index1 = left;
	int index2 = middle+1;
	int i = left;
	while((index1 <= middle)&&(index2 <= right))
	{
		if(TempArray[index1] < TempArray[index2])
			Array[i++] = TempArray[index1++];
		else
			Array[i++] = TempArray[index2++];
	}
	while(index1 <= middle)
		Array[i++] = TempArray[index1++];
	while(index2 <= right)
		Array[i++] = TempArray[index2++];
}

void Sort(double Array[],double TempArray[],int left,int right)
{
	if(left<right)
	{
		int middle = (left+right)/2;
		Sort(Array,TempArray,left,middle);
		Sort(Array,TempArray,middle+1,right);
		Merge(Array,TempArray,left,right,middle);
	}
}

void mainfun( double p[],double Arr[],double Temp[],int n1,int n2 )
{
	//double Arr[20] = {2,44,56,7,8, 77,35,12,46,43, 27,90,59,71,3, 60,82,11,38,49};
	//double Temp[20];
	Sort( Arr,Temp,n1,n2 );
	//*
	for(int i=0;i<=n2;i++)
		p[ i ] = Arr[ i ];
	//*/
}

void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
	int M,N;
	double *data;
	double *play,*Temp;
	data=mxGetPr(prhs[0]);//获得指向输入数据矩阵的指针
	M=mxGetM(prhs[0]); //获得矩阵的行数
	N=mxGetN(prhs[0]); //获得矩阵的列数
	plhs[ 0 ] = mxCreateDoubleMatrix(M,N,mxREAL);
	plhs[ 1 ] = mxCreateDoubleMatrix(M,N,mxREAL);
	play = mxGetPr( plhs[0] );
	Temp = mxGetPr( plhs[1] );
	mainfun( play,data,Temp,M-1,N-1 );
}