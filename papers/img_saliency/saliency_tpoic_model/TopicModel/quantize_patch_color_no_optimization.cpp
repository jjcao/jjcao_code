#include<mex.h>
#include<iostream>
#include<list>
#include<vector>
using namespace std;
#define dim 3

typedef struct N
{
	int  num;
	int  colorR;
	int  colorG;
	int  colorB;

bool operator < (N &N)
{
	return ( num >= N.num );
}
bool operator > (N &N)
{
	return ( num <= N.num );
}
} N; 


N*  Quantize( const double *Image, const int w,const float ratio,int &qcnum,const int binnum)
{
	vector< vector< vector<int> > > colorHistogram( binnum, vector< vector<int> >( binnum, vector<int>( binnum, 0 )));
	int index[3];
	float bin = 256.0/binnum*1.0;
	for(int i=0;i<w;++i)
	{
		index[0] = (int)Image[i]/bin;
		index[1] = (int)Image[i+w]/bin;
		index[2] = (int)Image[i+2*w]/bin;
		++colorHistogram[index[0]][index[1]][index[2]];
	}
    list<N> colorNum;
	N temp;
	
	for(int i=0;i<binnum;++i)
	for(int j=0;j<binnum;++j)
	for(int k=0;k<binnum;++k)
	{
		if(colorHistogram[i][j][k])
		{
			temp.num = colorHistogram[i][j][k];
			temp.colorR = i;
			temp.colorG = j;
			temp.colorB = k;
			colorNum.push_back(temp);
		}
	}
     colorNum.sort();
	 list<N>::iterator iter = colorNum.begin();
	 int Max = w*ratio;
	 int mark = 0;
	 int tmp = 0;
	 int tmp2 = colorNum.size();//共tmp2个颜色
	
	 while(mark<Max&&iter!=colorNum.end())
	 {
		 mark += iter->num;
		 ++iter;
		 ++tmp;
	 }   //第tmp个恰好取完 主要颜色index到tmp-1 
	 list<N>::iterator iter2 = colorNum.begin();
	 N *Histogram = new N [tmp2];
	 
	 for(int i=0;i<tmp2;++i)
	 {
		 Histogram[i].num = iter2->num;
		 Histogram[i].colorB = iter2->colorB;
		 Histogram[i].colorG = iter2->colorG;
		 Histogram[i].colorR = iter2->colorR;
		 ++iter2;
	 }
	 int distance ;
	 
	 for(int i = tmp;i<tmp2;++i)
	 { int index = 0, simVal = INT_MAX;
		 for(int j = 0; j<tmp ;++j)
		 {
			 distance = ( Histogram[i].colorB-Histogram[j].colorB)*( Histogram[i].colorB-Histogram[j].colorB)+( Histogram[i].colorG-Histogram[j].colorG)*( Histogram[i].colorG-Histogram[j].colorG)+( Histogram[i].colorR-Histogram[j].colorR)*( Histogram[i].colorR-Histogram[j].colorR);
			 if(distance<simVal)
				{simVal = distance;
				 index = j;}
		 }
		 ++Histogram[index].num;
		 Histogram[i].num = 0;
	 }//剩下5%放到前95%的颜色中  0~tmp-1 为直方图
	qcnum = tmp;
   
	return  Histogram;
	}

int* quantizePatches(const double *Image,const double *labels,const double *count,const int w,const int maxlable,const float ratio,vector< vector<N> > &colorArray,const int binnum)
{

	int bin = (256-0.1)/binnum;
	list<N> colorImage;
	N tmp;
	for(int i=0; i<w; ++i)
	{
		tmp.num = labels[i];
		tmp.colorR = Image[i];
		tmp.colorG = Image[i+w];
		tmp.colorB = Image[i+w*2];
		colorImage.push_back(tmp);
	}
	colorImage.sort();
	colorImage.reverse();
	 N *numAddImage = new N [w];//54321
	list<N>::iterator iter = colorImage.begin();
	for(int i=0;i<w;++i)
	{
		numAddImage[i].num = iter->num;
		numAddImage[i].colorB = iter->colorB;
		numAddImage[i].colorG = iter->colorG;
		numAddImage[i].colorR = iter->colorR;
		++iter;
	}
	
	 int k=0;
	for(int i=0; i<maxlable;++i)
	{
       
		double *temp =  new double [(int)count[i]*dim];
		memset(temp,0,(int)count[i]*dim*sizeof(double));
		for(int j=0;j<(int)count[i];++j)
		{
			temp[j] = (double)numAddImage[k+j].colorR;
			temp[j+(int)count[i]] = (double)numAddImage[k+j].colorG;
			temp[j+(int)count[i]*2] = (double)numAddImage[k+j].colorB;
		}
		k += (int)count[i];
		int qcnum = 0;
		N *Histogram = Quantize( temp, (int)count[i], ratio,qcnum,binnum);//69cuo
		delete[] temp;
		
		colorArray[i].resize(qcnum);
		for(int m = 0;m<qcnum;++m)
		{
			
			N tmpp;
			tmpp.num = Histogram[m].num;
			tmpp.colorB = Histogram[m].colorB;
			tmpp.colorG = Histogram[m].colorG;
			tmpp.colorR  =Histogram[m].colorR;
			colorArray[i][m] = tmpp;
		}
		delete[] Histogram;
			
	}
	int *patchQuantizedColor = new int [maxlable];
	for(int i=0;i<maxlable;++i)
	{
		patchQuantizedColor[i] = colorArray[i].size();
	}
	return patchQuantizedColor;
	//vector[i]为第i个patch中的颜色直方图
}
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
	double *img= mxGetPr(prhs[0]);
	double *labels = mxGetPr(prhs[1]);
	double *count = mxGetPr(prhs[2]);
	int w = mxGetN(prhs[0])/dim;
	double *maxlable = mxGetPr(prhs[3]);
	double *ratio = mxGetPr(prhs[4]);
	double *binnum = mxGetPr(prhs[5]);
	vector< vector<N> >  colorArray(maxlable[0],vector<N>(0));
	
	int *patchQuantizedColor = quantizePatches(img,labels,count, w,maxlable[0],ratio[0],colorArray,binnum[0]);
	int sum = 0;
	for(int i=0;i<maxlable[0];++i)
	{
		sum+=patchQuantizedColor[i];
	}

	plhs[0] = mxCreateDoubleMatrix(1,sum,mxREAL);//num
	plhs[1] = mxCreateDoubleMatrix(1,sum,mxREAL);//R
	plhs[2] = mxCreateDoubleMatrix(1,sum,mxREAL);//G
	plhs[3] = mxCreateDoubleMatrix(1,sum,mxREAL);//B
	plhs[4] = mxCreateDoubleMatrix(1,maxlable[0],mxREAL);//每一个patch颜色数目
	double *num,*R,*G,*B,*quantizeNum;
	num = mxGetPr(plhs[0]);
	R = mxGetPr(plhs[1]);
	G = mxGetPr(plhs[2]);
	B = mxGetPr(plhs[3]);
	quantizeNum = mxGetPr(plhs[4]);
	int l=0;
	int tmp1 = 0,tmp2 = 0;
		for(int i=0;i<maxlable[0];++i)
		{
			for(int j=0;j<colorArray[i].size();++j)
			{
			    
					num[l] = colorArray[i][j].num;
					R[l] = colorArray[i][j].colorR;
					G[l] = colorArray[i][j].colorG;
					B[l] = colorArray[i][j].colorB;
				    ++l;
					
			}
		}
		for(int i=0;i<maxlable[0];++i)
		{
			quantizeNum[i] = patchQuantizedColor[i];
		}
	
}