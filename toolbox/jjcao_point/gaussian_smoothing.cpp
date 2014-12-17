#include "KDTree.h"
#include "mex.h"
#include <numeric>
const double pi = 3.1415926;

void retrieve_tree( const mxArray* matptr, KDTree* & tree){
    // retrieve pointer from the MX form
    double* pointer0 = mxGetPr(matptr);
    // check that I actually received something
    if( pointer0 == NULL )
        mexErrMsgTxt("vararg{1} must be a valid k-D tree pointer\n");
    // convert it to "long" datatype (good for addresses)
    long pointer1 = (long) pointer0[0];
    // convert it to "KDTree"
    tree = (KDTree*) pointer1;
    // check that I actually received something
    if( tree == NULL )
        mexErrMsgTxt("vararg{1} must be a valid k-D tree pointer\n");
    if( tree -> ndims() <= 0 )
        mexErrMsgTxt("the k-D tree must have k>0"); 
}
//% function [nfunc,nneigh] = gaussian_smoothing(verts, func, neighDist, sigma2, kdtree)
//% smoothing a scalar or vector function func defined on verts, using
//% Gaussian with sigma^2 = sigma2, and find neighbor vertices within
//% neighDist via kdtree.
//%
//% verts: n*m matrix, which are n m-dimentional points, such as n*3 in 3d space
//% func: n*m1 matrix, which is a function defined on verts, you can set func=verts to smooth the point set
//% neighDist: n*1 distance vector, neighDist(i) is used to collect neighbors of vertex i
//% sigma2: n*1 vector, sigma2(i) is sigma^2 for the Gaussian filter of vertex i
//%
//% nneigh: numbers of neighbors selected for each vertex
//%
//% jjcao @ 2014
//%
void mexFunction(int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[])
{
	/////////////////////////////////////////////////// retrieve input arguments
	if( nrhs!=5 )
		mexErrMsgTxt("This function requires 5 arguments\n");
    
	if( !mxIsNumeric(prhs[0]) )
		mexErrMsgTxt("varargin{0} must be a n*m matrix, which are n m-dimentional points, such as n*3 in 3d space\n");
	if( !mxIsNumeric(prhs[1]) )
		mexErrMsgTxt("varargin{1} must be a n*m1 matrix, which is a function defined on verts, you can set func=verts to smooth the point set\n");
	if( !mxIsNumeric(prhs[2]) )
		mexErrMsgTxt("varargin{2} n*1 a distance vector, neighDist(i) is used to collect neighbors of vertex i\n");
	if( !mxIsNumeric(prhs[3]) )
		mexErrMsgTxt("varargin{3} must be a n*1 vector, sigma2(i) is sigma^2 for the Gaussian filter of vertex i\n");	
    if( !mxIsNumeric(prhs[4]) )
		mexErrMsgTxt("varargin{4} must be a valid kdtree pointer\n");

	int row_vertex=mxGetM(prhs[0]); //获得矩阵的行数
	int column_vertex=mxGetN(prhs[0]); //获得矩阵的列数
    double *vertex = mxGetPr( prhs[0] );

	//int row_func=mxGetM(prhs[1]); //获得矩阵的行数
	int column_func=mxGetN(prhs[1]); //获得矩阵的列数
	double *func = mxGetPr( prhs[1] );

	 // retrieve neighDist, a vector or scalar
	int row_neighDist=mxGetM(prhs[2]);
    double *neighDist = mxGetPr( prhs[2] );

	 // retrieve sigma2, a vector or scalar
	int row_sigma2 = mxGetM(prhs[3]);
    double *sigma2 = mxGetPr( prhs[3] );

	// retrieve the tree pointer
    KDTree* tree;
    retrieve_tree( prhs[4], tree ); 
    
    /////////////////////////////////////////////////// Gaussian smoothing
	plhs[0] = mxCreateDoubleMatrix(row_vertex, column_func, mxREAL);
	double* result =(double*)mxGetPr(plhs[0]);
	plhs[1] = mxCreateNumericMatrix(row_vertex, 1, mxINT32_CLASS, mxREAL);
	int* nneigh =(int*)mxGetPr(plhs[1]);
	for (int i = 0;i < row_vertex;++i)
	{
		double rval(0.0);
		vector<int> idxsInRange;
		vector<double> dists;
		vector<double> point(column_vertex,0);
		for (int j = 0;j < column_vertex; ++j)
		{
			point[j] = vertex[j*row_vertex+i];			
		}
		tree->ball_query( point, neighDist[i], idxsInRange, dists );		
		
		nneigh[i] = idxsInRange.size();
		std::vector<double> gauss(nneigh[i],0);
		
		for ( int j = 0; j < nneigh[i]; ++j)
		{
			gauss[j] = exp(-dists[j]*dists[j]/(2*sigma2[i])) / sqrt(2*pi*sigma2[i]);
		}
		for (int j = 0;j < column_func; ++j)
		{
			int offset = j*row_vertex;
			//result[i + offset] = gaussian_smoothing(func, idxsInRange, dists);
			std::vector<double> value(nneigh[i],0);
			for ( int k = 0; k < nneigh[i]; ++k)
			{
				value[k] = func[idxsInRange[k] + offset];
			}			
			double tmp1 = std::inner_product(gauss.begin(),gauss.end(),value.begin(),0.0);
			double tmp2 = std::accumulate(gauss.begin(), gauss.end(), 0.0);
			result[i + offset] = tmp1/tmp2;
		}
	}
}
//#endif
