#include "KDTree.h"
#include "mex.h"
#include <numeric>

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

void retrieve_delta( const mxArray* matptr, double& delta ){
    // check that I actually received something
    if( matptr == NULL )
        mexErrMsgTxt("vararg{3} must be a scalar\n");

    if( 1 != mxGetM(matptr) || !mxIsNumeric(matptr) || 1 != mxGetN(matptr) )
    	mexErrMsgTxt("vararg{3} must be a scalar\n");    
    
    // retrieve point
	delta = mxGetScalar(matptr);	
}

double gaussian_weighted_curcature(vector<double> &point, double qradi, KDTree* tree, double* curvature)
{
    double rval(0.0);
    vector<int> idxsInRange;
    vector<double> dists;
    tree->ball_query( point, qradi, idxsInRange, dists );

    ///////////////////////////////////////////////
    vector<double> temp(idxsInRange.size(),0);
    for(int k = 0;k < temp.size(); ++k)
    {
        temp[k] = -dists[k]*dists[k];
        temp[k] = temp[k]/(qradi*qradi*0.5);
        temp[k] = exp(temp[k]);
    }
    for(int k = 0;k < temp.size(); ++k)
    {
        rval += curvature[idxsInRange[k]]*temp[k];
    }
    rval /= std::accumulate(temp.begin(), temp.end(), 0);

	return rval;
}

//GaussWeighcurvature= compute_gaussian_weighted_curvature(vertex,Cmean,delta,tree)
void mexFunction(int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[])
{
	//check number of arguments
	if( nrhs!=4 )
		mexErrMsgTxt("This function requires 4 arguments\n");
    
	if( !mxIsNumeric(prhs[0]) )
		mexErrMsgTxt("varargin{0} must be a vertex matrix\n");
	if( !mxIsNumeric(prhs[1]) )
		mexErrMsgTxt("varargin{1} must be a mean curvature vector\n");
	if( !mxIsNumeric(prhs[2]) )
		mexErrMsgTxt("varargin{2} must be a double (delta)\n");
    if( !mxIsNumeric(prhs[3]) )
		mexErrMsgTxt("varargin{3} must be a valid kdtree pointer\n");

	int row_vertex=mxGetM(prhs[0]); //获得矩阵的行数
	int column_vertex=mxGetN(prhs[0]); //获得矩阵的列数
    double *vertex = mxGetPr( prhs[0] );

	int row_curvature=mxGetM(prhs[1]); //获得矩阵的行数
	int column_curvature=mxGetN(prhs[1]); //获得矩阵的列数
	double *curvature = mxGetPr( prhs[1] );

	 // retrieve the delta
    double delta = 0;
    retrieve_delta(prhs[2], delta);

	// retrieve the tree pointer
    KDTree* tree;
    retrieve_tree( prhs[3], tree ); 
    
    /////////////////////////////////////////////////// retrieve the query point  
    vector<double> gw1(row_vertex,0),gw2(row_vertex,0);
	double qradi = 2*delta;
    double  qradii = 2*qradi;
	for (int i = 0;i < row_vertex;++i)
	{
		vector<double> point(column_vertex,0);
		for (int j = 0;j < column_vertex; ++j)
		{
			point[j] = vertex[j*row_vertex+i];			
		}
        gw1[i] = gaussian_weighted_curcature(point, qradi, tree, curvature);
        gw2[i] = gaussian_weighted_curcature(point, qradii, tree, curvature);
	}
    
    plhs[0] = mxCreateDoubleMatrix(row_vertex, 2, mxREAL);
    double* Arr =(double*)mxGetPr(plhs[0]);

    for(int i = 0;i < row_vertex;++i)
    {
        Arr[i] = gw1[i];
        Arr[i + row_vertex] = gw2[i];
    }    
    //printf("function end\n");   
}
//#endif
