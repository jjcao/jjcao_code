/*=================================================================
* Transform vertices (n*3) with a 3D affine transform
* 
* according to David Legland's matlab toolbox: geom3d 
* 
* Junjie Cao, 2013, jjcao1231@gmail.com
*
=================================================================*/

#include <mex.h>
#include <Eigen/Dense>
#include <vector>

Eigen::MatrixXd transform_point3d(Eigen::MatrixXd& verts, Eigen::Matrix4d& trans)
{
	int NP  = verts.rows();
	Eigen::MatrixXd res(NP,4);
	res.col(0) = verts.col(0);
	res.col(1) = verts.col(1);
	res.col(2) = verts.col(2);
	res.col(3) = Eigen::VectorXd::Ones(NP,1);
	//std::cout << res << std::endl;

	//std::cout << trans * res.row(0).transpose() << std::endl;
	res = res * trans.transpose();
	//std::cout << res << std::endl;
	res.conservativeResize(NP, 3);
	 
	return res;
}
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray*prhs[])
{  
	////////////////// parse input
	if ( nrhs != 2)
		mexErrMsgTxt("2 arguments needed!");

	double *verts_, *trans_;
	verts_ = mxGetPr(prhs[0]);
	trans_ = mxGetPr(prhs[1]);
	int nverts = mxGetM(prhs[0]);
	if ( mxGetN(prhs[0]) != 3)
		mexErrMsgTxt("center must be n*3 matrix!");
	if ( mxGetN(prhs[1]) != 4 && mxGetM(prhs[1]) != 4)
		mexErrMsgTxt("vector must be 4*4 matrix!");
	
	Eigen::MatrixXd verts = Eigen::Map<Eigen::MatrixXd>(verts_, nverts, 3);
	Eigen::Matrix4d trans = Eigen::Map<Eigen::Matrix4d>(trans_);

	/////////////////////////////
	Eigen::MatrixXd newverts = transform_point3d(verts,trans);

	//////////////////////// output
	plhs[0] = mxCreateDoubleMatrix(nverts,3,mxREAL);
	double *newverts_ = mxGetPr(plhs[0]);
	Eigen::Map<Eigen::MatrixXd>(newverts_, nverts, 3) = newverts;  
}

