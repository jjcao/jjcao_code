/*=================================================================
* Create rotation around a line by an angle theta
* 
* according to David Legland's matlab toolbox: geom3d 
* 
* Junjie Cao, 2013, jjcao1231@gmail.com
*
=================================================================*/

#include <mex.h>
#include <Eigen/Dense>
#include <vector>

Eigen::Matrix4d create_translation3d(Eigen::Vector3d center)
{
	double dx(center.coeff(0)),dy(center.coeff(1)),dz(center.coeff(2));
	Eigen::Matrix4d trans;
	trans << 1,0,0,dx,
		    0,1,0,dy,
			0,0,1,dz,
			0,0,0,1;
	return trans;
}
Eigen::Matrix4d recenter_transform3d(Eigen::Matrix4d &transfo, Eigen::Vector3d& center)
{
	//% remove former translation part
	Eigen::Matrix4d res = Eigen::Matrix4d::Identity();
	res.block(0, 0, 3, 3) = transfo.block(0, 0, 3, 3);

	//% create translations
	Eigen::Matrix4d t1 = create_translation3d(-center);
	Eigen::Matrix4d t2 = create_translation3d(center);

	//% compute translated transform
	res = t2*res*t1;
	return res;
}
Eigen::Matrix4d create_rotation3d_line_angle(Eigen::Vector3d& center,Eigen::Vector3d& v, double theta)
{
	//% normalize vector
	v.normalize();

	//% compute projection matrix P and anti-projection matrix
	Eigen::Matrix3d P = v * v.transpose();
	Eigen::Matrix3d Q;
	Q << 0, -v.coeff(2), v.coeff(1),
		v.coeff(2), 0, -v.coeff(0),
		-v.coeff(1), v.coeff(0), 0;
	Eigen::Matrix3d I = Eigen::Matrix3d::Identity();

	//% compute vectorial part of the transform
	Eigen::Matrix4d mat = Eigen::Matrix4d::Identity();
	mat.block(0, 0, 3, 3) = P + (I - P)*cos(theta) + Q*sin(theta);

	//% add translation coefficient
	mat = recenter_transform3d(mat, center);
	
	return mat;
}
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray*prhs[])
{  
	////////////////// parse input
	if ( nrhs != 3)
		mexErrMsgTxt("3 arguments needed!");

	double *center_, *vector_, *theta_;
	center_ = mxGetPr(prhs[0]);
	vector_ = mxGetPr(prhs[1]);
	theta_ = mxGetPr(prhs[2]);

	if ( mxGetN(prhs[0]) != 3)
		mexErrMsgTxt("center must be 1*3 matrix!");
	if ( mxGetN(prhs[1]) != 3)
		mexErrMsgTxt("vector must be 1*3 matrix!");
	
	Eigen::Vector3d center = Eigen::Map<Eigen::Vector3d>(center_);
	Eigen::Vector3d v = Eigen::Map<Eigen::Vector3d>(vector_);

	/////////////////////////////
	Eigen::Matrix4d rot = create_rotation3d_line_angle(center,v,*theta_);

	//////////////////////// output
	plhs[0] = mxCreateDoubleMatrix(4,4,mxREAL);
	double *rot_ = mxGetPr(plhs[0]);
	Eigen::Map<Eigen::Matrix4d>(rot_, 4, 4) = rot;  // may contain mem error
}

