// usage of Eigen
// jjcao

#include <iostream>
#include <vector>
#include <Eigen/Dense>

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
void eigen_matrix2array()
{
	Eigen::Matrix4d mat;
	mat << 1, 2, 3, 4,
           5, 6, 7, 8,
           9, 10, 11, 12,
		   0, 0, 0, 1;
	std::cout << mat << std::endl;

	// method 1
	//double *ptr = new double[4*4];	
	//Eigen::Map<Eigen::Matrix4d>(ptr, 4, 4) = mat; 
	// method 2
	double *ptr = mat.data(); 
	for ( int i = 0; i < 4; ++i)
	{
		std::cout << *ptr << ", "; ++ptr;
		std::cout << *ptr << ", "; ++ptr;
		std::cout << *ptr << ", "; ++ptr;
		std::cout << *ptr << std::endl; ++ptr;
	}
	//delete []ptr; // why error?
}
void array2eigen_matrix()
{
	//int data[] = {1,2,3,4};
	double *data = new double[4];
	for ( int i = 0; i < 4; ++i)
	{
		data[i] = i;
	}
	Eigen::Matrix2i mat221(data);
	mat221.coeffRef(0,0) = 1;
	mat221.coeffRef(0,1) = 2;
	//mat221.coeffRef(1,0) = 3;
	//mat221.coeffRef(1,1) = 4;
	std::cout << mat221 << std::endl;
	Eigen::MatrixXi mat222 = Eigen::Map<Eigen::Matrix2i>(data);
	std::cout << mat222 << std::endl;
	Eigen::MatrixXi mat223 = Eigen::Map<Eigen::MatrixXi>(data, 2, 2);
	std::cout << mat223 << std::endl;

}
void vector_matrix()
{
	Eigen::MatrixXd verts(2,3);
	verts << 0.20822,      0.12405,      0.87662,
        0.207,      0.12458,      0.81992;
	std::cout << verts << std::endl;

	Eigen::Vector3d tmp = verts.row(1);
	std::cout << tmp << std::endl;

	verts.row(1) = Eigen::Vector3d(1,1,1);
	verts.row(1) = verts.row(1) * 2;
	std::cout << verts << std::endl;

}
int main()
{
	vector_matrix();
	array2eigen_matrix();
	eigen_matrix2array();

	double center_[] = {0.0029885,     -0.31246,      0.23711};
	double vector_[] = { 0,   2.2204e-16,            1};
	Eigen::Vector3d center(0,0,0);
	Eigen::Vector3d v = Eigen::Map<Eigen::Vector3d>(vector_);
	for ( int i = 0; i < 3; ++i)
	{
		center.coeffRef(i) = center_[i];
	}
	std::cout << center << std::endl;
	std::cout << v << std::endl;

	/////////////////////////
	Eigen::Matrix4d rot = create_rotation3d_line_angle(
		center,v,3.1415926*0.5);
	std::cout << rot << std::endl;

	Eigen::MatrixXd verts(2,3);
	verts << 0.20822,      0.12405,      0.87662,
        0.207,      0.12458,      0.81992;
	std::cout << verts << std::endl;

	Eigen::MatrixXd newverts = transform_point3d(verts,rot);
	std::cout << newverts << std::endl;
	return 0;
}
