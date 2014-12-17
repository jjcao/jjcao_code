#include <iostream>
#include <Eigen/Dense>
#include <Eigen/SVD>
using namespace Eigen;
using namespace std;

void construct()
{
	//////////////
	int row(3),col(5);
	MatrixXd m2(row,col);
	double* data = m2.data();
	for(int i = 0; i < row*col; ++i) data[i] = i;
	cout << "Column-major:\n" << m2 << endl;
	MatrixXd m3 = m2.block(0,0,3,4);
	cout << "after resize:\n" << m3 << endl;

	//////////
	MatrixXd m(2,2);
	m(0,0) = 3;
	m(1,0) = 2.5;
	m(0,1) = -1;
	m(1,1) = m(1,0) + m(0,1);
	std::cout << "Here is the matrix m:\n" << m << std::endl;
	VectorXd v(2);
	v(0) = 4;
	v(1) = v(0) - 1;
	std::cout << "Here is the vector v:\n" << v << std::endl;

	///////
	Matrix3f m1;
	m1 << 1, 2, 3,
	4, 5, 6,
	7, 8, 9;
	std::cout << m1;

}
void svd()
{	
	MatrixXf m = MatrixXf::Random(3,8);
	cout << "Here is the matrix m:" << endl << m << endl;
	JacobiSVD<MatrixXf> svd(m, ComputeThinU | ComputeThinV);
	cout << "Its singular values are:" << endl << svd.singularValues() << endl;
	cout << "Its left singular vectors are the columns of the thin U matrix:" << endl << svd.matrixU() << endl;
	cout << "Its right singular vectors are the columns of the thin V matrix:" << endl << svd.matrixV() << endl;
	Vector3f v = svd.matrixU().col(2);
	cout << v << endl;
	Vector3f rhs(1, 0, 0);
	cout << "Now consider this rhs vector:" << endl << rhs << endl;
	cout << "A least-squares solution of m*x = rhs is:" << endl << svd.solve(rhs) << endl;


	MatrixXf m2 = m * m.transpose();
	SelfAdjointEigenSolver<Matrix3f> eigensolver(m2);
	cout << "The eigenvalues of A are:\n" << eigensolver.eigenvalues() << endl;
	cout << "Here's a matrix whose columns are eigenvectors of A \n"
	<< "corresponding to these eigenvalues:\n"
	<< eigensolver.eigenvectors() << endl;
}

int main(int argc, char** args)
{
	construct();
	svd();
}