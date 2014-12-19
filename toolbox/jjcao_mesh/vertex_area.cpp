/*=================================================================
*
* todo!!!!!!!!
* compute "voronoi" area of each vertex
* reference: Discrete Differential Geometry Operators for triangulated 2-manifolds_02
* usage: 
		va = vertex_area(verts, faces);
* inputs:
		verts: 3*nverts
		faces: 3*faces
		bVoronoi: use fast approximation or voronoi area
*
* output:
*		va: nverts*1
*          
*
* JJCAO, 2013
*
*=================================================================*/

#include <mex.h>
#include <Eigen/Dense>

using namespace std;

void compute_vertex_area_fast(Eigen::MatrixXd& V, Eigen::MatrixXd& F, double* va)
{
	int n1,n2,n3;
	Eigen::Vector3d v1,v2,v3;
	Eigen::Vector3d e1,e2;
	Eigen::Vector3d tmpV;
	double tmp;
	double farea;
	for( int i = 0; i < F.rows(); ++i)
	{
		n1 = F.coeff(i,0)-1; n2 = F.coeff(i,1)-1; n3 = F.coeff(i,2)-1;
		v1 = V.row(n1); v2 = V.row(n2); v3 = V.row(n3);
		e1 = v2 - v1; e2 = v3 - v1;

		tmpV = e1.cross(e2);
		farea = tmpV.norm()*0.5;
		tmp = farea / 3.0;
		
		va[n1] += tmp;
		va[n2] += tmp;
		va[n3] += tmp;
	}
}
void compute_vertex_voronoi_area(Eigen::MatrixXd& V, Eigen::MatrixXd& F, double* va)
{
	mexErrMsgTxt("not implemented still !");
}
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray*prhs[])
{
	///////////// Error Check
	if ( nrhs < 2) 
		mexErrMsgTxt("Number of input should > 1!");

	///////////// input & output arguments	
	// input 0: verts: 3*nverts
	int nverts = mxGetM(prhs[0]);
	int col = mxGetN(prhs[0]);
	if(col != 3)
		mexErrMsgTxt("The mesh must be triangle mesh! it is excepted to be n*3");

	double *verts = mxGetPr(prhs[0]);

	// input 1: faces: 3*nfaces
	int nfaces = mxGetM(prhs[1]);
	col = mxGetN(prhs[1]);
	if(col != 3)
		mexErrMsgTxt("The mesh must be triangle mesh! it is excepted to be n*3");

	double* faces = mxGetPr(prhs[1]);

	
	bool bVoronoi(false);
	if ( nrhs > 2) 
	{
		double* tmp = mxGetPr(prhs[2]);
		if ( *tmp != 0)
			bVoronoi = true;
	}

	///////////////////////////////////////////////
	// output 0
	plhs[0] = mxCreateDoubleMatrix( nverts, 1, mxREAL);   
	double *va = mxGetPr(plhs[0]);
	for ( int i = 0; i < nverts; ++i)
	{
		va[i] = 0;
	}

	///////////////////////////////////////////////	
	// process
	Eigen::MatrixXd V = Eigen::Map<Eigen::MatrixXd>(verts, nverts, 3);
	Eigen::MatrixXd F = Eigen::Map<Eigen::MatrixXd>(faces, nfaces, 3);

	if ( bVoronoi)
		compute_vertex_voronoi_area(V,F,va);
	else
		compute_vertex_area_fast(V,F,va);
}