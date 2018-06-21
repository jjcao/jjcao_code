/*=================================================================
*
* compute various kinds of Laplacian weight
*
* usage: 
		L = perform_mesh_weight(verts, faces, type, options);
* inputs:
		verts: 3*nverts
		faces: 3*faces
		type: type is either 
			%       0: 'combinatorial' or 'graph': W(i,j)=1 is vertex i is conntected to vertex j.
			%       1: 'distance': W(i,j) = 1/d_ij^2 where d_ij is distance between vertex
			%           i and j.
			%       2: 'spring': W(i,j) = 1/d_ij where d_ij is distance between vertex
			%           i and j.
			%       3: 'conformal' or 'dcp': W(i,j) = (cot(alpha_ij)+cot(beta_ij))/2 where alpha_ij and
			%           beta_ij are the adjacent angle to edge (i,j). (do not offer W(i,j) = cot(alpha_ij)+cot(beta_ij) anymore)
			%           Refer to Computing discrete minimal surfaces and their conjugates_93, 
			%           Lemma 2 of On the convergence of metric and geometric properties of polyhedral surfaces_06 and 
			%           Characterizing Shape Using Conformal Factors_08.
			%           Refer to Skeleton Extraction by Mesh Extraction_08, and Intrinsic Parameterizations of Surface Meshes_02.
			%       4: 'Mean_curvature' or 'Laplace-Beltrami': W(i,j) = (1/area_i)*(cot(alpha_ij)+cot(beta_ij))/2 where alpha_ij and
			%           beta_ij are the adjacent angle to edge (i,j), area_i is the area of vertex i's Voroni vicinity. 
			%           Refer to Discrete Differential-Geometry Operators_for triangulated 2-manifolds_02
			%       5: 'Manifold-harmonic': W(i,j) = (1/sqrt(area_i*area_j))*(cot(alpha_ij)+cot(beta_ij))/2 where alpha_ij and
			%           beta_ij are the adjacent angle to edge (i,j). 
			%           Refer to Spectral Geometry Processing with Manifold Harmonics_08
			%       6: 'mvc': W(i,j) = [tan(/_kij/2)+tan(/_jil/2)]/d_ij where /_kij and /_jil are angles at i
		options.?:
*
*           4 or 5 can be built on 3 in matlab.          
*           just 0 & 3 are handled here. 1, 2 and 6 are left as todo.
*
* JJCAO, 2013
*
*=================================================================*/

#include <mex.h>
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <vector>

using namespace Eigen;
using namespace std;

typedef Triplet<double> T;
typedef SparseMatrix<double>::Index Index;
typedef SparseMatrix<double>::Scalar Scalar;


/// Return cotangent of (P,Q,R) corner (ie cotan of QP,QR angle).
double cotangent(const Vector3d& P,
                    const Vector3d& Q,
                    const Vector3d& R)
{

    Vector3d u = P - Q;
    Vector3d v = R - Q;
    // (u . v)/((u x v).len)
    double dot = u.dot(v);
    Vector3d cross_vector = u.cross(v);
    double cross_norm = std::sqrt(cross_vector.dot(cross_vector));
    if(cross_norm != 0.0)
        return (dot/cross_norm);
    else
        return 0.0; // undefined
}

////                                                  -> ->
///// Return tangent of (P,Q,R) corner (ie tangent of QP,QR angle).
//double tangent(const Point_3& P,
//                const Point_3& Q,
//                const Point_3& R)
//{
//    Vector_3 u = P - Q;
//    Vector_3 v = R - Q;
//    // (u . v)/((u x v).len)
//    double dot = (u*v);
//    CGAL_surface_mesh_parameterization_assertion(dot != 0.0);
//    Vector_3 cross_vector = CGAL::cross_product(u,v);
//    double cross_norm = std::sqrt(cross_vector*cross_vector);
//    if(dot != 0.0)
//        return (cross_norm/dot);
//    else
//        return 0.0; // undefined
//}

void compute_dcp_weight(double* verts, int* faces, int nverts, int nfaces, SparseMatrix<double>& sm)
{
	int i,j,m;
	double dtmp;
	std::vector<T> coef(nverts*6);	
	for (int k = 0; k < nfaces; ++k)
	{
		int tmp = 3*k;
		int face[3] = {faces[tmp],faces[tmp+1],faces[tmp+2]};			
		for (int l = 0; l < 3; ++l)
		{
			i = face[l]; j = face[(l+1)%3]; m = face[(l+2)%3];
			Vector3d p(verts[3*i],verts[3*i+1],verts[3*i+2]);
			Vector3d q(verts[3*m],verts[3*m+1],verts[3*m+2]);
			Vector3d r(verts[3*j],verts[3*j+1],verts[3*j+2]);
			dtmp = 0.5*cotangent(p, q, r);
			coef.push_back(T(i,j,dtmp));
		}			
	}
	{
		SparseMatrix<double> smtmp(nverts, nverts);
		smtmp.setFromTriplets(coef.begin(), coef.end());
		sm = SparseMatrix<double>(smtmp.transpose()) + smtmp;
	}
}

void perform_mesh_weight(double* verts, int nverts, int *faces, int nfaces, int type, double *vert_areas, SparseMatrix<double>& sm)
{	
	switch(type)
	{
	case 0: //'combinatorial' or 'graph'
		// ���ܴ�������
		{
		int i,j;	
		std::vector<T> coef(nverts*6);	
		for (int k = 0; k < nfaces; ++k)
		{
			int tmp = 3*k;
			int face[3] = {faces[tmp],faces[tmp+1],faces[tmp+2]};			
			for (int l = 0; l < 3; ++l)
			{
				i = face[l]; j = face[(l+1)%3];
				coef.push_back(T(i,j,1));
			}			
		}
		sm.setFromTriplets(coef.begin(), coef.end());
		//sm.coeffRef(1,1) = 2;sm.coeffRef(2,2) = 3;//sm.coeffRef(0,0) = 1;
		}
		break;
	case 3: //'conformal' or 'dcp'		
		compute_dcp_weight(verts, faces, nverts, nfaces, sm);
		break;
	//case 4: // 'Mean_curvature' or 'Laplace-Beltrami'	
	//	//if (vert_areas == 0)
	//	//{
	//	//	stringstream ss("options.vert_areas is not offered! ");	 
	//	//	mexErrMsgTxt(ss.str().c_str());
	//	//}
	//	compute_dcp_weight(verts, faces, nverts, nfaces, sm);
	//	//{
	//	//	for ( int i = 0; i < sm.rows(); ++i)
	//	//	{
	//	//		sm.row(i) = sm.row(i) * (1.0/vert_areas[i]);
	//	//	}
	//	//}
	//	break;
	default:
		stringstream ss("type: ");	 
		ss << type << " is not supported!";
		mexErrMsgTxt(ss.str().c_str());
	}

	sm.makeCompressed();
}

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray*prhs[])
{
	///////////// Error Check
	if ( nrhs < 3) 
		mexErrMsgTxt("Number of input should be > 2");
	if (1 != nlhs) 
		mexErrMsgTxt("Number of output should be 1");

	///////////// input & output arguments	
	// input 0: verts: 3*nverts
	int row = mxGetM(prhs[0]);
	int nverts = mxGetN(prhs[0]);
	if(row != 3)
		mexErrMsgTxt("The mesh must be triangle mesh! it is excepted to be 3*n");

	double *verts = mxGetPr(prhs[0]);

	// input 1: faces: 3*nfaces
	row = mxGetM(prhs[1]);
	int nfaces = mxGetN(prhs[1]);
	if(row != 3)
		mexErrMsgTxt("The mesh must be triangle mesh! it is excepted to be 3*n");

	double* dfaces = mxGetPr(prhs[1]);
    int* faces = new int[nfaces*3];
    for(int i = 0; i < nfaces*3; ++i)
    {
        faces[i] = int(dfaces[i]);
		--faces[i];
    }
    
	// input 3: type
	double* type = mxGetPr(prhs[2]);
	
	// input 4: options
	double* vert_areas(0);
	if ( nrhs > 3) 
	{
		mxArray* tmp;
		const mxArray *options = prhs[3];
		if ( mxSTRUCT_CLASS != mxGetClassID(options))
			mexErrMsgTxt("4th arguments is not a structure!");
		else
		{
			// options.vert_areas: 1*nverts
			tmp = mxGetField(options,0,"vert_areas");
			if (tmp)
				vert_areas = mxGetPr(tmp);// not used!
			//mexPrintf("%f, %f, %f\n", vert_areas[0],vert_areas[1],vert_areas[2]);	
		}
	}

	///////////////////////////////////////////////		
	SparseMatrix<double> sm(nverts, nverts);
	perform_mesh_weight(verts, nverts, faces, nfaces, *type, vert_areas, sm);

	SparseMatrix<double>::StorageIndex* innerInd = sm.innerIndexPtr();
	SparseMatrix<double>::StorageIndex* outerInd = sm.outerIndexPtr();
	Scalar* valuePtr = sm.valuePtr();

	///////////////////////////////////////////////
	// output 0
	//plhs[0] = mxCreateDoubleMatrix( nverts, nverts, mxREAL);
    mwSize nzmax= sm.nonZeros();
	plhs[0] = mxCreateSparse( nverts, nverts, nzmax, mxREAL);
	double *L = mxGetPr(plhs[0]);
	mwIndex *irs = mxGetIr(plhs[0]);//row index
    mwIndex *jcs = mxGetJc(plhs[0]);//column index	
	

	for (int k = 0; k<nzmax; ++k)
	{
		//double dtmp = valuePtr[k];
		//long ltmp = innerInd[k];
		L[k] = valuePtr[k];
		irs[k] = innerInd[k];
		//mexPrintf("%f %f\n", L[k], irs[k]);
	}

	for (int k = 0; k<nverts; ++k)
	{
		//long ltmp = outerInd[k];
		jcs[k] = outerInd[k];
	}
	//long ltmp = outerInd[nverts];
	jcs[nverts]=outerInd[nverts];
    
    delete[] faces;
}