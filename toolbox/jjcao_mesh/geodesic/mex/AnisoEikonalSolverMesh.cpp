//================================================================
//================================================================
// File: AnisoEikonalSolverMesh.cpp
// (C) 02/2010 by Fethallah Benmansour & Gabriel Peyrée
//================================================================
//================================================================

#include "AnisoEikonalSolverMesh.h"

void mexFunction(	int nlhs, mxArray *plhs[], 
				 int nrhs, const mxArray*prhs[] ) 
{
    //==================================================================
	/* retrive arguments */
	if( nrhs<7 ) 
		mexErrMsgTxt("at least 7 input arguments are required.");
	if( nlhs!=0  && nlhs!=2 ) 
		mexErrMsgTxt("0 or 2 output arguments are required.");
    //==================================================================
	// arg1 : Nb_calls: How many times this function has been called
    // This parameter avoids the recomputation of the mesh connectivity
	Nb_calls = (int) mxGetScalar(prhs[0]);
    //mexPrintf("Nb_calls = %d\n", Nb_calls);
    //NbMax_calls = (int) mxGetPr(prhs[0])[1];
    //==================================================================
	// arg2 : vertex
	vertex = mxGetPr(prhs[1]);
	nverts = mxGetN(prhs[1]); 
	if( mxGetM(prhs[1])!=3 )
		mexErrMsgTxt("vertex must be of size 3 x nverts."); 
	//==================================================================
    // arg3 : faces
	faces = mxGetPr(prhs[2]);
	nfaces= mxGetN(prhs[2]);
	if( mxGetM(prhs[2])!=3 )
		mexErrMsgTxt("face must be of size 3 x nfaces."); 
	//==================================================================
    // arg4 : Metric (should be symmetric definite positive on the tangent plan)
    // order of the 6 components for xx, yy, zz, xy, yz, zx
	if( (mxGetDimensions(prhs[3])[0] != 6) || (mxGetDimensions(prhs[3])[1] != nverts) )
        mexErrMsgTxt("T must be of same size as vertex with 6 components : 6xnverts.");
	T = mxGetPr(prhs[3]);
    //==================================================================
	// arg5 : start_points
	start_points = mxGetPr(prhs[4]);
	nstart = mxGetM(prhs[4]);
	//==================================================================	
	// arg 6 and 7 : if initial distance values at source points are given
    // In this case, provide the intial voronoi indices as well
	U_ini_seeds = mxGetPr(prhs[5]);
    V_ini_seeds = mxGetPr(prhs[6]);
	if( mxGetM(prhs[5])==0 && mxGetN(prhs[5])==0 )
		U_ini_seeds=NULL;
	if( mxGetM(prhs[6])==0 && mxGetN(prhs[6])==0 )
		V_ini_seeds=NULL;
	if( U_ini_seeds!=NULL && (mxGetM(prhs[5])!=nstart || mxGetN(prhs[5])!=1) )
        mexErrMsgTxt("values must be of size nb_start_points x 1."); 
	if( V_ini_seeds!=NULL && (mxGetM(prhs[6])!=nstart || mxGetN(prhs[6])!=1) )
		mexErrMsgTxt("values must be of size nb_start_points x 1.");
	//==================================================================
	// arg 8 : boolean array for region to update
    doUpdate = (bool*) mxGetPr(prhs[7]);
    //==================================================================
    if(nrhs == 8){
        // first ouput : geodesic distance
        plhs[0] = mxCreateNumericArray(1,&nverts, mxDOUBLE_CLASS, mxREAL ); 
    	U = (double*) mxGetPr(plhs[0]);
        // second output : voronoi
    	plhs[1] = mxCreateNumericArray(1,&nverts, mxINT16_CLASS, mxREAL ); 
        Vor = (short*) mxGetPr(plhs[1]);
        given_u = false;
    }
    else if(nrhs == 10){
        U   = mxGetPr(prhs[8]);
        Vor = (short*) mxGetPr(prhs[9]);
        given_u = true;
    }
	//==================================================================
    if(Nb_calls == 0){
        mexPrintf("Creating the Mesh\n");
        create_the_mesh();
	}
    //------------------------------------------------------------------
    InitializeArrays();
	//------------------------------------------------------------------
	InitializeQueue();
    //------------------------------------------------------------------
	GaussSiedelIterate();
    //------------------------------------------------------------------
    //Nb_calls++;
	//==================================================================
    DELETEARRAY(S);
    DELETEARRAY(Q);
    DELETEARRAY(ITER);
    DELETEARRAY(TAB);
    DELETEARRAY(U_n_D);
    //==================================================================
    //if(Nb_calls >= NbMax_calls){
    //    mexPrintf("Clening the geodesic Mesh\n");
    //    Mesh.~GW_Mesh();
	//}
    
};