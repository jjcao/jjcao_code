//================================================================
//================================================================
// File: AnisoEikonalSolverMatlabMesh.cpp
// (C) 02/2010 by Fethallah Benmansour
//================================================================
//================================================================

#include "AnisoEikonalSolverMatlabMesh.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    {
        //==================================================================
        /* retrive arguments */
        if (nrhs != 5 && nrhs != 7)
            mexErrMsgTxt("5 or 7 input arguments are required.");
        if (nlhs != 0 && nlhs != 2)
            mexErrMsgTxt("0 or 2 output arguments are required.");
        //==================================================================
        // arg1 : vertex
        vertex = mxGetPr(prhs[0]);
        nverts = mxGetN(prhs[0]);
        if (mxGetM(prhs[0]) != 3)
            mexErrMsgTxt("vertex must be of size 3 x nverts.");
        //==================================================================
        // arg2 : connectivity
        if ((int) mxGetN(prhs[1]) != nverts)
            mexErrMsgTxt("connectivity must have the same size as vertex.");
        Connectivity = (mxArray*) prhs[1];
        nb_neigh_field = mxGetFieldNumber(Connectivity, "nb_neighbors");
        neigh_idx_field = mxGetFieldNumber(Connectivity, "neighbours_idx");
        //==================================================================
        // arg3 : Metric (should be symmetric definite positive on the tangent plan)
        // order of the 6 components for xx, yy, zz, xy, yz, zx
        if (((int) mxGetDimensions(prhs[2])[0] != 6) || ((int) mxGetDimensions(prhs[2])[1] != nverts))
            mexErrMsgTxt("T must be of same size as vertex with 6 components : 6xnverts.");
        T = mxGetPr(prhs[2]);
        //==================================================================
        // arg4 : start_points
        start_points = mxGetPr(prhs[3]);
        nstart = mxGetM(prhs[3]);
        //==================================================================
        // arg5 : boolean array for region to update
        doUpdate = (bool*) mxGetPr(prhs[4]);
        //==================================================================
        if (nrhs == 5) {
            // first ouput : geodesic distance
            //      plhs[0] = mxCreateNumericArray(1,&nverts, mxDOUBLE_CLASS, mxREAL );
            //    	U = (double*) mxGetPr(plhs[0]);
            U = mxGetPr(plhs[0] = (mxArray*) mxCreateDoubleMatrix(nverts, 1, mxREAL));
            // second output : voronoi
            //plhs[1] = mxCreateNumericArray(1, &nverts, mxINT16_CLASS, mxREAL);
            //Vor = (short*) mxGetPr(plhs[1]);
            Vor = (short*) mxGetPr(plhs[1] = mxCreateNumericArray(1, (mwSize*) &nverts, mxINT16_CLASS, mxREAL));
            given_u = false;
        } else if (nrhs == 7) {
            U = mxGetPr(prhs[5]);
            Vor = (short*) mxGetPr(prhs[6]);
            given_u = true;
        }
        //------------------------------------------------------------------
        InitializeArrays();
        //------------------------------------------------------------------
        InitializeQueue();
        //------------------------------------------------------------------
        GaussSiedelIterate();
        //==================================================================
        DELETEARRAY(S);
        DELETEARRAY(Q);
        DELETEARRAY(ITER);
        DELETEARRAY(TAB);
        DELETEARRAY(U_n_D);
        DELETEARRAY(Utmp);
        DELETEARRAY(Vtmp);
        DELETEARRAY(heap);
        DELETEARRAY(heapIndex);
        //==================================================================
    }
}