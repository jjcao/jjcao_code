//================================================================
//================================================================
// File: ComputeMeshConnectivity.cpp
// (C) 02/2010 by Fethallah Benmansour & Gabriel Peyré
//================================================================
//================================================================


#include <math.h>
#include "config.h"
#include <algorithm>
#include <map>
#include <vector>
#include <list>
#include <string>
#include <iostream>
#include <fstream>
#include <string.h>
#include <queue>
using std::string;
using std::cerr;
using std::cout;
using std::endl;

#include "mex.h"
#include "gw/gw_core/GW_Config.h"
#include "gw/gw_core/GW_MathsWrapper.h"
#include "gw/gw_core/GW_Mesh.h"
#include "gw/gw_geodesic/GW_GeodesicMesh.h"
using namespace GW;
using namespace std;

double *vertex = NULL;
int     nverts;
double *faces  = NULL;
int     nfaces;
GW_GeodesicMesh Mesh;

#define faces_(k,i) faces[k+3*i]
#define vertex_(k,i) vertex[k+3*i]

//================================================================
void create_the_mesh()
//================================================================
{
    int i;
	//--------------------------------------------------------
    //Mesh.GW_Mesh();
    Mesh.SetNbrVertex(nverts);
    Mesh.SetNbrFace(nfaces);
    //--------------------------------------------------------
	for( i=0; i<nverts; ++i )
	{
		GW_GeodesicVertex& vert = (GW_GeodesicVertex&) Mesh.CreateNewVertex();
		vert.SetPosition( GW_Vector3D(vertex_(0,i),vertex_(1,i),vertex_(2,i)) );
		Mesh.SetVertex(i, &vert);
	}
    //--------------------------------------------------------
	for( i=0; i<nfaces; ++i )
	{
		GW_GeodesicFace& face = (GW_GeodesicFace&) Mesh.CreateNewFace();
		GW_Vertex* v1 = Mesh.GetVertex((int) faces_(0,i)); GW_ASSERT( v1!=NULL );
		GW_Vertex* v2 = Mesh.GetVertex((int) faces_(1,i)); GW_ASSERT( v2!=NULL );
		GW_Vertex* v3 = Mesh.GetVertex((int) faces_(2,i)); GW_ASSERT( v3!=NULL );
		face.SetVertex( *v1,*v2,*v3 );
		Mesh.SetFace(i, &face);
	}
    //--------------------------------------------------------
	Mesh.BuildConnectivity();
    //--------------------------------------------------------
};


void mexFunction(	int nlhs, mxArray *plhs[], 
				 int nrhs, const mxArray*prhs[] ) 
{
	// arg1 : vertex
	vertex = (double*) mxGetPr(prhs[0]);
	nverts = mxGetN(prhs[0]); 
	if( mxGetM(prhs[0])!=3 )
		mexErrMsgTxt("vertex must be of size 3 x nverts."); 
	//==================================================================
    // arg2 : faces
	faces = (double*) mxGetPr(prhs[1]);
	nfaces= mxGetN(prhs[1]);
	if( mxGetM(prhs[1])!=3 )
		mexErrMsgTxt("face must be of size 3 x nfaces."); 
	//==================================================================
    create_the_mesh();
    //==================================================================
    const char *field_names[] = {"nb_neighbors", "neighbours_idx"};
	int dims[2] = {1,nverts};
    int nb_fields = 2;
	plhs[0] = mxCreateStructArray(2,dims,nb_fields,field_names);
    int nb_neigh_field = mxGetFieldNumber(plhs[0],"nb_neighbors");
    int neigh_idx_field = mxGetFieldNumber(plhs[0],"neighbours_idx");
    //==================================================================
	double* nb_neighbors;
    double* neigh_idx;
    int nb_neigh;
    int point;
    //------------------------------------------------------------------
    for( point = 0; point < nverts; point++){
        nb_neigh = 0;
        mxArray *field_value1;
        field_value1 = mxCreateDoubleMatrix(1,1,mxREAL);
    	nb_neighbors = mxGetPr(field_value1);
        GW_GeodesicVertex* v = (GW_GeodesicVertex*) Mesh.GetVertex((GW_U32) point);
        for( GW_VertexIterator VertIt = v->BeginVertexIterator(); VertIt!=v->EndVertexIterator(); ++VertIt ){
        	GW_GeodesicVertex* pNewVert = (GW_GeodesicVertex*) *VertIt;
            GW_ASSERT( pNewVert!=NULL );
            nb_neigh++;
        }
        nb_neighbors[0] = double(nb_neigh);
        mxSetFieldByNumber(plhs[0],point,nb_neigh_field,field_value1);
        //--------------------------------------------------------------
        mxArray *field_value2;
        field_value2 = mxCreateDoubleMatrix(1,nb_neigh,mxREAL);
    	neigh_idx = mxGetPr(field_value2);
        nb_neigh = 0;
        for( GW_VertexIterator VertIt = v->BeginVertexIterator(); VertIt!=v->EndVertexIterator(); ++VertIt ){
        	GW_GeodesicVertex* pNewVert = (GW_GeodesicVertex*) *VertIt;
            GW_ASSERT( pNewVert!=NULL );
            neigh_idx[nb_neigh++] = double(pNewVert->GetID());
        }
        mxSetFieldByNumber(plhs[0],point,neigh_idx_field,field_value2);
    }    
    //------------------------------------------------------------------
};