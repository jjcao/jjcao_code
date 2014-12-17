//================================================================
//================================================================
// File: AnisoEikonalSolverMesh.h
// (C) Fethallah Benmansour 2010 -- CVLab-EPFL --
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

#define isEqual(a, b) (a==b ? 1 : 0);
#define kSeed -1 
#define kEstimated -2
#define kFar -3
#define kBorder -4
#define kEnqueued 0
#define kUnqueued 1
#define INFINITE 1e9
#define tol 1e-15
#define DELETEARRAY(p) {if (p!=NULL) delete [] p; p=NULL;}

queue<GW_U32> waiting_Points;

int Nb_calls;
//int NbMax_calls;
double* vertex = NULL;
int nverts = -1; 
double* faces = NULL;
int nfaces = -1;
GW_GeodesicMesh Mesh;
double* start_points = NULL;
int nstart = -1;
double* T = NULL;	// weight

double* U_ini_seeds = NULL;
double* V_ini_seeds = NULL;

bool given_u;
// outputs 
double *U        = NULL;	// distance
short  *Vor      = NULL;	// Voronoi
short  *S        = NULL;	// state
short  *Q        = NULL;
bool   *doUpdate = NULL;
int    *ITER     = NULL;
double *TAB      = NULL;
double *U_n_D    = NULL;
double *X1, *X2, *X12;
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

//================================================================
void InitializeArrays()
//================================================================
{
	//--------------------------------------------------------
    int x;
    //--------------------------------------------------------
    S      = new short[nverts];
    Q      = new short[nverts];
    ITER   = new int[nverts];
    TAB    = new double[3];
    U_n_D  = new double[2];
    X1     = new double[3];
    X2     = new double[3];
    X12    = new double[3];
    //--------------------------------------------------------
    if(given_u){
        for(x = 0; x < nverts; x++){
            Q[x]    = kUnqueued;
            S[x]    = kEstimated;
            ITER[x] = 0;
        }
    }
    else{
        for(x = 0; x < nverts; x++){
            U[x]    = INFINITE;
            Vor[x]  = kBorder;
            Q[x]    = kUnqueued;
            S[x]    = kBorder;
            ITER[x] = 0;
        }
    }
    //--------------------------------------------------------
};

//================================================================
void InitializeQueue()
//================================================================
{
	int i, point, x;
    GW_U32 npoint;
    short maxVor = -1;
    if(given_u){
        for(x = 0; x < nverts; x++)
            maxVor = MAX(maxVor, Vor[x]);
        if(U_ini_seeds == NULL){
            for( i=0; i<nstart; ++i){
                point = (int) start_points[i];
                if( point >= nverts)
                    mexErrMsgTxt("start_points should be in the domain.");
                //--------------------------------------------------------
                U[point] = 0.0;
                S[point] = kSeed;
                Vor[point] = i + 1 + maxVor;
                //--------------------------------------------------------
                GW_GeodesicVertex* v = (GW_GeodesicVertex*) Mesh.GetVertex((GW_U32) point);
                for( GW_VertexIterator VertIt = v->BeginVertexIterator(); VertIt!=v->EndVertexIterator(); ++VertIt ){
                    GW_GeodesicVertex* pNewVert = (GW_GeodesicVertex*) *VertIt;
                    GW_ASSERT( pNewVert!=NULL );
                    npoint = pNewVert->GetID();
                    if( (S[npoint] != kSeed) ){
                        waiting_Points.push(npoint);
                        Q[npoint] = kEnqueued;
                    }
                }
                //--------------------------------------------------------
            }
        }
        else{ // intialize using given intial values
            for( i=0; i<nstart; ++i){
                point = (int) start_points[i];
                if( point >= nverts)
                    mexErrMsgTxt("start_points should be in the domain.");
                //--------------------------------------------------------
                U[point]   = U_ini_seeds[i]; 
                S[point]   = kEstimated;// Or kSeed ? 
                Vor[point] = V_ini_seeds[i];
                //--------------------------------------------------------
                GW_GeodesicVertex* v = (GW_GeodesicVertex*) Mesh.GetVertex((GW_U32) point);
                for( GW_VertexIterator VertIt = v->BeginVertexIterator(); VertIt!=v->EndVertexIterator(); ++VertIt ){
                    GW_GeodesicVertex* pNewVert = (GW_GeodesicVertex*) *VertIt;
                    GW_ASSERT( pNewVert!=NULL );
                    npoint = pNewVert->GetID();
                    if( (S[npoint] != kSeed) ){
                        waiting_Points.push(npoint);
                        Q[npoint] = kEnqueued;
                    }
                }
                //--------------------------------------------------------
            }
        }
    }
    else{
        if(U_ini_seeds == NULL){
            for( i=0; i<nstart; ++i){
                point = (int) start_points[i];
                if( point >= nverts)
                    mexErrMsgTxt("start_points should be in the domain.");
                //--------------------------------------------------------
                U[point] = 0.0;
                S[point] = kSeed; 
                Vor[point] = i;
                //--------------------------------------------------------
                GW_GeodesicVertex* v = (GW_GeodesicVertex*) Mesh.GetVertex((GW_U32) point);
                for( GW_VertexIterator VertIt = v->BeginVertexIterator(); VertIt!=v->EndVertexIterator(); ++VertIt ){
                    GW_GeodesicVertex* pNewVert = (GW_GeodesicVertex*) *VertIt;
                    GW_ASSERT( pNewVert!=NULL );
                    npoint = pNewVert->GetID();
                    if( (S[npoint] != kSeed) ){
                        waiting_Points.push(npoint);
                        Q[npoint] = kEnqueued;
                    }
                }
                //--------------------------------------------------------
            }
        }
        else{ // intialize using given intial values
            for( i=0; i<nstart; ++i){
                point = (int) start_points[i];
                if( point >= nverts)
                    mexErrMsgTxt("start_points should be in the domain.");
                //--------------------------------------------------------
                U[point]   = U_ini_seeds[i]; 
                S[point]   = kEstimated; 
                Vor[point] = V_ini_seeds[i]; ;
                //--------------------------------------------------------
                GW_GeodesicVertex* v = (GW_GeodesicVertex*) Mesh.GetVertex((GW_U32) point);
                for( GW_VertexIterator VertIt = v->BeginVertexIterator(); VertIt!=v->EndVertexIterator(); ++VertIt ){
                    GW_GeodesicVertex* pNewVert = (GW_GeodesicVertex*) *VertIt;
                    GW_ASSERT( pNewVert!=NULL );
                    npoint = pNewVert->GetID();
                    if( (S[npoint] != kSeed) ){
                        waiting_Points.push(npoint);
                        Q[npoint] = kEnqueued;
                    }
                }
                //--------------------------------------------------------
            }
        }
    }
};


//================================================================
void DotProductMetric(double *res, double* M, double *V1, double *V2)
/*
 COMMENTS : 
 * compute V1^t M V2
 * M is symetric
*/
//================================================================
{
    *res    = M[0]*V1[0]*V2[0] + M[1]*V1[1]*V2[1] + M[2]*V1[2]*V2[2]
            + M[3]*( V1[0]*V2[1] + V2[0]*V1[1] )
            + M[4]*( V1[1]*V2[2] + V2[1]*V1[2] )
            + M[5]*( V1[2]*V2[0] + V2[2]*V1[0] );
};

//================================================================
void TsitsiklisOnePoint(double *res, double* M, double u, double* Z)
//================================================================
{
    double NormZ_M;
    DotProductMetric(&NormZ_M, M, Z, Z);
    *res = u + sqrt(NormZ_M);
}

//================================================================
void TsitsiklisTwoPoints(double* res, double* M, double  k, double u, double* z1, double* z2)
/*
COMMENTS : 
 * computes the minimum and the arg min of :
 * \alpha*k + u + \| \alpha* z_1 + z_2 \|_M ; \alpha \in [0, 1]
 * res is a double array of size 2 containing the argmin and the min
*/
//================================================================
{
    double r11, r22, r12, R;
    DotProductMetric(&r11, M, z1, z1);
    DotProductMetric(&r22, M, z2, z2);
    DotProductMetric(&r12, M, z1, z2);
    R = r11*r22-r12*r12;
    if(!(R > 0))
        mexErrMsgTxt("z1 and z2 should not be collinear !!");
    if( k >= sqrt(r11) ){
        res[0] = 0.0;
        res[1] = u + sqrt(r22);
    }
    else if(k <= -sqrt(r11) ){
        res[0] = 1.0;
        res[1] = k + u + sqrt(r11 + r22 + 2.0*r12);
    }
    else{
        if( r12>= -k*sqrt(R/(r11-k*k)) ){
            res[0] = 0.0;
            res[1] = u + sqrt(r22);
        }
        else if( r12 <= (-r11-k*sqrt(R/(r11-k*k)))  ){
            res[0] = 1.0;
            res[1] = k + u + sqrt(r11 + r22 + 2.0*r12);
        }
        else{
            res[0] = -(r12 + k*sqrt(R/(r11-k*k))) / r11;
            res[1] = res[0]*k + u + sqrt(R/(r11-k*k));
        }
    }
};

//================================================================
void TsitsiklisTriangle(double* res, GW_GeodesicVertex* pVert, GW_GeodesicFace* pFace)
//================================================================
{
    int	  point, npoint1, npoint2, i;
    double k, u;
    double* M;
    double U1, U2, V1, V2;
	double Ur;
    bool b_point1, b_point2;
    //--------------------------------------------------------------
    point = pVert->GetID();
    M = T + 6*point;
    Ur = U[point];
    //--------------------------------------------------------------
    res[1] = Ur+1.0;// assign bigger value
    //--------------------------------------------------------------
    GW_GeodesicVertex* pVert1 = (GW_GeodesicVertex*) pFace->GetNextVertex( *pVert );
	GW_ASSERT( pVert1!=NULL );
	GW_GeodesicVertex* pVert2 = (GW_GeodesicVertex*) pFace->GetNextVertex( *pVert1 );
	GW_ASSERT( pVert2!=NULL );
    //----------------------------------------------------------
    npoint1 = pVert1->GetID();
    npoint2 = pVert2->GetID();
    //----------------------------------------------------------
	b_point1 = ((S[npoint1]==kEstimated) || (S[npoint1]==kSeed)) && doUpdate[npoint1];
    b_point2 = ((S[npoint2]==kEstimated) || (S[npoint2]==kSeed)) && doUpdate[npoint2];
    if(b_point1){  U1 = U[npoint1]; V1 = Vor[npoint1];}
    if(b_point2){  U2 = U[npoint2]; V2 = Vor[npoint2];}
    //----------------------------------------------------------
    GW_Vector3D Edge1  = pVert->GetPosition() - pVert1->GetPosition();
	GW_Vector3D Edge2  = pVert->GetPosition() - pVert2->GetPosition();
	GW_Vector3D Edge12 = Edge1 - Edge2;
    //----------------------------------------------------------
    for(i = 0; i < 3; i++){
        X1[i]  = (double) Edge1[i];
        X2[i]  = (double) Edge2[i];
        X12[i] = (double) Edge12[i];
	}
    //----------------------------------------------------------
    if(b_point1 && b_point2){
        k = U1-U2; u = U2;
        TsitsiklisTwoPoints(res, M, k, u, X12, X2);
        res[2] = (res[0]>0.5 ? V1 : V2);
    }
    //----------------------------------------------------------
    else if(b_point1 && !b_point2 ){
        u = U1;
        TsitsiklisOnePoint(res, M, u, X1);
        res[1] = res[0];
        res[0] = 1.0;
        res[2] = V1;
    }
    //----------------------------------------------------------
    else if(!b_point1 && b_point2 ){
        u = U2;
        TsitsiklisOnePoint(res, M, u, X2);
        res[1] = res[0];
        res[0] = 0;
        res[2] = V2;
    }
};


//================================================================
void TsitsiklisUpdate(int point, double* res)
/*
COMMENTS : 
*/
//================================================================
{
	int i;
    bool is_updated = false;
    double Ur = U[point], Vr;
    //--------------------------------------------------------------
    GW_GeodesicVertex* pVert = (GW_GeodesicVertex*) Mesh.GetVertex((GW_U32) point);
    //--------------------------------------------------------------
    for( GW_FaceIterator FaceIt=pVert->BeginFaceIterator(); FaceIt!=pVert->EndFaceIterator(); ++FaceIt )
	{
        GW_GeodesicFace* pFace = (GW_GeodesicFace*) *FaceIt;
		//----------------------------------------------------------
        TsitsiklisTriangle(TAB, pVert, pFace);
        //----------------------------------------------------------
        if (TAB[1]<Ur){
            Ur   = TAB[1];
            Vr   = TAB[2];
            is_updated = true;
    	}
	}
	//--------------------------------------------------------------
	if (is_updated){
		res[0] = Ur;
        res[1] = Vr;
	}
    else{// copy old values
		res[0] = U[point];
        res[1] = Vor[point];
	}
    //--------------------------------------------------------------
};

//================================================================
void GaussSiedelIterate()
//================================================================
{
    int point,npoint,k;
    int iter = 0;
	//------------------------------------------------------------
    while( !waiting_Points.empty() ) {
        point = waiting_Points.front();
        waiting_Points.pop();
        TsitsiklisUpdate(point, U_n_D);
        S[point] = kEstimated;
        Q[point] = kUnqueued;
        ITER[point] = ITER[point] + 1;
        iter = MAX(iter, ITER[point]);
        if(fabs(U_n_D[0] - U[point]) > tol){
            //----------------------------------------------------
            U[point]   = U_n_D[0];
            Vor[point] = U_n_D[1];
            //----------------------------------------------------
            GW_GeodesicVertex* v = (GW_GeodesicVertex*) Mesh.GetVertex((GW_U32) point);
            for( GW_VertexIterator VertIt = v->BeginVertexIterator(); VertIt!=v->EndVertexIterator(); ++VertIt ){
                GW_GeodesicVertex* pNewVert = (GW_GeodesicVertex*) *VertIt;
                npoint = (int) pNewVert->GetID();
                if( (S[npoint] != kSeed) && (Q[npoint] != kEnqueued) && (doUpdate[point])){
                    waiting_Points.push(npoint);
                    Q[npoint] = kEnqueued;
                }
            }
        }
	}
};