//================================================================
//================================================================
// File: AnisoEikonalSolverMatlabMesh.h
// (C) Fethallah Benmansour 2010 -- CVLab-EPFL --
// and Thomas Satzger 2010 - TU München
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

#define FIRST	1						/* root of heap */
#define OFF		FIRST-1
#define QUEUE_NOT_EMPTY (last >= FIRST)
#define LEFT_SON(k) (2*k)				/* genealogy in the heap (first = 1) */
#define FATHER(k)  (k/2)				/* genealogy in the heap (first = 1) */

queue<int> waiting_Points;

double* vertex = NULL;
int nverts = -1;
mxArray* Connectivity = NULL;
int nb_neigh_field;
int neigh_idx_field;

double* start_points = NULL;
int nstart = -1;
double* T = NULL; // weight

bool given_u;
// outputs
double *U = NULL; // distance
short *Vor = NULL; // Voronoi
short *S = NULL; // state
short *Q = NULL;
bool *doUpdate = NULL;
int *ITER = NULL;
double *TAB = NULL;
double *U_n_D = NULL;
int *heap = NULL; //heap
int *heapIndex=NULL;
double *X1, *X2, *X12;
#define faces_(k,i) faces[k+3*i]
#define vertex_(k,i) vertex[k+3*i]
double* Utmp = NULL;
short* Vtmp = NULL;
int last=1; // needed by heap
int Hilfsupdates;

//================================================================

void getNb_Neighbors(int point, double* nb_neighbors)
//================================================================
{
    *nb_neighbors = mxGetPr(mxGetFieldByNumber(Connectivity, point, nb_neigh_field))[0];
}

//================================================================

void get_Neighbors(int point, double **Neighbors)
//================================================================
{
    *Neighbors = mxGetPr(mxGetFieldByNumber(Connectivity, point, neigh_idx_field));
}

//================================================================

bool areNeighbors(int point1, int point2)
//================================================================
{
    double nb_neigh1, nb_neigh2;
    double* neigh;
    int i;
    getNb_Neighbors(point1, &nb_neigh1);
    getNb_Neighbors(point2, &nb_neigh2);
    if (nb_neigh1 > nb_neigh2) {
        get_Neighbors(point2, &neigh);
        for (i = 0; i < nb_neigh2; i++) {
            if (point1 == neigh[i]) {
                return true;
            }
        }
    } else {
        get_Neighbors(point1, &neigh);
        for (i = 0; i < nb_neigh1; i++) {
            if (point2 == neigh[i]) {
                return true;
            }
        }
    }
    return false;
}

//================================================================

void InitializeArrays()
//================================================================
{
    //--------------------------------------------------------
    int x;
    //--------------------------------------------------------
    S = new short[nverts];
    Q = new short[nverts];
    ITER = new int [nverts];
    TAB = new double[3];
    U_n_D = new double[2];
    X1 = new double[3];
    X2 = new double[3];
    X12 = new double[3];
    Utmp = new double[nverts];
    Vtmp = new short [nverts];
    //--------------------------------------------------------
    if (given_u) {
        for (x = 0; x < nverts; x++) {
            Q[x] = kUnqueued;
            ITER[x] = 0;
            if (doUpdate[x])
                S[x] = kEstimated;
            if (U[x] > 1e6)
                S[x] = kFar;
        }
    } else {
        for (x = 0; x < nverts; x++) {
            U[x] = INFINITE;
            Vor[x] = kBorder;
            Q[x] = kUnqueued;
            S[x] = kBorder;
            ITER[x] = 0;
        }
    }
    //--------------------------------------------------------
}

//================================================================

void InitializeQueue()
//================================================================
{
    int i, j, point, x;
    double nb_neigh;
    double* neigh;
    int npoint;
    short maxVor = -1;
    if (given_u) {
        for (x = 0; x < nverts; x++)
            maxVor = MAX(maxVor, Vor[x]);
        for (i = 0; i < nstart; i++) {
            point = (int) start_points[i];
            if (point >= nverts)
                mexErrMsgTxt("start_points should be in the domain.");
            //--------------------------------------------------------
            U[point] = 0.0;
            S[point] = kSeed;
            Vor[point] = i + 1 + maxVor;
            //--------------------------------------------------------
            getNb_Neighbors(point, &nb_neigh);
            get_Neighbors(point, &neigh);
            for (j = 0; j < nb_neigh; j++) {
                npoint = (int) neigh[j];
                if ((S[npoint] != kSeed)) {
                    waiting_Points.push(npoint);
                    Q[npoint] = kEnqueued;
                }
            }
            //--------------------------------------------------------
        }
    } else {
        for (i = 0; i < nstart; i++) {
            point = (int) start_points[i];
            if (point >= nverts)
                mexErrMsgTxt("start_points should be in the domain.");
            //--------------------------------------------------------
            U[point] = 0.0;
            S[point] = kSeed;
            Vor[point] = i;
            //--------------------------------------------------------
            getNb_Neighbors(point, &nb_neigh);
            get_Neighbors(point, &neigh);
            for (j = 0; j < nb_neigh; j++) {
                npoint = neigh[j];
                if ((S[npoint] != kSeed)) {
                    waiting_Points.push(npoint);
                    Q[npoint] = kEnqueued;
                }
            }
            //--------------------------------------------------------
        }
    }
}


//================================================================

void DotProductMetric(double *res, double* M, double *V1, double *V2)
/*
 COMMENTS :
 * compute V1^t M V2
 * M is symetric
 */
//================================================================
{
    *res = M[0] * V1[0] * V2[0] + M[1] * V1[1] * V2[1] + M[2] * V1[2] * V2[2]
            + M[3]*(V1[0] * V2[1] + V2[0] * V1[1])
            + M[4]*(V1[1] * V2[2] + V2[1] * V1[2])
            + M[5]*(V1[2] * V2[0] + V2[2] * V1[0]);
}

//================================================================

void TsitsiklisOnePoint(double *res, double* M, double u, double* Z)
//================================================================
{
    double NormZ_M;
    DotProductMetric(&NormZ_M, M, Z, Z);
    *res = u + sqrt(NormZ_M);
}

//================================================================

void TsitsiklisTwoPoints(double* res, double* M, double k, double u, double* z1, double* z2)
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
    R = r11 * r22 - r12*r12;
    if (!(R > 0)) {
        mexPrintf("z1 : %f %f %f\n", z1[0], z1[1], z1[2]);
        mexPrintf("z2 : %f %f %f\n", z2[0], z2[1], z2[2]);
        mexErrMsgTxt("z1 and z2 should not be collinear !!");
    }
    if (k >= sqrt(r11)) {
        res[0] = 0.0;
        res[1] = u + sqrt(r22);
    } else if (k <= -sqrt(r11)) {
        res[0] = 1.0;
        res[1] = k + u + sqrt(r11 + r22 + 2.0 * r12);
    } else {
        if (r12 >= -k * sqrt(R / (r11 - k * k))) {
            res[0] = 0.0;
            res[1] = u + sqrt(r22);
        } else if (r12 <= (-r11 - k * sqrt(R / (r11 - k * k)))) {
            res[0] = 1.0;
            res[1] = k + u + sqrt(r11 + r22 + 2.0 * r12);
        } else {
            res[0] = -(r12 + k * sqrt(R / (r11 - k * k))) / r11;
            res[1] = res[0] * k + u + sqrt(R / (r11 - k * k));
        }
    }
}

//================================================================

void TsitsiklisTriangle(double* res, int point, int npoint1, int npoint2)
//================================================================
{
    int i;
    double k, u;
    double *M, *p, *p1, *p2;
    double U1, U2, V1, V2;
    double Ur;
    bool b_point1, b_point2;
    //--------------------------------------------------------------
    M = T + 6 * point;
    Ur = U[point];
    //--------------------------------------------------------------
    res[1] = Ur + 1.0; // assign bigger value
    //----------------------------------------------------------
    b_point1 = ((S[npoint1] == kEstimated) || (S[npoint1] == kSeed)) && doUpdate[npoint1];
    b_point2 = ((S[npoint2] == kEstimated) || (S[npoint2] == kSeed)) && doUpdate[npoint2];
    if (b_point1) {
        U1 = U[npoint1];
        V1 = Vor[npoint1];
    }
    if (b_point2) {
        U2 = U[npoint2];
        V2 = Vor[npoint2];
    }
    //----------------------------------------------------------
    p = vertex + 3 * point;
    p1 = vertex + 3 * npoint1;
    p2 = vertex + 3 * npoint2;
    //----------------------------------------------------------
    for (i = 0; i < 3; i++) {
        X1[i] = p[i] - p1[i];
        X2[i] = p[i] - p2[i];
        X12[i] = p2[i] - p1[i];
    }
    //----------------------------------------------------------
    if (b_point1 && b_point2) {
        k = U1 - U2;
        u = U2;
        TsitsiklisTwoPoints(res, M, k, u, X12, X2);
        res[2] = (res[0] > 0.5 ? V1 : V2);
    }//----------------------------------------------------------
    else if (b_point1 && !b_point2) {
        u = U1;
        TsitsiklisOnePoint(res, M, u, X1);
        res[1] = res[0];
        res[0] = 1.0;
        res[2] = V1;
    }//----------------------------------------------------------
    else if (!b_point1 && b_point2) {
        u = U2;
        TsitsiklisOnePoint(res, M, u, X2);
        res[1] = res[0];
        res[0] = 0;
        res[2] = V2;
    }
}


//================================================================

void TsitsiklisUpdate(int point, double* res)
/*
COMMENTS :
 */
//================================================================
{
    int i, npoint1, npoint2;
    bool is_updated = false;
    bool doComp;
    double Ur = U[point], Vr;
    double nb_neigh;
    double* neigh;
    //--------------------------------------------------------------
    getNb_Neighbors(point, &nb_neigh);
    get_Neighbors(point, &neigh);
    //--------------------------------------------------------------
    for (i = 0; i < nb_neigh; i++) {
        if (i < nb_neigh - 1) {
            npoint1 = neigh[i];
            npoint2 = neigh[i + 1];
        } else {
            npoint1 = neigh[i];
            npoint2 = neigh[0];
        }
        doComp = areNeighbors(npoint1, npoint2);
        if (doComp) {
            //------------------------------------------------------
            TsitsiklisTriangle(TAB, point, npoint1, npoint2);
            //------------------------------------------------------
            if (TAB[1] <= Ur) {
                Ur = TAB[1];
                Vr = TAB[2];
                is_updated = true;
            }
        }
    }
    //--------------------------------------------------------------
    if (is_updated) {
        res[0] = Ur;
        res[1] = Vr;
    } else {// copy old values
        res[0] = U[point];
        res[1] = Vor[point];
    }
    //--------------------------------------------------------------
}

//================================================================

double GaussSeidelIterateVoronoi_sequential()
//================================================================
{
    int point;
    double error = 0;
    for (point = 0; point < nverts; point++) {
        TsitsiklisUpdate(point, U_n_D);
        Utmp[point] = U_n_D[0];
        Vtmp[point] = (short) U_n_D[1];
        error = MAX(error, fabs(U_n_D[0] - U[point]));
    }
    for (point = 0; point < nverts; point++) {
        U[point] = Utmp[point];
        Vor[point] = Vtmp[point];
    }
    return error;
}

void heapInit() {
    int k = 0;

    heap = new int[nverts + 2];
    heapIndex = new int[nverts + 2];
    last = OFF;
    for (k = 0; k < nverts; k++) {
        heap[k] = 0;
    }
    heap[nverts] = 0;
}

void Store(int pos, int ind) {

    heap[ind] = pos;
    heapIndex[pos] = ind;
}

void DownHeap(int ind) {
    int k;
    int pos;

    k = LEFT_SON(ind);
    pos = heap[ind];
    while (k <= last) {
        if (k < last && U[heap[k + 1]] < U[heap[k]]) {
            k++;
        }

        if (U[pos] <= U[heap[k]]) {
            break;
        }
        Store(heap[k], ind);
        ind = k;
        k = LEFT_SON(ind);
    }
    Store(pos, ind);
}

void UpHeap(int ind) {
    int k;
    int pos;

    k = FATHER(ind);
    pos = heap[ind];
    while (k >= FIRST) {
        if (U[pos] >= U[heap[k]]) {
            break;
        }
        Store(heap[k], ind);
        ind = k;
        k = FATHER(ind);
    }
    Store(pos, ind);
}

int RemoveFirstFromQueue() {
    int pos = heap[FIRST];

    heap[FIRST] = heap[last--];
    DownHeap(FIRST);

    return pos;
}

void InsertToQueue(int pos) {
    heap[++last] = pos;
    UpHeap(last);
}

void DecreaseInQueue(int pos) {
    UpHeap(heapIndex[pos]);
}

void initVoronoi() {
    int k;
    int pos;

    // delete voronoi information
    for (k = 0; k < nverts; k++) {
        Vor[k] = -1;
    }

    // set voronoi for start points
    for (k = 0; k < nstart; k++) {
        pos = (int) start_points[k];
        Vor[pos] = k;
        //        printf("Punkt %i ist Startpunkt mit dem Voronoi-Index %i\n",pos,k);
    }

    // insert all not-starting points to heap
    for (k = 0; k < nverts; k++) {
        // test if point is not a starting point
        if (Vor[k] < 0) {
            InsertToQueue(k);
        }
    }
}

int searchMinNeighbor(int pos) {

    double nb_neigh;
    double* neigh;
    double tmin = INFINITE;
    int neighborIndex;
    int minIndex = -1;
    int k;

    getNb_Neighbors(pos, &nb_neigh);
    get_Neighbors(pos, &neigh);

    //loop over all neighbors of point pos
    for (k = 0; k < nb_neigh; k++) {
        // Index of neighbor point
        neighborIndex = (int) neigh[k];

        // Falls der Nachbar eine gültige Farbe hat
        if (Vor[neighborIndex] >= 0) {
            if (U[neighborIndex] < tmin) {
                tmin = U[neighborIndex];
                minIndex = neighborIndex;
            }
        }
    }
    return minIndex;
}

void setVoronoi() {
    int pos;
    int minNeighbor;
    int k;
    double nb_neigh;
    double* neigh;
    int VoronoiTmp;
    int firstColour = -1;
    int neighborIndex;
    bool setNeighborFalse = false;

    pos = RemoveFirstFromQueue();

    TsitsiklisUpdate(pos, U_n_D);
    //U[pos] = U_n_D[0];
    VoronoiTmp = U_n_D[1];

    // Test, ob die Charakteristische Richtung ein sinnvolles Ergebnis liefert.
    if (VoronoiTmp < 0) {
        Hilfsupdates++;
        // Suche einfach den nächstgelegenen Nachbarn
        minNeighbor = searchMinNeighbor(pos);
        if (minNeighbor >= 0) {
            Vor[pos] = Vor[minNeighbor];
        } else {
            Vor[pos] = -1;
        }

    } else {
        Vor[pos] = VoronoiTmp;
    }

    getNb_Neighbors(pos, &nb_neigh);
    get_Neighbors(pos, &neigh);

    //loop over all neighbors of point pos
    for (k = 0; k < nb_neigh; k++) {
        // Index of neighbor point
        neighborIndex = (int) neigh[k];
        // look weather there is a neighbor point with valid color.
        // Falls der Nachbar eine gültige Farbe hat
        if (Vor[neighborIndex] >= 0) {
            // schaue, ob es schon einen Nachbarn mit einer gültigen Farbe gibt.
            if (firstColour < 0) {
                firstColour = Vor[neighborIndex];
            } else {
                //Falls es Nachbarn mit verschiedenen gültigen Farben gibt.
                if (firstColour != Vor[neighborIndex]) {
                    Vor[pos] = -1;
                    // alle Nachbarn ungültig setzen.
                    setNeighborFalse = true;
                    break;
                }
            }
        }
    }

    if (setNeighborFalse) {
        Vor[pos] = -1;
        //loop over all neighbors of point pos
        for (k = 0; k < nb_neigh; k++) {
            // Index of neighbor point
            neighborIndex = (int) neigh[k];
            Vor[neighborIndex] = -1;
        }
    }
}

void transportVoronoi() {
    Hilfsupdates = 0;
    while ((QUEUE_NOT_EMPTY)) {
        setVoronoi();
    }
    printf("Insgesamt %i Hilfsupdates\n", Hilfsupdates);
}

void secondVoronoiIteration() {
    int k;
    int pos;
    int minNeighbor;
    int VoronoiTmp;

    // insert all unknown points to heap
    for (k = 0; k < nverts; k++) {
        // test if point is not a starting point
        if (Vor[k] < 0) {
            InsertToQueue(k);
        }
    }

    Hilfsupdates = 0;
    while ((QUEUE_NOT_EMPTY)) {
        pos = RemoveFirstFromQueue();
        TsitsiklisUpdate(pos, U_n_D);
        VoronoiTmp = U_n_D[1];

        // Test, ob die Charakteristische Richtung ein sinnvolles Ergebnis liefert.
        if (VoronoiTmp < 0) {
            Hilfsupdates++;
            // Suche einfach den nächstgelegenen Nachbarn
            minNeighbor = searchMinNeighbor(pos);
            if (minNeighbor >= 0) {
                Vor[pos] = Vor[minNeighbor];
            } else {
                Vor[pos]= -1;
            }

        } else {
            Vor[pos] = VoronoiTmp;
        }
    }
    printf("Insgesamt %i Hilfsupdates bei der zweiten Iteration\n", Hilfsupdates);
}

//================================================================

void GaussSiedelIterate()
//================================================================
{
    int point, npoint, i;
    double nb_neigh;
    double* neigh;
    int iter = 0;
    //------------------------------------------------------------
    while (!waiting_Points.empty()) {
        point = waiting_Points.front();
        waiting_Points.pop();
        TsitsiklisUpdate(point, U_n_D);
        S[point] = kEstimated;
        Q[point] = kUnqueued;
        ITER[point] = ITER[point] + 1;
        iter = MAX(iter, ITER[point]);
        if (fabs(U_n_D[0] - U[point]) > tol) {
            //----------------------------------------------------
            U[point] = U_n_D[0];
            Vor[point] = U_n_D[1];
            //----------------------------------------------------
            getNb_Neighbors(point, &nb_neigh);
            get_Neighbors(point, &neigh);
            for (i = 0; i < nb_neigh; i++) {
                npoint = neigh[i];
                if ((S[npoint] != kSeed) && (Q[npoint] != kEnqueued)
                        && (doUpdate[npoint])) {
                    waiting_Points.push(npoint);
                    Q[npoint] = kEnqueued;
                }
            }
        }
    }
    //------------------------------------------------------------
    // iterated sequential update until convergence
    //------------------------------------------------------------
    /*double error = GaussSeidelIterateVoronoi_sequential();
    if (error > tol)
        mexPrintf("Convergence condition not reached, Iterate until then !\n");
    while (error > tol) {
        error = GaussSeidelIterateVoronoi_sequential();
    }*/
    // Initialize Heap
    heapInit();
    initVoronoi();
    printf("Voroni initialisiert\n");
    transportVoronoi();
    printf("Voroni transportiert\n");
    secondVoronoiIteration();
}
