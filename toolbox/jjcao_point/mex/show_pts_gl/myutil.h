#ifndef MY_UTIL
#define MY_UTIL

#include "float.h"     // maximum float value
#include <cmath>       // basic math functions
#include <cassert>     // to use assert()
#include <queue>       // STL priority queue
#include <vector>      // STL vector
#include <list>      // STL list, added by jjcao
#include <string>	   // C++ strings
#include <fstream>   
#include "time.h"	   // Timing facilities

using namespace std;   // default namespace throughout the application

#define PI 3.14159265

/**
 * define the infinity as the maximum float value that can be expressed
 */
#define INF FLT_MAX

/////////////////////////////////////////////////////////////////////////////
//                                                                         //
//                                                                         //
//                          OUTPUT   UTILITIES                             //
//                                                                         //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////
#ifdef __MATLAB__
	#define toSysout(...) mexPrintf(__VA_ARGS__)
    #define toSyserr(...) mexErrMsgTxt(__VA_ARGS__)
#else
     #define toSysout(...) printf(__VA_ARGS__)
     #define toSyserr(...)                  \
     do {                                   \
             fprintf(stderr, "Error: ");    \
             fprintf(stderr, __VA_ARGS__ ); \
             exit(1)                        \
     } while(0)
#endif

/////////////////////////////////////////////////////////////////////////////
//                                                                         //
//                                                                         //
//                        TIMING UTILITIES                                 //
//                                                                         //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////
time_t start_time;
void tic(){
	start_time = clock();
}
void toc(){
	float delta_t = (float) (clock() - start_time ) / CLOCKS_PER_SEC;
	toSysout("time elapsed: %f seconds\n", delta_t);
}


/////////////////////////////////////////////////////////////////////////////
//                                                                         //
//                                                                         //
//                       COMPUTATION UTILITIES                             //
//                                                                         //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////
// 3D vector sum, vectors must be preallocated
void sum( const double (&A)[3], const double (&B)[3], double (&SUM)[3] ){
	SUM[ 0 ] = A[ 0 ] + B[ 0 ];
	SUM[ 1 ] = A[ 1 ] + B[ 1 ];
	SUM[ 2 ] = A[ 2 ] + B[ 2 ];
}
// 3D vector sum, vectors must be preallocated
void sum( double (&A)[3], const double (&B)[3] ){
	A[ 0 ] += B[ 0 ];
	A[ 1 ] += B[ 1 ];
	A[ 2 ] += B[ 2 ];
}
template<class T>
void sum( vector<T>& A, const vector<T>& B ){
	A[ 0 ] += B[ 0 ];
	A[ 1 ] += B[ 1 ];
	A[ 2 ] += B[ 2 ];
}
// added by jjcao
template<class T>
void sum(const vector<T>& A, const vector<T>& B, vector<T>& SUM ){
	SUM[ 0 ] = A[ 0 ] + B[ 0 ];
	SUM[ 1 ] = A[ 1 ] + B[ 1 ];
	SUM[ 2 ] = A[ 2 ] + B[ 2 ];
}

// 3D vector scale by constant, vectors must be preallocated
template<class T>
void scale( vector<T>& A, const double alpha ){
	A[ 0 ] = A[ 0 ] * alpha;
	A[ 1 ] = A[ 1 ] * alpha;
	A[ 2 ] = A[ 2 ] * alpha;
}
template<class T>
vector<T> scale1( const vector<T>& A, const double alpha ){
	vector<T> B(3,0);
	B[ 0 ] = A[ 0 ] * alpha;
	B[ 1 ] = A[ 1 ] * alpha;
	B[ 2 ] = A[ 2 ] * alpha;
	return B;
}
// 3D vector difference, vectors must be preallocated
void diff( const vector<double>& A, const vector<double>& B, vector<double>& DIFF ){
	DIFF[ 0 ] = A[ 0 ] - B[ 0 ];
	DIFF[ 1 ] = A[ 1 ] - B[ 1 ];
	DIFF[ 2 ] = A[ 2 ] - B[ 2 ];
}
// 3D vector difference A-=B, vectors must be preallocated
void diff( vector<double>& A, const vector<double>& B ){
	A[ 0 ] = A[ 0 ] - B[ 0 ];
	A[ 1 ] = A[ 1 ] - B[ 1 ];
	A[ 2 ] = A[ 2 ] - B[ 2 ];
}
// added by jjcao
vector<int> diff( vector<int>& A, const vector<int>& B ){
	vector<int> C(2);
	C[ 0 ] = A[ 0 ] - B[ 0 ];
	C[ 1 ] = A[ 1 ] - B[ 1 ];
	return C;
}
void diff( double* A, double* B, double* DIFF ){
    DIFF[ 0 ] = A[ 0 ] - B[ 0 ];
	DIFF[ 1 ] = A[ 1 ] - B[ 1 ];
	DIFF[ 2 ] = A[ 2 ] - B[ 2 ];
}

// @DEPRECATED
double length( const double (&A)[3] ){
	return sqrt( A[0]*A[0] + A[1]*A[1] + A[2]*A[2] );
}
// vector length
double length( const vector<double> &A ){
	return sqrt( A[0]*A[0] + A[1]*A[1] + A[2]*A[2] );
}
// added by jjcao
double length( const vector<int> &A ){
	double a(A[0]), b(A[1]);
	return sqrt( a*a + b*b );
}
/**
 * @brief normlizes a vector to a versor
 *
 * @param  A      the 3D vector to normalize
 *
 * @modifies A    with the normalized vector
 */
void normalize( vector<double> &A ){
	assert( A.size()==3 );

	double len = length( A );
	A[0] /= len;
	A[1] /= len;
	A[2] /= len;
}

/// @DEPRECATED
void normalize( double (&A)[3] ){
	double len = length( A );
	A[0] /= len;
	A[1] /= len;
	A[2] /= len;
}
/**
 * @brief calculates the distance between two vertices
 *
 * @param A the first  point
 * @param B the second point
 *
 * @return the distance between A and B
 */
double euclidean_distance( const vector<double>& A, const vector<double>& B ){
    vector<double> DIFF(3,0);
    diff( A, B, DIFF );
    return length( DIFF );
}
/**
 * @brief computes the square of the length of a vector A
 *
 * @param  A	  const reference to the vector
 *
 * @return scalar containing the squared length
 */
double lengthSquared( const double (&A)[3] ){
	return A[0]*A[0] + A[1]*A[1] + A[2]*A[2];
}
/**
 * @brief computes the dot product between two vectors dot = <A,B>
 *
 * @param  A	  const reference to the first vector
 * @param  B      const reference to the second vector
 *
 * @return scalar containing the dot product of the two vectors
 */
double dot( const vector<double> &A, const vector<double> &B ){
	return ( A[0]*B[0] + A[1]*B[1] + A[2]*B[2] );
}
/**
 * @brief computes the dot product between two vectors dot = <A,B>
 *
 * @param  A	  double* of size 3 pointing to the first vector
 * @param  B      double* of size 3 pointing to the second vector
 *
 * @return scalar containing the dot product of the two vectors
 */
double dot( double* A, double* B ){
	return A[0]*B[0] + A[1]*B[1] + A[2]*B[2];
}

/**
 * @brief computes the cross product between two vectors A X B
 *
 * @param  A	  const reference to the first vector
 * @param  B      const reference to the second vector
 *
 * @return CROSS  the resulting cross product A X B
 */
void cross( const vector<double>& A, const vector<double>& B, vector<double>& CROSS ){
	 CROSS[0] = A[1] * B[2] - A[2] * B[1];
     CROSS[1] = A[2] * B[0] - A[0] * B[2];
     CROSS[2] = A[0] * B[1] - A[1] * B[0];
}
/**
 * @brief computes the cross product between two vectors A X B
 * @deprecated you should use the stl::vector version!!
 *
 * @param  A	  const reference to the first vector
 * @param  B      const reference to the second vector
 *
 * @return CROSS  the resulting cross product A X B
 */
void cross( const double (&A)[3], const double (&B)[3], double (&CROSS)[3] ){
	 CROSS[0] = A[1] * B[2] - A[2] * B[1];
     CROSS[1] = A[2] * B[0] - A[0] * B[2];
     CROSS[2] = A[0] * B[1] - A[1] * B[0];
}
/**
 * @brief computes the minimum between two elements (template type) and returns it
 *
 * @param  a	  first element to compare
 * @param  b      second element to compare
 *
 * @return min(a,b)
 */
template<class T>
T min( T& a, T& b ){
	return (a>b)?b:a;
}
/**
 * @brief computes the minimum between two elements (template type) and returns it
 *
 * @param  a	  first element to compare
 * @param  b      second element to compare
 *
 * @return min(a,b)
 */
template<class T>
T max( T& a, T& b ){
	return (b>a)?b:a;
}

/**
 * @brief computes the normal of a triangle defined by three points (A,B,C) in the space
 *
 * @note vertices are assumed to be in clockwise order
 *
 * @param  A	  const reference to first triangle vertice
 * @param  B      const reference to second triangle vertice
 * @param  C	  const reference to third triangle vertice
 *
 * @return normal  the normal of the triangle
 */
void triangle_normal( const vector<double>& A, const vector<double>& B, const vector<double>& C, vector<double>& normal ){
	vector<double> v1(3,0);
	vector<double> v2(3,0);
	diff( C, A, v1 );
	diff( B, A, v2 );
	cross( v1, v2, normal );
	normalize( normal );
}

/**
 * @brief computes the area of a triangle defined by three points (A,B,C) in the space.
 *
 * @description   to compute the area, the famous property of the cross product is used
 *
 * @note vertices are assumed to be in clockwise order
 *
 * @param  A	  const reference to first triangle vertice
 * @param  B      const reference to second triangle vertice
 * @param  C	  const reference to third triangle vertice
 *
 * @return the area of the triangle
 */
double triangle_area( const vector<double> (&A), const vector<double> (&B), const vector<double> (&C) ){
	// Area = norm( cross3D( (A-B), (A-C) ) ) / 2;
	vector<double> v1(3,0);
	vector<double> v2(3,0);
	diff( A, B, v1 );
	diff( A, C, v2 );
	vector<double> v3(3,0);
	cross( v1, v2, v3 );
	return length( v3 )/2;
}

/////////////////////////////////////////////////////////////////////////////
//                                                                         //
//                                                                         //
//                           ACCESS  UTILITIES                             //
//                                                                         //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////

/**
 * @brief retrieves a vertex from the vertex matrix
 *
 * @note no boundary checking is done
 * @note i is in C-notation
 *
 * @param i the vertex offset in the vertex list [0..Nvertices-1]
 * @param the vertex coordinate matrix
 * @param the number of vertices in the coordinate matrix
 * @param a coordinate vector where to store the retrieved vertex coordinates
 */
void getVertex( int i, double* vertices, int Nvertices, double (&coordinates)[3] ){
	assert( i<Nvertices );
	coordinates[ 0 ] = vertices[ i + 0*Nvertices ];
	coordinates[ 1 ] = vertices[ i + 1*Nvertices ];
	coordinates[ 2 ] = vertices[ i + 2*Nvertices ];
}
/**
 * @brief retrieves a face from the faces matrix
 *
 * @note no boundary checking is done
 * @note i is in C-notation [0..Nfaces]
 * @note the returned vertex indexes are in C-notation [0..Nvertices]
 */
void getFace( int i, double* faces, int Nfaces, int (&vidxs)[3] ){
	assert( i<Nfaces );
	vidxs[ 0 ] = (int)( faces[ i + 0*Nfaces ] - 1 );
	vidxs[ 1 ] = (int)( faces[ i + 1*Nfaces ] - 1 );
	vidxs[ 2 ] = (int)( faces[ i + 2*Nfaces ] - 1 );
}

#endif


