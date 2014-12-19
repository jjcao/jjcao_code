#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <cassert>
using namespace std;

#include <CholmodWrapper/Cholmod_conversion.h>
#include <CholmodWrapper/Sparse_coordinate_matrix.h>
#include <CholmodWrapper/Cholmod_dense_matrix.h>
typedef Sparse_coordinate_matrix<double> Sparse_matrix;
typedef Cholmod_dense_matrix<double>	 Dense_matrix;

//"factor of At", i.e. factor of A, since A is unsymmetric.
cholmod_factor* factor_At(cholmod_sparse* At, cholmod_common& c, FILE* file)
{	
	c.final_asis = 0 ;
    c.final_super = 0 ;
    c.final_ll = 0 ;
    c.final_pack = 1 ;
    c.final_monotonic = 1 ;
	c.nmethods = 1 ;
	c.method [0].ordering = CHOLMOD_NATURAL ;
	c.postorder = 0 ;

	cholmod_factor* ldlA = cholmod_analyze(At, &c);	assert(ldlA != 0);
	cholmod_factorize(At, ldlA, &c);

	cholmod_factor * tmpFactor = cholmod_copy_factor(ldlA, &c);
	cholmod_sparse* tmp = cholmod_factor_to_sparse(tmpFactor, &c);
	cholmod_write_sparse(file, tmp, 0, 0, &c);
	cholmod_free_sparse(&tmp, &c);
	cholmod_free_factor(&tmpFactor, &c);

	return ldlA;
}

////////////////////////////////
//prepare data
//D = [ 4.5 0.0 3.2 0.0
//	  3.1 2.9 0.0 0.9
//	  0.0 1.7 3.0 0.0
//	  3.5 0.4 0.0 1.0 ];
//W = [0 2 0 0
//	 0 0 0 1];
//A=[D;W];

Sparse_matrix D(4, 4);
Sparse_matrix W(2, 4);
//Sparse_matrix W(1, 4);
Sparse_matrix A(6, 4);
void prepare_data()
{
	D.add_coef(0, 0, 4.5);D.add_coef(1, 0, 3.1);D.add_coef(3, 0, 3.5);
	D.add_coef(1, 1, 2.9);D.add_coef(2, 1, 1.7);D.add_coef(3, 1, 0.4); 
	D.add_coef(0, 2, 3.2);D.add_coef(2, 2, 3.0);
	D.add_coef(1, 3, 0.9);D.add_coef(3, 3, 1.0);

	W.add_coef(0, 1, 2.0);
	W.add_coef(1, 3, 1.0);	

	A.copy(&D);
	A.add_coef(4, 1, 2.0);
	A.add_coef(5, 3, 1.0);
}

int main(int argc, char* argv[])
{
	FILE* file;
	fopen_s(&file, "log.txt","w");
	prepare_data();
	/////////////////////////////////////
	/////////////////////////
	cholmod_common common;
	cholmod_start(&common);

	cholmod_sparse* sA = create_cholmod_sparse(A, &common);	assert(sA != 0);
	cholmod_sparse* sAt = cholmod_transpose(sA, 1 /* array transpose */, &common);	assert (sAt != 0);
	fprintf( file, "A\n");
	cholmod_write_sparse(file, sA, 0, 0, &common);// A

	Dense_matrix B(6,1);
	B.set(0,0,1);B.set(1,0,2);B.set(2,0,3);B.set(3,0,4);B.set(4,0,5);B.set(5,0,6);
	cholmod_dense* b = B.create_cholmod_dense(&common);
	fprintf( file, "B\n");
	cholmod_write_dense(file, b, 0, &common);// B
	cholmod_dense* Atb = cholmod_zeros(sAt->nrow, b->ncol, CHOLMOD_REAL, &common);	
	double alpha[2] = {1, 0}, beta[2] = {0, 0};
	cholmod_sdmult(sAt, 0, alpha, beta, b, Atb, &common);
	

	fprintf( file, "\n factor of At\n");
	cholmod_factor* ldlA = factor_At(sAt, common, file);
	cholmod_dense* x = cholmod_solve(CHOLMOD_LDLt, ldlA, Atb, &common);	assert (x != 0);//CHOLMOD_LDLt//CHOLMOD_A
	fprintf( file, "\n X using factorization of A'A:\n"); 
	cholmod_write_dense(file, x, 0, &common);// x


	//////////////////////////////
	//use factorizaiton of D and W
	cholmod_sparse* sD = create_cholmod_sparse(D, &common);	assert(sD != 0);
	cholmod_sparse* sDt = cholmod_transpose(sD, 1 /* array transpose */, &common);assert (sDt != 0);
	fprintf( file, "\n factor of Dt\n");
	cholmod_factor* ldlD = factor_At(sDt, common, file);

	cholmod_sparse* sW = create_cholmod_sparse(W, &common);	assert (sW != 0);
	cholmod_sparse* sWt = cholmod_transpose(sW, 1 /* array transpose */, &common);assert (sWt != 0);
	fprintf( file, "W\n");
	cholmod_write_sparse(file, sW, 0, 0, &common);// W

	cholmod_updown(1, sWt, ldlD, &common);
	cholmod_dense* y = cholmod_solve(CHOLMOD_LDLt, ldlD, Atb, &common);	assert (y != 0);//CHOLMOD_LDLt//CHOLMOD_A
	fprintf( file, "\n X using factorization of D'D+W'W:\n"); 
	cholmod_write_dense(file, y, 0, &common);// y
	//cholmod_updown_solve(1, sWt, ldlD, y, Atb, &common);assert (y != 0);
	

	////////////////////////////////
	// free memory
	cholmod_free_sparse(&sA, &common);
	cholmod_free_sparse(&sAt, &common);
	cholmod_free_factor(&ldlA, &common);
	cholmod_free_sparse(&sD, &common);
	cholmod_free_sparse(&sDt, &common);
	cholmod_free_sparse(&sW, &common);
	cholmod_free_sparse(&sWt, &common);	
	cholmod_free_factor(&ldlD, &common);
	cholmod_finish(&common);
	fclose(file);
	return 0;
}
