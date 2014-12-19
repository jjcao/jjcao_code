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

//factor of A, A is symmetric.
cholmod_factor* factor_A(cholmod_sparse* A, cholmod_common& c, FILE* file)
{	
	c.final_asis = 0 ;
    c.final_super = 0 ;
    c.final_ll = 0 ;
    c.final_pack = 1 ;
    c.final_monotonic = 1 ;
	c.nmethods = 1 ;
	c.method [0].ordering = CHOLMOD_NATURAL ;
	c.postorder = 0 ;

	cholmod_factor* ldlA = cholmod_analyze(A, &c);	assert(ldlA != 0);
	cholmod_factorize(A, ldlA, &c);

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
//      0.0 1.0 0.0 0.0
//      3.2 0.0 3.0 0.9
//      0.0 0.0 0.9 1.0 ];
//W = [1 2 0 3];
Sparse_matrix D(4, 4, 1);
Sparse_matrix W(1, 4);
void prepare_data()
{
	D.add_coef(0, 0, 4.5);D.add_coef(2, 0, 3.2);
	D.add_coef(1, 1, 1.0);
	D.add_coef(0, 2, 3.2);D.add_coef(2, 2, 3.0);D.add_coef(3, 2, 0.9);
	D.add_coef(2, 3, 0.9);D.add_coef(3, 3, 1.0);
	std::cout << D <<std::endl;

	W.add_coef(0, 0, 1.0);W.add_coef(0, 1, 2.0);W.add_coef(0, 3, 3.0);
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

	//////////////////////////////
	//factor D and add row W
	cholmod_sparse* sD = create_cholmod_sparse(D, &common);	assert(sD != 0);	
	cholmod_sparse* sDt = cholmod_transpose(sD, 1 /* array transpose */, &common);	assert (sDt != 0);
	fprintf( file, "D\n");
	cholmod_write_sparse(file, sD, 0, 0, &common);// D
	fprintf( file, "Dt\n");
	cholmod_write_sparse(file, sDt, 0, 0, &common);// Dt


	//////////////////////////////////////////////////////
	////////////////
	fprintf( file, "\n factor of D\n");
	cholmod_factor* ldlD = factor_A(sD, common, file);

	cholmod_sparse* sW = create_cholmod_sparse(W, &common);	assert (sW != 0);
	cholmod_sparse* sWt = cholmod_transpose(sW, 1 /* array transpose */, &common);assert (sWt != 0);
	fprintf( file, "W\n");
	cholmod_write_sparse(file, sW, 0, 0, &common);// W
	

	cholmod_rowadd(1, sWt, ldlD, &common);
	fprintf( file, "\n add a row to factor of D \n");
	cholmod_factor * tmpFactor = cholmod_copy_factor(ldlD, &common);
	cholmod_sparse* tmp = cholmod_factor_to_sparse(tmpFactor, &common);
	cholmod_write_sparse(file, tmp, 0, 0, &common);
	cholmod_free_sparse(&tmp, &common);	
	cholmod_free_factor(&tmpFactor, &common);

	////////////////////////////////
	// free memory
	cholmod_free_sparse(&sD, &common);
	cholmod_free_sparse(&sW, &common);
	cholmod_free_sparse(&sWt, &common);	
	cholmod_free_factor(&ldlD, &common);
	cholmod_finish(&common);
	fclose(file);
	return 0;
}