#include <string>
#include <iostream>
#include <assert.h>
using namespace std;

void init_text(string open_text, string& save_text, string& open_scalarField_text, string& save_scalarField_text, string& optimize_scalarField_text, string& optimize_debug_text)
{
	string::size_type idx = open_text.find_last_of(".smf");
	string basename = open_text.substr(0, idx-3);
	save_text = basename + "_mc.smf";
	open_scalarField_text = basename + "_sf_in.txt";
	save_scalarField_text = basename + "_sf_out.txt";
	optimize_scalarField_text = basename + "_sf_opt.txt";
	optimize_debug_text = basename + "_sf_deb.txt";
}
int main( int argc, char **argv )
{
	string  open_text("./mesh/sphere_4k.smf");
	string sace_text, open_scalarField_text, save_scalarField_text, optimize_scalarField_text, optimize_debug_text;
	/*char  save_text[] = "./mesh/sphere_4k_mc.smf";
	char  open_scalarField_text[] = "./mesh/shpere_4k_sf_in.txt";
	char  save_scalarField_text[] = "./mesh/sphere_4k_sf_out.txt";
	char  optimize_scalarField_text[] = "./mesh/sphere_4k_sf_opt.txt";
	char  optimize_debug_text[] = "./mesh/sphere_4k_sf_deb.txt";*/
	init_text(open_text, sace_text, open_scalarField_text, save_scalarField_text, optimize_scalarField_text, optimize_debug_text);

	cout << "open_text = " << open_text <<  endl;
	cout << "sace_text = " << sace_text << endl;
	cout << "open_scalarField_text = " << open_scalarField_text << endl;
	cout << "save_scalarField_text = " << save_scalarField_text << endl;
	cout << "open_optimize_scalarField_text = " << optimize_scalarField_text << endl;
	cout << "optimize_debug_text = " << optimize_debug_text << endl;

	assert( sace_text.compare("./mesh/sphere_4k_mc.smf")==0 );
	assert( open_scalarField_text.compare("./mesh/sphere_4k_sf_in.txt")==0 );
	assert( save_scalarField_text.compare("./mesh/sphere_4k_sf_out.txt")==0 );
	assert( optimize_scalarField_text.compare("./mesh/sphere_4k_sf_opt.txt")==0 );
	assert( optimize_debug_text.compare("./mesh/sphere_4k_sf_deb.txt")==0 );

	return 0;
}