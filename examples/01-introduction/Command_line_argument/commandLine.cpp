#include <iostream>
using namespace std;
int main(int argc, char** args){
	// Notice I start from i=1 not 0 because the args[0] is reserved 
    //        for the name of this program.
	for(int i = 1; i < argc; i++)
	{
	     cerr << i << "th argument is " << args[i] << "\n";
	}
	return 0;
}
