#include <vector>
#include <iostream>
#include <numeric>
#include <math.h>
using namespace std;
void main()
{
	int len(10);
	vector<int> vec;
	vec.reserve(len);
	for ( int i = 0; i < len; ++i)
		vec.push_back(len-i);

	for ( int i = 0; i < len; ++i)
		vec[i] = i;

	for ( vector<int>::iterator it = vec.begin(); it != vec.end(); ++it)
		cout << *it;
	

	return;
}