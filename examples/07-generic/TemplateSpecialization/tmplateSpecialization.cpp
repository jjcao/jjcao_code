#include <iostream>
#include <vector>

// fun1: compares two objects
template <typename T> int compare(const T&, const T&)
{
	std::cout << "template <typename T> int compare(const T&, const T&)" << std::endl;
	return 0;
}
// fun2: compares elements in two sequences, overloading of fun1
template <class U, class V> int compare(U, U, V)
{
	std::cout << "template <class U, class V> int compare(U, U, V)" << std::endl;
	return 0;
}
// fun3: plain functions to handle C-style character strings, overloading of fun1
int compare(const char*, const char*)
{
	std::cout << "int compare(const char*, const char*)" << std::endl;
	return 0;
}

int main()
{
	// calls compare(const T&, const T&) with T bound to int
	compare(1, 0);
	// calls compare(U, U, V), with U and V bound to vector<int>::iterator
	std::vector<int> ivec1(10), ivec2(20);
	compare(ivec1.begin(), ivec1.end(), ivec2.begin());
	int ia1[] = {0,1,2,3,4,5,6,7,8,9};
	// calls compare(U, U, V) with U bound to int* and V bound to vector<int>::iterator
	compare(ia1, ia1 + 10, ivec1.begin());
	// calls the ordinary function taking const char* parameters
	const char const_arr1[] = "world", const_arr2[] = "hi";
	compare(const_arr1, const_arr2);
	// calls the ordinary function taking const char* parameters
	char ch_arr1[] = "world", ch_arr2[] = "hi";
	compare(ch_arr1, ch_arr2);

	return 0;
}