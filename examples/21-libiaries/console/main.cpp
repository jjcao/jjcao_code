#include "mylib.h"
#include "mydll.h"

#include <iostream>
void main()
{
	// use a class in a static libary
	Point pt_static;
	pt_static.getX();
	// call a function in a static library
	fun_static();

	// use a class in a static libary
	PointArray pts_dynamic;
	pts_dynamic.push_back( pt_static);
	// call a function in a dynamic library
	print_dynamic();

	std::cout << "main" << std::endl;
}