#include "mydll.h"
#include <iostream>
void print_dynamic(){
	std::cout << "tt" << std::endl;

}

void PointArray::push_back(const Point &p){
	pts.push_back(p);
}