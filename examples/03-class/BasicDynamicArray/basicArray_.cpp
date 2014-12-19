#include "basicArray.h"
using namespace jj;

BasicArray::BasicArray(int size, double value)
{
}

BasicArray::BasicArray(const BasicArray& ba)
{
}

BasicArray& BasicArray::operator = (const BasicArray& array)
{
	return *this;
}

double	BasicArray::at(int ind)
{
	return -1.0;
}
double 	BasicArray::operator[] (int ind) const
{
	return -1.0;
}
double& BasicArray::operator[] (int ind)
{
	double rvalue;
	return rvalue;
}
int BasicArray::push_back(double elem)
{
	return -1;
}

int BasicArray::insert(int ind, double value)
{
	return -1;
}

void BasicArray::erase(int ind)
{
}
int BasicArray::resize(int num, double elem)
{
	return -1;
}
bool BasicArray::isValidateIndex(int ind)
{
	return false;
}
void BasicArray::print()
{
}

void jj::selectionSort(BasicArray& a)
{
}