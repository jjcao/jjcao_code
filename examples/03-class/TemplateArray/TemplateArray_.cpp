#include "efficientArray.h"
using namespace jj;

EfficientArray::EfficientArray(int size, double value)
{
}

EfficientArray::EfficientArray(const EfficientArray& ba)
{
}

EfficientArray& EfficientArray::operator = (const EfficientArray& array)
{
	return *this;
}

double	EfficientArray::at(int ind)
{
	return -1.0;
}
double 	EfficientArray::operator[] (int ind) const
{
	return -1.0;
}
double& EfficientArray::operator[] (int ind)
{
	double rvalue;
	return rvalue;
}
int EfficientArray::push_back(double elem)
{
	return -1;
}

int EfficientArray::insert(int ind, double value)
{
	return -1;
}

void EfficientArray::erase(int ind)
{
}
int EfficientArray::reserve(int num, double value)
{
	return -1;
}
bool EfficientArray::isValidateIndex(int ind)
{
	return false;
}
void EfficientArray::print()
{
}

void jj::insertionSort(EfficientArray& a)
{
}