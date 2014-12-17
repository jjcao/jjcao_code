#include "basicArray.h"
#include <numeric>
#include <iostream>
using namespace jj;

BasicArray::BasicArray(int size, double value):size_(size)
{
	data_ = new double[size_];
	for (int i = 0; i < size_; ++i)
	{
		data_[i] = value;
	}
}

BasicArray::BasicArray(const BasicArray& ba)
{
	size_ = ba.size();
	data_ = new double[size_];
	
	double *pos = ba.data_;
	for (int i = 0; i < size_; ++i,++pos)
	{
		data_[i] = *pos;
	}
}

BasicArray& BasicArray::operator = (const BasicArray& ba)
{
	if (this == &ba)  return *this;

	if(data_) delete [] data_;

	size_ = ba.size();
	data_ = new double[size_];
	
	double *pos = ba.data_;
	for (int i = 0; i < size_; ++i,++pos)
	{
		data_[i] = *pos;
	}
	return *this;
}

double	BasicArray::at(int ind)
{
	if (isValidateIndex(ind))
		return data_[ind];
	return -std::numeric_limits<double>::max();
}
double 	BasicArray::operator[] (int ind) const
{
	return data_[ind];
}
double& BasicArray::operator[] (int ind)
{
	return data_[ind];
}
int BasicArray::resize(int num, double elem)
{
	double* newdata = new double[num];

	if ( num < size_)
	{
		for ( int i = 0; i < num; ++i)
		{
			newdata[i] = data_[i];
		}
	}
	else
	{
		for ( int i = 0; i < size_; ++i)
		{
			newdata[i] = data_[i];
		}
		for ( int i = size_; i < num; ++i)
		{
			newdata[i] = elem;
		}
	}

	delete [] data_;
	data_ = newdata;
	size_ = num;
	return num;
}
int BasicArray::push_back(double elem)
{
	return resize(size_+1, elem);
}

int BasicArray::insert(int ind, double value)
{
	if ( isValidateIndex(ind) )
	{
		resize(size_+1);
		for ( int i = size_ -1; i > ind; --i)
		{
			data_[i] = data_[i-1];
		}

		data_[ind] = value;
		return size();
	}
	else
	{
		return -1;
	}
}

void BasicArray::erase(int ind)
{
	if ( isValidateIndex(ind) )
	{
		for ( int i = ind; i < size_-1; ++i){
			data_[i] = data_[i+1];
		}
		resize(size_-1);
	}
}

bool BasicArray::isValidateIndex(int ind)
{
	return ind>-1 && ind<size_;
}
std::ostream& jj::operator <<(std::ostream& os, const BasicArray &ba)
{
	os << "[";
	for ( int i = 0; i < ba.size_; ++i)
	{
		os << ba.data_[i] << ", ";
	}
	os << "]" <<  std::endl;

	return os;
}

//void jj::selectionSort(BasicArray& a){
//	int tmp(0);
//	for ( int i = 0; i < a.size(); ++i)	{
//		int k = i;
//		for (int j = i+1; j < a.size(); ++j){
//			if ( a[j] < a[k] ){
//				k = j;
//			}
//		}
//		tmp = a[k];
//		a[k] = a[i];
//		a[i] = tmp;
//	}
//}
void jj::selectionSort(BasicArray& a)
{
	int i,j,k;
	for(i=0;i<a.size();++i)
	{
		k=i;
		for(j=i+1;j<a.size();++j)
	    {
			if(a[j]<a[k])
				k=j;
		}
		//swap(a[i],a[k]);
		if ( i != k)
		{
			int tmp = a[k];
			a[k] = a[i];
			a[i] = tmp;
		}
		std::cout << a;
	}
}

//void BasicArray::bubbleSort()
//{
//	double tmp(0.0);
//	for (int i = 0; i < size_; ++i)
//	{
//		for (int j = size_-1; j>i; --j)
//		{
//			if ( data_[j] < data_[j-1])
//			{
//				tmp = data_[j-1];
//				data_[j-1] = data_[j];
//				data_[j] = tmp;
//			}
//		}
//	}
//}
