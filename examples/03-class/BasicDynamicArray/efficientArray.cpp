#include "efficientArray.h"
#include <numeric>
#include <iostream>
using namespace jj;

EfficientArray::EfficientArray(int size, double value):m_capacity(size),size_(size)
{
	data_ = new double[m_capacity];
	for (int i = 0; i < m_capacity; ++i)
	{
		data_[i] = value;
	}	
}

EfficientArray::EfficientArray(const EfficientArray& ba)
{
	size_ = ba.size();
	m_capacity = ba.m_capacity;
	data_ = new double[m_capacity];
	
	double *pos = ba.data_;
	for (int i = 0; i < m_capacity; ++i,++pos)
	{
		data_[i] = *pos;
	}
}

EfficientArray& EfficientArray::operator = (const EfficientArray& ba)
{
	if (this == &ba)  return *this;

	delete [] data_;
	size_ = ba.size();
	m_capacity = ba.m_capacity;
	data_ = new double[m_capacity];
	
	double *pos = ba.data_;
	for (int i = 0; i < m_capacity; ++i,++pos)
	{
		data_[i] = *pos;
	}
	return *this;
}

double	EfficientArray::at(int ind)
{
	if (isValidateIndex(ind))
		return data_[ind];
	return -std::numeric_limits<double>::max();
}
double 	EfficientArray::operator[] (int ind) const
{
	return data_[ind];
}
double& EfficientArray::operator[] (int ind)
{
	return data_[ind];
}
int EfficientArray::reserve(int num, double value)
{
	if ( num <= m_capacity) return m_capacity;

	m_capacity = num;
	double* newdata = new double[m_capacity];

	for ( int i = 0; i < size_; ++i)
	{
		newdata[i] = data_[i];
	}
	for ( int i = size_; i < m_capacity; ++i)
	{
		newdata[i] = value;
	}

	delete [] data_;
	data_ = newdata;
	return m_capacity;
}
int EfficientArray::push_back(double elem)
{
	if ( size_+1 > m_capacity )
	{
		reserve( (size_+1)*2, elem );	
	}
	else
	{
		data_[size_] = elem;
	}
	++size_;
	return size_;
}

int EfficientArray::insert(int ind, double value)
{
	if ( isValidateIndex(ind) )
	{
		if ( size_+1 > m_capacity )
		{
			reserve( (size_+1)*2 );	
		}

		for ( int i = size_; i >= ind; --i)
		{
			data_[i] = data_[i-1];
		}
		data_[ind] = value;
		++size_;
		return size_;
	}
	else
	{
		return -1;
	}
}

void EfficientArray::erase(int ind)
{
	if ( isValidateIndex(ind) )
	{
		for ( int i = ind; i < size_-1; ++i){
			data_[i] = data_[i+1];
		}
		--size_;
	}
}

bool EfficientArray::isValidateIndex(int ind)
{
	return ind>-1 && ind<size_;
}
void EfficientArray::print()
{
	std::cout << "[";
	for ( int i = 0; i < size_; ++i)
	{
		std::cout << data_[i] << ", ";
	}
	std::cout << "]" <<  std::endl;
}

void jj::insertionSort(EfficientArray& a)
{
	double key(0.0);
	int i(0);
	for ( int j = 1; j < a.size(); ++j)
	{
		key = a[j];
		i = j -1;
		while ( i>-1 && a[i]>key)
		{
			a[i+1] = a[i];
			--i;
		}
		a[i+1] = key;
	}
}