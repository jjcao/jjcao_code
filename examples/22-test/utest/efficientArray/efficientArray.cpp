#include "efficientArray.h"
#include <numeric>
#include <iostream>
using namespace jj;

EfficientArray::EfficientArray(int size, double value):m_capacity(size),m_size(size)
{
	m_data = new double[m_capacity];
	for (int i = 0; i < m_capacity; ++i)
	{
		m_data[i] = value;
	}	
}

EfficientArray::EfficientArray(const EfficientArray& ba)
{
	m_size = ba.size();
	m_capacity = ba.m_capacity;
	m_data = new double[m_capacity];
	
	double *pos = ba.m_data;
	for (int i = 0; i < m_capacity; ++i,++pos)
	{
		m_data[i] = *pos;
	}
}

EfficientArray& EfficientArray::operator = (const EfficientArray& ba)
{
	if (this == &ba)  return *this;

	delete [] m_data;
	m_size = ba.size();
	m_capacity = ba.m_capacity;
	m_data = new double[m_capacity];
	
	double *pos = ba.m_data;
	for (int i = 0; i < m_capacity; ++i,++pos)
	{
		m_data[i] = *pos;
	}
	return *this;
}

double	EfficientArray::at(int ind)
{
	if (isValidateIndex(ind))
		return m_data[ind];
	return -std::numeric_limits<double>::max();
}
double 	EfficientArray::operator[] (int ind) const
{
	return m_data[ind];
}
double& EfficientArray::operator[] (int ind)
{
	return m_data[ind];
}
int EfficientArray::reserve(int num, double value)
{
	if ( num <= m_capacity) return m_capacity;

	m_capacity = num;
	double* newdata = new double[m_capacity];

	for ( int i = 0; i < m_size; ++i)
	{
		newdata[i] = m_data[i];
	}
	for ( int i = m_size; i < m_capacity; ++i)
	{
		newdata[i] = value;
	}

	delete [] m_data;
	m_data = newdata;
	return m_capacity;
}
int EfficientArray::push_back(double elem)
{
	if ( m_size+1 > m_capacity )
	{
		reserve( (m_size+1)*2, elem );	
	}
	else
	{
		m_data[m_size] = elem;
	}
	++m_size;
	return m_size;
}

int EfficientArray::insert(int ind, double value)
{
	if ( isValidateIndex(ind) )
	{
		if ( m_size+1 > m_capacity )
		{
			reserve( (m_size+1)*2 );	
		}

		for ( int i = m_size; i >= ind; --i)
		{
			m_data[i] = m_data[i-1];
		}
		m_data[ind] = value;
		++m_size;
		return m_size;
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
		for ( int i = ind; i < m_size-1; ++i){
			m_data[i] = m_data[i+1];
		}
		--m_size;
	}
}

bool EfficientArray::isValidateIndex(int ind)
{
	return ind>-1 && ind<m_size;
}
void EfficientArray::print()
{
	std::cout << "[";
	for ( int i = 0; i < m_size; ++i)
	{
		std::cout << m_data[i] << ", ";
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