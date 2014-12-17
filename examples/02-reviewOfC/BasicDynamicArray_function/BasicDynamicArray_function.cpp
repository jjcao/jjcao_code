#include "BasicDynamicArray_function.h"
#include <iostream>
#include <assert.h>

int main()
{
	SetArraySize(2);
	assert( 2 == g_arraysize); // test-driven development (TDD)
	assert( 0 == p_arraydata[1]);

	SetValue( 0, 1.0);
	SetValue( 1, 2.0);
	assert( 2.0 == p_arraydata[1]);

	PrintArray();

	FreeArray();
	assert ( 0 == g_arraysize);
	assert ( 0 == p_arraydata);

    return 0;
}

int SetArraySize( int size)
{
	g_arraysize = size;
	p_arraydata = new double [size];
	if (NULL == p_arraydata) // Defensive Programming
	{
		std::cout << " no enough memory" << std::endl;
		return 0;
	}
	return 1;
}
int FreeArray()
{
	if (NULL != p_arraydata) // Defensive Programming
	{
		delete [] p_arraydata;
		p_arraydata = NULL;
		//g_arraysize = 0;
		return 0;
	}
	return -1;
}
int SetValue( int k, double value)
{
	if (NULL == p_arraydata) // Defensive Programming
		return 0;

	if ( k < 0 || k >= g_arraysize) // Defensive Programming
		return 0;

	p_arraydata[k] = value;
	return k;
}
void PrintArray()
{
	std::cout << "[";
	for (int i = 0; i < g_arraysize; ++i)
	{
		std::cout << p_arraydata[i] << ", ";
	}
	std::cout << std::endl;
}