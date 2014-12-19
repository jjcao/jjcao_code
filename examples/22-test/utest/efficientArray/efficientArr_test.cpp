#include "efficientArray.h"
#include <boost/timer.hpp>
#include <iostream>

typedef jj::EfficientArray EArray;
#include <assert.h>

void testEfficientArray()
{
		// test constructor & at()
	EArray a(1); 
	assert( a.at(0)==0 ); 

	 // test push_back and []
	assert( a.push_back(2) == 2 );// 0 2
	assert( a[1]==2 );
	assert( a.size()==2 ); 

	// test insert
	assert( a.insert(1, 3)==3 ); 
	assert( a[1]==3 );
	a.print(); //0 3 2 

	// test erase
	a.erase(1); // 0 2 
	assert( a[1]==2 );	

	// test copy constructor and size()
	EArray b = a;  //此处用到了拷贝构造函数
	assert( b.size() ==2 );
	EArray c(a);  //该语句等同于上面的语句，都是初始化
	assert( b.size() ==2 );
	b.insert(1,1);
	assert( a.size()==2 );
	assert( b.size()==3 );
	
	// test = operator
	c = b;				//此处用到了赋值操作符号"="的重载
	assert( c.size()==3 );
	c.print();// 0 1 2 1 1 

	// test sort
	c[0] = 5; c.push_back(4); c.push_back(3);
	c.print(); // 5 1 2 4 3
	jj::insertionSort(c);
	assert( c[0]==1);
	assert( c[4]==5);
	c.print();// 1 2 3 4 5 
}
// coding begins with testing.
int main()
{
	testEfficientArray();
	return 0;
}
