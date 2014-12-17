#include "templateArray.h"
#include "sort.h"
#include <iostream>
#include <assert.h>

// coding begins with testing.
int main()
{
	// test constructor & at()
	jj::TemplateArray<int> a(1); 
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
	//typedef jj::TemplateArray<double, BubbleSort> DBArray;
	typedef jj::TemplateArray<double> DBArray;
	DBArray b = a;  //此处用到了拷贝构造函数
	assert( b.size() ==2 );
	DBArray c(b);  //该语句等同于上面的语句，都是初始化
	assert( b.size() ==2 );
	b.insert(1,1);
	assert( a.size()==2 );
	assert( b.size()==3 );
	
	// test = operator
	c = b;				//此处用到了赋值操作符号"="的重载
	assert( c.size()==3 );
	c.print();// 0 1 2  

	// test sort
	c[0] = 5; c.push_back(4); c.push_back(3);
	c.print(); // 5 1 2 4 3
	bubbleSort(c);
	assert( c[0]==1);	assert( c[4]==5);	
	c.print();// 1 2 3 4 5 

	return 0;
}
