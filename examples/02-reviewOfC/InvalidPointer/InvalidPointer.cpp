#include <iostream>

void GetMemory(char *p, int num)
{
 p = (char*)malloc(sizeof(char) * num);
}

void basic()
{
	int x(10);
	int& dx = x;
	int* px = &x;
	px = &dx; //error C2440: '=' : cannot convert from 'int' to 'int *'
	std::cout << *px << std::endl;
	std::cout << px[0] << std::endl;

	int* parray = new int[2];
	parray[0] = 0;
	parray[1] = 1;
	std::cout << *parray << std::endl;
	std::cout << parray[0] << std::endl;
}
void convert()
{
	double x = 2.3; 
	int *p;
	p = (int*)(&x);
	std::cout << *p << std::endl;
	std::cout << *((double*)p) << std::endl;
}

void wildPointer()
{
	int  *p;
    //*p = 3; // Run-Time Check Failure #3 - The variable 'p' is being used without being initialized.
	//std::cout << *p << std::endl;
}
void constPointer()
{
	//char* const p1;　　// p1是指针常量 // error C2734: 'p1' : const object must be initialized if not extern
	char* const p1 = "hello";
	const char* p2 = "hello"; // p2是指向常量的指针或称常量指针
	const char* const p3 = "hello"; // p3是指向常量的指针常量

	//p1 = "world";
	p2 = "world";

	//int i(2),j(3);
	//int* const p1 = &i; int* const p3 = &j;
	//const int* p2 = &i;
	//const int ci = *p2;
	////ci = 3;
	////*p2 = 3;
	//*p1 = 3;
	////p1 = p2; // cannot convert from 'const int *' to 'int *const '
	////p1 = p3;//error C3892: 'p1' : you cannot assign to a variable that is const
	//p2 = p3;
}
void reference()
{
	int x[2]={0,1};
	//int& xr = x;//error C2440: 'initializing' : cannot convert from 'int [2]' to 'int &'
	int& xr = x[0];
}
int main()
{
	char* p1 = "hello";
	basic();
	convert();
	wildPointer();
	constPointer();
	reference();

	/************** part 2 *********************/
    char *str = NULL;
	delete str;

    GetMemory(str, 100);
    strcpy(str, "hello");
    return 0;
}


//#include <iostream>
//#include <algorithm>
//
//using namespace std;
//int *myFunc(){
//	int phantom = 4;
//	cout << "address of phantom: " << &phantom << endl;
//	return &phantom;
//}
//
//int main()
//{
//	int *ptr(0);
//	//cout << "uninitialized pointer ptr: " << ptr << endl;
//	ptr = myFunc();
//	cout << "after assigning an invalid address to ptr: " << ptr << endl;
//	cout << "is the address saved the phantom 4? " << *ptr << endl;
//
//	//ptr = new int(5);
//	delete ptr;
//	//delete ptr; // do not delete a pointer more than once! runtime error!
//	cout << "is the address saved the phantom 4? " << *ptr << endl;
//	return 0;
//}
//
////typedef char* Char;
////void setString(Char& strPtr){
////strPtr="negatrive";
////}
////void setString(Char* strPtr){
////*strPtr="negatrive";
////}
////int main(){
////Char str;
////setString(&str);
////cout << str << endl;
////return 0;
////}
