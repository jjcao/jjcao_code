#include<iostream>
#include<math.h>
using namespace std;
double a(1),b(1),c(1);
double x,m,n;
void fun1(double *x1, double& x2)
{
	*x1=10;
    x2=(-b-sqrt(b*b-4*a*c))/(2*a);
}
void main()
{
	int a;
	a = 1;
	
	double x1,x2;
	fun1(&x1, x2);
}

