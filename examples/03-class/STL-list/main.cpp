#include "Polynomial.h"
#include <list>
#include <iostream>

using namespace std;


int main()
{
    CPolynomial p1("P3.txt");
	CPolynomial p2("P4.txt");
    CPolynomial p3;
	p1.Print();
	p2.Print();

	p3=p1+p2;
	p3.Print();
    p3=p1-p2;
	p3.Print();

	p3=p1*p2;
	p3.Print();

	return 0;
}