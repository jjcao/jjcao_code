#include"dll.h"
#include<fstream>
#include<iostream>

int main()
{
	CPolynomial p1("P3.txt");
	CPolynomial p2("P4.txt");
	CPolynomial p3;
	p1.Print();
	p2.Print();

	p3=p1+p2;
	p3.Print();
	CPolynomial p4("P3.txt");
	p3=p4-p2;
	p3.Print();

	CPolynomial p5("P3.txt");
	CPolynomial p6("P4.txt");
	p3=p5*p6;
	p3.Print();

	return 0;
}