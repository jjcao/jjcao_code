// Polynomial.h: interface for the CPolynomial class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_POLYNOMIAL_H__DA256A82_70A5_4CCF_AD38_2278A66EFE5D__INCLUDED_)
#define AFX_POLYNOMIAL_H__DA256A82_70A5_4CCF_AD38_2278A66EFE5D__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include <list>
#include <vector>
#include <cstring>
using namespace std;

typedef struct Node
{
    double  cof;      // coefficient 
	int     deg;      // degree
} Node;               // the node of polynomial

class CPolynomial  
{
private:
    list<Node> m_Polynomial;	
public:
	CPolynomial();
	CPolynomial(string file);                     // initialization using file
	CPolynomial(double *cof,double *deg,int n);
	CPolynomial(vector<double> cof, vector<double> deg);
	virtual ~CPolynomial();

    // overload
	CPolynomial& operator+( CPolynomial &right );	//Overload operator +
	CPolynomial& operator-( CPolynomial &right );	//Overload operator -
	CPolynomial& operator*( CPolynomial &right );	//Overload operator *
	CPolynomial& operator=( CPolynomial &right );	//Overload operator =

	void Print();
private:
	void ReadFromFile(string file);  
    void AddOneTerm(Node term);   // add one term into m_Polynomial
	
};

#endif // !defined(AFX_POLYNOMIAL_H__DA256A82_70A5_4CCF_AD38_2278A66EFE5D__INCLUDED_)
