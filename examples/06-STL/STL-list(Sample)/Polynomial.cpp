// Polynomial.cpp: implementation of the CPolynomial class.
//
//////////////////////////////////////////////////////////////////////

#include "Polynomial.h"
#include <fstream>
#include <iostream>
#include <algorithm>
#include <string>
using namespace std;

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////
CPolynomial::CPolynomial()
{
	m_Polynomial.clear();
}

CPolynomial::CPolynomial(string file)   
{
	ReadFromFile(file);
}

CPolynomial::CPolynomial(vector<double> cof, vector<double> deg)
{
	m_Polynomial.clear();
	int n= cof.size();
	for (int i=0;i<n;i++)
	{
		Node nd;
		nd.cof=cof[i]; nd.deg=deg[i];
		AddOneTerm(nd);
	}
}

CPolynomial::CPolynomial(double *cof,double *deg,int n)
{
	m_Polynomial.clear();
	for (int i=0;i<n;i++)
	{
		Node nd;
		nd.cof=cof[i]; nd.deg=deg[i];
		AddOneTerm(nd);
	}
}

CPolynomial::~CPolynomial()
{
	m_Polynomial.clear();
}

//overload operator + function
CPolynomial& CPolynomial::operator+(CPolynomial &right)
{
	CPolynomial *tmpPoly = new CPolynomial();
	tmpPoly->m_Polynomial=this->m_Polynomial;
	list<Node> tp;
	tp = right.m_Polynomial;
	
	list<Node>::iterator it;
	for (it=tp.begin();it!=tp.end();it++)
	{
		tmpPoly->AddOneTerm(*it);
	}
	
	return *tmpPoly;
}

//overload operator - function
CPolynomial& CPolynomial::operator-(CPolynomial &right)
{
	CPolynomial *tmpPoly = new CPolynomial();
	tmpPoly->m_Polynomial=this->m_Polynomial;
	list<Node> tp;
	tp = right.m_Polynomial;
	
	list<Node>::iterator it;
	for (it=tp.begin();it!=tp.end();it++)
	{
		(*it).cof=-(*it).cof;
		tmpPoly->AddOneTerm(*it);
	}
	
	return *tmpPoly;
}

//overload operator * function
CPolynomial& CPolynomial::operator*(CPolynomial &right)
{
	CPolynomial *tmpPoly = new CPolynomial();
	
    list<Node> p1,p2;
	p1=this->m_Polynomial;
	p2=right.m_Polynomial;


    list<Node>::iterator it1,it2;
    for (it1=p1.begin();it1!=p1.end();it1++)
		for (it2=p2.begin();it2!=p2.end();it2++)
		{
			Node tn1,tn2;
			tn1 = *it1; tn2=*it2;
			tn1.cof*=tn2.cof;tn1.deg+=tn2.deg;
			tmpPoly->AddOneTerm(tn1);
		}
	
	return *tmpPoly;
}

//overload operator = function
CPolynomial& CPolynomial::operator=(CPolynomial &right)
{
	this->m_Polynomial = right.m_Polynomial;
	return *this;
}

void CPolynomial::ReadFromFile(string file)
{
	m_Polynomial.clear();
	ifstream inp;
    inp.open(file.c_str());
	char ch;
	int n;
	inp>>ch; 
	inp>>n;
	for (int i=0;i<n;i++)
	{
		Node nd;
		inp>>nd.deg; 
		inp>>nd.cof;

        AddOneTerm(nd);
	}
	inp.close();
}

void CPolynomial::AddOneTerm(Node term)  // insert one term into the list 
{
	int n=m_Polynomial.size();
	if(n==0) { m_Polynomial.push_back(term); return;}

	if(term.cof==0) return;
	if(term.deg<0) { cout<<"Degree can not be negative!"<<endl; return; }
    
	list<Node>::iterator it;  // iterator of list
	Node td;
	for (it=m_Polynomial.begin();it!=m_Polynomial.end();it++)
	{
		td = *it;
		if(term.deg<td.deg) break;
	}
	
	// add 
	it--;
	td=*it;

	if (term.deg==td.deg)
	{
		(*it).cof+=term.cof;
		if ((*it).cof==0)
		{
			m_Polynomial.erase(it);
		}
		return;
	}
    it++;
    
	// insert
    if (it==m_Polynomial.end())
    {
		m_Polynomial.push_back(term);
		return;
    }

	m_Polynomial.insert(it,term);
}

void CPolynomial::Print()    // print the polynomial
{
	list<Node>::iterator it;
	cout<<"The polynomial is:"<<endl;

	for (it=m_Polynomial.begin();it!=m_Polynomial.end();it++)
	{
		Node td = *it;
		
        if (it!=m_Polynomial.begin())
        {
			if (td.cof>0)
			{
				cout<<"+";
			}
		}
        cout<<td.cof<<" ";
        
		if (td.deg!=0)
        {
			cout<<"x^"<<td.deg;
        }
		cout<<" ";
	}

	cout<<endl;
}