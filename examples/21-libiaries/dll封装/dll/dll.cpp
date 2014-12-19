#include<iostream>
#include<stdlib.h>
#include<fstream>
#include<string.h>
#include<math.h>
#include"dll.h"
using namespace std;

CPolynomial::CPolynomial(string file)
{
	const char* p = file.c_str();
	ifstream fin(p);
	Node x;
	if(fin)
	{
		char a;
		int num,co,de;
		fin>>a>>num;
		for(int i=0;i<num;++i)
		{
			fin>>de>>co;
			x.cof = co;
			x.deg = de;
			AddOneTerm(x);
		}
	}
	fin.close();
}

CPolynomial::CPolynomial(double *cof,double *deg,int n)
{
	Node x;
	for(int i=0;i<n;++i)
	{
		x.deg = deg[i];
		x.cof = cof[i];
		AddOneTerm(x);
	}
}

CPolynomial::CPolynomial(vector<double> cof, vector<double> deg)
{
	Node x;
	for(int i=0;i<cof.size();++i)
	{   
		x.deg = deg[i];
		x.cof = cof[i];
		AddOneTerm(x);
	}
}

CPolynomial CPolynomial:: operator+( CPolynomial &right ) //Overload operator +
{
	CPolynomial temp;
	temp.m_Polynomial.assign(right.m_Polynomial.begin(),right.m_Polynomial.end());
	temp.m_Polynomial.merge(m_Polynomial);
	temp.m_Polynomial.sort();
	list<Node>::iterator i,j,k,p;
	for(i=temp.m_Polynomial.begin();i!=temp.m_Polynomial.end();++i)
	{	 
		p = i;
		j = ++p;
		while(j != temp.m_Polynomial.end())
		{
			if(i->deg == j->deg)
			{
				i->cof = i->cof + j->cof;
				temp.m_Polynomial.erase(j);
				p = i;
		        j = ++p;
			}
			else
				break;
		}	
	}
	for(i=temp.m_Polynomial.begin();i!=temp.m_Polynomial.end();++i)
	{
		if(i->cof == 0)
			temp.m_Polynomial.erase(i);
	}
	return temp;
}

CPolynomial CPolynomial::operator-( CPolynomial &right ) //Overload operator -
{
	CPolynomial tmp;
	tmp.m_Polynomial.assign(m_Polynomial.begin(),m_Polynomial.end());
	tmp.m_Polynomial.merge(right.m_Polynomial);
	tmp.m_Polynomial.sort();
	list<Node>::iterator i,j,k,p;
	for(i=tmp.m_Polynomial.begin();i!=tmp.m_Polynomial.end();++i)
	{	 
		p = i;
		j = ++p;
		while(j != tmp.m_Polynomial.end())
		{
			if(i->deg == j->deg)
			{
				i->cof = i->cof - j->cof;
				tmp.m_Polynomial.erase(j);
				p = i;
		        j = ++p;
			}
			else
				break;
		}	
	}
	for(i=tmp.m_Polynomial.begin();i!=tmp.m_Polynomial.end();++i)
	{
		if(i->cof == 0)
		{
			tmp.m_Polynomial.erase(i);
			i = tmp.m_Polynomial.begin();//”–Œ Ã‚£°£°£°£°
		}
	}
	return tmp;
}

CPolynomial CPolynomial::operator*( CPolynomial &right ) //Overload operator *
{
	CPolynomial tmp,tmp1;
    Node nod;
	for(list<Node>::iterator i=right.m_Polynomial.begin();i!=right.m_Polynomial.end();++i)
	{	
		for(list<Node>::iterator j=this->m_Polynomial.begin();j!=this->m_Polynomial.end();++j)
		{
			nod.cof=i->cof*j->cof;
			nod.deg=i->deg+j->deg;
			tmp1.m_Polynomial.push_back(nod);
		}
		tmp=tmp+tmp1;
		tmp1.m_Polynomial.clear();
	}
	return tmp;
}

CPolynomial CPolynomial::operator=( CPolynomial &right ) //Overload operator =
{
	if(this != &right)
	{
		this->m_Polynomial.clear();
		this->m_Polynomial.assign(right.m_Polynomial.begin(),right.m_Polynomial.end());	
	}
	return *this;
}

void CPolynomial::Print()
{
	cout<<"f(x)=";
	list<Node>::iterator i;
	for(i=m_Polynomial.begin();i!=m_Polynomial.end();++i)
	{
		if(i == m_Polynomial.begin() && i->deg ==0 && i->cof != 0)
			cout<<i->cof;
		else if(i == m_Polynomial.begin() && i->deg !=0 && i->cof != 0)
		     {
				if((abs(i->cof) - 1) > 0.0001)
					cout<<i->cof;
				if(i->deg == 1)
					cout<<"x";
				else
					cout<<"x^"<<i->deg;
		     }		     
			else
			{       
				if(i != m_Polynomial.begin() && i->cof > 0)
					cout<<"+";
				if((abs(i->cof) - 1) > 0.0001)
					cout<<i->cof;
				if(i->deg == 1)
					cout<<"x";
				else
					cout<<"x^"<<i->deg;
			}
	}
	cout<<endl;
}

void CPolynomial::ReadFromFile(string file)
{
	const char* p = file.c_str();
	ifstream fin(p);
	if(fin)
	{
		char a;
		int num;
		fin>>a>>num;
		char* pp = new char[num];
		cin.getline(pp,num);
	}
	fin.close();
}

void CPolynomial::AddOneTerm(Node term) // add one term into m_Polynomial
{
	list<Node>::iterator i;
	for(i=m_Polynomial.begin();i!=m_Polynomial.end();++i)
	{
		if(term.deg <= i->deg)
			break;
	}
	if(i == m_Polynomial.end())
		m_Polynomial.insert(i,term);
	else if(term.deg == i->deg)
		{
			i->cof = i->cof + term.cof;
			if(abs(i->cof) < 0.0001)
				m_Polynomial.erase(i);
		}
	else
		m_Polynomial.insert(i,term);
}