#include<vector>
#include<list>
using namespace std;

#ifdef DLL_EXPORTS
#define MY_DLL_API __declspec(dllexport)
#else               
#define MY_DLL_API __declspec(dllimport)
#endif 

MY_DLL_API void print_dynamic();

typedef struct Node
{
    double cof; // coefficient 
    int deg; // degree
} Node; // the node of polynomial

bool operator < (Node &l, Node &r)
{
	return l.deg<r.deg? 1:0;
}

class MY_DLL_API CPolynomial 
{
private:
	list<Node> m_Polynomial; 
public:
	CPolynomial(){};
	CPolynomial(string file); // initialization using file
	CPolynomial(double *cof,double *deg,int n);
	CPolynomial(vector<double> cof, vector<double> deg);
	virtual ~CPolynomial(){m_Polynomial.clear();};

	// overload
	CPolynomial operator+( CPolynomial &right ); //Overload operator +
	CPolynomial operator-( CPolynomial &right ); //Overload operator -
	CPolynomial operator*( CPolynomial &right ); //Overload operator *
	CPolynomial operator=( CPolynomial &right ); //Overload operator =

	void Print();
private:
	void ReadFromFile(string file); 
	void AddOneTerm(Node term); // add one term into m_Polynomial
};
