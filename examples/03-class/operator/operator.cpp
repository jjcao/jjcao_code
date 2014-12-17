#include <iostream>
using namespace std;

class MyClass{
public:
	friend istream& operator>>(istream& in, MyClass& a);
private:
    int data;
};

	 istream&
     operator>>(istream& in, MyClass& a)
     {
		 in >> a.data;
		 return in;
	 }

int main()
{
	MyClass a;
	cin >> a;

	char s2[20] = {'A', 'N', 'S', 'I', '\0', 'C', '+', '+'};
	cout << s2[8] << endl;

	char c2[20] = {'A', 'N', 'S', 'I', '\0', 'C', '+', '+'};
	cout << c2[8] << endl;

	return 0;
}