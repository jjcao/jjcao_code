#include <iostream>
#include <string>
using namespace std;

class Base{
public:
	int val;
	Base(int val=0):val(val){
		val = val;
		cout << "Base(int i)" << endl;
	}
};

class Derived: public Base{
public:	
	Derived(int i):Base(){ val = i;
    }
    Derived(char s, int i):Base(i){
    }
};

int main()
{
	Base b(3);
	Base bb=5;
    Derived *d = new Derived(10);
	 Derived dd(10);

	Base * bptr = static_cast<Base*>(d);
	cout << bptr->val << endl;
	return 0;
}