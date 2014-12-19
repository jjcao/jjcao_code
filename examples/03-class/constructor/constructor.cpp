#include <iostream>
#include <string>
using namespace std;

struct USCurrency {
	int dollars;
	int cents;
	//USCurrency(int a, int b){
	//}
};
class CurrencyC{
public:
	int dollars;
	int cents;
	char a;
	//CurrencyC(){}
};

class Base{
public:
	int val;
	Base(int val):val(val){
		val = val;
		cout << "Base(int i)" << endl;
	}
};

enum suit_t {CLUBS ,
DIAMONDS ,
HEARTS ,
SPADES };

const char * print_suit ( const suit_t suit )
{
const char * names [] ={" Clubs "," Diamonds "," Hearts "," Spades "};
return names [ suit ];
}

int main(){
	print_suit(CLUBS);
	// how to initialized a class instance
	Base b(3);
	Base bb=5;

	//
	// both struct and class can be initialized by {...}, if no constructor is defined. This is not the common usage. 
	USCurrency a = {1,4}; // ok
	// USCurrency aa(1,4); // error without an proper constructor is defined
	CurrencyC c={1,4,'a'};
	return 0;
}