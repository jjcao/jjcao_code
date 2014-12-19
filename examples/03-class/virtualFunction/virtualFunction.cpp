#include <iostream>
using namespace std;

class Vehicle
{
public:
	virtual const string getDesc()
	{
		return desc_;
	}
	int getID()
	{
		return id_;
	}
	Vehicle():desc_("vehicle"),id_(1){}
protected:
	string desc_;
	int id_;
private:
	string desc_private_;
};
class Car: public Vehicle
{
public:
	Car():Vehicle(){}
	virtual const string getDesc()
	{
		return "car";
	}
	int getID()
	{
		return ++id_;
	}
};

int main()
{
	Car c;
	Vehicle* vptr(&c);
	Vehicle &v(c);
	cout << vptr->getDesc().c_str() << endl;
	cout << v.getDesc().c_str() << endl;
	cout << vptr->getID() << endl;
	cout << v.getID() << endl;

	return 0;
}