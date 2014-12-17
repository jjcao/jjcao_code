#include <iostream>
//在C++中，explicit关键字用来修饰类的构造函数，被修饰的构造函数的类，不能发生相应的隐式类型转换
//explicit 关键字只能用于类内部的构造函数声明上。
//explicit 关键字作用于单个参数的构造函数。

class Number
{
public:
	int type;
	Number(): type(1){};
	explicit Number(short) : type(2){};
	Number(int) : type(3) { };
};
void Show(const Number n) { 
	std::cout << n.type << std::endl;
	//printf("%d",  n.type); 
}
void main()
{
	short s = 42;	
	Show(s);// 3
}