#include <iostream>

int test(char var[])
{
    return sizeof(var);
//4, char*
}

int main(){
	char s1[] ="\0" ;
	std::cout << sizeof (s1) << "\n";
	char var[10] = {0, '1', '2'};

	std::cout << var << std::endl;
	std::cout << test(var) << std::endl;
	std::cout << strlen(var) << std::endl;
	std::cout << var[4] << std::endl;

	int i;
	float f;double d;
	char c(0);
	std::cout << sizeof(c) << std::endl;
	std::cout << sizeof(i) << std::endl;
	std::cout << sizeof(f) << std::endl;
	std::cout << sizeof(d) << std::endl;
	//10
	return 0;
}
