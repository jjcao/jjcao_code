#include <iostream>
using namespace std;

char cn1 [] = {'a','6','.','0','9','\0'};
char* cn2 = "6.09";

void assign()
{
	cout << "cn1: " << cn1 << endl;
	cout << "cn2: " << cn2 << endl;

	cn1[0] = '7'; // ok
	//cn2[0] = '7'; // runtime error
	cn2 = "7.09"; // ok
	cout << "cn1: " << cn1 << endl;
	cout << "cn2: " << cn2 << endl;
}
void cpycat()
{
	char fragment1[] = "I'm a s";
	char fragment2[] = "string!";
	char fragment3[20];
	char finalStr[20] = "";

	strcpy(fragment3, fragment1);
	strcat(finalStr, fragment3);
	strcat(finalStr, fragment2);

	cout << finalStr << endl;
}

char* strA()
{
	char *str= "hell";
	cout << str << endl;
	return str;
}

void sizeofStr()
{
	char *str= "hellow";
	cout << strlen(str) << sizeof(str) << endl;
	char* pstr = new char[strlen(str)];	
	strcpy(pstr,str);
	cout << strlen(pstr) << sizeof(pstr) << endl;
	cout << pstr << endl;

	for(int k=0; *pstr != '\0'; ++pstr)
	{
		cout << *pstr;
	}

}
int main()
{
	sizeofStr();

	char* tmp = strA();	
	cout << *tmp << endl;
	cout << tmp << endl;


		//////////////////////////////////
std::string str("abcd"); 
const char *p = str.c_str(); 
for (int i = 0; i < 2; ++i)
{
	str.append("abcdefg");//what will happened after this
	std::cout << p << std::endl;
}

	///////////////////////////
	cpycat();

	char c = ' ';
	for (int i = 0; c != '\0'; c=cn1[i++])
	{
		if (isalpha(c))
			cout << c << endl;
		else if (isdigit(c)) //ispunct(c)
			cout << c << endl;
	}

	//////////////////////////////
	assign();

	return 0;
}