#include <iostream>
using namespace std;

int s1(int a, char b){
	cout << "s1, int a, char b" << endl;
	return a;
}
// // diff name of arguments does not lead to an overloading function
//int s1(int b, char a){
//	cout << "s1, int b, char a" << endl;
//	return b;
//}
int s1(char b, int a){
	cout << "s1, char b, int a" << endl;
	return a;
}
char s1(char a, float b)
{
	return a;
}
// // just diff return types does not lead to an overloading
//char s1(char b, int a){
//	cout << "char, s1, char b, int a" << endl;
//	return b;
//}

int main (){
	s1(1,'1');
	s1('1',1);
}