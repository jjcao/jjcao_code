#include <iostream>
using namespace std;

int s1(int a, char b){
	cout << "s1, int a, char b" << endl;
	return a;
}
//int s1(int b, char a){
//	cout << "s1, int b, char a" << endl;
//	return b;
//}
int s1(char b, int a){
	cout << "s1, char b, int a" << endl;
	return a;
}
//char s1(char b, int a){
//	cout << "char, s1, char b, int a" << endl;
//	return b;
//}

int main (){
	s1(1,'1');
	s1('1',1);
}