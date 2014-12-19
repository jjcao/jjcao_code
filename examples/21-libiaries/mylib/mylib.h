#ifndef MYLIB__H
#define MYLIB__H      

void fun_static();

class Point{
private:
	int x, y;
public:
	Point(int x=0, int y = 0):x(x),y(y){}
	int getX();
};
 
#endif