#ifdef MYDLL_EXPORTS
#define MY_DLL_API __declspec(dllexport)
#else               
#define MY_DLL_API __declspec(dllimport)
#endif 

#include "mylib.h"
#include <vector>

MY_DLL_API void print_dynamic();
class MY_DLL_API PointArray
{
private:
	std::vector<Point> pts;
	int size;
public:
	PointArray(){}
	void push_back(const Point &p);
};