#include <iostream>

float sum_elements(float a[], unsigned length) 
 {
 	float result = 0;

 	for (int i = 0; i <= length-1; ++i)
 		result += a[i];
 	return result;
}
#define MIN(A,B) ((A) <= (B) ? (A) : (B))  

int main()
{
	int tmp = sizeof "1";

	//////////////////////////
	int r1 = MIN(3,4);
	float r2 = MIN(3.0, 4.0);  
	float a[] = {1,2,3};
	float r = sum_elements(a, 0);
	return 0;
}
