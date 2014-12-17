#include <iostream>
using namespace std;

void test1()
{
	int v[2][10] = {{1, 2, 3, 4, 5, 6, 7, 8, 9, 10},
    {11, 12, 13, 14, 15, 16, 17, 18, 19, 20}};
    int (*a)[10] = v;
    cout << **a <<endl;
    //1
    cout << **(a+1) << endl;
    //11
    cout << *(*a+1) << endl;
    //2
    cout << *(a[0]+1) << endl;
    //2
    cout << *(a[1]) << endl;
    //11
}
//void funArray(int  arr[][4], int nrow, int ncol) //ok
//void funArray(int*  arr[4], int nrow, int ncol) // error
//void funArray(int**  arr, int nrow, int ncol) //error
void funArray(int (*arr)[4], int nrow, int ncol) //ok
{
	cout << arr[1][3] << endl;
}
int main(int argc, char* argv)
{
	test1();

	//int arr[8];
	int arr[2][4]={0,1,2,3,4,5,6,7};
	cout << arr[1][3] << endl;
	arr[1][3] = 8;	
	cout << arr[1][3] << endl;
	funArray(arr, 2,4);

	cout << arr[7] << endl;
	cin >> arr[0][0];
	return 0;
}