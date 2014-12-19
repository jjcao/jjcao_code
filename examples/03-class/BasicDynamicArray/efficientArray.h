#ifndef __JJ_EFFICIENTARRAY__
#define __JJ_EFFICIENTARRAY__

namespace jj{

// interfaces of Dynamic Array class EfficientArray
class EfficientArray
{
private:
	double* data_; // the pointer to the array memory
	int		size_; // the size of the array
	int		m_capacity; // the max memory of the array

public:
	EfficientArray():data_(0),m_capacity(0),size_(0){} // default constructor
	EfficientArray(int size, double value = 0);	// other constructor, set an array with default values. Of course, more constructors can be provided.
	EfficientArray(const EfficientArray& ba); // copy constructor (It is a best practice to provide copy constructor for all classes which contians members allocated dynamically)
	EfficientArray& operator = (const EfficientArray& array);  	// overload "=" operator
	~EfficientArray(){ delete[] data_;} // deconstructor

//public:	
	int size() const { return size_;}			// get the size of the array

	double	at(int ind);				// get an element at an index
	double 	operator[] (int ind) const; 	// overload "[]" operator, get an element, such as double tmp = a[k]
	double& operator[] (int ind) ;	// overload "[]" operator, set value of specified position, such as a[k]=3.14;
	
	int push_back(double elem);		// add a new element at the end of the array, return the size of the array.
	int insert(int ind, double value);	// insert a new element at some index, return the size of the array
	void erase(int ind);			// delete an element at specified index,

	void print(); // print all elements

private:
	int reserve(int num, double value=0.0);		// 为数组增加内存空间，不能比现在的空间小，新增的元素的值设置缺省值0即可，返回m_capacity。备注：size_不变
	inline bool isValidateIndex(int ind);	// judge the validate of an index
};

void insertionSort(EfficientArray& a); // sort all elements using the Insertion sort algorithm
}
#endif //__JJ_EFFICIENTARRAY__