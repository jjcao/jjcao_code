#ifndef __JJ_BASICARRAY__
#define __JJ_BASICARRAY__

namespace jj{

// interfaces of Dynamic Array class BasicArray
class BasicArray
{
private:
	double* m_data; // the pointer to the array memory
	int		m_size; // the size of the array
	//int		m_capacity; // the max memory of the array

public:
	BasicArray():m_data(0),m_size(0){} // default constructor
	BasicArray(int size, double value = 0);	// other constructor, set an array with default values. Of course, more constructors can be provided.
	BasicArray(const BasicArray& ba); // copy constructor (It is a best practice to provide copy constructor for all classes which contians members allocated dynamically)
	BasicArray& operator = (const BasicArray& array);  	// overload "=" operator
	~BasicArray(){ delete[] m_data;} // deconstructor

//public:	
	int size() const { return m_size;}			// get the size of the array
	int resize(int num, double elem=0.0);		// re-set the size of the array. 注：若m_size小于原数组大小，可截断取前m_size个元素作为新数组的元素；若m_size大于原数组大小，新增的元素的值设置缺省值0即可

	double	at(int ind);				// get an element at an index
	double 	operator[] (int ind) const; 	// overload "[]" operator, get an element, such as double tmp = a[k]
	double& operator[] (int ind) ;	// overload "[]" operator, set value of specified position, such as a[k]=3.14;
	
	int push_back(double elem);		// add a new element at the end of the array, return the size of the array.
	int insert(int ind, double value);	// insert a new element at some index, return the size of the array
	void erase(int ind);			// delete an element at specified index

	void print(); // print all elements

private:
	inline bool isValidateIndex(int ind);	// judge the validate of an index
};

void selectionSort(BasicArray& a); // sort all elements using the selection sort algorithm
}
#endif //__JJ_BASICARRAY__