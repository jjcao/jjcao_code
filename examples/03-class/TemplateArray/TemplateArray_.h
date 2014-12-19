#ifndef __JJ_TEMPLATEARRAY__
#define __JJ_TEMPLATEARRAY__

namespace jj{

// interfaces of Dynamic Array class TemplateArray
template<class T>
class TemplateArray
{
private:
	T* m_data; // the pointer to the array memory
	int m_size; // the size of the array
	int m_capacity; // the max memory of the array

public:
	typedef T Value_type;
	TemplateArray():m_data(0),m_capacity(0),m_size(0){} // default constructor
	TemplateArray(int size, T value = 0);	// other constructor, set an array with default values. Of course, more constructors can be provided.
	TemplateArray(const TemplateArray& ba)// copy constructor (It is a best practice to provide copy constructor for all classes which contians members allocated dynamically)
	{
	}
	template< class T1> TemplateArray(const TemplateArray<T1>& ba)// copy constructor (It is a best practice to provide copy constructor for all classes which contians members allocated dynamically)
	{
	}
	TemplateArray& operator = (const TemplateArray& ba) 	// overload "=" operator
	{
		return *this;
	}
	template< class T1> TemplateArray& operator = (const TemplateArray<T1>& ba) 	// overload "=" operator
	{
		return *this;
	}
	~TemplateArray(){ delete[] m_data;} // deconstructor

public:	
	int size() const { return m_size;}			// get the size of the array
	int capcity() const {return m_capacity;}

	T	at(int ind);				// get an element at an index
	T 	operator[] (int ind) const; 	// overload "[]" operator, get an element, such as T tmp = a[k]
	T& operator[] (int ind) ;	// overload "[]" operator, set value of specified position, such as a[k]=3.14;
	
	int push_back(T elem);		// add a new element at the end of the array, return the size of the array.
	int insert(int ind, T value);	// insert a new element at some index, return the size of the array
	void erase(int ind);			// delete an element at specified index,

	void print(); // print all elements

private:
	int reserve(int num, T value=0.0);		// 为数组增加内存空间，不能比现在的空间小，新增的元素的值设置缺省值0即可，返回m_capacity。备注：m_size不变
	inline bool isValidateIndex(int ind);	// judge the validate of an index
};

template<class T>
TemplateArray<T>::TemplateArray(int size, T value):m_capacity(size),m_size(size)
{
}

template<class T>
T TemplateArray<T>::at(int ind)
{
	T result;
	return result;
}

template<class T>
T 	TemplateArray<T>::operator[] (int ind) const
{
	T result;
	return result;
}

template<class T>
T& TemplateArray<T>::operator[] (int ind)
{
	T result;
	return result;
}
template<class T>
int TemplateArray<T>::reserve(int num, T value)
{
	return -1;
}
template<class T>
int TemplateArray<T>::push_back(T elem)
{
	return -1;
}

template<class T>
int TemplateArray<T>::insert(int ind, T value)
{
	return -1;
}

template<class T>
void TemplateArray<T>::erase(int ind)
{
}

template<class T>
bool TemplateArray<T>::isValidateIndex(int ind)
{
	return 0;
}

template<class T>
void TemplateArray<T>::print()
{
}

}

#endif //__JJ_TEMPLATEARRAY__