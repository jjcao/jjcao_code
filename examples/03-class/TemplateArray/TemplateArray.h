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
	template< class T1> TemplateArray(const TemplateArray<T1>& ba)// copy constructor (It is a best practice to provide copy constructor for all classes which contians members allocated dynamically)
	{
		m_size = ba.size();
		m_capacity = ba.capcity();
		m_data = new T[m_capacity];
	
		for (int i = 0; i < m_capacity; ++i)
		{
			m_data[i] = ba[i];
		}
	}
	TemplateArray(const TemplateArray& ba)// copy constructor (It is a best practice to provide copy constructor for all classes which contians members allocated dynamically)
	{
		m_size = ba.size();
		m_capacity = ba.capcity();
		m_data = new T[m_capacity];
	
		for (int i = 0; i < m_capacity; ++i)
		{
			m_data[i] = ba[i];
		}
	}
	template< class T1> TemplateArray& operator = (const TemplateArray<T1>& ba) 	// overload "=" operator
	{
		if (this == &ba)  return *this;

		delete [] m_data;
		m_size = ba.size();
		m_capacity = ba.m_capacity;
		m_data = new T[m_capacity];
	
		T *pos = ba.m_data;
		for (int i = 0; i < m_capacity; ++i,++pos)
		{
			m_data[i] = *pos;
		}
		return *this;
	}
	TemplateArray& operator = (const TemplateArray& ba) 	// overload "=" operator
	{
		if (this == &ba)  return *this;

		delete [] m_data;
		m_size = ba.size();
		m_capacity = ba.m_capacity;
		m_data = new T[m_capacity];
	
		T *pos = ba.m_data;
		for (int i = 0; i < m_capacity; ++i,++pos)
		{
			m_data[i] = *pos;
		}
		return *this;
	}
	~TemplateArray(){ 
		if ( m_data)  
			delete[] m_data;
	} // deconstructor

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
	if ( size < 1)
	{
		m_data = 0; m_capacity = 0; m_size=0;
	}
	else
	{
		m_data = new T[m_capacity];
		for (int i = 0; i < m_capacity; ++i)
		{
			m_data[i] = value;
		}	
	}
}

template<class T>
T TemplateArray<T>::at(int ind)
{
	if (isValidateIndex(ind))
		return m_data[ind];
	return -std::numeric_limits<T>::max();
}

template<class T>
T 	TemplateArray<T>::operator[] (int ind) const
{
	return m_data[ind];
}

template<class T>
T& TemplateArray<T>::operator[] (int ind)
{
	return m_data[ind];
}
template<class T>
int TemplateArray<T>::reserve(int num, T value)
{
	if ( num <= m_capacity) return m_capacity;

	m_capacity = num;
	T* newdata = new T[m_capacity];

	for ( int i = 0; i < m_size; ++i)
	{
		newdata[i] = m_data[i];
	}
	for ( int i = m_size; i < m_capacity; ++i)
	{
		newdata[i] = value;
	}

	delete [] m_data;
	m_data = newdata;
	return m_capacity;
}
template<class T>
int TemplateArray<T>::push_back(T elem)
{
	if ( m_size+1 > m_capacity )
	{
		reserve( (m_size+1)*2, elem );	
	}
	else
	{
		m_data[m_size] = elem;
	}
	++m_size;
	return m_size;
}

template<class T>
int TemplateArray<T>::insert(int ind, T value)
{
	if ( isValidateIndex(ind) )
	{
		if ( m_size+1 > m_capacity )
		{
			reserve( (m_size+1)*2 );	
		}

		for ( int i = m_size; i >= ind; --i)
		{
			m_data[i] = m_data[i-1];
		}
		m_data[ind] = value;
		++m_size;
		return m_size;
	}
	else
	{
		return -1;
	}
}

template<class T>
void TemplateArray<T>::erase(int ind)
{
	if ( isValidateIndex(ind) )
	{
		for ( int i = ind; i < m_size-1; ++i){
			m_data[i] = m_data[i+1];
		}
		--m_size;
	}
}

template<class T>
bool TemplateArray<T>::isValidateIndex(int ind)
{
	return ind>-1 && ind<m_size;
}

template<class T>
void TemplateArray<T>::print()
{
	std::cout << "[";
	for ( int i = 0; i < m_size; ++i)
	{
		std::cout << m_data[i] << ", ";
	}
	std::cout << "]" <<  std::endl;
}

}
#endif //__JJ_TEMPLATEARRAY__