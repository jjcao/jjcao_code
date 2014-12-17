#ifndef __JJ_SORT__
#define __JJ_SORT__

namespace jj{
	template<class T>
	void selectionSort(T& a)
	{
		int tmp(0);
		for ( int i = 0; i < a.size(); ++i)
		{
			int k = i;
			for (int j = i+1; j < a.size(); ++j)
			{
				if ( a[j] < a[k] )
				{
					k = j;
				}
			}
			tmp = a[k];
			a[k] = a[i];
			a[i] = tmp;
		}
	}

	template<class T>
	void insertionSort(T& a)
	{
		T::Value_type key(0.0);
		int i(0);
		for ( int j = 1; j < a.size(); ++j)
		{
			key = a[j];
			i = j -1;
			while ( i>-1 && a[i]>key)
			{
				a[i+1] = a[i];
				--i;
			}
			a[i+1] = key;
		}
	}

	template<class T>
	void bubbleSort(T& a)
	{		
		for(int i=0;i<a.size();i++)
		{
			int flag(1);
			for(int j=0;j<a.size()-i-1;++j)
			{
				if(a[j]>a[j+1])
				{
					swap(a[j],a[j+1]);
					flag=0;
				}
			}
			if(flag)break;
		}	
	}
} // end of namespace jj

#endif //__JJ_SORT__