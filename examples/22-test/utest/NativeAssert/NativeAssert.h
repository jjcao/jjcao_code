#ifdef __NATIVE__ASSERT__H__
#include <string>

class NativeAssert
{
   public:
      static void AreEqual(int a,int b);
      static void AreEqual(std::string a,std::string b);
}

#endif //__NATIVE__ASSERT__H__