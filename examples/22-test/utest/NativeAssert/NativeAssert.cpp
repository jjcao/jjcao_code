#include <CppUnitTest.h>

using namespace System;
using namespace System::Text;
using namespace System::Collections::Generic;
using namespace	Microsoft::VisualStudio::TestTools::UnitTesting;
#include "NativeAssert.h"
void  NativeAssert::AreEqual(int a,int b)
{
    Assert::AreEqual(a,b);
}
void  NativeAssert::AreEqual(std::string a,std::string b)
{
    Assert::AreEqual<System::String^>(gcnew System::String(a.c_str()), gcnew System::String(b.c_str()));
}