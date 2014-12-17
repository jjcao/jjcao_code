#include "stdafx.h"
#include <memory>
//#include "efficientArray.h"
//typedef jj::EfficientArray EArray;

using namespace System;
using namespace System::Text;
using namespace System::Collections::Generic;
using namespace	Microsoft::VisualStudio::TestTools::UnitTesting;

namespace UnitTest
{
	[TestClass]
	public ref class UnitTest
	{
	public: 
		[TestMethod]
		void TestMethod1()
		{
			std::unique_ptr<jj::EfficientArray> pClass(new jj::EfficientArray());
			Assert::AreEqual<int>(1, pClass->size());
		}
	};
}
