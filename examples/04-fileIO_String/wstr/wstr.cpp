#include <string>
#include <iostream>
#include <Windows.h>
using namespace std;

std::wstring stow(const std::string& s)
{
std::wstring temp(s.length(),L' ');
std::copy(s.begin(), s.end(), temp.begin());
return temp; 
}


std::string wtos(const std::wstring& s)
{
std::string temp(s.length(), ' ');
std::copy(s.begin(), s.end(), temp.begin());
return temp; 
}

//Converting a WChar string to a Ansi string
std::string ws2s(LPCWSTR pwszSrc)
{
int nLen = WideCharToMultiByte(CP_ACP, 0, pwszSrc, -1, NULL, 0, NULL, NULL);
if (nLen<= 0) return std::string("");
char* pszDst = new char[nLen];
if (NULL == pszDst) return std::string("");
WideCharToMultiByte(CP_ACP, 0, pwszSrc, -1, pszDst, nLen, NULL, NULL);
pszDst[nLen -1] = 0;
std::string strTemp(pszDst);
delete [] pszDst;
return strTemp;
}
//Converting a Ansi string to WChar string
std::wstring s2ws(const std::string& s)
{
    int len;
    int slength = (int)s.length() + 1;
    len = MultiByteToWideChar(CP_ACP, 0, s.c_str(), slength, 0, 0); 
    wchar_t* buf = new wchar_t[len];
    MultiByteToWideChar(CP_ACP, 0, s.c_str(), slength, buf, len);
    std::wstring r(buf);
    delete[] buf;
    return r;
}

//std::wstring stemp = s2ws(myString);
//LPCWSTR result = stemp.c_str();

int main( int argc, char **argv ){
	// 1
	std::string str ="中国ab";
	std::wstring wstr =L"中国ab";
	cout << "print string: " << str << endl;
	//cout << wstr.c_str() << endl;//compile error
	//wcout << L"print wstring: " << wstr << endl;// without output

	str = wtos(wstr.c_str());
	cout << "print string: " << str << endl;

	str = ws2s(wstr.c_str());
	cout << "print string: " << str << endl;
	return 0;
}