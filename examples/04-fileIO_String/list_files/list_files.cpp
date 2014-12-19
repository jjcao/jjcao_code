#define _WIN32_WINNT 0x0601

#ifndef		_WINDOWS_
#include	<windows.h>
#endif

#include <string>
#include <iostream>
#include <fstream>
using namespace std;

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
//converting char[] to WChar string

//std::wstring stemp = s2ws(myString);
//LPCWSTR result = stemp.c_str();

int main( int argc, char **argv ){
	WIN32_FIND_DATA FindData;
	string strWildcard= "*.dll";	
	string a_strPath = "plugins";
	
	//verify
	string filename = a_strPath + "/test.dll";
	ifstream ifs(filename);
	if ( !ifs.is_open() ){
		cout << "open the input file: " << filename << " failed!" << endl;
		return -1;
	}
	ifs.close();

	//do	
	//::SetCurrentDirectory( LPCWSTR( a_strPath.c_str()) );
	//HANDLE hDirectory = ::FindFirstFileEx( LPCWSTR(strWildcard.c_str()), FindExInfoStandard, 	&FindData, FindExSearchNameMatch, NULL, 0 );
	wstring wa_strPath = s2ws(a_strPath);
	wstring wstrWildcard = s2ws( strWildcard);	
	::SetCurrentDirectory(  wa_strPath.c_str() );
	HANDLE hDirectory = ::FindFirstFileEx( wstrWildcard.c_str(), FindExInfoStandard, 	&FindData, FindExSearchNameMatch, NULL, 0 );

	if( hDirectory != INVALID_HANDLE_VALUE )
	{
		do{
			if( wcscmp( FindData.cFileName, L"." ) != 0 )
			{
				cout << ws2s(FindData.cFileName) << endl;
			}else
			{
				cout << ws2s(FindData.cFileName) << endl;
			}
		} while( ::FindNextFile( hDirectory, &FindData ) );
	}

	::FindClose( hDirectory );

	char szUUID[ 41 ] = "abcde";
	string tmp(szUUID);
	cout << tmp << endl;
	return 0;
}