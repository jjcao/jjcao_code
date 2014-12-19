// mShow.h : main header file for the mShow DLL
//

#pragma once

#ifndef __AFXWIN_H__
	#error "include 'stdafx.h' before including this file for PCH"
#endif

#include "resource.h"		// main symbols


// CmShowApp
// See mShow.cpp for the implementation of this class
//

class CmShowApp : public CWinApp
{
public:
	CmShowApp();

// Overrides
public:
	virtual BOOL InitInstance();

	DECLARE_MESSAGE_MAP()
};
