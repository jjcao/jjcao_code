// movePoint.h : main header file for the movePoint DLL
//

#pragma once

#ifndef __AFXWIN_H__
	#error "include 'stdafx.h' before including this file for PCH"
#endif

#include "resource.h"		// main symbols


// CmovePointApp
// See movePoint.cpp for the implementation of this class
//

class CmovePointApp : public CWinApp
{
public:
	CmovePointApp();

// Overrides
public:
	virtual BOOL InitInstance();

	DECLARE_MESSAGE_MAP()
};
