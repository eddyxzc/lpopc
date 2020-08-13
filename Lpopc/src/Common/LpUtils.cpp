// Copyright (C) 2014-2015 Xue Zhichen, Wang Yujie,Wang Na 
// All Rights Reserved.
// This file is a part of LPOPC , published under the Eclipse Public License.
// Author:Xue Zhichen 7/10 2015   17:50
// Email:eddy_lpopc@163.com
#include "LpUtils.h"
#include<cstdarg>
#include<cstdio>
int Snprintf( char* str, long size, const char* format, ... )
{
	va_list ap;
	va_start(ap, format);
	int ret;\
	ret = vsprintf(str, format, ap);
	va_end(ap);
	return ret;
}

