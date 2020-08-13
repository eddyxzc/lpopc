// Copyright (C) 2014-2015 Xue Zhichen, Wang Yujie,Wang Na 
// All Rights Reserved.
// This file is a part of LPOPC , published under the Eclipse Public License.
// Author:Xue Zhichen 7/10 2015   17:50
// Email:eddy_lpopc@163.com
#ifndef  LPCONF_H
#define LPCONF_H

#include"armadillo"//all matrix calculation rely on armadillo
#include<memory>

#define LPOPC_VERSION "1.0"//     2015/7/12


#ifdef WIN32
#ifdef _DEBUG

#define LPOPC_REPORT_DBG_MSG ///print dbg msg
#define LPOPC_DEBUG_CHECK    ///check input args
#endif
#endif // WIN32

//#define LPOPC_REPORT_DBG_MSG /// uncomment when debug
//#define LPOPC_DEBUG_CHECK    /// uncomment when debug
#include <cassert>
namespace Lpopc
{


    typedef  arma::uword lp_index;
	using std::shared_ptr;
	using namespace arma;

}
#endif
