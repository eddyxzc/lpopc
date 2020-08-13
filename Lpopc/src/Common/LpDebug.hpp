// Copyright (C) 2014-2015 Xue Zhichen, Wang Yujie,Wang Na 
// All Rights Reserved.
// This file is a part of LPOPC , published under the Eclipse Public License.
// Author:Xue Zhichen 7/10 2015   17:50
// Email:eddy_lpopc@163.com
#ifndef LPOPC_DEBUG_HPP
#define LPOPC_DEBUG_HPP
#include "LpConf.h"
#ifndef LPOPC_REPORT_DBG_MSG 
///NO Debug Msg Print
#define LP_DBG_START_FUN(__fun_name)
#define LP_DBG_START_OBJ(__fun_name)
# define LP_DBG_PRINT(__args) 
#else
/// Print debug Msg
#include "LpReporter.hpp"
#include <string>
namespace Lpopc
{
	class DebugMsgPrinter
	{
	public:
		DebugMsgPrinter(std::string fun_name);
		DebugMsgPrinter(std::string fun_name, const void* const method_owner);
		~DebugMsgPrinter();
		const shared_ptr<spdlog::logger> DbgMsgPrint()const;;
		static void SetReporter();
	private:
		static shared_ptr<spdlog::logger> dbg_logger_; 
		DebugMsgPrinter();
		DebugMsgPrinter(const DebugMsgPrinter&);
		DebugMsgPrinter& operator=(const DebugMsgPrinter&);
		std::string fun_name_;
		const void* method_owner_;
	};
#define LP_DBG_START_FUN(__fun_name)\
	DebugMsgPrinter lpopc_dbg_printer((__fun_name));

#define LP_DBG_START_OBJ(__fun_name)\
	DebugMsgPrinter lpopc_dbg_printer((__fun_name),this);

# define LP_DBG_PRINT \
  lpopc_dbg_printer.DbgMsgPrint();
}// namespace Lpopc
#endif

#endif