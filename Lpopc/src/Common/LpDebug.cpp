// Copyright (C) 2014-2015 Xue Zhichen, Wang Yujie,Wang Na 
// All Rights Reserved.
// This file is a part of LPOPC , published under the Eclipse Public License.
// Author:Xue Zhichen 7/10 2015   13:39
// Email:eddy_lpopc@163.com
#include "LpDebug.hpp"
#ifdef LPOPC_REPORT_DBG_MSG 
namespace Lpopc
{
	shared_ptr<spdlog::logger> DebugMsgPrinter::dbg_logger_ = shared_ptr<spdlog::logger>();
	DebugMsgPrinter::DebugMsgPrinter(std::string fun_name, const void* const method_owner)
		:fun_name_(fun_name),
		method_owner_(method_owner)
	{
		if (dbg_logger_ == NULL)
		{
			return;
		}
		DbgMsgPrint()->trace("-> Calling to: %s in obj: 0x%x\n", fun_name_.c_str(), method_owner_);
			
	}
	DebugMsgPrinter::DebugMsgPrinter(std::string fun_name)
		:fun_name_(fun_name),
		method_owner_(NULL)
	{
		if (dbg_logger_ == NULL)
		{
			return;
		}
		DbgMsgPrint()->trace("-> Calling to: {} \n", fun_name_.c_str());

	}
	DebugMsgPrinter::~DebugMsgPrinter()
	{
		if (dbg_logger_)
		{
			if (method_owner_ == NULL) {
				DbgMsgPrint()->trace("<- Returning from : {}\n", fun_name_.c_str());
			}
			else {
				DbgMsgPrint()->trace("<- Returning from : {} in obj: {}\n",
					fun_name_.c_str(), method_owner_);
			}
		}
	}

	const shared_ptr<spdlog::logger> DebugMsgPrinter::DbgMsgPrint() const
	{
		return dbg_logger_;
	}

	void DebugMsgPrinter::SetReporter()
	{
		try{
			dbg_logger_ = spdlog::get("lpopc_debug_logger");
		}
		catch (const spdlog::spdlog_ex& ex)
		{
			std::cout << "Debug Logger hasn't been initialized!  " << ex.what() << std::endl;
		}
	}


	

	

}//namespace lpopc
#endif

