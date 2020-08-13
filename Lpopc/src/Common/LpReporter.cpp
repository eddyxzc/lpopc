// Copyright (C) 2014-2015 Xue Zhichen, Wang Yujie,Wang Na 
// All Rights Reserved.
// This file is a part of LPOPC , published under the Eclipse Public License.
// Author:Xue Zhichen 7/10 2015   17:51
// Email:eddy_lpopc@163.com
#include <LpReporter.hpp>
#include<LpConf.h>
#include<assert.h>
namespace Lpopc
{
	


	LpReporter::LpReporter()
	{
		try{
			logger_ = spdlog::get("lpopc_main_logger");
		}
		catch (const spdlog::spdlog_ex& ex)
		{
			std::cout << "Main Logger hasn't been initialized!  " << ex.what() << std::endl;
		}
	}

	void LpReporter::Flush()
	{
		logger_->flush();
	}

}