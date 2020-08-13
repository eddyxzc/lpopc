// Copyright (C) 2014-2015 Xue Zhichen, Wang Yujie,Wang Na 
// All Rights Reserved.
// This file is a part of LPOPC , published under the Eclipse Public License.
// Author:Xue Zhichen 7/10 2015   17:51
// Email:eddy_lpopc@163.com
#ifndef LPREPORTER_HPP
#define LPREPORTER_HPP
#include"LpConf.h"
#include<cstdarg>
#include<cstdio>
#include<string>
#include<vector>
#include "spdlog/spdlog.h"
namespace Lpopc
{
	/**This class modified from IPOPT's Journalist class
	* I replace the io with spdlog,which is more faster,using modern C++(C++11 needed)
	* Also,I modify the source code of spdlog.
	*/
	
	/*spdlog::level::level_enum

	typedef enum
	{
	trace    = 0,
	debug    = 1,
	info     = 2,
	notice   = 3,
	warn     = 4,
	err      = 5,
	critical = 6,
	alert    = 7,
	emerg    = 8,
	off      = 9
	} level_enum;
	*/
	typedef spdlog::level::level_enum WriteLevel;


	class LpReporter
	{
	public:
		LpReporter();

		virtual~LpReporter(){}

		//reporter.log()->info("Hello {} {} !!", "param1", 123.4);
		shared_ptr<spdlog::logger> Printf() const
		{
			return logger_;
		}
	
		void Flush();
	private:
		LpReporter(const LpReporter&);

		void operator=(const LpReporter&);

		/*std::vector<shared_ptr<LpWriter>> writers_;*/

		std::shared_ptr<spdlog::logger> logger_;
	};

	
	

}//namespace Lpopc
#endif

