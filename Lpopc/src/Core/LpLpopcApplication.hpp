// Copyright (C) 2014-2015 Xue Zhichen, Wang Yujie,Wang Na 
// All Rights Reserved.
// This file is a part of LPOPC , published under the Eclipse Public License.
// Author:Xue Zhichen 7/11 2015   14:28
// Email:eddy_lpopc@163.com
#ifndef LPOPC_APPLICATION_HPP
#define  LPOPC_APPLICATION_HPP

#include "LpConf.h"
#include "LpDebug.hpp"
#include "LpException.hpp"
#include "LpReporter.hpp"
#include "LpOption.hpp"
#include "LpOptionList.hpp"

#include "LpOptimalProblem.hpp"
#include "LpCalculateData.hpp"
#include "LpLpopcAlgorithm.hpp"
#include <iostream>
namespace Lpopc
{
	LP_DECLARE_EXCEPTION(LPOPC_APPLICATION_ERROR); 
	enum console_print_not_print
	{
		console_not_print = 0,
		console_print = 1
	};
	class LpopcApplication
	{
	public:
		
		LpopcApplication(console_print_not_print if_console_print);
		void SetOptimalControlProblem(shared_ptr<OptimalProblem>& user_optimal_control_problem);
		void SolveOptimalProblem();
		shared_ptr<LpOptionsList> Options();
		~LpopcApplication();;
	private:
		LpopcApplication(const LpopcApplication&);
		LpopcApplication& operator=(const LpopcApplication&);
		void RegisterAllLpopcOptions(shared_ptr<LpRegisteredOptions>&);
		shared_ptr<LpReporter> msg_reporter_;
		shared_ptr<OptimalProblem> optpro_;
		shared_ptr<LpopcAlgorithm> algorithm_;
		shared_ptr<LpCalculateData> data_;
		shared_ptr<LpCauculateHistoryData> history_data_;

		shared_ptr<LpOptionsList> optionlist_;
		shared_ptr<LpRegisteredOptions> registeredoption_;

	};
#endif

}//namespace lpopc