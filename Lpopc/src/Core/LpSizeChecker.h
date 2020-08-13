// Copyright (C) 2014-2015 Xue Zhichen, Wang Yujie,Wang Na 
// All Rights Reserved.
// This file is a part of LPOPC , published under the Eclipse Public License.
// Author:Xue Zhichen 7/11 2015   0:35
// Email:eddy_lpopc@163.com
#ifndef LPOPC_LPSIZECHECKER_H
#define LPOPC_LPSIZECHECKER_H
#include "LpConf.h"
#include "LpOptimalProblem.hpp"
#include "LpCalculateData.hpp"
namespace Lpopc
{
	LP_DECLARE_EXCEPTION(LP_SIZECHECK_ERROR);
	class LpSizeChecker
	{
		/**Get sizes of optimal control problem 
		*Exam the input date and get the number of state,control,parameter,path,event
		*number of phase,number of linkage
		*/
	public:
		LpSizeChecker(){};
		void  GetSize(shared_ptr<OptimalProblem> opt_problem,shared_ptr<LpCalculateData>& calcu_data);
		~LpSizeChecker(){};

	private:
		LpSizeChecker(const LpSizeChecker&);
		LpSizeChecker& operator = (const LpSizeChecker&);
	};
}//namespace lpopc

#endif