// Copyright (C) 2014-2015 Xue Zhichen, Wang Yujie,Wang Na 
// All Rights Reserved.
// This file is a part of LPOPC , published under the Eclipse Public License.
// Author:Xue Zhichen 7/11 2015   18:09
// Email:eddy_lpopc@163.com
#ifndef LPOPC_BOUNDS_CHECKER_H
#define LPOPC_BOUNDS_CHECKER_H
#include "LpOptimalProblem.hpp"
#include "LpCalculateData.hpp"
#include "LpConf.h"
#include "LpDebug.hpp"
namespace Lpopc
{
	LP_DECLARE_EXCEPTION(LP_BOUNDSCHECKER_ERROR);
	//Check the input bounds,and get the bounds used in NLPSolver IPOPT
	class LpBoundsChecker
	{
	public:
		LpBoundsChecker(){}
		void GetBounds(shared_ptr<OptimalProblem> opt_problem, shared_ptr<LpCalculateData>& calcu_data);
		~LpBoundsChecker(){}

	private:
		LpBoundsChecker(const LpBoundsChecker&);
		LpBoundsChecker& operator = (const LpBoundsChecker&);
	};

}//end of namespace Lpopc
#endif // !LPOPC_BOUNDS_CHECKER_H
