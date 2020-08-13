// Copyright (C) 2014-2015 Xue Zhichen, Wang Yujie,Wang Na 
// All Rights Reserved.
// This file is a part of LPOPC , published under the Eclipse Public License.
// Author:Xue Zhichen 7/12 2015   19:11
// Email:eddy_lpopc@163.com
#ifndef LPOPC_NLPSOLVER_H
#define LPOPC_NLPSOLVER_H
#include "LpNLPWrapper.hpp"
#include "LpConf.h"
#include "LpOptionList.hpp"
#include "LpDebug.hpp"
#include "LpException.hpp"
namespace Lpopc
{
	LP_DECLARE_EXCEPTION(LPOPC_NLPSOLVER_ERROR);
	class NLPSolver
	{
	public:
		NLPSolver(shared_ptr<LpOptionsList> optionlist)
			:optionlist_(optionlist)
		{};
		~NLPSolver(){};
		void SolveNlp(shared_ptr<NLPWrapper> nlpwrapper, shared_ptr<FunctionWrapper> Funs,
			shared_ptr<LpCalculateData> Data, shared_ptr<OptimalProblem> optpro);
	private:
		NLPSolver();
		NLPSolver(const NLPSolver&);
		NLPSolver& operator=(const NLPSolver&);

		shared_ptr<LpOptionsList> optionlist_;
	};
}// namespace Lpopc
#endif // !LPOPC_NLPSOLVER_H
