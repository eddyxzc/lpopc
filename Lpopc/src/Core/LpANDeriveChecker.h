// Copyright (C) 2014-2015 Xue Zhichen, Wang Yujie,Wang Na 
// All Rights Reserved.
// This file is a part of LPOPC , published under the Eclipse Public License.
// Author:Xue Zhichen 8/4 2015   21:40
// Email:eddy_lpopc@163.com
#ifndef LPOPC_ANDERIVECHECKER_H
#define LPOPC_ANDERIVECHECKER_H
#include "LpConf.h"
#include "LpDebug.hpp"
#include "LpReporter.hpp"
#include "LpAnalyticDerive.hpp"
#include "LpFiniteDifferenceDerive.hpp"
#include "LpCalculateData.hpp"
#include "LpOptimalProblem.hpp"
namespace Lpopc
{
	///check user specified analytic derive using finite difference
	// make sure the function  is correct before you check its derive
	class LpANDeriveChecker
	{
	public:
		LpANDeriveChecker(shared_ptr<LpReporter> reporter, 
			shared_ptr<OptDerive> anderive,shared_ptr<FunctionWrapper> fun,shared_ptr<LpCalculateData>Data,
			shared_ptr<OptimalProblem> optpro,double tol)
			:reporter_(reporter), ANderive_(anderive), Data_(Data), optpro_(optpro), FDderive_(new LpFDderive(fun,tol))
		{
			tol_ = tol;
		};
		~LpANDeriveChecker(){};
		void CheckeAnlyticlDerive();
	private:
		LpANDeriveChecker(const LpANDeriveChecker&);
		LpANDeriveChecker& operator=(const LpANDeriveChecker&);
		double tol_;
		shared_ptr<LpReporter> reporter_;
		shared_ptr<OptDerive> ANderive_;
		shared_ptr<LpFDderive> FDderive_;
		shared_ptr<LpCalculateData> Data_;
		shared_ptr<OptimalProblem> optpro_;
	};
}



#endif // !LPOPC_ANDERIVECHECKER_H
