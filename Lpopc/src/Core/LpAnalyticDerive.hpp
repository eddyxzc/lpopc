// Copyright (C) 2014-2015 Xue Zhichen, Wang Yujie,Wang Na 
// All Rights Reserved.
// This file is a part of LPOPC , published under the Eclipse Public License.
// Author:Xue Zhichen 8/4 2015   21:40
// Email:eddy_lpopc@163.com
#ifndef LPOPC_ANALYTICDERIVE_HPP_
#define LPOPC_ANALYTICDERIVE_HPP_
#include "LpOptDerive.hpp"
#include "LpFunctionWrapper.h"
#include "LpConf.h"
#include "LpDebug.hpp"
namespace Lpopc
{
	
	class LpAnalyticDerive :public OptDerive
	{
	public:
		LpAnalyticDerive(shared_ptr<FunctionWrapper> fun):fun_(fun){};
		virtual  void DerivMayer(SolCost&mySolcost, rowvec& deriv_mayer)
		{
			LP_DBG_START_FUN("LpAnalyticDerive::DerivMayer()")
			fun_->DerivMayer(mySolcost, deriv_mayer);
		};

		virtual  void DerivLagrange(SolCost& mySolcost, mat& deriv_langrange)
		{
			LP_DBG_START_FUN("LpAnalyticDerive::DerivLagrange()")
			fun_->DerivLagrange(mySolcost, deriv_langrange);
		};

		virtual void DerivDae(SolDae& mySolDae, mat& deriv_state, mat& deriv_path)
		{
			LP_DBG_START_FUN("LpAnalyticDerive::DerivDae()")
			fun_->DerivDae(mySolDae, deriv_state, deriv_path);
		};

		virtual void DerivEvent(SolEvent& mySolEvent, mat& deriv_event)
		{
			LP_DBG_START_FUN("LpAnalyticDerive::DerivEvent()")
			fun_->DerivEvent(mySolEvent, deriv_event);
		};

		virtual void DerivLink(SolLink& mySolLink, mat& derive_link)
		{
			LP_DBG_START_FUN("LpAnalyticDerive::DerivLink()")
			fun_->DerivLink(mySolLink, derive_link);
		};
		virtual ~LpAnalyticDerive(){};

	private:
		shared_ptr<FunctionWrapper> fun_;
		LpAnalyticDerive(const LpAnalyticDerive&);
		LpAnalyticDerive& operator=(const LpAnalyticDerive&);
	};
}


#endif // !LPOPC_ANALYTICDERIVE_HPP_


