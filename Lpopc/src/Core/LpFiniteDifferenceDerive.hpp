// Copyright (C) 2014-2015 Xue Zhichen, Wang Yujie,Wang Na 
// All Rights Reserved.
// This file is a part of LPOPC , published under the Eclipse Public License.
// Author:Xue Zhichen 7/11 2015   23:55
// Email:eddy_lpopc@163.com
// Email:eddy_xuezhichen@qq.com
#ifndef LPOPC_FINITEDIFFERENCEDERIVE_HPP
#define LPOPC_FINITEDIFFERENCEDERIVE_HPP
#include "LpOptDerive.hpp"
#include "LpConf.h"
#include "LpDebug.hpp"
#include "LpFunctionWrapper.h"
#include "LpCalculateData.hpp"
#include "LpOptimalProblem.hpp"
#include "LpOption.hpp"
namespace Lpopc
{
	class LpFDderive :public OptDerive
	{
	public:
		LpFDderive(){};
		LpFDderive(const shared_ptr<FunctionWrapper>& userfun,double tol)
			:fun_(userfun)
		{
			tol_ = tol;
		};
		virtual  void DerivMayer(SolCost&mySolcost, rowvec& deriv_mayer);

		virtual  void DerivLagrange(SolCost& mySolcost, mat& deriv_langrange);

		virtual void DerivDae(SolDae& mySolDae, mat& deriv_state, mat& deriv_path);

		virtual void DerivEvent(SolEvent& mySolEvent, mat& deriv_event);

		virtual void DerivLink(SolLink& mySolLink, mat& derive_link);
		virtual ~LpFDderive(){}
		
	private:
		shared_ptr<FunctionWrapper> fun_;
		double tol_;
		LpFDderive(const LpFDderive&);
		LpFDderive& operator=(const LpFDderive&);
	};
}
#endif // !LPOPC_FINITEDIFFERENCEDERIVE_HPP
