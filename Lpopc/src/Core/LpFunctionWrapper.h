// Copyright (C) 2014-2015 Xue Zhichen, Wang Yujie,Wang Na 
// All Rights Reserved.
// This file is a part of LPOPC , published under the Eclipse Public License.
// Author:Xue Zhichen 7/11 2015   0:39
// Email:eddy_lpopc@163.com
#ifndef LPFUNCTIONWRAPPER_H
#define LPFUNCTIONWRAPPER_H
#include "LpConf.h"
namespace Lpopc
{
	
	struct SolCost
	{
		int phase_num_;
		double initial_time_;
		vec initial_state_;
		double terminal_time_;
		vec terminal_state_;
		vec time_;
		mat state_;
		mat control_;
		mat parameter_;
	};
	struct SolDae{
		int phase_num_;
		vec time_;
		mat state_;
		mat contol_;
		mat parameter_;
	};
	struct SolEvent
	{
		int phase_num_;
		double initial_time_;
		double terminal_time_;
		vec initial_state_;
		vec terminal_state_;
		vec parameter_;
	};
	struct SolLink
	{
		int left_phase_num_;
		int right_phase_num_;
		size_t ipair;
		vec left_state_;
		vec right_state_;
		mat left_parameter_;
		mat right_parameter_;
	};
	class FunctionWrapper
	{
	public:
		FunctionWrapper(){}
		virtual void MayerCost(SolCost&mySolcost,double& mayer){}
		virtual  void DerivMayer(SolCost&mySolcost,rowvec& deriv_mayer){};

		virtual void LagrangeCost(SolCost& mySolcost,vec& langrange){}
		virtual  void DerivLagrange(SolCost& mySolcost,mat& deriv_langrange){};

		virtual void DaeFunction(SolDae& mySolDae,mat& stateout,mat& pathout){}
		virtual void DerivDae(SolDae& mySolDae,mat& deriv_state,mat& deriv_path){}

		virtual void EventFunction(SolEvent& mySolEvent,vec& eventout){}
		virtual void DerivEvent(SolEvent& mySolEvent,mat& deriv_event){}

		virtual void LinkFunction(SolLink& mySolLink,vec& linkageout){}
		virtual void DerivLink(SolLink& mySolLink,mat& derive_link){}
		virtual ~FunctionWrapper(){}
	};
}
#endif // !LPFUNCTIONWRAPPER_H
