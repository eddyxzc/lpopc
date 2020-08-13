// Copyright (C) 2014-2015 Xue Zhichen, Wang Yujie,Wang Na 
// All Rights Reserved.
// This file is a part of LPOPC , published under the Eclipse Public License.
// Author:Xue Zhichen 7/15 2015   11:11
// Email:eddy_lpopc@163.com
#ifndef LPOPC_BRYSON_DENHAM_HPP
#define LPOPC_BRYSON_DENHAM_HPP
#include "LpFunctionWrapper.h"
using namespace Lpopc;
class BrysonDenhamFunction :public FunctionWrapper
{
public:
	BrysonDenhamFunction(){}
	virtual void MayerCost(SolCost&mySolcost, double& mayer);
	virtual  void DerivMayer(SolCost&mySolcost, rowvec& deriv_mayer);

	virtual void LagrangeCost(SolCost& mySolcost, vec& langrange);
	virtual  void DerivLagrange(SolCost& mySolcost, mat& deriv_langrange);

	virtual void DaeFunction(SolDae& mySolDae, mat& stateout, mat& pathout);
	virtual void DerivDae(SolDae& mySolDae, mat& deriv_state, mat& deriv_path);

	virtual void EventFunction(SolEvent& mySolEvent, vec& eventout);
	virtual void DerivEvent(SolEvent& mySolEvent, mat& deriv_event);

	virtual void LinkFunction(SolLink& mySolLink, vec& linkageout);
	virtual void DerivLink(SolLink& mySolLink, mat& derive_link);
	virtual ~BrysonDenhamFunction(){}
};

#endif // !LPOPC_BRYSON_DENHAM_HPP
