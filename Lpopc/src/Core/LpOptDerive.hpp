// Copyright (C) 2014-2015 Xue Zhichen, Wang Yujie,Wang Na 
// All Rights Reserved.
// This file is a part of LPOPC , published under the Eclipse Public License.
// Author:Xue Zhichen 7/11 2015   23:02
// Email:eddy_lpopc@163.com
#ifndef LPOPC_OPT_DERIVE_HPP
#define  LPOPC_OPT_DERIVE_HPP
#include "LpConf.h"
#include "LpOption.hpp"
#include "LpFunctionWrapper.h"
namespace Lpopc
{
	class OptDerive
	{
	public:
		OptDerive(){};
		virtual~OptDerive(){};
		virtual  void DerivMayer(SolCost&mySolcost, rowvec& deriv_mayer){};

		virtual  void DerivLagrange(SolCost& mySolcost, mat& deriv_langrange){};

		virtual void DerivDae(SolDae& mySolDae, mat& deriv_state, mat& deriv_path){}

		virtual void DerivEvent(SolEvent& mySolEvent, mat& deriv_event){}

		virtual void DerivLink(SolLink& mySolLink, mat& derive_link){}
		static void SetAllOptions(shared_ptr<LpRegisteredOptions>& registeredOptions)
		{
			registeredOptions->AddNumberOption("finite-difference-tol", "the perturbation h used in finite difference", 1e-6);
			registeredOptions->AddStringOption2("first-derive", "which first derive is used by LPOPC", "finite-difference",
				"finite-difference", "use finite-difference method of LPOPC",
				"analytic", "use user specified derive");

			registeredOptions->AddStringOption2("analytic-derive-check", "When Set to Yes ,LPOPC check user specified derive using finite difference method (consume a lot time)"
				, "no", "yes", "Check specified derive", "no", "Don't check specified derive");
			registeredOptions->AddNumberOption("analytic-derive-check-tol", "the perturbation h used in finite difference methods when check user specified derive", 1e-7);
		}
	private:
		OptDerive(const OptDerive&);
		OptDerive& operator=(const OptDerive&);
	};


}// namespace Lpopc
#endif // LPOPC_OPT_DERIVE_HPP
