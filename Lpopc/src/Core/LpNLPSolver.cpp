// Copyright (C) 2014-2015 Xue Zhichen, Wang Yujie,Wang Na 
// All Rights Reserved.
// This file is a part of LPOPC , published under the Eclipse Public License.
// Author:Xue Zhichen 7/12 2015   19:11
// Email:eddy_lpopc@163.com
#include "LpNLPSolver.h"
#include "IpIpoptApplication.hpp"
#include "LpopcIpopt.h"
namespace Lpopc
{


	void NLPSolver::SolveNlp(shared_ptr<NLPWrapper> nlpwrapper, shared_ptr<FunctionWrapper> Funs, 
		shared_ptr<LpCalculateData> Data, shared_ptr<OptimalProblem> optpro)
	{
		SmartPtr<TNLP> lpopc_nlp = new LpopcIpopt(nlpwrapper,Funs,Data,optpro);
		// Create a new instance of IpoptApplication
		//  (use a SmartPtr, not raw)
		// We are using the factory, since this allows us to compile this
		// example with an Ipopt Windows DLL
		SmartPtr<IpoptApplication> app = IpoptApplicationFactory();
		app->RethrowNonIpoptException(true);
		double ipopt_tol;
		std::stringstream ipopt_out_put_name;
		ipopt_out_put_name << "grid-" << Data->current_grid << "Ipopt-out.txt";
		optionlist_->GetNumericValue("Ipopt-tol", ipopt_tol, "");
		app->Options()->SetNumericValue("tol", ipopt_tol);
		app->Options()->SetStringValue("mu_strategy", "adaptive");
		app->Options()->SetStringValue("output_file", ipopt_out_put_name.str());
		//app->Options()->SetStringValue("hessian_approximation", "limited-memory");
		std::string hessian_option;
		optionlist_->GetStringValue("hessian-approximation", hessian_option, "");
		app->Options()->SetStringValue("hessian_approximation", hessian_option);
		//!<Debug Used
		//app->Options()->SetStringValue("derivative_test", "first-order");

		ApplicationReturnStatus status;
		status = app->Initialize();
		if (status != Solve_Succeeded) {
			printf("\n\n*** Error during initialization!\n");
			LP_THROW_EXCEPTION(NLPWRAPPER_ERROR, "NLP problem solve falied")
		}

		// Ask Ipopt to solve the problem
		status = app->OptimizeTNLP(lpopc_nlp);

		if (status == Solve_Succeeded) {
			printf("\n\n*** The problem solved!\n");
		}
		else {
			printf("\n\n*** The problem FAILED!\n");
			LP_THROW_EXCEPTION(LPOPC_NLPSOLVER_ERROR,"Ipopt repoters NLP prolem failed!")
		}
	}

}//namespace Lpopc