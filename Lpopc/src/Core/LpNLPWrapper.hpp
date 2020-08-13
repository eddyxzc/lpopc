// Copyright (C) 2014-2015 Xue Zhichen, Wang Yujie,Wang Na 
// All Rights Reserved.
// This file is a part of LPOPC , published under the Eclipse Public License.
// Author:Xue Zhichen 7/12 2015   20:15
// Email:eddy_lpopc@163.com
#ifndef LPNLPWRAPPER_HPP
#define LPNLPWRAPPER_HPP
#include "LpConf.h"
#include "LpOptimalProblem.hpp"
#include "LpFunctionWrapper.h"
#include"LpCalculateData.hpp"
#include "LpAnalyticDerive.hpp"
#include "LpHessian.h"
#include"LpException.hpp"
#include "LpOptDerive.hpp"
#include "LpSacleOCP.hpp"
#include "LpOption.hpp"
namespace Lpopc
{
	LP_DECLARE_EXCEPTION(NLPWRAPPER_ERROR);
	class NLPWrapper
	{
	public:
		NLPWrapper(shared_ptr<FunctionWrapper> Funs, 
			shared_ptr<LpCalculateData> Data,
			shared_ptr<OptimalProblem> optpro, 
			shared_ptr<OptDerive> First_deriver,
			shared_ptr<LpHessianCalculator>hessain,
			shared_ptr<LpScaleOCP> ocpScalor,
			shared_ptr<RPMGenerator> rpm
			)
			: optimalFunction_(Funs), 
			calculateData_(Data), 
			optpro_(optpro),
			derive_(First_deriver),
			hessian_(hessain),
			ocpScalor_(ocpScalor),
			rpmGenerator_(rpm)


		{
		}
		~NLPWrapper(){};

		void GetPhaseJacbi(int iphase, const vec& x_all, mat& idependencies, vec& Sjac_V, vec& Sconstant_V);
		void GetPhaseScale(size_t iphase, mat& idependencies, vec& Nconstant_scale, vec& constant_scale);
		void GetPhaseSparsity(int iphase, mat& idependencies, vec& Sjac_I, vec& Sjac_J, vec &Sconstant_I, vec &Sconstant_J);

		void GetWholeJacbi(const vec&x, vec& Sjac_V, vec& Sconstant_V);
		void GetWholeScale(vec& Nconstant_scale, vec& constant_scale);
		void GetWholeSparsity(vec& Sjac_I, vec& Sjac_J, vec &Sconstant_I, vec &Sconstant_J);

		void GetConsJacbi(const vec& x, vec& Sjac_V);
		void GetConsScale(vec& Sjac_Scale);
		void GetConsSparsity(vec& Sjac_I, vec& Sjac_J);

		void GetObjGrad(const vec& x, vec& grad_f);
		void GetObjGradScale(vec& grad_scale);

		double GetObjFun(const vec& x);
		void  GetConsFun(const vec& x, vec& Cons);//Not including linear constrains about tf-t0>0.
		void GetAllCons(vec& x, vec& Cons);//All Constraints including linear and nonlinear.

		void GetHessainValue(const vec& x, const double sigma, const vec&lambda, vec & H_V);

		void GetHessainSparsity(vec& H_I, vec& H_J);
		shared_ptr<LpCalculateData>& GetCalculateData(){ return calculateData_; }

		static void SetAllOptions(shared_ptr<LpRegisteredOptions>& registeredOptions)
		{
			registeredOptions->AddStringOption2("hessian-approximation", "How to obtain NLP Hessian ", "limited-memory", "limited-memory", "using limited-memory methods provided by Ipopt",
				"exact", "Using Sparse difference methods provided by Lpopc");
			registeredOptions->AddNumberOption("Ipopt-tol", "The tolerance used by Ipopt ", 1e-6);
			registeredOptions->AddStringOption2("auto-scale", "wether of not use autoscle method in Lpopc", "no",
				"yes", "turn on autoscale", "no", "turn off autoscale");
		}
		static int nnz(mat&);
		void GetScales();

		/************************************
		// Call this Function every time calling Ipopt.New NLP problem has new sparsity
		// Now it's called in Algorithm::GetSize();
		// Method:    RefreshSparsity
		// FullName:  Lpopc::NLPWrapper::RefreshSparsity
		// Access:    public 
		// Returns:   void
		// Qualifier: 
		//************************************/
		inline void RefreshSparsity()
		{
			jac_sparsity_need_refresh = true;
			hessain_sparsity_need_refresh = true;
		}
	private:
		
		NLPWrapper();
		NLPWrapper(const NLPWrapper&);
		NLPWrapper& operator=(const NLPWrapper&);
		

		vec conscale_;
		vec objscale_;
		shared_ptr<FunctionWrapper> optimalFunction_;
		shared_ptr<OptimalProblem> optpro_;
		shared_ptr<LpCalculateData> calculateData_;
		shared_ptr<RPMGenerator> rpmGenerator_;
		shared_ptr<LpHessianCalculator> hessian_;
		shared_ptr<OptDerive> derive_;
		shared_ptr<LpScaleOCP> ocpScalor_;
		bool jac_sparsity_need_refresh;//!<Set to true everytime call Ipopt
		bool  hessain_sparsity_need_refresh;
	};

}// lpopc
#endif // !LPNLPSOLVER_HPP
