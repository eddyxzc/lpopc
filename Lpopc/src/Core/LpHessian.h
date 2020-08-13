// Copyright (C) 2014-2015 Xue Zhichen, Wang Yujie,Wang Na 
// All Rights Reserved.
// This file is a part of LPOPC , published under the Eclipse Public License.
// Author:Xue Zhichen 7/14 2015   13:04
// Email:eddy_lpopc@163.com
#ifndef LPOPC_HESSIAN_H
#define LPOPC_HESSIAN_H
#include "LpConf.h"
#include "LpFunctionWrapper.h"
#include "LpCalculateData.hpp"
#include "LpOptimalProblem.hpp"
#include "LpFiniteDifferenceDerive.hpp"
#include "LpAnalyticDerive.hpp"
#include "LpDebug.hpp"
namespace Lpopc
{
	struct PhaseHessain
	{
		field<mat> hDae_state_state;
		field<mat> hPath_state_state;
		field<mat> hDae_control_state;
		field<mat> hPath_control_state;
		field<mat> hDae_control_control;
		field<mat> hPath_control_control;
		field<mat> hDae_time_state;
		field<mat> hPath_time_state;
		field<mat> hDae_time_control;
		field<mat> hPath_time_control;
		mat hdae_time_time;
		mat hPath_time_time;
		field<mat> hDae_para_state;
		field<mat> hPath_para_state;
		field<mat> hDae_para_control;
		field<mat> hPath_para_control;
		field<mat> hDae_para_time;
		field<mat> hPath_para_time;
		field<mat> hDae_para_para;
		field<mat> hPath_para_para;
		field<mat> hEvents_x0_x0;
		field<mat> hEvents_x0_xf;
		field<mat> hEvents_xf_x0;
		field<mat> hEvents_xf_xf;
		field<mat> hEvents_t0_x0;
		field<mat> hEvents_t0_xf;
		mat hEvents_t0_t0;
		field<mat> hEvents_tf_x0;
		field<mat> hEvents_tf_xf;
		mat hEvents_tf_t0;
		mat hEvents_tf_tf;
		field<mat>hEvents_para_x0;
		field<mat>hEvents_para_xf;
		field<mat>hEvents_para_t0;
		field<mat>hEvents_para_tf;
		field<mat>hEvents_para_para;
		mat hMayer_x0_x0;
		mat hMayer_x0_xf;
		mat hMayer_xf_x0;
		mat hMayer_xf_xf;
		vec hMayer_t0_x0;
		vec hMayer_t0_xf;
		double hMayer_t0_t0;
		vec hMayer_tf_x0;
		vec hMayer_tf_xf;
		double hMayer_tf_t0;
		double hMayer_tf_tf;
		mat hMayer_para_x0;
		mat hMayer_para_xf;
		vec hMayer_para_t0;
		vec hMayer_para_tf;
		mat hMayer_para_para;
		field<vec> hLagrange_x_x;
		field<vec> hLagrange_u_x;
		field<vec> hLagrange_u_u;
		field<vec> hLagrange_t_x;
		field<vec> hLagrange_t_u;
		vec        hLagrange_t_t;
		field<vec> hLagrange_p_x;
		field<vec> hLagrange_p_u;
		field<vec> hLagrange_p_t;
		field<vec> hLagrange_p_p;
	};
	struct LinkHessian
	{
		field<vec> hLink_xfL_xfL;
		field<vec> hLink_xfL_x0R;
		field<vec>hLink_paraL_xfL;
		field<vec>hLink_paraL_paraL;
		field<vec>hLink_paraL_x0R;
		field<vec> hLink_x0R_x0R;
		field<vec>hLink_paraR_xfL;
		field<vec>hLink_paraR_paraL;
		field<vec>hLink_paraR_x0R;
		field<vec>hLink_paraR_paraR;
	};
	class LpHessianCalculator
	{
		
	public:
		
		LpHessianCalculator(const shared_ptr<FunctionWrapper>& userfun, const shared_ptr<OptDerive> derive_fun, const shared_ptr<LpCalculateData>& userdata,
			const shared_ptr<OptimalProblem>& optpro, double tolerance)
			:fun_(userfun), Data_(userdata), optpro_(optpro), derive_fun_(derive_fun)
		{
			tol = tolerance;
		};
		
		void GetPhaseHessian(int iphase,double sigma,const vec& lambada, const vec& x_all, umat& idependencies, vec& Hessian_V);
		void GetPhaseHessianSparsity(int iphase, umat& idependencies, vec &Hessian_I, vec &Hessian_J);
		void GetHessian(double sigma,const vec&x,const vec& lambada, vec& hessain_V);
		void GetHessianSparsity(vec& Hessian_I, vec& Hessian_J);
		void GetLinkHessian(int ipair, const std::vector<SolCost>solcCostTotal,const vec& lambada,  vec& linkHessain_V);
		void GetLinkHessianSparsity(int ipair, vec& linkHessain_I, vec& linkHessain_J);
		void CalculatePhaseHessian(size_t iphase, const vec& x_all, PhaseHessain& iphase_hessian);
		void CalculateLinkHessain(size_t ipair, const std::vector<SolCost>solcCostTotal, LinkHessian& ilink_hessian);
		~LpHessianCalculator(){};

	private:
		LpHessianCalculator(){};
		LpHessianCalculator(const LpHessianCalculator&);
		LpHessianCalculator& operator=(const LpHessianCalculator&);

		shared_ptr<FunctionWrapper> fun_;
		shared_ptr<LpCalculateData> Data_;
		shared_ptr<OptimalProblem> optpro_;
		shared_ptr<OptDerive>      derive_fun_;
		double tol;
		std::vector<PhaseHessain> phase_hessian_;
		std::vector<LinkHessian>  link_hessian_;
		static int nnz(umat& op_mat)
		{

			int nonZeroNumber = 0;
			for (int i = 0; i < op_mat.n_elem; i++)
			{
				if (op_mat(i) != 0)
				{
					nonZeroNumber++;
				}
			}
			return nonZeroNumber;
		}

	};

	

}// namespace Lpopc
#endif // !LPOPC_HESSIAN_H
