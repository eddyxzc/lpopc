// Copyright (C) 2014-2015 Xue Zhichen, Wang Yujie,Wang Na 
// All Rights Reserved.
// This file is a part of LPOPC , published under the Eclipse Public License.
// Author:Xue Zhichen 8/4 2015   21:24
// Email:eddy_lpopc@163.com
#ifndef LPOPC_SCALE_OCP_HPP
#define LPOPC_SCALE_OCP_HPP
#include "LpConf.h"
#include "LpOptimalProblem.hpp"
#include "LpCalculateData.hpp"
#include "LpOptDerive.hpp"
namespace Lpopc
{
	/**auto-scale in based on Bett's book
	*and still underdevelopment
	*auto-scale is much solwer than manual  scale.
	*So the best choice is manual scale
	*/
	struct PhaseScaleShift
	{
		//var scale and shift
		vec statescale;
		vec controlscale;
		double t0scale;
		double tfscale;
		double tscale;
		vec parameterscale;

		vec stateshift;
		vec controlshift;
		double t0shift;
		double tfshift;
		double tshift;
		vec parametershift;

	};
	struct PhaseFunScale
	{
		// fun scale
		vec ode_scale;
		vec path_scale;
		vec event_scale;
	};
	class LpScaleOCP
	{
	public:
		struct  PhaseVarMaxMin
		{
			double tmin;
			double tmax;
			double t0min;
			double t0max;
			double tfmin;
			double tfmax;
			vec x_min;
			vec x_max;
			vec x0_min;
			vec x0_max;
			vec xf_min;
			vec xf_max;
			vec u_min;
			vec u_max;
			vec p_min;
			vec p_max;
		};
		LpScaleOCP(shared_ptr<OptimalProblem>& optpro,
			shared_ptr<LpCalculateData>& data,
			shared_ptr<OptDerive>& firstderive){
			optpro_ = optpro;
			data_ = data;
			firstderive_ = firstderive;
			big_num_as_inf_ = 1e20;
			nsamples_ = 90;
		};
		~LpScaleOCP(){};
		void GetOCPScaleAndShift();
		void CalculateVarScale();
		void CalculateFunScaleFromRand();
		void CalculateFunScaleFromGuess();
		std::vector<PhaseScaleShift>allScaleShift_;
		std::vector<PhaseFunScale>  allPhaseFunscale_;
		std::vector<vec> allLinkFunscale_;
	private:
		shared_ptr<OptimalProblem> optpro_;
		shared_ptr<LpCalculateData> data_;
		shared_ptr<OptDerive> firstderive_;
		double big_num_as_inf_;
		lp_index nsamples_;
		std::vector<PhaseVarMaxMin> allVarMaxMin_;
		std::vector<double> objscales_;
		void GetScaleAndShift(double val_max, double val_min, double& scale, double& shift);
		inline void CheckInf(double& val_);
		LpScaleOCP(const LpScaleOCP&);
		LpScaleOCP& operator= (const LpScaleOCP&);
	};
}//namespace Lpopc
#endif // !LPOPC_SCALE_OCP_HPP
