// Copyright (C) 2014-2015 Xue Zhichen, Wang Yujie,Wang Na 
// All Rights Reserved.
// This file is a part of LPOPC , published under the Eclipse Public License.
// Author:Xue Zhichen 8/4 2015   21:25
// Email:eddy_lpopc@163.com
#include "LpSacleOCP.hpp"
namespace Lpopc
{
	/*
	% Compute the Column Scales Using Maximum and Minimum Values %
	% of Time, State, Control, and Parameters Provided by User   %
	*/
	void LpScaleOCP::CalculateVarScale()
	{
		std::vector<PhaseScaleShift> allScaleShift(data_->numphases_);
		std::vector<PhaseVarMaxMin> allVarMaxMin(data_->numphases_);
		for (size_t iphase = 0; iphase < data_->numphases_; iphase++)
		{
			size_t nstates = data_->SIZES_[iphase][0];
			size_t ncontrols = data_->SIZES_[iphase][1];
			size_t nparameters = data_->SIZES_[iphase][2];
			size_t npaths = data_->SIZES_[iphase][3];
			size_t nevents = data_->SIZES_[iphase][4];
			shared_ptr<Phase> currentphase(optpro_->GetPhase(iphase));
			/*
			Get the lower and upper limits on the initial and terminal
			% time in the current phase
			*/
			double t0_min, tf_min, t0_max, tf_max, t_min, t_max;
			currentphase->GetTimeMin(t0_min, tf_min);
			CheckInf(t0_min);
			CheckInf(tf_min);
			currentphase->GetTimeMax(t0_max, tf_max);
			CheckInf(t0_max);
			CheckInf(tf_max);

			t_min = t0_min;//Fix me!
			t_max = tf_max;

			double iphase_t0scale, iphase_t0shift, iphase_tfscale, iphase_tfshift;
			double iphase_tscale, iphase_tshift;

			GetScaleAndShift(t0_max, t0_min, iphase_t0scale, iphase_t0shift);
			GetScaleAndShift(tf_max, tf_min, iphase_tfscale, iphase_tfshift);

			GetScaleAndShift(t_max, t_min, iphase_tscale, iphase_tshift);
			/*Get the lower and upper limits on the states in the current phase*/
			int nodes = (currentphase->GetTotalNodes());
			size_t statenum = currentphase->GetstateMin().size();

			double minstate0 = 0, minstate = 0, minstatef = 0;
			double maxstate0 = 0, maxstate = 0, maxstatef = 0;
			vec iphase_statescale(nstates), iphase_stateshift(nstates);
			vec iphase_xmin(nstates), iphase_xmax(nstates);
			vec iphase_x0min(nstates), iphase_x0max(nstates);
			vec iphase_xfmin(nstates), iphase_xfmax(nstates);
			for (int istate = 0; istate < currentphase->GetstateMin().size(); istate++){
				currentphase->GetstateMin()[istate]->GetLimit(minstate0, minstate, minstatef);
				currentphase->GetstateMax()[istate]->GetLimit(maxstate0, maxstate, maxstatef);
				CheckInf(maxstate);
				CheckInf(minstate);
				iphase_xmin(istate) = minstate;
				iphase_xmax(istate) = maxstate;

				CheckInf(maxstate0);
				CheckInf(minstate0);
				iphase_x0min(istate) = minstate0;
				iphase_x0max(istate) = maxstate0;

				CheckInf(maxstatef);
				CheckInf(minstatef);
				iphase_xfmin(istate) = minstatef;
				iphase_xfmax(istate) = maxstatef;
				GetScaleAndShift(maxstate, minstate, iphase_statescale(istate), iphase_stateshift(istate));
			}
			//Get the lower and upper limits on the controls in the current phase
			double controlmin = 0.0, controlmax = 0.0;
			vec iphase_controlscale(ncontrols), iphase_controlshift(ncontrols);
			vec iphase_umin(ncontrols), iphase_umax(ncontrols);
			for (int icontrol = 0; icontrol < currentphase->GetcontrolMin().size(); icontrol++){
				controlmin = currentphase->GetcontrolMin()[icontrol];
				controlmax = currentphase->GetcontrolMax()[icontrol];
				CheckInf(controlmin);
				CheckInf(controlmax);
				iphase_umin(icontrol) = controlmin;
				iphase_umax(icontrol) = controlmax;
				GetScaleAndShift(controlmax, controlmin, iphase_controlscale(icontrol), iphase_controlshift(icontrol));

			}//end for int icontrol*/

			//Get the lower and upper limits on the static parameters in the current phase
			double parametermin = 0.0, parametermax = 0.0;
			vec iphase_parascale(nparameters), iphase_parashift(nparameters);
			vec iphase_pmin(nparameters), iphase_pmax(nparameters);
			for (int ipara = 0; ipara < currentphase->GetparameterMin().size(); ipara++){
				parametermin = currentphase->GetparameterMin()[ipara];
				parametermax = currentphase->GetparameterMax()[ipara];
				CheckInf(parametermin);
				CheckInf(parametermax);
				iphase_pmin(ipara) = parametermin;
				iphase_pmax(ipara) = parametermax;
				GetScaleAndShift(parametermax, parametermin, iphase_parascale(ipara), iphase_parashift(ipara));
			}//end for int ipara*/
			allScaleShift[iphase].t0scale = iphase_t0scale;
			allScaleShift[iphase].tfscale = iphase_tfscale;
			allScaleShift[iphase].tscale = iphase_tscale;
			allScaleShift[iphase].statescale = iphase_statescale;
			allScaleShift[iphase].controlscale = iphase_controlscale;
			allScaleShift[iphase].parameterscale = iphase_parascale;

			allScaleShift[iphase].t0shift = iphase_t0shift;
			allScaleShift[iphase].tfshift = iphase_tfshift;
			allScaleShift[iphase].tshift = iphase_tshift;
			allScaleShift[iphase].stateshift = iphase_stateshift;
			allScaleShift[iphase].controlshift = iphase_controlshift;
			allScaleShift[iphase].parametershift = iphase_parashift;

			allVarMaxMin[iphase].tmin = t_min;
			allVarMaxMin[iphase].tmax = t_max;
			allVarMaxMin[iphase].t0min = t0_min;
			allVarMaxMin[iphase].t0max = t0_max;
			allVarMaxMin[iphase].tfmin = tf_min;
			allVarMaxMin[iphase].tfmax = tf_max;
			allVarMaxMin[iphase].x_min = iphase_xmin;
			allVarMaxMin[iphase].x_max = iphase_xmax;
			allVarMaxMin[iphase].x0_min = iphase_x0min;
			allVarMaxMin[iphase].x0_max = iphase_x0max;
			allVarMaxMin[iphase].xf_min = iphase_xfmin;
			allVarMaxMin[iphase].xf_max = iphase_xfmax;
			allVarMaxMin[iphase].u_min = iphase_umin;
			allVarMaxMin[iphase].u_max = iphase_umax;
			allVarMaxMin[iphase].p_min = iphase_pmin;
			allVarMaxMin[iphase].p_max = iphase_pmax;
		}
		allScaleShift_ = allScaleShift;
		allVarMaxMin_ = allVarMaxMin;
	}
	void LpScaleOCP::CalculateFunScaleFromRand()
	{
		std::vector<std::vector<SolCost>> solCostTotal(data_->numphases_);
		PhaseFunScale iphasefunscale;
		std::vector<double> allphaseobjscale(data_->numphases_);
		objscales_ = allphaseobjscale;
		std::vector<PhaseFunScale> allPhaseFunscale(data_->numphases_);
		allPhaseFunscale_ = allPhaseFunscale;
		std::vector<vec> allLinkFunscale(data_->numlinkpairs_);
		allLinkFunscale_ = allLinkFunscale;
		for (size_t iphase = 0; iphase < data_->numphases_; iphase++)
		{
			std::vector<SolCost> isolCostTotal(nsamples_);
			size_t nstates = data_->SIZES_[iphase][0];
			size_t ncontrols = data_->SIZES_[iphase][1];
			size_t nparameters = data_->SIZES_[iphase][2];
			size_t npaths = data_->SIZES_[iphase][3];
			size_t nevents = data_->SIZES_[iphase][4];
			// get samples t,x,u,p
			double  tupp = allVarMaxMin_[iphase].tmax;
			double  tlow = allVarMaxMin_[iphase].tmin;

			double  t0upp = allVarMaxMin_[iphase].t0max;
			double  t0low = allVarMaxMin_[iphase].t0min;

			double  tfupp = allVarMaxMin_[iphase].tfmax;
			double  tflow = allVarMaxMin_[iphase].tfmin;
			vec randtime = randu(nsamples_, 1);
			vec trand = randtime*tupp + (1 - randtime)*tlow;
			vec t0rand = randtime*t0upp + (1 - randtime)*t0low;
			vec tfrand = randtime*tfupp + (1 - randtime)*tflow;
			mat xrand, x0rand, xfrand;
			if (nstates > 0)
			{
				vec xupp = allVarMaxMin_[iphase].x_max;
				vec xlow = allVarMaxMin_[iphase].x_min;

				vec x0upp = allVarMaxMin_[iphase].x0_max;
				vec x0low = allVarMaxMin_[iphase].x0_min;

				vec xfupp = allVarMaxMin_[iphase].xf_max;
				vec xflow = allVarMaxMin_[iphase].xf_min;
				xrand = randu(nsamples_, nstates);
				x0rand = randu(nsamples_, nstates);
				xfrand = randu(nsamples_, nstates);
				for (size_t istate = 0; istate < nstates; istate++)
				{
					xrand.col(istate) = xrand.col(istate)*xupp(istate) + (1 - xrand.col(istate))*xlow(istate);
					x0rand.col(istate) = x0rand.col(istate)*x0upp(istate) + (1 - x0rand.col(istate))*x0low(istate);
					xfrand.col(istate) = xfrand.col(istate)*xfupp(istate) + (1 - xfrand.col(istate))*xflow(istate);
				}

			}
			mat urand;
			if (ncontrols > 0)
			{
				vec uupp = allVarMaxMin_[iphase].u_max;
				vec ulow = allVarMaxMin_[iphase].u_min;
				urand = randu(nsamples_, ncontrols);
				for (size_t icontrol = 0; icontrol < ncontrols; icontrol++)
				{
					urand.col(icontrol) = urand.col(icontrol)*uupp(icontrol) + (1 - urand.col(icontrol))*ulow(icontrol);
				}
			}
			mat prand;
			if (nparameters > 0)
			{
				vec pupp = allVarMaxMin_[iphase].p_max;
				vec plow = allVarMaxMin_[iphase].p_min;
				prand = randu(nsamples_, nparameters);
				for (size_t ipara = 0; ipara < nparameters; ipara++)
				{
					prand.col(ipara) = prand.col(ipara)*pupp(ipara) + (1 - prand.col(ipara))*plow(ipara);
				}
			}
			//GetPathEventObjGrad from sample points
			SolDae testdae;
			SolEvent testevent;
			SolCost testcost;
			mat testeventgrad, dummdiffjac, testpathjac;

			vec obj_norm(nsamples_);
			mat path_norm(npaths, nsamples_);
			mat event_norm(nevents, nsamples_);

			for (size_t isample = 0; isample < nsamples_; isample++)
			{
				//Get Obj sample grad
				testcost.time_ = trand(isample);
				testcost.initial_time_ = t0rand(isample);
				testcost.terminal_time_ = tfrand(isample);
				testcost.state_ = testdae.state_ = xrand.row(isample);
				testcost.initial_state_ = trans(x0rand.row(isample));
				testcost.terminal_state_ = trans(xfrand.row(isample));
				testcost.phase_num_ = iphase + 1;
				if (ncontrols > 0)
				{
					testcost.control_ = urand.row(isample);
				}
				if (nparameters > 0)
				{
					testcost.parameter_ = trans(prand.row(isample));
				}
				isolCostTotal[isample] = testcost;
				rowvec testmayergrad;
				rowvec    testlagrangegrad;

				firstderive_->DerivMayer(testcost, testmayergrad);
				firstderive_->DerivLagrange(testcost, testlagrangegrad);
				obj_norm(isample) = sqrt(sum(testlagrangegrad%testlagrangegrad) + sum(testmayergrad%testmayergrad));
				/*///////////////////////////////////////////////*/
				// Get path sample jacobi
				if (npaths > 0)
				{
					testdae.phase_num_ = iphase + 1;
					testdae.time_ = trand(isample);
					testdae.state_ = xrand.row(isample);
					if (ncontrols > 0)
					{
						testdae.contol_ = urand.row(isample);
					}
					if (nparameters > 0)
					{
						testdae.parameter_ = trans(prand.row(isample));
					}
					firstderive_->DerivDae(testdae, dummdiffjac, testpathjac);//pathjac npath*(2*nstate+2*npath+t+npara);
					path_norm.col(isample) = sqrt(sum(testpathjac%testpathjac, 1));
				}//if npath>0

				//Get event sample jacobi
				if (nevents > 0)
				{
					testevent.phase_num_ = iphase + 1;
					testevent.initial_time_ = t0rand(isample);
					testevent.terminal_time_ = tfrand(isample);
					testevent.initial_state_ = trans(x0rand.row(isample));
					testevent.terminal_state_ = trans(xfrand.row(isample));
					if (nparameters > 0)
					{
						testevent.parameter_ = trans(prand.row(isample));
					}
					firstderive_->DerivEvent(testevent, testeventgrad);//eventgrad nevent*(2*nstate+2+npara)
					event_norm.col(isample) = sqrt(sum(testeventgrad%testeventgrad, 1));
				}//if nevent>0

			}// for isample
			
			double obj_average_norm;
			obj_average_norm = mean(obj_norm);
			if (true)// if obj_avnorm<eps
			{
				obj_average_norm = 1.0;//FIX
			}
			objscales_[iphase] = 1.0 / obj_average_norm;
			iphasefunscale.ode_scale = allScaleShift_[iphase].statescale;
			vec path_average_norm, event_average_norm;
			if (npaths > 0)
			{
				path_average_norm = mean(path_norm, 1);
				path_average_norm(find(path_average_norm < datum::eps)).fill(1.0);

				iphasefunscale.path_scale = 1.0 / path_average_norm;
			}
			if (nevents > 0)
			{
				event_average_norm = mean(event_norm, 1);
				event_average_norm(find(event_average_norm < datum::eps)).fill(1.0);
				iphasefunscale.event_scale = 1.0 / event_average_norm;
			}
			allPhaseFunscale_[iphase] = iphasefunscale;
			solCostTotal[iphase] = isolCostTotal;
		}// for iphase
		//get link_fun scales
		vec linkfunscale;
		for (size_t ipair = 0; ipair < data_->numlinkpairs_; ipair++)
		{
			int left_index = optpro_->GetLinkage(ipair)->LeftPhase();
			int right_index = optpro_->GetLinkage(ipair)->RightPhase();
			int nlinks = optpro_->GetLinkage(ipair)->GetLinkageMin().size();
			std::vector<SolCost> left_sol_cost = solCostTotal[left_index];
			std::vector<SolCost> right_sol_cost = solCostTotal[right_index];
			SolLink testlink;
			testlink.ipair = ipair + 1;
			testlink.left_phase_num_ = left_index + 1;
			testlink.right_phase_num_ = right_index + 1;
			mat link_norm(nlinks, nsamples_);
			mat testlinkjac;
			for (size_t isamples = 0; isamples < nsamples_; isamples++)
			{
				testlink.left_state_ = left_sol_cost[isamples].terminal_state_;
				testlink.left_parameter_ = left_sol_cost[isamples].parameter_;

				testlink.right_state_ = right_sol_cost[isamples].initial_state_;
				testlink.right_parameter_ = right_sol_cost[isamples].parameter_;
				firstderive_->DerivLink(testlink, testlinkjac);
				link_norm.col(isamples) = sqrt(sum(testlinkjac%testlinkjac, 1));
			}
			vec link_average_norm = mean(link_norm, 1);
			link_average_norm(find(link_average_norm < datum::eps)).fill(1.0);
			linkfunscale = 1.0 / link_average_norm;
			allLinkFunscale_[ipair] = linkfunscale;
		}// for ipair
		data_->objScale = objscales_;
	}

	void LpScaleOCP::GetScaleAndShift(double val_max, double val_min, double& scale, double& shift)
	{
		double istatediff = val_max - val_min;
		double istatesum = val_max + val_min;
		if (istatediff != 0)
		{
			scale = 1.0 / istatediff;
			shift = -(istatesum / istatediff) / 2.0;

		}
		else if (val_max == 0)
		{
			scale = 1.0;
			shift = 0.0;
		}
		else
		{
			scale = 1.0 / abs(val_max);
			shift = -val_max / abs(val_max);//-sign(val_max)
		}
	}

	void LpScaleOCP::CheckInf(double& val_)
	{
		if (abs(val_) == arma::datum::inf)
		{
			val_ = val_ / abs(val_)*big_num_as_inf_;
		}
	}

	void LpScaleOCP::GetOCPScaleAndShift()
	{
		int variable_offset = 0, constraint_offset = 0;
		std::vector<int> nodes = data_->totalnodes_perphase;

		vec var_scale(data_->varbounds_min.size());
		vec var_shift(data_->varbounds_min.size());
		vec fun_scale(data_->conbounds_min.size() + data_->linmax.n_elem);
		lp_index funindex_shift = 0;
		for (lp_index iphase = 0; iphase < data_->numphases_; iphase++)
		{
			PhaseScaleShift ivarscales = allScaleShift_[iphase];
			PhaseFunScale   ifunscales = allPhaseFunscale_[iphase];
			int inodes = nodes[iphase];
			size_t nstates = data_->SIZES_[iphase][0];
			size_t ncontrols = data_->SIZES_[iphase][1];
			size_t nparameters = data_->SIZES_[iphase][2];
			size_t npaths = data_->SIZES_[iphase][3];
			size_t nevents = data_->SIZES_[iphase][4];

			lp_index stateindex_start = data_->phase_indices[iphase]->state[0] - 1;
			lp_index stateindex_end = data_->phase_indices[iphase]->state[data_->phase_indices[iphase]->state.size() - 1] - 1;
			lp_index controlindex_start = data_->phase_indices[iphase]->control[0] - 1;
			lp_index controlindex_end = data_->phase_indices[iphase]->control[data_->phase_indices[iphase]->control.size() - 1] - 1;
			lp_index t0_index = data_->phase_indices[iphase]->time[0] - 1;
			lp_index tf_index = data_->phase_indices[iphase]->time[1] - 1;
			lp_index para_index_start, para_index_end;
			if (nparameters>0)
			{
				para_index_start = data_->phase_indices[iphase]->parameter[0] - 1;
				para_index_end = data_->phase_indices[iphase]->parameter[data_->phase_indices[iphase]->parameter.size() - 1] - 1;
			}
			

			/*
			Get the lower and upper limits on the initial and terminal
			% time in the current phase
			*/
			/*Get the lower and upper limits on the states in the current phase*/
			lp_index stateindex_shift = stateindex_start;
			for (lp_index istate = 0; istate < nstates; istate++){
				var_scale.subvec(stateindex_shift, stateindex_shift + inodes).fill(ivarscales.statescale(istate));
				var_shift.subvec(stateindex_shift, stateindex_shift + inodes).fill(ivarscales.stateshift(istate));
				stateindex_shift += inodes + 1;
			}
			//Get the lower and upper limits on the controls in the current phase
			lp_index controlindex_shift = controlindex_start;
			for (lp_index icontrol = 0; icontrol < ncontrols; icontrol++){
				var_scale.subvec(controlindex_shift, controlindex_shift + inodes - 1).fill(ivarscales.controlscale(icontrol));
				var_shift.subvec(controlindex_shift, controlindex_shift + inodes - 1).fill(ivarscales.controlshift(icontrol));
				controlindex_shift += inodes;
			}//end for int j*/
			var_scale(t0_index) = ivarscales.t0scale;
			var_shift(t0_index) = ivarscales.t0shift;

			var_scale(tf_index) = ivarscales.tfscale;
			var_shift(tf_index) = ivarscales.tfshift;
			//Get the lower and upper limits on the static parameters in the current phase
			if (nparameters > 0)
			{
				var_scale.subvec(para_index_start, para_index_end) = ivarscales.parameterscale;
				var_shift.subvec(para_index_start, para_index_end) = ivarscales.parametershift;
			}
			//Get the scale of diff fun in current phase
			funindex_shift = data_->constraint_indices[iphase][0] - 1;
			for (size_t istate = 0; istate < nstates; istate++)
			{
				fun_scale.subvec(funindex_shift, funindex_shift + inodes).fill(ifunscales.ode_scale(istate));
				funindex_shift += inodes;
			}
			//Get the scale of the path constraints in the current phase
			for (size_t ipath = 0; ipath < npaths; ipath++){
				fun_scale.subvec(funindex_shift, funindex_shift + inodes).fill(ifunscales.path_scale(ipath));
				funindex_shift += inodes;
			}//end for int j*/
			//Get the scale of  the event constraints in the current phase
			if (nevents > 0)
			{
				fun_scale.subvec(funindex_shift, funindex_shift + nevents-1) = ifunscales.event_scale;
				funindex_shift += nevents;
			}
		}//for iphase
		for (int ipair = 0; ipair < data_->numlinkpairs_; ipair++)
		{
			size_t nlinkfun = data_->link_indices[ipair].size();
			fun_scale.subvec(funindex_shift, funindex_shift + nlinkfun - 1) = allLinkFunscale_[ipair];
			funindex_shift += nlinkfun;
		}// for ipair

		//Set Scales on Linear Constraints equal to 1
		lp_index numlinearcons = data_->linmin.n_elem;//never be zero
		fun_scale.subvec(funindex_shift, funindex_shift + numlinearcons - 1).fill(1.0);
		data_->varscale = var_scale;
		data_->varshift = var_shift;
		data_->funscale = fun_scale;
// 		var_scale.save("var_scale.txt", raw_ascii);
// 		var_shift.save("var_shift.txt", raw_ascii);
// 		fun_scale.save("fun_scale.txt", raw_ascii);
	}

	void LpScaleOCP::CalculateFunScaleFromGuess()
	{
		// under test
		data_->nlpGuess;
		std::vector<std::vector<SolCost>> solCostTotal(data_->numphases_);
		PhaseFunScale iphasefunscale;
		std::vector<double> allphaseobjscale(data_->numphases_);
		objscales_ = allphaseobjscale;
		std::vector<PhaseFunScale> allPhaseFunscale(data_->numphases_);
		allPhaseFunscale_ = allPhaseFunscale;
		std::vector<vec> allLinkFunscale(data_->numlinkpairs_);
		allLinkFunscale_ = allLinkFunscale;
		for (size_t iphase = 0; iphase < data_->numphases_; iphase++)
		{
			std::vector<SolCost> isolCostTotal(nsamples_);
			size_t nstates = data_->SIZES_[iphase][0];
			size_t ncontrols = data_->SIZES_[iphase][1];
			size_t nparameters = data_->SIZES_[iphase][2];
			size_t npaths = data_->SIZES_[iphase][3];
			size_t nevents = data_->SIZES_[iphase][4];
			// get samples t,x,u,p
			double  tupp = allVarMaxMin_[iphase].tmax;
			double  tlow = allVarMaxMin_[iphase].tmin;

			double  t0upp = allVarMaxMin_[iphase].t0max;
			double  t0low = allVarMaxMin_[iphase].t0min;

			double  tfupp = allVarMaxMin_[iphase].tfmax;
			double  tflow = allVarMaxMin_[iphase].tfmin;
			vec randtime = randu(nsamples_, 1);
			vec trand = randtime*tupp + (1 - randtime)*tlow;
			vec t0rand = randtime*t0upp + (1 - randtime)*t0low;
			vec tfrand = randtime*tfupp + (1 - randtime)*tflow;
			mat xrand, x0rand, xfrand;
			if (nstates > 0)
			{
				vec xupp = allVarMaxMin_[iphase].x_max;
				vec xlow = allVarMaxMin_[iphase].x_min;

				vec x0upp = allVarMaxMin_[iphase].x0_max;
				vec x0low = allVarMaxMin_[iphase].x0_min;

				vec xfupp = allVarMaxMin_[iphase].xf_max;
				vec xflow = allVarMaxMin_[iphase].xf_min;
				xrand = randu(nsamples_, nstates);
				x0rand = randu(nsamples_, nstates);
				xfrand = randu(nsamples_, nstates);
				for (size_t istate = 0; istate < nstates; istate++)
				{
					xrand.col(istate) = xrand.col(istate)*xupp(istate) + (1 - xrand.col(istate))*xlow(istate);
					x0rand.col(istate) = x0rand.col(istate)*x0upp(istate) + (1 - x0rand.col(istate))*x0low(istate);
					xfrand.col(istate) = xfrand.col(istate)*xfupp(istate) + (1 - xfrand.col(istate))*xflow(istate);
				}

			}
			mat urand;
			if (ncontrols > 0)
			{
				vec uupp = allVarMaxMin_[iphase].u_max;
				vec ulow = allVarMaxMin_[iphase].u_min;
				urand = randu(nsamples_, ncontrols);
				for (size_t icontrol = 0; icontrol < ncontrols; icontrol++)
				{
					urand.col(icontrol) = urand.col(icontrol)*uupp(icontrol) + (1 - urand.col(icontrol))*ulow(icontrol);
				}
			}
			mat prand;
			if (nparameters > 0)
			{
				vec pupp = allVarMaxMin_[iphase].p_max;
				vec plow = allVarMaxMin_[iphase].p_min;
				prand = randu(nsamples_, nparameters);
				for (size_t ipara = 0; ipara < nparameters; ipara++)
				{
					prand.col(ipara) = prand.col(ipara)*pupp(ipara) + (1 - prand.col(ipara))*plow(ipara);
				}
			}
			//GetPathEventObjGrad from sample points
			SolDae testdae;
			SolEvent testevent;
			SolCost testcost;
			mat testeventgrad, dummdiffjac, testpathjac;

			vec obj_norm(nsamples_);
			mat path_norm(npaths, nsamples_);
			mat event_norm(nevents, nsamples_);

			for (size_t isample = 0; isample < nsamples_; isample++)
			{
				//Get Obj sample grad
				testcost.time_ = trand(isample);
				testcost.initial_time_ = t0rand(isample);
				testcost.terminal_time_ = tfrand(isample);
				testcost.state_ = testdae.state_ = xrand.row(isample);
				testcost.initial_state_ = trans(x0rand.row(isample));
				testcost.terminal_state_ = trans(xfrand.row(isample));
				testcost.phase_num_ = iphase + 1;
				if (ncontrols > 0)
				{
					testcost.control_ = urand.row(isample);
				}
				if (nparameters > 0)
				{
					testcost.parameter_ = trans(prand.row(isample));
				}
				isolCostTotal[isample] = testcost;
				rowvec testmayergrad;
				rowvec    testlagrangegrad;

				firstderive_->DerivMayer(testcost, testmayergrad);
				firstderive_->DerivLagrange(testcost, testlagrangegrad);
				obj_norm(isample) = sqrt(sum(testlagrangegrad%testlagrangegrad) + sum(testmayergrad%testmayergrad));
				/*///////////////////////////////////////////////*/
				// Get path sample jacobi
				if (npaths > 0)
				{
					testdae.phase_num_ = iphase + 1;
					testdae.time_ = trand(isample);
					testdae.state_ = xrand.row(isample);
					if (ncontrols > 0)
					{
						testdae.contol_ = urand.row(isample);
					}
					if (nparameters > 0)
					{
						testdae.parameter_ = trans(prand.row(isample));
					}
					firstderive_->DerivDae(testdae, dummdiffjac, testpathjac);//pathjac npath*(2*nstate+2*npath+t+npara);
					path_norm.col(isample) = sqrt(sum(testpathjac%testpathjac, 1));
				}//if npath>0

				//Get event sample jacobi
				if (nevents > 0)
				{
					testevent.phase_num_ = iphase + 1;
					testevent.initial_time_ = t0rand(isample);
					testevent.terminal_time_ = tfrand(isample);
					testevent.initial_state_ = trans(x0rand.row(isample));
					testevent.terminal_state_ = trans(xfrand.row(isample));
					if (nparameters > 0)
					{
						testevent.parameter_ = trans(prand.row(isample));
					}
					firstderive_->DerivEvent(testevent, testeventgrad);//eventgrad nevent*(2*nstate+2+npara)
					event_norm.col(isample) = sqrt(sum(testeventgrad%testeventgrad, 1));
				}//if nevent>0

			}// for isample

			double obj_average_norm;
			obj_average_norm = mean(obj_norm);
			if (true)// if obj_avnorm<eps
			{
				obj_average_norm = 1.0;//FIX
			}
			objscales_[iphase] = 1.0 / obj_average_norm;
			iphasefunscale.ode_scale = allScaleShift_[iphase].statescale;
			vec path_average_norm, event_average_norm;
			if (npaths > 0)
			{
				path_average_norm = mean(path_norm, 1);
				path_average_norm(find(path_average_norm < datum::eps)).fill(1.0);

				iphasefunscale.path_scale = 1.0 / path_average_norm;
			}
			if (nevents > 0)
			{
				event_average_norm = mean(event_norm, 1);
				event_average_norm(find(event_average_norm < datum::eps)).fill(1.0);
				iphasefunscale.event_scale = 1.0 / event_average_norm;
			}
			allPhaseFunscale_[iphase] = iphasefunscale;
			solCostTotal[iphase] = isolCostTotal;
		}// for iphase
		//get link_fun scales
		vec linkfunscale;
		for (size_t ipair = 0; ipair < data_->numlinkpairs_; ipair++)
		{
			int left_index = optpro_->GetLinkage(ipair)->LeftPhase();
			int right_index = optpro_->GetLinkage(ipair)->RightPhase();
			int nlinks = optpro_->GetLinkage(ipair)->GetLinkageMin().size();
			std::vector<SolCost> left_sol_cost = solCostTotal[left_index];
			std::vector<SolCost> right_sol_cost = solCostTotal[right_index];
			SolLink testlink;
			testlink.ipair = ipair + 1;
			testlink.left_phase_num_ = left_index + 1;
			testlink.right_phase_num_ = right_index + 1;
			mat link_norm(nlinks, nsamples_);
			mat testlinkjac;
			for (size_t isamples = 0; isamples < nsamples_; isamples++)
			{
				testlink.left_state_ = left_sol_cost[isamples].terminal_state_;
				testlink.left_parameter_ = left_sol_cost[isamples].parameter_;

				testlink.right_state_ = right_sol_cost[isamples].initial_state_;
				testlink.right_parameter_ = right_sol_cost[isamples].parameter_;
				firstderive_->DerivLink(testlink, testlinkjac);
				link_norm.col(isamples) = sqrt(sum(testlinkjac%testlinkjac, 1));
			}
			vec link_average_norm = mean(link_norm, 1);
			link_average_norm(find(link_average_norm < datum::eps)).fill(1.0);
			linkfunscale = 1.0 / link_average_norm;
			allLinkFunscale_[ipair] = linkfunscale;
		}// for ipair
		data_->objScale = objscales_;
	}



}//namespace Lpopc