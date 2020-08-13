// Copyright (C) 2014-2015 Xue Zhichen, Wang Yujie,Wang Na 
// All Rights Reserved.
// This file is a part of LPOPC , published under the Eclipse Public License.
// Author:Xue Zhichen 7/14 2015   13:04
// Email:eddy_lpopc@163.com
#include "LpHessian.h"
namespace Lpopc
{



	void LpHessianCalculator::GetPhaseHessian(int iphase, double sigma, const vec& lambada, const vec& x_all, umat& idependencies, vec& Hessian_V)
	{
		LP_DBG_START_FUN("LpHessianCalculator::GetPhaseHessian")
		PhaseHessain iphase_hessain;
		CalculatePhaseHessian(iphase, x_all, iphase_hessain);
		int sumnodes = optpro_->GetPhase(iphase)->GetTotalNodes();
		size_t nstates = Data_->SIZES_[iphase][0];
		size_t ncontrols = Data_->SIZES_[iphase][1];
		size_t nparameters = Data_->SIZES_[iphase][2];
		size_t npaths = Data_->SIZES_[iphase][3];
		size_t nevents = Data_->SIZES_[iphase][4];
		int stateindex_start = Data_->phase_indices[iphase]->state[0] - 1;
		int stateindex_end = Data_->phase_indices[iphase]->state[Data_->phase_indices[iphase]->state.size() - 1] - 1;
		int controlindex_start = Data_->phase_indices[iphase]->control[0] - 1;
		int controlindex_end = Data_->phase_indices[iphase]->control[Data_->phase_indices[iphase]->control.size() - 1] - 1;
		vec state_vector = x_all.subvec(stateindex_start, stateindex_end);
		vec control_vector = x_all.subvec(controlindex_start, controlindex_end);
		int t0_index = Data_->phase_indices[iphase]->time[0] - 1;
		int tf_index = Data_->phase_indices[iphase]->time[1] - 1;
		double t0 = x_all(t0_index);
		double tf = x_all(tf_index);
		double tspan = tf - t0;
		vec t_radau = (Data_->PS[iphase]->Points + 1)*(tspan / 2.0) + t0;
		//t_radau.print("tradau");
		mat state_matrix = reshape(state_vector, sumnodes + 1, nstates);
		//state_matrix.print("state_matrix");
		mat state_radau = state_matrix.rows(0, state_matrix.n_rows - 2);//index start form 0

		vec x0 = trans(state_matrix.row(0));
		vec xf = trans(state_matrix.row(state_matrix.n_rows - 1));
		mat control_radau = reshape(control_vector, sumnodes, ncontrols);
		mat parameters;
		if (nparameters > 0)
		{
			int para_index_start = Data_->phase_indices[iphase]->parameter[0] - 1;
			int para_index_end = Data_->phase_indices[iphase]->parameter[Data_->phase_indices[iphase]->parameter.size() - 1] - 1;
			parameters = x_all.rows(para_index_start, para_index_end);
		}
		//get derivatives of DAE function
		SolDae mySolDae;
		mySolDae.time_ = t_radau;
		mySolDae.state_ = state_radau;
		mySolDae.contol_ = control_radau;
		mySolDae.parameter_ = parameters;
		mySolDae.phase_num_ = iphase + 1;
		SolCost mysolcost;
		mysolcost.initial_time_ = t0;
		mysolcost.initial_state_ = x0;
		mysolcost.terminal_time_ = tf;
		mysolcost.terminal_state_ = xf;
		mysolcost.time_ = t_radau;
		mysolcost.state_ = state_radau;
		mysolcost.control_ = control_radau;
		mysolcost.parameter_ = parameters;
		mysolcost.phase_num_ = iphase + 1;

		//allocate the hessian ,obtain the struct and the value

		size_t con_index_start = Data_->constraint_indices[iphase][0] - 1;
		size_t con_index_end = Data_->constraint_indices[iphase][Data_->constraint_indices[iphase].size()-1] - 1;

		size_t diff_index_start = 0;
		size_t diff_index_end = diff_index_start + sumnodes*nstates - 1;
		size_t path_index_start = 0, path_index_end = 0, event_index_start = 0, event_index_end = 0;
		if (npaths > 0)
		{
			path_index_start = sumnodes*nstates;
			path_index_end = sumnodes*nstates + npaths*sumnodes - 1;
		}
		if (nevents > 0)
		{
			event_index_start = sumnodes*(nstates + npaths);
			event_index_end = sumnodes*(nstates + npaths) + nevents - 1;
		}
		//get lambda of differential constrains , path constrains and events constrains
		vec phase_lambda = lambada.subvec(con_index_start, con_index_end);
		mat diff_lambda(sumnodes, nstates);
		mat path_lambda(sumnodes, npaths);
		vec event_lambda(nevents);
		size_t lambda_shift = 0;
		for (size_t idiff = 0; idiff < nstates; idiff++)
		{
			diff_lambda.col(idiff) = phase_lambda.subvec(lambda_shift, lambda_shift + sumnodes - 1);
			lambda_shift += sumnodes;
		}
		for (size_t ipath = 0; ipath < npaths; ipath++)
		{
			path_lambda.col(ipath) = phase_lambda.subvec(lambda_shift, lambda_shift + sumnodes - 1);
			lambda_shift += sumnodes;
		}
		for (size_t ievent = 0; ievent < nevents; ievent++)
		{
			event_lambda(ievent) = phase_lambda(lambda_shift);
			lambda_shift++;
		}
		//////////////////////////////////////////////////////////////////////////
		/*dd_LI=
		*dd_xx
		*dd_ux, dd_uu
		*dd_t0x,dd_t0u,dd_t0t0
		*dd_tfx,dd_tfu,dd_tft0,dd_tftf
		*dd_px, dd_pu, dd_pt0 ,dd_ptf, dd_pp
		*/
		//get dd_xx
		vec rpm_w = Data_->PS[iphase]->Weights;
		vec rpm_tau = Data_->PS[iphase]->Points;
		field<vec> hLI_xx(nstates, nstates);
		for (size_t istate = 0; istate < nstates; istate++)
		{
			for (size_t jstate = 0; jstate <= istate; jstate++)
			{
				vec sum_lambda_plus_hdae_xx = sum(diff_lambda % iphase_hessain.hDae_state_state(istate, jstate), 1);
				vec sum_lambda_plus_hpath_xx = zeros<vec>(sumnodes);
				if (npaths>0)sum_lambda_plus_hpath_xx = sum(path_lambda % iphase_hessain.hPath_state_state(istate, jstate), 1);
				vec sigma_plus_hL_xx = sigma*rpm_w%iphase_hessain.hLagrange_x_x(istate, jstate);
				hLI_xx(istate, jstate) = (tf - t0) / 2.0 * (sigma_plus_hL_xx - sum_lambda_plus_hdae_xx) + sum_lambda_plus_hpath_xx;
			}//for jstate
		}//for istate
		//get dd_ux
		field<vec> hLI_ux(ncontrols, nstates);
		for (size_t icontrol = 0; icontrol < ncontrols; icontrol++)
		{
			for (size_t jstate = 0; jstate < nstates; jstate++)
			{
				vec sum_lambda_plus_hdae_ux = sum(diff_lambda % iphase_hessain.hDae_control_state(icontrol, jstate), 1);
				vec sum_lambda_plus_hpath_ux = zeros<vec>(sumnodes);
				if (npaths>0) sum_lambda_plus_hpath_ux = sum(path_lambda %iphase_hessain.hPath_control_state(icontrol, jstate), 1);
				vec sigma_plus_hL_ux = sigma*rpm_w%iphase_hessain.hLagrange_u_x(icontrol, jstate);
				hLI_ux(icontrol, jstate) = (tf - t0) / 2.0 * (sigma_plus_hL_ux - sum_lambda_plus_hdae_ux) + sum_lambda_plus_hpath_ux;
			}//for jstate
		}// for istate
		//get dd_uu
		field<vec> hLI_uu(ncontrols, ncontrols);
		for (size_t icontrol = 0; icontrol < ncontrols; icontrol++)
		{
			for (size_t jcontrol = 0; jcontrol <= icontrol; jcontrol++)
			{
				vec sum_lambda_plus_hdae_uu = sum(diff_lambda % iphase_hessain.hDae_control_control(icontrol, jcontrol), 1);
				vec sum_lambda_plus_hpath_uu = zeros<vec>(sumnodes);
				if (npaths>0)sum_lambda_plus_hpath_uu = sum(path_lambda % iphase_hessain.hPath_control_control(icontrol, jcontrol), 1);
				vec sigma_plus_hL_uu = sigma*rpm_w%iphase_hessain.hLagrange_u_u(icontrol, jcontrol);
				hLI_uu(icontrol, jcontrol) = (tf - t0) / 2.0 * (sigma_plus_hL_uu - sum_lambda_plus_hdae_uu) + sum_lambda_plus_hpath_uu;
			}//for jstate
		}// for istate
		// get dd_t0x
		field<vec> hLI_t0x(nstates);
		//get dd_tfx
		field<vec> hLI_tfx(nstates);
		mat dstate, dpath;
		derive_fun_->DerivDae(mySolDae, dstate, dpath);
		mat dstate_x = dstate.cols(0, nstates - 1);
		mat dstate_u;
		if (ncontrols > 0)
		{
			dstate_u = dstate.cols(nstates, nstates + ncontrols - 1);
		}
		vec dstate_t = dstate.col(nstates + ncontrols);
		mat dstate_parameter;
		if (nparameters > 0)
		{
			dstate_parameter = dstate.tail_cols(nparameters);
		}
		mat dLagrange;
		derive_fun_->DerivLagrange(mysolcost, dLagrange);
																					
		mat talpha = (1 - rpm_tau) / 2.0, tbeta = (1 + rpm_tau) / 2.0;
		for (size_t istate = 0; istate < nstates; istate++)
		{
			vec sum_lambda_plus_hdae_tx = sum(diff_lambda % iphase_hessain.hDae_time_state(istate), 1);
			vec sum_lambda_plus_hpath_tx = zeros<vec>(sumnodes);
			if (npaths>0)sum_lambda_plus_hpath_tx = sum(path_lambda % iphase_hessain.hPath_time_state(istate), 1);
			vec sigma_plus_hL_tx = sigma*rpm_w % iphase_hessain.hLagrange_t_x(istate);
			vec sum_lambda_plus_ddae_x = sum(diff_lambda%reshape(dstate_x.col(istate),sumnodes,nstates), 1);
			vec sigma_plus_dLagrange_x = sigma*rpm_w%dLagrange.col(istate);
			hLI_t0x(istate) = 0.5*(sum_lambda_plus_ddae_x - sigma_plus_dLagrange_x) +
				talpha % ((tf - t0) / 2.0 * (sigma_plus_hL_tx - sum_lambda_plus_hdae_tx) + sum_lambda_plus_hpath_tx);
			hLI_tfx(istate) = -0.5*(sum_lambda_plus_ddae_x - sigma_plus_dLagrange_x) +
				tbeta % ((tf - t0) / 2.0 * (sigma_plus_hL_tx - sum_lambda_plus_hdae_tx) + sum_lambda_plus_hpath_tx);
		}

		field<vec> hLI_t0u(ncontrols);
		field<vec> hLI_tfu(ncontrols);
		for (size_t icontrol = 0; icontrol < ncontrols; icontrol++)
		{
			vec sum_lambda_plus_hdae_tu = sum(diff_lambda % iphase_hessain.hDae_time_control(icontrol), 1);
			vec sum_lambda_plus_hpath_tu = zeros<vec>(sumnodes);
			if (npaths>0)sum_lambda_plus_hpath_tu = sum(path_lambda % iphase_hessain.hPath_time_control(icontrol), 1);
			vec sigma_plus_hL_tu = sigma*rpm_w % iphase_hessain.hLagrange_t_u(icontrol);
			vec sum_lambda_plus_ddae_u = sum(diff_lambda%reshape(dstate_u.col(icontrol),sumnodes,nstates), 1);
			vec sigma_plus_dLagrange_u = sigma*rpm_w%dLagrange.col(nstates + icontrol);
			hLI_t0u(icontrol) = 0.5*(sum_lambda_plus_ddae_u - sigma_plus_dLagrange_u) +
				talpha % ((tf - t0) / 2.0 * (sigma_plus_hL_tu - sum_lambda_plus_hdae_tu) + sum_lambda_plus_hpath_tu);
			hLI_tfu(icontrol) = -0.5*(sum_lambda_plus_ddae_u - sigma_plus_dLagrange_u) +
				tbeta % ((tf - t0) / 2.0 * (sigma_plus_hL_tu - sum_lambda_plus_hdae_tu) + sum_lambda_plus_hpath_tu);
		}
		vec hLI_t0t0;
		vec hLI_tft0;
		vec hLI_tftf;
		vec sigma_plus_hL_tt = sigma*rpm_w % iphase_hessain.hLagrange_t_t;
		vec sum_lambda_plus_hdae_tt = sum(diff_lambda % iphase_hessain.hdae_time_time, 1);
		vec sum_lambda_plus_hpath_tt = zeros < vec >(sumnodes);
		if (npaths > 0)sum_lambda_plus_hpath_tt = sum(path_lambda % iphase_hessain.hPath_time_time, 1);
		vec sigma_plus_dLagrange_t = sigma*rpm_w%dLagrange.col(nstates + ncontrols);
		vec sum_lambda_plus_ddae_t = sum(diff_lambda%reshape(dstate_t, sumnodes, nstates), 1);
		hLI_t0t0 = trans(talpha)*((sum_lambda_plus_ddae_t - sigma_plus_dLagrange_t) + talpha % ((tf - t0) / 2.0 * (sigma_plus_hL_tt - sum_lambda_plus_hdae_tt) + sum_lambda_plus_hpath_tt));
		hLI_tftf = trans(tbeta)*(-(sum_lambda_plus_ddae_t - sigma_plus_dLagrange_t) + tbeta % ((tf - t0) / 2.0 * (sigma_plus_hL_tt - sum_lambda_plus_hdae_tt) + sum_lambda_plus_hpath_tt));
		hLI_tft0 = 0.5*(trans(tbeta - talpha))*(sum_lambda_plus_ddae_t - sigma_plus_dLagrange_t) + trans(talpha)*(tbeta % ((tf - t0) / 2.0 * (sigma_plus_hL_tt - sum_lambda_plus_hdae_tt) + sum_lambda_plus_hpath_tt));

		field<vec> hLI_px(nparameters, nstates);
		for (size_t ipara = 0; ipara < nparameters; ipara++)
		{
			for (size_t jstate = 0; jstate < nstates; jstate++)
			{
				vec sum_lambda_plus_hdae_px = sum(diff_lambda %iphase_hessain.hDae_para_state(ipara, jstate), 1);
				vec sum_lambda_plus_hpath_px = zeros<vec>(sumnodes);
				if (npaths>0)sum_lambda_plus_hpath_px = sum(path_lambda %iphase_hessain.hPath_para_state(ipara, jstate), 1);
				vec sigma_plus_hL_px = sigma*rpm_w%iphase_hessain.hLagrange_p_x(ipara, jstate);
				hLI_px(ipara, jstate) = (tf - t0) / 2.0 * (sigma_plus_hL_px - sum_lambda_plus_hdae_px) + sum_lambda_plus_hpath_px;
			}//for jstate
		}//for istate
		field<vec> hLI_pu(nparameters, ncontrols);
		for (size_t ipara = 0; ipara < nparameters; ipara++)
		{
			for (size_t jcontrol = 0; jcontrol < ncontrols; jcontrol++)
			{
				vec sum_lambda_plus_hdae_pu = sum(diff_lambda % iphase_hessain.hDae_para_control(ipara, jcontrol), 1);
				vec sum_lambda_plus_hpath_pu = zeros<vec>(sumnodes);
				if (npaths>0)sum_lambda_plus_hpath_pu = sum(path_lambda % iphase_hessain.hPath_control_control(ipara, jcontrol), 1);
				vec sigma_plus_hL_uu = sigma*rpm_w%iphase_hessain.hLagrange_p_u(ipara, jcontrol);
				hLI_pu(ipara, jcontrol) = (tf - t0) / 2.0 * (sigma_plus_hL_uu - sum_lambda_plus_hdae_pu) + sum_lambda_plus_hpath_pu;
			}//for jstate
		}// for ipara

		vec hLI_pt0(nparameters);
		vec hLI_ptf(nparameters);
		for (size_t ipara = 0; ipara < nparameters; ipara++)
		{
			vec sum_lambda_plus_hdae_pt = sum(diff_lambda % iphase_hessain.hDae_para_time(ipara), 1);
			vec sum_lambda_plus_hpath_pt = zeros<vec>(sumnodes);
			if (npaths>0)sum_lambda_plus_hpath_pt = sum(path_lambda %iphase_hessain.hPath_para_time(ipara), 1);
			vec sigma_plus_hL_pt = sigma*rpm_w %iphase_hessain.hLagrange_p_t(ipara);
			double sum_lambda_plus_ddae_p = accu(diff_lambda%reshape(dstate_parameter.col(ipara),sumnodes,nstates));
			double sigma_plus_dLagrange_p = accu(sigma*rpm_w%dLagrange.col(nstates + ncontrols + 1 + ipara));
			hLI_pt0(ipara) = 0.5*(sum_lambda_plus_ddae_p - sigma_plus_dLagrange_p) +
				accu(talpha % ((tf - t0) / 2.0 * (sigma_plus_hL_pt - sum_lambda_plus_hdae_pt) + sum_lambda_plus_hpath_pt));
			hLI_ptf(ipara) = -0.5*(sum_lambda_plus_ddae_p - sigma_plus_dLagrange_p) +
				accu(tbeta % ((tf - t0) / 2.0 * (sigma_plus_hL_pt - sum_lambda_plus_hdae_pt) + sum_lambda_plus_hpath_pt));
		}

		mat hLI_pp(nparameters, nparameters);
		for (size_t ipara = 0; ipara < nparameters; ipara++)
		{
			for (size_t jpara = 0; jpara < ipara; jpara++)
			{
				double sum_lambda_plus_hdae_pp = accu(diff_lambda% iphase_hessain.hDae_para_para(ipara, jpara));
				double sum_lambda_plus_hpath_pp = 0;
				if (npaths>0)sum_lambda_plus_hpath_pp = accu(diff_lambda %iphase_hessain.hPath_para_para(ipara, jpara));
				double sigma_plus_hL_pp = accu(sigma*rpm_w%iphase_hessain.hLagrange_p_p(ipara, jpara));
				hLI_pp(ipara, jpara) = (tf - t0) / 2.0 * (sigma_plus_hL_pp - sum_lambda_plus_hdae_pp) + sum_lambda_plus_hpath_pp;
			}//for jpara
		}//for ipara

		//////////////////////////////////////////////////////////////////////////
		/*dd_LE=
		*dd_xx
		*dd_ux(0),dd_uu (0)
		*dd_t0x,  dd_t0u(0), dd_t0t0
		*dd_tfx,  dd_tfu(0), dd_tft0,dd_tftf
		*dd_px,   dd_pu (0), dd_pt0 ,dd_ptf, dd_pp
		*/
		mat hE_x0x0(nstates, nstates);
		mat hE_x0xf(nstates, nstates);
		mat hE_xfx0(nstates, nstates);
		mat hE_xfxf(nstates, nstates);

		for (size_t istate = 0; istate < nstates; istate++)
		{
			for (size_t jstate = 0; jstate < nstates; jstate++)
			{
				double lambda_sum_event=0.0;

				if(nevents>0) lambda_sum_event = accu(iphase_hessain.hEvents_x0_x0(istate, jstate) % event_lambda);
				double sigma_plus_mayer = sigma*iphase_hessain.hMayer_x0_x0(istate, jstate);
				hE_x0x0(istate, jstate) = sigma_plus_mayer + lambda_sum_event;

				if (nevents>0) lambda_sum_event = accu(iphase_hessain.hEvents_x0_xf(istate, jstate) % event_lambda);
				sigma_plus_mayer = sigma*iphase_hessain.hMayer_x0_xf(istate, jstate);
				hE_x0xf(istate, jstate) = sigma_plus_mayer + lambda_sum_event;

				if (nevents>0) lambda_sum_event = accu(iphase_hessain.hEvents_xf_x0(istate, jstate) % event_lambda);
				sigma_plus_mayer = sigma*iphase_hessain.hMayer_xf_x0(istate, jstate);
				hE_xfx0(istate, jstate) = sigma_plus_mayer + lambda_sum_event;

				if (nevents>0)  lambda_sum_event = accu(iphase_hessain.hEvents_xf_xf(istate, jstate) % event_lambda);
				sigma_plus_mayer = sigma*iphase_hessain.hMayer_xf_xf(istate, jstate);
				hE_xfxf(istate, jstate) = sigma_plus_mayer + lambda_sum_event;
			}//for jstate
		}//for istate
		vec hE_t0x0(nstates);
		vec hE_t0xf(nstates);
		vec hE_tfx0(nstates);
		vec hE_tfxf(nstates);
		for (size_t istate = 0; istate < nstates; istate++)
		{
			double lambda_sum_event = 0.0;

			if(nevents>0) lambda_sum_event = accu(iphase_hessain.hEvents_t0_x0(istate) % event_lambda);
			double sigma_plus_mayer = sigma*iphase_hessain.hMayer_t0_x0(istate);
			hE_t0x0(istate) = sigma_plus_mayer + lambda_sum_event;

			if (nevents>0)lambda_sum_event = accu(iphase_hessain.hEvents_t0_xf(istate) % event_lambda);
			sigma_plus_mayer = sigma*iphase_hessain.hMayer_t0_xf(istate);
			hE_t0xf(istate) = sigma_plus_mayer + lambda_sum_event;

			if (nevents>0)lambda_sum_event = accu(iphase_hessain.hEvents_tf_x0(istate) % event_lambda);
			sigma_plus_mayer = sigma*iphase_hessain.hMayer_tf_x0(istate);
			hE_tfx0(istate) = sigma_plus_mayer + lambda_sum_event;

			if (nevents>0)lambda_sum_event = accu(iphase_hessain.hEvents_tf_xf(istate) % event_lambda);
			sigma_plus_mayer = sigma*iphase_hessain.hMayer_tf_xf(istate);
			hE_tfxf(istate) = sigma_plus_mayer + lambda_sum_event;
		}

		double hE_t0t0, hE_tft0, hE_tftf;
		double lambda_sum_event = 0.0;
		if(nevents>0) lambda_sum_event = accu(iphase_hessain.hEvents_t0_t0 % event_lambda);
		double sigma_plus_mayer = sigma*iphase_hessain.hMayer_t0_t0;
		hE_t0t0 = sigma_plus_mayer + lambda_sum_event;

		if (nevents>0)lambda_sum_event = accu(iphase_hessain.hEvents_tf_t0 % event_lambda);
		sigma_plus_mayer = sigma*iphase_hessain.hMayer_tf_t0;
		hE_tft0 = sigma_plus_mayer + lambda_sum_event;

		if (nevents>0)lambda_sum_event = accu(iphase_hessain.hEvents_tf_tf % event_lambda);
		sigma_plus_mayer = sigma*iphase_hessain.hMayer_tf_tf;
		hE_tftf = sigma_plus_mayer + lambda_sum_event;

		mat hE_px0(nparameters, nstates);
		mat hE_pxf(nparameters, nstates);
		for (size_t ipara = 0; ipara < nparameters; ipara++)
		{
			for (size_t jstate = 0; jstate < nstates; jstate++)
			{
				double lambda_sum_event = 0;
				if(nevents>0) lambda_sum_event = accu(iphase_hessain.hEvents_para_x0(ipara, jstate) % event_lambda);
				double sigma_plus_mayer = sigma*iphase_hessain.hMayer_para_x0(ipara, jstate);
				hE_px0(ipara, jstate) = sigma_plus_mayer + lambda_sum_event;

				if (nevents>0) lambda_sum_event = accu(iphase_hessain.hEvents_para_xf(ipara, jstate) % event_lambda);
				sigma_plus_mayer = sigma*iphase_hessain.hMayer_para_xf(ipara, jstate);
				hE_pxf(ipara, jstate) = sigma_plus_mayer + lambda_sum_event;
			}// for jstate
		}//for ipara

		vec hE_pt0(nparameters);
		vec hE_ptf(nparameters);
		for (size_t ipara = 0; ipara < nparameters; ipara++)
		{
			double lambda_sum_event = 0;
			if (nevents>0) lambda_sum_event = accu(iphase_hessain.hEvents_para_t0(ipara) % event_lambda);
			double sigma_plus_mayer = sigma*iphase_hessain.hMayer_para_t0(ipara);
			hE_pt0(ipara) = sigma_plus_mayer + lambda_sum_event;

			if (nevents>0) lambda_sum_event = accu(iphase_hessain.hEvents_para_tf(ipara) % event_lambda);
			sigma_plus_mayer = sigma*iphase_hessain.hMayer_para_tf(ipara);
			hE_ptf(ipara) = sigma_plus_mayer + lambda_sum_event;
		}
		mat hE_pp(nparameters, nparameters);
		for (size_t ipara = 0; ipara < nparameters; ipara++)
		{
			for (size_t jpara = 0; jpara < nparameters; jpara++)
			{
				double lambda_sum_event = 0;
				if (nevents>0) lambda_sum_event = accu(iphase_hessain.hEvents_para_para(ipara, jpara) % event_lambda);
				double sigma_plus_mayer = sigma*iphase_hessain.hMayer_para_para(ipara, jpara);
				hE_pp(ipara, jpara) = sigma_plus_mayer + lambda_sum_event;
			}//for jpara
		}//for ipara

		int nDependancies = nnz(idependencies);
		//allocate memory for Hessain

		//get the number of nonzeros of the Lagrangian that are functions of collocation points functions
		size_t hessainI_nonzeros = sumnodes*((nDependancies - nstates - ncontrols) / 2 + nstates + ncontrols) +(2 + nparameters)*(nstates + ncontrols)*sumnodes
			+ (2 + nparameters)*(2 + nparameters - 1) / 2 + nparameters + 2;
		vec ShessainI_V(hessainI_nonzeros);

		//get the number of nonzeros of the Lagrangian that are functions of endpoints functions
		size_t hessainE_nonzeros = (nstates + nstates)*(nstates + nstates - 1) / 2 + nstates + nstates +
			2 * (2 + nparameters)*(nstates ) + (2 + nparameters)*(2 + nparameters - 1) / 2 + nparameters + 2;
		vec ShessainE_V(hessainE_nonzeros);

		//differential equations
		int rowstart, colstart;
		int ShessI_rowShift = 0, ShessE_rowShift = 0;
		vec indexvector = linspace<vec>(0, sumnodes - 1, sumnodes);
		//inset hI_xx
		for (int i = 0; i < nstates; i++)
		{
			rowstart = i*(sumnodes + 1);
			for (int j = 0; j <= i; j++)
			{
				colstart = j*(sumnodes + 1);
				if (idependencies(i, j))
				{
					
					ShessainI_V.subvec(ShessI_rowShift, ShessI_rowShift + sumnodes - 1) = hLI_xx(i, j);
					ShessI_rowShift += sumnodes;
				}
				ShessainE_V(ShessE_rowShift) = hE_x0x0(i, j);//xi0,xj0
				ShessE_rowShift++;
				if (i != j)
				{
					ShessainE_V(ShessE_rowShift) = hE_x0xf(i, j);
					ShessE_rowShift++;
				}
				ShessainE_V(ShessE_rowShift) = hE_xfx0(i, j);//xif,xj0
				ShessE_rowShift++;
				ShessainE_V(ShessE_rowShift) = hE_xfxf(i, j);//xif,xjf
				ShessE_rowShift++;
			}//for j
		}//for i
		//insert hI_ux.hI_uu
		int rowshift = nstates*(sumnodes + 1);
		for (int i = 0; i < ncontrols; i++)
		{
			rowstart = rowshift + i*sumnodes;
			for (int j = 0; j < nstates; j++)//ux
			{
				colstart = j*(sumnodes + 1);
				if (idependencies(i + nstates, j) )
				{
					//Insert hI_ux
					ShessainI_V.subvec(ShessI_rowShift, ShessI_rowShift + sumnodes - 1) = hLI_ux(i, j);
					ShessI_rowShift += sumnodes;
				}
			}//for j

			int colshift = nstates*(sumnodes + 1);
			for (int j = 0; j <= i; j++)
			{
				colstart = colshift + j*sumnodes;
				if (idependencies(i + nstates, j + nstates))
				{
					//insert hI_uu
					ShessainI_V.subvec(ShessI_rowShift, ShessI_rowShift + sumnodes - 1) = hLI_uu(i, j);
					ShessI_rowShift += sumnodes;
				}
			}//for j
		}//for i
		//////////////////////////////////////////////////////////////////////////
		//insert hI_t0x,t0u,t0t0;
		rowshift = nstates*(sumnodes + 1) + ncontrols*sumnodes;
		rowstart = rowshift;

		// hI_t0x
		for (size_t i = 0; i < nstates; i++)
		{
			colstart = i*(sumnodes + 1);
			ShessainI_V.subvec(ShessI_rowShift, ShessI_rowShift + sumnodes - 1) = (hLI_t0x(i));
			ShessI_rowShift += sumnodes;

			ShessainE_V(ShessE_rowShift) = hE_t0x0(i);//t0,x0
			ShessE_rowShift++;

			ShessainE_V(ShessE_rowShift) = hE_t0xf(i);//t0,xf
			ShessE_rowShift++;
		}
		//hI_t0u
		for (size_t i = 0; i < ncontrols; i++)
		{
			colstart = i*(sumnodes)+nstates*(sumnodes + 1);
			ShessainI_V.subvec(ShessI_rowShift, ShessI_rowShift + sumnodes - 1) = (hLI_t0u(i));
			ShessI_rowShift += sumnodes;
		}
		//hI_t0t0
		colstart = nstates*(sumnodes + 1) + ncontrols*sumnodes;
		ShessainI_V(ShessI_rowShift) = as_scalar(hLI_t0t0);
		ShessI_rowShift++;

		ShessainE_V(ShessE_rowShift) = hE_t0t0;//t0,xf
		ShessE_rowShift++;

		//////////////////////////////////////////////////////////////////////////
		////insert hI_tfx,tfu,tft0.tf,tf;

		rowstart++;
		// hI_tfx
		for (size_t i = 0; i < nstates; i++)
		{
			colstart = i*(sumnodes + 1);
			ShessainI_V.subvec(ShessI_rowShift, ShessI_rowShift + sumnodes - 1) = (hLI_tfx(i));
			ShessI_rowShift += sumnodes;

			ShessainE_V(ShessE_rowShift) = hE_tfx0(i);//tf,x0
			ShessE_rowShift++;

			ShessainE_V(ShessE_rowShift) = hE_tfxf(i);//tf,xf
			ShessE_rowShift++;
		}
		//hI_tfu
		for (size_t i = 0; i < ncontrols; i++)
		{
			colstart = i*(sumnodes)+nstates*(sumnodes + 1);
			ShessainI_V.subvec(ShessI_rowShift, ShessI_rowShift + sumnodes - 1) = (hLI_tfu(i));
			ShessI_rowShift += sumnodes;
		}
		//hI_tft0
		colstart = nstates*(sumnodes + 1) + ncontrols*sumnodes;
		
		ShessainI_V(ShessI_rowShift) = as_scalar(hLI_tft0);
		ShessI_rowShift++;

		ShessainE_V(ShessE_rowShift) = hE_tft0;//tf,t0
		ShessE_rowShift++;

		//hI_tftf
		colstart = nstates*(sumnodes + 1) + ncontrols*sumnodes + 1;
		
		ShessainI_V(ShessI_rowShift) = as_scalar(hLI_tftf);
		ShessI_rowShift++;

		ShessainE_V(ShessE_rowShift) = hE_tftf;//tf,tf
		ShessE_rowShift++;

		//////////////////////////////////////////////////////////////////////////
		// insert hI_px,pu,pt0,ptu,pp
		rowstart++;
		rowshift = nstates*(sumnodes + 1) + ncontrols*sumnodes + 2;
		for (size_t i = 0; i < nparameters; i++)
		{
			rowstart = rowstart + i;
			for (size_t j = 0; j < nstates; j++)
			{
				colstart = i*(sumnodes + 1);
				ShessainI_V.subvec(ShessI_rowShift, ShessI_rowShift + sumnodes - 1) = (hLI_px(i, j));
				ShessI_rowShift += sumnodes;

				ShessainE_V(ShessE_rowShift) = hE_px0(i, j);//p,x0
				ShessE_rowShift++;

				ShessainE_V(ShessE_rowShift) = hE_px0(i, j);//p,xf
				ShessE_rowShift++;
			}

			for (size_t j = 0; j < ncontrols; j++)
			{
				colstart = nstates*(sumnodes + 1) + sumnodes*j;
				ShessainI_V.subvec(ShessI_rowShift, ShessI_rowShift + sumnodes - 1) = (hLI_px(i, j));
				ShessI_rowShift += sumnodes;
			}
			//hI_pit0
			colstart = nstates*(sumnodes + 1) + ncontrols*sumnodes;
			
			ShessainI_V(ShessI_rowShift) = hLI_pt0(i);
			ShessI_rowShift++;

			
			ShessainE_V(ShessE_rowShift) = hE_pt0(i);//p,t0
			ShessE_rowShift++;

			//hI_pitf
			colstart = nstates*(sumnodes + 1) + ncontrols*sumnodes + 1;
			
			ShessainI_V(ShessI_rowShift) = hLI_ptf(i);
			ShessI_rowShift++;

			
			ShessainE_V(ShessE_rowShift) = hE_ptf(i);//p,xf
			ShessE_rowShift++;

			//hI_pipj
			for (size_t j = 0; j <= i; j++)
			{
				colstart = nstates*(sumnodes + 1) + ncontrols*sumnodes + 2 + j;
				
				ShessainI_V(ShessI_rowShift) = hLI_pp(i, j);
				ShessI_rowShift++;

				
				ShessainE_V(ShessE_rowShift) = hE_pp(i, j);//p,xf
				ShessE_rowShift++;
			}
		}//for i

		
		Hessian_V = join_vert(ShessainI_V, ShessainE_V);
	}

	void LpHessianCalculator::GetPhaseHessianSparsity(int iphase, umat& idependencies, vec &Hessian_I, vec &Hessian_J)
	{
		LP_DBG_START_FUN("LpHessianCalculator::GetPhaseHessianSparsity")
		int sumnodes = optpro_->GetPhase(iphase)->GetTotalNodes();
		size_t nstates = Data_->SIZES_[iphase][0];
		size_t ncontrols = Data_->SIZES_[iphase][1];
		size_t nparameters = Data_->SIZES_[iphase][2];
		size_t npaths = Data_->SIZES_[iphase][3];
		size_t nevents = Data_->SIZES_[iphase][4];

		int stateindex_start = Data_->phase_indices[iphase]->state[0] - 1;
		int stateindex_end = Data_->phase_indices[iphase]->state[Data_->phase_indices[iphase]->state.size() - 1] - 1;
		int controlindex_start = Data_->phase_indices[iphase]->control[0] - 1;
		int controlindex_end = Data_->phase_indices[iphase]->control[Data_->phase_indices[iphase]->control.size() - 1] - 1;
		int t0_index = Data_->phase_indices[iphase]->time[0] - 1;
		int tf_index = Data_->phase_indices[iphase]->time[1] - 1;
		
		mat parameters;
		if (nparameters > 0)
		{
			int para_index_start = Data_->phase_indices[iphase]->parameter[0] - 1;
			int para_index_end = Data_->phase_indices[iphase]->parameter[Data_->phase_indices[iphase]->parameter.size() - 1] - 1;
		}
		//allocate the hessian ,obtain the struct and the value

		size_t con_index_start = Data_->constraint_indices[iphase][0] - 1;
		size_t con_index_end = Data_->constraint_indices[iphase][Data_->constraint_indices[iphase].size() - 1] - 1;

		size_t diff_index_start = 0;
		size_t diff_index_end = diff_index_start + sumnodes*nstates - 1;
		size_t path_index_start = 0, path_index_end = 0, event_index_start = 0, event_index_end = 0;
		if (npaths > 0)
		{
			path_index_start = sumnodes*nstates;
			path_index_end = sumnodes*nstates + npaths*sumnodes - 1;
		}
		if (nevents > 0)
		{
			event_index_start = sumnodes*(nstates + npaths);
			event_index_end = sumnodes*(nstates + npaths) + nevents - 1;
		}
		
		int nDependancies = nnz(idependencies);
		//allocate memory for Hessain

		//get the number of nonzeros of the Lagrangian that are functions of collocation points functions
		size_t hessainI_nonzeros = sumnodes*((nDependancies - nstates - ncontrols) / 2 + nstates + ncontrols) + (2 + nparameters)*(nstates + ncontrols)*sumnodes
			+ (2 + nparameters)*(2 + nparameters - 1) / 2 + nparameters + 2;
		vec ShessainI_I(hessainI_nonzeros), ShessainI_J(hessainI_nonzeros);

		//get the number of nonzeros of the Lagrangian that are functions of endpoints functions
		size_t hessainE_nonzeros = (nstates + nstates)*(nstates + nstates - 1) / 2 + nstates + nstates +
			2 * (2 + nparameters)*(nstates)+(2 + nparameters)*(2 + nparameters - 1) / 2 + nparameters + 2;
		vec ShessainE_I(hessainE_nonzeros), ShessainE_J(hessainE_nonzeros);
		vec ShessainE_V(hessainE_nonzeros);

		//differential equations
		int rowstart, colstart;
		int ShessI_rowShift = 0, ShessE_rowShift = 0;
		vec indexvector = linspace<vec>(0, sumnodes - 1, sumnodes);
		//inset hI_xx
		for (int i = 0; i < nstates; i++)
		{
			rowstart = i*(sumnodes + 1);
			for (int j = 0; j <= i; j++)
			{
				colstart = j*(sumnodes + 1);
				if (idependencies(i, j))
				{
					// 
					ShessainI_I.subvec(ShessI_rowShift, ShessI_rowShift + sumnodes - 1) = indexvector + rowstart;
					ShessainI_J.subvec(ShessI_rowShift, ShessI_rowShift + sumnodes - 1) = indexvector + colstart;
					ShessI_rowShift += sumnodes;
				}
				ShessainE_I(ShessE_rowShift) = rowstart;//xi0,xj0
				ShessainE_J(ShessE_rowShift) = colstart;
				ShessE_rowShift++;
				if (i != j)
				{
					ShessainE_I(ShessE_rowShift) = rowstart;//xi0,xjf
					ShessainE_J(ShessE_rowShift) = colstart + sumnodes;
					ShessE_rowShift++;
				}
				ShessainE_I(ShessE_rowShift) = rowstart + sumnodes;//xif,xj0
				ShessainE_J(ShessE_rowShift) = colstart;
				ShessE_rowShift++;

				ShessainE_I(ShessE_rowShift) = rowstart + sumnodes;//xif,xjf
				ShessainE_J(ShessE_rowShift) = colstart + sumnodes;
				ShessE_rowShift++;
			}//for j
		}//for i
		//insert hI_ux.hI_uu
		int rowshift = nstates*(sumnodes + 1);
		for (int i = 0; i < ncontrols; i++)
		{
			rowstart = rowshift + i*sumnodes;
			for (int j = 0; j < nstates; j++)//ux
			{
				colstart = j*(sumnodes + 1);
				if (idependencies(i + nstates, j))
				{
					//Insert hI_ux
					ShessainI_I.subvec(ShessI_rowShift, ShessI_rowShift + sumnodes - 1) = indexvector + rowstart;
					ShessainI_J.subvec(ShessI_rowShift, ShessI_rowShift + sumnodes - 1) = indexvector + colstart;
					ShessI_rowShift += sumnodes;
				}
			}//for j

			int colshift = nstates*(sumnodes + 1);
			for (int j = 0; j <= i; j++)
			{
				colstart = colshift + j*sumnodes;
				if (idependencies(i + nstates, j + nstates))
				{
					//insert hI_uu
					ShessainI_I.subvec(ShessI_rowShift, ShessI_rowShift + sumnodes - 1) = indexvector + rowstart;
					ShessainI_J.subvec(ShessI_rowShift, ShessI_rowShift + sumnodes - 1) = indexvector + colstart;
					ShessI_rowShift += sumnodes;
				}
			}//for j
		}//for i
		//////////////////////////////////////////////////////////////////////////
		//insert hI_t0x,t0u,t0t0;
		rowshift = nstates*(sumnodes + 1) + ncontrols*sumnodes;
		rowstart = rowshift;

		// hI_t0x
		for (size_t i = 0; i < nstates; i++)
		{
			colstart = i*(sumnodes + 1);
			ShessainI_I.subvec(ShessI_rowShift, ShessI_rowShift + sumnodes - 1).fill(rowstart);
			ShessainI_J.subvec(ShessI_rowShift, ShessI_rowShift + sumnodes - 1) = indexvector + colstart;
			ShessI_rowShift += sumnodes;

			ShessainE_I(ShessE_rowShift) = rowstart;//t0,x0
			ShessainE_J(ShessE_rowShift) = colstart;
			ShessE_rowShift++;

			ShessainE_I(ShessE_rowShift) = rowstart;//t0,xf
			ShessainE_J(ShessE_rowShift) = colstart + sumnodes;
			ShessE_rowShift++;
		}
		//hI_t0u
		for (size_t i = 0; i < ncontrols; i++)
		{
			colstart = i*(sumnodes)+nstates*(sumnodes + 1);
			ShessainI_I.subvec(ShessI_rowShift, ShessI_rowShift + sumnodes - 1).fill(rowstart);
			ShessainI_J.subvec(ShessI_rowShift, ShessI_rowShift + sumnodes - 1) = indexvector + colstart;
			ShessI_rowShift += sumnodes;
		}
		//hI_t0t0
		colstart = nstates*(sumnodes + 1) + ncontrols*sumnodes;
		ShessainI_I(ShessI_rowShift) = rowstart;
		ShessainI_J(ShessI_rowShift) = colstart;
		ShessI_rowShift++;

		ShessainE_I(ShessE_rowShift) = rowstart;//t0,xf
		ShessainE_J(ShessE_rowShift) = colstart;
		ShessE_rowShift++;

		//////////////////////////////////////////////////////////////////////////
		////insert hI_tfx,tfu,tft0.tf,tf;

		rowstart++;
		// hI_tfx
		for (size_t i = 0; i < nstates; i++)
		{
			colstart = i*(sumnodes + 1);
			ShessainI_I.subvec(ShessI_rowShift, ShessI_rowShift + sumnodes - 1).fill(rowstart);
			ShessainI_J.subvec(ShessI_rowShift, ShessI_rowShift + sumnodes - 1) = indexvector + colstart;
			ShessI_rowShift += sumnodes;

			ShessainE_I(ShessE_rowShift) = rowstart;//tf,x0
			ShessainE_J(ShessE_rowShift) = colstart;
			ShessE_rowShift++;

			ShessainE_I(ShessE_rowShift) = rowstart;//tf,xf
			ShessainE_J(ShessE_rowShift) = colstart + sumnodes;
			ShessE_rowShift++;
		}
		//hI_tfu
		for (size_t i = 0; i < ncontrols; i++)
		{
			colstart = i*(sumnodes)+nstates*(sumnodes + 1);
			ShessainI_I.subvec(ShessI_rowShift, ShessI_rowShift + sumnodes - 1).fill(rowstart);
			ShessainI_J.subvec(ShessI_rowShift, ShessI_rowShift + sumnodes - 1) = indexvector + colstart;
			ShessI_rowShift += sumnodes;
		}
		//hI_tft0
		colstart = nstates*(sumnodes + 1) + ncontrols*sumnodes;
		ShessainI_I(ShessI_rowShift) = rowstart;
		ShessainI_J(ShessI_rowShift) = colstart;
		ShessI_rowShift++;

		ShessainE_I(ShessE_rowShift) = rowstart;//tf,t0
		ShessainE_J(ShessE_rowShift) = colstart;
		ShessE_rowShift++;

		//hI_tftf
		colstart = nstates*(sumnodes + 1) + ncontrols*sumnodes + 1;
		ShessainI_I(ShessI_rowShift) = rowstart;
		ShessainI_J(ShessI_rowShift) = colstart;
		ShessI_rowShift++;

		ShessainE_I(ShessE_rowShift) = rowstart;//tf,tf
		ShessainE_J(ShessE_rowShift) = colstart;
		ShessE_rowShift++;

		//////////////////////////////////////////////////////////////////////////
		// insert hI_px,pu,pt0,ptu,pp
		rowstart++;
		rowshift = nstates*(sumnodes + 1) + ncontrols*sumnodes + 2;
		for (size_t i = 0; i < nparameters; i++)
		{
			rowstart = rowstart + i;
			for (size_t j = 0; j < nstates; j++)
			{
				colstart = i*(sumnodes + 1);
				ShessainI_I.subvec(ShessI_rowShift, ShessI_rowShift + sumnodes - 1).fill(rowstart);
				ShessainI_J.subvec(ShessI_rowShift, ShessI_rowShift + sumnodes - 1) = indexvector + colstart;
				ShessI_rowShift += sumnodes;

				ShessainE_I(ShessE_rowShift) = rowstart;//p,x0
				ShessainE_J(ShessE_rowShift) = colstart;
				ShessE_rowShift++;

				ShessainE_I(ShessE_rowShift) = rowstart;//p,xf
				ShessainE_J(ShessE_rowShift) = colstart + sumnodes;
				ShessE_rowShift++;
			}

			for (size_t j = 0; j < ncontrols; j++)
			{
				colstart = nstates*(sumnodes + 1) + sumnodes*j;
				ShessainI_I.subvec(ShessI_rowShift, ShessI_rowShift + sumnodes - 1).fill(rowstart);
				ShessainI_J.subvec(ShessI_rowShift, ShessI_rowShift + sumnodes - 1) = indexvector + colstart;
				ShessI_rowShift += sumnodes;
			}
			//hI_pit0
			colstart = nstates*(sumnodes + 1) + ncontrols*sumnodes;
			ShessainI_I(ShessI_rowShift) = rowstart;
			ShessainI_J(ShessI_rowShift) = colstart;
			ShessI_rowShift++;

			ShessainE_I(ShessE_rowShift) = rowstart;//p,t0
			ShessainE_J(ShessE_rowShift) = colstart;
			ShessE_rowShift++;

			//hI_pitf
			colstart = nstates*(sumnodes + 1) + ncontrols*sumnodes + 1;
			ShessainI_I(ShessI_rowShift) = rowstart;
			ShessainI_J(ShessI_rowShift) = colstart;
			ShessI_rowShift++;

			ShessainE_I(ShessE_rowShift) = rowstart;//p,xf
			ShessainE_J(ShessE_rowShift) = colstart;
			ShessE_rowShift++;

			//hI_pipj
			for (size_t j = 0; j <= i; j++)
			{
				colstart = nstates*(sumnodes + 1) + ncontrols*sumnodes + 2 + j;
				ShessainI_I(ShessI_rowShift) = rowstart;
				ShessainI_J(ShessI_rowShift) = colstart;
				ShessI_rowShift++;

				ShessainE_I(ShessE_rowShift) = rowstart;//p,xf
				ShessainE_J(ShessE_rowShift) = colstart;
				ShessE_rowShift++;
			}
		}//for i

		Hessian_I = join_vert(ShessainI_I, ShessainE_I);
		Hessian_J = join_vert(ShessainI_J, ShessainE_J);
	}

	void LpHessianCalculator::GetHessian(double sigma, const vec&x, const vec& lambada, vec& hessain_V)
	{
		LP_DBG_START_FUN("LpHessianCalculator::GetHessian")
		//!< Get dependencies before calling Ipopt
		std::vector<umat> dependencies = Data_->allPhaseDependencies;


		//find total non-zero elements in hessian
		// note that the hessian matrix stored in lower triangular form
		int  nonZerosHessain = 0;


		for (int iphase = 0; iphase < optpro_->GetPhaseNum(); iphase++)
		{
			size_t nstates = Data_->SIZES_[iphase][0];
			size_t ncontrols = Data_->SIZES_[iphase][1];
			size_t nparameters = Data_->SIZES_[iphase][2];
			size_t npaths = Data_->SIZES_[iphase][3];
			size_t nevents = Data_->SIZES_[iphase][4];
			size_t sumnodes = Data_->PS[iphase]->Points.n_elem;
			//!< find number of non-zero elements in dependencies (always including diagonal)

			umat temJacdependencies =trans( dependencies[iphase]);
			umat temDependencies = temJacdependencies*dependencies[iphase];
			temDependencies.diag().fill( 1);
			dependencies[iphase] = temDependencies;
			int nDependancies = nnz(temDependencies);

			
			nonZerosHessain += sumnodes*((nDependancies - nstates - ncontrols) / 2 + nstates + ncontrols) + (2 + nparameters)*(nstates + ncontrols)*sumnodes
				+ (2 + nparameters)*(2 + nparameters - 1) / 2 + nparameters + 2;
			nonZerosHessain += (nstates + nstates)*(nstates + nstates - 1) / 2 + nstates + nstates +
				2 * (2 + nparameters)*nstates + (2 + nparameters)*(2 + nparameters - 1) / 2 + nparameters + 2;
		}// for iphase

		for (int ipair = 0; ipair < optpro_->GetLinkageNum(); ipair++)
		{
			shared_ptr<Linkage> temLink = optpro_->GetLinkage(ipair);
			int nstates_left = 0, ncontrol_left = 0, nparameter_left = 0, npath_left = 0, nevent_left = 0;
			optpro_->GetPhase(temLink->LeftPhase())->
				get_optimal_info(nstates_left, ncontrol_left, nparameter_left, npath_left, nevent_left);
			int nstates_right = 0, ncontrol_right = 0, nparameter_right = 0, npath_right = 0, nevent_right = 0;
			optpro_->GetPhase(temLink->LeftPhase())->
				get_optimal_info(nstates_right, ncontrol_left, nparameter_right, npath_right, nevent_right);

			int numlinks = temLink->GetLinkageMin().size();//FIX ME!!!
			size_t nonzeros_link_hessian = (nstates_left + nstates_right)*(nstates_left + nstates_right - 1) / 2 + (nstates_left + nstates_right);
			nonzeros_link_hessian += nparameter_left*(nstates_right + nstates_left) + nparameter_left*(nparameter_left - 1) + nparameter_left;
			nonzeros_link_hessian += nparameter_right*(nstates_left + nstates_right + nparameter_left) + nparameter_right*(nparameter_right - 1) + nparameter_right;
			nonZerosHessain += nonzeros_link_hessian;
		}//for ipair
		hessain_V = zeros<vec>(nonZerosHessain);
		size_t rowshift = 0;
		size_t colshift = 0;
		size_t hessian_shift = 0;
		//////////////////////////////////////////////////////////////////////////
		//Get Jac in each Phase
		for (int i = 0; i < optpro_->GetPhaseNum(); i++)
		{
			int i_nstate = optpro_->GetPhase(i)->GetstateMin().size();
			int i_ncontrol = optpro_->GetPhase(i)->GetcontrolMin().size();
			int i_nparameters = optpro_->GetPhase(i)->GetparameterMin().size();
			int i_npaths = optpro_->GetPhase(i)->GetpathMin().size();
			int i_nevents = optpro_->GetPhase(i)->GeteventMin().size();
			int i_totnodes = optpro_->GetPhase(i)->GetTotalNodes();
			vec hessian_phase_V;
			GetPhaseHessian(i, sigma, lambada, x, dependencies[i], hessian_phase_V);
			hessain_V.subvec(hessian_shift, hessian_shift + hessian_phase_V.n_elem - 1) = hessian_phase_V;
			hessian_shift += hessian_phase_V.n_elem;

			size_t	numcons = i_nstate*(i_totnodes + 1) + i_ncontrol*i_totnodes + i_nparameters + 2;
			size_t numvars = i_nstate*(i_totnodes + 1) + i_ncontrol*i_totnodes + i_nparameters + 2;
			rowshift += numcons;
			colshift += numvars;
		}
		//derive of linkage and it's Sparsity Pattern
		//int sjacRowsStart=calculateData_->conbounds_min.size()-optpro_->GetPhaseNum()-1;//count form zero
		std::vector<SolCost> solcCostTotal(optpro_->GetPhaseNum());
		for (int iphase = 0; iphase < optpro_->GetPhaseNum(); iphase++)
		{
			int sumnodes = optpro_->GetPhase(iphase)->GetTotalNodes();
			int nstates = optpro_->GetPhase(iphase)->GetstateMin().size();
			int ncontrols = optpro_->GetPhase(iphase)->GetcontrolMin().size();
			int nparameters = optpro_->GetPhase(iphase)->GetparameterMin().size();
			int npaths = optpro_->GetPhase(iphase)->GetpathMin().size();
			int nevents = optpro_->GetPhase(iphase)->GeteventMin().size();
			int stateindex_start = Data_->phase_indices[iphase]->state[0] - 1;
			int stateindex_end = Data_->phase_indices[iphase]->state[Data_->phase_indices[iphase]->state.size() - 1] - 1;
			vec state_vector = x.subvec(stateindex_start, stateindex_end);
			int controlindex_start, controlindex_end;
			vec control_vector;
			if (ncontrols > 0)
			{
				controlindex_start = Data_->phase_indices[iphase]->control[0] - 1;
				controlindex_end = Data_->phase_indices[iphase]->control[Data_->phase_indices[iphase]->control.size() - 1] - 1;
				control_vector = x.subvec(controlindex_start, controlindex_end);
			}
			int paralindex_start, paraindex_end;
			mat parameters;
			if (nparameters > 0)
			{
				paralindex_start = Data_->phase_indices[iphase]->parameter[0] - 1;
				paraindex_end = Data_->phase_indices[iphase]->parameter[Data_->phase_indices[iphase]->parameter.size() - 1] - 1;
				parameters = x.subvec(paralindex_start, paraindex_end);
			}
			int t0_index = Data_->phase_indices[iphase]->time[0] - 1;
			int tf_index = Data_->phase_indices[iphase]->time[1] - 1;
			double t0 = x(t0_index, 0);
			double tf = x(tf_index, 0);
			double tspan = tf - t0;
			vec t_radau = (Data_->PS[iphase]->Points + 1)*(tspan / 2.0) + t0;
			mat state_matrix = reshape(state_vector, sumnodes + 1, nstates);
			mat state_radau = state_matrix.rows(0, state_matrix.n_rows - 2);//index start form 0
			vec x0 = trans(state_matrix.row(0));
			vec xf = trans(state_matrix.row(state_matrix.n_rows - 1));
			mat control_radau = reshape(control_vector, sumnodes, ncontrols);

			SolCost mySolCost;
			mySolCost.initial_time_ = t0;
			mySolCost.initial_state_ = x0;
			mySolCost.terminal_time_ = tf;
			mySolCost.terminal_state_ = xf;
			mySolCost.time_ = t_radau;
			mySolCost.state_ = state_radau;
			mySolCost.control_ = control_radau;
			mySolCost.parameter_ = parameters;
			mySolCost.phase_num_ = iphase;
			solcCostTotal[iphase] = mySolCost;
		}
		rowshift = 0;
		colshift = 0;
		for (size_t ipair = 0; ipair < optpro_->GetLinkageNum(); ipair++)
		{

			vec hessian_ipair_V;
			GetLinkHessian(ipair, solcCostTotal, lambada, hessian_ipair_V);
			hessain_V.subvec(hessian_shift, hessian_shift + hessian_ipair_V.n_elem - 1) = hessian_ipair_V;
			hessian_shift += hessian_ipair_V.n_elem;
		}

	}

	void LpHessianCalculator::GetLinkHessian(int ipair, const std::vector<SolCost>solcCostTotal, const vec& lambada, vec& linkHessain_V)

	{
		LP_DBG_START_FUN("LpHessianCalculator::GetLinkHessian")
		LinkHessian ilink_hessain ;
		CalculateLinkHessain(ipair, solcCostTotal, ilink_hessain);
		int rowstart = 0, colstart = 0;
		int rowshift = 0, colshift = 0;
		int hessindexshift = 0;
		
		int nlinks = optpro_->GetLinkage(ipair)->GetLinkageMin().size();
		vec link_lambda;
		size_t link_start = Data_->link_indices[ipair][0] -1;
		size_t link_end = Data_->link_indices[ipair][Data_->link_indices[ipair].size() - 1] - 1;
		link_lambda = lambada.subvec(link_start, link_end);
		int left_index = optpro_->GetLinkage(ipair)->LeftPhase();
		int right_index = optpro_->GetLinkage(ipair)->RightPhase();
		vec xf_left = solcCostTotal[left_index].terminal_state_;
		mat p_left = solcCostTotal[left_index].parameter_;

		vec x0_right = solcCostTotal[right_index].initial_state_;
		mat p_right = solcCostTotal[right_index].parameter_;

		SolLink mySolLink;
		mySolLink.left_state_ = xf_left;
		mySolLink.left_parameter_ = p_left;
		mySolLink.left_phase_num_ = left_index + 1;

		mySolLink.right_state_ = x0_right;
		mySolLink.right_parameter_ = p_left;
		mySolLink.right_phase_num_ = right_index + 1;

		mySolLink.ipair = ipair + 1;
		//////////////////////////////////////////////////////////////////////////
		//Get Linkage Hessian
		vec Linkout;
		fun_->LinkFunction(mySolLink, Linkout);
		SolLink isolink, jsollink, ijsollink;

		int nstateLeft = mySolLink.left_state_.n_elem;
		int nparametersLeft = mySolLink.left_parameter_.n_elem;
		int nnodesLeft = optpro_->GetPhase(left_index)->GetTotalNodes();

		int nstateRight = mySolLink.right_state_.n_elem;
		int nparametersRight = mySolLink.right_parameter_.n_elem;
		int nnodesRight = optpro_->GetPhase(right_index)->GetTotalNodes();


		//////////////////////////////////////////////////////////////////////////
		size_t nonzeros_link_hessian = (nstateLeft + nstateRight)*(nstateLeft + nstateRight - 1) / 2 + (nstateLeft + nstateRight);
		nonzeros_link_hessian += nparametersLeft*(nstateRight + nstateLeft) + nparametersLeft*(nparametersLeft - 1)/2 + nparametersLeft;
		nonzeros_link_hessian += nparametersRight*(nstateLeft + nstateRight + nparametersLeft) + nparametersRight*(nparametersRight - 1) + nparametersRight;
		linkHessain_V = zeros(nonzeros_link_hessian);

		int nstatesLeft, ncontrolsLeft, npathLeft, neventLeft;
		int nstatesRight, ncontrolsRight, npathRight, neventRight;
		optpro_->GetPhase(left_index)->get_optimal_info(
			nstatesLeft, ncontrolsLeft, nparametersLeft, npathLeft, neventLeft);
		optpro_->GetPhase(right_index)->get_optimal_info(
			nstatesRight, ncontrolsRight, nparametersRight, npathRight, neventRight);
		int disc_left = nnodesLeft + 1;
		int disc_right = nnodesRight + 1;


		int stateindexstartL = Data_->phase_indices[left_index]->state[0] - 1;

		for (size_t istate = 0; istate < nstatesLeft; istate++)
		{
			rowstart = stateindexstartL + (nnodesLeft + 1)*(istate+1) - 1;
			for (size_t jstate = 0; jstate <= istate; jstate++)
			{
				colstart = stateindexstartL + (nnodesLeft + 1)*(jstate + 1) - 1;
				
				linkHessain_V(hessindexshift) = accu(ilink_hessain.hLink_xfL_xfL(istate, jstate) % link_lambda);
				hessindexshift++;
			}
		}

		//////////////////////////////////////////////////////////////////////////
		int paraindexstartL = 0;
		if (nparametersLeft > 0)
		{
			paraindexstartL = Data_->phase_indices[left_index]->parameter[0] - 1;
			for (size_t ipara = 0; ipara < nparametersLeft; ipara++)
			{
				rowstart = paraindexstartL + ipara;
				for (size_t jstate = 0; jstate < nstatesLeft; jstate++)
				{
					colstart = stateindexstartL + (nnodesLeft + 1)*(jstate + 1) - 1;
					
					linkHessain_V(hessindexshift) = accu(ilink_hessain.hLink_paraL_xfL(ipara, jstate) % link_lambda);
					hessindexshift++;
				}
				for (size_t jpara = 0; jpara <= ipara; jpara++)
				{
					colstart = paraindexstartL + jpara;
					
					linkHessain_V(hessindexshift) = accu(ilink_hessain.hLink_paraL_paraL(ipara, jpara) % link_lambda);
					hessindexshift++;
				}//for jpara
			}
		}// if nparameters
		//////////////////////////////////////////////////////////////////////////
		int stateindexstartR = Data_->phase_indices[right_index]->state[0] - 1;
		for (size_t istate = 0; istate < nstatesRight; istate++)
		{
			rowstart = stateindexstartR + (nnodesRight + 1)*istate;
			for (size_t jstate = 0; jstate < nstatesLeft; jstate++)
			{
				colstart = stateindexstartL + (nnodesLeft + 1)*(jstate + 1) - 1;
			
				linkHessain_V(hessindexshift) = accu(ilink_hessain.hLink_xfL_x0R(jstate, istate) % link_lambda);
				hessindexshift++;
			}
			for (size_t jpara = 0; jpara < nparametersLeft; jpara++)
			{
				colstart = paraindexstartL + jpara;
			
				linkHessain_V(hessindexshift) = accu(ilink_hessain.hLink_paraL_x0R(jpara, istate) % link_lambda);
				hessindexshift++;
			}
			for (size_t jstate = 0; jstate <= istate; jstate++)
			{
				colstart = stateindexstartR + (nnodesLeft + 1)*(jstate);
			
				linkHessain_V(hessindexshift) = accu(ilink_hessain.hLink_x0R_x0R(istate, jstate) % link_lambda);
				hessindexshift++;
			}
		}// for i< nstateRight
		//////////////////////////////////////////////////////////////////////////
		int paraindexstartR = 0;
		if (nparametersRight > 0)
		{
			paraindexstartR = Data_->phase_indices[right_index]->parameter[0] - 1;
			for (size_t ipara = 0; ipara < nparametersRight; ipara++)
			{
				rowstart = paraindexstartR + ipara;
				for (size_t jstate = 0; jstate < nstatesLeft; jstate++)
				{
					colstart = stateindexstartL + (nnodesLeft + 1)*(jstate + 1) - 1;
				
					linkHessain_V(hessindexshift) = accu(ilink_hessain.hLink_paraR_xfL(ipara, jstate) % link_lambda);
					hessindexshift++;
				}// xf_left
				for (size_t jpara = 0; jpara < nparametersLeft; jpara++)
				{
					colstart = paraindexstartL + jpara;
					
					linkHessain_V(hessindexshift) = accu(ilink_hessain.hLink_paraR_paraL(ipara, jpara) % link_lambda);
					hessindexshift++;
				}//for para_left

				for (size_t jstate = 0; jstate < nstatesLeft; jstate++)
				{
					colstart = stateindexstartR + (nnodesLeft + 1)*(jstate);
					
					linkHessain_V(hessindexshift) = accu(ilink_hessain.hLink_paraR_x0R(ipara, jstate) % link_lambda);
					hessindexshift++;
				}// x0_right
				for (size_t jpara = 0; jpara <= ipara; jpara++)
				{
					colstart = paraindexstartR + jpara;
					
					linkHessain_V(hessindexshift) = accu(ilink_hessain.hLink_paraR_paraR(ipara, jpara) % link_lambda);
					hessindexshift++;
				}//for para_right
			}// right_para
		}// if nparameters

	}


	void LpHessianCalculator::CalculatePhaseHessian(size_t iphase, const vec& x_all, PhaseHessain& iphase_hessian)
	{
		LP_DBG_START_FUN("LpHessianCalculator::CalculatePhaseHessian")
		int sumnodes = Data_->PS[iphase]->Points.n_elem;
		size_t nstates = Data_->SIZES_[iphase][0];
		size_t ncontrols = Data_->SIZES_[iphase][1];
		size_t nparameters = Data_->SIZES_[iphase][2];
		size_t npaths = Data_->SIZES_[iphase][3];
		size_t nevents = Data_->SIZES_[iphase][4];
		int stateindex_start = Data_->phase_indices[iphase]->state[0] - 1;
		int stateindex_end = Data_->phase_indices[iphase]->state[Data_->phase_indices[iphase]->state.size() - 1] - 1;
		int controlindex_start = Data_->phase_indices[iphase]->control[0] - 1;
		int controlindex_end = Data_->phase_indices[iphase]->control[Data_->phase_indices[iphase]->control.size() - 1] - 1;
		vec state_vector = x_all.subvec(stateindex_start, stateindex_end);
		vec control_vector = x_all.subvec(controlindex_start, controlindex_end);
		int t0_index = Data_->phase_indices[iphase]->time[0] - 1;
		int tf_index = Data_->phase_indices[iphase]->time[1] - 1;
		double t0 = x_all(t0_index);
		double tf = x_all(tf_index);
		double tspan = tf - t0;
		vec t_radau = (Data_->PS[iphase]->Points + 1)*(tspan / 2.0) + t0;
		//t_radau.print("tradau");
		mat state_matrix = reshape(state_vector, sumnodes + 1, nstates);
		//state_matrix.print("state_matrix");
		mat state_radau = state_matrix.rows(0, state_matrix.n_rows - 2);//index start form 0

		vec x0 = trans(state_matrix.row(0));
		vec xf = trans(state_matrix.row(state_matrix.n_rows - 1));
		mat control_radau = reshape(control_vector, sumnodes, ncontrols);
		//control_radau.print("control_radau");
		mat parameters;
		if (nparameters > 0)
		{
			int para_index_start = Data_->phase_indices[iphase]->parameter[0] - 1;
			int para_index_end = Data_->phase_indices[iphase]->parameter[Data_->phase_indices[iphase]->parameter.size() - 1] - 1;
			parameters = x_all.rows(para_index_start, para_index_end);
		}
		//get derivatives of DAE function
		SolDae mySolDae;
		mySolDae.time_ = t_radau;
		mySolDae.state_ = state_radau;
		mySolDae.contol_ = control_radau;
		mySolDae.parameter_ = parameters;
		mySolDae.phase_num_ = iphase + 1;

		mat stateOut, pathOut;
		fun_->DaeFunction(mySolDae, stateOut, pathOut);


		
		vec pertTime = tol*(1 + abs(mySolDae.time_));
		mat pertState = tol*(1 + abs(mySolDae.state_));
		mat pertControl, perParameter;
		if (ncontrols > 0) pertControl = tol*(1 + abs(mySolDae.contol_));
		if (nparameters > 0) perParameter = tol*(1 + abs(mySolDae.parameter_));

		
		vec tRadauPert = mySolDae.time_ + pertTime;
		mat stateRadauPert = mySolDae.state_ + pertState;
		mat controlRadauPert, parameterPert;
		if (ncontrols > 0)
		{
			controlRadauPert = mySolDae.contol_ + pertControl;
		}
		if (nparameters > 0)
		{
			parameterPert = mySolDae.parameter_ + perParameter;
		}

		mat  denominatorij;

		mat stateouti, stateoutj, stateoutij, pathouti, pathoutj, pathoutij;
		// get df/(dx_i,*dx_j) & dpath/(dx_i,dx_j)
		field<mat> hDae_state_state(nstates, nstates);
		field<mat> hPath_state_state(nstates, nstates);
		for (size_t istate = 0; istate < nstates; istate++)
		{
			SolDae isoldae = mySolDae;
			isoldae.state_.col(istate) = stateRadauPert.col(istate);
			fun_->DaeFunction(isoldae, stateouti, pathouti);//f(x+h_i)
			SolDae ijsoldae = isoldae;
			for (size_t jstate = 0; jstate < nstates; jstate++)
			{
				ijsoldae = isoldae;
				SolDae jsoldae = mySolDae;
				jsoldae.state_.col(jstate) = stateRadauPert.col(jstate);
				fun_->DaeFunction(jsoldae, stateoutj, pathoutj);//f(x+h_j)
				ijsoldae.state_.col(jstate) += pertState.col(jstate);
				fun_->DaeFunction(ijsoldae, stateoutij, pathoutij);//f(x+h_i+h_j)
				denominatorij = repmat(pertState.col(istate) % pertState.col(jstate), 1, nstates + npaths);
				hDae_state_state(istate, jstate) = (stateoutij  - stateouti - stateoutj+ stateOut) / denominatorij.head_cols(nstates);
				if (npaths > 0)
				{
					hPath_state_state(istate, jstate) = (pathoutij  - pathouti - pathoutj+ pathOut) / denominatorij.tail_cols(npaths);

				}

			}//jstate

		}//istate
		iphase_hessian.hDae_state_state = hDae_state_state;
		iphase_hessian.hPath_state_state = hPath_state_state;

			field<mat> hDae_control_state(ncontrols, nstates);
			field<mat> hPath_control_state(ncontrols, nstates);
			field<mat> hDae_control_control(ncontrols, ncontrols);
			field<mat> hPath_control_control(ncontrols, ncontrols);

			iphase_hessian.hDae_control_state = hDae_control_state;
			iphase_hessian.hPath_control_state = hPath_control_state;
			iphase_hessian.hDae_control_control = hDae_control_control;
			iphase_hessian.hPath_control_control = hPath_control_control;
		// get df/(du_i,*dx_j) & dpath/(du_i,dx_j)
		for (size_t icontrol = 0; icontrol < ncontrols; icontrol++)
		{
			SolDae isoldae = mySolDae;
			isoldae.contol_.col(icontrol) = controlRadauPert.col(icontrol);
			fun_->DaeFunction(isoldae, stateouti, pathouti);//f(x+h_i)
			SolDae ijsoldae = isoldae;
			for (size_t jstate = 0; jstate < nstates; jstate++)
			{
				ijsoldae = isoldae;
				SolDae jsoldae = mySolDae;
				jsoldae.state_.col(jstate) = stateRadauPert.col(jstate);
				fun_->DaeFunction(jsoldae, stateoutj, pathoutj);//f(x+h_j)
				ijsoldae.state_.col(jstate) += pertState.col(jstate);
				fun_->DaeFunction(ijsoldae, stateoutij, pathoutij);//f(x+h_i+h_j)
				denominatorij = repmat(pertControl.col(icontrol) % pertState.col(jstate), 1, nstates + npaths);
				///////////////////////////////////////////
				// tem = (stateoutij + stateOut - stateouti - stateoutj);
				// //////////////////////////////////////
				iphase_hessian.hDae_control_state(icontrol, jstate) = (stateoutij  - stateouti - stateoutj+ stateOut) / denominatorij.head_cols(nstates);
				if (npaths > 0)
				{
					iphase_hessian.hPath_control_state(icontrol, jstate) = (pathoutij  - pathouti - pathoutj+ pathOut) / denominatorij.tail_cols(npaths);
				}

			}//jstate
			
			for (size_t jcontrol = 0; jcontrol < ncontrols; jcontrol++)
			{
				ijsoldae = isoldae;
				SolDae jsoldae = mySolDae;
				jsoldae.contol_.col(jcontrol) = controlRadauPert.col(jcontrol);
				fun_->DaeFunction(jsoldae, stateoutj, pathoutj);//f(x+h_j)
				ijsoldae.contol_.col(jcontrol) += pertControl.col(jcontrol);
				fun_->DaeFunction(ijsoldae, stateoutij, pathoutij);//f(x+h_i+h_j)
				denominatorij = repmat(pertControl.col(icontrol) % pertControl.col(jcontrol), 1, nstates + npaths);
				iphase_hessian.hDae_control_control(icontrol, jcontrol) = (stateoutij  - stateouti - stateoutj+ stateOut) / denominatorij.head_cols(nstates);
				if (npaths > 0)
				{
					iphase_hessian.hPath_control_control(icontrol, jcontrol) = (pathoutij  - pathouti - pathoutj+ pathOut) / denominatorij.tail_cols(npaths);
				}
			}
		}//jcontrol
		

		field<mat> hDae_time_state(nstates);
		field<mat> hPath_time_state(nstates);
		field<mat> hDae_time_control(ncontrols);
		field<mat> hPath_time_control(ncontrols);

		iphase_hessian.hDae_time_state = hDae_time_state;
		iphase_hessian.hPath_time_state = hPath_time_state;
		iphase_hessian.hDae_time_control = hDae_time_control;
		iphase_hessian.hPath_time_control = hPath_time_control;
		SolDae tsoldae = mySolDae;
		tsoldae.time_ = tRadauPert;
		fun_->DaeFunction(tsoldae, stateouti, pathouti);
		SolDae tjsoldae = tsoldae;
		//get df/(dt*dx)
		for (size_t jstate = 0; jstate < nstates; jstate++)
		{
			tjsoldae = tsoldae;
			SolDae jsoldae = mySolDae;
			jsoldae.state_.col(jstate) = stateRadauPert.col(jstate);
			fun_->DaeFunction(jsoldae, stateoutj, pathoutj);//f(x+h_j)
			tjsoldae.state_.col(jstate) += pertState.col(jstate);
			fun_->DaeFunction(tjsoldae, stateoutij, pathoutij);//f(x+h_i+h_j)
			denominatorij = repmat(pertTime % pertState.col(jstate), 1, nstates + npaths);
			iphase_hessian.hDae_time_state(jstate) = (stateoutij  - stateouti - stateoutj+ stateOut) / denominatorij.head_cols(nstates);
			if (npaths > 0)
			{
				iphase_hessian.hPath_time_state(jstate) = (pathoutij  - pathouti - pathoutj+ pathOut) / denominatorij.tail_cols(npaths);
			}

		}//jstate
		//get df/(dt*du)
		
		for (size_t jcontrol = 0; jcontrol < ncontrols; jcontrol++)
		{
			tjsoldae = tsoldae;
			SolDae jsoldae = mySolDae;
			jsoldae.contol_.col(jcontrol) = controlRadauPert.col(jcontrol);
			fun_->DaeFunction(jsoldae, stateoutj, pathoutj);//f(x+h_j)
			tjsoldae.contol_.col(jcontrol) += pertControl.col(jcontrol);
			fun_->DaeFunction(tjsoldae, stateoutij, pathoutij);//f(x+h_i+h_j)
			denominatorij = repmat(pertTime % pertControl.col(jcontrol), 1, nstates + npaths);
			iphase_hessian.hDae_time_control(jcontrol) = (stateoutij  - stateouti - stateoutj+ stateOut) / denominatorij.head_cols(nstates);
			if (npaths > 0)
			{
				iphase_hessian.hPath_time_control(jcontrol) = (pathoutij  - pathouti - pathoutj+ pathOut) / denominatorij.tail_cols(npaths);
			}

		}//jstate
		//get df/(dt*dt)
		mat hdae_time_time, hPath_time_time;
		tjsoldae = tsoldae;
		tjsoldae.time_ += pertTime;
		fun_->DaeFunction(tjsoldae, stateoutij, pathoutij);//f(x+h_i+h_j)
		denominatorij = repmat(pertTime % pertTime, 1, nstates + npaths);
		stateoutj = stateouti;
		hdae_time_time = (stateoutij  - stateouti - stateoutj+ stateOut) / denominatorij.head_cols(nstates);
		if (npaths > 0)
		{
			pathoutj = pathouti;
			hPath_time_time = (pathoutij  - pathouti - pathoutj+ pathOut) / denominatorij.tail_cols(npaths);
		}
		iphase_hessian.hdae_time_time = hdae_time_time;
		iphase_hessian.hPath_time_time = hPath_time_time;
		//get df/(dpi*dxj) &get df/(dpi*duj)
		field<mat> hDae_para_state(nparameters, nstates);
		field<mat> hPath_para_state(nparameters, nstates);
		field<mat> hDae_para_control(nparameters, ncontrols);
		field<mat> hPath_para_control(nparameters, ncontrols);
		field<mat> hDae_para_time(nparameters);
		field<mat> hPath_para_time(nparameters);
		field<mat> hDae_para_para(nparameters, nparameters);
		field<mat> hPath_para_para(nparameters, nparameters);
		iphase_hessian.hDae_para_state = hDae_para_state;
		iphase_hessian.hPath_para_state = hPath_para_state;
		iphase_hessian.hDae_para_control = hDae_para_control;
		iphase_hessian.hPath_para_control = hPath_para_control;
		iphase_hessian.hDae_para_time = hDae_para_time;
		iphase_hessian.hPath_para_time = hPath_para_time;
		iphase_hessian.hDae_para_para = hDae_para_para;
		iphase_hessian.hPath_para_para = hPath_para_para;
		for (size_t iparameter = 0; iparameter < nparameters; iparameter++)
		{
			SolDae isoldae = mySolDae;
			isoldae.parameter_.row(iparameter) = parameterPert.row(iparameter);
			fun_->DaeFunction(isoldae, stateouti, pathouti);//f(x+h_i)
			SolDae ijsoldae = isoldae;
			for (size_t jstate = 0; jstate < nstates; jstate++)
			{
				ijsoldae = isoldae;
				SolDae jsoldae = mySolDae;
				jsoldae.state_.col(jstate) = stateRadauPert.col(jstate);
				fun_->DaeFunction(jsoldae, stateoutj, pathoutj);//f(x+h_j)
				ijsoldae.state_.col(jstate) += pertState.col(jstate);
				fun_->DaeFunction(ijsoldae, stateoutij, pathoutij);//f(x+h_i+h_j)
				denominatorij = repmat(perParameter(iparameter) * pertState.col(jstate), 1, nstates + npaths);

				iphase_hessian.hDae_para_state(iparameter, jstate) = (stateoutij  - stateouti - stateoutj+ stateOut) / denominatorij.head_cols(nstates);
				if (npaths > 0)
				{
					iphase_hessian.hPath_para_state(iparameter, jstate) = (pathoutij  - pathouti - pathoutj+ pathOut) / denominatorij.tail_cols(npaths);
				}

			}//jstate
			
			for (size_t jcontrol = 0; jcontrol < ncontrols; jcontrol++)
			{
				ijsoldae = isoldae;
				SolDae jsoldae = mySolDae;
				jsoldae.contol_.col(jcontrol) = controlRadauPert.col(jcontrol);
				fun_->DaeFunction(jsoldae, stateoutj, pathoutj);//f(x+h_j)
				ijsoldae.contol_.col(jcontrol) += pertControl.col(jcontrol);
				fun_->DaeFunction(ijsoldae, stateoutij, pathoutij);//f(x+h_i+h_j)
				denominatorij = repmat(perParameter(iparameter)* pertControl.col(jcontrol), 1, nstates + npaths);
				iphase_hessian.hDae_para_control(iparameter, jcontrol) = (stateoutij  - stateouti - stateoutj+ stateOut) / denominatorij.head_cols(nstates);
				if (npaths > 0)
				{
					iphase_hessian.hPath_para_control(iparameter, jcontrol) = (pathoutij  - pathouti - pathoutj+ pathOut) / denominatorij.tail_cols(npaths);
				}
			}//jcontrol
			//get df/(dpi*dt)
			ijsoldae = isoldae;
			SolDae jsoldae = mySolDae;
			jsoldae.time_ = tRadauPert;
			fun_->DaeFunction(jsoldae, stateoutj, pathoutj);//f(x+h_j)
			ijsoldae.time_ += pertTime;
			fun_->DaeFunction(ijsoldae, stateoutij, pathoutij);//f(x+h_i+h_j)
			denominatorij = repmat(pertTime * perParameter(iparameter), 1, nstates + npaths);
			iphase_hessian.hDae_para_time(iparameter) = (stateoutij  - stateouti - stateoutj+ stateOut) / denominatorij.head_cols(nstates);
			if (npaths > 0)
			{
				iphase_hessian.hPath_para_time(iparameter) = (pathoutij  - pathouti - pathoutj+ pathOut) / denominatorij.tail_cols(npaths);
			}
			//get df/(dpi*dpj)
			
			for (size_t jparameter = 0; jparameter < nparameters; jparameter++)
			{
				ijsoldae = isoldae;
				SolDae jsoldae = mySolDae;
				jsoldae.parameter_.row(jparameter) = parameterPert.row(jparameter);
				fun_->DaeFunction(jsoldae, stateoutj, pathoutj);//f(x+h_j)
				ijsoldae.parameter_.row(jparameter) += perParameter.row(jparameter);
				fun_->DaeFunction(ijsoldae, stateoutij, pathoutij);//f(x+h_i+h_j)
				denominatorij = perParameter(iparameter)*perParameter(jparameter)*ones(sumnodes, nstates + npaths);
				iphase_hessian.hDae_para_para(iparameter, jparameter) = (stateoutij  - stateouti - stateoutj+ stateOut) / denominatorij.head_cols(nstates);
				if (npaths > 0)
				{
					iphase_hessian.hPath_para_para(iparameter, jparameter) = (pathoutij - pathouti - pathoutj + pathOut) / denominatorij.tail_cols(npaths);
				}
			}

		}//for ipara


		//Get hessian of event function
		vec events, eventsi, eventsj, eventsij;


		if (nevents > 0)
		{
			field<mat> hEvents_x0_x0(nstates,nstates);
			field<mat> hEvents_x0_xf(nstates,nstates);
			field<mat> hEvents_xf_x0(nstates,nstates);
			field<mat> hEvents_xf_xf(nstates,nstates);
			field<mat> hEvents_t0_x0(nstates);
			field<mat> hEvents_t0_xf(nstates);
			mat hEvents_t0_t0;
			field<mat> hEvents_tf_x0(nstates);
			field<mat> hEvents_tf_xf(nstates);
			mat hEvents_tf_t0;
			mat hEvents_tf_tf;
			field<mat>hEvents_para_x0(nparameters,nstates);
			field<mat>hEvents_para_xf(nparameters,nstates);
			field<mat>hEvents_para_t0(nparameters);
			field<mat>hEvents_para_tf(nparameters);
			field<mat>hEvents_para_para(nparameters,nparameters);

			iphase_hessian.hEvents_x0_x0 = hEvents_x0_x0;
			iphase_hessian.hEvents_x0_xf = hEvents_x0_xf;
			iphase_hessian.hEvents_xf_x0 = hEvents_xf_x0;
			iphase_hessian.hEvents_xf_xf = hEvents_xf_xf;
			iphase_hessian.hEvents_t0_x0 = hEvents_t0_x0;
			iphase_hessian.hEvents_t0_xf = hEvents_t0_xf;
			iphase_hessian.hEvents_t0_t0 = hEvents_t0_t0;
			iphase_hessian.hEvents_tf_x0 = hEvents_tf_x0;
			iphase_hessian.hEvents_tf_xf = hEvents_tf_xf;
			iphase_hessian.hEvents_tf_t0 = hEvents_tf_t0;
			iphase_hessian.hEvents_tf_tf = hEvents_tf_tf;
			iphase_hessian.hEvents_tf_t0 = hEvents_tf_t0;
			iphase_hessian.hEvents_tf_tf = hEvents_tf_tf;
			iphase_hessian.hEvents_para_x0 = hEvents_para_x0;
			iphase_hessian.hEvents_para_xf = hEvents_para_xf;
			iphase_hessian.hEvents_para_t0 = hEvents_para_t0;
			iphase_hessian.hEvents_para_tf = hEvents_para_tf;
			iphase_hessian.hEvents_para_para = hEvents_para_para;
			SolEvent mySolEvent, isolevent, jsolevent, ijsolevent;
			mySolEvent.initial_time_ = t0;
			mySolEvent.initial_state_ = x0;
			mySolEvent.terminal_time_ = tf;
			mySolEvent.terminal_state_ = xf;
			mySolEvent.parameter_ = parameters;
			mySolEvent.phase_num_ = iphase + 1;
			events = zeros(nevents, 1);
			fun_->EventFunction(mySolEvent, events);

			double pert0 = tol*(1+abs(mySolEvent.initial_time_));
			double pertf = tol*(1 + abs(mySolEvent.terminal_time_));
			mat pertx0 = tol*(1 + abs(mySolEvent.initial_state_));
			mat pertxf = tol*(1 + abs(mySolEvent.terminal_state_));

			double t0Pert = mySolEvent.initial_time_ + pert0;
			vec x0Pert = mySolEvent.initial_state_ + pertx0;
			double tfPert = mySolEvent.terminal_time_ + pertf;
			vec xfPert = mySolEvent.terminal_state_ + pertxf;
			vec x0 = mySolEvent.initial_state_;
			vec xf = mySolEvent.terminal_state_;


			for (size_t istate = 0; istate < nstates; istate++)
			{
				isolevent = mySolEvent;
				isolevent.initial_state_(istate) = x0Pert(istate);
				fun_->EventFunction(isolevent, eventsi);

				for (size_t jstate = 0; jstate < nstates; jstate++)
				{
					ijsolevent = isolevent;
					jsolevent = mySolEvent;
					jsolevent.initial_state_(jstate) = x0Pert(jstate);
					fun_->EventFunction(jsolevent, eventsj);
					ijsolevent.initial_state_(jstate) += pertx0(jstate);
					fun_->EventFunction(ijsolevent, eventsij);
					iphase_hessian.hEvents_x0_x0(istate, jstate) = (eventsij  - eventsi - eventsj+ events) / (pertx0(istate)*pertx0(jstate)*ones<vec>(nevents));

					ijsolevent = isolevent;
					jsolevent = mySolEvent;
					jsolevent.terminal_state_(jstate) = xfPert(jstate);
					fun_->EventFunction(jsolevent, eventsj);
					ijsolevent.terminal_state_(jstate) += pertxf(jstate);
					fun_->EventFunction(ijsolevent, eventsij);
					iphase_hessian.hEvents_x0_xf(istate, jstate) = (eventsij  - eventsi - eventsj+ events) / (pertx0(istate)*pertxf(istate)*ones<vec>(nevents));
				}



				isolevent = mySolEvent;
				isolevent.terminal_state_(istate) = xfPert(istate);
				fun_->EventFunction(isolevent, eventsi);
				for (size_t jstate = 0; jstate < nstates; jstate++)
				{
					ijsolevent = isolevent;
					jsolevent = mySolEvent;
					jsolevent.initial_state_(jstate) = x0Pert(jstate);
					fun_->EventFunction(jsolevent, eventsj);
					ijsolevent.initial_state_(jstate) += pertx0(jstate);
					fun_->EventFunction(ijsolevent, eventsij);
					iphase_hessian.hEvents_xf_x0(istate, jstate) = (eventsij  - eventsi - eventsj+ events) / (pertxf(istate)*pertx0(jstate)*ones<vec>(nevents));

					ijsolevent = isolevent;
					jsolevent = mySolEvent;
					jsolevent.terminal_state_(jstate) = xfPert(jstate);
					fun_->EventFunction(jsolevent, eventsj);
					ijsolevent.terminal_state_(jstate) += pertxf(jstate);
					fun_->EventFunction(ijsolevent, eventsij);
					iphase_hessian.hEvents_xf_xf(istate, jstate) = (eventsij  - eventsi - eventsj+ events) / (pertxf(istate)*pertxf(istate)*ones<vec>(nevents));
				}
			}//istate
			//Get hessian of Event with Respect to Initial Time
			isolevent = mySolEvent;
			isolevent.initial_time_ = t0Pert;
			fun_->EventFunction(isolevent, eventsi);
			for (size_t jstate = 0; jstate < nstates; jstate++)
			{
				ijsolevent = isolevent;
				jsolevent = mySolEvent;
				jsolevent.initial_state_(jstate) = x0Pert(jstate);
				fun_->EventFunction(jsolevent, eventsj);
				ijsolevent.initial_state_(jstate) += pertx0(jstate);
				fun_->EventFunction(ijsolevent, eventsij);
				iphase_hessian.hEvents_t0_x0(jstate) = (eventsij  - eventsi - eventsj+ events) / (pert0*pertx0(jstate)*ones<vec>(nevents));

				ijsolevent = isolevent;
				jsolevent = mySolEvent;
				jsolevent.terminal_state_(jstate) = xfPert(jstate);
				fun_->EventFunction(jsolevent, eventsj);
				ijsolevent.terminal_state_(jstate) += pertxf(jstate);
				fun_->EventFunction(ijsolevent, eventsij);
				iphase_hessian.hEvents_t0_xf(jstate) = (eventsij  - eventsi - eventsj+ events) / (pert0*pertxf(jstate)*ones<vec>(nevents));
			}

			isolevent = mySolEvent;
			isolevent.initial_time_ = t0Pert;
			fun_->EventFunction(isolevent, eventsi);
			ijsolevent = isolevent;
			jsolevent = mySolEvent;
			jsolevent.initial_time_ = t0Pert;
			fun_->EventFunction(isolevent, eventsj);
			ijsolevent.initial_time_ += pert0;
			fun_->EventFunction(ijsolevent, eventsij);
			iphase_hessian.hEvents_t0_t0 = (eventsij  - eventsi - eventsj+ events) / (pert0*pert0*ones<vec>(nevents));

			//Get hessian of Event with Respect to Terminal Time
			isolevent = mySolEvent;
			isolevent.terminal_time_ = tfPert;
			fun_->EventFunction(isolevent, eventsi);
			for (size_t jstate = 0; jstate < nstates; jstate++)
			{
				ijsolevent = isolevent;
				jsolevent = mySolEvent;
				jsolevent.initial_state_(jstate) = x0Pert(jstate);
				fun_->EventFunction(jsolevent, eventsj);
				ijsolevent.initial_state_(jstate) += pertx0(jstate);
				fun_->EventFunction(ijsolevent, eventsij);
				iphase_hessian.hEvents_tf_x0(jstate) = (eventsij  - eventsi - eventsj+ events) / (pertf*pertx0(jstate)*ones<vec>(nevents));

				ijsolevent = isolevent;
				jsolevent = mySolEvent;
				jsolevent.terminal_state_(jstate) = xfPert(jstate);
				fun_->EventFunction(jsolevent, eventsj);
				ijsolevent.terminal_state_(jstate) += pertxf(jstate);
				fun_->EventFunction(ijsolevent, eventsij);
				iphase_hessian.hEvents_tf_xf(jstate) = (eventsij  - eventsi - eventsj+ events) / (pertf*pertxf(jstate)*ones<vec>(nevents));
			}

			isolevent = mySolEvent;
			isolevent.terminal_time_ = tfPert;
			fun_->EventFunction(isolevent, eventsi);
			ijsolevent = isolevent;
			jsolevent = mySolEvent;
			jsolevent.initial_time_ = t0Pert;
			fun_->EventFunction(jsolevent, eventsj);
			ijsolevent.initial_time_ += pert0;
			fun_->EventFunction(ijsolevent, eventsij);
			iphase_hessian.hEvents_tf_t0 = (eventsij  - eventsi - eventsj+ events) / (pertf*pert0*ones<vec>(nevents));

			isolevent = mySolEvent;
			isolevent.terminal_time_ = tfPert;
			fun_->EventFunction(isolevent, eventsi);
			ijsolevent = isolevent;
			jsolevent = mySolEvent;
			jsolevent.terminal_time_ = tfPert;
			fun_->EventFunction(jsolevent, eventsj);
			ijsolevent.terminal_time_ += pertf;
			fun_->EventFunction(ijsolevent, eventsij);
			iphase_hessian.hEvents_tf_tf = (eventsij  - eventsi - eventsj+ events) / (pertf*pertf*ones<vec>(nevents));

			//Get hessian of Event with Respect to parameters
			for (size_t ipara = 0; ipara < nparameters; ipara++)
			{
				isolevent = mySolEvent;
				isolevent.parameter_(ipara) = parameterPert(ipara);
				fun_->EventFunction(isolevent, eventsi);
				for (size_t jstate = 0; jstate < nstates; jstate++)
				{
					ijsolevent = isolevent;
					jsolevent = mySolEvent;
					jsolevent.initial_state_(jstate) = x0Pert(jstate);
					fun_->EventFunction(jsolevent, eventsj);
					ijsolevent.initial_state_(jstate) += pertx0(jstate);
					fun_->EventFunction(ijsolevent, eventsij);
					iphase_hessian.hEvents_para_x0(ipara, jstate) = (eventsij  - eventsi - eventsj+ events) / (perParameter(ipara)*pertx0(jstate)*ones<vec>(nevents));

					ijsolevent = isolevent;
					jsolevent = mySolEvent;
					jsolevent.terminal_state_(jstate) = xfPert(jstate);
					fun_->EventFunction(jsolevent, eventsj);
					ijsolevent.terminal_state_(jstate) += pertxf(jstate);
					fun_->EventFunction(ijsolevent, eventsij);
					iphase_hessian.hEvents_para_xf(ipara, jstate) = (eventsij  - eventsi - eventsj+ events) / (perParameter(ipara)*pertxf(jstate)*ones<vec>(nevents));
				}
				//devent/dpara*dt0

				ijsolevent = isolevent;
				jsolevent = mySolEvent;
				jsolevent.initial_time_ = t0Pert;
				fun_->EventFunction(jsolevent, eventsj);
				ijsolevent.initial_time_ += pert0;
				fun_->EventFunction(ijsolevent, eventsij);
				iphase_hessian.hEvents_para_t0(ipara) = (eventsij  - eventsi - eventsj+ events) / (perParameter(ipara)*pert0*ones<vec>(nevents));
				//devent/dpara*dtf

				ijsolevent = isolevent;
				jsolevent = mySolEvent;
				jsolevent.terminal_time_ = tfPert;
				fun_->EventFunction(jsolevent, eventsj);
				ijsolevent.terminal_time_ += pertf;
				fun_->EventFunction(ijsolevent, eventsij);
				iphase_hessian.hEvents_para_tf(ipara) = (eventsij  - eventsi - eventsj+ events) / (perParameter(ipara)*pertf*ones<vec>(nevents));
				//devent/dpara*dpara
				for (size_t jpara = 0; jpara < nparameters; jpara++)
				{
					ijsolevent = isolevent;
					jsolevent = mySolEvent;
					jsolevent.parameter_(jpara) = parameterPert(jpara);
					fun_->EventFunction(jsolevent, eventsj);
					ijsolevent.parameter_ (jpara)+= perParameter(jpara);
					fun_->EventFunction(ijsolevent, eventsij);
					iphase_hessian.hEvents_para_para(ipara,jpara) = (eventsij  - eventsi - eventsj+ events) / (perParameter(ipara)*perParameter(jpara)*ones<vec>(nevents));
				}
			}//for inpara
		}//if nevent

		//Get hessian of Mayer 
		SolCost mysolcost, isolcost, jsolcost, ijsolcost;
		mysolcost.initial_time_ = t0;
		mysolcost.initial_state_ = x0;
		mysolcost.terminal_time_ = tf;
		mysolcost.terminal_state_ = xf;
		mysolcost.time_ = t_radau;
		mysolcost.state_ = state_radau;
		mysolcost.control_ = control_radau;
		mysolcost.parameter_ = parameters;
		mysolcost.phase_num_ = iphase + 1;
		double mayer, mayeri, mayerj, mayerij;
		fun_->MayerCost(mysolcost, mayer);

		//Get hessain of Mayer with Respect to Initial and Terminal State
		mat hMayer_x0_x0(nstates,nstates);
		mat hMayer_x0_xf(nstates,nstates);
		mat hMayer_xf_x0(nstates,nstates);
		mat hMayer_xf_xf(nstates,nstates);
		vec hMayer_t0_x0(nstates);
		vec hMayer_t0_xf(nstates);

		vec hMayer_tf_x0(nstates);
		vec hMayer_tf_xf(nstates);

		mat hMayer_para_x0(nparameters,nstates);
		mat hMayer_para_xf(nparameters,nstates);
		vec hMayer_para_t0(nparameters);
		vec hMayer_para_tf(nparameters);
		mat hMayer_para_para(nparameters,nparameters);

		iphase_hessian.hMayer_x0_x0 = hMayer_x0_x0;
		iphase_hessian.hMayer_x0_xf = hMayer_x0_xf;
		iphase_hessian.hMayer_xf_x0 = hMayer_xf_x0;
		iphase_hessian.hMayer_xf_xf = hMayer_xf_xf;
		iphase_hessian.hMayer_t0_x0 = hMayer_t0_x0;
		iphase_hessian.hMayer_t0_xf = hMayer_t0_xf;
		iphase_hessian.hMayer_tf_x0 = hMayer_tf_x0;
		iphase_hessian.hMayer_tf_xf = hMayer_tf_xf;
		iphase_hessian.hMayer_para_x0 = hMayer_para_x0;
		iphase_hessian.hMayer_para_xf = hMayer_para_xf;
		iphase_hessian.hMayer_para_t0 = hMayer_para_t0;
		iphase_hessian.hMayer_para_tf = hMayer_para_tf;
		iphase_hessian.hMayer_para_para = hMayer_para_para;

		double pert0 = tol*(1 + abs(mysolcost.initial_time_));
		double pertf = tol*(1 + abs(mysolcost.terminal_time_));
		mat pertx0 = tol*(1 + abs(mysolcost.initial_state_));
		mat pertxf = tol*(1 + abs(mysolcost.terminal_state_));

		double t0Pert = mysolcost.initial_time_ + pert0;
		vec x0Pert = mysolcost.initial_state_ + pertx0;
		double tfPert = mysolcost.terminal_time_ + pertf;
		vec xfPert = mysolcost.terminal_state_ + pertxf;
		for (size_t istate = 0; istate < nstates; istate++)
		{
			isolcost = mysolcost;
			isolcost.initial_state_(istate) = x0Pert(istate);
			fun_->MayerCost(isolcost, mayeri);

			for (size_t jstate = 0; jstate < nstates; jstate++)
			{
				ijsolcost = isolcost;
				jsolcost = mysolcost;
				jsolcost.initial_state_(jstate) = x0Pert(jstate);
				fun_->MayerCost(jsolcost, mayerj);
				ijsolcost.initial_state_(jstate) += pertx0(jstate);
				fun_->MayerCost(ijsolcost, mayerij);
				iphase_hessian.hMayer_x0_x0(istate, jstate) = (mayerij  - mayeri - mayerj+ mayer) / (pertx0(istate)*pertx0(jstate));

				ijsolcost = isolcost;
				jsolcost = mysolcost;
				jsolcost.terminal_state_(jstate) = xfPert(jstate);
				fun_->MayerCost(jsolcost, mayerj);
				ijsolcost.terminal_state_(jstate) += pertxf(jstate);
				fun_->MayerCost(ijsolcost, mayerij);
				iphase_hessian.hMayer_x0_xf(istate, jstate) = (mayerij  - mayeri - mayerj+ mayer) / (pertx0(istate)*pertxf(istate));
			}



			isolcost = mysolcost;
			isolcost.terminal_state_(istate) = xfPert(istate);
			fun_->MayerCost(isolcost, mayeri);
			for (size_t jstate = 0; jstate < nstates; jstate++)
			{
				ijsolcost = isolcost;
				jsolcost = mysolcost;
				jsolcost.initial_state_(jstate) = x0Pert(jstate);
				fun_->MayerCost(jsolcost, mayerj);
				ijsolcost.initial_state_(jstate) += pertx0(jstate);
				fun_->MayerCost(ijsolcost, mayerij);
				iphase_hessian.hMayer_xf_x0(istate, jstate) = (mayerij  - mayeri - mayerj+ mayer) / (pertxf(istate)*pertx0(jstate));

				ijsolcost = isolcost;
				jsolcost = mysolcost;
				jsolcost.terminal_state_(jstate) = xfPert(jstate);
				fun_->MayerCost(jsolcost, mayerj);
				ijsolcost.terminal_state_(jstate) += pertxf(jstate);
				fun_->MayerCost(ijsolcost, mayerij);
				iphase_hessian.hMayer_xf_xf(istate, jstate) = (mayerij  - mayeri - mayerj+ mayer) / (pertxf(istate)*pertxf(istate));
			}
		}//istate
		//Get hessian of Mayer with Respect to Initial Time
		isolcost = mysolcost;
		isolcost.initial_time_ = t0Pert;
		fun_->MayerCost(isolcost, mayeri);
		for (size_t jstate = 0; jstate < nstates; jstate++)
		{
			ijsolcost = isolcost;
			jsolcost = mysolcost;
			jsolcost.initial_state_(jstate) = x0Pert(jstate);
			fun_->MayerCost(jsolcost, mayerj);
			ijsolcost.initial_state_(jstate) += pertx0(jstate);
			fun_->MayerCost(ijsolcost, mayerij);
			iphase_hessian.hMayer_t0_x0(jstate) = (mayerij  - mayeri - mayerj+ mayer) / (pert0*pertx0(jstate));

			ijsolcost = isolcost;
			jsolcost = mysolcost;
			jsolcost.terminal_state_(jstate) = xfPert(jstate);
			fun_->MayerCost(jsolcost, mayerj);
			ijsolcost.terminal_state_(jstate) += pertxf(jstate);
			fun_->MayerCost(ijsolcost, mayerij);
			iphase_hessian.hMayer_t0_xf(jstate) = (mayerij  - mayeri - mayerj+ mayer) / (pert0*pertxf(jstate));
		}

		isolcost = mysolcost;
		isolcost.initial_time_ = t0Pert;
		fun_->MayerCost(isolcost, mayeri);
		ijsolcost = isolcost;
		jsolcost = mysolcost;
		jsolcost.initial_time_ = t0Pert;
		fun_->MayerCost(jsolcost, mayerj);
		ijsolcost.initial_time_ += pert0;
		fun_->MayerCost(ijsolcost, mayerij);
		iphase_hessian.hMayer_t0_t0 = (mayerij  - mayeri - mayerj+ mayer) / (pert0*pert0);

		//Get hessian of Mayer with Respect to Terminal Time
		isolcost = mysolcost;
		isolcost.terminal_time_ = tfPert;
		fun_->MayerCost(isolcost, mayeri);
		for (size_t jstate = 0; jstate < nstates; jstate++)
		{
			ijsolcost = isolcost;
			jsolcost = mysolcost;
			jsolcost.initial_state_(jstate) = x0Pert(jstate);
			fun_->MayerCost(jsolcost, mayerj);
			ijsolcost.initial_state_(jstate) += pertx0(jstate);
			fun_->MayerCost(ijsolcost, mayerij);
			iphase_hessian.hMayer_tf_x0(jstate) = (mayerij  - mayeri - mayerj+ mayer) / (pertf*pertx0(jstate));

			ijsolcost = isolcost;
			jsolcost = mysolcost;
			jsolcost.terminal_state_(jstate) = xfPert(jstate);
			fun_->MayerCost(jsolcost, mayerj);
			ijsolcost.terminal_state_(jstate) += pertxf(jstate);
			fun_->MayerCost(ijsolcost, mayerij);
			iphase_hessian.hMayer_tf_xf(jstate) = (mayerij  - mayeri - mayerj+ mayer) / (pertf*pertxf(jstate));
		}

		isolcost = mysolcost;
		isolcost.terminal_time_ = tfPert;
		fun_->MayerCost(isolcost, mayeri);
		ijsolcost = isolcost;
		jsolcost = mysolcost;
		jsolcost.initial_time_ = t0Pert;
		fun_->MayerCost(jsolcost, mayerj);
		ijsolcost.initial_time_ += pert0;
		fun_->MayerCost(ijsolcost, mayerij);
		iphase_hessian.hMayer_tf_t0 = (mayerij  - mayeri - mayerj+ mayer) / (pertf*pert0);

		isolcost = mysolcost;
		isolcost.terminal_time_ = tfPert;
		fun_->MayerCost(isolcost, mayeri);
		ijsolcost = isolcost;
		jsolcost = mysolcost;
		jsolcost.terminal_time_ = tfPert;
		fun_->MayerCost(jsolcost, mayerj);
		ijsolcost.terminal_time_ += pertf;
		fun_->MayerCost(ijsolcost, mayerij);
		iphase_hessian.hMayer_tf_tf = (mayerij  - mayeri - mayerj+ mayer) / (pertf*pertf);
		//Get hessian of Mayer with Respect to parameters
		for (size_t ipara = 0; ipara < nparameters; ipara++)
		{
			isolcost = mysolcost;
			isolcost.parameter_(ipara) = parameterPert(ipara);
			fun_->MayerCost(isolcost, mayeri);
			for (size_t jstate = 0; jstate < nstates; jstate++)
			{
				ijsolcost = isolcost;
				jsolcost = mysolcost;
				jsolcost.initial_state_(jstate) = x0Pert(jstate);
				fun_->MayerCost(jsolcost, mayerj);
				ijsolcost.initial_state_(jstate) += pertx0(jstate);
				fun_->MayerCost(ijsolcost, mayerij);
				iphase_hessian.hMayer_para_x0(ipara, jstate) = (mayerij  - mayeri - mayerj+ mayer) / (perParameter(ipara)*pertx0(jstate));

				ijsolcost = isolcost;
				jsolcost = mysolcost;
				jsolcost.terminal_state_(jstate) = xfPert(jstate);
				fun_->MayerCost(jsolcost, mayerj);
				ijsolcost.terminal_state_(jstate) += pertxf(jstate);
				fun_->MayerCost(ijsolcost, mayerij);
				iphase_hessian.hMayer_para_xf(ipara, jstate) = (mayerij  - mayeri - mayerj+ mayer) / (perParameter(ipara)*pertxf(jstate));
			}
			//devent/dpara*dt0

			ijsolcost = isolcost;
			jsolcost = mysolcost;
			jsolcost.initial_time_ = t0Pert;
			fun_->MayerCost(jsolcost, mayerj);
			ijsolcost.initial_time_ += pert0;
			fun_->MayerCost(ijsolcost, mayerij);
			iphase_hessian.hMayer_para_t0(ipara) = (mayerij  - mayeri - mayerj+ mayer) / (perParameter(ipara)*pert0);
			//devent/dpara*dtf

			ijsolcost = isolcost;
			jsolcost = mysolcost;
			jsolcost.terminal_time_ = tfPert;
			fun_->MayerCost(jsolcost, mayerj);
			ijsolcost.terminal_time_ += pertf;
			fun_->MayerCost(ijsolcost, mayerij);
			iphase_hessian.hMayer_para_tf(ipara) = (mayerij  - mayeri - mayerj+ mayer) / (perParameter(ipara)*pertf);
			//devent/dpara*dpara
			for (size_t jpara = 0; jpara < nparameters; jpara++)
			{
				ijsolcost = isolcost;
				jsolcost = mysolcost;
				jsolcost.parameter_(jpara) = parameterPert(jpara);
				fun_->MayerCost(jsolcost, mayerj);
				ijsolcost.parameter_ (jpara)+= perParameter(jpara);
				fun_->MayerCost(ijsolcost, mayerij);
				iphase_hessian.hMayer_para_para(ipara,jpara) = (mayerij  - mayeri - mayerj+ mayer) / (perParameter(ipara)*perParameter(jpara));
			}
		}//for ipara
		// Get hessian of Lagrange

		field<vec> hLagrange_x_x(nstates,nstates);
		field<vec> hLagrange_u_x(ncontrols,nstates);
		field<vec> hLagrange_u_u(ncontrols,ncontrols);
		field<vec> hLagrange_t_x(nstates);
		field<vec> hLagrange_t_u(ncontrols);
		vec        hLagrange_t_t;
		field<vec> hLagrange_p_x(nparameters,nstates);
		field<vec> hLagrange_p_u(nparameters,ncontrols);
		field<vec> hLagrange_p_t(nparameters);
		field<vec> hLagrange_p_p(nparameters,nparameters);

		iphase_hessian.hLagrange_x_x = hLagrange_x_x;
		iphase_hessian.hLagrange_u_x = hLagrange_u_x;
		iphase_hessian.hLagrange_u_u = hLagrange_u_u;
		iphase_hessian.hLagrange_t_x = hLagrange_t_x;
		iphase_hessian.hLagrange_t_u = hLagrange_t_u;
		iphase_hessian.hLagrange_t_t = hLagrange_t_t;
		iphase_hessian.hLagrange_p_x = hLagrange_p_x;
		iphase_hessian.hLagrange_p_u = hLagrange_p_u;
		iphase_hessian.hLagrange_p_t = hLagrange_p_t;
		iphase_hessian.hLagrange_p_p = hLagrange_p_p;
		vec lagrange, lagrangei, lagrangej, lagrangeij;
		fun_->LagrangeCost(mysolcost, lagrange);

		// get dL/(dx_i,*dx_j) 
		for (size_t istate = 0; istate < nstates; istate++)
		{
			SolCost isolcost = mysolcost;
			isolcost.state_.col(istate) = stateRadauPert.col(istate);
			fun_->LagrangeCost(isolcost, lagrangei);//f(x+h_i)
			SolCost ijsolcost = isolcost;
			for (size_t jstate = 0; jstate < nstates; jstate++)
			{
				ijsolcost = isolcost;
				SolCost jsolcost = mysolcost;
				jsolcost.state_.col(jstate) = stateRadauPert.col(jstate);
				fun_->LagrangeCost(jsolcost, lagrangej);//f(x+h_j)
				ijsolcost.state_.col(jstate) += pertState.col(jstate);
				fun_->LagrangeCost(ijsolcost, lagrangeij);//f(x+h_i+h_j)
				denominatorij = pertState.col(istate) % pertState.col(jstate);
				iphase_hessian.hLagrange_x_x(istate, jstate) = (lagrangeij  - lagrangei - lagrangej+ lagrange) / denominatorij;
			}//jstate

		}//istate
		// get dL/(du_i,*dx_j) & dL/(du_i,du_j)
		for (size_t icontrol = 0; icontrol < ncontrols; icontrol++)
		{
			SolCost isolcost = mysolcost;
			isolcost.control_.col(icontrol) = controlRadauPert.col(icontrol);
			fun_->LagrangeCost(isolcost, lagrangei);//f(x+h_i)
			SolCost ijsolcost = isolcost;
			for (size_t jstate = 0; jstate < nstates; jstate++)
			{
				ijsolcost = isolcost;
				SolCost jsolcost = mysolcost;
				jsolcost.state_.col(jstate) = stateRadauPert.col(jstate);
				fun_->LagrangeCost(jsolcost, lagrangej);//f(x+h_j)
				ijsolcost.state_.col(jstate) += pertState.col(jstate);
				fun_->LagrangeCost(ijsolcost, lagrangeij);//f(x+h_i+h_j)
				denominatorij = pertControl.col(icontrol) % pertState.col(jstate);
				iphase_hessian.hLagrange_u_x(icontrol, jstate) = (lagrangeij  - lagrangei - lagrangej+ lagrange) / denominatorij;
			}//jstate
			
			for (size_t jcontrol = 0; jcontrol < ncontrols; jcontrol++)
			{
				ijsolcost = isolcost;
				SolCost jsolcost = mysolcost;
				jsolcost.control_.col(jcontrol) = controlRadauPert.col(jcontrol);
				fun_->LagrangeCost(jsolcost, lagrangej);//f(x+h_j)
				ijsolcost.control_.col(jcontrol) += pertControl.col(jcontrol);
				fun_->LagrangeCost(ijsolcost, lagrangeij);//f(x+h_i+h_j)
				denominatorij = pertControl.col(icontrol) % pertControl.col(jcontrol);
				iphase_hessian.hLagrange_u_u(icontrol, jcontrol) = (lagrangeij  - lagrangei - lagrangej+ lagrange) / denominatorij;

			}
		}//jcontrol

		SolCost tsolcost = mysolcost;
		tsolcost.time_ = tRadauPert;
		fun_->LagrangeCost(tsolcost, lagrangei);
		SolCost tjsolcost = tsolcost;
		//get dL/(dt*dx)
		for (size_t jstate = 0; jstate < nstates; jstate++)
		{
			tjsolcost = tsolcost;
			SolCost jsolcost = mysolcost;
			jsolcost.state_.col(jstate) = stateRadauPert.col(jstate);
			fun_->LagrangeCost(jsolcost, lagrangej);//f(x+h_j)
			tjsolcost.state_.col(jstate) += pertState.col(jstate);
			fun_->LagrangeCost(tjsolcost, lagrangeij);//f(x+h_i+h_j)
			denominatorij = pertTime % pertState.col(jstate);
			iphase_hessian.hLagrange_t_x(jstate) = (lagrangeij  - lagrangei - lagrangej+ lagrange) / denominatorij;

		}//jstate
		//get dL/(dt*du)
		
		for (size_t jcontrol = 0; jcontrol < ncontrols; jcontrol++)
		{
			tjsolcost = tsolcost;
			SolCost jsolcost = mysolcost;
			jsolcost.control_.col(jcontrol) = controlRadauPert.col(jcontrol);
			fun_->LagrangeCost(jsolcost, lagrangej);//f(x+h_j)
			tjsolcost.control_.col(jcontrol) += pertControl.col(jcontrol);
			fun_->LagrangeCost(tjsolcost, lagrangeij);//f(x+h_i+h_j)
			denominatorij = pertTime % pertControl.col(jcontrol);
			iphase_hessian.hLagrange_t_u(jcontrol) = (lagrangeij  - lagrangei - lagrangej+ lagrange) / denominatorij;
		}//jstate
		//get dL/(dt*dt)
		tjsolcost = tsolcost;
		tjsolcost.time_ += pertTime;
		lagrangej = lagrangei;
		fun_->LagrangeCost(tjsolcost, lagrangeij);//f(x+h_i+h_j)
		denominatorij = pertTime % pertTime;

		iphase_hessian.hLagrange_t_t = (lagrangeij  - lagrangei - lagrangej+ lagrange) / denominatorij;

		//get dL/(dpi*dxj) &get dL/(dpi*duj)
		for (size_t iparameter = 0; iparameter < nparameters; iparameter++)
		{
			SolCost isolcost = mysolcost;
			isolcost.parameter_(iparameter) = parameterPert(iparameter);
			fun_->LagrangeCost(isolcost, lagrangei);//f(x+h_i)
			SolCost ijsolcost = isolcost;
			for (size_t jstate = 0; jstate < nstates; jstate++)
			{
				ijsolcost = isolcost;
				SolCost jsolcost = mysolcost;
				jsolcost.state_.col(jstate) = stateRadauPert.col(jstate);
				fun_->LagrangeCost(jsolcost, lagrangej);//f(x+h_j)
				ijsolcost.state_.col(jstate) += pertState.col(jstate);
				fun_->LagrangeCost(ijsolcost, lagrangeij);//f(x+h_i+h_j)
				denominatorij = (perParameter(iparameter) * pertState.col(jstate));

				iphase_hessian.hLagrange_p_x(iparameter, jstate) = (lagrangeij  - lagrangei - lagrangej+ lagrange) / denominatorij;

			}//jstate
			
			for (size_t jcontrol = 0; jcontrol < ncontrols; jcontrol++)
			{
				ijsolcost = isolcost;
				SolCost jsolcost = mysolcost;
				jsolcost.control_.col(jcontrol) = controlRadauPert.col(jcontrol);
				fun_->LagrangeCost(jsolcost, lagrangej);//f(x+h_j)
				ijsolcost.control_.col(jcontrol) += pertControl.col(jcontrol);
				fun_->LagrangeCost(ijsolcost, lagrangeij);//f(x+h_i+h_j)
				denominatorij = (perParameter(iparameter)* pertControl.col(jcontrol));
				iphase_hessian.hLagrange_p_u(iparameter, jcontrol) = (lagrangeij  - lagrangei - lagrangej+ lagrange) / denominatorij;

			}//jcontrol
			//get dL/(dpi*dt)
			ijsolcost = isolcost;
			SolCost jsolcost = mysolcost;
			jsolcost.time_ = tRadauPert;
			fun_->LagrangeCost(jsolcost, lagrangej);//f(x+h_j)
			ijsolcost.time_ += pertTime;
			fun_->LagrangeCost(ijsolcost, lagrangeij);//f(x+h_i+h_j)
			denominatorij = (pertTime * perParameter(iparameter));
			iphase_hessian.hLagrange_p_t(iparameter) = (lagrangeij  - lagrangei - lagrangej+ lagrange) / denominatorij;

			//get dL/(dpi*dpj)
			
			for (size_t jparameter = 0; jparameter < nparameters; jparameter++)
			{
				ijsolcost = isolcost;
				SolCost jsolcost = mysolcost;
				jsolcost.parameter_(jparameter) = parameterPert(jparameter);
				fun_->LagrangeCost(jsolcost, lagrangej);//f(x+h_j)
				ijsolcost.parameter_(jparameter) += perParameter(jparameter);
				fun_->LagrangeCost(ijsolcost, lagrangeij);//f(x+h_i+h_j)
				denominatorij = perParameter(iparameter)*perParameter(jparameter)*ones(sumnodes, 1);
				iphase_hessian.hLagrange_p_p(iparameter, jparameter) = (lagrangeij  - lagrangei - lagrangej+ lagrange) / denominatorij;
			}

		}//for ipara
	}

	void LpHessianCalculator::CalculateLinkHessain(size_t ipair, const std::vector<SolCost>solcCostTotal, LinkHessian& ilink_hessian)
	{
		LP_DBG_START_FUN("LpHessianCalculator::CalculateLinkHessain")
		int nlinks = optpro_->GetLinkage(ipair)->GetLinkageMin().size();
		int left_index = optpro_->GetLinkage(ipair)->LeftPhase();
		int right_index = optpro_->GetLinkage(ipair)->RightPhase();
		vec xf_left = solcCostTotal[left_index].terminal_state_;
		mat p_left = solcCostTotal[left_index].parameter_;

		vec x0_right = solcCostTotal[right_index].initial_state_;
		mat p_right = solcCostTotal[right_index].parameter_;

		SolLink mySolLink;
		mySolLink.left_state_ = xf_left;
		mySolLink.left_parameter_ = p_left;
		mySolLink.left_phase_num_ = left_index + 1;

		mySolLink.right_state_ = x0_right;
		mySolLink.right_parameter_ = p_right;
		mySolLink.right_phase_num_ = right_index + 1;

		mySolLink.ipair = ipair + 1;
		//////////////////////////////////////////////////////////////////////////
		//Get Linkage Hessian
		vec Linkout;
		fun_->LinkFunction(mySolLink, Linkout);
		SolLink isolink, jsollink, ijsollink;

		int nstateLeft = mySolLink.left_state_.n_elem;
		int nparametersLeft = mySolLink.left_parameter_.n_elem;

		int nstateRight = mySolLink.right_state_.n_elem;
		int nparametersRight = mySolLink.right_parameter_.n_elem;

		vec pertxfLeft = tol*(1 + abs(xf_left));
		vec pertpLeft = tol*(1 + abs(p_left));
		mat pertx0Right = tol*(1 + abs(x0_right));
		mat pertpRight = tol*(1 + abs(p_right));

		mat xfLeftPert = xf_left + pertxfLeft;
		mat pLeftPert;
		if (nparametersLeft>0)pLeftPert= p_left + pertpLeft;
		mat x0RightPert = x0_right + pertx0Right;
		mat pRightPert;
		if (nparametersRight>0)pRightPert = p_right + pertpRight;
		mat perLinkouti, perLinkoutj, perLinkoutij;

		field<vec> hLink_xfL_xfL(nstateLeft,nstateLeft);
		field<vec> hLink_xfL_x0R(nstateLeft,nstateRight);
		field<vec>hLink_paraL_xfL(nparametersLeft,nstateLeft);
		field<vec>hLink_paraL_paraL(nparametersLeft,nparametersLeft);
		field<vec>hLink_paraL_x0R(nparametersLeft,nstateRight);
		field<vec> hLink_x0R_x0R(nstateRight,nstateRight);
		field<vec>hLink_paraR_xfL(nparametersRight,nstateLeft);
		field<vec>hLink_paraR_paraL(nparametersRight,nparametersLeft);
		field<vec>hLink_paraR_x0R(nparametersRight,nstateRight);
		field<vec>hLink_paraR_paraR(nparametersRight,nparametersRight);
		ilink_hessian.hLink_xfL_xfL = hLink_xfL_xfL;
		ilink_hessian.hLink_xfL_x0R = hLink_xfL_x0R;
		ilink_hessian.hLink_x0R_x0R = hLink_x0R_x0R;
		ilink_hessian.hLink_paraL_xfL = hLink_paraL_xfL;
		ilink_hessian.hLink_paraL_x0R = hLink_paraL_x0R;
		ilink_hessian.hLink_paraL_paraL = hLink_paraL_paraL;
		ilink_hessian.hLink_paraR_xfL = hLink_paraR_xfL;
		ilink_hessian.hLink_paraR_x0R = hLink_paraR_x0R;
		ilink_hessian.hLink_paraR_paraL = hLink_paraR_paraL;
		ilink_hessian.hLink_paraR_paraR = hLink_paraR_paraR;
		//Get hLink_xfl_xfl,xfL_x0R;
		vec linkouti, linkoutj, linkoutij;
		for (size_t istate = 0; istate < nstateLeft; istate++)
		{
			isolink = mySolLink;
			isolink.left_state_(istate)=(xfLeftPert(istate));
			fun_->LinkFunction(isolink, linkouti);
			for (size_t jstate = 0; jstate < nstateLeft; jstate++)
			{
				ijsollink = isolink;
				jsollink = mySolLink;
				jsollink.left_state_(jstate)=(xfLeftPert(jstate));
				fun_->LinkFunction(jsollink, linkoutj);
				ijsollink.left_state_(jstate) += pertxfLeft(jstate);
				fun_->LinkFunction(ijsollink, linkoutij);
				ilink_hessian.hLink_xfL_xfL(istate, jstate) = (linkoutij  - linkouti - linkoutj+ Linkout) / (pertxfLeft(istate)*pertxfLeft(jstate));
			}

			for (size_t jstate = 0; jstate < nstateRight; jstate++)
			{
				ijsollink = isolink;
				jsollink = mySolLink;
				jsollink.right_state_(jstate)=(x0RightPert(jstate));
				fun_->LinkFunction(jsollink, linkoutj);
				ijsollink.right_state_(jstate) += pertx0Right(jstate);
				fun_->LinkFunction(ijsollink, linkoutij);
				ilink_hessian.hLink_xfL_x0R(istate, jstate) = (linkoutij  - linkouti - linkoutj+ Linkout) / (pertxfLeft(istate)*pertx0Right(jstate));
			}
		}// for istate
		//////////////////////////////////////////////////////////////////////////
		// get hLink_paraL_xfL,hLink_paraL_x0R,hLink_paraL_paraL
		for (size_t ipara = 0; ipara < nparametersLeft; ipara++)
		{
			isolink = mySolLink;
			isolink.left_parameter_(ipara) = (pLeftPert(ipara));
			fun_->LinkFunction(isolink, linkouti);
			for (size_t jstate = 0; jstate < nstateLeft; jstate++)
			{
				ijsollink = isolink;
				jsollink = mySolLink;
				jsollink.left_state_(jstate) = xfLeftPert(jstate);
				fun_->LinkFunction(jsollink, linkoutj);
				ijsollink.left_state_(jstate) += pertxfLeft(jstate);
				fun_->LinkFunction(ijsollink, linkoutij);
				ilink_hessian.hLink_paraL_xfL(ipara, jstate) = (linkoutij  - linkouti - linkoutj+ Linkout) / (pertpLeft(ipara)*pertxfLeft(jstate));
			}
			for (size_t jstate = 0; jstate < nstateRight; jstate++)
			{
				ijsollink = isolink;
				jsollink = mySolLink;
				jsollink.right_state_(jstate) = x0RightPert(jstate);
				fun_->LinkFunction(jsollink, linkoutj);
				ijsollink.right_state_(jstate) += pertx0Right(jstate);
				fun_->LinkFunction(ijsollink, linkoutij);
				ilink_hessian.hLink_paraL_x0R(ipara, jstate) = (linkoutij  - linkouti - linkoutj+ Linkout) / (pertpLeft(ipara)*pertx0Right(jstate));
			}
			for (size_t jpara = 0; jpara < nparametersLeft; jpara++)
			{
				ijsollink = isolink;
				jsollink = mySolLink;
				jsollink.left_parameter_(jpara) = pLeftPert(jpara);
				fun_->LinkFunction(jsollink, linkoutj);
				ijsollink.left_parameter_(jpara) += pertpLeft(jpara);
				fun_->LinkFunction(ijsollink, linkoutij);
				ilink_hessian.hLink_paraL_paraL(ipara, jpara) = (linkoutij  - linkouti - linkoutj+ Linkout) / (pertpLeft(ipara)*pertpLeft(jpara));
			}

		}
		//////////////////////////////////////////////////////////////////////////
		//get hLink_x0R_x0R
		for (size_t istate = 0; istate < nstateRight; istate++)
		{
			isolink = mySolLink;
			isolink.right_state_(istate)=(x0RightPert(istate));
			fun_->LinkFunction(isolink, linkouti);
			for (size_t jstate = 0; jstate < nstateRight; jstate++)
			{
				ijsollink = isolink;
				jsollink = mySolLink;
				jsollink.right_state_(jstate)=(x0RightPert(jstate));
				fun_->LinkFunction(jsollink, linkoutj);
				ijsollink.right_state_(jstate) += pertx0Right(jstate);
				fun_->LinkFunction(ijsollink, linkoutij);
				ilink_hessian.hLink_x0R_x0R(istate, jstate) = (linkoutij  - linkouti - linkoutj+ Linkout) / (pertx0Right(istate)*pertx0Right(jstate));
			}
		}// for istate
		//////////////////////////////////////////////////////////////////////////
		// get hLink_paraR_xfL,hLink_paraR_x0R,hLink_paraR_paraL,hLink_paraR_paraR

		for (size_t ipara = 0; ipara < nparametersRight; ipara++)
		{
			isolink = mySolLink;
			isolink.right_parameter_(ipara) = (pRightPert(ipara));
			fun_->LinkFunction(isolink, linkouti);
			for (size_t jstate = 0; jstate < nstateLeft; jstate++)
			{
				ijsollink = isolink;
				jsollink = mySolLink;
				jsollink.left_state_(jstate) = xfLeftPert(jstate);
				fun_->LinkFunction(jsollink, linkoutj);
				ijsollink.left_state_(jstate) += pertxfLeft(jstate);
				fun_->LinkFunction(ijsollink, linkoutij);
				ilink_hessian.hLink_paraR_xfL(ipara, jstate) = (linkoutij  - linkouti - linkoutj+ Linkout) / (pertpRight(ipara)*pertxfLeft(jstate));
			}
			for (size_t jstate = 0; jstate < nstateRight; jstate++)
			{
				ijsollink = isolink;
				jsollink = mySolLink;
				jsollink.right_state_(jstate) = x0RightPert(jstate);
				fun_->LinkFunction(jsollink, linkoutj);
				ijsollink.right_state_(jstate) += pertx0Right(jstate);
				fun_->LinkFunction(ijsollink, linkoutij);
				ilink_hessian.hLink_paraR_x0R(ipara, jstate) = (linkoutij  - linkouti - linkoutj+ Linkout) / (pertpRight(ipara)*pertx0Right(jstate));
			}
			for (size_t jpara = 0; jpara < nparametersLeft; jpara++)
			{
				ijsollink = isolink;
				jsollink = mySolLink;
				jsollink.left_parameter_(jpara) = pLeftPert(jpara);
				fun_->LinkFunction(jsollink, linkoutj);
				ijsollink.left_parameter_(jpara) += pertpLeft(jpara);
				fun_->LinkFunction(ijsollink, linkoutij);
				ilink_hessian.hLink_paraR_paraL(ipara, jpara) = (linkoutij  - linkouti - linkoutj+ Linkout) / (pertpRight(ipara)*pertpLeft(jpara));
			}
			for (size_t jpara = 0; jpara < nparametersRight; jpara++)
			{
				ijsollink = isolink;
				jsollink = mySolLink;
				jsollink.right_parameter_(jpara) = pRightPert(jpara);
				fun_->LinkFunction(jsollink, linkoutj);
				ijsollink.right_parameter_(jpara) += pertpRight(jpara);
				fun_->LinkFunction(ijsollink, linkoutij);
				ilink_hessian.hLink_paraR_paraR(ipara, jpara) = (linkoutij  - linkouti - linkoutj+ Linkout) / (pertpRight(ipara)*pertpRight(jpara));
			}


		}
	}

	void LpHessianCalculator::GetLinkHessianSparsity(int ipair, vec& linkHessain_I, vec& linkHessain_J)
	{
		LP_DBG_START_FUN("LpHessianCalculator::GetLinkHessianSparsity")
		int rowstart = 0, colstart = 0;
		int rowshift = 0, colshift = 0;
		int hessindexshift = 0;

		int nlinks = optpro_->GetLinkage(ipair)->GetLinkageMin().size();
		vec link_lambda;
		size_t link_start = Data_->link_indices[ipair][0] - 1;
		size_t link_end = Data_->link_indices[ipair][Data_->link_indices[ipair].size() - 1] - 1;
		int left_index = optpro_->GetLinkage(ipair)->LeftPhase();
		int right_index = optpro_->GetLinkage(ipair)->RightPhase();


		int nstatesLeft, ncontrolsLeft, npathLeft, neventLeft, nparali, nparametersLeft;
		int nstatesRight, ncontrolsRight, npathRight, neventRight, nparametersRight;
		optpro_->GetPhase(left_index)->get_optimal_info(
			nstatesLeft, ncontrolsLeft, nparametersLeft, npathLeft, neventLeft);
		optpro_->GetPhase(right_index)->get_optimal_info(
			nstatesRight, ncontrolsRight, nparametersRight, npathRight, neventRight);
		
		int nnodesLeft = optpro_->GetPhase(left_index)->GetTotalNodes();
		int nnodesRight = optpro_->GetPhase(right_index)->GetTotalNodes();
		int disc_left = nnodesLeft + 1;
		int disc_right = nnodesRight + 1;

		//////////////////////////////////////////////////////////////////////////
		size_t nonzeros_link_hessian = (nstatesLeft + nstatesRight)*(nstatesLeft + nstatesRight - 1) / 2 + (nstatesLeft + nstatesRight);
		nonzeros_link_hessian += nparametersLeft*(nstatesRight + nstatesLeft) + nparametersLeft*(nparametersLeft - 1) / 2 + nparametersLeft;
		nonzeros_link_hessian += nparametersRight*(nstatesLeft + nstatesRight + nparametersLeft) + nparametersRight*(nparametersRight - 1) + nparametersRight;
		linkHessain_I = zeros(nonzeros_link_hessian);
		linkHessain_J = zeros(nonzeros_link_hessian);

		int stateindexstartL = Data_->phase_indices[left_index]->state[0] - 1;

		for (size_t istate = 0; istate < nstatesLeft; istate++)
		{
			rowstart = stateindexstartL + (nnodesLeft + 1)*(istate + 1) - 1;
			for (size_t jstate = 0; jstate <= istate; jstate++)
			{
				colstart = stateindexstartL + (nnodesLeft + 1)*(jstate + 1) - 1;
				linkHessain_I(hessindexshift) = rowstart;
				linkHessain_J(hessindexshift) = colstart;
				hessindexshift++;
			}
		}

		//////////////////////////////////////////////////////////////////////////
		int paraindexstartL = 0;
		if (nparametersLeft > 0)
		{
			paraindexstartL = Data_->phase_indices[left_index]->parameter[0] - 1;
			for (size_t ipara = 0; ipara < nparametersLeft; ipara++)
			{
				rowstart = paraindexstartL + ipara;
				for (size_t jstate = 0; jstate < nstatesLeft; jstate++)
				{
					colstart = stateindexstartL + (nnodesLeft + 1)*(jstate + 1) - 1;
					linkHessain_I(hessindexshift) = rowstart;
					linkHessain_J(hessindexshift) = colstart;
					hessindexshift++;
				}
				for (size_t jpara = 0; jpara <= ipara; jpara++)
				{
					colstart = paraindexstartL + jpara;
					linkHessain_I(hessindexshift) = rowstart;
					linkHessain_J(hessindexshift) = colstart;
					hessindexshift++;
				}//for jpara
			}
		}// if nparameters
		//////////////////////////////////////////////////////////////////////////
		int stateindexstartR = Data_->phase_indices[right_index]->state[0] - 1;
		for (size_t istate = 0; istate < nstatesRight; istate++)
		{
			rowstart = stateindexstartR + (nnodesRight + 1)*istate;
			for (size_t jstate = 0; jstate < nstatesLeft; jstate++)
			{
				colstart = stateindexstartL + (nnodesLeft + 1)*(jstate + 1) - 1;
				linkHessain_I(hessindexshift) = rowstart;
				linkHessain_J(hessindexshift) = colstart;
				hessindexshift++;
			}
			for (size_t jpara = 0; jpara < nparametersLeft; jpara++)
			{
				colstart = paraindexstartL + jpara;
				linkHessain_I(hessindexshift) = rowstart;
				linkHessain_J(hessindexshift) = colstart;

				hessindexshift++;
			}
			for (size_t jstate = 0; jstate <= istate; jstate++)
			{
				colstart = stateindexstartR + (nnodesLeft + 1)*(jstate);
				linkHessain_I(hessindexshift) = rowstart;
				linkHessain_J(hessindexshift) = colstart;
				hessindexshift++;
			}
		}// for i< nstateRight
		//////////////////////////////////////////////////////////////////////////
		int paraindexstartR = 0;
		if (nparametersRight > 0)
		{
			paraindexstartR = Data_->phase_indices[right_index]->parameter[0] - 1;
			for (size_t ipara = 0; ipara < nparametersRight; ipara++)
			{
				rowstart = paraindexstartR + ipara;
				for (size_t jstate = 0; jstate < nstatesLeft; jstate++)
				{
					colstart = stateindexstartL + (nnodesLeft + 1)*(jstate + 1) - 1;
					linkHessain_I(hessindexshift) = rowstart;
					linkHessain_J(hessindexshift) = colstart;
					hessindexshift++;
				}// xf_left
				for (size_t jpara = 0; jpara < nparametersLeft; jpara++)
				{
					colstart = paraindexstartL + jpara;
					linkHessain_I(hessindexshift) = rowstart;
					linkHessain_J(hessindexshift) = colstart;
					hessindexshift++;
				}//for para_left

				for (size_t jstate = 0; jstate < nstatesLeft; jstate++)
				{
					colstart = stateindexstartR + (nnodesLeft + 1)*(jstate);
					linkHessain_I(hessindexshift) = rowstart;
					linkHessain_J(hessindexshift) = colstart;
					hessindexshift++;
				}// x0_right
				for (size_t jpara = 0; jpara <= ipara; jpara++)
				{
					colstart = paraindexstartR + jpara;
					linkHessain_I(hessindexshift) = rowstart;
					linkHessain_J(hessindexshift) = colstart;
					hessindexshift++;
				}//for para_right
			}// right_para
		}// if nparameters
	}

	void LpHessianCalculator::GetHessianSparsity(vec& hessain_I, vec& hessain_J)
	{
		LP_DBG_START_FUN("LpHessianCalculator::GetHessianSparsity")
		//!< Get dependencies before calling Ipopt
		std::vector<umat> dependencies = Data_->allPhaseDependencies;


		//find total non-zero elements in hessian
		// note that the hessian matrix stored in lower triangular form
		int  nonZerosHessain = 0;


		for (int iphase = 0; iphase < optpro_->GetPhaseNum(); iphase++)
		{
			size_t nstates = Data_->SIZES_[iphase][0];
			size_t ncontrols = Data_->SIZES_[iphase][1];
			size_t nparameters = Data_->SIZES_[iphase][2];
			size_t npaths = Data_->SIZES_[iphase][3];
			size_t nevents = Data_->SIZES_[iphase][4];
			size_t sumnodes = Data_->PS[iphase]->Points.n_elem;
			//!< find number of non-zero elements in dependencies (always including diagonal)

			umat temJacdependencies = trans(dependencies[iphase]);
			umat temDependencies = temJacdependencies*dependencies[iphase];
			temDependencies.diag().fill(1);
			dependencies[iphase] = temDependencies;
			int nDependancies = nnz(temDependencies);


			nonZerosHessain += sumnodes*((nDependancies - nstates - ncontrols) / 2 + nstates + ncontrols) + (2 + nparameters)*(nstates + ncontrols)*sumnodes
				+ (2 + nparameters)*(2 + nparameters - 1) / 2 + nparameters + 2;
			nonZerosHessain += (nstates + nstates)*(nstates + nstates - 1) / 2 + nstates + nstates +
				2 * (2 + nparameters)*nstates + (2 + nparameters)*(2 + nparameters - 1) / 2 + nparameters + 2;
		}// for iphase

		for (int ipair = 0; ipair < optpro_->GetLinkageNum(); ipair++)
		{
			shared_ptr<Linkage> temLink = optpro_->GetLinkage(ipair);
			int nstates_left = 0, ncontrol_left = 0, nparameter_left = 0, npath_left = 0, nevent_left = 0;
			optpro_->GetPhase(temLink->LeftPhase())->
				get_optimal_info(nstates_left, ncontrol_left, nparameter_left, npath_left, nevent_left);
			int nstates_right = 0, ncontrol_right = 0, nparameter_right = 0, npath_right = 0, nevent_right = 0;
			optpro_->GetPhase(temLink->LeftPhase())->
				get_optimal_info(nstates_right, ncontrol_left, nparameter_right, npath_right, nevent_right);

			int numlinks = temLink->GetLinkageMin().size();//FIX ME!!!
			size_t nonzeros_link_hessian = (nstates_left + nstates_right)*(nstates_left + nstates_right - 1) / 2 + (nstates_left + nstates_right);
			nonzeros_link_hessian += nparameter_left*(nstates_right + nstates_left) + nparameter_left*(nparameter_left - 1) + nparameter_left;
			nonzeros_link_hessian += nparameter_right*(nstates_left + nstates_right + nparameter_left) + nparameter_right*(nparameter_right - 1) + nparameter_right;
			nonZerosHessain += nonzeros_link_hessian;
		}//for ipair
		hessain_I = zeros<vec>(nonZerosHessain);
		hessain_J = zeros<vec>(nonZerosHessain);

		size_t rowshift = 0;
		size_t colshift = 0;
		size_t hessian_shift = 0;
		//////////////////////////////////////////////////////////////////////////
		//Get Jac in each Phase
		for (int i = 0; i < optpro_->GetPhaseNum(); i++)
		{
			int i_nstate = optpro_->GetPhase(i)->GetstateMin().size();
			int i_ncontrol = optpro_->GetPhase(i)->GetcontrolMin().size();
			int i_nparameters = optpro_->GetPhase(i)->GetparameterMin().size();
			int i_npaths = optpro_->GetPhase(i)->GetpathMin().size();
			int i_nevents = optpro_->GetPhase(i)->GeteventMin().size();
			int i_totnodes = optpro_->GetPhase(i)->GetTotalNodes();
			vec hessian_phase_V, hessian_phase_I, hessian_phase_J;
			GetPhaseHessianSparsity(i,  dependencies[i], hessian_phase_I, hessian_phase_J);
			hessain_I.subvec(hessian_shift, hessian_shift + hessian_phase_I.n_elem - 1) = hessian_phase_I + rowshift;
			hessain_J.subvec(hessian_shift, hessian_shift + hessian_phase_J.n_elem - 1) = hessian_phase_J + colshift;
			hessian_shift += hessian_phase_I.n_elem;

			size_t	numcons = i_nstate*(i_totnodes + 1) + i_ncontrol*i_totnodes + i_nparameters + 2;
			size_t numvars = i_nstate*(i_totnodes + 1) + i_ncontrol*i_totnodes + i_nparameters + 2;
			rowshift += numcons;
			colshift += numvars;
		}
		//!<derive of linkage and it's Sparsity Pattern
		
		for (size_t ipair = 0; ipair < optpro_->GetLinkageNum(); ipair++)
		{

			vec hessian_ipair_V, hessian_ipair_I, hessian_ipair_J;
			GetLinkHessianSparsity(ipair,  hessian_ipair_I, hessian_ipair_J);
			hessain_I.subvec(hessian_shift, hessian_shift + hessian_ipair_I.n_elem - 1) = hessian_ipair_I;
			hessain_J.subvec(hessian_shift, hessian_shift + hessian_ipair_J.n_elem - 1) = hessian_ipair_J;
			hessian_shift += hessian_ipair_I.n_elem;
		}
	}

}//namespace Lpopc