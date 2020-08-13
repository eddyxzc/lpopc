// Copyright (C) 2014-2015 Xue Zhichen, Wang Yujie,Wang Na 
// All Rights Reserved.
// This file is a part of LPOPC , published under the Eclipse Public License.
// Author:Xue Zhichen 8/4 2015   21:54
// Email:eddy_lpopc@163.com
#include "LpFiniteDifferenceDerive.hpp"
namespace Lpopc
{


	void LpFDderive::DerivMayer(SolCost&mySolcost, rowvec& deriv_mayer)
	{
		LP_DBG_START_FUN("LpFDderive::DerivMayer")
		size_t iphase = mySolcost.phase_num_ - 1;
		size_t nstate = mySolcost.state_.n_cols;
		size_t ncontrols = mySolcost.control_.n_cols;
		size_t nparameters = mySolcost.parameter_.n_elem;
		size_t nnodes = mySolcost.state_.n_rows;


		double pert0 = tol_*(1 + abs(mySolcost.initial_time_));
		double pertf = tol_*(1 + abs(mySolcost.terminal_time_));
		vec pertx0 = tol_*(1 + abs(mySolcost.initial_state_));
		vec pertxf = tol_*(1 + abs(mySolcost.terminal_state_));
		mat perParameter;
		if (nparameters > 0) perParameter = tol_*(1 + abs(mySolcost.parameter_));
		double t0Pert = mySolcost.initial_time_ + pert0;
		vec x0Pert = mySolcost.initial_state_ + pertx0;
		double tfPert = mySolcost.terminal_time_ + pertf;
		vec xfPert = mySolcost.terminal_state_ + pertxf;

		double DMayer_t0 = 0, DMayer_tf = 0;
		rowvec DMayer_x0 = zeros(1, nstate);
		rowvec DMayer_xf = zeros(1, nstate);
		rowvec DMyayer_parameter;
		if (nparameters > 0) DMyayer_parameter = zeros<rowvec>(nparameters);
		double mayerout = 0;
		fun_->MayerCost(mySolcost, mayerout);
		double perMayerout = 0;
		//Get Derivatives of Mayer Cost with Respect to Initial Time
		double t0 = mySolcost.initial_time_;
		mySolcost.initial_time_ = t0Pert;
		fun_->MayerCost(mySolcost, perMayerout);
		DMayer_t0 = (perMayerout - mayerout) / pert0;
		mySolcost.initial_time_ = t0;
		//Get Derivatives of Mayer Cost with Respect to Terminal Time
		double tf = mySolcost.terminal_time_;
		mySolcost.terminal_time_ = tfPert;
		fun_->MayerCost(mySolcost, perMayerout);
		DMayer_tf = (perMayerout - mayerout) / pertf;
		mySolcost.terminal_time_ = tf;
		//Get Derivatives of Mayer Cost with Respect to Initial and Terminal State
		vec x0 = mySolcost.initial_state_;
		vec xf = mySolcost.terminal_state_;
		for (size_t istate = 0; istate < nstate; istate++)
		{   //dmayer/dx0
			mySolcost.initial_state_(istate) = x0Pert(istate);
			fun_->MayerCost(mySolcost, perMayerout);
			DMayer_x0(istate) = (perMayerout - mayerout) / pertx0(istate);
			mySolcost.initial_state_(istate) = x0(istate);
			//dmayer/dxf
			mySolcost.terminal_state_(istate) = xfPert(istate);
			fun_->MayerCost(mySolcost, perMayerout);
			DMayer_xf(istate) = (perMayerout - mayerout) / pertxf(istate);
			mySolcost.terminal_state_(istate) = xf(istate);

		}
		//Get Derivatives of Mayer Cost with Respect to Parameters
		vec parameters = mySolcost.parameter_;
		vec parameterPert;
		if (nparameters > 0)
		{
			parameterPert = mySolcost.parameter_ + perParameter;
		}
		for (size_t ipara = 0; ipara < nparameters; ipara++)
		{
			mySolcost.parameter_(ipara) = parameterPert(ipara);
			fun_->MayerCost(mySolcost, perMayerout);
			DMyayer_parameter(ipara) = (perMayerout - mayerout) / perParameter(ipara);
			mySolcost.parameter_(ipara) = parameters(ipara);
		}
		//insert derive
		//deriv_mayer = [DMayer_x0,DMayer_t0,DMayer_xf,DMayer_tf,DMyayer_parameter]
		deriv_mayer = zeros(1, nstate + 1 + nstate + 1 + nparameters);
		size_t colshift = 0;
		deriv_mayer.subvec(0, nstate - 1) = DMayer_x0;
		colshift += nstate;
		deriv_mayer(colshift) = DMayer_t0;
		colshift++;
		deriv_mayer.subvec(colshift, colshift + nstate - 1) = DMayer_xf;
		colshift += nstate;
		deriv_mayer(colshift) = DMayer_tf;
		colshift++;
		if (nparameters > 0)
		{
			deriv_mayer.subvec(colshift, colshift + nparameters - 1) = DMyayer_parameter;
		}
	}

	void LpFDderive::DerivLagrange(SolCost& mySolcost, mat& deriv_langrange)
	{
		LP_DBG_START_FUN("LpFDderive::DerivLagrange")
		size_t iphase = mySolcost.phase_num_ - 1;
		size_t nstate = mySolcost.state_.n_cols;
		size_t ncontrols = mySolcost.control_.n_cols;
		size_t nparameters = mySolcost.parameter_.n_elem;
		size_t nnodes = mySolcost.state_.n_rows;

		vec pertTime = tol_*(1 + abs(mySolcost.time_));
		mat pertState = tol_*(1 + abs(mySolcost.state_));
		mat pertControl, perParameter;
		if (ncontrols > 0) pertControl = tol_*(1 + abs(mySolcost.control_));
		if (nparameters > 0) perParameter = tol_*(1 + abs(mySolcost.parameter_));

		vec DLagrange_time = zeros(nnodes, 1);
		mat DLagrange_state = zeros(nnodes, nstate);
		mat DLagrange_control, DLagrange_parameter;
		if (ncontrols > 0)DLagrange_control = zeros(nnodes, ncontrols);
		if (nparameters > 0)DLagrange_parameter = zeros(nnodes, nparameters);
		vec tRadauPert = mySolcost.time_ + pertTime;
		mat stateRadauPert = mySolcost.state_ + pertState;
		mat controlRadauPert, parameterPert;

		if (ncontrols > 0)
		{
			controlRadauPert = mySolcost.control_ + pertControl;
			DLagrange_control = zeros(nnodes, ncontrols);
		}
		if (nparameters > 0)
		{
			parameterPert = mySolcost.parameter_ + perParameter;
			DLagrange_parameter = zeros(nnodes, nparameters);
		}
		vec lagrangeOut, perlagrangeOut;
		fun_->LagrangeCost(mySolcost, lagrangeOut);
		//Get Derivatives of Lagrange Cost with Respect to  Time
		vec t_radau = mySolcost.time_;
		mySolcost.time_ = tRadauPert;
		fun_->LagrangeCost(mySolcost, perlagrangeOut);
		DLagrange_time = (perlagrangeOut - lagrangeOut) / pertTime;
		mySolcost.time_ = t_radau;

		//Get Derivatives of Lagrange Cost with Respect to State
		mat state_radau = mySolcost.state_;
		for (size_t istate = 0; istate < nstate; istate++)
		{
			mySolcost.state_.col(istate) = stateRadauPert.col(istate);
			fun_->LagrangeCost(mySolcost, perlagrangeOut);
			DLagrange_state.col(istate) = (perlagrangeOut - lagrangeOut) / pertState.col(istate);
			mySolcost.state_.col(istate) = state_radau.col(istate);
		}

		//Get Derivatives of Lagrange Cost with Respect to control
		mat control_radau = mySolcost.control_;
		for (size_t icontrol = 0; icontrol < ncontrols; icontrol++)
		{
			mySolcost.control_.col(icontrol) = controlRadauPert.col(icontrol);
			fun_->LagrangeCost(mySolcost, perlagrangeOut);
			DLagrange_control.col(icontrol) = (perlagrangeOut - lagrangeOut) / pertControl.col(icontrol);
			mySolcost.control_.col(icontrol) = control_radau.col(icontrol);
		}
		//Get Derivatives of Lagrange Cost with Respect to Parameters
		vec parameter = mySolcost.parameter_;
		for (size_t ipara = 0; ipara < nparameters; ipara++)
		{
			mySolcost.parameter_(ipara) = parameterPert(ipara);
			fun_->LagrangeCost(mySolcost, perlagrangeOut);
			DLagrange_parameter.col(ipara) = (perlagrangeOut - lagrangeOut) / (perParameter(ipara)*ones<vec>(lagrangeOut.n_elem));
			mySolcost.parameter_(ipara) = parameter(ipara);
		}
		//insert derive
		//DLagrange=[dL/dx,dL/du,dL/dt,dl/dp]
		deriv_langrange = zeros(nnodes, nstate + ncontrols + 1 + nparameters);
		size_t colshift = 0;
		for (size_t istate = 0; istate < nstate; istate++)
		{
			deriv_langrange.col(colshift) = DLagrange_state.col(istate);
			colshift++;
		}

		for (size_t icontrol = 0; icontrol < ncontrols; icontrol++)
		{
			deriv_langrange.col(colshift) = DLagrange_control.col(icontrol);
			colshift++;
		}
		deriv_langrange.col(nstate + ncontrols) = DLagrange_time;
		colshift = nstate + ncontrols + 1;
		for (size_t ipara = 0; ipara < nparameters; ipara++)
		{
			deriv_langrange.col(colshift + ipara) = DLagrange_parameter.col(ipara);
		}
	}

	void LpFDderive::DerivDae(SolDae& mySolDae, mat& deriv_state, mat& deriv_path)
	{
		LP_DBG_START_FUN("LpFDderive::DerivDae")
		size_t iphase = mySolDae.phase_num_ - 1;
		size_t nstate = mySolDae.state_.n_cols;
		size_t ncontrols = mySolDae.contol_.n_cols;
		size_t nparameters = mySolDae.parameter_.n_elem;

		size_t nnodes = mySolDae.state_.n_rows;
		//////////////////////////////////////////////////////////////////////////
		//Get user dae function
		mat daeout, pathout;
		fun_->DaeFunction(mySolDae, daeout, pathout);
		size_t npaths = pathout.n_cols;
		vec pertTime = tol_*(1 + abs(mySolDae.time_));
		mat pertState = tol_*(1 + abs(mySolDae.state_));
		mat pertControl, perParameter;
		if (ncontrols > 0) pertControl = tol_*(1 + abs(mySolDae.contol_));
		if (nparameters > 0) perParameter = tol_*(1 + abs(mySolDae.parameter_));
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

		//Compute the Derivative of DAEs with Respect to Time
		mat dDae_time = zeros(nnodes, nstate + npaths);
		vec t_radau = mySolDae.time_;
		mySolDae.time_ = tRadauPert;
		mat denominator = repmat(pertTime, 1, nstate + npaths);
		mat perstateout, perpathout;
		fun_->DaeFunction(mySolDae, perstateout, perpathout);
		mySolDae.time_ = t_radau;

		if (pathout.n_elem > 0)
		{
			dDae_time = (join_horiz(perstateout, perpathout) - join_horiz(daeout, pathout)) / (denominator);
		}
		else
		{
			dDae_time = (perstateout - daeout) / denominator;
		}
		//Compute the Derivative of DAEs with Respect to State
		mat state_radau = mySolDae.state_;
		std::vector<mat> dDae_state(nstate);
		for (size_t istate = 0; istate < nstate; istate++)
		{
			mySolDae.state_.col(istate) = stateRadauPert.col(istate);
			denominator = repmat(pertState.col(istate), 1, nstate + npaths);
			fun_->DaeFunction(mySolDae, perstateout, perpathout);
			if (pathout.n_elem > 0)
			{
				dDae_state[istate] = (join_horiz(perstateout, perpathout) - join_horiz(daeout, pathout)) / (denominator);
			}
			else
			{
				dDae_state[istate] = (perstateout - daeout) / denominator;
			}
			mySolDae.state_.col(istate) = state_radau.col(istate);
		}
		//Compute the Derivative of DAEs with Respect to Control
		mat control_radau = mySolDae.contol_;
		std::vector<mat> dDae_control(ncontrols);
		for (size_t icontrol = 0; icontrol < ncontrols; icontrol++)
		{
			mySolDae.contol_.col(icontrol) = controlRadauPert.col(icontrol);
			denominator = repmat(pertControl.col(icontrol), 1, nstate + npaths);
			fun_->DaeFunction(mySolDae, perstateout, perpathout);
			if (pathout.n_elem > 0)
			{
				dDae_control[icontrol] = (join_horiz(perstateout, perpathout) - join_horiz(daeout, pathout)) / (denominator);
			}
			else
			{

				dDae_control[icontrol] = (perstateout - daeout) / denominator;
			}
			mySolDae.contol_.col(icontrol) = control_radau.col(icontrol);
		}
		//Compute the Derivative of DAEs with Respect to Parameter.
		mat parameters = mySolDae.parameter_;
		std::vector<mat> dDae_parameter(nparameters);
		for (size_t iparameter = 0; iparameter < nparameters; iparameter++)
		{
			mySolDae.parameter_.row(iparameter) = parameterPert.row(iparameter);
			denominator = perParameter(iparameter)*ones(nnodes, nstate + npaths);
			fun_->DaeFunction(mySolDae, perstateout, perpathout);
			if (pathout.n_elem > 0)
			{
				dDae_parameter[iparameter] = (join_horiz(perstateout, perpathout) - join_horiz(daeout, pathout)) / (denominator);
			}
			else
			{

				dDae_parameter[iparameter] = (perstateout - daeout) / denominator;
			}
			mySolDae.parameter_(iparameter) = parameters(iparameter);
		}
		//insert derive
		mat derive_dae(nnodes*(nstate + npaths), nstate + ncontrols + 1 + nparameters);
		size_t colshift = 0;
		for (size_t istate = 0; istate < nstate; istate++)
		{
			derive_dae.col(colshift) = reshape(dDae_state[istate], nnodes*(nstate + npaths), 1);
			colshift++;
		}

		for (size_t icontrol = 0; icontrol < ncontrols; icontrol++)
		{
			derive_dae.col(colshift) = reshape(dDae_control[icontrol], nnodes*(nstate + npaths), 1);
			colshift++;
		}
		derive_dae.col(nstate + ncontrols) = reshape(dDae_time, derive_dae.n_rows, 1);
		colshift = nstate + ncontrols + 1;
		for (size_t ipara = 0; ipara < nparameters; ipara++)
		{
			derive_dae.col(colshift + ipara) = reshape(dDae_parameter[ipara], derive_dae.n_rows, 1);
		}
		//deriv_state = zeros(nnodes*nstate, nstate + ncontrols + 1 + nparameters);
		deriv_state = derive_dae.rows(0, nnodes*nstate - 1);
		if (npaths > 0)
		{
			deriv_path = derive_dae.rows(nnodes*nstate, nnodes*(nstate + npaths) - 1);
		}
	}

	void LpFDderive::DerivEvent(SolEvent& mySolEvent, mat& deriv_event)
	{
		LP_DBG_START_FUN("LpFDderive::DerivEvent")
		size_t iphase = mySolEvent.phase_num_ - 1;
		size_t nstate = mySolEvent.initial_state_.n_elem;
		size_t nparameters = mySolEvent.parameter_.n_elem;

		double t0Scales = 1.0;
		double tfScales = 1.0;
		double pert0 = tol_*(1 + abs(mySolEvent.initial_time_));
		double pertf = tol_*(1 + abs(mySolEvent.terminal_time_));
		vec pertx0 = tol_*(abs(mySolEvent.initial_state_) + 1);
		vec pertxf = tol_*(abs(mySolEvent.terminal_state_) + 1);
		mat perParameter;
		if (nparameters > 0) perParameter = tol_*(abs(mySolEvent.parameter_) + 1);
		vec eventout;
		fun_->EventFunction(mySolEvent, eventout);
		size_t nevents = eventout.n_elem;
		double t0Pert = mySolEvent.initial_time_ + pert0;
		vec x0Pert = mySolEvent.initial_state_ + pertx0;
		double tfPert = mySolEvent.terminal_time_ + pertf;
		vec xfPert = mySolEvent.terminal_state_ + pertxf;

		vec perEventOut;
		//Get Derivatives of Event with Respect to Initial Time
		double t0 = mySolEvent.initial_time_;
		mySolEvent.initial_time_ = t0Pert;
		fun_->EventFunction(mySolEvent, perEventOut);
		vec DEvent_t0 = (perEventOut - eventout) / pert0;
		mySolEvent.initial_time_ = t0;

		//Get Derivatives of Event with Respect to Terminal Time
		double tf = mySolEvent.terminal_time_;
		mySolEvent.terminal_time_ = tfPert;
		fun_->EventFunction(mySolEvent, perEventOut);
		vec DEvent_tf = (perEventOut - eventout) / pertf;
		mySolEvent.terminal_time_ = tf;

		//Get Derivatives of Event with Respect to Initial and Terminal State
		mat DEvent_x0(nevents, nstate), DEvent_xf(nevents, nstate);
		vec x0 = mySolEvent.initial_state_;
		vec xf = mySolEvent.terminal_state_;

		for (size_t istate = 0; istate < nstate; istate++)
		{
			mySolEvent.initial_state_(istate) = x0Pert(istate);
			fun_->EventFunction(mySolEvent, perEventOut);
			DEvent_x0.col(istate) = (perEventOut - eventout) / (pertx0(istate)*ones<vec>(eventout.n_elem));
			mySolEvent.initial_state_(istate) = x0(istate);

			mySolEvent.terminal_state_(istate) = xfPert(istate);
			fun_->EventFunction(mySolEvent, perEventOut);
			DEvent_xf.col(istate) = (perEventOut - eventout) / (pertxf(istate)*ones<vec>(eventout.n_elem));
			mySolEvent.terminal_state_(istate) = xf(istate);
		}
		//Get Derivatives of Event with Respect to Parameters
		mat DEvent_para;
		if (nparameters>0)DEvent_para = zeros(nevents, nparameters);
		vec paramemter = mySolEvent.parameter_;
		for (size_t ipara = 0; ipara < nparameters; ipara++)
		{
			mySolEvent.parameter_(ipara) += perParameter(ipara);
			fun_->EventFunction(mySolEvent, perEventOut);
			DEvent_para.col(ipara) = (perEventOut - eventout) / (perParameter(ipara)*ones<vec>(nevents));
			mySolEvent.parameter_(ipara) = paramemter(ipara);
		}

		//insert Derive 
		//deriv_event = [DEvent_x0,DEvent_t0,DEvent_xf,DEvent_tf,DEvent_para]
		deriv_event = zeros(nevents, nstate + 1 + nstate + 1 + nparameters);
		size_t colshift = 0;
		deriv_event.cols(0, nstate - 1) = DEvent_x0;
		colshift += nstate;
		deriv_event.col(colshift) = DEvent_t0;
		colshift++;
		deriv_event.cols(colshift, colshift + nstate - 1) = DEvent_xf;
		colshift += nstate;
		deriv_event.col(colshift) = DEvent_tf;
		colshift++;
		if (nparameters > 0)
		{
			deriv_event.cols(colshift, colshift + nparameters - 1) = DEvent_para;
		}
	}

	void LpFDderive::DerivLink(SolLink& mySolLink, mat& derive_link)
	{
		LP_DBG_START_FUN("LpFDderive::DerivLink")
		vec Linkout;
		fun_->LinkFunction(mySolLink, Linkout);
		size_t left_iphase = mySolLink.left_phase_num_-1;
		size_t right_iphase = mySolLink.right_phase_num_-1;

		size_t nstateLeft =mySolLink.left_state_.n_elem;
		size_t nparametersLeft = mySolLink.left_parameter_.n_elem;
		size_t nnodesLeft = mySolLink.left_state_.n_elem;

		size_t nstateRight = mySolLink.right_state_.n_elem;
		size_t nparametersRight = mySolLink.right_parameter_.n_elem;
		size_t nnodesRight = mySolLink.right_state_.n_elem;

		mat xf_left = mySolLink.left_state_;
		mat x0_right = mySolLink.right_state_;

		mat p_left = mySolLink.left_parameter_;
		mat p_right = mySolLink.right_parameter_;

		size_t ipair = mySolLink.ipair - 1;
		size_t nlinks = Linkout.n_rows;
		mat DLink_xf_left = zeros(nlinks, nstateLeft);
		mat DLink_x0_right = zeros(nlinks, nnodesRight);
		mat DLink_p_left, DLink_p_right;
		if (nparametersLeft>0)
		{
			DLink_p_left = zeros(nlinks, nparametersLeft);
		}
		if (nparametersRight>0)
		{
			DLink_p_right = zeros(nlinks, nparametersRight);
		}
		mat pertxfLeft = tol_*(1 + abs(xf_left));
		mat pertpLeft = tol_*(1 + abs(p_left));
		mat pertx0Right = tol_*(1 + abs(x0_right));
		mat pertpRight = tol_*(1 + abs(p_right));

		mat xfLeftPert = xf_left + pertxfLeft;
		mat pLeftPert = p_left + pertpLeft;
		mat x0RightPert = x0_right + pertx0Right;
		mat pRightPert = p_right + pertpRight;
		vec perLinkout;
		//Get Derivatives of Linkage with Respect to xf_left
		for (size_t istate = 0; istate < nstateLeft; istate++)
		{
			mySolLink.left_state_(istate)=(xfLeftPert(istate));
			fun_->LinkFunction(mySolLink, perLinkout);
			DLink_xf_left.col(istate) = (perLinkout - Linkout) / (ones<vec>(Linkout.n_rows)*pertxfLeft(istate));
			mySolLink.left_state_(istate)=(xf_left(istate));
		}
		//Get Derivatives of Linkage with Respect to p_left
		for (size_t ipara = 0; ipara < nparametersLeft; ipara++)
		{
			mySolLink.left_parameter_(ipara)=(pLeftPert(ipara));
			fun_->LinkFunction(mySolLink, perLinkout);
			DLink_p_left.col(ipara) = (perLinkout - Linkout) / (ones<vec>(Linkout.n_rows)*pertpLeft(ipara));
			mySolLink.left_parameter_(ipara)=(p_left(ipara));
		}
		//Get Derivatives of Linkage with Respect to x0_right
		for (size_t istate = 0; istate < nstateRight; istate++)
		{
			mySolLink.right_state_(istate)=(x0RightPert(istate));
			fun_->LinkFunction(mySolLink, perLinkout);
			DLink_x0_right.col(istate) = (perLinkout - Linkout) / (ones<vec>(Linkout.n_rows)*pertx0Right(istate));
			mySolLink.right_state_(istate)=(x0_right(istate));
		}
		//Get Derivatives of Linkage with Respect to p_right
		for (size_t ipara = 0; ipara < nparametersRight; ipara++)
		{
			mySolLink.right_parameter_(ipara)=(pRightPert(ipara));
			fun_->LinkFunction(mySolLink, perLinkout);
			DLink_p_right.col(ipara) = (perLinkout - Linkout) / (ones<vec>(Linkout.n_rows)*pertpLeft(ipara));
			mySolLink.right_parameter_(ipara)=(p_right(ipara));
		}

		//insert derive
		//derive_link=[DLink_xf_left,DLink_p_left,DLink_x0_right,DLink_p_left]
		derive_link = zeros(nlinks, nstateLeft + nparametersLeft + nstateRight + nparametersRight);
		derive_link.cols(0, nstateLeft - 1) = DLink_xf_left;
		if (nparametersLeft>0)
		{
			derive_link.cols(nstateLeft, nstateLeft + nparametersLeft - 1) = DLink_p_left;
		}
		derive_link.cols(nstateLeft + nparametersLeft, nstateLeft + nparametersLeft + nstateRight - 1) = DLink_x0_right;
		if (nparametersRight>0)
		{
			derive_link.cols(nstateLeft + nparametersLeft + nstateRight, derive_link.n_cols-1) = DLink_p_left;
		}
	}

}//namespace Lpopc