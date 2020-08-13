// Copyright (C) 2014-2015 Xue Zhichen, Wang Yujie,Wang Na 
// All Rights Reserved.
// This file is a part of LPOPC , published under the Eclipse Public License.
// Author:Xue Zhichen 7/12 2015   20:16
// Email:eddy_lpopc@163.com
#include "LpNLPWrapper.hpp"
#include "LpFunctionWrapper.h"
#include "LpSizeChecker.h"
#include "LpBoundsChecker.hpp"
#include "LpGuessChecker.h"
#include "IpIpoptApplication.hpp"
#include "LpopcIpopt.h"
#include "LpSizeChecker.h"
#include "LpGuessChecker.h"
#include <strstream>
#include <vector>

using std::vector;
namespace Lpopc
{

	int NLPWrapper::nnz(mat& op_mat)
	{
		int nonZeroNumber = 0;
		for (int i = 0; i < op_mat.n_elem; i++)
		{
			if (op_mat[i] != 0)
			{
				nonZeroNumber++;
			}
		}
		return nonZeroNumber;
	}
	void NLPWrapper::GetAllCons(vec& x, vec& Cons)
	{
		LP_DBG_START_FUN("NLPWrapper::GetAllCons()")
		vec y = x;
		if (calculateData_->autoscale)
		{
			y = (y - calculateData_->varshift) / calculateData_->varscale;
		}
		
		vec nonLinearCons, linearCons;
		GetConsFun(y, nonLinearCons);
		linearCons = calculateData_->AlinearMatrix*y;
		Cons = zeros(nonLinearCons.n_elem + linearCons.n_elem, 1);
		Cons.subvec(0, nonLinearCons.n_elem - 1) = nonLinearCons;
		Cons.subvec(nonLinearCons.n_elem, Cons.n_elem - 1) = linearCons;
		if (calculateData_->autoscale)
		{
			Cons %= calculateData_->funscale;
		}
	}

	void NLPWrapper::GetConsFun(const vec& x, vec& Cons)
	{
		LP_DBG_START_FUN("NLPWrapper::GetConsFun()")
		vector<vec> AllCons(calculateData_->numphases_);
		vector<SolCost> solTotal(calculateData_->numphases_);
		int total_constraints = 0;
		for (int i = 0; i < calculateData_->numphases_; i++)
		{
			size_t nstate = calculateData_->SIZES_[i][0];
			size_t ncontrols = calculateData_->SIZES_[i][1];
			size_t nparameters = calculateData_->SIZES_[i][2];
			size_t npaths = calculateData_->SIZES_[i][3];
			size_t nevents = calculateData_->SIZES_[i][4];

			int state_index_start = calculateData_->phase_indices[i]->state[0] - 1;
			int state_index_end = calculateData_->phase_indices[i]->state[calculateData_->phase_indices[i]->state.size() - 1] - 1;
			vec state_vector = x.subvec(state_index_start, state_index_end);

			int control_index_start = calculateData_->phase_indices[i]->control[0] - 1;
			int control_index_end = calculateData_->phase_indices[i]->control[calculateData_->phase_indices[i]->control.size() - 1] - 1;
			vec control_vector = x.subvec(control_index_start, control_index_end);

			double t0 = x(calculateData_->phase_indices[i]->time[0] - 1);
			double tf = x(calculateData_->phase_indices[i]->time[1] - 1);
			double tspan = tf - t0;
			vec t_radau = (calculateData_->PS[i]->Points + 1)*(tspan / 2.0) + t0;

			mat state_matrix = reshape(state_vector,
				optpro_->GetPhase(i)->GetTotalNodes() + 1, nstate);
			mat state_radau = state_matrix(span(0, state_matrix.n_rows - 2), span::all);

			vec x0 = trans(state_matrix.row(0));
			vec xf = trans(state_matrix.row(state_matrix.n_rows - 1));
			mat control_radau = reshape(control_vector,
				optpro_->GetPhase(i)->GetTotalNodes(), ncontrols);
			vec parameter_vec;
			if (nparameters > 0)
			{
				int para_index_start = calculateData_->phase_indices[i]->parameter[0] - 1;
				int para_index_end = calculateData_->phase_indices[i]->parameter[calculateData_->phase_indices[i]->parameter.size() - 1] - 1;
				parameter_vec = x.subvec(para_index_start, para_index_end);
			}

			//Get DAE
			SolDae mySolDae;
			mySolDae.time_ = t_radau;
			//t_radau.print("t_radau");
			mySolDae.state_ = state_radau;
			//state_radau.print("state_radau");
			mySolDae.contol_ = control_radau;
			//control_radau.print("control_radau");
			mySolDae.parameter_ = parameter_vec;
			mySolDae.phase_num_ = i + 1;
			mat stateout(optpro_->GetPhase(i)->GetTotalNodes(), nstate);
			mat pathout(optpro_->GetPhase(i)->GetTotalNodes(), ncontrols);
			optimalFunction_->DaeFunction(mySolDae, stateout, pathout);
			mat odeleft = calculateData_->PS[i]->D*state_matrix;
			//odeleft.print("odeleft");
			mat oderight = stateout*(tspan / 2.0);
			//stateout.print("stateout");
			//pathout.print("pathout");
			////////////////////////////////////////////////////////
			//odeleft.print("odeleft");
			SolDae temSolDae = mySolDae;
			//temSolDae.state_.()
			////////////////////////////////////////////////////
			//oderight.print("oderight");
			mat defects = odeleft - oderight;
			//defects.print("defects");
			vec events;
			if (nevents > 0)
			{
				SolEvent mySolEvent;
				mySolEvent.initial_time_ = t0;
				mySolEvent.initial_state_ = x0;
				mySolEvent.terminal_time_ = tf;
				mySolEvent.terminal_state_ = xf;
				mySolEvent.parameter_ = parameter_vec;
				mySolEvent.phase_num_ = i + 1;
				events = zeros(nevents, 1);
				optimalFunction_->EventFunction(mySolEvent, events);
			}

			if (npaths > 0 && nevents > 0)
			{
				vec consi(defects.n_elem + pathout.n_elem + events.n_elem, 1);
				consi.subvec(0, defects.n_elem - 1) = reshape(defects, defects.n_elem, 1);
				consi.subvec(defects.n_elem, defects.n_elem + pathout.n_elem - 1) = reshape(pathout, pathout.n_elem, 1);
				consi.subvec(defects.n_elem + pathout.n_elem,
					defects.n_elem + pathout.n_elem + events.n_elem - 1) = reshape(events, events.n_elem, 1);
				AllCons[i] = consi;
			}
			else if (npaths > 0)
			{
				mat consi(defects.n_rows, defects.n_cols + pathout.n_cols);
				consi.cols(0, defects.n_cols - 1) = defects;
				consi.cols(defects.n_cols, consi.n_cols - 1) = pathout;
				AllCons[i] = reshape(consi, consi.n_elem, 1);
			}
			else if (nevents > 0)
			{
				vec consi(defects.n_elem + events.n_elem, 1);
				consi.subvec(0, defects.n_elem - 1) = reshape(defects, defects.n_elem, 1);
				consi.subvec(defects.n_elem, defects.n_elem + events.n_elem - 1) = reshape(events, events.n_elem, 1);
				AllCons[i] = consi;
			}
			else
			{
				AllCons[i] = reshape(defects, defects.n_elem, 1);
			}
			//AllCons[i].print("AllCons[i]");
			total_constraints += AllCons[i].n_elem;
			SolCost mySolCost;
			mySolCost.initial_time_ = t0;
			mySolCost.initial_state_ = x0;
			mySolCost.terminal_time_ = tf;
			mySolCost.terminal_state_ = xf;
			mySolCost.time_ = t_radau;
			mySolCost.state_ = state_radau;
			mySolCost.control_ = control_radau;
			mySolCost.parameter_ = parameter_vec;
			mySolCost.phase_num_ = i + 1;
			solTotal[i] = mySolCost;
		}//end for

		vector<vec> AllLink(optpro_->GetLinkageNum());
		int total_links = 0;
		if (calculateData_->numlinks_ > 0)
		{
			for (int ipair = 0; ipair < optpro_->GetLinkageNum(); ipair++)
			{
				int left_index = optpro_->GetLinkage(ipair)->LeftPhase();
				int right_index = optpro_->GetLinkage(ipair)->RightPhase();
				vec xf_left = solTotal[left_index].terminal_state_;
				mat p_left = solTotal[left_index].parameter_;

				vec x0_right = solTotal[right_index].initial_state_;
				mat p_right = solTotal[right_index].parameter_;

				SolLink mySolLink;
				mySolLink.left_state_ = xf_left;
				//xf_left.print("xf_left");
				//x0_right.print("x0_right");
				mySolLink.left_parameter_ = p_left;
				mySolLink.left_phase_num_ = left_index;

				mySolLink.right_state_ = x0_right;
				mySolLink.right_parameter_ = p_left;
				mySolLink.right_phase_num_ = right_index;
				vec linkout(x0_right.n_elem, 1);
				optimalFunction_->LinkFunction(mySolLink, linkout);
				AllLink[ipair] = linkout;
				//linkout.print("linkout");
				total_links += AllLink[ipair].n_elem;
			}

		}//end for ipair

		vec ConsTemp(total_links + total_constraints, 1);
		int rowshift = 0;
		for (int i = 0; i < AllCons.size(); i++)
		{
			ConsTemp.subvec(rowshift, rowshift + (AllCons[i].n_elem) - 1) =
				AllCons[i];
			rowshift += AllCons[i].n_elem;
		}

		for (int ipair = 0; ipair < AllLink.size(); ipair++)
		{
			ConsTemp.subvec(rowshift, rowshift + AllLink[ipair].n_elem - 1) =
				AllLink[ipair];
			rowshift += AllLink[ipair].n_elem;
		}
		Cons = ConsTemp;
	}
	void NLPWrapper::GetConsJacbi(const vec& scaled_x, vec& Sjac_V)
	{
		   LP_DBG_START_FUN("NLPWrapper::GetConsJacbi()");
			vec x = scaled_x;
			if (calculateData_->autoscale)
			{
				x = (scaled_x - calculateData_->varshift) / calculateData_->varscale;
			}
			vec  NL_V, L_I, L_J, L_V, C_V;

			GetWholeJacbi(x, NL_V, C_V);

		dsmatrix::Find(calculateData_->AlinearMatrix, L_I, L_J, L_V);

			int nonZerosNum = NL_V.n_elem + L_V.n_elem + C_V.n_elem;
			int linearConsRowStart = calculateData_->conbounds_min.size();//Count From 0
			Sjac_V = zeros(nonZerosNum, 1);
			int rowshift = 0;
			Sjac_V.rows(0, NL_V.n_elem - 1) = NL_V;
			rowshift += NL_V.n_elem;
			Sjac_V.rows(NL_V.n_elem, rowshift + L_V.n_elem - 1) = L_V;
			rowshift += L_V.n_elem;
			Sjac_V.rows(rowshift, rowshift + C_V.n_elem - 1) = C_V;

			if (calculateData_->autoscale)
			{
				Sjac_V %= conscale_;
			}

	}
	void NLPWrapper::GetWholeJacbi(const vec&x, vec& Sjac_V, vec& Sconstant_V)
	{
		LP_DBG_START_FUN("NLPWrapper::GetWholeJacbi()")
		vector<mat> dependencies(optpro_->GetPhaseNum());

		for (int itor = 0; itor < dependencies.size(); itor++)
		{
			int i_nstate = optpro_->GetPhase(itor)->GetstateMin().size();
			int i_ncontrol = optpro_->GetPhase(itor)->GetcontrolMin().size();
			int i_npaths = optpro_->GetPhase(itor)->GetpathMin().size();
			dependencies[itor] = zeros(i_nstate + i_npaths, i_nstate + i_ncontrol);
		}

		//find total non-zero elements in Jacobian and constant derivatives
		int  nonZerosSjac = 0;
		int  nonZerosSconstant = 0;
		int  nonZerosSnonconstant = 0;

		for (int i = 0; i < optpro_->GetPhaseNum(); i++)
		{
			int i_nstate = optpro_->GetPhase(i)->GetstateMin().size();
			int i_ncontrol = optpro_->GetPhase(i)->GetcontrolMin().size();
			int i_nparameters = optpro_->GetPhase(i)->GetparameterMin().size();
			int i_npaths = optpro_->GetPhase(i)->GetpathMin().size();
			int i_nevents = optpro_->GetPhase(i)->GeteventMin().size();
			int i_totnodes = optpro_->GetPhase(i)->GetTotalNodes();
			// find number of non-zero elements in dependencies (always including diagonal)


			//////////////////////////////////////////////////////////////////////////
			//only when use Analytical Derives 
			dependencies[i].fill(1.0);
			////////////////////////////////////////////////////////////////////////
			mat temDependencies = dependencies[i];
			for (int j = 0; j < i_nstate; j++)
			{
				temDependencies(j, j) = 1.0;
			}

			int nDependancies = nnz(temDependencies);

			nonZerosSjac = nonZerosSjac + nDependancies*i_totnodes + 2 * (i_nstate + i_npaths)*i_totnodes +
				(i_nstate + i_npaths)*i_nparameters*i_totnodes + i_nevents*(2 * i_nstate + i_nparameters + 2);
			vec DiffMatOffDiag_I, DiffMatOffDiag_J, DiffMatOffDiag_V;
			dsmatrix::Find(calculateData_->PS[i]->Doffdiag, DiffMatOffDiag_I, DiffMatOffDiag_J, DiffMatOffDiag_V);
			int nonZerosDiffMat = DiffMatOffDiag_I.n_elem;
			nonZerosSconstant = nonZerosSconstant + nonZerosDiffMat*i_nstate;
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
			nonZerosSjac += numlinks*(nstates_left + nparameter_left + nstates_right + nparameter_right);
		}//for ipair
		//allocate Memory for Jacobian
		Sjac_V = zeros(nonZerosSjac, 1);
		Sconstant_V = zeros(nonZerosSconstant, 1);
		int Sjac_rowsShift = 0, Sconstant_rowsShift = 0;
		int rowshift = 0;
		int colshift = 0;
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
			vec  Sjac_phase_V, Sconstant_phase_V;
			GetPhaseJacbi(i, x, dependencies[i], Sjac_phase_V, Sconstant_phase_V);
			Sjac_V.rows(Sjac_rowsShift, Sjac_rowsShift + Sjac_phase_V.n_elem - 1) =
				Sjac_phase_V;
			Sjac_rowsShift += Sjac_phase_V.n_elem;
			Sconstant_V.rows(Sconstant_rowsShift, Sconstant_rowsShift + Sconstant_phase_V.n_elem - 1) =
				Sconstant_phase_V;
			Sconstant_rowsShift += Sconstant_phase_V.n_elem;
			int	numcons = i_nstate*i_totnodes + i_npaths*i_totnodes + i_nevents;
			int numvars = i_nstate*(i_totnodes + 1) + i_ncontrol*i_totnodes + i_nparameters + 2;
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
			int stateindex_start = calculateData_->phase_indices[iphase]->state[0] - 1;
			int stateindex_end = calculateData_->phase_indices[iphase]->state[calculateData_->phase_indices[iphase]->state.size() - 1] - 1;
			vec state_vector = x.subvec(stateindex_start, stateindex_end);
			int controlindex_start, controlindex_end;
			vec control_vector;
			if (ncontrols > 0)
			{
				controlindex_start = calculateData_->phase_indices[iphase]->control[0] - 1;
				controlindex_end = calculateData_->phase_indices[iphase]->control[calculateData_->phase_indices[iphase]->control.size() - 1] - 1;
				control_vector = x.subvec(controlindex_start, controlindex_end);
			}
			int paralindex_start, paraindex_end;
			mat parameters;
			if (nparameters > 0)
			{
				paralindex_start = calculateData_->phase_indices[iphase]->parameter[0] - 1;
				paraindex_end = calculateData_->phase_indices[iphase]->parameter[calculateData_->phase_indices[iphase]->parameter.size() - 1] - 1;
				parameters = x.subvec(paralindex_start, paraindex_end);
			}
			int t0_index = calculateData_->phase_indices[iphase]->time[0] - 1;
			int tf_index = calculateData_->phase_indices[iphase]->time[1] - 1;
			double t0 = x(t0_index, 0);
			double tf = x(tf_index, 0);
			double tspan = tf - t0;
			vec t_radau = (calculateData_->PS[iphase]->Points + 1)*(tspan / 2.0) + t0;
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
		std::vector<std::vector<lp_index>> variable_indices = calculateData_->variable_indices;
		int linkrow = rowshift;//index start from zero
		for (int ipair = 0; ipair < optpro_->GetLinkageNum(); ipair++)
		{
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
			mySolLink.right_parameter_ = p_left;
			mySolLink.right_phase_num_ = right_index + 1;

			mySolLink.ipair = ipair + 1;
			//////////////////////////////////////////////////////////////////////////
			//Get Linkage Derivatives
			mat dLinkOut;
			derive_->DerivLink(mySolLink, dLinkOut);
			int nstatesLeft, ncontrolsLeft, nparametersLeft, npathLeft, neventLeft;
			int nstatesRight, ncontrolsRight, nparametersRight, npathRight, neventRight;
			optpro_->GetPhase(left_index)->get_optimal_info(
				nstatesLeft, ncontrolsLeft, nparametersLeft, npathLeft, neventLeft);
			optpro_->GetPhase(right_index)->get_optimal_info(
				nstatesRight, ncontrolsRight, nparametersRight, npathRight, neventRight);
			int nnodesLeft = optpro_->GetPhase(left_index)->GetTotalNodes();
			int nnodesRight = optpro_->GetPhase(right_index)->GetTotalNodes();
			int disc_left = nnodesLeft + 1;
			int disc_right = nnodesRight + 1;
			mat dLinkLeft = dLinkOut.cols(0, nstatesLeft + nparametersLeft - 1);
			mat DLink_xf_left = dLinkLeft.cols(0, nstatesLeft - 1);
			//DLink_xf_left.print("DLink_xf_left");
			mat DLink_p_left;
			if (nparametersLeft > 0)
			{
				DLink_p_left = dLinkLeft.cols(nstatesLeft, nstatesLeft + nparametersLeft - 1);
			}
			mat dLinkRight = dLinkOut.cols(dLinkLeft.n_cols, dLinkOut.n_cols - 1);
			mat DLink_x0_Right = dLinkRight.cols(0, nstatesRight - 1);
			//DLink_x0_Right.print("Dlink_x0_Right");
			mat DLink_p_Right;
			if (nparametersRight > 0)
			{
				DLink_p_Right = dLinkRight.cols(nstatesRight, nstatesRight + nparametersRight - 1);
			}

			int stateindexstart = calculateData_->phase_indices[left_index]->state[0] - 1;
			int statindexend = calculateData_->phase_indices[left_index]->state[calculateData_->phase_indices[left_index]->state.size() - 1] - 1;

			for (int jcol = 0; jcol < DLink_xf_left.n_cols; jcol++)
			{
				for (int irow = 0; irow < DLink_xf_left.n_rows; irow++)
				{
					double retvalue = DLink_xf_left(irow, jcol);
					Sjac_V.row(Sjac_rowsShift) =
						retvalue;
					Sjac_rowsShift++;

				}//for irow
			}//for jcol
			//////////////////////////////////////////////////////////////////////////
			int paraindexstart = 0, paraindexend = 0;
			if (nparametersLeft > 0)
			{
				paraindexstart = calculateData_->phase_indices[left_index]->parameter[0] - 1;
				paraindexend = calculateData_->phase_indices[left_index]->parameter[calculateData_->phase_indices[left_index]->parameter.size() - 1] - 1;
				for (int jcol = 0; jcol < DLink_p_left.n_cols; jcol++)
				{
					for (int irow = 0; irow < DLink_p_left.n_rows; irow++)
					{
						double retvalue = DLink_p_left(irow, jcol);
						Sjac_V.row(Sjac_rowsShift) =
							retvalue;
						Sjac_rowsShift++;
					}//for irow
				}//for jcol
			}

			stateindexstart = calculateData_->phase_indices[right_index]->state[0] - 1;
			statindexend = calculateData_->phase_indices[right_index]->state[calculateData_->phase_indices[right_index]->state.size() - 1] - 1;
			for (int jcol = 0; jcol < DLink_x0_Right.n_cols; jcol++)
			{
				for (int irow = 0; irow < DLink_x0_Right.n_rows; irow++)
				{
					double retvalue = DLink_x0_Right(irow, jcol);
					Sjac_V.row(Sjac_rowsShift) = retvalue;
					Sjac_rowsShift++;

				}//for irow
			}//for jcol
			//////////////////////////////////////////////////////////////////////////
			paraindexstart = 0; paraindexend = 0;
			if (nparametersRight > 0)
			{
				paraindexstart = calculateData_->phase_indices[right_index]->parameter[0] - 1;
				paraindexend = calculateData_->phase_indices[right_index]->parameter[calculateData_->phase_indices[right_index]->parameter.size() - 1] - 1;
				for (int jcol = 0; jcol < DLink_p_Right.n_cols; jcol++)
				{
					for (int irow = 0; irow < DLink_p_Right.n_rows; irow++)
					{
						double retvalue = DLink_p_Right(irow, jcol);
						Sjac_V.row(Sjac_rowsShift) =
							retvalue;
						Sjac_rowsShift++;

					}//for irow
				}//for jcol
			}
			calculateData_->phase_indices[right_index]->state[0];
			linkrow += DLink_x0_Right.n_rows;
		}// for ipair
	}
	void NLPWrapper::GetPhaseJacbi(int iphase, const vec& x_all, mat& idependencies, vec& Sjac_V, vec& Sconstant_V)
	{
		LP_DBG_START_FUN("NLPWrapper::GetPhaseJacbi()")
		int sumnodes = optpro_->GetPhase(iphase)->GetTotalNodes();
		size_t nstates = calculateData_->SIZES_[iphase][0];
		size_t ncontrols = calculateData_->SIZES_[iphase][1];
		size_t nparameters = calculateData_->SIZES_[iphase][2];
		size_t npaths = calculateData_->SIZES_[iphase][3];
		size_t nevents = calculateData_->SIZES_[iphase][4];
		int stateindex_start = calculateData_->phase_indices[iphase]->state[0] - 1;
		int stateindex_end = calculateData_->phase_indices[iphase]->state[calculateData_->phase_indices[iphase]->state.size() - 1] - 1;
		int controlindex_start = calculateData_->phase_indices[iphase]->control[0] - 1;
		int controlindex_end = calculateData_->phase_indices[iphase]->control[calculateData_->phase_indices[iphase]->control.size() - 1] - 1;
		vec state_vector = x_all.subvec(stateindex_start, stateindex_end);
		vec control_vector = x_all.subvec(controlindex_start, controlindex_end);
		int t0_index = calculateData_->phase_indices[iphase]->time[0] - 1;
		int tf_index = calculateData_->phase_indices[iphase]->time[1] - 1;
		double t0 = x_all(t0_index);
		double tf = x_all(tf_index);
		double tspan = tf - t0;
		vec t_radau = (calculateData_->PS[iphase]->Points + 1)*(tspan / 2.0) + t0;
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
			int para_index_start = calculateData_->phase_indices[iphase]->parameter[0] - 1;
			int para_index_end = calculateData_->phase_indices[iphase]->parameter[calculateData_->phase_indices[iphase]->parameter.size() - 1] - 1;
			parameters = x_all.rows(para_index_start, para_index_end);
		}
		//get derivatives of DAE function
		SolDae mySolDae;
		mySolDae.time_ = t_radau;
		mySolDae.state_ = state_radau;
		mySolDae.contol_ = control_radau;
		mySolDae.parameter_ = parameters;
		mySolDae.phase_num_ = iphase + 1;
		mat dDaeOut, dPathOut;
		derive_->DerivDae(mySolDae, dDaeOut, dPathOut);
		mat daeOut(sumnodes, nstates);
		mat pathOut(sumnodes, ncontrols);
		optimalFunction_->DaeFunction(mySolDae, daeOut, pathOut);
		mat dDae_time = zeros(sumnodes, nstates);
		std::vector<mat> dDae_state(nstates);
		std::vector<mat> dDae_control(ncontrols);
		std::vector<mat> dDae_parameter(nparameters);

		mat dPath_time;
		if (npaths > 0) dPath_time = zeros(sumnodes, npaths);
		std::vector<mat> dPath_state(nstates);
		std::vector<mat> dPath_control(ncontrols);
		std::vector<mat> dPath_parameter(nparameters);
		//Compute the Derivative  Respect to State
		for (int istate = 0; istate < nstates; istate++)
		{
			dDae_state[istate] = reshape(dDaeOut.col(istate), sumnodes, nstates);
			//dDae_state[istate].print("dDae_state[istate]");
			if (npaths > 0)
			{
				dPath_state[istate] = reshape(dPathOut.col(istate), sumnodes, npaths);
				//dPath_state[istate].print("dpath_state[istate]");
			}
		}
		//Compute the Derivative  Respect to Control
		for (int icontrol = 0; icontrol < ncontrols; icontrol++)
		{
			dDae_control[icontrol] = reshape(dDaeOut.col(nstates + icontrol), sumnodes, nstates);
			//dDae_control[icontrol].print("dDae_control[icontrol]");
			if (npaths > 0)
			{
				dPath_control[icontrol] = reshape(dPathOut.col(nstates + icontrol), sumnodes, npaths);
				//dPath_control[icontrol].print("dPath_control[icontrol]");
			}
		}
		//Compute the Derivative  Respect to time
		dDae_time = reshape(dDaeOut.col(nstates + ncontrols), sumnodes, nstates);
		//dDae_time.print("dDae_time");
		if (npaths > 0)
		{
			dPath_time = reshape(dPathOut.col(nstates + ncontrols), sumnodes, npaths);
			//dPath_time.print("dPath_time");
		}

		if (nparameters > 0)
		{
			for (int iparameter = 0; iparameter < nparameters; iparameter++)
			{
				dDae_parameter[iparameter] = reshape(dDaeOut.col(iparameter + nstates + ncontrols), sumnodes, nstates);
				if (npaths > 0)
				{
					dPath_parameter[iparameter] = reshape(dPathOut.col(iparameter + nstates + ncontrols), sumnodes, npaths);
				}
			}
		}//end if nparameters
		//Get derivatives of event function

		vec dEvent_t0;
		vec dEvent_tf;
		if (nevents > 0)
		{
			dEvent_t0 = zeros(nevents, 1);
			dEvent_tf = zeros(nevents, 1);
		}

		std::vector<vec> dEvent_x0(nstates);
		std::vector<vec> dEvent_xf(nstates);
		std::vector<vec> dEvent_parameter(nparameters);
		if (nevents > 0)
		{
			SolEvent mySolEvents;
			mySolEvents.initial_time_ = t0;
			mySolEvents.terminal_time_ = tf;
			mySolEvents.initial_state_ = x0;
			mySolEvents.terminal_state_ = xf;
			mySolEvents.parameter_ = parameters;
			mySolEvents.phase_num_ = iphase + 1;
			mat dEventOut;
			derive_->DerivEvent(mySolEvents, dEventOut);
			//dEventOut.print("dEventOut");
			//Get Derivatives of Event with Respect to Initial Time
			dEvent_t0 = dEventOut.col(nstates);
			//dEvent_t0.print("devent_t0");
			//Get Derivatives of Event with Respect to Initial Time
			dEvent_tf = dEventOut.col(2 * nstates + 1);
			//dEvent_tf.print("devent_tf");
			//Get Derivatives of Event with Respect to Initial and Terminal State
			for (int istate = 0; istate < nstates; istate++)
			{
				dEvent_x0[istate] = dEventOut.col(istate);
				//dEvent_x0[istate].print("devent_x0");
				dEvent_xf[istate] = dEventOut.col(istate + nstates + 1);
				//dEvent_xf[istate].print("devent_xf");
			}
			//Get Derivatives of Event with Respect to Parameters
			for (int iparameter = 0; iparameter < nparameters; iparameter++)
			{
				dEvent_parameter[iparameter] = dEventOut.col(iparameter + 2 * (nstates + 1));
			}
		}//if event>0
		//////////////////////////////////////////////////////////////////////////
		//Get Sparsity And Insert Derives of constraints
		int ndiffeqs = nstates;
		int disc_pts = sumnodes + 1;
		int numcons = ndiffeqs*(disc_pts - 1) + npaths*sumnodes + nevents;
		//find number of non-zero elements in dependencies (always including diagonal)
		for (int i = 0; i < ndiffeqs; i++)
		{
			idependencies(i, i) = 1.0;
		}
		int nDependancies = nnz(idependencies);
		//allocate memory for Jacobian
		Sjac_V = zeros(nDependancies*sumnodes + 2 * (ndiffeqs + npaths)*sumnodes +
			(ndiffeqs + npaths)*nparameters*sumnodes + nevents*(2 * nstates + nparameters + 2), 1);
		//find non-zeros in off diagonal D matrix
		vec DiffMatOffDiag_I, DiffMatOffDiag_J, DiffMatOffDiag_V;
		dsmatrix::Find(calculateData_->PS[iphase]->Doffdiag, DiffMatOffDiag_I, DiffMatOffDiag_J, DiffMatOffDiag_V);
		int nonZerosDiffMat = DiffMatOffDiag_I.n_elem;
		//allocate memory for constant derivatives
		Sconstant_V = zeros(ndiffeqs*nonZerosDiffMat, 1);

		//differential equations
		int rowstart, colstart;
		int Sjac_rowShift = 0, Sconstant_rowShift = 0;
		vec indexvector = linspace(0, sumnodes - 1, sumnodes);
		for (int i = 0; i < ndiffeqs; i++)
		{
			rowstart = i*sumnodes;
			for (int j = 0; j < nstates; j++)//inset df/dx
			{
				colstart = j*disc_pts;
				if (i == j)
				{
					// DdiagQuadBlock
					int Ddiag_ele_num = calculateData_->PS[iphase]->Diag.GetnRows();
					double*dDiagValues = new double[Ddiag_ele_num];
					double*PsDiag = calculateData_->PS[iphase]->Diag.GetSArrayConst()->GetDataPtr();
					for (int ii = 0; ii < Ddiag_ele_num; ii++)
					{
						dDiagValues[ii] = PsDiag[ii];
					}
					vec Ddiag(dDiagValues, Ddiag_ele_num, false);
					vec ret = Ddiag - dDae_state[j].col(i)*(tf - t0) / 2.0;
					Sjac_V.subvec(Sjac_rowShift, Sjac_rowShift + sumnodes - 1) = ret;
					Sjac_rowShift += sumnodes;
					//DiffMatOffDiag
					//dmatrix Sconstant_rows=diffmatrixindex_vector+i*nonZerosDiffMat;
					Sconstant_rowShift = i*nonZerosDiffMat;
					Sconstant_V.subvec(Sconstant_rowShift, Sconstant_rowShift + nonZerosDiffMat - 1) = DiffMatOffDiag_V;
				}
				else
				{
					if (idependencies(i, j) == 1.0)
					{
						//OffDiag_State
						vec ret = dDae_state[j].col(i)*(tf - t0) / 2.0;
						Sjac_V.rows(Sjac_rowShift, Sjac_rowShift + sumnodes - 1) = -ret;
						Sjac_rowShift += sumnodes;
					}
				}
			}//for j
			//insert df/du
			int colshift = nstates*disc_pts;
			for (int j = 0; j < ncontrols; j++)
			{
				colstart = colshift + j*sumnodes;
				if (idependencies(i, j + nstates) == 1.0)
				{
					//Diffeq_Block_Control
					vec ret = dDae_control[j].col(i)*(tf - t0) / 2.0;
					Sjac_V.rows(Sjac_rowShift, Sjac_rowShift + sumnodes - 1) = -ret;
					Sjac_rowShift += sumnodes;
				}
			}//for j
			//////////////////////////////////////////////////////////////////////////
			//insert df/dt0
			colshift += ncontrols*sumnodes;
			colstart = colshift;
			vec ret = daeOut.col(i)*(0.5);
			vec ret2 = -(calculateData_->PS[iphase]->Points*0.5) + 0.5;
			ret -= ret2 % (dDae_time.col(i)*(tf - t0) / 2.0);
			Sjac_V.rows(Sjac_rowShift, Sjac_rowShift + sumnodes - 1) = ret;
			Sjac_rowShift += sumnodes;
			colshift++;
			// insert df/dtf
			colstart = colshift;
			ret = -daeOut.col(i)*(0.5);
			ret2 = (calculateData_->PS[iphase]->Points*0.5) + 0.5;
			ret = ret + ret2 % (dDae_time.col(i)*(tf - t0) / 2.0);
			Sjac_V.rows(Sjac_rowShift, Sjac_rowShift + sumnodes - 1) = ret;
			Sjac_rowShift += sumnodes;
			colshift++;
			//insert df/dpara
			for (int j = 0; j < nparameters; j++)
			{
				colstart = colshift + j;
				vec ret = -dDae_parameter[j].col(i)*(tf - t0) / 2.0;
				Sjac_V.rows(Sjac_rowShift, Sjac_rowShift + sumnodes - 1) = ret;
				Sjac_rowShift += sumnodes;
			}//for j
		}//for i
		//path contraints
		int rowshift = ndiffeqs*(sumnodes);
		for (int i = 0; i < npaths; i++)
		{
			rowstart = rowshift + i*sumnodes;
			for (int j = 0; j < nstates; j++)
			{
				colstart = j*disc_pts;
				if (idependencies(i + ndiffeqs, j) == 1.0)
				{
					//Insert dc/dx
					Sjac_V.rows(Sjac_rowShift, Sjac_rowShift + sumnodes - 1) = dPath_state[j].col(i);
					Sjac_rowShift += sumnodes;
				}
			}//for j
			int colshift = nstates*disc_pts;
			for (int j = 0; j < ncontrols; j++)
			{
				colstart = colshift + j*sumnodes;
				if (idependencies(i + ndiffeqs, j + nstates))
				{
					//insert dc/du
					Sjac_V.rows(Sjac_rowShift, Sjac_rowShift + sumnodes - 1) = dPath_control[j].col(i);
					Sjac_rowShift += sumnodes;
				}
			}//for j
			//////////////////////////////////////////////////////////////////////////
			//insert dc/dt0
			colshift += ncontrols*sumnodes;
			colstart = colshift;
			vec ret = dPath_time.col(i);
			vec ret2 = -(calculateData_->PS[iphase]->Points*0.5) + 0.5;
			Sjac_V.rows(Sjac_rowShift, Sjac_rowShift + sumnodes - 1) = ret2 % (ret);
			Sjac_rowShift += sumnodes;
			colshift++;
			// insert dc/dtf
			colstart = colshift;
			ret = dPath_time.col(i);
			ret2 = (calculateData_->PS[iphase]->Points*0.5) + 0.5;
			Sjac_V.rows(Sjac_rowShift, Sjac_rowShift + sumnodes - 1) = ret2 % (ret);
			Sjac_rowShift += sumnodes;
			colshift++;
			//insert dc/dpara
			for (int j = 0; j < nparameters; j++)
			{
				colstart = colshift + j;
				Sjac_V.rows(Sjac_rowShift, Sjac_rowShift + sumnodes - 1) = dDae_parameter[j].col(i);
				Sjac_rowShift += sumnodes;
			}//for j
		}//for i
		rowshift += npaths*sumnodes;
		//Event Constraints
		rowshift = ndiffeqs*(sumnodes)+npaths*sumnodes;
		vec event_indces_init = linspace(0, nstates - 1, nstates)*(sumnodes + 1);
		vec event_indces_term = event_indces_init + (sumnodes);
		vec event_indeces(nstates * 2);
		if (nevents > 0)
		{
			event_indeces.rows(0, nstates - 1) = event_indces_init;
			event_indeces.rows(nstates, 2 * nstates - 1) = event_indces_term;
		}

		for (int i = 0; i < nevents; i++)
		{
			int row = rowshift + i;
			//devent/dx
			for (int j = 0; j < nstates; j++)
			{
				int col0 = sumnodes*j + j;
				int colf = sumnodes*(j + 1) + j;
				Sjac_V.row(Sjac_rowShift) = dEvent_x0[j](i);
				Sjac_rowShift++;
				Sjac_V.row(Sjac_rowShift) = dEvent_xf[j](i);
				Sjac_rowShift++;
			}
			//insert devent/dt0
			int cols = nstates*(sumnodes + 1) + ncontrols*sumnodes;
			Sjac_V.row(Sjac_rowShift) = dEvent_t0(i, 0);
			Sjac_rowShift++;
			//insert devent/df
			cols++;
			Sjac_V.row(Sjac_rowShift) = dEvent_tf(i, 0);
			Sjac_rowShift++;
			//insert dEvent/dparam
			for (int j = 0; j < nparameters; j++)
			{
				cols++;
				Sjac_V.row(Sjac_rowShift) = dEvent_parameter[j](i);
				Sjac_rowShift++;
			}
		}//for ievent
	}
	double NLPWrapper::GetObjFun(const vec& scaled_x)
	{
		LP_DBG_START_FUN("NLPWrapper::GetObjFun()")
		vec x = scaled_x;
		if (calculateData_->autoscale)
		{
			x = (scaled_x - calculateData_->varshift) / calculateData_->varscale;
		}
		double cost = 0.0;
		for (int i = 0; i < calculateData_->numphases_; i++)
		{
			int nstate = calculateData_->SIZES_[i][0];
			int ncontrols = calculateData_->SIZES_[i][1];
			int nparameters = calculateData_->SIZES_[i][2];
			int npaths = calculateData_->SIZES_[i][3];
			int nevents = calculateData_->SIZES_[i][4];

			int state_index_start = calculateData_->phase_indices[i]->state[0] - 1;
			int state_index_end = calculateData_->phase_indices[i]->state[calculateData_->phase_indices[i]->state.size() - 1] - 1;
			vec state_vector = x.rows(state_index_start, state_index_end);

			int control_index_start = calculateData_->phase_indices[i]->control[0] - 1;
			int control_index_end = calculateData_->phase_indices[i]->control[calculateData_->phase_indices[i]->control.size() - 1] - 1;
			vec control_vector = x.rows(control_index_start, control_index_end);

			double t0 = x(calculateData_->phase_indices[i]->time[0] - 1, 0);
			double tf = x(calculateData_->phase_indices[i]->time[1] - 1, 0);
			double tspan = tf - t0;
			vec t_radau = (calculateData_->PS[i]->Points + 1)*(tspan / 2.0) + t0;

			mat state_matrix = reshape(state_vector,
				optpro_->GetPhase(i)->GetTotalNodes() + 1, nstate);
			mat state_radau = state_matrix.rows(0, state_matrix.n_rows - 2);

			vec x0 = trans(state_matrix.row(0));
			vec xf = trans(state_matrix.row(state_matrix.n_rows - 1));
			mat control_radau = reshape(control_vector,
				optpro_->GetPhase(i)->GetTotalNodes(), ncontrols);
			mat parameter_vec;
			if (nparameters > 0)
			{
				int para_index_start = calculateData_->phase_indices[i]->parameter[0] - 1;
				int para_index_end = calculateData_->phase_indices[i]->parameter[calculateData_->phase_indices[i]->parameter.size() - 1] - 1;
				parameter_vec = x.rows(para_index_start, para_index_end);
			}

			//Get Cost;
			SolCost mySolCost;
			mySolCost.initial_time_ = t0;
			mySolCost.initial_state_ = x0;
			//x0.print("x0");
			mySolCost.terminal_time_ = tf;
			mySolCost.terminal_state_ = xf;
			//xf.print("xf");
			mySolCost.time_ = t_radau;
			//t_radau.print("tradau");
			mySolCost.state_ = state_radau;
			//state_radau.print("state_radau");
			mySolCost.control_ = control_radau;
			//	control_radau.print("control_radau");
			mySolCost.parameter_ = parameter_vec;
			mySolCost.phase_num_ = i + 1;

			double temMayer = 0.0;
			vec 	temLagrange;
			optimalFunction_->MayerCost(mySolCost, temMayer);
			cost += temMayer;
			optimalFunction_->LagrangeCost(mySolCost, temLagrange);
			vec integrand = trans(calculateData_->PS[i]->Weights)*temLagrange*(tspan / 2.0);
			cost += integrand(0);
			if (calculateData_->autoscale)
			{
				cost *= calculateData_->objScale[i];
			}
		}//end for
		return cost;
	}
	void NLPWrapper::GetObjGrad(const vec& scaled_x, vec& grad_f)
	{
		LP_DBG_START_FUN("NLPWrapper::GetObjGrad()")
		vec x = scaled_x;
		if (calculateData_->autoscale)
		{
			x = (scaled_x - calculateData_->varshift) / calculateData_->varscale;
		}
		 grad_f = zeros(calculateData_->varbounds_max.size(), 1);
		int grad_shift = 0;
		for (int iphase = 0; iphase < optpro_->GetPhaseNum(); iphase++)
		{
			int sumnodes = optpro_->GetPhase(iphase)->GetTotalNodes();
			int nstates = optpro_->GetPhase(iphase)->GetstateMin().size();
			int ncontrols = optpro_->GetPhase(iphase)->GetcontrolMin().size();
			int nparameters = optpro_->GetPhase(iphase)->GetparameterMin().size();
			int npaths = optpro_->GetPhase(iphase)->GetpathMin().size();
			int nevents = optpro_->GetPhase(iphase)->GeteventMin().size();
			int stateindex_start = calculateData_->phase_indices[iphase]->state[0] - 1;
			int stateindex_end = calculateData_->phase_indices[iphase]->state[calculateData_->phase_indices[iphase]->state.size() - 1] - 1;
			int controlindex_start = calculateData_->phase_indices[iphase]->control[0] - 1;
			int controlindex_end = calculateData_->phase_indices[iphase]->control[calculateData_->phase_indices[iphase]->control.size() - 1] - 1;
			vec state_vector = x.subvec(stateindex_start, stateindex_end);
			vec control_vector = x.subvec(controlindex_start, controlindex_end);
			int t0_index = calculateData_->phase_indices[iphase]->time[0] - 1;
			int tf_index = calculateData_->phase_indices[iphase]->time[1] - 1;
			double t0 = x(t0_index, 0);
			double tf = x(tf_index, 0);
			double tspan = tf - t0;
			vec t_radau = (calculateData_->PS[iphase]->Points + 1)*(tspan / 2.0) + t0;
			mat state_matrix = reshape(state_vector, sumnodes + 1, nstates);
			mat state_radau = state_matrix.rows(0, state_matrix.n_rows - 2);//index start form 0
			vec x0 = trans(state_matrix.row(0));
			vec xf = trans(state_matrix.row(state_matrix.n_rows - 1));
			mat control_radau = reshape(control_vector, sumnodes, ncontrols);
			mat parameters;
			SolCost mySolCost;
			mySolCost.initial_time_ = t0;
			mySolCost.initial_state_ = x0;
			mySolCost.terminal_time_ = tf;
			mySolCost.terminal_state_ = xf;
			mySolCost.time_ = t_radau;
			mySolCost.state_ = state_radau;
			mySolCost.control_ = control_radau;
			mySolCost.parameter_ = parameters;
			mySolCost.phase_num_ = iphase + 1;
			rowvec dmayerOut;
			mat dLagrangeOut;
			vec	LagrangeOut;
			optimalFunction_->LagrangeCost(mySolCost, LagrangeOut);
			derive_->DerivMayer(mySolCost, dmayerOut);
			derive_->DerivLagrange(mySolCost, dLagrangeOut);
			double dMayer_t0 = 0.0;
			double dMayer_tf = 0.0;
			rowvec dMayer_x0 = zeros(1, nstates);
			rowvec dMayer_xf = zeros(1, nstates);
			rowvec dMayer_para;
			if (nparameters > 0)
			{
				dMayer_para = zeros(1, nparameters);
			}
			mat dLagrange_state = zeros(sumnodes, nstates);
			vec dLagrange_time = zeros(sumnodes, 1);
			mat dLagrange_control, dLagrange_para;
			if (ncontrols > 0)dLagrange_control = zeros(sumnodes, ncontrols);
			if (nparameters > 0)dLagrange_para = zeros(sumnodes, nparameters);
			//Get Derivatives of mayer cost 
			if (dmayerOut.n_elem > 0)
			{

				dMayer_t0 = dmayerOut(0, nstates);
				dMayer_tf = dmayerOut(0, 2 * nstates + 1);
				for (int istate = 0; istate < nstates; istate++)
				{
					dMayer_x0(istate) = dmayerOut(istate);
					dMayer_xf(istate) = dmayerOut(istate + nstates + 1);
				}
				for (int ipara = 0; ipara < nparameters; ipara++)
				{
					dMayer_para(ipara) = dmayerOut(ipara + 2 * (nstates + 1));
				}
			}//if dmayerOut.length>0
			//Get Derivatives of Lagrange cost 
			if (dLagrangeOut.n_elem > 0)
			{
				dLagrange_time = dLagrangeOut.col(dLagrangeOut.n_cols - 1);
				//dLagrange_time.print("L_time");
				for (int istate = 0; istate < nstates; istate++)
				{
					dLagrange_state.col(istate) = dLagrangeOut.col(istate);
				}
				for (int icontrol = 0; icontrol < ncontrols; icontrol++)
				{
					dLagrange_control.col(icontrol) = dLagrangeOut.col(icontrol + nstates);
				}
				for (int ipara = 0; ipara < nparameters; ipara++)
				{
					dLagrange_para.col(ipara) = dLagrangeOut.col(ipara + nstates + ncontrols);
				}
			}//if dlagrangeout.length>0

			//////////////////////////////////////////////////////////////////////////
			//insert Cost derivatives
			int colshift = 0;
			rowvec Jcost = zeros(1, nstates*(sumnodes + 1) + ncontrols*sumnodes + nparameters + 2);
			for (int j = 0; j < nstates; j++)
			{
				int col0 = sumnodes*(j)+j;
				int colf = sumnodes*(j + 1) + j;
				vec index_vector = linspace(0, sumnodes - 1, sumnodes);
				Jcost(0, col0) = dMayer_x0(0, j);
				vec ret = calculateData_->PS[iphase]->Weights*tspan / 2.0;
				ret %= (dLagrange_state.col(j));
				Jcost.subvec(col0, colf - 1) = trans(ret);
				Jcost(colf) = dMayer_xf(j);
			}
			//dCost/du
			colshift += nstates*(sumnodes + 1);
			for (int j = 0; j < ncontrols; j++)
			{
				int colstart = colshift + j*sumnodes;
				vec ret = calculateData_->PS[iphase]->Weights*tspan / 2.0;
				ret %= (dLagrange_control.col(j));
				Jcost.subvec(colstart, colstart + sumnodes - 1) = trans(ret);
			}
			colshift += ncontrols*sumnodes;
			int t0_col = colshift;

			//dCost/dt0
			mat ret = calculateData_->PS[iphase]->Weights*(-0.5);
			ret = trans(ret)*(LagrangeOut);
			//ret.print("ret");
			mat ret2 = calculateData_->PS[iphase]->Weights*(tspan / 2.0);
			ret2 = trans(ret2)*diagmat(dLagrange_time);
			//ret2.print("ret2");
			//dLagrange_time.diag().print("diag");
			mat ret3 = calculateData_->PS[iphase]->Points*(-0.5) + 0.5;
			ret3 = ret2*ret3 + dMayer_t0 + ret(0);
			Jcost(t0_col) = ret3(0);
			int tf_col = colshift + 1;
			//dCost/dtf
			ret = calculateData_->PS[iphase]->Weights*(0.5);
			ret = trans(ret)*(LagrangeOut);
			ret2 = calculateData_->PS[iphase]->Weights*(tspan / 2.0);
			ret2 = trans(ret2)*diagmat(dLagrange_time);
			ret3 = calculateData_->PS[iphase]->Points*(0.5) + 0.5;
			ret3 *= ret2;
			Jcost(tf_col) = dMayer_tf + ret(0) + ret3(0);
			colshift = colshift + 1;
			//dCost/dpara
			for (int j = 0; j < nparameters; j++)
			{
				int paracol = colshift + j;
				ret = calculateData_->PS[iphase]->Weights*(tspan / 2.0)*dLagrange_para.col(j);
				Jcost(paracol) = dMayer_para(j) + ret(0);
			}
			grad_f.subvec(grad_shift, grad_shift + Jcost.n_elem - 1) = trans(Jcost);
			grad_shift += Jcost.n_elem;
			//Jcost.print("Jcost");
		}//for iphase
		if (calculateData_->autoscale)
		{
			grad_f %= objscale_;
		}
	}

	void NLPWrapper::GetPhaseSparsity(int iphase, mat& idependencies, vec& Sjac_I, vec& Sjac_J, vec &Sconstant_I, vec &Sconstant_J)
	{
		LP_DBG_START_FUN("NLPWrapper::GetPhaseSparsity()")
		int sumnodes = optpro_->GetPhase(iphase)->GetTotalNodes();
		int nstates = optpro_->GetPhase(iphase)->GetstateMin().size();
		int ncontrols = optpro_->GetPhase(iphase)->GetcontrolMin().size();
		int nparameters = optpro_->GetPhase(iphase)->GetparameterMin().size();
		int npaths = optpro_->GetPhase(iphase)->GetpathMin().size();
		int nevents = optpro_->GetPhase(iphase)->GeteventMin().size();
		int stateindex_start = calculateData_->phase_indices[iphase]->state[0] - 1;
		int stateindex_end = calculateData_->phase_indices[iphase]->state[calculateData_->phase_indices[iphase]->state.size() - 1] - 1;
		int controlindex_start = calculateData_->phase_indices[iphase]->control[0] - 1;
		int controlindex_end = calculateData_->phase_indices[iphase]->control[calculateData_->phase_indices[iphase]->control.size() - 1] - 1;
		int t0_index = calculateData_->phase_indices[iphase]->time[0] - 1;
		int tf_index = calculateData_->phase_indices[iphase]->time[1] - 1;


		//Get Sparsity And Insert Derives of constraints
		int ndiffeqs = nstates;
		int disc_pts = sumnodes + 1;
		int numcons = ndiffeqs*(disc_pts - 1) + npaths*sumnodes + nevents;
		//find number of non-zero elements in dependencies (always including diagonal)
		for (int i = 0; i < ndiffeqs; i++)
		{
			idependencies(i, i) = 1.0;
		}
		int nDependancies = nnz(idependencies);
		//allocate memory for Jacobian
		Sjac_I = zeros(nDependancies*sumnodes + 2 * (ndiffeqs + npaths)*sumnodes +
			(ndiffeqs + npaths)*nparameters*sumnodes + nevents*(2 * nstates + nparameters + 2), 1);
		Sjac_J = zeros(nDependancies*sumnodes + 2 * (ndiffeqs + npaths)*sumnodes +
			(ndiffeqs + npaths)*nparameters*sumnodes + nevents*(2 * nstates + nparameters + 2), 1);
		//find non-zeros in off diagonal D matrix
		vec DiffMatOffDiag_I, DiffMatOffDiag_J, DiffMatOffDiag_V;
		dsmatrix::Find(calculateData_->PS[iphase]->Doffdiag, DiffMatOffDiag_I, DiffMatOffDiag_J, DiffMatOffDiag_V);
		int nonZerosDiffMat = DiffMatOffDiag_I.n_elem;
		//allocate memory for constant derivatives
		Sconstant_I = zeros(ndiffeqs*nonZerosDiffMat, 1);
		Sconstant_J = zeros(ndiffeqs*nonZerosDiffMat, 1);

		//differential equations
		int rowstart, colstart;
		int Sjac_rowShift = 0, Sconstant_rowShift = 0;
		vec indexvector = linspace(0, sumnodes - 1, sumnodes);
		for (int i = 0; i < ndiffeqs; i++)
		{
			rowstart = i*sumnodes;
			for (int j = 0; j < nstates; j++)//inset df/dx
			{
				colstart = j*disc_pts;
				if (i == j)
				{
					// DdiagQuadBlock
					Sjac_I.subvec(Sjac_rowShift, Sjac_rowShift + sumnodes - 1) = indexvector + rowstart;
					Sjac_J.subvec(Sjac_rowShift, Sjac_rowShift + sumnodes - 1) = indexvector + colstart;
					Sjac_rowShift += sumnodes;
					//DiffMatOffDiag
					//dmatrix Sconstant_rows=diffmatrixindex_vector+i*nonZerosDiffMat;
					Sconstant_rowShift = i*nonZerosDiffMat;
					Sconstant_I.subvec(Sconstant_rowShift, Sconstant_rowShift + nonZerosDiffMat - 1) = DiffMatOffDiag_I + rowstart;
					Sconstant_J.subvec(Sconstant_rowShift, Sconstant_rowShift + nonZerosDiffMat - 1) = DiffMatOffDiag_J + colstart;
				}
				else
				{
					if (idependencies(i, j) == 1.0)
					{
						//OffDiag_State
						Sjac_I.rows(Sjac_rowShift, Sjac_rowShift + sumnodes - 1) = indexvector + rowstart;
						Sjac_J.rows(Sjac_rowShift, Sjac_rowShift + sumnodes - 1) = indexvector + colstart;
						Sjac_rowShift += sumnodes;
					}
				}
			}//for j
			//insert df/du
			int colshift = nstates*disc_pts;
			for (int j = 0; j < ncontrols; j++)
			{
				colstart = colshift + j*sumnodes;
				if (idependencies(i, j + nstates) == 1.0)
				{
					//Diffeq_Block_Control
					Sjac_I.rows(Sjac_rowShift, Sjac_rowShift + sumnodes - 1) = indexvector + rowstart;
					Sjac_J.rows(Sjac_rowShift, Sjac_rowShift + sumnodes - 1) = indexvector + colstart;
					Sjac_rowShift += sumnodes;
				}
			}//for j
			//////////////////////////////////////////////////////////////////////////
			//insert df/dt0
			colshift += ncontrols*sumnodes;
			colstart = colshift;
			Sjac_I.rows(Sjac_rowShift, Sjac_rowShift + sumnodes - 1) = indexvector + rowstart;
			Sjac_J.rows(Sjac_rowShift, Sjac_rowShift + sumnodes - 1).fill(colstart);
			Sjac_rowShift += sumnodes;
			colshift++;
			// insert df/dtf
			colstart = colshift;
			Sjac_I.rows(Sjac_rowShift, Sjac_rowShift + sumnodes - 1) = indexvector + rowstart;
			Sjac_J.rows(Sjac_rowShift, Sjac_rowShift + sumnodes - 1).fill(colstart);
			Sjac_rowShift += sumnodes;
			colshift++;
			//insert df/dpara
			for (int j = 0; j < nparameters; j++)
			{
				colstart = colshift + j;
				Sjac_I.rows(Sjac_rowShift, Sjac_rowShift + sumnodes - 1) = indexvector + rowstart;
				Sjac_J.rows(Sjac_rowShift, Sjac_rowShift + sumnodes - 1) = indexvector + colstart;
				Sjac_rowShift += sumnodes;
			}//for j
		}//for i
		//path contraints
		int rowshift = ndiffeqs*(sumnodes);
		for (int i = 0; i < npaths; i++)
		{
			rowstart = rowshift + i*sumnodes;
			for (int j = 0; j < nstates; j++)
			{
				colstart = j*disc_pts;
				if (idependencies(i + ndiffeqs, j) == 1.0)
				{
					//Insert dc/dx
					Sjac_I.rows(Sjac_rowShift, Sjac_rowShift + sumnodes - 1) = indexvector + rowstart;
					Sjac_J.rows(Sjac_rowShift, Sjac_rowShift + sumnodes - 1) = indexvector + colstart;
					Sjac_rowShift += sumnodes;
				}
			}//for j
			int colshift = nstates*disc_pts;
			for (int j = 0; j < ncontrols; j++)
			{
				colstart = colshift + j*sumnodes;
				if (idependencies(i + ndiffeqs, j + nstates))
				{
					//insert dc/du
					Sjac_I.rows(Sjac_rowShift, Sjac_rowShift + sumnodes - 1) = indexvector + rowstart;
					Sjac_J.rows(Sjac_rowShift, Sjac_rowShift + sumnodes - 1) = indexvector + colstart;
					Sjac_rowShift += sumnodes;
				}
			}//for j
			//////////////////////////////////////////////////////////////////////////
			//insert dc/dt0
			colshift += ncontrols*sumnodes;
			colstart = colshift;
			Sjac_I.rows(Sjac_rowShift, Sjac_rowShift + sumnodes - 1) = indexvector + rowstart;
			Sjac_J.rows(Sjac_rowShift, Sjac_rowShift + sumnodes - 1).fill(colstart);
			Sjac_rowShift += sumnodes;
			colshift++;
			// insert dc/dtf
			colstart = colshift;
			Sjac_I.rows(Sjac_rowShift, Sjac_rowShift + sumnodes - 1) = indexvector + rowstart;
			Sjac_J.rows(Sjac_rowShift, Sjac_rowShift + sumnodes - 1).fill(colstart);
			Sjac_rowShift += sumnodes;
			colshift++;
			//insert dc/dpara
			for (int j = 0; j < nparameters; j++)
			{
				colstart = colshift + j;
				Sjac_I.rows(Sjac_rowShift, Sjac_rowShift + sumnodes - 1) = indexvector + rowstart;
				Sjac_J.rows(Sjac_rowShift, Sjac_rowShift + sumnodes - 1) = indexvector + colstart;
				Sjac_rowShift += sumnodes;
			}//for j
		}//for i
		rowshift += npaths*sumnodes;
		//Event Constraints
		rowshift = ndiffeqs*(sumnodes)+npaths*sumnodes;
		vec event_indces_init = linspace(0, nstates - 1, nstates)*(sumnodes + 1);
		vec event_indces_term = event_indces_init + (sumnodes);
		vec event_indeces(nstates * 2);
		if (nevents > 0)
		{
			event_indeces.rows(0, nstates - 1) = event_indces_init;
			event_indeces.rows(nstates, 2 * nstates - 1) = event_indces_term;
		}

		for (int i = 0; i < nevents; i++)
		{
			int row = rowshift + i;
			//devent/dx
			for (int j = 0; j < nstates; j++)
			{
				int col0 = sumnodes*j + j;
				int colf = sumnodes*(j + 1) + j;
				Sjac_I.row(Sjac_rowShift) = row;
				Sjac_J.row(Sjac_rowShift) = col0;
				Sjac_rowShift++;
				Sjac_I.row(Sjac_rowShift) = row;
				Sjac_J.row(Sjac_rowShift) = colf;
				Sjac_rowShift++;
			}
			//insert devent/dt0
			int cols = nstates*(sumnodes + 1) + ncontrols*sumnodes;
			Sjac_I.row(Sjac_rowShift) = row;
			Sjac_J.row(Sjac_rowShift) = cols;
			Sjac_rowShift++;
			//insert devent/df
			cols++;
			Sjac_I.row(Sjac_rowShift) = row;
			Sjac_J.row(Sjac_rowShift) = cols;
			Sjac_rowShift++;
			//insert dEvent/dparam
			for (int j = 0; j < nparameters; j++)
			{
				cols++;
				Sjac_I.row(Sjac_rowShift) = row;
				Sjac_J.row(Sjac_rowShift) = cols;
				Sjac_rowShift++;
			}
		}//for ievent
	}

	void NLPWrapper::GetWholeSparsity(vec& Sjac_I, vec& Sjac_J, vec &Sconstant_I, vec &Sconstant_J)
	{
		LP_DBG_START_FUN("NLPWrapper::GetWholeSparsity()")
		vector<mat> dependencies(optpro_->GetPhaseNum());

		for (int itor = 0; itor < dependencies.size(); itor++)
		{
			int i_nstate = optpro_->GetPhase(itor)->GetstateMin().size();
			int i_ncontrol = optpro_->GetPhase(itor)->GetcontrolMin().size();
			int i_npaths = optpro_->GetPhase(itor)->GetpathMin().size();
			dependencies[itor] = zeros(i_nstate + i_npaths, i_nstate + i_ncontrol);
		}

		//find total non-zero elements in Jacobian and constant derivatives
		int  nonZerosSjac = 0;
		int  nonZerosSconstant = 0;
		int  nonZerosSnonconstant = 0;

		for (int i = 0; i < optpro_->GetPhaseNum(); i++)
		{
			int i_nstate = optpro_->GetPhase(i)->GetstateMin().size();
			int i_ncontrol = optpro_->GetPhase(i)->GetcontrolMin().size();
			int i_nparameters = optpro_->GetPhase(i)->GetparameterMin().size();
			int i_npaths = optpro_->GetPhase(i)->GetpathMin().size();
			int i_nevents = optpro_->GetPhase(i)->GeteventMin().size();
			int i_totnodes = optpro_->GetPhase(i)->GetTotalNodes();
			// find number of non-zero elements in dependencies (always including diagonal)


			//////////////////////////////////////////////////////////////////////////
			//only when use Analytical Derives 
			dependencies[i].fill(1.0);
			////////////////////////////////////////////////////////////////////////
			mat temDependencies = dependencies[i];
			for (int j = 0; j < i_nstate; j++)
			{
				temDependencies(j, j) = 1.0;
			}

			int nDependancies = nnz(temDependencies);

			nonZerosSjac = nonZerosSjac + nDependancies*i_totnodes + 2 * (i_nstate + i_npaths)*i_totnodes +
				(i_nstate + i_npaths)*i_nparameters*i_totnodes + i_nevents*(2 * i_nstate + i_nparameters + 2);
			vec DiffMatOffDiag_I, DiffMatOffDiag_J, DiffMatOffDiag_V;
			dsmatrix::Find(calculateData_->PS[i]->Doffdiag, DiffMatOffDiag_I, DiffMatOffDiag_J, DiffMatOffDiag_V);
			int nonZerosDiffMat = DiffMatOffDiag_I.n_elem;
			nonZerosSconstant = nonZerosSconstant + nonZerosDiffMat*i_nstate;
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
			nonZerosSjac += numlinks*(nstates_left + nparameter_left + nstates_right + nparameter_right);
		}//for ipair
		//allocate Memory for Jacobian
		Sjac_I = zeros(nonZerosSjac, 1);
		Sjac_J = zeros(nonZerosSjac, 1);
		Sconstant_I = zeros(nonZerosSconstant, 1);
		Sconstant_J = zeros(nonZerosSconstant, 1);
		int Sjac_rowsShift = 0, Sconstant_rowsShift = 0;
		int rowshift = 0;
		int colshift = 0;
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
			vec Sjac_phase_I, Sjac_phase_J, Sconstant_phase_I, Sconstant_phase_J;
			GetPhaseSparsity(i, dependencies[i], Sjac_phase_I, Sjac_phase_J, Sconstant_phase_I, Sconstant_phase_J);
			Sjac_I.rows(Sjac_rowsShift, Sjac_rowsShift + Sjac_phase_I.n_elem - 1) =
				Sjac_phase_I + rowshift;
			Sjac_J.rows(Sjac_rowsShift, Sjac_rowsShift + Sjac_phase_I.n_elem - 1) =
				Sjac_phase_J + colshift;

			Sjac_rowsShift += Sjac_phase_I.n_elem;
			Sconstant_I.rows(Sconstant_rowsShift, Sconstant_rowsShift + Sconstant_phase_I.n_elem - 1) =
				Sconstant_phase_I + rowshift;
			Sconstant_J.rows(Sconstant_rowsShift, Sconstant_rowsShift + Sconstant_phase_I.n_elem - 1) =
				Sconstant_phase_J + colshift;
			Sconstant_rowsShift += Sconstant_phase_I.n_elem;
			int	numcons = i_nstate*i_totnodes + i_npaths*i_totnodes + i_nevents;
			int numvars = i_nstate*(i_totnodes + 1) + i_ncontrol*i_totnodes + i_nparameters + 2;
			rowshift += numcons;
			colshift += numvars;
		}
		//derive of linkage and it's Sparsity Pattern
		//int sjacRowsStart=calculateData_->conbounds_min.size()-optpro_->GetPhaseNum()-1;//count form zero
		std::vector<SolCost> solcCostTotal(optpro_->GetPhaseNum());
		/*	for (int iphase = 0; iphase < optpro_->GetPhaseNum(); iphase++)
			{
			int sumnodes = optpro_->GetPhase(iphase)->GetTotalNodes();
			int nstates = optpro_->GetPhase(iphase)->GetstateMin().size();
			int ncontrols = optpro_->GetPhase(iphase)->GetcontrolMin().size();
			int nparameters = optpro_->GetPhase(iphase)->GetparameterMin().size();
			int npaths = optpro_->GetPhase(iphase)->GetpathMin().size();
			int nevents = optpro_->GetPhase(iphase)->GeteventMin().size();
			int stateindex_start = calculateData_->phase_indices[iphase]->state[0] - 1;
			int stateindex_end = calculateData_->phase_indices[iphase]->state[calculateData_->phase_indices[iphase]->state.size() - 1] - 1;

			int t0_index = calculateData_->phase_indices[iphase]->time[0] - 1;
			int tf_index = calculateData_->phase_indices[iphase]->time[1] - 1;

			}*/
		std::vector<std::vector<lp_index>> variable_indices = calculateData_->variable_indices;
		int linkrow = rowshift;//index start from zero
		for (int ipair = 0; ipair < optpro_->GetLinkageNum(); ipair++)
		{
			int nlinks = optpro_->GetLinkage(ipair)->GetLinkageMin().size();
			int left_index = optpro_->GetLinkage(ipair)->LeftPhase();
			int right_index = optpro_->GetLinkage(ipair)->RightPhase();

			//////////////////////////////////////////////////////////////////////////
			//Get Linkage Derivatives

			int nstatesLeft, ncontrolsLeft, nparametersLeft, npathLeft, neventLeft;
			int nstatesRight, ncontrolsRight, nparametersRight, npathRight, neventRight;
			optpro_->GetPhase(left_index)->get_optimal_info(
				nstatesLeft, ncontrolsLeft, nparametersLeft, npathLeft, neventLeft);
			optpro_->GetPhase(right_index)->get_optimal_info(
				nstatesRight, ncontrolsRight, nparametersRight, npathRight, neventRight);
			mat dLinkOut;
			//dLinkOut = zeros(xf_left.n_elem, xf_left.n_elem + x0_right.n_elem);
			dLinkOut = zeros(nlinks, nstatesLeft + nstatesRight);
			int nnodesLeft = optpro_->GetPhase(left_index)->GetTotalNodes();
			int nnodesRight = optpro_->GetPhase(right_index)->GetTotalNodes();
			int disc_left = nnodesLeft + 1;
			int disc_right = nnodesRight + 1;
			mat dLinkLeft = dLinkOut.cols(0, nstatesLeft + nparametersLeft - 1);
			mat DLink_xf_left = dLinkLeft.cols(0, nstatesLeft - 1);
			//DLink_xf_left.print("DLink_xf_left");
			mat DLink_p_left;
			if (nparametersLeft > 0)
			{
				DLink_p_left = dLinkLeft.cols(nstatesLeft, nstatesLeft + nparametersLeft - 1);
			}
			mat dLinkRight = dLinkOut.cols(dLinkLeft.n_cols, dLinkOut.n_cols - 1);
			mat DLink_x0_Right = dLinkRight.cols(0, nstatesRight - 1);
			//DLink_x0_Right.print("Dlink_x0_Right");
			mat DLink_p_Right;
			if (nparametersRight > 0)
			{
				DLink_p_Right = dLinkRight.cols(nstatesRight, nstatesRight + nparametersRight - 1);
			}

			int stateindexstart = calculateData_->phase_indices[left_index]->state[0] - 1;
			int statindexend = calculateData_->phase_indices[left_index]->state[calculateData_->phase_indices[left_index]->state.size() - 1] - 1;
			/*	dmatrix cols=Linespace(stateindexstart,statindexend,nstatesLeft);
			dmatrix indexVector=Linespace(0,nlinks-1,nlinks);
			Sjac_I.SetRow(Sjac_rowsShift,Sjac_rowsShift+nstatesLeft-1,
			indexVector+linkrow);
			Sjac_J.SetRow(Sjac_rowsShift,Sjac_rowsShift+nstatesLeft-1,
			cols);
			Sjac_V.SetRow(Sjac_rowsShift,Sjac_rowsShift+nstatesLeft-1,
			DLink_xf_left);
			Sjac_rowsShift+=DLink_xf_left.GetLength();*/
			for (int jcol = 0; jcol < DLink_xf_left.n_cols; jcol++)
			{
				for (int irow = 0; irow < DLink_xf_left.n_rows; irow++)
				{

					Sjac_I.row(Sjac_rowsShift) =
						irow + linkrow;
					Sjac_J.row(Sjac_rowsShift) =
						(jcol + 1)*(nnodesLeft)+jcol + stateindexstart;
					Sjac_rowsShift++;

				}//for irow
			}//for jcol
			//////////////////////////////////////////////////////////////////////////
			int paraindexstart = 0, paraindexend = 0;
			if (nparametersLeft > 0)
			{
				paraindexstart = calculateData_->phase_indices[left_index]->parameter[0] - 1;
				paraindexend = calculateData_->phase_indices[left_index]->parameter[calculateData_->phase_indices[left_index]->parameter.size() - 1] - 1;
				for (int jcol = 0; jcol < DLink_p_left.n_cols; jcol++)
				{
					for (int irow = 0; irow < DLink_p_left.n_rows; irow++)
					{
						Sjac_I.row(Sjac_rowsShift) =
							irow + linkrow;
						Sjac_J.row(Sjac_rowsShift) =
							jcol + paraindexstart;
						Sjac_rowsShift++;
					}//for irow
				}//for jcol
			}

			stateindexstart = calculateData_->phase_indices[right_index]->state[0] - 1;
			statindexend = calculateData_->phase_indices[right_index]->state[calculateData_->phase_indices[right_index]->state.size() - 1] - 1;
			for (int jcol = 0; jcol < DLink_x0_Right.n_cols; jcol++)
			{
				for (int irow = 0; irow < DLink_x0_Right.n_rows; irow++)
				{

					Sjac_I.row(Sjac_rowsShift) = irow + linkrow;
					Sjac_J.row(Sjac_rowsShift) = jcol*(nnodesRight + 1) + stateindexstart;
					Sjac_rowsShift++;

				}//for irow
			}//for jcol
			//////////////////////////////////////////////////////////////////////////
			paraindexstart = 0; paraindexend = 0;
			if (nparametersRight > 0)
			{
				paraindexstart = calculateData_->phase_indices[right_index]->parameter[0] - 1;
				paraindexend = calculateData_->phase_indices[right_index]->parameter[calculateData_->phase_indices[right_index]->parameter.size() - 1] - 1;
				for (int jcol = 0; jcol < DLink_p_Right.n_cols; jcol++)
				{
					for (int irow = 0; irow < DLink_p_Right.n_rows; irow++)
					{
						Sjac_I.row(Sjac_rowsShift) =
							irow + linkrow;
						Sjac_J.row(Sjac_rowsShift) =
							jcol + paraindexstart;
						Sjac_rowsShift++;

					}//for irow
				}//for jcol
			}
			calculateData_->phase_indices[right_index]->state[0];
			linkrow += DLink_x0_Right.n_rows;
		}// for ipair
	}

	void NLPWrapper::GetConsSparsity(vec& Sjac_I_return, vec& Sjac_J_return)
	{
		LP_DBG_START_FUN("NLPWrapper::GetConsSparsity()")
		static vec Sjac_I;
		static vec Sjac_J;
		if (jac_sparsity_need_refresh)
		{
			vec NL_I, NL_J, L_I, L_J, L_V, C_I, C_J;
			GetWholeSparsity(NL_I, NL_J, C_I, C_J);
			dsmatrix::Find(calculateData_->AlinearMatrix, L_I, L_J, L_V);
			int nonZerosNum = NL_I.n_elem + L_I.n_elem + C_I.n_elem;
			int linearConsRowStart = calculateData_->conbounds_min.size();//Count From 0

			Sjac_I = zeros(nonZerosNum, 1);
			Sjac_J = zeros(nonZerosNum, 1);
			int rowshift = 0;
			Sjac_I.rows(0, NL_I.n_elem - 1) = NL_I;
			Sjac_J.rows(0, NL_J.n_elem - 1) = NL_J;
			rowshift += NL_I.n_elem;
			Sjac_I.rows(NL_I.n_elem, rowshift + L_I.n_elem - 1) = L_I + linearConsRowStart;
			Sjac_J.rows(NL_J.n_elem, rowshift + L_J.n_elem - 1) = L_J;
			rowshift += L_I.n_elem;
			Sjac_I.rows(rowshift, rowshift + C_I.n_elem - 1) = C_I;
			Sjac_J.rows(rowshift, rowshift + C_J.n_elem - 1) = C_J;
			jac_sparsity_need_refresh = false;
		}
		Sjac_I_return = Sjac_I;
		Sjac_J_return = Sjac_J;
	}


	void NLPWrapper::GetPhaseScale(size_t iphase, mat& idependencies, vec& Nconstant_scale, vec& constant_scale)
	{
		LP_DBG_START_FUN("NLPWrapper::GetPhaseScale()")
		PhaseScaleShift var_scale=ocpScalor_->allScaleShift_[iphase];
		PhaseFunScale   fun_scale=ocpScalor_->allPhaseFunscale_[iphase];


		int sumnodes = optpro_->GetPhase(iphase)->GetTotalNodes();
		size_t nstates = calculateData_->SIZES_[iphase][0];
		size_t ncontrols = calculateData_->SIZES_[iphase][1];
		size_t nparameters = calculateData_->SIZES_[iphase][2];
		size_t npaths = calculateData_->SIZES_[iphase][3];
		size_t nevents = calculateData_->SIZES_[iphase][4];
		//get derivatives of DAE function
		mat dDaeScale = ones(nstates, nstates + ncontrols + 1 + nparameters);//allocate mem
		dDaeScale.each_col() %= fun_scale.ode_scale;//scale fun
		vec allvarscale = ones(nstates + ncontrols + 1 + nparameters);
		allvarscale.subvec(0, nstates - 1) = var_scale.statescale;//x
		if (ncontrols > 0)allvarscale.subvec(nstates, nstates + ncontrols - 1) = var_scale.controlscale;//u
		allvarscale(nstates + ncontrols) = var_scale.tscale;//t
		if (nparameters > 0)allvarscale.subvec(nstates + ncontrols + 1, nstates + ncontrols + 1 + nparameters - 1) = var_scale.parameterscale;//p
		dDaeScale.each_row() /= trans(allvarscale);//scale var
		mat dPathScale;
		if (npaths > 0)
		{
			dPathScale = ones(npaths, nstates + ncontrols + 1 + nparameters);//allocate men
			dPathScale.each_col() %= fun_scale.path_scale;//scale fun
			dPathScale.each_row() /= trans(allvarscale);//scale var
		}

		//mat dDae_time = zeros(sumnodes, nstates);
		mat dDae_time = zeros(sumnodes, nstates);
		std::vector<mat> dDae_state(nstates);
		std::vector<mat> dDae_control(ncontrols);
		std::vector<mat> dDae_parameter(nparameters);

		mat dPath_time;
		if (npaths > 0) dPath_time = zeros(sumnodes, npaths);
		std::vector<mat> dPath_state(nstates);
		std::vector<mat> dPath_control(ncontrols);
		std::vector<mat> dPath_parameter(nparameters);
		//Compute the Derivative  Respect to State
		for (int istate = 0; istate < nstates; istate++)
		{
			//dDae_state[istate] = reshape(dDaeScale.col(istate), sumnodes, ndiffeqs);
			//dDae_state[istate] = ones(sumnodes, nstates);
			dDae_state[istate] = repmat(trans(dDaeScale.col(istate)), sumnodes, 1);
			//dDae_state[istate].print("dDae_state[istate]");
			if (npaths > 0)
			{
				dPath_state[istate] = repmat(trans(dPathScale.col(istate)), sumnodes, 1);
				//dPath_state[istate] = reshape(dPathScale.col(istate), sumnodes, npaths);
				//dPath_state[istate].print("dpath_state[istate]");
			}
		}
		//Compute the Derivative  Respect to Control
		for (int icontrol = 0; icontrol < ncontrols; icontrol++)
		{
			dDae_control[icontrol] = repmat(trans(dDaeScale.col(nstates + icontrol)), sumnodes, 1);
			//dDae_control[icontrol].print("dDae_control[icontrol]");
			if (npaths > 0)
			{
				dPath_control[icontrol] = repmat(trans(dPathScale.col(nstates + icontrol)), sumnodes, 1);
				//dPath_control[icontrol].print("dPath_control[icontrol]");
			}
		}
		//Compute the Derivative  Respect to time
		// 
		dDae_time = repmat(trans(dDaeScale.col(nstates + ncontrols)), sumnodes, 1);
		//dDae_time.print("dDae_time");
		if (npaths > 0)
		{
			dPath_time = repmat(trans(dPathScale.col(nstates + ncontrols)), sumnodes, npaths);
			//dPath_time.print("dPath_time");
		}

		if (nparameters > 0)
		{
			for (int iparameter = 0; iparameter < nparameters; iparameter++)
			{
				dDae_parameter[iparameter] = repmat(trans(dDaeScale.col(iparameter + nstates + ncontrols)), sumnodes, 1);
				if (npaths > 0)
				{
					dPath_parameter[iparameter] = repmat(trans(dPathScale.col(iparameter + nstates + ncontrols)), sumnodes, 1);
				}
			}
		}//end if nparameters
		//Get derivatives of event function

		vec dEvent_t0;
		vec dEvent_tf;
		if (nevents > 0)
		{
			dEvent_t0 = zeros(nevents, 1);
			dEvent_tf = zeros(nevents, 1);
		}

		std::vector<vec> dEvent_x0(nstates);
		std::vector<vec> dEvent_xf(nstates);
		std::vector<vec> dEvent_parameter(nparameters);
		if (nevents > 0)
		{
			mat dEventScale = ones(nevents, nstates + 1 + nstates + 1 + nparameters);
			dEventScale.each_col() %= fun_scale.event_scale;//scale fun
			vec xandtscale = ones(nstates * 2 + 2 + nparameters);
			xandtscale.subvec(0, nstates - 1) = var_scale.statescale;
			xandtscale(nstates) = var_scale.t0scale;
			xandtscale.subvec(nstates + 1, nstates + nstates) = var_scale.statescale;
			xandtscale(nstates + nstates + 1) = var_scale.tfscale;
			if (nparameters > 0)
			{
				xandtscale.tail_rows(nparameters) = var_scale.parameterscale;
			}
			dEventScale.each_row() /= trans(xandtscale);//scale var
			//derive_->DerivEvent(mySolEvents, dEventOut);
			//dEventOut.print("dEventOut");
			//Get Derivatives of Event with Respect to Initial Time
			dEvent_t0 = dEventScale.col(nstates);
			//dEvent_t0.print("devent_t0");
			//Get Derivatives of Event with Respect to Initial Time
			dEvent_tf = dEventScale.col(2 * nstates + 1);
			//dEvent_tf.print("devent_tf");
			//Get Derivatives of Event with Respect to Initial and Terminal State
			for (int istate = 0; istate < nstates; istate++)
			{
				dEvent_x0[istate] = dEventScale.col(istate);
				//dEvent_x0[istate].print("devent_x0");
				dEvent_xf[istate] = dEventScale.col(istate + nstates + 1);
				//dEvent_xf[istate].print("devent_xf");
			}
			//Get Derivatives of Event with Respect to Parameters
			for (int iparameter = 0; iparameter < nparameters; iparameter++)
			{
				dEvent_parameter[iparameter] = dEventScale.col(iparameter + 2 * (nstates + 1));
			}
		}//if event>0
		//////////////////////////////////////////////////////////////////////////
		//Get Sparsity And Insert Derives of constraints
		int ndiffeqs = nstates;
		int disc_pts = sumnodes + 1;
		int numcons = ndiffeqs*(disc_pts - 1) + npaths*sumnodes + nevents;
		//find number of non-zero elements in dependencies (always including diagonal)
		for (int i = 0; i < ndiffeqs; i++)
		{
			idependencies(i, i) = 1.0;
		}
		int nDependancies = nnz(idependencies);
		//allocate memory for Jacobian scale
		Nconstant_scale = zeros(nDependancies*sumnodes + 2 * (ndiffeqs + npaths)*sumnodes +
			(ndiffeqs + npaths)*nparameters*sumnodes + nevents*(2 * nstates + nparameters + 2), 1);
		//find non-zeros in off diagonal D matrix
		vec DiffMatOffDiag_I, DiffMatOffDiag_J, DiffMatOffDiag_V;
		dsmatrix::Find(calculateData_->PS[iphase]->Doffdiag, DiffMatOffDiag_I, DiffMatOffDiag_J, DiffMatOffDiag_V);
		lp_index nonZerosDiffMat = DiffMatOffDiag_I.n_elem;
		//allocate memory for constant derivatives scales
		constant_scale = zeros(ndiffeqs*nonZerosDiffMat, 1);
		//differential equations
		lp_index rowstart, colstart;
		lp_index Sjac_rowShift = 0, Sconstant_rowShift = 0;
		vec indexvector = linspace(0, sumnodes - 1, sumnodes);
		for (int i = 0; i < ndiffeqs; i++)
		{
			rowstart = i*sumnodes;
			for (int j = 0; j < nstates; j++)//inset df/dx
			{
				colstart = j*disc_pts;
				if (i == j)
				{
					// DdiagQuadBlock
					// 					int Ddiag_ele_num = calculateData_->PS[iphase]->Diag.GetnRows();
					// 					double*dDiagValues = new double[Ddiag_ele_num];
					// 					double*PsDiag = calculateData_->PS[iphase]->Diag.GetSArrayConst()->GetDataPtr();
					// 					for (int ii = 0; ii < Ddiag_ele_num; ii++)
					// 					{
					// 						dDiagValues[ii] = PsDiag[ii];
					// 					}
					//vec Ddiag(dDiagValues, Ddiag_ele_num, false);
					//vec ret = Ddiag - dDae_state[j].col(i)*(tf - t0) / 2.0;
					Nconstant_scale.subvec(Sjac_rowShift, Sjac_rowShift + sumnodes - 1) = dDae_state[j].col(i);
					Sjac_rowShift += sumnodes;
					//DiffMatOffDiag
					//dmatrix Sconstant_rows=diffmatrixindex_vector+i*nonZerosDiffMat;
					Sconstant_rowShift = i*nonZerosDiffMat;
					constant_scale.subvec(Sconstant_rowShift, Sconstant_rowShift + nonZerosDiffMat - 1).
						fill(fun_scale.ode_scale(i)*var_scale.statescale(j));
				}
				else
				{
					if (idependencies(i, j) == 1.0)
					{
						//OffDiag_State
						Nconstant_scale.rows(Sjac_rowShift, Sjac_rowShift + sumnodes - 1) = dDae_state[j].col(i);
						Sjac_rowShift += sumnodes;
					}
				}
			}//for j
			//insert df/du
			lp_index colshift = nstates*disc_pts;
			for (int j = 0; j < ncontrols; j++)
			{
				colstart = colshift + j*sumnodes;
				if (idependencies(i, j + nstates) == 1.0)
				{
					//Diffeq_Block_Control
					Nconstant_scale.rows(Sjac_rowShift, Sjac_rowShift + sumnodes - 1) = dDae_control[j].col(i);
					Sjac_rowShift += sumnodes;
				}
			}//for j
			//////////////////////////////////////////////////////////////////////////
			//insert df/dt0
			colshift += ncontrols*sumnodes;
			colstart = colshift;
			Nconstant_scale.rows(Sjac_rowShift, Sjac_rowShift + sumnodes - 1) = dDae_time.col(i);
			Sjac_rowShift += sumnodes;
			colshift++;
			// insert df/dtf
			colstart = colshift;
			Nconstant_scale.rows(Sjac_rowShift, Sjac_rowShift + sumnodes - 1) = dDae_time.col(i);
			Sjac_rowShift += sumnodes;
			colshift++;
			//insert df/dpara
			for (int j = 0; j < nparameters; j++)
			{
				colstart = colshift + j;
				Nconstant_scale.rows(Sjac_rowShift, Sjac_rowShift + sumnodes - 1) = dDae_parameter[j].col(i);
				Sjac_rowShift += sumnodes;
			}//for j
		}//for i
		//path contraints
		lp_index rowshift = ndiffeqs*(sumnodes);
		for (size_t i = 0; i < npaths; i++)
		{
			rowstart = rowshift + i*sumnodes;
			for (int j = 0; j < nstates; j++)
			{
				colstart = j*disc_pts;
				if (idependencies(i + ndiffeqs, j) == 1.0)
				{
					//Insert dc/dx
					Nconstant_scale.rows(Sjac_rowShift, Sjac_rowShift + sumnodes - 1) = dPath_state[j].col(i);
					Sjac_rowShift += sumnodes;
				}
			}//for j
			lp_index colshift = nstates*disc_pts;
			for (size_t j = 0; j < ncontrols; j++)
			{
				colstart = colshift + j*sumnodes;
				if (idependencies(i + ndiffeqs, j + nstates))
				{
					//insert dc/du
					Nconstant_scale.rows(Sjac_rowShift, Sjac_rowShift + sumnodes - 1) = dPath_control[j].col(i);
					Sjac_rowShift += sumnodes;
				}
			}//for j
			//////////////////////////////////////////////////////////////////////////
			//insert dc/dt0
			colshift += ncontrols*sumnodes;
			colstart = colshift;
			Nconstant_scale.rows(Sjac_rowShift, Sjac_rowShift + sumnodes - 1) = dPath_time.col(i);
			Sjac_rowShift += sumnodes;
			colshift++;
			// insert dc/dtf
			colstart = colshift;
			Nconstant_scale.rows(Sjac_rowShift, Sjac_rowShift + sumnodes - 1) = dPath_time.col(i);
			Sjac_rowShift += sumnodes;
			colshift++;
			//insert dc/dpara
			for (int j = 0; j < nparameters; j++)
			{
				colstart = colshift + j;
				Nconstant_scale.rows(Sjac_rowShift, Sjac_rowShift + sumnodes - 1) = dDae_parameter[j].col(i);
				Sjac_rowShift += sumnodes;
			}//for j
		}//for i
		rowshift += npaths*sumnodes;
		//Event Constraints
		rowshift = ndiffeqs*(sumnodes)+npaths*sumnodes;

		for (size_t i = 0; i < nevents; i++)
		{
			lp_index row = rowshift + i;
			//devent/dx
			for (size_t j = 0; j < nstates; j++)
			{
				Nconstant_scale(Sjac_rowShift) = dEvent_x0[j](i);
				Sjac_rowShift++;
				Nconstant_scale(Sjac_rowShift) = dEvent_xf[j](i);
				Sjac_rowShift++;
			}
			//insert devent/dt0
			lp_index cols = nstates*(sumnodes + 1) + ncontrols*sumnodes;
			Nconstant_scale(Sjac_rowShift) = dEvent_t0(i, 0);
			Sjac_rowShift++;
			//insert devent/df
			cols++;
			Nconstant_scale(Sjac_rowShift) = dEvent_tf(i, 0);
			Sjac_rowShift++;
			//insert dEvent/dparam
			for (size_t j = 0; j < nparameters; j++)
			{
				cols++;
				Nconstant_scale(Sjac_rowShift) = dEvent_parameter[j](i);
				Sjac_rowShift++;
			}
		}//for ievent

	}

	void NLPWrapper::GetWholeScale(vec& Nconstant_scale, vec& constant_scale)
	{
		LP_DBG_START_FUN("NLPWrapper::GetWholeScale()")
		vector<mat> dependencies(optpro_->GetPhaseNum());

		for (int itor = 0; itor < calculateData_->numphases_; itor++)
		{
			int i_nstate = optpro_->GetPhase(itor)->GetstateMin().size();
			int i_ncontrol = optpro_->GetPhase(itor)->GetcontrolMin().size();
			int i_npaths = optpro_->GetPhase(itor)->GetpathMin().size();
			dependencies[itor] = zeros(i_nstate + i_npaths, i_nstate + i_ncontrol);
			dependencies[itor].fill(1.0);
		}
		std::vector<PhaseScaleShift> allphasevar_scale=ocpScalor_->allScaleShift_;
		std::vector<PhaseFunScale>   allphasefun_scale=ocpScalor_->allPhaseFunscale_;
		std::vector<vec> alllinkfun_scale=ocpScalor_->allLinkFunscale_;
		//find total non-zero elements in Jacobian and constant derivatives
		lp_index  nonZerosSjac = 0;
		lp_index  nonZerosSconstant = 0;
		lp_index  nonZerosSnonconstant = 0;

		for (int i = 0; i < optpro_->GetPhaseNum(); i++)
		{
			PhaseScaleShift var_scale = allphasevar_scale[i];
			PhaseFunScale   fun_scale = allphasefun_scale[i];


			size_t i_totnodes = optpro_->GetPhase(i)->GetTotalNodes();
			size_t i_nstate = calculateData_->SIZES_[i][0];
			size_t ncontrols = calculateData_->SIZES_[i][1];
			size_t i_nparameters = calculateData_->SIZES_[i][2];
			size_t i_npaths = calculateData_->SIZES_[i][3];
			size_t i_nevents = calculateData_->SIZES_[i][4];
			// find number of non-zero elements in dependencies (always including diagonal)


			//////////////////////////////////////////////////////////////////////////
			//only when use Analytical Derives 
			//dependencies[i].fill(1.0);
			////////////////////////////////////////////////////////////////////////
			mat temDependencies = dependencies[i];
			for (size_t j = 0; j < i_nstate; j++)
			{
				temDependencies(j, j) = 1.0;
			}

			int nDependancies = nnz(temDependencies);

			nonZerosSjac = nonZerosSjac + nDependancies*i_totnodes + 2 * (i_nstate + i_npaths)*i_totnodes +
				(i_nstate + i_npaths)*i_nparameters*i_totnodes + i_nevents*(2 * i_nstate + i_nparameters + 2);
			vec DiffMatOffDiag_I, DiffMatOffDiag_J, DiffMatOffDiag_V;
			dsmatrix::Find(calculateData_->PS[i]->Doffdiag, DiffMatOffDiag_I, DiffMatOffDiag_J, DiffMatOffDiag_V);
			int nonZerosDiffMat = DiffMatOffDiag_I.n_elem;
			nonZerosSconstant = nonZerosSconstant + nonZerosDiffMat*i_nstate;
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
			nonZerosSjac += numlinks*(nstates_left + nparameter_left + nstates_right + nparameter_right);
		}//for ipair
		//allocate Memory for Jacobian
		Nconstant_scale = zeros(nonZerosSjac, 1);
		constant_scale = zeros(nonZerosSconstant, 1);
		lp_index Sjac_rowsShift = 0, Sconstant_rowsShift = 0;
		lp_index rowshift = 0;
		lp_index colshift = 0;
		//////////////////////////////////////////////////////////////////////////
		//Get Jac in each Phase
		for (int i = 0; i < optpro_->GetPhaseNum(); i++)
		{
			size_t i_nstate = calculateData_->SIZES_[i][0];
			size_t i_ncontrol = calculateData_->SIZES_[i][1];
			size_t i_nparameters = calculateData_->SIZES_[i][2];
			size_t i_npaths = calculateData_->SIZES_[i][3];
			size_t i_nevents = calculateData_->SIZES_[i][4];
			size_t i_totnodes = optpro_->GetPhase(i)->GetTotalNodes();
			vec  Nconstant_phase_scale, constant_phase_scale;
			GetPhaseScale(i, dependencies[i], Nconstant_phase_scale, constant_phase_scale);
			Nconstant_scale.rows(Sjac_rowsShift, Sjac_rowsShift + Nconstant_phase_scale.n_elem - 1) =
				Nconstant_phase_scale;
			Sjac_rowsShift += Nconstant_phase_scale.n_elem;
			constant_scale.rows(Sconstant_rowsShift, Sconstant_rowsShift + constant_phase_scale.n_elem - 1) =
				constant_phase_scale;
			Sconstant_rowsShift += constant_phase_scale.n_elem;
			size_t	numcons = i_nstate*i_totnodes + i_npaths*i_totnodes + i_nevents;
			size_t numvars = i_nstate*(i_totnodes + 1) + i_ncontrol*i_totnodes + i_nparameters + 2;
			rowshift += numcons;
			colshift += numvars;
		}
		//derive of linkage and it's Sparsity Pattern
		//int sjacRowsStart=calculateData_->conbounds_min.size()-optpro_->GetPhaseNum()-1;//count form zero
		int linkrow = rowshift;//index start from zero
		for (int ipair = 0; ipair < optpro_->GetLinkageNum(); ipair++)
		{
			int nlinks = optpro_->GetLinkage(ipair)->GetLinkageMin().size();
			int left_index = optpro_->GetLinkage(ipair)->LeftPhase();
			int right_index = optpro_->GetLinkage(ipair)->RightPhase();

			//////////////////////////////////////////////////////////////////////////


			int nstatesLeft, ncontrolsLeft, nparametersLeft, npathLeft, neventLeft;
			int nstatesRight, ncontrolsRight, nparametersRight, npathRight, neventRight;
			optpro_->GetPhase(left_index)->get_optimal_info(
				nstatesLeft, ncontrolsLeft, nparametersLeft, npathLeft, neventLeft);
			optpro_->GetPhase(right_index)->get_optimal_info(
				nstatesRight, ncontrolsRight, nparametersRight, npathRight, neventRight);

			//Get Linkage Derivatives
			mat dLinkScale = zeros(nlinks, nstatesLeft + nparametersLeft + nstatesRight + nparametersRight);
			dLinkScale.each_col() %= alllinkfun_scale[ipair];//scale fun
			vec xlxrplprscall = zeros(nstatesLeft + nparametersLeft + nstatesRight + nparametersRight);
			xlxrplprscall.subvec(0, nstatesLeft - 1) = allphasevar_scale[left_index].statescale;
			if (nparametersLeft > 0)
			{
				xlxrplprscall.subvec(nstatesLeft, nstatesLeft + nparametersLeft - 1) = allphasevar_scale[left_index].parameterscale;
			}

			xlxrplprscall.subvec(nstatesLeft + nparametersLeft, nstatesLeft + nparametersLeft + nstatesRight - 1)
				= allphasevar_scale[right_index].statescale;
			if (nparametersRight > 0)
			{
				xlxrplprscall.tail_rows(nparametersRight) = allphasevar_scale[right_index].parameterscale;
			}
			dLinkScale.each_row() /= trans(xlxrplprscall);//scale var
			mat dLinkLeft = dLinkScale.cols(0, nstatesLeft + nparametersLeft - 1);
			mat DLink_xf_left = dLinkLeft.cols(0, nstatesLeft - 1);
			//DLink_xf_left.print("DLink_xf_left");
			mat DLink_p_left;
			if (nparametersLeft > 0)
			{
				DLink_p_left = dLinkLeft.cols(nstatesLeft, nstatesLeft + nparametersLeft - 1);
			}
			mat dLinkRight = dLinkScale.cols(dLinkLeft.n_cols, dLinkScale.n_cols - 1);
			mat DLink_x0_Right = dLinkRight.cols(0, nstatesRight - 1);
			//DLink_x0_Right.print("Dlink_x0_Right");
			mat DLink_p_Right;
			if (nparametersRight > 0)
			{
				DLink_p_Right = dLinkRight.cols(nstatesRight, nstatesRight + nparametersRight - 1);
			}
			for (int jcol = 0; jcol < DLink_xf_left.n_cols; jcol++)
			{
				for (int irow = 0; irow < DLink_xf_left.n_rows; irow++)
				{
					double retvalue = DLink_xf_left(irow, jcol);
					Nconstant_scale.row(Sjac_rowsShift) =
						retvalue;
					Sjac_rowsShift++;

				}//for irow
			}//for jcol
			//////////////////////////////////////////////////////////////////////////

			if (nparametersLeft > 0)
			{
				for (int jcol = 0; jcol < DLink_p_left.n_cols; jcol++)
				{
					for (int irow = 0; irow < DLink_p_left.n_rows; irow++)
					{
						double retvalue = DLink_p_left(irow, jcol);
						Nconstant_scale.row(Sjac_rowsShift) =
							retvalue;
						Sjac_rowsShift++;
					}//for irow
				}//for jcol
			}

			for (int jcol = 0; jcol < DLink_x0_Right.n_cols; jcol++)
			{
				for (int irow = 0; irow < DLink_x0_Right.n_rows; irow++)
				{
					double retvalue = DLink_x0_Right(irow, jcol);
					Nconstant_scale.row(Sjac_rowsShift) = retvalue;
					Sjac_rowsShift++;

				}//for irow
			}//for jcol
			//////////////////////////////////////////////////////////////////////////
			if (nparametersRight > 0)
			{
				for (int jcol = 0; jcol < DLink_p_Right.n_cols; jcol++)
				{
					for (int irow = 0; irow < DLink_p_Right.n_rows; irow++)
					{
						double retvalue = DLink_p_Right(irow, jcol);
						Nconstant_scale.row(Sjac_rowsShift) =
							retvalue;
						Sjac_rowsShift++;

					}//for irow
				}//for jcol
			}
			linkrow += DLink_x0_Right.n_rows;
		}// for ipair
	}

	void NLPWrapper::GetConsScale(vec& Sjac_Scale)
	{
		LP_DBG_START_FUN("NLPWrapper::GetConsScale()")
		std::vector<PhaseScaleShift> allphasevar_scale=ocpScalor_->allScaleShift_;
		vec  NL_V, L_I, L_J, L_V, C_V;
		GetWholeScale(NL_V, C_V);
		//Get Bounds on Linear Constraints
		//int numvars = calcu_data->varbounds_min.size();

		int alinrowshift = 0;

		vec Alinear_V = zeros(2 * (calculateData_->numphases_ + calculateData_->numlinkpairs_), 1);
		//Part 1:  Monotonicity of Independent Variable
		for (int i = 0; i < calculateData_->numphases_; i++){
			PhaseScaleShift var_scale = allphasevar_scale[i];
			shared_ptr<Phase> currentphase(optpro_->GetPhase(i));
			int nstate = currentphase->GetstateMin().size();
			int ncontrol = currentphase->GetcontrolMin().size();
			size_t nodes = optpro_->GetPhase(i)->GetTotalNodes();
			int ishift = 0;
			if (i != 0){
				ishift = calculateData_->variable_indices[i - 1][calculateData_->variable_indices[i - 1].size() - 1];
			}

			int t0_index = ishift + ((nodes + 1)*nstate)
				+ nodes * ncontrol + 1;
			int tf_index = t0_index + 1;

			Alinear_V(alinrowshift) = var_scale.t0scale;
			alinrowshift++;
			Alinear_V(alinrowshift) = var_scale.tfscale;
			alinrowshift++;

		}//for int i

		int istart = calculateData_->numphases_;
		//Part 2:  Linkage of Time Across Phases
		for (int i = 0; i < calculateData_->numlinkpairs_; i++){

			shared_ptr<Linkage> currentlinkage(optpro_->GetLinkage(i));

			int left_phase = currentlinkage->LeftPhase();
			int right_phase = currentlinkage->RightPhase();
			PhaseScaleShift left_var_scale = allphasevar_scale[left_phase];
			PhaseScaleShift right_var_scale = allphasevar_scale[right_phase];
			int nparameters_left = calculateData_->SIZES_[left_phase][2];
			int nparameters_right = calculateData_->SIZES_[right_phase][2];

			int tf_index_left = calculateData_->variable_indices[left_phase][calculateData_->variable_indices[left_phase].size() - nparameters_left - 1];
			int t0_index_right = calculateData_->variable_indices[right_phase][calculateData_->variable_indices[right_phase].size() - nparameters_right - 2];
			
			Alinear_V(alinrowshift) = left_var_scale.tfscale;
			alinrowshift++;
			
			Alinear_V(alinrowshift) = right_var_scale.t0scale;
			alinrowshift++;
		}


		//////////////////////////////////////////////////////////////////////////
		L_V = 1.0/Alinear_V;

		//////////////////////////////////////////////////////////////////////////
		int nonZerosNum = NL_V.n_elem + L_V.n_elem + C_V.n_elem;
		int linearConsRowStart = calculateData_->conbounds_min.size();//Count From 0
		Sjac_Scale = zeros(nonZerosNum, 1);
		int rowshift = 0;
		Sjac_Scale.rows(0, NL_V.n_elem - 1) = NL_V;
		rowshift += NL_V.n_elem;
		Sjac_Scale.rows(NL_V.n_elem, rowshift + L_V.n_elem - 1).fill(1.0);
		rowshift += L_V.n_elem;
		Sjac_Scale.rows(rowshift, rowshift + C_V.n_elem - 1) = C_V;
	}

	void NLPWrapper::GetObjGradScale(vec& grad_scale)
	{
		LP_DBG_START_FUN("NLPWrapper::GetObjGradScale()")
		vec allphaseobjscales=calculateData_->objScale;
		std::vector<PhaseScaleShift> allphasevar_scale=ocpScalor_->allScaleShift_;
		std::vector<PhaseFunScale>   allphasefun_scale=ocpScalor_->allPhaseFunscale_;
		//////////////////////////////////////////////////////////////////////////
		grad_scale = zeros(calculateData_->varbounds_max.size(), 1);
		int grad_shift = 0;
		for (int iphase = 0; iphase < optpro_->GetPhaseNum(); iphase++)
		{
			PhaseScaleShift var_scale = allphasevar_scale[iphase];
			PhaseFunScale   fun_scale = allphasefun_scale[iphase];

			int sumnodes = optpro_->GetPhase(iphase)->GetTotalNodes();
			int nstates = optpro_->GetPhase(iphase)->GetstateMin().size();
			int ncontrols = optpro_->GetPhase(iphase)->GetcontrolMin().size();
			int nparameters = optpro_->GetPhase(iphase)->GetparameterMin().size();
			int npaths = optpro_->GetPhase(iphase)->GetpathMin().size();
			int nevents = optpro_->GetPhase(iphase)->GeteventMin().size();
			int stateindex_start = calculateData_->phase_indices[iphase]->state[0] - 1;
			int stateindex_end = calculateData_->phase_indices[iphase]->state[calculateData_->phase_indices[iphase]->state.size() - 1] - 1;
			int controlindex_start = calculateData_->phase_indices[iphase]->control[0] - 1;
			int controlindex_end = calculateData_->phase_indices[iphase]->control[calculateData_->phase_indices[iphase]->control.size() - 1] - 1;

			int t0_index = calculateData_->phase_indices[iphase]->time[0] - 1;
			int tf_index = calculateData_->phase_indices[iphase]->time[1] - 1;

			rowvec dmayerscale = ones<rowvec>(nstates * 2 + 2 + nparameters);
			dmayerscale *= allphaseobjscales(iphase);//scale fun
			vec xandtscale = ones(nstates * 2 + 2 + nparameters);
			xandtscale.subvec(0, nstates - 1) = var_scale.statescale;
			xandtscale(nstates) = var_scale.t0scale;
			xandtscale.subvec(nstates + 1, nstates + nstates) = var_scale.statescale;
			xandtscale(nstates + nstates + 1) = var_scale.tfscale;
			if (nparameters > 0)
			{
				xandtscale.tail_rows(nparameters) = var_scale.parameterscale;
			}
			dmayerscale /= trans(xandtscale);//scale var



			mat dLagrangeScale = ones(sumnodes,nstates + ncontrols + 1 + nparameters);
			dLagrangeScale *= allphaseobjscales(iphase);//scale fun
			vec allvarscale = ones(nstates + ncontrols + 1 + nparameters);
			allvarscale.subvec(0, nstates - 1) = var_scale.statescale;//x
			if (ncontrols > 0)allvarscale.subvec(nstates, nstates + ncontrols - 1) = var_scale.controlscale;//u
			allvarscale(nstates + ncontrols) = var_scale.tscale;//t
			if (nparameters > 0)allvarscale.subvec(nstates + ncontrols + 1, nstates + ncontrols + 1 + nparameters - 1) = var_scale.parameterscale;//p
			dLagrangeScale.each_row() /= trans(allvarscale);//scale var


			double dMayer_t0 = 0.0;
			double dMayer_tf = 0.0;
			rowvec dMayer_x0 = zeros(1, nstates);
			rowvec dMayer_xf = zeros(1, nstates);
			rowvec dMayer_para;
			if (nparameters > 0)
			{
				dMayer_para = zeros(1, nparameters);
			}
			mat dLagrange_state = zeros(sumnodes, nstates);
			vec dLagrange_time = zeros(sumnodes, 1);
			mat dLagrange_control, dLagrange_para;
			if (ncontrols > 0)dLagrange_control = zeros(sumnodes, ncontrols);
			if (nparameters > 0)dLagrange_para = zeros(sumnodes, nparameters);
			//Get Derivatives of mayer cost 
			if (dmayerscale.n_elem > 0)
			{

				dMayer_t0 = dmayerscale(0, nstates);
				dMayer_tf = dmayerscale(0, 2 * nstates + 1);
				for (int istate = 0; istate < nstates; istate++)
				{
					dMayer_x0(istate) = dmayerscale(istate);
					dMayer_xf(istate) = dmayerscale(istate + nstates + 1);
				}
				for (int ipara = 0; ipara < nparameters; ipara++)
				{
					dMayer_para(ipara) = dmayerscale(ipara + 2 * (nstates + 1));
				}
			}//if dmayerOut.length>0
			//Get Derivatives of Lagrange cost 
			if (dLagrangeScale.n_elem > 0)
			{
				dLagrange_time = dLagrangeScale.col(dLagrangeScale.n_cols - 1);
				//dLagrange_time.print("L_time");
				for (int istate = 0; istate < nstates; istate++)
				{
					dLagrange_state.col(istate) = dLagrangeScale.col(istate);
				}
				for (int icontrol = 0; icontrol < ncontrols; icontrol++)
				{
					dLagrange_control.col(icontrol) = dLagrangeScale.col(icontrol + nstates);
				}
				for (int ipara = 0; ipara < nparameters; ipara++)
				{
					dLagrange_para.col(ipara) = dLagrangeScale.col(ipara + nstates + ncontrols);
				}
			}//if dlagrangeout.length>0

			//////////////////////////////////////////////////////////////////////////
			//insert Cost derivatives
			int colshift = 0;
			rowvec Jscale = zeros(1, nstates*(sumnodes + 1) + ncontrols*sumnodes + nparameters + 2);
			for (int j = 0; j < nstates; j++)
			{
				int col0 = sumnodes*(j)+j;
				int colf = sumnodes*(j + 1) + j;
				vec index_vector = linspace(0, sumnodes - 1, sumnodes);
				Jscale(0, col0) = dMayer_x0(0, j);
				Jscale.subvec(col0, colf - 1) = trans(dLagrange_state.col(j));
				Jscale(colf) = dMayer_xf(j);
			}
			//dCost/du
			colshift += nstates*(sumnodes + 1);
			for (int j = 0; j < ncontrols; j++)
			{
				int colstart = colshift + j*sumnodes;
				Jscale.subvec(colstart, colstart + sumnodes - 1) = trans(dLagrange_control.col(j));
			}
			colshift += ncontrols*sumnodes;
			int t0_col = colshift;

			//dCost/dt0
			//ret.print("ret");
			mat ret2 = calculateData_->PS[iphase]->Weights;
			ret2.fill(1.0);
			ret2 = trans(ret2)*diagmat(dLagrange_time);
			//ret2.print("ret2");
			//dLagrange_time.diag().print("diag");
			mat ret3 = calculateData_->PS[iphase]->Points;
			ret3.fill(1.0);
			ret3 = ret2*ret3 + dMayer_t0;
			Jscale(t0_col) = ret3(0);
			int tf_col = colshift + 1;
			//dCost/dtf

			ret2 = calculateData_->PS[iphase]->Weights;
			ret2.fill(1.0);
			ret2 = trans(ret2)*diagmat(dLagrange_time);
			ret3 = calculateData_->PS[iphase]->Points;
			ret3.fill(1.0);
			ret3 *= ret2;
			Jscale(tf_col) = dMayer_tf + ret3(0);
			colshift = colshift + 1;
			//dCost/dpara
			for (int j = 0; j < nparameters; j++)
			{
				int paracol = colshift + j;
				vec ret = ones(calculateData_->PS[iphase]->Weights.n_elem);
				ret *= dLagrange_para.col(j);
				Jscale(paracol) = dMayer_para(j) + as_scalar(ret);
			}
			grad_scale.subvec(grad_shift, grad_shift + Jscale.n_elem - 1) = trans(Jscale);
			grad_shift += Jscale.n_elem;
			//Jcost.print("Jcost");
		}//for iphase
	}

	void NLPWrapper::GetScales()
	{
		LP_DBG_START_FUN("NLPWrapper::GetScales()")
		ocpScalor_->CalculateVarScale();
		ocpScalor_->CalculateFunScaleFromRand();
		ocpScalor_->GetOCPScaleAndShift();
		//scale guess
		calculateData_->scaled_nlpGuessVector = calculateData_->nlpGuessVector%calculateData_->varscale + calculateData_->varshift;
		//scale var bounds
		vec unscaled_var_min = calculateData_->varbounds_min;
		vec unscaled_var_max = calculateData_->varbounds_max;
		calculateData_->scaled_varbounds_min = unscaled_var_min%calculateData_->varscale + calculateData_->varshift;
		calculateData_->scaled_varbounds_max = unscaled_var_max%calculateData_->varscale + calculateData_->varshift;
		//scale fun bounds
		vec unscaled_fun_min = calculateData_->conbounds_min;
		vec unscaled_fun_max = calculateData_->conbounds_max;
		calculateData_->scaled_conbounds_min = unscaled_fun_min%calculateData_->funscale.head_rows(unscaled_fun_min.n_elem);//exclude Alinear constrains
		calculateData_->scaled_conbounds_max = unscaled_fun_max%calculateData_->funscale.head_rows(unscaled_fun_min.n_elem);//exclude Alinear constrains

		//scale jacobi
		GetConsScale(conscale_);
		//scale grad
		GetObjGradScale(objscale_);
	}

	void NLPWrapper::GetHessainValue(const vec& x, const double sigma, const vec&lambda, vec & H_V)
	{
		LP_DBG_START_FUN("NLPWrapper::GetHessainValue()")
		hessian_->GetHessian(sigma, x, lambda, H_V);
	}

	void NLPWrapper::GetHessainSparsity(vec& H_I_return, vec& H_J_return)
	{
		LP_DBG_START_FUN("NLPWrapper::GetHessainSparsity()")
		if (!hessian_.get())
		{
			return;
		}
		static vec H_I;
		static vec H_J;

		if (hessain_sparsity_need_refresh)
		{
			hessian_->GetHessianSparsity(H_I, H_J);
			hessain_sparsity_need_refresh = false;
		}
		H_I_return = H_I;
		H_J_return = H_J;
	}

}//Lpopc