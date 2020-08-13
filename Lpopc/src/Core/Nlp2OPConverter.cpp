// Copyright (C) 2014-2015 Xue Zhichen, Wang Yujie,Wang Na 
// All Rights Reserved.
// This file is a part of LPOPC , published under the Eclipse Public License.
// Author:Xue Zhichen 7/12 2015   19:11
// Email:eddy_lpopc@163.com
#include "Nlp2OPConverter.h"
#include "LpGuessChecker.h"
#include "LpOptimalProblem.hpp"
#include <cstdio>
#include<strstream>
namespace Lpopc
{
	void Nlp2OpConverter::Nlp2OpControl(shared_ptr<FunctionWrapper>& Funs_, shared_ptr<LpCalculateData>& Data_, shared_ptr<OptimalProblem>& optpro_)
	{

		std::vector<shared_ptr<SolutionData>> resultData(Data_->numphases_);
		vec x, lamba;
		if (Data_->autoscale)
		{
			Data_->scaled_nlpreturn_x = (Data_->nlpreturn_x-Data_->varshift) / Data_->varscale;
			x = Data_->scaled_nlpreturn_x;
			Data_->scaled_nlpreturn_lambda = (Data_->nlpreturn_lambda) / Data_->funscale;
			lamba = Data_->scaled_nlpreturn_lambda;
		}
		else
		{
			x = Data_->nlpreturn_x;
			lamba = Data_->nlpreturn_lambda;
		}
		//Data_->scaled_nlpreturn_x.save("scaledx", raw_ascii);
		size_t varshift = 0, conshift = 0;
		Data_->optcontrol_cost = 0;
		vec xcurrent = x;//start from 0
		for (size_t i = 0; i < Data_->numphases_; i++)
		{

			uvec xindex(Data_->variable_indices[i]);

			size_t nodes = optpro_->GetPhase(i)->GetTotalNodes();
			int nstate = Data_->SIZES_[i][0];
			int ncontrols = Data_->SIZES_[i][1];
			int nparameters = Data_->SIZES_[i][2];
			int npaths = Data_->SIZES_[i][3];
			int nevents = Data_->SIZES_[i][4];
			std::vector<lp_index> state_index, control_index, time_index, paramemter_index;
			state_index = Data_->phase_indices[i]->state;
			control_index = Data_->phase_indices[i]->control;
			time_index = Data_->phase_indices[i]->time;
			paramemter_index = Data_->phase_indices[i]->parameter;
			uvec uvec_state_index(state_index), uvec_control_index(control_index), uvec_paramemter_index(paramemter_index);

			double t0 = xcurrent(time_index[0] - 1);
			double tf = xcurrent(time_index[1] - 1);
			vec state_vec = xcurrent(uvec_state_index - 1);
			vec control_vec = xcurrent(uvec_control_index - 1);
			vec tau_all = Data_->PS[i]->Points;
			tau_all = join_vert(tau_all, ones(1, 1));
			vec t_all = (tf - t0)*(tau_all + 1) / 2 + t0;
			mat state_matrix = reshape(state_vec, nodes + 1, nstate);
			mat control_matrix, control_matrixTotal;
			if (ncontrols > 0)
			{
				control_matrix = reshape(control_vec, nodes, ncontrols);
				rowvec control_tf = zeros(1, ncontrols);
				for (size_t j = 0; j < ncontrols; j++)
				{
					vec xdata = tau_all(span(0, tau_all.n_elem - 2));
					vec ydata = control_matrix.col(j);
					LpGuessChecker::spline_interpolation(control_tf.memptr() + j, 1.0, xdata, ydata);
				}
				control_matrixTotal = join_vert(control_matrix, control_tf);
			}
			vec parameter;
			if (nparameters > 0)
			{
				parameter = xcurrent(uvec_paramemter_index - 1);
			}
			/************************************************************************/
			/*Next  Contruct the costates                                                                     */
			/************************************************************************/
			vec lambacurrrent_vec = lamba.subvec(conshift, conshift + nodes*nstate - 1);
			mat lambacurrent = reshape(lambacurrrent_vec, nodes, nstate);
			mat inverse_weight_Matrix = diagmat(1 / Data_->PS[i]->Weights);
			mat costate_gauss = inverse_weight_Matrix*lambacurrent;
			sp_mat D = dsmatrix::GetArmaSP_Mat(Data_->PS[i]->D);
			mat costatef = trans(D.col(D.n_cols - 1))*lambacurrent;
			mat costate = -join_vert(costate_gauss, costatef);//The Ipopt 's returned lambda is different the book's,we need add a minus
			mat pathmultTotal;
			if (npaths > 0)
			{
				mat pathmult = Data_->nlpreturn_lambda.subvec(nodes*nstate, nodes*nstate + npaths*nodes - 1);
				pathmult = reshape(pathmult, nodes, npaths);
				for (size_t j = 0; j < npaths; j++)
				{
					pathmult.col(j) = 2 * inverse_weight_Matrix*pathmult.col(j) / (tf - t0);
				}
				if (ncontrols > 0)
				{
					rowvec pathmulF = zeros(1, npaths);
					for (size_t k = 0; k < npaths; k++)
					{
						vec xdata = tau_all(span(0, tau_all.n_elem - 2));
						vec ydata = pathmult.col(k);
						LpGuessChecker::spline_interpolation(pathmulF.memptr() + k, 1.0, xdata, ydata);
					}
					pathmultTotal = join_vert(pathmult, pathmulF);
				}
				else{
					rowvec pathmulF = zeros(1, npaths);
					for (size_t k = 0; k < npaths; k++)
					{
						vec xdata = tau_all(span(0, tau_all.n_elem - 2));
						vec ydata = pathmult.col(k);
						LpGuessChecker::spline_interpolation(pathmulF.memptr() + k, 1.0, xdata, ydata);
					}
					pathmultTotal = join_vert(pathmult, pathmulF);

				}
			}
			SolCost mycost;
			mycost.initial_state_ = t_all(0);
			mycost.initial_state_ = trans(state_matrix.row(0));
			mycost.terminal_time_ = t_all(t_all.n_elem - 1);
			mycost.terminal_state_ = trans(state_matrix.row(state_matrix.n_rows - 1));
			mycost.time_ = t_all;
			mycost.state_ = state_matrix;
			mycost.control_ = control_matrixTotal;
			mycost.parameter_ = parameter;
			mycost.phase_num_ = i + 1;
			vec lagrange;
			double mayer;
			Funs_->MayerCost(mycost, mayer);
			Funs_->LagrangeCost(mycost, lagrange);
			double mayer_cost = mayer;
			mat lagrange_cost = (tf - t0)* trans(Data_->PS[i]->Weights)*lagrange.subvec(0, lagrange.n_elem - 2) / 2.0;
			Data_->optcontrol_cost += mayer_cost + lagrange_cost(0);
			SolDae mydae;
			mydae.time_ = t_all;
			mydae.state_ = state_matrix;
			mydae.contol_ = control_matrixTotal;
			mydae.parameter_ = parameter;
			mydae.phase_num_ = i + 1;
			mat dae, path;
			Funs_->DaeFunction(mydae, dae, path);
			dae = join_horiz(dae, path);
			vec Hamiltonian = lagrange + sum(costate%dae(span::all, span(0, nstate - 1)), 1);
			varshift += Data_->variable_indices[i].size();
			conshift += Data_->constraint_indices[i].size();
			shared_ptr<SolutionData> result_data_i(new SolutionData());
			result_data_i->time = mydae.time_;
			result_data_i->state = mydae.state_;
			result_data_i->control = control_matrixTotal;
			result_data_i->parameter = mydae.parameter_;
			result_data_i->costate = costate;
			result_data_i->pathmult = pathmultTotal;
			result_data_i->Hamiltonian = Hamiltonian;
			result_data_i->mayerCost = mayer_cost;
			result_data_i->lagrangeCost = lagrange_cost(0);
			resultData[i] = (result_data_i);
			std::vector<double> vtimeguess(mydae.time_.n_elem);
			for (size_t j = 0; j < mydae.time_.n_elem; j++)
			{
				vtimeguess[j] = mydae.time_(j);
			}
			std::vector<std::vector<double>> vstateguess(nstate);
			for (size_t istate = 0; istate < nstate; istate++)
			{
				std::vector<double> jstate(mydae.state_.n_rows);
				for (size_t j = 0; j < jstate.size(); j++)
				{
					jstate[j] = mydae.state_(j, istate);
				}
				vstateguess[istate] = jstate;
			}
			std::vector<std::vector<double>> vcontrolguess(ncontrols);
			for (size_t icontrol = 0; icontrol < ncontrols; icontrol++)
			{
				std::vector<double> jcontrol(mydae.contol_.n_rows);
				for (size_t j = 0; j < jcontrol.size(); j++)
				{
					jcontrol[j] = mydae.contol_(j, icontrol);
				}
				vcontrolguess[icontrol] = jcontrol;
			}
			std::vector<double> vparameterguess(nparameters);
			for (size_t j = 0; j < mydae.parameter_.n_elem; j++)
			{
				vparameterguess[j] = mydae.parameter_[j];
			}
			optpro_->GetPhase(i)->GetStateGuess() = vstateguess;
			optpro_->GetPhase(i)->GetControlGuess() = vcontrolguess;
			optpro_->GetPhase(i)->GetTimeGuess() = vtimeguess;
			optpro_->GetPhase(i)->GetparameterGuess() = vparameterguess;
		}//for iphase
		Data_->result = resultData;
	};

	void Nlp2OpConverter::FinalResultSave(shared_ptr<LpCalculateData>& Data_)
	{
		char filename[20];
		for (size_t iphase = 0; iphase < Data_->numphases_; iphase++)
		{
			shared_ptr<SolutionData> result_data_i(Data_->result[iphase]);

			sprintf(filename, "time%u", iphase+1);
			result_data_i->time.save(filename, raw_ascii);

			sprintf(filename, "state%u", iphase + 1);
		    result_data_i->state.save(filename, raw_ascii);

			sprintf(filename, "control%u", iphase + 1);
			result_data_i->control.save(filename, raw_ascii);
			sprintf(filename, "parameter%u", iphase + 1);
			result_data_i->parameter.save(filename, raw_ascii);

			sprintf(filename, "costate%u", iphase + 1);
			result_data_i->costate.save(filename, raw_ascii);

			sprintf(filename, "Hamiltonian%u", iphase + 1);
			result_data_i->Hamiltonian.save(filename, raw_ascii);
		}

	}

}//namespace Lpopc