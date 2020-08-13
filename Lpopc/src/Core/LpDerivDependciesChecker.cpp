// Copyright (C) 2014-2015 Xue Zhichen, Wang Yujie,Wang Na 
// All Rights Reserved.
// This file is a part of LPOPC , published under the Eclipse Public License.
// Author:Xue Zhichen 7/12 2015   21:27
// Email:eddy_lpopc@163.com
#include "LpDerivDependciesChecker.h"
namespace Lpopc
{

	void DeriveDependicieshecker::GetDependiciesForJacobiInEveryPhase()
	{
		LP_DBG_START_FUN("DeriveDependicieshecker::GetDependiciesForJacobiInEveryPhase")
		std::vector<umat> dependices(Data_->numphases_);
		Jocobidependices = dependices;
		for (size_t iphase = 0; iphase < Data_->numphases_; iphase++)
		{
			int nstate = Data_->SIZES_[iphase][0];
			int ncontrols = Data_->SIZES_[iphase][1];
			int nparameters = Data_->SIZES_[iphase][2];
			int npaths = Data_->SIZES_[iphase][3];
			int sumnodes = Data_->PS[iphase]->Points.n_elem;
			umat iphasedependicies(nstate + npaths, ncontrols + nstate);

			int stateindex_start = Data_->phase_indices[iphase]->state[0] - 1;
			int stateindex_end = Data_->phase_indices[iphase]->state[Data_->phase_indices[iphase]->state.size() - 1] - 1;
			int controlindex_start = Data_->phase_indices[iphase]->control[0] - 1;
			int controlindex_end = Data_->phase_indices[iphase]->control[Data_->phase_indices[iphase]->control.size() - 1] - 1;
			vec state_vector = Data_->nlpGuessVector.subvec(stateindex_start, stateindex_end);
			vec control_vector;
			if (ncontrols > 0)
			{
				control_vector = Data_->nlpGuessVector.subvec(controlindex_start, controlindex_end);
			}
			int t0_index = Data_->phase_indices[iphase]->time[0] - 1;
			int tf_index = Data_->phase_indices[iphase]->time[1] - 1;
			double t0 = Data_->nlpGuessVector(t0_index);
			double tf = Data_->nlpGuessVector(tf_index);
			double tspan = tf - t0;
			vec t_radau = (Data_->PS[iphase]->Points + 1)*(tspan / 2.0) + t0;
			//t_radau.print("tradau");
			mat state_matrix = reshape(state_vector, sumnodes + 1, nstate);
			//state_matrix.print("state_matrix");
			mat state_radau = state_matrix.rows(0, state_matrix.n_rows - 2);//index start form 0

			vec x0 = trans(state_matrix.row(0));
			vec xf = trans(state_matrix.row(state_matrix.n_rows - 1));
			mat control_radau;
			if (ncontrols > 0)
			{
				control_radau = reshape(control_vector, sumnodes, ncontrols);
			}
			
			mat parameters;
			if (nparameters > 0)
			{
				int para_index_start = Data_->phase_indices[iphase]->parameter[0] - 1;
				int para_index_end = Data_->phase_indices[iphase]->parameter[Data_->phase_indices[iphase]->parameter.size() - 1] - 1;
				parameters = Data_->nlpGuessVector.rows(para_index_start, para_index_end);
			}
			//get derivatives of DAE function
			SolDae mydae;
			mydae.time_ = t_radau.row(1);
			mydae.state_ = state_radau.row(1);
			mydae.contol_ = control_radau.row(1);
			mydae.parameter_ = parameters;
			mydae.phase_num_ = iphase + 1;
			mat dae, path;
			for (size_t istate = 0; istate < nstate; istate++)
			{
					mydae.state_.col(istate) .fill( datum::nan);
					funs_->DaeFunction(mydae, dae, path);
					urowvec tem=zeros<urowvec>(iphasedependicies.n_rows);
					tem(find_nonfinite(join_horiz(dae, path))).fill(1);
					iphasedependicies.col(istate) = trans(tem);
					mydae.state_.col(istate).fill( state_radau(1,istate));
			}
			for (size_t icontrol = 0; icontrol < ncontrols; icontrol++)
			{
				mydae.contol_(0, icontrol) = datum::nan;
				funs_->DaeFunction(mydae, dae, path);
				urowvec tem = zeros<urowvec>(iphasedependicies.n_rows);
				
				tem(find_nonfinite(join_horiz(dae, path))).fill(1);
				iphasedependicies.col(nstate+icontrol) =trans( tem);
				mydae.contol_(0, icontrol) =control_radau(1,icontrol);
			}
			//////////////////////////////////////////////////////////////////////////
			//Xue Zhichen :
			// Since the  finite-difference methods approximating the hessian doesn't work very well
			// I turn off obtaining dependices using sparse nan,which need more testing.
			//iphasedependicies.fill(1.0);//Fix Me
			Jocobidependices[iphase] = iphasedependicies;
		}//for iphase
		

	}

}