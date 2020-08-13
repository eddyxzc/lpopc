// Copyright (C) 2014-2015 Xue Zhichen, Wang Yujie,Wang Na 
// All Rights Reserved.
// This file is a part of LPOPC , published under the Eclipse Public License.
// Author:Xue Zhichen 7/11 2015   22:57
// Email:eddy_lpopc@163.com
#include "LpSolutionError.h"
namespace Lpopc
{
	
	void SolutionErrorChecker::BarLagrangeInterp(const vec& data_x, const vec& data_y, const vec& x, vec& y)
	{
		/************************************************************************/
		/* data_x is interpolate nodes,data_y is function values of nodes,x is your nodes,y is the interpolated data                                                                  */
		/************************************************************************/
		LP_ASSERT_EXCEPTION(data_x.n_elem == data_y.n_elem, LPOPC_SOLUTION_ERROR,
			"The size of data_x is different from the size of data_y in BarLagrangeInterp ");

		size_t M = data_x.n_elem;
		size_t N = x.n_elem;
		y = zeros(N);
		// Compute the barycentric weights
		mat X = repmat(data_x, 1, M);
		// matrix of weights
		mat W = repmat(1/ prod(X - trans(X)+eye(M,M),0),N,1);
		// Get distances between nodes and interpolation points
		mat xdist = repmat(x, 1, M) - repmat(trans(data_x),N,1);

		// Find all of the elements where the interpolation point is on a node
		uvec fixindex=find(xdist == 0);
		uvec fixi(fixindex.n_elem), fixj(fixindex.n_elem);
		size_t k = 0;
		for (size_t i = 0; i < xdist.n_cols; i++)
		{
			for (size_t j = 0; j < xdist.n_rows; j++)
			{
				if (xdist(j, i) == 0)
				{
					fixi[k] = i;
					fixj[k] = j;
					k++;
				}
			}
		}
		// Use NaNs as a place - holder
		xdist(fixindex).fill( datum::nan);
		mat H = W/ xdist;
		// Compute the interpolated polynomial
		mat data = join_horiz(data_x, data_y);
		y = (H*data_y)/sum(H, 1);
		// Replace NaNs with the given exact values.
		y(fixj) = data_y(fixi);
	}

	void SolutionErrorChecker::SolutionInterpolation(const int iphase,vec& vec_time, mat& mat_state, mat& mat_control)const
	{
		uvec nodesCum(optpro_->GetPhase(iphase)->GetNodesPerInterval());
		uvec indices = join_vert(zeros<uvec>(1),cumsum(nodesCum) );
		vec meshPoints(optpro_->GetPhase(iphase)->GetMeshPoints());
		vec tau = join_vert(Data_->PS[iphase]->Points, ones<vec>(1));
		size_t nstate = Data_->SIZES_[iphase][0];
		size_t ncontrol = Data_->SIZES_[iphase][1];
		for (size_t segindex = 0; segindex < nodesCum.n_elem; segindex++)
		{
			size_t n = optpro_->GetPhase(iphase)->GetNodesPerInterval()[segindex];

			size_t istart = indices(segindex);
			size_t ifinish = indices(segindex + 1);
			double time0 = tau(istart), timef = tau(ifinish);
			vec timecurrent = tau.subvec(istart, ifinish);
			mat statecurrent = Data_->result[iphase]->state.rows(istart, ifinish);
			mat controlcurrent;
			if (ncontrol>0)
			{
				controlcurrent = Data_->result[iphase]->control.rows(istart, ifinish-1);
			}
			vec rpmpoints, rpmweights;
			RPMGenerator::GetLGRPoints(n + 1, rpmpoints, rpmweights);
			vec ttem = join_vert((rpmpoints + 1)*(timef - time0) / 2 + time0, meshPoints.row(segindex + 1));
			mat stem = zeros(ttem.n_elem, nstate);
			mat ctem;
			for (size_t i = 0; i < nstate; i++)
			{
				vec stemcoli;
				BarLagrangeInterp(timecurrent, statecurrent.col(i), ttem, stemcoli);
				stem.col(i) = stemcoli;
			}
			if (ncontrol>0)
			{
				ctem = zeros(ttem.n_elem-1, ncontrol);
				for (size_t i = 0; i < ncontrol; i++)
				{
					vec ctemcoli;
					BarLagrangeInterp(timecurrent.subvec(0,timecurrent.n_elem-2), controlcurrent.col(i), ttem.subvec(0,ttem.n_rows-2), ctemcoli);
					ctem.col(i) = ctemcoli;
				}
			}//if ncontrol
			if (segindex==0)
			{
				vec_time = ttem.subvec(0, ttem.n_elem - 2);
				mat_state = stem.rows(0, stem.n_rows - 2);
				mat_control = ctem;
			}else{
				vec_time = join_vert(vec_time, ttem.subvec(0, ttem.n_elem - 2));
				mat_state = join_vert(mat_state, stem.rows(0, stem.n_rows - 2));
				if (ncontrol>0)mat_control = join_vert(mat_control, ctem);
			}
		}//for segindex
		vec_time = join_vert(vec_time, ones<vec>(1));
		mat_state = join_vert(mat_state, Data_->result[iphase]->state.tail_rows(1));
	}

	void SolutionErrorChecker::CheckSolutionDiffError(int iphase, mat& relative_error)const
	{
		vec temTime; mat temState; mat temControl;
		/************************************************************************/
		/* ------------------------------------------------------------------
		 Interpolate the solution on a mesh that consists of one more LGR 
		 point in each mesh interval than was used to solve the NLP       
		/************************************************************************/
		SolutionInterpolation(iphase, temTime, temState, temControl);
		size_t nstate = Data_->SIZES_[iphase][0];
		size_t ncontrol = Data_->SIZES_[iphase][1];
		vec phase_time = Data_->result[iphase]->time;
		double t0 = phase_time(0);
		double tf = phase_time(phase_time.n_elem - 1);

		SolDae mydae;
		mydae.time_ = (tf - t0) / 2 * temTime.subvec(0, temTime.n_elem - 2) + (tf - t0) / 2;
		mydae.state_ = temState.rows(0, temState.n_rows - 2);
		mydae.phase_num_ = iphase+1;
		mydae.contol_ = temControl;
		mydae.parameter_ = Data_->result[iphase]->parameter;
		mat daeout, pathout;
		Funs_->DaeFunction(mydae, daeout, pathout);
		daeout *= (tf - t0) / 2.0;
		vec meshPoints(optpro_->GetPhase(iphase)->GetMeshPoints());
		uvec nodesPerInterval(optpro_->GetPhase(iphase)->GetNodesPerInterval());
		nodesPerInterval += 1;
		size_t numMeshIntervals = nodesPerInterval.n_elem;
		shared_ptr<RPMGenerator> RPM(new RPMGenerator());
		std::vector<lp_index> stdvec_nodesPerInterval(nodesPerInterval.n_elem);
		std::vector<double>stdvec_meshpoint(meshPoints.n_elem);
		for (size_t i = 0; i < nodesPerInterval.n_elem; i++)
		{
			stdvec_nodesPerInterval[i] = nodesPerInterval[i];
		}
		for (size_t i = 0; i < meshPoints.n_elem; i++)
		{
			stdvec_meshpoint[i] = meshPoints[i];
		}
		//daeout.save("daeot", raw_ascii);
		RPM->initialize(numMeshIntervals, stdvec_meshpoint, stdvec_nodesPerInterval);
		mat integratedRHS = join_vert(temState.row(0),RPM->UnityMatrix()*temState + RPM->IntegrationMatrix()*daeout);
		mat absolute_error = abs(integratedRHS - temState);
		// absolute_error=(nLGR+1)*nstate
		//
		//
		//RPM->IntegrationMatrix().print("intermatrix");
		//integratedRHS.save("integratedRHS", raw_ascii);
		//temState.save("temstate", raw_ascii);
		
		rowvec maxofeachcol = 1 + max(temState, 0);
		relative_error = absolute_error;
		for (size_t i = 0; i < relative_error.n_cols; i++)
		{
			relative_error.col(i) /= maxofeachcol(i);
		}
		//relative_error.save("relateerror", raw_ascii);
	}

	void SolutionErrorChecker::CheckSolutionDiffError(int iphase, mat& relative_error, mat& state_with_one_more_point) const
	{
		vec temTime; mat temState; mat temControl;
		/************************************************************************/
		/* ------------------------------------------------------------------
		Interpolate the solution on a mesh that consists of one more LGR
		point in each mesh interval than was used to solve the NLP
		/************************************************************************/
		SolutionInterpolation(iphase, temTime, temState, temControl);
		size_t nstate = Data_->SIZES_[iphase][0];
		size_t ncontrol = Data_->SIZES_[iphase][1];
		vec phase_time = Data_->result[iphase]->time;
		double t0 = phase_time(0);
		double tf = phase_time(phase_time.n_elem - 1);

		SolDae mydae;
		mydae.time_ = (tf - t0) / 2 * temTime.subvec(0, temTime.n_elem - 2) + (tf - t0) / 2;
		mydae.state_ = temState.rows(0, temState.n_rows - 2);
		mydae.phase_num_ = iphase + 1;
		mydae.contol_ = temControl;
		mydae.parameter_ = Data_->result[iphase]->parameter;
		mat daeout, pathout;
		Funs_->DaeFunction(mydae, daeout, pathout);
		daeout *= (tf - t0) / 2.0;
		vec meshPoints(optpro_->GetPhase(iphase)->GetMeshPoints());
		uvec nodesPerInterval(optpro_->GetPhase(iphase)->GetNodesPerInterval());
		nodesPerInterval += 1;
		size_t numMeshIntervals = nodesPerInterval.n_elem;
		shared_ptr<RPMGenerator> RPM(new RPMGenerator());
		std::vector<lp_index> stdvec_nodesPerInterval(nodesPerInterval.n_elem);
		std::vector<double>stdvec_meshpoint(meshPoints.n_elem);
		for (size_t i = 0; i < nodesPerInterval.n_elem; i++)
		{
			stdvec_nodesPerInterval[i] = nodesPerInterval[i];
		}
		for (size_t i = 0; i < meshPoints.n_elem; i++)
		{
			stdvec_meshpoint[i] = meshPoints[i];
		}
		//daeout.save("daeot", raw_ascii);
		RPM->initialize(numMeshIntervals, stdvec_meshpoint, stdvec_nodesPerInterval);
		mat integratedRHS = join_vert(temState.row(0), RPM->UnityMatrix()*temState + RPM->IntegrationMatrix()*daeout);
		mat absolute_error = abs(integratedRHS - temState);
		// absolute_error=(nLGR+1)*nstate
		//
		//
		//RPM->IntegrationMatrix().print("intermatrix");
		//integratedRHS.save("integratedRHS", raw_ascii);
		//temState.save("temstate", raw_ascii);

		rowvec maxofeachcol = 1 + max(temState, 0);
		relative_error = absolute_error;
		for (size_t i = 0; i < relative_error.n_cols; i++)
		{
			relative_error.col(i) /= maxofeachcol(i);
		}
		//relative_error.save("relateerror", raw_ascii);
		state_with_one_more_point = temState;
	}

}//namespace Lpopc


