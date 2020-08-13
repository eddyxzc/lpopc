#include "LpLiuHpMeshRefineAlg.hpp"
#include "RPMGenerator.hpp"
#include "LpSolutionError.h"
#include <utility>
#include <cmath>
#include <algorithm>
#include "LpOptimalProblem.hpp"
namespace Lpopc
{
	std::map<int, mat> LiuHpMeshRefineAlg::alj_map_ = std::map<int, mat>();

	bool LiuHpMeshRefineAlg::RefineMesh(shared_ptr<OptimalProblem>optpro, std::vector<shared_ptr<LpMesh>>& allMeshVector)
	{
		size_t nphases = optpro->GetPhaseNum();
		shared_ptr<SolutionErrorChecker> error_checker(new SolutionErrorChecker(Funs_, Data_, optpro));
		bool NoMoreRefine = true;
		allMeshVector = std::vector<shared_ptr<LpMesh>>(nphases);
		std::vector<std::vector<MeshOperation>> refine_tag(nphases);
		std::vector<mat> phase_state(nphases);
		std::vector<vec> phase_mesh_point(nphases);
		std::vector<shared_ptr<LpMesh>> beforemeshInfo(nphases);
		if (mesh_index_ == 0)
		{//stor user set fist mesh onformation
			for (size_t iphase = 0; iphase < nphases; iphase++)
			{
				beforemeshInfo[iphase].reset(new LpMesh());
				beforemeshInfo[iphase]->nodesPerInterval = optpro->GetPhase(iphase)->GetNodesPerInterval();
				beforemeshInfo[iphase]->meshpoints = optpro->GetPhase(iphase)->GetMeshPoints();
				beforemeshInfo[iphase]->e_k = zeros<vec>(optpro->GetPhase(iphase)->GetNodesPerInterval().size());
				//No e_k of User set
			}
			mesh_history_.push_back(beforemeshInfo);
		}

		beforemeshInfo = mesh_history_[mesh_history_.size() - 1];
		for (size_t iphase = 0; iphase < nphases; iphase++)
		{

			mat relative_error, temstate;//state with one more LGR point.
			error_checker->CheckSolutionDiffError(iphase, relative_error, temstate);
			//relative_error=(nLGR+1)*nstate
			//
			//
			//relative_error.save("relative_error", raw_ascii);
			size_t numseg = optpro->GetPhase(iphase)->GetNodesPerInterval().size();
			//for one more point
			uvec nodesPerInterval(optpro->GetPhase(iphase)->GetNodesPerInterval());
			nodesPerInterval += 1;
			uvec indices = join_vert(zeros<uvec>(1), cumsum(nodesPerInterval));

			//the actual index before inserting one more point
			uvec actual_nodesPerInterval(optpro->GetPhase(iphase)->GetNodesPerInterval());
			uvec state_indices = join_vert(zeros<uvec>(1), cumsum(actual_nodesPerInterval));

			vec meshpoints = optpro->GetPhase(iphase)->GetMeshPoints();

			std::vector<shared_ptr<LpMesh>> tempMeshVector(numseg);
			std::vector<uword> newnodesPerInterval;
			std::vector<double>   new_e_k;
			std::vector<double> newmeshpoints;
			newmeshpoints.push_back(-1);
			mat intervalState = Data_->result[iphase]->state;
			phase_state[iphase] = intervalState;
			std::vector<MeshOperation> phase_tag(numseg);

			std::vector<mat> segment_state(numseg);
			for (size_t iseg = 0; iseg < numseg; iseg++)
			{
				size_t istart = indices(iseg);
				size_t ifinish = indices(iseg + 1);
				tempMeshVector[iseg].reset(new LpMesh());

				mat intervalError = relative_error.rows(istart, ifinish);
				rowvec  maxofstates = max(intervalError, 0);
				bool errorsatisfied = maxofstates.max() <= mesh_tol_;
				beforemeshInfo[iphase]->e_k(iseg) = max(maxofstates.max());

				
				if (errorsatisfied)
				{
					// satisfied.
					//1.reducing the degree of polynomial approximation
					//2.merging mesh intervals(after 1)
					// If all mesh satisfied,and neither 1 nor 2 satisfied,then mesh refinement is completed.
					reporter_->Printf()->info("Error Tolerance Satisfied in Segment {} of Phase {}\n", iseg + 1, iphase + 1);
					double  meshpoint0 = optpro->GetPhase(iphase)->GetMeshPoints()[iseg];
					double  meshpointf = optpro->GetPhase(iphase)->GetMeshPoints()[iseg + 1];



					//1.reducing the degree of polynomial approximation
					//void Reducing_N(size_t iphase, size_t segindex, uvec state_index, vec betai, uword return_N, shared_ptr<OptimalProblem> optpro);
					rowvec betai = 1 + max(intervalState, 0);
					uword Nneed = 0;
					uword Nactual = actual_nodesPerInterval(iseg);
					Reducing_N(iphase, iseg, state_indices, betai, Nneed, optpro);
					tempMeshVector[iseg]->nodesPerInterval = Nneed*ones<uvec>(1);
					tempMeshVector[iseg]->meshpoints = linspace(meshpoint0, meshpointf, 2);

					if (Nneed==Nactual)
					{
						phase_tag[iseg] = SATISFIED;
					}
					else
					{
						phase_tag[iseg] = REDUCED;
						NoMoreRefine = false;
					}
					//tag?///////////////////////
					
				}
				else
				{
					// Not all satisfied ,modify segment 
					//1.increasing the degree of the polynomial approximation
					//2.dividing a mesh interval
					reporter_->Printf()->info("Error Tolerance Not Satisfied in Segment {} of Phase {}\n", iseg + 1, iphase + 1);

					double  meshpoint0 = optpro->GetPhase(iphase)->GetMeshPoints()[iseg];
					double  meshpointf = optpro->GetPhase(iphase)->GetMeshPoints()[iseg + 1];
					if (mesh_index_ == 0)
					{
						//Create Second Mesh
						//Just add three collocation points
						int nodes_current = optpro->GetPhase(iphase)->GetNodesPerInterval()[iseg];
						int newnodes = nodes_current + 3;


						tempMeshVector[iseg]->nodesPerInterval = newnodes*ones<uvec>(1);
						tempMeshVector[iseg]->meshpoints = linspace(meshpoint0, meshpointf, 2);
					}
					else
					{
						//using method 1 or method 2

						bool ifincreaseN = CanWeIncreaseN(iphase, iseg, state_indices);
						bool ifdivide = !ifincreaseN;
						if (ifincreaseN)
						{
							//increase N
							uword Nneed = Increasing_N(iphase, iseg, state_indices, maxofstates.max(), optpro);
							if (Nneed > Nmax_)
							{
								ifdivide = true;
							}
							else
							{
								tempMeshVector[iseg]->nodesPerInterval = Nneed*ones<uvec>(1);
								tempMeshVector[iseg]->meshpoints = linspace(meshpoint0, meshpointf, 2);
							}
						}
						if (ifdivide)
						{
							//divide mesh
							uword diviedH = Dividing_mesh(iphase, iseg, state_indices, maxofstates.max(), optpro);
							uword  nodes_current = optpro->GetPhase(iphase)->GetNodesPerInterval()[iseg];
							tempMeshVector[iseg]->nodesPerInterval = nodes_current*ones<uvec>(diviedH);
							tempMeshVector[iseg]->meshpoints = linspace(meshpoint0, meshpointf, diviedH + 1);
						}
					}

					phase_tag[iseg] = NOT_SATISFIED;
					NoMoreRefine = false;

				}
				

			}//for iseg

			auto reservedTempMesh = tempMeshVector;
			if (NoMoreRefine)continue;//all good No need merge
			
			uword refinedsegment_index = 0;
			uvec refined_state_indices;
			vec refined_meshpoints;
			auto itormesh = tempMeshVector.begin();
			auto itorstate = segment_state.begin();
			auto itortag = phase_tag.begin();
			uword seg_x_start =0;
			uword seg_x_finish=0 ;

			for (size_t iseg = 0; iseg < numseg; iseg++)//start with second mesh
			{
				 seg_x_start = state_indices(iseg);
				 seg_x_finish = state_indices(iseg + 1);

				segment_state[refinedsegment_index] = intervalState.rows(seg_x_start, seg_x_finish);

				if (iseg>0 &&(phase_tag[refinedsegment_index] != NOT_SATISFIED) && (phase_tag[refinedsegment_index - 1] != NOT_SATISFIED))
				{
					

					uword N_k = tempMeshVector[refinedsegment_index]->nodesPerInterval(0);
					uword N_km1 = tempMeshVector[refinedsegment_index -1]->nodesPerInterval(0);
					if (N_k == N_km1)
					{
						
						bool merged = Merging_mesh(iphase, refinedsegment_index, tempMeshVector,segment_state, N_k, N_km1);

						NoMoreRefine = false;
// 						auto itorpoint = newmeshpoints.begin();
// 						newmeshpoints.erase(itorpoint + refinedsegment_index);
// 						auto itornodes = newnodesPerInterval.begin();
// 						newnodesPerInterval.erase(itornodes + refinedsegment_index);
		//Number of meshpoints in a segment tagged SATISFIED mush be 2
						tempMeshVector[refinedsegment_index - 1]->meshpoints(1) = tempMeshVector[refinedsegment_index]->meshpoints(1);
						tempMeshVector.erase(itormesh + refinedsegment_index);

						seg_x_start = state_indices(iseg-1);
						seg_x_finish = state_indices(iseg + 1);

						segment_state[refinedsegment_index-1] = intervalState.rows(seg_x_start, seg_x_finish);
						segment_state.erase(itorstate + refinedsegment_index);

						phase_tag.erase(itortag+ refinedsegment_index);
						phase_tag[refinedsegment_index - 1] = MERGED;
						
					}
					else {
						// needn't merge
						refinedsegment_index++;
					}
				}
				else
				{
					refinedsegment_index++;
				}//if all satisfied
			}
			for (size_t iseg = 0; iseg < tempMeshVector.size(); iseg++)
			{
				for (size_t i = 1; i < tempMeshVector[iseg]->meshpoints.n_elem; i++)
				{
					newmeshpoints.push_back(tempMeshVector[iseg]->meshpoints(i));
				}
				for (size_t i = 0; i < tempMeshVector[iseg]->nodesPerInterval.n_elem; i++)
				{
					newnodesPerInterval.push_back(tempMeshVector[iseg]->nodesPerInterval(i));

				}
			}
			optpro->GetPhase(iphase)->GetMeshPoints() = newmeshpoints;
			optpro->GetPhase(iphase)->GetNodesPerInterval() = newnodesPerInterval;
			allMeshVector[iphase].reset(new LpMesh());
			allMeshVector[iphase]->meshpoints = newmeshpoints;
			allMeshVector[iphase]->nodesPerInterval = newnodesPerInterval;
			allMeshVector[iphase]->e_k = zeros<vec>(newnodesPerInterval.size());
			refine_tag[iphase] = phase_tag;
			phase_mesh_point[iphase] = newmeshpoints;
		}//for iphase
		mesh_history_.push_back(allMeshVector);


		tag_history_.push_back(refine_tag);

		state_history_.push_back(phase_state);
		mesh_points_history_.push_back(phase_mesh_point);

		mesh_index_++;
		return NoMoreRefine;
	}


	void LiuHpMeshRefineAlg::GetLagrangeInterpPowerCoefficients(int N, mat& alj)
	{
		std::map<int, mat>::iterator aitor = alj_map_.find(N);

		if ((aitor != alj_map_.end()))
		{

			alj = aitor->second;
		}
		else
		{
			///Calculate new coefficients
			GetLagrangeInterpCoefficientsImpl(N, alj);
			auto ret1 = alj_map_.insert(std::make_pair(N, alj));
			if (!ret1.second)
			{
				LP_THROW_EXCEPTION(LIU_HP_ERROR, "Insert lagrange cofficients failed!");
			}

		}
	}

	void LiuHpMeshRefineAlg::CalculateDi(vec x, vec& D_i)
	{
		// Based on
		//������, ����. �������ղ�ֵ��ʽ����ȫչ��[J]. ͨ��ʦ��ѧԺѧ��, 
		//2007, 28:10-12. DOI:doi:10.3969/j.issn.1008-7974.2007.02.003.
		lp_index N = x.n_elem;
		mat T = zeros(N, N);
		T.row(0) = trans(x);
		for (size_t i = 1; i < N; i++)
		{
			for (int j = N - 2; j >= 0; j--)
			{
				T(i, j) = T(i, j + 1) + T(i - 1, j + 1);
			}
			T.row(i) = T.row(i) % trans(x);
		}
		D_i = ones(N + 1);
		D_i.subvec(1, N) = sum(T, 1);
	}

	void LiuHpMeshRefineAlg::GetLagrangeInterpCoefficientsImpl(int N, mat& alj)
	{
		vec x, w, temx, Di;
		alj = zeros(N + 1, N + 1);
		RPMGenerator::GetLGRPoints(N, x, w);
		x = join_vert(x, ones(1, 1));
		for (size_t i = 0; i < N + 1; i++)
		{
			temx = x;
			temx(i) = datum::nan;
			//remove xi;
			temx = temx(find_finite(temx));
			CalculateDi(-temx, Di);
			vec powxi = ones(N + 1);
			for (int k = N - 1; k >= 0; k--)
			{
				powxi(k) = powxi(k + 1)*x(i);
			}
			alj.col(i) = Di / sum(powxi%Di);
		}

	}

	void LiuHpMeshRefineAlg::ModifySegment(shared_ptr<OptimalProblem> optpro, size_t iphase, size_t segindex, mat seg_error, uvec state_index, shared_ptr<LpMesh> newmesh) const
	{
		//Find which state has the max error
		rowvec  maxofstates = max(seg_error, 0);
		uword which_state_max;
		maxofstates.max(which_state_max);

		uword istart = state_index(segindex);
		uword iend = state_index(segindex + 1);

		vec state_max_error = Data_->result[iphase]->state(span(istart, iend), which_state_max);


	}





	uword LiuHpMeshRefineAlg::Dividing_mesh(size_t iphase, size_t segindex, uvec state_index, double seg_error, shared_ptr<OptimalProblem> optpro)
	{

		////////////////////////////Find q////////////////////////////////////

		vec meshpointbefore = mesh_history_[mesh_history_.size() - 2][iphase]->meshpoints;
		double t0 = optpro->GetPhase(iphase)->GetMeshPoints()[segindex];
		double tf = optpro->GetPhase(iphase)->GetMeshPoints()[segindex + 1];

		double e_k = seg_error;
		//double e_k_b = mesh_history_[mesh_history_.size() - 2][iphase]->e_k(segindex);

		double h = tf - t0;
		uword itmin_b;
		uword itmax_b;

		//Here tf<>t0
		itmin_b = max(find(meshpointbefore <= t0));
		itmax_b = min(find(meshpointbefore >= tf));



		double h_b = meshpointbefore(itmax_b) - meshpointbefore(itmin_b);

		uword N = optpro->GetPhase(iphase)->GetNodesPerInterval()[segindex];
		uvec  nodePerIntervalbefore = mesh_history_[mesh_history_.size() - 2][iphase]->nodesPerInterval;

		uword N_b = 0;
		double e_k_b = max(mesh_history_[mesh_history_.size() - 2][iphase]->e_k.subvec(itmin_b, itmax_b - 1));
		for (size_t i = 0; i < itmax_b - itmin_b; i++)
		{
			N_b += nodePerIntervalbefore[itmin_b + i];
		}

		//Set some var
		double fN = (double)N / (double)N_b;
		double fh = h / h_b;
		double fe = e_k / e_k_b;
		//calculate q
		double q = ceil(log((fe / pow(N, 5.0 / 2.0))) / log(fh / fN));
		//////////////////////////////////////////////////////////////////////////
		uword divided_H = (uword)ceil(pow(e_k / mesh_tol_, 1 / q));

		uword divided_H_max = (uword)ceil(log(e_k / mesh_tol_) / log(N));

		uword S = std::min(divided_H, divided_H_max);
		S = std::max(S, (uword)2);
		return S;
	}

	uword LiuHpMeshRefineAlg::Increasing_N(size_t iphase, size_t segindex, uvec state_index, double seg_error, shared_ptr<OptimalProblem> optpro)
	{
		////////////////////////////Find q////////////////////////////////////

		vec meshpointbefore = mesh_history_[mesh_history_.size() - 2][iphase]->meshpoints;
		double t0 = optpro->GetPhase(iphase)->GetMeshPoints()[segindex];
		double tf = optpro->GetPhase(iphase)->GetMeshPoints()[segindex + 1];

		double e_k = seg_error;
		//double e_k_b = mesh_history_[mesh_history_.size() - 2][iphase]->e_k(segindex);

		double h = tf - t0;
		uword itmin_b;
		uword itmax_b;

		//Here tf<>t0
		itmin_b = max(find(meshpointbefore <= t0));
		itmax_b = min(find(meshpointbefore >= tf));



		double h_b = meshpointbefore(itmax_b) - meshpointbefore(itmin_b);

		uword N = optpro->GetPhase(iphase)->GetNodesPerInterval()[segindex];
		uvec  nodePerIntervalbefore = mesh_history_[mesh_history_.size() - 2][iphase]->nodesPerInterval;

		uword N_b = 0;
		double e_k_b = max(mesh_history_[mesh_history_.size() - 2][iphase]->e_k.subvec(itmin_b, itmax_b - 1));
		for (size_t i = 0; i < itmax_b - itmin_b; i++)
		{
			N_b += nodePerIntervalbefore[itmin_b + i];
		}

		//Set some var
		double fN = (double)N / (double)N_b;
		double fh = h / h_b;
		double fe = e_k / e_k_b;
		//calculate q
		double q = ceil(log((fe / pow(N, 5.0 / 2.0))) / log(fh / fN));
		//////////////////////////////////////////////////////////////////////////

		uword Nk_next = static_cast<uword>(ceil(N*(pow(e_k / mesh_tol_, 1.0 / (q - 5.0 / 2.0)))));
		return Nk_next;
		//note Nmax!!!!!!!!!!!!!!!!!!!!!!1
	}

	void LiuHpMeshRefineAlg::Reducing_N(size_t iphase, size_t segindex, uvec state_index, rowvec betai, uword& return_N, shared_ptr<OptimalProblem> optpro)
	{
		uword nodes_current = optpro->GetPhase(iphase)->GetNodesPerInterval()[segindex];
		mat alj;
		GetLagrangeInterpPowerCoefficients(nodes_current, alj);
		//alj=(N+1)*(N+1);
		uword istart = state_index(segindex);
		uword iend = state_index(segindex + 1);
		mat segment_state = Data_->result[iphase]->state.rows(istart, iend);

		mat bil = alj*segment_state;
		/**
		 bij.
					state1,state2 ....
		power high
		 ...
		power low
		*/
		double beta_of_x = 0.;
		mat normalized_bil = bil;
		uvec less_than_eps;
		uvec maxN_of_state(segment_state.n_cols);
		uword N_interp_points = bil.n_rows;
		for (size_t i = 0; i < segment_state.n_cols; i++)
		{
			beta_of_x = betai(i);
			normalized_bil.col(i) /= beta_of_x;
			less_than_eps = find(normalized_bil.col(i) > mesh_tol_);
			if (less_than_eps.n_elem == 0)
			{
				maxN_of_state(i) = 1;
			}
			else
			{
				maxN_of_state(i) = N_interp_points - 1 - less_than_eps.min();
			}
		}
		return_N = std::max(uword(2), maxN_of_state.max());
	}

	bool LiuHpMeshRefineAlg::Merging_mesh(size_t iphase, size_t segindex,
		std::vector<shared_ptr<LpMesh>>& tempmesh,std::vector<mat> tempstate, uword N_k, uword N_km1)
	{
		//uvec state_index;////////
		//vec meshpoins;
		//Find pair adjacent mesh interval;begin with second segment


			//check N_{k}=N_{K-1}i
			//uword N_k= ;//N_K
			//uword N_km1;//N_{k-1}
		if (N_k == N_km1)
		{

			double hk = tempmesh[segindex]->meshpoints(1) - tempmesh[segindex]->meshpoints(0);
			double hk_b = tempmesh[segindex-1]->meshpoints(1) - tempmesh[segindex-1]->meshpoints(0);
			double hm = std::min(hk, hk_b);

			std::vector<uword> refinednodesPerInterval;
			for (size_t iseg = 0; iseg < tempmesh.size(); iseg++)
			{
				for (size_t i = 0; i < tempmesh[iseg]->nodesPerInterval.n_elem; i++)
				{
					refinednodesPerInterval.push_back(tempmesh[iseg]->nodesPerInterval(i));
				}
			}
			uvec state_index = refinednodesPerInterval;
			state_index =join_vert(zeros<uvec>(1), cumsum(state_index));
			mat alj;
			GetLagrangeInterpPowerCoefficients(N_k, alj);
			//alj=(N+1)*(N+1);

			
			mat original_segment_state = tempstate[segindex];
			mat segment_state;
			uword nstate = original_segment_state.n_cols;
			if (original_segment_state.n_rows==N_k+1)
			{
				segment_state = original_segment_state;
			}
			else 
			{
				vec changed_meshpoints, temp_w, return_interp_col;
				segment_state = zeros(N_k+1, nstate);
				RPMGenerator::GetLGRPoints(N_k, changed_meshpoints, temp_w);
				vec original_segment_time;
				RPMGenerator::GetLGRPoints(original_segment_state.n_rows-1, original_segment_time, temp_w);
				original_segment_time = join_vert(original_segment_time, ones<vec>(1));

				changed_meshpoints = join_vert(changed_meshpoints, ones<vec>(1));
				for (size_t istate = 0; istate < original_segment_state.n_cols; istate++)
				{
					SolutionErrorChecker::BarLagrangeInterp(original_segment_time, original_segment_state.col(istate), changed_meshpoints, return_interp_col);
					segment_state.col(istate) = return_interp_col;
				}
			}

			mat bil = alj*segment_state;
			vec hm_d_hk = ones(bil.n_rows);
			for (size_t ielem = 0; ielem < hm_d_hk.n_elem; ielem++)
			{
				hm_d_hk(ielem) *= 2 * (hm / hk);
			}

			bil.each_col() %= hm_d_hk;
			///////////////////////get bil_b///////////////////////////////////

			mat alj_b = alj;
			//GetLagrangeInterpPowerCoefficients(nodes_current_b, alj_b);
			//alj=(N+1)*(N+1);
			uword istart_b = state_index(segindex - 1);
			uword iend_b = state_index(segindex);
			mat original_segment_state_b = tempstate[segindex-1];

			//////////////////////////////////////////////////////////////////////////
			mat segment_state_b;
			if (original_segment_state_b.n_rows == N_k+1)
			{
				segment_state_b = original_segment_state_b;
			}
			else
			{
				vec changed_meshpoints, temp_w, return_interp_col;
				segment_state_b = zeros(N_k+1, nstate);
				RPMGenerator::GetLGRPoints(N_k, changed_meshpoints, temp_w);
				vec original_segment_time_b;
				RPMGenerator::GetLGRPoints(original_segment_state_b.n_rows - 1, original_segment_time_b, temp_w);
				original_segment_time_b = join_vert(original_segment_time_b, ones<vec>(1));

				changed_meshpoints = join_vert(changed_meshpoints, ones<vec>(1));
				for (size_t istate = 0; istate < original_segment_state.n_cols; istate++)
				{
					SolutionErrorChecker::BarLagrangeInterp(original_segment_time_b, original_segment_state_b.col(istate), changed_meshpoints, return_interp_col);
					segment_state_b.col(istate) = return_interp_col;
				}
			}
			/////////////////////////////////////////////////////////////////////////
			mat bil_b = alj_b*segment_state_b;
			vec hmb_d_hk = ones(bil_b.n_rows);
			for (size_t ielem = 0; ielem < hmb_d_hk.n_elem; ielem++)
			{
				hmb_d_hk(ielem) *= 2 * (hm / hk_b);
			}

			bil_b.each_col() %= hmb_d_hk;
			///////////////////////////////////////////////////////
			///////////////////Find beta//////////////////////////
			mat intervalState = Data_->result[iphase]->state;
			rowvec betai = 1 + max(intervalState, 0);
			/////////////////////////////////////////////////////////
			mat diff_of_b = abs(bil_b - bil);
			rowvec sum_diff = sum(diff_of_b);

			uvec less_index = find(sum_diff < (betai*mesh_tol_));
			if (less_index.n_elem == sum_diff.n_elem) return true;

		}// if N_{k}=N_{K-1}

		return false;
	}

	bool LiuHpMeshRefineAlg::CanWeIncreaseN(size_t iphase, size_t segindex, uvec state_index)
	{
		uword istart = state_index(segindex);
		uword iend = state_index(segindex + 1);
		mat segment_state = Data_->result[iphase]->state.rows(istart, iend);
		vec tau = join_vert(Data_->PS[iphase]->Points, ones<vec>(1.0));
		vec tau_segment = tau.subvec(istart, iend);
		mat derive2nd; vec interp_t;
		calculate2nd_derive(tau_segment, segment_state, interp_t, derive2nd);
		mat absderive2nd_M = abs(derive2nd);
		//absderive2nd_M.save("2nd", raw_ascii);
		uword nstate = segment_state.n_cols;
		uvec t_maxindex_in_every_2(nstate);//store the indexes of time of max 2nd derives in each state.
		rowvec Pij_M(nstate);
		for (size_t i = 0; i < nstate; i++)
		{
			Pij_M(i) = absderive2nd_M.col(i).max(t_maxindex_in_every_2(i));
		}
		//t,corresponding to the local maxima 
		vec tmax_in_erery_2 = interp_t(t_maxindex_in_every_2);
		double mintime = tmax_in_erery_2.min();
		double maxtime = tmax_in_erery_2.max();

		//Next, we need to find Pij_(M-1).
		//b means before
		//find the interval that contains mintime and maxtime.
		uvec nodesPerInterval = mesh_history_[mesh_history_.size() - 1][iphase]->nodesPerInterval;
		vec meshpointbefore = mesh_history_[mesh_history_.size() - 1][iphase]->meshpoints;
		uword itmin_b;
		uword itmax_b;

		if (mintime == maxtime)
		{
			if (mintime == tau_segment(0))
			{
				itmin_b = max(find(meshpointbefore <= mintime));
				itmax_b = min(find(meshpointbefore > maxtime));
			}
			else if (maxtime == tau_segment(tau_segment.n_elem - 1))
			{
				itmin_b = max(find(meshpointbefore < mintime));
				itmax_b = min(find(meshpointbefore >= maxtime));
			}
			else
			{
				itmin_b = max(find(meshpointbefore < mintime));
				itmax_b = min(find(meshpointbefore > maxtime));
			}
		}
		else
		{
			itmin_b = max(find(meshpointbefore <= mintime));
			itmax_b = min(find(meshpointbefore >= maxtime));
		}

		uvec indices_b = join_vert(zeros<uvec>(1), cumsum(nodesPerInterval));

		vec time_b = mesh_points_history_[mesh_points_history_.size() - 1][iphase];
		mat state_b = state_history_[state_history_.size() - 1][iphase];

		vec seg_t_b = time_b.subvec(itmin_b, itmax_b);
		mat seg_x_b = state_b.rows(itmin_b, itmax_b);

		mat derive2nd_b; vec interp_t_b;
		calculate2nd_derive(seg_t_b, seg_x_b, interp_t_b, derive2nd_b);
		mat absderive2nd_M_b = abs(derive2nd_b);

		uvec t_maxindex_in_every_2_b(nstate);//store the indexes of time of max 2nd derives in each state.
		rowvec Pij_M_b(nstate);
		for (size_t i = 0; i < nstate; i++)
		{
			Pij_M_b(i) = absderive2nd_M_b.col(i).max(t_maxindex_in_every_2_b(i));
		}

		rowvec Rij = Pij_M / Pij_M_b;
		bool need_dividing_mesh = max(Rij) > ratio_R_;
		return !need_dividing_mesh;
	}

	void LiuHpMeshRefineAlg::calculate2nd_derive(vec t, mat x, vec & interp_t, mat &derive2nd)
	{
		//forward difference
		uword ntime = t.n_rows;
		uword nstates = x.n_cols;
		double tf = t(ntime - 1);
		double t0 = t(0);
		vec tau = 2.0 * (t - t0) / (tf - t0) - 1.0;

		double taustep = 2.0 / 500.0;
		vec tpert = linspace(t0, tf, 501);
		mat xpert(tpert.n_rows, x.n_cols);
		vec temcol;
		for (size_t istate = 0; istate < nstates; istate++)
		{
			SolutionErrorChecker::BarLagrangeInterp(tau, x.col(istate), tpert, temcol);
			xpert.col(istate) = temcol;
		}
		uword npert_t = xpert.n_rows;
		mat xdderiv = xpert.rows(2, npert_t - 1) - 2 * xpert.rows(1, npert_t - 2)
			+ xpert.rows(0, npert_t - 3);
		xdderiv /= (taustep*taustep);
		derive2nd = xdderiv;
		interp_t = tpert.subvec(0, npert_t - 3);
	}


}//namespace lpopc