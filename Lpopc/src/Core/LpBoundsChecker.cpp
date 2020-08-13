// Copyright (C) 2014-2015 Xue Zhichen, Wang Yujie,Wang Na 
// All Rights Reserved.
// This file is a part of LPOPC , published under the Eclipse Public License.
// Author:Xue Zhichen 7/11 2015   18:09
// Email:eddy_lpopc@163.com
#include "LpBoundsChecker.hpp"
#include <strstream>
#include <vector>
namespace Lpopc
{


	void LpBoundsChecker::GetBounds(shared_ptr<OptimalProblem> opt_problem, shared_ptr<LpCalculateData>& calcu_data)
	{
		LP_DBG_START_FUN("LpBoundsChecker::GetBounds()")
		int variable_offset = 0, constraint_offset = 0;
		std::vector<int> nodes(calcu_data->numphases_);
		for (int i = 0; i<nodes.size(); i++){
			nodes[i] = (opt_problem->GetPhase(i)->GetTotalNodes());
		};
		calcu_data->totalnodes_perphase = nodes;

		std::vector<double> var_min;
		std::vector<double> var_max;
		std::vector<double> con_min;
		std::vector<double> con_max;
		std::vector<int> opt_varnum(opt_problem->GetPhaseNum());
		std::vector<int> opt_connum(opt_problem->GetPhaseNum());
		vector<shared_ptr<indices>> opt_phase_indices(opt_problem->GetPhaseNum());
		vector<std::vector<lp_index>> opt_variable_indices(opt_problem->GetPhaseNum());
		vector<std::vector<lp_index>> opt_constraint_indices(opt_problem->GetPhaseNum());
		for (int i = 0; i<calcu_data->numphases_; i++){
			int varnum = 0;//the number of  nlp variable
			int connum = 0;///the number of nlp nonlinear constrains
			shared_ptr<Phase> currentphase(opt_problem->GetPhase(i));
			/*
			Get the lower and upper limits on the initial and terminal
			% time in the current phase
			*/
			double t0_min, tf_min, t0_max, tf_max;
			currentphase->GetTimeMin(t0_min, tf_min);
			currentphase->GetTimeMax(t0_max, tf_max);
			/*Get the lower and upper limits on the states in the current phase*/
			int nodes = (currentphase->GetTotalNodes());
			int statenum = currentphase->GetstateMin().size();


			double minstate0 = 0, minstate = 0, minstatef = 0;
			double maxstate0 = 0, maxstate = 0, maxstatef = 0;

			for (int j = 0; j<currentphase->GetstateMin().size(); j++){
				currentphase->GetstateMin()[j]->GetLimit(minstate0, minstate, minstatef);
				currentphase->GetstateMax()[j]->GetLimit(maxstate0, maxstate, maxstatef);

				bool checkstate0 = (minstate0 <= maxstate0);//make sure statemax>=statemin
				bool checkstate = (minstate <= maxstate);
				bool checkstatef = (minstatef <= maxstatef);
				if (checkstate0&&checkstate&&checkstatef){
					var_min.push_back(minstate0);
					var_max.push_back(maxstate0);
					con_min.push_back(0);
					con_max.push_back(0);
					connum++;
					varnum++;
					int k = 1;
					for (k = 1; k<nodes; k++){
						var_min.push_back(minstate);
						var_max.push_back(maxstate);
						varnum++;
						con_min.push_back(0);
						con_max.push_back(0);
						connum++;
					}
					var_min.push_back(minstatef);
					var_max.push_back(maxstatef);
					varnum++;

				}
				else{
					std::string errmsg1, errmsg2;
					errmsg1 = "Bounds on State are Inconsistent (i.e. max < min) in Phase:";
					std::strstream ss;
					ss << i + 1; ss >> errmsg2;
					LP_THROW_EXCEPTION(LP_BOUNDSCHECKER_ERROR, errmsg1 + errmsg2);
				}

			}


			//Get the lower and upper limits on the controls in the current phase
			double controlmin = 0.0, controlmax = 0.0;
			for (int j = 0; j<currentphase->GetcontrolMin().size(); j++){
				controlmin = currentphase->GetcontrolMin()[j];
				controlmax = currentphase->GetcontrolMax()[j];
				bool checkscontrol = (controlmin <= controlmax);//make sure controlmax>=controlmin
				if (checkscontrol){
					for (int k = 0; k<nodes; k++){
						varnum++;
						var_min.push_back(controlmin);
						var_max.push_back(controlmax);
					}

				}
				else{
					std::string errmsg1, errmsg2;
					errmsg1 = "Bounds on Control are Inconsistent (i.e. max < min) in Phase:";
					std::strstream ss;
					ss << i + 1; ss >> errmsg2;
					LP_THROW_EXCEPTION(LP_BOUNDSCHECKER_ERROR, errmsg1 + errmsg2);
				}//if
			}//end for int j*/
			var_min.push_back(t0_min);
			var_min.push_back(tf_min);
			varnum++;
			varnum++;
			var_max.push_back(t0_max);
			var_max.push_back(tf_max);
			//Get the lower and upper limits on the static parameters in the current phase
			double parametermin = 0.0, parametermax = 0.0;
			for (int j = 0; j<currentphase->GetparameterMin().size(); j++){
				parametermin = currentphase->GetparameterMin()[j];
				parametermax = currentphase->GetparameterMax()[j];
				bool checksparameter = (parametermin <= parametermax);//make sure parametermax>=parametermin
				if (checksparameter){

						var_min.push_back(parametermin);
						var_max.push_back(parametermax);
						varnum++;

				}
				else{
					std::string errmsg1, errmsg2;
					errmsg1 = "Bounds on parameter are Inconsistent (i.e. max < min) in Phase:";
					std::strstream ss;
					ss << i + 1; ss >> errmsg2;
					LP_THROW_EXCEPTION(LP_BOUNDSCHECKER_ERROR, errmsg1 + errmsg2);
				}//if
			}//end for int j*/

			//Get the lower and upper limits on the path constraints in the current phase
			double pathmin = 0.0, pathmax = 0.0;
			for (int j = 0; j<currentphase->GetpathMin().size(); j++){
				pathmin = currentphase->GetpathMin()[j];
				pathmax = currentphase->GetpathMax()[j];
				bool checkspath = (pathmin <= pathmax);//make sure pathmax>=pathmin
				if (checkspath){

					for (int k = 0; k<nodes; k++){
						con_min.push_back(pathmin);
						con_max.push_back(pathmax);
						connum++;
					}

				}
				else{
					std::string errmsg1, errmsg2;
					errmsg1 = "Bounds on path are Inconsistent (i.e. max < min) in Phase:";
					std::strstream ss;
					ss << i + 1; ss >> errmsg2;
					LP_THROW_EXCEPTION(LP_BOUNDSCHECKER_ERROR, errmsg1 + errmsg2);
				}//if
			}//end for int j*/
			//Get the lower and upper limits on the event constraints in the current phase
			double eventmin = 0.0, eventmax = 0.0;
			for (int j = 0; j<currentphase->GeteventMin().size(); j++){
				eventmin = currentphase->GeteventMin()[j];
				eventmax = currentphase->GeteventMax()[j];
				bool checksevent = (eventmin <= eventmax);//make sure eventmax>=eventmin
				if (checksevent){


					con_min.push_back(eventmin);
					con_max.push_back(eventmax);
					connum++;

				}
				else{
					std::string errmsg1, errmsg2;
					errmsg1 = "Bounds on event are Inconsistent (i.e. max < min) in Phase:";
					std::strstream ss;
					ss << i + 1; ss >> errmsg2;
					LP_THROW_EXCEPTION(LP_BOUNDSCHECKER_ERROR, errmsg1 + errmsg2);
				}//if
			}//end for int j*/
			opt_varnum[i]=(varnum);
			opt_connum[i]=(connum);
			
			std::vector<lp_index> var_indices(varnum);
			for (int itor = 0; itor<varnum; itor++){
				var_indices[itor] = variable_offset + itor + 1;
			}
			std::vector<lp_index> con_indices(connum);
			for (int itor = 0; itor<connum; itor++){
				con_indices[itor] = constraint_offset + itor + 1;
			}
			opt_variable_indices[i] = (var_indices);
			opt_constraint_indices[i] = con_indices;

			shared_ptr<indices> var_index(new indices());

			for (int itor = 0; itor<(nodes + 1)*currentphase->GetstateMin().size(); itor++){
				var_index->state.push_back(variable_offset + itor + 1);
			}
			int state_index = var_index->state[var_index->state.size() - 1];
			for (int itor = 0; itor<nodes*currentphase->GetcontrolMin().size(); itor++){
				var_index->control.push_back(state_index + itor + 1);
			}
			int t0_index = 0, tf_index = 0;
			if (var_index->control.size() != 0){
				t0_index = var_index->control[var_index->control.size() - 1] + 1;
			}
			else{
				t0_index = var_index->state[var_index->state.size() - 1] + 1;
			}
			tf_index = t0_index + 1;
			var_index->time.push_back(t0_index);
			var_index->time.push_back(tf_index);
			for (int itor = 0; itor<currentphase->GetparameterMin().size(); itor++){
				var_index->parameter.push_back(tf_index + 1 + itor);
			}
			opt_phase_indices[i] = var_index;

			variable_offset = variable_offset + opt_varnum[i];
			constraint_offset = constraint_offset + opt_connum[i];
		}//for
		//Get linkage bounds
		double linkmin = 0.0;
		double linkmax = 0.0;
		std::vector<std::vector<lp_index>> link_index(calcu_data->numlinkpairs_);
		for (int i = 0; i<calcu_data->numlinkpairs_; i++)
		{
			shared_ptr<Linkage> currentlinkage(opt_problem->GetLinkage(i));

			for (int j = 0; j<currentlinkage->GetLinkageMin().size(); j++){
				linkmin = currentlinkage->GetLinkageMin()[j];
				linkmax = currentlinkage->GetLinkageMax()[j];
				bool checksevent = (linkmin <= linkmax);//make sure eventmax>=eventmin
				if (checksevent){


					con_min.push_back(linkmin);
					con_max.push_back(linkmax);
					link_index[i].push_back(constraint_offset + j + 1);
				}
				else{
					std::string errmsg1, errmsg2;
					errmsg1 = "Bounds on event are Inconsistent (i.e. max < min) in Phase:";
					std::strstream ss;
					ss << i + 1; ss >> errmsg2;
					LP_THROW_EXCEPTION(LP_BOUNDSCHECKER_ERROR, errmsg1 + errmsg2);
				}//if
			}//i
		}
		calcu_data->variables = opt_varnum;
		calcu_data->constraints = opt_connum;
		calcu_data->variable_indices=opt_variable_indices;
		calcu_data->constraint_indices=opt_constraint_indices;
		calcu_data->link_indices = link_index;
		calcu_data->phase_indices = opt_phase_indices;;
		calcu_data->varbounds_min = var_min;
		calcu_data->varbounds_max = var_max;
		calcu_data->conbounds_min = con_min;
		calcu_data->conbounds_max = con_max;

		//Get Bounds on Linear Constraints
		int numvars = calcu_data->varbounds_min.size();

		int alinrowshift = 0;
		vec Alinear_I = zeros(2 * (calcu_data->numphases_ + calcu_data->numlinkpairs_), 1);
		vec Alinear_J = zeros(2 * (calcu_data->numphases_ + calcu_data->numlinkpairs_), 1);
		vec Alinear_V = zeros(2 * (calcu_data->numphases_ + calcu_data->numlinkpairs_), 1);

		vec Alinmin = zeros(calcu_data->numphases_ + calcu_data->numlinkpairs_, 1);
		vec Alinmax = zeros(calcu_data->numphases_ + calcu_data->numlinkpairs_, 1);
		//Part 1:  Monotonicity of Independent Variable
		for (int i = 0; i<calcu_data->numphases_; i++){
			shared_ptr<Phase> currentphase(opt_problem->GetPhase(i));
			int nstate = currentphase->GetstateMin().size();
			int ncontrol = currentphase->GetcontrolMin().size();
			int ishift = 0;
			if (i != 0){
				ishift = calcu_data->variable_indices[i - 1][calcu_data->variable_indices[i - 1].size() - 1];
			}

			int t0_index = ishift + ((nodes[i] + 1)*currentphase->GetstateMin().size())
				+ nodes[i] * currentphase->GetcontrolMin().size() + 1;
			int tf_index = t0_index + 1;
			Alinear_I(alinrowshift) = i;
			Alinear_J(alinrowshift) = t0_index - 1;
			Alinear_V(alinrowshift) = -1;
			alinrowshift++;
			Alinear_I(alinrowshift) = i;
			Alinear_J(alinrowshift) = tf_index - 1;
			Alinear_V(alinrowshift) = 1;
			alinrowshift++;
			// Check if only ONE of the lower and upper bounds on the phase duration are specified
			if (currentphase->HasDuration()){
				double duramin, duramax;
				currentphase->Getduration(duramin, duramax);
				if (duramin <= duramax){
					Alinmin(i) = duramin;
					Alinmax(i) = duramax;
				}
				else{
					std::string errmsg1, errmsg2;
					errmsg1 = "Bounds on duration are Inconsistent (i.e. max < min) in Phase:";
					std::strstream ss;
					ss << i + 1; ss >> errmsg2;
					LP_THROW_EXCEPTION(LP_BOUNDSCHECKER_ERROR, errmsg1 + errmsg2);
				}
			}
			else{
				Alinmin(i) = 0;
				Alinmax(i) = datum::inf;
			}

		}//for int i

		int istart = calcu_data->numphases_;
		//Part 2:  Linkage of Time Across Phases
		for (int i = 0; i<calcu_data->numlinkpairs_; i++){
			shared_ptr<Linkage> currentlinkage(opt_problem->GetLinkage(i));

			int left_phase = currentlinkage->LeftPhase();
			int right_phase = currentlinkage->RightPhase();
			int nparameters_left = calcu_data->SIZES_[left_phase][2];
			int nparameters_right = calcu_data->SIZES_[right_phase][2];

			int tf_index_left = calcu_data->variable_indices[left_phase][calcu_data->variable_indices[left_phase].size() - nparameters_left - 1];
			int t0_index_right = calcu_data->variable_indices[right_phase][calcu_data->variable_indices[right_phase].size() - nparameters_right - 2];
			Alinear_I(alinrowshift) = istart + i;
			Alinear_J(alinrowshift) = tf_index_left - 1;
			Alinear_V(alinrowshift) = -1;
			alinrowshift++;
			Alinear_I(alinrowshift) = istart + i;
			Alinear_J(alinrowshift) = t0_index_right - 1;
			Alinear_V(alinrowshift) = 1;
			alinrowshift++;
			Alinmin(istart + i) = 0;
			Alinmax(istart + i) = 0;
		}
		
		calcu_data->linmin = Alinmin;
		calcu_data->linmax = Alinmax;
		calcu_data->AlinearMatrix = dsmatrix::Sparse(Alinear_I, Alinear_J, Alinear_V,
		calcu_data->numphases_ + calcu_data->numlinkpairs_, numvars);

	}

}//end of namespace  Lpopc