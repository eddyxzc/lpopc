// Copyright (C) 2014-2015 Xue Zhichen, Wang Yujie,Wang Na 
// All Rights Reserved.
// This file is a part of LPOPC , published under the Eclipse Public License.
// Author:Xue Zhichen 7/11 2015   0:34
// Email:eddy_lpopc@163.com
#include "LpSizeChecker.h"
#include <strstream>
#include <vector>
namespace Lpopc
{


	void LpSizeChecker::GetSize(shared_ptr<OptimalProblem> opt_problem, shared_ptr<LpCalculateData>& calcu_data)
	{
		std::vector<std::vector<lp_index>> opt_size(opt_problem->GetPhaseNum());
		for (int i = 0; i < opt_problem->GetPhaseNum(); i++){
			int phasetotalnodes = 0;
			for (int interval = 0; interval < opt_problem->GetPhase(i)->GetNodesPerInterval().size(); interval++){

				phasetotalnodes += opt_problem->GetPhase(i)->GetNodesPerInterval()[interval];
			}
			opt_problem->GetPhase(i)->SetTotalNodes(phasetotalnodes);
			// Determine the number of states in the current phase 
			std::vector<shared_ptr<Limit>> statemin = opt_problem->GetPhase(i)->GetstateMin();
			std::vector<shared_ptr<Limit>> statemax = opt_problem->GetPhase(i)->GetstateMax();
			int nstate = 0;
			if ((statemin.size() != 0) && (statemax.size() != 0)){
				if (statemin.size() == statemax.size()){
					nstate = statemin.size();
				}
				else{
					std::string errmsg1, errmsg2;
					errmsg1 = "State upper & lower bound  MUST be same size in phase";
					std::strstream ss;
					ss << i + 1; ss >> errmsg2;
					LP_THROW_EXCEPTION(LP_SIZECHECK_ERROR, errmsg1 + errmsg2);
				}
			}
			else{
				nstate = 0;
			}

			// Determine the number of controls in the current phase
			std::vector<double> controlmin = opt_problem->GetPhase(i)->GetcontrolMin();
			std::vector<double> controlmax = opt_problem->GetPhase(i)->GetcontrolMax();
			int ncontrol = 0;
			if ((controlmin.size() != 0) && (controlmax.size() != 0)){
				if (controlmin.size() == controlmax.size()){
					ncontrol = controlmin.size();
				}
				else{
					std::string errmsg1, errmsg2;
					errmsg1 = "Control upper & lower bound  MUST be same size in phase";
					std::strstream ss;
					ss << i + 1; ss >> errmsg2;
					LP_THROW_EXCEPTION(LP_SIZECHECK_ERROR, errmsg1 + errmsg2);
				}
			}
			else{
				ncontrol = 0;
			}

			//Determine the number of static parameters in the current phase
			std::vector<double> parametermin = opt_problem->GetPhase(i)->GetparameterMin();
			std::vector<double> parametermax = opt_problem->GetPhase(i)->GetparameterMax();
			int nparameter = 0;
			if ((parametermin.size() != 0) && (parametermax.size() != 0)){
				if (parametermin.size() == parametermax.size()){
					nparameter = parametermin.size();
				}
				else{
					std::string errmsg1, errmsg2;
					errmsg1 = "Parameter upper & lower bound  MUST be same size in phase";
					std::strstream ss;
					ss << i + 1; ss >> errmsg2;
					LP_THROW_EXCEPTION(LP_SIZECHECK_ERROR, errmsg1 + errmsg2);
				}
			}
			else{
				nparameter = 0;
			}
			// Number of path constraints in the current phase
			std::vector<double> pathmin = opt_problem->GetPhase(i)->GetpathMin();
			std::vector<double> pathmax = opt_problem->GetPhase(i)->GetpathMax();
			int npath = 0;
			if ((pathmin.size() != 0) && (pathmax.size() != 0)){
				if (pathmin.size() == pathmax.size()){
					npath = pathmin.size();
				}
				else{
					std::string errmsg1, errmsg2;
					errmsg1 = "Path upper & lower bound  MUST be same size in phase";
					std::strstream ss;
					ss << i + 1; ss >> errmsg2;
					LP_THROW_EXCEPTION(LP_SIZECHECK_ERROR, errmsg1 + errmsg2);
				}
			}
			else{
				npath = 0;
			}
			//Number of event constraints in the current phase
			std::vector<double> eventmin = opt_problem->GetPhase(i)->GeteventMin();
			std::vector<double> eventmax = opt_problem->GetPhase(i)->GeteventMax();
			int nevent = 0;
			if ((eventmin.size() != 0) && (eventmax.size() != 0)){
				if (eventmin.size() == eventmax.size()){
					nevent = eventmin.size();
				}
				else{
					std::string errmsg1, errmsg2;
					errmsg1 = "Event upper & lower bound  MUST be same size in phase";
					std::strstream ss;
					ss << i + 1; ss >> errmsg2;
					LP_THROW_EXCEPTION(LP_SIZECHECK_ERROR, errmsg1 + errmsg2);
				}
			}
			else{
				nevent = 0;
			}
			std::vector<lp_index> iphase_size(5);
			iphase_size[0] = nstate;
			iphase_size[1] = ncontrol;
			iphase_size[2] = nparameter;
			iphase_size[3] = npath;
			iphase_size[4] = nevent;
			opt_size[i]=iphase_size;

		}//end for
		calcu_data->SIZES_= opt_size;

		int numlinks = 0;
		for (int i = 0; i < opt_problem->GetLinkageNum(); i++){
			//Check the sizes of the lower and upper limits
			std::vector<double> linkagemin = opt_problem->GetLinkage(i)->GetLinkageMin();
			std::vector<double> linkagemax = opt_problem->GetLinkage(i)->GetLinkageMax();
			if (linkagemin.size() != linkagemax.size())
			{
				std::string errmsg1, errmsg2;
				errmsg1 = "'Linkage Upper and Lower Bound Vector Must Be Same Size in linkage";
				std::strstream ss;
				ss << i + 1; ss >> errmsg2;
				LP_THROW_EXCEPTION(LP_SIZECHECK_ERROR, errmsg1 + errmsg2);
			}
			else{
				numlinks = numlinks + linkagemin.size();

			}
		}//end for
		calcu_data->numphases_ = opt_problem->GetPhaseNum();
		calcu_data->numlinkpairs_ = opt_problem->GetLinkageNum();
		calcu_data->numlinks_ = numlinks;
	}

}//end of namespace Lpopc