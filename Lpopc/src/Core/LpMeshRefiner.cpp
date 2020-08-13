// Copyright (C) 2014-2015 Xue Zhichen, Wang Yujie,Wang Na 
// All Rights Reserved.
// This file is a part of LPOPC , published under the Eclipse Public License.
// Author:Xue Zhichen 7/11 2015   22:57
// Email:eddy_lpopc@163.com
#include "LpMeshRefiner.h"
#include <math.h>
namespace Lpopc
{ 
	void Lpopc::MeshRefiner::SetAndCheckMesh(shared_ptr<OptimalProblem>& optpro)
	{
		LP_DBG_START_FUN("MeshRefiner::SetAndCheckMesh()")
		for (int i = 0; i < optpro->GetPhaseNum(); i++){
			shared_ptr<Phase> temPhase(optpro->GetPhase(i));
			std::vector<double> meshpoints = temPhase->GetMeshPoints();
			std::vector<lp_index> nodesperinterval = temPhase->GetNodesPerInterval();
			if (meshpoints.size() == 0){
				reporter_->Printf()->info("Start with default mesh points!");
				if (nodesperinterval.size() != 0)
				{
					vec temVec = linspace(-1, 1, nodesperinterval.size() + 1);
					for (int j = 0; j < temVec.n_elem; j++){
						temPhase->GetMeshPoints().push_back(temVec(i));
					};
				}
				else{
					temPhase->GetMeshPoints().push_back(-1);
					temPhase->GetMeshPoints().push_back(1);
				};

			}
			else if (meshpoints.size() == 1){
				std::string errmsg1, errmsg2;
				errmsg1 = "MeshRefinement need at least two  meshPoints,but there's only one in phase";
				std::strstream ss;
				ss << i + 1; ss >> errmsg2;
				LP_THROW_EXCEPTION(MESHREFINE_ERROR, errmsg1 + errmsg2);

			}
			else if ((meshpoints[0] != -1) || (meshpoints[meshpoints.size() - 1] != 1)){
				std::string errmsg1, errmsg2;
				errmsg1 = "meshPoints must span -1 to +1 in phase";
				std::strstream ss;
				ss << i + 1; ss >> errmsg2;
				LP_THROW_EXCEPTION(MESHREFINE_ERROR, errmsg1 + errmsg2);
			};
			if (nodesperinterval.size() == 0){
				reporter_->Printf()->info("Start with default nodes per interval (10)!");
				for (int k = 0; k < (temPhase->GetMeshPoints().size() - 1); k++){
					temPhase->GetNodesPerInterval().push_back(20);
				}
			}
			else if (meshpoints.size() != nodesperinterval.size() + 1)
			{
				std::string errmsg2, errmsg1 = "Number of nodesPerInterval must match number of mesh intervals in phase";
				std::strstream ss;
				ss << i + 1; ss >> errmsg2;
				LP_THROW_EXCEPTION(MESHREFINE_ERROR, errmsg1 + errmsg2);
			}
		}//for end
	}


	bool MeshRefiner::RefineMesh(shared_ptr<OptimalProblem> optpro)
	{
		LP_DBG_START_FUN("MeshRefiner::RefineMesh()")
			if (grid_ > gridmaxnum_)
			{
				LP_THROW_EXCEPTION(MESHREFINE_ERROR, "The problem reach the max number of refine grid,but hasn't reach the derized error tolrance");
			}
		if (grid_ == 0)
		{
			std::vector<shared_ptr<LpMesh>> firstLpMesh(optpro->GetPhaseNum());
			for (size_t iphase = 0; iphase < optpro->GetPhaseNum(); iphase++)
			{
				firstLpMesh[iphase].reset(new LpMesh());
				firstLpMesh[iphase]->meshpoints = optpro->GetPhase(iphase)->GetMeshPoints();
				firstLpMesh[iphase]->nodesPerInterval = optpro->GetPhase(iphase)->GetNodesPerInterval();
			}
			meshhistory.push_back(firstLpMesh);
		}
		std::vector<shared_ptr<LpMesh>> newmesh;
		reporter_->Printf()->info("Start refine grid") << grid_ + 1;
		bool no_more_refine = meshrefine_alg_->RefineMesh(optpro, newmesh);
		if (!no_more_refine)
		{
			grid_++;
			meshhistory.push_back(newmesh);
		}

		return no_more_refine;
	}

}

