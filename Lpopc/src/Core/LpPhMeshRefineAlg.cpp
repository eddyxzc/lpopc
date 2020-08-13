// Copyright (C) 2014-2015 Xue Zhichen, Wang Yujie,Wang Na 
// All Rights Reserved.
// This file is a part of LPOPC , published under the Eclipse Public License.
// Author:Xue Zhichen 8/4 2015   21:26
// Email:eddy_lpopc@163.com
#include "LpPhMeshRefineAlg.hpp"
#include "LpSolutionError.h"
namespace Lpopc
{


	bool PhMeshRefineAlg::RefineMesh(shared_ptr<OptimalProblem>optpro, std::vector<shared_ptr<LpMesh>>& allMeshVector) 
	{
		size_t nphases = optpro->GetPhaseNum();
		shared_ptr<SolutionErrorChecker> error_checker(new SolutionErrorChecker(Funs_, Data_, optpro));
		bool NoMoreRefine = true;
		allMeshVector = std::vector<shared_ptr<LpMesh>>(nphases);
		for (size_t iphase = 0; iphase < nphases; iphase++)
		{
			mat relative_error;
			error_checker->CheckSolutionDiffError(iphase, relative_error);
			//relative_error.save("relative_error", raw_ascii);
			size_t numseg = optpro->GetPhase(iphase)->GetNodesPerInterval().size();
			uvec nodesPerInterval(optpro->GetPhase(iphase)->GetNodesPerInterval());
			nodesPerInterval += 1;

			uvec indices = join_vert(zeros<uvec>(1), cumsum(nodesPerInterval));
			std::vector<shared_ptr<LpMesh>> tempMeshVector(numseg);
			std::vector<lp_index> newnodesPerInterval;
			std::vector<double> newmeshpoints;
			newmeshpoints.push_back(-1);
			for (size_t iseg = 0; iseg < numseg; iseg++)
			{
				size_t istart = indices(iseg);
				size_t ifinish = indices(iseg + 1);
				mat intervalError = relative_error.rows(istart, ifinish);
				double maxtem = intervalError.max();
				bool errorsatisfied = intervalError.max() <= tol_;

				if (errorsatisfied)
				{
					//all satisfied
					reporter_->Printf()->info( "Error Tolerance Satisfied in Segment {} of Phase {}\n", iseg + 1, iphase + 1);
					double  meshpoint0 = optpro->GetPhase(iphase)->GetMeshPoints()[iseg];
					double  meshpointf = optpro->GetPhase(iphase)->GetMeshPoints()[iseg + 1];
					tempMeshVector[iseg].reset(new LpMesh());
					tempMeshVector[iseg]->meshpoints = linspace(meshpoint0, meshpointf, 2);
					tempMeshVector[iseg]->nodesPerInterval = optpro->GetPhase(iphase)->GetNodesPerInterval()[iseg];
				}
				else
				{
					// not all satisfied ,modify segment 
					reporter_->Printf()->info ("Error Tolerance Not Satisfied in Segment {} of Phase {}\n", iseg + 1, iphase + 1);
					tempMeshVector[iseg].reset(new LpMesh());
					ModifySegment(optpro, iphase, iseg, intervalError.max(), tol_, tempMeshVector[iseg]);
					NoMoreRefine = false;
				}
				for (size_t i = 1; i < tempMeshVector[iseg]->meshpoints.n_elem; i++)
				{
					newmeshpoints.push_back(tempMeshVector[iseg]->meshpoints(i));
				}
				for (size_t i = 0; i < tempMeshVector[iseg]->nodesPerInterval.n_elem; i++)
				{
					newnodesPerInterval.push_back(tempMeshVector[iseg]->nodesPerInterval(i));
				}

			}//for iseg
			optpro->GetPhase(iphase)->GetMeshPoints() = newmeshpoints;
			optpro->GetPhase(iphase)->GetNodesPerInterval() = newnodesPerInterval;
			allMeshVector[iphase].reset(new LpMesh());
			allMeshVector[iphase]->meshpoints = newmeshpoints;
			allMeshVector[iphase]->nodesPerInterval = newnodesPerInterval;
		}//for iphase

		return NoMoreRefine;
	}

	void PhMeshRefineAlg::ModifySegment(shared_ptr<OptimalProblem> optpro, int iphase, int segindex, double emax, double tol, shared_ptr<LpMesh> newmesh)const
	{
		int nodes_current = optpro->GetPhase(iphase)->GetNodesPerInterval()[segindex];
		int  Pq = static_cast<int>(log(emax / tol) / log(nodes_current));
		int newnodes = nodes_current + Pq;
		double  meshpoint0 = optpro->GetPhase(iphase)->GetMeshPoints()[segindex];
		double  meshpointf = optpro->GetPhase(iphase)->GetMeshPoints()[segindex + 1];
		if (newnodes <= Nmax_)
		{
			// p methods :only increase the number of collacation points in this mesh
			newmesh->nodesPerInterval = newnodes*ones<uvec>(1);
			newmesh->meshpoints = linspace(meshpoint0, meshpointf, 2);
		}
		else
		{
			//N methods :this mesh interval is divided into Bq subintervals,each newly created subinterval
			//will contain Nmin colllocation points
			int Bq = static_cast<int>(std::max(ceil(double(newnodes) / double(Nmin_)), 2.0));
			newmesh->nodesPerInterval = Nmin_*ones<uvec>(Bq);
			newmesh->meshpoints = linspace(meshpoint0, meshpointf, Bq + 1);
		}
	}

}// namespace Lpopc