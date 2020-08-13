// Copyright (C) 2014-2015 Xue Zhichen, Wang Yujie,Wang Na 
// All Rights Reserved.
// This file is a part of LPOPC , published under the Eclipse Public License.
// Author:Xue Zhichen 8/4 2015   21:25
// Email:eddy_lpopc@163.com
#ifndef LPOPC_PH_MESHREFINE_ALG_HPP
#define LPOPC_PH_MESHREFINE_ALG_HPP

#include "LpConf.h"
#include "LpMeshRefineImpletation.hpp"
#include "LpNLPWrapper.hpp"
#include "LpCalculateData.hpp"
namespace Lpopc
{
	/**This class implement  adaptive-ph mesh refinement
	*/
	class PhMeshRefineAlg :public MeshRefineImpl
	{
	public:
		PhMeshRefineAlg(int Nmax,int Nmin,double tol, shared_ptr<FunctionWrapper> Funs, 
			 shared_ptr<LpCalculateData> Data,
			 shared_ptr<LpReporter>reporter)
			 :Nmax_(Nmax), Nmin_(Nmin), tol_(tol),
			 Funs_(Funs), Data_(Data), reporter_(reporter)
		{}
		virtual ~PhMeshRefineAlg(){}
		virtual bool RefineMesh(shared_ptr<OptimalProblem>optpro, std::vector<shared_ptr<LpMesh>>& tempMeshVector) ;
		///tempMeshVector contains new mesh information
	private:
 		PhMeshRefineAlg(const PhMeshRefineAlg&);
 		void operator=(const PhMeshRefineAlg&);
		void ModifySegment(shared_ptr<OptimalProblem> optpro, int iphase, int segindex, double emax, double tol, shared_ptr<LpMesh> newmesh)const;

		double tol_;
		int Nmin_;
		int Nmax_;
		shared_ptr<FunctionWrapper> Funs_;
		shared_ptr<LpCalculateData> Data_;
		shared_ptr<LpReporter> reporter_;
	};
}
#endif // !LPOPC_PH_MESHREFINE_ALG_HPP

