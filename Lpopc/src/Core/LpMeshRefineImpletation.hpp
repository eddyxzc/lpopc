// Copyright (C) 2014-2015 Xue Zhichen, Wang Yujie,Wang Na 
// All Rights Reserved.
// This file is a part of LPOPC , published under the Eclipse Public License.
// Author:Xue Zhichen 8/4 2015   21:55
// Email:eddy_lpopc@163.com
#ifndef LPOPC_MESHREFINE_IMPLETATION_HPP
#define LPOPC_MESHREFINE_IMPLETATION_HPP
#include "LpOptimalProblem.hpp"
#include "LpConf.h"
namespace Lpopc
{
	struct LpMesh
	{
		uvec nodesPerInterval;
		vec meshpoints;
		vec e_k;//only used in Liu-hp-method
	};
	class MeshRefineImpl
	{
	public:
		MeshRefineImpl(){};
		virtual ~MeshRefineImpl(){};
		virtual bool RefineMesh(shared_ptr<OptimalProblem>optpro,std::vector<shared_ptr<LpMesh>>& tempMeshVector) = 0 ;
	private:
		MeshRefineImpl(const MeshRefineImpl&);
		MeshRefineImpl& operator=(const MeshRefineImpl&);
	};
}
#endif // !LPOPC_MESHREFINE_IMPLETATION_HPP
