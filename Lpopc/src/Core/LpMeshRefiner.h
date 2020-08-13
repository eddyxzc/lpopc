// Copyright (C) 2014-2015 Xue Zhichen, Wang Yujie,Wang Na 
// All Rights Reserved.
// This file is a part of LPOPC , published under the Eclipse Public License.
// Author:Xue Zhichen 7/11 2015   22:57
// Email:eddy_lpopc@163.com
#ifndef LPMESHREFINER_H
#define LPMESHREFINER_H
#include "LpOptimalProblem.hpp"
#include "LpException.hpp"
#include"LpReporter.hpp"
#include "LpCalculateData.hpp"
#include "LpSolutionError.h"
#include "LpOption.hpp"
#include "LpOptionList.hpp"
#include "LpMeshRefineImpletation.hpp"
#include "LpPhMeshRefineAlg.hpp"
#include "LpLiuHpMeshRefineAlg.hpp"
#include "LpDebug.hpp"
#include <strstream> 
#include<vector>
#include <memory>
namespace Lpopc{

	LP_DECLARE_EXCEPTION(MESHREFINE_ERROR);
	

	class MeshRefiner{
		
	public:
		 MeshRefiner(shared_ptr<LpReporter> repoter,shared_ptr<LpOptionsList> optionlist)
		 :reporter_(repoter),optionlist_(optionlist){}
		virtual ~MeshRefiner(){};
		
		/*return true if all error satisfied,false in contrast */
		bool RefineMesh(shared_ptr<OptimalProblem> optpro);
		void Initialize(shared_ptr<FunctionWrapper> Funs, shared_ptr<LpCalculateData> Data,
			shared_ptr<OptimalProblem> optpro)
		{
			LP_DBG_START_FUN("MeshRefiner::Initialize()")
			std::string mesh_refine_alg;
			optionlist_->GetNumericValue("desired-relative-error", tol_, "");
			optionlist_->GetStringValue("mesh-refine-methods", mesh_refine_alg, "");
			 grid_ = 0;
			if (mesh_refine_alg == "ph")
			{
				int Nmin, Nmax;
				optionlist_->GetIntegerValue("Nmax", Nmax, "");
				optionlist_->GetIntegerValue("Nmin", Nmin, "");
				optionlist_->GetIntegerValue("max-grid-num", gridmaxnum_, "");
// 				shared_ptr<PhMeshRefineAlg> ph_refine(new PhMeshRefineAlg(Nmax, Nmin, tol_, Funs, Data, reporter_));
// 				// ph_refine = std::make_shared<PhMeshRefineAlg>();
				meshrefine_alg_ = shared_ptr<PhMeshRefineAlg>(new PhMeshRefineAlg(Nmax, Nmin, tol_, Funs, Data, reporter_));
			}
			else if (mesh_refine_alg=="hp-Liu")
			{
				int  Nmax; double R;
				optionlist_->GetIntegerValue("Nmax", Nmax, "");
				optionlist_->GetNumericValue("R", R, "");
				optionlist_->GetIntegerValue("max-grid-num", gridmaxnum_, "");
				meshrefine_alg_ = std::make_shared<LiuHpMeshRefineAlg>(Nmax, R, tol_, Funs, Data, reporter_);
			}
			else
			{
				LP_THROW_EXCEPTION(MESHREFINE_ERROR,"Unknown Refine Methods")
			}
		};
		static void SetAllOptions(shared_ptr<LpRegisteredOptions>& registeredOptions)
		{
			registeredOptions->AddStringOption2("mesh-refine-methods", "mesh refine methods", "ph",
				"ph", "ph methods",
				"hp-Liu", "hp-Liu methods");
			registeredOptions->AddIntegerOption("max-grid-num", "the max times of grid refinement", 10);
			registeredOptions->AddNumberOption("desired-relative-error", "the desired relative error in mesh refinement.", 1e-6, 
				"the desired relative error in mesh refinement.It should be equal of a litter bigger than the tolerance of finite difference");

			registeredOptions->AddIntegerOption("Nmax", "max number of LGL points", 16);
			registeredOptions->AddIntegerOption("Nmin", "min number of LGL points", 4);
			registeredOptions->AddNumberOption("R", "R in hp-Liu second derive",1.2);

		}
		 void SetAndCheckMesh(shared_ptr<OptimalProblem>& optpro);///Set the first grid 
		  inline size_t CurrentGrid()const{ return grid_; };
	private:
		void operator=(const MeshRefiner&);
		MeshRefiner(const MeshRefiner&);

		bool reachgridMax(){ return grid_ >= gridmaxnum_; };
		double tol_;
		int gridmaxnum_;
		size_t grid_;
		shared_ptr<LpReporter> reporter_;
		shared_ptr<MeshRefineImpl> meshrefine_alg_;
		shared_ptr<LpOptionsList> optionlist_;
		std::vector<std::vector<shared_ptr<LpMesh>>> meshhistory;
		
	};//MeshRefiner
}//Lpopc
#endif // !LPMESHREFINER_H


