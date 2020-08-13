// Copyright (C) 2014-2015 Xue Zhichen, Wang Yujie,Wang Na 
// All Rights Reserved.
// This file is a part of LPOPC , published under the Eclipse Public License.
// Author:Xue Zhichen 7/11 2015   15:36
// Email:eddy_lpopc@163.com

#ifndef LPOPC_ALGORITHM_HPP
#define LPOPC_ALGORITHM_HPP
#include "LpConf.h"
#include "LpDebug.hpp"
#include "LpOptionList.hpp"
#include "LpNLPSolver.h"
#include "LpNLPWrapper.hpp"
#include "LpAnalyticDerive.hpp"
#include "LpANDeriveChecker.h"
#include "LpException.hpp"
#include "LpMeshRefiner.h"
#include "LpCalculateData.hpp"
#include "LpOption.hpp"
#include "LpSizeChecker.h"
#include "LpBoundsChecker.hpp"
#include "LpGuessChecker.h"
#include "RPMGenerator.hpp"
#include "LpSacleOCP.hpp"
#include "Nlp2OPConverter.h"
#include "LpDerivDependciesChecker.h"
#include "LpFunctionWrapper.h"
#include <string>
#include <vector>
namespace Lpopc
{
	LP_DECLARE_EXCEPTION(LPOPC_ALGORITHM_ERROR);
	enum ocpstatus
	{
		solved=0,
		reach_max_grid,
		nlp_failed,
		matrix_operation_failed,
		lpopc_failed
	};
	class LpopcAlgorithm
	{
	public:
		LpopcAlgorithm(shared_ptr<OptimalProblem>optpro, shared_ptr<LpCalculateData> cd_data,
			shared_ptr<LpOptionsList> optionlist, 
			shared_ptr<LpReporter>    msg_reporter);
		~LpopcAlgorithm();;

		void Initialized();
		void SetFirstMesh();;
		//do
		void GetSizes();;
		void GetBounds();;
		void GetGuess();
		void GetDependecies();
		void GetScales();
		void SolveNlp();
		void Nlp2OpControl();
		bool RefineMesh();
		//while;
		ocpstatus SolveOptimalControlProblem();
		shared_ptr<LpOptionsList>Options();
	private:
		LpopcAlgorithm();
		LpopcAlgorithm(const LpopcAlgorithm&);
		LpopcAlgorithm& operator=(const LpopcAlgorithm&);
		void UpdateGrid();
		void OutputProblemInfo();
		shared_ptr<LpReporter>    msg_reporter_;
		shared_ptr<LpOptionsList> optionlist_;
		shared_ptr<LpCalculateData> cd_data_;
		shared_ptr<LpSizeChecker> size_checker_;
		shared_ptr<LpBoundsChecker> bounds_checker_;
		shared_ptr<LpGuessChecker> guess_checker_;
		shared_ptr<OptimalProblem> optpro_;
		shared_ptr<MeshRefiner> meshrefiner_;
		shared_ptr<RPMGenerator> rpm_;
		shared_ptr<LpScaleOCP>   ocp_scale_;
		shared_ptr<NLPSolver> nlp_splver_;
		shared_ptr<NLPWrapper> nlp_wrapper_;
		shared_ptr<Nlp2OpConverter> nlp2op_converter_;
		shared_ptr<DeriveDependicieshecker> dependicies_checker_;
		shared_ptr<FunctionWrapper> funs_;
		shared_ptr<OptDerive> first_derive_;
		shared_ptr<LpHessianCalculator> hessian_;
        
		shared_ptr<LpANDeriveChecker> an_derive_checker_;
	};
}//namespace Lpopc
#endif
