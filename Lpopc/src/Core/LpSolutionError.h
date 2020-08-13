// Copyright (C) 2014-2015 Xue Zhichen, Wang Yujie,Wang Na 
// All Rights Reserved.
// This file is a part of LPOPC , published under the Eclipse Public License.
// Author:Xue Zhichen 7/11 2015   22:57
// Email:eddy_lpopc@163.com
#ifndef LPOPC_SOLUTION_ERROR_CHECKER_H
#define LPOPC_SOLUTION_ERROR_CHECKER_H
#include "LpConf.h"
#include "LpException.hpp"
#include "LpOptimalProblem.hpp"
#include "LpFunctionWrapper.h"
#include"LpCalculateData.hpp"
namespace Lpopc
{
	LP_DECLARE_EXCEPTION(LPOPC_SOLUTION_ERROR);
	class SolutionErrorChecker
	{
	public:
		SolutionErrorChecker(shared_ptr<FunctionWrapper> Funs, shared_ptr<LpCalculateData> Data, shared_ptr<OptimalProblem> optpro)
			: Funs_(Funs), Data_(Data), optpro_(optpro){};
		~SolutionErrorChecker(){};
		///Calculate error in differential equations,using one more points
		void CheckSolutionDiffError( int iphase,mat& relative_error)const ;
		///Same with above one,return state interpolated  with one more points,Liu-hp-methods needs.
		void CheckSolutionDiffError(int iphase, mat& relative_error, mat& temstate)const;
		/************************************************************************/
		/* % Interpolates the given data using the Barycentric
		Lagrange Interpolation formula. Vectorized to remove all loops

		data_x contains the
        nodes and data_y column two contains the function value
        at the nodes
		y - interpolated data. 

		Reference:
		(1) Jean-Paul Berrut & Lloyd N. Trefethen, "Barycentric Lagrange
		Interpolation"
		 http://web.comlab.ox.ac.uk/oucl/work/nick.trefethen/berrut.ps.gz
		(2) Walter Gaustschi, "Numerical Analysis, An Introduction" (1997) pp. 94-95                                                                     */
		/************************************************************************/
		static void BarLagrangeInterp(const vec& data_x, const vec& data_y, const vec& x, vec& y);
	private:
		SolutionErrorChecker(const SolutionErrorChecker&);
		SolutionErrorChecker& operator=(const SolutionErrorChecker&);
		void SolutionInterpolation(const int iphase, vec& vec_time, mat& mat_state, mat& mat_control)const;
		shared_ptr<FunctionWrapper> Funs_;
		shared_ptr<LpCalculateData> Data_;
		shared_ptr<OptimalProblem> optpro_;
	};
}//namespace Lpopc
#endif // !LPOPC_SOLUTION_ERROR_CHECKER_H
