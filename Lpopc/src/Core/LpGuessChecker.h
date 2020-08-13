// Copyright (C) 2014-2015 Xue Zhichen, Wang Yujie,Wang Na 
// All Rights Reserved.
// This file is a part of LPOPC , published under the Eclipse Public License.
// Author:Xue Zhichen 7/11 2015   18:09
// Email:eddy_lpopc@163.com
#ifndef LPOPC_GUESS_CHECKER_H
#define LPOPC_GUESS_CHECKER_H
#include "LpOptimalProblem.hpp"
#include "LpOptimalProblem.hpp"
#include "LpCalculateData.hpp"
#include "RPMGenerator.hpp"
#include "LpConf.h"
#include "LpDebug.hpp"
namespace Lpopc
{
	LP_DECLARE_EXCEPTION(LP_GUESSCHECKER_ERROR);
	LP_DECLARE_EXCEPTION(SPLINE_ERROR);
	struct PhaseGuess
	{
		vec tauGuess;
		mat xtauGuess;
		mat utauGuess;
		vec ptauGuess;
	};
	class LpGuessChecker
	{
	public:
		LpGuessChecker(){}
		void GetGuess(shared_ptr<OptimalProblem> opt_problem, shared_ptr<RPMGenerator>& rpm, shared_ptr<LpCalculateData>& calcul_data);
		~LpGuessChecker(){}
        //   Given the arrays xdata[i] and ydata[i], i = 0,...n-1, which tabulate a function, with xdata[i] < xdata[i+1],
		//   and given the array d2y[i], which is the output from function spline_second_derivative(),
		//   and given a value of x, this function returns the interpolated value y using (natural) cubic-spline interpolation
		//
		//   Reference: Burden and Faires (2005) "Numerical Analysis". Thompson.
		static void spline_interpolation(double* y, double& x, double* xdata, double* ydata, int n);
		static void spline_interpolation(double* y, double x, vec& xdata, vec& ydata);
		static void spline_second_derivative(double *x, double *y, int n, double *d2y);
		static void spline_interpolation(vec& y, vec& x, vec& xdata, vec& ydata);
	private:
		LpGuessChecker(const LpGuessChecker&);
		LpGuessChecker& operator= (const LpGuessChecker&);

		
	};

}// namespace Lpopc
#endif // !LPOPC_GUESS_CHECKER_H
