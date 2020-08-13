// Copyright (C) 2014-2015 Xue Zhichen, Wang Yujie,Wang Na 
// All Rights Reserved.
// This file is a part of LPOPC , published under the Eclipse Public License.
// Author:Xue Zhichen 7/11 2015   0:07
// Email:eddy_lpopc@163.com
#ifndef LPCALCULATEDATA_HPP_
#define LPCALCULATEDATA_HPP_
#include "LpConf.h"
#include "LpSparseMatrix.h"
#include "RPMGenerator.hpp"
#include "LpReporter.hpp"
#include <vector>
namespace Lpopc
{
	using  std::vector;
	struct SolutionData
	{
		vec time;
		mat state;
		mat control;
		mat parameter;
		mat costate;
		mat pathmult;
		mat Hamiltonian;
		double mayerCost;
		double lagrangeCost;
		
	};
	struct indices{
		std::vector<lp_index> state;
		std::vector<lp_index> control;
		std::vector<lp_index> time;
		std::vector<lp_index> parameter;
	};
	struct ps {
		vec Points;
		vec Weights;
		dsmatrix  D;
		dsmatrix  Diag;
		dsmatrix  Doffdiag;
	};

	struct LpCalculateData
	{
/*		LpCalculateData();*/
		bool autoscale; ///whether or not using autoscale.Wraning!The best choice is scaling the ocp problem your self
		std::vector<std::vector<lp_index>> SIZES_;
		/* SIZES
		       inner std::vector is <nstate ncontrol nparameters npaths nevents> 
			   outer std::vector is every phase
			   do not change this outside sizechecker!!!!
		*/
		int numphases_;
		int numlinkpairs_;
		int numlinks_;

		vector<int> variables;
		vector<int> constraints;

		vec linmin;
		vec linmax;
		vector<double> varbounds_min;
		vector<double> varbounds_max;
		vector<double> conbounds_min;
		vector<double> conbounds_max;

		vec scaled_varbounds_min;
		vec scaled_varbounds_max;
		vec scaled_conbounds_min;
		vec scaled_conbounds_max;

		vector<int> totalnodes_perphase;
		vector<shared_ptr<indices>> phase_indices; 
		vector<std::vector<lp_index>> variable_indices;
		vector<std::vector<lp_index>> constraint_indices;
		vector<std::vector<lp_index>> link_indices;
		vec nlpGuessVector;
		vec scaled_nlpGuessVector;
		//the sparse structure of jacobi matrix
		vec jac_I;
		vec jac_J;

		//the sparse structure of hessian matrix
		vec hessian_I;
		vec hessian_J;

		/*!<jacobi dependencies in every phase  (nstate+npath>*(nstate+ncontrol+nparamemter)
		!*/
		std::vector<umat> allPhaseDependencies;

		vec varshift;
		vec varscale;
		vec funscale;
		vec objScale;
		dsmatrix AlinearMatrix;//Linear function,about duration and phase by phase time link
		vector<vec> nlpGuess;
		vector<shared_ptr<ps>> PS;
		vec scaled_nlpreturn_x;
		vec scaled_nlpreturn_lambda;

		vec nlpreturn_x;
		vec nlpreturn_lambda;
		double return_obj;
		double optcontrol_cost;
		vector<shared_ptr<SolutionData>> result;///the result of every grid
		size_t current_grid;///Set in MeshRefiner,current grid in Lpopc calculation
// 	private:
// 		LpCalculateData(const LpCalculateData&);
// 		LpCalculateData& operator=(const LpCalculateData&);
// 		//void AllPrint(shared_ptr< LpReporter>&);
	};

	struct LpCauculateHistoryData
	{

	};
}

#endif