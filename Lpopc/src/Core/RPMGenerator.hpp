// Copyright (C) 2014-2015 Xue Zhichen, Wang Yujie,Wang Na 
// All Rights Reserved.
// This file is a part of LPOPC , published under the Eclipse Public License.
// Author:Xue Zhichen 7/10 2015   17:54
// Email:eddy_lpopc@163.com
#ifndef RPMGENERATOR_HPP
#define RPMGENERATOR_HPP
#include"LpConf.h"
#include "LpReporter.hpp"
#include "LpSparseMatrix.h"
#include "LpException.hpp"
#include<map>
#include <vector>
namespace Lpopc
{
	LP_DECLARE_EXCEPTION(RPMGENERATOR_ERROR);
	class RPMGenerator
	{
	/*<This class calculate the LGR points and weights,
	* and return the differential matrix, intergration matrix
	*/
	public:
		RPMGenerator();
		static void GetLGRPoints(const int N,vec &x,vec& w);
		static void GetLGRPointsImp(const int N, vec &x, vec& w);
		void initialize(int sections,
			std::vector<double>& mesh_points,
			std::vector<lp_index>& node_per_interval);
		vec& RPMpoints(){return RPM_points_;}
		vec& RPMweights(){return RPM_weights_;}
		dsmatrix&  DifferentiationMatrix(){return RPM_Differentiation_matrix_;}
		dsmatrix&  DifferentiationMatrixDiag(){return RPM_Differentiation_matrix_diag_;}
		dsmatrix&  DifferentiationMatrixOffDiag(){return RPM_Differentiation_matrix_off_diag_;}
		dsmatrix&  IntegrationMatrix(){ return RPM_Integration_matrix_; };
		dsmatrix&  UnityMatrix(){ return RPM_Unity_matrix_; }

	private:
		shared_ptr<LpReporter> reporter_;
		void operator=(const RPMGenerator&);
		RPMGenerator(const RPMGenerator&);
		int sections_;
		std::vector<double> mesh_points_;
		std::vector<lp_index> node_per_interval_;

		vec RPM_points_;
		vec RPM_weights_;
		dsmatrix RPM_Differentiation_matrix_;
		dsmatrix RPM_Differentiation_matrix_diag_;
		dsmatrix RPM_Differentiation_matrix_off_diag_;
		dsmatrix RPM_Integration_matrix_;
		dsmatrix RPM_Unity_matrix_;
		static std::map<int,vec> LGR_points_;///calculate points and weights onece
		static std::map<int,vec> LGR_wights_;
		void CompositeD(std::vector<mat>& Dsect,dsmatrix& sparse_matrix);
		void CompositeA(std::vector<mat>& Asect, dsmatrix& sparse_matrix);
		void CollocD(vec&x,mat&D,mat&Dd,mat& Do )const;
		bool CheckDataReady(int sections,std::vector<double>& mesh_points,std::vector<lp_index>& node_per_interval);
		
	};
}
#endif