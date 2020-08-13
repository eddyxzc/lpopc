// Copyright (C) 2014-2015 Xue Zhichen, Wang Yujie,Wang Na 
// All Rights Reserved.
// This file is a part of LPOPC , published under the Eclipse Public License.
// Author:Xue Zhichen 7/12 2015   20:19
// Email:eddy_lpopc@163.com
#include "LpopcIpopt.h"
namespace Lpopc
{


	bool LpopcIpopt::get_nlp_info(Index& n, Index& m, Index& nnz_jac_g, Index& nnz_h_lag, IndexStyleEnum& index_style)
	{
		n = Data_->varbounds_min.size();
		m = Data_->conbounds_min.size() + Data_->linmin.n_elem;
		vec I, J;
		NlpWrapper_->GetConsSparsity(I, J);
		nnz_jac_g = I.n_elem;
		//////////////////////////////////////////////////////////////////////////
		NlpWrapper_->GetHessainSparsity(I, J);
		nnz_h_lag = I.n_elem;
		//////////////////////////////////////////////////////////////////////////
		index_style = TNLP::C_STYLE;
		return true;
	}

	bool LpopcIpopt::get_bounds_info(Index n, Number* x_l, Number* x_u, Index m, Number* g_l, Number* g_u)
	{
		//Fix me exec once !
		if (Data_->autoscale)// autoscale
		{
			for (size_t i = 0; i < Data_->scaled_varbounds_min.n_elem; i++)
			{
				x_l[i] = Data_->scaled_varbounds_min(i);
			}
			for (size_t i = 0; i < Data_->scaled_varbounds_max.n_elem; i++)
			{
				x_u[i] = Data_->scaled_varbounds_max(i);
			}

			for (size_t i = 0; i < Data_->scaled_conbounds_min.n_elem; i++)
			{
				g_l[i] = Data_->scaled_conbounds_min(i);
			}
			for (size_t i = 0; i < Data_->scaled_conbounds_max.n_elem; i++)
			{
				g_u[i] = Data_->scaled_conbounds_max(i);
			}
		}
		else
		{
			for (size_t i = 0; i < Data_->varbounds_min.size(); i++)
			{
				x_l[i] = Data_->varbounds_min[i];
			}
			for (size_t i = 0; i < Data_->varbounds_max.size(); i++)
			{
				x_u[i] = Data_->varbounds_max[i];
			}

			for (size_t i = 0; i < Data_->conbounds_min.size(); i++)
			{
				g_l[i] = Data_->conbounds_min[i];
			}
			for (size_t i = 0; i < Data_->conbounds_max.size(); i++)
			{
				g_u[i] = Data_->conbounds_max[i];
			}
		}
		

		//AlinearMatrix
		for (size_t i = Data_->conbounds_min.size(); i < m; i++)
		{
			g_l[i] = Data_->linmin(i - Data_->conbounds_min.size());
		}
		for (size_t i = Data_->conbounds_max.size(); i < m; i++)
		{
			g_u[i] = Data_->linmax(i - Data_->conbounds_min.size());
		}

		return true;
	}

	bool LpopcIpopt::get_starting_point(Index n, bool init_x, Number* x, bool init_z, Number* z_L, Number* z_U, Index m, bool init_lambda, Number* lambda)
	{
		assert(init_x == true);
		assert(init_z == false);
		assert(init_lambda == false);
		if (Data_->autoscale)// autoscale
		{
			for (lp_index i = 0; i < Data_->varbounds_min.size(); i++)
			{
				x[i] = Data_->scaled_nlpGuessVector(i);
			}
		}
		else{
			for (lp_index i = 0; i < Data_->varbounds_min.size(); i++)
			{
				x[i] = Data_->nlpGuessVector(i);
				//std::cout<<"x start"<<x[i]<<std::endl;
			}
		}
		return true;
	}

	bool LpopcIpopt::eval_f(Index n, const Number* x, bool new_x, Number& obj_value)
	{
		double* x_all = new double[n];
		for (int i = 0; i<n; i++)
		{
			x_all[i] = x[i];
		}
		vec x_vec(x_all, n, false);
		obj_value = NlpWrapper_->GetObjFun(x_vec);
		return true;
	}

	bool LpopcIpopt::eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f)
	{
		double* x_all = new double[n];
		for (int i = 0; i<n; i++)
		{
			x_all[i] = x[i];
		}
		vec x_vec(x_all, n, false);
		vec grad_fvec;
		NlpWrapper_->GetObjGrad(x_vec, grad_fvec);
		for (int i = 0; i < grad_fvec.n_elem; i++)
		{
			grad_f[i] = grad_fvec(i);
		}
		return true;
	}

	bool LpopcIpopt::eval_g(Index n, const Number* x, bool new_x, Index m, Number* g)
	{
		double* x_all = new double[n];
		for (int i = 0; i<n; i++)
		{
			x_all[i] = x[i];
		}
		vec x_vec(x_all, n,  false);
		vec cons_vec;
		NlpWrapper_->GetAllCons(x_vec, cons_vec);
		for (int i = 0; i < cons_vec.n_elem; i++)
		{
			g[i] = cons_vec(i);
		}
		return true;
	}

	bool LpopcIpopt::eval_jac_g(Index n, const Number* x, bool new_x, Index m, Index nele_jac, Index* iRow, Index *jCol, Number* values)
	{
		vec I, J, V;
		
		if (values == NULL)
		{
			NlpWrapper_->GetConsSparsity(I, J);
			for (int i = 0; i < I.n_elem; i++)
			{
				iRow[i] =(Index) I(i);
				jCol[i] = (Index)J(i);
			}
		}
		else{
			double* x_all = new double[n];
			for (int i = 0; i<n; i++)
			{
				x_all[i] = x[i];
				//	std::cout<<"x "<<x[i]<<std::endl;
			}
			vec x_vec(x_all, n,  false);
			NlpWrapper_->GetConsJacbi(x_vec,V);
			int tem_len = V.n_elem;
			for (int j = 0; j <tem_len; j++)
			{
				values[j] = V(j);
			}
		}
		return true;
	}

	bool LpopcIpopt::eval_h(Index n, const Number* x, bool new_x, Number obj_factor, Index m, const Number* lambda, bool new_lambda, Index nele_hess, Index* iRow, Index* jCol, Number* values)
	{
		vec I, J, V;
			
		if (values == NULL)
		{
			NlpWrapper_->GetHessainSparsity(I, J);
			for (int i = 0; i < I.n_elem; i++)
			{
				iRow[i] = (Index)I(i);
				jCol[i] = (Index)J(i);
			}
		}
		else{
			double* x_all = new double[n];
			double* lambda_all = new double[m];
			for (int i = 0; i < n; i++)
			{
				x_all[i] = x[i];
				//	std::cout<<"x "<<x[i]<<std::endl;
			}
			vec x_vec(x_all, n, false);
			for (int i = 0; i < m - 1; i++)
			{
				lambda_all[i] = lambda[i];
			}
			vec lambda_vec(lambda_all, m , false);
			NlpWrapper_->GetHessainValue( x_vec,obj_factor, lambda_vec, V);
			size_t tem_len = V.n_elem;
			for (int j = 0; j < tem_len; j++)
			{
				values[j] = V(j);
			}
		}
		return true;
	}

	void LpopcIpopt::finalize_solution(SolverReturn status, Index n, const Number* x, const Number* z_L, const Number* z_U, Index m, const Number* g, const Number* lambda, Number obj_value, const IpoptData* ip_data, IpoptCalculatedQuantities* ip_cq)
	{
		/*	std::cout << std::endl << std::endl << "Solution of the primal variables, x" << std::endl;
		for (Index i=0; i<n; i++) {
		std::cout << "x[" << i << "] = " << x[i] << std::endl;
		}*/
		/*
		std::cout << std::endl << std::endl << "Solution of the bound multipliers, z_L and z_U" << std::endl;
		for (Index i=0; i<n; i++) {
		std::cout << "z_L[" << i << "] = " << z_L[i] << std::endl;
		}
		for (Index i=0; i<n; i++) {
		std::cout << "z_U[" << i << "] = " << z_U[i] << std::endl;
		}
		*/
		std::cout << std::endl << std::endl << "Objective value" << std::endl;
		std::cout << "f(x*) = " << obj_value << std::endl;
		vec return_x(x, n);
		Data_->nlpreturn_x = return_x;

		vec return_lambda(lambda, m);
		Data_->nlpreturn_lambda = return_lambda;

		Data_->return_obj = obj_value;
		/*	std::cout << std::endl << "Final value of the constraints:" << std::endl;
		for (Index i=0; i<m ;i++) {
		std::cout << "g(" << i << ") = " << g[i] << std::endl;
		}*/
	}

}