// Copyright (C) 2014-2015 Xue Zhichen, Wang Yujie,Wang Na 
// All Rights Reserved.
// This file is a part of LPOPC , published under the Eclipse Public License.
// Author:Xue Zhichen 7/10 2015   17:54
// Email:eddy_lpopc@163.com
// 
#include "RPMGenerator.hpp"
#include <limits>
#include <utility>
namespace Lpopc
{
	RPMGenerator::RPMGenerator()
	{}
	    //Reference
	    //C.Canuto, M.Y.Hussaini, A.Quarteroni, T.A.Tang, "Spectral Methods
		//   in Fluid Dynamics," Section 2.3. Springer-Verlag 1987
	void RPMGenerator::GetLGRPoints( const int iniN,vec &x,vec &w )
	{
		
		std::map<int, vec>::iterator xitor = LGR_points_.find(iniN);
		std::map<int, vec>::iterator witor = LGR_wights_.find(iniN);
		if ((xitor!=LGR_points_.end() )&& (witor!=LGR_wights_.end()))
		{
			
			x = xitor->second;
			w = witor->second;
		}
		else
		{
			///Calculate new points and wights
			GetLGRPointsImp(iniN, x, w);
			auto ret1 = LGR_points_.insert(std::make_pair(iniN,x));
			auto ret2 = LGR_wights_.insert(std::make_pair(iniN, w));

			if (!ret1.second || !ret2.second)
			{
				LP_THROW_EXCEPTION(RPMGENERATOR_ERROR, "Insert LGR points or weights failed!");
			}

		}
	}//getLGRPoints

	void RPMGenerator::initialize( int sections, std::vector<double>& mesh_points,std::vector<lp_index>& node_per_interval )
	{
		sections_=sections;
		mesh_points_=mesh_points;
		node_per_interval_=node_per_interval;
		std::vector<vec> sSeg(sections_);
		std::vector<vec> wscaled(sections_);
		std::vector<mat> Dsect(sections_);
		std::vector<mat> Ddsect(sections_);
		std::vector<mat> Dosect(sections_);
		std::vector<mat> Asect(sections_);
		std::vector<mat> Bsect(sections_);
		int tau_size_all=0;
		for(int i=0;i<sections_;i++)
		{	
			vec x(node_per_interval_[i]),w(node_per_interval_[i]);
			GetLGRPoints(node_per_interval_[i],x,w);
			double tspan =mesh_points_[i+1]-mesh_points_[i];
			vec sSegi(x.n_rows);
			vec wscaledi(x.n_rows);
			vec sall(sSegi.n_elem+1);
			mat D2(node_per_interval_[i],node_per_interval_[i]+1);
			mat Dd2(node_per_interval_[i],node_per_interval_[i]+1);
			mat Do2(node_per_interval_[i],node_per_interval_[i]+1);
			sSegi=x+1;
			sSegi*=tspan/2.0;
			sSegi+=mesh_points_[i];
			sall.subvec(0,sSegi.n_elem-1)=sSegi;
			sall(sSegi.n_elem,0)=mesh_points_[i+1];
			wscaledi=w/2;
			wscaledi*=tspan;

			CollocD(sall,D2,Dd2,Do2);
			/*
			% Compose the points, weights, and matrices computed in each mesh %
			% interval into single quantities defined across the entire phase %
			*/
			sSeg[i]=sSegi;
			wscaled[i]=wscaledi;
			Dsect[i]= D2;
			Ddsect[i]= Dd2;
			Dosect[i]= Do2;
			Asect[i] = inv(D2(span::all, span(1, D2.n_cols - 1)));
			Bsect[i] = zeros(Dsect[i].n_rows,Dsect[i].n_cols);
			Bsect[i].col(0) += 1 ;
			tau_size_all+=sSeg[i].n_elem;
		}//for i

		int itor=0;
		RPM_points_=zeros(tau_size_all,1);
		RPM_weights_=zeros(tau_size_all,1);
		for (int i=0;i<sSeg.size();i++)
		{	
			RPM_points_.subvec(itor,itor+sSeg[i].n_elem-1)=sSeg[i];
			RPM_weights_.subvec(itor,itor+wscaled[i].n_elem-1)=wscaled[i];
			itor=itor+sSeg[i].n_elem;
		}
		CompositeD(Dsect,RPM_Differentiation_matrix_);
		CompositeD(Ddsect,RPM_Differentiation_matrix_diag_);
		CompositeD(Dosect,RPM_Differentiation_matrix_off_diag_);
		CompositeA(Asect, RPM_Integration_matrix_);
		CompositeD(Bsect, RPM_Unity_matrix_);
	}

	void RPMGenerator::CollocD( vec&x,mat&D,mat&Dd,mat& Do ) const
	{
		int N1=x.n_elem-1;
		int M=x.n_elem,M1=M+1,M2=M*M;
		mat Y=repmat(x,1,M);
		mat Ydiff=eye(M,M);
		Ydiff=Ydiff+Y-trans(Y);
		mat ww = repmat(1 / (prod(Ydiff, 1)), 1, M);
		D = ww / (trans(ww) % (Ydiff));
	    mat ret2 =1 -sum(D);
		int j=0;
		for (int i=0;i<M2;i+=M1)
		{
			D(j,j)=ret2(0,j);
			j++;
		}
		ret2=-trans(D);
		D=ret2.rows(0,ret2.n_rows-2);
		Dd=zeros(N1,N1+1);
		mat ret=D.cols(0,M-2);
		ret=ret.diag();
		Dd.cols(0,N1-1)=diagmat(ret);// Diag diff matrix
		Do=D-Dd;//off-diag diff matrix*/
	}

	void RPMGenerator::CompositeD( std::vector<mat>& Dsect,dsmatrix& sparse_matrix )
	{
		int sections = Dsect.size();
		int totalrows = 0;
		int totalcols = 0;
		int nonzerros = 0;
		std::vector<lp_index> rows(sections);
		std::vector<lp_index> cols(sections);
		std::vector<vec> irows(sections);
		std::vector<vec> jcols(sections);
		std::vector<vec> values(sections);
		int nodes = 0;
		for (int i = 0; i < sections; i++)
		{
			nodes += Dsect[i].n_rows;
		}
		int rowshift = 0;
		int colshift = 0;
		for (int i = 0; i < sections; i++){
			rows[i] = Dsect[i].n_rows;
			cols[i] = rows[i] + 1;
			vec irow(rows[i] * cols[i], 1);
			vec	jcol(rows[i] * cols[i], 1);
			vec valu(rows[i] * cols[i], 1);
			dsmatrix::GeneratRowColValue(Dsect[i], irow, jcol, valu, rowshift, colshift);
			nonzerros += Dsect[i].n_elem;
			rowshift += rows[i];
			colshift += cols[i]-1;
			totalrows = totalrows + rows[i];
			totalcols = totalcols + cols[i];
			irows[i] = irow;
			jcols[i] = jcol;
			values[i] = valu;
		}
		vec vec_i(nonzerros, 1);
		vec vec_j(nonzerros, 1);
		vec vec_v(nonzerros, 1);
		int k = 0;
		for (int i = 0; i < sections; i++){
			for (int j = 0; j < rows[i] * cols[i]; j++){
				vec_i(k) = irows[i](j);
				vec_j(k) = jcols[i](j);
				vec_v(k) = values[i](j);
				k++;
			}
		}
		sparse_matrix = dsmatrix::Sparse(vec_i, vec_j, vec_v, nodes, nodes + 1);
		dsmatrix::Find(sparse_matrix, vec_i, vec_j, vec_v);
		sparse_matrix = dsmatrix::Sparse(vec_i, vec_j, vec_v, nodes, nodes + 1);
	}

	bool RPMGenerator::CheckDataReady(int sections, std::vector<double>& mesh_points,std::vector<lp_index>& node_per_interval)
	{
		bool testresult = false;
		if (sections == sections_)
		{
			if (mesh_points == mesh_points_)
			{
				if (node_per_interval==node_per_interval_)
				{
					testresult = true;
				}
			}

		}
		return testresult;
	}

	void RPMGenerator::CompositeA(std::vector<mat>& Asect, dsmatrix& sparse_matrix)
	{
		int sections = Asect.size();
		int totalrows = 0;
		int totalcols = 0;
		int nonzerros = 0;
		std::vector<int> rows(sections);
		std::vector<int> cols(sections);
		std::vector<vec> irows(sections);
		std::vector<vec> jcols(sections);
		std::vector<vec> values(sections);
		int nodes = 0;
		for (int i = 0; i < sections; i++)
		{
			nodes += Asect[i].n_rows;
		}
		int rowshift = 0;
		int colshift = 0;
		for (int i = 0; i < sections; i++){
			//Asect[i].print("[asect]");
			rows[i] = Asect[i].n_rows;
			cols[i] = rows[i];
			vec irow(rows[i] * cols[i], 1);
			vec	jcol(rows[i] * cols[i], 1);
			vec valu(rows[i] * cols[i], 1);
			dsmatrix::GeneratRowColValue(Asect[i], irow, jcol, valu, rowshift, colshift);
			nonzerros += Asect[i].n_elem;
			rowshift += rows[i];
			colshift += cols[i];
			totalrows = totalrows + rows[i];
			totalcols = totalcols + cols[i];
			irows[i] = irow;
			jcols[i] = jcol;
			values[i] = valu;
		}
		vec vec_i(nonzerros, 1);
		
		vec vec_j(nonzerros, 1);
		vec vec_v(nonzerros, 1);
		int k = 0;
		for (int i = 0; i < sections; i++){
			for (int j = 0; j < rows[i] * cols[i]; j++){
				vec_i(k) = irows[i](j);
				vec_j(k) = jcols[i](j);
				vec_v(k) = values[i](j);
				k++;
			}
		}
		sparse_matrix = dsmatrix::Sparse(vec_i, vec_j, vec_v, nodes, nodes );
		dsmatrix::Find(sparse_matrix, vec_i, vec_j, vec_v);
		sparse_matrix = dsmatrix::Sparse(vec_i, vec_j, vec_v, nodes, nodes );
	}

	void RPMGenerator::GetLGRPointsImp(const int iniN, vec &x, vec& w)
	{
		int N = iniN - 1, N1 = N + 1;
		double eps = datum::eps;
		// Initial guess for LGR nodes
		vec ret_matrix = linspace(0, N, N + 1)*((2 * datum::pi) / (2 * N + 1));
		x = -1 * cos(ret_matrix);
		// The Legendre Vandermonde Matrix
		mat P = zeros(N1, N1 + 1);
		//dmatrix ret_matrix1=x-2.0;
		vec xold = x;
		xold.fill(2);
		//dmatrix temp;

		while (max(abs(x - xold)) > eps)
		{
			xold = x;
			P.col(0) = ones(P.n_rows, 1);
			P.col(1) = x;
			for (int k = 1; k < N1; k++)
			{
				vec ret1 = x*(2 * k + 1) % (P.col(k)) - (P.col(k - 1)*k);
				P.col(k + 1) = ret1 / (k + 1);
			}
			vec ret2 = ones(x.n_rows, 1) - xold;
			ret2 /= N1;
			vec temp = P.col(N1 - 1) + P.col(N1);
			ret2 = ret2 % (temp);
			ret2 = xold - (ret2 / (P.col(N1 - 1) - P.col(N1)));
			ret2 = ret2.subvec(1, N);
			x.rows(1, N) = ret2;

		}
		w = zeros(N1, 1);
		w(0, 0) = 2.0 / (N1*N1);
		vec ret4 = (P.col(N)*N1);
		mat temp2 = (1 - x) / (ret4 % (ret4));
		w.subvec(1, N) = temp2.rows(1, N);
	}

	std::map<int, vec> RPMGenerator::LGR_wights_ = std::map<int, vec>();

	std::map<int, vec> RPMGenerator::LGR_points_ = std::map<int, vec>();

	

}//namespace lpopc