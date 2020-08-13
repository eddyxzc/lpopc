// Copyright (C) 2014-2015 Xue Zhichen, Wang Yujie,Wang Na 
// All Rights Reserved.
// This file is a part of LPOPC , published under the Eclipse Public License.
// Author:Xue Zhichen 7/10 2015   17:50
// Email:eddy_lpopc@163.com
#ifndef LPSPARSEMATRIX_HPP
#define LPSPARSEMATRIX_HPP
#include "LpSparseArray.h"
#include "LpConf.h"
#include"armadillo"
#include "LpException.hpp"
namespace Lpopc
{
	 LP_DECLARE_EXCEPTION(DSMATRIX_ERROR);
	 /**This class is an sparse matrix ,just for storing sparse matrix in triplet form
	 *Modified from GNU OCATVE .
	 * this class will be replaced with arma::sp_mat in the future.
	 */
	class dsmatrix:public SArray
	{
	public:
		dsmatrix():m_(0),n_(0),SArray() {}
		dsmatrix(int nRow,int nCol,int nNzeros)
			:SArray(nNzeros),m_(nRow),n_(nCol){}
		dsmatrix(int* rows,int* cols,double* values,int m,int n,int nonzeros,bool has_allocated)
			:SArray(rows,cols,values,nonzeros,has_allocated),m_(m),n_(n) {}

		dsmatrix(const dsmatrix& copy_matrix)
			:m_(copy_matrix.GetnRows()),n_(copy_matrix.GetnCols()),SArray(copy_matrix){}
		dsmatrix& operator=(const dsmatrix&equal_matrix)
		{
			m_=equal_matrix.GetnRows();
			n_=equal_matrix.GetnCols();
			(*this).SArray::operator=(equal_matrix);
			return *this;
		}
		virtual ~dsmatrix(){}
		int GetnRows()const {return m_;}
		int GetnCols()const{return n_;}
		void print( const shared_ptr< LpReporter>& reporter, 
			std::string matrix_name )const;
		void print (std::string matrix_name) const;
		double operator()(int nrow,int ncol)const;
		mat operator*(mat& op_matrix);//note !sorted by column,slow!
		static dsmatrix Sparse( vec row_vector,vec col_vector,vec value_vector, int nrow,int ncol);
		static void  Find(dsmatrix op_smatrix,vec& imatrix,vec& jmatrix,vec& vmatrix);
		static void  GeneratRowColValue(mat& dataMatrix,vec& irow,vec& jcol,vec& values,int rowShift,int colShift);
		static sp_mat  GetArmaSP_Mat(dsmatrix op_dsmatrix);
	private:
		int m_;
		int n_;
	};
}
#endif // !LPSPARSEMATRIX_HPP
