// Copyright (C) 2014-2015 Xue Zhichen, Wang Yujie,Wang Na 
// All Rights Reserved.
// This file is a part of LPOPC , published under the Eclipse Public License.
// Author:Xue Zhichen 7/10 2015   17:50
// Email:eddy_lpopc@163.com
#include "LpSparseMatrix.h"
#include<iostream>
namespace Lpopc{

void Lpopc::dsmatrix::print( const shared_ptr< LpReporter>& reporter,std::string matrix_name ) const
{
	reporter->Printf()->debug("\n");
	reporter->Printf()->debug(
		"dmatrix \"{}\" with {}rows and {} columns:\n",matrix_name.c_str(),m_,n_);
	if (CheckAllocated()&&CheckInitialized())
	{
		int *   rows_ptr=GetSArrayConst()->GetRowsPtr();
		int * cols_ptr=GetSArrayConst()->GetColsPtr();
		double* data_ptr=GetSArrayConst()->GetDataPtr();
		for (int i=0; i<GetLength(); i++) {
			
				reporter->Printf()->info(
					"{}[{5d},{5d}]={16f}\n",
					matrix_name.c_str(), rows_ptr[i], cols_ptr[i], data_ptr[i]);
		}
	} 
	else
	{
		reporter->Printf()->warn("The matrix has not yet been initialized!\n");
	}
}

void Lpopc::dsmatrix::print( std::string matrix_name ) const
{

	printf("\ndmatrix \"%s\" with %d rows and %d columns:\n",matrix_name.c_str(),m_,n_);
	if (CheckAllocated()&&CheckInitialized())
	{
		int *   rows_ptr=GetSArrayConst()->GetRowsPtr();
		int * cols_ptr=GetSArrayConst()->GetColsPtr();
		double* data_ptr=GetSArrayConst()->GetDataPtr();
		for (int i=0; i<GetLength(); i++) {

			printf("%s[%5d,%5d]=%23.16e\n",
				matrix_name.c_str(), rows_ptr[i], cols_ptr[i], data_ptr[i]);
		}
	} 
	else
	{
		printf("The matrix has not yet been initialized!\n");
	}
}
dsmatrix dsmatrix:: Sparse( vec row_vector,vec col_vector,vec value_vector, int nrow,int ncol )
{
	LP_ASSERT_EXCEPTION(row_vector.n_elem==col_vector.n_elem,
		DSMATRIX_ERROR,"the number of rows don't match the number of columns elements ")
		LP_ASSERT_EXCEPTION(row_vector.n_elem==value_vector.n_elem,
		DSMATRIX_ERROR,"the number of rows don't match the number of values ")
		LP_ASSERT_EXCEPTION(col_vector.n_elem==value_vector.n_elem,
		DSMATRIX_ERROR,"the number of columns don't match the number of columns elements ")
		int nonzeros=row_vector.n_elem;
	double *row_data=row_vector.memptr();
	double *col_data=col_vector.memptr();
	double *value_dat=value_vector.memptr();

	int *rows=new int[nonzeros];
	int *cols=new int  [nonzeros];
	double* values=new double[nonzeros];
	for (int i = 0; i < nonzeros; i++)
	{
		rows[i]=int(row_data[i]);
		cols[i]=int(col_data[i]);
		values[i]=double(value_dat[i]);
	}

	dsmatrix result_smatrix(rows,cols,values,nrow,ncol,nonzeros,true);
	return result_smatrix;
}

void  dsmatrix::GeneratRowColValue(mat& dataMatrix, vec& irow, vec& jcol, vec& values, int rowShift, int colShift)
{
	LP_ASSERT_EXCEPTION(dataMatrix.n_elem!=0,
		DSMATRIX_ERROR, "dataMatrix has no ele in GeneratRowColValue ")
	irow=zeros(dataMatrix.n_elem,1);
	jcol = zeros(dataMatrix.n_elem, 1);
	values = zeros(dataMatrix.n_elem, 1);
	int nirow = dataMatrix.n_rows;
	int njcol = dataMatrix.n_cols;
	int k=0;
	for(int j=0;j<njcol;j++){
		for (int i=0; i<nirow; i++) {
			irow[k]= i+rowShift;
			jcol[k]=j+colShift;
			k++;
		}
	}

	double* value_data=dataMatrix.memptr();
	for (int i = 0; i < values.n_elem; i++)
	{
		values[i]=value_data[i];
	}
}


double dsmatrix::operator()( int nrow,int ncol ) const
{
	LP_ASSERT_EXCEPTION(CheckInitialized(),
		DSMATRIX_ERROR,"The dsmatrix hasn't been initialized!")
		LP_ASSERT_EXCEPTION(nrow>-1&&ncol>-1&&nrow<m_&&ncol<n_,
		DSMATRIX_ERROR,"Out of range index in dsmatrix")
		int* rowindex= GetSArrayConst()->GetRowsPtr();
	    int* colindex= GetSArrayConst()->GetColsPtr();
		double* values=GetSArrayConst()->GetDataPtr();
		double retval=datum::eps;
		for (int i = 0; i < GetLength(); i++)
		{
			if (rowindex[i]==nrow && colindex[i]==ncol)
			{
				retval = values[i];
				break;
			}
		}
		return retval;
}

mat dsmatrix::operator*( mat& op_matrix )
{
	LP_ASSERT_EXCEPTION(op_matrix.n_elem!=0,
		DSMATRIX_ERROR,"op_matrix  hasn't non elem")
		LP_ASSERT_EXCEPTION(CheckInitialized(),
		DSMATRIX_ERROR,"The dsmatrix hasn't been initialized!")
	LP_ASSERT_EXCEPTION(GetnCols()==op_matrix.n_rows,
	DSMATRIX_ERROR,"dsmatrix's cols didn't match the op_matrix's rows")
	int* rowindex= GetSArrayConst()->GetRowsPtr();
	int* colindex= GetSArrayConst()->GetColsPtr();
	double* values=GetSArrayConst()->GetDataPtr();

	mat result_matrix=zeros(m_,op_matrix.n_cols);
	double *new_values=result_matrix.memptr();
	int i,j;
	for (int icol = 0; icol < op_matrix.n_cols; icol++)
	{
		vec temcol=op_matrix.col(icol);
		
		int tem=0;
		for (int k = 0; k < this->GetLength(); k++)
		{
			j=colindex[k];
			i=rowindex[k];
			new_values[i+icol*m_]+=values[k]*temcol[j];
		}
	}
	return result_matrix;

	 
	//The code below is different from above.The below will return an dsmatrix.
	//but I didn't test it.Test it before you use ////By Xue Zhichen
	/*
	int new_nonzero=0;
	int has_nonzeros_rows=0;
	int* rowindex= GetSArrayConst()->GetRowsPtr();
	int* colindex= GetSArrayConst()->GetColsPtr();
	double* values=GetSArrayConst()->GetDataPtr();
	
	int* nonzeros_in_row=new int[m_];

		for (int irow=0;irow<m_;irow++)
		{
			nonzeros_in_row[irow]=0.0;
			for (int i = 0; i < GetLength(); i++)
			{
				if (rowindex[i]==irow)
				{
					nonzeros_in_row[irow]=1.0;
					has_nonzeros_rows++;
					break;
				}
				
			}
		}//end for

		new_nonzero=has_nonzeros_rows*op_matrix.GetnCols();

   int  *which_row_has_nonzero=new int[has_nonzeros_rows];

   int last_has_nonzero_row=0;
   for (int i = 0; i < has_nonzeros_rows; i++)
   {
	   for (int j = last_has_nonzero_row; j < m_; j++)
	   {
		   if (nonzeros_in_row[j]==1.0)
		   {
			   which_row_has_nonzero[i]=j;
			   last_has_nonzero_row=j+1;
			   break;
		   }
		   
	   }

   }

	
	int*new_rows =new int[new_nonzero];
	int*new_cols  =new int[new_nonzero];
	double* new_values=new double[new_nonzero];
	for (int i = 0; i < new_nonzero; i++)
	{
		new_values[i]=0.0;
	}

	int new_index=0;
	int i,j;
	for (int icol = 0; icol < op_matrix.GetnCols(); icol++)
	{
		dmatrix temcol=op_matrix.Column(icol);
		int tem=0;
		for (int k = 0; k < this->GetLength(); k++)
		{
			j=colindex[k];
			i=rowindex[k];
			new_values[i+icol*has_nonzeros_rows]+=values[k]*temcol.GetElem(j,0);
		}
		for (int q = 0; q < has_nonzeros_rows; q++)
		{
			new_cols[q+icol*has_nonzeros_rows]=icol;
			new_rows[q+icol*has_nonzeros_rows]=which_row_has_nonzero[q];
		}
	}


	dmatrix result_dsmatrix(new_rows,new_cols,new_values,m_,op_matrix.GetnCols(),new_nonzero,true);
	//delete[] which_row_has_nonzero;
	//delete[] nonzeros_in_row;*/
//	return result_dsmatrix;
}



void dsmatrix::Find( dsmatrix op_smatrix,vec& imatrix,vec& jmatrix,vec& vmatrix )
{
	LP_ASSERT_EXCEPTION(op_smatrix.CheckInitialized(),
		DSMATRIX_ERROR,"dsmatrix  hasn't been initialized!")
	int nonzeros=0;
	int* data_i=op_smatrix.GetSArrayConst()->GetRowsPtr();
	int* data_j=op_smatrix.GetSArrayConst()->GetColsPtr();
	double*data_v=op_smatrix.GetSArrayConst()->GetDataPtr();
	for (int i = 0; i < op_smatrix.GetLength(); i++)
	{
		if (data_v[i]!=0.0) nonzeros++;
	}
	double* new_i=new double[nonzeros];
	double* new_j=new double[nonzeros];
	double* new_v=new double[nonzeros];
	int incx=0;
	for (int i = 0; i < op_smatrix.GetLength(); i++)
	{
		if (data_v[i]!=0.0)
		{
			new_i[incx]=(double)data_i[i];
			new_j[incx]=(double)data_j[i];
			new_v[incx]=(double)data_v[i];
			incx++;
		}
	}
	vec ret_i(new_i,nonzeros,false);
	vec ret_j(new_j,nonzeros,false);
	vec ret_v(new_v,nonzeros,false);
	imatrix=ret_i;
	jmatrix=ret_j;
	vmatrix=ret_v;
}
 sp_mat dsmatrix::GetArmaSP_Mat(dsmatrix op_dsmatrix)
{
	urowvec S_I(op_dsmatrix.GetLength()), S_J(op_dsmatrix.GetLength());
	vec S_V(op_dsmatrix.GetLength());
	int* i_ptr = op_dsmatrix.GetSArrayConst()->GetRowsPtr();
	int* j_ptr = op_dsmatrix.GetSArrayConst()->GetColsPtr();
	double* v_ptr = op_dsmatrix.GetSArrayConst()->GetDataPtr();
	for (size_t i = 0; i < op_dsmatrix.GetLength(); i++)
	{
		S_I[i] = (lp_index)i_ptr[i];
		S_J[i] = (lp_index)j_ptr[i];
		S_V[i] = v_ptr[i];
	}
	sp_mat return_spmat(join_vert(S_I, S_J), S_V, op_dsmatrix.GetnRows(), op_dsmatrix.GetnCols());
	return return_spmat;
}
}//namespace Lpopc











