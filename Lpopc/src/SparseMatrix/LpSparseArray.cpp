// Copyright (C) 2014-2015 Xue Zhichen, Wang Yujie,Wang Na 
// All Rights Reserved.
// This file is a part of LPOPC , published under the Eclipse Public License.
// Author:Xue Zhichen 7/10 2015   17:50
// Email:eddy_lpopc@163.com
#include "LpSparseArray.h"

namespace Lpopc
{

	
	SArray::SArray()
		:len_(0),sarry_(0),allocated_(false),initialized_(false)
	{
			SetSArray(new SparseArrayRep());
	}

	SArray::SArray( int data_len )
		:len_(data_len),sarry_(0),allocated_(true),initialized_(false)
	{
		SetSArray(new SparseArrayRep(data_len));
	}

	SArray::SArray( int* rows,int* cols,double* values,int data_len,bool has_initialized_ )
		:len_(0),sarry_(0),allocated_(true),initialized_(has_initialized_)
	{
		SetSArray(new SparseArrayRep(rows,cols,values,data_len));
		len_=data_len;
	}

	SArray::SArray( const SArray& copy )
		:sarry_(0),len_(copy.GetLength())
	{
		if (NotEqual(copy) )
		{
			ReleaseSArray();
			sarry_=copy.GetSArrayMutable();
			sarry_->AddRef(this);
			allocated_=copy.CheckAllocated();
			initialized_=copy.CheckInitialized();
			len_=copy.GetLength();
		}
		len_=copy.GetLength();
	}

	SArray& SArray::operator=( const SArray& equal_sarray )
	{
		
		if (NotEqual(equal_sarray) )
		{
			ReleaseSArray();
			sarry_=equal_sarray.GetSArrayMutable();
			sarry_->AddRef(this);
			allocated_=equal_sarray.CheckAllocated();
			initialized_=equal_sarray.CheckInitialized();
            len_=equal_sarray.GetLength();
		}
		
		return *this;
	}

	void SArray::MakeUnique()
	{
		LP_ASSERT_EXCEPTION(CheckAllocated(),SARRAY_ERROR,"has't been allocated");
		if (GetSArrayConst()->RefCount(this)>1)
		{
			bool initialized_old=initialized_;
			bool allocated_old=allocated_;
			int len_old=len_;
			SetSArray(new SparseArrayRep(*(this->GetSArrayConst())));
			initialized_=initialized_old;
			allocated_=allocated_old;
			len_=len_old;
		}
	}

	void SArray::Copy( const SArray& copy_sarray )
	{
		SetSArray(new SparseArrayRep(*(copy_sarray.GetSArrayMutable())));
		len_=copy_sarray.GetLength();
		initialized_=copy_sarray.CheckInitialized();
		allocated_=copy_sarray.CheckAllocated();
	}

	void SArray::ReleaseSArray()
	{
		if(sarry_){
			if (sarry_->GetDataPtr()) {
				sarry_->ReleaseRef(this);
				if (sarry_->RefCount(this) == 0 ) {
					delete sarry_;
				}
				sarry_ = 0;
				allocated_=false;
				initialized_=false;
				len_=0;
			}
		}
	}

	void SArray::SetSArray( SparseArrayRep* new_sarray_ )
	{
		ReleaseSArray();//Release any old ArrayRep
		new_sarray_->AddRef(this);
		sarry_=new_sarray_;
	}

	bool SArray::NotEqual( const SArray& compare_sarray )
	{
		if ((allocated_!=compare_sarray.CheckAllocated())||
			(initialized_!=compare_sarray.CheckInitialized())
			||(sarry_!=compare_sarray.GetSArrayConst()))
		{
			return true;
		}

		if (!HasSameSArray(compare_sarray))
		{
			return true;
		}

		return false;
	}

	bool SArray::HasSameSArray( const SArray& compare_sarray )
	{
		if (GetSArrayConst()->length()!=compare_sarray.GetSArrayConst()->length())
		{
			return false;
		}

		int nelem=GetSArrayConst()->length();
		int* row_lhs=GetSArrayConst()->GetRowsPtr();
		int* row_rhs=compare_sarray.GetSArrayConst()->GetRowsPtr();
		
		for(long i=0;i<nelem;i++)
		{
			if (row_lhs[i]!=row_rhs[i]){return false;}
		}

		int* col_lhs=GetSArrayConst()->GetColsPtr();
		int* col_rhs=compare_sarray.GetSArrayConst()->GetColsPtr();

		for(long i=0;i<nelem;i++)
		{
			if (col_lhs[i]!=col_rhs[i]){return false;}
		}

		double* data_lhs=GetSArrayConst()->GetDataPtr();
		double* data_rhs=compare_sarray.GetSArrayConst()->GetDataPtr();

		for(long i=0;i<nelem;i++)
		{
			if (data_lhs[i]!=data_rhs[i]){return false;}
		}
		return true;
	}

}//namespace Lpopx