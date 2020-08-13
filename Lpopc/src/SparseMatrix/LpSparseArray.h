// Copyright (C) 2014-2015 Xue Zhichen, Wang Yujie,Wang Na 
// All Rights Reserved.
// This file is a part of LPOPC , published under the Eclipse Public License.
// Author:Xue Zhichen 7/10 2015   17:50
// Email:eddy_lpopc@163.com
#ifndef LPSPARSEARRAY_H
#define LPSPARSEARRAY_H
#include "LpException.hpp"
#include<cassert>
#include<string>
namespace Lpopc
{
	 LP_DECLARE_EXCEPTION(SARRAY_ERROR);
	class SArray
	{
		class SparseArrayRep
		{
		public:

			SparseArrayRep(int* rows,int * cols ,double *d,int l) : irow_(rows),jcol_(cols),data_(d), len (l), count (0) {}

			SparseArrayRep() : irow_(0),jcol_(0),data_(0), len (0),count (0) { }

			SparseArrayRep (int l) : irow_(new int[l]),jcol_(new int[l]),data_ (new double [l]), len (l),count (0) { }

			SparseArrayRep (const SparseArrayRep& copy)
				:irow_(new int[copy.len]),jcol_(new int[copy.len]), data_ (new double [copy.len]), len (copy.len),count (0)
			{
				memcpy(irow_,copy.GetRowsPtr(),sizeof(lp_index)*copy.length() );
				memcpy(jcol_,copy.GetColsPtr(),sizeof(lp_index)*copy.length());
				memcpy(data_,copy.GetDataPtr(),sizeof(double)*copy.length());
			}
			virtual ~SparseArrayRep()
			{
				delete [] data_;
				delete [] irow_;
				delete []jcol_;
			}
			int length () const { return len; }
			int * GetRowsPtr()const{return irow_;};
			int * GetColsPtr()const {return jcol_;};
			double* GetDataPtr ()const { return data_; }

			int RefCount(const SArray*dm)const{return count; }
			void AddRef(const SArray*dm)const {count++;}
			void ReleaseRef(const SArray*dm)const {count--;}
		private:
			SparseArrayRep& operator = (const SparseArrayRep&);
			int *irow_;
			int *jcol_;
			double *data_;
			int len;
			mutable int count;
		};// End of SparseArrayRep

	public:
		SArray();
		SArray(int data_len);
		SArray(int* rows,int* cols,double* values,int data_len,bool has_allocated);
		SArray(const SArray& copy);
		SArray& operator=(const SArray& equal_arry);
		virtual ~SArray(){ReleaseSArray();}
		int GetLength()const {return len_;}
		void MakeUnique();
		//completely copy another SArray , other than  sharing  An SArrayRep with copy_array
		void Copy(const SArray& copy_sarray);
		bool CheckAllocated()const {return allocated_;}

		bool CheckInitialized()const {return initialized_;}

		void ChageToAllocated() {allocated_=true;}

		void ChageToInitialized() {initialized_=true;}
		bool NotEqual(const SArray& compare_sarray);
		const SparseArrayRep* GetSArrayConst()const{return sarry_;}

		SparseArrayRep* GetSArrayMutable()const{
			return sarry_;
		}
	private:
		void ReleaseSArray();
		void SetSArray(SparseArrayRep* new_sarray_);
		bool HasSameSArray(const SArray& compare_sarray);
		bool allocated_;
		bool initialized_;
		int len_;
		SparseArrayRep* sarry_;
	};

}//namespace Lpopc
#endif