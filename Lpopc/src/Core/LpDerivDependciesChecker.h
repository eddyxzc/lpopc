// Copyright (C) 2014-2015 Xue Zhichen, Wang Yujie,Wang Na 
// All Rights Reserved.
// This file is a part of LPOPC , published under the Eclipse Public License.
// Author:Xue Zhichen 7/12 2015   21:26
// Email:eddy_lpopc@163.com
#ifndef LPOPC_DERIVDEPENDCIESCHECKER_H
#define LPOPC_DERIVDEPENDCIESCHECKER_H
#include "LpConf.h"
#include "LpFunctionWrapper.h"
#include "LpException.hpp"
#include "LpCalculateData.hpp"
#include "LpDebug.hpp"
#include <vector>
namespace Lpopc
{
	LP_DECLARE_EXCEPTION(DERIVECHECKER_ERROR);
	///Check jacobi dependicies use finite-nan
	// considering a F(x,y,z)=0,we set x=nan,if F=nan,then F dependes on x,df/dx should be calculate.
	// The jacobi of RPM problem is sparse,most of its elements are zero,we only tell the NLP solver which
	// of them are nonzero,this will accelerate the calculation
	class DeriveDependicieshecker
	{
	public:

		DeriveDependicieshecker(shared_ptr<FunctionWrapper> funs, shared_ptr<LpCalculateData> datas)
			:funs_(funs), Data_(datas)
		{
			GetDependiciesForJacobiInEveryPhase();

		};
		void GetDependiciesForJacobi( std::vector<umat>& dependicies)const
		{
			LP_DBG_START_FUN("DeriveDependicieshecker::GetDependiciesForJacobi()")
			dependicies = Jocobidependices;
		};

		~DeriveDependicieshecker(){};

	private:
		DeriveDependicieshecker(){};
		DeriveDependicieshecker& operator=(const DeriveDependicieshecker&);
		DeriveDependicieshecker(const DeriveDependicieshecker&);

		void GetDependiciesForJacobiInEveryPhase();
		std::vector<umat> Jocobidependices;
		shared_ptr<FunctionWrapper> funs_;
		shared_ptr<LpCalculateData> Data_;
	};
}



#endif // !LPOPC_DERIVDEPENDCIESCHECKER_H
