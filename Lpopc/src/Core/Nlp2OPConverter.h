// Copyright (C) 2014-2015 Xue Zhichen, Wang Yujie,Wang Na 
// All Rights Reserved.
// This file is a part of LPOPC , published under the Eclipse Public License.
// Author:Xue Zhichen 7/12 2015   19:11
// Email:eddy_lpopc@163.com
#ifndef LPOPC_NLP2OPCONVERTER_H
#define LPOPC_NLP2OPCONVERTER_H
#include "LpConf.h"
#include "LpOptimalProblem.hpp"
#include "LpFunctionWrapper.h"
#include"LpCalculateData.hpp"
namespace Lpopc
{
	class Nlp2OpConverter
	{
		/** This class converts the data Ipopt returned to the
		 * original optimal problem
		 */
	public:
		Nlp2OpConverter(){
			/*Data_ = Data;
			optpro_ = optpro;*/
		};
		~Nlp2OpConverter(){};
		void Nlp2OpControl(shared_ptr<FunctionWrapper>& Funs_, shared_ptr<LpCalculateData>& Data_, shared_ptr<OptimalProblem>& optpro_);
		void FinalResultSave(shared_ptr<LpCalculateData>& Data_);
		///this function  save the optimal problem result finally.state1 means the state matrix in phase 1.
		///Every column of state matrix is a state.
	private:
		Nlp2OpConverter(const Nlp2OpConverter&);
		Nlp2OpConverter& operator=(const Nlp2OpConverter&);
	};

}//namespace Lpopc
#endif // !LPOPC_NLP2OPCONVERTER_H
