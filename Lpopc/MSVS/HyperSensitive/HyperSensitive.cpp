// Copyright (C) 2014-2015 Xue Zhichen, Wang Yujie,Wang Na 
// All Rights Reserved.
// This file is a part of LPOPC , published under the Eclipse Public License.
// Author:Xue Zhichen 7/14 2015   13:16
// Email:eddy_lpopc@163.com

///only need including two headfils
#include "hypersensitive.h"
#include "LpLpopcApplication.hpp"


//////////////////////////////////////////////////////////////////////////
using namespace Lpopc;
int main(void)
{
	
	double t0 = 0.0, tf = 5000, x0 = 1.5, xf = 1;
	double xmin = -10, xmax = 10, umin = -10, umax = 10;


	//Step 1
	shared_ptr<LpopcApplication> app(new LpopcApplication(console_print));

	//Step2
	shared_ptr<Phase> Phase1(new Phase(1, 1, 1, 0, 0, 0));
	//Set bounds
	Phase1->SetTimeMin(t0, tf);
	Phase1->SetTimeMax(t0, tf);
	Phase1->SetStateMin(x0, xmin, xf);
	Phase1->SetStateMax(x0, xmax, xf);
	Phase1->SetcontrolMin(umin);
	Phase1->SetcontrolMax(umax);

	///Set guess
	Phase1->SetTimeGuess(t0);
	Phase1->SetTimeGuess(tf);
	Phase1->SetStateGuess(1, x0);//count from 1
	Phase1->SetStateGuess(1, xf);//count from 1
	Phase1->SetControlGuess(1, -1);
	Phase1->SetControlGuess(1, 1);

	//Step3
	shared_ptr<FunctionWrapper> userfun(new HyperSensitiveFunction());//declare optimal problem function
	shared_ptr<OptimalProblem> optpro(new OptimalProblem(1, 0, userfun));
	optpro->AddPhase(Phase1);

	
	LP_DBG_START_FUN("main")
	try
	{
		//Step 4
		app->SetOptimalControlProblem(optpro);
		app->Options()->SetStringValue("hessian-approximation", "exact"); //use sparse-finite hessian
		app->Options()->SetStringValue("first-derive", "analytic");//use analytic derive
 
		app->Options()->SetIntegerValue("max-grid-num", 20);  //maximum of mesh refinement

		app->SolveOptimalProblem();// Solve the problem
	}
	catch (LpopcException & e)
	{
		std::cout << "got an error" << std::endl;
	}
	catch (...)
	{
		std::cout << "got an error" << std::endl;
	}
	
	return 0;
}


void HyperSensitiveFunction::MayerCost(SolCost&mySolcost, double& mayer)
{
	vec x0 = mySolcost.initial_state_;
	vec xf = mySolcost.terminal_state_;
	double t0 = mySolcost.initial_time_;
	double tf = mySolcost.terminal_time_;
	mat p = mySolcost.parameter_;
	int iphase = mySolcost.phase_num_;
	assert(iphase == 1 );

	mayer= 0.0;
}

void HyperSensitiveFunction::DerivMayer(SolCost&mySolcost, rowvec& deriv_mayer)
{
	double t0 = mySolcost.initial_time_;
	double tf = mySolcost.terminal_time_;
	vec x0 = mySolcost.initial_state_;
	vec xf = mySolcost.terminal_state_;
	mat p = mySolcost.parameter_;
	size_t iphase = mySolcost.phase_num_;
	assert(iphase == 1 );
	deriv_mayer = zeros(1, x0.n_elem + 1 + xf.n_elem + 1 + p.n_elem);
}

void HyperSensitiveFunction::LagrangeCost(SolCost& mySolcost, vec& langrange)
{
	
	mat x = mySolcost.state_;
	mat u = mySolcost.control_;
	vec t = mySolcost.time_;
	mat p = mySolcost.parameter_;
	size_t iphase = mySolcost.phase_num_;
	assert(iphase == 1 );
	langrange =0.5*( x%(x) + u%(u));

}

void HyperSensitiveFunction::DerivLagrange(SolCost& mySolcost, mat& deriv_langrange)
{
	mat x = mySolcost.state_;
	mat u = mySolcost.control_;
	vec t = mySolcost.time_;
	mat p = mySolcost.parameter_;
	size_t iphase = mySolcost.phase_num_;
	deriv_langrange = zeros(x.n_rows, x.n_cols + u.n_cols+t.n_cols+p.n_elem);
	deriv_langrange.cols(0, x.n_cols-1) = x;
	deriv_langrange.cols(x.n_cols, x.n_cols + u.n_cols - 1) = u;
}

void HyperSensitiveFunction::DaeFunction(SolDae& mySolDae, mat& stateout, mat& pathout)
{
	vec t = mySolDae.time_;
	mat x = mySolDae.state_;
	mat u = mySolDae.contol_;
	mat p = mySolDae.parameter_;
	size_t iphase = mySolDae.phase_num_;
	assert(iphase == 1 );
	stateout = -x%x%x + u;
		
}

void HyperSensitiveFunction::DerivDae(SolDae& mySolDae, mat& deriv_state, mat& deriv_path)
{
	vec t = mySolDae.time_;
	mat x = mySolDae.state_;
	mat u = mySolDae.contol_;
	mat p=mySolDae.parameter_;

	mat df_dx = -3*(x%x);
	mat df_du = ones(u.n_rows, u.n_cols);
	mat df_dt = zeros(t.n_rows, t.n_cols);
	
	deriv_state = zeros(df_dx.n_rows, df_dx.n_cols + df_du.n_cols + df_dt.n_cols);
	deriv_state.cols(0, df_dx.n_cols - 1) = df_dx;
	deriv_state.cols(df_dx.n_cols, df_dx.n_cols + df_du.n_cols - 1) = df_du;
	deriv_state.cols(df_dx.n_cols + df_du.n_cols, deriv_state.n_cols - 1) = df_dt;
}
void HyperSensitiveFunction::EventFunction(SolEvent& mySolEvent, vec& eventout)
{
}

void HyperSensitiveFunction::DerivEvent(SolEvent& mySolEvent, mat& deriv_event)
{
	
}

void HyperSensitiveFunction::LinkFunction(SolLink& mySolLink, vec& linkageout)
{
}

void HyperSensitiveFunction::DerivLink(SolLink& mySolLink, mat& derive_link)
{
}
