// Copyright (C) 2014-2015 Xue Zhichen, Wang Yujie,Wang Na 
// All Rights Reserved.
// This file is a part of LPOPC , published under the Eclipse Public License.
// Author:Xue Zhichen 7/15 2015   11:18
// Email:eddy_lpopc@163.com

#include "BrysonDenham.h"
#include "LpLpopcApplication.hpp"
int main()
{
	double x10 = 0;
	double x20 = -10;
	double x30 = 0;
	double x1f = 0;
	double x2f = -1;
	double x1min = 0;
	double x1max = 1.0 / 9;
	double x2min = -10;
	double x2max = 10;
	double x3min = -10;
	double x3max = 10;
	//Step1
	shared_ptr<LpopcApplication> app(new LpopcApplication(console_print));
	//Step2
	shared_ptr<Phase> phase1(new Phase(1, 3, 1, 0, 0, 5));
	phase1->SetTimeMin(0.0, 0.0);
	phase1->SetTimeMax(0, 50);

	phase1->SetStateMin(0, 0, 0);
	phase1->SetStateMax(1.0/9.0, 1.0/9.0, 1.0/9.0);

	phase1->SetStateMin(-10, -10, -10);
	phase1->SetStateMax(10, 10, 10);

	phase1->SetStateMin(-10, -10, -10);
	phase1->SetStateMax(10, 10, 10);

	phase1->SetcontrolMin(-10);
	phase1->SetcontrolMax(10);

	phase1->SetTimeGuess(0.0);
	phase1->SetTimeGuess(1.0);

	phase1->SeteventMin(0);
	phase1->SeteventMin(1);
	phase1->SeteventMin(0);
	phase1->SeteventMin(0);
	phase1->SeteventMin(-1);

	phase1->SeteventMax(0);
	phase1->SeteventMax(1);
	phase1->SeteventMax(0);
	phase1->SeteventMax(0);
	phase1->SeteventMax(-1);

	phase1->SetStateGuess(1, 0);// count from 1
	phase1->SetStateGuess(1, 0);

	phase1->SetStateGuess(2, 1.0);
	phase1->SetStateGuess(2, -1.0);

	phase1->SetStateGuess(3, 0.0);
	phase1->SetStateGuess(3, 0.0);

	phase1->SetControlGuess(1, 0.0);
	phase1->SetControlGuess(1, 0.0);

	//Step3
	shared_ptr<FunctionWrapper> userfun(new BrysonDenhamFunction());
	shared_ptr<OptimalProblem> optpro(new OptimalProblem(1, 0, userfun));
	optpro->AddPhase(phase1);
	
		try
	{
		//Step4
		app->SetOptimalControlProblem(optpro);
		app->Options()->SetStringValue("hessian-approximation", "exact");
		app->Options()->SetIntegerValue("max-grid-num", 20);
		app->SolveOptimalProblem();
	}
	catch (LpopcException & e)
	{
		std::cout << "got an error" << std::endl;
	}
	catch (std::logic_error& e)
	{
		std::cout << "got an error" << std::endl;
	}
	catch (std::runtime_error& e)
	{
		std::cout << "got an error" << std::endl;
	}
	catch (std::bad_alloc& e)
	{
		std::cout << "got an error" << std::endl;
	}
	return 0;
}

void BrysonDenhamFunction::MayerCost(SolCost&mySolcost, double& mayer)
{
	vec xf = mySolcost.terminal_state_;
	mayer = xf(2);
}

void BrysonDenhamFunction::LagrangeCost(SolCost& mySolcost, vec& langrange)
{
	vec t = mySolcost.time_;
	langrange = zeros(t.n_elem);

}

void BrysonDenhamFunction::DaeFunction(SolDae& mySolDae, mat& stateout, mat& pathout)
{
	vec t = mySolDae.time_;
	mat x = mySolDae.state_;
	mat u = mySolDae.contol_;
	mat p = mySolDae.parameter_;
	size_t iphase = mySolDae.phase_num_;
	assert(iphase == 1);
	stateout = zeros(x.n_rows, x.n_cols);
	stateout.col(0) = x.col(1);
	stateout.col(1) = u;
	stateout.col(2) = 0.5*(u%u);

}

void BrysonDenhamFunction::DerivMayer(SolCost&mySolcost, rowvec& deriv_mayer)
{
}

void BrysonDenhamFunction::DerivLagrange(SolCost& mySolcost, mat& deriv_langrange)
{
}

void BrysonDenhamFunction::DerivDae(SolDae& mySolDae, mat& deriv_state, mat& deriv_path)
{
}

void BrysonDenhamFunction::EventFunction(SolEvent& mySolEvent, vec& eventout)
{
	vec x0 = mySolEvent.initial_state_;
	vec xf = mySolEvent.terminal_state_;
	eventout = zeros(5);
	double x10 = x0(0);
	double x20 = x0(1);
	double x30 = x0(2);
	double x1f = xf(0);
	double x2f = xf(1);
	eventout(0) = x10;
	eventout(1) = x20;
	eventout(2) = x30;
	eventout(3) = x1f;
	eventout(4) = x2f;
}

void BrysonDenhamFunction::DerivEvent(SolEvent& mySolEvent, mat& deriv_event)
{
}

void BrysonDenhamFunction::LinkFunction(SolLink& mySolLink, vec& linkageout)
{
}

void BrysonDenhamFunction::DerivLink(SolLink& mySolLink, mat& derive_link)
{
}
