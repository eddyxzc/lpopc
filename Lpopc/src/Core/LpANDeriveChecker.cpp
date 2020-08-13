// Copyright (C) 2014-2015 Xue Zhichen, Wang Yujie,Wang Na 
// All Rights Reserved.
// This file is a part of LPOPC , published under the Eclipse Public License.
// Author:Xue Zhichen 8/4 2015   21:40
// Email:eddy_lpopc@163.com
#include "LpANDeriveChecker.h"
#include <string>
#include <sstream>
namespace Lpopc
{


	void LpANDeriveChecker::CheckeAnlyticlDerive()
	{
		LP_DBG_START_FUN("LpANDeriveChecker::CheckeAnlyticlDerive()")
		reporter_->Printf()->info("\n");
		reporter_->Printf()->info( "Checking user defined analytic derivatives against finite difference \n");

		std::vector<SolCost> solcCostTotal(optpro_->GetPhaseNum());
		vec x_all = Data_->nlpGuessVector;
		for (size_t iphase = 0; iphase < Data_->numphases_; iphase++)
		{
			//Get the guess in each phase of the problem
			int sumnodes = Data_->PS[iphase]->Points.n_elem;
			size_t nstates = Data_->SIZES_[iphase][0];
			size_t ncontrols = Data_->SIZES_[iphase][1];
			size_t nparameters = Data_->SIZES_[iphase][2];
			size_t npaths = Data_->SIZES_[iphase][3];
			size_t nevents = Data_->SIZES_[iphase][4];
			int stateindex_start = Data_->phase_indices[iphase]->state[0] - 1;
			int stateindex_end = Data_->phase_indices[iphase]->state[Data_->phase_indices[iphase]->state.size() - 1] - 1;
			int controlindex_start = Data_->phase_indices[iphase]->control[0] - 1;
			int controlindex_end = Data_->phase_indices[iphase]->control[Data_->phase_indices[iphase]->control.size() - 1] - 1;
			vec state_vector = x_all.subvec(stateindex_start, stateindex_end);
			vec control_vector = x_all.subvec(controlindex_start, controlindex_end);
			int t0_index = Data_->phase_indices[iphase]->time[0] - 1;
			int tf_index = Data_->phase_indices[iphase]->time[1] - 1;
			double t0 = x_all(t0_index);
			double tf = x_all(tf_index);
			double tspan = tf - t0;
			vec t_radau = (Data_->PS[iphase]->Points + 1)*(tspan / 2.0) + t0;
			mat state_matrix = reshape(state_vector, sumnodes + 1, nstates);
			mat state_radau = state_matrix.rows(0, state_matrix.n_rows - 2);//index start form 0
			vec x0 = trans(state_matrix.row(0));
			vec xf = trans(state_matrix.row(state_matrix.n_rows - 1));
			mat control_radau = reshape(control_vector, sumnodes, ncontrols);
			mat parameters;
			if (nparameters > 0)
			{
				int para_index_start = Data_->phase_indices[iphase]->parameter[0] - 1;
				int para_index_end = Data_->phase_indices[iphase]->parameter[Data_->phase_indices[iphase]->parameter.size() - 1] - 1;
				parameters = x_all.rows(para_index_start, para_index_end);
			}
			
			SolCost mysolcost;
			mysolcost.initial_time_ = t0;
			mysolcost.initial_state_ = x0;
			mysolcost.terminal_time_ = tf;
			mysolcost.terminal_state_ = xf;
			mysolcost.time_ = t_radau;
			mysolcost.state_ = state_radau;
			mysolcost.control_ = control_radau;
			mysolcost.parameter_ = parameters;
			mysolcost.phase_num_ = iphase + 1;
			//////////////////////////////////////////////////////////////////////////
			//check user cost function
			rowvec dANmayer;
			ANderive_->DerivMayer(mysolcost, dANmayer);
			if (dANmayer.n_cols!=nstates*2+nparameters+2)
			{
				std::string errormsg = "User's Analytical derivative of Mayer Cost returned invalid size of row vector in Phase%u";
				reporter_->Printf()->warn( errormsg.c_str(),iphase+1);
			}
			mat dANlagrange;
			ANderive_->DerivLagrange(mysolcost, dANlagrange);
			if (dANlagrange.n_rows != t_radau.n_rows)
			{
				std::string errormsg = "User's Analytical derivative of Lagrange Cost returned invalid number of rows  in Phase";
				errormsg += iphase + 1;
				errormsg += ",which should equal to the number of elements of time";
				reporter_->Printf()->warn(errormsg.c_str());

			}
			if (dANlagrange.n_cols!=nstates+ncontrols+nparameters+1)
			{
				std::string errormsg = "User's Analytical derivative of Lagrange Cost returned invalid number of columns in Phase%u";
				reporter_->Printf()->warn(errormsg.c_str(), iphase + 1);
			}
			//estimate cost derive using finite difference
			rowvec dFDmayer;
			FDderive_->DerivMayer(mysolcost, dFDmayer);
			rowvec mayer_error=abs(dFDmayer - dANmayer);
			reporter_->Printf()->warn("Phase %u, Derive Mayer Cost:\n", iphase + 1);
			for (size_t istate = 0; istate < nstates; istate++)
			{
				double abs_error = mayer_error(istate);
				if (abs_error>tol_)
				{
					std::string msg = "dMayer/dx0(%u): \tuser=%16f,\t finite difference =%16f,\t error=%16f";
					reporter_->Printf()->warn(msg.c_str(),
						istate, dANmayer(istate), dFDmayer(istate),mayer_error(istate));
				}
				
			}
			if (mayer_error(nstates)>tol_)
			{
				reporter_->Printf()->warn("dMayer/dt0: \tuser=%16f,\t finite difference =%16f,\t error=%16f\n",
					dANmayer(nstates), dFDmayer(nstates), mayer_error(nstates));
			}
			
			for (size_t istate = 0; istate < nstates; istate++)
			{
				if (mayer_error(istate + nstates + 1)>tol_)
				{
					std::string msg = "dMayer/dxf(%u): \tuser=%16f,\t finite difference =%16f\t error=%16f\n";
					reporter_->Printf()->warn(msg.c_str(),
						istate, dANmayer(istate + nstates + 1), dFDmayer(istate+1+nstates), mayer_error(istate + nstates + 1));
				}
				
			}
			if (mayer_error(2 * nstates + 1)>tol_)
			{
				reporter_->Printf()->warn("dMayer/dtf: \tuser=%16f,\t finite difference =%16f,\t error=%16f\n",
					dANmayer(2 * nstates + 1), dFDmayer(2 * nstates + 1), mayer_error(2 * nstates + 1));
			}
			
			if (nparameters>0)
			{
				for (size_t ipara = 0; ipara < nparameters; ipara++)
				{
					if (mayer_error(ipara + 2 * nstates + 2)>tol_)
					{
						std::string msg = "dMayer/dpara(%u): \tuser=%16f,\t finite difference =%16f,\t error=%16f\n";
						reporter_->Printf()->warn(msg.c_str(),
							ipara, dANmayer(ipara + 2 * nstates + 2), dFDmayer(ipara + 2 * nstates + 2), mayer_error(ipara + 2 * nstates + 2));
					}
					
				}
			}
			//////////////////////////////////////////////////////////////////////////
			// check lagrange cost
			reporter_->Printf()->warn("Phase %u, Derive Lagrange Cost:\n", iphase + 1);
			mat dFDlagrange;
			FDderive_->DerivLagrange(mysolcost, dFDlagrange);
			mat lagrange_error = abs(dANlagrange - dFDlagrange);
			double abs_error=0.0;
			for (size_t istate = 0; istate < nstates; istate++)
			{
				for (size_t j = 0; j < sumnodes; j++)
				{
					abs_error = lagrange_error(j, istate);// major column
					if (abs_error>tol_)
					{
						std::string msg = "dLagrange/dx(%u) elem(%u): \tuser=%16f,\t finite difference =%16f,\t error=%16f\n";
						reporter_->Printf()->warn(msg.c_str(),
							istate,j,dANlagrange(j,istate), dFDlagrange(j,istate), abs_error);
					}
				}
			}

			for (size_t icontrol = 0; icontrol < ncontrols; icontrol++)
			{
				for (size_t j = 0; j < sumnodes; j++)
				{
					abs_error = lagrange_error(j, nstates+icontrol);// major column
					if (abs_error > tol_)
					{
						std::string msg = "dLagrange/du(%u) elem(%u): \tuser=%16f,\t finite difference =%16f,\t error=%16f\n";
						reporter_->Printf()->warn(msg.c_str(),
							icontrol, j, dANlagrange(j, nstates+icontrol), dFDlagrange(j, nstates+icontrol), abs_error);
					}
				}
			}// for icontrol

			for (size_t j = 0; j < sumnodes; j++)
			{
				abs_error = lagrange_error(j, nstates + ncontrols);// major column
				if (abs_error > tol_)
				{
					std::string msg = "dLagrange/dt elem(%u): \tuser=%16f,\t finite difference =%16f,\t error=%16f\n";
					reporter_->Printf()->warn(msg.c_str(),
						 j, dANlagrange(j, nstates + ncontrols), dFDlagrange(j, nstates + ncontrols), abs_error);
				}
			}
			for (size_t ipara = 0; ipara < nparameters; ipara++)
			{
				for (size_t j = 0; j < sumnodes; j++)
				{
					abs_error = lagrange_error(j, nstates + ncontrols+1+ipara);// major column
					if (abs_error > tol_)
					{
						std::string msg = "dLagrange/dp(%u) elem(%u): \tuser=%16f,\t finite difference =%16f,\t error=%16f\n";
						reporter_->Printf()->warn(msg.c_str(),
							ipara, j, dANlagrange(j, nstates + ncontrols+1+ipara), dFDlagrange(j, nstates + ncontrols+1+ipara), abs_error);
					}
				}
			}// for icontrol

			//////////////////////////////////////////////////////////////////////////
			//check dae derive,check path derive
			SolDae mySolDae;
			mySolDae.time_ = t_radau;
			mySolDae.state_ = state_radau;
			mySolDae.contol_ = control_radau;
			mySolDae.parameter_ = parameters;
			mySolDae.phase_num_ = iphase + 1;
			mat dANdiff, dANpath;
			reporter_->Printf()->warn("Phase %u,  check difference equation sizes:\n", iphase + 1);
			ANderive_->DerivDae(mySolDae, dANdiff, dANpath);
			
			if (dANdiff.n_rows!=sumnodes*nstates)
			{
				std::string errormsg = "User's Analytical derivative of different constrains returned invalid number of rows  in Phase";
				errormsg += iphase + 1;
				reporter_->Printf()->warn(errormsg.c_str());
			}
			if (dANdiff.n_cols!=nstates+ncontrols+nparameters+1)
			{
				std::string errormsg = "User's Analytical derivative of different constrains returned invalid number of columns  in Phase";
				errormsg += iphase + 1;
				reporter_->Printf()->warn(errormsg.c_str());
			}
			if (npaths>0)
			{
				if (dANpath.n_rows != sumnodes*npaths)
				{
					std::string errormsg = "User's Analytical derivative of path constrains returned invalid number of rows  in Phase %u";
					reporter_->Printf()->warn(errormsg.c_str(), iphase + 1);
				}
				if (dANpath.n_cols != nstates + ncontrols + nparameters + 1)
				{
					std::string errormsg = "User's Analytical derivative of path constrains returned invalid number of columns  in Phase";
					errormsg += iphase + 1;
					reporter_->Printf()->warn(errormsg.c_str());
				}
			}
			reporter_->Printf()->warn("Phase %u,   difference equation sizes check down \n", iphase + 1);
			reporter_->Printf()->warn("Phase %u,  check difference equation value:\n", iphase + 1);
			mat dFDdiff, dFDpath;
			FDderive_->DerivDae(mySolDae, dFDdiff, dFDpath);
			mat diff_error = abs(dANdiff - dFDdiff);
			mat path_error = abs(dANpath - dFDpath);
			for (size_t jstate = 0; jstate < nstates; jstate++)
			{
				for (size_t idiff = 0; idiff < nstates; idiff++)
				{
					for (size_t inodes = 0; inodes < sumnodes; inodes++)
					{
						abs_error = diff_error(idiff*nstates + inodes, jstate);
						if (abs_error > tol_)
						{
							std::string msg = "dDae(%u)/dx(%u) elem(%u): \tuser=%16f,\t finite difference =%16f,\t error=%16f\n";
							reporter_->Printf()->warn(msg.c_str(),
								idiff, jstate, inodes, dANdiff(idiff*nstates + inodes, jstate), dFDdiff(idiff*nstates + inodes, jstate), abs_error);
						}
					}
				}
			}// for jstate
			for (size_t jcontrol = 0; jcontrol < ncontrols; jcontrol++)
			{
				for (size_t idiff = 0; idiff < nstates; idiff++)
				{
					for (size_t inodes = 0; inodes < sumnodes; inodes++)
					{
						abs_error = diff_error(idiff*nstates + inodes, jcontrol+nstates);
						if (abs_error > tol_)
						{
							std::string msg = "dDae(%u)/du(%u) elem(%u): \tuser=%16f,\t finite difference =%16f,\t error=%16f\n";
							reporter_->Printf()->warn(msg.c_str(),
								idiff, jcontrol, inodes, dANdiff(idiff*nstates + inodes, jcontrol + nstates), dFDdiff(idiff*nstates + inodes, jcontrol + nstates), abs_error);
						}
					}
				}
			}// for jcontrol
			for (size_t idiff = 0; idiff < nstates; idiff++)
			{
				for (size_t inodes = 0; inodes < sumnodes; inodes++)
				{
					abs_error = diff_error(idiff*nstates + inodes, ncontrols + nstates);
					if (abs_error > tol_)
					{
						std::string msg = "dDae(%u)/dt elem(%u): \tuser=%16f,\t finite difference =%16f,\t error=%16f\n";
						reporter_->Printf()->warn(msg.c_str(),
							idiff,  inodes, dANdiff(idiff*nstates + inodes, ncontrols + nstates ), dFDdiff(idiff*nstates + inodes, ncontrols + nstates ), abs_error);
					}
				}
			}
			for (size_t jpara = 0; jpara < nparameters; jpara++)
			{
				for (size_t idiff = 0; idiff < nstates; idiff++)
				{
					for (size_t inodes = 0; inodes < sumnodes; inodes++)
					{
						abs_error = diff_error(idiff*nstates + inodes, nstates+ncontrols+1+jpara);
						if (abs_error > tol_)
						{
							std::string msg = "dDae(%u)/dp(%u) elem(%u): \tuser=%16f,\t finite difference =%16f,\t error=%16f\n";
							reporter_->Printf()->warn(msg.c_str(),
								idiff, jpara, inodes, dANdiff(idiff*nstates + inodes, nstates + ncontrols + 1 + jpara), dFDdiff(idiff*nstates + inodes, nstates + ncontrols + 1 + jpara), abs_error);
						}
					}
				}
			}// for jpara
			reporter_->Printf()->warn("Phase {},   difference equation values check down \n", iphase + 1);
			//////////////////////////////////////////////////////////////////////////
			if (npaths>0)
			{
				reporter_->Printf()->warn("Phase {},   path constraints values checking... \n", iphase + 1);
				for (size_t jstate = 0; jstate < nstates; jstate++)
				{
					for (size_t ipath = 0; ipath < npaths; ipath++)
					{
						for (size_t inodes = 0; inodes < sumnodes; inodes++)
						{
							abs_error = path_error(ipath*nstates + inodes, jstate);
							if (abs_error > tol_)
							{
								std::string msg = "dPath(%u)/dx(%u) elem(%u): \tuser=%16f,\t finite difference =%16f,\t error=%16f\n";
								reporter_->Printf()->warn(msg.c_str(),
									ipath, jstate, inodes, dANpath(ipath*nstates + inodes, jstate), dFDpath(ipath*nstates + inodes, jstate), abs_error);
							}
						}
					}
				}// for jstate
				for (size_t jcontrol = 0; jcontrol < ncontrols; jcontrol++)
				{
					for (size_t ipath = 0; ipath < npaths; ipath++)
					{
						for (size_t inodes = 0; inodes < sumnodes; inodes++)
						{
							abs_error = path_error(ipath*nstates + inodes, jcontrol + nstates);
							if (abs_error > tol_)
							{
								std::string msg = "dPath(%u)/du(%u) elem(%u): \tuser=%16f,\t finite difference =%16f,\t error=%16f\n";
								reporter_->Printf()->warn(msg.c_str(),
									ipath, jcontrol, inodes, dANpath(ipath*nstates + inodes, jcontrol + nstates), dFDpath(ipath*nstates + inodes, jcontrol + nstates), abs_error);
							}
						}
					}
				}// for jcontrol
				for (size_t ipath = 0; ipath < npaths; ipath++)
				{
					for (size_t inodes = 0; inodes < sumnodes; inodes++)
					{
						abs_error = path_error(ipath*nstates + inodes, ncontrols + nstates );
						if (abs_error > tol_)
						{
							std::string msg = "dPath(%u)/dt elem(%u): \tuser=%16f,\t finite difference =%16f,\t error=%16f\n";
							reporter_->Printf()->warn(msg.c_str(),
								ipath, inodes, dANpath(ipath*nstates + inodes, ncontrols + nstates ), dFDpath(ipath*nstates + inodes, ncontrols + nstates ), abs_error);
						}
					}
				}
				for (size_t jpara = 0; jpara < nparameters; jpara++)
				{
					for (size_t ipath = 0; ipath < npaths; ipath++)
					{
						for (size_t inodes = 0; inodes < sumnodes; inodes++)
						{
							abs_error = path_error(ipath*nstates + inodes, nstates + ncontrols + 1 + jpara);
							if (abs_error > tol_)
							{
								std::string msg = "dPath(%u)/dp(%u) elem(%u): \tuser=%16f,\t finite difference =%16f,\t error=%16f\n";
								reporter_->Printf()->warn(msg.c_str(),
									ipath, jpara, inodes, dANpath(ipath*nstates + inodes, nstates + ncontrols + 1 + jpara), dFDpath(ipath*nstates + inodes, nstates + ncontrols + 1 + jpara));
							}
						}
					}
				}// for jpara
				reporter_->Printf()->warn("Phase %u,   path constrains values check down \n", iphase + 1);
			}
			
			//////////////////////////////////////////////////////////////////////////
			// check event derive
			mat event_error;
			if (nevents>0)
			{
				SolEvent mySolEvent, isolevent, jsolevent, ijsolevent;
				mySolEvent.initial_time_ = t0;
				mySolEvent.initial_state_ = x0;
				mySolEvent.terminal_time_ = tf;
				mySolEvent.terminal_state_ = xf;
				mySolEvent.parameter_ = parameters;
				mySolEvent.phase_num_ = iphase + 1;
				mat dANevent;
				reporter_->Printf()->warn("Phase %u,  check event equation sizes:\n", iphase + 1);
				ANderive_->DerivEvent(mySolEvent, dANevent);
				if (dANevent.n_rows!=nevents)
				{
					std::string errormsg = "User's Analytical derivative of Event constrains returned invalid number of rows  in Phase";
					errormsg += iphase + 1;
					reporter_->Printf()->warn(errormsg.c_str());
				}
				if (dANevent.n_cols != nstates*2+nparameters+2)
				{
					std::string errormsg = "User's Analytical derivative of Event constrains returned invalid number of cols  in Phase";
					errormsg += iphase + 1;
					reporter_->Printf()->warn(errormsg.c_str());
				}
				reporter_->Printf()->warn("Phase %u,   event equation sizes check down\n", iphase + 1);
				reporter_->Printf()->warn("Phase %u,  checking event equation values...\n", iphase + 1);

				mat dFDevent;
				FDderive_->DerivEvent(mySolEvent, dFDevent); 
				event_error = abs(dANevent - dFDevent);
				reporter_->Printf()->warn("Phase %u, Derive Event Constrains:\n", iphase + 1);
				for (size_t ievent = 0; ievent < nevents; ievent++)
				{
					for (size_t jstate = 0; jstate < nstates; jstate++)
					{
						double abs_error = event_error(ievent, jstate);
						if (abs_error > tol_)
						{
							std::string msg = "dEvent(%u)/dx0(%u): \tuser=%16f,\t finite difference =%16f,\t error=%16f\n";
							reporter_->Printf()->warn(msg.c_str(),
								ievent, jstate, dANevent(ievent, jstate), dFDevent(ievent, jstate), abs_error);
						}

					}


					if (event_error(ievent,nstates) > tol_)
					{
						reporter_->Printf()->warn("dEvent(%u)/dt0: \tuser=%16f,\t finite difference =%16f,\t error=%16f\n",
							ievent,dANevent(ievent,nstates), dFDevent(ievent,nstates), event_error(ievent,nstates));
					}
					for (size_t jstate = 0; jstate < nstates; jstate++)
					{
						abs_error = event_error(ievent, nstates + jstate + 1);
						if (abs_error>tol_)
						{
							std::string msg = "dEvent(%u)/dxf(%u): \tuser=%16f,\t finite difference =%16f\t error=%16f\n";
							reporter_->Printf()->warn(msg.c_str(),
								ievent,jstate, dANevent(ievent,jstate + nstates + 1), dFDevent(ievent,jstate + 1 + nstates), abs_error);
						}

					}
					if (event_error(2 * nstates + 1) > tol_)
					{
						reporter_->Printf()->warn("dEvent(%u)/dtf: \tuser=%16f,\t finite difference =%16f,\t error=%16f\n",
							ievent,dANevent(ievent,2 * nstates + 1), dFDevent(ievent,2 * nstates + 1), event_error(2 * nstates + 1));
					}

					if (nparameters > 0)
					{
						for (size_t ipara = 0; ipara < nparameters; ipara++)
						{
							abs_error = event_error(ievent, nstates * 2 + 2 + ipara);
							if (abs_error > tol_)
							{
								std::string msg = "dEvent(%u)/dpara(%u): \tuser=%16f,\t finite difference =%16f,\t error=%16f\n";
								reporter_->Printf()->warn( msg.c_str(),
									ievent,ipara, dANevent(ievent,ipara + 2 * nstates + 2), dFDevent(ievent,ipara + 2 * nstates + 2), abs_error);
							}

						}
					}
				}// for ievent
				reporter_->Printf()->warn("Phase %u,   event function values check down \n", iphase + 1);
			}// if events
			solcCostTotal[iphase] = mysolcost;
		}// for iphase
		//////////////////////////////////////////////////////////////////////////
		//check linkage derive
		size_t numlinkpairs = Data_->numlinkpairs_;
		for (size_t ipair = 0; ipair < numlinkpairs; ipair++)
		{
			size_t nlinks = optpro_->GetLinkage(ipair)->GetLinkageMin().size();
			int left_index = optpro_->GetLinkage(ipair)->LeftPhase();
			int right_index = optpro_->GetLinkage(ipair)->RightPhase();
			vec xf_left = solcCostTotal[left_index].terminal_state_;
			mat p_left = solcCostTotal[left_index].parameter_;
			
			vec x0_right = solcCostTotal[right_index].initial_state_;
			mat p_right = solcCostTotal[right_index].parameter_;

			SolLink mySolLink;
			mySolLink.left_state_ = xf_left;
			mySolLink.left_parameter_ = p_left;
			mySolLink.left_phase_num_ = left_index + 1;

			mySolLink.right_state_ = x0_right;
			mySolLink.right_parameter_ = p_left;
			mySolLink.right_phase_num_ = right_index + 1;

			mySolLink.ipair = ipair + 1;

			int nstateLeft = mySolLink.left_state_.n_elem;
			int nparametersLeft = mySolLink.left_parameter_.n_elem;
			int nstateRight = mySolLink.right_state_.n_elem;
			int nparametersRight = mySolLink.right_parameter_.n_elem;

			//////////////////////////////////////////////////////////////////////////
			reporter_->Printf()->warn("Linkage %u,  check linkage equation sizes:\n", ipair + 1);
			mat dANlink;
			ANderive_->DerivLink(mySolLink, dANlink);
			if (dANlink.n_rows!=nlinks)
			{
				std::string errormsg = "User's Analytical derivative of Link Function returned invalid number of rows  in Linkage %u";
				reporter_->Printf()->warn(errormsg.c_str(), ipair + 1);
			}
			if (dANlink.n_cols != nstateLeft + nstateRight + nparametersLeft + nparametersLeft)
			{
				std::string errormsg = "User's Analytical derivative of Link Function returned invalid number of columns  in Linkage %u";
				reporter_->Printf()->warn(errormsg.c_str(), ipair + 1);
			}
			reporter_->Printf()->warn("Linkage %u,   linkage equation sizes check down\n", ipair + 1);
			reporter_->Printf()->warn("Linkage %u,  check linkage equation values:\n", ipair + 1);
			mat dFDlink;
			FDderive_->DerivLink(mySolLink, dFDlink);
			mat link_error = abs(dANlink - dFDlink);
			reporter_->Printf()->warn("Linkage %u, Linkage Constrains:\n", ipair);
			double abs_error = 0;
			for (size_t istate = 0; istate < nstateLeft; istate++)
			{
				for (size_t ilink = 0; ilink < nlinks; ilink++)
				{
					abs_error = link_error(ilink, istate);
					if (abs_error>tol_)
					{
						std::string msg = "dLink(%u)/dxfL(%u): \tuser=%16f,\t finite difference =%16f,\t error=%16f\n";
						reporter_->Printf()->warn(msg.c_str(),
							ilink, istate, dANlink(ilink, istate), dFDlink(ilink, istate), abs_error);
					}
				}
			}
			for (size_t ipara = 0; ipara < nparametersLeft; ipara++)
			{
				for (size_t ilink = 0; ilink < nlinks; ilink++)
				{
					abs_error = link_error(ilink, ipara+nstateLeft);
					if (abs_error > tol_)
					{
						std::string msg = "dLink(%u)/dpL(%u): \tuser=%16f,\t finite difference =%16f,\t error=%16f\n";
						reporter_->Printf()->warn(msg.c_str(),
							ilink, ipara, dANlink(ilink, ipara), dFDlink(ilink, ipara), abs_error);
					}
				}
			}
			//////////////////////////////////////////////////////////////////////////
			for (size_t istate = 0; istate < nstateRight; istate++)
			{
				for (size_t ilink = 0; ilink < nlinks; ilink++)
				{
					abs_error = link_error(ilink, istate+nstateLeft+nparametersLeft);
					if (abs_error > tol_)
					{
						std::string msg = "dLink(%u)/dxR0(%u): \tuser=%16f,\t finite difference =%16f,\t error=%16f\n";
						reporter_->Printf()->warn(msg.c_str(),
							ilink, istate, dANlink(ilink, istate+nstateLeft+nparametersLeft), dFDlink(ilink, istate+nstateLeft+nparametersLeft), abs_error);
					}
				}
			}
			for (size_t ipara = 0; ipara < nparametersRight; ipara++)
			{
				for (size_t ilink = 0; ilink < nlinks; ilink++)
				{
					abs_error = link_error(ilink, ipara + nstateLeft+nstateRight+nparametersLeft);
					if (abs_error > tol_)
					{
						std::string msg = "dLink(%u)/dpR(%u): \tuser=%16f,\t finite difference =%16f,\t error=%16f\n";
						reporter_->Printf()->warn(msg.c_str(),
							ilink, ipara, dANlink(ilink, ipara + nstateLeft + nstateRight + nparametersLeft), dFDlink(ilink, ipara + nstateLeft+nstateRight+nparametersLeft), abs_error);
					}
				}
			}
			reporter_->Printf()->warn("Linkage %u,   linkage equation values check down\n", ipair + 1);
		}//for ipair
		reporter_->Printf()->warn("Analytic Derive Check Down! Error tol =%f\n", tol_);
	}

}// namespace Lpopc