// Copyright (C) 2014-2015 Xue Zhichen, Wang Yujie,Wang Na 
// All Rights Reserved.
// This file is a part of LPOPC , published under the Eclipse Public License.
// Author:Xue Zhichen 7/11 2015   17:54
// Email:eddy_lpopc@163.com
#include "LpLpopcAlgorithm.hpp"
#include <sstream>
namespace Lpopc
{


	shared_ptr<LpOptionsList> LpopcAlgorithm::Options()
	{
		return optionlist_;
	}

	Lpopc::ocpstatus LpopcAlgorithm::SolveOptimalControlProblem()
	{
		LP_DBG_START_FUN("LpopcAlgorithm::SolveOptimalControlProblem()")
			try
		{
			SetFirstMesh();
			GetSizes();
			GetBounds();
			GetGuess();
			GetScales();

			GetDependecies();
			if (an_derive_checker_.get())
			{
				an_derive_checker_->CheckeAnlyticlDerive();
			}
			OutputProblemInfo();//Print info about optimal problem such as bounds,options
			SolveNlp();
			Nlp2OpControl();
			while (!RefineMesh())
			{
				UpdateGrid();
				GetSizes();
				GetBounds();
				GetGuess();
				GetScales();
				SolveNlp();
				Nlp2OpControl();
			}
			nlp2op_converter_->FinalResultSave(cd_data_);
			msg_reporter_->Printf()->info("Optimal Problem Solved!Lpopc Exited!");
		}
		catch (MESHREFINE_ERROR&e)
		{
			e.ReportException(*msg_reporter_);
		}
		catch (OPTION_INVALID& e)
		{
			e.ReportException(*msg_reporter_);
		}
		catch (LpopcException& e)
		{
			e.ReportException(*msg_reporter_);
			LP_THROW_EXCEPTION(LPOPC_ALGORITHM_ERROR, "Caught  Lpopc exception");
		}
		catch (std::logic_error& e)
		{
			msg_reporter_->Printf()->error(e.what());
			LP_THROW_EXCEPTION(LPOPC_ALGORITHM_ERROR, "Caught Arma Matrix Lib exception");
		}
		catch (std::runtime_error& e)
		{
			msg_reporter_->Printf()->error(e.what());
			LP_THROW_EXCEPTION(LPOPC_ALGORITHM_ERROR, "Caught Arma Matrix Lib exception");
		}
		catch (std::bad_alloc& e)
		{
			msg_reporter_->Printf()->error("\nEXIT: Not enough memory.\n");
			LP_THROW_EXCEPTION(LPOPC_ALGORITHM_ERROR,
				"Not enough memory");
		}
		catch (...)
		{
			LpopcException exc("Unknown Exception caught in LPOPC", "Unknown File", -1);
			exc.ReportException(*msg_reporter_);
			LP_THROW_EXCEPTION(LPOPC_ALGORITHM_ERROR,
				"Caught unknown exception");
		}
		return solved;
	}

	bool LpopcAlgorithm::RefineMesh()
	{
		LP_DBG_START_FUN("LpopcAlgorithm::RefineMesh()")
			return meshrefiner_->RefineMesh(optpro_);
	}

	void LpopcAlgorithm::Nlp2OpControl()
	{
		LP_DBG_START_FUN("LpopcAlgorithm::Nlp2OpControl()")
			nlp2op_converter_->Nlp2OpControl(optpro_->GetOpimalProblemFuns(), cd_data_, optpro_);
	}

	void LpopcAlgorithm::SolveNlp()
	{
		LP_DBG_START_FUN("LpopcAlgorithm::SolveNlp()")
			msg_reporter_->Printf()->info("Call Ipopt to salve NLP...");
		try{
			nlp_splver_->SolveNlp(nlp_wrapper_, optpro_->GetOpimalProblemFuns(), cd_data_, optpro_);
			msg_reporter_->Printf()->info("NLP problem solved!");
		}
		catch (LPOPC_NLPSOLVER_ERROR e){
			msg_reporter_->Printf()->error(e.Message());
			throw e;
		}
	}

	void LpopcAlgorithm::GetScales()
	{
		LP_DBG_START_FUN("LpopcAlgorithm::GetScales()")
			if (cd_data_->autoscale)
			{
				nlp_wrapper_->GetScales();
			}
	}

	void LpopcAlgorithm::GetDependecies()
	{
		//!Set NLP dependencies 
		LP_DBG_START_FUN("LpopcAlgorithm::GetDependecies()")
			dependicies_checker_.reset(new DeriveDependicieshecker(funs_, cd_data_));
		dependicies_checker_->GetDependiciesForJacobi(cd_data_->allPhaseDependencies);
	}

	void LpopcAlgorithm::GetGuess()
	{
		LP_DBG_START_FUN("LpopcAlgorithm::GetGuess()")
			guess_checker_->GetGuess(optpro_, rpm_, cd_data_);
	}

	void LpopcAlgorithm::GetBounds()
	{
		LP_DBG_START_FUN("LpopcAlgorithm::GetBounds()")
			bounds_checker_->GetBounds(optpro_, cd_data_);
	}

	void LpopcAlgorithm::GetSizes()
	{
		LP_DBG_START_FUN("LpopcAlgorithm::GetSizes()")
			nlp_wrapper_->RefreshSparsity();
		size_checker_->GetSize(optpro_, cd_data_);
	}

	void LpopcAlgorithm::SetFirstMesh()
	{
		LP_DBG_START_FUN("LpopcAlgorithm::SetFirstMesh()")
			meshrefiner_->SetAndCheckMesh(optpro_);
		UpdateGrid();
	}

	void LpopcAlgorithm::Initialized()
	{
		LP_DBG_START_FUN("LpopcAlgorithm::Initialized()");
		std::string first_derive_type, hessain_approximation_type;
		double fd_tol;
		optionlist_->GetStringValue("first-derive", first_derive_type, "");
		optionlist_->GetNumericValue("finite-difference-tol", fd_tol, "");
		optionlist_->GetStringValue("hessian-approximation", hessain_approximation_type, "");

		/*!<
		 *Set first deriver
		 >*/
		if (first_derive_type == "finite-difference")
		{
			first_derive_ = shared_ptr<LpFDderive>(new LpFDderive(funs_, fd_tol));
		}
		else if (first_derive_type == "analytic")
		{
			first_derive_ = shared_ptr<LpAnalyticDerive>(new LpAnalyticDerive(funs_));
			//!< check user input derive
			std::string user_derive_check;
			double check_tol;
			optionlist_->GetStringValue("analytic-derive-check", user_derive_check, "");
			if (user_derive_check == "yes")
			{
				optionlist_->GetNumericValue("analytic-derive-check-tol", check_tol, "");
				//!< Set analytic derive checker
				an_derive_checker_.reset(new LpANDeriveChecker(msg_reporter_, first_derive_, funs_, cd_data_, optpro_, check_tol));
			}
		}

		else
		{
			LP_THROW_EXCEPTION(LPOPC_ALGORITHM_ERROR, "Error first deriver detected!");
		}
		//!<Set second deriver
		if (hessain_approximation_type == "exact")
		{
			hessian_.reset(new LpHessianCalculator(funs_, first_derive_, cd_data_, optpro_, fd_tol));
		}
		else if (hessain_approximation_type == "limited-memory")
		{
		}
		else
		{
			LP_THROW_EXCEPTION(LPOPC_ALGORITHM_ERROR, "Error second deriver detected!");
		}
		//!<check auto scale
		std::string auto_scale_setting;
		optionlist_->GetStringValue("auto-scale", auto_scale_setting, "");
		if (auto_scale_setting == "yes")
		{
			cd_data_->autoscale = true;
			ocp_scale_.reset(new LpScaleOCP(optpro_, cd_data_, first_derive_));
		}
		else
		{
			cd_data_->autoscale = false;
		}
		//!< Set Size Checker

		size_checker_.reset(new LpSizeChecker());

		//!< Set Bounds Checker

		bounds_checker_.reset(new LpBoundsChecker());
		//!< Set LGL generator


		rpm_.reset(new RPMGenerator());
		//!< Set Guess Checker
		guess_checker_.reset(new LpGuessChecker());


		//!< Set NLP Wrapper
		nlp_wrapper_.reset(new NLPWrapper(funs_, cd_data_, optpro_, first_derive_, hessian_, ocp_scale_, rpm_));


		//dependicies_checker_->GetDependiciesForJacobi(cd_data_->allPhaseDependencies);

		//!< Set NLP Solver
		nlp_splver_.reset(new NLPSolver(optionlist_));

		//!<Set NLP to OCP converter
		nlp2op_converter_.reset(new Nlp2OpConverter());

		//!< Set mesh refiner
		meshrefiner_.reset(new MeshRefiner(msg_reporter_, optionlist_));
		meshrefiner_->Initialize(funs_, cd_data_, optpro_);
	}

	LpopcAlgorithm::~LpopcAlgorithm()
	{

	}

	LpopcAlgorithm::LpopcAlgorithm(shared_ptr<OptimalProblem>optpro, shared_ptr<LpCalculateData> cd_data, shared_ptr<LpOptionsList> optionlist, shared_ptr<LpReporter> msg_reporter) :optpro_(optpro), cd_data_(cd_data), optionlist_(optionlist), msg_reporter_(msg_reporter)
	{
		LP_DBG_START_FUN("LpopcAlgorithm::LpopcAlgorithm()");
		funs_ = optpro->GetOpimalProblemFuns();
	}

	void LpopcAlgorithm::UpdateGrid()
	{
		cd_data_->current_grid = meshrefiner_->CurrentGrid();
	}

	void LpopcAlgorithm::OutputProblemInfo()
	{
		std::stringstream outputmsg;
		outputmsg.setf(ios::scientific);
		outputmsg
			<< "<<<<<<<<<< Problem Information>>>>>>>>>>" << std::endl
			<< cd_data_->numphases_ << " phases" << std::endl
			<< cd_data_->numlinks_ << " linkages" << std::endl;
			//out put every phase
		for (size_t iphase = 0; iphase < cd_data_->numphases_; iphase++)
		{
			size_t nstates = cd_data_->SIZES_[iphase][0];
			size_t ncontrols = cd_data_->SIZES_[iphase][1];
			size_t nparameters = cd_data_->SIZES_[iphase][2];
			size_t npaths = cd_data_->SIZES_[iphase][3];
			size_t nevents = cd_data_->SIZES_[iphase][4];
			outputmsg << "  --------phase" << iphase << "--------" << std::endl;
			outputmsg << " nstates----" << nstates
				<< " ncontrols----" << ncontrols
				<< " nparameters----" << nparameters
				<< " npaths----" << npaths
				<< " nevents----" << nevents << std::endl;
		}

			outputmsg<< "<<<<<<<<<<Options>>>>>>>>> " << std::endl;
		std::string auto_scale;
		optionlist_->GetStringValue("auto-scale", auto_scale, "");
		outputmsg << " auto-scale     -->" << auto_scale <<std::endl;

		if (auto_scale=="yes")
		{
			outputmsg << "  Warning,auto-scale may fail!" << endl;
		}
		double ipopt_tol;
		optionlist_->GetNumericValue("Ipopt-tol", ipopt_tol, "");
		outputmsg << " Ipopt-tol     -->" << ipopt_tol << std::endl;

		std::string first_derive;
		optionlist_->GetStringValue("first-derive", first_derive,"");
       outputmsg << " first-derive     -->" << first_derive << std::endl;

		std::string hessian_approximation;
		optionlist_->GetStringValue("hessian-approximation", hessian_approximation, "");
		outputmsg << " hessian_approximation     -->" << hessian_approximation << std::endl;
        
		 if (first_derive == "analytic")
		{
			std::string analytic_derive_check;
			optionlist_->GetStringValue("analytic-derive-check", analytic_derive_check, "");
			outputmsg << " analytic-derive-check     -->" << analytic_derive_check << std::endl;
			if (analytic_derive_check=="yes")
			{
				std::string analytic_derive_check_tol;
				optionlist_->GetStringValue("analytic-derive-check-tol", analytic_derive_check_tol, "");
				outputmsg << " analytic-derive-check-tol     -->" << analytic_derive_check_tol << std::endl;
			}
		}
		 double finite_difference_tol;
		 optionlist_->GetNumericValue("finite-difference-tol", finite_difference_tol, "");
		 outputmsg << " finite-difference-tol     -->" << finite_difference_tol << std::endl;

		std::string mesh_refine_methods;
		optionlist_->GetStringValue("mesh-refine-methods", mesh_refine_methods, "");
		outputmsg << " mesh-refine-methods     -->" << mesh_refine_methods << std::endl;

		int max_grid_num;
		optionlist_->GetIntegerValue("max-grid-num", max_grid_num, "");
		outputmsg << " max-grid-num     -->" << max_grid_num << std::endl;
		
		if (mesh_refine_methods=="ph")
		{
			int Nmax, Nmin;
			optionlist_->GetIntegerValue("Nmin", Nmin, "");
			outputmsg << " Nmin(ph)     -->" << Nmin << std::endl;
			optionlist_->GetIntegerValue("Nmax", Nmax, "");
			outputmsg << " Nmax(ph)     -->" << Nmax << std::endl;
		}
		double desired_relative_error;
		optionlist_->GetNumericValue("desired-relative-error", desired_relative_error, "");
		outputmsg << " desired-relative-error-->" << desired_relative_error << std::endl;
		msg_reporter_->Printf()->info(outputmsg.str());
	}

}//namespaace Lpopc