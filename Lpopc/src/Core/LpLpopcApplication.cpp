// Copyright (C) 2014-2015 Xue Zhichen, Wang Yujie,Wang Na 
// All Rights Reserved.
// This file is a part of LPOPC , published under the Eclipse Public License.
// Author:Xue Zhichen 8/4 2015   21:54
// Email:eddy_lpopc@163.com
#include "LpLpopcApplication.hpp"
#include  <sstream>
namespace Lpopc
{


	void LpopcApplication::RegisterAllLpopcOptions(shared_ptr<LpRegisteredOptions>&)
	{
		registeredoption_->SetRegisteringCategory("derive");
		OptDerive::SetAllOptions(registeredoption_);
		registeredoption_->SetRegisteringCategory("mesh-refine");
		MeshRefiner::SetAllOptions(registeredoption_);
		registeredoption_->SetRegisteringCategory("NLP-wrapper");
		NLPWrapper::SetAllOptions(registeredoption_);

		optionlist_->SetRegisteredOptions(registeredoption_);
	}

	LpopcApplication::~LpopcApplication()
	{

	}

	shared_ptr<LpOptionsList> LpopcApplication::Options()
	{
		return optionlist_;
	}

	void LpopcApplication::SolveOptimalProblem()
	{
		LP_DBG_START_FUN("LpopcApplication::SolveOptimalProblem()")
			try
		{
			algorithm_.reset(new LpopcAlgorithm(optpro_, data_, optionlist_, msg_reporter_));
			algorithm_->Initialized();
			algorithm_->SolveOptimalControlProblem();
		}
		catch (OPTION_INVALID& e)
		{
			e.ReportException(*msg_reporter_);
		}
		catch (LpopcException& e)
		{
			e.ReportException(*msg_reporter_);
			//	LP_THROW_EXCEPTION(LPOPC_APPLICATION_ERROR, "Caught unknown Lpopc exception");
		}
		catch (std::logic_error& e)
		{
			msg_reporter_->Printf()->error(e.what());
			//	LP_THROW_EXCEPTION(LPOPC_APPLICATION_ERROR, "Caught Arma Matrix Lib exception");
		}
		catch (std::runtime_error& e)
		{
			msg_reporter_->Printf()->error(e.what());
			//	LP_THROW_EXCEPTION(LPOPC_APPLICATION_ERROR, "Caught Arma Matrix Lib exception");
		}
		catch (std::bad_alloc& e)
		{
			msg_reporter_->Printf()->error("\nEXIT: Not enough memory.\n");
			//	LP_THROW_EXCEPTION(LPOPC_APPLICATION_ERROR,
			//		"Not enough memory");
		}
		catch (...)
		{
			LpopcException exc("Unknown Exception caught in LPOPC", "Unknown File", -1);
			exc.ReportException(*msg_reporter_);
			//	LP_THROW_EXCEPTION(LPOPC_APPLICATION_ERROR,
			//		"Caught unknown exception");
		}
	}

	void LpopcApplication::SetOptimalControlProblem(shared_ptr<OptimalProblem>& user_optimal_control_problem)
	{
		optpro_ = user_optimal_control_problem;
	}

	LpopcApplication::LpopcApplication(console_print_not_print if_console_print)
	{
		LP_DBG_START_FUN("LpopcApplication::LpopcApplication(console_print_not_print)")
			optionlist_.reset(new LpOptionsList());

		data_.reset(new LpCalculateData());
		try
		{
#ifdef LPOPC_REPORT_DBG_MSG


			std::vector<spdlog::sink_ptr> dbg_sinks;
			//dbg_sinks.push_back(std::make_shared<spdlog::sinks::stdout_sink_st>());
			dbg_sinks.push_back(std::make_shared<spdlog::sinks::simple_file_sink_st>("lpopc-dbg-msg.txt", true));
			auto dbgcombined_logger = std::make_shared<spdlog::logger>("lpopc_debug_logger", begin(dbg_sinks), end(dbg_sinks));
			dbgcombined_logger->set_level(WriteLevel::trace);
			spdlog::register_logger(dbgcombined_logger);
			DebugMsgPrinter::SetReporter();
#endif 
			std::vector<spdlog::sink_ptr> sinks;
			sinks.push_back(std::make_shared<spdlog::sinks::stdout_sink_st>());
			sinks.push_back(std::make_shared<spdlog::sinks::simple_file_sink_st>("lpopc-main-msg.txt", true));
			auto combined_logger = std::make_shared<spdlog::logger>("lpopc_main_logger", begin(sinks), end(sinks));
			combined_logger->set_level(WriteLevel::info);
			combined_logger->set_pattern("%v");
			spdlog::register_logger(combined_logger);

			msg_reporter_.reset(new LpReporter());
			registeredoption_.reset(new LpRegisteredOptions());
			RegisterAllLpopcOptions(registeredoption_);
			optionlist_->SetReporter(msg_reporter_);
			optionlist_->SetRegisteredOptions(registeredoption_);
			std::stringstream Lpopc_copyright_msg;
			Lpopc_copyright_msg
				<< "********************************************************************************" << std::endl
				<< "                LPOPC Optimal Control problem solver using RPM                                 " << std::endl
				<< "     Copyright (C) 2014-2015 Xue Zhichen, Wang Yujie, Wang Na. NJUST CHINA       " << std::endl
				<< "                LPOPC Version " << LPOPC_VERSION << ",Published with EPL license            " << std::endl
				<< "                           Email:eddy_lpopc@163.com                                                                     " << std::endl
				<< "********************************************************************************" << std::endl;
			msg_reporter_->Printf()->info(Lpopc_copyright_msg.str());
		}
		catch (OPTION_INVALID& e)
		{
			e.ReportException(*msg_reporter_);
		}
		catch (LpopcException& e)
		{
			e.ReportException(*msg_reporter_);
			//	LP_THROW_EXCEPTION(LPOPC_APPLICATION_ERROR, "Caught unknown Lpopc exception");
		}
		catch (std::logic_error& e)
		{
			msg_reporter_->Printf()->info(e.what());
			//	LP_THROW_EXCEPTION(LPOPC_APPLICATION_ERROR, "Caught Arma Matrix Lib exception");
		}
		catch (std::runtime_error& e)
		{
			msg_reporter_->Printf()->info(e.what());
			//	LP_THROW_EXCEPTION(LPOPC_APPLICATION_ERROR, "Caught Arma Matrix Lib exception");
		}
		catch (std::bad_alloc& e)
		{
			msg_reporter_->Printf()->error("\nEXIT: Not enough memory.\n");
			//	LP_THROW_EXCEPTION(LPOPC_APPLICATION_ERROR,
			//		"Not enough memory");
		}
		catch (...)
		{
			LpopcException exc("Unknown Exception caught in LPOPC", "Unknown File", -1);
			exc.ReportException(*msg_reporter_);
			//	LP_THROW_EXCEPTION(LPOPC_APPLICATION_ERROR,
			//		"Caught unknown exception");
		}
	}

}