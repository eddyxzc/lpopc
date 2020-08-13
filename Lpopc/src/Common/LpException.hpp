// Copyright (C) 2014-2015 Xue Zhichen, Wang Yujie,Wang Na 
// All Rights Reserved.
// This file is a part of LPOPC , published under the Eclipse Public License.
// Author:Xue Zhichen 7/10 2015   13:27
// Email:eddy_lpopc@163.com
#ifndef  LPEXCEPTION_HPP
#define LPEXCEPTION_HPP

#include"LpReporter.hpp"
#include<string>
#include <exception>
namespace Lpopc
{
	class LpopcException
	{
	public:
		/**@name Constructors/Destructors */
		//@{
		/** Constructor */
		LpopcException(std::string msg, std::string file_name, int line_number, std::string type="LpopcException")
			:
			msg_(msg),
			file_name_(file_name),
			line_number_(line_number),
			type_(type)
		{}

		/** Copy Constructor */
		LpopcException(const LpopcException& copy)
			:
			msg_(copy.msg_),
			file_name_(copy.file_name_),
			line_number_(copy.line_number_),
			type_(copy.type_)
		{}

		/** Default destructor */
		virtual ~LpopcException()
		{}
		//@}

		/** Method to report the exception to a journalist */
		void ReportException(const LpReporter& lreporter) const
		{
			lreporter.Printf()->error(
				"Exception of type: {} in file \"{}\" at line {}:\n Exception message: {}\n",
				type_.c_str(), file_name_.c_str(),  line_number_, msg_.c_str());
		}

		const std::string& Message() const
		{
			return msg_;
		}

	private:
		/**@name Default Compiler Generated Methods
		* (Hidden to avoid implicit creation/calling).
		* These methods are not implemented and 
		* we do not want the compiler to implement
		* them for us, so we declare them private
		* and do not define them. This ensures that
		* they will not be implicitly created/called. */
		//@{
		/** Default Constructor */
		LpopcException();

		/** Overloaded Equals Operator */
		void operator=(const LpopcException&);
		//@}

		std::string msg_;
		std::string file_name_;
		int line_number_;
		std::string type_;
	};
}//namespaceLPOPC

#define LP_THROW_EXCEPTION(__except_type, __msg) \
	throw __except_type( (__msg), (__FILE__), (__LINE__) );

#define LP_DECLARE_EXCEPTION(__except_type) \
class __except_type : public Lpopc::LpopcException \
{ \
public: \
	__except_type(std::string msg, std::string fname, int line) \
	: Lpopc::LpopcException(msg,fname,line, #__except_type) {} \
	__except_type(const __except_type& copy) \
	: Lpopc::LpopcException(copy) {} \
private: \
	__except_type(); \
	void operator=(const __except_type&); \
}

#define LP_ASSERT_EXCEPTION(__condition, __except_type, __msg) \
	if (! (__condition) ) { \
	std::string newmsg = #__condition; \
	newmsg += " evaluated false: "; \
	newmsg += __msg; \
	throw __except_type( (newmsg), (__FILE__), (__LINE__) ); \
	}
#endif