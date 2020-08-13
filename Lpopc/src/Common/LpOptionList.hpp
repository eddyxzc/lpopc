// Copyright (C) 2014-2015 Xue Zhichen, Wang Yujie,Wang Na 
// All Rights Reserved.
// This file is a part of LPOPC , published under the Eclipse Public License.
// Author:Xue Zhichen 8/4 2015   21:43
// Email:eddy_lpopc@163.com
#ifndef  LPOPTIONLIST_HPP
#define LPOPTIONLIST_HPP

#include"LpConf.h"
#include"LpException.hpp"
#include"LpOption.hpp"
#include<iostream>
#include<map>
namespace Lpopc
{
	LP_DECLARE_EXCEPTION(OPTION_INVALID);
class LpOptionsList 
  {
	/** Class for storing the value and counter for each option in
	 *  OptionsList. */
	class OptionValue
	{
	public:
	  /**@name Constructors/Destructors */
	  //@{
	  /** Default constructor (needed for the map) */
	  OptionValue()
		  :
		  initialized_(false)
	  {}

	  /** Constructor given the value */
	  OptionValue(std::string value, bool allow_clobber, bool dont_print)
		  :
		  value_(value),
		  counter_(0),
		  initialized_(true),
		  allow_clobber_(allow_clobber),
		  dont_print_(dont_print)
	  {}

	  /** Copy Constructor */
	  OptionValue(const OptionValue& copy)
		  :
		  value_(copy.value_),
		  counter_(copy.counter_),
		  initialized_(copy.initialized_),
		  allow_clobber_(copy.allow_clobber_),
		  dont_print_(copy.dont_print_)
	  {}

	  /** Equals operator */
	  void operator=(const OptionValue& copy)
	  {
		value_=copy.value_;
		counter_=copy.counter_;
		initialized_=copy.initialized_;
		allow_clobber_=copy.allow_clobber_;
		dont_print_=copy.dont_print_;
	  }

	  /** Default Destructor */
	  ~OptionValue()
	  {}
	  //@}

	  /** Method for retrieving the value of an option.  Calling this
	   *  method will increase the counter by one. */
	  std::string GetValue() const
	  {
		assert(initialized_);
		counter_++;
		return value_;
	  }

	  /** Method for retrieving the value without increasing the
	   *  counter */
	  std::string Value() const
	  {
		assert(initialized_);
		return value_;
	  }

	  /** Method for accessing current value of the request counter */
	  int Counter() const
	  {
		assert(initialized_);
		return counter_;
	  }

	  /** True if the option can be overwritten */
	  bool AllowClobber() const
	  {
		assert(initialized_);
		return allow_clobber_;
	  }

	  /** True if this option is not to show up in the
	   *  print_user_options output */
	  bool DontPrint() const
	  {
		assert(initialized_);
		return dont_print_;
	  }

	private:
	  /** Value for this option */
	  std::string value_;

	  /** Counter for requests */
	  mutable int counter_;

	  /** for debugging */
	  bool initialized_;

	  /** True if the option can be overwritten */
	  bool allow_clobber_;

	  /** True if this option is not to show up in the
	   *  print_user_options output */
	  bool dont_print_;
	};

  public:
	/**@name Constructors/Destructors */
	//@{
	LpOptionsList(shared_ptr<Lpopc::LpRegisteredOptions> reg_options, shared_ptr<LpReporter> reporter)
		: reg_options_(reg_options), reporter_(reporter)
	{}
	LpOptionsList(){};



	/** Copy Constructor */
	LpOptionsList(const LpOptionsList& copy)
	{
	  // copy all the option strings and values
	  options_ = copy.options_;
	  // copy the registered options pointer
	  reg_options_ = copy.reg_options_;
	}

	/** Default destructor */
	virtual ~LpOptionsList()
	{}

	/** Overloaded Equals Operator */
	virtual void operator=(const LpOptionsList& source)
	{
	  options_ = source.options_;
	  reg_options_ = source.reg_options_;
	  reporter_ = source.reporter_;
	}
	//@}

	/** Method for clearing all previously set options */
	virtual void clear()
	{
	  options_.clear();
	}

	/** @name Get / Set Methods */
	//@{
	virtual void SetRegisteredOptions(const shared_ptr<LpRegisteredOptions> reg_options)
	{
	  reg_options_ = reg_options;
	}
	virtual void SetReporter(const shared_ptr<LpReporter> reporter)
	{
	  reporter_ = reporter;
	}
	//@}
	/** @name Methods for setting options */
	//@{
	virtual bool SetStringValue(const std::string& tag, const std::string& value,
								bool allow_clobber = true, bool dont_print = false);
	virtual bool SetNumericValue(const std::string& tag, double value,
								 bool allow_clobber = true, bool dont_print = false);
	virtual bool SetIntegerValue(const std::string& tag, int value,
								 bool allow_clobber = true, bool dont_print = false);
	//@}

	/** @name Methods for setting options only if they have not been
	 *  set before*/
	//@{
	virtual bool SetStringValueIfUnset(const std::string& tag, const std::string& value,
									   bool allow_clobber = true, bool dont_print = false);
	virtual bool SetNumericValueIfUnset(const std::string& tag, double value,
										bool allow_clobber = true, bool dont_print = false);
	virtual bool SetIntegerValueIfUnset(const std::string& tag, int value,
										bool allow_clobber = true, bool dont_print = false);
	//@}

	/** @name Methods for retrieving values from the options list.  If
	 *  a tag is not found, the methods return false, and value is set
	 *  to the default value defined in the registered options. */
	//@{
	virtual bool GetStringValue(const std::string& tag, std::string& value,
								const std::string& prefix) const;
	virtual bool GetEnumValue(const std::string& tag, int& value,
							  const std::string& prefix) const;
	virtual bool GetBoolValue(const std::string& tag, bool& value,
							  const std::string& prefix) const;
	virtual bool GetNumericValue(const std::string& tag, double& value,
								 const std::string& prefix) const;
	virtual bool GetIntegerValue(const std::string& tag, int& value,
								 const std::string& prefix) const;
	//@}

	/** Get a string with the list of all options (tag, value, counter) */
	virtual void PrintList(std::string& list) const;

	/** Get a string with the list of all options set by the user
	 *  (tag, value, use/notused).  Here, options with dont_print flag
	 *  set to true are not printed. */
	virtual void PrintUserOptions(std::string& list) const;

	/** Read options from the stream is.  Returns false if
	 *  an error was encountered. */
	virtual bool ReadFromStream(const LpReporter& reporter, std::istream& is);

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
	//    OptionsList();

	//@}

	/** map for storing the options */
	std::map< std::string, OptionValue > options_;

	/** list of all the registered options to validate against */
	shared_ptr<LpRegisteredOptions> reg_options_;

	/** Journalist for writing error messages, etc. */
	shared_ptr<LpReporter> reporter_;

	/** auxilliary method for converting sting to all lower-case
	 *  letters */
	const std::string& lowercase(const std::string tag) const;

	/** auxilliary method for finding the value for a tag in the
	 *  options list.  This method first looks for the concatenated
	 *  string prefix+tag (if prefix is not ""), and if this is not
	 *  found, it looks for tag.  The return value is true iff
	 *  prefix+tag or tag is found.  In that case, the corresponding
	 *  string value is copied into value. */
	bool find_tag(const std::string& tag, const std::string& prefix,
				  std::string& value) const;

	/** tells whether or not we can clobber a particular option.
	 *  returns true if the option does not already exist, or if
	 *  the option exists but is set to allow_clobber
	 */
	bool will_allow_clobber(const std::string& tag) const;

	/** read the next token from stream is.  Returns false, if EOF was
	 *  reached before a tokens was ecountered. */
	bool readnexttoken(std::istream& is, std::string& token);

	/** auxilliary string set by lowercase method */
	mutable std::string lowercase_buffer_;
  };
}//namespace Lpopc
#endif