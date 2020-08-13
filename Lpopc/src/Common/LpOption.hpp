// Copyright (C) 2014-2015 Xue Zhichen, Wang Yujie,Wang Na 
// All Rights Reserved.
// This file is a part of LPOPC , published under the Eclipse Public License.
// Author:Xue Zhichen 8/4 2015   21:43
// Email:eddy_lpopc@163.com
#ifndef LPOPTION_HPP
#define LPOPTION_HPP
#include"LpException.hpp"
#include"LpReporter.hpp"
#include<string>
#include<cassert>
#include<map>
#include<list>
namespace Lpopc
{
	enum OptionType
	{
		Number_type,
		Integer_type,
		String_type,
		Unknown_type
	};

	class LpOption
	{
  public:
	/** class to hold the valid string settings for a string option */
	class string_entry
	{
	public:
	  string_entry(const std::string& value, const std::string& description)
		  : value_(value), description_(description)
	  {}
	  std::string value_;
	  std::string description_;
	};

	/** Constructors / Destructors */
	//@{
	LpOption(int counter)
		:
		type_(Unknown_type),
		has_lower_(false),
		has_upper_(false),
		counter_(counter)
	{}

	LpOption(const std::string& name,
					 const std::string& short_description,
					 const std::string& long_description,
					 const std::string& registering_category,
					 int counter)
		:
		name_(name),
		short_description_(short_description),
		long_description_(long_description),
		registering_category_(registering_category),
		type_(Unknown_type),
		has_lower_(false),
		has_upper_(false),
		counter_(counter)
	{}

	LpOption(const LpOption& copy)
		:
		name_(copy.name_),
		short_description_(copy.short_description_),
		long_description_(copy.long_description_),
		registering_category_(copy.registering_category_),
		type_(copy.type_),
		has_lower_(copy.has_lower_),
		lower_(copy.lower_),
		has_upper_(copy.has_upper_),
		upper_(copy.upper_),
		valid_strings_(copy.valid_strings_),
		counter_(copy.counter_)
	{}

	 ~LpOption()
	{}
	//@}

	LP_DECLARE_EXCEPTION(ERROR_CONVERTING_STRING_TO_ENUM);

	/** Standard Get / Set Methods */
	//@{
	/** Get the option's name (tag in the input file) */
	 const std::string& Name() const
	{
	  return name_;
	}
	/** Set the option's name (tag in the input file) */
	 void SetName(const std::string& name)
	{
	  name_ = name;
	}
	/** Get the short description */
	const std::string& ShortDescription() const
	{
	  return short_description_;
	}
	/** Get the long description */
	 const std::string& LongDescription() const
	{
	  return long_description_;
	}
	/** Set the short description */
	 void SetShortDescription(const std::string& short_description)
	{
	  short_description_ = short_description;
	}
	/** Set the long description */
	 void SetLongDescription(const std::string& long_description)
	{
	  long_description_ = long_description;
	}
	/** Get the registering class */
	 const std::string& RegisteringCategory() const
	{
	  return registering_category_;
	}
	/** Set the registering class */
	 void SetRegisteringCategory(const std::string& registering_category)
	{
	  registering_category_ = registering_category;
	}
	/** Get the Option's type */
	 const OptionType& Type() const
	{
	  return type_;
	}
	/** Get the Option's type */
	 void SetType(const OptionType& type)
	{
	  type_ = type;
	}
	/** Counter */
	 int Counter() const
	{
	  return counter_;
	}
	//@}

	/** @name Get / Set methods valid for specific types - NOTE: the Type
	 *  must be set before calling these methods.
	 */
	//@{
	/** check if the option has a lower bound - can be called for
	 *  Number_type & Integer_type*/
	 const bool& HasLower() const
	{
	  assert(type_ == Number_type || type_ == Integer_type);
	  return has_lower_;
	}
	/** check if the lower bound is strict - can be called for
	Number_type */
	 const bool& LowerStrict() const
	{
	  assert(type_ == Number_type && has_lower_ == true);
	  return lower_strict_;
	}
	/** get the Number version of the lower bound - can be called for
	 *  Number_type */
	 double LowerNumber() const
	{
	  assert(has_lower_ == true && type_ == Number_type);
	  return lower_;
	}
	/** set the Number version of the lower bound - can be called for
	 *  Number_type */
	 void SetLowerNumber(const double& lower, const bool& strict)
	{
	  assert(type_ == Number_type);
	  lower_ = lower;
	  lower_strict_ = strict;
	  has_lower_ = true;
	}
	/** get the Integer version of the lower bound can be called for
	 *  Integer_type*/
	 int LowerInteger() const
	{
	  assert(has_lower_ == true && type_ == Integer_type);
	  return (int)lower_;
	}
	/** set the Integer version of the lower bound - can be called for
	 *  Integer_type */
	 void SetLowerInteger(const int& lower)
	{
	  assert(type_ == Integer_type);
	  lower_ = (double)lower;
	  has_lower_ = true;
	}
	/** check if the option has an upper bound - can be called for
	 *  Number_type & Integer_type*/
	 const bool& HasUpper() const
	{
	  assert(type_ == Number_type || type_ == Integer_type);
	  return has_upper_;
	}
	/** check if the upper bound is strict - can be called for
	 *  Number_type */
	 const bool& UpperStrict() const
	{
	  assert(type_ == Number_type && has_upper_ == true);
	  return upper_strict_;
	}
	/** get the Number version of the upper bound - can be called for
	 *  Number_type */
	 double UpperNumber() const
	{
	  assert(has_upper_ == true && type_ == Number_type);
	  return upper_;
	}
	/** set the Number version of the upper bound - can be called for
	 *  Number_type */
	 void SetUpperNumber(const double& upper, const bool& strict)
	{
	  assert(type_ == Number_type);
	  upper_ = upper;
	  upper_strict_ = strict;
	  has_upper_ = true;
	}
	/** get the Integer version of the upper bound - can be called for
	 *  Integer_type*/
	 int UpperInteger() const
	{
	  assert(has_upper_ == true && type_ == Integer_type);
	  return (int)upper_;
	}
	/** set the Integer version of the upper bound - can be called for
	 *  Integer_type */
	void SetUpperInteger(const int& upper)
	{
	  assert(type_ == Integer_type);
	  upper_ = (double)upper;
	  has_upper_ = true;
	}
	/** method to add valid string entries - can be called for
	 *  String-type */
	 void AddValidStringSetting(const std::string value,
									   const std::string description)
	{
	  assert(type_ == String_type);
	  valid_strings_.push_back(string_entry(value, description));
	}
	/** get the default as a Number - can be called for Number_type */
	 double DefaultNumber() const
	{
	  assert(type_ == Number_type);
	  return default_number_;
	}
	/** Set the default as a Number - can be called for Number_type */
	 void SetDefaultNumber(const double& default_value)
	{
	  assert(type_ == Number_type);
	  default_number_ = default_value;
	}
	/** get the default as an Integer - can be called for Integer_type*/
	 int DefaultInteger() const
	{
	  assert(type_ == Integer_type);
	  return (int)default_number_;
	}
	/** Set the default as an Integer - can be called for
	Integer_type */
	 void SetDefaultInteger(const int& default_value)
	{
	  assert(type_ == Integer_type);
	  default_number_ = (double)default_value;
	}
	/** get the default as a string - can be called for String_type */
	 std::string DefaultString() const
	{
	  assert(type_ == String_type);
	  return default_string_;
	}
	/** get the default as a string, but as the index of the string in
	 *  the list - helps map from a string to an enum- can be called
	 *  for String_type */
	 int DefaultStringAsEnum() const
	{
	  assert(type_ == String_type);
	  return MapStringSettingToEnum(default_string_);
	}
	/** Set the default as a string - can be called for String_type */
	  void SetDefaultString(const std::string& default_value)
	{
	  assert(type_ == String_type);
	  default_string_ = default_value;
	}
	/** get the valid string settings - can be called for String_type */
	  std::vector<string_entry> GetValidStrings() const
	{
	  assert(type_ == String_type);
	  return valid_strings_;
	}
	/** Check if the Number value is a valid setting - can be called
	 *  for Number_type */
	  bool IsValidNumberSetting(const double& value) const
	{
	  assert(type_ == Number_type);
	  if (has_lower_ && ((lower_strict_ == true && value <= lower_) ||
						 (lower_strict_ == false && value < lower_))) {
		return false;
	  }
	  if (has_upper_ && ((upper_strict_ == true && value >= upper_) ||
						 (upper_strict_ == false && value > upper_))) {
		return false;
	  }
	  return true;
	}
	/** Check if the Integer value is a valid setting - can be called
	 *  for Integer_type */
	 bool IsValidIntegerSetting(const int& value) const
	{
	  assert(type_ == Integer_type);
	  if (has_lower_ && value < lower_) {
		return false;
	  }
	  if (has_upper_ && value > upper_) {
		return false;
	  }
	  return true;
	}
	/** Check if the String value is a valid setting - can be called
	 *  for String_type */
	 bool IsValidStringSetting(const std::string& value) const;

	/** Map a user setting (allowing any case) to the case used when
	 *  the setting was registered.
	 */
	 std::string MapStringSetting(const std::string& value) const;

	/** Map a user setting (allowing any case) to the index of the
	 *  matched setting in the list of string settings. Helps map a
	 *  string setting to an enumeration.
	 */
	 int MapStringSettingToEnum(const std::string& value) const;
	//@}

	/** output a description of the option */
	 void OutputDescription(const LpReporter& reporter) const;
	/** output a more concise version */
	 void OutputShortDescription(const LpReporter& reporter) const;


  private:
	std::string name_;
	std::string short_description_;
	std::string long_description_;
	std::string registering_category_;
	OptionType type_;

	bool has_lower_;
	bool lower_strict_;
	double lower_;
	bool has_upper_;
	bool upper_strict_;
	double upper_;
	double default_number_;

	LpOption& operator=(const LpOption&);
	/** Compare two strings and return true if they are equal (case
	insensitive comparison) */
	bool string_equal_insensitive(const std::string& s1,
								  const std::string& s2) const;

	std::vector<string_entry> valid_strings_;
	std::string default_string_;

	/** Has the information as how many-th option this one was
	 *  registered. */
	const int counter_;
  };

  /** Class for storing registered options. Used for validation and
   *  documentation.
   */
  class LpRegisteredOptions 
  {
  public:
	/** Constructors / Destructors */
	//@{
	/** Standard Constructor */
	LpRegisteredOptions()
		:
		next_counter_(0),
		current_registering_category_("Uncategorized")
	{}

	/** Standard Destructor */
	 ~LpRegisteredOptions()
	{}
	//@}

	LP_DECLARE_EXCEPTION(OPTION_ALREADY_REGISTERED);

	/** Methods to interact with registered options */
	//@{
	/** set the registering class. All subsequent options will be
	 *  added with the registered class */
	 void SetRegisteringCategory(const std::string& registering_category)
	{
	  current_registering_category_ = registering_category;
	}

	/** retrieve the value of the current registering category */
	 std::string RegisteringCategory()
	{
	  return current_registering_category_;
	}

	/** Add a Number option (with no restrictions) */
	 void AddNumberOption(const std::string& name,
								 const std::string& short_description,
								 double default_value,
								 const std::string& long_description="");
	/** Add a Number option (with a lower bound) */
	 void AddLowerBoundedNumberOption(const std::string& name,
		const std::string& short_description,
		double lower, bool strict,
		double default_value,
		const std::string& long_description="");
	/** Add a Number option (with a upper bound) */
	 void AddUpperBoundedNumberOption(const std::string& name,
		const std::string& short_description,
		double upper, bool strict,
		double default_value,
		const std::string& long_description="");
	/** Add a Number option (with a both bounds) */
	 void AddBoundedNumberOption(const std::string& name,
										const std::string& short_description,
										double lower, bool lower_strict,
										double upper, bool upper_strict,
										double default_value,
										const std::string& long_description="");
	/** Add a Integer option (with no restrictions) */
	 void AddIntegerOption(const std::string& name,
								  const std::string& short_description,
								  int default_value,
								  const std::string& long_description="");
	/** Add a Integer option (with a lower bound) */
	 void AddLowerBoundedIntegerOption(const std::string& name,
		const std::string& short_description,
		int lower, int default_value,
		const std::string& long_description="");
	/** Add a Integer option (with a upper bound) */
	 void AddUpperBoundedIntegerOption(const std::string& name,
		const std::string& short_description,
		int upper, int default_value,
		const std::string& long_description="");
	/** Add a Integer option (with a both bounds) */
	 void AddBoundedIntegerOption(const std::string& name,
										 const std::string& short_description,
										 int lower, int upper,
										 int default_value,
										 const std::string& long_description="");

	/** Add a String option (with no restrictions) */
	 void AddStringOption(const std::string& name,
								 const std::string& short_description,
								 const std::string& default_value,
								 const std::vector<std::string>& settings,
								 const std::vector<std::string>& descriptions,
								 const std::string& long_description="");
	/** Methods that make adding string options with only a few
	 *  entries easier */
	 void AddStringOption1(const std::string& name,
								  const std::string& short_description,
								  const std::string& default_value,
								  const std::string& setting1,
								  const std::string& description1,
								  const std::string& long_description="");
	 void AddStringOption2(const std::string& name,
								  const std::string& short_description,
								  const std::string& default_value,
								  const std::string& setting1,
								  const std::string& description1,
								  const std::string& setting2,
								  const std::string& description2,
								  const std::string& long_description="");
	 void AddStringOption3(const std::string& name,
								  const std::string& short_description,
								  const std::string& default_value,
								  const std::string& setting1,
								  const std::string& description1,
								  const std::string& setting2,
								  const std::string& description2,
								  const std::string& setting3,
								  const std::string& description3,
								  const std::string& long_description="");
	 void AddStringOption4(const std::string& name,
								  const std::string& short_description,
								  const std::string& default_value,
								  const std::string& setting1,
								  const std::string& description1,
								  const std::string& setting2,
								  const std::string& description2,
								  const std::string& setting3,
								  const std::string& description3,
								  const std::string& setting4,
								  const std::string& description4,
								  const std::string& long_description="");
	 void AddStringOption5(const std::string& name,
								  const std::string& short_description,
								  const std::string& default_value,
								  const std::string& setting1,
								  const std::string& description1,
								  const std::string& setting2,
								  const std::string& description2,
								  const std::string& setting3,
								  const std::string& description3,
								  const std::string& setting4,
								  const std::string& description4,
								  const std::string& setting5,
								  const std::string& description5,
								  const std::string& long_description="");
	 void AddStringOption6(const std::string& name,
								  const std::string& short_description,
								  const std::string& default_value,
								  const std::string& setting1,
								  const std::string& description1,
								  const std::string& setting2,
								  const std::string& description2,
								  const std::string& setting3,
								  const std::string& description3,
								  const std::string& setting4,
								  const std::string& description4,
								  const std::string& setting5,
								  const std::string& description5,
								  const std::string& setting6,
								  const std::string& description6,
								  const std::string& long_description="");
	 void AddStringOption7(const std::string& name,
								  const std::string& short_description,
								  const std::string& default_value,
								  const std::string& setting1,
								  const std::string& description1,
								  const std::string& setting2,
								  const std::string& description2,
								  const std::string& setting3,
								  const std::string& description3,
								  const std::string& setting4,
								  const std::string& description4,
								  const std::string& setting5,
								  const std::string& description5,
								  const std::string& setting6,
								  const std::string& description6,
								  const std::string& setting7,
								  const std::string& description7,
								  const std::string& long_description="");
	 void AddStringOption8(const std::string& name,
								  const std::string& short_description,
								  const std::string& default_value,
								  const std::string& setting1,
								  const std::string& description1,
								  const std::string& setting2,
								  const std::string& description2,
								  const std::string& setting3,
								  const std::string& description3,
								  const std::string& setting4,
								  const std::string& description4,
								  const std::string& setting5,
								  const std::string& description5,
								  const std::string& setting6,
								  const std::string& description6,
								  const std::string& setting7,
								  const std::string& description7,
								  const std::string& setting8,
								  const std::string& description8,
								  const std::string& long_description="");
	 void AddStringOption9(const std::string& name,
								  const std::string& short_description,
								  const std::string& default_value,
								  const std::string& setting1,
								  const std::string& description1,
								  const std::string& setting2,
								  const std::string& description2,
								  const std::string& setting3,
								  const std::string& description3,
								  const std::string& setting4,
								  const std::string& description4,
								  const std::string& setting5,
								  const std::string& description5,
								  const std::string& setting6,
								  const std::string& description6,
								  const std::string& setting7,
								  const std::string& description7,
								  const std::string& setting8,
								  const std::string& description8,
								  const std::string& setting9,
								  const std::string& description9,
								  const std::string& long_description="");
	 void AddStringOption10(const std::string& name,
								   const std::string& short_description,
								   const std::string& default_value,
								   const std::string& setting1,
								   const std::string& description1,
								   const std::string& setting2,
								   const std::string& description2,
								   const std::string& setting3,
								   const std::string& description3,
								   const std::string& setting4,
								   const std::string& description4,
								   const std::string& setting5,
								   const std::string& description5,
								   const std::string& setting6,
								   const std::string& description6,
								   const std::string& setting7,
								   const std::string& description7,
								   const std::string& setting8,
								   const std::string& description8,
								   const std::string& setting9,
								   const std::string& description9,
								   const std::string& setting10,
								   const std::string& description10,
								   const std::string& long_description="");

	/** Get a registered option - this will return NULL if the option
	 *  does not exist */
	 shared_ptr<const LpOption> GetOption(const std::string& name);

	/** Output documentation for the options - gives a description,
	 *  etc. */
	 void OutputOptionDocumentation(const LpReporter& reporter, std::list<std::string>& categories);


	typedef std::map<std::string, shared_ptr<LpOption> > RegOptionsList;

	/** Giving access to iteratable representation of the registered
	 *  options */
	 const RegOptionsList& RegisteredOptionsList () const
	{
	  return registered_options_;
	}

  private:
	int next_counter_;
	std::string current_registering_category_;
	std::map<std::string, shared_ptr<LpOption> > registered_options_;
  };

}//namespace Lpopc
#endif