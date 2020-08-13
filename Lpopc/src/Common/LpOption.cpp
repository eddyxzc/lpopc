// Copyright (C) 2014-2015 Xue Zhichen, Wang Yujie,Wang Na 
// All Rights Reserved.
// This file is a part of LPOPC , published under the Eclipse Public License.
// Author:Xue Zhichen 8/4 2015   21:43
// Email:eddy_lpopc@163.com
#include "LpOption.hpp"
#include <set>
#include <cstdio>
#include<cctype>
#include "LpException.hpp"

namespace Lpopc
{
	void LpOption::OutputDescription(const LpReporter& reporter) const
	{
		std::string type_str = "Unknown";
		if (type_ ==Number_type) {
			type_str = "Real Number";
		}
		else if (type_ ==Integer_type) {
			type_str = "Integer";
		}
		else if (type_ ==String_type) {
			type_str = "String";
		}

		reporter.Printf()->info(
			"\n### {} ({}) ###\nCategory: {}\nDescription: {}\n",
			name_.c_str(), type_str.c_str(),
			registering_category_.c_str(), short_description_.c_str());

		if (type_ ==Number_type) {
			if (has_lower_) {
				reporter.Printf()->info("({g})",lower_)  ;
			}
			else {
				reporter.Printf()->info()<< "-inf";
			}

			if (lower_strict_) {
				reporter.Printf()->info(" < ");
			}
			else {
				reporter.Printf()->info(" <= ");
			}

			reporter.Printf()->info("({g})", default_number_);

			if (has_upper_ && upper_strict_) {
				reporter.Printf()->info( " < ");
			}
			else {
				reporter.Printf()->info( " <= ");
			}

			if (has_upper_) {
				reporter.Printf()->info( "{g}\n", upper_);
			}
			else {
				reporter.Printf()->info( "+inf\n");
			}
		}
		else if (type_ ==Integer_type) {
			if (has_lower_) {
				reporter.Printf()->info( "{}", (int)lower_);
			}
			else {
				reporter.Printf()->info( "-inf");
			}

			reporter.Printf()->info( " <= ({d}) <= ", (int)default_number_);

			if (has_upper_) {
				reporter.Printf()->info( "{d}\n", (int)upper_);
			}
			else {
				reporter.Printf()->info( "+inf\n");
			}
		}
		else if (type_ ==String_type) {
			std::vector<string_entry>::const_iterator i;
			reporter.Printf()->info( "Valid Settings:\n");
			for (i = valid_strings_.begin(); i != valid_strings_.end(); i++) {
				reporter.Printf()->info("\t{} ({})\n",
					(*i).value_.c_str(), (*i).description_.c_str());
			}
			reporter.Printf()->info( "Default: \"{}\"\n",
				default_string_.c_str());
		}
	}

	
	

	void LpOption::OutputShortDescription(const LpReporter& reporter) const
	{
		reporter.Printf()->info( "{}",  name_.c_str());


		if (type_ == Number_type) {
			if (has_lower_) {
				reporter.Printf()->info("{g}", lower_);
			}
			else {
				reporter.Printf()->info( "{}", "-inf");
			}

			if (has_lower_ && !lower_strict_) {
				reporter.Printf()->info("{} "," <= ");
			}
			else {
				reporter.Printf()->info("{}", " <  ");
			}

			reporter.Printf()->info( "{g}", default_number_);

			if (has_upper_ && !upper_strict_) {
				reporter.Printf()->info({}, " <= ");
			}
			else {
				reporter.Printf()->info({}, " <  ");
			}

			if (has_upper_) {
				reporter.Printf()->info( "{g}\n", upper_);
			}
			else {
				reporter.Printf()->info( "{}\n", "+inf");
			}
		}
		else if (type_ == Integer_type) {
			if (has_lower_) {
				reporter.Printf()->info( "{} <= ", (int)lower_);
			}
			else {
				reporter.Printf()->info( "{} <  ", "-inf");
			}

			reporter.Printf()->info( "({})",
				(int)default_number_);

			if (has_upper_) {
				reporter.Printf()->info( " <= {}\n", (int)upper_);
			}
			else {
				reporter.Printf()->info( " <  {}\n", "+inf");
			}
		}
		else if (type_ == String_type) {
			reporter.Printf()->info( "(\"{}\")\n",
				default_string_.c_str());
		}
		reporter.Printf()->info( "   ");
		reporter.Printf()->info(
			short_description_.c_str());
		if (long_description_ != "") {
			reporter.Printf()->info( "\n     ");
			reporter.Printf()->info(
				long_description_.c_str());
		}
		if (type_ == String_type) {
			reporter.Printf()->info( "\n   Possible values:\n");
			for (std::vector<string_entry>::const_iterator
				i = valid_strings_.begin();
				i != valid_strings_.end(); i++) {
				reporter.Printf()->info( "    {s}",
						(*i).value_.c_str());

					if( (*i).description_.length() > 0 ) {
						reporter.Printf()->info(" [");
						reporter.Printf()->info(
							(*i).description_.c_str());
						reporter.Printf()->info("]");
					}
					reporter.Printf()->info("\n");
			}
		}
		else {
			reporter.Printf()->info( "\n");
		}
		reporter.Printf()->info("\n");
	}

	bool LpOption::IsValidStringSetting(const std::string& value) const
	{
		assert(type_ == String_type);

		std::vector<string_entry>::const_iterator i;
		for (i = valid_strings_.begin(); i != valid_strings_.end(); i++) {
			if (i->value_ == "*" || string_equal_insensitive(i->value_, value)) {
				return true;
			}
		}
		return false;
	}

	std::string
		LpOption::MapStringSetting(const std::string& value) const
	{
		assert(type_ == String_type);

		std::string matched_setting = "";

		std::vector<string_entry>::const_iterator i;
		for (i = valid_strings_.begin(); i != valid_strings_.end(); i++) {
			if (i->value_ == "*") {
				matched_setting = value;
			}
			else if (string_equal_insensitive(i->value_, value)) {
				matched_setting = i->value_;
			}
		}
		return matched_setting;
	}

	int
		LpOption::MapStringSettingToEnum(const std::string& value) const
	{
		assert(type_ == String_type);

		int matched_setting = -1;

		int cnt = 0;
		std::vector<string_entry>::const_iterator i;
		for (i = valid_strings_.begin(); i != valid_strings_.end(); i++) {
			LP_ASSERT_EXCEPTION(i->value_ != "*", LpopcException,
				"Cannot map a wildcard setting to an enumeration");
			if (string_equal_insensitive(i->value_, value)) {
				matched_setting = cnt;
				break;
			}
			cnt++;
		}

		LP_ASSERT_EXCEPTION(matched_setting != -1, ERROR_CONVERTING_STRING_TO_ENUM,
			std::string("Could not find a match for setting ") + value +
			" in option: " + name_);
		return matched_setting;
	}

	bool
		LpOption::string_equal_insensitive(const std::string& s1,
		const std::string& s2) const
	{
		using namespace std;

		if (s1.size()!=s2.size())
			return false;

		string::const_iterator i1 = s1.begin();
		string::const_iterator i2 = s2.begin();

		while (i1!=s1.end()) {
			if (toupper(*i1)!=toupper(*i2))
				return false;
			i1++;
			i2++;
		}
		return true;
	}

	void
		LpRegisteredOptions::AddNumberOption(const std::string& name,
		const std::string& short_description,
		double default_value,
		const std::string& long_description)
	{
		shared_ptr<LpOption> option (
			new LpOption(name, short_description, long_description,
			current_registering_category_, next_counter_++));
		option->SetType(Number_type);
		option->SetDefaultNumber(default_value);
		LP_ASSERT_EXCEPTION(registered_options_.find(name) == registered_options_.end(),
			OPTION_ALREADY_REGISTERED,
			std::string("The option: ") + option->Name() + " has already been registered by someone else");
		registered_options_[name] = option;
	}

	void
		LpRegisteredOptions::AddLowerBoundedNumberOption(const std::string& name,
		const std::string& short_description,
		double lower, bool strict,
		double default_value,
		const std::string& long_description)
	{
		shared_ptr<LpOption> option (
			new LpOption(name, short_description, long_description,
			current_registering_category_, next_counter_++));
		option->SetType(Number_type);
		option->SetDefaultNumber(default_value);
		option->SetLowerNumber(lower, strict);
		LP_ASSERT_EXCEPTION(registered_options_.find(name) == registered_options_.end(), OPTION_ALREADY_REGISTERED,
			std::string("The option: ") + option->Name() + " has already been registered by someone else");
		registered_options_[name] = option;
	}

	void
		LpRegisteredOptions::AddUpperBoundedNumberOption(const std::string& name,
		const std::string& short_description,
		double upper, bool strict,
		double default_value,
		const std::string& long_description)
	{
		shared_ptr<LpOption> option (
			new LpOption(name, short_description, long_description,
			current_registering_category_, next_counter_++));
		option->SetType(Number_type);
		option->SetDefaultNumber(default_value);
		option->SetUpperNumber(upper, strict);
		LP_ASSERT_EXCEPTION(registered_options_.find(name) == registered_options_.end(), OPTION_ALREADY_REGISTERED,
			std::string("The option: ") + option->Name() + " has already been registered by someone else");
		registered_options_[name] = option;
	}

	void
		LpRegisteredOptions::AddBoundedNumberOption(const std::string& name,
		const std::string& short_description,
		double lower, bool lower_strict,
		double upper, bool upper_strict,
		double default_value,
		const std::string& long_description)
	{
		shared_ptr<LpOption> option (
			new LpOption(name, short_description, long_description,
			current_registering_category_, next_counter_++));
		option->SetType(Number_type);
		option->SetDefaultNumber(default_value);
		option->SetLowerNumber(lower, lower_strict);
		option->SetUpperNumber(upper, upper_strict);
		LP_ASSERT_EXCEPTION(registered_options_.find(name) == registered_options_.end(), OPTION_ALREADY_REGISTERED,
			std::string("The option: ") + option->Name() + " has already been registered by someone else");
		registered_options_[name] = option;
	}

	void
		LpRegisteredOptions::AddIntegerOption(const std::string& name,
		const std::string& short_description,
		int default_value,
		const std::string& long_description)
	{
		shared_ptr<LpOption> option (
			new LpOption(name, short_description, long_description,
			current_registering_category_, next_counter_++));
		option->SetType(Integer_type);
		option->SetDefaultInteger(default_value);
		LP_ASSERT_EXCEPTION(registered_options_.find(name) == registered_options_.end(), OPTION_ALREADY_REGISTERED,
			std::string("The option: ") + option->Name() + " has already been registered by someone else");
		registered_options_[name] = option;
	}

	void
		LpRegisteredOptions::AddLowerBoundedIntegerOption(const std::string& name,
		const std::string& short_description,
		int lower, int default_value,
		const std::string& long_description)
	{
		shared_ptr<LpOption> option (
			new LpOption(name, short_description, long_description,
			current_registering_category_, next_counter_++));
		option->SetType(Integer_type);
		option->SetDefaultInteger(default_value);
		option->SetLowerInteger(lower);
		LP_ASSERT_EXCEPTION(registered_options_.find(name) == registered_options_.end(), OPTION_ALREADY_REGISTERED,
			std::string("The option: ") + option->Name() + " has already been registered by someone else");
		registered_options_[name] = option;
	}

	void
		LpRegisteredOptions::AddUpperBoundedIntegerOption(const std::string& name,
		const std::string& short_description,
		int upper, int default_value,
		const std::string& long_description)
	{
		shared_ptr<LpOption> option (
			new LpOption(name, short_description, long_description,
			current_registering_category_, next_counter_++));
		option->SetType(Integer_type);
		option->SetDefaultInteger(default_value);
		option->SetUpperInteger(upper);
		LP_ASSERT_EXCEPTION(registered_options_.find(name) == registered_options_.end(), OPTION_ALREADY_REGISTERED,
			std::string("The option: ") + option->Name() + " has already been registered by someone else");
		registered_options_[name] = option;
	}

	void
		LpRegisteredOptions::AddBoundedIntegerOption(const std::string& name,
		const std::string& short_description,
		int lower, int upper,
		int default_value,
		const std::string& long_description)
	{
		shared_ptr<LpOption> option (
			new LpOption(name, short_description, long_description,
			current_registering_category_, next_counter_++));
		option->SetType(Integer_type);
		option->SetDefaultInteger(default_value);
		option->SetLowerInteger(lower);
		option->SetUpperInteger(upper);
		LP_ASSERT_EXCEPTION(registered_options_.find(name) == registered_options_.end(), OPTION_ALREADY_REGISTERED,
			std::string("The option: ") + option->Name() + " has already been registered by someone else");
		registered_options_[name] = option;
	}

	void
		LpRegisteredOptions::AddStringOption(const std::string& name,
		const std::string& short_description,
		const std::string& default_value,
		const std::vector<std::string>& settings,
		const std::vector<std::string>& descriptions,
		const std::string& long_description)
	{
		shared_ptr<LpOption> option (
			new LpOption(name, short_description, long_description,
			current_registering_category_, next_counter_++));
		option->SetType(String_type);
		option->SetDefaultString(default_value);
		assert(settings.size() == descriptions.size());
		for (int i=0; i<(int)settings.size(); i++) {
			option->AddValidStringSetting(settings[i], descriptions[i]);
		}
		LP_ASSERT_EXCEPTION(registered_options_.find(name) == registered_options_.end(), OPTION_ALREADY_REGISTERED,
			std::string("The option: ") + option->Name() + " has already been registered by someone else");
		registered_options_[name] = option;
	}

	void
		LpRegisteredOptions::AddStringOption1(const std::string& name,
		const std::string& short_description,
		const std::string& default_value,
		const std::string& setting1,
		const std::string& description1,
		const std::string& long_description)
	{
		shared_ptr<LpOption> option (
			new LpOption(name, short_description, long_description,
			current_registering_category_, next_counter_++));
		option->SetType(String_type);
		option->SetDefaultString(default_value);
		option->AddValidStringSetting(setting1, description1);
		LP_ASSERT_EXCEPTION(registered_options_.find(name) == registered_options_.end(), OPTION_ALREADY_REGISTERED,
			std::string("The option: ") + option->Name() + " has already been registered by someone else");
		registered_options_[name] = option;
	}

	void
		LpRegisteredOptions::AddStringOption2(const std::string& name,
		const std::string& short_description,
		const std::string& default_value,
		const std::string& setting1,
		const std::string& description1,
		const std::string& setting2,
		const std::string& description2,
		const std::string& long_description)
	{
		shared_ptr<LpOption> option (
			new LpOption(name, short_description, long_description,
			current_registering_category_, next_counter_++));
		option->SetType(String_type);
		option->SetDefaultString(default_value);
		option->AddValidStringSetting(setting1, description1);
		option->AddValidStringSetting(setting2, description2);
		LP_ASSERT_EXCEPTION(registered_options_.find(name) == registered_options_.end(), OPTION_ALREADY_REGISTERED,
			std::string("The option: ") + option->Name() + " has already been registered by someone else");
		registered_options_[name] = option;
	}

	void
		LpRegisteredOptions::AddStringOption3(const std::string& name,
		const std::string& short_description,
		const std::string& default_value,
		const std::string& setting1,
		const std::string& description1,
		const std::string& setting2,
		const std::string& description2,
		const std::string& setting3,
		const std::string& description3,
		const std::string& long_description)
	{
		shared_ptr<LpOption> option (
			new LpOption(name, short_description, long_description,
			current_registering_category_, next_counter_++));
		option->SetType(String_type);
		option->SetDefaultString(default_value);
		option->AddValidStringSetting(setting1, description1);
		option->AddValidStringSetting(setting2, description2);
		option->AddValidStringSetting(setting3, description3);
		LP_ASSERT_EXCEPTION(registered_options_.find(name) == registered_options_.end(), OPTION_ALREADY_REGISTERED,
			std::string("The option: ") + option->Name() + " has already been registered by someone else");
		registered_options_[name] = option;
	}

	void
		LpRegisteredOptions::AddStringOption4(const std::string& name,
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
		const std::string& long_description)
	{
		shared_ptr<LpOption> option (
			new LpOption(name, short_description, long_description,
			current_registering_category_, next_counter_++));
		option->SetType(String_type);
		option->SetDefaultString(default_value);
		option->AddValidStringSetting(setting1, description1);
		option->AddValidStringSetting(setting2, description2);
		option->AddValidStringSetting(setting3, description3);
		option->AddValidStringSetting(setting4, description4);
		LP_ASSERT_EXCEPTION(registered_options_.find(name) == registered_options_.end(), OPTION_ALREADY_REGISTERED,
			std::string("The option: ") + option->Name() + " has already been registered by someone else");
		registered_options_[name] = option;
	}

	void
		LpRegisteredOptions::AddStringOption5(const std::string& name,
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
		const std::string& long_description)
	{
		shared_ptr<LpOption> option (
			new LpOption(name, short_description, long_description,
			current_registering_category_, next_counter_++));
		option->SetType(String_type);
		option->SetDefaultString(default_value);
		option->AddValidStringSetting(setting1, description1);
		option->AddValidStringSetting(setting2, description2);
		option->AddValidStringSetting(setting3, description3);
		option->AddValidStringSetting(setting4, description4);
		option->AddValidStringSetting(setting5, description5);
		LP_ASSERT_EXCEPTION(registered_options_.find(name) == registered_options_.end(), OPTION_ALREADY_REGISTERED,
			std::string("The option: ") + option->Name() + " has already been registered by someone else");
		registered_options_[name] = option;
	}

	void
		LpRegisteredOptions::AddStringOption6(const std::string& name,
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
		const std::string& long_description)
	{
		shared_ptr<LpOption> option (
			new LpOption(name, short_description, long_description,
			current_registering_category_, next_counter_++));
		option->SetType(String_type);
		option->SetDefaultString(default_value);
		option->AddValidStringSetting(setting1, description1);
		option->AddValidStringSetting(setting2, description2);
		option->AddValidStringSetting(setting3, description3);
		option->AddValidStringSetting(setting4, description4);
		option->AddValidStringSetting(setting5, description5);
		option->AddValidStringSetting(setting6, description6);
		LP_ASSERT_EXCEPTION(registered_options_.find(name) == registered_options_.end(), OPTION_ALREADY_REGISTERED,
			std::string("The option: ") + option->Name() + " has already been registered by someone else");
		registered_options_[name] = option;
	}

	void
		LpRegisteredOptions::AddStringOption7(const std::string& name,
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
		const std::string& long_description)
	{
		shared_ptr<LpOption> option (
			new LpOption(name, short_description, long_description,
			current_registering_category_, next_counter_++));
		option->SetType(String_type);
		option->SetDefaultString(default_value);
		option->AddValidStringSetting(setting1, description1);
		option->AddValidStringSetting(setting2, description2);
		option->AddValidStringSetting(setting3, description3);
		option->AddValidStringSetting(setting4, description4);
		option->AddValidStringSetting(setting5, description5);
		option->AddValidStringSetting(setting6, description6);
		option->AddValidStringSetting(setting7, description7);
		LP_ASSERT_EXCEPTION(registered_options_.find(name) == registered_options_.end(), OPTION_ALREADY_REGISTERED,
			std::string("The option: ") + option->Name() + " has already been registered by someone else");
		registered_options_[name] = option;
	}

	void
		LpRegisteredOptions::AddStringOption8(const std::string& name,
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
		const std::string& long_description)
	{
		shared_ptr<LpOption> option (
			new LpOption(name, short_description, long_description,
			current_registering_category_, next_counter_++));
		option->SetType(String_type);
		option->SetDefaultString(default_value);
		option->AddValidStringSetting(setting1, description1);
		option->AddValidStringSetting(setting2, description2);
		option->AddValidStringSetting(setting3, description3);
		option->AddValidStringSetting(setting4, description4);
		option->AddValidStringSetting(setting5, description5);
		option->AddValidStringSetting(setting6, description6);
		option->AddValidStringSetting(setting7, description7);
		option->AddValidStringSetting(setting8, description8);
		LP_ASSERT_EXCEPTION(registered_options_.find(name) == registered_options_.end(), OPTION_ALREADY_REGISTERED,
			std::string("The option: ") + option->Name() + " has already been registered by someone else");
		registered_options_[name] = option;
	}

	void
		LpRegisteredOptions::AddStringOption9(const std::string& name,
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
		const std::string& long_description)
	{
		shared_ptr<LpOption> option(
			new LpOption(name, short_description, long_description,
			current_registering_category_, next_counter_++));
		option->SetType(String_type);
		option->SetDefaultString(default_value);
		option->AddValidStringSetting(setting1, description1);
		option->AddValidStringSetting(setting2, description2);
		option->AddValidStringSetting(setting3, description3);
		option->AddValidStringSetting(setting4, description4);
		option->AddValidStringSetting(setting5, description5);
		option->AddValidStringSetting(setting6, description6);
		option->AddValidStringSetting(setting7, description7);
		option->AddValidStringSetting(setting8, description8);
		option->AddValidStringSetting(setting9, description9);
		LP_ASSERT_EXCEPTION(registered_options_.find(name) == registered_options_.end(), OPTION_ALREADY_REGISTERED,
			std::string("The option: ") + option->Name() + " has already been registered by someone else");
		registered_options_[name] = option;
	}

	void
		LpRegisteredOptions::AddStringOption10(const std::string& name,
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
		const std::string& long_description)
	{
		shared_ptr<LpOption> option (
			new LpOption(name, short_description, long_description,
			current_registering_category_, next_counter_++)  );
		option->SetType(String_type);
		option->SetDefaultString(default_value);
		option->AddValidStringSetting(setting1, description1);
		option->AddValidStringSetting(setting2, description2);
		option->AddValidStringSetting(setting3, description3);
		option->AddValidStringSetting(setting4, description4);
		option->AddValidStringSetting(setting5, description5);
		option->AddValidStringSetting(setting6, description6);
		option->AddValidStringSetting(setting7, description7);
		option->AddValidStringSetting(setting8, description8);
		option->AddValidStringSetting(setting9, description9);
		option->AddValidStringSetting(setting10, description10);
		LP_ASSERT_EXCEPTION(registered_options_.find(name) == registered_options_.end(), OPTION_ALREADY_REGISTERED,
			std::string("The option: ") + option->Name() + " has already been registered by someone else");
		registered_options_[name] = option;
	}

	shared_ptr<const LpOption> LpRegisteredOptions::GetOption(const std::string& name)
	{
		std::string tag_only = name;
		std::string::size_type pos = name.rfind(".", name.length());
		if (pos != std::string::npos) {
			tag_only = name.substr(pos+1, name.length()-pos);
		}
		shared_ptr<const LpOption> option;
		std::map< std::string, shared_ptr<LpOption> >::iterator reg_option = registered_options_.find(tag_only);

		if (reg_option == registered_options_.end()) {
		}
		else {
			option =(reg_option->second);
		}

		return option;
	}

	void LpRegisteredOptions::OutputOptionDocumentation(const LpReporter& reporter, std::list<std::string>& categories)
	{
		// create a set to print sorted output
		//     std::set
		//       <std::string> classes;
		//     std::map <std::string, shared_ptr<LpOption> >::iterator option;
		//     for (option = registered_options_.begin(); option != registered_options_.end(); option++) {
		//       classes.insert(option->second->RegisteringCategory());
		//     }

		std::list<std::string>::iterator i;
		for( i = categories.begin(); i != categories.end(); i++) {
			reporter.Printf()->info("\n### {} ###\n\n", (*i).c_str());
			std::map<int, shared_ptr<LpOption> > class_options;
			std::map <std::string, shared_ptr<LpOption> >::iterator option;
			for (option = registered_options_.begin();
				option != registered_options_.end(); option++) {
					if (option->second->RegisteringCategory() == (*i)) {

						class_options[option->second->Counter()] = option->second;
					}
			}
			std::map<int, shared_ptr<LpOption> >::const_iterator co;
			for (co = class_options.begin(); co != class_options.end(); co++) {
				co->second->OutputShortDescription(reporter);
			}
			reporter.Printf()->info("/n");
		}
	}


}//namespace lpopc