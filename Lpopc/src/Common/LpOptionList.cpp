// Copyright (C) 2014-2015 Xue Zhichen, Wang Yujie,Wang Na 
// All Rights Reserved.
// This file is a part of LPOPC , published under the Eclipse Public License.
// Author:Xue Zhichen 8/4 2015   21:43
// Email:eddy_lpopc@163.com
#include "LpOptionList.hpp"
#include"LpException.hpp"
#include"LpOption.hpp"
#include"LpUtils.h"
#include"LpConf.h"
#include<cctype>
#include<cstdio>
#include <cstdlib>
#include <cstring>

#ifndef _CRT_SECURE_NO_WARNINGS
#define _CRT_SECURE_NO_WARNINGS
#endif
namespace Lpopc
{
	bool LpOptionsList::SetStringValue(const std::string& tag,
		const std::string& value,
		bool allow_clobber, /* = true */
		bool dont_print /* = false */)
	{
		if (reg_options_.get()) {
			shared_ptr<const LpOption>option = reg_options_->GetOption(tag);

			if (!option.get()) {
				std::string msg;
				if (reporter_.get()) {
					msg = "Tried to set Option: " + tag;
					msg += ". It is not a valid option. Please check the list of available options.\n";
					reporter_->Printf()->error(msg) ;
				}
				LP_THROW_EXCEPTION(OPTION_INVALID, msg);
				return false;
			}

			if (option->Type() != String_type) {
				std::string msg;
				if (reporter_.get()) {
					msg= "Tried to set Option: " + tag;
					msg += ". It is a valid option, but it is of type ";
					if (option->Type() == Number_type) {
						msg += " double";
					}
					else if (option->Type() == Integer_type) {
						msg += " Integer";
					}
					else {
						msg += " Unknown";
					}
					msg += ", not of type String. Please check the documentation for options.\n";
					reporter_->Printf()->error(msg) ;
					option->OutputDescription(*reporter_);
				}
				LP_THROW_EXCEPTION(OPTION_INVALID, msg);
				return false;
			}

			if (!option->IsValidStringSetting(value)) {
				std::string msg;
				if (reporter_.get()) {
					msg = "Setting: \"" + value;
					msg += "\" is not a valid setting for Option: ";
					msg += tag;
					msg += ". Check the option documentation.\n";
					reporter_->Printf()->error(msg);
					option->OutputDescription(*reporter_);
				}
				LP_THROW_EXCEPTION(OPTION_INVALID, msg);
				return false;
			}
		}

		if (!will_allow_clobber(tag)) {
			if (reporter_.get()) {
				std::string msg = "WARNING: Tried to set option \"" + tag;
				msg += "\" to a value of \"" + value;
				msg += "\",\n         but the previous value is set to disallow clobbering.\n";
				msg += "         The setting will remain as: \"" + tag;
				msg += " " + options_[lowercase(tag)].GetValue();
				msg += "\"\n";
				reporter_->Printf()->warn(msg);
			}
		}
		else {
			//    if (will_allow_clobber(tag)) {
			LpOptionsList::OptionValue optval(value, allow_clobber, dont_print);
			options_[lowercase(tag)] = optval;
		}
		return true;

		//     std::string msg = "Option: \"" + tag;
		//     msg += " ";
		//     msg += value;
		//     msg += "\" not taken because a value of \n\"" ;
		//     msg += options_[lowercase(tag)].GetValue();
		//     msg += "\" already exists and is set to disallow clobbering.\n\n";
		//     reporter_->Printf(L_ERROR, C_MAIN, msg.c_str());
		//     return false;
	}

	bool LpOptionsList::SetNumericValue(const std::string& tag, double value,
		bool allow_clobber, /* = true */
		bool dont_print /* = false */)
	{
		char buffer[256];
		Snprintf(buffer, 255, "%g", value);

		if (reg_options_.get()) {
			shared_ptr<const LpOption> option (reg_options_->GetOption(tag));

			if (!option.get()) {
				if (reporter_.get()) {
					std::string msg = "Tried to set Option: " + tag;
					msg += ". It is not a valid option. Please check the list of available options.\n";
					reporter_->Printf()->error(msg.c_str());
				}
				//LP_THROW_EXCEPTION(OPTION_INVALID, msg);
				return false;
			}

			if (option->Type() != Number_type) {
				if (reporter_.get()) {
					std::string msg = "Tried to set Option: " + tag;
					msg += ". It is a valid option, but it is of type ";
					if (option->Type() == String_type) {
						msg += " String";
					}
					else if (option->Type() == Integer_type) {
						msg += " Integer";
					}
					else {
						msg += " Unknown";
					}
					msg += ", not of type double. Please check the documentation for options.\n";
					reporter_->Printf()->error( msg.c_str());
					option->OutputDescription(*reporter_);
				}
				//LP_THROW_EXCEPTION(OPTION_INVALID, msg);
				return false;
			}

			if (!option->IsValidNumberSetting(value)) {
				if (reporter_.get()) {
					std::string msg = "Setting: \"";
					msg += buffer;
					msg += "\" is not a valid setting for Option: ";
					msg += tag;
					msg += ". Check the option documentation.\n";
					reporter_->Printf()->error(msg) ;
					option->OutputDescription(*reporter_);
				}
				//LP_THROW_EXCEPTION(OPTION_INVALID, msg);
				return false;
			}
		}

		if (!will_allow_clobber(tag)) {
			if (reporter_.get()) {
				std::string msg = "WARNING: Tried to set option \"" + tag;
				msg += "\" to a value of \"";
				msg += buffer;
				msg += "\",\n         but the previous value is set to disallow clobbering.\n";
				msg += "         The setting will remain as: \"" + tag;
				msg += " " + options_[lowercase(tag)].GetValue();
				msg += "\"\n";
				reporter_->Printf()->warn(msg);
			}
		}
		else {
			LpOptionsList::OptionValue optval(buffer, allow_clobber, dont_print);
			options_[lowercase(tag)] = optval;
		}
		return true;
	}

	bool LpOptionsList::SetIntegerValue(const std::string& tag, int value,
		bool allow_clobber, /* = true */
		bool dont_print /* = false */)
	{
		char buffer[256];
		Snprintf(buffer, 255, "%d", value);

		if (reg_options_.get()) {
			shared_ptr<const LpOption> option = reg_options_->GetOption(tag);

			if (!option.get()) {
				std::string msg = "Tried to set Option: " + tag;
				msg += ". It is not a valid option. Please check the list of available options.\n";
				if (reporter_.get()) {
					reporter_->Printf()->error(msg) ;
				}
				LP_THROW_EXCEPTION(OPTION_INVALID, msg);
				return false;
			}

			if (option->Type() != Integer_type) {
				if (reporter_.get()) {
					std::string msg = "Tried to set Option: " + tag;
					msg += ". It is a valid option, but it is of type ";
					if (option->Type() == String_type) {
						msg += " String";
					}
					else if (option->Type() == Number_type) {
						msg += " double";
					}
					else {
						msg += " Unknown";
					}
					msg += ", not of type Integer. Please check the documentation for options.\n";
					reporter_->Printf()->error(msg);
					option->OutputDescription(*reporter_);
				
				LP_THROW_EXCEPTION(OPTION_INVALID, msg);
				}
				return false;
			}

			if (!option->IsValidIntegerSetting(value)) {
				if (reporter_.get()) {
					std::string msg = "Setting: \"";
					msg += buffer;
					msg += "\" is not a valid setting for Option: ";
					msg += tag;
					msg += ". Check the option documentation.\n";
					reporter_->Printf()->error(msg);
					option->OutputDescription(*reporter_);
				
				LP_THROW_EXCEPTION(OPTION_INVALID, msg);
				}
				return false;
			}
		}

		if (!will_allow_clobber(tag)) {
			if (reporter_.get()) {
				std::string msg = "WARNING: Tried to set option \"" + tag;
				msg += "\" to a value of \"";
				msg += buffer;
				msg += "\",\n         but the previous value is set to disallow clobbering.\n";
				msg += "         The setting will remain as: \"" + tag;
				msg += " " + options_[lowercase(tag)].GetValue();
				msg += "\"\n";
				reporter_->Printf()->warn(msg) ;
			}
		}
		else {
			//    if (will_allow_clobber(tag)) {
			LpOptionsList::OptionValue optval(buffer, allow_clobber, dont_print);
			options_[lowercase(tag)] = optval;
		}
		return true;
	}

	bool LpOptionsList::SetStringValueIfUnset(const std::string& tag,
		const std::string& value,
		bool allow_clobber, /* = true */
		bool dont_print /* = false */)
	{
		std::string val;
		bool found = GetStringValue(tag, val, "");
		if (!found) {
			return SetStringValue(tag, value, allow_clobber, dont_print);
		}
		return true;
	}

	bool LpOptionsList::SetNumericValueIfUnset(const std::string& tag,
		double value,
		bool allow_clobber, /* = true */
		bool dont_print /* = false */)
	{
		double val;
		bool found = GetNumericValue(tag, val, "");
		if (!found) {
			return SetNumericValue(tag, value, allow_clobber, dont_print);
		}
		return true;
	}

	bool LpOptionsList::SetIntegerValueIfUnset(const std::string& tag,
		int value,
		bool allow_clobber, /* = true */
		bool dont_print /* = false */)
	{
		int val;
		bool found = GetIntegerValue(tag, val, "");
		if (!found) {
			return SetIntegerValue(tag, value, allow_clobber, dont_print);
		}
		return true;
	}

	bool LpOptionsList::GetStringValue(const std::string& tag, std::string& value,
		const std::string& prefix) const
	{
		shared_ptr<const LpOption> option;

		bool found = find_tag(tag, prefix, value);

		if (reg_options_.get()) {
			option = reg_options_->GetOption(tag);
			if (!option.get()) {
				std::string msg = "IPOPT tried to get the value of Option: " + tag;
				msg += ". It is not a valid registered option.";
				LP_THROW_EXCEPTION(OPTION_INVALID, msg);
			}

			if (option->Type() != String_type) {
				std::string msg = "IPOPT tried to get the value of Option: " + tag;
				msg += ". It is a valid option, but it is of type ";
				if (option->Type() == Integer_type) {
					msg += " Integer";
				}
				else if (option->Type() == Number_type) {
					msg += " double";
				}
				else {
					msg += " Unknown";
				}
				msg += ", not of type String. Please check the documentation for options.";
				if (reporter_.get()) {
					option->OutputDescription(*reporter_);
				}
				LP_THROW_EXCEPTION(OPTION_INVALID, msg);
			}

			if (found) {
				value = option->MapStringSetting(value);
			}
			else {
				value = option->DefaultString();
			}
		}

		return found;
	}

	bool LpOptionsList::GetEnumValue(const std::string& tag, int& value,
		const std::string& prefix) const
	{
		std::string str;
		shared_ptr<const LpOption> option ;

		bool found = find_tag(tag, prefix, str);

		if (reg_options_.get()) {
			option = reg_options_->GetOption(tag);
			if (!option.get()) {
				std::string msg = "IPOPT tried to get the value of Option: " + tag;
				msg += ". It is not a valid registered option.";
				LP_THROW_EXCEPTION(OPTION_INVALID, msg);
			}

			if (option->Type() != String_type) {
				std::string msg = "IPOPT tried to get the value of Option: " + tag;
				msg += ". It is a valid option, but it is of type ";
				if (option->Type() == Integer_type) {
					msg += " Integer";
				}
				else if (option->Type() == Number_type) {
					msg += " double";
				}
				else {
					msg += " Unknown";
				}
				msg += ", not of type String. Please check the documentation for options.";
				if (reporter_.get()) {
					option->OutputDescription(*reporter_);
				}
				LP_THROW_EXCEPTION(OPTION_INVALID, msg);
			}

			if (found) {
				value = option->MapStringSettingToEnum(str);
			}
			else {
				value = option->DefaultStringAsEnum();
			}
		}

		return found;
	}

	bool LpOptionsList::GetBoolValue(const std::string& tag, bool& value,
		const std::string& prefix) const
	{
		std::string str;
		bool ret = GetStringValue(tag, str, prefix);
		if (str == "no" || str == "false" || str == "off") {
			value = false;
		}
		else if (str == "yes" || str == "true" || str == "on") {
			value = true;
		}
		else {
			LP_THROW_EXCEPTION(OPTION_INVALID, "Tried to get a boolean from an option and failed.");
		}

		return ret;
	}

	bool LpOptionsList::GetNumericValue(const std::string& tag, double& value,
		const std::string& prefix) const
	{
		shared_ptr<const LpOption> option ;

		if (reg_options_.get()) {
			option = reg_options_->GetOption(tag);
			if (!option.get()) {
				std::string msg = "IPOPT tried to get the value of Option: " + tag;
				msg += ". It is not a valid registered option.";
				LP_THROW_EXCEPTION(OPTION_INVALID, msg);
			}

			if (option->Type() != Number_type) {
				std::string msg = "IPOPT tried to get the value of Option: " + tag;
				msg += ". It is a valid option, but it is of type ";
				if (option->Type() == Integer_type) {
					msg += " Integer";
				}
				else if (option->Type() == String_type) {
					msg += " String";
				}
				else {
					msg += " Unknown";
				}
				msg += ", not of type double. Please check the documentation for options.";
				if (reporter_.get()) {
					option->OutputDescription(*reporter_);
				}
				LP_THROW_EXCEPTION(OPTION_INVALID, msg);
			}
		}

		std::string strvalue;
		if (find_tag(tag, prefix, strvalue)) {
			// Some people like to use 'd' instead of 'e' in floating point
			// doubles.  Therefore, we change a 'd' to an 'e'
			char* buffer = new char[strvalue.length()+1];
			strcpy(buffer, strvalue.c_str());
			for (int i=0; i<(int)strvalue.length(); ++i) {
				if (buffer[i]=='d' || buffer[i]=='D') {
					buffer[i] = 'e';
				}
			}
			char* p_end;
			double retval = strtod(buffer, &p_end);
			if (*p_end!='\0' && !isspace(*p_end)) {
				delete [] buffer;
				std::string msg = "Option \"" + tag +
					"\": Double value expected, but non-numeric value \"" +
					strvalue+"\" found.\n";
				LP_THROW_EXCEPTION(OPTION_INVALID, msg);
			}
			delete [] buffer;
			value = retval;
			return true;
		}
		else if (option.get()) {
			value = option->DefaultNumber();
			return false;
		}
		return false;
	}

	bool LpOptionsList::GetIntegerValue(const std::string& tag, int& value,
		const std::string& prefix) const
	{
		shared_ptr<const LpOption> option ;

		if (reg_options_.get()) {
			option = reg_options_->GetOption(tag);
			if (!option.get()) {
				std::string msg = "IPOPT tried to get the value of Option: " + tag;
				msg += ". It is not a valid registered option.";
				LP_THROW_EXCEPTION(OPTION_INVALID, msg);
			}

			if (option->Type() != Integer_type) {
				std::string msg = "IPOPT tried to get the value of Option: " + tag;
				msg += ". It is a valid option, but it is of type ";
				if (option->Type() == Number_type) {
					msg += " double";
				}
				else if (option->Type() == String_type) {
					msg += " String";
				}
				else {
					msg += " Unknown";
				}
				msg += ", not of type Integer. Please check the documentation for options.";
				if (reporter_.get()) {
					option->OutputDescription(*reporter_);
				}
				LP_THROW_EXCEPTION(OPTION_INVALID, msg);
			}
		}

		std::string strvalue;
		if (find_tag(tag, prefix, strvalue)) {
			char* p_end;
			size_t retval = strtol(strvalue.c_str(), &p_end, 10);
			if (*p_end!='\0' && !isspace(*p_end)) {
				std::string msg = "Option \"" + tag +
					"\": Integer value expected, but non-integer value \"" +
					strvalue+"\" found.\n";
				LP_THROW_EXCEPTION(OPTION_INVALID, msg);
			}
			value = static_cast<int>(retval);
			return true;
		}
		else if (option.get()) {
			value = option->DefaultInteger();
			return false;
		}

		return false;
	}

	const std::string& LpOptionsList::lowercase(const std::string tag) const
	{
		lowercase_buffer_ = tag;
		for (int i=0; i<(int)tag.length(); i++) {
			lowercase_buffer_[i] = (char)tolower(tag[i]);
		}
		return lowercase_buffer_;
	}

	void LpOptionsList::PrintList(std::string& list) const
	{
		list.erase();
		char buffer[256];
		Snprintf(buffer, 255, "%40s   %-20s %s\n", "Name", "Value", "# times used");
		list += buffer;
		for (std::map< std::string, OptionValue >::const_iterator p = options_.begin();
			p != options_.end();
			p++ ) {
				Snprintf(buffer, 255, "%40s = %-20s %6d\n", p->first.c_str(),
					p->second.Value().c_str(), p->second.Counter());
				list += buffer;
		}
	}

	void LpOptionsList::PrintUserOptions(std::string& list) const
	{
		list.erase();
		char buffer[256];
		Snprintf(buffer, 255, "%40s   %-20s %s\n", "Name", "Value", "used");
		list += buffer;
		for (std::map< std::string, OptionValue >::const_iterator p = options_.begin();
			p != options_.end();
			p++ ) {
				if (!p->second.DontPrint()) {
					const char yes[] = "yes";
					const char no[] = "no";
					const char* used;
					if (p->second.Counter()>0) {
						used = yes;
					}
					else {
						used = no;
					}
					Snprintf(buffer, 255, "%40s = %-20s %4s\n", p->first.c_str(),
						p->second.Value().c_str(), used);
					list += buffer;
				}
		}
	}

	bool LpOptionsList::ReadFromStream(const LpReporter& reporter,
		std::istream& is)
	{
		reporter.Printf()->trace( "Start reading options from stream.\n");

		while (true) {
			std::string tag;
			std::string value;

			if (!readnexttoken(is, tag)) {
				// That's it - end of file reached.
				reporter.Printf()->trace(
					"Finished reading options from file.\n");
				return true;
			}

			if (!readnexttoken(is, value)) {
				// Can't read value for a given tag
				reporter.Printf()->error(
					"Error reading value for tag {} from file.\n",
					tag.c_str());
				return false;
			}

			// Now add the value for the options list
			reporter.Printf()->trace(
				"Adding option \"{}\" with value \"{}\" to LpOptionsList.\n",
				tag.c_str(), value.c_str());

			if (reg_options_.get()) {
				shared_ptr<const LpOption> option = reg_options_->GetOption(tag);
				if (option.get()) {
					std::string msg = "Read Option: \"";
					msg += tag;
					msg += "\". It is not a valid option. Check the list of available options.";
					LP_THROW_EXCEPTION(OPTION_INVALID, msg);
				}

				if (option->Type() == String_type) {
					bool result = SetStringValue(tag, value, false);
					LP_ASSERT_EXCEPTION(result, OPTION_INVALID,
						"Error setting string value read from option file.");
				}
				else if (option->Type() == Number_type) {
					// Some people like to use 'd' instead of 'e' in floating
					// point doubles.  Therefore, we change a 'd' to an 'e'
					char* buffer = new char[value.length()+1];
					strcpy(buffer, value.c_str());
					for (int i=0; i<(int)value.length(); ++i) {
						if (buffer[i]=='d' || buffer[i]=='D') {
							buffer[i] = 'e';
						}
					}
					char* p_end;
					double retval = strtod(buffer, &p_end);
					if (*p_end!='\0' && !isspace(*p_end)) {
						delete [] buffer;
						std::string msg = "Option \"" + tag +
							"\": Double value expected, but non-numeric option value \"" +
							value + "\" found.\n";
						LP_THROW_EXCEPTION(OPTION_INVALID, msg);
					}
					delete [] buffer;
					bool result = SetNumericValue(tag, retval, false);
					LP_ASSERT_EXCEPTION(result, OPTION_INVALID,
						"Error setting numeric value read from file.");
				}
				else if (option->Type() == Integer_type) {
					char* p_end;
					size_t retval = strtol(value.c_str(), &p_end, 10);
					if (*p_end!='\0' && !isspace(*p_end)) {
						std::string msg = "Option \"" + tag +
							"\": Integer value expected, but non-integer option value \"" +
							value + "\" found.\n";
						if (reporter_.get()) {
							option->OutputDescription(*reporter_);
						}
						LP_THROW_EXCEPTION(OPTION_INVALID, msg);
					}
					bool result = SetIntegerValue(tag, static_cast<int>(retval), false);
					LP_ASSERT_EXCEPTION(result, OPTION_INVALID,
						"Error setting integer value read from option file.");
				}
				else {
					assert(false && "Option Type: Unknown");
				}
			}
			else {
				bool result = SetStringValue(tag, value, false);
				LP_ASSERT_EXCEPTION(result, OPTION_INVALID,
					"Error setting value read from option file.");
			}
		}
	}

	bool LpOptionsList::find_tag(const std::string& tag,
		const std::string& prefix,
		std::string& value) const
	{
		bool found=false;
		std::map< std::string, OptionValue >::const_iterator p;

		if (prefix != "") {
			p = options_.find(lowercase(prefix+tag));
			if (p != options_.end()) {
				found = true;
			}
		}

		if (!found) {
			p = options_.find(lowercase(tag));
			if (p != options_.end()) {
				found = true;
			}
		}

		if (found) {
			value = p->second.GetValue();
		}

		return found;
	}

	bool LpOptionsList::will_allow_clobber(const std::string& tag) const
	{
		bool allow_clobber=true;
		std::map< std::string, OptionValue >::const_iterator p;

		p = options_.find(lowercase(tag));
		if (p != options_.end()) {
			allow_clobber = p->second.AllowClobber();
		}

		return allow_clobber;
	}

	bool LpOptionsList::readnexttoken(std::istream& is, std::string& token)
	{
		token.erase();
		int c = is.get();

		// First get rid of all comments and white spaces
		while (!is.eof() && (isspace(c) || c=='#') ) {
			if (c=='#') {
				is.ignore(10000000, '\n');
			}
			c=is.get();
		}

		bool inside_quotes = (c=='"');
		if (inside_quotes) {
			if (is.eof()) return false; // eof after quotation symbol
			c=is.get();
		}

		if( is.eof() )
			return false;

		// Now read the token
		while (!is.eof() && (inside_quotes || !isspace(c))) {
			token += (char)c;
			c = is.get();
			if (inside_quotes && (c=='"')) {
				inside_quotes = false;
				if (!is.eof())
					c = is.get();
			}
		}

		return !inside_quotes;
	}
}//namespace Lpopc