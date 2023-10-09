/* 
 * File:   OptionParser.cpp
 * Author: John
 * 
 * Created on April 27, 2015, 6:07 PM
 */
#include <pch.h>
#include <iostream>
#include <sstream>
#include "options/OptionParser.h"

using namespace std;
ArgvParser OptionParser::cmd;

void OptionParser::addOption(const string& option, const string& desc){
    cmd.defineOption(option, desc, ArgvParser::NoOptionAttribute);
}

void OptionParser::setup(){
    //define error codes
    cmd.addErrorCode(0, "Success");
    cmd.addErrorCode(1, "Error");
}

void OptionParser::parseOptions(int argc, char** argv){
    int result = cmd.parse(argc, argv);
    
    if (result != ArgvParser::NoParserError) {
      cout << "Parse ERROR!!!: " << cmd.parseErrorDescription(result) << endl;
      exit(1);
    }
}

string OptionParser::optionValue(const string& option){
    return cmd.optionValue(option);
}

bool OptionParser::foundOption(const string& option){
    return cmd.foundOption(option);
}

OptionParser::OptionParser() {
	this->myCmd.addErrorCode(0, "Success");
	this->myCmd.addErrorCode(1, "Error");
}

void OptionParser::add(const string& option, const string& desc) {
	myCmd.defineOption(option, desc, ArgvParser::NoOptionAttribute);
}

int OptionParser::parse(string optionsText) {
	vector<char *> args;
	istringstream iss(optionsText);
	string token;
	int optionsCount = 0;
	while (iss >> token) {
		char *arg = new char[token.size() + 1];
		copy(token.begin(), token.end(), arg);
		arg[token.size()] = '\0';
		args.push_back(arg);
		optionsCount++;
	}
	args.push_back(0);
	int result = myCmd.parse(optionsCount, &args[0]);

	if (result != ArgvParser::NoParserError) {
		cout << "Parse ERROR!!!: " << myCmd.parseErrorDescription(result) << endl;
		return 0;
	}
	return 1;
}
string OptionParser::value(const string& option) {
	return myCmd.optionValue(option);
}
bool OptionParser::has(const string& option) {
	return myCmd.foundOption(option);
}

double OptionParser::parseDouble(const string& option) {
	return parseDouble_(cmd.optionValue(option));
}

int OptionParser::parseInt(const string& option) {
	return parseInt_(cmd.optionValue(option));
}

int OptionParser::parseInt_(const string& str) {
	int aux;
	sscanf(str.c_str(), "%d", &aux);
	return aux;
}

double OptionParser::parseDouble_(const string& str) {
	double aux;
	sscanf(str.c_str(), "%lf", &aux);
	return aux;
}

vector<IndexedValue> OptionParser::parseIndexedValues(const string& option) {
	vector<IndexedValue> values;
	vector<string> tokens;
	CommandLineProcessing::splitString(tokens, cmd.optionValue(option), "/");

	for (unsigned int i = 0; i < tokens.size(); i++) {
		string pair_str = tokens[i];
		vector<string> subtokens;
		CommandLineProcessing::splitString(subtokens, pair_str, "_");
		if (subtokens.size() != 2) {
			cout << "Parse ERROR!!!: Could not parse indexed value." << endl;
			return values;
		}
		IndexedValue value;
		value.index = parseInt_(subtokens[0]);
		value.value = parseDouble_(subtokens[1]);
		values.push_back(value);
	}

	return values;
}
