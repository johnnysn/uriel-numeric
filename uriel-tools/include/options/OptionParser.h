/* 
 * File:   OptionParser.h
 * Author: John
 *
 * Created on April 27, 2015, 6:07 PM
 */

#ifndef OPTIONPARSER_H
#define	OPTIONPARSER_H

#include <string>
#include "options/argvparser.h"

using namespace std;
using namespace CommandLineProcessing;

class OptionParser {
public:
	OptionParser();

    static void setup();
    static void addOption(const string& option, const string& desc);
    static void parseOptions(int argc, char** argv);
    static string optionValue(const string& option);
    static bool foundOption(const string& option);

	void add(const string& option, const string& desc);
	int parse(string optionsText);
	string value(const string& option);
	bool has(const string& option);

	static double parseDouble(const string& option);
	static int parseInt(const string& option);
private:
    static ArgvParser cmd;
	ArgvParser myCmd;
};

#endif	/* OPTIONPARSER_H */

