/*
 * logger.h
 *
 *  Created on: Feb 11, 2021
 *      Author: Milos Prokop
 */

#ifndef SRC_LOGGER_H_
#define SRC_LOGGER_H_

#include <string>
#include <sstream>
#include <iostream>
#include "indicators/termcolor.hpp"
#include "indicators/progress_bar.hpp"

using namespace indicators;

const Color colors[9] = {Color::red, Color::green, Color::yellow, Color::blue, Color::magenta, Color::cyan};

/*#define bar_opts(counter, num_lattices, lattice_name, optionsName)   option::BarWidth{50},\
					option::Start{"["},\
					option::Fill{"="},\
					option::Lead{">"},\
					option::Remainder{" "},\
					option::End{"]"},\
					option::PrefixText{std::to_string(counter+1) + "/" + std::to_string(num_lattices) + " " + lattice_name},\
					option::ForegroundColor{colors[counter % 7]},\
					option::ShowElapsedTime{true},\
					option::ShowRemainingTime{true},\
					option::FontStyles{std::vector<FontStyle>{FontStyle::bold}},\
					option::MaxProgress{optionsName->max_iters}*/

#define LOG_LEVEL 0 //0 - debug, 1 - info, 2 - warning, 3 - error

inline std::string getCurrentDateTime(){

    time_t now = time(0);
    struct tm  tstruct;
    char  buf[80];
    tstruct = *localtime(&now);
    strftime(buf, sizeof(buf), "%Y-%m-%d %X", &tstruct);

    return std::string(buf);

};
template<typename Functor>
inline void printLog(std::string logType, Functor termColor, std::string logMsg){

	std::string now = getCurrentDateTime();
	std::cerr << termColor << "[[" << now << "]][" << logType << "] " << logMsg << termcolor::reset << std::endl;

}

inline void logd( std::string logMsg, int log_level = LOG_LEVEL ){

	if(log_level == 0)
		printLog("DEBUG", termcolor::cyan, logMsg);

}

inline void logi( std::string logMsg, int log_level = LOG_LEVEL ){

	if(log_level <= 1)
		printLog("INFO", termcolor::white, logMsg);

}

inline void logw( std::string logMsg, int log_level = LOG_LEVEL ){

	if(log_level <= 2)
		printLog("WARNING", termcolor::yellow, logMsg);

}

inline void loge( std::string logMsg ){

	printLog("ERROR", termcolor::red, logMsg);

}

inline void throw_runtime_error( std::string logMsg ){

	printLog("ERROR", termcolor::red, logMsg);
	throw std::runtime_error(logMsg);

}


#endif /* SRC_LOGGER_H_ */
