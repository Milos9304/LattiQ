/*
 * lattiq.h
 *
 *  Created on: Apr 19, 2021
 *      Author: Milos Prokop
 */

#ifndef LATTIQ_LATTIQ_H_
#define LATTIQ_LATTIQ_H_
/*
#include "io/logger.h"
#include "io/indicators/progress_bar.hpp"

using namespace indicators;

const Color colors[9] = {Color::red, Color::green, Color::yellow, Color::blue, Color::magenta, Color::cyan};

#define bar_opts(counter, num_lattices, lattice_name, optionsName)   option::BarWidth{50},\
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
														option::MaxProgress{optionsName->max_iters}


std::string process(std::string const& s);*/

#endif /* LATTIQ_LATTIQ_H_ */
