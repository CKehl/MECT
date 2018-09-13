/*
 * main.hpp
 *
 *  Created on: Apr 5, 2018
 *      Author: christian
 */

#ifndef MAIN_HPP_
#define MAIN_HPP_

#include <boost/regex.hpp>

#include <boost/program_options.hpp>
#include <boost/program_options/cmdline.hpp>
#include <boost/program_options/errors.hpp>
#include <boost/program_options/option.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/config.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/positional_options.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/program_options/version.hpp>
#include <boost/filesystem.hpp>

#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/algorithm/string/classification.hpp>


using namespace boost::program_options;

#include "Utilities/imaging/src/image.hpp"
#include "Utilities/imaging/src/image3D.hpp"
#include "Utilities/imaging/src/image4D.hpp"
#include "Utilities/imaging/src/mhdImage.hpp"
//#include "Utilities/io/src/MHD.hpp"
#include "Utilities/io/src/io_commons.hpp"
#include "Utilities/common/src/std_typedefs.h"
#include "GraphCut4D.hpp"
#include <iostream>
#include <fstream>
#include <vector>
#include <time.h>
#include <cmath>
#include <limits>



#endif /* MAIN_HPP_ */
