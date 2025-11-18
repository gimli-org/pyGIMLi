/******************************************************************************
 *   Copyright (C) 2007-2024 by the GIMLi development team                    *
 *   Carsten Rücker carsten@resistivity.net                                   *
 *                                                                            *
 *   Licensed under the Apache License, Version 2.0 (the "License");          *
 *   you may not use this file except in compliance with the License.         *
 *   You may obtain a copy of the License at                                  *
 *                                                                            *
 *       http://www.apache.org/licenses/LICENSE-2.0                           *
 *                                                                            *
 *   Unless required by applicable law or agreed to in writing, software      *
 *   distributed under the License is distributed on an "AS IS" BASIS,        *
 *   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. *
 *   See the License for the specific language governing permissions and      *
 *   limitations under the License.                                           *
 *                                                                            *
 ******************************************************************************/

#ifndef _GIMLI_EXITCODES__H
#define _GIMLI_EXITCODES__H

// http://gcc.gnu.org/onlinedocs/libstdc++/latest-doxygen/a00653.html
// here we will define some std::exception instead if exit codes

#include <stdexcept>

namespace GIMLI {

//#define USE_EXIT_CODES 1


inline void throwToImplement(const std::string & errString){
#ifndef USE_EXIT_CODES
    throw std::length_error(errString);
#else
    std::cerr << errString << std::endl;
#endif
}

inline void throwRangeError(const std::string & errString, int idx, int low, int high){
    std::stringstream str(errString);
    str << " " << idx << " [" << low << ".." << high << ")" << std::endl;
#ifndef USE_EXIT_CODES
    throw std::out_of_range(str.str());
#else
    std::cerr << str.str() << std::endl;
#endif
}

inline void throwLengthError(const std::string & errString){
#ifndef USE_EXIT_CODES
    throw std::length_error(errString);
#else
    std::cerr << errString << std::endl;
#endif
}

DLLEXPORT bool debug();

inline void throwError(const std::string & errString){
#ifndef USE_EXIT_CODES
    if (debug()){
        std::cerr << "Debug: " << errString << std::endl;
    } 
    throw std::length_error(errString);
#else
    std::cerr << errString << std::endl;
#endif
}

} // namespace GIMLI

#endif // _GIMLI_EXITCODES__H
