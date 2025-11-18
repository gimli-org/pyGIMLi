/******************************************************************************
 *   Copyright (C) 2006-2024 by the GIMLi development team                    *
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

#ifndef _GIMLI_PLATFORM__H
#define _GIMLI_PLATFORM__H

#if defined(WINDOWS) || defined(_WIN32) || defined(WIN32)
    //#error checkme

    #define PATHSEPARATOR "\\"

	#define WIN32_LEAN_AND_MEAN
    #define BOOST_USE_WINDOWS_H
    //#define BOOST_PLAT_WINDOWS_RUNTIME
    #ifndef DLLEXPORT
        #if defined(DLL_EXPORT) || defined(gimli_EXPORTS)
            //* We are building this library
            //#warning (GIMLI We are building this library)
            #define DLLEXPORT __declspec(dllexport)
        #else
            //* We are using this library
            //#warning (GIMLI We are using this library)
            #define DLLEXPORT __declspec(dllimport)
        #endif
    #endif

    // Don't let win32api windef.h define min and max as macros
    // if included after c++config.h.
    #ifndef NOMINMAX
        #define NOMINMAX
    #endif

    #include <stddef.h>

    #ifndef PYGIMLI_CAST
        #include <windows.h>
    #endif

    #undef near
    #undef far

#else /* ELSE NO WINDOWS */
    #define PATHSEPARATOR "/"
    #define DLLEXPORT
#endif /* NO WINDOWS */

/*!
 * Return the numbers of virtual CPUS on this system
 */

#include <cmath>

namespace GIMLI{

/*!Return the number of available CPUs. 1 if unknown.*/
DLLEXPORT long numberOfCPU();

/*!Return the number of the currently used CPU from the scheduler.*/
DLLEXPORT int schedGetCPU();

// Microsoft Visual C++ 10 does not provide some C99 functions
#if defined(_MSC_VER)
template< typename T > T rint( T x ){ return std::floor(x + 0.5); }
template< typename T > inline bool isinf(T value){
    return value >= std::numeric_limits<T>::min() &&
	       value <= std::numeric_limits<T>::max();}
template< typename T > inline bool isnan(T value){return (value) != (value);}
#elif defined(_WIN32) // mingw
template< typename T > inline bool isinf(T value){return std::isinf(value);}
template< typename T > inline bool isnan(T value){return std::isnan(value);}
#else
template< typename T > inline bool isinf(T value){return std::isinf(value);}
template< typename T > inline bool isnan(T value){return std::isnan(value);}
#endif


}

#if defined ( __APPLE__ ) || ( defined (__SVR4) && defined (__sun) )
// #include <ieeefp.h>
// inline int isinf( double x ){ return !finite(x) && x==x; }
                              // /*
#include <unistd.h>
#include <sys/types.h>

//   Apple (OS X) and Sun systems declare getopt in unistd.h,
//   other systems (Linux) use getopt.h
// */
// #if defined ( __APPLE__ ) || ( defined (__SVR4) && defined (__sun) )
// #include <unistd.h>
// struct option {
// # if defined __STDC__ && __STDC__
//   const char *name;
// # else
//   char *name;
// # endif
//   /* has_arg can't be an enum because some compilers complain about
//      type mismatches in all the code that assumes it is an int.  */
//   int has_arg;
//   int *flag;
//   int val;
// };
// #else

#endif

#endif // _GIMLI_PLATFORM__H
