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

#include "platform.h"

#ifdef WIN32_LEAN_AND_MEAN

#else
	#include <unistd.h>
#endif

#include <iostream>
#include <cerrno>
#include <cstring>

namespace GIMLI{

long numberOfCPU(){
    long nprocs = -1;
    long nprocs_max = -1;
#if defined(WINDOWS) || defined(_WIN32)
    #ifndef _SC_NPROCESSORS_ONLN
        SYSTEM_INFO info;
        GetSystemInfo(&info);
    #define sysconf(a) info.dwNumberOfProcessors
    #define _SC_NPROCESSORS_ONLN
    #endif
#endif // no windows

#ifdef _SC_NPROCESSORS_ONLN
    nprocs = sysconf(_SC_NPROCESSORS_ONLN);

    if ( nprocs < 1 ) {
        std::cerr << "Could not determine number of CPUs online:" << std::strerror(errno) << std::endl;
    } else {
        nprocs = 1;
    }
    nprocs_max = sysconf(_SC_NPROCESSORS_CONF);

    if ( nprocs_max < 1 ) {
        std::cerr << "Could not determine number of CPUs configured:" << std::strerror(errno) << std::endl;
    }
#else
    std::cerr << "Could not determine number of CPUs (!_SC_NPROCESSORS_ONLN)" << std::endl;
    nprocs = 1;
    return 1;
#endif
    return nprocs_max;
}

int schedGetCPU(){
#if defined(WINDOWS) || defined(_WIN32)
    return -1;
#elif defined(__APPLE__)
    return -1;
#else
    return sched_getcpu();
#endif
}

} // namespace GIMLI{
