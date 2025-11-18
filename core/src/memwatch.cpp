/******************************************************************************
 *   Copyright (C) 2012-2024 by the GIMLi development team                    *
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

#include "memwatch.h"
#include "stopwatch.h"

#include <iostream>

#ifdef WIN32_LEAN_AND_MEAN
    #include <psapi.h>
#elif READPROC_FOUND || defined(HAVE_PROC_READPROC)
    #include <proc/readproc.h>
    #define USE_PROC_READPROC TRUE
    #error readproc
#else
    #define USE_PROC_READPROC 0
#endif

#include <mutex>
static std::mutex __readproc__mutex__;

#ifdef WIN32_LEAN_AND_MEAN
    #include <windows.h>
#else
    #include <unistd.h>
#endif

namespace GIMLI {

// Global static pointer used to ensure a single instance of the class.
template < > DLLEXPORT MemWatch * Singleton < MemWatch>::pInstance_ = NULL;

MemWatch::MemWatch(){
    last_ = inUse();
    swatchAll_ = new Stopwatch(true);
    swatchDur_ = new Stopwatch(true);
}

MemWatch::~MemWatch(){
    delete swatchAll_; swatchAll_ = NULL;
    delete swatchDur_; swatchDur_ = NULL;
}

double MemWatch::current(){
    double ret = inUse() - last_;
    last_ = inUse();
    return ret;
}

double MemWatch::inUse() {
#if USE_BOOST_THREAD
    boost::mutex::scoped_lock lock(__readproc__mutex__);
#else
    std::unique_lock < std::mutex > lock(__readproc__mutex__);
        //std::lock_guard< std::mutex > lock(__readproc__mutex__);
#endif

#ifdef WIN32_LEAN_AND_MEAN

    PROCESS_MEMORY_COUNTERS pmc;
    if (GetProcessMemoryInfo(GetCurrentProcess(), &pmc, sizeof(pmc))){
        double ret = MByte(pmc.WorkingSetSize);
        return ret;
    } else { return 0; }
#else // no WINDOWS

#if USE_PROC_READPROC
    struct proc_t usage;
    look_up_our_self(& usage);
    // vsize .. number of pages of virtual memory ...
     //__MS("vsize: " << usage.vsize/1024) // virtual memory
//      __MS("size: " << usage.size/1024)//   total # of pages of memory
//      __MS("resident: " << usage.resident/1024)//         number of resident set (non-swapped) pages (4k)
//      __MS("share: " << usage.share/1024)
//      __MS("rss: " << usage.rss/1024)
    double ret = MByte(usage.vsize);
    return ret;
#else
    // no windows and no libproc
#endif


#endif // no WINDOWS
    return 0;
}

void MemWatch::info(const std::string & str){
    if (debug()){
#if defined(WIN32_LEAN_AND_MEAN) || USE_PROC_READPROC
        std::cout << "\t" << str << " Memory "
#if USE_BOOST_THREAD
                    << "(mt)"
#else
                    << "(st)"
#endif // HAVE_BOOST_THREAD_HPP
                    << " in use: abs: " << inUse() << " rel: "
                    << current() << " MByte. t = "
                    << swatchAll_->duration() << "/" << swatchDur_->duration(true) << " s " <<  std::endl;
    #else
        std::cout << "\t" << str << " Memory no info"  <<  std::endl;
    #endif
    }
}

long maxMem(){
#if defined(WIN32_LEAN_AND_MEAN)
    MEMORYSTATUSEX status;
    status.dwLength = sizeof(status);
    GlobalMemoryStatusEx(&status);
    return status.ullTotalPhys;
#else
    long pages = sysconf(_SC_PHYS_PAGES);
    long page_size = sysconf(_SC_PAGE_SIZE);
    return pages * page_size;
#endif
}

} // namespace GIMLI
