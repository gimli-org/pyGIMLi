/******************************************************************************
 *   Copyright (C) 2005-2024 by the GIMLi development team                    *
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

#ifndef _GIMLI_CALCULATE_MULTI_THREAD__H
#define _GIMLI_CALCULATE_MULTI_THREAD__H

#include "gimli.h"

#ifdef USE_BOOST_THREAD
    #include <boost/thread.hpp>
#else
    #include <mutex>
    #include <thread>
#endif

namespace GIMLI{

class BaseCalcMT{
public:
    BaseCalcMT(bool verbose=false)
        : verbose_(verbose), start_(0), end_(0), _threadNumber(0){
    }

    virtual ~BaseCalcMT(){ }

    void operator () () { calc(); }

    void setRange(Index start, Index end, Index threadNumber=0){
        start_ = start;
        end_ = end;
        _threadNumber = threadNumber;
    }

    virtual void calc()=0;

    Index start() const { return start_;}
    Index end() const { return end_;}
protected:
    bool verbose_;
    Index start_;
    Index end_;
    Index _threadNumber;
};
template < class T > void distributeCalc(T calc, uint nCalcs, uint nThreads, bool verbose=false){
    log(Debug, "Create distributed calculation of " + str(nCalcs) + " jobs on "
        + str(nThreads) + " threads for " + str(numberOfCPU()) + " CPU");

    if (nThreads == 1){
        calc.setRange(0, nCalcs);
        Stopwatch swatch(true);
        calc();
        log(Debug, "time: " + str(swatch.duration()) + "s");
    } else {
        uint singleCalcCount = (uint)ceil((double)nCalcs / (double)nThreads);

        std::vector < T > calcObjs;
        for (uint i = 0; i < nThreads; i ++){
            calcObjs.push_back(calc);
            Index start = singleCalcCount * i;
            Index end   = min(singleCalcCount * (i + 1), nCalcs);
            log(Debug, "Threaded calculation: #" + str(i) + ": " + str(start)  +" " + str(end));
            calcObjs.back().setRange(start, end, i);
            if (end >= nCalcs) break;
        }
#if USE_BOOST_THREAD
        boost::thread_group threads;
        for (uint i = 0; i < calcObjs.size(); i++) {
            threads.create_thread(calcObjs[i]);
        }
        threads.join_all();
#else

        std::mutex iomutex;
        std::vector<std::thread> threads(calcObjs.size());

        for (uint i = 0; i < calcObjs.size(); i++) {
            //threads.emplace_back(calcObjs[i]);
            threads[i] = std::thread( [&iomutex, i, &calcObjs] {
                // Stopwatch swatch(true);
                {
                    std::lock_guard<std::mutex> iolock(iomutex);
                    //log(Debug, "Thread #" + str(i) + ": on CPU "
                    //+ str(schedGetCPU()) + " slice " + str(calcObjs[i].start()) + ":" + str(calcObjs[i].end()));
                }
                calcObjs[i]();
                {
                    std::lock_guard<std::mutex> iolock(iomutex);
                    // log(Debug, "time: #" + str(i) + " " + str(swatch.duration()) + "s");
                }

            });
        }

        for (auto & t: threads) if (t.joinable()) t.join();

#endif
    }
}

} // namespace GIMLI{

#endif //_GIMLI_IPC_CLIENT__H
