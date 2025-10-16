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

#ifndef _GIMLI_STOPWATCH__H
#define _GIMLI_STOPWATCH__H

#include "gimli.h"

#include <sys/timeb.h>
#include <chrono>

#if defined(__i386__)
static __inline__ size_t rdtsc__(void){
    size_t x;
    __asm__ volatile (".byte 0x0f, 0x31" : "=A" (x));
    return x;
}
#elif defined(__x86_64__)
static __inline__ size_t rdtsc__(void){
    unsigned hi, lo;
    __asm__ __volatile__ ("rdtsc" : "=a"(lo), "=d"(hi));
    return ((size_t)lo)|(((size_t)hi)<<32);
}
#else
static inline size_t rdtsc__(void){
    return 0;
}
#endif

namespace GIMLI{

class DLLEXPORT CycleCounter{
public:
    CycleCounter() : var_(0) {}

    ~CycleCounter(){}

    inline void tic(){ var_ = rdtsc__(); }

    inline size_t toc() const { return (rdtsc__() - var_); }

protected:

    size_t var_;
};

class DLLEXPORT Stopwatch {
public:
    Stopwatch(bool start=true);

    ~Stopwatch();

    void start();

    void stop(bool verbose=false);

    /*! Restart the stopwatch.*/
    void restart();

    /*! Reset the stopwatch, and delete stored values.*/
    void reset();

    /*! Returns the current duration in seconds. Optional you can restart the stopwatch.*/
    double duration(bool restart=false);

    /*! Returns the cpu cycles. Optional you can restart the stopwatch.*/
    size_t cycles(bool restart=false);

    const CycleCounter & cycleCounter() const { return _cCounter; }

    /*!Save current duration value into store.
    Reset current time if needed to store relative times.*/
    void store(bool restart=false);

    /*!Get stored duration values.*/
    const RVector & stored() const { return *this->_store;}

protected:
    enum watchstate {undefined,halted,running} _state;

    std::chrono::time_point<std::chrono::high_resolution_clock> _start, _stop;
    RVector *_store;

    CycleCounter _cCounter;
};

#define TIC__ std::cout.precision(12); GIMLI::Stopwatch __swatch__(true);
#define TOC__ std::cout << __swatch__.duration(true) << std::endl;
#define toc__ __swatch__.duration()

// little debugging and testing fun
DLLEXPORT void waitms(Index us, Index count=1);
DLLEXPORT void waitmsOMP(Index us, Index count=1);

/*! Microseconds are rounded somewhat to system clock of 1000Hz.
Also there seems to be an offset of 0.05 ms. */
DLLEXPORT void waitus(Index us, Index count=1);
DLLEXPORT void waitusOMP(Index us, Index count=1);


class DLLEXPORT Swatches: public Singleton< Swatches >{
public:
    friend class Singleton< Swatches >;

    Stopwatch & operator[](const std::string & key);

    std::vector < std::string > keys();

    std::vector < const Stopwatch * > vals();

    std::vector < std::pair< std::string, Stopwatch * > > items();

    void remove(const std::string & key, bool isRoot=false);

    void setTrace(const std::string & trace) { _trace = trace; }

    const std::string & trace() const { return _trace; }

protected:
    std::map < std::string, Stopwatch * > _sw;

private:
    /*! Private so it can not be called */
    Swatches();
    /*! Private so it can not be called */
    virtual ~Swatches();
    /*! Copy constructor is private, so don't use it */
    Swatches(const Swatches &){};
    /*! Assignment operator is private, so don't use it */
    void operator = (const Swatches &){ };

    std::string _trace;
};

class DLLEXPORT TicToc{
public:
    TicToc(const std::string & name, bool reset=false);

    ~TicToc();

protected:

    Stopwatch *_sw;
    std::string _parentTrace;

};

#define WITH_TICTOC(name) TicToc tictoc_name(name);


class DLLEXPORT PickleTest{
public:
    PickleTest();

    ~PickleTest();

    std::string serialize();

    void deserialize(const std::string & s);

    void setName(const std::string &n );
    std::string name() const { return this->_name; }

protected:
    std::string _name;
};

} // namespace GIMLI

#endif // _GIMLI_STOPWATCH__H