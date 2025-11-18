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

#include "gimli.h"

#include <cerrno>
#include <cstring>
#include <fstream>
#include <cmath>
#include <cstdio>
#include <cstdlib>

//#include <omp.h> need -lgomp of -fopenmp

#if OPENBLAS_CBLAS_FOUND
    #include <cblas.h>
#endif

#if USE_BOOST_THREAD
    #include <boost/thread.hpp>
    boost::mutex __GIMLILogWriteMutex__;
#else
    #include <mutex>
    std::mutex __GIMLILogWriteMutex__;
#endif

#include <omp.h>

namespace GIMLI{

static bool __SAVE_PYTHON_GIL__ = false;
static bool __GIMLI_DEBUG__ = false;
static int __GIMLI_DEEP_DEBUG__ = 0;

Index __setTC__(){
    // to close to number of cpu can bring significent shedular overhead
    int omp = getEnvironment("OMP_NUM_THREADS", -1, false);
    if (omp == -1){
        omp_set_num_threads(min(8, numberOfCPU()-2));
    }

    int tc = getEnvironment("OPENBLAS_NUM_THREADS", -1, false);

    if (tc == -1){
        int tc = getEnvironment("GIMLI_NUM_THREADS", -1, false);

    }

// __MS(tc)
    tc = Index(min(8, numberOfCPU()-2));
// __MS(Index(tc))
    setThreadCount(Index(tc));
    return Index(tc);
}

static Index __GIMLI_THREADCOUNT__ = __setTC__();

void setThreadCount(Index nThreads){

    nThreads = max(1, nThreads);
    log(Debug, "Set amount of threads to " + str(nThreads));
    //log(Debug, "omp_get_max_threads: " + str(omp_get_max_threads()));
    //omp_set_num_threads
#if OPENBLAS_CBLAS_FOUND
    openblas_set_num_threads(nThreads);
    //omp_set_num_threads
#else
    log(Debug, "can't set openblas thread count. ");
#endif

    __GIMLI_THREADCOUNT__ = nThreads;
    // __MS(__GIMLI_THREADCOUNT__)
}

Index threadCount(){
    return __GIMLI_THREADCOUNT__;
}


bool __setOMP__(){

    int tc = getEnvironment("PG_USE_OMP", -1, false);
    if (tc == -1) {
        // default is false
        tc = false;
    }
    setUseOMP(tc);
    return (bool)tc;
}

static bool __GIMLI_USE_OMP__ = __setOMP__();

void setUseOMP(bool o){
    __GIMLI_USE_OMP__ = o;
}
bool useOMP(){
    return __GIMLI_USE_OMP__;
}

// //** end forward declaration
// // static here gives every .cpp its own static bool
// extern static bool __SAVE_PYTHON_GIL__;
// extern static bool __GIMLI_DEBUG__;

std::string versionStr(){
    std::string vers(str(PACKAGE_NAME) + "-" + PACKAGE_VERSION);
    return vers;
}

void savePythonGIL(bool s){ __SAVE_PYTHON_GIL__ = s; }
bool pythonGIL(){ return __SAVE_PYTHON_GIL__; }

void setDebug(bool s){ __GIMLI_DEBUG__ = s; }
bool debug(){ return __GIMLI_DEBUG__;}

void setDeepDebug(int level){__GIMLI_DEEP_DEBUG__ = level;}
int deepDebug(){ return __GIMLI_DEEP_DEBUG__; }

void PythonGILSave::save() {
    if (!saved_) {
#ifdef PYGIMLI
//#warning  "save_ = PyEval_SaveThread();"
        if (__SAVE_PYTHON_GIL__) save_ =  PyEval_SaveThread();
#endif
        saved_ = true;
    }
}
void PythonGILSave::restore() { if (saved_) {
#ifdef PYGIMLI
    if (__SAVE_PYTHON_GIL__) PyEval_RestoreThread(save_);
#endif
        saved_ = false;
    }
}

void showSizes(){
    std::cout << "size_t: " << sizeof(size_t) << std::endl;
    std::cout << "ssize_t: " << sizeof(ssize_t) << std::endl;
    std::cout << "Index: " << sizeof(GIMLI::Index) << std::endl;
    std::cout << "Sindex: " << sizeof(GIMLI::SIndex) << std::endl;

    std::cout << "int: " << sizeof(int) << std::endl;
    std::cout << "long: " << sizeof(long) << std::endl;
    std::cout << "long long int: " << sizeof(long long int) << std::endl;

    std::cout << "int8: " << sizeof(GIMLI::int8) << std::endl;
    std::cout << "int16: " << sizeof(GIMLI::int16) << std::endl;
    std::cout << "int32: " << sizeof(GIMLI::int32) << std::endl;
    std::cout << "int64: " << sizeof(GIMLI::int64) << std::endl;

    std::cout << "uint8: " << sizeof(GIMLI::uint8) << std::endl;
    std::cout << "uint16: " << sizeof(GIMLI::uint16) << std::endl;
    std::cout << "uint32: " << sizeof(GIMLI::uint32) << std::endl;
    std::cout << "uint64: " << sizeof(GIMLI::uint64) << std::endl;

    std::cout << "float: " << sizeof(float) << std::endl;
    std::cout << "double: " << sizeof(double) << std::endl;
}

std::string authors(){
    return (std::string)("bugs and suggestions to:") + PACKAGE_AUTHORS;
}

int openFile(const std::string & fname, std::fstream * file,
             std::ios_base::openmode farg, bool terminate){

    file->open(fname.c_str(), farg);
    if (!*file){
        if (terminate) {
            throwError(WHERE_AM_I + " '" + fname + "': " +strerror(errno) + str(errno));
        } else {
			std::cerr << fname << ": " << strerror(errno) << " " << errno << std::endl;
		}
        return 0;
    }
    return 1;
}

bool fileExist(const std::string & filename){
    bool result = false;
    std::ifstream file; file.open(filename.c_str());
    if (file) {
        result = true;
        file.close();
    }
    return result;
}

uint fileLength(std::fstream & file){
    std::streamoff oldPos = file.tellg();
    file.seekg(0, std::ios::end);
    std::streamoff length = file.tellg();
    file.seekg(oldPos);
    return (uint)length;
}

uint countColumnsInFile(const std::string & fname){
    uint columnsCount = 0;
    return countColumnsInFile(fname, columnsCount);
}

uint countColumnsInFile(const std::string & fname, uint & columnsCount){
    columnsCount = 0;
    std::fstream file; if (!openInFile(fname, & file, false)) { return 0; }

    std::vector < std::string > subStrings;
    std::string str, tmp;

    while (!file.eof()){
        getline(file, str);
        if (str.find('#') != std::string::npos || str.empty()) {
            columnsCount++;
        } else {
            file.close();
            return getSubstrings(str).size();
        }
    }

    file.close();
    return 0;
}

uint countRowsInFile(const std::string & fname){
    std::fstream file; openInFile(fname, & file);
    CERR_TO_IMPL;
    file.close();
    return 0;
}

std::vector < std::string > getRowSubstrings(std::fstream & file, char comment){
  std::vector < std::string > subStrings;
  std::string str, tmp;
  getline(file, str);

  std::istringstream is(str.substr(0, str.find(comment)));
  while (is >> tmp) subStrings.push_back(tmp);
  return subStrings;
}

std::vector < std::string > getNonEmptyRow(std::fstream & file, char comment){
  std::vector < std::string > row;
  while((row = getRowSubstrings(file, comment)).empty() && !file.eof());
  return row;
}

std::vector < std::string > getCommentLine(std::fstream & file, char comment){
    std::vector< std::string > row;
    std::string str;
    getline(file, str);
    row = getSubstrings(str.substr(str.find(comment), std::string::npos));
    return row;
}

std::vector < std::string > getSubstrings(const std::string & str){
  std::vector < std::string > subStrings;
  std::istringstream is(str);
  std::string tmp;
  while (is >> tmp) subStrings.push_back(tmp);
  return subStrings;
}

std::vector < std::string > split(const std::string & str, char delimiter){
    std::vector < std::string > subStrings;
    size_t pos = 0;
    size_t lastPos = 0;
    //std::cout << str << std::endl;
    while ((pos = str.find(delimiter, lastPos)) != std::string::npos){
        //std::cout << pos << " " << lastPos << " " << str.substr(lastPos, pos -lastPos) << std::endl;
        subStrings.push_back(str.substr(lastPos, pos - lastPos));
        lastPos = pos + 1;
    }
    //std::cout << pos << " " << lastPos << " " << str.substr(lastPos, pos -lastPos) << std::endl;
    subStrings.push_back(str.substr(lastPos));
    return subStrings;
}

std::string replace(const std::string & str, const std::string & from,
                    const std::string & to){
    std::string ret(str);
    if (ret.find(from) != std::string::npos){
        ret.replace(ret.find(from), from.length(), to);
    }
    return ret;
}
std::string replace(const std::string & str, char from, char to){
    std::string ret(str);
    std::replace(ret.begin(), ret.end(), from, to);
    return ret;
}
std::string lower(const std::string & str){
    std::string lo(str);
    std::transform(lo.begin(), lo.end(), lo.begin(), tolower);
    return lo;
}


std::map < float, float > loadFloatMap(const std::string & filename){
    std::map < float, float > aMap;

    std::fstream file; if (!openInFile(filename, & file)){
        std::cerr << WHERE_AM_I << " Map not found: " << filename << std::endl;
        return aMap;
    }

    std::vector < std::string > row;
    while (!file.eof()) {
        row = getNonEmptyRow(file);
        if (row.size() == 2){
            aMap[toFloat(row[0])] = toFloat(row[1]);
        } else {
            if (aMap.size() == 0){
                throwError("no proper format found for map <float, Complex> in " + filename  + " " + str(row.size()) );
            }
        }
    }

    file.close();
    return aMap;
}

std::map < float, Complex > loadCFloatMap(const std::string & filename){
    std::map < float, Complex > aMap;

    std::fstream file; if (!openInFile(filename, & file)){
        std::cerr << WHERE_AM_I << " Map not found: " << filename << std::endl;
        return aMap;
    }

    std::vector < std::string > row;
    while (!file.eof()) {
        row = getNonEmptyRow(file);
        if (row.size() == 3){
            aMap[toFloat(row[0])] = Complex(toDouble(row[1]),
                                            toDouble(row[2]));
        } else {
            if (aMap.size() == 0){
                throwError("no proper format found for map <float, Complex> in " + filename  + " " + str(row.size()) );
            }
        }
    }

    file.close();
    return aMap;
}

std::map < int, int > loadIntMap(const std::string & filename){
    std::map < int, int > aMap;

    std::fstream file; if (!openInFile(filename, & file)){
        std::cerr << WHERE_AM_I << " Map not found: " << filename << std::endl;
        return aMap;
    }

    std::vector < std::string > row;
    while (! file.eof()) {
        row = getNonEmptyRow(file);
        if (row.size() == 2) aMap[toInt(row[0])] = toInt(row[1]);
    }

    file.close();
    return aMap;
}

std::string logStrShort_(LogType type){
    switch (type){
        case Verbose: return "Verbose"; break;
        case Info: return "info"; break;
        case Warning: return "warning";  break;
        case Error: return "error";  break;
        case Debug: return "debug";  break;
        case Critical: return "critical";  break;
        default: return str(type) + "-unknown";
    }
}

std::string logStr_(LogType type){
    switch (type){
        case Verbose: return "Verbose"; break;
        case Info: return "Info"; break;
        case Warning: return "Warning";  break;
        case Debug: return "Debug";  break;
        case Error: return "Error!";  break;
        case Critical: return "!!! Critical !!!";  break;
        default: return str(type) + "-unknown";
    }
}

#ifdef PYTHON_FOUND
    #define HAVE_ROUND 1
    #include <Python.h>
#endif


void log(LogType type, const std::string & msg){
#if USE_BOOST_THREAD
    #ifndef PYGIMLI_CAST // fails because of boost threads and clang problems
        boost::mutex::scoped_lock lock(__GIMLILogWriteMutex__);
    #endif
#else
    std::unique_lock < std::mutex > lock(__GIMLILogWriteMutex__);
#endif

#if defined(PYTHON_FOUND) && not defined(WIN32)
    //#if defined(RUN_FROM_PY) && defined(PYTHON_FOUND) && not defined(WIN32)
    static PyObject *logging = NULL;
    static PyObject *logger = NULL;
    static PyObject *string = NULL;

    if (Py_IsInitialized()){ // only use python logger if called from python runtime

        logging = PyImport_ImportModule("logging");
        //logging = PyImport_ImportModuleNoBlock("logging");

        if (logging != NULL){
            logger = PyObject_CallMethod(logging, "getLogger", "s", "Core");
            string = Py_BuildValue("s", msg.c_str());
            PyObject_CallMethod(logger, logStrShort_(type).c_str(), "O", string);
            Py_DECREF(string);
            return;
        }
    }
#endif
    if (type == Debug && !debug()) return;
    if (type == Critical) throwError(logStr_(type) + ": " + msg);
    std::cout << logStr_(type) << ": " << msg << std::endl;
}

} // namespace GIMLI

