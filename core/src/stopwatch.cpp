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

#include "stopwatch.h"
#include "vector.h"

#include <iostream>

namespace GIMLI{

Stopwatch::Stopwatch(bool start) : _store(nullptr) {
    //newPtr(_store);
    _state = undefined;
    if (start) this->start();
    this->_store = new RVector();
}

Stopwatch::~Stopwatch() {
    // deletePtr(this->_store);
    delete this->_store;
}

void Stopwatch::start(){
    this->_start = std::chrono::high_resolution_clock::now();
    _state = running;
    _cCounter.tic();
}

void Stopwatch::stop(bool verbose){
    this->_stop = std::chrono::high_resolution_clock::now();
    _state = halted;
    if (verbose) std::cout << "time: " << duration() << "s" << std::endl;
}

void Stopwatch::restart(){
    stop();
    start();
}

void Stopwatch::reset(){
    restart();
    this->_store->clear();
}

double Stopwatch::duration(bool restart){
    std::chrono::time_point<std::chrono::high_resolution_clock> now;
    
    if (_state == undefined) {
        log(Error, "Stopwatch not started!");
    }
        
    if (_state == running) {
        //ftime(&stoptime);
        // double t = (stoptime.time - starttime.time) + 
        //     double(stoptime.millitm - starttime.millitm) / 1000.0;
            
        now = std::chrono::high_resolution_clock::now();
    } else {
        now = this->_stop;
    }
    
    std::chrono::duration< double > t = now - this->_start;

    if (restart) this->restart();
    return t.count();
}

void Stopwatch::store(bool restart){
    this->_store->push_back(this->duration(restart));
}

size_t Stopwatch::cycles(bool res){
    if (_state == undefined) {
        log(Error, "Stopwatch not started!");
    }
    size_t t = 0;
    if (_state == running) t = _cCounter.toc();
    if (res) restart();
    return t;
}

} // namespace GIMLI