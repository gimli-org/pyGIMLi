/******************************************************************************
 *   Copyright (C) 2011-2024 by the resistivity.net development team          *
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

#include "bertDataContainer.h"

#include <map>

namespace GIMLI{

DataContainerERT::DataContainerERT()
    : DataContainer(){
    init();
}

DataContainerERT::DataContainerERT(const std::string & fileName, bool removeInvalid)
    : DataContainer(){
    init();
    load(fileName, true, removeInvalid);
}

DataContainerERT::DataContainerERT(const DataContainerERT & data)
    : DataContainer(){
    init();
    this->copy_(data);
}

DataContainerERT::~DataContainerERT(){
}

void DataContainerERT::init(){
    fillCounter_ = 0;
    this->registerSensorIndex("a");
    this->registerSensorIndex("b");
    this->registerSensorIndex("m");
    this->registerSensorIndex("n");

    dataMap_["i"] = RVector(0);
    dataMap_["u"] = RVector(0);
    dataMap_["r"] = RVector(0);
    dataMap_["err"] = RVector(0);
    dataMap_["k"] = RVector(0);
    dataMap_["rhoa"] = RVector(0);
    dataMap_["ip"] = RVector(0);
    dataMap_["iperr"] = RVector(0);
    sensorIndexOnFileFromOne_ = true;
    initTokenTranslator();
}

void DataContainerERT::initTokenTranslator(){
    DataContainer::initTokenTranslator();
    std::map< std::string, std::string > l;

    l["a"] = "a c1"; //** s and g not needed in ERT data container
    l["b"] = "b c2";
    l["m"] = "m p1";
    l["n"] = "n p2";
    l["rhoa"] = "rhoa rho_a ra rs rhos rhoa/Ohmm rhoa(Ohmm)"; //** apparent resistivity
    l["eca"] = "ECa EC_a";
    l["r"] = "r rho r(Ohm) imp z u/i"; //** rho is confusing!
// done in base l["err"] = "err std error err/%";
    l["ip"] = "ip ip/mrad ip/° phase phase/mrad phase/° phi phi/mrad phi/°";
    l["iperr"] = "iperr iperr/mrad iperr/° phierr phierr/mrad phierr/°";
    l["u"] = "u u/V u(V) u/mV u(mV) v v/V v(V) v/mV v(mV)";
    l["i"] = "i i/A i(A) i/mA i(mA)";
    l["k"] = "k";
// done in base     l["t"] = "t t/s t/ms t/(ms)";

    for (std::map< std::string, std::string >::iterator it = l.begin(); it != l.end(); it ++){
        std::vector< std::string> row(getSubstrings(it->second));
        for (Index i = 0; i < row.size(); i ++){
            tT_.insert(std::pair< std::string, std::string >(row[i], it->first));
        }
    }
}

void DataContainerERT::checkDataValidityLocal(){

    if (size() > 0){
        if (this->haveData("rhoa") || this->haveData("r") || this->haveData("u")){
            //** no shm only here

            if (!allNonZero("rhoa")) { //** there are zero rhoa

                //** we are looking for impedance
                if (!allNonZero("r")) { //** there are zero r

                    //** can we build impedance, myself
                    if (allNonZero("i") && allNonZero("u")) {
                        __DS("setting r = u/i")
                        set("r", get("u") / get("i"));
                    } else {
                        std::cerr << WHERE_AM_I << " In datafile zero current or voltage found." << std::endl;
                    }
                }
            }
        }

        if (this->haveData("rhoa")){
            IndexArray idx(find((*this)("rhoa") < TOLERANCE));

            if (debug()){
                if (idx.size()){
                    __DS("rhoa < 0")
                    for (Index i = 0; i < idx.size(); i++){
                        __DS(str(idx[i]) + " " + str((*this)("rhoa")[idx[i]]))
                    }
                }
            }
            this->markInvalid(idx);
        }

        if (this->haveData("r")){
            IndexArray idx(find(abs((*this)("r")) < TOLERANCE));
            if (debug()){
                if (idx.size()){
                    __DS("|r| == 0")
                    for (Index i = 0; i < idx.size(); i++){
                        __DS(str(idx[i]) + " " + str((*this)("r")[idx[i]]))
                    }
                }
            }
            this->markInvalid(idx);
        }

        this->markInvalid(find(((abs((*this)("a") - (*this)("b")) < TOLERANCE) & ((*this)("a") > 0)) |
                                ((abs((*this)("a") - (*this)("m")) < TOLERANCE) & ((*this)("a") > 0)) |
                                ((abs((*this)("a") - (*this)("n")) < TOLERANCE) & ((*this)("a") > 0)) |
                                ((abs((*this)("b") - (*this)("m")) < TOLERANCE) & ((*this)("b") > 0)) |
                                ((abs((*this)("b") - (*this)("n")) < TOLERANCE) & ((*this)("b") > 0)) |
                                ((abs((*this)("m") - (*this)("n")) < TOLERANCE) & ((*this)("m") > 0))));
    }
}

Index DataContainerERT::electrodeToCurrentPattern(Index a, Index b) const {
    //std::cout << a << " " << b << " " << electrodes_.size() << std::endl;
    /*! +1 ensures that -1 (pol) electrodes are considered. */
    return (a + 1) * (sensorCount() + 1) + b + 1; //** better take min(a,b) and max(a,b)?
}

CurrentPattern DataContainerERT::currentPatternToElectrode(Index pattern){
    /*! +1 ensures that -1 (pol) electrodes are considered. */

    Index nElecs = sensorCount() + 1;
    Index a = (SIndex)ceil(pattern / nElecs);
    Index b = pattern - a * nElecs;
    a-=1; b-=1;
    //std::cout << pattern << " " << a << " " << b << std::endl;
    return CurrentPattern(a, b);
}

std::set < Index > DataContainerERT::currentPattern(bool reciprocity){
    std::set < Index > pattern;

    const RVector & valid = dataMap_["valid"];
    const RVector & a = dataMap_["a"];
    const RVector & b = dataMap_["b"];
    const RVector & m = dataMap_["m"];
    const RVector & n = dataMap_["n"];

    for (Index i = 0; i < this->size(); i ++){
        if (valid[i]){
            pattern.insert(electrodeToCurrentPattern(a[i], b[i]));

            if (reciprocity){
                pattern.insert(electrodeToCurrentPattern(m[i], n[i]));
            }
        }
    }
    return pattern;
}

Index DataContainerERT::addFourPointData(SIndex a, SIndex b, SIndex m, SIndex n){
    return this->createFourPointData(this->size(), a, b, m, n);
}

Index DataContainerERT::createFourPointData(Index i, SIndex a, SIndex b, SIndex m, SIndex n){
    if (this->size() <= i) {
        resize(max(i+1,1));
        // memory reservation is vectors job so this should be ok
        //resize(max(i,1) * 2);
    }
    fillCounter_ = max(fillCounter_, i+1);

    if ((a == b && a != -1) ||
         (a == m && a != -1) ||
         (a == n && a != -1) ||
         (b == m && b != -1) ||
         (b == n && b != -1) ||
         (m == n && m != -1)) {
            std::cerr << "WARNING! " << WHERE_AM_I << " Error:  Electrode config wrong data: "   << i
                << " ; a = " << a << " b = " << b << " m = " << m << " n = " << n << std::endl;
//             warn(EXIT_DATACONTAINER_ELECS, str.str());
        dataMap_["valid"][i] = 0;
    } else if (this->size() > i){
        dataMap_["a"][i] = a;
        dataMap_["b"][i] = b;
        dataMap_["m"][i] = m;
        dataMap_["n"][i] = n;
        dataMap_["valid"][i] = 1;
    } else {
        throwError(WHERE_AM_I + " index out of size, resize data first:" + str(this->size()) + " " + str(i));
    }
    return i;
}

void DataContainerERT::fitFillSize(){
    if (fillCounter_ > 0 && this->size() != fillCounter_) resize(fillCounter_);
}

void DataContainerERT::averageDuplicateData(bool verbose){
    DataContainerERT origData(*this);

    /*! +1 ensures that -1 (pol) electrodes are considered. */
    Index N = this->sensorCount() + 1;
    Index N2 = N * N;
    Index N3 = N * N * N;

    std::map< SIndex, IndexArray > allDataIdxMap;

    std::map< SIndex, IndexArray >::iterator it;

    for (Index i = 0; i < this->size(); i ++){
        Index a = dataMap_["a"][i] + 1;
        Index b = dataMap_["b"][i] + 1;
        Index m = dataMap_["m"][i] + 1;
        Index n = dataMap_["n"][i] + 1;

//        SIndex dataIdx = a * N3 + b * N2 + m * N + n;
        Index dataIdx = min(a, b) * N3 + max(a, b) * N2 + min(m, n) * N + max(m, n);

        it = allDataIdxMap.find(dataIdx);
        if (it == allDataIdxMap.end()){
            it = allDataIdxMap.insert(std::pair< SIndex, IndexArray >(dataIdx, IndexArray())).first;
        }
        it->second.push_back(i);
    }

    this->resize(allDataIdxMap.size());

    Index i = 0;
    for (it = allDataIdxMap.begin(); it != allDataIdxMap.end(); it ++, i++){

        SIndex dataIdx = it->first;

        SIndex a = (SIndex)ceil(dataIdx / N3);
        SIndex b = (SIndex)ceil((dataIdx - (a * N3)) / N2);
        SIndex m = (SIndex)ceil((dataIdx - (a * N3 + b * N2)) / N);
        SIndex n = (SIndex)ceil((dataIdx - (a * N3 + b * N2 + m * N)));

        dataMap_["a"][i] = a - 1;
        dataMap_["b"][i] = b - 1;
        dataMap_["m"][i] = m - 1;
        dataMap_["n"][i] = n - 1;

        for (std::map< std::string, RVector >::iterator itD = dataMap_.begin(); itD != dataMap_.end(); itD ++){

            if (!isSensorIndex(itD->first)){
                //** merge data fields
                itD->second.setVal(mean(origData.get(itD->first)(it->second)), i);
            }
        }
    }

    if (verbose){
        std::cout << "Merged " << origData.size() - this->size() << " data points." << std::endl;
    }

}


} // namespace BERT
