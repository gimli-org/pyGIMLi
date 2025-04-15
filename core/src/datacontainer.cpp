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

#include <list>

#include "datacontainer.h"
#include "pos.h"
#include "numericbase.h"
#include "vectortemplates.h"


namespace GIMLI{

DataContainer::DataContainer(){
    initDefaults();
    //std::cout << "DataContainer(){" << std::endl;
}

DataContainer::DataContainer(const std::string & fileName,
                             bool sensorIndicesFromOne,
                             bool removeInvalid){
    initDefaults();
    this->load(fileName, sensorIndicesFromOne, removeInvalid);
}

DataContainer::DataContainer(const std::string & fileName,
                             const std::string & sensorTokens,
                             bool sensorIndicesFromOne,
                             bool removeInvalid){
    initDefaults();

    std::vector < std::string > tokenList = getSubstrings(sensorTokens);
    for (Index i=0; i < tokenList.size(); i++) registerSensorIndex(tokenList[i]);
    this->load(fileName, sensorIndicesFromOne, removeInvalid);
    //std::cout << "DataContainer(const std::string & fileName){" << std::endl;
}

DataContainer::DataContainer(const DataContainer & data){
    initDefaults();
    this->copy_(data);
    //std::cout << "DataContainer(const DataContainer & data){" << std::endl;
}

DataContainer::~DataContainer(){
    initDefaults();
    clear();
}

DataContainer & DataContainer::operator = (const DataContainer & data){
    if (this != & data){
        this->copy_(data);
    }
    //std::cout << "DataContainer = " << std::endl;
    return * this;
}

void DataContainer::initDefaults(){
    dataMap_["valid"] = RVector(0);
    sensorIndexOnFileFromOne_ = true;
    init();
    initTokenTranslator();
}

void DataContainer::init(){
}

void DataContainer::initTokenTranslator(){
    std::map< std::string, std::string > l;

    l["err"] = "err std error err/%";
    l["t"] = "t t/s t/ms t/(ms)";

    for (std::map< std::string, std::string >::iterator it = l.begin(); it != l.end(); it ++){
        std::vector< std::string> row(getSubstrings(it->second));
        for (Index i = 0; i < row.size(); i ++){
            tT_.insert(std::pair< std::string, std::string >(row[i], it->first));
        }
    }
}

bool DataContainer::haveTranslationForAlias(const std::string & alias) const {
    std::string test(alias);
    test[0] = tolower(test[0]);
    return tT_.count(test) != 0;
}

std::string DataContainer::translateAlias(const std::string & alias) const {
    std::string test(alias);
    test[0] = tolower(test[0]);
    if (tT_.find(test) != tT_.end())
        return tT_.find(test)->second;
    return "no alias";
}

void DataContainer::clear(){
    topoPoints_.clear();
    sensorPoints_.clear();
    // sensor token list should not be cleared
    // dataSensorIdx_.clear();
    dataMap_.clear();
    initDefaults();
}

void DataContainer::copy_(const DataContainer & data){
    clear();

    topoPoints_ = data.additionalPoints();
    sensorPoints_= data.sensorPositions();

    this->resize(data.size());

    inputFormatString_ = data.inputFormatString();
    inputFormatStringSensors_ = data.formatStringSensors();
    dataSensorIdx_  = data.sensorIdx();
    dataMap_ = data.dataMap();

    dataDescription_ = data.dataDescription();

    tT_ = data.tokenTranslator();

    sensorIndexOnFileFromOne_ = data.sensorIndexOnFileFromOne();
}

void DataContainer::add(const DataContainer & data, double snap){

    Index start = this->size();
    this->resize(this->size() + data.size());

    IndexArray perm(data.sensorCount(), 0);

    //** merge sensor data
    for (Index i = 0; i < data.sensorCount(); i ++){
        perm[i] = this->createSensor(data.sensorPositions()[i], snap);
    }

    for (std::map< std::string, RVector >::iterator it = dataMap_.begin(); it!= dataMap_.end(); it ++){
        //** permutate sensorindices
        if (isSensorIndex(it->first)){
            RVector tmp(data.get(it->first));
            for (Index i = 0; i < tmp.size(); i ++){
                SIndex id = (SIndex)tmp[i];
                if (id > -1 && id < (SIndex)perm.size()) {
                    it->second[start + i] = perm[id];
                } else {
                    it->second[start + i] = -1;
                }
            }
        } else {
        //** merge data fields
            it->second.setVal(data.get(it->first), start, -1);
        }
    }
}

SIndex DataContainer::createSensor(const RVector3 & pos, double tolerance){

    SIndex ret = -1;
    for (Index i = 0; i < sensorPoints_.size(); i ++){
        if (pos.distance(sensorPoints_[i]) < tolerance){
            ret = i;
        }
    }

    if (ret == -1){
        ret = sensorPoints_.size();
        sensorPoints_.push_back(pos);
    }
    return ret;
}

void DataContainer::registerSensorIndex(const std::string & token) {
    dataSensorIdx_.insert(token);
    if (!this->exists(token)){
        this->set(token, RVector(this->size(), -1.0));
    }
}

bool DataContainer::isSensorIndex(const std::string & token) const {
    return dataSensorIdx_.find(token) != dataSensorIdx_.end();
}

void DataContainer::setSensorPositions(const RVector & sensors){
    std::vector< RVector3 > s;
    for (Index i = 0; i < sensors.size(); i ++ ) {
        s.push_back(RVector3(sensors[i], 0.0, 0.0));
    }
    return this->setSensorPositions(s);
}

void DataContainer::setSensorPositions(const PosVector & sensors) {
    sensorPoints_ = sensors;
}

int DataContainer::load(const std::string & fileName,
                        bool sensorIndicesFromOne,
                        bool removeInvalid){

    // __MS(std::locale::locale(LC_NUMERIC));
    // __MS(std::localeconv()->decimal_point[0]));
    //putenv("LC_NUMERIC=C");
    setlocale(LC_NUMERIC, "C"); // no std::setlocale on appleclang
    //std::setlocale(LC_NUMERIC, "C");
    // __MS(std::locale(LC_NUMERIC, "C"));
    // __MS(std::localeconv()->decimal_point[0]));

    clear();
    setSensorIndexOnFileFromOne(sensorIndicesFromOne);

	std::fstream file; openInFile(fileName, & file, true);

    std::vector < std::string > row(getNonEmptyRow(file));

    if (row.size() != 1){
        throwError(WHERE_AM_I + " cannot determine data format. " + str(row.size()));
    }

    //** read number of electrodes
    int nSensors = toInt(row[0]);
    if (nSensors < 1){
        throwError(" cannot determine sensor count " + row[0]);
    }
    RVector x(nSensors, 0.0), y(nSensors, 0.0), z(nSensors, 0.0);

    //** read electrodes format
    //** if no electrodes format is given (no x after comment symbol) take defaults
    std::string sensorFormatDefault("x y z");
    std::vector < std::string > format(getSubstrings(sensorFormatDefault));

    char c; file.get(c);
    if (c == '#') {
        format = getRowSubstrings(file);
//         if (format[0][0] != 'x' && format[0][0] != 'X'){
//             format = getSubstrings(elecsFormatDefault);
//         }
    } else {
        format = getSubstrings(sensorFormatDefault);
    }
    file.unget();

    inputFormatStringSensors_.clear();
    for (Index i = 0; i < format.size(); i ++)
        inputFormatStringSensors_ += format[i] + " ";

    //** read sensor
    for (int i = 0; i < nSensors; i ++){
        row = getNonEmptyRow(file);

        if (row.empty()){
            throwError(
                       WHERE_AM_I + "To few sensor data. " +
                       str(nSensors) + " Sensors expected but " +
                       str(i) + " found.");
        }

        for (Index j = 0; j < row.size(); j ++){

            if (j == format.size()) break; // no or to few format defined, ignore

            if (format[j] == "x" || format[j] == "X" ||  format[j] == "x/m" || format[j] == "X/m") x[i] = toDouble(row[j]);
            else if (format[j] == "x/mm" || format[j] == "X/mm") x[i] = toDouble(row[j]) / 1000.0;
            else if (format[j] == "y" || format[j] == "Y" ||  format[j] == "y/m" || format[j] == "Y/m" ) y[i] = toDouble(row[j]);
            else if (format[j] == "y/mm" || format[j] == "Y/mm") y[i] = toDouble(row[j]) / 1000.0;
            else if (format[j] == "z" || format[j] == "Z" ||  format[j] == "z/m" || format[j] == "Z/m" ) z[i] = toDouble(row[j]);
            else if (format[j] == "z/mm" || format[j] == "z/mm") z[i] = toDouble(row[j]) / 1000.0;
            else {
                std::cerr << WHERE_AM_I << " Warning! format description unknown: format[" << j << "] = " << format[j] << " column ignored." << std::endl;
            }
        }
    }

    std::list < std::pair< Index, Index > > duplicatedSensors;
    for (SIndex i = 0; i < nSensors; i ++) {
        Index oldCount = this->sensorCount();
        Index s = createSensor(RVector3(x[i], y[i], z[i]).round(1e-12));

        if (this->sensorCount() != oldCount+1){
            duplicatedSensors.insert(duplicatedSensors.end(),
                                     std::pair< Index, Index >(this->sensorCount(), s));
            log(Warning, "Duplicated sensor position found at: " + str(RVector3(x[i], y[i], z[i])));
            //** we add them temporary and delete them later as unused
            sensorPoints_.push_back(RVector3(x[i], y[i], z[i]).round(1e-12));
        }
    }
    //****************************** Start read the data;
    row = getNonEmptyRow(file);
    if (row.size() != 1) {
        for (Index i = 0; i < row.size(); i ++){
            std::cerr << row[i] << " ";
        }
        std::cerr << std::endl;

        throwError(WHERE_AM_I + " cannot determine data size. " + str(row.size()));
    }

    int nData = toInt(row[0]);

    //bool fromOne    = true;
    //bool schemeOnly = false;

    if (nData > 0){
        this->resize(nData);

        //** looking for # symbol which start format description section
        file.get(c);
        if (c == '#') {
            format = getRowSubstrings(file);
            if (format.size() == 0){
                throwError(WHERE_AM_I + "Can not determine data format.");
            }
//            if (format.size() == 4) schemeOnly = true;
        }
        file.unget();

    }

    std::map< std::string, RVector > tmpMap;

    for (int data = 0; data < nData; data ++){
        row = getNonEmptyRow(file);

        if (row.empty()){
            throwError(
                       WHERE_AM_I + " To few data. " + str(nData) +
                       " data expected and " + str(data) + " data found.");
        }

        for (Index j = 0; j < row.size(); j ++){
            if (j == format.size()) break;

            if (!tmpMap.count(format[j])){
                tmpMap.insert(std::pair< std::string, RVector > (format[j], RVector(nData, 0.0))) ;
            }
// 	    std::cout << format[j] << " " << tmpMap[format[j]].size() << std::endl;
// 	    std::cout << row[j] << " " << toDouble(row[j]) << std::endl;
            tmpMap[format[j]][data] = toDouble(row[j]);
        }
    }

    //** renaming formats with the token translator (tT_):
    for (std::map< std::string, RVector >::iterator it = tmpMap.begin();
         it != tmpMap.end(); it++){
        if (haveTranslationForAlias(it->first)){

            double scale = 1.0;
            // should only used in specialization
            if (it->first.rfind("mV") != std::string::npos ||
                it->first.rfind("mA") != std::string::npos ||
                it->first.rfind("ms") != std::string::npos){
                scale = 1.0 / 1000.0;
            }
            if (it->first.rfind("%") != std::string::npos ){
                scale = 1.0 / 100.0;
            }

// 	    std::cout << min(it->second) <<  " " << max(it->second) << std::endl;

//std::cout << "Insert translated:-" << tT_[it->first] << "-" << it->second.size() << std::endl;
            dataMap_[translateAlias(it->first)] = it->second * scale;
        } else {
//std::cout << "Insert:-" <<  it->first << "-" << it->second.size() << std::endl;
            dataMap_[it->first] = it->second;
        }
    }

    inputFormatString_.clear();
    for (Index i = 0; i < format.size(); i ++) {
        if (haveTranslationForAlias(format[i])){
            inputFormatString_ += translateAlias(format[i]) + " ";
        } else {
            inputFormatString_ += format[i] + " ";
        }
    }

    if (sensorIndexOnFileFromOne_){
        for (std::map< std::string, RVector >::iterator it = dataMap_.begin();
             it!= dataMap_.end(); it ++){
            if (isSensorIndex(it->first) && size() > 0){
                it->second -= RVector(size(), 1.);
            }
        }
    }

    dataMap_["valid"] = 1;

    //** Fix duplicated sensors
    if (!duplicatedSensors.empty()){
        for (auto p: duplicatedSensors){
            // print(p.first, p.second);

            for (auto &it: dataMap_){
                if (isSensorIndex(it.first) && size() > 0){
                    RVector &v = it.second;
                    for (Index i = 0; i < v.size(); i ++){
                        if (v[i] == p.first) v[i] = p.second;
                    }
                }
            }
        }
        this->removeUnusedSensors();
    }


    // validity check should only used in specialization
    this->checkDataValidity(removeInvalid);

    //** start read topography;
    row = getNonEmptyRow(file);

    if (row.size() == 1) {
        //** we found topography
        Index nTopoPoints = toInt(row[0]);

        if (nTopoPoints > 0){

            RVector xt(nTopoPoints), yt(nTopoPoints), zt(nTopoPoints);

            std::string topoFormatDefault("x y z");

            file.get(c);
            if (c == '#') {
                format = getRowSubstrings(file);
                if (format[0] != "x" && format[0] != "X"){
                    //** if no electrodes format is given (no x after comment symbol) take defaults;
                    format = getSubstrings(topoFormatDefault);
                }
            } else {
                file.unget();
            }

            //** read topography points;
            for (Index i = 0; i < nTopoPoints; i ++){
                row = getNonEmptyRow(file);

                if (row.empty()) {
                    throwError(WHERE_AM_I
                            + "To few topo data. " + str(nTopoPoints)
                            + " Topopoints expected and " + str(i) + " found.");
                }

                for (Index j = 0; j < row.size(); j ++){
                    if (j == format.size()) break; // no or to few format defined, ignore
        //__MS(j << " "<< format[j]<< " " << row[j] )
                    if (     format[j] == "x" || format[j] == "X") xt[i] = toDouble(row[j]);
                    else if (format[j] == "y" || format[j] == "Y") yt[i] = toDouble(row[j]);
                    else if (format[j] == "z" || format[j] == "Z") zt[i] = toDouble(row[j]);
                    else {
                        std::stringstream str;
                        str << " Warning! format description unknown: topo electrode format["
                            << j << "] = " << format[j] << " column ignored." << std::endl;
                        throwError(str.str());
                    }
                }
            }

            for (Index i = 0 ; i < nTopoPoints; i ++) {
                topoPoints_.push_back(RVector3(xt[i], yt[i], zt[i]).round(1e-12));
            }
        } // if nTopo > 0
    } // if topo

    file.close();
    return 1;
}

void DataContainer::checkDataValidity(bool remove){
    //** mark inf and nans as invalid
    int nInvalid0 = find(get("valid") < 1).size();
    for (std::map< std::string, RVector >::iterator it = dataMap_.begin();
            it!= dataMap_.end(); it ++){
        this->markInvalid(find(isInfNaN(it->second)));
    }
    int nInvalidNaN = find(get("valid") < 1).size() - nInvalid0;
    if (nInvalidNaN > 0) {
        std::cout << "Warning: removed " << nInvalidNaN << " values due to NaN/Inf!" << std::endl;
    }

    //** check sensor indices < -1 and >= sensorCount()
    for (std::map< std::string, RVector >::iterator it = dataMap_.begin(); it!= dataMap_.end(); it ++){
        if (isSensorIndex(it->first)){
//           std::cout << it->first << " "<< sensorCount()<< " " << find((it->second < -1) | (it->second >= sensorCount())) << std::endl;

            this->markInvalid(find((it->second < -1) | (it->second >= double(sensorCount()))));
        }
    }
    int nInvalidIdx = find(get("valid") < 1).size() - nInvalidNaN - nInvalid0;
    if (nInvalidIdx > 0) {
        std::cout << "Warning: removed " << nInvalidIdx << " values due to sensor index out of bounds!" << std::endl;
    }

    //** call local specialization if any
    checkDataValidityLocal();

    if (find(get("valid") < 1).size() > 0){
        std::cout << "Data validity check: found " << find(get("valid") < 1).size() << " invalid data. " << std::endl;
        GIMLI::save(RVector(find(get("valid") < 1)), "invalid.data");
        if (remove) {
            std::cout << "Data validity check: remove invalid data." << std::endl;
            this->removeInvalid();
        }
    }
}

int DataContainer::save(const std::string & fileName,
                        const std::string & formatData,
                        const std::string & formatSensor,
                        bool noFilter,
                        bool verbose) const {
    std::fstream file; if (!openOutFile(fileName, & file)) return 0;
    return write(file, formatData, formatSensor, noFilter, verbose);
}

int DataContainer::write(std::fstream & file,
                        const std::string & formatData,
                        const std::string & formatSensor,
                        bool noFilter,
                        bool verbose) const {

    file.precision(14);

    //** START write sensor data
    file << sensorPoints_.size() << std::endl;
    std::string sensorsString(formatSensor);
    file << "# " << sensorsString << std::endl;

    std::vector < std::string > format(getSubstrings(sensorsString));
    for (Index i = 0, imax = sensorPoints_.size(); i < imax; i ++){
        for (Index j = 0; j < format.size(); j ++){
            if (format[j] == "x" || format[j] == "X") file << sensorPoints_[i].x();
            if (format[j] == "y" || format[j] == "Y") file << sensorPoints_[i].y();
            if (format[j] == "z" || format[j] == "Z") file << sensorPoints_[i].z();
            if (j < format.size() -1) file << "\t";
        }
        file << std::endl;
    }

    //** START write data map
    std::vector < const RVector * > outVec;
    std::vector < bool > outInt;

    IndexArray toSaveIdx;
    std::string formatString;

    if (lower(formatData) == "all" ){
        toSaveIdx = find(get("valid") > -1);
        formatString = this->tokenList(false);
    } else if (noFilter){
        toSaveIdx = find(get("valid") > -1);
        formatString = formatData;
    } else {
        toSaveIdx = find(get("valid") == 1);
        formatString = formatData;
    }

    int count = toSaveIdx.size();
    file << count << std::endl;
    file << "# " << formatString << std::endl;

    std::vector < std::string > token(getSubstrings(formatString));

    for (Index i = 0; i < token.size(); i ++){
        //  std::cout << token[i] << std::endl;

        std::string valName;
        if (haveTranslationForAlias(token[i])){
            valName = translateAlias(token[i]);
        } else {
            valName = token[i];
        }

        if (dataMap_.count(valName)){
            outVec.push_back(&dataMap_.find(valName)->second );

            if (isSensorIndex(valName)){
                outInt.push_back(true);
            } else {
                outInt.push_back(false);
            }
        } else {
            throwError(WHERE_AM_I + " no such data: " + valName);
        }
    }

    for (Index i = 0; i < toSaveIdx.size(); i ++){
        file.setf(std::ios::scientific, std::ios::floatfield);
        file.precision(14);

        for (Index j = 0; j < token.size(); j ++){

            if (outInt[j]){
                file << int((*outVec[j])[toSaveIdx[i]]) + sensorIndexOnFileFromOne_;
            } else {
                 if (token[j] == "valid"){
                    file << int((*outVec[j])[toSaveIdx[i]]);
                 } else {
                    file << (*outVec[j])[toSaveIdx[i]];
                 }
            }
            if (j < token.size() -1) file << "\t";
        }
        file << std::endl;
    }

    //** START write additional points

    if (topoPoints_.size() > 0){
        std::cout << "Additional (topo) points are unhandled and will not be saved." << std::endl;
    }
    file << 0 << std::endl;

//     file << topoPoints_.size() << std::endl;
//     for (Index i = 0; i < topoPoints_.size(); i ++){
//         std::cout   << topoPoints_[i].x() << "\t"
//                     << topoPoints_[i].y() << "\t"
//                     << topoPoints_[i].z() << std::endl;
//     }

    file.close();
    return 1;
}

std::string DataContainer::tokenList(bool withAnnotation) const {
    std::string tokenList;
    if (withAnnotation) tokenList += "SensorIdx: ";

    for (std::map< std::string, RVector >::const_iterator it = dataMap_.begin(); it!= dataMap_.end(); it ++){
        if (isSensorIndex(it->first)){
            tokenList += it->first;
            tokenList += " ";
        }
    }
    if (withAnnotation) tokenList += " Data: ";
    for (std::map< std::string, RVector >::const_iterator it = dataMap_.begin(); it!= dataMap_.end(); it ++){
        if (!isSensorIndex(it->first)){
            tokenList += it->first;
            tokenList += " ";
        }
    }
    return tokenList;
}

void DataContainer::showInfos() const {
    std::cout << "Sensors: " << this->sensorCount() << ", Data: " << this->size();
    if (topoPoints_.size() > 0){
        std::cout << " Topopoints: " << topoPoints_.size();
    }
    std::cout << std::endl << tokenList() << std::endl;
}

Index DataContainer::hash() const{
    Index seed = GIMLI::hash(this->sensorPoints_, this->topoPoints_,
                             this->dataMap_);
    return seed;
}


void DataContainer::add(const std::string & token, const RVector & data,
                        const std::string & description){
    this->set(token, data);
    this->setDataDescription(token, description);
//     if (data.size() == this->size()) {
//         dataMap_.insert(make_pair(token, data));
//         this->setDataDescription(token, description);
//     } else {
//         throwError(WHERE_AM_I + " wrong data size: " + str(this->size()) + " " + str(data.size()));
//     }
}

void DataContainer::set(const std::string & token, const RVector & data){
    if (data.size() == this->size()){
        dataMap_[token] = data;
    } else {
        throwError(WHERE_AM_I + " wrong data size: " + str(this->size()) + " " + str(data.size()));
    }
}

const RVector & DataContainer::get(const std::string & token) const {
    if (dataMap_.count(token)) {
        return dataMap_.find(token)->second;
    }

    throwError(WHERE_AM_I + " unknown token data for get: " + token + " available are: " + tokenList() );
    return *new RVector(0);
}

const IndexArray DataContainer::id(const std::string & token) const {
    if (!dataMap_.count(token)) {
        throwError(WHERE_AM_I + " unknown token data for get: " + token + " available are: " + tokenList() );
    }

    if (!isSensorIndex(token)) {
        throwError(WHERE_AM_I + " token: " + token + " is not an index list: " + tokenList() );
    }

    //RVector idx = dataMap_.find(token)->second;
    IndexArray ret(dataMap_.find(token)->second.size());

    for (Index i = 0; i < ret.size(); i ++){
        ret[i] = Index(dataMap_.find(token)->second[i]);
    }

    return ret;
}

RVector * DataContainer::ref(const std::string & token){
    if (dataMap_.count(token)) {
        return &dataMap_.find(token)->second;
    }

    throwError(WHERE_AM_I + " unknown token data for ref: " + token + " available are: " + tokenList() );
    return NULL;
}

void DataContainer::setDataDescription(const std::string & token, const std::string & description){
//    std::cout << "(this->exists(token)) " << this->exists(token) << std::endl;
    if (this->exists(token)){
        dataDescription_[token] = description;
    }
}

std::string DataContainer::dataDescription(const std::string & token) const {
    if (this->exists(token) && (dataDescription_.count(token))) {
        return dataDescription_.find(token)->second;
    }
    return "";
}

void DataContainer::resize(Index size) {
    for (std::map< std::string, RVector >::iterator it = dataMap_.begin();
            it!= dataMap_.end(); it ++){

        if (isSensorIndex(it->first)){
            // if the data field represent a sensorIDX fill with -1
            it->second.resize(size, -1.0);
        } else {
            // else pure data, fill with 0.0
            it->second.resize(size, 0.0);
        }
    }
}

void DataContainer::removeInvalid(){
    IndexArray validIdx(find(get("valid") == 1));

    for (std::map< std::string, RVector >::iterator it = dataMap_.begin(); it!= dataMap_.end(); it ++){
        dataMap_[it->first] = it->second(validIdx);
    }
}

void DataContainer::remove(const IndexArray & idx){
    this->markValid(idx, false);
    this->removeInvalid();
}

IndexArray DataContainer::findSensorIndex(const RVector & d) const{
    IndexArray ret(d.size());
    for (Index i = 0; i < d.size(); i ++){
        if (d[i] > -1 && d[i] < sensorCount()) {
            ret[i] = Index(d[i]);
        } else {
            throwError(WHERE_AM_I + " Sensor index value out of range 0--"
                        + str(sensorCount()) + " " + str(Index(d[i])));

        }
    }
    return ret;
}

// START Sensor related section
void DataContainer::removeSensorIdx(Index idx){
    IndexArray i(1, idx);
    this->removeSensorIdx(i);
}

void DataContainer::removeSensorIdx(const IndexArray & idx){
    for (std::map< std::string, RVector >::iterator it = dataMap_.begin(); it!= dataMap_.end(); it ++){
        if (isSensorIndex(it->first)){
            for (IndexArray::iterator id = idx.begin(); id != idx.end(); id ++){
                this->markValid(find(it->second == double(*id)), false);
            }
        }
    }
    this->removeInvalid();
    this->removeUnusedSensors();
}

void DataContainer::removeUnusedSensors(bool verbose){

    Index oldSensorCount = this->sensorCount();

    BVector activeSensors(this->sensorCount(), false);

    for (std::map< std::string, RVector >::iterator it = dataMap_.begin(); it!= dataMap_.end(); it ++){
        if (isSensorIndex(it->first)){
            for (Index i = 0; i < it->second.size(); i ++){
                SIndex id = (SIndex)it->second[i];
                if (id > -1 && id < (SIndex)this->sensorCount()) activeSensors[id] = true;
            }
        }
    }

    IndexArray perm(this->sensorCount(), 0);

    R3Vector tmp(sensorPoints_);
    sensorPoints_.clear();

    for (size_t i = 0; i < activeSensors.size(); i ++) {
        if (activeSensors[i]){
            sensorPoints_.push_back(tmp[i]);
        }
        perm[i] = sensorPoints_.size() -1;
    }

    for (std::map< std::string, RVector >::iterator it = dataMap_.begin(); it!= dataMap_.end(); it ++){
        if (isSensorIndex(it->first)){
            for (Index i = 0; i < it->second.size(); i ++){
                SIndex id = (SIndex)it->second[i];
                if (id > -1 && id < (SIndex)perm.size()) it->second[i] = perm[id];
            }
        }
    }

    if (verbose){
        std::cout << "Removed " << oldSensorCount - this->sensorCount() << " sensors." << std::endl;
    }
}

void DataContainer::setSensorPosition(Index i, const RVector3 & pos) {
    if (i >= sensorPoints_.size()) {
        // std::cout << "Warning! .. sensor count was " << sensorCount() << " resize to " << i+1 << std::endl;
        sensorPoints_.resize((i+1));
    }
    sensorPoints_[i] = pos;
}

bool idPosLesserX(const std::pair < RVector3, Index > & a, const std::pair < RVector3, Index > & b){
    return posLesserX(a.first, b.first);
}
bool idPosLesserXrY(const std::pair < RVector3, Index > & a, const std::pair < RVector3, Index > & b){
    return posLesserXrY(a.first, b.first);
}
bool idPosLesserXYrZ(const std::pair < RVector3, Index > & a, const std::pair < RVector3, Index > & b){
    return posLesserXYrZ(a.first, b.first);
}

void DataContainer::sortSensorsX(bool incX, bool incY, bool incZ){
    std::vector < std::pair < RVector3, Index > > permSens(sensorCount());
    for (Index i = 0; i < sensorCount(); i ++) permSens[i] = std::pair< RVector3, Index >(sensorPoints_[i], i);

    if (incX && incY && incZ) {
        std::sort(permSens.begin(), permSens.end(), idPosLesserX);
    } else if (incX && !incY) {
        std::sort(permSens.begin(), permSens.end(), idPosLesserXrY);
    } else if (incX && incY && !incZ) {
        std::sort(permSens.begin(), permSens.end(), idPosLesserXYrZ);
    } else {
        THROW_TO_IMPL
    }

    IndexArray perm(sensorCount());
    for (Index i = 0; i < perm.size(); i ++){
        sensorPoints_[i] = permSens[i].first;
        perm[permSens[i].second] = i;
    }

    for (std::map< std::string, RVector >::iterator it = dataMap_.begin(); it!= dataMap_.end(); it ++){
        if (isSensorIndex(it->first)){
            for (Index i = 0; i < it->second.size(); i ++){
                SIndex id = (SIndex)it->second[i];
                if (id > -1 && id < (SIndex)perm.size()) it->second[i] = perm[id];
            }
        }
    }
}

IndexArray DataContainer::dataIndex(){
    IndexArray ids(this->size());
    Index nSensorsIdx = dataSensorIdx_.size();
    for (Index i = 0; i < this->size(); i ++) {
        Index sensorUniqueID = 0;
        Index count = nSensorsIdx;
        for (std::set< std::string >::iterator it = dataSensorIdx_.begin(); it!= dataSensorIdx_.end(); it ++){
            count --;
            sensorUniqueID += ((Index)dataMap_[*it][i] + 1) * (Index)powInt(this->sensorCount(), count);
        }
        ids[i] = sensorUniqueID;
    }
    return ids;
}

IndexArray DataContainer::sortSensorsIndex(){

    IndexArray ids(this->dataIndex());

    // use IndexArray instead of std::vector would be nice but need some more advanced iterator definition -> TODO
    std::vector< Index > perm(this->size());
    std::iota(perm.begin(), perm.end(), 0);

    auto idxComp = [&ids](Index i1, Index i2) { return ids[i1] < ids[i2]; };
    std::sort(perm.begin(), perm.end(), idxComp);

    for (std::map< std::string, RVector >::iterator it = dataMap_.begin(); it!= dataMap_.end(); it ++){
        it->second = it->second(perm);
    }
    return perm;
}

void DataContainer::markInvalidSensorIndices(){
    for (std::map< std::string, RVector >::iterator it = dataMap_.begin(); it!= dataMap_.end(); it ++){
        if (isSensorIndex(it->first)){
            this->markValid(find(it->second >= double(this->sensorCount())), false);
        }
    }
}

void DataContainer::translate(const RVector3 & trans){
    for (Index i = 0; i < sensorPoints_.size(); i ++) {
        sensorPoints_[i].translate(trans);
    }
}

void DataContainer::scale(const RVector3 & scale){
    for (Index i = 0; i < sensorPoints_.size(); i ++) {
        sensorPoints_[i].scale(scale);
    }
}

// END Sensor related section

} // namespace GIMLI{
