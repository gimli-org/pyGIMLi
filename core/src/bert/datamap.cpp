/******************************************************************************
 *   Copyright (C) 2005-2024 by the resistivity.net development team          *
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

#include "datamap.h"

#include "bertDataContainer.h"
#include "electrode.h"

#include <interpolate.h>
#include <node.h>
#include <meshentities.h>

#include <fstream>

namespace GIMLI{

DataMap::DataMap(){
    complex_ = false;
}

DataMap::DataMap(const std::string & filename){
    this->load(filename);
}

DataMap::DataMap(const DataMap & map){
    this->copy_(map);
}

DataMap & DataMap::operator = (const DataMap & map){
    if (this != & map){
        this->copy_(map);
    }
    return *this;
}

void DataMap::copy_(const DataMap & map){
    elecs_ = map.electrodes();
    map_ = map.map();
    complex_ = map.isComplex();
}

void DataMap::collect(const std::vector < ElectrodeShape * > & electrodes,
                      const RMatrix & sol, bool isCEM){
    elecs_.clear();
    map_.clear();
    for (uint i = 0; i < electrodes.size(); i ++) {
        if (electrodes[i]){
            if (electrodes[i]->id() > -1){
                elecs_.push_back(electrodes[i]->pos());
            }
        } else {
            std::cerr << WHERE_AM_I << " no position found for electrode " << i << std::endl;
        }
    }

    if (sol.rows() == 0){
        throwLengthError(WHERE_AM_I + " sol.size() == 0 ");
    }

    if ((sol.rows()     == elecs_.size() && sol.cols() == elecs_.size()) ||
        (sol.rows()     == elecs_.size() && sol.cols() == elecs_.size() + 1) ||
        (sol.rows() + 1 == elecs_.size() && sol.cols() == elecs_.size()) ||
        isCEM ){

        map_ = sol;
    } else if((sol.rows()     == 2. * elecs_.size() && sol.cols() == elecs_.size()) ){
        map_ = sol;
        complex_ = true;
    } else {
        RVector tmp(elecs_.size());
        for (size_t row = 0; row < sol.rows(); row ++) {
//             __MS(row)
//             __MS(sol.rows())
            for (size_t i = 0; i < tmp.size(); i ++) {
//                 __MS(i)
//                 __MS(sol.rows())
//                 __MS(sol.cols())
                tmp[i] = electrodes[i]->pot(sol[row]);
            }
            map_.push_back(tmp);
        }
    }
}

int DataMap::save(const std::string & filename){
    std::fstream outfile; if (!openOutFile(filename, & outfile)){ return -1; }
    uint nElecs = elecs_.size();

    outfile << nElecs << std::endl;
    for (unsigned int i = 0; i < elecs_.size(); i ++){
        outfile << elecs_[i] << std::endl;
    }

    outfile.setf(std::ios::scientific, std::ios::floatfield);
    outfile.precision(14);

    for (size_t i = 0; i < map_.rows(); i ++){
        for (size_t j = 0; j < map_[i].size(); j ++){
            outfile << map_[i][j] << "\t";
        }
        outfile << std::endl;
    }
    outfile.close();
    return 1;
}

int DataMap::load(const std::string & filename){
    std::fstream file; if (!openInFile(filename, & file)){ return -1; }

    std::vector < std::string > row; row = getNonEmptyRow(file);
    size_t nElecs = toInt(row[0]);
    for (size_t i = 0; i < nElecs; i ++){
        row = getNonEmptyRow(file);

        switch (row.size()){
            case 2:
                elecs_.push_back(RVector3(toDouble(row[0]), toDouble(row[1]), 0.0));
                break;
            case 3:
                elecs_.push_back(RVector3(toDouble(row[0]),
                                            toDouble(row[1]),
                                            toDouble(row[2])));
                break;
            default:
                std::cerr << WHERE_AM_I << " WARNING: rowsize unknown " << row.size() << std::endl;
                THROW_TO_IMPL
                break;
        }
    }

    while (!file.eof()){
        row = getNonEmptyRow(file);
        if (row.size() > 0){
            RVector tmp(row.size());
            for (size_t i = 0; i < row.size(); i ++) tmp[i] = toDouble(row[i ]);
            map_.push_back(tmp);
        }
    }

    file.close();
    return 1;
}

RVector DataMap::data(const DataContainerERT & data, bool reciprocity, bool imag){
    size_t nData = data.size();
    size_t nElecs = elecs_.size();

    Index offSet = 0;
    if (imag){
        if (map_.rows() == 2 * nElecs){
            offSet = nElecs;
        } else {
            throwError(WHERE_AM_I + " imaginary values requested but not calculated");
        }
    }

    RVector uVec(nData);
    for (size_t i = 0; i < nData; i ++){
        int a = data("a")[i];
        int b = data("b")[i];
        int m = data("m")[i];
        int n = data("n")[i];

        if (reciprocity){
            std::swap(a, m);
            std::swap(b, n);
        }

        if ((a > (int)nElecs - 1) || (a < -1) || (b > (int)nElecs - 1) || (b < -1) ||
             (m > (int)nElecs - 1) || (m < -1) || (n > (int)nElecs - 1) || (n < -1)) {
            std::stringstream str1;
            str1 << WHERE_AM_I << " Collect matrix to small. data: "
                << i << " Number of electrodes = "
                << nElecs
                << "; a = " << a+1  << " b = " << b+1 << " m = " << m+1 << " n = " << n+1 << std::endl;
            //dat.save("tmp.data");
            throwLengthError(str1.str());
        }

        double uAM = 0.0, uAN = 0.0, uBM = 0.0, uBN = 0.0;

        if (map_.rows() == nElecs || map_.rows() == 2*nElecs){
            if (a != -1 && m != -1) uAM = map_[a + offSet][m];
            if (a != -1 && n != -1) uAN = map_[a + offSet][n];
            if (b != -1 && m != -1) uBM = map_[b + offSet][m];
            if (b != -1 && n != -1) uBN = map_[b + offSet][n];
        } else {
            if (map_.rows() + 1 == nElecs && map_.cols() == nElecs){
                //!** probably measured agains last electrode as current and power reference
                if (a == (int)nElecs -1) a = -1;
                if (b == (int)nElecs -1) b = -1;

                if (a != -1 && m != -1) uAM = map_[a][m];
                if (a != -1 && n != -1) uAN = map_[a][n];
                if (b != -1 && m != -1) uBM = map_[b][m];
                if (b != -1 && n != -1) uBN = map_[b][n];
            }
        }

        uVec[i] = (uAM - uAN) - (uBM - uBN);

        if (std::fabs(uVec[i]) < TOLERANCE/1000. || std::isnan(uVec[i])){
            if (!imag){
                
                std::stringstream str1;
                str1 << WHERE_AM_I << std::endl << " a = " << a
                    << " b = " << b << " m = " << m << " n = " << n 	<< std::endl
                    <<  " " << uAM <<  " " << uAN <<  " " << uBM <<  " " << uBN
                    << " U = " << uVec[i] << std::endl;
                log(Warning, str1.str());
            }
        }
    }

    return uVec;
}

} //namespace GIMLI
