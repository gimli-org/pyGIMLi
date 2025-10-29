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

#include "pos.h"
#include "matrix.h"
#include "elementmatrixmap.h"
#include "meshentities.h"
//#include "vectortemplates.h"

namespace GIMLI{

std::vector < Pos > loadRVector3(const std::string & fileName){
    std::vector < Pos > l;
    std::fstream file; openInFile(fileName, & file, true);
    std::vector < std::string > row;

    while (! file.eof()){
        row = getNonEmptyRow(file);
        switch (row.size()){
            case 1 : l.push_back(RVector3(toDouble(row[0]), 0.0, 0.0)); break;
            case 2 : l.push_back(RVector3(toDouble(row[0]), toDouble(row[1]), 0.0)); break;
            case 3 : l.push_back(RVector3(toDouble(row[0]), toDouble(row[1]),
                                    toDouble(row[2]))); break;
        }
    }
    file.close();
    return l;
}

R3Vector r3(const RVector & x, const RVector & y, const RVector & z){
    ASSERT_EQUAL_SIZE(x, y)
    ASSERT_EQUAL_SIZE(x, z)
    R3Vector ret(x.size());
    for (Index i = 0; i < x.size(); i ++) {
        ret[i] = RVector3(x[i], y[i], z[i]);
    }
    return ret;
}

void vectorizePosVectorList(const std::vector < PosVector > & v, PosVector & r){
    r.clear();
    for (Index i = 0; i < v.size(); i ++) {
        for (Index j = 0; j < v[i].size(); j ++) {
            r.push_back(v[i][j]);
        }
    }

    // Index count = 0;
    // for (Index i = 0; i < v.size(); i ++) {
    //     count += v[i].size();
    // }
    // r.resize(count);
    // Indes start = 0;
    // for (Index i = 0; i < v.size(); i ++) {
    //     r.setVal(v[count += v[i].size();
    // }
}


void deVectorizeRVectorToPosVectorList(std::vector < RVector > & ret,
                                       const RVector & r,
                                       const std::vector < PosVector > & v){
    ret.clear();
    Index start = 0;
    Index end = 0;
    for (Index i = 0; i < v.size(); i ++) {
        end = start + v[i].size();
        ret.push_back(r.getVal(start, end));
        //#std::copy(&data_[se.first], &data_[se.second], &v[0]);
        start = end;
    }
}


void vectorizePosVectorList(const std::vector < PosVector > & v, PosVector & r,
                            int marker, const ElementMatrixMap & eMap){
    ASSERT_EQUAL_SIZE(v, eMap)
    r.clear();
    for (Index i = 0; i < v.size(); i ++) {
        if (eMap.mats()[i].entity()->marker() == marker){
            for (Index j = 0; j < v[i].size(); j ++) {
                r.push_back(v[i][j]);
            }
        }
    }
}


void deVectorizeRVectorToPosVectorList(std::vector < RVector > & ret,
                                       const RVector & r,
                                       const std::vector < PosVector > & v,
                                       int marker,
                                       const ElementMatrixMap & eMap){
    if (ret.size() != eMap.size()){
        ret.resize(eMap.size(), RVector(0));
    }

    Index start = 0;
    Index end = 0;
    for (Index i = 0; i < v.size(); i ++) {
        if (eMap.mats()[i].entity()->marker() == marker){
            end = start + v[i].size();
            ret[i] = r.getVal(start, end);
            //#std::copy(&data_[se.first], &data_[se.second], &v[0]);
            start = end;
        }
    }
}

void dot(const PosVector & a,
         const PosVector & b,
         RVector & r){

    ASSERT_EQUAL_SIZES(a, b)
    r.resize(a.size());
    for (Index i = 0; i < a.size(); i ++ ){
        r[i] = a[i].dot(b[i]);
    }
}

void dot(const std::vector < PosVector > & a,
         const std::vector < PosVector > & b,
         std::vector < RVector > & r){

    ASSERT_EQUAL_SIZES(a, b)
    r.resize(a.size());
    for (Index i = 0; i < a.size(); i ++ ){
        dot(a[i], b[i], r[i]);
    }
}

void dot(const stdVectorRMatrix & a,
         const stdVectorRMatrix & b,
         RVector & r){

    ASSERT_EQUAL_SIZES(a, b)
    r.resize(a.size());
    for (Index i = 0; i < a.size(); i ++ ){
        r[i] = (a[i]*b[i]).sum();
    }
}

void dot(const stdVectorMatrixVector & a,
         const stdVectorMatrixVector & b,
         std::vector < RVector > & r){

    ASSERT_EQUAL_SIZES(a, b)
    r.resize(a.size());

    for (Index i = 0; i < a.size(); i ++ ){
        dot(a[i], b[i], r[i]);
    }
}

template < class Val > void _sum_vec_T(const Val & a, RVector & r){
    r.resize(a.size());
    for (Index i = 0; i < a.size(); i ++ ){
        r[i] = a[i].sum();
    }
}
void sum(const PosVector & a, RVector & r){
    _sum_vec_T(a, r);
}
void sum(const stdVectorRMatrix & a, RVector & r){
    _sum_vec_T(a, r);
}

template < class Val > void _sum_vec_vec_T(const std::vector < Val > & a,
                                           std::vector < RVector > & r){
    r.resize(a.size());
    for (Index i = 0; i < a.size(); i ++ ){
        sum(a[i], r[i]);
    }
}
void sum(const std::vector < PosVector > & a,
         std::vector < RVector > & r){
    _sum_vec_vec_T(a, r);
}

void sum(const std::vector < std::vector< GIMLI::RMatrix > > & a,
                   std::vector < RVector > & r){
    _sum_vec_vec_T(a, r);
}

void trace(const stdVectorRMatrix & a, RVector & r){
    r.resize(a.size());
    for (Index i = 0; i < a.size(); i ++ ){
        r[i] = a[i].trace();
    }
}

RVector trace(const stdVectorRMatrix & a){
    RVector ret;
    trace(a, ret);
    return ret;
}

void trace(const stdVectorMatrixVector & a, stdVectorRVector & r){
    r.resize(a.size());
    for (Index i = 0; i < a.size(); i ++ ){
        trace(a[i], r[i]);
    }
}

stdVectorRVector trace(const stdVectorMatrixVector & a){
    stdVectorRVector ret;
    trace(a, ret);
    return ret;
}



void saveRVector3(const std::vector < Pos > l, const std::string & fileName){
    std::fstream file; openOutFile(fileName, & file);
    for (uint i = 0; i < l.size(); i ++) file << l[i] << std::endl;
    file.close();
}

RVector3 center(const R3Vector & vPos){
  Pos pos(0.0, 0.0, 0.0);

  //  std::cout << WHERE_AM_I << vPos.size() << std::endl;

  if (vPos.size() == 0) {
    pos.setValid(false);
    return pos;
  }

  for (uint i = 0; i < vPos.size(); i ++) pos += vPos[i];
  return pos /= vPos.size();
}

R3Vector normalise(const R3Vector & vPos){
    R3Vector ret(vPos);
    for (uint i = 0; i < vPos.size(); i ++) ret[i] = ret[i]/ret[i].abs();
    return ret;
}

double jacobianDetXY(const Pos & p1, const Pos & p2, const Pos & p3){
  double x1 = p1.x(), x2 = p2.x(), x3 = p3.x();
  double y1 = p1.y(), y2 = p2.y(), y3 = p3.y();

  return (x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1);
}

double angle(const Pos & p1, const Pos & p2, const Pos & p3){
  double angle = p2.angle(p1, p3);
  if (jacobianDetXY(p1, p2, p3) > 0) angle = 2.0 * PI - angle;
  return angle;
}

bool xVari(const R3Vector & electrodeList){
    if (electrodeList.empty()) return false;
    double start = electrodeList[0].x();
    for (size_t i = 1; i < electrodeList.size(); i ++){
        if (std::fabs(electrodeList[i].x() - start) > TOLERANCE) return true;
    }
    return false;
}
bool yVari(const R3Vector & electrodeList){
    if (electrodeList.empty()) return false;
    double start = electrodeList[0].y();
    for (size_t i = 1; i < electrodeList.size(); i ++){
        if (std::fabs(electrodeList[i].y() - start) > TOLERANCE) return true;
    }
    return false;
}
bool zVari(const R3Vector & electrodeList){
    if (electrodeList.empty()) return false;
    double start = electrodeList[0].z();
    for (size_t i = 1; i < electrodeList.size(); i ++){
        if (std::fabs(electrodeList[i].z() - start) > TOLERANCE) return true;
    }
    return false;
}

RVector x(const R3Vector & rv){
    RVector t(rv.size());
    //std::copy(rv.begin(), rv.end(), &t[0], std::mem_fun_ref(&Pos< double >::x));
    for (uint i = 0, imax = rv.size(); i < imax; i ++) t[i] = rv[i].x();
    return t;
}

RVector y(const R3Vector & rv){
    RVector t(rv.size());
    for (uint i = 0, imax = rv.size(); i < imax; i ++) t[i] = rv[i].y();
    return t;
}

RVector z(const R3Vector & rv){
    RVector t(rv.size());
    for (uint i = 0, imax = rv.size(); i < imax; i ++) t[i] = rv[i].z();
    return t;
}
template < class ValueType > void swap(ValueType & v1, ValueType & v2){
    ValueType tmp = v1; v1 = v2; v2 = tmp;
}

RVector absR3(const R3Vector & vPos){
    RVector r(vPos.size());
    for (Index i = 0; i < vPos.size(); i ++ ) r[i] = vPos[i].abs();
    return r;
}

void swapXY(R3Vector & rv){ for (uint i = 0, imax = rv.size(); i < imax; i ++) swap(rv[i][0], rv[i][1]); }
void swapXZ(R3Vector & rv){ for (uint i = 0, imax = rv.size(); i < imax; i ++) swap(rv[i][0], rv[i][2]); }
void swapYZ(R3Vector & rv){ for (uint i = 0, imax = rv.size(); i < imax; i ++) swap(rv[i][1], rv[i][2]); }

RVector toArray(const R3Vector & vec){
    RVector ret(vec.size() * 3);
    for (Index i = 0; i < vec.size(); i ++) {
        ret[i*3] = vec[i][0];
        ret[i*3 + 1] = vec[i][1];
        ret[i*3 + 2] = vec[i][2];
    }
    return ret;
}

RMatrix toMatrix(const R3Vector & vec){
    RMatrix ret(vec.size(), 3);
    for (Index i = 0; i < vec.size(); i ++) {
        std::copy(&vec[i][0], &vec[i][2], &ret[i][0]);
    }
    return ret;
}

R3Vector stdVectorR3VectorToR3Vector(const std::vector < Pos > & rv){
    R3Vector ret(rv.size());
    for (Index i = 0; i < rv.size(); i ++) ret[i] = rv[i];
    return ret;
}

std::vector < Pos > R3VectorTostdVectorR3Vector(const R3Vector & rv){
    std::vector < Pos > ret(rv.size());
    for (Index i = 0; i < rv.size(); i ++) ret[i] = rv[i];
    return ret;
}

RVector3 RVector3::cross(const Pos & p) const{
//     r[0] = (a2 * b3 - a3 * b2);
//     r[1] = (a3 * b1 - a1 * b3);
//     r[2] = (a1 * b2 - a2 * b1);
  return (RVector3((*this)[1] * p[2] - (*this)[2] * p[1],
                   (*this)[2] * p[0] - (*this)[0] * p[2],
                   (*this)[0] * p[1] - (*this)[1] * p[0]));
}

RVector3 RVector3::norm(const Pos & p1, const Pos & p2) const {
    Pos a(p1 - (*this));
    Pos b(p2 - (*this));
    Pos r(a.cross(b));
    return r.norm();
}

RVector3 RVector3::normXY(const Pos & p) const {
    Pos dir(*this - p);
    return Pos(-dir[1], dir[0]).normalize();

    // Pos result((*this + p) / 2.0);
    // result.setZ(1.0);
    // result = result.norm(*this, p);
    // result.setZ(0.0);
    // return result;
}

double Pos::angle(const Pos & p1, const Pos & p3) const {
    Pos a(p1 - (*this));
    Pos b(p3 - (*this));
    return (a).angle(b);
}

double Pos::angle(const Pos & p) const {
    double result = acos(this->dot(p) / (this->abs() * p.abs()));
    if (isnan(result) || isinf(result)){
        result = 0.0;
    }
    return result;
}

// template <> double RVector3::angle(const Pos & p) const {
//   double result = acos(this->dot(p) / (this->abs() * p.abs()));
//   if (std::isnan(result) || std::isinf(result)) {
//     result = 0.0;
//   }
//   return result;
// }

// template <> double RVector3::angle(const Pos & p1, const Pos & p3) const {
//   Pos a(p1 - (*this));
//   Pos b(p3 - (*this));
//   return (a).angle(b);
// }

} // namespace GIMLI
