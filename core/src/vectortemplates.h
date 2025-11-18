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

#ifndef _GIMLI_VECTORTEMPLATES__H
#define _GIMLI_VECTORTEMPLATES__H

#include "gimli.h"
#include <cmath>
#include <cstdlib>
#include <cerrno>
#include <cstring>
#include <string>
#include <vector>
#include <fstream>
#include <algorithm>
#include <numeric>

namespace GIMLI{

template < typename T, class Iter, template < typename, class > class Vec >
// template < typename T, template < typename, class > class Vec > 
IndexArray ids(const Vec < T, Iter > & e){
    IndexArray id(e.size());
    for (Index i = 0; i < e.size(); i ++ ) id[i] = e[i]->id();
    return id;
}

template < class ValueType > bool save(std::vector < ValueType > & a, const std::string & filename, IOFormat format = Ascii){
    return saveVec(a, filename, format);
}

template < class ValueType > bool load(std::vector < ValueType > & a, const std::string & filename, IOFormat format = Ascii,
                                  bool verbose = true){
    return loadVec(a, filename, format, verbose);
}

template < class Vec > bool saveVec(const Vec & a, const std::string & filename,
                                     IOFormat format, bool verbose = true){

    if (filename.rfind(VECTORASCSUFFIX) != std::string::npos) format = Ascii;
    else if (filename.rfind(VECTORBINSUFFIX) != std::string::npos) format = Binary;
    std::string fname(filename);

    if (format == Ascii){
        if (fname.rfind(".") == std::string::npos) fname += VECTORASCSUFFIX;

        std::ofstream file; file.open(fname.c_str());
        if (!file) {
            std::cerr << filename << ": " << strerror(errno) << " " << errno << std::endl;
            return false;
        }

        file.setf(std::ios::scientific, std::ios::floatfield);
        file.precision(14);

        for (uint i = 0, imax = a.size(); i < imax; i ++) file << a[i] << std::endl;
        file.close();
    } else {
        if (fname.rfind(".") == std::string::npos) fname += VECTORBINSUFFIX;
    // so haett ich das gern //
//     std::ofstream file(filename.c_str(), std::ofstream::binary);
//     std::copy(&a[0], &a[a.size()-1], ostream_iterator< double >(&file));
//     file.close();
        FILE *file; file = fopen(fname.c_str(), "w+b");
        if (!file) {
            if (verbose) std::cerr << filename << ": " << strerror(errno) << " " << errno << std::endl;
            return false;
        }

        int count = a.size();
        uint ret = 0; ret = fwrite((char*)&count, sizeof(int), 1, file);
        if (ret){
			for (uint i = 0; i < a.size(); i++) ret = fwrite((char*)&a[i], sizeof(double), 1, file);
		}
        fclose(file);
    }
    return true;
}

template < class Vec > bool loadVec(Vec & a, const std::string & filename,
                                     IOFormat format, bool verbose = true){

    if (filename.rfind(VECTORASCSUFFIX) != std::string::npos) format = Ascii;
    else if (filename.rfind(VECTORBINSUFFIX) != std::string::npos) format = Binary;

    if (!fileExist(filename)){
        if (fileExist(filename + VECTORBINSUFFIX))
            return loadVec(a, filename + VECTORBINSUFFIX, Binary);
        if (fileExist(filename + VECTORASCSUFFIX))
            return loadVec(a, filename + VECTORASCSUFFIX, Ascii);
    }

    if (format == Ascii){

        a.clear();
        std::fstream file; openInFile(filename.c_str(), &file);
        double val; while(file >> val) a.push_back(val);

    //so haett ich das gern
//     std::ifstream file(filename.c_str());
//     std::copy( std::istream_iterator<double>(file),
//                 std::istream_iterator<double>(),
//                 std::back_inserter(tmp));

//std::back_inserter< double > (tmp));
    //std::copy(file.begin(), file.end(), back_inserter< double >(& tmp[0]));

//         a.resize(tmp.size());
//         std::copy(tmp.begin(), tmp.end(), &a[0]);
        file.close();

  } else {
    FILE *file;
    file = fopen(filename.c_str(), "r+b");
    if (!file) {
        if (verbose) std::cerr << filename << ": " << strerror(errno) << " " << errno << std::endl;
        return false;
    }
    uint ret = 0;
    int size; ret = fread(&size, sizeof(int), 1, file);
    if (ret){
        a.resize(size);
        ret = fread(&a[0], sizeof(double), size, file);
    }else {

    }
    fclose(file);
  }
  return true;
}


//** START std::vector shortcuts
template < class T > std::vector< T > sort(const std::vector < T > & a){
    std::vector < T > t(a);
    sort(t.begin(), t.end());
    return t;
}

template < class T > std::vector< T > unique(const std::vector < T > & a){
    std::vector < T > t;
    std::unique_copy(a.begin(), a.end(), back_inserter(t));
    return t;
}

//** END std::vector shortcuts

// template < class Vec > Vec abs(const Vec & v) {
//   Vec tmp(v.size());
//   for (uint i = 0; i < v.size(); i ++) tmp[i] = abs(v[i]);
//   return tmp;
// }

template < class Vec > void clear(Vec & a) {
    for (int i = 0; i < (int)a.size(); i ++) a[i] = 0.0;
}

template < typename T, class Iter, template < typename, class > class Vec >
T min(const Vec< T, Iter > & v){
    return *std::min_element(&v[0], &v[0] + v.size());
}

template < typename T, class Iter, template < typename, class > class Vec >
T max(const Vec< T, Iter > & v){
    return *std::max_element(v.begin(), v.end());
}
//
// template < class T, class Iter, template < class T, class Iter > class Vec > T sum(const Vec< T, Iter > & v){
//     return std::accumulate(v.begin(), v.end(), 0.0);
// //     T sum = 0.0;
// //     for (size_t i = v.size(); i--;) sum += v[i];
// //     return sum;
// }
//
// template < class ValueType, template < class T > class Vec > T dot(const Vec< T > & a, const Vec < T > & b) {
//     return sum(a * b);
// }

template < class Vec > void echoMinMax(const Vec & vec, const std::string & name){
    if (vec.size() > 0){
        std::cout << "min " << name << " = " << min(vec)
              << " max " << name << " = " << max(vec) << " (" << vec.size() << ")" << std::endl;
    } else {
        std::cout << "min " << name << " = ndef."
              << " max " << name << " = ndef." << std::endl;
    }
}

template < class Vec > double median(const Vec & a) {
    Index dim = a.size();
    if (dim == 1) return a[0];
    if (dim > 1){
        Vec tmp(sort(a));
        if (std::fabs(dim / 2.0 - rint(dim / 2.0)) < 1e-12){ // even
            return (tmp[dim / 2] + tmp[dim / 2 - 1]) / 2.0;
        } else { // odd
            return tmp[(dim - 1)/ 2];
        }
    }
    return 0.0;
}

// template < class Vec > double mean(const Vec & a) { return sum(a) / a.size(); }

// template < class Vec > double stdDev(const Vec & a) {
//     return std::sqrt(sum(square(a - mean(a))) / (a.size() - 1));
// }

template < class Vec > double arithmeticMean(const Vec & a) { return mean(a); }

template < class Vec > double geometricMean(const Vec & a) {
    int dim = a.size();
    double result = 0.0;
    for (int i = 0; i < dim; i ++) result += std::log(a[i]);
    result /= (double)dim;
    return std::exp(result);
}

template < class Vec > double harmonicMean(const Vec & a) {
    int dim = a.size();
    double result = 1.0 / a[0];

    for (int i = 1; i < dim; i ++) result += 1.0 / a[i];
    result /= dim;
    return 1.0 / result;
}

template < class Vec > double rms(const Vec & a) { return std::sqrt(mean(square(a))); }

template < class Vec > double rms(const Vec & a, const Vec & b) { return rms(a - b); }

template < class Vec > double rrms(const Vec & a, const Vec & b) { return rms((a - b) / a); }

template < class Vec > double normlp(const Vec & a, int p) {
  return std::pow(sum(pow(abs(a), p)), 1.0/(double)p);
}

template < class Vec > double norml1(const Vec & a) {
  //http://mathworld.wolfram.com/L1-Norm.html
  return normlp(a, 1);
}

template < class Vec > double norml2(const Vec & a) {
  // vector norm \ell^2 nicht L^2
  return normlp(a, 2);
}
template < class Vec > double normlInfinity(const Vec & a) {
  return max(abs(a));
}
template < class Vec > double euclideanNorm(const Vec & a) {
  return norml2(a);
}

template < class Vec > double norm(const Vec & a) {
  //** sqrt(a_0^2 + a_i^2 + a_n^2) ; 0 < i < a.size()
  return norml2(a);
}

template < class Vec > double chiQuad(const Vec & a, const Vec & b, const Vec & err) {
    Vec tmp((a - b) / err);
    return dot(tmp, tmp) / a.size();
}

template < class Vec > double chiQuad(const Vec & a, const Vec & b, const Vec & err, bool isLog) {
  double chiq = 0.0;

  Vec tmp(a.size());

  if (isLog){
    tmp = ((log(b) - log(a)) / log(err + 1.0));
    chiq = mean(tmp * tmp);
  } else {
    chiq = mean(((b / a - 1.0) * (b / a - 1.0)) / (err * err));
  }
  return chiq;
}

template < class ValueType >
void rand(Vector < ValueType > & vec, ValueType min = 0.0, ValueType max = 1.0){
    for (int i = 0, imax = vec.size(); i < imax; i ++){
        vec[i] = (ValueType)std::rand() * ((max -min)/ RAND_MAX) + min;
    }
}

template < class ValueType >
void randn(Vector< ValueType > & vec){
    //   vec = sqrt(log(1.0 - vec) * -2.0) * cos(vec * 2.0 * PI_);
    //** keine Ahnung was passiert aber es funktioniert. vgl. tests/vector/baseio.cpp
    for (uint i = 0, imax = vec.size(); i < imax; i ++){
        double sum = 0.0;
        for (int j = 0; j < 16; j++) sum += std::rand() & 0xfff;
        vec[i] = (sum - 0x8000) * (1.0 / 4729.7);
    }
}

/*!Create a array of len n with normal distributed randomized values.*/
inline RVector randn(Index n){
    RVector r(n);
    randn(r);
    return r;
}


// //*!Position of the minimum Value position as IndexArray */
// template < class ValueType >
// IndexArray minIdx(Vector< ValueType > & vec){
// minId

} // namespace GIMLI{

#endif // _GIMLI_VECTORTEMPLATES__H
