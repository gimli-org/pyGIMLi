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

#ifndef _BERT_DATAMAP__H
#define _BERT_DATAMAP__H

#include "bert.h"

#include <matrix.h>
#include <pos.h>

namespace GIMLI{

/*! Potential matrix for BERT.
Stores the potential values for every current injection at every electrode positions. Filtering a specific data configuration is done by the \ref data call. In generall, the map size is equal to the amount of electrodes. There are exceptions: Dipol-Pattern sources, Dipol-Sources against a reference Electrode. */
class DLLEXPORT DataMap {
public:
    /*! Default contructor, builds an empty map */
    DataMap();

    /*! Contructs a map and load content from a file */
    DataMap(const std::string & filename);

    /*! Copyconstructor */
    DataMap(const DataMap & map);

    /*! Assignment operator */
    DataMap & operator = (const DataMap & map);

    /*! Complex valued entries.*/
    bool isComplex() const { return complex_; }

    /*! Save the collect matrix */
    int save(const std::string & filename);

    /*! Load a collect matrix */
    int load(const std::string & filename);

    /*! Fill the map. Collect potential values from complete mesh related solution matrix
        and a vector of assiciated electrodes. */
    void collect(const std::vector < ElectrodeShape * > & electrodes,
                 const RMatrix & sol, bool isCEM=false);

    /*! Returns \ref RVector of power values corresponding to the \ref DataConatiner,
     * Return reciprocity values if reciprocity set to True.
     * In the case of reciprocity, a <-> m and b <-> n are swapped. */
    RVector data(const DataContainerERT & dat, bool reciprocity=false, bool imag=false);

    /*! Set the electrode positions */
    inline void setElectrodes(const std::vector < RVector3 > & elecs){ elecs_ = elecs; }

    /*! Return a reference to the electrode positions */
    inline const std::vector < RVector3 > & electrodes() const { return elecs_; }

    /*! Return a reference to the potential matrix */
    inline const RMatrix & map() const { return map_; }

protected:
    /*! Internal copy function */
    void copy_(const DataMap & map);

    /*! Hold electrode positions */
    std::vector < RVector3 > elecs_;

    /*! Hold potential map */
    RMatrix map_;

    bool complex_;
};

} //namespace BERT

#endif // _BERT_DATAMAP__H
