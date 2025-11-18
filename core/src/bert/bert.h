/******************************************************************************
 *   Copyright (C) 2006-2024 by the resistivity.net development team          *
 *   Carsten R�cker carsten@resistivity.net                                   *
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

#ifndef _BERT_BERT__H
#define _BERT_BERT__H

#ifdef bert_EXPORTS
	#define gimli_EXPORTS
#endif

#include <gimli.h>

namespace GIMLI{

static const int MARKER_NODE_ELECTRODE = -99;
static const int MARKER_NODE_REFERENCEELECTRODE = -999;
static const int MARKER_NODE_CALIBRATION = -1000;
static const int MARKER_BOUND_ELECTRODE = -10000;

class ElectrodeShape;
class DataContainerERT;
class DataMap;

}

#include "bertMisc.h"
#include "bertJacobian.h"
#include "datamap.h"
#include "dcfemmodelling.h"
#include "bertDataContainer.h"
#include "electrode.h"

#endif // _BERT_BERT__H

