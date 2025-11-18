/******************************************************************************
 *   Copyright (C) 2012-2024 by the GIMLi development team                    *
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

#ifndef _GIMLI_MEMWATCH__H
#define _GIMLI_MEMWATCH__H

#include "gimli.h"

#define MEMINFO GIMLI::MemWatch::instance().info(WHERE);

namespace GIMLI{

/*! Convert byte into KByte. */
inline double KByte(long byte){ return double(byte / 1024.0); }

/*! Convert byte into MByte */
inline double MByte(long byte){ return double(byte / (1024.0 * 1024.0)); }

/*! Convert byte into GByte */
inline double GByte(long byte){ return double(byte / (1024.0 * 1024.0 * 1024)); }

//! Memory watch.
/*! Class that might help debugging memory usage.
 * Informations are taken from /proc system, so only available for linux systems.
 * This is a singleton class to ensure a single instance.
 * To call it use e.g.: MemWatch::instance().info(WHERE);*/
class DLLEXPORT MemWatch : public Singleton< MemWatch > {
public:
    friend class Singleton< MemWatch >;

    /*! Return the current memory usage of the process. Values are in MByte. */
    double inUse();

    /*! Return the current memory usage relative to the last call of this method. Values are in MByte. */
    double current();

    /*! Shows the current and the relative (to the last call) memory usage. */
    void info(const std::string & str="");

protected:
    double last_;

private:
    /*! Private so that it can not be called */
    MemWatch();
    /*! Private so that it can not be called */
    virtual ~MemWatch();
    /*! Copy constructor is private, so don't use it */
    MemWatch(const MemWatch &){};
    /*! Assignment operator is private, so don't use it */
    void operator = (const MemWatch &){ };

    Stopwatch * swatchAll_;
    Stopwatch * swatchDur_;

};

/*! Current amount of memory in use for the current process in MByte. */
inline double memoryInUse(){
    return GIMLI::MemWatch::instance().inUse();
}

DLLEXPORT long maxMem();

} // namespace GIMLI

#endif
