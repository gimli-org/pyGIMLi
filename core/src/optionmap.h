/******************************************************************************
 *   Copyright (C) 2009-2024 by the GIMLi development team                    *
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

#ifndef GIMLI_LONGOPTIONS__H
#define GIMLI_LONGOPTIONS__H

#include <list>
#include <map>
#include <string>
#include <iostream>
#include <algorithm>
#include <cstdio>

#ifdef _MSC_VER
	#include <wingetopt.h>
#else
	#include <getopt.h>
#endif

#include "gimli.h"

namespace GIMLI{

class DLLEXPORT OptionBase{
public:
    virtual ~OptionBase(){}

    virtual void assign_(char * opt=NULL) = 0;

    virtual std::string typname() = 0;

    virtual std::string defautToString() = 0;

    inline const char key() const { return key_; }

    void assign(char * opt){
        if (hasArg_ == no_argument) assign_(NULL);
        if (hasArg_ == required_argument) assign_(opt);
    }

    void            * var_;
    char              key_;
    std::string       name_;
    std::string       help_;
    int               hasArg_;
};

inline std::ostringstream & operator << (std::ostringstream & s, const std::vector < std::string > & str){
    for (uint i = 0; i < str.size(); i ++) s << str[ i ] << " ";
    return s;
}

template < class T > class Option : public OptionBase{
public:
    Option(T & var){
        var_ = &var;
    }
    Option(T & var, const T & defaultVar){
        var_ = &var;
        *(T*)(var_) = defaultVar;
        defaultVar_ = defaultVar;
    }

    virtual ~Option(){}

    virtual void assign_(char * opt = NULL){ convert(*(T*)(var_), opt); }

    virtual std::string typname(){ return type(*(T*)(var_)); }

    T value() { return *(T*)(var_); }

    virtual std::string defautToString() {
        std::ostringstream streamOut;
        streamOut << defaultVar_;
        return streamOut.str();
    }

    T defaultVar_;
};

//! Simplified command line parser
/*! Simplified command line parser. The default option are -h [--help]. --version */
class DLLEXPORT OptionMap{
public:
    OptionMap(const std::string & description = "");

    virtual ~OptionMap();

    typedef std::map< char,         OptionBase * >    CharOptionMap;
    typedef std::map< std::string,  OptionBase * >    LongOptionMap;

    template < class T > void addLastArg(T & var, const std::string & lastArgString){
        lastArgString_ = lastArgString;
        lastOption_ = new Option< T >(var);
        lastOption_->hasArg_ = required_argument;
    }

    template < class T > void add(T & var, const std::string & key, const std::string & longOpt,
                                    const std::string & help){
        add(var, key, longOpt, help, T(var));
    }

    template < class T > void add(T & var, const std::string & key, const std::string & longOpt,
                                    const std::string & help, const T & defaultVar){

        if (key[ 0 ] != ':') allKeys_ += key;

        OptionBase * o = new Option< T >(var, defaultVar);
        o->name_ = longOpt;
        if ((*--key.end()) == ':' || (*--longOpt.end() ==':')) o->hasArg_ = required_argument;
        else o->hasArg_ = no_argument;

        std::string k( key.substr(0, key.rfind(':')));
        std::string lo(longOpt.substr(0, longOpt.rfind(':')));

        if (k.length() > 0) {
            o->key_ = k[ 0 ];
            cMap_.insert(std::pair< char, OptionBase * >(k[ 0 ], o));
        }

        lMap_.insert(std::pair< std::string, OptionBase * >(lo, o));
        options_.push_back(o);

        o->help_ = help;
    }

    void parse(int argc, char * argv[]);

    void buildLongOptions();

    inline void setDescription(const std::string & description){ descriptionString_ = description; }

    void printHelp(const std::string & main);

    OptionBase * findOption(const char & key) const {
        log(Warning, "command line parsing from core is obsolete and will be replaced.");
        THROW_TO_IMPL
        return NULL;
    }

protected:
    struct option               * opts_;
    std::string                   allKeys_;
    std::string                   lastArgString_;
    std::string                   descriptionString_;

    OptionBase                  * lastOption_;
    std::list< OptionBase * >     options_;
    LongOptionMap                 lMap_;
    CharOptionMap                 cMap_;
    bool                          showVersion_;
    bool                          debug_;
};

} // namespace GIMLI
#endif //GIMLI_LONGOPTIONS__H
