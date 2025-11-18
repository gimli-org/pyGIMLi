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

#include "optionmap.h"
#include <cstdio>

namespace GIMLI{

OptionMap::OptionMap(const std::string & description)
    : descriptionString_(description) {
    opts_           = NULL;
    lastArgString_  = "";
    lastOption_     = NULL;
    bool help;
    showVersion_    = false;
    debug_          = false;
    add(help,              "h" , "help"    , "Show this help.");
    add(showVersion_,      ""  , "version" , "Show library version and exit programm.");
    add(debug_,            ""  , "debug"   , "Enable global debug mode.");
}

OptionMap::~OptionMap(){
    for (std::list < OptionBase * >::iterator it = options_.begin(), itmax = options_.end(); it != itmax; it++) delete *it;
    delete lastOption_;
    if (opts_) delete[] opts_;
}

void OptionMap::parse(int argc, char * argv[]){

    int option_char = 0, option_index = 0;

    buildLongOptions();
    //std::cout << allKeys_ << std::endl;
    while ((option_char = getopt_long(argc, argv, allKeys_.c_str()
              , opts_, & option_index)) != EOF){

        if (option_char == 'h') {
            printHelp(argv[0]);
            exit(EXIT_SUCCESS);
        }

        if (option_char == 0){
            LongOptionMap::iterator it = lMap_.find(opts_[option_index].name);
            if (it != lMap_.end()){
                if (it->second->name_ == "help") {
                    printHelp(argv[0]);
                    exit(EXIT_SUCCESS);
                }
                it->second->assign(optarg);
                //std::cout << "long: "<< it->second->name_ << " = "  <<  GIMLI::__GIMLI_DEBUG__<< std::endl;
            } else {
                //std::cout << "\t not in longmap" << std::endl;
                printHelp(argv[0]);
                exit(EXIT_SUCCESS);
            }
        } else {
            CharOptionMap::iterator it = cMap_.find((char)option_char);
            if (it != cMap_.end()){
            //std::cout << "short: "<< it->second->name_ << " = "<<  0<< std::endl;
               it->second->assign(optarg);
            } else {
                //std::cout << "\t not in charmap" << std::endl;
                printHelp(argv[0]);
                exit(EXIT_SUCCESS);
            }
        }
    }

    if (showVersion_){
        std::cout << versionStr() << std::endl;
        exit(EXIT_SUCCESS);
    }
    GIMLI::setDebug(debug_);

    std::string lastString;
    for (int index = optind; index < argc; index++){
        lastString += std::string(argv[index]);
        if (index < argc -1) lastString += " " ;
    }
    if (lastOption_){
        if (lastString.empty()){
            printHelp(argv[0]);
            exit(EXIT_SUCCESS);
        } else {
      //      std::cout << "lastoption"<< lastOption_->name_ << " = "<< std::endl;
            lastOption_->assign(&lastString[0]);
        }
    } else {
        std::cout << "Non-option argument: " << lastString << std::endl;
    }
}

void OptionMap::buildLongOptions(){
    if (opts_) delete[] opts_;

    opts_ = new option[lMap_.size() + 1];
    int count = 0;

    for (LongOptionMap::iterator it = lMap_.begin(); it != lMap_.end(); it++){
        opts_[count].name     = it->first.c_str();
        opts_[count].has_arg  = it->second->hasArg_;
        opts_[count].val      = 0; //(*it).key();
        opts_[count].flag     = 0;
        count ++;
    }

    opts_[count].name     = 0;
    opts_[count].has_arg  = 0;
    opts_[count].val      = 0;
    opts_[count].flag     = 0;
}

void OptionMap::printHelp(const std::string & main){
    std::cout << "Usage: " << main << " [options] " << lastArgString_ << std::endl;
    std::cout << "Description: " << descriptionString_ << std::endl;
    std::cout << "Options:" << std::endl;
    std::cout << "  -shortcut [--longname] type (defaultvalue) " << "\t\t\t: General option description" << std::endl<< std::endl;

    for (std::list < OptionBase * >::iterator it = options_.begin(), itmax = options_.end(); it != itmax; it++) {
        std::cout << "  -" << (char)(*it)->key_;
        if ((*it)->name_.size() > 0){
            std::cout << " [--" << (*it)->name_  << "]";
        }
        if ((*it)->hasArg_ == required_argument ) {
            std::cout << " " << (*it)->typname() << " (" << (*it)->defautToString() << ")";
        } else if ((*it)->typname() == "int"){
            std::cout << " incremental (" << (*it)->defautToString() << ")";
        }
        std::cout << "\t\t\t: " << (*it)->help_ << std::endl;
    }
    std::cout << std::endl;
}

} // namespace GIMLI
