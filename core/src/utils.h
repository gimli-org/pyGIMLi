
/******************************************************************************
 *   Copyright (C) 2006-2024 by the GIMLi development team                    *
 *   Carsten Rücker carsten@gimli.org                                         *
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
#pragma once

#include "gimli.h"

#include <fstream>

namespace GIMLI{

class ByteBuffer;

/*! Simple hexdump.
borrowed from from https://stackoverflow.com/questions/29242/off-the-shelf-c-hex-dump-code
*/
DLLEXPORT void hexdump(void *ptr, Index buflen);


class ByteBuffer{
    public:
    ByteBuffer(){

    }
    void push_back(unsigned char c){
        _buf.push_back(c);
    }
    inline char * begin() const { return (char *)&_buf[0]; }
    inline char * end() const { return begin() + size(); }
    inline Index size() const { return _buf.size(); }

    #ifndef PYGIMLI_CAST
    std::vector < uint8 > & buf() { return _buf; }
    #endif
    protected:
    std::vector < uint8 > _buf;
};

struct OutBuffer : public std::streambuf{
    ByteBuffer outbuffer;
    virtual int_type overflow (int_type c) {
        outbuffer.push_back(c);
        return c;
    }
};

struct InBuffer : public std::streambuf{

    InBuffer(const ByteBuffer & s) : std::streambuf() {
        char * begin = s.begin();
        char * end = s.end();
        this->setg(begin, begin, end);
    }

#ifndef PYGIMLI_CAST
    pos_type seekoff(off_type off, std::ios_base::seekdir dir,
                    std::ios_base::openmode which = std::ios_base::in) override{
        if (dir == std::ios_base::cur){
            gbump(off);
        } else if (dir == std::ios_base::end){
            setg(eback(), egptr() + off, egptr());
        } else if (dir == std::ios_base::beg){
            setg(eback(), eback() + off, egptr());
        } return gptr() - eback();
    }
#endif
    pos_type seekpos(pos_type sp, std::ios_base::openmode which) override {
        return seekoff(sp - pos_type(off_type(0)), std::ios_base::beg, which);
    }
};


/*! Serialize an obj into a ByteBuffer buff.
The obj need have a method `writeToStream` */
template < class T > ByteBuffer serialize(const T & obj) {
    OutBuffer buf;
    std::ostream out(&buf);
    if (!out){
        throw std::system_error(errno, std::system_category(),
                                "failed to serialize.");
    }

    obj.writeToStream(out);
    return buf.outbuffer;

}
/*! Deserialize ByteBuffer buff into obj.
That need have a method `readFromStream` */
template < class T > void deserialize(const ByteBuffer & buf, T & obj) {
    InBuffer ib(buf);
    std::istream in(&ib);
    if (!in){
        throw std::system_error(errno, std::system_category(),
                                "failed to deserialize.");
    }
    obj.readFromStream(in);
}

/*! Serialize an amount of values to a stream */
template < class ValueType > void writeStream(std::ostream & out,
                                              const ValueType & v, int count=1){
    out.write((const char*)&v, sizeof(ValueType)*count);
}

/*! Deserialize an amount of values to from a stream into a reserved
array. */
template < class ValueType > void readStream(std::istream & in,
                                             ValueType & v, int count=1){
    in.read((char*)&v, sizeof(ValueType)*count);
}

/*! Save object from a binary file.
The object need have a method `writeToStream`.
Returns the name of the filename. */
template < class T > std::string saveBinary(const std::string & fbody,
                                            const T & obj,
                                            const std::string & suffix="") {
    std::string fileName(fbody.substr(0, fbody.rfind(suffix)) + suffix);

    std::ofstream out(fileName, std::ios::binary | std::ios::out);
    if (!out){
        throw std::system_error(errno, std::system_category(),
                                "failed to open: " + fileName);
    }
    obj.writeToStream(out);
    out.flush();
    out.close();
    return fileName;
}

/*! Load object from a binary file.
The object need have a method `readFromStream` */
template < class T > void loadBinary(const std::string & fbody,
                                     T & obj,
                                     const std::string & suffix="") {
    std::string fileName(fbody.substr(0, fbody.rfind(suffix)) + suffix);
    std::ifstream in(fileName, std::ios::binary | std::ios::in);
    if (!in){
        throw std::system_error(errno, std::system_category(),
                                "failed to open: " + fileName);
    }

    obj.readFromStream(in);
    in.close();
}

}