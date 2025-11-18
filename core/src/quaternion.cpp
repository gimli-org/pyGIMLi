/******************************************************************************
 *   Copyright (C) 2008-2024 by the GIMLi development team                    *
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

#include "quaternion.h"
#include "matrix.h"

namespace GIMLI{

std::ostream & operator << (std::ostream & str, const RQuaternion & q){
    str << "Quaternion(" << q[ 0 ] << ", " << q[ 1 ] << ", "
                            << q[ 2 ] << ", " << q[ 3 ] << ") ";
    return str;
}

RMatrix getRotation(const Pos & src, const Pos & dest){
// Based on Stan Melax's article in Game Programming Gems
    RQuaternion q;
    Pos v0(src);
    Pos v1(dest);
    if (v0.abs() < TOLERANCE || v1.abs() < TOLERANCE) {
        q = RQuaternion(1.0, 0.0, 0.0, 0.0);
    } else {
        v0.normalise();
        v1.normalise();

        double d = v0.dot(v1);

        if (::fabs((d - 1.0)) < TOLERANCE) { //** v1 == v2
            q = RQuaternion(1.0, 0.0, 0.0, 0.0);
        } else if (::fabs((d + 1.0)) < TOLERANCE) { //** v1 == -v2
            Pos a(RVector3(1.0, 0.0, 0.0).cross(v0));
            if (a.length() < TOLERANCE){
                a = RVector3(0.0, 1.0, 0.0).cross(v0);
            }
            a.normalise();
            q.createFromAxisAngle(a, PI);
        } else {
            double s = std::sqrt((1.0 + d) * 2.0);
            Pos c = v0.cross(v1) / s;
            q = RQuaternion(s * 0.5, c);
            q.normalise();
        }
    }

    RMatrix rot(3, 3);
    q.rotMatrix(rot);
    return rot;
}


} // namespace GIMLI
