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

#ifndef GIMLI_LINE__H
#define GIMLI_LINE__H

#include "gimli.h"
#include "pos.h"

namespace GIMLI{

//! A line
/*! A straight line segment in the geometrically point of view. */

class DLLEXPORT Line {
public:
    /*! Default constructor. Constructs an invalid straight line. */
    Line();

    /*! Constructor, construct a line between p0 and the orign [0.0, 0.0, 0.0 ]. */
    Line(const RVector3 & p0);

    /*! Constructor construct a line betweem p0, and p1.*/
    Line(const RVector3 & p0, const RVector3 & p1);

    /*! Copycontructor */
    Line(const Line & line);

    /*! Default destructor. Destroys the object.*/
    ~Line();

    /*! Assignment operator*/
    Line & operator = (const Line & line);

    /*! Functor to the line function. ret = _p0 + (_p1 - _p0) * t; */
    inline const RVector3 operator()(double t) const { return this->at(t); }

    /*! Equal_to operator */
    inline bool operator == (const Line & line) const { return this->compare(line); }

    /*! Not_equal_to operator */
    inline bool operator != (const Line & line) const { return !(*this == line); }

    /*! Check validility of the straight line. False if [ p0-p1] < tol. */
    bool checkValidity(double tol=TOLERANCE);

    /*! Return the validity of this line.*/
    inline bool valid() const { return valid_; }

    /*! Compare two lines with a given tolerance. Check if both ref points from line touch this line..
    distance(this, line.p1)< tol && distance(this, line.p2)< tol */
    bool compare(const Line & line, double tol=TOLERANCE) const;

    /*! Return const reference to ref. point p0 */
    const RVector3 & p0() const { return _p0; }

    /*! Return const reference to ref. point p1 */
    const RVector3 & p1() const { return _p1; }

    /*! Return true if this line segment intersect with a the ray given by
     start and direction dir. Fill pos with the intersection position.
     Returns also true if the ray is parallel and the distance between
     Line and ray is smaller than the tolerance but pos is not valid.*/
    bool intersectRay(const RVector3 & start, const RVector3 & dir,
                      RVector3 & pos, double tol=TOLERANCE) const;

    /*! Return insection position of this line segment with a the ray given by
     start and direction dir. Pos is invalid for no intersection.*/
    RVector3 intersect(const RVector3 & start, const RVector3 & dir, double tol=TOLERANCE) const;

    /*! Returns the intersection point between line and this line.
        Returns an invalid RVector3 of both do not touch, are equal or parallel.
        Return is invalid if they don't intersect.
    */
    RVector3 intersect(const Line & line, double tol=TOLERANCE) const;

    /*! Returns line parameter t that represents the nearest position to the
     *line for the given point p. Feed it to at() and you get the nearest point. */
    double nearest(const RVector3 & p) const;

    /*! Apply the line as function. ret = _p0 + (_p1 - _p0) * t*/
    inline RVector3 at(double t) const { return _p0 + (_p1 - _p0) * t; }

    /*!  Return the distance between this line and pos,
    https://mathworld.wolfram.com/Point-LineDistance3-Dimensional.html;
        | x2-x1 cross x1-x0 | / | x2 -x1 | */
    double distance(const RVector3 & pos) const;

    /*! Return position at line.
        solves: pos = _p0 + (_p1 - _p0) * ret.
        ret = (pos - _p0) / (_p1 - _p0)
        Throws an exception if pos is not on Line. */
    double t(const RVector3 & pos, double tol=TOLERANCE) const ;

    /*! Return true if pos touches the line. pFunIdx gives an identifier which shows the dependency between this Line and the RVector3 pos. Possible return values are: \n
    -1 -- this straight line don't touch the Position pos \n
    1 -- this straight line touch the Position pos and pos lies before _p0 \n
    2 -- this straight line touch the Position pos and pos lies at the same position like pos _p0\n
    3 -- this straight line touch the Position pos and pos lies within the definition positions _p0 and _p1\n
    4 -- this straight line touch the Position pos and pos lies at the same position like pos _p1\n
    5 -- this straight line touch the Position pos and pos lies behind _p1\n  */
    bool touch1(const RVector3 & pos, int & pFunIdx, double tol=TOLERANCE) const ;

    /*! Return true if pos touches the line. */
    bool touch(const RVector3 & pos, double tol=TOLERANCE) const ;

    /*! Return the length if the line.*/
    double length() const { return _p0.dist(_p1); }

protected:

    /*! internal: copy content of line into this line */
    void copy_(const Line & line);

    RVector3 _p0;
    RVector3 _p1;
    bool valid_;
};

std::ostream & operator << (std::ostream & str, const Line & l);

} // namespace GIMLI;

#endif // GIMLI_LINE__H
