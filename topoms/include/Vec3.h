/*
Copyright (c) 2018, Lawrence Livermore National Security, LLC.
Produced at the Lawrence Livermore National Laboratory.
Written by Harsh Bhatia (hbhatia@llnl.gov) and Attila G Gyulassy
(jediati@sci.utah.edu).
LLNL-CODE-745278. All rights reserved.

This file is part of TopoMS, Version 1.0. For details, see
https://github.com/LLNL/TopoMS. Please also read this link – Additional BSD
Notice.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

• Redistributions of source code must retain the above copyright notice, this
list of conditions and the disclaimer below.
• Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the disclaimer (as noted below) in the
documentation and/or other materials provided with the distribution.
• Neither the name of the LLNS/LLNL nor the names of its contributors may be
used to endorse or promote products derived from this software without specific
prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY, LLC, THE
U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

Additional BSD Notice

1. This notice is required to be provided under our contract with the U.S.
Department of Energy (DOE). This work was produced at Lawrence Livermore
National Laboratory under Contract No.  DE-AC52-07NA27344 with the DOE.

2. Neither the United States Government nor Lawrence Livermore National
Security, LLC nor any of their employees, makes any warranty, express or
implied, or assumes any liability or responsibility for the accuracy,
completeness, or usefulness of any information, apparatus, product, or process
disclosed, or represents that its use would not infringe privately-owned
rights.

3. Also, reference herein to any specific commercial products, process, or
services by trade name, trademark, manufacturer or otherwise does not
necessarily constitute or imply its endorsement, recommendation, or favoring by
the United States Government or Lawrence Livermore National Security, LLC. The
views and opinions of authors expressed herein do not necessarily state or
reflect those of the United States Government or Lawrence Livermore National
Security, LLC, and shall not be used for advertising or product endorsement
purposes.
*/

/**
 *  @file    Vec3.h
 *  @author  Harsh Bhatia (hbhatia@llnl.gov)
 *  @date    10/01/2017
 *
 *  @brief This class handles 3-dimensional vector objects
 *
 *  @section DESCRIPTION
 *
 *  This class handles 3-dimensional vector objects
 *
 */

#ifndef _VEC3_H_
#define _VEC3_H_

#include <cmath>
#include <cstdio>

template <typename T>
class Vec3 {

private:
    T v[3];

public:
    Vec3() {                v[0] = v[1] = v[2] = T(0);      }
    Vec3(T x, T y, T z) {   v[0] = x; v[1] = y; v[2] = z;   }

    T& operator[](int i) {              return v[i]; }
    const T& operator[](int i) const {  return v[i]; }

    inline T magnSq() const {    return (v[0] * v[0]) + (v[1] * v[1]) + (v[2] * v[2]);   }
    inline T magn() const {      return sqrt(magnSq());                                   }

    void operator=(const Vec3& other){
        this->v[0] = other.v[0];
        this->v[1] = other.v[1];
        this->v[2] = other.v[2];
    }
    void operator+=(const Vec3& other) {
        this->v[0] += other.v[0];
        this->v[1] += other.v[1];
        this->v[2] += other.v[2];
    }
    void operator-=(const Vec3& other) {
        this->v[0] -= other.v[0];
        this->v[1] -= other.v[1];
        this->v[2] -= other.v[2];
    }
    Vec3 operator+(const Vec3  &a) const {
        Vec3 res(*this);
        res += a;
        return res;
    }
    Vec3 operator-(const Vec3  &a) const {
        Vec3 res(*this);
        res -= a;
        return res;
    }
    Vec3 operator+(const T  val) const {
        Vec3 res(*this);
        res.v[0] += val; res.v[1] += val; res.v[2] += val;
        return res;
    }
    Vec3 operator*(const T  val) const {
        Vec3 res(this->v[0], this->v[1], this->v[2]);// = (*this);
        res.v[0] *= val; res.v[1] *= val; res.v[2] *= val;
        return res;
    }
    bool operator==(const Vec3& other) const {
        return memcmp(this, &other, sizeof(Vec3)) == 0;
    }

    // THIS GIVES REMAINDER, NOT MODULO
    Vec3 operator%(const Vec3& other) const{
        Vec3 res;
        res.v[0] = this->v[0] % other.v[0];
        res.v[1] = this->v[1] % other.v[1];
        res.v[2] = this->v[2] % other.v[2];
        return res;
    }

    operator Vec3<int>() const {    return Vec3<int>(this->v[0], this->v[1], this->v[2]); }
    operator Vec3<double>() const { return Vec3<double>(this->v[0], this->v[1], this->v[2]); }

    // linear interpolation betwee 2 vectors: returns the vector a * (1-t) + b * t
    static Vec3<double> lerp(const Vec3<double>& a, const Vec3<double>& b, const double t)  {
        return (a * (1.0 - t)) + (b * (t));
    }

    void print_vi() const {
        printf("(%d, %d, %d)\n", v[0], v[1], v[2]);
    }
    void print_vf() const {
        printf("(%f, %f, %f)\n", v[0], v[1], v[2]);
    }

    // this is a hack that rounds towards -infinity for "small" negative values
    Vec3<int> int_floor() const {
        Vec3 dres = *this;
        dres = dres + 1000.0;
        Vec3<int> ires = dres;
        ires = ires + -1000;
        return ires;
    }
};


typedef Vec3<int> Vec3i;
typedef Vec3<bool> Vec3b;
typedef Vec3<float> Vec3f;
typedef Vec3<double> Vec3d;

#endif
