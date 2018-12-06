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
 *  @date    10/01/2018
 *
 *  @brief This class handles 3x3 matrix objects
 *
 *  @section DESCRIPTION
 *
 *  This class handles 3x3 matrix objects
 *
 */

#ifndef _MAT3_H_
#define _MAT3_H_

#include <cmath>
#include <vector>
#include <iostream>

#include "Vec3.h"

template <typename T=double>
class Mat3 {

public:
    T v[3][3];

public:
    Mat3() {
        for(uint8_t i = 0; i < 3; i++)
        for(uint8_t j = 0; j < 3; j++)
            v[i][j] = T(0);
    }

    void eye() {
        for(uint8_t i = 0; i < 3; i++)
        for(uint8_t j = 0; j < 3; j++)
            v[i][j] = T(0);

        for(uint8_t i = 0; i < 3; i++)
            v[i][i] = T(1);
    }

    bool is_eye() const {

        Mat3<T> m;  m.eye();
        for(uint8_t i = 0; i < 3; i++)
        for(uint8_t j = 0; j < 3; j++){
            if (fabs(v[i][j] - m.v[i][j]) > 0.0000001)
                return false;
        }
        return true;
    }

    bool is_cuboid() const {

        for(uint8_t i = 0; i < 3; i++)
        for(uint8_t j = 0; j < 3; j++){
            if (i == j) {
                continue;
            }
            if (fabs(v[i][j]) > 0.0000001)
                return false;
        }
        return true;
    }
    T determinant() const {
        return   v[0][0]*(v[1][1]*v[2][2] - v[1][2]*v[2][1])
               - v[0][1]*(v[1][0]*v[2][2] - v[1][2]*v[2][0])
               + v[0][2]*(v[1][0]*v[2][1] - v[1][1]*v[2][1]);
    }
    Mat3 inverse() const {

        Mat3 inv;
        double det = determinant();
        if (fabs(det) < 0.000001) {
            std::cerr << "Cannot compute inverse for singular matrix!\n";
            inv.eye();
            return inv;
        }

        det = 1.0/det;
        inv.v[0][0] = det*(v[1][1]*v[2][2] - v[1][2]*v[2][1]);
        inv.v[0][1] = det*(v[0][2]*v[2][1] - v[0][1]*v[2][2]);
        inv.v[0][2] = det*(v[0][1]*v[1][2] - v[0][2]*v[1][1]);
        inv.v[1][0] = det*(v[1][2]*v[2][0] - v[1][0]*v[2][2]);
        inv.v[1][1] = det*(v[0][0]*v[2][2] - v[0][2]*v[2][0]);
        inv.v[1][2] = det*(v[0][2]*v[1][0] - v[0][0]*v[1][2]);
        inv.v[2][0] = det*(v[1][0]*v[2][1] - v[1][1]*v[2][0]);
        inv.v[2][1] = det*(v[0][1]*v[2][0] - v[0][0]*v[2][1]);
        inv.v[2][2] = det*(v[0][0]*v[1][1] - v[0][1]*v[1][0]);
        return inv;
    }

    void transform(const T in[3], T out[3]) const {
        for(uint8_t d = 0; d < 3; d++){
            out[d] = in[0]*v[0][d] + in[1]*v[1][d] + in[2]*v[2][d];
        }
    }

    Vec3<T> transform(const T in[3]) const {
        Vec3<T> out;
        for(uint8_t d = 0; d < 3; d++){
            out[d] = in[0]*v[0][d] + in[1]*v[1][d] + in[2]*v[2][d];
        }
        return out;
    }

    Vec3<T> transform(const Vec3<T> &in) const {
        Vec3<T> out;
        for(uint8_t d = 0; d < 3; d++){
            out[d] = in[0]*v[0][d] + in[1]*v[1][d] + in[2]*v[2][d];
        }
        return out;
    }

    std::vector<T> linearize(bool homogenous = false) const {

        const size_t sz = (homogenous) ? 4 : 3;
        std::vector<T> mat (sz*sz, T(0));

        for(uint8_t r = 0; r < 3; r++)
        for(uint8_t c = 0; c < 3; c++)
            mat[sz*c + r] = v[r][c];

        if (homogenous) {   mat[15] = T(1); }
        return mat;
    }
};

#endif
