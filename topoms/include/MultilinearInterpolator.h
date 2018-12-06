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
 *  @file    MultilinearInterpolator.h
 *  @author  Harsh Bhatia (hbhatia@llnl.gov)
 *  @date    05/12/2017
 *
 *  @brief This file provides a templated interpolator
 *
 *  @section DESCRIPTION
 *
 * This file provides a templated interpolator
 *
 */

#ifndef _MLINEARINTERP_H_
#define _MLINEARINTERP_H_

#include <cmath>
#include <vector>
#include <cstdlib>

#ifdef USE_VTK
    #include <vtkDataArray.h>
    #include <vtkDoubleArray.h>
    #include <vtkPointData.h>
    #include <vtkImageData.h>
#endif

// lerp based on lambda
class MultilinearInterpolator {
public:
    template <typename T, typename T2=double>
    static inline T lerp (T2 l, T f0, T f1) {
        return (1.0 - l)*f0 + l*f1;
    }
    template <typename T, typename T2=double>
    static inline T lerp2(T2 lx, T2 ly, T f00, T f01, T f10, T f11) {
        return lerp(ly,
                        lerp(lx, f00, f01),
                        lerp(lx, f10, f11));
    }
    template <typename T, typename T2=double>
    static inline T lerp3(T2 lx, T2 ly, T2 lz,
                                     T f000, T f010, T f100, T f110,
                                     T f001, T f011, T f101, T f111) {
        return lerp(lz,
                        lerp2(lx, ly, f000, f010, f100, f110),
                        lerp2(lx, ly, f001, f011, f101, f111));
    }

    template <typename T>
    static inline T trilinear_interpolation(const T pos[3], const T *data, const size_t dims[3]) {

        // get a corner pixel and gather the barycentric
        size_t o[3];
        double l[3];
        for(unsigned d = 0; d < 3; d++){
            l[d] = pos[d] - std::floor(pos[d]);
            o[d] = std::floor(pos[d]);
        }

        // gather the eight corners
        std::vector<T> vals;
        vals.reserve(8);
        for(uint8_t iz = 0; iz <= 1; iz++){
        for(uint8_t iy = 0; iy <= 1; iy++){
        for(uint8_t ix = 0; ix <= 1; ix++){

            // check perioudic boundary
            size_t z = iz+o[2];
            size_t y = iy+o[1];
            size_t x = ix+o[0];

            if (z >= dims[2])   z -= dims[2];
            if (y >= dims[1])   y -= dims[1];
            if (x >= dims[0])   x -= dims[0];

            size_t idx = z*dims[1]*dims[0] + y*dims[0] + x;
            vals.push_back(data[idx]);
        }
        }
        }

        return lerp3(l[0], l[1], l[2], vals[0], vals[1], vals[2], vals[3], vals[4], vals[5], vals[6], vals[7]);
    }

#ifdef USE_VTK
    template <typename T>
    static inline T trilinear_interpolation(const T pos[3], vtkImageData *volume) {

        int *dims = volume->GetDimensions();

        // get a corner pixel and gather the barycentric
        size_t o[3];
        double l[3];
        for(unsigned d = 0; d < 3; d++){
            l[d] = pos[d] - std::floor(pos[d]);
            o[d] = std::floor(pos[d]);
        }

        // gather the eight corners
        std::vector<T> vals;
        vals.reserve(8);
        for(uint8_t iz = 0; iz <= 1; iz++){
        for(uint8_t iy = 0; iy <= 1; iy++){
        for(uint8_t ix = 0; ix <= 1; ix++){

            // check perioudic boundary
            size_t z = iz+o[2];
            size_t y = iy+o[1];
            size_t x = ix+o[0];

            if (z >= dims[2])   z -= dims[2];
            if (y >= dims[1])   y -= dims[1];
            if (x >= dims[0])   x -= dims[0];

            size_t idx = z*dims[1]*dims[0] + y*dims[0] + x;
            vals.push_back(volume->GetPointData()->GetScalars()->GetComponent(idx, 0));
        }
        }
        }

        return lerp3(l[0], l[1], l[2], vals[0], vals[1], vals[2], vals[3], vals[4], vals[5], vals[6], vals[7]);
    }
#endif
};
#endif
