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
 *  @file    MSCBond.h
 *  @author  Harsh Bhatia (hbhatia@llnl.gov)
 *  @date    10/01/2018
 *
 *  @brief This file provides a class for analyzing a bond.
 *
 *  @section DESCRIPTION
 *
 *  This file provides a class for analyzing a bond.
 *
 */

#ifndef _MSCBOND_H_
#define _MSCBOND_H_

#include <vector>
#include <utility>

#include "basic_types.h"
#include "vectors.h"

namespace MS {
class SystemInfo;
}

/**
  *  @brief This struct provides a colleciton of utilities to handle MSC bonds
  */
struct MSCBond {

    // saddle node id and MSC coordinates
    INT_TYPE saddle;
    MSC::Vec3d scoords;

    // extrema ids and MSC coordiantes
    std::vector<INT_TYPE> extrema;
    std::vector<MSC::Vec3d> ecoords;

    // the corresponding atoms (starting with 1)
    std::vector<size_t> atomIds;

    // geometric representation of paths
    std::vector<std::vector<MSC::Vec3d>> paths;

    // parameterization of the path
    std::vector<std::pair<float, MSC::Vec3d>> parameterization;

    // integrated charge and area (at the saddle)
    double ichg, iarea;

    // integrated charge and area along the path (indexed into the parameterization)
    std::vector<double> bichg, biarea;

    // charge at the bond path
    std::vector<double> bchg;

    // -------------------------------------------------------------------------
    static void fix_periodic(MSC::Vec3d &p, const MSC::Vec3d &orig, const size_t dims[3]);

    // -------------------------------------------------------------------------
    MSCBond(INT_TYPE _saddle=-1) : saddle(_saddle), ichg(-1), iarea(-1) {}

    void print() const;
    bool check2() const;
    void parameterize(const MS::SystemInfo &metadata);
    void get_points(MSC::Vec3d &origin, std::vector<MSC::Vec3d> &nbrs, const size_t dims[], float p) const;
    void get_points_idx(MSC::Vec3d &origin, std::vector<MSC::Vec3d> &nbrs, const size_t dims[], int pidx) const;
};

#endif
