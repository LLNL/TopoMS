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
 *  @file    MSCBond.cpp
 *  @author  Harsh Bhatia (hbhatia@llnl.gov)
 *  @date    05/12/2017
 *
 *  @brief This file provides a collection of utilities to handle MSC bonds
 *
 *  @section DESCRIPTION
 *
 * This file provides a collection of utilities to handle MSC bonds
 *
 */

#include <map>
#include <utility>
#include <algorithm>

#include "MSCBond.h"
#include "MultilinearInterpolator.h"

// TODO: this should go in utils
void MSCBond::fix_periodic(MSC::Vec3d &p, const MSC::Vec3d &orig, const size_t dims[3]) {

    for(uint8_t d = 0; d < 3; d++) {

        if (p[d] - orig[d] > size_t(0.5*dims[d]))
            p[d] -= dims[d];

        else if (orig[d] - p[d] > size_t(0.5*dims[d]))
            p[d] += dims[d];
    }
}

void MSCBond::get_points_idx(MSC::Vec3d &origin, std::vector<MSC::Vec3d> &nbrs, const size_t dims[3], int pidx) const {

    if (pidx == -1) {
        origin = scoords;    nbrs.push_back(ecoords[0]);    nbrs.push_back(ecoords[1]);
    }
    else {
        origin = parameterization[pidx].second;
        if (pidx != 0) {
            nbrs.push_back(parameterization[pidx-1].second);
        }
        if (pidx != parameterization.size()-1) {
            nbrs.push_back(parameterization[pidx+1].second);
        }

        // if either was not found
        if (nbrs.size() != 2) {
            nbrs.push_back(origin - (nbrs.front()-origin));
        }
    }

    // fix periodic
    for(uint8_t i = 0; i < nbrs.size(); i++) {
        fix_periodic(nbrs[i], origin, dims);
    }
}

void MSCBond::get_points(MSC::Vec3d &origin, std::vector<MSC::Vec3d> &nbrs, const size_t dims[3], float p) const {

    // TODO: float comparison should move to utils
    if (fabs(p) < 0.000001) {
        return get_points(origin, nbrs, dims, -1);
    }

    if (p < 0)          p = -1.0 * p * this->parameterization.front().first;
    else if (p > 0)     p =        p * this->parameterization.back().first;


    std::multimap<float, unsigned int> dmap;
    for(unsigned int i = 0; i < parameterization.size(); i++) {
        dmap.insert(std::make_pair(fabs(p-parameterization[i].first), i));
    }

    return get_points_idx(origin, nbrs, dims, int(dmap.begin()->second));
}

void MSCBond::parameterize(const size_t dims[3]) {

    if (!this->check2()){
        return;
    }

    parameterization.clear();
    for(unsigned int k = 0; k < 2; k++) {

        const std::vector<MSC::Vec3d>& path = this->paths[k];
        for(unsigned int i = 0; i < path.size(); i++) {

            MSC::Vec3d pp = path[i];

            // fix periodic
            fix_periodic(pp, scoords, dims);

            float p = (pp-scoords).Mag();
            if (k == 0) p *= -1;

            if(parameterization.size() > 0) {
                // TODO: float comparison should move to utils
                if (fabs(p-parameterization.back().first) < 0.00001)
                    continue;
            }
            parameterization.push_back(std::pair<float, MSC::Vec3d>(p, path[i]));
        }

        if (k == 0)
            std::reverse(parameterization.begin(), parameterization.end());
    }
    return;
    for(int k = 0; k < parameterization.size(); k++)
        printf(" %f, (%f %f %f)\n", parameterization[k].first, parameterization[k].second[0], parameterization[k].second[1], parameterization[k].second[2]);
}

void MSCBond::study_value(const double *func, const size_t dims[], std::vector<std::pair<float, float>> &vals) {

    vals.clear();
    vals.resize(this->parameterization.size());
    for(size_t i = 0; i < this->parameterization.size(); i++) {

        float x = this->parameterization[i].first;

        const MSC::Vec3d &p = this->parameterization[i].second;
        double pos[3] = {p[0], p[1], p[2]};

        double y = MultilinearInterpolator::trilinear_interpolation(pos, func, dims);
        vals[i] = std::make_pair(x,y);
    }
}
