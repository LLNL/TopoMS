/*
 * Copyright (c) 2017 University of Utah 
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 * 3. Neither the name of the copyright holder nor the names of its
 *    contributors may be used to endorse or promote products derived from
 *    this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */

#ifndef TOPOLOGICAL_REGULAR_GRID_RESTRICTED_H
#define TOPOLOGICAL_REGULAR_GRID_RESTRICTED_H


#include "topological_regular_grid.h"
#include "labeling.h"


namespace MSC {

class TopologicalRegularGridRestricted : virtual public TopologicalRegularGrid {

protected:
    const DenseLabeling<char> *m_restriction;

public:
    TopologicalRegularGridRestricted(RegularGrid* base_grid) :
        TopologicalRegularGrid(base_grid), m_restriction(nullptr) {
        //printf(" -- Created TopologicalRegularGridRestricted \n");
    }

    void set_restriction(const DenseLabeling<char> *restriction) {
        m_restriction = restriction;
    }
    char restrictionLabel(INDEX_TYPE cellid) const {

        if (m_restriction == nullptr) {
            std::cerr << " restriction not set for TopologicalRegularGridRestricted!\n";
        }
        return m_restriction->GetLabel(cellid);
    }

    BOUNDARY_TYPE boundaryValue(INDEX_TYPE cellid) const {

        //printf("TopologicalRegularGridRestricted::boundaryValue()\n");

        // TopologicalRegularGrid::boundaryValue returns 0,1,2,3
        // m_restriction has labels 0 or 1.. so this function returns (0,1,2,3) or (4,5,6,7)

        // Harsh added the following conditional on 07.28.2018 to make this class
        // also work for unrestricted meshes
        BOUNDARY_TYPE bval = TopologicalRegularGrid::boundaryValue (cellid);
        if (m_restriction == nullptr) { return bval;                                                        }
        else {                          return bval + (m_restriction->GetLabel (cellid) * (maxDim()+1));    }
    }
};
}   // end of namespace
#endif
