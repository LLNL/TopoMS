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

#ifndef REGULAR_GRID_H
#define REGULAR_GRID_H

#include "basic_types.h"
#include "vectors.h"


namespace MSC {
    // just a class to do some index stuff on a grid
    // e.g. get the 6 neighbors of a point
    // or get the 8 vertices surrounding a sample location.
    // x moves fastest then y then z

    class RegularGrid {
    protected:



        const Vec3l m_xyz;
        const Vec3b m_periodic;
        static const Vec3l kNeighborOffsets6[6];
    public:

        RegularGrid(Vec3l xyz, Vec3b p) : m_xyz(xyz), m_periodic(p) {
            printf(" -- Created RegularGrid [%d %d %d] with periodicity [%d %d %d]\n",
                   m_xyz[0], m_xyz[1], m_xyz[2], m_periodic[0], m_periodic[1], m_periodic[2]);
        }

        inline const Vec3l& XYZ() const { return m_xyz; }
        inline const Vec3b& Periodic() const { return m_periodic; }
        inline INDEX_TYPE Index3d(const Vec3l& v) const {
            return v[0] + v[1] * m_xyz[0] + v[2] * m_xyz[0] * m_xyz[1];
        }
        inline Vec3l XYZ3d(INDEX_TYPE id) const {
            Vec3l res(id % m_xyz[0], (id / m_xyz[0]) % m_xyz[1], id / (m_xyz[0] * m_xyz[1]));
            return res;
        }
        static inline INDEX_TYPE PositiveModulo(INDEX_TYPE a, INDEX_TYPE b) {
            return (a % b + b) % b;
        }
        static inline Vec3l PositiveModulo(const Vec3l& a, const Vec3l& b) {
            return (a % b + b) % b;
        }
        // return Inbounds version of the vertex - in case of negative indices
        inline Vec3l Inbounds(const Vec3l& v) const {
            Vec3l res = v;
            for (int i = 0; i < 3; i++) {
                if (m_periodic[i]) {
                    res[i] = PositiveModulo(v[i], m_xyz[i]);
                }
                else {
                    if (res[i] < 0) res[i] = 0;
                    if (res[i] >= m_xyz[i]) res[i] = m_xyz[i] - 1;
                }
            }
            return res;
        }



        // return Inbounds version of the vector - find the base point - towards -infinity,
        inline Vec3d Inbounds(const Vec3d& v) const {


            Vec3d res = v;
            for (int i = 0; i < 3; i++) {
                if (m_periodic[i]) {
                    if (res[i] < 0) res[i] += m_xyz[i];
                    if (res[i] >= m_xyz[i]) res[i] -= m_xyz[i];
                }
                else {
                    if (res[i] < 0) res[i] = 0;
                    if (res[i] >= m_xyz[i]) res[i] = m_xyz[i] - 1;
                }
            }
            return res;

            //Vec3l fl = v.IntFloor();
            //Vec3d diff = fl;
            //diff = v - fl;
            //Vec3l res = Inbounds(fl);
            //Vec3d resd = res;
            //resd = resd + diff;
            //return resd;
        }

        INDEX_TYPE NumElements() const {
            return m_xyz[0] * m_xyz[1] * m_xyz[2];
        }

        static void PrintVector(const Vec3d& v) {
            printf("(%f, %f, %f)\n", v[0], v[1], v[2]);
        }
        static void PrintVector(const Vec3l& v) {
            printf("(%d, %d, %d)\n", v[0], v[1], v[2]);
        }

        int GatherExistingNeighbors6(const Vec3l& s, Vec3l* results) const {
            int numsofar = 0;
            if (m_periodic[0]) {
                results[numsofar++] = PositiveModulo(s + kNeighborOffsets6[0], m_xyz);
                results[numsofar++] = PositiveModulo(s + kNeighborOffsets6[1], m_xyz);
            }
            else {
                if (s[0] < m_xyz[0] - 1) results[numsofar++] = s + kNeighborOffsets6[0];
                if (s[0] > 0) results[numsofar++] = s + kNeighborOffsets6[1];
            }
            if (m_periodic[1]) {
                results[numsofar++] = PositiveModulo(s + kNeighborOffsets6[2], m_xyz);
                results[numsofar++] = PositiveModulo(s + kNeighborOffsets6[3], m_xyz);
            }
            else {
                if (s[1] < m_xyz[1] - 1) results[numsofar++] = s + kNeighborOffsets6[2];
                if (s[1] > 0) results[numsofar++] = s + kNeighborOffsets6[3];
            }
            if (m_periodic[2]) {
                results[numsofar++] = PositiveModulo(s + kNeighborOffsets6[4], m_xyz);
                results[numsofar++] = PositiveModulo(s + kNeighborOffsets6[5], m_xyz);
            }
            else {
                if (s[2] < m_xyz[2] - 1) results[numsofar++] = s + kNeighborOffsets6[4];
                if (s[2] > 0) results[numsofar++] = s + kNeighborOffsets6[5];
            }
            return numsofar;
        }

        // get neighbors along chosen axis, where kernelsize is number in each direction
        // results must be big enough to contain, = 2 * kernelsize + 1 (we include original point)
        // for non-m_periodic grids, the boundary point is replicated, e.g. the 2-kernel around point 0 is 0 0 0 1 2
        void Gather1DNeighborhood(const Vec3l& s, int axis, int kernelsize, Vec3l* results) const {
            int current = 0;
            if (m_periodic[axis]) {
                for (int i = -kernelsize; i <= kernelsize; i++) {
                    Vec3l offset = s;
                    offset[axis] = PositiveModulo(offset[axis] + i, m_xyz[axis]);
                    results[i + kernelsize] = offset;
                }
            }
            else {
                for (int i = kernelsize; i > 0; i--) {
                    Vec3l offset = s;
                    offset[axis] =offset[axis] - i ;
                    if (offset[axis] < 0) {
                        offset[axis] = 0;
                    }
                    results[kernelsize - i] = offset;
                }
                results[kernelsize] = s;
                for (int i = 0; i < kernelsize; i++) {
                    Vec3l offset = s;
                    offset[axis] = i + 1 + offset[axis];
                    if (offset[axis] >= m_xyz[axis]) {
                        offset[axis] = m_xyz[axis] - 1;
                    }
                    results[kernelsize + i + 1] = offset;
                }
            }
        }


        int GatherSurrounding(const Vec3d& s, Vec3l* results) const {
            Vec3l basep = s.IntFloor();
            return GatherSurrounding(basep, results);
        }

        // this gives minimum manhattan distance to a boundary
        int DistToBoundary(const Vec3l& a) const {
            Vec3l b = (m_xyz - 1) - a; // b stores distance to extents boundary
            // get closest distance to 0,0,0 planes
            INT_TYPE min01a = (a[0] < a[1] ? a[0] : a[1]);
            INT_TYPE min12a = (a[1] < a[2] ? a[1] : a[2]);
            // get closest distance to extents boundary
            INT_TYPE min01b = (b[0] < b[1] ? b[0] : b[1]);
            INT_TYPE min12b = (b[1] < b[2] ? b[1] : b[2]);
            // are we closer to which planes?
            INT_TYPE mina = (min01a < min12a ? min01a : min12a);
            INT_TYPE minb = (min01b < min12b ? min01b : min12b);
            return (mina < minb ? mina : minb);
        }

        // this function should only be used with extreme care!
        // all bets are off if the queried point lies withing 1 cell of the boundary
        int GatherSurroundingNoBoundaryCheck(const Vec3d& s, Vec3l* results) const {
            Vec3l basep = s;
            return GatherSurroundingNoBoundaryCheck(basep, results);
        }
        int GatherSurroundingNoBoundaryCheck(const Vec3l& basep, Vec3l* results) const {

            //printf("GatherSurrounding::input = "); PrintVector(s);
            //PrintVector(basep);

            int xvecs[2];
            xvecs[0] = basep[0];
            xvecs[1] = basep[0] + 1;
            int yvecs[2];
            yvecs[0] = basep[1];
            yvecs[1] = basep[1] + 1;

            int zvecs[2];
            zvecs[0] = basep[2];
            zvecs[1] = basep[2] + 1;

            int numsofar = 0;
            for (int k = 0; k < 2; k++)
                for (int j = 0; j < 2; j++)
                    for (int i = 0; i < 2; i++) {
                results[numsofar++] = Vec3l(xvecs[i], yvecs[j], zvecs[k]);
                    }
            return 8;
        }


        // THIS ASSUMES that we are REASONABLE, i.e. basep is not out of bounds
        INT_TYPE GatherSurrounding(const Vec3l& basep, Vec3l* results) const {

            //printf("GatherSurrounding::input = "); PrintVector(s);
            //PrintVector(basep);
            for (int i = 0; i < 3; i++) {
                if (basep[i] < 0 || basep[i] > m_xyz[i]-1) {
                    printf("basep="); basep.PrintInt();;
                }
            }

            INT_TYPE xvecs[2];
            if (m_periodic[0]) {
                xvecs[0] = basep[0];
                xvecs[1] = (basep[0] == m_xyz[0]-1 ? 0 : basep[0]+1 ); // logic here: if we are on periodic boundary, wrap around, else increment
            }
            else {
                xvecs[0] = basep[0]; xvecs[1] = 1 + basep[0];
                if (basep[0] < 0) xvecs[0] = 0;
                if (basep[0] >= m_xyz[0] - 1) { xvecs[0] = xvecs[1] = m_xyz[0] - 1; }
                else if (basep[0] >= m_xyz[0]) { xvecs[1] = m_xyz[0] - 1; }
            }
            INT_TYPE yvecs[2];
            if (m_periodic[1]) {
                yvecs[0] = basep[1];
                yvecs[1] = (basep[1] == m_xyz[1]-1 ? 0 : basep[1] + 1);
            }
            else {
                yvecs[0] = basep[1]; yvecs[1] = 1 + basep[1];
                if (basep[1] < 0) yvecs[0] = 0;
                if (basep[1] >= m_xyz[1] - 1) { yvecs[0] = yvecs[1] = m_xyz[1] - 1; }
                else if (basep[1] >= m_xyz[1]) { yvecs[1] = m_xyz[1] - 1; }
            }

            INT_TYPE zvecs[2];
            if (m_periodic[2]) {
                zvecs[0] = basep[2];
                zvecs[1] = (basep[2] == m_xyz[2]-1 ? 0 : basep[2] + 1);
            }
            else {
                zvecs[0] = basep[2]; zvecs[1] = 1 + basep[2];
                if (basep[2] < 0) zvecs[0] = 0;
                if (basep[2] >= m_xyz[2] - 1) { zvecs[0] = zvecs[1] = m_xyz[2] - 1; }
                else if (basep[2] >= m_xyz[2]) { zvecs[1] = m_xyz[2] - 1; }
            }

            INT_TYPE numsofar = 0;
            for (int k = 0; k < 2; k++)
                for (int j = 0; j < 2; j++)
                    for (int i = 0; i < 2; i++) {
                results[numsofar++] = Vec3l(xvecs[i], yvecs[j], zvecs[k]);
                    }
            return 8;
        }
        ////int GatherSurrounding(const Vec3d& const s, INDEX_TYPE* results, double* factors) {
        ////	Vec3l surrounding[8];
        ////	GatherSurrounding(s, surrounding);

        ////	return 8;
        ////}
    };


}

#endif
