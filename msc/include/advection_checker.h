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

#ifndef ADVECTION_CHECKER_H
#define ADVECTION_CHECKER_H

#include <unordered_set>

#include "basic_types.h"
#include "vectors.h"
#include "regular_grid.h"
#include "labeling.h"
#include "advection_events.h"



namespace MSC {
    class AdvectionChecker {
    public:

        virtual ADVECTION_EVENT CheckAndDoStuff(const Vec3d& old_point, const Vec3d& new_point, bool nearboundary){ return NONE; }
    };

    // no early termination during advection interior to a voxel
    class NoTermination : public AdvectionChecker {
    public:
        virtual ADVECTION_EVENT CheckAndDoStuff(const Vec3d& old_point, const Vec3d& new_point, bool nearboundary) {
            return NONE;
        }
    };


    // this variant
    class TerminateNearCrits : public AdvectionChecker {
    protected:
        DenseLabeling<INDEX_TYPE>* m_labeling;
        RegularGrid* m_grid;
    public:
        TerminateNearCrits(DenseLabeling<INDEX_TYPE>* labeling, RegularGrid*grid) : m_labeling(labeling), m_grid(grid) {}
        virtual ADVECTION_EVENT CheckAndDoStuff(const Vec3d& old_point, const Vec3d& new_point, bool nearboundary) {
            if (!nearboundary) {
                Vec3l closest_vertex = new_point + 0.5;
                INDEX_TYPE id = m_grid->Index3d(closest_vertex);
                if (m_labeling->GetLabel(id) == id) { // THIS IS WRONG!!!
                    return HIT_EXTREMUM; // reached a critical point
                }
            }
            else {
                Vec3l closest_vertex = m_grid->Inbounds(new_point + 0.5);
                INDEX_TYPE id = m_grid->Index3d(closest_vertex);
                if (m_labeling->GetLabel(id) == id) { // THIS IS WRONG!!!
                    return HIT_EXTREMUM; // reached a critical point
                }

            }
            return NONE;
        }

    };

    class TerminateNearAssigned : public AdvectionChecker {
    protected:
        DenseLabeling<INDEX_TYPE>* m_labeling;
        RegularGrid* m_grid;
    public:
        TerminateNearAssigned(DenseLabeling<INDEX_TYPE>* labeling, RegularGrid*grid) : m_labeling(labeling), m_grid(grid) {}
        virtual ADVECTION_EVENT CheckAndDoStuff(const Vec3d& old_point, const Vec3d& new_point, bool nearboundary) {
            if (!nearboundary) {
                Vec3l closest_vertex = new_point + 0.5;
                INDEX_TYPE id = m_grid->Index3d(closest_vertex);
                if (m_labeling->GetLabel(id) != -1) {
                    return HIT_PREASSIGNED; // reached a critical point
                }
            }
            else {
                Vec3l closest_vertex = m_grid->Inbounds(new_point + 0.5);
                INDEX_TYPE id = m_grid->Index3d(closest_vertex);
                if (m_labeling->GetLabel(id) != -1) {
                    return HIT_PREASSIGNED; // reached a critical point
                }

            }
            return NONE;
        }

    };
    class TerminateNearCertain : public AdvectionChecker {
    protected:
        DenseLabeling<int>* m_labeling;
        RegularGrid* m_grid;
    public:
        TerminateNearCertain(DenseLabeling<int>* labeling, RegularGrid*grid) : m_labeling(labeling), m_grid(grid) {}
        virtual ADVECTION_EVENT CheckAndDoStuff(const Vec3d& old_point, const Vec3d& new_point, bool nearboundary) {
            if (!nearboundary) {
                Vec3l closest_vertex = new_point + 0.5;
                INDEX_TYPE id = m_grid->Index3d(closest_vertex);
                if (m_labeling->GetLabel(id) != -1) {
                    return HIT_PREASSIGNED; // reached a critical point
                }
            }
            else {
                Vec3l closest_vertex = m_grid->Inbounds(new_point + 0.5);
                INDEX_TYPE id = m_grid->Index3d(closest_vertex);
                if (m_labeling->GetLabel(id) != -1) {
                    return HIT_PREASSIGNED; // reached a critical point
                }

            }
            return NONE;
        }

    };

    class TerminateNearPathCompressedRegion : public AdvectionChecker {
    protected:
        DenseLabeling<DestType>* m_labeling;
        RegularGrid* m_grid;
    public:
        TerminateNearPathCompressedRegion(DenseLabeling<DestType>* labeling, RegularGrid*grid) : m_labeling(labeling), m_grid(grid) {}
        virtual ADVECTION_EVENT CheckAndDoStuff(const Vec3d& old_point, const Vec3d& new_point, bool nearboundary) {
            if (!nearboundary) {
                Vec3l closest_vertex = new_point + 0.5;
                INDEX_TYPE id = m_grid->Index3d(closest_vertex);
                if (m_labeling->GetLabel(id) == DestType::ASSIGNED || m_labeling->GetLabel(id) == DestType::CERTAIN_TERMINAL) {
                    return HIT_PREASSIGNED; // reached a path-compressed region
                }
            }
            else {
                Vec3l closest_vertex = m_grid->Inbounds(new_point + 0.5);
                INDEX_TYPE id = m_grid->Index3d(closest_vertex);
                if (m_labeling->GetLabel(id) == DestType::ASSIGNED || m_labeling->GetLabel(id) == DestType::CERTAIN_TERMINAL) {
                    return HIT_PREASSIGNED; // reached a path-compressed region
                }

            }
            return NONE;
        }

    };

    class TerminateNearOriginalCertain : public AdvectionChecker {
    protected:
        DenseLabeling<DestType>* m_labeling;
        RegularGrid* m_grid;
    public:
        TerminateNearOriginalCertain(DenseLabeling<DestType>* labeling, RegularGrid*grid) : m_labeling(labeling), m_grid(grid) {}
        virtual ADVECTION_EVENT CheckAndDoStuff(const Vec3d& old_point, const Vec3d& new_point, bool nearboundary) {
            if (!nearboundary) {
                Vec3l closest_vertex = new_point + 0.5;
                INDEX_TYPE id = m_grid->Index3d(closest_vertex);
                if (m_labeling->GetLabel(id) == DestType::CERTAIN_TERMINAL) {
                    return HIT_PREASSIGNED;
                }
            }
            else {
                Vec3l closest_vertex = m_grid->Inbounds(new_point + 0.5);
                INDEX_TYPE id = m_grid->Index3d(closest_vertex);
                if (m_labeling->GetLabel(id) == DestType::CERTAIN_TERMINAL) {
                    return HIT_PREASSIGNED;
                }

            }
            return NONE;
        }

    };

    // THIS IS NOT FINISHED!!!
    class TerminateNearAssignedAndUpdate : public AdvectionChecker {
    protected:
        DenseLabeling<bool>* m_labeling;
        RegularGrid* m_grid;
    public:
        virtual ADVECTION_EVENT CheckAndDoStuff(const Vec3d& old_point, const Vec3d& new_point, bool nearboundary) {
            if (!nearboundary) {
                Vec3l closest_vertex = new_point + 0.5;
                if (m_labeling->GetLabel(m_grid->Index3d(closest_vertex))) {
                    return HIT_EXTREMUM; // reached a critical point
                }
            }
            return NONE;
        }

    };




    class AdvectionChecker2D {
    public:

        virtual ADVECTION_EVENT CheckAndDoStuff(const Vec2d& old_point, const Vec2d& new_point, bool nearboundary){ return NONE; }
    };

    // no early termination during advection interior to a voxel
    class NoTermination2D : public AdvectionChecker2D {
    public:
        virtual ADVECTION_EVENT CheckAndDoStuff(const Vec2d& old_point, const Vec2d& new_point, bool nearboundary) {
            return NONE;
        }
    };



    class TerminateNearExtrema : public AdvectionChecker {
    protected:
        RegularGrid* m_grid;
    public:
        std::unordered_set<INDEX_TYPE> m_extrema;
        TerminateNearExtrema(const std::unordered_set<INDEX_TYPE>& extrema, RegularGrid* grid) : m_grid(grid), m_extrema(extrema) {}
        TerminateNearExtrema(const std::vector<INDEX_TYPE>& extrema, RegularGrid* grid) : m_grid(grid) {
            for(size_t i = 0; i < extrema.size(); i++)
                m_extrema.insert(extrema[i]);
        }
        void Printstuff() {
            printf("Termination minima are: \n");
            for (auto it = m_extrema.begin(); it != m_extrema.end(); it++) {
                Vec3i co = m_grid->XYZ3d(*it);
                printf("id %d = (%d, %d, %d)\n", *it, co[0], co[1], co[2]);
            }
        }

        INDEX_TYPE WhichExtremum(Vec3d point) {
            Vec3l closest_vertex = m_grid->Inbounds(point + 0.5); // dont know if our point is near boundary!
            //printf("(%f, %f, %f)->(%d, %d, %d)\n", point[0], point[1], point[2], closest_vertex[0], closest_vertex[1], closest_vertex[2]);
            return m_grid->Index3d(closest_vertex);
        }

        //void AddExtremum(INDEX_TYPE id) { m_extrema.insert(id); }
        virtual ADVECTION_EVENT CheckAndDoStuff(const Vec3d& old_point, const Vec3d& new_point, bool nearboundary) {
            if (!nearboundary) {
                Vec3l closest_vertex = new_point + 0.5;
                INDEX_TYPE id = m_grid->Index3d(closest_vertex);
                if (m_extrema.count(id) != 0) {
                    return HIT_PREASSIGNED; // reached a critical point
                }
            }
            else {
                Vec3l closest_vertex = m_grid->Inbounds(new_point + 0.5);
                INDEX_TYPE id = m_grid->Index3d(closest_vertex);
                if (m_extrema.count(id) != 0) {
                    return HIT_PREASSIGNED; // reached a critical point
                }

            }
            return NONE;
        }

    };
}
#endif
