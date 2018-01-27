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

#ifndef NUMERIC_INTEGRATOR_EXPANDING_REGION_STOP_H
#define NUMERIC_INTEGRATOR_EXPANDING_REGION_STOP_H

#include <set>
#include <queue>
#include "basic_types.h"
#include "vectors.h"
#include "labeling.h"
#include "regular_grid.h"
#include "regular_grid_trilinear_function.h"
#include "adaptive_euler_advector.h"
#include "timing.h"
//#include "union_find_labeling.h"

//#define OUTPUTINTERMEDIATE
namespace MSC {


    class IndexCompareLessThan {
    protected:
        RegularGridTrilinearFunction* m_func;
    public:
        IndexCompareLessThan(RegularGridTrilinearFunction* func) : m_func(func) {}
        inline bool Compare(INDEX_TYPE a, INDEX_TYPE b) const {
            return m_func->IsGreater(b, a);
        }
        inline bool CompareF(FLOATTYPE a, FLOATTYPE b) const {
            return a < b;
        }
        bool operator()(INDEX_TYPE a, INDEX_TYPE b) const {
            return m_func->IsGreater(a, b);
        }
    };
    class IndexCompareGreaterThan {
    protected:
        RegularGridTrilinearFunction* m_func;
    public:
        IndexCompareGreaterThan(RegularGridTrilinearFunction* func) : m_func(func) {}
        inline bool Compare(INDEX_TYPE a, INDEX_TYPE b) const {
            return m_func->IsGreater(a, b);
        }
        inline bool CompareF(FLOATTYPE a, FLOATTYPE b) const {
            return a > b;
        }
        bool operator()(INDEX_TYPE a, INDEX_TYPE b) const {
            return m_func->IsGreater(a, b);
        }
    };


    template< class Advector, class Comparer>
    class NumericIntegratorExpandingRegionStop {

    protected:

        //struct bridge {
        //	INDEX_TYPE labela;
        //	INDEX_TYPE labelb;
        //	float value;
        //};
        //

        DenseLabeling<INDEX_TYPE>* m_destinations;
        DenseLabeling<int>* m_certains;
        RegularGrid* m_grid;
        RegularGridTrilinearFunction* m_func;
        int m_num_iterations_left;

        Vec3i m_xyz;
        Vec3b m_periodic;
        FLOATTYPE m_error_threshold;
        FLOATTYPE m_gradient_threshold;

        char* m_fname;

        //struct cmp{
        //	bool operator()(INDEX_TYPE a, INDEX_TYPE b) {
        //		return mCompare->Compare(a, b);
        //	}
        //};

        bool IsExtremeVertexIn6Neighborhood(INDEX_TYPE id) const {
            Vec3l t_neighbors[6];
            Vec3l t_coords = m_grid->XYZ3d(id);
            int t_num_neighbors = m_grid->GatherExistingNeighbors6(t_coords, t_neighbors);

            INDEX_TYPE t_current_lowest = id;
            for (int i = 0; i < t_num_neighbors; i++) {
                INDEX_TYPE t_neighbor_vertex = m_grid->Index3d(t_neighbors[i]);
                if (mCompare->Compare(t_neighbor_vertex, t_current_lowest)) {
                    return false;
                }
            }
            return true;
        }
        void Enqueue_Later_Neighbors(Vec3l xyz, std::priority_queue<INDEX_TYPE, std::vector<INDEX_TYPE>, Comparer > &expansion, std::set<INDEX_TYPE>&seen) {
            INDEX_TYPE tid = m_grid->Index3d(xyz);

            Vec3l neighbors[6];
            int nn = m_grid->GatherExistingNeighbors6(xyz, neighbors);
            for (int i = 0; i < nn; i++) {
                INDEX_TYPE tneg = m_grid->Index3d(neighbors[i]);
                if (m_certains->GetLabel(tneg) == -1 && mCompare->Compare(tid, tneg) && seen.count(tneg) == 0) {
                    seen.insert(tneg);
                    expansion.push(tneg);
                }

            }
        }

        int Inspect_Higher_Certains(INDEX_TYPE tid) {
            INDEX_TYPE tneg;
            int extremal_certain = m_certains->GetLabel(tid);
            bool has_extremal = false;

            Vec3l neighbors[6];
            int nn = m_grid->GatherExistingNeighbors6(m_grid->XYZ3d(tid), neighbors);
            for (int i = 0; i < nn; i++) {
                INDEX_TYPE tneg = m_grid->Index3d(neighbors[i]);
                if (mCompare->Compare(tneg, tid)) {
                    if (m_certains->GetLabel(tneg) < 0) return -1; // if a extremal one is uncertain, we are uncertain
                    if (!has_extremal) {
                        extremal_certain = m_certains->GetLabel(tneg);
                        has_extremal = true;
                    }
                    else {
                        if (extremal_certain != m_certains->GetLabel(tneg)) return -1;
                    }
                }
            }

            if (!has_extremal) {
                printf("ERROR should never get here\n");
                return -1;
            }
            return extremal_certain;

        }

        void Expand_Lower_Neighborhood(INDEX_TYPE startid) {
            Vec3l xyz = m_grid->XYZ3d(startid);
            std::set<INDEX_TYPE> seen;

            INDEX_TYPE tid = startid;
            // the natural ordering using the < operator on pairs will give us the highest
            // element first, simulating region growing from high to low
            std::priority_queue<INDEX_TYPE, std::vector<INDEX_TYPE>, Comparer > growing_front(*mCompare);
            seen.insert(startid);
            Enqueue_Later_Neighbors(xyz, growing_front, seen);

            while (!growing_front.empty()) {

                INDEX_TYPE currid = growing_front.top();
                growing_front.pop();

                int cellvale = Inspect_Higher_Certains(currid);
                // find extremals

                // cellvalue >=0 indicates that there is certainty here, so lets expand
                if (cellvale >= 0) {
                    m_certains->SetLabel(currid, cellvale);
                    m_destinations->SetLabel(currid, startid);
                    Enqueue_Later_Neighbors(m_grid->XYZ3d(currid), growing_front, seen);
                }

            }
        }


        Comparer* mCompare;
        std::unordered_set<INDEX_TYPE> m_extrema;
    public:

        bool mDoHackCutoff;

        NumericIntegratorExpandingRegionStop(RegularGridTrilinearFunction* func, RegularGrid* grid,  FLOATTYPE error_threshold, FLOATTYPE gradient_threshold, int interation_limit) :
            m_num_iterations_left(interation_limit), m_xyz(func->GetGrid()->XYZ()), m_periodic(func->GetGrid()->Periodic()), m_func(func), m_grid(grid),
            m_gradient_threshold(gradient_threshold), m_error_threshold(error_threshold) {
            mCompare = new Comparer(func);
        }

        ~NumericIntegratorExpandingRegionStop() {
            delete m_destinations;
            delete m_certains;
            delete mCompare;
        }


        DenseLabeling<INDEX_TYPE>* GetOutputLabels() { return m_destinations; }
        RegularGrid* GetGrid() { return m_grid; }
        RegularGridTrilinearFunction* GetFunction() { return m_func; }
        const std::unordered_set<INDEX_TYPE>& GetExtrema() const { return m_extrema; }
        void BeginIntegration() {
            printf(" -- Performing numeric integration for volume assignment...\n");

            //m_func->ComputeGradFromImage(m_rkindex);

            m_destinations = new DenseLabeling<INDEX_TYPE>(m_grid->NumElements());
            m_certains = new DenseLabeling<int>(m_grid->NumElements());
            const INDEX_TYPE t_num_vertices = m_grid->NumElements();
            AdvectionChecker* inside_voxel_critical_advection_checker = new TerminateNearCertain(m_certains, m_grid);
            AdvectionChecker* no_check = new NoTermination();//AdvectionChecker* inside_voxel_advection_checker = new TerminateNearAssigned(m_destinations, m_grid);

            printf("   -- finding extrema to terminate integral lines...");
            fflush(stdout);

            std::vector<INDEX_TYPE> extrema;

            // set all potential extrema, so we terminate near them
#pragma omp parallel for
            for (INDEX_TYPE i = 0; i < t_num_vertices; i++) {
                m_certains->SetLabel(i, -1);
                if (IsExtremeVertexIn6Neighborhood(i)) {
                    m_destinations->SetLabel(i, i);
#pragma omp critical
                    {
                        extrema.push_back(i);
                        m_extrema.insert(i);
                    }
                }
                else {
                    m_destinations->SetLabel(i, -1);
                }

            }
            printf(" done! found %d extrema!\n", m_extrema.size());
            printf("   -- expanding extrema certain regions...");
            fflush(stdout);

            int num_extrema = extrema.size();
#pragma omp parallel shared(extrema)
            {
#pragma omp for schedule(dynamic)  nowait
                for (int m = 0; m < num_extrema; m++) {
                    INDEX_TYPE maximum = extrema[m];
                    m_certains->SetLabel(maximum, m);
                    Expand_Lower_Neighborhood(maximum);
                    m_func->SetGradExplicit(maximum, Vec3d(0, 0, 0));
                }
            }
#ifdef OUTPUTINTERMEDIATE
            m_certains->OutputToFile("certain_expansion.raw");
            TopologicalRegularGridRestricted* ttgrid = new TopologicalRegularGridRestricted(m_grid);
            VertexLabelingToBoundaryLabeling<int>* tedge = new VertexLabelingToBoundaryLabeling<int>(m_certains, ttgrid);
            tedge->ComputeBoundary();
            tedge->OutputEdgesToFile("certain_edges.txt");
            FILE* fout = fopen("Linesout.txt", "w");
#endif
            printf(" expansion done!\n", extrema.size());
            printf("   -- doing numerical integration first pass with path compression...");
            fflush(stdout);

            //m_destinations->OutputToIntFile("dest.raw");
            int t1, t2; t1 = t2 = 0;
#pragma omp parallel
    {
                Advector t_advector(m_grid, m_func, m_gradient_threshold, m_error_threshold, inside_voxel_critical_advection_checker);
            std::vector<INDEX_TYPE> t_path;
            t_path.reserve(100);
#pragma omp for schedule(guided)  nowait
            for (INDEX_TYPE i = 0; i < t_num_vertices; i++) {
                // early skip if this is already a maximum
                if (m_certains->GetLabel(i) != -1 || m_destinations->GetLabel(i) != -1) {
                    continue;
                }


                t_path.clear();
                t_path.push_back(i);
                Vec3l t_coords = m_grid->XYZ3d(i); // get the coordinates of the poitn
                //if (t_coords[0] == 0 && t_coords[1] == 0) printf("doing %d\n", t_coords[2]);
                Vec3d t_current_point = t_coords;
                int t_num_iterations_left = m_num_iterations_left;

                bool t_continue = true;
#ifdef OUTPUTINTERMEDIATE
                std::vector<Vec3d> line_soup;
                line_soup.push_back(t_current_point);
#endif
                while (t_continue) {
                    Vec3d t_next_point;
                    ADVECTION_EVENT t_return_code;
                    if (m_grid->DistToBoundary(t_coords) <= 1) {
                        t_return_code = t_advector.AdvectThroughVoxelNearBoundary(t_current_point, t_num_iterations_left);
                        t_coords = m_grid->Inbounds(t_current_point + 0.5); // get nearest integer voxel
                        t1++;
                    }
                    else {
                        t_return_code = t_advector.AdvectThroughVoxelNoCheck(t_current_point, t_num_iterations_left);
                        t_coords = (t_current_point + 0.5);
                        t2++;
                    }
                    INDEX_TYPE t_next_id = m_grid->Index3d(t_coords);
                    t_path.push_back(t_next_id);
#ifdef OUTPUTINTERMEDIATE
                    line_soup.push_back(t_current_point);
#endif					// if we terminated or hit a critical point, then we are done
                    if (t_return_code == ADVECTION_EVENT::LOW_GRADIENT ||
                        t_return_code == ADVECTION_EVENT::HIT_EXTREMUM ||
                        t_return_code == ADVECTION_EVENT::HIT_PREASSIGNED ||
                        t_return_code == ADVECTION_EVENT::OVER_MAX_ITERATIONS) {
                        INDEX_TYPE t_dest_label = m_destinations->GetLabel(t_next_id);
                        int t_certain_label = m_certains->GetLabel(t_next_id);
//#pragma omp critical
                        {
                            for (int j = 0; j < t_path.size(); j++) {
                                m_destinations->SetLabel(t_path[j], t_dest_label);
                                //m_certains->SetLabel(t_path[j], t_certain_label);
                            }
#ifdef OUTPUTINTERMEDIATE
#pragma omp critical
                            {
                                int tn = omp_get_thread_num();
                                fprintf(fout, "%d %d %d %d\n", i, tn, line_soup.size(), t_dest_label);
                                for (int j = 0; j < line_soup.size(); j++) {
                                    fprintf(fout, "%f %f %f\n", line_soup[j][0], line_soup[j][1], line_soup[j][2]);
                                }
                            }
#endif
                        }
                        t_continue = false;
                    }
                }
            }
#pragma omp single
            {
#ifdef OUTPUTINTERMEDIATE
                fclose(fout);
                m_destinations->OutputToIntFile("first_integration.raw");
#endif
                printf(" done!\n");
                printf("   -- checking unambiguous voxels...");
                fflush(stdout);
            }
#pragma omp for schedule(static)
            for (INDEX_TYPE i = 0; i < t_num_vertices; i++) {
                // early skip if this is already a maximum
                if (m_certains->GetLabel(i) != -1) {
                    continue;
                }

                Vec3l t_coords = m_grid->XYZ3d(i); // get the coordinates of the poitn
                //if (t_coords[0] == 0 && t_coords[1] == 0) printf("doing %d\n", t_coords[2]);
                Vec3l negs[6];

                int nn = m_grid->GatherExistingNeighbors6(t_coords, negs);

                bool needsupdate = false;
                INDEX_TYPE t_dest_label = m_destinations->GetLabel(i);
                for (int j = 0; j < nn; j++)  {
                    if (m_destinations->GetLabel(m_grid->Index3d(negs[j])) != t_dest_label) {
                        needsupdate = true;
                        break;
                    }
                }

                if (needsupdate) continue;

                m_certains->SetLabel(i, 1);
            }
#pragma omp single
            {
                printf(" done!\n");
                printf("   -- re-integrating boundaries that may have ended too soon...");
                fflush(stdout);
            }
#pragma omp for schedule(guided)  nowait
                for (INDEX_TYPE i = 0; i < t_num_vertices; i++) {
                    // early skip if this is already a maximum
                    if (m_certains->GetLabel(i) != -1) {
                        continue;
                    }

                    Vec3l t_coords = m_grid->XYZ3d(i); // get the coordinates of the poitn
                    //if (t_coords[0] == 0 && t_coords[1] == 0) printf("doing %d\n", t_coords[2]);



                    Vec3d t_current_point = t_coords;
                    int t_num_iterations_left = 500;
                    bool t_continue = true;

                    while (t_continue) {
                        Vec3d t_next_point;
                        ADVECTION_EVENT t_return_code;
                        if (m_grid->DistToBoundary(t_coords) <= 1) {
                            t_return_code = t_advector.AdvectThroughVoxelNearBoundary(t_current_point, t_num_iterations_left);
                            t_coords = m_grid->Inbounds(t_current_point + 0.5); // get nearest integer voxel
                            t1++;
                        }
                        else {
                            t_return_code = t_advector.AdvectThroughVoxelNoCheck(t_current_point, t_num_iterations_left);
                            t_coords = (t_current_point + 0.5);
                            t2++;
                        }
                        INDEX_TYPE t_next_id = m_grid->Index3d(t_coords);
                        t_path.push_back(t_next_id);
                        // if we terminated or hit a critical point, then we are done
                        if (t_return_code == ADVECTION_EVENT::LOW_GRADIENT ||
                            t_return_code == ADVECTION_EVENT::HIT_EXTREMUM ||
                            t_return_code == ADVECTION_EVENT::HIT_PREASSIGNED ||
                            t_return_code == ADVECTION_EVENT::OVER_MAX_ITERATIONS) {
                            m_destinations->SetLabel(i, m_destinations->GetLabel(t_next_id));
                            t_continue = false;
                        }
                    }
                }

                //if (m_grid->DistToBoundary()) {}
                //if (IsLowestVertexIn6Neighborhood(i)) {
                //	m_destinations->SetLabel(i, i);
                //}
                //else {
                //	m_destinations->SetLabel(i, -1);
                //}



            //printf("%d boundary, %d noboundary\n", t1, t2);
            }
#ifdef OUTPUTINTERMEDIATE
            m_destinations->OutputToIntFile("reintegration.raw");
#endif
            printf(" done!\n");
            printf(" -- Done numeric integration of 3d volume!\n");

    }

    void BeginIntegration(ThreadedTimer* timer) {
        enum TimingType {
            I_FINDMINS,
            I_EXPAND_CERTAINS,
            I_INTEGRATION,
            I_REINTEGRATION,
            I_DONE
        };
        printf(" -- Performing numeric integration for volume assignment...\n");

        //m_func->ComputeGradFromImage(m_rkindex);

        m_destinations = new DenseLabeling<INDEX_TYPE>(m_grid->NumElements());
        m_certains = new DenseLabeling<int>(m_grid->NumElements());
        const INDEX_TYPE t_num_vertices = m_grid->NumElements();
        AdvectionChecker* inside_voxel_critical_advection_checker = new TerminateNearCertain(m_certains, m_grid);
        AdvectionChecker* no_check = new NoTermination();//AdvectionChecker* inside_voxel_advection_checker = new TerminateNearAssigned(m_destinations, m_grid);

        printf("   -- finding extrema to terminate integral lines...");
        fflush(stdout);

        std::vector<INDEX_TYPE> extrema;
#pragma omp parallel
        {
            int thread_num = omp_get_thread_num();
            timer->RecordThreadActivity(thread_num, I_FINDMINS);

        }
        // set all potential extrema, so we terminate near them
#pragma omp parallel for
        for (INDEX_TYPE i = 0; i < t_num_vertices; i++) {
            m_certains->SetLabel(i, -1);
            if (IsExtremeVertexIn6Neighborhood(i)) {
                m_destinations->SetLabel(i, i);
#pragma omp critical
                {
                    extrema.push_back(i);
                    m_extrema.insert(i);
                }
            }
            else {
                m_destinations->SetLabel(i, -1);
            }

        }
        printf(" found %d extrema!\n", extrema.size());
        printf("   -- expanding extrema certain regions...");
        fflush(stdout);
#pragma omp parallel
        {
            int thread_num = omp_get_thread_num();
            timer->RecordThreadActivity(thread_num, I_EXPAND_CERTAINS);

        }
        int num_extrema = extrema.size();
#pragma omp parallel shared(extrema)
        {
#pragma omp for schedule(dynamic)  nowait
            for (int m = 0; m < num_extrema; m++) {
                INDEX_TYPE maximum = extrema[m];
                m_certains->SetLabel(maximum, m);
                Expand_Lower_Neighborhood(maximum);
                m_func->SetGradExplicit(maximum, Vec3d(0, 0, 0));
            }
        }

        printf(" expansion done!\n", extrema.size());
        printf("   -- doing numerical integration first pass with path compression...");
        fflush(stdout);
#pragma omp parallel
        {
            int thread_num = omp_get_thread_num();
            timer->RecordThreadActivity(thread_num, I_INTEGRATION);

        }
        //m_destinations->OutputToIntFile("dest.raw");
        int t1, t2; t1 = t2 = 0;
#pragma omp parallel
        {
            Advector t_advector(m_grid, m_func, m_gradient_threshold, m_error_threshold, inside_voxel_critical_advection_checker);
            std::vector<INDEX_TYPE> t_path;
            t_path.reserve(100);
#pragma omp for schedule(guided)  nowait
            for (INDEX_TYPE i = 0; i < t_num_vertices; i++) {
                // early skip if this is already a maximum
                if (m_certains->GetLabel(i) != -1 || m_destinations->GetLabel(i) != -1) {
                    continue;
                }

                t_path.clear();
                t_path.push_back(i);
                Vec3l t_coords = m_grid->XYZ3d(i); // get the coordinates of the poitn
                //if (t_coords[0] == 0 && t_coords[1] == 0) printf("doing %d\n", t_coords[2]);
                Vec3d t_current_point = t_coords;
                int t_num_iterations_left = m_num_iterations_left;

                bool t_continue = true;

                while (t_continue) {
                    Vec3d t_next_point;
                    ADVECTION_EVENT t_return_code;
                    if (m_grid->DistToBoundary(t_coords) <= 1) {
                        t_return_code = t_advector.AdvectThroughVoxelNearBoundary(t_current_point, t_num_iterations_left);
                        t_coords = m_grid->Inbounds(t_current_point + 0.5); // get nearest integer voxel
                        t1++;
                    }
                    else {
                        t_return_code = t_advector.AdvectThroughVoxelNoCheck(t_current_point, t_num_iterations_left);
                        t_coords = (t_current_point + 0.5);
                        t2++;
                    }
                    INDEX_TYPE t_next_id = m_grid->Index3d(t_coords);
                    t_path.push_back(t_next_id);
                    // if we terminated or hit a critical point, then we are done
                    if (t_return_code == ADVECTION_EVENT::LOW_GRADIENT ||
                        t_return_code == ADVECTION_EVENT::HIT_EXTREMUM ||
                        t_return_code == ADVECTION_EVENT::HIT_PREASSIGNED ||
                        t_return_code == ADVECTION_EVENT::OVER_MAX_ITERATIONS) {
                        INDEX_TYPE t_dest_label = m_destinations->GetLabel(t_next_id);
                        int t_certain_label = m_certains->GetLabel(t_next_id);
                        //#pragma omp critical
                        {
                            for (int j = 0; j < t_path.size(); j++) {
                                m_destinations->SetLabel(t_path[j], t_dest_label);
                                //m_certains->SetLabel(t_path[j], t_certain_label);
                            }
                        }
                        t_continue = false;
                    }
                }
            }
#pragma omp single
            {
                printf(" done!\n");
                printf("   -- checking unambiguous voxels...");
                fflush(stdout);
            }

            int thread_num = omp_get_thread_num();
            timer->RecordThreadActivity(thread_num, I_REINTEGRATION);


#pragma omp for schedule(static)
            for (INDEX_TYPE i = 0; i < t_num_vertices; i++) {
                // early skip if this is already a maximum
                if (m_certains->GetLabel(i) != -1) {
                    continue;
                }

                Vec3l t_coords = m_grid->XYZ3d(i); // get the coordinates of the poitn
                //if (t_coords[0] == 0 && t_coords[1] == 0) printf("doing %d\n", t_coords[2]);
                Vec3l negs[6];

                int nn = m_grid->GatherExistingNeighbors6(t_coords, negs);

                bool needsupdate = false;
                INDEX_TYPE t_dest_label = m_destinations->GetLabel(i);
                for (int j = 0; j < nn; j++)  {
                    if (m_destinations->GetLabel(m_grid->Index3d(negs[j])) != t_dest_label) {
                        needsupdate = true;
                        break;
                    }
                }

                if (needsupdate) continue;

                m_certains->SetLabel(i, 1);
            }
#pragma omp single
            {
                printf(" done!\n");
                printf("   -- re-integrating boundaries that may have ended too soon...");
                fflush(stdout);
            }
#pragma omp for schedule(guided)  nowait
            for (INDEX_TYPE i = 0; i < t_num_vertices; i++) {
                // early skip if this is already a maximum
                if (m_certains->GetLabel(i) != -1) {
                    continue;
                }

                Vec3l t_coords = m_grid->XYZ3d(i); // get the coordinates of the poitn
                //if (t_coords[0] == 0 && t_coords[1] == 0) printf("doing %d\n", t_coords[2]);



                Vec3d t_current_point = t_coords;
                int t_num_iterations_left = 500;
                bool t_continue = true;

                while (t_continue) {
                    Vec3d t_next_point;
                    ADVECTION_EVENT t_return_code;
                    if (m_grid->DistToBoundary(t_coords) <= 1) {
                        t_return_code = t_advector.AdvectThroughVoxelNearBoundary(t_current_point, t_num_iterations_left);
                        t_coords = m_grid->Inbounds(t_current_point + 0.5); // get nearest integer voxel
                        t1++;
                    }
                    else {
                        t_return_code = t_advector.AdvectThroughVoxelNoCheck(t_current_point, t_num_iterations_left);
                        t_coords = (t_current_point + 0.5);
                        t2++;
                    }
                    INDEX_TYPE t_next_id = m_grid->Index3d(t_coords);
                    t_path.push_back(t_next_id);
                    // if we terminated or hit a critical point, then we are done
                    if (t_return_code == ADVECTION_EVENT::LOW_GRADIENT ||
                        t_return_code == ADVECTION_EVENT::HIT_EXTREMUM ||
                        t_return_code == ADVECTION_EVENT::HIT_PREASSIGNED ||
                        t_return_code == ADVECTION_EVENT::OVER_MAX_ITERATIONS) {
                        m_destinations->SetLabel(i, m_destinations->GetLabel(t_next_id));
                        t_continue = false;
                    }
                }
            }

            //if (m_grid->DistToBoundary()) {}
            //if (IsLowestVertexIn6Neighborhood(i)) {
            //	m_destinations->SetLabel(i, i);
            //}
            //else {
            //	m_destinations->SetLabel(i, -1);
            //}



            //printf("%d boundary, %d noboundary\n", t1, t2);
        }
        printf(" done!\n");
        printf(" -- Done numeric integration of 3d volume!\n");
#pragma omp parallel
        {
            int thread_num = omp_get_thread_num();
            timer->RecordThreadActivity(thread_num, I_DONE);

        }

    }

    };




}
#endif
