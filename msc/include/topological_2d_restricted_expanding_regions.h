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

#ifndef TOPOLOGICAL_2D_RESTRICTED_EXPANDING_REGIONS_H
#define TOPOLOGICAL_2D_RESTRICTED_EXPANDING_REGIONS_H

//#include <algorithm>
#include <stdio.h>
#include <set>
#include <map>
#include <queue>
#include <stack>
#include "basic_types.h"
#include "vectors.h"
#include "labeling.h"
#include "regular_grid.h"
#include "regular_grid_trilinear_function.h"
#include "adaptive_euler_advector.h"
#include "topological_regular_grid.h"
#include "array_index_partition.h"
#include "topological_explicit_mesh_function.h"
//#include "topological_explicit_mesh_function.h"

namespace MSC {



    // some naming conventions:
    //
    // the word *primal* refers to anything to do with the "original" mesh, i.e. the mesh where
    // image values were provided on vertices. there are primal_vertices, primal_edges, primal_quads,
    // and primal_voxels
    //
    // the word *dual* refers to anything to do with the dual of the mesh. i.e. a dual_vertex corresponds
    // to a primal_voxel, etc. The correspondence is based on index in the topological mesh. a cell and its
    // dual will also share the same centroid, which is the same as the coordiante location given by the
    // topological mesh (except the topological mesh multiplies by 2 in each direction)
    //
    // any index that refers to the topological grid has *topo* id in front of it, and any that refers to a vertex
    // in the space of the function has *grid* in front of it. e.g. topo_id = 10 is the 5th vertex, where grid_id=10
    // is the 10th vertex
    //
    // i call an object that is not excluded by the restriction labels as "permissable"
    class Topological2DRestrictedExpandingRegions {

    public:

        // will need some other kind of tuple and not dense?
        DenseLabeling<char>* m_restriction_labels; // 1 if we can, 0 if we can;t
        TopologicalExplicitDenseMeshFunction* m_func;
        TopologicalRegularGrid* m_topo_grid;

        // consistent edge comparisons
        inline bool IsGreater(INDEX_TYPE a, INDEX_TYPE b) const {
            return m_func->lessThan(b, a);
        }


        bool IsHighestPrimalEdgeInUNASSIGNEDRestrictedNeighborhood(INDEX_TYPE primal_edge_topo_id) const {

            // return false if the edge is not part of our permissable set
            if (m_restriction_labels->GetLabel(primal_edge_topo_id) != 1 || TEST_PRIOR_STUFF->GetLabel(primal_edge_topo_id) != 0) {
                return false;
            }

            // to get surrounding edges, look at surrounding quads, and all their edges
            TopologicalRegularGrid::CofacetsIterator t_adjacent_primal_quad_iterator(m_topo_grid);
            TopologicalRegularGrid::FacetsIterator t_neighboring_primal_edge_iterator(m_topo_grid);
            for (t_adjacent_primal_quad_iterator.begin(primal_edge_topo_id);
                t_adjacent_primal_quad_iterator.valid();
                t_adjacent_primal_quad_iterator.advance()){
                INDEX_TYPE t_adjacent_primal_quad_topo_id = t_adjacent_primal_quad_iterator.value();
                for (t_neighboring_primal_edge_iterator.begin(t_adjacent_primal_quad_topo_id);
                    t_neighboring_primal_edge_iterator.valid();
                    t_neighboring_primal_edge_iterator.advance()){
                    INDEX_TYPE t_neighboring_primal_edge_topo_id = t_neighboring_primal_edge_iterator.value();
                    // ignore an edge if it is not permissable

                    if (m_restriction_labels->GetLabel(t_neighboring_primal_edge_topo_id) != 1 || TEST_PRIOR_STUFF->GetLabel(t_neighboring_primal_edge_topo_id) != 0) continue;

                    // else compare the values;
                    if (IsGreater(t_neighboring_primal_edge_topo_id, primal_edge_topo_id)) return false;
                }
            }
            return true;
        }


        bool IsHighestPrimalEdgeInRestrictedNeighborhood(INDEX_TYPE primal_edge_topo_id) const {

            // return false if the edge is not part of our permissable set
            if (m_restriction_labels->GetLabel(primal_edge_topo_id) != 1) {
                return false;
            }

            // to get surrounding edges, look at surrounding quads, and all their edges
            TopologicalRegularGrid::CofacetsIterator t_adjacent_primal_quad_iterator(m_topo_grid);
            TopologicalRegularGrid::FacetsIterator t_neighboring_primal_edge_iterator(m_topo_grid);
            for (t_adjacent_primal_quad_iterator.begin(primal_edge_topo_id);
                t_adjacent_primal_quad_iterator.valid();
                t_adjacent_primal_quad_iterator.advance()){
                INDEX_TYPE t_adjacent_primal_quad_topo_id = t_adjacent_primal_quad_iterator.value();
                for (t_neighboring_primal_edge_iterator.begin(t_adjacent_primal_quad_topo_id);
                    t_neighboring_primal_edge_iterator.valid();
                    t_neighboring_primal_edge_iterator.advance()){
                    INDEX_TYPE t_neighboring_primal_edge_topo_id = t_neighboring_primal_edge_iterator.value();
                    // ignore an edge if it is not permissable

                    if (m_restriction_labels->GetLabel(t_neighboring_primal_edge_topo_id) != 1) continue;

                    // else compare the values;
                    if (IsGreater(t_neighboring_primal_edge_topo_id, primal_edge_topo_id)) return false;
                }
            }
            return true;
        }

        void Enqueue_Lower_Neighbors(INDEX_TYPE primal_edge_topo_id,
            FLOATTYPE t_original_edge_value,
            std::priority_queue<std::pair<FLOATTYPE, INDEX_TYPE> > &expansion,
            std::set<INDEX_TYPE>&narrow_band) const {
            // return if the edge is not part of our permissable set
            if (m_restriction_labels->GetLabel(primal_edge_topo_id) != 1) {
                return;
            }
            // to get surrounding edges, look at surrounding quads, and all their edges
            TopologicalRegularGrid::CofacetsIterator t_adjacent_primal_quad_iterator(m_topo_grid);
            TopologicalRegularGrid::FacetsIterator t_neighboring_primal_edge_iterator(m_topo_grid);
            for (t_adjacent_primal_quad_iterator.begin(primal_edge_topo_id);
                t_adjacent_primal_quad_iterator.valid();
                t_adjacent_primal_quad_iterator.advance()){
                INDEX_TYPE t_adjacent_primal_quad_topo_id = t_adjacent_primal_quad_iterator.value();
                for (t_neighboring_primal_edge_iterator.begin(t_adjacent_primal_quad_topo_id);
                    t_neighboring_primal_edge_iterator.valid();
                    t_neighboring_primal_edge_iterator.advance()){
                    INDEX_TYPE t_neighboring_primal_edge_topo_id = t_neighboring_primal_edge_iterator.value();

                    // make sure it is permissable and not already in our traversal

                    if (m_restriction_labels->GetLabel(t_neighboring_primal_edge_topo_id) == 1 &&
                        TEST_PRIOR_STUFF->GetLabel(t_neighboring_primal_edge_topo_id) == 0 && // inserted --
                        narrow_band.count(t_neighboring_primal_edge_topo_id) == 0) {
                        // i think we can ignore the value test? - no we can't don't want to "spill over"
                        // into another basin - only insert strictly lower stuff!
                        FLOATTYPE t_neg_value = m_func->cellValue(t_neighboring_primal_edge_topo_id);
                        if (IsGreater(primal_edge_topo_id, t_neighboring_primal_edge_topo_id)) {
                            narrow_band.insert(t_neighboring_primal_edge_topo_id);
                            expansion.push(std::pair<FLOATTYPE, INDEX_TYPE>(t_neg_value, t_neighboring_primal_edge_topo_id));
                        }
                    }

                }
            }
        }

        // returns true if all the higher neighbors are "certain"
        bool Inspect_Higher_Certains(INDEX_TYPE primal_edge_topo_id, FLOATTYPE original_edge_value, std::set<INDEX_TYPE>& processed_band) const {
            INDEX_TYPE tneg;
            bool has_higher = false;

            TopologicalRegularGrid::CofacetsIterator t_adjacent_primal_quad_iterator(m_topo_grid);
            TopologicalRegularGrid::FacetsIterator t_neighboring_primal_edge_iterator(m_topo_grid);
            for (t_adjacent_primal_quad_iterator.begin(primal_edge_topo_id);
                t_adjacent_primal_quad_iterator.valid();
                t_adjacent_primal_quad_iterator.advance()){
                INDEX_TYPE t_adjacent_primal_quad_topo_id = t_adjacent_primal_quad_iterator.value();
                for (t_neighboring_primal_edge_iterator.begin(t_adjacent_primal_quad_topo_id);
                    t_neighboring_primal_edge_iterator.valid();
                    t_neighboring_primal_edge_iterator.advance()){
                    INDEX_TYPE t_neighboring_primal_edge_topo_id = t_neighboring_primal_edge_iterator.value();

                    // skip any non-permissable edges
                    if (m_restriction_labels->GetLabel(t_neighboring_primal_edge_topo_id) != 1
                        || TEST_PRIOR_STUFF->GetLabel(t_neighboring_primal_edge_topo_id) != 0) // inserted --
                        continue;

                    // if our neighbor is highter, check it is certain
                    if (IsGreater(t_neighboring_primal_edge_topo_id, primal_edge_topo_id)) {
                        if (processed_band.count(t_neighboring_primal_edge_topo_id) == 0) {
                            // then its not certain so we can't be certain
                            return false;
                        }
                    }
                }
            }
            return true;
        }


        //keep only growing front?

        struct RegionState {
            INDEX_TYPE region_id;
            INDEX_TYPE bottleneck_id;
            std::priority_queue<std::pair<FLOATTYPE, INDEX_TYPE> > unprocessed_band;
            std::set<INDEX_TYPE> narrow_band;
            std::set<INDEX_TYPE> processed_band;
        };

        std::map<INDEX_TYPE, RegionState*> m_bottlenecks_map;
        std::queue<std::pair<RegionState*, RegionState*> > m_ready_to_merge_list;




        // warning destroys a region and updates the other, returns updated one
        RegionState* MergeRegions(RegionState* a, RegionState* b) const {
            RegionState* aa;
            RegionState* bb;
            if (a->narrow_band.size() > b->narrow_band.size()) {
                aa = a;
                bb = b;
            }
            else  {
                aa = b;
                bb = a;
            }

            aa->processed_band.insert(bb->processed_band.begin(), bb->processed_band.end());
            aa->narrow_band.insert(bb->narrow_band.begin(), bb->narrow_band.end());
            while (!bb->unprocessed_band.empty()) {
                aa->unprocessed_band.push(bb->unprocessed_band.top());
                bb->unprocessed_band.pop();
            }

            aa->region_id = (aa->region_id > bb->region_id ? aa->region_id + 1 : bb->region_id + 1);
            //printf("regionid = %d\n", a->region_id);
            return aa;
        }


        RegionState* Expand_MergeRegions(RegionState* region)  {

            RegionState* local_region = region;

            INDEX_TYPE primal_edge_topo_id = region->unprocessed_band.top().second;
            // the natural ordering using the < operator on pairs will give us the highest
            // element first, simulating region growing from high to low

            std::priority_queue<std::pair<FLOATTYPE, INDEX_TYPE> > growing_front(region->unprocessed_band);
            region->unprocessed_band = std::priority_queue<std::pair<FLOATTYPE, INDEX_TYPE> >(); // clear the old delayed front
            //std::priority_queue<std::pair<FLOATTYPE, INDEX_TYPE> >  new_delayed_front (growing_front);

            while (!growing_front.empty()) {

                INDEX_TYPE currid = growing_front.top().second;
                FLOATTYPE currid_value = growing_front.top().first;
                growing_front.pop();

                // cellvalue >=0 indicates that there is certainty here, so lets expand
                if (Inspect_Higher_Certains(currid, currid_value, local_region->processed_band)) {
                    local_region->processed_band.insert(currid);
                    TEST_OUTPUT->SetLabel(currid, omp_get_thread_num() + 1);
                    m_thread_work_count[omp_get_thread_num()]++;
                    Enqueue_Lower_Neighbors(currid, currid_value, growing_front, local_region->narrow_band);
                }
                else {
                    region->unprocessed_band.push(std::pair<FLOATTYPE, INDEX_TYPE>(currid_value, currid));
                }

            }
            if (!local_region->unprocessed_band.empty()) {
                local_region->bottleneck_id = region->unprocessed_band.top().second;
            }
            //// REGION FOR DEBUG ONLY
            //#pragma omp critical
            //			{
            //				if (m_regions.count(primal_edge_topo_id) != 0)
            //					printf("something went wrong\n");
            //				m_regions[primal_edge_topo_id] = local_region;
            //
            //				if (!local_region->unprocessed_band.empty()) {
            //					INDEX_TYPE tmp = local_region->unprocessed_band.top().second;
            //					if (bottleneckcount.count(tmp) == 0) {
            //						bottleneckcount[tmp] = 1;
            //						bottleneckregions[tmp].push_back(primal_edge_topo_id);
            //					}
            //					else {
            //						bottleneckcount[tmp]++;
            //						bottleneckregions[tmp].push_back(primal_edge_topo_id);
            //					}
            //				}
            //			}
            //// END REGION FOR DEBUG
            return local_region;
        }


        RegionState* Expand_Lower_Neighborhood2(INDEX_TYPE primal_edge_topo_id)  {

            RegionState* local_region = new RegionState();

            local_region->processed_band.clear();
            INDEX_TYPE t_id = primal_edge_topo_id;
            // the natural ordering using the < operator on pairs will give us the highest
            // element first, simulating region growing from high to low

            local_region->narrow_band.insert(primal_edge_topo_id);
            local_region->processed_band.insert(primal_edge_topo_id);

            TEST_OUTPUT->SetLabel(primal_edge_topo_id, omp_get_thread_num() + 1);
            m_thread_work_count[omp_get_thread_num()]++;

            FLOATTYPE original_edge_value = m_func->cellValue(primal_edge_topo_id);
            std::priority_queue<std::pair<FLOATTYPE, INDEX_TYPE> > growing_front;
            Enqueue_Lower_Neighbors(primal_edge_topo_id, original_edge_value, growing_front, local_region->narrow_band);

            while (!growing_front.empty()) {

                INDEX_TYPE currid = growing_front.top().second;
                FLOATTYPE currid_value = growing_front.top().first;
                growing_front.pop();

                // cellvalue >=0 indicates that there is certainty here, so lets expand
                if (Inspect_Higher_Certains(currid, currid_value, local_region->processed_band)) {
                    local_region->processed_band.insert(currid);
                    TEST_OUTPUT->SetLabel(currid, omp_get_thread_num() + 1);
                    m_thread_work_count[omp_get_thread_num()]++;

                    Enqueue_Lower_Neighbors(currid, currid_value, growing_front, local_region->narrow_band);
                }
                else {
                    local_region->unprocessed_band.push(std::pair<FLOATTYPE, INDEX_TYPE>(currid_value, currid));
                }

            }
            if (!local_region->unprocessed_band.empty()) {
                local_region->bottleneck_id = local_region->unprocessed_band.top().second;
            }
            //// THIS WHOLE SECTION IS JUST FOR DEBUGGING
            //#pragma omp critical
            //			{
            //				if (m_regions.count(primal_edge_topo_id) != 0)
            //					printf("something went wrong\n");
            //				m_regions[primal_edge_topo_id] = local_region;
            //
            //				if (!local_region->unprocessed_band.empty()) {
            //					INDEX_TYPE tmp = local_region->unprocessed_band.top().second;
            //					if (bottleneckcount.count(tmp) == 0) {
            //						bottleneckcount[tmp] = 1;
            //						bottleneckregions[tmp].push_back(primal_edge_topo_id);
            //					}
            //					else {
            //						bottleneckcount[tmp]++;
            //						bottleneckregions[tmp].push_back(primal_edge_topo_id);
            //					}
            //				}
            //			}
            local_region->region_id = 1;
            //// END DEBUGGING SECTION

            return local_region;
        }

    public:


        Topological2DRestrictedExpandingRegions(
            TopologicalExplicitDenseMeshFunction* func,
            TopologicalRegularGrid* topo_grid,
            DenseLabeling<char>* restriction_labels) :
            m_func(func), m_restriction_labels(restriction_labels),
            m_topo_grid(topo_grid){
        }

        //DenseLabeling<INDEX_TYPE>* GetOutputLabels() { return m_destinations; }

        TopologicalExplicitDenseMeshFunction* GetFunction() { return m_func; }


        DenseLabeling<int>* TEST_OUTPUT; // debug
        DenseLabeling<char>* TEST_PRIOR_STUFF; // debug
        std::vector<INDEX_TYPE> delayed_cells; // debug
        std::vector<INDEX_TYPE> critical_primal_edges;// debug

        enum WorkType { MERGE_AND_EXPAND, EXPAND_ONLY, NONE };
        struct WorkRequest {
            WorkType worktype;
            INDEX_TYPE start_cell;
            RegionState* region1;
            RegionState* region2;
        };

        std::stack<WorkRequest> m_work_to_do;

        //std::map<INDEX_TYPE, INDEX_TYPE> bottleneckcount;
        std::map<INDEX_TYPE, std::vector<INDEX_TYPE> > bottleneckregions;

        int* m_thread_work_count;

        std::map<INDEX_TYPE, RegionState*> m_bottleneckpoint_to_region_map;
        void BeginIntegration() {


            TEST_OUTPUT = new DenseLabeling<int>(m_topo_grid->numCells());
            TEST_OUTPUT->SetAll(0);
            TEST_PRIOR_STUFF = new DenseLabeling<char>(m_topo_grid->numCells());
            TEST_PRIOR_STUFF->SetAll(0);

            //const INDEX_TYPE t_num_vertices = m_grid->NumElements();

            //AdvectionChecker* inside_voxel_critical_advection_checker = new TerminateNearCertain(m_certains, m_grid);
            //AdvectionChecker* no_check = new NoTermination();//AdvectionChecker* inside_voxel_advection_checker = new TerminateNearAssigned(m_destinations, m_grid);

            m_thread_work_count = new int[128];
            memset(m_thread_work_count, 0, 128 * sizeof(int));
            // set all potential critical_primal_edges, so we terminate near them
            std::vector<INDEX_TYPE> topo_index_partition;
            int num_threads;
#pragma omp parallel
            {
#pragma omp single
                {
                    num_threads = omp_get_num_threads();
                    ArrayIndexPartitioner::EvenChunkSplit(m_topo_grid->numCells(), num_threads, topo_index_partition);
                }

                int thread_num = omp_get_thread_num();
                TopologicalRegularGrid::DCellsIterator primal_edge_iterator(m_topo_grid, 1, topo_index_partition[thread_num], topo_index_partition[thread_num + 1]);
                for (primal_edge_iterator.begin(); primal_edge_iterator.valid(); primal_edge_iterator.advance()) {
                    INDEX_TYPE primal_edge_topo_id = primal_edge_iterator.value();

                    if (IsHighestPrimalEdgeInRestrictedNeighborhood(primal_edge_topo_id)) {
#pragma omp critical
                        {
                            m_work_to_do.push(WorkRequest{ EXPAND_ONLY, primal_edge_topo_id, NULL, NULL });
                            critical_primal_edges.push_back(primal_edge_topo_id);
                            TEST_OUTPUT->SetLabel(primal_edge_topo_id, thread_num + 1);
                        }
                    }
                }

            }
            printf("found %d critical_primal_edges!\nExpanding critical_primal_edges certain regions\n", critical_primal_edges.size());

            FILE* fmerge = fopen("merge.txt", "w");
            std::vector< std::pair<std::chrono::steady_clock::time_point, int> > timings;
            timings.push_back(std::pair<std::chrono::steady_clock::time_point, int>(std::chrono::steady_clock::now(), m_work_to_do.size()));
#pragma omp parallel
            {
                bool keep_looping = true;
                WorkRequest work;
                while (keep_looping) { // begin while
#pragma omp critical
                        { // begin critical
                            if (m_work_to_do.empty()) {
                                keep_looping = false;
                            }
                            else {
                                work = m_work_to_do.top();
                                m_work_to_do.pop();
                                timings.push_back(std::pair<std::chrono::steady_clock::time_point, int>(std::chrono::steady_clock::now(), m_work_to_do.size()));
                            }
                        } // end critical

                    if (keep_looping) {
                        if (work.worktype == EXPAND_ONLY) {

                            std::set<INDEX_TYPE> local_certains;														// debug
                            int thread_num = omp_get_thread_num();														// debug

                            RegionState* local_region = Expand_Lower_Neighborhood2(work.start_cell);

                            if (!local_region->unprocessed_band.empty()) {
#pragma omp critical
                                { // begin critical update to merge list and work to do
                                    if (m_bottleneckpoint_to_region_map.count(local_region->bottleneck_id) > 0) {
                                        RegionState* other_region = m_bottleneckpoint_to_region_map[local_region->bottleneck_id];
                                        m_bottleneckpoint_to_region_map.erase(local_region->bottleneck_id);
                                        m_work_to_do.push(WorkRequest{ MERGE_AND_EXPAND, 0, local_region, other_region });
                                        timings.push_back(std::pair<std::chrono::steady_clock::time_point, int>(std::chrono::steady_clock::now(), m_work_to_do.size()));

                                    }
                                    else {
                                        m_bottleneckpoint_to_region_map[local_region->bottleneck_id] = local_region;
                                    }

                                    // BEGIN DEBUG REGION
                                    std::priority_queue<std::pair<FLOATTYPE, INDEX_TYPE> > tmp_front(local_region->unprocessed_band);
                                    while (!tmp_front.empty()) {
                                        delayed_cells.push_back(tmp_front.top().second);
                                        tmp_front.pop();
                                    }
                                    // END DEBUG REGION

                                    //while (!local_region->unprocessed_band.empty()) {											// debug
                                    //	delayed_cells.push_back(local_region->unprocessed_band.top().second);					// debug
                                    //	local_region->unprocessed_band.pop();													// debug
                                    //}																						// debug
                                } // end critical update
                            }

                        }
                        else if (work.worktype == MERGE_AND_EXPAND) {
                            RegionState* merged = MergeRegions(work.region1, work.region2);
                            RegionState* modified = Expand_MergeRegions(merged);
                            // only add more work if there IS work to be added - could end at a max!
                            if (!modified->unprocessed_band.empty()){
#pragma omp critical
                                { // begin critical update to merge list and work to do
                                    if (m_bottleneckpoint_to_region_map.count(modified->bottleneck_id)) {
                                        RegionState* other_region = m_bottleneckpoint_to_region_map[modified->bottleneck_id];
                                        m_bottleneckpoint_to_region_map.erase(modified->bottleneck_id);
                                        if (m_work_to_do.size() < num_threads) {
                                            keep_looping = false;

                                        }
                                        else {
                                            m_work_to_do.push(WorkRequest{ MERGE_AND_EXPAND, 0, modified, other_region });
                                            timings.push_back(std::pair<std::chrono::steady_clock::time_point, int>(std::chrono::steady_clock::now(), m_work_to_do.size()));
                                            fprintf(fmerge, "pushing merge point %d %d %d %d %d\n", (int)modified->bottleneck_id, modified->narrow_band.size(), other_region->narrow_band.size(), m_work_to_do.size(), omp_get_thread_num());
                                        }
                                    }
                                    else {
                                        m_bottleneckpoint_to_region_map[modified->bottleneck_id] = modified;
                                    }

                                    // debug
                                    std::priority_queue<std::pair<FLOATTYPE, INDEX_TYPE> > tmp_front(modified->unprocessed_band);
                                    while (!tmp_front.empty()) {
                                        delayed_cells.push_back(tmp_front.top().second);
                                        tmp_front.pop();
                                    }
                                } // end critical update
                            } // end if delayed front not empty

                        }
                    } // end if(keep_looping)
                } // end while (keep_looping)
                timings.push_back(std::pair<std::chrono::steady_clock::time_point, int>(std::chrono::steady_clock::now(), m_work_to_do.size()));
            } // end parallel region

            for (int i = 0; i < num_threads; i++) {
                printf("thread %d did %d work\n", i, m_thread_work_count[i]);
            }

            int countassigned = 0;
            int countunassigned = 0;
            for (INDEX_TYPE i = 0; i < m_topo_grid->numCells(); i++) {

                if (m_topo_grid->dimension(i) == 1 && m_restriction_labels->GetLabel(i) == 1) {
                    if (TEST_OUTPUT->GetLabel(i) == 0) {
                        countunassigned++;
                        //printf("ERROR: edge that is part of boundar is not certain!\n");
                    }
                    else {
                        countassigned++;
                    }
                }
            }
            printf("%d assigned, %d unassigned\n", countassigned, countunassigned);



            for (int rounds = 0; rounds < 10; rounds++) {

#pragma omp parallel
            {

                int thread_num = omp_get_thread_num();
                TopologicalRegularGrid::DCellsIterator primal_edge_iterator(m_topo_grid, 1, topo_index_partition[thread_num], topo_index_partition[thread_num + 1]);
                for (primal_edge_iterator.begin(); primal_edge_iterator.valid(); primal_edge_iterator.advance()) {
                    INDEX_TYPE primal_edge_topo_id = primal_edge_iterator.value();

                    if (TEST_OUTPUT->GetLabel(primal_edge_topo_id) > 0) {
                        TEST_PRIOR_STUFF->SetLabel(primal_edge_topo_id, 1);
                    }
                }

            }

                critical_primal_edges.clear();
#pragma omp parallel
                {

                    int thread_num = omp_get_thread_num();
                    TopologicalRegularGrid::DCellsIterator primal_edge_iterator(m_topo_grid, 1, topo_index_partition[thread_num], topo_index_partition[thread_num + 1]);
                    for (primal_edge_iterator.begin(); primal_edge_iterator.valid(); primal_edge_iterator.advance()) {
                        INDEX_TYPE primal_edge_topo_id = primal_edge_iterator.value();

                        if (IsHighestPrimalEdgeInUNASSIGNEDRestrictedNeighborhood(primal_edge_topo_id)) {
#pragma omp critical
                        {
                            m_work_to_do.push(WorkRequest{ EXPAND_ONLY, primal_edge_topo_id, NULL, NULL });
                            //printf("found critical2\n");
                            critical_primal_edges.push_back(primal_edge_topo_id);
                            //TEST_OUTPUT->SetLabel(primal_edge_topo_id, thread_num + 1);
                        }
                        }
                    }

                }
                printf("found %d critical_primal_edges!\nExpanding critical_primal_edges certain regions\n", critical_primal_edges.size());







#pragma omp parallel
                {
                    bool keep_looping = true;
                    WorkRequest work;
                    while (keep_looping) { // begin while
#pragma omp critical
                        { // begin critical
                            if (m_work_to_do.empty()) {
                                keep_looping = false;
                            }
                            else {
                                work = m_work_to_do.top();
                                m_work_to_do.pop();
                                timings.push_back(std::pair<std::chrono::steady_clock::time_point, int>(std::chrono::steady_clock::now(), m_work_to_do.size()));
                            }
                        } // end critical

                        if (keep_looping) {
                            if (work.worktype == EXPAND_ONLY) {

                                std::set<INDEX_TYPE> local_certains;														// debug
                                int thread_num = omp_get_thread_num();														// debug

                                RegionState* local_region = Expand_Lower_Neighborhood2(work.start_cell);

                                if (!local_region->unprocessed_band.empty()) {
#pragma omp critical
                                { // begin critical update to merge list and work to do
                                    if (m_bottleneckpoint_to_region_map.count(local_region->bottleneck_id) > 0) {
                                        RegionState* other_region = m_bottleneckpoint_to_region_map[local_region->bottleneck_id];
                                        m_bottleneckpoint_to_region_map.erase(local_region->bottleneck_id);
                                        m_work_to_do.push(WorkRequest{ MERGE_AND_EXPAND, 0, local_region, other_region });
                                        timings.push_back(std::pair<std::chrono::steady_clock::time_point, int>(std::chrono::steady_clock::now(), m_work_to_do.size()));

                                    }
                                    else {
                                        m_bottleneckpoint_to_region_map[local_region->bottleneck_id] = local_region;
                                    }

                                    // BEGIN DEBUG REGION
                                    std::priority_queue<std::pair<FLOATTYPE, INDEX_TYPE> > tmp_front(local_region->unprocessed_band);
                                    while (!tmp_front.empty()) {
                                        delayed_cells.push_back(tmp_front.top().second);
                                        tmp_front.pop();
                                    }
                                    // END DEBUG REGION

                                    //while (!local_region->unprocessed_band.empty()) {											// debug
                                    //	delayed_cells.push_back(local_region->unprocessed_band.top().second);					// debug
                                    //	local_region->unprocessed_band.pop();													// debug
                                    //}																						// debug
                                } // end critical update
                                }

                            }
                            else if (work.worktype == MERGE_AND_EXPAND) {
                                RegionState* merged = MergeRegions(work.region1, work.region2);
                                RegionState* modified = Expand_MergeRegions(merged);
                                // only add more work if there IS work to be added - could end at a max!
                                if (!modified->unprocessed_band.empty()){
#pragma omp critical
                                { // begin critical update to merge list and work to do
                                    if (m_bottleneckpoint_to_region_map.count(modified->bottleneck_id)) {
                                        RegionState* other_region = m_bottleneckpoint_to_region_map[modified->bottleneck_id];
                                        m_bottleneckpoint_to_region_map.erase(modified->bottleneck_id);
                                        if (m_work_to_do.size() < num_threads) {
                                            keep_looping = false;

                                        }
                                        else {
                                            m_work_to_do.push(WorkRequest{ MERGE_AND_EXPAND, 0, modified, other_region });
                                            timings.push_back(std::pair<std::chrono::steady_clock::time_point, int>(std::chrono::steady_clock::now(), m_work_to_do.size()));
                                            fprintf(fmerge, "pushing merge point %d %d %d %d %d\n", (int)modified->bottleneck_id, modified->narrow_band.size(), other_region->narrow_band.size(), m_work_to_do.size(), omp_get_thread_num());
                                        }
                                    }
                                    else {
                                        m_bottleneckpoint_to_region_map[modified->bottleneck_id] = modified;
                                    }

                                    // debug
                                    std::priority_queue<std::pair<FLOATTYPE, INDEX_TYPE> > tmp_front(modified->unprocessed_band);
                                    while (!tmp_front.empty()) {
                                        delayed_cells.push_back(tmp_front.top().second);
                                        tmp_front.pop();
                                    }
                                } // end critical update
                                } // end if delayed front not empty

                            }
                        } // end if(keep_looping)
                    } // end while (keep_looping)
                    timings.push_back(std::pair<std::chrono::steady_clock::time_point, int>(std::chrono::steady_clock::now(), m_work_to_do.size()));
                } // end parallel region


                countassigned = 0;
                countunassigned = 0;
                for (INDEX_TYPE i = 0; i < m_topo_grid->numCells(); i++) {

                    if (m_topo_grid->dimension(i) == 1 && m_restriction_labels->GetLabel(i) == 1) {
                        if (TEST_OUTPUT->GetLabel(i) == 0) {
                            countunassigned++;
                            //printf("ERROR: edge that is part of boundar is not certain!\n");
                        }
                        else {
                            countassigned++;
                        }
                    }
                }
                printf("round%d: %d assigned, %d unassigned\n", rounds, countassigned, countunassigned);

                for (int i = 0; i < num_threads; i++) {
                    printf("thread %d did %d work\n", i, m_thread_work_count[i]);
                }

                if (countunassigned == 0) rounds = 100;
            }
            fclose(fmerge);
            ////			// THIS WILL ALL CHANGE -> one more indirection needed
            ////			std::map<INDEX_TYPE, INDEX_TYPE>* pmaps =
            ////				new std::map<INDEX_TYPE, INDEX_TYPE>[topo_index_partition.size() - 1];
            ////
            ////			INDEX_TYPE num_critical_primal_edges = critical_primal_edges.size();
            ////#pragma omp parallel for schedule(dynamic)
            ////			for (INDEX_TYPE i = 0; i < num_critical_primal_edges; i++) {
            ////				std::set<INDEX_TYPE> local_certains;
            ////				int thread_num = omp_get_thread_num();
            ////				RegionState* local_region = Expand_Lower_Neighborhood2(critical_primal_edges[i]);
            ////
            ////				for (auto it = local_region->processed_band.begin(); it != local_region->processed_band.end(); it++) {
            ////
            ////					pmaps[thread_num][*it] = critical_primal_edges[i];
            ////					TEST_OUTPUT->SetLabel(*it, TEST_OUTPUT->GetLabel(*it) + 1);
            ////				}
            ////#pragma omp critical
            ////				{
            ////					while (!local_region->unprocessed_band.empty()){
            ////						delayed_cells.push_back(local_region->unprocessed_band.top().second);
            ////						local_region->unprocessed_band.pop();
            ////					}
            ////				}
            ////			}
        }



    };




}
#endif
