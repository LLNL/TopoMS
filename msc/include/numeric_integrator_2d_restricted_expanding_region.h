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

#ifndef NUMERIC_INTEGRATOR_2D_RESTRICTED_EXPANDING_REGION_H
#define NUMERIC_INTEGRATOR_2D_RESTRICTED_EXPANDING_REGION_H

//#include <algorithm>
#include <stdio.h>
#include <set>
#include <queue>
#include "basic_types.h"
#include "vectors.h"
#include "labeling.h"
#include "regular_grid.h"
#include "regular_grid_trilinear_function.h"
#include "adaptive_euler_advector.h"
#include "topological_regular_grid.h"
#include "array_index_partition.h"
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
    class NumericIntegrator2DRestrictedExpandingRegion {

    protected:

        // will need some other kind of tuple and not dense?
        DenseLabeling<char>* m_restriction_labels; // 1 if we can, 0 if we can;t
        RegularGrid* m_grid;
        RegularGridTrilinearFunction* m_func;
        TopologicalRegularGrid* m_topo_grid;
        static const Vec3l kNeighborOffsets[3][4];

        Vec3i m_xyz;
        Vec3b m_periodic;
        FLOATTYPE m_error_threshold;
        FLOATTYPE m_gradient_threshold;
        int m_iteration_limit;

        char* m_fname;

        // this is where we can change how edges get their values. - just interpolate for now?
        inline FLOATTYPE PrimalEdgeValue(INDEX_TYPE primal_edge_topo_id) const {
            // this will be slow but whatever! i do what i want!
            TopologicalRegularGrid::FacetsIterator t_adjacent_primal_vertex_iterator(m_topo_grid);
            int t_vertex_count = 0;
            FLOATTYPE t_value_sum = 0.0;
            for (t_adjacent_primal_vertex_iterator.begin(primal_edge_topo_id);
                t_adjacent_primal_vertex_iterator.valid();
                t_adjacent_primal_vertex_iterator.advance()) {
                INDEX_TYPE t_primal_vertex_topo_id = t_adjacent_primal_vertex_iterator.value();
                INDEX_TYPE t_primal_vertex_grid_id = m_topo_grid->VertexNumberFromCellID(t_primal_vertex_topo_id);
                t_value_sum += m_func->SampleImage(t_primal_vertex_grid_id);
                t_vertex_count++;
            }
            return t_value_sum / ((FLOATTYPE)t_vertex_count);
        }

        // consistent edge comparisons
        inline bool IsPrimalEdgeGreater(INDEX_TYPE a, FLOATTYPE a_value, INDEX_TYPE b, FLOATTYPE b_value) const {
            if (a_value != b_value) return a_value > b_value;
            return a > b;
        }

        bool IsHighestPrimalEdgeInRestrictedNeighborhood(INDEX_TYPE primal_edge_topo_id) const {

            // return false if the edge is not part of our permissable set
            if (m_restriction_labels->GetLabel(primal_edge_topo_id) != 1) {
                return false;
            }

            FLOATTYPE t_original_edge_value = PrimalEdgeValue(primal_edge_topo_id);
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
                    FLOATTYPE t_npe_value = PrimalEdgeValue(t_neighboring_primal_edge_topo_id);
                    if (IsPrimalEdgeGreater(t_neighboring_primal_edge_topo_id, t_npe_value,
                        primal_edge_topo_id, t_original_edge_value)) return false;
                }
            }
            return true;
        }

        void Enqueue_Lower_Neighbors(INDEX_TYPE primal_edge_topo_id,
            FLOATTYPE t_original_edge_value,
            std::priority_queue<std::pair<FLOATTYPE, INDEX_TYPE> > &expansion,
            std::set<INDEX_TYPE>&seen) const {
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
                        seen.count(t_neighboring_primal_edge_topo_id) == 0) {
                        // i think we can ignore the value test? - no we can't don't want to "spill over"
                        // into another basin - only insert strictly lower stuff!
                        FLOATTYPE t_neg_value = PrimalEdgeValue(t_neighboring_primal_edge_topo_id);
                        if (IsPrimalEdgeGreater(primal_edge_topo_id, t_original_edge_value,
                            t_neighboring_primal_edge_topo_id, t_neg_value)) {
                            seen.insert(t_neighboring_primal_edge_topo_id);
                            expansion.push(std::pair<FLOATTYPE, INDEX_TYPE>(t_neg_value, t_neighboring_primal_edge_topo_id));
                        }
                    }

                }
            }
        }

        // returns true if all the higher neighbors are "certain"
        bool Inspect_Higher_Certains(INDEX_TYPE primal_edge_topo_id, FLOATTYPE original_edge_value, std::set<INDEX_TYPE>& certains) const {
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
                    if (m_restriction_labels->GetLabel(t_neighboring_primal_edge_topo_id) != 1) continue;

                    FLOATTYPE t_npe_value = PrimalEdgeValue(t_neighboring_primal_edge_topo_id);
                    // if our neighbor is highter, check it is certain
                    if (IsPrimalEdgeGreater(t_neighboring_primal_edge_topo_id, t_npe_value,
                        primal_edge_topo_id, original_edge_value)) {
                        if (certains.count(t_neighboring_primal_edge_topo_id) == 0) {
                            // then its not certain so we can't be certain
                            return false;
                        }
                    }
                }
            }
            return true;
        }

        void Expand_Lower_Neighborhood(INDEX_TYPE primal_edge_topo_id, std::set<INDEX_TYPE>& certains) const {

            certains.clear();
            std::set<INDEX_TYPE> seen;
            INDEX_TYPE t_id = primal_edge_topo_id;
            // the natural ordering using the < operator on pairs will give us the highest
            // element first, simulating region growing from high to low
            std::priority_queue<std::pair<FLOATTYPE, INDEX_TYPE> > growing_front;
            seen.insert(primal_edge_topo_id);
            certains.insert(primal_edge_topo_id);
            FLOATTYPE original_edge_value = PrimalEdgeValue(primal_edge_topo_id);
            Enqueue_Lower_Neighbors(primal_edge_topo_id, original_edge_value, growing_front, seen);

            while (!growing_front.empty()) {

                INDEX_TYPE currid = growing_front.top().second;
                FLOATTYPE currid_value = growing_front.top().first;
                growing_front.pop();

                // cellvalue >=0 indicates that there is certainty here, so lets expand
                if (Inspect_Higher_Certains(currid, currid_value, certains)) {
                    certains.insert(currid);
                    Enqueue_Lower_Neighbors(currid, currid_value, growing_front, seen);
                }

            }
        }

        // only works on the primal grid
        // the corners gives the
        // return value is success or failure
        inline bool GetDualQuadFromPrimalEdgeWithGrad(INDEX_TYPE primal_edge_topo_id, int& proj_axis,
            Vec3d* corners_grid_coords, Vec2d* grads_2d) const {


            Vec3l ec;
            m_topo_grid->cellid2Coords(primal_edge_topo_id, ec);
            for (int i = 0; i < 3; i++) if (ec[i] % 2 == 1) proj_axis = i; // gets the odd occrdinate

            if (m_topo_grid->MemoryBlockBoundaryValue(ec)) { return false; } // ignore domain boundaries for now

            //int num_quads = 0;
            //INDEX_TYPE quadlist[4];
            //TopologicalRegularGrid::CofacetsIterator primal_quads_iter(m_topo_grid);
            //for (primal_quads_iter.begin(primal_edge_topo_id);
            //	primal_quads_iter.valid();
            //	primal_quads_iter.advance()) {
            //	INDEX_TYPE primal_quad = primal_quads_iter.value();
            //	quadlist[num_quads++] = primal_quad;
            //}

            //// ignore boundary cases for now
            //if (num_quads == 2) return;

            // er ignore quadlist, and simply compute voxels
            // (x, y, [z]) -> (x-1, y-1,) (x+1, y-1, )...

            for (int i = 0; i < 4; i++) {
                corners_grid_coords[i] = ((Vec3d) ec + kNeighborOffsets[proj_axis][i])*0.5;
            }

            for (int i = 0; i < 4; i++) {
                //Vec3d corners_grid_coords = corners_topo_coords[i];
                //corners_grid_coords *= 0.5; // grid coordinates of func
                Vec3d t_grad = m_func->TriLinInterpGrad(corners_grid_coords[i]);
                grads_2d[i] = t_grad.subset(proj_axis); // this projects the 3d gradient to the 2d dual quad plane
            }
            return true;
        }

        inline Vec2d Put3dPointOntoDualQuad(Vec3d& point_grid_coords, int proj_axis, Vec3d* quad_corners_grid_coords) const {
            //it is assumed that the point is co-planar (pretty much) with the quad
            // also the FIRST dual quad corner in the list is always the one with the LOWEST
            // coordinate - EVEN in the case where it is periodic and we wrap around - this
            // is the one to use to subtract!
            Vec2d point_proj = point_grid_coords.subset(proj_axis);
            Vec2d base_proj = quad_corners_grid_coords[0].subset(proj_axis);
            Vec2d parameterized_point = point_proj - base_proj; // if this is not between 0-1 we screwed up somewhere - but maybe not!
            // now FORCE it to be inside
            if (parameterized_point[0] <= 0) {
                parameterized_point[0] = 0.00001; // my epsilon
            }
            else if (parameterized_point[0] >= 1.0) {
                parameterized_point[0] = 1 - 0.00001; // my epsilon
            }

            if (parameterized_point[1] <= 0) {
                parameterized_point[1] = 0.00001; // my epsilon
            }
            else if (parameterized_point[1] >= 1.0) {
                parameterized_point[1] = 1 - 0.00001; // my epsilon
            }
            return parameterized_point;
        }


        // the current_primal_edge_topo_id will become the quad that we advect through
        void RecursiveIntegrateUpwardsFromPoint(INDEX_TYPE current_primal_edge_topo_id, Vec3d current_point_grid_coords,
            std::set<INDEX_TYPE>& seen_dual_quads, std::set<INDEX_TYPE>& destinations, int iterations_left, AdaptiveInQuadEulerAdvector* advector2d,
            std::vector<Vec3d>& line_soup /* LINE_SOUP */){

            // mark that this quad has been visited
            seen_dual_quads.insert(current_primal_edge_topo_id);
            Vec3l primal_edge_topo_coords;
            m_topo_grid->cellid2Coords(current_primal_edge_topo_id, primal_edge_topo_coords);


            line_soup.push_back(current_point_grid_coords); // LINE SOUP

            int proj_axis;
            Vec3d corners_grid_coords[4];
            Vec2d grads_2d[4];

            // set up problem
            int temp_iter_left = m_iteration_limit;

            if (!GetDualQuadFromPrimalEdgeWithGrad(current_primal_edge_topo_id, proj_axis, corners_grid_coords, grads_2d)){
                // the we failed, might as well continue
                printf("put on quad failed - boundary\n");
                line_soup.push_back(current_point_grid_coords);
                return;
            }

            Vec2d current_2d_point = Put3dPointOntoDualQuad(current_point_grid_coords, proj_axis, corners_grid_coords);

            // DELETE THIS WHole test!! ONCE LINE_SOUP is not used
            Vec3d newpoint_test = corners_grid_coords[0];
            int tmp_test = 0;
            for (int i = 0; i < 3; i++) {
                if (i == proj_axis) continue;
                newpoint_test[i] += current_2d_point[tmp_test++];
            }
            line_soup.push_back(newpoint_test); // LINE_SOUP
            line_soup.push_back(newpoint_test); // LINE_SOUP
            // END PART TO delete LINE SOUP

            AdaptiveInQuadEulerAdvector::PLANEKIND exitplane;
            ADVECTION_EVENT event = advector2d->AdvectThroughQuadNoCheck(current_2d_point, grads_2d, temp_iter_left, exitplane);


            // MOVE THIS BACK BELOW EXIT to reinflate!! ONCE LINE_SOUP is not used
            Vec3d newpoint = corners_grid_coords[0];
            int tmp = 0;
            for (int i = 0; i < 3; i++) {
                if (i == proj_axis) continue;
                newpoint[i] += current_2d_point[tmp++];
            }
            line_soup.push_back(newpoint); // LINE_SOUP
            // END PART TO PUT DOWN LINE SOUP

            if (exitplane == AdaptiveInQuadEulerAdvector::NONE ||
                exitplane == AdaptiveInQuadEulerAdvector::PLANE_ERROR) {
                // return
                printf("ended in kind = %d\n", exitplane);

                return;
            }

            // reinflate!!


            Vec3l next_primal_quad_quad_coords = primal_edge_topo_coords;
            // get the primal quad
            if (exitplane == AdaptiveInQuadEulerAdvector::X0PLANE) {
                // then first of 2d coordinates is < 0
                int pos_of_x = (proj_axis == 0 ? 1 : 0);
                next_primal_quad_quad_coords[pos_of_x] -= 1;
            }
            else if (exitplane == AdaptiveInQuadEulerAdvector::X1PLANE) {
                // then first of 2d coordinates is < 0
                int pos_of_x = (proj_axis == 0 ? 1 : 0);
                next_primal_quad_quad_coords[pos_of_x] += 1;
            }
            else if (exitplane == AdaptiveInQuadEulerAdvector::Y0PLANE) {
                // then first of 2d coordinates is < 0
                int pos_of_y = (proj_axis <= 1 ? 2 : 1);
                next_primal_quad_quad_coords[pos_of_y] -= 1;
            }
            else if (exitplane == AdaptiveInQuadEulerAdvector::Y1PLANE) {
                // then first of 2d coordinates is < 0
                int pos_of_y = (proj_axis <= 1 ? 2 : 1);
                next_primal_quad_quad_coords[pos_of_y] += 1;
            }

            //now we have the quad!
            INDEX_TYPE next_primal_quad_topo_id = m_topo_grid->coords2Cellid(next_primal_quad_quad_coords);

            seen_dual_quads.insert(next_primal_quad_topo_id);

            // start expanding on all downwards things
            //Vec3d next_quad_grid_coords = ((Vec3d)next_primal_quad_quad_coords) * 0.5;
            //Vec3d next_quad_grad = m_func->TriLinInterpGrad(next_quad_grid_coords);
            int countsplits = 0;
            TopologicalRegularGrid::FacetsIterator next_edge_iter(m_topo_grid);
            for (next_edge_iter.begin(next_primal_quad_topo_id); next_edge_iter.valid(); next_edge_iter.advance()) {

                INDEX_TYPE next_primal_edge_topo_id = next_edge_iter.value();
                if (next_primal_edge_topo_id == current_primal_edge_topo_id) continue; // don't go back to same edge
                if (seen_dual_quads.count(next_primal_edge_topo_id) != 0) continue; // it's seen so exit
                if (m_restriction_labels->GetLabel(next_primal_edge_topo_id) != 1) continue; // stay in restriction
                // TEST CRITICALITY


                //Vec3l next_primal_edge_topo_coords;
                //m_topo_grid->cellid2Coords(next_primal_edge_topo_id, next_primal_edge_topo_coords);

                // recurse on ALL

                //Vec3d next_edge_grid_coords = ((Vec3d)next_primal_edge_topo_coords) * 0.5;

                //if (next_quad_grad.Dot(next_edge_grid_coords - next_quad_grid_coords) > 0) {
                    // // then the next edge is in the direction of the gradient so lets go!. NONONONON
                RecursiveIntegrateUpwardsFromPoint(next_primal_edge_topo_id, newpoint, seen_dual_quads, destinations, temp_iter_left, advector2d, line_soup /* LINE_SOUP */);
                //}
                    countsplits++;


            }

            if (countsplits > 1)
                printf("countsplit = %d\n", countsplits);



        }


    public:


        NumericIntegrator2DRestrictedExpandingRegion(
            RegularGridTrilinearFunction* func,
            RegularGrid* grid,
            TopologicalRegularGrid* topo_grid,
            DenseLabeling<char>* restriction_labels,
            FLOATTYPE error_threshold, FLOATTYPE gradient_threshold, int interation_limit) :
            m_xyz(func->GetGrid()->XYZ()), m_periodic(func->GetGrid()->Periodic()), m_func(func), m_grid(grid), m_restriction_labels(restriction_labels),
            m_topo_grid(topo_grid), m_gradient_threshold(gradient_threshold), m_error_threshold(error_threshold), m_iteration_limit(interation_limit) {
        }

        //DenseLabeling<INDEX_TYPE>* GetOutputLabels() { return m_destinations; }

        RegularGrid* GetGrid() { return m_grid; }
        RegularGridTrilinearFunction* GetFunction() { return m_func; }


        DenseLabeling<char>* TEST_OUTPUT;
        DenseLabeling<int>* TEST_LINEOUT;

        void BeginIntegration() {


            TEST_OUTPUT = new DenseLabeling<char>(m_topo_grid->numCells());
            TEST_OUTPUT->SetAll(0);
            TEST_LINEOUT = new DenseLabeling<int>(m_topo_grid->numCells());
            TEST_LINEOUT->SetAll(0);
            //const INDEX_TYPE t_num_vertices = m_grid->NumElements();

            //AdvectionChecker* inside_voxel_critical_advection_checker = new TerminateNearCertain(m_certains, m_grid);
            //AdvectionChecker* no_check = new NoTermination();//AdvectionChecker* inside_voxel_advection_checker = new TerminateNearAssigned(m_destinations, m_grid);

            std::vector<INDEX_TYPE> critical_primal_edges;

            // set all potential critical_primal_edges, so we terminate near them
            std::vector<INDEX_TYPE> topo_index_partition;

#pragma omp parallel
            {
#pragma omp single
                {
                    int num_threads = omp_get_num_threads();
                    ArrayIndexPartitioner::EvenChunkSplit(m_topo_grid->numCells(), num_threads, topo_index_partition);
                }

                int thread_num = omp_get_thread_num();
                TopologicalRegularGrid::DCellsIterator primal_edge_iterator(m_topo_grid, 1, topo_index_partition[thread_num], topo_index_partition[thread_num + 1]);
                for (primal_edge_iterator.begin(); primal_edge_iterator.valid(); primal_edge_iterator.advance()) {
                    INDEX_TYPE primal_edge_topo_id = primal_edge_iterator.value();

                    if (IsHighestPrimalEdgeInRestrictedNeighborhood(primal_edge_topo_id)) {
#pragma omp critical
                        {
                            critical_primal_edges.push_back(primal_edge_topo_id);
                            TEST_OUTPUT->SetLabel(primal_edge_topo_id, 1);
                        }
                    }
                }

            }
            printf("found %d critical_primal_edges!\nExpanding critical_primal_edges certain regions\n", critical_primal_edges.size());


            // THIS WILL ALL CHANGE -> one more indirection needed
            std::map<INDEX_TYPE, INDEX_TYPE>* pmaps =
                new std::map<INDEX_TYPE, INDEX_TYPE>[topo_index_partition.size() - 1];

            INDEX_TYPE num_critical_primal_edges = critical_primal_edges.size();
#pragma omp parallel for schedule(dynamic)
            for (INDEX_TYPE i = 0; i < num_critical_primal_edges; i++) {
                std::set<INDEX_TYPE> local_certains;
                int thread_num = omp_get_thread_num();
                Expand_Lower_Neighborhood(critical_primal_edges[i], local_certains);

                for (auto it = local_certains.begin(); it != local_certains.end(); it++) {

                    pmaps[thread_num][*it] = critical_primal_edges[i];
                    TEST_OUTPUT->SetLabel(*it, TEST_OUTPUT->GetLabel(*it) + 1);
                }
            }

            NoTermination2D* no_op = new NoTermination2D();

            FILE* fout = fopen("Linesout.txt", "w");

#pragma omp parallel
            {
                int thread_num = omp_get_thread_num();
                AdaptiveInQuadEulerAdvector* advector2d = new AdaptiveInQuadEulerAdvector(m_grid, m_func, m_gradient_threshold, m_error_threshold, no_op);
                TopologicalRegularGrid::DCellsIterator primal_edge_iterator(m_topo_grid, 1, topo_index_partition[thread_num], topo_index_partition[thread_num + 1]);
                for (primal_edge_iterator.begin(); primal_edge_iterator.valid(); primal_edge_iterator.advance()) {
                    INDEX_TYPE primal_edge_topo_id = primal_edge_iterator.value();

                    //skip anything that is not a boundary edge or is certain
                    if (m_restriction_labels->GetLabel(primal_edge_topo_id) != 1 ||
                        pmaps[thread_num].count(primal_edge_topo_id) != 0) continue;



                    INDEX_TYPE current_primal_edge_topo_id = primal_edge_topo_id;
                    Vec3l primal_edge_topo_coords;
                    m_topo_grid->cellid2Coords(primal_edge_topo_id, primal_edge_topo_coords);
                    Vec3d current_point_grid_coords = ((Vec3d)primal_edge_topo_coords) * 0.5;

                    std::set<INDEX_TYPE> seen_quads;
                    std::set<INDEX_TYPE> destinations;
                    std::vector<Vec3d> line_soup; // LINE_SOUP
                    RecursiveIntegrateUpwardsFromPoint(primal_edge_topo_id, current_point_grid_coords,
                        seen_quads, destinations, m_iteration_limit, advector2d, line_soup /*LINE_SOUP*/);

#pragma omp critical
                    {
                            int val = primal_edge_topo_id;
                        for (auto it = seen_quads.begin(); it != seen_quads.end(); it++) {
                            TEST_LINEOUT->SetLabel(*it, val);

                        }
                        // LINE SOUP
                        fprintf(fout, "%d %d\n", val, line_soup.size());
                        for (int i = 0; i < line_soup.size(); i++) {
                            fprintf(fout, "%f %f %f\n", line_soup[i][0], line_soup[i][1], line_soup[i][2]);
                        }
                    }
                    // start a downward integration





                }
            }


            fclose(fout); //LINE_SOUP








            // PURELY TESTING
#pragma omp parallel
            {
                int thread_num = omp_get_thread_num();
                TopologicalRegularGrid::DCellsIterator edges(m_topo_grid, 2, topo_index_partition[thread_num], topo_index_partition[thread_num + 1]);
                for (edges.begin(); edges.valid(); edges.advance()) {
                    INDEX_TYPE quad = edges.value();
                    if (m_restriction_labels->GetLabel(quad) != 1) continue;
                    TopologicalRegularGrid::FacetsIterator eit(m_topo_grid);
                    bool isgood = true;
                    for (eit.begin(quad); eit.valid(); eit.advance()) {
                        INDEX_TYPE e = eit.value();
                        if (m_restriction_labels->GetLabel(e) == 1 && TEST_OUTPUT->GetLabel(e) == 0) isgood = false;
                    }
                    if (isgood) TEST_OUTPUT->SetLabel(quad, 1);
                }
            }
#pragma omp parallel
            {
                int thread_num = omp_get_thread_num();
                TopologicalRegularGrid::DCellsIterator edges(m_topo_grid, 3, topo_index_partition[thread_num], topo_index_partition[thread_num + 1]);
                for (edges.begin(); edges.valid(); edges.advance()) {
                    INDEX_TYPE quad = edges.value();
                    if (m_restriction_labels->GetLabel(quad) != 1) continue;
                    TopologicalRegularGrid::FacetsIterator eit(m_topo_grid);
                    bool isgood = true;
                    for (eit.begin(quad); eit.valid(); eit.advance()) {
                        INDEX_TYPE e = eit.value();
                        if (m_restriction_labels->GetLabel(e) == 1 && TEST_OUTPUT->GetLabel(e) == 0) isgood = false;
                    }
                    if (isgood) TEST_OUTPUT->SetLabel(quad, 1);
                }
            }

            INDEX_TYPE ttnun = m_topo_grid->numCells();
#pragma omp parallel for schedule(static)
            for (INDEX_TYPE i = 0; i < ttnun; i++) {
                TEST_OUTPUT->SetLabel(i, m_restriction_labels->GetLabel(i) + TEST_OUTPUT->GetLabel(i));
            }

            for (INDEX_TYPE i = 0; i < num_critical_primal_edges; i++) {
                INDEX_TYPE ii = critical_primal_edges[i];
                TEST_OUTPUT->SetLabel(ii, TEST_OUTPUT->GetLabel(ii) + 2);

            }
        }


        //			int num_maxima = critical_primal_edges.size();
        //#pragma omp parallel shared(critical_primal_edges)
        //			{
        //#pragma omp for schedule(dynamic)  nowait
        //				for (int m = 0; m < num_maxima; m++) {
        //					INDEX_TYPE maximum = critical_primal_edges[m];
        //					m_certains->SetLabel(maximum, m);
        //					Expand_Lower_Neighborhood(maximum);
        //				}
        //			}
        //
        //
        //			int t1, t2; t1 = t2 = 0;
        //#pragma omp parallel
        //	{
        //				AdaptiveEulerAdvector t_advector(m_grid, m_func, m_gradient_threshold, m_error_threshold, inside_voxel_critical_advection_checker);
        //			std::vector<INDEX_TYPE> t_path;
        //			t_path.reserve(100);
        //#pragma omp for schedule(guided)  nowait
        //			for (INDEX_TYPE i = 0; i < t_num_vertices; i++) {
        //				// early skip if this is already a maximum
        //				if (m_certains->GetLabel(i) != -1 || m_destinations->GetLabel(i) != -1) {
        //					continue;
        //				}
        //
        //				t_path.clear();
        //				t_path.push_back(i);
        //				Vec3l t_coords = m_grid->XYZ3d(i); // get the coordinates of the poitn
        //				if (t_coords[0] == 0 && t_coords[1] == 0) printf("doing %d\n", t_coords[2]);
        //				Vec3d t_current_point = t_coords;
        //				int t_num_iterations_left = 500;
        //				bool t_continue = true;
        //
        //				while (t_continue) {
        //					Vec3d t_next_point;
        //					ADVECTION_EVENT t_return_code;
        //					if (m_grid->DistToBoundary(t_coords) <= 1) {
        //						t_return_code = t_advector.AdvectThroughVoxelNearBoundary(t_current_point, t_num_iterations_left);
        //						t_coords = m_grid->Inbounds(t_current_point + 0.5); // get nearest integer voxel
        //						t1++;
        //					}
        //					else {
        //						t_return_code = t_advector.AdvectThroughVoxelNoCheck(t_current_point, t_num_iterations_left);
        //						t_coords = (t_current_point + 0.5);
        //						t2++;
        //					}
        //					INDEX_TYPE t_next_id = m_grid->Index3d(t_coords);
        //					t_path.push_back(t_next_id);
        //					// if we terminated or hit a critical point, then we are done
        //					if (t_return_code == ADVECTION_EVENT::LOW_GRADIENT ||
        //						t_return_code == ADVECTION_EVENT::HIT_EXTREMUM ||
        //						t_return_code == ADVECTION_EVENT::HIT_PREASSIGNED ||
        //						t_return_code == ADVECTION_EVENT::OVER_MAX_ITERATIONS) {
        //						INDEX_TYPE t_dest_label = m_destinations->GetLabel(t_next_id);
        //						int t_certain_label = m_certains->GetLabel(t_next_id);
        //						for (int j = 0; j < t_path.size(); j++) {
        //							m_destinations->SetLabel(t_path[j], t_dest_label);
        //							//m_certains->SetLabel(t_path[j], t_certain_label);
        //						}
        //						t_continue = false;
        //					}
        //				}
        //			}
        //#pragma omp for schedule(static)
        //			for (INDEX_TYPE i = 0; i < t_num_vertices; i++) {
        //				// early skip if this is already a maximum
        //				if (m_certains->GetLabel(i) != -1) {
        //					continue;
        //				}
        //
        //				Vec3l t_coords = m_grid->XYZ3d(i); // get the coordinates of the poitn
        //				if (t_coords[0] == 0 && t_coords[1] == 0) printf("doing %d\n", t_coords[2]);
        //				Vec3l negs[6];
        //
        //				int nn = m_grid->GatherExistingNeighbors6(t_coords, negs);
        //
        //				bool needsupdate = false;
        //				INDEX_TYPE t_dest_label = m_destinations->GetLabel(i);
        //				for (int j = 0; j < nn; j++)  {
        //					if (m_destinations->GetLabel(m_grid->Index3d(negs[j])) != t_dest_label) {
        //						needsupdate = true;
        //						break;
        //					}
        //				}
        //
        //				if (needsupdate) continue;
        //
        //				m_certains->SetLabel(i, 1);
        //			}
        //#pragma omp for schedule(guided)  nowait
        //				for (INDEX_TYPE i = 0; i < t_num_vertices; i++) {
        //					// early skip if this is already a maximum
        //					if (m_certains->GetLabel(i) != -1) {
        //						continue;
        //					}
        //
        //					Vec3l t_coords = m_grid->XYZ3d(i); // get the coordinates of the poitn
        //					if (t_coords[0] == 0 && t_coords[1] == 0) printf("doing %d\n", t_coords[2]);
        //
        //
        //
        //					Vec3d t_current_point = t_coords;
        //					int t_num_iterations_left = 500;
        //					bool t_continue = true;
        //
        //					while (t_continue) {
        //						Vec3d t_next_point;
        //						ADVECTION_EVENT t_return_code;
        //						if (m_grid->DistToBoundary(t_coords) <= 1) {
        //							t_return_code = t_advector.AdvectThroughVoxelNearBoundary(t_current_point, t_num_iterations_left);
        //							t_coords = m_grid->Inbounds(t_current_point + 0.5); // get nearest integer voxel
        //							t1++;
        //						}
        //						else {
        //							t_return_code = t_advector.AdvectThroughVoxelNoCheck(t_current_point, t_num_iterations_left);
        //							t_coords = (t_current_point + 0.5);
        //							t2++;
        //						}
        //						INDEX_TYPE t_next_id = m_grid->Index3d(t_coords);
        //						t_path.push_back(t_next_id);
        //						// if we terminated or hit a critical point, then we are done
        //						if (t_return_code == ADVECTION_EVENT::LOW_GRADIENT ||
        //							t_return_code == ADVECTION_EVENT::HIT_EXTREMUM ||
        //							t_return_code == ADVECTION_EVENT::HIT_PREASSIGNED ||
        //							t_return_code == ADVECTION_EVENT::OVER_MAX_ITERATIONS) {
        //							m_destinations->SetLabel(i, m_destinations->GetLabel(t_next_id));
        //							t_continue = false;
        //						}
        //					}
        //				}
        //
        //				//if (m_grid->DistToBoundary()) {}
        //				//if (IsLowestVertexIn6Neighborhood(i)) {
        //				//	m_destinations->SetLabel(i, i);
        //				//}
        //				//else {
        //				//	m_destinations->SetLabel(i, -1);
        //				//}
        //
        //
        //
        //			printf("%d boundary, %d noboundary\n", t1, t2);
        //			}

        //	}



    };

/*
    const Vec3l NumericIntegrator2DRestrictedExpandingRegion::kNeighborOffsets[3][4] = {
            { Vec3l(0, -1, -1), Vec3l(0, 1, -1), Vec3l(0, -1, 1), Vec3l(0, 1, 1) },
            { Vec3l(-1, 0, -1), Vec3l(1, 0, -1), Vec3l(-1, 0, 1), Vec3l(1, 0, 1) },
            { Vec3l(-1, -1, 0), Vec3l(1, -1, 0), Vec3l(-1, 1, 0), Vec3l(1, 1, 0) }
    };
*/
}
#endif
