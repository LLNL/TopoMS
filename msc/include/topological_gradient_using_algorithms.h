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

#ifndef TOPOLOGICAL_GRADIENT_USING_ALGORITHMS_H
#define TOPOLOGICAL_GRADIENT_USING_ALGORITHMS_H

#include "basic_types.h"
#include "discrete_gradient_labeling.h"
#include "topological_explicit_mesh_function.h"
#include "topological_regular_grid.h"
#include "topological_regular_masked_restricted_grid.h"

#include <vector>
#include <queue>
#include <map>
#include <set>




namespace MSC {

    template <class GridType, class MeshFunction, class LabelingType>
    class TopologicalGradientUsingAlgorithms {
    protected:

        MeshFunction* my_mesh_function;
        GridType* my_mesh_handler;
        LabelingType* my_grad_field;

    public:

        TopologicalGradientUsingAlgorithms(
            MeshFunction* mesh_function,
            GridType* mesh_handler,
            LabelingType* grad_field) :
            my_mesh_function(mesh_function),
            my_mesh_handler(mesh_handler),
            my_grad_field(grad_field) {
        }

        // trace "down" in gradient and fill in the result std::vector
        // with all cells that are found

        void CheckGradientConsistency() {
            typename GridType::AllCellsIterator allit(my_mesh_handler);
            for (allit.begin(); allit.valid(); allit.advance()) {
                INDEX_TYPE cellid = allit.value();
                if (my_grad_field->getAssigned(cellid) == 0) {
                    printf("CheckGradientConsistency(): error: cellid %d is not assigned\n",cellid);
                }
                if (my_grad_field->getCritical(cellid)) {
                }
                else {
                    INDEX_TYPE pairid = my_grad_field->getPair(cellid);
                    if (my_grad_field->getCritical(pairid)) {
                        printf("CheckGradientConsistency(): error: cell %d is paired with critical cell %d\n", cellid, pairid);
                    }
                    else {
                        INDEX_TYPE pairpair = my_grad_field->getPair(pairid);
                        if (pairpair != cellid) {
                            printf("CheckGradientConsistency(): error: pair pair is not cellid (%d -> %d -> %d)\n", cellid, pairid, pairpair);
                        }
                    }
                    if (my_mesh_handler->dimension(pairid) != my_mesh_handler->dimension(cellid) - 1 &&
                        my_mesh_handler->dimension(pairid) != my_mesh_handler->dimension(cellid) + 1) {
                        printf("CheckGradientConsistency(): error: dimensions of cell (%d) and pair (%d) dont match\n",
                               my_mesh_handler->dimension(cellid), my_mesh_handler->dimension(pairid));
                    }
                }
            }
        }


        virtual void count_critical_points(int dim) {
            int* counts = new int[dim];
            for (int i = 0; i < dim; i++) counts[i] = 0;

            for (INDEX_TYPE i = 0; i < my_mesh_handler->numCells(); i++) {
                if (my_grad_field->getCritical(i))
                    counts[my_mesh_handler->dimension(i)]++;
            }

            for (int i = 0; i < dim; i++)
                printf("index-%d=%d\n", i, counts[i]);
        }

        virtual void trace_down_cells(const INDEX_TYPE& cellid,
            std::vector<INDEX_TYPE>& result) {

            std::queue<INDEX_TYPE> cell_queue;
            cell_queue.push(cellid);

            result.clear();
            std::set<INDEX_TYPE> cell_visited;

            while (!cell_queue.empty()) {
                INDEX_TYPE current = cell_queue.front();
                cell_queue.pop();

                cell_visited.insert(current);
                result.push_back(current);

                typename GridType::FacetsIterator fit(my_mesh_handler);

                for (fit.begin(current); fit.valid(); fit.advance()) {
                    INDEX_TYPE temp_id = fit.value();

                    if (my_grad_field->getCritical(temp_id) &&
                        cell_visited.count(temp_id) == 0) {
                        result.push_back(temp_id);
                        cell_visited.insert(temp_id);
                    }
                    else if (cell_visited.count(temp_id) == 0) {
                        INDEX_TYPE pair = my_grad_field->getPair(temp_id);
                        result.push_back(temp_id);
                        result.push_back(pair);
                        cell_visited.insert(temp_id);
                        cell_visited.insert(pair);
                        cell_queue.push(pair);
                    }
                }
            }

        }

        virtual void trace_up_cells(const INDEX_TYPE& cellid,
            std::vector<INDEX_TYPE>& result) {

            std::queue<INDEX_TYPE> cell_queue;
            cell_queue.push(cellid);

            result.clear();
            std::set<INDEX_TYPE> cell_visited;

            while (!cell_queue.empty()) {
                INDEX_TYPE current = cell_queue.front();
                cell_queue.pop();

                cell_visited.insert(current);
                result.push_back(current);

                typename GridType::CofacetsIterator cofacets(my_mesh_handler);
                for (cofacets.begin(current); cofacets.valid(); cofacets.advance()) {
                    INDEX_TYPE temp_id = cofacets.value();

                    if (my_grad_field->getCritical(temp_id) &&
                        cell_visited.count(temp_id) == 0) {
                        result.push_back(temp_id);
                        cell_visited.insert(temp_id);
                    }
                    else if (cell_visited.count(temp_id) == 0) {
                        INDEX_TYPE pair = my_grad_field->getPair(temp_id);
                        result.push_back(temp_id);
                        result.push_back(pair);
                        cell_visited.insert(temp_id);
                        cell_visited.insert(pair);
                        cell_queue.push(pair);
                    }
                }
            }

        }


        virtual void trace_down_cells_restricted(const INDEX_TYPE& cellid,
            std::vector<INDEX_TYPE>& result) {

            std::queue<INDEX_TYPE> cell_queue;
            cell_queue.push(cellid);

            DIM_TYPE temp_dim = my_grad_field->getDimAscMan(cellid) + 1;
            result.clear();
            std::set<INDEX_TYPE> cell_visited;

            while (!cell_queue.empty()) {
                INDEX_TYPE current = cell_queue.front();
                cell_queue.pop();

                cell_visited.insert(current);
                result.push_back(current);

                typename GridType::FacetsIterator fit(my_mesh_handler);
                for (fit.begin(current); fit.valid(); fit.advance()) {
                    INDEX_TYPE temp_id = fit.value();

                    if (my_grad_field->getCritical(temp_id) &&
                        cell_visited.count(temp_id) == 0) {
                        result.push_back(temp_id);
                        cell_visited.insert(temp_id);
                    }
                    else if (cell_visited.count(temp_id) == 0 &&
                        my_grad_field->getDimAscMan(temp_id) == temp_dim) {
                        INDEX_TYPE pair = my_grad_field->getPair(temp_id);
                        result.push_back(temp_id);
                        result.push_back(pair);
                        cell_visited.insert(temp_id);
                        cell_visited.insert(pair);
                        cell_queue.push(pair);
                    }
                }
            }

        }
        virtual void trace_down_cells_restricted_counting(const INDEX_TYPE& cellid,
            std::vector<INDEX_TYPE>& result, std::vector<int>& counts) {

            std::queue<INDEX_TYPE> cell_queue;
            cell_queue.push(cellid);

            DIM_TYPE temp_dim = my_grad_field->getDimAscMan(cellid) + 1;
            result.clear();
            counts.clear();
            std::set<INDEX_TYPE> cell_visited;

            // build the graph
            std::map<INDEX_TYPE, std::set<INDEX_TYPE> > node_graph;
            std::map<INDEX_TYPE, int > visit_counts;

            while (!cell_queue.empty()) {
                INDEX_TYPE current = cell_queue.front();
                cell_queue.pop();

                std::set<INDEX_TYPE> neighbors;

                cell_visited.insert(current);

                typename GridType::FacetsIterator fit(my_mesh_handler);
                for (fit.begin(current); fit.valid(); fit.advance()) {
                    INDEX_TYPE temp_id = fit.value();

                    if (my_grad_field->getCritical(temp_id)) {
                        neighbors.insert(temp_id);
                        if (visit_counts.count(temp_id) == 0) {
                            visit_counts[temp_id] = 1;
                        }
                        else {
                            visit_counts[temp_id]++;
                        }

                        cell_visited.insert(temp_id);
                    }
                    else if (my_grad_field->getDimAscMan(temp_id) == temp_dim) {
                        INDEX_TYPE pair = my_grad_field->getPair(temp_id);
                        if (current == pair) continue;

                        neighbors.insert(pair);
                        if (visit_counts.count(pair) == 0) {
                            visit_counts[pair] = 1;
                        }
                        else {
                            visit_counts[pair]++;
                        }
                        if (cell_visited.count(pair) == 0) {
                            cell_queue.push(pair);
                        }
                        cell_visited.insert(temp_id);
                        cell_visited.insert(pair);

                    }
                }
                node_graph[current].insert(neighbors.begin(), neighbors.end());
            }
            //print graph
            printf("\ngraph of %d:\n", cellid);
            for (std::map<INDEX_TYPE, std::set<INDEX_TYPE> >::iterator mit = node_graph.begin();
                mit != node_graph.end(); mit++) {
                INDEX_TYPE tempid = (*mit).first;
                printf(" n=%d\n", tempid);
                for (std::set<INDEX_TYPE>::iterator sit = (*mit).second.begin();
                    sit != (*mit).second.end(); sit++)
                    printf("  -->%d\n", *sit);
            }
            // traverse graph from root
            cell_queue.push(cellid);
            while (!cell_queue.empty()) {
                INDEX_TYPE current = cell_queue.front();
                cell_queue.pop();
                result.push_back(current);
                counts.push_back(0);

                for (std::set<INDEX_TYPE>::iterator it = node_graph[current].begin();
                    it != node_graph[current].end(); it++) {
                    INDEX_TYPE tempid = *it;
                    visit_counts[tempid]--;
                    if (visit_counts[tempid] == 0) {
                        cell_queue.push(tempid);
                    }
                }
            }

            // the base case, 1 path from cell to itself
            visit_counts[cellid] = 1;
            for (int i = 0; i < result.size(); i++) {
                INDEX_TYPE current = result[i];
                int temp_count = visit_counts[current];
                counts[i] = temp_count;
                for (std::set<INDEX_TYPE>::iterator it = node_graph[current].begin();
                    it != node_graph[current].end(); it++) {
                    INDEX_TYPE tempid = *it;
                    visit_counts[tempid] += temp_count;
                }
            }
        }


        void rec_man_trace_up(INDEX_TYPE& cellid, std::set<INDEX_TYPE>& res) {
            res.insert(cellid);
            INDEX_TYPE current = cellid;
            DIM_TYPE cdim = this->my_mesh_handler->dimension(cellid);
            typename GridType::CofacetsIterator cofacets(my_mesh_handler);
            for (cofacets.begin(current); cofacets.valid(); cofacets.advance()) {
                INDEX_TYPE temp_id = cofacets.value();

                if (this->my_grad_field->getCritical(temp_id) || !my_grad_field->getAssigned(temp_id)) continue;

                INDEX_TYPE temp_pair = my_grad_field->getPair(temp_id);

                if (temp_pair == cellid) continue;

                if (my_mesh_handler->dimension(temp_pair) != cdim) continue;

                rec_man_trace_up(temp_pair, res);
            }
        }


protected:
        void rec_man_trace_up_marking(INDEX_TYPE& cellid, DIM_TYPE value) {
            if (my_grad_field->getDimAscMan(cellid) == value) return;
            my_grad_field->setDimAscMan(cellid, value);

            INDEX_TYPE current = cellid;
            DIM_TYPE cdim = this->my_mesh_handler->dimension(cellid);
            typename GridType::CofacetsIterator cofacets(my_mesh_handler);
            for (cofacets.begin(current); cofacets.valid(); cofacets.advance()) {
                INDEX_TYPE temp_id = cofacets.value();

                if (this->my_grad_field->getCritical(temp_id) || !my_grad_field->getAssigned(temp_id)) continue;

                INDEX_TYPE temp_pair = my_grad_field->getPair(temp_id);

                if (temp_pair == cellid) continue;

                if (my_mesh_handler->dimension(temp_pair) != cdim) continue;
                if (my_grad_field->getDimAscMan(temp_id) == value) return;
                my_grad_field->setDimAscMan(temp_id, value);

                rec_man_trace_up_marking(temp_pair, value);
            }
        }

public:


        void setAscendingManifoldDimensions() {

            std::vector<INDEX_TYPE> criticals[4];
            std::vector<INDEX_TYPE> topo_index_partition;
            int num_threads;
#pragma omp parallel
            {
#pragma omp single
                {
                    num_threads = omp_get_num_threads();
                    ArrayIndexPartitioner::EvenChunkSplit(my_mesh_handler->numCells(), num_threads, topo_index_partition);
                }

                int thread_num = omp_get_thread_num();
                typename GridType::AllCellsIterator all_cells_iterator(my_mesh_handler,  topo_index_partition[thread_num], topo_index_partition[thread_num + 1]);
                for (all_cells_iterator.begin(); all_cells_iterator.valid(); all_cells_iterator.advance()) {
                    INDEX_TYPE cell_id = all_cells_iterator.value();
                    //my_grad_field->setMark(cell_id, 0);
                    my_grad_field->setDimAscMan(cell_id, 3);
                    if (my_grad_field->getCritical(cell_id)) {
                        DIM_TYPE tdim = my_mesh_handler->dimension(cell_id);
#pragma omp critical
                        {
                            criticals[tdim].push_back(cell_id);
                        }
                    }
                }
            }
            // no now every cell is assigned to 3-manifold, and have list of critical points of each dimension
            //printf("found %d %d %d %d crits\n", criticals[0].size(), criticals[1].size(), criticals[2].size(), criticals[3].size());
            INDEX_TYPE num_1s = criticals[1].size();
//#pragma omp parallel for schedule(dynamic)
            for (INDEX_TYPE vid = 0; vid < num_1s; vid++) {
                INDEX_TYPE cid = criticals[1][vid];
                rec_man_trace_up_marking(cid, 2);
            }
            INDEX_TYPE num_2s = criticals[2].size();
//#pragma omp parallel for schedule(dynamic)
            for (INDEX_TYPE vid = 0; vid < num_2s; vid++) {
                INDEX_TYPE cid = criticals[2][vid];
                rec_man_trace_up_marking(cid, 1);
            }
            INDEX_TYPE num_3s = criticals[3].size();
//#pragma omp parallel for schedule(static)
            for (INDEX_TYPE vid = 0; vid < num_3s; vid++) {
                INDEX_TYPE cid = criticals[3][vid];
                my_grad_field->setDimAscMan(cid, 0);
            }
        }

    };
}
#endif
