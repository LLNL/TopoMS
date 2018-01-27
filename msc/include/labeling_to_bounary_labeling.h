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

#ifndef VERTEX_LABELING_TO_BOUNDARY_LABELING_H
#define VERTEX_LABELING_TO_BOUNDARY_LABELING_H

#include <omp.h>
#include "basic_types.h"
#include "regular_grid.h"
#include "topological_regular_grid.h"
#include "labeling.h"
#include "array_index_partition.h"

namespace MSC {




    template<typename INPUT_LABEL_TYPE>
    class VertexLabelingToBoundaryLabeling {
    protected:
        DenseLabeling<INPUT_LABEL_TYPE>* m_input_labels;
        DenseLabeling<char>* m_output_labels;

        TopologicalRegularGrid* m_topological_grid;
    public:

        VertexLabelingToBoundaryLabeling(DenseLabeling<INPUT_LABEL_TYPE>* input_labels, TopologicalRegularGrid* topological_grid) :
        m_input_labels(input_labels), m_topological_grid(topological_grid) {
            m_output_labels = new DenseLabeling<char>(m_topological_grid->numCells());
            m_output_labels->SetAll(0);
        }


        void OutputEdgesToFile(const char* filename) {
            FILE* fout = fopen(filename, "wb");

            TopologicalRegularGrid::DCellsIterator edges(m_topological_grid, 1);
            for (edges.begin(); edges.valid(); edges.advance()) {
                INDEX_TYPE edge = edges.value();
                if (m_output_labels->GetLabel(edge) == 1)
                    fwrite(&edge, sizeof(INDEX_TYPE), 1, fout);
            }
            fclose(fout);

        }


        DenseLabeling<char>* GetOutputLabels() { return m_output_labels; }

        DenseLabeling<char>* ComputeBoundary() {
#pragma omp parallel
            {
                int num_threads = omp_get_num_threads();
                int thread_num = omp_get_thread_num();
                std::vector<INDEX_TYPE> partition;
                ArrayIndexPartitioner::EvenChunkSplit(m_topological_grid->numCells(), num_threads, partition);

                TopologicalRegularGrid::DCellsIterator edges(m_topological_grid, 1, partition[thread_num], partition[thread_num + 1]);
                for (edges.begin(); edges.valid(); edges.advance()) {
                    TopologicalRegularGrid::FacetsIterator vertices(m_topological_grid);
                    INDEX_TYPE edge = edges.value();
                    vertices.begin(edge);
                    INDEX_TYPE vertex1 = vertices.value();
                    vertices.advance();
                    INDEX_TYPE vertex2 = vertices.value();

                    INDEX_TYPE vertex_number1 = m_topological_grid->VertexNumberFromCellID(vertex1);
                    INDEX_TYPE vertex_number2 = m_topological_grid->VertexNumberFromCellID(vertex2);

                    if ((*m_input_labels)[vertex_number1] != (*m_input_labels)[vertex_number2]) {
                        (*m_output_labels)[edge] = 1;
                    }

                }
#pragma omp barrier
                TopologicalRegularGrid::DCellsIterator quads(m_topological_grid, 2, partition[thread_num], partition[thread_num + 1]);
                for (quads.begin(); quads.valid(); quads.advance()) {
                    INDEX_TYPE quad = quads.value();
                    TopologicalRegularGrid::FacetsIterator quadedges(m_topological_grid);

                    for (quadedges.begin(quad); quadedges.valid(); quadedges.advance()) {
                        if ((*m_output_labels)[quadedges.value()] == 1) {
                            (*m_output_labels)[quad] = 1;
                            break;
                        }
                    }


                }
#pragma omp barrier
                TopologicalRegularGrid::DCellsIterator voxels(m_topological_grid, 3, partition[thread_num], partition[thread_num + 1]);
                for (voxels.begin(); voxels.valid(); voxels.advance()) {
                    INDEX_TYPE voxel = voxels.value();
                    TopologicalRegularGrid::FacetsIterator voxelquads(m_topological_grid);
                    for (voxelquads.begin(voxel); voxelquads.valid(); voxelquads.advance()) {
                        if ((*m_output_labels)[voxelquads.value()] == 1) {
                            (*m_output_labels)[voxel] = 1;
                            break;
                        }
                    }

                }


            }
            return m_output_labels;
        }



        DenseLabeling<char>* ComputeBoundaryHACK() {
#pragma omp parallel
            {
                int num_threads = omp_get_num_threads();
                int thread_num = omp_get_thread_num();
                std::vector<INDEX_TYPE> partition;
                ArrayIndexPartitioner::EvenChunkSplit(m_topological_grid->numCells(), num_threads, partition);

                TopologicalRegularGrid::DCellsIterator edges(m_topological_grid, 1, partition[thread_num], partition[thread_num + 1]);
                for (edges.begin(); edges.valid(); edges.advance()) {
                    TopologicalRegularGrid::FacetsIterator vertices(m_topological_grid);
                    INDEX_TYPE edge = edges.value();
                    vertices.begin(edge);
                    INDEX_TYPE vertex1 = vertices.value();
                    vertices.advance();
                    INDEX_TYPE vertex2 = vertices.value();

                    INDEX_TYPE vertex_number1 = m_topological_grid->VertexNumberFromCellID(vertex1);
                    INDEX_TYPE vertex_number2 = m_topological_grid->VertexNumberFromCellID(vertex2);

                    if ((*m_input_labels)[vertex_number1] != (*m_input_labels)[vertex_number2]) {
                        INDEX_TYPE lv = ((*m_input_labels)[vertex_number1] > (*m_input_labels)[vertex_number2] ? (*m_input_labels)[vertex_number1] : (*m_input_labels)[vertex_number2]);
                        (*m_output_labels)[edge] = (lv % 126 + 1);
                    }

                }
#pragma omp barrier
                TopologicalRegularGrid::DCellsIterator quads(m_topological_grid, 2, partition[thread_num], partition[thread_num + 1]);
                for (quads.begin(); quads.valid(); quads.advance()) {
                    INDEX_TYPE quad = quads.value();
                    TopologicalRegularGrid::FacetsIterator quadedges(m_topological_grid);

                    for (quadedges.begin(quad); quadedges.valid(); quadedges.advance()) {
                        if ((*m_output_labels)[quadedges.value()] == 1) {
                            (*m_output_labels)[quad] = 1;
                            break;
                        }
                    }


                }
#pragma omp barrier
                TopologicalRegularGrid::DCellsIterator voxels(m_topological_grid, 3, partition[thread_num], partition[thread_num + 1]);
                for (voxels.begin(); voxels.valid(); voxels.advance()) {
                    INDEX_TYPE voxel = voxels.value();
                    TopologicalRegularGrid::FacetsIterator voxelquads(m_topological_grid);
                    for (voxelquads.begin(voxel); voxelquads.valid(); voxelquads.advance()) {
                        if ((*m_output_labels)[voxelquads.value()] == 1) {
                            (*m_output_labels)[voxel] = 1;
                            break;
                        }
                    }

                }


            }
            return m_output_labels;
        }

    };




}

#endif
