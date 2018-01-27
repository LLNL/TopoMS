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

#ifndef ISOLATED_REGION_REMOVER_MASKED_H
#define ISOLATED_REGION_REMOVER_MASKED_H


#include "basic_types.h"
#include "vectors.h"
#include "labeling.h"
#include "regular_grid.h"
#include "regular_grid_trilinear_function.h"
#include "adaptive_euler_advector.h"
#include "timing.h"

#include <map>

namespace MSC {
    template< class Comparer>
    class IsolatedRegionRemoverMasked {

    protected:
        const DenseLabeling<int>* m_input;
        std::map<INDEX_TYPE, int> m_inversemask;

        DenseLabeling<int>* m_destinations_unmasked;

        DenseLabeling<INDEX_TYPE>* m_destinations;
        const RegularGrid* m_grid;
        RegularGridTrilinearFunction* m_func;


        // if the vertex is a max, return own index
        // if the vertex has a higher neighbor with same label, return that (or highest such)
        // if the vertex does not have higher with same label, return highest id
        INDEX_TYPE GetHighestUpVertexIn6Neighborhood(INDEX_TYPE id) const {
            Vec3l t_neighbors[6];
            Vec3l t_coords = m_grid->XYZ3d(id);
            int t_num_neighbors = m_grid->GatherExistingNeighbors6(t_coords, t_neighbors);

            bool has_samelabel_higher = false;
            bool is_max = true;
            INDEX_TYPE t_overall_highest = id;
            INDEX_TYPE t_highest_same_label = id;

            for (int i = 0; i < t_num_neighbors; i++) {
                INDEX_TYPE t_neighbor_vertex = m_grid->Index3d(t_neighbors[i]);
                if (mCompare->Compare(t_neighbor_vertex, t_overall_highest)) {
                    t_overall_highest = t_neighbor_vertex;
                    is_max = false;
                }
                if (m_input->GetLabel(id) == m_input->GetLabel(t_neighbor_vertex) &&
                    mCompare->Compare(t_neighbor_vertex, t_highest_same_label)) {
                    t_highest_same_label = t_neighbor_vertex;
                    has_samelabel_higher = true;
                }
            }
            if (is_max) return id;
            if (has_samelabel_higher) return t_highest_same_label;
            return t_overall_highest;
        }


        INDEX_TYPE PathCompressFind(INDEX_TYPE id) {
            if (m_destinations->GetLabel(id) == id) return id;
            INDEX_TYPE retval = PathCompressFind(m_destinations->GetLabel(id));
            INDEX_TYPE idref = m_destinations->operator[](id);
#pragma omp atomic
            m_destinations->operator[](id) += retval - idref; // this is stupid - openmp 2.0 supports only binops= for atomics
            return retval;
        }

        Comparer* mCompare;

    public:

        IsolatedRegionRemoverMasked(RegularGridTrilinearFunction* func, const DenseLabeling<int> *input, const std::vector<INDEX_TYPE> *mask) :
            m_input(input), m_func(func) {
            m_grid = func->GetGrid();
            mCompare = new Comparer(func);

            for(int i = 0; i < mask->size(); i++) {
                m_inversemask.insert( std::pair<INDEX_TYPE, int> (mask->at(i), i) );
            }
        }

        ~IsolatedRegionRemoverMasked(){
            delete mCompare;
            delete m_destinations;
        }

        DenseLabeling<INDEX_TYPE>* GetOutputLabels() { return m_destinations; }
        DenseLabeling<int>* GetOutputLabelsUnmasked() { return m_destinations_unmasked; }
        const RegularGrid* GetGrid() { return m_grid; }
        RegularGridTrilinearFunction* GetFunction() { return m_func; }
        void ComputeOutput(bool verbose = false) {


            ThreadedTimer gtimer(1);
            gtimer.StartGlobal();

            if(verbose){
                printf(" -- Cleaning noise in integration...");
                fflush(stdout);
            }

            const INDEX_TYPE t_num_vertices = m_grid->NumElements();

            m_destinations = new DenseLabeling<INDEX_TYPE>(t_num_vertices);


            // set all potential extrema, so we terminate near them
#pragma omp parallel for
            for (INDEX_TYPE i = 0; i < t_num_vertices; i++) {
                INDEX_TYPE nid = GetHighestUpVertexIn6Neighborhood(i);
                m_destinations->SetLabel(i, nid);
            }

//#pragma omp parallel for
            for (INDEX_TYPE i = 0; i < t_num_vertices; i++) {
                PathCompressFind(i);
            }


            ThreadedTimer ltimer(1);
            ltimer.StartGlobal();

            // create unmasked versions!
            m_destinations_unmasked = new DenseLabeling<int>(t_num_vertices);

#pragma omp parallel for
            for(int i = 0; i < t_num_vertices; i++) {

                if (m_input->GetLabel(i) < 0){
                    m_destinations->SetLabel(i, m_input->GetLabel(i));
                    continue;
                }
                m_destinations_unmasked->SetLabel(i, m_inversemask.at( m_destinations->GetLabel(i) ));
            }

            ltimer.EndGlobal();
            gtimer.EndGlobal ();
            if (verbose){
                printf(" done!");   gtimer.PrintAll();
                //ltimer.PrintAll();
            }
        }
    };
#if 0
    class PathCompressor {

    protected:
        DenseLabeling<INDEX_TYPE>* m_input;
        DenseLabeling<INDEX_TYPE>* m_destinations;
        const RegularGrid* m_grid;




        // if the vertex is a max, return own index
        // if the vertex has a higher neighbor with same label, return that (or highest such)
        // if the vertex does not have higher with same label, return highest id

        INDEX_TYPE PathCompressFind(INDEX_TYPE id) {
            if (m_destinations->GetLabel(id) == id) return id;
            INDEX_TYPE retval = PathCompressFind(m_destinations->GetLabel(id));
//			INDEX_TYPE idref = m_destinations->operator[](id);
//#pragma omp atomic
//			m_destinations->operator[](id) += retval - idref; // this is stupid - openmp 2.0 supports only binops= for atomics
            m_destinations->SetLabel(id, retval);
            return retval;
        }

        //Comparer* mCompare;

    public:
        PathCompressor(RegularGridTrilinearFunction* func, DenseLabeling<INDEX_TYPE>* input) :
            m_input(input) {
            m_grid = func->GetGrid();
        }


        DenseLabeling<INDEX_TYPE>* GetOutputLabels() { return m_destinations; }
        const RegularGrid* GetGrid() { return m_grid; }
        void ComputeOutput() {


            m_destinations = new DenseLabeling<INDEX_TYPE>(m_grid->NumElements());
            const INDEX_TYPE t_num_vertices = m_grid->NumElements();
            for (INDEX_TYPE i = 0; i < t_num_vertices; i++) {
                m_destinations->SetLabel(i, m_input->GetLabel(i));
            }

            for (INDEX_TYPE i = 0; i < t_num_vertices; i++) {
                PathCompressFind(i);
            }

        }



    };

#endif

}
#endif
