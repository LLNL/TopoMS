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

#ifndef TOPOLOGICAL_EXPLICIT_DENSE_MESH_FUNCTION
#define TOPOLOGICAL_EXPLICIT_DENSE_MESH_FUNCTION

#include <algorithm>
#include <omp.h>
#include "basic_types.h"
#include "topological_regular_grid.h"
#include "array_index_partition.h"
#include "regular_grid_trilinear_function.h"

namespace MSC {

    typedef TopologicalRegularGrid Mesh;
    typedef float dtype;

    class TopologicalExplicitDenseMeshFunction {
    protected:
        Mesh* mMesh;
        dtype* mVals;
    public:

        //// assume all facet values are correct!!!
        dtype maxOfFacets(INDEX_TYPE cellid) {
            dtype result;
            Mesh::FacetsIterator fit(mMesh);
            fit.begin(cellid);
            result = this->cellValue(fit.value());
            fit.advance();
            for (; fit.valid(); fit.advance()) {
                dtype tmpother = cellValue(fit.value());
                result = (result > tmpother ? result : tmpother);
            }
            return result;
        }

        INDEX_TYPE maxIdOfVerts(INDEX_TYPE cellid) {
            INDEX_TYPE result;
            Mesh::CellVerticesIterator vit(mMesh);
            vit.begin(cellid);
            result = vit.value();
            vit.advance();
            for (; vit.valid(); vit.advance()) {
                INDEX_TYPE tmpotherid = vit.value();
                if (this->greaterThan(tmpotherid, result)) result = tmpotherid;
            }
            return result;
        }



        dtype minOfFacets(INDEX_TYPE cellid) {
            dtype result;
            Mesh::FacetsIterator fit(mMesh);
            fit.begin(cellid);
            result = this->cellValue(fit.value());
            fit.advance();
            for (; fit.valid(); fit.advance()) {
                dtype tmpother = cellValue(fit.value());
                result = (result < tmpother ? result : tmpother);
            }
            return result;
        }
        //   // assume all facet values are correct!!!
        //dtype ave_of_facets(CELL_INDEX_TYPE cellid) {
        //
        // dtype result = 0;
        // dtype count = 0;
        //
        // cellIterator it;
        // iteratorOperator& ito = my_mesh_handler->facets(cellid, it);
        // ito.begin(it);
        //
        // while (ito.valid(it)) {
        //  result += cell_value(ito.value(it));
        //  count += 1;
        // 	   ito.advance(it);
        // }
        // return result / count;
        //}
        //// first do 0's. Assume that the i'th 0-cell in the mesh is the i'th vertex in data
        //void set_vertex_values() {
        // mscBasicArray<dtype>& values_r = *(my_values);

        // cellIterator it;
        // iteratorOperator& ito = my_mesh_handler->d_cells_iterator(0, it);
        // ito.begin(it);
        // CELL_INDEX_TYPE temp_grad_2_data = 0;
        //
        // while (ito.valid(it)) {
        //  values_r[ito.value(it)] = my_data_handler->value(temp_grad_2_data);
        //  temp_grad_2_data++;
        //  ito.advance(it);
        // }
        //}

        //void set_cell_values(DIM_TYPE dim) {
        // mscBasicArray<dtype>& values_r = *(my_values);

        // cellIterator it;
        // iteratorOperator& ito = my_mesh_handler->d_cells_iterator(dim, it);
        // ito.begin(it);
        //
        // while (ito.valid(it)) {
        //  CELL_INDEX_TYPE temp_id = ito.value(it);
        //  values_r[temp_id] = max_of_facets(temp_id);
        //  ito.advance(it);
        // }
        //}

    public:
        TopologicalExplicitDenseMeshFunction()  {
            mVals = NULL;
            mMesh = NULL;
        }

        ~TopologicalExplicitDenseMeshFunction()  {
            if (mVals != NULL) delete[] mVals;
        }

        void setMeshAndAllocate(Mesh* m) {
            mMesh = m;
            if (mVals != NULL) delete[] mVals;
            mVals = new dtype[m->numCells()];
        }

        dtype cellValue(INDEX_TYPE cellid) const {
            return mVals[cellid];
        }

        void setCellValue(INDEX_TYPE cellid, dtype val) {
            mVals[cellid] = val;
        }

        void setCellValuesMaxOfVerts() {
            for (int i = 1; i <= mMesh->maxDim(); i++) {
#pragma omp parallel
            {
                int num_threads = omp_get_num_threads();
                int thread_num = omp_get_thread_num();
                std::vector<INDEX_TYPE> partition;
                ArrayIndexPartitioner::EvenChunkSplit(mMesh->numCells(), num_threads, partition);
                Mesh::DCellsIterator dcells(mMesh, i, partition[thread_num], partition[thread_num + 1]);
                for (dcells.begin(); dcells.valid(); dcells.advance()) {
                    INDEX_TYPE cellid = dcells.value();
                    setCellValue(cellid, maxOfFacets(cellid));
                }
            }
            }
        }

        void setCellValuesMinOfVerts() {
            for (int i = 1; i <= mMesh->maxDim(); i++) {
#pragma omp parallel
            {
                int num_threads = omp_get_num_threads();
                int thread_num = omp_get_thread_num();
                std::vector<INDEX_TYPE> partition;
                ArrayIndexPartitioner::EvenChunkSplit(mMesh->numCells(), num_threads, partition);
                Mesh::DCellsIterator dcells(mMesh, i, partition[thread_num], partition[thread_num + 1]);
                for (dcells.begin(); dcells.valid(); dcells.advance()) {
                    INDEX_TYPE cellid = dcells.value();
                    setCellValue(cellid, minOfFacets(cellid));
                }
            }
            }
        }

        const Mesh* mesh() const {
            return this->mMesh;
        }
        bool lessThan(INDEX_TYPE a, INDEX_TYPE b) const {
            dtype av = cellValue(a);
            dtype bv = cellValue(b);
            if (av < bv) return true;
            if (bv < av) return false;
            return a < b;
        }
        bool greaterThan(INDEX_TYPE a, INDEX_TYPE b) const {
            return lessThan (b, a);
        }

        bool copyVertexValuesFromGridFunction(RegularGridTrilinearFunction* func) {
#pragma omp parallel
            {
                int num_threads = omp_get_num_threads();
                int thread_num = omp_get_thread_num();
                std::vector<INDEX_TYPE> partition;
                ArrayIndexPartitioner::EvenChunkSplit(mMesh->numCells(), num_threads, partition);
                Mesh::DCellsIterator vertices(mMesh, 0, partition[thread_num], partition[thread_num + 1]);
                for (vertices.begin(); vertices.valid(); vertices.advance()) {
                    INDEX_TYPE cellid = vertices.value();
                    setCellValue(cellid,func->SampleImage(this->mMesh->VertexNumberFromCellID(cellid)) );
                }

            }
            return true;
        }


        // assumes file is binary, with contiguous dense vertex values, where the ith value corresponds to the
        // value of the vertex after the vertex iterator has been advanced i times.
        // must come after the set mesh and allocate call, because it uses the mesh iterator
        // to go through vertices
        // return value is bool representing success or failure
        bool loadVerticesFromDenseFile(const char* filename) {
            // check that mesh has been set, return if not
            if (mMesh == NULL) {
                printf("Error: load called before mesh was set and allocated\n");
                return false;
            }

            FILE* fdat = fopen(filename, "rb");
            if (fdat == NULL) {
                printf("ERROR: file %s not found\n", filename);
                return false;
            }

#define TBUFFSIZE 4096
            dtype tbuff[TBUFFSIZE];

            Mesh::DCellsIterator zeros(mMesh, 0);
            INDEX_TYPE dataposition = 0;
            INDEX_TYPE maxposition = mMesh->numCells(0);

            for (zeros.begin(); zeros.valid(); zeros.advance()) {

                if (dataposition % TBUFFSIZE == 0) {
                    int numtoread = std::min((INDEX_TYPE) TBUFFSIZE, maxposition - dataposition);
                    int numbytesread = fread(tbuff, sizeof(dtype), numtoread, fdat);
                    if (numtoread != numbytesread) {
                        printf("Error: read misalignment, %d, %d, %d\n", numtoread, numbytesread, dataposition);
                        return false;
                    }
                }

                INDEX_TYPE vert = zeros.value();
                this->setCellValue(vert, tbuff[dataposition % TBUFFSIZE]);
                dataposition++;
            }
            fclose(fdat);
            return true;
        }

    };
}
#endif
