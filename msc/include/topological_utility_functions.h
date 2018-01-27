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

#ifndef TOPOLOGICAL_UTILITY_FUNCTIONS_H
#define TOPOLOGICAL_UTILITY_FUNCTIONS_H

#include "basic_types.h"
#include "vectors.h"
#include "regular_grid.h"
#include "regular_grid_trilinear_function.h"
#include "topological_regular_grid.h"


namespace MSC {

    void RunMeshConsistencyChecks(TopologicalRegularGrid* mesh) {

        TopologicalRegularGrid::AllCellsIterator allcells(mesh);

        for (allcells.begin(); allcells.valid(); allcells.advance()) {

            INDEX_TYPE id = allcells.value();
            Vec3l coords;
            mesh->cellid2Coords(id, coords);
            DIM_TYPE dim = mesh->dimension(coords);
            // check facets/cofacets
            {
                TopologicalRegularGrid::FacetsIterator fit(mesh);
                for (fit.begin(id); fit.valid(); fit.advance()) {
                    INDEX_TYPE facetid = fit.value();

                    bool has_reflex = false;
                    TopologicalRegularGrid::CofacetsIterator cfit(mesh);
                    for (cfit.begin(facetid); cfit.valid(); cfit.advance()) {
                        INDEX_TYPE cofacetid = cfit.value();
                        if (cofacetid == id) has_reflex = true;
                    }
                    if (!has_reflex) {
                        printf("Error 0: facet of cell (%llu, %llu, %llu) has no cofacet with same id\n", coords[0], coords[1], coords[2]);
                    }

                    int dir = mesh->Compress6NeighborOffsetToByte(id, facetid);
                    INDEX_TYPE test = mesh->UncompressByteTo6NeighborOffset(id, dir);
                    if (test != facetid) {
                        printf("Error 1:\n");
                    }
                    int dir2 = mesh->Compress6NeighborOffsetToByte(facetid, id);
                    INDEX_TYPE test2 = mesh->UncompressByteTo6NeighborOffset(facetid, dir2);
                    if (test2 != id) {
                        printf("Error 2:\n");
                    }

                }
            }
            {
                TopologicalRegularGrid::CofacetsIterator fit(mesh);
                for (fit.begin(id); fit.valid(); fit.advance()) {
                    INDEX_TYPE facetid = fit.value();

                    bool has_reflex = false;
                    TopologicalRegularGrid::FacetsIterator cfit(mesh);
                    for (cfit.begin(facetid); cfit.valid(); cfit.advance()) {
                        INDEX_TYPE cofacetid = cfit.value();
                        if (cofacetid == id) has_reflex = true;
                    }
                    if (!has_reflex) {
                        printf("Error 3: cofacet of cell (%llu, %llu, %llu) has no facet with same id\n", coords[0], coords[1], coords[2]);
                    }
                    int dir = mesh->Compress6NeighborOffsetToByte(id, facetid);
                    INDEX_TYPE test = mesh->UncompressByteTo6NeighborOffset(id, dir);
                    if (test != facetid) {
                        printf("Error 4:\n");
                    }
                    int dir2 = mesh->Compress6NeighborOffsetToByte(facetid, id);
                    INDEX_TYPE test2 = mesh->UncompressByteTo6NeighborOffset(facetid, dir2);
                    if (test2 != id) {
                        printf("Error 5:\n");
                    }
                }
            }

            // ---------------------------------------------------------------
            // check adjacency
            {
                // Harsh changed this on 01.18.17
                int boundary = mesh->boundaryValue(coords);
                //int boundary = mesh->MemoryBlockBoundaryValue(coords);

                //printf(" id = %d, coords = (%d %d %d), boundary = %d\n", id, coords[0], coords[1], coords[2], boundary);

                int expected = 0;
                int actual = 0;

                if (boundary == 0) {            expected = 27;      }
                else if (boundary == 1) {       expected = 27 - 9;  }
                else if (boundary == 2) {       expected = 12;      }
                else if (boundary == 3) {       expected = 8;       }
                else {
                    printf("Error 6: Invalid value of boundary for cell %d\n", id);
                }


                TopologicalRegularGrid::AdjacentCellsIterator ajit(mesh);
                for (ajit.begin(coords); ajit.valid(); ajit.advance()) {
                    INDEX_TYPE ajid = ajit.value();
                    actual++;

                    // check reflexive relation
                    bool has_reflex = false;
                    TopologicalRegularGrid::AdjacentCellsIterator ajit2(mesh);
                    for (ajit2.begin(ajid); ajit2.valid(); ajit2.advance()) {
                        if (ajit2.value() == id){ has_reflex = true;    break;  }
                    }
                    if (!has_reflex) {
                        printf("Error 7: Reflexive relation failed for cell %d\n", id);
                        exit(1);
                    }

                    // check compression and decompression
                    BYTE_TYPE dir = mesh->CompressAdjacentOffsetToByte(id, ajid, ajit);
                    INDEX_TYPE test = mesh->UncompressByteToAdjacentOffset(id, dir);
                    if (test != ajid) {
                        printf("Error 8: Compresssion and decompression failed for adjCell %d of cell %d!\n", ajid, id);
                        exit(1);
                    }
                }

                if (expected != actual) {
                    printf("Error 9: Expected %d adjacent cells for cell %d, but found %d\n", expected, id, actual);
                    exit(1);
                }
            }
            // ---------------------------------------------------------------

            // check cell verts
            {

                int countverts = 0;
                TopologicalRegularGrid::CellVerticesIterator cvit(mesh);

                for (cvit.begin(coords); cvit.valid(); cvit.advance()) {

                    INDEX_TYPE cvid = cvit.value();
                    countverts++;

                    // ----------------------------------
                    // added by Harsh..
                    if(cvid >= mesh->numCells()){
                        printf("Error 10: Found cell %d as a vertex of cell %d, but max cells = %d\n",
                               cvid, id, mesh->numCells());
                        exit(1);
                    }
                    // ----------------------------------

                    // check compression and decompression
                    BYTE_TYPE dir = mesh->CompressVertexOffsetToByte(id, cvid, cvit);
                    INDEX_TYPE test = mesh->UncompressByteToVertexOffset(id, dir);
                    if (test != cvid) {
                        printf("Error 11: Compresssion and decompression failed for vert %d of cell %d!\n", cvid, id);
                        exit(1);
                    }


                }
                if (countverts != 1 << dim) {
                    printf("Error 12: Number of verts for cell %d does not match the dimension %d\n");
                    exit(1);
                }
            }
            // ---------------------------------------------------------------
            //exit(1);
        }

        printf("\n\n ==== RunMeshConsistencyChecks()... Done!\n");
    }


    //INT_TYPE Inspect_Higher_Certains(INT_TYPE tid) {
    //	INT_TYPE tneg;
    //	int higher_certain = certains[tid];
    //	bool has_higher = false;

    //	Vec3l neighbors[6];
    //	int nn = GlobalGrid->gather_existing_neighbors6(GlobalGrid->xyz3d(tid), neighbors);
    //	for (int i = 0; i < nn; i++) {
    //		INT_TYPE tneg = GlobalGrid->index3d(neighbors[i]);
    //		if (GlobalGrid->is_greater(tneg, tid)) {
    //			if (certains[tneg] < 0) return -1; // if a higher one is uncertain, we are uncertain
    //			if (!has_higher) {
    //				higher_certain = certains[tneg];
    //				has_higher = true;
    //			}
    //			else {
    //				if (higher_certain != certains[tneg]) return -1;
    //			}
    //		}
    //	}

    //	if (!has_higher) {
    //		printf("ERROR should never get here\n");
    //		return -1;
    //	}
    //	return higher_certain;

    //}


    //void Enqueue_Lower_Neighbors(Vec3l xyz, std::priority_queue<std::pair<FLOATTYPE, INT_TYPE> > &expansion, std::set<INT_TYPE>&seen) {
    //	INT_TYPE tid = GlobalGrid->index3d(xyz);

    //	Vec3l neighbors[6];
    //	int nn = GlobalGrid->gather_existing_neighbors6(xyz, neighbors);
    //	for (int i = 0; i < nn; i++) {
    //		INT_TYPE tneg = GlobalGrid->index3d(neighbors[i]);
    //		if (certains[tneg] == -1 && GlobalGrid->is_greater(tid, tneg) && seen.count(tneg) == 0) { seen.insert(tneg); expansion.push(std::pair<FLOATTYPE, INT_TYPE>(GlobalGrid->sampleI(tneg), tneg)); }

    //	}
    //}

    //// takes a maximum, and starts expanding downwards
    //// computes set of cells from which all monotonically increasing paths terminate at the
    //// extremum
    //void Expand_Maximum_Neighborhood(INT_TYPE startid, const RegularGridTrilinearFunction* func) {
    //	Vec3l xyz = func->GetGrid()->XYZ3d(startid);
    //	std::set<INT_TYPE> seen;

    //	INT_TYPE tid = startid;
    //	// the natural ordering using the < operator on pairs will give us the highest
    //	// element first, simulating region growing from high to low
    //	std::priority_queue<std::pair<FLOATTYPE, INT_TYPE> > growing_front; // this is a hack! - we aught to be using the same simulation of simplicity as elsewhere
    //	seen.insert(startid);
    //	Enqueue_Lower_Neighbors(xyz, growing_front, seen);

    //	while (!growing_front.empty()) {

    //		INT_TYPE currid = growing_front.top().second;
    //		growing_front.pop();

    //		int cellvale = Inspect_Higher_Certains(currid);
    //		// find highers

    //		// cellvalue >=0 indicates that there is certainty here, so lets expand
    //		if (cellvale >= 0) {
    //			certains[currid] = cellvale;
    //			dests[currid] = startid;
    //			Enqueue_Lower_Neighbors(func->GetGrid()->XYZ3d(currid), growing_front, seen);
    //		}

    //	}
    //}


}


#endif
