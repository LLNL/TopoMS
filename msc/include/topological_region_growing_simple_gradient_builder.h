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

#ifndef TOPOLOCICAL_REGION_GROWING_SIMPLE_GRADIENT_BUILDER_H
#define TOPOLOCICAL_REGION_GROWING_SIMPLE_GRADIENT_BUILDER_H

#include "basic_types.h"
#include "topological_regular_grid.h"
#include "topological_explicit_mesh_function.h"
#include "discrete_gradient_labeling.h"

#include <vector>
#include <queue>

namespace MSC {

    typedef TopologicalRegularGrid MeshType;
    typedef TopologicalExplicitDenseMeshFunction FuncType;
    typedef DiscreteGradientLabeling GradType;

    //template <typename dtype>
    class TopologicalRegionGrowingSimpleGradientBuilder {
    protected:

        MeshType* mMesh;
        FuncType* mFunc;
        GradType* mGrad;
        //mscBasicMeshFunction<dtype>* mFunc;
        //   mscBasicMeshHandler* mMesh;
        //   mscBasicGradientField* mGrad;
        //mscArrayFactory* my_array_factory;

        // set all cells to unassigned
        virtual void initAssigned() {

            MeshType::AllCellsIterator allit(mMesh);
            for (allit.begin(); allit.valid(); allit.advance()) {
                INDEX_TYPE cellid = allit.value();
                mGrad->setAssigned(cellid, 0);
                mGrad->setMark(cellid, 0);
            }
        }
        // set number of unpaired facets for all cells
        virtual void initNumberUnpairedFacets() {
            MeshType::AllCellsIterator allit(mMesh);
            for (allit.begin(); allit.valid(); allit.advance()) {
                INDEX_TYPE cellid = allit.value();
                INT_TYPE temp_num = 0;
                MeshType::FacetsIterator facets(mMesh);
                for (facets.begin(cellid); facets.valid(); facets.advance()) {
                    temp_num++;
                }
                //DIM_TYPE dd = mMesh->dimension(cellid);
                //if (temp_num != dd * 2)
                //	printf("whoathere\n");

                mGrad->setNumUnpairedFacets(cellid, temp_num);
            }
        }

        INDEX_TYPE doublecount;

        virtual void initAll() {
            initAssigned();
            initNumberUnpairedFacets();
        }

        // comparison object for queue
        struct comparison_element {
            INDEX_TYPE cellid;
            INDEX_TYPE insertion_time;
            dtype value;
            BOUNDARY_TYPE boundary;
        };

        // comparator of objects
        struct element_comparator {
            // return true if b comes first, then a
            bool operator() (const comparison_element& a, const comparison_element& b) {

                // then values
                if (a.value < b.value) return false;
                if (a.value > b.value) return true;

                // first boundaries
                if (a.boundary < b.boundary) return true;
                if (a.boundary > b.boundary) return false;

                // then insertion time, later first
                if (a.insertion_time > b.insertion_time) return true;
                if (a.insertion_time < b.insertion_time) return false;

                return a.cellid > b.cellid;
            }
        };

        struct oneleft_element {
            INDEX_TYPE cellid;
            BOUNDARY_TYPE boundary;
            DIM_TYPE dim;
            dtype value;
            INT_TYPE weight;
        };

        struct oneleft_comparator {
            bool operator() (const oneleft_element& a, const oneleft_element& b) {
                // first boundaries
                if (a.boundary < b.boundary) return true;
                if (a.boundary > b.boundary) return false;

                // then dim
                if (a.dim < b.dim) return false;
                if (a.dim > b.dim) return true;

                // then values
                if (a.value < b.value) return false;
                if (a.value > b.value) return true;

                // then weight, lower first
                if (a.weight < b.weight) return false;
                if (a.weight > b.weight) return true;

                return a.cellid > b.cellid;

            }
        };

        //// simple simulation of simplicity comparator
        //virtual bool lessThan(const INDEX_TYPE& a, const INDEX_TYPE& b) {
        //	if (mFunc->cellValue(a) < mFunc->cellValue(b)) return true;
        //	if (mFunc->cellValue(a) > mFunc->cellValue(b)) return false;
        //	//printf("INDEX BREAKING TIES\n");
        //	return a < b;
        //}

        //virtual bool fGreaterThan(const INDEX_TYPE& a, const INDEX_TYPE& b) {
        //	return mFunc->cellValue(a) > mFunc->cellValue(b);
        //}

        virtual bool lessThanAllUnassignedCofacets(const INDEX_TYPE& cellid) {

            BOUNDARY_TYPE boundary = mMesh->boundaryValue(cellid);
            MeshType::CofacetsIterator cofacets(mMesh);
            MeshType::FacetsIterator facets(mMesh);
            cofacets.begin(cellid);
            while (cofacets.valid()) {
                INDEX_TYPE cofacetid = cofacets.value();
                if (mGrad->getAssigned(cofacetid) == 0 &&
                    boundary == mMesh->boundaryValue(cofacetid)) {
                    //if (fGreaterThan(cofacetid, cellid)) {
                    //	cofacets.advance(it);
                    //	continue;
                    //} else {

                    facets.begin(cofacetid);
                    while (facets.valid()) {
                        INDEX_TYPE facetid = facets.value();
                        if (facetid != cellid &&
                            mGrad->getAssigned(facetid) == 0 &&
                            !this->mFunc->lessThan(cellid, facetid))
                            return false;
                        facets.advance();
                    }
                    //}
                }
                cofacets.advance();
            }
            return true;
        }

        std::priority_queue<comparison_element, std::vector< comparison_element >, element_comparator> mSortedCellQueue;
        std::priority_queue<oneleft_element, std::vector< oneleft_element >, oneleft_comparator> mOneleftCellQueue;

        virtual void enqueueOneleftElement(const INDEX_TYPE& cellid) {
            oneleft_element to_insert;
            // basic optimization: 1-cells will never have to be "zipped up"
            to_insert.dim = mMesh->dimension(cellid);
            if (to_insert.dim == 1) return;

            // compute weight
            int weight = 0;
            int SANITY_COUNT = 0;
            MeshType::FacetsIterator facets(mMesh);
            for (facets.begin(cellid); facets.valid(); facets.advance()) {
                INDEX_TYPE temp_id = facets.value();
                if (mGrad->getAssigned(temp_id)) {
                    if (mGrad->getCritical(temp_id)) {
                        weight += mMesh->dimension(temp_id);
                    }
                    else {
                        weight +=
                            mMesh->dimension(mGrad->getPair(temp_id));
                    }
                }
                else {
                    SANITY_COUNT++;
                }
            }

            //if (SANITY_COUNT != 1) printf("INSANE IN THE MAINFRAME\n");

            to_insert.cellid = cellid;
            to_insert.boundary = mMesh->boundaryValue(cellid);
            to_insert.dim = mMesh->dimension(cellid);
            to_insert.value = mFunc->cellValue(cellid);
            to_insert.weight = weight;
            mOneleftCellQueue.push(to_insert);

        }

        INDEX_TYPE insert_time;
        virtual void enqueueSortedElement(const INDEX_TYPE& cellid, INDEX_TYPE itime, unsigned char mark) {
            comparison_element to_insert;
            to_insert.cellid = cellid;
            to_insert.insertion_time = itime++;
            to_insert.value = mFunc->cellValue(cellid);
            to_insert.boundary = mMesh->boundaryValue(cellid);
            //printf("insert %d -> fval=%.2f\n", cellid, to_insert.value);
            mGrad->setMark(cellid, mark);
            mSortedCellQueue.push(to_insert);
        }

        // seed my queue
        virtual void seedQueueWithDMinima(const DIM_TYPE& dim) {
            int counter = 0;
            MeshType::DCellsIterator d_cells(mMesh, dim);
            for (d_cells.begin(); d_cells.valid(); d_cells.advance()) {
                INDEX_TYPE cellid = d_cells.value();

                if (!mGrad->getAssigned(cellid) &&
                    lessThanAllUnassignedCofacets(cellid)) {
                    // potential critical point, so enqueue
                    //printf("potential minimum %d\n", cellid);
                    enqueueSortedElement(cellid, doublecount, 0);
                    counter++;
                }
            }
            printf("seeding queue with %d %d-cells\n", counter, dim);
        }

        // we don't need to add d+1 cells, only d+2 cells
        virtual void decrementCofacetsNumUnpairedFacets(const INDEX_TYPE& cellid,
            const bool& add_to_oneleft) {
            MeshType::CofacetsIterator cofacets(mMesh);
            for (cofacets.begin(cellid); cofacets.valid(); cofacets.advance()) {
                INDEX_TYPE temp_cell = cofacets.value();
                if (mGrad->getAssigned(temp_cell) == false) {
                    INDEX_TYPE num_unpaired =
                        mGrad->getNumUnpairedFacets(temp_cell);
                    //if (num_unpaired == 0) continue;
                    num_unpaired -= 1;
                    //if (num_unpaired < 0 || num_unpaired > 7) printf("UNP=%lld\n", num_unpaired);
                    mGrad->setNumUnpairedFacets(temp_cell, num_unpaired);
                    if (num_unpaired == 1 && add_to_oneleft) {
                        enqueueOneleftElement(temp_cell);
                    }
                }
            }
        }


        virtual void makeCritical(const INDEX_TYPE& cellid,
            const DIM_TYPE& dim) {
            //printf("making critical %d %d\n", cellid, dim);
            mGrad->setAssigned(cellid, true);
            mGrad->setCritical(cellid, true);
            mGrad->setDimAscMan(cellid, mMesh->maxDim() - dim);

            decrementCofacetsNumUnpairedFacets(cellid, false);
        }


        virtual void pair(const INDEX_TYPE& tail,
            const INDEX_TYPE& head,
            const bool& add_tail_to_oneleft) {
            //printf("p(%d->%d)\n", tail, head);
            mGrad->setAssigned(tail, true);
            mGrad->setAssigned(head, true);

            DIM_TYPE mindim = mMesh->maxDim();
            MeshType::FacetsIterator facets(mMesh);
            for (facets.begin(head); facets.valid(); facets.advance()) {
                INDEX_TYPE temp_id = facets.value();
                if (temp_id == tail) continue;
                DIM_TYPE otherdim = mGrad->getDimAscMan(temp_id);
                if (otherdim < mindim) mindim = otherdim;
            }

            mGrad->setCritical(head, false);
            mGrad->setCritical(tail, false);
            mGrad->setDimAscMan(head, mindim);
            mGrad->setDimAscMan(tail, mindim);
            mGrad->setPair(head, tail);
            mGrad->setPair(tail, head);

            decrementCofacetsNumUnpairedFacets(tail, add_tail_to_oneleft);
            decrementCofacetsNumUnpairedFacets(head, true);

        }
        ////////////////////////
        /////// ADD NEIGHBORS - facets of co-facets that are not assigned and not marked
        ////////////

        virtual void addNeighborsToSort(const INDEX_TYPE& cellid) {

            MeshType::CofacetsIterator cofacets(mMesh);
            MeshType::FacetsIterator facets(mMesh);
            for (cofacets.begin(cellid); cofacets.valid(); cofacets.advance()) {
                INDEX_TYPE temp_cell = cofacets.value();

                for (facets.begin(temp_cell); facets.valid(); facets.advance()) {
                    INDEX_TYPE temp_neg = facets.value();

                    if (temp_neg != cellid &&
                        !mGrad->getAssigned(temp_neg) &&
                        !mGrad->getMark(temp_neg)
                        ) {
                        //printf("enqueuing %d\n", temp_neg);
                        enqueueSortedElement(temp_neg, ++insert_time, 1);
                    }
                }
            }
        }



        // return index in candidates of pair
        virtual int pickFromCandidates(const INDEX_TYPE& cellid,
            const std::vector<INDEX_TYPE>& candidates) {
            int minloc = 0;
            dtype minval = lowestFacetValue(candidates[minloc]);
            //printf("%d=%.2f\n", candidates[0], minval);
            for (int i = 1; i < candidates.size(); i++) {

                dtype otherval = lowestFacetValue(candidates[i]);
                //printf("%d=%.2f\n", candidates[i], otherval);
                if (otherval < minval) {
                    minval = otherval;
                    minloc = i;
                }
            }
            return minloc;
        }

        virtual void pickAndPair(const comparison_element& element,
            const DIM_TYPE& dim) {
            // assume it's unassigned
            //if (mGrad->getAssigned(element.cellid))
            //	printf("WHOA assigned already!!\n");

            std::vector<INDEX_TYPE> candidates;
            MeshType::CofacetsIterator cofacets(mMesh);
            for (cofacets.begin(element.cellid); cofacets.valid(); cofacets.advance()) {
                INDEX_TYPE temp_cell = cofacets.value();
                // it's not assigned and is "lower"
                if (mGrad->getAssigned(temp_cell) == false &&
                    mGrad->getNumUnpairedFacets(temp_cell) == 1 &&
                    element.value >= mFunc->cellValue(temp_cell) &&
                    element.boundary == mMesh->boundaryValue(temp_cell)) {
                    candidates.push_back(temp_cell);
                }
            }

            addNeighborsToSort(element.cellid);

            // if no candidates, we have a critical point
            if (candidates.size() == 0) {

                //MeshType::CofacetsIterator cofacets2(mMesh);
                //for (cofacets2.begin(element.cellid); cofacets2.valid(); cofacets2.advance()) {
                //	INDEX_TYPE temp_cell = cofacets2.value();
                //	// it's not assigned and is "lower"
                //	DIM_TYPE temp_dim = mMesh->dimension(temp_cell);
                //	ASSIGNED_TYPE temp_assigned = mGrad->getAssigned(temp_cell);
                //	int temp_numm = mGrad->getNumUnpairedFacets(temp_cell);
                //	BOUNDARY_TYPE temp_boundary = mMesh->boundaryValue(temp_cell);
                //	float fval = mFunc->cellValue(temp_cell);
                //
                //}

                makeCritical(element.cellid, dim);
                return;
            }

            // so we have candidates
            if (candidates.size() == 1) {
                pair(element.cellid, candidates[0], false);

                return;
            }

            // else we have a choice.
            // for now pick lowest value
            //printf("have coice:\n");
            int minloc = pickFromCandidates(element.cellid, candidates);
            pair(element.cellid, candidates[minloc], false);

        }

        virtual void assignDArrows(const DIM_TYPE& dim) {
            // initialize queue
            seedQueueWithDMinima(dim);

            while (!mSortedCellQueue.empty()) {
                comparison_element top_element = mSortedCellQueue.top();
                mSortedCellQueue.pop();

                //if (! mGrad->getAssigned(top_element.cellid) &&
                //	! mGrad->getMark(top_element.cellid)) {
                //		printf("WHOA THERE NELLY\n");
                //}
                //intf("considering %d\n", top_element.cellid);
                if (!mGrad->getAssigned(top_element.cellid)) {
                    //intf("pairing->%d\n", top_element.cellid);
                    pickAndPair(top_element, dim);
                }
                //else {
                //	printf("whoa\n");
                //}
            }
        }

        virtual void findOnlyAndPair(const oneleft_element& element) {

            int SANITY_COUNT = 0;
            MeshType::FacetsIterator facets(mMesh);
            for (facets.begin(element.cellid); facets.valid(); facets.advance()) {
                INDEX_TYPE temp_id = facets.value();

                if (mGrad->getAssigned(temp_id)) continue;

                SANITY_COUNT++;
                // now we have the unassigned one
                // TEST FOR DECREASING
                if (mFunc->cellValue(temp_id) ==
                    mFunc->cellValue(element.cellid) &&
                    element.boundary == mMesh->boundaryValue(temp_id))
                    pair(temp_id, element.cellid, true);

            }
            //if (SANITY_COUNT != 1) printf("UUUUBER INSAAAANE! %d %d\n", SANITY_COUNT,
            //	mMesh->dimension(element.cellid));


        }

        virtual void zipUp() {

            while (!mOneleftCellQueue.empty()) {
                oneleft_element top_element = mOneleftCellQueue.top();
                mOneleftCellQueue.pop();

                if (mGrad->getAssigned(top_element.cellid)) continue;
                findOnlyAndPair(top_element);
            }

        }

    public:

        TopologicalRegionGrowingSimpleGradientBuilder(
            FuncType* mesh_function,
            MeshType* mesh_handler,
            GradType* grad_field) :
            mFunc(mesh_function),
            mMesh(mesh_handler),
            mGrad(grad_field) {
            insert_time = 0;
            doublecount = mesh_handler->numCells() * 2;

        }

        virtual void computeGradient() {
            initAll();
            for (DIM_TYPE i = 0; i <= mMesh->maxDim(); i++) {
                assignDArrows(i);
                zipUp();
            }
        }

        virtual dtype lowestFacetValue(const INDEX_TYPE& cellid) {

            if (mMesh->dimension(cellid) == 0) return mFunc->cellValue(cellid);

            bool init = false;
            dtype value;

            MeshType::FacetsIterator facets(mMesh);

            for (facets.begin(cellid); facets.valid(); facets.advance()) {
                INDEX_TYPE temp_neg = facets.value();
                //dtype nvalue = mFunc->cellValue(temp_neg);
                dtype nvalue = lowestFacetValue(temp_neg);
                if (!init) {
                    value = nvalue;	init = true;
                }
                else if (nvalue < value) {
                    value = nvalue;
                }
            }
            return value;
        }
    };

}


#endif
