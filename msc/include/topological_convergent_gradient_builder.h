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

#ifndef TOPOLOGICAL_CONVERGENT_GRADIENT_BUILDER_H
#define TOPOLOGICAL_CONVERGENT_GRADIENT_BUILDER_H

#include <chrono>
#include "basic_types.h"
#include "topological_region_growing_simple_gradient_builder.h"

#include <vector>
#include <queue>
#include <map>
#ifndef WIN32
#include <cmath>
#endif

namespace MSC {

    struct idfpair {
        INDEX_TYPE id;
        float prob;
        bool picked;
    };
    struct MemberDist {
        INT_TYPE count;
        INDEX_TYPE mydest;
        std::vector<idfpair> pairs;
    };

    class TopologicalConverventGradientBuilder : public TopologicalRegionGrowingSimpleGradientBuilder  {
    public:
        //struct MemberDist;
    protected:


        bool my_erase;

#ifdef FANCY_PROB_OUTPUT
        std::vector<idfpair> maxvals;
#endif

        virtual bool mdDecrement(const INDEX_TYPE& cellid, const DIM_TYPE& dim) {

            DIM_TYPE mydim = this->mMesh->maxDim() - dim;

            MeshType::FacetsIterator facets(this->mMesh);
            for (facets.begin(cellid); facets.valid(); facets.advance()) {
                INDEX_TYPE temp_cell = facets.value();
                // it's not assigned and is "lower"
                if (this->mGrad->getAssigned(temp_cell) == true &&
                    this->mGrad->getDimAscMan(temp_cell) == mydim) {
                    // now decrement
                    if (my_dists[temp_cell].count == 1 && my_erase) {
                        //printf("erasing %d\n", candidates[i]);
                        my_dists.erase(temp_cell);

                    }
                    else {
                        my_dists[temp_cell].count--;
                    }
                }
            }

            return true;
        }

        virtual bool mdCombineAllFacetsAndDecrement(const INDEX_TYPE& cellid,
            const DIM_TYPE& dim, MemberDist& result) {

            std::vector<INDEX_TYPE> candidates;
            DIM_TYPE mydim = this->mMesh->maxDim() - dim;

            MeshType::FacetsIterator facets(this->mMesh);
            for (facets.begin(cellid); facets.valid(); facets.advance()) {
                INDEX_TYPE temp_cell = facets.value();
                // it's not assigned and is "lower"
                if (this->mGrad->getAssigned(temp_cell) == true &&
                    this->mGrad->getDimAscMan(temp_cell) == mydim) {
                    candidates.push_back(temp_cell);
                }
            }

            if (candidates.size() == 0) {
                //printf("ERROR: mdCombineAllFacets candidates.size == 0, %d, %d\n", dim, mydim);
                return false;
            }

            ////if (candidates.size() == 1) {
            ////	if (my_dists.count(candidates[0]) == 0) {
            ////		printf("ERROR: mdCombineAllFacets, my_dists does not have candidate\n");
            ////	}
            ////	return my_dists[candidates[0]];
            ////}

            float oneoversize = 1.0f / (float)candidates.size();
            for (int i = 0; i < candidates.size(); i++) {
                if (my_dists.count(candidates[i]) == 0) {
                    printf("ERROR: mdCombineAllFacets, my_dists does not have candidate %d\n", i);
                }
                mdCombine(oneoversize, my_dists[candidates[i]], result);

                // now decrement
                if (my_dists[candidates[i]].count == 1 && my_erase) {
                    //printf("erasing %d\n", candidates[i]);
                    my_dists.erase(candidates[i]);

                }
                else {
                    my_dists[candidates[i]].count--;
                }

            }
            return true;
        }




        std::map<INDEX_TYPE, MemberDist> my_dists;





        // set all cells to unassigned
        virtual void initAssigned() {
            MeshType::AllCellsIterator allit(this->mMesh);
            for (allit.begin(); allit.valid(); allit.advance()) {
                INDEX_TYPE cellid = allit.value();
                this->mGrad->setAssigned(cellid, 0);
                this->mGrad->setMark(cellid, 0);
            }
        }

        virtual void initAll() {
            initAssigned();
            this->initNumberUnpairedFacets();
        }









        // seed my queue
        virtual void seedQueueWithDMinima(const DIM_TYPE& dim) {
            int counter = 0;
            MeshType::DCellsIterator d_cells(this->mMesh, dim);
            for (d_cells.begin(); d_cells.valid(); d_cells.advance()) {
                INDEX_TYPE cellid = d_cells.value();

                if (!this->mGrad->getAssigned(cellid) &&
                    this->lessThanAllUnassignedCofacets(cellid)) {
                    // potential critical point, so enqueue
                    //if (dim ==0)printf("potential minimum %d\n", cellid);
                    this->enqueueSortedElement(cellid, this->doublecount, 0);
                    counter++;
                }
            }
            printf("seeding queue with %d %d-cells\n", counter, dim);
        }

        ////// we don't need to add d+1 cells, only d+2 cells
        ////virtual void decrementCofacetsNumUnpairedFacets(const INDEX_TYPE& cellid,
        ////	const bool& add_to_oneleft, int& count) {
        ////		count = 0;
        ////	cellIterator it;
        ////	iteratorOperator& cofacets = mMesh->cofacets(cellid, it);
        ////	for (cofacets.begin(it); cofacets.valid(it); cofacets.advance(it)) {
        ////		INDEX_TYPE temp_cell = cofacets.value(it);
        ////		if (this->mGrad->getAssigned(temp_cell) == false) {
        ////			count++;
        ////			INDEX_TYPE num_unpaired =
        ////				this->mGrad->getNumUnpairedFacets(temp_cell);
        ////			//if (num_unpaired == 0) continue;
        ////			num_unpaired -= 1;
        ////			//if (num_unpaired < 0 || num_unpaired > 7) printf("UNP=%lld\n", num_unpaired);
        ////			this->mGrad->set_num_unpaired_facets(temp_cell, num_unpaired);
        ////			if (num_unpaired == 1 && add_to_oneleft) {
        ////				enqueue_oneleft_element(temp_cell);
        ////			}
        ////		}
        ////	}
        ////}

        // number of cofacets that will look up this distribution via a call to
        // mdCombineAllFacetsAndDecrement
        int countUnpairedFacets(const INDEX_TYPE& cellid) {
            int counter = 0;
            //cellIterator it;
            MeshType::FacetsIterator facets(this->mMesh);
            for (facets.begin(cellid); facets.valid(); facets.advance()) {
                INDEX_TYPE temp_cell = facets.value();
                // it's not assigned and is "lower"
                if (this->mGrad->getAssigned(temp_cell) == false) {
                    counter++;
                }
            }

            return counter;
        }

        virtual int countMDCofacets(const INDEX_TYPE& cellid) {
            int count = 0;
            MeshType::CofacetsIterator cofacets(this->mMesh);
            for (cofacets.begin(cellid); cofacets.valid(); cofacets.advance()) {
                INDEX_TYPE temp_cell = cofacets.value();
                //if (this->mGrad->getAssigned(temp_cell) == false &&
                //	count_unpaired_facets(temp_cell) > 1) {
                //		count++;
                //}
                if (this->mGrad->getAssigned(temp_cell) == false &&
                    this->mGrad->getNumUnpairedFacets(temp_cell) > 1 /*&&
                    this->mFunc->cell_value(cellid) >= this->mFunc->cell_value(temp_cell)*/) {
                    count++;
                }
            }
            return count;
        }

        virtual void makeCritical(const INDEX_TYPE& cellid,
            const DIM_TYPE& dim) {
            //printf("making critical %d %d\n", cellid, dim);
            this->mGrad->setAssigned(cellid, true);
            this->mGrad->setCritical(cellid, true);
            this->mGrad->setDimAscMan(cellid, this->mMesh->maxDim() - dim);

            int count = countMDCofacets(cellid);

            this->decrementCofacetsNumUnpairedFacets(cellid, false);

            idfpair p; p.id = cellid; p.prob = 1.0f; p.picked = true;
#ifdef FANCY_PROB_OUTPUT
            maxvals.push_back(p);
#endif
            //now add the distribution!!
            if (count > 0) {
                MemberDist md;
                md.count = count;
                md.mydest = cellid;
                md.pairs.push_back(p);
                my_dists[cellid] = md;
            }
        }


        // return index in candidates of pair
        virtual int pickFromCandidates(const INDEX_TYPE& cellid,
            const std::vector<INDEX_TYPE>& candidates,
            const std::vector<MemberDist>& mdcand) {

            // find weights for computing my probability
            dtype cell_value = this->mFunc->cellValue(cellid);

            int result = 0;
            std::vector<float> values;
            float temp_sum = 0;

            for (int i = 0; i < candidates.size(); i++) {
                float diff_val = (float)(cell_value - this->lowestFacetValue(candidates[i]));
                temp_sum += diff_val;
                values.push_back(diff_val);
            }
            if (temp_sum == 0.0f) {
                for (int i = 0; i < values.size(); i++) values[i] = 1.0f / (float)values.size();
            }
            else {
                for (int i = 0; i < values.size(); i++) values[i] = values[i] / temp_sum;
            }

            // compute my local probabilities
            MemberDist my_dist;

            for (int i = 0; i < mdcand.size(); i++) {
                this->mdCombine(values[i], mdcand[i], my_dist);
            }
            // yay, have my distribution now!
            my_dist.count = countMDCofacets(cellid);

#ifdef DEBUG_SIZE_OF_MD
            this->mdcounts[my_dist.pairs.size()]++;
#endif

#ifdef FANCY_PROB_OUTPUT
            idfpair p; p.id = cellid;
            p.prob = this->mdMax(my_dist).prob;
            maxvals.push_back(p);
#endif
#ifdef DEBUG_SIZE_OF_MD
            if (mdcounter % 200 == 0) {
                mdsizes.push_back(my_dists.size());
            }
#endif
            // compute my local weights

            int maxloc = 0;
            std::vector<float> tress; // the similarity of candidate[i]
            float maxval = mdDot3(mdcand[0], my_dist);
            tress.push_back(maxval);
            //printf("%d=%.2f\n", candidates[0], minval);
            for (int i = 1; i < mdcand.size(); i++) {

                float otherval = mdDot3(mdcand[i], my_dist);
                //printf("%d=%.2f\n", candidates[i], otherval);
                tress.push_back(otherval);
                if (otherval > maxval) {
                    maxval = otherval;
                    maxloc = i;
                }
            }

            // check if lowest
            for (int i = 0; i < mdcand.size(); i++) {
                if (i == maxloc) continue;
                if (tress[maxloc] == tress[i]) {
                    // pick lower value
                    if (this->lowestFacetValue(candidates[i]) <= this->lowestFacetValue(candidates[maxloc])){
                        //printf("EHHEHEHE\n");
                        maxloc = i;
                    }
                }
            }


            // maxloc has my pair
            for (int i = 0; i < my_dist.pairs.size(); i++) {
                my_dist.pairs[i].picked = false;
            }
            for (int i = 0; i < my_dist.pairs.size(); i++) {
                for (int j = 0; j < mdcand[maxloc].pairs.size(); j++) {
                    if (my_dist.pairs[i].id == mdcand[maxloc].pairs[j].id) {
                        my_dist.pairs[i].picked |= mdcand[maxloc].pairs[j].picked;
                    }
                }
            }

            if (my_dist.count > 0) this->my_dists[cellid] = my_dist;
            return maxloc;


        }

        virtual void pickAndPair(const TopologicalRegionGrowingSimpleGradientBuilder::comparison_element& element,
            const DIM_TYPE& dim) {
            // assume it's unassigned
            //if (this->mGrad->getAssigned(element.cellid))
            //	printf("WHOA assigned already!!\n");

            std::vector<INDEX_TYPE> candidates;
            std::vector<MemberDist> mdcand;
            MeshType::CofacetsIterator cofacets(this->mMesh);
            for (cofacets.begin(element.cellid); cofacets.valid(); cofacets.advance()) {
                INDEX_TYPE temp_cell = cofacets.value();
                // it's not assigned and is "lower"
                if (this->mGrad->getAssigned(temp_cell) == false &&
                    this->mGrad->getNumUnpairedFacets(temp_cell) == 1) {
                    if (element.value >= this->mFunc->cellValue(temp_cell) &&
                        element.boundary == this->mMesh->boundaryValue(temp_cell)) {
                        MemberDist res;
                        if (mdCombineAllFacetsAndDecrement(temp_cell, dim, res)) {
                            mdcand.push_back(res);
                            candidates.push_back(temp_cell);
                        }
                    }
                    else {
                        mdDecrement(temp_cell, dim);
                    }
                }
            }

            this->addNeighborsToSort(element.cellid);

            // if no candidates, we have a critical point
            if (candidates.size() == 0) {
                makeCritical(element.cellid, dim);
                return;
            }

            // so we have candidates.
            // for now pick lowest value
            //printf("have coice:\n");
            int minloc = pickFromCandidates(element.cellid, candidates, mdcand);

            this->pair(element.cellid, candidates[minloc], false);

        }

#ifdef DEBUG_SIZE_OF_MD
        std::vector<int> mdcounts;
        std::vector<int> mdsizes;
        int mdcounter;
#endif

    public:

        void set_eraser(bool val) { my_erase = val; }

        TopologicalConverventGradientBuilder(
            FuncType* mesh_function,
            MeshType* mesh_handler,
            GradType* grad_field) :
            TopologicalRegionGrowingSimpleGradientBuilder(mesh_function, mesh_handler, grad_field) {
            my_erase = true;


#ifdef DEBUG_SIZE_OF_MD
            mdcounts.resize(1024, 0);
#endif
        }

        void checkLdirVsUafacetcount(int stage0, int stage, int stage2) {
            MeshType::AllCellsIterator allit(this->mMesh);
            for (allit.begin(); allit.valid(); allit.advance()) {
                INDEX_TYPE cellid = allit.value();

                if (!this->mGrad->getAssigned(cellid)) {
                    // check values
                    int actual = 0;
                    MeshType::FacetsIterator facets(this->mMesh);
                    for (facets.begin(cellid); facets.valid(); facets.advance()) {
                        INDEX_TYPE temp_id = facets.value();
                        if (!this->mGrad->getAssigned(temp_id)) actual++;
                    }

                    int recorded = this->mGrad->getNumUnpairedFacets(cellid);

                    if (recorded != actual)
                        printf("WHOA-%d-%d-%d, %d's recorded = %d, actual %d num unpaired facets dim = %d\n",
                        stage0, stage, stage2, (int)cellid, recorded, actual,
                        this->mMesh->dimension(cellid));
                }
            }
        }


        std::vector<std::chrono::steady_clock::time_point > m_timing_vector;
#ifdef DEBUG_SIZE_OF_MD
        virtual void computeGradient(int stage=0) {
            mdcounter = 0;
#else
        virtual void computeGradient() {
#endif
            m_timing_vector.push_back(std::chrono::steady_clock::now());

            initAll();
            for (DIM_TYPE i = 0; i <= this->mMesh->maxDim(); i++) {
                std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
                this->my_dists.clear();
                this->assignDArrows(i);
                printf("Have %d remaining dists\n", my_dists.size());
#ifdef DEBUG_SIZE_OF_MD
                //check_ldir_vs_uafacetcount(stage,i,0);


                char fname[1024];
                sprintf(fname, "mdcount_%d_%d.txt", stage, i);
                FILE* fout = fopen(fname, "w");
                for (int i = 0; i < mdcounts.size(); i++) {
                    fprintf(fout, "%d %d\n", i, mdcounts[i]);
                    mdcounts[i] = 0;
                }
                fclose(fout);
                sprintf(fname, "mdsizes_%d_%d.txt", stage, i);
                fout = fopen(fname, "w");
                for (int i = 0; i < mdsizes.size(); i++) {
                    fprintf(fout, "%d %d\n", i, mdsizes[i]);
                }
                fclose(fout);
                mdsizes.clear();
                mdcounter = 0;
#endif
                this->my_dists.clear();
                this->zipUp();
                std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
                printf("took %d ms\n", std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count());
                m_timing_vector.push_back(end);

#ifdef DEBUG_SIZE_OF_MD
                //check_ldir_vs_uafacetcount(stage, i, 1);
#endif
            }
        }

#ifdef FANCY_PROB_OUTPUT
        std::vector<idfpair>& getmaxvals() {
            return maxvals;
        }
#endif





        virtual float mdDot(const MemberDist& a, const MemberDist& b) {
            float result = 0.0f;
            for (int i = 0; i < a.pairs.size(); i++) {
                for (int j = 0; j < b.pairs.size(); j++) {
                    if (a.pairs[i].id == b.pairs[j].id) {
                        result += a.pairs[i].prob * b.pairs[j].prob;
                        //break;
                    }
                }
            }
            return result;
        }

        virtual float mdDot3(const MemberDist& a, const MemberDist& b) {
            float result = 0.0f;
            for (int i = 0; i < a.pairs.size(); i++) {
                for (int j = 0; j < b.pairs.size(); j++) {
                    if (a.pairs[i].id == b.pairs[j].id) {
                        result += a.pairs[i].prob * b.pairs[j].prob;
                        //break;
                    }
                }
            }
            return result;


            //float result = 0.0f;
            //for (int i = 0; i < a.pairs.size(); i++) {
            //	for (int j = 0; j < b.pairs.size(); j++) {
            //		if (a.pairs[i].id == b.pairs[j].id && a.pairs[i].picked) {
            //			result += /*a.pairs[i].prob **/ b.pairs[j].prob;
            //			//break;
            //		}
            //	}
            //}
            //return result;
        }

        virtual float mdDot2(const MemberDist& a, const MemberDist& b) {
            int tempi = 0;
            float mv = 0;
            //float sum = a.pairs[tempi].prob;
            for (int i = 0; i < a.pairs.size(); i++) {
                for (int j = 0; j < b.pairs.size(); j++) {
                    if (a.pairs[i].id == b.pairs[j].id) {
                        float v = a.pairs[i].prob * b.pairs[j].prob;
                        if (v > mv) { tempi = i; mv = v; }
                        //sum += a.pairs[i].prob;
                        //if (a.pairs[i].prob > a.pairs[tempi].prob) tempi = i;
                    }
                }
            }
            //float result = a.pairs[tempi].prob;
            return mv; //result;
        }
        virtual float mdDist(const MemberDist& a, const MemberDist& b) {
            int tempi = 0;
            float mv = 0;
            //float sum = a.pairs[tempi].prob;
            for (int i = 0; i < a.pairs.size(); i++) {
                for (int j = 0; j < b.pairs.size(); j++) {
                    if (a.pairs[i].id == b.pairs[j].id) {
                        mv += (a.pairs[i].prob - b.pairs[j].prob)*(a.pairs[i].prob - b.pairs[j].prob);
                        //float v = a.pairs[i].prob * b.pairs[j].prob;
                        //if (v > mv) { tempi = i; mv = v; }
                        //sum += a.pairs[i].prob;
                        //if (a.pairs[i].prob > a.pairs[tempi].prob) tempi = i;
                    }
                }
            }
            //float result = a.pairs[tempi].prob;
            return sqrt(mv); //result;
        }
        virtual idfpair mdMax(MemberDist& a) {
            int size = a.pairs.size();
            if (size == 0) printf("ERRROROROROROR: mdMax has no elements\n");
            int tempi = 0;
            float sum = a.pairs[tempi].prob;
            for (int i = 1; i < size; i++) {
                sum += a.pairs[i].prob;
                if (a.pairs[i].prob > a.pairs[tempi].prob) tempi = i;
            }
            if (sum < .99f || sum > 1.001f) printf("WHOA Sum of probs = %f\n", sum);
            return a.pairs[tempi];
        }

        // res += scale*a
        virtual void mdCombine(float scale, const MemberDist& a, MemberDist& res) {
            for (int i = 0; i < a.pairs.size(); i++) {
                bool has = false;
                for (int j = 0; j < res.pairs.size(); j++) {
                    if (a.pairs[i].id == res.pairs[j].id) {
                        has = true;
                        res.pairs[j].prob += a.pairs[i].prob * scale;
                        res.pairs[j].picked = res.pairs[j].picked || a.pairs[i].picked;
                        break;
                    }
                }
                if (!has) {
                    idfpair p;
                    p.id = a.pairs[i].id;
                    p.prob = scale * a.pairs[i].prob;
                    p.picked = a.pairs[i].picked;
                    res.pairs.push_back(p);
                }
            }
        }

    };

}


#endif
