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

#ifndef GI_MORSE_SMALE_COMPLEX_RESTRICTED
#define GI_MORSE_SMALE_COMPLEX_RESTRICTED


#include <stdio.h>
#include <vector>
#include <map>
#include <queue>
#include <set>
#include <unordered_map>

#include "basic_types.h"
#include "topological_explicit_mesh_function.h"
#include "discrete_gradient_labeling.h"
#include "topological_regular_grid_restricted.h"
#include "morse_smale_complex_basic.h"
#include "msc_selectors.h"


using namespace std;

// a class that computes an MSC (in parallel) with NO arc geometry, and a hierarchy (in serial) storing the merging of manifolds
namespace MSC {

    template<typename SCALAR_TYPE, class MESH_TYPE, class FUNC_TYPE, class GRAD_TYPE>
    class MorseSmaleComplexRestricted : public MorseSmaleComplexBasic<SCALAR_TYPE, MESH_TYPE, FUNC_TYPE, GRAD_TYPE> {

    public:
        MorseSmaleComplexRestricted(GRAD_TYPE* grad, MESH_TYPE* mesh, FUNC_TYPE* func) :
            MorseSmaleComplexBasic<SCALAR_TYPE, MESH_TYPE, FUNC_TYPE, GRAD_TYPE>(grad, mesh, func)
            //MorseSmaleComplexBasic(grad, mesh, func)
        {}


    protected:

        virtual bool isValid(INT_TYPE a, arc<SCALAR_TYPE>& ap) {
            // test for boundary
            if (this->nodes[ap.lower].boundary !=
                this->nodes[ap.upper].boundary) return false;

            // 2 endpoints must be connected by exactly one arc
            if (countMultiplicity(ap, this->num_cancelled) != 1) return false;
            // test for inversions?

            return true;


        }

        virtual void InsertArcIntoSimplification(arc<SCALAR_TYPE>& a, INT_TYPE id) {

            //if (na.persistence <= gPersThreshold) {
                typename MorseSmaleComplexBasic<SCALAR_TYPE, MESH_TYPE, FUNC_TYPE, GRAD_TYPE>::sortedEdge se;

                se.persistence = a.persistence;
                se.countweight = 0; edgeCountWeight(a);
                se.ep = id;
                this->edges_to_cancel.push(se);
            //}
            //if (arcs.size() == arcs.capacity()) {


        }
        bool get_next_to_cancel_restricted(INT_TYPE& a) {
            ///printf("getnext to cancel called\n");
            while (!this->edges_to_cancel.empty()) {
                //sortedEdge se = edges_to_cancel.top();
                typename MorseSmaleComplexBasic<SCALAR_TYPE, MESH_TYPE, FUNC_TYPE, GRAD_TYPE>::sortedEdge se = this->edges_to_cancel.top();
                this->edges_to_cancel.pop();
                a = se.ep;
                arc<SCALAR_TYPE>& ap = this->arcs[a];
                ///printf("%u ", a);
                // is it alive in the current context
                if (!isAlive(ap, this->num_cancelled)) {
                    //printf("adsf1\n");
                    //printf("skip1->");
                    ///printArc(a);
                    continue;
                }

                //test if it's a valid cancellation
                if (!isValid(a, ap)) {
                    //printf("adsf2\n");
                    continue;
                }

                return true;

            }
            return false;
        }

        bool get_next_to_cancel(INT_TYPE& a) {
            ///printf("getnext to cancel called\n");
            while (!this->edges_to_cancel.empty()) {
                typename MorseSmaleComplexBasic<SCALAR_TYPE, MESH_TYPE, FUNC_TYPE, GRAD_TYPE>::sortedEdge se = this->edges_to_cancel.top();
                this->edges_to_cancel.pop();
                a = se.ep;
                arc<SCALAR_TYPE>& ap = this->arcs[a];
                ///printf("%u ", a);
                // is it alive in the current context
                if (!isAlive(ap, this->num_cancelled)) {
                    //printf("adsf1\n");
                    //printf("skip1->");
                    ///printArc(a);
                    continue;
                }

                //test if it's a valid cancellation
                if (!isValid(a, ap)) {
                    //printf("adsf2\n");
                    continue;
                }

                int newcountweight = edgeCountWeight(ap);
                if (newcountweight > se.countweight) {
                    //printf("adsf3\n");
                    se.countweight = newcountweight;
                    this->edges_to_cancel.push(se);
                    continue;
                }

                if (newcountweight > 1500) {
                    se.persistence += 1;
                    if (se.persistence <= this->gPersThreshold)
                        this->edges_to_cancel.push(se);
                    continue;
                }

                //printf("getnext to cancel returned true\n");

                return true;


            }
            //printf("getnext to cancel returned false\n");
            return false;
        }


    public:


        //void test_grow() {

        //	set<INT_TYPE> ntc;
        //	for (INT_TYPE i = 0; i < arcs.size(); i++) {
        //		arc<SCALAR_TYPE>& a = arcs[i];
        //		if (a.created != 0) continue;
        //		node<SCALAR_TYPE>& l = nodes[a.lower];
        //		node<SCALAR_TYPE>& u = nodes[a.upper];

        //		if (l.boundary)

        //	}

        //}

        virtual void ComputeHierarchy(SCALAR_TYPE pers_limit){
            this->cancel_num_to_pers.clear();
            printf(" Computing instance\n");
            this->gPersThreshold = pers_limit;
            this->max_pers_so_far = 0;
            this->num_cancelled = 0;
            // insert every arc to cancel list
            printf("  - Adding arcs to sorter...");
            INT_TYPE mysize = (INT_TYPE) this->arcs.size();
            INT_TYPE a;
            int specialcount = 0;
            for (INT_TYPE i = 0; i < mysize; i++) {
                // test if it passes the hierarchy test
                //if (!hierarchy_test->testme(i, this)) continue;
                arc<SCALAR_TYPE> &a = this->arcs[i];

                //if (a.persistence > gPersThreshold) continue;

                node<SCALAR_TYPE>& lower = getNode(a.lower);
                node<SCALAR_TYPE>& upper = getNode(a.upper);

                if (lower.boundary != upper.boundary) continue;

                typename MorseSmaleComplexBasic<SCALAR_TYPE, MESH_TYPE, FUNC_TYPE, GRAD_TYPE>::sortedEdge se;

                se.persistence = a.persistence;
                se.countweight = 0; // edgeCountWeight(a);
                se.ep = i;
                this->edges_to_cancel.push(se);
                specialcount++;
            }
            printf("found %d internal arcs\n", specialcount);

            while (get_next_to_cancel_restricted(a)) {
                //if (arcs[a].persistence > maxv) maxv = arcs[a].persistence;
                if (this->num_cancelled % 1000 == 0) {
                    printf("\r    - Cancelling: %u val=%f", this->num_cancelled, (float)this->max_pers_so_far);
                }
                //printf("\n\ncancelling: %d", num_cancelled);
                //printArc(a);

                this->cancel(a);
                this->cancel_num_to_pers.push_back(this->max_pers_so_far);

                //printf("Done\n");
            }

            printf("\r    - Cancelling finished: %u val=%f\n", this->num_cancelled, (float)this->max_pers_so_far);
            return;

            for (INT_TYPE i = 0; i < mysize; i++) {
                // test if it passes the hierarchy test
                //if (!hierarchy_test->testme(i, this)) continue;
                arc<SCALAR_TYPE> &a = this->arcs[i];

                if (a.persistence > this->gPersThreshold) continue;

                typename MorseSmaleComplexBasic<SCALAR_TYPE, MESH_TYPE, FUNC_TYPE, GRAD_TYPE>::sortedEdge se;

                se.persistence = a.persistence;
                se.countweight = edgeCountWeight(a);
                se.ep = i;
                this->edges_to_cancel.push(se);
            }

            printf("Done!\n  - Cancelling:");

            //INT_TYPE a;
            //float maxv = 0;

            while (get_next_to_cancel(a) && this->arcs[a].persistence <= this->gPersThreshold) {
                //if (arcs[a].persistence > maxv) maxv = arcs[a].persistence;
                if (this->num_cancelled % 1000 == 0) {
                    printf("\r    - Cancelling: %u val=%f", this->num_cancelled, (float)this->max_pers_so_far);
                }
                //printf("\n\ncancelling: %d", num_cancelled);
                //printArc(a);

                this->cancel(a);
                this->cancel_num_to_pers.push_back(this->max_pers_so_far);

                //printf("Done\n");
            }
            printf("\r    - Cancelling finished: %u val=%f\n", this->num_cancelled, (float)this->max_pers_so_far);



        }


    };




};
#endif


