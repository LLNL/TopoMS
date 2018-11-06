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

#ifndef MC_LIGHT_GEOM_MSC
#define MC_LIGHT_GEOM_MSC


#include <stdio.h>
#include <vector>
#include <map>
#include <queue>
#include <set>
#include <unordered_map>

#ifdef __APPLE__
#include <tr1/unordered_map>
#endif

#include "basic_types.h"
#include "topological_explicit_mesh_function.h"
#include "discrete_gradient_labeling.h"
#include "topological_regular_grid_restricted.h"




#define MAX_CELLID_VALUE 999999999
//#define WRITEREMAP( x ) (x==-1?-1:x+1)
//#define datap unsigned long long
//#define datap unsigned long long
#define UINT_INFTY -1
#define INT_INFTY 4294967295 /2
#define NULLID 4294967295/2

using namespace std;

// a class that computes an MSC (in parallel) with NO arc geometry, and a hierarchy (in serial) storing the merging of manifolds
namespace MSC {

    //typedef DiscreteGradientLabeling GRAD_TYPE;
    //typedef TopologicalRegularGrid MESH_TYPE;
    //typedef TopologicalExplicitDenseMeshFunction FUNC_TYPE;
    //typedef float SCALAR_TYPE;

    struct msbitfield
    {
        //// rehash this
        unsigned char dim : 3;
        unsigned char boundary : 1;
        unsigned char f1 : 1;
        unsigned char f2 : 1;
        unsigned char f3 : 2;
    };


    // nodes start out with -1 merged_manifolds - means that you use discrete grad to fill in the merged_manifold geometry
    // as cancellations go, we add to global merged_manifold arrays
    //
    struct merged_manifold {
        INT_TYPE merged[2]; // the merged_manifold ids
        INT_TYPE basenode; // the node id
        INT_TYPE mergetime;
    };


    // original arcs get their geometry explicitly, while new arcs from merging get different ones
    struct arc_base_geometry {
        vector<INDEX_TYPE> geometry;
    };
    struct arc_merged_geometry {
        INT_TYPE fields[3];
    };

    template<typename SCALAR_TYPE>
    struct node
    {
        INDEX_TYPE cellindex;
        INT_TYPE firstarc;
        INT_TYPE destroyed; // the time this node is cancelled.
        INT_TYPE amanifoldid; // set to -1 if this is the base
        INT_TYPE dmanifoldid; // set to -1 if this is the base
        unsigned short numarcs;
        unsigned short numlower;
        SCALAR_TYPE value;
        DIM_TYPE dim;
        BOUNDARY_TYPE boundary;
    };


    template<typename SCALAR_TYPE>
    struct arc
    {
        INT_TYPE lower; // node
        INT_TYPE lower_next; //arc - INT_INFTY is the null arc!
        INT_TYPE upper; // node
        INT_TYPE upper_next; //arc
        INT_TYPE created;
        INT_TYPE destroyed;
        INT_TYPE geom; // if created == 0, this is an original arc, and its geometry will be found in base_geom list of MSC, else in the merged_geom list
        SCALAR_TYPE persistence;
        DIM_TYPE dim; // the .dim stores the dim of the lower endpoint
        BOUNDARY_TYPE boundary;
    };

    template<typename SCALAR_TYPE, class MESH_TYPE, class FUNC_TYPE, class GRAD_TYPE>
    class MorseSmaleComplexBasic {
    public:


        typedef SCALAR_TYPE ScalarType;
        typedef MESH_TYPE MeshType;
        typedef FUNC_TYPE FuncType;
        typedef GRAD_TYPE GradType;

    protected:
        struct sortedEdge {
            SCALAR_TYPE persistence;
            int countweight;
            INT_TYPE ep;

            bool operator()(const sortedEdge& _Left, const sortedEdge& _Right) {
                if (_Left.persistence < _Right.persistence) return false;
                if (_Left.persistence > _Right.persistence) return true;
                if (_Left.countweight < _Right.countweight) return false;
                if (_Left.countweight > _Right.countweight) return true;
                return _Left.ep > _Right.ep;
            }
        };

        priority_queue<sortedEdge, vector< sortedEdge>, sortedEdge > edges_to_cancel;

        struct cancellation_record {
            int index;
            SCALAR_TYPE persistence;
            SCALAR_TYPE lval;
            SCALAR_TYPE uval;
            SCALAR_TYPE persPerc;
            INT_TYPE arcp;
            int boundary;
        };

        vector<cancellation_record> mCRecords;
        bool output_cancellation_records;

        const char* mCRecordFName;
        float m_temp_perc_pers;
        //INDEX_TYPE select_persistence;
        INDEX_TYPE num_destroyed;
        GRAD_TYPE* mGrad;
        MESH_TYPE* mMesh;
        FUNC_TYPE* mFunc;
        Vec3b mBuildArcGeometry;
        Vec3b mBuildArcAtAll;

        vector<node<SCALAR_TYPE>> nodes;
        vector<arc<SCALAR_TYPE>> arcs;
        vector<arc_base_geometry> arc_base_geoms;
        vector<arc_merged_geometry> arc_merge_geoms;
        vector<merged_manifold> mans;

#ifdef WIN32
        std::unordered_map<INDEX_TYPE, INT_TYPE> nodemap;
#else
#ifdef __APPLE__
        //__gnu_cxx::unordered_map<INDEX_TYPE, INT_TYPE> nodemap;
        tr1::unordered_map<INDEX_TYPE, INT_TYPE> nodemap;
#else
        std::unordered_map<INDEX_TYPE, INT_TYPE> nodemap;
#endif
#endif

        // this is NOT thread safe - uses references to items in a vector that could be reallocated
        inline void connectArc(INT_TYPE arcID, arc<SCALAR_TYPE>& a, node<SCALAR_TYPE>& lower, node<SCALAR_TYPE>& upper) {
            a.lower_next = lower.firstarc;
            lower.firstarc = arcID;
            a.upper_next = upper.firstarc;
            upper.firstarc = arcID;

            lower.numarcs++;
            upper.numarcs++;
            upper.numlower++;

        }

        void connectArc(INT_TYPE arcID) {

            arc<SCALAR_TYPE>& a = arcs[arcID];
            node<SCALAR_TYPE>& lower = nodes[a.lower];
            node<SCALAR_TYPE>& upper = nodes[a.upper];
            connectArc(arcID, a, lower, upper);
        }
        merged_manifold& getManifold(INT_TYPE m) {
            return mans[m];
        }

    public:

        INT_TYPE numArcs() {
            return (INT_TYPE)arcs.size();
        }
        INT_TYPE numNodes() {
            return (INT_TYPE)nodes.size();
        }
        node<SCALAR_TYPE>& getNode(INT_TYPE e) {
            return nodes[e];
        }

        arc<SCALAR_TYPE>& getArc(INT_TYPE e) {
            return arcs[e];
        }

        MESH_TYPE* const GetMesh() const {
            return mMesh;
        }

        void set_output_cancellation_records(const char* fname) {
            output_cancellation_records = true;
            mCRecordFName = fname;
        }


        MorseSmaleComplexBasic(GRAD_TYPE* grad,
            MESH_TYPE* mesh,
            FUNC_TYPE* func) :
            mGrad(grad), mMesh(mesh), mFunc(func),
            mBuildArcAtAll(Vec3b(true, true, true)), mBuildArcGeometry(Vec3b(true, true, true))
        {
            num_destroyed = 0;
            select_persistence = 0;
            output_cancellation_records = false;
        }

        // must be calleb before computefromgrad
        void SetBuildArcs(Vec3b v) { mBuildArcAtAll = v; }
        void SetBuildArcGeometry(Vec3b v) { mBuildArcGeometry = v; }

        protected:

        INT_TYPE createManifold(INT_TYPE baseId, INT_TYPE baseManId, INT_TYPE mergeManId, INT_TYPE mergetime) {
            INT_TYPE manid = mans.size();
            merged_manifold m;
            m.basenode = baseId;
            m.merged[0] = baseManId;
            m.merged[1] = mergeManId;
            m.mergetime = mergetime;
            mans.push_back(m);
            return manid;
        }

        INT_TYPE createNode(INDEX_TYPE cellID) {
            node<SCALAR_TYPE> tn;
            INT_TYPE listID = nodes.size();
            //if (nodes.size() == nodes.capacity()) printf("about to expand capacity - nodes createnode!\b");
            nodes.push_back(tn);
            nodemap[cellID] = listID;

            //node<SCALAR_TYPE> &n = nodes->operator [](listID);
            node<SCALAR_TYPE> &n = nodes[listID];
            n.cellindex = cellID;

            n.destroyed = INT_INFTY;
            n.firstarc = NULLID;
            n.boundary = mMesh->boundaryValue(cellID);

             //if (n.cellindex == 23436 || n.cellindex == 23652)
             //printf("n.cellindex = %d, n.boundary = %d\n", n.cellindex, n.boundary);
             
            //if (n.boundary) {
            //	INDEX_TYPE coords[3];
            //	sgg->getCoords(coords, cellID);
            //	printf("(%d, %d, %d) = %d\n",
            //	   (int) coords[0], (int) coords[1], (int) coords[2], n.boundary);
            //}
            n.dim = mMesh->dimension(cellID);
            n.numarcs = 0;
            n.numlower = 0;
            n.value = mFunc->cellValue(cellID);

            n.amanifoldid = createManifold(listID, -1, -1, 0);
            n.dmanifoldid = createManifold(listID, -1, -1, 0);

            return listID;
        }

        //// create an arc connecting lower to upper <- cell id's from gradient
        //INT_TYPE createArc(arc<SCALAR_TYPE> &a){
        //	INT_TYPE listID = (INT_TYPE)arcs.push_back(a);
        //	//printf("%f - %f\n", (float) a.persistence, arcs[listID].persistence);
        //	return listID;
        //}

        INT_TYPE createArc(INDEX_TYPE lowerCellID, INDEX_TYPE upperCellID, std::vector<INDEX_TYPE>& geometry) {

            arc<SCALAR_TYPE> at;
            INT_TYPE listID = arcs.size();
            arcs.push_back(at);
            INT_TYPE geomID = this->arc_base_geoms.size();
            arc_base_geoms.push_back(arc_base_geometry());

            arc<SCALAR_TYPE>& a = arcs[listID];
            arc_base_geometry& ge = arc_base_geoms[geomID];
            ge.geometry.insert(ge.geometry.begin(), geometry.begin(), geometry.end());



            a.created = 0;
            a.destroyed = INT_INFTY;
            a.lower = nodemap[lowerCellID];
            a.upper = nodemap[upperCellID];
            a.geom = geomID;
            if (nodes[a.upper].value < nodes[a.lower].value)
                printf("ERROR: upper (%f, %d) - lower (%f, %d)\n",
                float(nodes[a.upper].value), nodes[a.upper].dim,
                float(nodes[a.lower].value), nodes[a.lower].dim);

            a.persistence = nodes[a.upper].value - nodes[a.lower].value;
            a.dim = nodes[a.lower].dim;
            connectArc(listID);

            return listID;
        }

        INT_TYPE createArc(INDEX_TYPE lowerCellID, INDEX_TYPE upperCellID) {
            std::vector<INDEX_TYPE> a;
            return createArc(lowerCellID, upperCellID, a);
        }
        virtual void InsertArcIntoSimplification(arc<SCALAR_TYPE>& a, INT_TYPE id) {

            if (a.persistence <= gPersThreshold) {
                sortedEdge se;

                se.persistence = a.persistence;
                se.countweight = edgeCountWeight(a);
                se.ep = id;
                edges_to_cancel.push(se);
            }
            //if (arcs.size() == arcs.capacity()) {


        }


        INT_TYPE createArc(INT_TYPE luap, INT_TYPE ma, INT_TYPE ulap, INT_TYPE ctime) {
            arc<SCALAR_TYPE>& lua = arcs[luap];
            arc<SCALAR_TYPE>& ula = arcs[ulap];
            INT_TYPE na_id = arcs.size();
            INT_TYPE na_geom_id = arc_merge_geoms.size();
            arc_merge_geoms.push_back(arc_merged_geometry());
            arc_merged_geometry& na_geom = arc_merge_geoms[na_geom_id];
            na_geom.fields[0] = luap;
            na_geom.fields[1] = ma;
            na_geom.fields[2] = ulap;

            //this->arcs.push_back(ma); // copy setting from old
            arc<SCALAR_TYPE> na;// = this->arcs[na_id];

            na.created = ctime;
            na.destroyed = INT_INFTY;
            na.lower = ula.lower;
            na.upper = lua.upper;
            na.geom = na_geom_id;

            node<SCALAR_TYPE>& nup = this->nodes[na.upper];
            node<SCALAR_TYPE>& nlo = this->nodes[na.lower];

            this->connectArc(na_id, na, nlo, nup);
            na.persistence = nup.value - nlo.value;
            na.boundary = nlo.boundary + nup.boundary;
            na.dim = nlo.dim;

            if (na.persistence < 0 && nlo.boundary == nup.boundary){
                //printf("creatinga inversion %f, %d-%d, b%d-b%d\n", (float)na.persistence, nlo.dim, nup.dim, nlo.boundary, nup.boundary);
            }

            // now insert into global sort if it has a chance of being cancelled
            InsertArcIntoSimplification(na, na_id);
            //	printf("about to expand capacity2 - arcs createnode!\n");
            //}
            this->arcs.push_back(na); // this could destroy references so put at end;
            return na_id;

        }
        //void rec_tdcr_no_geom(const INDEX_TYPE& cellid, DIM_TYPE& temp_dim, const INDEX_TYPE start) {
        //	INDEX_TYPE current = cellid;
        //	MESH_TYPE::FacetsIterator facets(mMesh);
        //	for (facets.begin(current); facets.valid(); facets.advance()) {
        //		INDEX_TYPE temp_id = facets.value();
        //		if (mGrad->getCritical(temp_id)) {
        //			//printf("adding arc: %llu\n", temp_id);
        //			createArc(temp_id, start);

        //		}
        //		else if (mGrad->getDimAscMan(temp_id) == temp_dim) {
        //			INDEX_TYPE pair = mGrad->getPair(temp_id);
        //			if (pair != cellid && mMesh->dimension(pair) == mMesh->dimension(cellid)) {

        //				//result.push_back(pair);
        //				rec_tdcr_no_geom(pair, temp_dim, start);

        //			}
        //		}

        //	}
        //}
        //void trace_down_cells_restricted_no_geom(const INDEX_TYPE& cellid) {

        //	DIM_TYPE temp_dim = mGrad->getDimAscMan(cellid) + 1;
        //	rec_tdcr_no_geom(cellid, temp_dim, cellid);
        //}

        void rec_tdcr(const INDEX_TYPE& cellid, DIM_TYPE& temp_dim, const INDEX_TYPE start, vector<INDEX_TYPE>& geom) {
            INDEX_TYPE current = cellid;

            geom.push_back(current);
            typename MESH_TYPE::FacetsIterator facets(mMesh);
            for (facets.begin(current); facets.valid(); facets.advance()) {
                INDEX_TYPE temp_id = facets.value();
                geom.push_back(temp_id);
                if (mGrad->getCritical(temp_id)) {
#pragma omp critical
                    {
                        createArc(temp_id, start, geom);
                    }

                }
                else if (mGrad->getDimAscMan(temp_id) == temp_dim) {
                    INDEX_TYPE pair = mGrad->getPair(temp_id);
                    if (pair != cellid && mMesh->dimension(pair) == mMesh->dimension(cellid)) {

                        rec_tdcr(pair, temp_dim, start, geom);

                    }
                }
                geom.pop_back();

            }
            geom.pop_back();
        }

        void rec_tdcr(const INDEX_TYPE& cellid, DIM_TYPE& temp_dim, const INDEX_TYPE start) {
            INDEX_TYPE current = cellid;

            typename MESH_TYPE::FacetsIterator facets(mMesh);
            for (facets.begin(current); facets.valid(); facets.advance()) {
                INDEX_TYPE temp_id = facets.value();
                if (mGrad->getCritical(temp_id)) {
#pragma omp critical
                    {
                        createArc(temp_id, start);
                    }

                }
                else if (mGrad->getDimAscMan(temp_id) == temp_dim) {
                    INDEX_TYPE pair = mGrad->getPair(temp_id);
                    if (pair != cellid && mMesh->dimension(pair) == mMesh->dimension(cellid)) {

                        rec_tdcr(pair, temp_dim, start);

                    }
                }

            }
        }

        void trace_down_cells_restricted(INDEX_TYPE cellid) {
            vector<INDEX_TYPE> geom;
            DIM_TYPE temp_dim = mGrad->getDimAscMan(cellid) + 1;
            rec_tdcr(cellid, temp_dim, cellid, geom);
        }
        void trace_down_cells_restricted_nogeom(INDEX_TYPE cellid) {
            DIM_TYPE temp_dim = mGrad->getDimAscMan(cellid) + 1;
            rec_tdcr(cellid, temp_dim, cellid);
        }

        void AddArcs(INT_TYPE nodeId) {
            node<SCALAR_TYPE>& n = getNode(nodeId);
            if (n.dim == 0) return;
            if (mBuildArcAtAll[n.dim - 1]) {
                if (mBuildArcGeometry[n.dim - 1]) {
                    trace_down_cells_restricted(nodes[nodeId].cellindex);
                }
                else {
                    trace_down_cells_restricted_nogeom(nodes[nodeId].cellindex);
                }
            }
        }

        public:
        void ComputeFromGrad(bool restricted = false) {
            typename MESH_TYPE::AllCellsIterator t_cells(mMesh);
            printf("   -- Finding critical points...");
            fflush(stdout);
            int count = 0;
            int counts[4] = {0,0,0,0};
            for (t_cells.begin(); t_cells.valid(); t_cells.advance()) {
                INDEX_TYPE t_id = t_cells.value();
                if (mGrad->getCritical(t_id)) {
                    createNode(t_id);
                    count++;

                    counts[ mMesh->dimension(t_id) ]++;
                }
            }
            printf(" done! found %d critical points (%d %d %d %d)\n", count, counts[0], counts[1], counts[2], counts[3]);
            printf("   -- Finding arcs...");
            fflush(stdout);
            // then add arcs
            INT_TYPE numnodes = nodes.size();

#pragma omp parallel for
            for (int j = 0; j < numnodes; j++) {
                AddArcs(j);
            }
            printf(" done! found %d arcs\n", this->arcs.size());
        }



    protected:

        //////////// NOW CANCELLATION STUFF

        int countUpperArcs(INT_TYPE id) {
            return this->nodes[id].numarcs - this->nodes[id].numlower;
        }
        int countLowerArcs(INT_TYPE id) {
            return this->nodes[id].numlower;
        }

        int edgeCountWeight(arc<SCALAR_TYPE> &a) {
            INT_TYPE lower = a.lower;
            INT_TYPE upper = a.upper;
            int nl = countUpperArcs(lower);
            int nu = countLowerArcs(upper);

            ////first option: return the number created
            return (nl - 1) * (nu - 1);

            //2nd option: return the change in the number of arcs
            return (nl - 1) * (nu - 1) - (nl + nu - 1);

            //return 1;
        }

        inline INT_TYPE current_pers() { return select_persistence; }
        inline bool isAlive(arc<SCALAR_TYPE>&a, INT_TYPE place) {
            return a.created <= place && a.destroyed > place;
        }
    public:
        bool isNodeAlive(INT_TYPE n) {
            return isAlive(this->nodes[n], select_persistence);
        }
        bool isArcAlive(INT_TYPE a) {
            bool res = isAlive(this->arcs[a], select_persistence);
            if (res && !isNodeAlive(this->arcs[a].lower)) printf("ERROR LOWER NODE NOT ALIVE\n");
            if (res && !isNodeAlive(this->arcs[a].upper)) printf("ERROR UPPER NODE NOT ALIVE\n");
            return isAlive(this->arcs[a], select_persistence);
        }

    protected:
        inline bool isAlive(node<SCALAR_TYPE>&n, INT_TYPE place) {
            //printf("n.destroyed = %d, place = %d\n", n.destroyed, place);
            return  n.destroyed > place;
        }
        inline bool isAlive(INT_TYPE n, INT_TYPE place) {
            return isAlive(this->nodes[n], place);
        }


    public:
        void SetSelectPersAbs(SCALAR_TYPE value) {
            //printf("mcLightGeomMSC::mcSuperLightMSC::SetSelectPersAbs -> %f\n", value);
            int offset = 1;// cancel_num_to_pers.size() - 1;
            for (int i = 0; i < cancel_num_to_pers.size(); i += offset) {
                if (cancel_num_to_pers[i] > value) {
                    select_persistence = i;

                    return;
                }
                //while (i + offset > cancel_num_to_pers.size() - 1) offset = offset / 2;
                //while (offset > 1 && cancel_num_to_pers[i + offset] > value) offset = offset / 2;
            }
            select_persistence = cancel_num_to_pers.size() - 1;
        }

        INT_TYPE nextArc(arc<SCALAR_TYPE>& ap, INT_TYPE n) {
            if (ap.lower == n) return ap.lower_next;
            if (ap.upper == n) return ap.upper_next;
            return 0;
        }
    protected:
        int cancel(INT_TYPE a) {

            int createcounter = 0;
            int deletecounter = 0;

            //int tmp1 = (int) edges_to_cancel.size();

            arc<SCALAR_TYPE>* ap = &(this->arcs[a]);
            node<SCALAR_TYPE>& lower = this->nodes[ap->lower];
            node<SCALAR_TYPE>& upper = this->nodes[ap->upper];

            //if (output_cancellation_records) {
            cancellation_record cr;
            cr.index = lower.dim;
            cr.lval = lower.value;
            cr.uval = upper.value;
            cr.persistence = ap->persistence;
            //SCALAR_TYPE diff = this->sgg->sgd->maxval - this->sgg->sgd->minval;
            //cr.persPerc = 100.0f * ap->persistence / diff;
            cr.arcp = a;
            cr.boundary = lower.boundary + upper.boundary;
            this->mCRecords.push_back(cr);
            //}



            if (lower.destroyed != INT_INFTY) printf("Error: MorseSmaleComplexBasic::cancel(%d) - lower.destroyed != INT_INFTY\n", a);
            if (upper.destroyed != INT_INFTY) printf("Error: MorseSmaleComplexBasic::cancel(%d) - upper.destroyed != INT_INFTY\n", a);
            //printf("\n");
            //printNode(ap->lower);
            //printNode(ap->upper);

            int initialguess = (lower.numarcs - lower.numlower - 1) * (upper.numlower - 1);
            int init2 = lower.numarcs + upper.numarcs - 1;

            if (ap->persistence > max_pers_so_far) max_pers_so_far = ap->persistence;
            num_cancelled++;


            // for each upwards arc connected to the lower node,
            // for each downwards arc connected to the upper node,
            // create a new connection from the other endpoints
            INT_TYPE la = lower.firstarc;
            while (la != NULLID) {
                arc<SCALAR_TYPE>& lap = this->arcs[la];

                // skip the arc itself
                if (la == a) {
                    la = lap.lower_next; // we're guarantee this is the next
                    continue;
                }

                //test if they are the same kind
                if (ap->lower != lap.lower) {
                    la = lap.upper_next;
                    continue;
                }
                if (!isAlive(lap, num_cancelled)) {
                    la = lap.lower_next;
                    continue;
                }

                INT_TYPE ua = upper.firstarc;
                while (ua != NULLID) {
                    arc<SCALAR_TYPE>& uap = this->arcs[ua];

                    // skip the arc itself
                    if (ua == a) {
                        ua = uap.upper_next; // we're guarantee this is the next
                        continue;
                    }

                    //test if they are the same kind
                    if (ap->upper != uap.upper) {
                        ua = uap.lower_next;
                        continue;
                    }
                    if (!isAlive(uap, num_cancelled)) {
                        ua = uap.upper_next;
                        continue;
                    }

                    // create the arc here!
                    INT_TYPE newarc = createArc(la, a, ua, num_cancelled);
                    //printf("ha1\n");
                    ap = &(this->arcs[a]);
                    //printf("ha\n");
                    createcounter++;

                    ua = this->arcs[ua].upper_next;
                    //printf("asdf\n");
                }

                la = this->arcs[la].lower_next;

            }
            //printf("hahan\n");
            // go through and set the delete time on the arcs that are to be removed,
            // and update the arc counters at the other endpoints

            // following comments are for 0-1 cancellations
            // first look in neighborhood of minimum = lower
            la = lower.firstarc;							// pick first arc attached to min
            while (la != NULLID) {							// while there are more arcs in neighborhood, keep looking
                arc<SCALAR_TYPE>& lap = this->arcs[la];			// get the reference to the arc attached to the min

                // skip the arc itself
                if (la == a) {
                    la = lap.lower_next; // we're guarantee this is the next
                    continue;
                }

                //test if they are the same kind
                if (ap->lower == lap.lower) {				// we do not want arcs that are not 0-1, so do this check
                    if (isAlive(lap, num_cancelled)) {		// make sure we are only working with living arcs
                        node<SCALAR_TYPE>& n = this->nodes[lap.upper]; // now "n" is the 1-saddle attached to minimum

                        n.numarcs--;						// if we remove this arc, we need to decrement the number of living arcs attached to it
                        n.numlower--;						// this cancellation will also remove a "lower" arc of the 1-saddle

                        // create merged_manifold for merging
                        // extend the descending merged_manifold of the 1-saddle by merging it with the descending merged_manifold of the cancelled 1-saddle
                        INT_TYPE nmanid =
                            createManifold(lap.upper, n.dmanifoldid, upper.dmanifoldid, num_cancelled);
                        n.dmanifoldid = nmanid;				// change the merged_manifold reference of the 1-saddle to reference this new decending merged_manifold

                        lap.destroyed = num_cancelled;		// we remove the old arc to the deleted min, so mark it destroyed with the num_cancelld
                        deletecounter++;					// keep track of deleted arc cound
                    }
                    la = lap.lower_next;					// continue on the next arc around the minimum
                }
                else if (ap->lower == lap.upper) {          // if this is in fact an i-1 - i arc, we just remove the arcs around it - no merging happens
                    if (isAlive(lap, num_cancelled)) {
                        node<SCALAR_TYPE>& n = this->nodes[lap.lower];
                        n.numarcs--;
                        lap.destroyed = num_cancelled;
                        deletecounter++;
                    }
                    la = lap.upper_next;
                }
                else  {
                    printf("ERROR SHOULD NEVER GET HERE\n");
                }



            }
            INT_TYPE ua = upper.firstarc;
            while (ua != NULLID) {
                arc<SCALAR_TYPE>& uap = this->arcs[ua];

                // skip the arc itself
                if (ua == a) {
                    ua = uap.upper_next; // we're guarantee this is the next
                    continue;
                }
                //test if they are the same kind
                if (ap->upper == uap.upper) {
                    if (isAlive(uap, num_cancelled)) {

                        node<SCALAR_TYPE>& n = this->nodes[uap.lower]; // pick the other

                        n.numarcs--;

                        // create merged_manifold for merging
                        INT_TYPE nmanid = createManifold(uap.lower, n.amanifoldid, lower.amanifoldid, num_cancelled);
                        n.amanifoldid = nmanid;

                        uap.destroyed = num_cancelled;
                        deletecounter++;
                    }
                    ua = uap.upper_next;
                }
                else if (ap->upper == uap.lower) {
                    if (isAlive(uap, num_cancelled)) {
                        node<SCALAR_TYPE>& n = this->nodes[uap.upper];
                        n.numarcs--;
                        n.numlower--;
                        uap.destroyed = num_cancelled;
                        deletecounter++;
                    }
                    ua = uap.lower_next;
                }
                else  {
                    printf("ERROR SHOULD NEVER GET HERE\n");
                }


            }

            ap->destroyed = num_cancelled;
            lower.destroyed = num_cancelled;
            upper.destroyed = num_cancelled;
            deletecounter++;

            //int tmp2 = (int) edges_to_cancel.size();

            //printf("initialguess: %d actual:%d del:%d td:%d tmp1:%d tmp2:%d \n", initialguess, createcounter, deletecounter,init2
            //	,tmp1, tmp2);
            //printf("gothere\n");
            return 1;
        }




        int countMultiplicity(arc<SCALAR_TYPE>& ap, INT_TYPE ctime) {
            INT_TYPE nu = ap.upper;
            INT_TYPE nl = ap.lower;
            INT_TYPE a = this->nodes[nu].firstarc;
            int counter = 0;
            while (a != NULLID){
                arc<SCALAR_TYPE>& nap = this->arcs[a];
                if (isAlive(nap, ctime) && nap.lower == nl && nap.upper == nu) {
                    counter++;
                }
                a = nextArc(nap, nu);
            }
            return counter;
        }

        virtual bool isValid(INT_TYPE a, arc<SCALAR_TYPE>& ap) {
            // test for boundary
            if (this->nodes[ap.lower].boundary !=
                this->nodes[ap.upper].boundary) return false;

            // 2 endpoints must be connected by exactly one arc
            if (countMultiplicity(ap, num_cancelled) != 1) return false;
            // test for inversions?

            return true;


        }

        // THIS WILL HAVE TO WORK LATER, but for now only need 0 and maxdim manifolds
        //void recCounLeaftManifolds(INT_TYPE mId, map< INT_TYPE, int >& counter) {
        //	merged_manifold &m = this->getManifold(mId);
        //	if (m->merge[0] != NULL) {
        //		recCounLeaftManifolds(m->merge[0], counter);
        //		recCounLeaftManifolds(m->merge[1], counter);
        //	}
        //	else {
        //		if (counter.count(m) == 0) {
        //			counter[m] = 1;
        //		}
        //		else {
        //			counter[m]++;
        //		}
        //	}
        //}

        void printmanifold(INT_TYPE man) {
            merged_manifold& m = getManifold(man);
            printf("man=%d, man.base=%d, man.mergetime=%d, man.merge[0]=%d, man.merge[1]=%d\n",
                man, m.basenode, m.mergetime, m.merged[0], m.merged[1]);
        }

        INT_TYPE getActiveMan(INT_TYPE man) {
            //printmanifold(man);
            while (getManifold(man).mergetime > this->current_pers()) {
                man = getManifold(man).merged[0];
                //printf("  --");  printmanifold(man);
            }
            //printf("ret %d\n", man);
            return man;
        }

        void GatherNodes(INT_TYPE nodeID, set<INT_TYPE>& res, bool ascending) {

            if (!isAlive(nodeID, this->select_persistence)) {
                return;
            }
            node<SCALAR_TYPE>& n = this->getNode(nodeID);
            INT_TYPE man;
            if (ascending) {
                man = n.amanifoldid;
            }
            else {
                man = n.dmanifoldid;
            }
            INT_TYPE manID = getActiveMan(man);
            recGatherNodes(manID, res);

        }

    protected:
        void recGatherNodes(INT_TYPE mid, set<INT_TYPE>& res) {
            vector<INT_TYPE> expand;
            expand.push_back(mid);
            while (expand.size() > 0) {

                INT_TYPE curr = expand.back(); expand.pop_back();

                merged_manifold& m = getManifold(curr);
                if (m.merged[0] != -1) {
                    expand.push_back(m.merged[0]);
                    expand.push_back(m.merged[1]);
                }
                res.insert(m.basenode);


            }


        }

        void rec_man_trace_up(INDEX_TYPE& cellid, set<INDEX_TYPE>& res) {
            res.insert(cellid);
            INDEX_TYPE current = cellid;
            DIM_TYPE cdim = mMesh->dimension(cellid);
            typename MESH_TYPE::CofacetsIterator cofacets(mMesh);
            for (cofacets.begin(current); cofacets.valid(); cofacets.advance()) {
                INDEX_TYPE temp_id = cofacets.value();
                if (mGrad->getCritical(temp_id)) continue;

                INDEX_TYPE temp_pair = mGrad->getPair(temp_id);

                if (temp_pair == cellid) continue;

                if (mMesh->dimension(temp_pair) != cdim) continue;

                rec_man_trace_up(temp_pair, res);
            }
        }
        void rec_man_trace_down(INDEX_TYPE& cellid, set<INDEX_TYPE>& res) {
            res.insert(cellid);
            INDEX_TYPE current = cellid;
            DIM_TYPE cdim = mMesh->dimension(cellid);
            typename MESH_TYPE::FacetsIterator facets(mMesh);
            for (facets.begin(current); facets.valid(); facets.advance()) {
                INDEX_TYPE temp_id = facets.value();
                if (mGrad->getCritical(temp_id)) continue;

                INDEX_TYPE temp_pair = mGrad->getPair(temp_id);

                if (temp_pair == cellid) continue;

                if (mMesh->dimension(temp_pair) != cdim) continue;

                rec_man_trace_down(temp_pair, res);
            }
        }

        void fillUnsimplifiedGeometry(INT_TYPE nodeID, set<INDEX_TYPE>& res, bool ascending) {

            node<SCALAR_TYPE>& n = this->getNode(nodeID);

            if (ascending){
                rec_man_trace_up(n.cellindex, res);
            }
            else  {
                rec_man_trace_down(n.cellindex, res);
            }

        }

    public:
        void fillGeometry(INT_TYPE nodeID, set<INDEX_TYPE>& res, bool ascending) {

            if (!isNodeAlive(nodeID)) return;

            set<INT_TYPE> nodeset;
            //printf("gothere1\n");
            GatherNodes(nodeID, nodeset, ascending);
            //printf("gothere2\n");

            for (set<INT_TYPE>::iterator it = nodeset.begin(); it != nodeset.end(); it++) {
                node<SCALAR_TYPE>& n = this->getNode(*it);

                if (ascending){
                    rec_man_trace_up(n.cellindex, res);
                }
                else  {
                    rec_man_trace_down(n.cellindex, res);
                }

            }
            //printf("gothere3\n");



        }




    protected:

        bool get_next_to_cancel(INT_TYPE& a) {
            ///printf("getnext to cancel called\n");
            while (!edges_to_cancel.empty()) {
                sortedEdge se = edges_to_cancel.top();
                edges_to_cancel.pop();
                a = se.ep;
                arc<SCALAR_TYPE>& ap = this->arcs[a];
                ///printf("%u ", a);
                // is it alive in the current context
                if (!isAlive(ap, num_cancelled)) {
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
                    edges_to_cancel.push(se);
                    continue;
                }

                if (newcountweight > 1500) {
                    se.persistence += 1;
                    if (se.persistence <= gPersThreshold)
                        edges_to_cancel.push(se);
                    continue;
                }

                //printf("getnext to cancel returned true\n");

                return true;


            }
            //printf("getnext to cancel returned false\n");
            return false;
        }

        vector<SCALAR_TYPE> cancel_num_to_pers;
        INT_TYPE select_persistence; // the persistence to select in hierarchy
        INT_TYPE num_cancelled;

        SCALAR_TYPE max_pers_so_far;
        SCALAR_TYPE gPersThreshold;


    public:

        virtual void ComputeHierarchy(SCALAR_TYPE pers_limit){
            cancel_num_to_pers.clear();
            printf(" -- Performing cancellation to %f...\n", pers_limit);
            gPersThreshold = pers_limit;
            max_pers_so_far = 0;
            num_cancelled = 0;
            // insert every arc to cancel list
            printf("  -- Adding arcs to sorter...");
            INT_TYPE mysize = (INT_TYPE) this->arcs.size();
            for (INT_TYPE i = 0; i < mysize; i++) {
                // test if it passes the hierarchy test
                //if (!hierarchy_test->testme(i, this)) continue;
                arc<SCALAR_TYPE> &a = this->arcs[i];

                if (a.persistence > gPersThreshold) continue;

                sortedEdge se;

                se.persistence = a.persistence;
                se.countweight = edgeCountWeight(a);
                se.ep = i;
                edges_to_cancel.push(se);
            }

            printf("Done!\n");
            printf("  -- Cancelling:");

            INT_TYPE a;
            //float maxv = 0;

            while (get_next_to_cancel(a) && this->arcs[a].persistence <= gPersThreshold) {
                //if (arcs[a].persistence > maxv) maxv = arcs[a].persistence;
                if (num_cancelled % 1000 == 0) {
                    printf("\r  -- Cancelling: %u val=%f", num_cancelled, (float)max_pers_so_far);
                }
                //printf("\n\ncancelling: %d", num_cancelled);
                //printArc(a);

                cancel(a);
                cancel_num_to_pers.push_back(max_pers_so_far);

                //printf("Done\n");
            }
            printf("\r  -- Cancelling finished. num_cancelled = %u, max_persistence = %f\n", num_cancelled, (float)max_pers_so_far);


            //_evercomputed = true;
            if (output_cancellation_records) {
                int mincount = 0;// this->countIDs[0];
                int maxcount = 0;//this->countIDs[3];
                int interiormincount = 0;// this->countInteriorIDs[0];
                FILE* fout = fopen(this->mCRecordFName, "w");
                for (int i = 0; i < this->mCRecords.size(); i++) {
                    cancellation_record& cr = this->mCRecords[i];
                    fprintf(fout, "%f %f %f %f %d %d %d %d %d\n", cr.persistence,
                        cr.lval, cr.uval, cr.persPerc, cr.index, mincount, maxcount,
                        interiormincount, cr.boundary);
                    if (cr.index == 0) {
                        mincount--;
                        node<SCALAR_TYPE>& n =
                            this->getNode(this->getArc(cr.arcp).lower);
                        if (!n.boundary) {
                            interiormincount--;
                        }
                    }
                    if (cr.index == 2) maxcount--;
                }
                fclose(fout);
            }
            // now fill in the lists?
        }
        void Output1SaddleRecord(const char* fname) {

            FILE* fsadrec = fopen(fname, "w");

            for (int i = 0; i < nodes.size(); i++) {
                node<SCALAR_TYPE>& n = nodes[i];
                if (n.dim > 1) continue;
                float persval =
                    (n.destroyed < cancel_num_to_pers.size() + 1 ?
                    cancel_num_to_pers[n.destroyed - 1] :
                    cancel_num_to_pers[cancel_num_to_pers.size() - 1]);
                fprintf(fsadrec, "%d %d %llu %f %d %f\n", n.dim, n.boundary, n.cellindex, n.value, n.destroyed, persval);


            }
            fclose(fsadrec);


        }
    protected:
        void _fillArcGeometry(INT_TYPE aid, vector<INDEX_TYPE>& v, bool direction) {
            arc<SCALAR_TYPE>& a = this->arcs[aid];
            //printf("gothere1\n");
            if (a.created == 0) {
                // this is base
                arc_base_geometry& base = this->arc_base_geoms[a.geom];
                if (direction) {
                    for (int i = 0; i < base.geometry.size(); i++) {
                        if (v.size() > 0 && v[v.size() - 1] == base.geometry[i]){
                            v.pop_back();
                        }
                        else{
                            v.push_back(base.geometry[i]);
                        }
                    }
                }
                else {
                    for (int i = base.geometry.size() - 1; i >= 0; i--) {
                        if (v.size() > 0 && v[v.size() - 1] == base.geometry[i]){
                            v.pop_back();
                        }
                        else{
                            v.push_back(base.geometry[i]);
                        }
                    }
                }
                return;
            }
            //printf("gothere2\n");
            // recurse on children
            arc_merged_geometry& m = arc_merge_geoms[a.geom];
            if (direction) {
                _fillArcGeometry(m.fields[0], v, true);
                _fillArcGeometry(m.fields[1], v, false);
                _fillArcGeometry(m.fields[2], v, true);
            }
            else {
                _fillArcGeometry(m.fields[2], v, false);
                _fillArcGeometry(m.fields[1], v, true);
                _fillArcGeometry(m.fields[0], v, false);
            }


        }
    public:
        INT_TYPE arcLowerNode(INT_TYPE aid) {
            return getArc(aid).lower;
        }
        INT_TYPE arcUpperNode(INT_TYPE aid) {
            return getArc(aid).upper;
        }
        INT_TYPE nextIncidentLivingArc(INT_TYPE aid, INT_TYPE nid) {
            //printf("%d-\n", aid);
            INT_TYPE naid = nextArc(getArc(aid), nid);
            //printf("  -%d\n", naid);
            while (naid != NULLID) {
                if (isArcAlive(naid)) return naid;
                naid = nextArc(getArc(naid), nid);
                //printf("  -%d\n", naid);
            }
            return naid;
        }
        bool isValidArcId(INT_TYPE aid) {
            return aid != NULLID;
        }

        INT_TYPE firstIncidentLivingArc(INT_TYPE nid) {
            node<SCALAR_TYPE>& n = getNode(nid);
            INT_TYPE aid = n.firstarc;
            if (!isArcAlive(aid)) return nextIncidentLivingArc(aid, nid);
            return aid;
        }

        void fillArcGeometry(INT_TYPE aid, vector<INDEX_TYPE>& v) {
            //v.push_back(nodes[arcs[aid].upper].cellindex);
            v.clear();
            _fillArcGeometry(aid, v, true);
            //v.push_back(nodes[arcs[aid].lower].cellindex);
        }
    };




}   // end of namespace
#endif
