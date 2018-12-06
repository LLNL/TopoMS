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

#ifndef ROBIN_LABELING_H
#define ROBIN_LABELING_H

#pragma once

#include "basic_types.h"
#include "regular_grid.h"
#include "labeling.h"
#include "topological_regular_grid.h"
#include "discrete_gradient_labeling.h"
#include "regular_grid_trilinear_function.h"
#include "topological_explicit_mesh_function.h"
#include "topological_regular_grid_restricted.h"

#include "timing.h"

#include <queue>
#include <vector>
#include <algorithm>
#include <fstream>



// -------------------------------------------------------------------------------

#define MAX_VERTS_PER_CELL 8
#define MAX_FACES_PER_CELL 6
#define MAX_COFACES_PER_CELL 8

#define MAX_FACES_PER_STAR 3
#define MAX_CELLS_PER_STAR 27


// -------------------------------------------------------------------------------
namespace MSC {

    // -------------------------------------------------------------------------------
    // comparison of cells based on the max vertex
    template <class MeshFunction>
    struct Cell_comparator {

        static MeshFunction const* ms_topofunc;

        static INDEX_TYPE vid(const INDEX_TYPE &a) { return ms_topofunc->mesh()->VertexNumberFromCellID(a); }
        static float value(const INDEX_TYPE &a) { return ms_topofunc->cellValue(a); }

        static bool is_lower(const INDEX_TYPE &a, const INDEX_TYPE &b) {
            return ms_topofunc->lessThan(a, b);
        }
        static bool is_greater(const INDEX_TYPE &a, const INDEX_TYPE &b) {
            return ms_topofunc->lessThan(b, a);
        }
    };

    // need to deal with this better!
    template<class MeshFunction>
    MeshFunction const* Cell_comparator<MeshFunction>::ms_topofunc = 0;

    // -------------------------------------------------------------------------------
    // comparision of pointers based on their values (used for priority queue)
    template <class T>
    struct Pointer_Comparator {

        static bool equals(const T* a, const T* b) { return *a == *b; }
        static bool greater_than(const T* a, const T* b) { return *a > *b; }
        static bool less_than(const T* a, const T* b) { return *a < *b; }

        // default behavior
        bool operator()(const T* a, const T* b) const { return greater_than(a, b); }
    };

    // -------------------------------------------------------------------------------
    // priority queue of pointers
    template <class T>
    class PriorityQueue : public std::priority_queue<const T*, std::vector<const T*>, Pointer_Comparator<T> > {
    public:

        const T* pop_front() {

            const T* t = this->top();
            this->pop();
            return t;
        }

        void push_unique(const T* x) {
#ifdef DEBUG_MODE
            if (exists(x)) {
                printf(" push2 error: element already exists :"); x->print("x"); this->print("Queue");
                exit(1);
            }
#endif
            this->push(x);
        }
        void remove(const T* x) {

            std::vector<const T*> backup;
            backup.reserve(this->size());

            bool found = false;

            // back the queue up until you find what you need
            while (!this->empty()) {

                const T* t = this->pop_front();
                if (Pointer_Comparator<T>::equals(t, x)) {
                    found = true;
                    break;
                }
                backup.push_back(t);
            }

            // now, put all backup elements back in the queue
            for (int i = 0; i < backup.size(); i++) {
                this->push(backup[i]);
            }
#ifdef DEBUG_MODE
            if (!found) {
                printf(" remove error: element not found: ");   x->print("x");
                exit(1);
            }
#endif
        }
        bool exists(const T* x) const {

            PriorityQueue<T> p = *this;

            while (!p.empty()) {

                const T* t = p.pop_front();
                if (Pointer_Comparator<T>::equals(t, x)) {
                    return true;
                }
            }
            return false;
        }

        bool validate() const {

            PriorityQueue<T> p = *this;

            int btype = -1;
            while (!p.empty()) {

                const T* t = p.pop_front();
                if (btype == -1) {
                    btype = t->boundary();
                }
                else if (btype != t->boundary()) {
                    return false;
                }
            }
            return true;
        }

        void print(const char *queuename) const {

            PriorityQueue<T> p = *this;

            printf(" Displaying queue %s of size %d\n", queuename, p.size());
            while (!p.empty()) {

                const T *element = p.pop_front();
                element->print("C");
            }
        }
    };

    // -------------------------------------------------------------------------------

    // store the lower star of a vertex
    template <class Mesh, class MeshFunction>
    class LowerStar;

    template <class Mesh, class MeshFunction>
    class MCell {

    private:
        INDEX_TYPE m_cellid;
        Vec3l m_cell_coords;
        DIM_TYPE m_dim;
        BOUNDARY_TYPE m_bdy;

        // store all vertices of this cell
        DIM_TYPE m_n0faces;
        INDEX_TYPE m_0faces[MAX_VERTS_PER_CELL];

        // store only the faces that are in the lower star
        DIM_TYPE m_nfaces_in_lstar;
        const MCell<Mesh, MeshFunction>* m_faces_in_lstar[MAX_FACES_PER_CELL];

    public:

        MCell() : m_cellid(0), m_dim(0) {}

        void initialize(INDEX_TYPE cellid, const MeshFunction *topofunc,
            const Mesh *tgrid) {

            m_cellid = cellid;
            m_dim = 0;

            tgrid->cellid2Coords(cellid, m_cell_coords);

            // G ordering of cells: will be done based upon the function values of corners
            m_n0faces = 0;
            typename Mesh::CellVerticesIterator cit(tgrid);
            for (cit.begin(m_cell_coords); cit.valid(); cit.advance()) {
                m_0faces[m_n0faces++] = cit.value();
            }
            std::sort(m_0faces, m_0faces + m_n0faces, Cell_comparator<MeshFunction>::is_greater);
            for (int i = m_n0faces; i < MAX_VERTS_PER_CELL; i++) {
                m_0faces[i++] = -1;
            }

            size_t sz = m_n0faces;
            while (sz >>= 1) { m_dim++; }

            m_bdy = tgrid->boundaryValue(m_cellid);
        }

        void populate_faces(const LowerStar<Mesh, MeshFunction> &lstar);

        inline DIM_TYPE dim() const { return m_dim; }
        inline INDEX_TYPE cellid() const { return m_cellid; }

        inline DIM_TYPE nfaces() const { return m_nfaces_in_lstar; }
        inline const MCell<Mesh, MeshFunction>* face(DIM_TYPE i) const { return m_faces_in_lstar[i]; }

        inline bool is_cofaceOf(const MCell<Mesh, MeshFunction> *other) const {
            return std::find(m_faces_in_lstar, m_faces_in_lstar + m_nfaces_in_lstar, other) != m_faces_in_lstar + m_nfaces_in_lstar;
        }

        inline BOUNDARY_TYPE boundary() const { return this->m_bdy; }
        inline bool boundary_match(const MCell<Mesh, MeshFunction> *other) const {
            return (this->boundary() == other->boundary());

        }

        inline bool is_candidate_for_pairing(const MCell<Mesh, MeshFunction> *other) const {


            return (this->is_cofaceOf(other) && this->boundary_match(other));

        }
        bool operator == (const MCell<Mesh, MeshFunction> &other) const { return (this->m_cellid == other.m_cellid); }
        bool operator > (const MCell<Mesh, MeshFunction> &other) const { return !(*this < other) && !(*this == other); }

        bool operator < (const MCell<Mesh, MeshFunction> &other) const {

            if (*this == other)        return false;
            /*
            // vertex always comes first
            if( this->m_dim == 0 )      return true;
            if( other.m_dim == 0 )      return false;

            // sort on boundary
            if(this->m_bdy != other.m_bdy) {
            return this->m_bdy < other.m_bdy;
            }
            */
            // normal Robin's sorting
            int n = std::min(m_n0faces, other.m_n0faces);
            for (int i = 0; i < n; i++) {

                // can't decide based on the same cell!
                if (m_0faces[i] == other.m_0faces[i])
                    continue;

                if (Cell_comparator<MeshFunction>::is_lower(m_0faces[i], other.m_0faces[i]))     return true;
                if (Cell_comparator<MeshFunction>::is_greater(m_0faces[i], other.m_0faces[i]))   return false;
            }

            // cell is ranked after its faces
            if (m_n0faces != other.m_n0faces)
                return (m_n0faces < other.m_n0faces);

            printf(" GOrdering.operator < (%d, %d) failed!\n", this->m_cellid, other.m_cellid);
            this->print("this");
            other.print("other");
            exit(0);
        }

        void print(const char* tag) const {

            int type = 2;

            printf(" %s(%d:%d %6d) = ", tag, this->m_dim, this->m_bdy, this->m_cellid);

            if (type == 1) {
                //for(int i = 0; i < m_0faces.size (); i++) {
                for (int i = 0; i < m_n0faces; i++) {
                    //printf("%.4f ", 100.0 * Vertex_cellID_comparator::value (m_0faces[i]) );
                    printf("%e ", Cell_comparator<MeshFunction>::value(m_0faces[i]));
                }
                printf("}\n");
            }
            else if (type == 2) {
                printf("[ ");
                /*for(int i = 0; i < m_n0faces; i++) {
                printf("%d ", m_0faces[i]);
                }*/
                for (int i = 0; i < m_nfaces_in_lstar; i++) {
                    printf("%6d ", m_faces_in_lstar[i]->cellid());
                }
                printf("]\n");

                /*printf("] --  [ ");
                for(int i = 0; i < m_nfaces_in_lstar; i++) {
                printf("<%d> ", m_faces_in_lstar[i]);
                }*/
                /*
                for(int i = 0; i < m_0faces.size (); i++) {
                //printf("%f ", 1000000.0 * Vertex_cellID_comparator::value (m_0faces[i]) );
                printf("%e ", Vertex_cellID_comparator::value (m_0faces[i]) );
                }
                printf("}\n");*/
            }
            else if (type == 3) {
                for (int i = 0; i < m_nfaces_in_lstar; i++) {
                    printf("%d ", m_faces_in_lstar[i]);
                }
                printf("] -- { ");
                for (int i = 0; i < m_n0faces; i++) {
                    printf("%f ", Cell_comparator<MeshFunction>::value(m_0faces[i]));
                    //printf("%e ", Vertex_cellID_comparator::value (m_0faces[i]) );
                }
            }
            else {

                printf(" %s(%d, %d) = [\n ", tag, this->m_cellid, this->m_dim);//Vertex_cellID_comparator::ms_mesh->dimension (this->m_cellid));
                for (int i = 0; i < m_n0faces; i++) {
                    printf("\t cell %d (vertex %d) = %e\n", m_0faces[i], Cell_comparator<MeshFunction>::vid(m_0faces[i]), Cell_comparator<MeshFunction>::value(m_0faces[i]));
                }
                printf(" ]\n");
            }

        }
    };

    // -------------------------------------------------------------------------------
    template <class Mesh, class MeshFunction>
    class LowerStar {

    private:
        INDEX_TYPE m_vertex_cellId;
        Vec3l m_vertex_coords;
        const MeshFunction *m_topofunc;
        const Mesh *m_tgrid;

        DIM_TYPE m_ncells;
        MCell<Mesh, MeshFunction> m_cells[MAX_CELLS_PER_STAR];

    public:
        LowerStar(INDEX_TYPE vertex_cellId, const MeshFunction *topofunc,
            const Mesh *tgrid,
            const DenseLabeling<DIM_TYPE> &max_Vlabeling) :
            m_vertex_cellId(vertex_cellId), m_topofunc(topofunc), m_tgrid(tgrid) {

            //TopologicalRegularGrid::CellVerticesIterator cviter(m_tgrid);
            m_tgrid->cellid2Coords(m_vertex_cellId, m_vertex_coords);

            m_ncells = 0;
            typename Mesh::AdjacentCellsIterator aciter(m_tgrid);
            for (aciter.begin(m_vertex_coords); aciter.valid(); aciter.advance()) {

                INDEX_TYPE adjCell = aciter.value();

                //printf(" checking if %d (%d) is in lower star of %d\n", adjCell, m_tgrid->dimension (adjCell), vertex_cellId);

                // get the max vertex for this adjacent cell
                INDEX_TYPE maxV = m_tgrid->UncompressByteToVertexOffset(adjCell, max_Vlabeling.GetLabel(adjCell));

                // if I am not the max 0face of this cell, this cell is not in my lower star
                if (maxV != m_vertex_cellId)
                    continue;

                m_cells[m_ncells++].initialize(adjCell, m_topofunc, m_tgrid);
            }
            std::sort(m_cells, m_cells + m_ncells);

            for (int i = 0; i < m_ncells; i++) {
                m_cells[i].populate_faces(*this);
            }
        }

        const Mesh* mesh() const { return m_tgrid; }

        size_t size() const { return m_ncells; }
        const MCell<Mesh, MeshFunction>* cell(unsigned int k) const { return &m_cells[k]; }

        bool contains(INDEX_TYPE cellid) const { return (find(cellid) != 0); }

        // have to do linear search since array is sorted on G-ordering
        const MCell<Mesh, MeshFunction>* find(INDEX_TYPE cellid) const {

            for (int i = 0; i < m_ncells; i++) {
                const MCell<Mesh, MeshFunction>* g = cell(i);
                if (g->cellid() == cellid)
                    return g;
            }
            return 0;
        }


        // -------------------------------------------------------------------
        // find the steepest edge in the star
        const MCell<Mesh, MeshFunction>* first_edge() const {

            for (DIM_TYPE k = 1; k < m_ncells; k++) {

                const MCell<Mesh, MeshFunction> *ocell = cell(k);
                if (ocell->dim() == 1) {
                    return ocell;
                }
            }
#ifdef DEBUG_MODE
            printf(" could not find an edge in lower star\n");
            this->print(1);
#endif
            return 0;
        }


        // find the steepest edge that can be paired with the vertex!
        const MCell<Mesh, MeshFunction>* find_vertex_pairing() const {

            for (DIM_TYPE k = 1; k < m_ncells; k++) {

                const MCell<Mesh, MeshFunction> *ocell = cell(k);
                if (ocell->dim() == 1) {
                    if (cell(0)->boundary_match(ocell))
                        return ocell;
                }
            }
            return 0;
        }
        // -------------------------------------------------------------------

        int num_unassigned(DiscreteGradientLabeling const*const labeling) const {

            int nunassigned = 0;
            for (int i = 0; i < m_ncells; i++) {

                if (1 != labeling->getAssigned(cell(i)->cellid())) {
                    nunassigned++;
                }
            }
            return nunassigned;
        }


        // -------------------------------------------------------------------
        void print(bool detailed) const {

            printf(" lower star of %d contains %d elements!\n", m_vertex_cellId, m_ncells);

            if (!detailed)
                return;

            for (int i = 0; i < m_ncells; i++) {
                cell(i)->print("C");
            }
        }

        void print_pairing(DiscreteGradientLabeling const*const labeling) const {

            int npaired = 0;

            // print assigned
            printf(" Assigned cells: \n");
            for (int i = 0; i < m_ncells; i++) {

                INDEX_TYPE cellid = cell(i)->cellid();
                if (1 == labeling->getAssigned(cellid)) {

                    if (labeling->getCritical(cellid)) {
                        printf("\t[%d <--> %d] critical\n", cellid, cellid);
                    }
                    else {
                        printf("\t[%d <--> %d]\n", cellid, labeling->getPair(cellid));
                    }
                    npaired++;
                }
            }

            if (npaired == m_ncells)
                return;

            printf(" Unssigned cells (%d): \n", (m_ncells - npaired));
            for (int i = 0; i < m_ncells; i++) {

                INDEX_TYPE cellid = cell(i)->cellid();
                if (1 != labeling->getAssigned(cellid)) {
                    printf("\t[%d --> ?]\n", cellid);
                }
            }
        }
    };

    // -------------------------------------------------------------------------------

    template<class Mesh, class MeshFunction>
    void MCell<Mesh, MeshFunction>::populate_faces(const LowerStar<Mesh, MeshFunction> &lstar) {

        // also need faces for quick queries later
        m_nfaces_in_lstar = 0;
        typename Mesh::FacetsIterator fit(lstar.mesh());
        for (fit.begin(m_cell_coords); fit.valid(); fit.advance()) {

            const MCell<Mesh, MeshFunction> *face = lstar.find(fit.value());
            if (face != 0) {
                m_faces_in_lstar[m_nfaces_in_lstar++] = face;
            }
        }
        //std::sort(m_faces_in_lstar, m_faces_in_lstar+m_nfaces_in_lstar);
        for (int i = m_nfaces_in_lstar; i < MAX_FACES_PER_CELL; i++) {
            m_faces_in_lstar[i++] = 0;
        }
    }

    // ===============================================================================================
    // ===============================================================================================

    template <class Mesh, class MeshFunction>
    class RobinsLabelingAlgorithm {

        MeshFunction *const m_topofunc;
        DiscreteGradientLabeling *const m_labeling;
        Mesh *m_tgrid;

        DenseLabeling<DIM_TYPE> *maxV_labeling;



    public:
        RobinsLabelingAlgorithm(MeshFunction *const topofunc, Mesh *const tgrid,
                                DiscreteGradientLabeling *const labeling) :

            m_topofunc(topofunc), m_tgrid(tgrid), m_labeling(labeling), maxV_labeling(0)
        {
            // need to handle the static pointers better!

            Cell_comparator<MeshFunction>::ms_topofunc = this->m_topofunc;
        }

        // -----------------------------------------------------------
        DIM_TYPE idxOf_maxVertex(INDEX_TYPE cellid) const {

            typename Mesh::CellVerticesIterator cviter(m_tgrid);
            cviter.begin(cellid);

            INDEX_TYPE maxV = cviter.value();
            DIM_TYPE test_dir = m_tgrid->CompressVertexOffsetToByte(cellid, maxV, cviter);

            for (cviter.advance(); cviter.valid(); cviter.advance()) {

                if (Cell_comparator<MeshFunction>::is_greater(cviter.value(), maxV)) {
                    maxV = cviter.value();
                    test_dir = m_tgrid->CompressVertexOffsetToByte(cellid, maxV, cviter);
                }
            }
            return test_dir;
        }

        void create_maxVertex_labeling() {

            ThreadedTimer timer(1);
            timer.StartGlobal();

            printf(" -- Creating maxV_labeling...");
            fflush(stdout);

            int num_cells = m_tgrid->numCells();
            maxV_labeling = new DenseLabeling<DIM_TYPE>(num_cells);

#pragma omp parallel
            {
                int num_threads = omp_get_num_threads();
                int thread_num = omp_get_thread_num();

                std::vector<INDEX_TYPE> partition;
                ArrayIndexPartitioner::EvenChunkSplit(m_tgrid->numCells(), num_threads, partition);

                for (INDEX_TYPE cellid = partition[thread_num]; cellid < partition[thread_num + 1]; cellid++) {
                    maxV_labeling->SetLabel(cellid, idxOf_maxVertex(cellid));
                }
            }
            timer.EndGlobal();
            printf(" Done! ");
            timer.PrintAll();
        }

        // -----------------------------------------------------------

        // find the single unpaired face of alpha
        const MCell<Mesh, MeshFunction>* single_unpaired_face(const MCell<Mesh, MeshFunction> *alpha) {

            DIM_TYPE nfaces = alpha->nfaces();
            for (DIM_TYPE i = 0; i < nfaces; i++) {

                const MCell<Mesh, MeshFunction> *face = alpha->face(i);
                if (!is_assigned(face->cellid())) {


                    if (alpha->boundary_match(face))

                        return face;
                }
            }
            return 0;
        }

        // number of unpaired faces restricted to a face's boundary
        DIM_TYPE num_unpaired_faces(const MCell<Mesh, MeshFunction> *alpha, const MCell<Mesh, MeshFunction> *restriction) {

            DIM_TYPE cnt = 0;

            DIM_TYPE nfaces = alpha->nfaces();
            for (DIM_TYPE i = 0; i < nfaces; i++) {

                const MCell<Mesh, MeshFunction> *face = alpha->face(i);
                if (!is_assigned(face->cellid())) {


                    if (alpha->boundary_match(restriction))

                        cnt++;
                }
            }
            return cnt;
        }

        // -----------------------------------------------------------
        // these functions interact with labeling

        bool is_critical(INDEX_TYPE cellid) const { return m_labeling->getCritical(cellid); }
        bool is_assigned(INDEX_TYPE cellid) const { return (1 == m_labeling->getAssigned(cellid)); }

        void set_critical(INDEX_TYPE cellid, int dim, bool verbose = false) {

#ifdef DEBUG_MODE
            if (verbose)
                printf("\n == set_critical (%d, %d)\n", cellid, dim);

            if (is_critical(cellid)) {
                printf(" : set_critical (%d, %d) : already critical!\n", cellid, dim);
                exit(1);
            }
            if (is_assigned(cellid)) {
                printf(" : set_critical (%d, %d) : already assigned!\n", cellid, dim);
                exit(1);
            }
#endif
            m_labeling->setCritical(cellid, true);
            m_labeling->setAssigned(cellid, 1);
        }

        void set_pair(INDEX_TYPE cellid, INDEX_TYPE pairid, bool verbose = false) {

#ifdef DEBUG_MODE
            if (verbose)
                printf("\n == set_pair (%d, %d)\n", cellid, pairid);


            if (!boundary_match(cellid, pairid)) {
                printf("\n\n ===== attempting to pair across boundary... %d (%d) and %d (%d) \n\n =====\n", cellid, m_tgrid->boundaryValue(cellid),
                    pairid, m_tgrid->boundaryValue(pairid));
                exit(1);
            }

            if (is_assigned(cellid) || is_assigned(pairid)) {
                printf(" : set_pair (%d, %d) : already assigned! (%d %d)\n", is_assigned(cellid), is_assigned(pairid));
                exit(1);
            }

            int dim1 = this->m_tgrid->dimension(cellid);
            int dim2 = this->m_tgrid->dimension(pairid);
            if (dim2 != dim1 + 1) {
                printf(" : set_pair (%d, %d) : wrong dimensions! (%d %d)\n", cellid, pairid, dim1, dim2);
                exit(1);
            }
#endif
            m_labeling->setPair(cellid, pairid);
            m_labeling->setPair(pairid, cellid);
            m_labeling->setAssigned(cellid, 1);
            m_labeling->setAssigned(pairid, 1);
        }

        bool boundary_match(INDEX_TYPE a, INDEX_TYPE b) {

            return (m_tgrid->boundaryValue(a) == m_tgrid->boundaryValue(b));

            printf(" a = %d, b = %d\n", this->m_tgrid->boundaryValue(a),
                this->m_tgrid->boundaryValue(b)
            );
        }

        // --------------------------------------------------------------

        void process_lowerStar(INDEX_TYPE vid) {


            bool debug = vid == -1;//325;//40215;//806;

#ifdef DEBUG_MODE
            if (0)
                printf(" process_lowerStar(%d)\n", vid);
#endif
            PriorityQueue<MCell<Mesh, MeshFunction> > PQ0, PQ1;

            INDEX_TYPE vert_cellId = m_tgrid->CellIDFromVertexNumber(vid);

            // compute the lower star
            LowerStar<Mesh, MeshFunction> lstar(vert_cellId, this->m_topofunc, this->m_tgrid, *maxV_labeling);

#ifdef DEBUG_MODE
            if (debug) {
                printf("\n ----------------------------\n");
                lstar.print(true);
                printf(" ----------------------------\n");
            }
#endif
            // critical vertex
            if (lstar.size() == 1) {
                set_critical(vert_cellId, 0);
#ifdef DEBUG_MODE
                if (debug) {
                    lstar.cell(0)->print("->critical v");
                    lstar.print(true);
                }
#endif
                return;
            }

            // ----------------------------------------------------------------------------
            // pass 1: go over all cells in lower star
            // (1) pair the current vertex to the steepest edge
            // (2) fill up PQ0 and PQ1

#ifdef DEBUG_MODE
            if (debug)
                printf("\n * initial assignment\n");
#endif
            const MCell<Mesh, MeshFunction> *del = lstar.find_vertex_pairing();

            // vertex pairing was found
            if (del != 0) {
#ifdef DEBUG_MODE
                if (debug) {
                    lstar.cell(0)->print("v");
                    del->print("del");
                }
#endif
                set_pair(vert_cellId, del->cellid(), debug);
            }

            // vertex pairing could not be found
            // start with the steepest edge
            else {
                set_critical(vert_cellId, 0);
                del = lstar.first_edge();
#ifdef DEBUG_MODE
                if (debug) {
                    lstar.cell(0)->print("->critical v");
                    del->print("steepest edge");
                    lstar.print(true);
                }
#endif
            }

            // iterate from the 2nd element
            for (int k = 1; k < lstar.size(); k++) {

                const MCell<Mesh, MeshFunction> *alpha = lstar.cell(k);
                if (is_assigned(alpha->cellid()))
                    continue;

                // number of unpaired faces restricted to del
                DIM_TYPE nalpha = num_unpaired_faces(alpha, del);

                // add all faces with no unpaired face
                if (nalpha == 0
                    && alpha->boundary_match(del)
                    ) {
                    PQ0.push_unique(alpha);
#ifdef DEBUG_MODE
                    if (debug)
                        alpha->print("pushed to PQ0");
#endif
                }

                // add all cofaces of del which have a single unpaired face left
                else if (nalpha == 1 && alpha->is_candidate_for_pairing(del)) {

                    PQ1.push_unique(alpha);
#ifdef DEBUG_MODE
                    if (debug)
                        alpha->print("pushed to PQ1");
#endif
                }
            }


            // ----------------------------------------------------------------------------
            // pass 2: assign all cells in the lower star

            while (!PQ0.empty() || !PQ1.empty()
                || lstar.num_unassigned(m_labeling) > 0
                ) {

#ifdef DEBUG_MODE
                if (debug) {
                    printf("\n * while there remains an unpaired cell in lstar! [%d %d %d]\n",
                        PQ0.size(), PQ1.size(), lstar.num_unassigned(m_labeling));
                    PQ0.print("\nPQ0");
                    PQ1.print("\nPQ1");
                }

                if (!PQ0.validate()) { printf(" PQ0 is invalid!"); PQ0.print("PQ0");   exit(1); }
                if (!PQ1.validate()) { printf(" PQ1 is invalid!"); PQ0.print("PQ1");   exit(1); }
#endif
                // homotopic expansion
                while (!PQ1.empty()) {

                    const MCell<Mesh, MeshFunction> *alpha = PQ1.pop_front();

#ifdef DEBUG_MODE
                    if (debug) {
                        printf("\n * homotopic expansion...\n");
                        alpha->print("popped from PQ1");
                    }
#endif
                    // no restriction needed!
                    if (num_unpaired_faces(alpha, alpha) == 0) {
                        PQ0.push_unique(alpha);
#ifdef DEBUG_MODE
                        if (debug)
                            alpha->print("pushed to PQ0");
#endif
                    }
                    else {

                        const MCell<Mesh, MeshFunction> *pairAlpha = single_unpaired_face(alpha);

                        if (pairAlpha == 0) {
                            set_critical(alpha->cellid(), alpha->dim());
#ifdef DEBUG_MODE
                            if (debug)
                                alpha->print("critical alpha");
#endif
                        }
                        else {
                            set_pair(pairAlpha->cellid(), alpha->cellid(), debug);
                            PQ0.remove(pairAlpha);
#ifdef DEBUG_MODE
                            if (debug)
                                pairAlpha->print("popped from PQ0");
#endif
                        }


                        // iterate from the 2nd element
                        for (int k = 1; k < lstar.size(); k++) {

                            const MCell<Mesh, MeshFunction> *beta = lstar.cell(k);
                            if (is_assigned(beta->cellid()))
                                continue;

                            if ((num_unpaired_faces(beta, alpha) == 1) &&
                                (beta->is_candidate_for_pairing(alpha) || beta->is_candidate_for_pairing(pairAlpha))) {

                                PQ1.push_unique(beta);
#ifdef DEBUG_MODE
                                if (debug)
                                    beta->print("pushed to PQ1");
#endif
                            }
                        }
                        //}
                    }
                }

                // pick the smallest unpaired edge to initiate another homotopic expansion
                if (!PQ0.empty()) {

#ifdef DEBUG_MODE
                    if (debug)
                        printf("\n * starting with a new edge...\n");
#endif
                    const MCell<Mesh, MeshFunction> *gamma = PQ0.pop_front();
                    if (is_assigned(gamma->cellid())) {
                        continue;
                    }
#ifdef DEBUG_MODE
                    if (debug)
                        gamma->print("popped from PQ0");
#endif
                    set_critical(gamma->cellid(), gamma->dim(), debug);

                    // iterate from the 2nd element
                    for (int k = 1; k < lstar.size(); k++) {

                        const MCell<Mesh, MeshFunction> *alpha = lstar.cell(k);
                        if (is_assigned(alpha->cellid()))
                            continue;

                        if ((num_unpaired_faces(alpha, gamma) == 1) && alpha->is_candidate_for_pairing(gamma)) {

                            PQ1.push_unique(alpha);
#ifdef DEBUG_MODE
                            if (debug)
                                alpha->print("pushed to PQ1");
#endif
                        }
                    }
                }



                // if PQ0 and PQ1 are empty, normal algorithm will exit.
                // but when boundary constraints are applied, we need to reset
                // the algorithm with a face on a different boundary (restriction class)
                // we pick the first face in order of G, and look for its pairable faces only
                if (lstar.num_unassigned(m_labeling) != 0 && PQ0.empty() && PQ1.empty()) {

                    // which boundary type to pick next
                    int btype = -1;
#ifdef DEBUG_MODE
                    if (debug)
                        printf("\n * another seeding\n");
#endif
                    for (int k = 1; k < lstar.size(); k++) {

                        const MCell<Mesh, MeshFunction> *alpha = lstar.cell(k);
                        if (is_assigned(alpha->cellid()))
                            continue;

                        if ((num_unpaired_faces(alpha, alpha) == 0)) {

                            if (btype == -1) { btype = alpha->boundary(); }
                            if (btype == alpha->boundary()) { PQ0.push_unique(alpha); }
                        }

                        else if ((num_unpaired_faces(alpha, alpha) == 1)) {

                            if (btype == -1) { btype = alpha->boundary(); }
                            if (btype == alpha->boundary()) { PQ1.push_unique(alpha); }
                        }
                    }
                }

            }

            // ----------------------------------------------------------------------------
            // sanity check
#ifdef DEBUG_MODE
            if (lstar.num_unassigned(m_labeling) != 0) {
                printf(" processing failed... for cell %d\n", vert_cellId);
                lstar.print_pairing(m_labeling);
                exit(1);
            }
#endif

#ifdef DEBUG_MODE
            if (debug) {

                lstar.print(true);
                lstar.print_pairing(m_labeling);
                exit(1);
            }
#endif
        }


    void compute_output() {

        create_maxVertex_labeling();

        //maxV_labeling->OutputToIntFile("maxlabeling.raw");
        //setlocale(LC_NUMERIC, "");

        printf(" -- Creating Robin's Discrete Gradient...");
        fflush(stdout);

        ThreadedTimer timer(1);
        timer.StartGlobal();

        size_t num_verts = m_tgrid->numCells(0);

        #pragma omp parallel
        {
            int thread_num = omp_get_thread_num();
            int num_threads = omp_get_num_threads();

            std::vector<INDEX_TYPE> partition;
            ArrayIndexPartitioner::EvenChunkSplit(num_verts, num_threads, partition);

            for (INDEX_TYPE vid = partition[thread_num]; vid < partition[thread_num + 1]; vid++) {
                process_lowerStar(vid);
            }
        }
        timer.EndGlobal();
        printf(" Done! ");
        timer.PrintAll();
    }

    void summarize() {

        size_t ncells = m_tgrid->numCells();
        size_t total[4] = { 0,0,0,0 };
        size_t critical[4] = { 0,0,0,0 };
        size_t assigned[4] = { 0,0,0,0 };
        size_t unassigned[4] = { 0,0,0,0 };

        for (size_t i = 0; i < ncells; i++) {
            int dim = this->m_tgrid->dimension(i);
            total[dim]++;

            if (this->m_labeling->getCritical(i)) { critical[dim]++;    }
            if (this->is_assigned(i)) {             assigned[dim]++;    }
            else {                                  unassigned[dim]++;  }
        }

        printf("   -- found %d critical cells (%d %d %d %d)\n",
                    (critical[0] + critical[1] + critical[2] + critical[3]), critical[0], critical[1], critical[2], critical[3]);

        if (unassigned[0] > 0 || unassigned[1] > 0 || unassigned[2] > 0 || unassigned[3] > 0) {
           std::cerr << " ERROR: Discrete gradient algorithm could not assign all cells: (" <<
                        unassigned[0] << ", " << unassigned[1] << ", " << unassigned[2] << ", " << unassigned[3] << ")\n";
        }
        /*
        printf("\n total cells %d\n", ncells);
        for(int k = 0; k < 4; k++)
            printf(" dim %d : total %d, assigned %d, critical %d, unassigned %d\n", k, total[k], assigned[k], critical[k], unassigned[k]);
        */
    }

    void write(std::string filename, const std::vector<INDEX_TYPE> &cp) {
        std::ofstream out(filename.c_str());
        for (size_t i = 0; i < cp.size(); i++) {
            out << cp[i] << std::endl;
        }
        out.close();
    }
};
}   // end of namespace
#endif
