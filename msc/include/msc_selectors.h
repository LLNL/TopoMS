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

#ifndef MSC_SELECTORS_H
#define MSC_SELECTORS_H

#include <cassert>
#include "morse_smale_complex_basic.h"
#include "topological_regular_grid.h"
#include "msc_iterators.h"

namespace MSC {

// -----------------------------------------------------------------------------
// Base case for selectors
// -----------------------------------------------------------------------------
class MSCSelector {

protected:

    vector<MSCSelector*> parents;
    vector<MSCSelector*> children;

    bool output_valid;

    // when this is called parents have been computed
    virtual void SelectorAction() {}

    size_t add_child(MSCSelector* s) {
        children.push_back(s);
        return children.size();
    }

public:

    set<INT_TYPE> output;

    MSCSelector() : output_valid(false) {}

    virtual void invalidate() {
        if (!output_valid) return; // to prevent infinite looping by accident
        output_valid = false;
        for (auto it = children.begin(); it != children.end(); it++) (*it)->invalidate();
    }

    size_t add_parent(MSCSelector* s) {
        s->add_child(this);
        parents.push_back(s);
        invalidate();
        return parents.size();
    }

    void compute_output() {

        if (output_valid) return;
        output.clear();
        for (auto it = parents.begin(); it != parents.end(); it++) {
            (*it)->compute_output();
        }
        SelectorAction();

        output_valid = true;
    }
};

// -----------------------------------------------------------------------------
// Living Nodes selector
// i.e., select nodes still living under current persistence
// -----------------------------------------------------------------------------
template<class MSCType>
class MSCSelectorLivingNodes : public MSCSelector {

protected:
    MSCType* mMsc;
    virtual void SelectorAction() {

        // Harsh added this assertion, since it is assumed
        assert(this->parents.empty());
        for (size_t i = 0; i < mMsc->numNodes(); i++) {
            if (mMsc->isNodeAlive(i)) this->output.insert(i);
        }
    }
public:
    MSCSelectorLivingNodes(MSCType* msc) : mMsc(msc) {}
};

// -----------------------------------------------------------------------------
// NodeIndex selector
// i.e., select nodes of a certain index
// -----------------------------------------------------------------------------
template<class MSCType>
class MSCSelectorNodeIndex : public MSCSelector {

protected:
    MSCType* mMsc;
    DIM_TYPE mIndex;

    virtual void SelectorAction() {
        for (auto pit = parents.begin(); pit != parents.end(); pit++) {
        for (auto it = (*pit)->output.begin(); it != (*pit)->output.end(); it++) {
            if (mMsc->getNode(*it).dim == mIndex) this->output.insert(*it);
        }
        }
    }
public:
    MSCSelectorNodeIndex(MSCType* msc, DIM_TYPE index) : mMsc(msc), mIndex(index) {}
};

// -----------------------------------------------------------------------------
// NearestNode selector
// i.e., select the nearest node to a given position
// -----------------------------------------------------------------------------
template<class MSCType>
class MSCSelectorNearestNode : public MSCSelector {

protected:
    MSCType* mMsc;
    MSC::Vec3d mPos;

    virtual void SelectorAction() {
        std::multimap<FLOATTYPE, INDEX_TYPE> distmap;
        for (auto pit = parents.begin(); pit != parents.end(); pit++) {
        for (auto it = (*pit)->output.begin(); it != (*pit)->output.end(); it++) {

            MSC::Vec3l lcoord;
            mMsc->GetMesh()->cellid2Coords(mMsc->getNode(*it).cellindex, lcoord);

            MSC::Vec3d dcoord(lcoord);
            distmap.insert(std::make_pair((mPos-dcoord).MagSq(), *it));
        }
        }
        this->output.clear();
        this->output.insert(distmap.begin()->second);
    }
public:
    // this expects the position in grid coordinates
    MSCSelectorNearestNode(MSCType* msc, const MSC::Vec3d &gPos) : mMsc(msc), mPos(gPos) {
        // convert from (regular) grid to topological grid
        mPos *= 2.0;
    }
};

// -----------------------------------------------------------------------------
// IncidentArcs selector
// i.e., select the arcs connected to given nodes
// -----------------------------------------------------------------------------
template<class MSCType>
class MSCSelectorIncidentArcs : public MSCSelector {

protected:
    MSCType* mMsc;
    int mIdx;

    virtual void SelectorAction() {
        this->output.clear();
        MSC::MSCIteratorSurroundingArcs<MSCType> sit(mMsc);
        for (auto pit = parents.begin(); pit != parents.end(); pit++) {                 // for each parent selector
        for (auto it = (*pit)->output.begin(); it != (*pit)->output.end(); it++) {      // for each node
        for (sit.begin(*it); sit.valid(); sit.advance()) {                              // each arc

            INDEX_TYPE arcId = sit.value();
            if (mIdx == -1 || mIdx == mMsc->getArc(arcId).dim)
                this->output.insert(arcId);
        }
        }
        }
    }
public:
    // this expects the position in grid coordinates
    MSCSelectorIncidentArcs(MSCType* msc, int idx=-1) : mMsc(msc), mIdx(idx) {}
};


template<class MSCType>
class MSCSelectorIncidentArcs2 : public MSCSelector {

protected:
    MSCType* mMsc;
    int mIdx;

    virtual void SelectorAction() {
        this->output.clear();

        std::set<INDEX_TYPE> inodes;
        for (auto pit = parents.begin(); pit != parents.end(); pit++) {                 // for each parent selector
        for (auto it = (*pit)->output.begin(); it != (*pit)->output.end(); it++) {      // for each node
            inodes.insert(*it);
        }
        }

        printf(" selected incident arcs to %d nodes!\n", inodes.size());
        size_t narcs = mMsc->numArcs();
        for(size_t i = 0; i < narcs; i++) {

            if (!mMsc->isArcAlive(i))                       continue;

            MSC::arc<FLOATTYPE> &a = mMsc->getArc(i);

            if (mIdx != -1 && mIdx != mMsc->getArc(i).dim)  continue;

            bool lowerMatch = inodes.count(a.lower) > 0;  // lower node of this arc exists in inodes
            bool upperMatch = inodes.count(a.upper) > 0;  // upper node of this arc exists in inodes
            if (!lowerMatch && !upperMatch)                 continue;

            this->output.insert(i);
        }
    }
public:
    // this expects the position in grid coordinates
    MSCSelectorIncidentArcs2(MSCType* msc, int idx=-1) : mMsc(msc), mIdx(idx) {}
};

// -----------------------------------------------------------------------------
// IncidentNodes selector
// i.e., select the end points of given arcs
// -----------------------------------------------------------------------------
template<class MSCType>
class MSCSelectorIncidentNodes : public MSCSelector {

protected:
    MSCType* mMsc;
    int mIdx;

    virtual void SelectorAction() {
        this->output.clear();
        for (auto pit = parents.begin(); pit != parents.end(); pit++) {                 // for each parent selector
        for (auto it = (*pit)->output.begin(); it != (*pit)->output.end(); it++) {      // for each arc

            arc<FLOATTYPE> &a = mMsc->getArc(*it);
            if(mIdx != 1)   this->output.insert(a.lower);
            if(mIdx != 0)   this->output.insert(a.upper);
        }
        }
    }
public:
    // this expects the position in grid coordinates
    MSCSelectorIncidentNodes(MSCType* msc, int idx=-1) : mMsc(msc), mIdx(idx) {}
};

// -----------------------------------------------------------------------------
// Representative1Saddle selector
// -----------------------------------------------------------------------------
template<class MSCType>
class MSCSelectorRepresentative1Saddle : public MSCSelector {

protected:
    MSCType* mMsc;
    DIM_TYPE mIndex;

    bool lessthan(INT_TYPE a, INT_TYPE b) {

        node<FLOATTYPE>& na = mMsc->getNode(a);
        node<FLOATTYPE>& nb = mMsc->getNode(b);
        if (na.value == nb.value) return a < b;
        return na.value < nb.value;
    }

    pair<INT_TYPE, INT_TYPE> getExtremumPair(INT_TYPE saddle) {

        MSCIteratorSurroundingLivingArcs<MSCType> sit(mMsc);
        int numext = 0;
        INT_TYPE extrema[2];

        for (sit.begin(saddle); sit.valid(); sit.advance()) {
            arc<FLOATTYPE>& a = mMsc->getArc(sit.value());
            if (a.upper == saddle) extrema[numext++] = a.lower;
        }
        if (numext == 1) {                  return pair<INT_TYPE, INT_TYPE>(extrema[0], extrema[0]);    }
        else if (extrema[0] < extrema[1]) { return pair<INT_TYPE, INT_TYPE>(extrema[0], extrema[1]);    }
        else {                              return pair<INT_TYPE, INT_TYPE>(extrema[1], extrema[0]);    }
    }


    virtual void SelectorAction() {
        map<pair<INT_TYPE, INT_TYPE>, INT_TYPE> saddlemap;

        for (auto pit = parents.begin(); pit != parents.end(); pit++) {
        for (auto it = (*pit)->output.begin(); it != (*pit)->output.end(); it++) {

            INT_TYPE saddle = *it;
            pair<INT_TYPE, INT_TYPE> p = getExtremumPair(saddle);

            if (saddlemap.count(p) == 0) {              saddlemap[p] = saddle;  }
            else {
                INT_TYPE othersaddle = saddlemap[p];
                if (lessthan(saddle, othersaddle)) {    saddlemap[p] = saddle;  }
            }
        }
        }

        for (auto it = saddlemap.begin(); it != saddlemap.end(); it++) {
            if ((*it).first.first != (*it).first.second)
                output.insert((*it).second);
            }
        }

public:
    MSCSelectorRepresentative1Saddle(MSCType* msc) : mMsc(msc) {}
};

}   // end of namespace
#endif
