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

#ifndef MSC_ITERATORS_H
#define MSC_ITERATORS_H

#include <cassert>
#include "morse_smale_complex_basic.h"
#include "topological_regular_grid.h"

namespace MSC {

// -----------------------------------------------------------------------------
// node and arc iterators
// -----------------------------------------------------------------------------

template<class MSCType>
class MSCIteratorNodes {
protected:
    MSCType* mMSC;
    INT_TYPE currid;

public:
    MSCIteratorNodes(MSCType* msc) : mMSC(msc){}
    void begin() {      currid = 0;     }
    void advance() {    currid++;       }
    INT_TYPE value() {  return currid;  }
    bool valid() {      return currid < mMSC->nodes.size(); }
};

// -----------------------------------------------------------------------------
template<class MSCType>
class MSCIteratorLivingNodes : public MSCIteratorNodes<MSCType> {
protected:
    void advance_until_alive() {
        this->currid++;
        while (this->valid() && !this->mMSC->isNodeAlive(this->currid)) {
            this->currid++;
        }
    }
public:
    MSCIteratorLivingNodes(MSCType* msc) : MSCIteratorNodes<MSCType>(msc) {}
    void begin() {
        MSCIteratorNodes<MSCType>::begin();
        if (this->valid() && !this->mMSC->isNodeAlive(this->currid))
          advance_until_alive();
    }
    void advance() {
        advance_until_alive();
    }
};

// -----------------------------------------------------------------------------
template<class MSCType>
class MSCIteratorArcs {
protected:
    MSCType* mMSC;
    INT_TYPE currid;

public:
    MSCIteratorArcs(MSCType* msc) : mMSC(msc){}
    void begin() {      currid = 0;     }
    void advance() {    currid++;       }
    INT_TYPE value() {  return currid;  }
    bool valid() {      return currid < mMSC->arcs.size();  }
};

// -----------------------------------------------------------------------------
template<class MSCType>
class MSCIteratorLivingArcs : public MSCIteratorArcs<MSCType> {
protected:
    void advance_until_alive() {
        this->currid++;
        while (this->valid() && !this->mMSC->isArcAlive(this->currid)) {
            this->currid++;
        }
    }
public:
    MSCIteratorLivingArcs(MSCType* msc) : MSCIteratorArcs<MSCType>(msc) {}
    void begin() {
        MSCIteratorArcs<MSCType>::begin();
        if (this->valid() && !this->mMSC->isNodeAlive(this->currid))
          advance_until_alive();
    }
    void advance() {
        advance_until_alive();
    }
};

// -----------------------------------------------------------------------------
template<class MSCType>
class MSCIteratorSurroundingArcs {
protected:
    MSCType* mMSC;
    INT_TYPE mNID;
    INT_TYPE currarc;
    INT_TYPE next_arc(INT_TYPE arcid) {
        return mMSC->nextArc(mMSC->getArc(arcid), mNID);
    }

public:
    MSCIteratorSurroundingArcs(MSCType* msc) : mMSC(msc){}

    bool valid() {      return currarc != NULLID;   }
    INT_TYPE value() {  return currarc;             }

    void begin(INT_TYPE nid) {
        mNID = nid;
        node<FLOATTYPE>& n = mMSC->getNode(mNID);
        currarc = n.firstarc;
    }
    void advance() {
        if (currarc == NULLID) return;
        currarc  = next_arc(currarc);
    }
};

// -----------------------------------------------------------------------------
template<class MSCType>
class MSCIteratorSurroundingLivingArcs : public MSCIteratorSurroundingArcs<MSCType> {
protected:
    bool advance_until_alive() {
        this->currarc = this->next_arc(this->currarc);
        if (this->currarc == NULLID) return false;
        while (!this->mMSC->isArcAlive(this->currarc)) {
            this->currarc = next_arc(this->currarc);
            if (this->currarc == NULLID) return false;
        }
        return true;
}
public:
    MSCIteratorSurroundingLivingArcs(MSCType* msc) : MSCIteratorSurroundingArcs<MSCType>(msc) {}
    void begin(INT_TYPE nid) {
        MSCIteratorSurroundingArcs<MSCType>::begin(nid);
        if (!this->mMSC->isArcAlive(this->currarc)) advance_until_alive();
    }
    void advance() {
        if (this->currarc == NULLID) return;
        advance_until_alive();
    }
};

}   // end of namespace
#endif
