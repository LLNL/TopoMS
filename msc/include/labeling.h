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

#ifndef VERTEX_LABELING_H
#define VERTEX_LABELING_H

#include "basic_types.h"

namespace MSC {
    template<typename LABEL_TYPE>
    class DenseLabeling {
    protected:
        LABEL_TYPE* m_labels;
        INDEX_TYPE m_num_labels;
    public:

        DenseLabeling(INDEX_TYPE num_labels) : m_num_labels(num_labels) {
            m_labels = new LABEL_TYPE[num_labels];
            //printf(" allocated space for %d values of size %d\n", num_labels, sizeof(LABEL_TYPE));
        }
        ~DenseLabeling() {
            delete[] m_labels;
        }

        void SetLabel(INDEX_TYPE id, LABEL_TYPE label) {
            m_labels[id] = label;
        }

        LABEL_TYPE GetLabel(INDEX_TYPE id) const {
            return m_labels[id];
        }

        INDEX_TYPE GetNumLabels() const {
            return m_num_labels;
        }

        LABEL_TYPE& operator[](const INDEX_TYPE id) { return m_labels[id]; }
        const LABEL_TYPE& operator[](const INDEX_TYPE id) const { return m_labels[id]; }

        void SetAll(LABEL_TYPE label){
#pragma omp parallel for schedule(static)
            for (int i = 0; i < m_num_labels; i++) {
                m_labels[i] = label;
            }
        }

        void ReadFromFile(const char* filename) {
            FILE* fout = fopen(filename, "rb");
            fread(m_labels, sizeof(LABEL_TYPE), m_num_labels, fout);
            fclose(fout);
        }

        void OutputToFile(const char* filename) const {
            FILE* fout = fopen(filename, "wb");
            //printf("Sizeof label type: %d\n", sizeof(LABEL_TYPE));
            fwrite(m_labels, sizeof(LABEL_TYPE), m_num_labels, fout);
            fclose(fout);
        }
        void OutputToIntFile(const char* filename) const {
            FILE* fout = fopen(filename, "wb");
            //printf("Sizeof int type: %d\n", sizeof(int));
            for (INDEX_TYPE i = 0; i < m_num_labels; i++) {
                int tval = (int)m_labels[i];
                fwrite(&tval, sizeof(int), 1, fout);
            }
            fclose(fout);
        }
        void OutputToFloatFile(const char* filename) const {
            FILE* fout = fopen(filename, "wb");
            //printf("Sizeof float type: %d\n", sizeof(float));
            for (INDEX_TYPE i = 0; i < m_num_labels; i++) {
                float tval = (float)m_labels[i];
                fwrite(&tval, sizeof(float), 1, fout);
            }
            fclose(fout);
        }

        LABEL_TYPE* LabelArray() { return m_labels; }
    };



}

#endif
