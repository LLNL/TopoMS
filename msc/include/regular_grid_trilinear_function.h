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

#ifndef REGULAR_GRID_TRILINEAR_FUNCTION
#define REGULAR_GRID_TRILINEAR_FUNCTION

#include "vectors.h"
#include "regular_grid.h"
#include <algorithm>


namespace MSC {
    // store image and gradient tied to a 3d grid.
    class RegularGridTrilinearFunction {
    protected:
        RegularGrid * m_grid;
        Vec3d* m_grad;

        FLOATTYPE* m_image;

        bool m_i_made_gradient;
        bool m_i_made_image;
    public:
        RegularGridTrilinearFunction(RegularGrid* grid, FLOATTYPE *image = 0) : m_grid(grid) {
            m_i_made_image = false;
            m_i_made_gradient = false;
            m_image = NULL;
            m_grad = NULL;
            // use the function if it is passed, otherwise simply allocate memory
            if(image != 0) {    m_image = image;                                }

            //m_grad = new Vec3d[m_grid->NumElements()];
        }
        ~RegularGridTrilinearFunction() {
            if (m_i_made_gradient) delete[] m_grad;
            if (m_i_made_image) delete[] m_image;
        }

        // return pointer to underlying mesh and function
        const RegularGrid* GetGrid() const {    return m_grid;  }
        const FLOATTYPE* GetImage() const {         return m_image; }

        // sample the image at integral location
        FLOATTYPE SampleImage(const Vec3l& p) const {
            return m_image[m_grid->Index3d(p)];
        }

        // sample the image at integral location
        FLOATTYPE SampleImage(const INDEX_TYPE id) const {
            return m_image[id];
        }

        // sample the gradient at integral location
        const Vec3d& SampleGrad(const Vec3l& p) const {
            return m_grad[m_grid->Index3d(p)];
        }


        FLOATTYPE TriLinInterpValue(const Vec3d& s) const {
            Vec3l n[8]; // for 8 vertices around s - some may be repeated based on boundary cond.
            m_grid->GatherSurrounding(s, n);
            Vec3d b = n[0];
            //s.print_vf();
            //b.print_vf();
            Vec3d factors = s - b;

            FLOATTYPE x0 = (1 - factors[0]) * SampleImage(n[0])  +  SampleImage(n[1]) * factors[0];
            FLOATTYPE x1 = (1 - factors[0]) * SampleImage(n[2]) + SampleImage(n[3]) * factors[0];
            FLOATTYPE x2 = (1 - factors[0]) * SampleImage(n[4]) + SampleImage(n[5]) * factors[0];
            FLOATTYPE x3 = (1 - factors[0]) * SampleImage(n[6]) + SampleImage(n[7]) * factors[0];

            FLOATTYPE y0 = (1 - factors[1]) *x0 + x1 * factors[1];
            FLOATTYPE y1 = (1 - factors[1]) *x2 + x3 * factors[1];

            return (1 - factors[2]) *y0 + y1 * factors[2];
        }

        // return trilinearly interpolated value
        Vec3d TriLinInterpGrad(const Vec3d& s) const {
            Vec3l n[8]; // for 8 vertices around s - some may be repeated based on boundary cond.
            m_grid->GatherSurrounding(s, n);
            Vec3d b = n[0];
            //s.print_vf();
            //b.print_vf();
            Vec3d factors = s - b;

            Vec3d x0 = Vec3d::Lerp(SampleGrad(n[0]), SampleGrad(n[1]), factors[0]);
            Vec3d x1 = Vec3d::Lerp(SampleGrad(n[2]), SampleGrad(n[3]), factors[0]);
            Vec3d x2 = Vec3d::Lerp(SampleGrad(n[4]), SampleGrad(n[5]), factors[0]);
            Vec3d x3 = Vec3d::Lerp(SampleGrad(n[6]), SampleGrad(n[7]), factors[0]);

            Vec3d y0 = Vec3d::Lerp(x0, x1, factors[1]);
            Vec3d y1 = Vec3d::Lerp(x2, x3, factors[1]);

            return Vec3d::Lerp(y0, y1, factors[2]);
        }

        void SetGradExplicit(INDEX_TYPE id, Vec3d vec) {
            this->m_grad[id] = vec;
        }

        // fill in vals with the 8 values of hte gradient around sample poitn
        void GetGradSurrounding(const Vec3d& s, Vec3d* vals) const {
            Vec3l n[8]; // for 8 vertices around s - some may be repeated based on boundary cond.
            m_grid->GatherSurrounding(s, n);
            for (int i = 0; i < 8; i++) vals[i] = SampleGrad(n[i]);
        }
        void GetGradSurrounding(const Vec3l& s, Vec3d* vals) const {
            Vec3l n[8]; // for 8 vertices around s - some may be repeated based on boundary cond.
            m_grid->GatherSurrounding(s, n);
            for (int i = 0; i < 8; i++) vals[i] = SampleGrad(n[i]);
        }

        // use with extreme care - no boundary checks, only do on really interior poitns
        void GetGradSurroundingNoBoundaryCheck(const Vec3d& s, Vec3d* vals) const {
            Vec3l n[8]; // for 8 vertices around s - some may be repeated based on boundary cond.
            m_grid->GatherSurroundingNoBoundaryCheck(s, n);
            for (int i = 0; i < 8; i++) vals[i] = SampleGrad(n[i]);
        }

        FLOATTYPE InterpolatedValue(const Vec3d& s) const {
            return TriLinInterpValue(s);
        }

        Vec3d InterpolatedGrad(const Vec3d& s) const {
            return TriLinInterpGrad(s);
        }

        // allow reuse of sampled gradient - the assumption that vals has the gradient arrows around s
        Vec3d TriLinInterpGrad(const Vec3d& s, const Vec3l& int_base, Vec3d* vals) const {

            //if (!(s.IntFloor() == int_base)) {
            //	printf("s=");  s.PrintFloat(); printf("d="); int_base.PrintFloat();
            //}
            //
            //Vec3d d = int_base.IntFloor();
            Vec3d factors = s - int_base;

            Vec3d x0 = Vec3d::Lerp(vals[0], vals[1], factors[0]);
            Vec3d x1 = Vec3d::Lerp(vals[2], vals[3], factors[0]);
            Vec3d x2 = Vec3d::Lerp(vals[4], vals[5], factors[0]);
            Vec3d x3 = Vec3d::Lerp(vals[6], vals[7], factors[0]);

            Vec3d y0 = Vec3d::Lerp(x0, x1, factors[1]);
            Vec3d y1 = Vec3d::Lerp(x2, x3, factors[1]);

            return Vec3d::Lerp(y0, y1, factors[2]);
        }
        void LoadImageFromFloatFile(const char* fname) {

            size_t image_size = m_grid->NumElements();

            // fill in image
            m_image = new FLOATTYPE[image_size];    m_i_made_image = true;


            FILE* fin = fopen(fname, "rb");
            for (INDEX_TYPE i = 0; i < image_size; i++) {

                float tval = 0;
                fread(&tval, sizeof(float), 1, fin);
                m_image[i] = tval;
            }

            fclose(fin);
            printf("min = %e, max = %e\n", *std::min_element(m_image, m_image+image_size), *std::max_element(m_image, m_image+image_size));
        }
        void LoadImageFromFile(const char* fname) {

            size_t image_size = m_grid->NumElements();

            // fill in image
            m_image = new FLOATTYPE[image_size]; m_i_made_image = true;

            FILE* fin = fopen(fname, "rb");
            fread(m_image, sizeof(FLOATTYPE), image_size, fin);
            fclose(fin);

            printf("min = %e, max = %e\n", *std::min_element(m_image, m_image+image_size), *std::max_element(m_image, m_image+image_size));
        }

        void ShallowCopyImage(FLOATTYPE *image) {

            m_image = image;

            INDEX_TYPE image_size = m_grid->NumElements();
            FLOATTYPE t_minval = *std::min_element(m_image, m_image+image_size);
            FLOATTYPE t_maxval = *std::max_element(m_image, m_image+image_size);
            printf(" function range = [%f, %f]\n", t_minval, t_maxval);
        }

        void DeepCopyImage(const FLOATTYPE *image) {

            m_image = new FLOATTYPE[m_grid->NumElements()];  m_i_made_image = true;

            INDEX_TYPE image_size = m_grid->NumElements();
            memcpy(m_image, image, image_size*sizeof(FLOATTYPE));

            FLOATTYPE t_minval = *std::min_element(m_image, m_image+image_size);
            FLOATTYPE t_maxval = *std::max_element(m_image, m_image+image_size);
            printf(" function range = [%f, %f]\n", t_minval, t_maxval);
        }

        static const FLOATTYPE kRKCoefficients[5][9];

        Vec3d GradientFromImage(const Vec3l& p, int rklevel) {
            int nume = rklevel * 2 + 1; // number of entries to average
            Vec3l negs[9]; // don't support more than 4th order - cmon. would be ridiculous

            double res_x = 0.0;
            m_grid->Gather1DNeighborhood(p, 0, rklevel, negs);
            for (int i = 0; i < nume; i++) {
                res_x += kRKCoefficients[rklevel][i] * SampleImage(negs[i]);
            }
            double res_y = 0.0;
            m_grid->Gather1DNeighborhood(p, 1, rklevel, negs);
            for (int i = 0; i < nume; i++) {
                res_y += kRKCoefficients[rklevel][i] * SampleImage(negs[i]);
            }
            double res_z = 0.0;
            m_grid->Gather1DNeighborhood(p, 2, rklevel, negs);
            for (int i = 0; i < nume; i++) {
                res_z += kRKCoefficients[rklevel][i] * SampleImage(negs[i]);
            }
            return Vec3d(res_x, res_y, res_z);
        }

        inline bool IsGreater(INDEX_TYPE a, INDEX_TYPE b) const {
            if (m_image[a] > m_image[b]) return true;
            if (m_image[b] > m_image[a]) return false;
            //if (a == b) printf("WHOA THERE NELLY\n");
            return a > b;
        }
        //Vec3d IStep(const Vec3d& p, const Vec3d& grad, const FLOATTYPE h) const {
        //	return m_grid->Inbounds(p + (grad * h));
        //}
        //Vec3d IStepNoBoundaryCheck(const Vec3d& p, const Vec3d& grad, const FLOATTYPE h) const {
        //	return p + (grad * h);
        //}


        // add in block structure




        void ComputeGradFromImage(int rklevel) {
            m_grad = new Vec3d[m_grid->NumElements()];
            m_i_made_gradient = true;
#pragma omp parallel for
            for (int i = 0; i < m_grid->XYZ()[0]; i++) {
                for (int j = 0; j < m_grid->XYZ()[1]; j++) {
                    for (int k = 0; k < m_grid->XYZ()[2]; k++) {
                        Vec3l p(i, j, k);
                        m_grad[m_grid->Index3d(p)] = GradientFromImage(p, rklevel);
                    }
                }
            }

        }



        void Negate() {
            if (m_grad != NULL) {
#pragma omp parallel for schedule(static)
                for (INDEX_TYPE i = 0; i < m_grid->NumElements(); i++) {
                    this->m_image[i] *= -1;
                    this->m_grad[i] *= -1.0;
                }

            }
            else {
#pragma omp parallel for schedule(static)
                for (INDEX_TYPE i = 0; i < m_grid->NumElements(); i++) {
                    this->m_image[i] *= -1;
                }
            }
        }

    };
#if 1

#endif

};

#endif
