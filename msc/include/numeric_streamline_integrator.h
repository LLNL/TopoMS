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

#ifndef GI_NUMERICAL_STREAMLINE_INTEGRATOR_H
#define GI_NUMERICAL_STREAMLINE_INTEGRATOR_H

#include <set>
#include <queue>
#include "basic_types.h"
#include "vectors.h"
#include "labeling.h"
#include "regular_grid.h"
#include "regular_grid_trilinear_function.h"
#include "adaptive_euler_advector.h"
//#include "adaptive_euler_advector_reverse.h"

namespace MSC {

    class StepperBase {

    protected:
        int m_direction;
        FLOATTYPE m_gradient_threshold;

        StepperBase(FLOATTYPE gradient_threshold, int direction = 1) :
            m_gradient_threshold(gradient_threshold), m_direction(direction) {}

    public:

        static Vec3d lerp3(const Vec3d &p, const Vec3l &base_p, const Vec3d *vals) {

            Vec3d weights = (p - base_p);

            Vec3d x0 = Vec3d::Lerp(vals[0], vals[1], weights[0]);
            Vec3d x1 = Vec3d::Lerp(vals[2], vals[3], weights[0]);
            Vec3d x2 = Vec3d::Lerp(vals[4], vals[5], weights[0]);
            Vec3d x3 = Vec3d::Lerp(vals[6], vals[7], weights[0]);

            Vec3d y0 = Vec3d::Lerp(x0, x1, weights[1]);
            Vec3d y1 = Vec3d::Lerp(x2, x3, weights[1]);

            return Vec3d::Lerp(y0, y1, weights[2]);
        }

        virtual bool step(Vec3d &p, const Vec3l &baseVtx, const Vec3d *gradients) { return false; }
    };

    class EulerStepper : public StepperBase {

        FLOATTYPE m_step_size;

    public:
        EulerStepper(FLOATTYPE step_size, FLOATTYPE gradient_threshold, int direction = 1) :
            m_step_size(step_size), StepperBase(gradient_threshold, direction) {}

        bool step(Vec3d &p, const Vec3l &baseVtx, const Vec3d *gradients)  {

            Vec3d grad = StepperBase::lerp3(p, baseVtx, gradients);
            double grad_mgn = grad.Mag();

            // check flatness! return if we are too flat
            if (grad_mgn < m_gradient_threshold) {
                //printf(" EulerStepper::step() stopping due to low gradient.. %f %f\n", grad_mgn, m_gradient_threshold);
                return false;
            }

            // else, take a step
            p = p + grad * (m_step_size*m_direction);
            return true;
        }

    };

    class AdaptiveEulerStepper : public StepperBase  {

        FLOATTYPE m_step_size;
        FLOATTYPE m_error_threshold;

    public:
        AdaptiveEulerStepper(FLOATTYPE step_size, FLOATTYPE error_threshold, FLOATTYPE gradient_threshold, int direction = 1) :
            m_step_size(step_size), m_error_threshold(error_threshold), StepperBase(gradient_threshold, direction) {}

        bool step(Vec3d &p, const Vec3l &baseVtx, const Vec3d *gradients)  {

            Vec3d grad = StepperBase::lerp3(p, baseVtx, gradients);
            double grad_mgn = grad.Mag();

            // check flatness! return if we are too flat
            if (grad_mgn < m_gradient_threshold) {
                //printf(" AdaptiveEulerStepper::step() stopping due to low gradient.. %f %f\n", grad_mgn, m_gradient_threshold);
                return false;
            }

            FLOATTYPE adaptive_step = m_step_size / grad_mgn;

            while (true) {

                //printf(" taking adaptive step... %f %f\n", grad_mgn, adaptive_step);

                // full-step displacement
                Vec3d step_full = grad * (adaptive_step*m_direction);

                // now do 2 steps with half sized steps to see if we are in error tolerance
                Vec3d step_half = grad * (0.5*adaptive_step*m_direction);
                Vec3d grad_half_step = StepperBase::lerp3((p + step_half), baseVtx, gradients);

                Vec3d step_2halves = grad_half_step * (0.5*adaptive_step*m_direction) + step_half;

                // now compare the points
                FLOATTYPE local_error = (step_2halves - step_full).MagSq();

                // if error is tolerable, accept the step!
                if (local_error < m_error_threshold) {
                    p = p + step_2halves;
                    m_step_size = min(m_step_size * 1.8, MAX_STEP_IN_GRID_CELLS / grad_mgn);
                    return true;
                }

                // the error is too high! reduce stepsize and try again
                adaptive_step *= 0.5;
                //printf(" ------ reduce step size!\n");
            }

            return false;
        }
    };

    class NumericIntegratorNew {

        std::vector<Vec3d> m_seeds;

        bool m_save_geometry;

        FLOATTYPE m_error_threshold;
        FLOATTYPE m_gradient_threshold;
        int m_max_steps;
        int m_direction;

        RegularGrid* m_grid;
        RegularGridTrilinearFunction* m_func;
        bool in_bounds(const Vec3d &x) const {
            return (0.0 <= x[0] && x[0] < 1.0 &&
                0.0 <= x[1] && x[1] < 1.0 &&
                0.0 <= x[2] && x[2] < 1.0);
        }
    public:
        std::vector<std::vector<Vec3d> > plines;

        NumericIntegratorNew(RegularGrid* grid, RegularGridTrilinearFunction* func,
            FLOATTYPE error_threshold, FLOATTYPE gradient_threshold, int max_steps, int direction) :
            m_grid(grid), m_func(func), m_save_geometry(0), m_direction(direction),
            m_gradient_threshold(gradient_threshold), m_error_threshold(error_threshold), m_max_steps(max_steps) {}

        void need_geometry(bool _ = true) { m_save_geometry = _; }
        void add_seed(const Vec3d &_) { m_seeds.push_back(_); }
        void add_seeds(const std::vector<Vec3d> &_) {
            m_seeds.insert(m_seeds.end(), _.begin(), _.end());
        }
        void set_seeds(const std::vector<Vec3d> &_) {
            m_seeds.clear();
            m_seeds.insert(m_seeds.end(), _.begin(), _.end());
        }

        MSC::ADVECTION_EVENT IntegrateStreamline(const Vec3d &seed, std::vector<Vec3d> &sline) const {
            static const FLOATTYPE stepsize = 0.25;
            static  MSC::AdaptiveEulerStepper intg(stepsize, m_error_threshold, m_gradient_threshold, m_direction);

            // ------------------------------------------------------
            MSC::ADVECTION_EVENT terminationCause = MSC::NONE;

            Vec3d currentP = seed;
            Vec3l baseVtx;

            Vec3d surrounding_grads[8];

            if (m_grid->DistToBoundary(currentP) > 1) {
                baseVtx = currentP;
                m_func->GetGradSurroundingNoBoundaryCheck(baseVtx, surrounding_grads);  // fill in local gradients
            }
            else {
                baseVtx = m_grid->Inbounds(currentP);
                m_func->GetGradSurrounding(baseVtx, surrounding_grads);                 // fill in local gradients
            }

            int nsteps = 0;
            while (true) {

                sline.push_back(currentP);

                if (++nsteps >= m_max_steps) {
                    return MSC::OVER_MAX_ITERATIONS;
                }

                Vec3d factors = currentP - baseVtx;

                // change the voxel if the current point is out of bounds!
                if (!in_bounds(factors)){

                    if (m_grid->DistToBoundary(currentP) > 1) {
                        baseVtx = currentP;
                        m_func->GetGradSurroundingNoBoundaryCheck(baseVtx, surrounding_grads);  // fill in local gradients
                    }
                    else {
                        baseVtx = m_grid->Inbounds(currentP);
                        m_func->GetGradSurrounding(baseVtx, surrounding_grads);                 // fill in local gradients
                    }
                }

                bool flag = intg.step(currentP, baseVtx, surrounding_grads);
                if (flag == false) {
                    return MSC::LOW_GRADIENT;
                }
            }
        }


        void BeginIntegration() {
            size_t nseeds = m_seeds.size();

            printf("NumericIntegratorNew::BeginIntegration() for %d seeds\n", nseeds);

            if (m_save_geometry) {
                plines.clear();
                plines.resize(m_seeds.size());
            }

#pragma omp parallel for
            for (int i = 0; i < nseeds; i++) {

                IntegrateStreamline(m_seeds.at(i), plines.at(i));
#if 0
                GInt::ADVECTION_EVENT terminationCause = GInt::NONE;

                Vec3d currentP = m_seeds.at(i);
                Vec3l baseVtx;

                Vec3d surrounding_grads[8];

                if (m_grid->DistToBoundary(currentP) > 1) {
                    baseVtx = currentP;
                    m_func->GetGradSurroundingNoBoundaryCheck(baseVtx, surrounding_grads);  // fill in local gradients
                }
                else {
                    baseVtx = m_grid->Inbounds(currentP);
                    m_func->GetGradSurrounding(baseVtx, surrounding_grads);                 // fill in local gradients
                }


                int nsteps = 0;
                while (true) {

                    if (m_save_geometry) {
                        plines.at(i).push_back(currentP);
                    }

                    if (++nsteps >= m_max_steps) {
                        terminationCause = GInt::OVER_MAX_ITERATIONS;
                        break;
                    }

                    Vec3d factors = currentP - baseVtx;
                    /*
                    printf("\n [%d] : currentP = (%.1f %.1f %.1f) : [%d %d %d]\n", count, currentP[0], currentP[1], currentP[2],
                    baseVtx[0], baseVtx[1], baseVtx[2]);
                    printf(" factors = (%f %f %f)\n", factors[0], factors[1], factors[2]);
                    */

                    //TODO --- handle periodic boundary!
                    // change the voxel if the current point is out of bounds!
                    if (!in_bounds(factors)){

                        if (m_grid->DistToBoundary(currentP) > 1) {
                            baseVtx = currentP;
                            m_func->GetGradSurroundingNoBoundaryCheck(baseVtx, surrounding_grads);  // fill in local gradients
                        }
                        else {
                            baseVtx = m_grid->Inbounds(currentP);
                            m_func->GetGradSurrounding(baseVtx, surrounding_grads);                 // fill in local gradients
                        }
                    }

                    bool flag = intg.step(currentP, baseVtx, surrounding_grads);
                    if (flag == false) {
                        terminationCause = GInt::LOW_GRADIENT;
                        break;
                    }
                }
#endif
            }
        }


    };


    template< class Advector>
    class NumericStreamlineIntegrator {


        FLOATTYPE m_error_threshold;
        FLOATTYPE m_gradient_threshold;
        int m_max_steps;

        RegularGrid* m_grid;
        RegularGridTrilinearFunction* m_func;
        bool in_bounds(const Vec3d &x) const {
            return (0.0 <= x[0] && x[0] < 1.0 &&
                0.0 <= x[1] && x[1] < 1.0 &&
                0.0 <= x[2] && x[2] < 1.0);
        }

        AdvectionChecker* m_checker;
    public:


        void SetAdvectionChecker(AdvectionChecker* checker) { m_checker = checker; }
        NumericStreamlineIntegrator(RegularGrid* grid, RegularGridTrilinearFunction* func,
            FLOATTYPE error_threshold, FLOATTYPE gradient_threshold, int max_steps) :
            m_grid(grid), m_func(func),
            m_gradient_threshold(gradient_threshold), m_error_threshold(error_threshold), m_max_steps(max_steps) {
            m_checker = new NoTermination();
        }


        MSC::ADVECTION_EVENT IntegrateStreamline(const Vec3d &seed, std::vector<Vec3d> &sline) const {

            Vec3l start_coords = (seed + 0.5); // since casting rounds down to integer... this is starting vertex
            INDEX_TYPE start_id = m_grid->Index3d(start_coords); // vertex id of start point
            int t1, t2; t1 = t2 = 0;
            Advector t_advector(m_grid, m_func, m_gradient_threshold, m_error_threshold, m_checker);
            sline.reserve(100);

            sline.push_back(seed);

            Vec3l t_coords = start_coords; // get the coordinates of the poitn
            //if (t_coords[0] == 0 && t_coords[1] == 0) printf("doing %d\n", t_coords[2]);

            //BUGFIX. starting point should be the seed
            //Vec3d t_current_point = t_coords;
            Vec3d t_current_point = seed;
            int t_num_iterations_left = m_max_steps;
            bool t_continue = true;

            while (t_continue) {
                Vec3d t_next_point;
                ADVECTION_EVENT t_return_code;
                if (m_grid->DistToBoundary(t_coords) <= 1) {
                    t_return_code = t_advector.AdvectThroughVoxelNearBoundary(t_current_point, t_num_iterations_left);
                    t_coords = m_grid->Inbounds(t_current_point + 0.5); // get nearest integer voxel
                    t1++;
                }
                else {
                    t_return_code = t_advector.AdvectThroughVoxelNoCheck(t_current_point, t_num_iterations_left);
                    t_coords = (t_current_point + 0.5);
                    t2++;
                }
                sline.push_back(t_current_point);
                INDEX_TYPE t_next_id = m_grid->Index3d(t_coords);

                // if we terminated or hit a critical point, then we are done
                if (t_return_code == ADVECTION_EVENT::LOW_GRADIENT ||
                    t_return_code == ADVECTION_EVENT::HIT_EXTREMUM ||
                    t_return_code == ADVECTION_EVENT::HIT_PREASSIGNED ||
                    t_return_code == ADVECTION_EVENT::OVER_MAX_ITERATIONS) {


                    return t_return_code;
                    t_continue = false;
                }
            }
        }



    };




}
#endif
