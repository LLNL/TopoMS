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

#ifndef ADAPTIVE_IN_QUAD_EULER_ADVECTOR_H
#define ADAPTIVE_IN_QUAD_EULER_ADVECTOR_H

#include "basic_types.h"
#include "vectors.h"
#include "regular_grid.h"
#include "regular_grid_trilinear_function.h"
#include "topological_regular_grid.h"
#include "advection_checker.h"
#include "advection_events.h"

#ifndef MAX_STEP_IN_GRID_CELLS
#define MAX_STEP_IN_GRID_CELLS 0.25
#endif
namespace MSC {

    class AdaptiveInQuadEulerAdvector {
    protected:
        //RegularGrid* m_grid;
        //
        //RegularGridTrilinearFunction* m_func;
        FLOATTYPE m_gradient_threshold;
        FLOATTYPE m_error_threshold;

        AdvectionChecker2D* m_inside_voxel_advection_checker;

        inline FLOATTYPE min(FLOATTYPE a, FLOATTYPE b) const {
            if (a < b) return a;
            return b;
        }


        //int m_num_noboundary;
        //int m_num_boundary;

        // this assumes p is already 0-1
        Vec2d BilinearInterp(Vec2d& p, Vec2d* surrounding) {
            Vec2d t0 = (surrounding[1] * p[0]) + (surrounding[0] * (1 - p[0]));
            Vec2d t1 = (surrounding[3] * p[0]) + (surrounding[2] * (1 - p[0]));
            return (t1 * p[1]) + (t0 * (1 - p[1]));

        }

        // tests intersection between line segment and coodinate plane, and if it does, returns other
        // coordinate of that intersection
        inline bool Intersects(const Vec2d& oldpoint, const Vec2d& currentpoint, int plane, float value) const {
            return (oldpoint[plane] <= value && currentpoint[plane] >= value) ||
                (oldpoint[plane] >= value && currentpoint[plane] <= value);
        }
        inline FLOATTYPE IntersectionPoint(const Vec2d& oldpoint, const Vec2d& currentpoint, int plane, float value) const {
            int other_plane = (plane + 3) % 2; // plane is 1 or 0, 1 + 3 = 4, 4 % 2 = 0, 0+3 = 3, 3%2 = 1
            FLOATTYPE plane_length =  currentpoint[plane]  - oldpoint[plane];
            // if collinear with plane, then simply return old points other coordinate
            // but this aught not to happen
            if (plane_length == 0)
                return oldpoint[other_plane];
            FLOATTYPE t = (value - oldpoint[plane]) / plane_length;
            FLOATTYPE ret = t * (currentpoint[other_plane] - oldpoint[other_plane]) + oldpoint[other_plane];
            if (ret < 0 || ret > 1)
                printf("whoa\n");
            return ret;
        }
    public:
        enum PLANEKIND { X0PLANE, X1PLANE, Y0PLANE, Y1PLANE, NONE, PLANE_ERROR };
    protected:

#define HACK_EPSILON 0.01
        PLANEKIND WhichPlane(Vec2d& oldpoint, Vec2d& currentpoint) {
            bool cross_x0 = false;
            bool cross_x1 = false;
            bool cross_y0 = false;
            bool cross_y1 = false;

            // test if it crosses a plane

            FLOATTYPE crosspx0 = -1.0;
            if(Intersects(oldpoint, currentpoint, 0, 0.0)) {
                cross_x0 = true;
                crosspx0 = IntersectionPoint(oldpoint, currentpoint, 0, 0.0);
                if (crosspx0 >= 0.0 - HACK_EPSILON && crosspx0 <= 1.0 + HACK_EPSILON) return X0PLANE;
            }
            FLOATTYPE crosspx1 = -1.0;
            if (Intersects(oldpoint, currentpoint, 0, 1.0)) {
                cross_x1 = true;
                crosspx1 = IntersectionPoint(oldpoint, currentpoint, 0, 1.0);
                if (crosspx1 >= 0.0 - HACK_EPSILON && crosspx1 <= 1.0 + HACK_EPSILON) return X1PLANE;
            }
            FLOATTYPE crosspy0 = -1.0;
            if (Intersects(oldpoint, currentpoint, 1, 0.0)) {
                cross_y0 = true;
                crosspy0 = IntersectionPoint(oldpoint, currentpoint, 1, 0.0);
                if (crosspy0 >= 0.0 - HACK_EPSILON && crosspy0 <= 1.0 + HACK_EPSILON) return Y0PLANE;
            }
            FLOATTYPE crosspy1 = -1.0;
            if (Intersects(oldpoint, currentpoint, 1, 1.0)) {
                cross_y1 = true;
                crosspy1 = IntersectionPoint(oldpoint, currentpoint, 1, 1.0);
                if (crosspy1 >= 0.0 - HACK_EPSILON && crosspy1 <= 1.0 + HACK_EPSILON) return Y1PLANE;
            }

            if (cross_x0 || cross_x1 || cross_y0 || cross_y1) {
                printf("ERROR: detected crossing but NONE from inside??\n");
                return PLANE_ERROR;
            }

            return NONE;
        }

    public:


        AdaptiveInQuadEulerAdvector(RegularGrid* grid, RegularGridTrilinearFunction* func, FLOATTYPE grad_thresh,
            FLOATTYPE error_thresh, AdvectionChecker2D* inside_voxel_advection_checker) :
            m_gradient_threshold(grad_thresh), m_error_threshold(error_thresh),
            m_inside_voxel_advection_checker(inside_voxel_advection_checker){
            //m_num_boundary = 0;
            //m_num_noboundary = 0;
        }

        void SetInteriorChecker(AdvectionChecker2D* checker) { m_inside_voxel_advection_checker = checker; }

        // starts at a point, identifies which grid voxel it is in, and advects the point until it either exits the voxel,
        // or hits an early termination condition. sets the output point to the new location, or the last location
        // before early termination. also sets the reason for return - either exiting voxel or early termination
        // exit plane is filled in with a PLANEKIND which says which plane the integration exited through.
        ADVECTION_EVENT AdvectThroughQuadNoCheck(Vec2d& current_point, Vec2d surrounding_grad[4],
            int& iterations_left, PLANEKIND& exit_plane) {
            Vec2d old_point = current_point; // set the last poitn to current one since we initialize cache
            Vec2d half_step_point_1; // these points are used to in the 2-step euler

            // declare variables to store gradient, and initialize with gradient at start position
            Vec2d gradient = BilinearInterp(current_point, surrounding_grad);

            FLOATTYPE gradient_magnitude = gradient.Mag(); // cache old gradient magnitude
            FLOATTYPE stepsize = MAX_STEP_IN_GRID_CELLS / gradient_magnitude; // estimate stepsize based on gradient magnitude- we dont wnat to step too far!
            while (iterations_left-- > 0) {

                // do we need to resample the gradient?
                if (!(old_point == current_point)) {
                    // INSERT TERMINATING CONDITION HERE - WE MOVED POINT -- RECOMMENT -> NO, CHECK HAS HAPPENED IF POINT WAS MOVED!!!!!
                    // yes! we need to sample the gradient since it is not cached

                    if (current_point[0] < 0.0 || current_point[0] > 1.0 ||
                        current_point[1] < 0.0 || current_point[1] > 1.0 ) {
                        exit_plane = WhichPlane(old_point, current_point);
                        return OUT_OF_VOXEL;
                    }

                    gradient = BilinearInterp(current_point, surrounding_grad);
                    gradient_magnitude = gradient.Mag();
                    // also adjust stepsize
                    if (gradient_magnitude * stepsize > MAX_STEP_IN_GRID_CELLS) stepsize = MAX_STEP_IN_GRID_CELLS / gradient_magnitude; // restrict stepsize if too big

                    old_point = current_point;
                }

                // check flatness! return if we are too flat
                if (gradient_magnitude < m_gradient_threshold) {
                    exit_plane = NONE;
                    printf("low_grad\n");
                    return LOW_GRADIENT;

                }

                // do simple 1 step euler
                // get the displacement of 1 euler step
                //full_step_point = current_point + (gradient * stepsize);

                Vec2d euler_displacement = gradient * stepsize;

                // now do 2 steps with half sized steps to see if we are in error tolerance
                Vec2d euler_2step_displacement;

                //// incorrect method: just re-use gradient and trilininterp
                {
                    // get the displacement we have for the first half step
                    Vec2d half_step_displacement = (gradient * (stepsize * 0.5));
                    // compute the half-step point
                    half_step_point_1 = current_point + half_step_displacement;
                    // do second step NOT CARING if we exited the voxel, will just do strange stuff with trilininterp

                    Vec2d half_step_gradient = BilinearInterp(half_step_point_1, surrounding_grad);

                    // get displacement for the two half-steps
                    euler_2step_displacement = (half_step_gradient * (stepsize * 0.5)) + half_step_displacement;
                }
                // now compare the points
                FLOATTYPE local_error = (euler_2step_displacement - euler_displacement).MagSq();
                if (local_error < m_error_threshold) {
                    // we are within tolerance so accept the new location, and increase stepsize
                    current_point = current_point + euler_2step_displacement;
                    stepsize = min(stepsize * 1.8, MAX_STEP_IN_GRID_CELLS / gradient_magnitude);

                    // INSERT ADVECTION CHECK HERE
                    ADVECTION_EVENT t_inside_event = m_inside_voxel_advection_checker->CheckAndDoStuff(old_point, current_point, false);
                    if (t_inside_event != ADVECTION_EVENT::NONE) {
                        exit_plane = NONE;
                        printf("advection vent\n");
                        return t_inside_event;
                    }

                }
                else {
                    // the error is too high! reduce stepsize and try again
                    stepsize = stepsize * 0.5;
                }

            }
            // only get here if we ran out of iterations inside the block
            exit_plane = NONE;
            printf("over_max_iter\n");
            return OVER_MAX_ITERATIONS;


        }



    };

}

#endif
