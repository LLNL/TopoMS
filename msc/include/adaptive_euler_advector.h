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

#ifndef ADAPTIVE_EULER_ADVECTOR_H
#define ADAPTIVE_EULER_ADVECTOR_H

#include "basic_types.h"
#include "vectors.h"
#include "regular_grid.h"
#include "regular_grid_trilinear_function.h"
#include "advection_checker.h"
#include "advection_events.h"

#ifndef MAX_STEP_IN_GRID_CELLS
#define MAX_STEP_IN_GRID_CELLS 0.5
#endif
namespace MSC {



    template< int StepMultiplier>
    class AdaptiveEulerAdvector {
    protected:
        RegularGrid* m_grid;
        RegularGridTrilinearFunction* m_func;
        FLOATTYPE m_gradient_threshold;
        FLOATTYPE m_error_threshold;

        AdvectionChecker* m_inside_voxel_advection_checker;

#undef min
        inline FLOATTYPE min(FLOATTYPE a, FLOATTYPE b) const {
            if (a < b) return a;
            return b;
        }


        //int m_num_noboundary;
        //int m_num_boundary;


    public:



        AdaptiveEulerAdvector(RegularGrid* grid, RegularGridTrilinearFunction* func, FLOATTYPE grad_thresh, FLOATTYPE error_thresh, AdvectionChecker* inside_voxel_advection_checker) : m_grid(grid), m_func(func), m_gradient_threshold(grad_thresh), m_error_threshold(error_thresh), m_inside_voxel_advection_checker(inside_voxel_advection_checker){
            //m_num_boundary = 0;
            //m_num_noboundary = 0;
        }

        // starts at a point, identifies which grid voxel it is in, and advects the point until it either exits the voxel,
        // or hits an early termination condition. sets the output point to the new location, or the last location
        // before early termination. also sets the reason for return - either exiting voxel or early termination
        ADVECTION_EVENT AdvectThroughVoxelNoCheck(Vec3d& current_point, int& iterations_left) {
            //m_num_boundary++;
            Vec3d surrounding_grad[8]; // used to cache surrounding gradient
            Vec3d old_point = current_point; // set the last poitn to current one since we initialize cache
            Vec3d half_step_point_1; // these points are used to in the 2-step euler

            // declare variables to store gradient, and initialize with gradient at start position
            Vec3l old_point_int_base = current_point; // figure out what the closest vertex is we are in
            m_func->GetGradSurroundingNoBoundaryCheck(old_point_int_base, surrounding_grad); // fill in local gradients
            Vec3d gradient = m_func->TriLinInterpGrad(current_point, old_point_int_base, surrounding_grad); // cache old gradient for old_point - used when we change stepsize but not location

            Vec3l half_step_point_1_int_base(-1, -1, -1); // used to figure out if cached half-step gradients can be reused
            Vec3d old_half_step_point_1_surrounding[8]; // cache surrounding gradeints for half-step

            FLOATTYPE gradient_magnitude = gradient.Mag(); // cache old gradient magnitude
            FLOATTYPE stepsize = MAX_STEP_IN_GRID_CELLS / gradient_magnitude; // estimate stepsize based on gradient magnitude- we dont wnat to step too far!
            while (iterations_left-- > 0) {

                // do we need to resample the gradient?
                if (!(old_point == current_point)) {
                    // INSERT TERMINATING CONDITION HERE - WE MOVED POINT -- RECOMMENT -> NO, CHECK HAS HAPPENED IF POINT WAS MOVED!!!!!
                    // yes! we need to sample the gradient since it is not cached

                    // THIS IS WRONG!! WILL GIVE WRONG RESULT IF WE MOVE ACROSS 0 PLANE! - THE INT CONVERSION
                    // TAKES US TOWARD ZERO, NOT NEGATIVE INFINITY!@!!!!!!!
                    if (((int)current_point[0]) != ((int)old_point[0]) ||
                        ((int)current_point[1]) != ((int)old_point[1]) ||
                        ((int)current_point[2]) != ((int)old_point[2]))
                        return OUT_OF_VOXEL;

                    //Vec3l new_base = current_point.IntFloor();
                    //if (!(new_base == old_point_int_base)) {
                    //	// we have exited the current voxel
                    //	return OUT_OF_VOXEL;
                    //}

                    gradient = m_func->TriLinInterpGrad(current_point, old_point_int_base, surrounding_grad);
                    gradient_magnitude = gradient.Mag();
                    // also adjust stepsize
                    if (gradient_magnitude * stepsize > MAX_STEP_IN_GRID_CELLS) stepsize = MAX_STEP_IN_GRID_CELLS / gradient_magnitude; // restrict stepsize if too big

                    old_point = current_point;
                }

                // check flatness! return if we are too flat
                if (gradient_magnitude < m_gradient_threshold) {
                    return LOW_GRADIENT;

                }

                // do simple 1 step euler
                // get the displacement of 1 euler step
                //full_step_point = current_point + (gradient * stepsize);

                Vec3d euler_displacement = gradient * (StepMultiplier * stepsize);

                // now do 2 steps with half sized steps to see if we are in error tolerance
                Vec3d euler_2step_displacement;

                //// incorrect method: just re-use gradient and trilininterp
                {
                    // get the displacement we have for the first half step
                    Vec3d half_step_displacement = (gradient * (stepsize * (StepMultiplier * 0.5)));
                    // compute the half-step point
                    half_step_point_1 = current_point + half_step_displacement;
                    // do second step NOT CARING if we exited the voxel, will just do strange stuff with trilininterp

                    Vec3d half_step_gradient = m_func->TriLinInterpGrad(half_step_point_1, old_point_int_base, surrounding_grad);




                    // get displacement for the two half-steps
                    euler_2step_displacement = (half_step_gradient * (stepsize * (StepMultiplier * 0.5))) + half_step_displacement;
                }

                //{
                //	// get the displacement we have for the first half step
                //	Vec3d half_step_displacement = (gradient * (stepsize * 0.5));
                //	half_step_point_1 = current_point + half_step_displacement;
                //	// to do the second step, we may have exited the voxel, so re-check gradient
                //	Vec3l new_base = half_step_point_1.IntFloor();
                //	Vec3d half_step_gradient;
                //	if (!(new_base == old_point_int_base)) {
                //		if (!(new_base == half_step_point_1_int_base)) {
                //			// cache intermediate gradient
                //			m_func->GetGradSurroundingNoBoundaryCheck(new_base, old_half_step_point_1_surrounding);
                //			half_step_point_1_int_base = new_base;
                //		}
                //		// get the intermediate gradient
                //		half_step_gradient = m_func->TriLinInterpGrad(half_step_point_1, new_base, old_half_step_point_1_surrounding);
                //	}
                //	else {
                //		// get the intermediate gradient using the same neighborhood as the initial one
                //		half_step_gradient = m_func->TriLinInterpGrad(half_step_point_1, old_point_int_base, surrounding_grad);
                //	}
                //
                //	// get displacement for the two half-steps
                //	euler_2step_displacement = (half_step_gradient * (stepsize * 0.5)) + half_step_displacement;
                //
                //}

                // now compare the points
                FLOATTYPE local_error = (euler_2step_displacement - euler_displacement).MagSq();
                if (local_error < m_error_threshold) {
                    // we are within tolerance so accept the new location, and increase stepsize
                    current_point = current_point + euler_2step_displacement;
                    stepsize = min(stepsize * 1.8, MAX_STEP_IN_GRID_CELLS / gradient_magnitude);

                    // INSERT ADVECTION CHECK HERE
                    ADVECTION_EVENT t_inside_event = m_inside_voxel_advection_checker->CheckAndDoStuff(old_point, current_point, false);
                    if (t_inside_event != NONE)
                        return t_inside_event;

                }
                else {
                    // the error is too high! reduce stepsize and try again
                    stepsize = stepsize * 0.5;
                }

            }
            // only get here if we ran out of iterations inside the block
            return OVER_MAX_ITERATIONS;


        }


        // starts at a point, identifies which grid voxel it is in, and advects the point until it either exits the voxel,
        // or hits an early termination condition. sets the output point to the new location, or the last location
        // before early termination. also sets the reason for return - either exiting voxel or early termination
        ADVECTION_EVENT AdvectThroughVoxelNearBoundary(Vec3d& current_point, int& iterations_left) {
            //m_num_noboundary++;
            Vec3d surrounding_grad[8]; // used to cache surrounding gradient
            Vec3d old_point = current_point; // set the last poitn to current one since we initialize cache
            Vec3d half_step_point_1; // these points are used to check the error of the current stepsize

            // declare variables to store gradient, and initialize with gradient at start position
            Vec3l old_point_int_base = m_grid->Inbounds(current_point); // figure out what the closest vertex is we are in
            m_func->GetGradSurrounding(old_point_int_base, surrounding_grad); // fill in local gradients
            Vec3d gradient = m_func->TriLinInterpGrad(current_point, old_point_int_base, surrounding_grad); // cache old gradient for old_point - used when we change stepsize but not location

            Vec3l half_step_point_1_int_base(-1, -1, -1); // used to figure out if cached half-step gradients can be reused
            Vec3d old_half_step_point_1_surrounding[8]; // cache surrounding gradeints for half-step

            FLOATTYPE gradient_magnitude = gradient.Mag(); // cache old gradient magnitude
            FLOATTYPE stepsize = MAX_STEP_IN_GRID_CELLS / gradient_magnitude; // estimate stepsize based on gradient magnitude- we dont wnat to step too far!
            while (iterations_left-- > 0) {

                // do we need to resample the gradient?
                if (!(old_point == current_point)) {
                    // yes! we need to sample the gradient since it is not cached
                    if (((int)current_point[0]) != ((int)old_point[0]) ||
                        ((int)current_point[1]) != ((int)old_point[1]) ||
                        ((int)current_point[2]) != ((int)old_point[2]))
                        return OUT_OF_VOXEL;
                    //Vec3l new_base = current_point.IntFloor();
                    //if (!(new_base == old_point_int_base)) {
                    //	// we have exited the current voxel
                    //	return OUT_OF_VOXEL;
                    //}

                    gradient = m_func->TriLinInterpGrad(current_point, old_point_int_base, surrounding_grad);
                    gradient_magnitude = gradient.Mag();
                    // also adjust stepsize
                    if (gradient_magnitude * stepsize > MAX_STEP_IN_GRID_CELLS) stepsize = MAX_STEP_IN_GRID_CELLS / gradient_magnitude; // restrict stepsize if too big

                    old_point = current_point;
                }

                // check flatness! return if we are too flat
                if (gradient_magnitude < m_gradient_threshold) {
                    return LOW_GRADIENT;

                }

                // get the displacement of 1 euler step
                Vec3d euler_displacement = gradient * (StepMultiplier * stepsize);

                // now do 2 steps with half sized steps to see if we are in error tolerance
                Vec3d euler_2step_displacement;


                //// incorrect method: just re-use gradient and trilininterp
                {
                    // get the displacement we have for the first half step
                    Vec3d half_step_displacement = (gradient * (stepsize * (StepMultiplier * 0.5)));
                    // compute the half-step point
                    half_step_point_1 =current_point + half_step_displacement;
                    // do second step NOT CARING if we exited the voxel, will just do strange stuff with trilininterp

                    Vec3d half_step_gradient = m_func->TriLinInterpGrad(half_step_point_1, old_point_int_base, surrounding_grad);




                    // get displacement for the two half-steps
                    euler_2step_displacement = (half_step_gradient * (stepsize * (StepMultiplier * 0.5))) + half_step_displacement;
                }

                ////// correct method:
                //{
                //	// get the displacement we have for the first half step
                //	Vec3d half_step_displacement = (gradient * (stepsize * 0.5));
                //	// compute the half-step point
                //	half_step_point_1 = m_grid->Inbounds(current_point + half_step_displacement);
                //	// to do the second step, we may have exited the voxel, so re-check gradient
                //
                //	Vec3l new_base = half_step_point_1.IntFloor();
                //	Vec3d half_step_gradient;
                //	if (!(new_base == old_point_int_base)) {
                //		if (!(new_base == half_step_point_1_int_base)) {
                //			// cache intermediate gradient
                //			m_func->GetGradSurrounding(new_base, old_half_step_point_1_surrounding);
                //			half_step_point_1_int_base = new_base;
                //		}
                //		// get the intermediate gradient
                //		half_step_gradient = m_func->TriLinInterpGrad(half_step_point_1, new_base, old_half_step_point_1_surrounding);
                //	}
                //	else {
                //		// get the intermediate gradient using the same neighborhood as the initial one
                //		half_step_gradient = m_func->TriLinInterpGrad(half_step_point_1, new_base, surrounding_grad);
                //	}
                //
                //
                //
                //	// get displacement for the two half-steps
                //	euler_2step_displacement = (half_step_gradient * (stepsize * 0.5)) + half_step_displacement;
                //}

                // now compare the points
                FLOATTYPE local_error = (euler_2step_displacement - euler_displacement).MagSq();
                if (local_error < m_error_threshold) {
                    // we are within tolerance so accept the new location, and increase stepsize
                    current_point = m_grid->Inbounds(current_point + euler_2step_displacement);
                    stepsize = min(stepsize * 1.8, MAX_STEP_IN_GRID_CELLS / gradient_magnitude);

                    //  now check if moving the point internally has any effect
                    ADVECTION_EVENT t_inside_event = m_inside_voxel_advection_checker->CheckAndDoStuff(old_point, current_point, true);
                    if (t_inside_event != NONE)
                        return t_inside_event;


                }
                else {
                    // the error is too high! reduce stepsize and try again
                    stepsize = stepsize * 0.5;
                }

            }
            // only get here if we ran out of iterations inside the block

            return OVER_MAX_ITERATIONS;

        }

    };

}

#endif
