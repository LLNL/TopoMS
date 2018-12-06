/*
Copyright (c) 2018, Lawrence Livermore National Security, LLC.
Produced at the Lawrence Livermore National Laboratory.
Written by Harsh Bhatia (hbhatia@llnl.gov) and Attila G Gyulassy
(jediati@sci.utah.edu).
LLNL-CODE-745278. All rights reserved.

This file is part of TopoMS, Version 1.0. For details, see
https://github.com/LLNL/TopoMS. Please also read this link – Additional BSD
Notice.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

• Redistributions of source code must retain the above copyright notice, this
list of conditions and the disclaimer below.
• Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the disclaimer (as noted below) in the
documentation and/or other materials provided with the distribution.
• Neither the name of the LLNS/LLNL nor the names of its contributors may be
used to endorse or promote products derived from this software without specific
prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY, LLC, THE
U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

Additional BSD Notice

1. This notice is required to be provided under our contract with the U.S.
Department of Energy (DOE). This work was produced at Lawrence Livermore
National Laboratory under Contract No.  DE-AC52-07NA27344 with the DOE.

2. Neither the United States Government nor Lawrence Livermore National
Security, LLC nor any of their employees, makes any warranty, express or
implied, or assumes any liability or responsibility for the accuracy,
completeness, or usefulness of any information, apparatus, product, or process
disclosed, or represents that its use would not infringe privately-owned
rights.

3. Also, reference herein to any specific commercial products, process, or
services by trade name, trademark, manufacturer or otherwise does not
necessarily constitute or imply its endorsement, recommendation, or favoring by
the United States Government or Lawrence Livermore National Security, LLC. The
views and opinions of authors expressed herein do not necessarily state or
reflect those of the United States Government or Lawrence Livermore National
Security, LLC, and shall not be used for advertising or product endorsement
purposes.
*/

/**
 *  @file    MolecularSystem.h
 *  @author  Harsh Bhatia (hbhatia@llnl.gov)
 *  @date    10/01/2017
 *
 *  @brief This class handles the meta data about the system under investigration
 *
 *  @section DESCRIPTION
 *
 *  This class handles the meta data about the system under investigration
 *
 */

#ifndef _MOLECULAR_SYSTEM_H_
#define _MOLECULAR_SYSTEM_H_

#include <cmath>
#include <string>
#include <vector>
#include <cstdint>
#include <iostream>

#include "Vec3.h"
#include "Mat3.h"

namespace MS {

class SystemInfo {

    public:
        // name of the system
        std::string m_sysname;
        std::vector<std::pair<std::string, unsigned int> > m_materials;

        // units of physical length (Bohr or Angs)
        double m_length_file2Angs;
        std::string m_length_unit;

        // units of charge (hartree or eV)
        double m_charge_file2electrons;
        std::string m_charge_unit;

        // information about spatial domain
        double m_scaling;
        Vec3<double> m_lattice_origin;
        Mat3<double> m_lattice;
        Mat3<double> m_lattice_inv;
        Vec3<size_t> m_grid_dims;

        std::string m_coordinate_type;    // direct, world, grid ?   (useful only for VASP)

        // information about time
        size_t m_num_tsteps;
        double m_time_step;
        std::string m_time_unit;
        unsigned int m_time_factor;

        // other information
        float m_avg_atomic_vol;
        float m_temperature;

        // constructor
        SystemInfo() : m_sysname("unknown"),
                        m_coordinate_type("unknown"), m_length_unit("unknown"), m_charge_unit("unknown"),
                        m_charge_file2electrons(1.0), m_length_file2Angs(1.0), m_scaling(1.0),
                        m_num_tsteps(1), m_time_step(1.0), m_time_unit("unknown"), m_time_factor(0),
                        m_avg_atomic_vol(0), m_temperature(0) {
            m_lattice.eye();
            m_lattice_inv.eye();
        }

        // some simple utils for metadata
        void print() const {
            std::cout << "    System name = <" << m_sysname <<">\n";
            std::cout << "    Lattice origin = (" << m_lattice_origin[0] << ", " << m_lattice_origin[1] << ", " << m_lattice_origin[2] << ") " << m_length_unit << "\n"
                      << "    Lattice matrix = (" << m_lattice.v[0][0] << ", " << m_lattice.v[0][1] << ", " << m_lattice.v[0][2] << ") " << m_length_unit << "\n"
                      << "                     (" << m_lattice.v[1][0] << ", " << m_lattice.v[1][1] << ", " << m_lattice.v[1][2] << ") " << m_length_unit << "\n"
                      << "                     (" << m_lattice.v[2][0] << ", " << m_lattice.v[2][1] << ", " << m_lattice.v[2][2] << ") " << m_length_unit << "\n"
                      << "    Lattice volume = " << this->volume_box() << " " << m_length_unit << "^3\n"
                      << "    Grid           = ["<< m_grid_dims[0] <<" x " << m_grid_dims[1] <<" x " << m_grid_dims[2] <<"]\n"
                      << "    Voxel   volume = " << this->volume_voxel() << " " << m_length_unit << "^3\n";
            print_time();
        }
        void print_time() const {
            std::cout << "    # time_steps = " << m_num_tsteps << ", time_step = " << m_time_step <<" " << m_time_unit << std::endl;
        }

        bool operator == (const SystemInfo &m) const {

            static const float eps = 0.00001;

            if (m_sysname.compare(m.m_sysname) != 0)                                return false;
            if (m_length_unit.compare(m.m_length_unit) != 0)                        return false;
            if (m_charge_unit.compare(m.m_charge_unit) != 0)                        return false;
            if (m_time_unit.compare(m.m_time_unit) != 0)                            return false;
            if (m_coordinate_type.compare(m.m_coordinate_type) != 0)                return false;

            if (fabs(m_length_file2Angs - m.m_length_file2Angs) > eps)              return false;
            if (fabs(m_charge_file2electrons - m.m_charge_file2electrons) > eps)    return false;

            if (fabs(m_scaling - m.m_scaling) > eps)                                return false;
            if (fabs(m_time_step - m.m_time_step) > eps)                            return false;

            if (fabs(m_temperature - m.m_temperature) > eps)                        return false;
            if (fabs(m_avg_atomic_vol - m.m_avg_atomic_vol) > eps)                  return false;

            if (m_num_tsteps != m.m_num_tsteps)                                     return false;
            if (m_time_factor != m.m_time_factor)                                   return false;

            for(uint8_t d = 0; d < 3; d++) {

                if (m_grid_dims[d] != m.m_grid_dims[d])                             return false;
                if (fabs(m_lattice_origin[d] - m.m_lattice_origin[d]) > eps)        return false;
                if (fabs(m_lattice.v[d][0] - m.m_lattice.v[d][0]) > eps)            return false;
                if (fabs(m_lattice.v[d][1] - m.m_lattice.v[d][1]) > eps)            return false;
                if (fabs(m_lattice.v[d][2] - m.m_lattice.v[d][2]) > eps)            return false;
            }
            return true;
        }
        // -----------------------------------------------------------------------------
        // -----------------------------------------------------------------------------
        // conversion between various types of spatial coordinates!
            // direct :: in range [0, 1]
            // grid   :: in range [0, grid_dims]
            // world  :: in range [lattice_origins, lattice_origins+lattice_vectors]

        void invert_lattice() { m_lattice_inv = m_lattice.inverse();    }

        /* Direct <--> Grid */
        inline void direct_to_grid(const float in_directc[3], float out_gridc[3]) const {
            for(uint8_t d = 0; d < 3; d++)      out_gridc[d] = in_directc[d] * m_grid_dims[d];
        }
        inline void grid_to_direct(const float in_gridc[3], float out_directc[3]) const {
            for(uint8_t d = 0; d < 3; d++)      out_directc[d] = in_gridc[d] / m_grid_dims[d];
        }

        /* Direct <--> World */
        inline void direct_to_world(const float in_directc[3], float out_worldc[3]) const {
            for(uint8_t d = 0; d < 3; d++){
                out_worldc[d] = in_directc[0]*m_lattice.v[0][d] +
                                in_directc[1]*m_lattice.v[1][d] +
                                in_directc[2]*m_lattice.v[2][d] +
                                                                  m_lattice_origin[d];
            }
        }
        inline void world_to_direct(const float in_worldc[3], float out_directc[3]) const {
            for(uint8_t d = 0; d < 3; d++){
                out_directc[d] = in_worldc[0]*m_lattice_inv.v[0][d] +
                                 in_worldc[1]*m_lattice_inv.v[1][d] +
                                 in_worldc[2]*m_lattice_inv.v[2][d] -
                                                                      m_lattice_origin[d];
            }
        }

        /* Grid <--> World */
        inline void grid_to_world(const float in_gridc[3], float out_worldc[3]) const {
            float directc[3];
            grid_to_direct(in_gridc, directc);
            direct_to_world(directc, out_worldc);
        }
        inline void world_to_grid(const float in_worldc[3], float out_gridc[3]) const {
            float directc[3];
            world_to_direct(in_worldc, directc);
            direct_to_grid(directc, out_gridc);
        }

        inline void world_to_grid_cuboid_lattice(float out_grid[3]) const {
            for(uint8_t i = 0; i < 3; i++)
                out_grid[i] = m_lattice.v[i][i] / float(m_grid_dims[i]);
        }

        inline void world_dims(float wdims[3]) const {
            float ddims[3] = {1,1,1};
            direct_to_world(ddims, wdims);
        }

        inline size_t grid_sz() const {     return m_grid_dims[0]*m_grid_dims[1]*m_grid_dims[2];    }

        // conversion factor of volume/charge from file units to Angstroms/electrons
        inline bool is_lunit_Angstrom() const {    return m_length_unit.compare("Ang") == 0;        }
        inline bool is_cunit_electron() const {    return m_charge_unit.compare("eV") == 0;         }

        inline double volume_file2Angs() const {
            static double f2a = -1.0;
            if (f2a < 0)
                f2a = 1.0 / std::pow(m_length_file2Angs, 3.0);
            return f2a;
        }
        inline double charge_file2electrons() const {
            static double f2e = -1.0;
            if (f2e < 0)
                f2e = 1.0 / m_charge_file2electrons;
            return f2e;
        }

        inline double volume_box() const {
            static double vol = -1.0;
            if (vol < 0)
                vol = fabs(m_lattice.determinant());
            return vol;
        }

        inline double volume_voxel() const {
            static double vol = -1.0;
            if (vol < 0) {
                const float gridc[3] = {1.0,1.0,1.0};
                float worldc[3] = {1.0,1.0,1.0};
                grid_to_world(gridc, worldc);
                vol = worldc[0]*worldc[1]*worldc[2];
            }
            return vol;
        }

        inline double volume_box_in_Angs() const {
            static double vol = volume_box();
            static double f2a = volume_file2Angs();
            return vol*f2a;
        }

        inline double volume_voxel_in_Angs() const {
            static double vol = volume_voxel();
            static double f2a = volume_file2Angs();
            return vol*f2a;
        }

        // -----------------------------------------------------------------------------
        // conversion between xyz and grid idx
        // -----------------------------------------------------------------------------

        template<typename T>
        inline void idx_to_grid(size_t idx, T out_gridc[3]) const {

            out_gridc[2] = idx / (m_grid_dims[0]*m_grid_dims[1]);
            int xy = idx % (m_grid_dims[0]*m_grid_dims[1]);

            out_gridc[1] = xy / m_grid_dims[0];
            out_gridc[0] = xy % m_grid_dims[0];
        }
        template<typename T>
        inline size_t grid_to_idx(T x, T y, T z) const {
            return z*m_grid_dims[0]*m_grid_dims[1] + y*m_grid_dims[0] + x;
        }

        // -----------------------------------------------------------------------------
        // conversion between time and tidx
        // -----------------------------------------------------------------------------

        inline double tidx_to_time(size_t tidx) const {     return double(tidx)*m_time_step;    }
        inline size_t time_to_tidx(double time) const {     return size_t(time/m_time_step);    }

        // -----------------------------------------------------------------------------
        // -----------------------------------------------------------------------------
    };
}
#endif
