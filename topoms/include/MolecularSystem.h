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
 *  @version 1.0
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

#include <vector>
#include <string>
#include <cmath>
#include <iostream>

namespace MS {

    class SystemInfo {

    public:
        // name of the system
        std::string m_sysname;
        std::vector<std::pair<std::string, unsigned int> > m_materials;

        // units of physical length (Bohr or Angs)
        double m_l_unit;
        std::string m_coordinate_unit;

        // units of charge (hartree or eV)
        double m_e_unit;
        std::string m_charge_unit;

        // information about spatial domain
        double m_scaling;
        double m_lattice_vectors[3][3];
        double m_lattice_origin[3];
        size_t m_grid_dims[3];

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
                        m_coordinate_type("unknown"), m_coordinate_unit("unknown"), m_charge_unit("unknown"),
                        m_e_unit(1.0), m_l_unit(1.0), m_scaling(1.0),
                        m_num_tsteps(0), m_time_step(1.0), m_time_unit("unknown"), m_time_factor(0),
                        m_avg_atomic_vol(0), m_temperature(0) {

            for(unsigned int i = 0; i < 3; i++){
                m_grid_dims[i] = 0;
                m_lattice_origin[i] = 0.0;
                for(unsigned int j = 0; j < 3; j++){
                    m_lattice_vectors[i][j] = 1.0;
                    m_lattice_vectors[i][(j+1)%3] = 0.0;
                    m_lattice_vectors[i][(j+2)%3] = 0.0;
                }
            }
        }

        // some simple utils for metadata
        void print() const {
            std::cout << "    System name = <" << m_sysname <<">\n";
            std::cout << "    Lattice vector x = (" << m_lattice_vectors[0][0] << ", " << m_lattice_vectors[0][1] << ", " << m_lattice_vectors[0][2] << ") " << m_coordinate_unit << "\n"
                      << "    Lattice vector y = (" << m_lattice_vectors[1][0] << ", " << m_lattice_vectors[1][1] << ", " << m_lattice_vectors[1][2] << ") " << m_coordinate_unit << "\n"
                      << "    Lattice vector z = (" << m_lattice_vectors[2][0] << ", " << m_lattice_vectors[2][1] << ", " << m_lattice_vectors[2][2] << ") " << m_coordinate_unit << "\n"
                      << "    Lattice origin = (" << m_lattice_origin[0] << ", " << m_lattice_origin[1] << ", " << m_lattice_origin[2] << ")\n"
                      << "    Lattice scaling = " << m_scaling << "\n"
                      << "    Lattice volume = " << this->volume() << " " << m_coordinate_unit << "^3\n"
                      << "    Grid = " << grid_sz() << " ["<< m_grid_dims[0] <<" x " << m_grid_dims[1] <<" x " << m_grid_dims[2] <<"]\n";
            print_time();
        }
        void print_time() const {
            std::cout << "    # time_steps = " << m_num_tsteps << ", time_step = " << m_time_step <<" " << m_time_unit << std::endl;
        }

        // -----------------------------------------------------------------------------
        // -----------------------------------------------------------------------------
        // conversion between various types of spatial coordinates!
            // direct :: in range [0, 1]
            // grid   :: in range [0, grid_dims]
            // world  :: in range [lattice_origins, lattice_origins+lattice_vectors]

        inline float direct_to_grid(float in_directc, unsigned int d = 0) const {   return in_directc * m_grid_dims[d];   }
        inline float grid_to_direct(float in_gridc, unsigned int d = 0) const {     return in_gridc / m_grid_dims[d];     }

        inline float direct_to_world(float in_directc, unsigned int d = 0) const {  return in_directc * m_lattice_vectors[d][d] + m_lattice_origin[d] ;  }
        inline float world_to_direct(float in_worldc, unsigned int d = 0) const {   return (in_worldc - m_lattice_origin[d]) / m_lattice_vectors[d][d];  }

        inline float grid_to_world(float in_gridc, unsigned int d = 0) const {      return direct_to_world( grid_to_direct(in_gridc, d), d );            }
        inline float world_to_grid(float in_worldc, unsigned int d = 0) const {     return direct_to_grid( world_to_direct(in_worldc, d), d );           }

        inline void direct_to_grid(const float in_directc[3], float out_gridc[3]) const {
            for(unsigned int d=0; d<3; d++)     out_gridc[d] = direct_to_grid(in_directc[d], d);
        }
        inline void grid_to_direct(const float in_gridc[3], float out_directc[3]) const {
            for(unsigned int d=0; d<3; d++)     out_directc[d] = grid_to_direct(in_gridc[d], d);
        }

        inline void direct_to_world(const float in_directc[3], float out_worldc[3]) const {
            for(unsigned int d=0; d<3; d++)     out_worldc[d] = direct_to_world(in_directc[d], d);
        }
        inline void world_to_direct(const float in_worldc[3], float out_directc[3]) const {
            for(unsigned int d=0; d<3; d++)     out_directc[d] = world_to_direct(in_worldc[d], d);
        }

        inline void grid_to_world(const float in_gridc[3], float out_worldc[3]) const {
            for(unsigned int d=0; d<3; d++)     out_worldc[d] = grid_to_world(in_gridc[d], d);
        }
        inline void world_to_grid(const float in_worldc[3], float out_gridc[3]) const {
            for(unsigned int d=0; d<3; d++)     out_gridc[d] = world_to_grid(in_worldc[d], d);
        }
        inline void world_to_grid(const float in_worldc[3], double out_gridc[3]) const {
            for(unsigned int d=0; d<3; d++)     out_gridc[d] = world_to_grid(in_worldc[d], d);
        }

        inline void world_dims(float worldc[3]) const {
            for(unsigned int d=0; d<3; d++)     worldc[d] = grid_to_world(m_grid_dims[d], d);
        }

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

        inline size_t grid_sz() const {     return m_grid_dims[0]*m_grid_dims[1]*m_grid_dims[2];    }

        inline double volume() const {
            return fabs(m_lattice_vectors[0][0]*(m_lattice_vectors[1][1]*m_lattice_vectors[2][2] - m_lattice_vectors[1][2]*m_lattice_vectors[2][1])
                      - m_lattice_vectors[0][1]*(m_lattice_vectors[1][0]*m_lattice_vectors[2][2] - m_lattice_vectors[1][2]*m_lattice_vectors[2][0])
                      + m_lattice_vectors[0][2]*(m_lattice_vectors[1][0]*m_lattice_vectors[2][1] - m_lattice_vectors[1][1]*m_lattice_vectors[2][1]));
        }

        // -----------------------------------------------------------------------------
        // -----------------------------------------------------------------------------
        // conversion between time and tidx

        inline double tidx_to_time(size_t tidx) const {     return double(tidx)*m_time_step;    }
        inline size_t time_to_tidx(double time) const {     return size_t(time/m_time_step);    }

        // -----------------------------------------------------------------------------
        // -----------------------------------------------------------------------------
    };
}
#endif
