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
 *  @file    Material.h
 *  @author  Harsh Bhatia (hbhatia@llnl.gov)
 *  @date    10/01/2017
 *
 *  @brief This class handles the properties of different types of materials
 *
 *  @section DESCRIPTION
 *
 *  This class handles the properties of different types of materials
 *
 */

#ifndef _MATERIAL_H_
#define _MATERIAL_H_

#include <cmath>
#include <vector>
#include <string>
#include <iostream>

class Material {

public:
    std::string m_symbol;
    unsigned int m_atom_number;
    float m_charge;
    float *m_pos;

    Material(unsigned int n, float chg, float x, float y, float z) {
        m_symbol = "";
        m_atom_number = n;
        m_charge = chg;
        m_pos = new float[3];
        m_pos[0] = x;   m_pos[1] = y;   m_pos[2] = z;
    }

    Material(std::string s, float x, float y, float z) {
        m_symbol = s;
        m_atom_number = 0;
        m_charge = -1;
        m_pos = new float[3];
        m_pos[0] = x;   m_pos[1] = y;   m_pos[2] = z;
    }

    float nelectrons() const {      return (float) m_atom_number;   }

    float radius() const {

        // get radius in Angstroms in Angstrom
        // http://www.periodictable.com/Properties/A/AtomicRadius.v.wt.html

        if(m_symbol == "F")             return 0.42;
        if(m_symbol == "O")             return 0.48;
        if(m_symbol == "H")             return 0.53;
        if(m_symbol == "C")             return 0.67;
        if(m_symbol == "Cl")            return 0.79;
        if(m_symbol == "B")             return 0.87;
        if(m_symbol == "P")             return 0.98;
        if(m_symbol == "Li")            return 1.67;
        if(m_symbol == "Ti")            return 1.76;
        if(m_symbol == "Ca")            return 1.94;

        return 1.0;
    }

    bool operator == (const Material &m) const {

        static const float eps = 0.00001;

        if (m_symbol.compare(m.m_symbol) != 0)      return false;
        if (m_atom_number != m.m_atom_number)       return false;
        if (fabs(m_charge - m.m_charge) > eps)      return false;

        if (fabs(m_pos[0] - m.m_pos[0]) > eps)      return false;
        if (fabs(m_pos[1] - m.m_pos[1]) > eps)      return false;
        if (fabs(m_pos[2] - m.m_pos[2]) > eps)      return false;
        return true;
    }

    void print(int id = -1) const {

        std::cout << " Atom ";
        if(id != -1) {              std::cout << id << ": ";                          }
        if(m_symbol != "") {        std::cout << m_symbol << " ";                     }
        if(m_atom_number != 0) {    std::cout << "[" << int(m_atom_number) << "] ";   }
        std::cout << "at (" << m_pos[0]<<", "<<m_pos[1]<<", "<<m_pos[2]<<")" << std::endl;
    }
};

#endif
