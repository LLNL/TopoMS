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

    static std::string atomic_number_to_symbol(const unsigned int &_) {
        switch(_) {
        case   1  : return   "H";
        case   2  : return   "He";
        case   3  : return   "Li";
        case   4  : return   "Be";
        case   5  : return   "B";
        case   6  : return   "C";
        case   7  : return   "N";
        case   8  : return   "O";
        case   9  : return   "F";
        case   10  : return   "Ne";
        case   11  : return   "Na";
        case   12  : return   "Mg";
        case   13  : return   "Al";
        case   14  : return   "Si";
        case   15  : return   "P";
        case   16  : return   "S";
        case   17  : return   "Cl";
        case   18  : return   "Ar";
        case   19  : return   "K";
        case   20  : return   "Ca";
        case   21  : return   "Sc";
        case   22  : return   "Ti";
        case   23  : return   "V";
        case   24  : return   "Cr";
        case   25  : return   "Mn";
        case   26  : return   "Fe";
        case   27  : return   "Co";
        case   28  : return   "Ni";
        case   29  : return   "Cu";
        case   30  : return   "Zn";
        case   31  : return   "Ga";
        case   32  : return   "Ge";
        case   33  : return   "As";
        case   34  : return   "Se";
        case   35  : return   "Br";
        case   36  : return   "Kr";
        case   37  : return   "Rb";
        case   38  : return   "Sr";
        case   39  : return   "Y";
        case   40  : return   "Zr";
        case   41  : return   "Nb";
        case   42  : return   "Mo";
        case   43  : return   "Tc";
        case   44  : return   "Ru";
        case   45  : return   "Rh";
        case   46  : return   "Pd";
        case   47  : return   "Ag";
        case   48  : return   "Cd";
        case   49  : return   "In";
        case   50  : return   "Sn";
        case   51  : return   "Sb";
        case   52  : return   "Te";
        case   53  : return   "I";
        case   54  : return   "Xe";
        case   55  : return   "Cs";
        case   56  : return   "Ba";
        case   57  : return   "La";
        case   58  : return   "Ce";
        case   59  : return   "Pr";
        case   60  : return   "Nd";
        case   61  : return   "Pm";
        case   62  : return   "Sm";
        case   63  : return   "Eu";
        case   64  : return   "Gd";
        case   65  : return   "Tb";
        case   66  : return   "Dy";
        case   67  : return   "Ho";
        case   68  : return   "Er";
        case   69  : return   "Tm";
        case   70  : return   "Yb";
        case   71  : return   "Lu";
        case   72  : return   "Hf";
        case   73  : return   "Ta";
        case   74  : return   "W";
        case   75  : return   "Re";
        case   76  : return   "Os";
        case   77  : return   "Ir";
        case   78  : return   "Pt";
        case   79  : return   "Au";
        case   80  : return   "Hg";
        case   81  : return   "Tl";
        case   82  : return   "Pb";
        case   83  : return   "Bi";
        case   84  : return   "Po";
        case   85  : return   "At";
        case   86  : return   "Rn";
        case   87  : return   "Fr";
        case   88  : return   "Ra";
        case   89  : return   "Ac";
        case   90  : return   "Th";
        case   91  : return   "Pa";
        case   92  : return   "U";
        case   93  : return   "Np";
        case   94  : return   "Pu";
        case   95  : return   "Am";
        case   96  : return   "Cm";
        case   97  : return   "Bk";
        case   98  : return   "Cf";
        case   99  : return   "Es";
        case   100  : return   "Fm";
        case   101  : return   "Md";
        case   102  : return   "No";
        case   103  : return   "Lr";
        case   104  : return   "Rf";
        case   105  : return   "Db";
        case   106  : return   "Sg";
        case   107  : return   "Bh";
        case   108  : return   "Hs";
        case   109  : return   "Mt";
        case   110  : return   "Ds";
        case   111  : return   "Rg";
        case   112  : return   "Uub";
        case   113  : return   "Uut";
        case   114  : return   "Uuq";
        case   115  : return   "Uup";
        case   116  : return   "Uuh";
        case   117  : return   "Uus";
        case   118  : return   "Uuo";
        }
        return "unknown";
    }
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
        m_symbol = atomic_number_to_symbol(n);
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
