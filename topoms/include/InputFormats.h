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
 *  @file    InputFormats.h
 *  @author  Harsh Bhatia (hbhatia@llnl.gov)
 *  @date    10/01/2017
 *
 *  @brief This header provides I/O functionality for different MD file formats
 *
 *  @section DESCRIPTION
 *
 *  This header provides I/O functionality for different MD file formats
 *
 */

#ifndef _INPUT_FORMATS_H_
#define _INPUT_FORMATS_H_

#include <cmath>
#include <vector>
#include <string>
#include <cctype>
#include <fstream>
#include <iostream>

#include "Utils.h"
#include "Material.h"
#include "MolecularSystem.h"

namespace MS {

    namespace VASP {

        /**
          *   @brief  Read a VASP CAR file: CHGCAR, AECCAR, LOCPOT
          *
          *   @param  filename is the input filename
          *   @param  mdata is the metadata about the system
          *   @param  atoms is the list of atoms in the system
          *
          *   @return array of function values
          */
        template <typename T>
        T* read_CAR(const std::string &filename, MS::SystemInfo &mdata, std::vector<Material> &atoms) {

            std::ifstream infile(filename.c_str());
            if(!infile.is_open()){
                std::cerr << " MD::VASP::read_CAR -- Could not open file " << filename << std::endl;
                exit(1);
            }

            // -----------------------------------------------------
            std::cout << " Reading VASP CAR file (" << filename << ")...";
            fflush(stdout);

            std::string line;

            mdata.m_length_unit = "Ang";   mdata.m_length_file2Angs = 1.0;
            mdata.m_charge_unit = "eV";    mdata.m_charge_file2electrons = 1.0;

            // -----------------------------------------------------
            // CAR header

            // read lattice
            std::getline(infile, mdata.m_sysname);
            Utils::rtrim(mdata.m_sysname);

            std::getline(infile, line);
            sscanf(line.c_str(), "%lf", &mdata.m_scaling);

            for(uint8_t i = 0; i < 3; i++){
                std::getline(infile, line);
                sscanf(line.c_str(),"%lf %lf %lf", &mdata.m_lattice.v[i][0], &mdata.m_lattice.v[i][1], &mdata.m_lattice.v[i][2]);

                mdata.m_lattice.v[i][0] *= mdata.m_scaling;
                mdata.m_lattice.v[i][1] *= mdata.m_scaling;
                mdata.m_lattice.v[i][2] *= mdata.m_scaling;
            }

            mdata.invert_lattice();

            // read materials
            std::getline(infile, line);
            std::vector<std::string> tokens = Utils::tokenize(line);

            size_t num_atoms = 0;
            mdata.m_materials.resize(tokens.size(), std::pair<std::string, unsigned int> ("", 0));

            // check if the files provided symbols!
            if( !std::isdigit(tokens.front().at(0))) {
                for(size_t i = 0; i < tokens.size(); i++){
                    mdata.m_materials[i].first = tokens[i];
                }

                std::getline(infile, line);
                tokens = Utils::tokenize(line);
            }

            // read the counts
            for(size_t i = 0; i < tokens.size(); i++){
                mdata.m_materials[i].second = atoi(tokens[i].c_str());
                num_atoms += mdata.m_materials[i].second;
            }

            // -----------------------------------------------------
            // read positions
            std::getline(infile, mdata.m_coordinate_type);     // read "Direct"
            Utils::rtrim(mdata.m_coordinate_type);

            atoms.reserve(num_atoms);

            double x, y, z;
            size_t prev = 0;
            for(size_t i = 0; i < mdata.m_materials.size(); i++) {

                size_t curr = mdata.m_materials[i].second;
                std::string symb = mdata.m_materials[i].first;

                for(size_t j = prev; j < curr; j++) {

                    std::getline(infile, line);
                    sscanf(line.c_str(), "%lf %lf %lf", &x, &y, &z);

                    atoms.push_back( Material(symb, x, y, z) );
                }
                curr = prev;
            }

            if(mdata.m_coordinate_type.compare("Direct") == 0){
                for(size_t i = 0; i < atoms.size(); i++){
                    mdata.direct_to_world(atoms[i].m_pos, atoms[i].m_pos);
                }
                mdata.m_coordinate_type = std::string("World");
            }

            // read grid dimensions
            while(true) {
               std::getline(infile, line);
                Utils::rtrim(line);
                if(line.length() >= 2)
                    break;
            }

            sscanf(line.c_str(),"%d %d %d", &mdata.m_grid_dims[0], &mdata.m_grid_dims[1], &mdata.m_grid_dims[2]);

            // VASP write 5 values in each line, and rounds off the last line with 0
            const size_t nvalues = mdata.grid_sz();

            // read values
            size_t idx = 0;
            T *values = new T[nvalues];

            while(true) {

                if (idx == nvalues) break;

                std::getline(infile, line);
                if(line.empty())    break;

                std::vector<std::string> tokens = Utils::tokenize(line);
                for(size_t i = 0; i < tokens.size(); i++){
                    values[idx++] = atof(tokens[i].c_str());
                    if (idx == nvalues) break;
                }
            }

            infile.close();
            std::cout <<" Done!\n";
            return values;
        }

        /**
          *   @brief  Write a VASP CAR file: CHGCAR, AECCAR, LOCPOT
          *
          *   @param  filename is the input filename
          *   @param  mdata is the metadata about the system
          *   @param  atoms is the list of atoms in the system
          *
          *   @return success flag
          */
        template <typename T>
        bool write_CAR(const std::string &filename, const MS::SystemInfo &mdata, const std::vector<Material> &atoms, const T *vals = nullptr) {

            FILE *outfile = fopen(filename.c_str(), "w");
            if(outfile == NULL){
                std::cerr <<" Could not open file " << filename << std::endl;
                return false;
            }

            std::cout << " - Writing " << filename << "...";

            // write lattice
            fprintf(outfile, "%s\n", mdata.m_sysname.c_str());
            fprintf(outfile, "%16f\n", 1.0);

            for(uint8_t i = 0; i < 3; i++){
                fprintf(outfile, "%14.6f %14.6f %14.6f\n", mdata.m_lattice.v[i][0], mdata.m_lattice.v[i][1], mdata.m_lattice.v[i][2]);
            }

            // write materials
            for(size_t i = 0; i < mdata.m_materials.size(); i++){
                fprintf(outfile, "%6s", mdata.m_materials[i].first.c_str());
            }
            fprintf(outfile, "\n");
            for(size_t i = 0; i < mdata.m_materials.size(); i++){
                fprintf(outfile, "%7d", mdata.m_materials[i].second);
            }
            fprintf(outfile, "\n");

            // write positions
            fprintf(outfile, "%s\n", mdata.m_coordinate_type.c_str());
            for(size_t i = 0; i < atoms.size(); i++){
                fprintf(outfile, "%13.7f %13.7f %13.7f\n", atoms[i].m_pos[0], atoms[i].m_pos[1], atoms[i].m_pos[2]);
            }

            // write grid
            fprintf(outfile, "\n%d %d %d\n", mdata.m_grid_dims[0], mdata.m_grid_dims[1], mdata.m_grid_dims[2]);

            if(vals != nullptr){
                for(size_t i = 0; i < mdata.grid_sz(); i+=5){
                    fprintf(outfile, "%E %E %E %E %E\n", vals[i], vals[i+1], vals[i+2], vals[i+3], vals[i+4]);
                }
            }

            printf(" Done!\n");
            fclose(outfile);
        }
    }

    namespace Cube {

        /**
          *   @brief  Read a Cube file
          *
          *   @param  filename is the input filename
          *   @param  mdata is the metadata about the system
          *   @param  atoms is the list of atoms in the system
          *
          *   @return array of function values
          */
        template <typename T>
        T* read(const std::string &filename, MS::SystemInfo &mdata, std::vector<Material> &atoms) {

            std::ifstream infile(filename.c_str());
            if(!infile.is_open()){
                std::cerr << " MD::Cube::read -- Could not open file " << filename << std::endl;
                exit(1);
            }

            // -----------------------------------------------------
            std::cout << " Reading Cube file " << filename << "...";
            fflush(stdout);

            std::string line;

            std::getline(infile, mdata.m_sysname);      // use first line as system name
            std::getline(infile, line);                 // ignore second line

            // number of atoms and origin of the volume
            int num_atoms;
            std::getline(infile, line);
            sscanf(line.c_str(),"%d %lf %lf %lf", &num_atoms, &mdata.m_lattice_origin[0], &mdata.m_lattice_origin[1], &mdata.m_lattice_origin[2]);

            int n;
            double dx, dy, dz;

            // read gridsize and lattice vectors
            for(uint8_t i = 0; i < 3; i++){

                std::getline(infile, line);
                sscanf(line.c_str(),"%d %lf %lf %lf", &n, &dx, &dy, &dz);

                // these numbers are N, dx, dy, dz in each dimension
                mdata.m_grid_dims[i] = abs(n);

                // for i = 0, N determines the unit of charge and spatial coordinates
                if (i == 0) {
                    // if N is positive, then the lattice vector is on Bohr
                    // if N is negative, then the lattice vector is in Angs
                    if (n > 0) {
                        mdata.m_length_unit = "Bohr";       mdata.m_length_file2Angs = 1.889725992;
                        mdata.m_charge_unit = "hartree";    mdata.m_charge_file2electrons = 0.036749309;
                    }
                    else {
                        mdata.m_length_unit = "Ang";        mdata.m_length_file2Angs = 1.0;
                        mdata.m_charge_unit = "eV";         mdata.m_charge_file2electrons = 1.0;
                    }
                }

                mdata.m_lattice.v[i][0] = (double)mdata.m_grid_dims[i] * dx;
                mdata.m_lattice.v[i][1] = (double)mdata.m_grid_dims[i] * dy;
                mdata.m_lattice.v[i][2] = (double)mdata.m_grid_dims[i] * dz;
            }

            mdata.m_coordinate_type = "World";
            mdata.invert_lattice();

            // now read atom locations
            atoms.reserve(num_atoms);

            double x, y, z;
            double chg;
            for(unsigned int i = 0; i < num_atoms; i++) {

                std::getline(infile, line);
                sscanf(line.c_str(),"%d %lf %lf %lf %lf", &n, &chg, &x, &y, &z);

                if(x >= (mdata.m_lattice_origin[0] + mdata.m_lattice.v[0][0]))
                    x -= mdata.m_lattice.v[0][0];

                if(y >= (mdata.m_lattice_origin[1] + mdata.m_lattice.v[1][1]))
                    y -= mdata.m_lattice.v[1][1];

                if(z > (mdata.m_lattice_origin[2] + mdata.m_lattice.v[2][2])){
                    z -= mdata.m_lattice.v[2][2];
                }
                atoms.push_back(Material(n, chg, x, y, z));
            }

            // read values
            // the ordering of values in CUBE files is fortran order
                // http://paulbourke.net/dataformats/cube/
            size_t cnt = 0;
            T *values = new T[mdata.grid_sz()];

            const size_t X = mdata.m_grid_dims[0];
            const size_t XY = mdata.m_grid_dims[0]*mdata.m_grid_dims[1];

            while(cnt < mdata.grid_sz()){

                std::getline(infile, line);
                if(line.empty()){
                    break;
                }

                std::vector<std::string> tokens = Utils::tokenize(line);
                for(size_t i = 0; i < tokens.size(); i++, cnt++){

                    // convert from column major to row major
                    size_t z = cnt / XY;
                    size_t x = (cnt % XY) % X;
                    size_t idx = cnt + int(XY - 1)*int(x - z);

                    values[idx] = atof(tokens[i].c_str());
                }
            }
            infile.close();
            std::cout <<" Done!\n";
            return values;
        }
    }
}
#endif
