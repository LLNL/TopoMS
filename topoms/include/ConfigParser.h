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
 *  @file    ConfigParser.h
 *  @author  Harsh Bhatia (hbhatia@llnl.gov)
 *  @date    10/01/2017
 *  @version 1.0
 *
 *  @brief This class handles the parsing of the input configuration file to TopoMS
 *
 *  @section DESCRIPTION
 *
 *  This class handles the parsing of the input configuration file to TopoMS.
 */

#ifndef _CONFIGPARSER_H
#define _CONFIGPARSER_H

#include <string>
#include <fstream>
#include <iostream>

#include "Utils.h"

class Config {

    const std::string configname;

public:

    std::string infiletype;
    std::string infilename;
    bool is_periodic[3];
    size_t grid_dims[3];

    double grad_threshold;
    double error_threshold;
    int numiter;
    int rkindex;

    bool bader;
    double vacuum_threshold;

    bool msc;
    double sval_threshold;
    double fval_threshold;

    std::string tfpath;   // optional

public:
    /**
      *   @brief  Initialize the constructor for a config file.
      *
      *   @param  fname is the name of the config file
      */
    Config(const std::string fname) : configname(fname),
        infilename(""), infiletype(""), tfpath(""),
        grad_threshold(-1), error_threshold(-1), numiter(1), rkindex(1),
        bader(false), vacuum_threshold(-1),
        msc(false), sval_threshold(-1), fval_threshold(-1)
    {
        for(int i = 0; i < 3; i++){
            is_periodic[i] = 0;
            grid_dims[i] = 0;
        }
    }

    /**
      *   @brief  Parses the config file specified to the constructor.
      *
      *   @return void
      */
    void parse() {

        bool need_newline = false;

        std::cout << " Parsing config file " << configname << "...";
        fflush(stdout);

        std::ifstream infile( configname.c_str() );
        std::string line;

        while (std::getline(infile, line)) {

            line = Utils::trim(line);
            if (line.empty() || line.at(0) == '#')
                continue;

            std::vector<std::string> toks = Utils::tokenize(line, '=');
            std::string &p = toks[0];
            std::string &v = toks[1];

            Utils::trim(p);
            Utils::trim(v);

            Utils::toupper(p);

            if(p == "INFILE"){            infilename = v;               continue;   }
            if(p == "INFILETYPE"){        infiletype = Utils::toupper(v);   continue;   }
            if(p == "THRESHOLD_GRAD"){    grad_threshold = stof(v);     continue;   }
            if(p == "THRESHOLD_ERROR"){   error_threshold = stof(v);    continue;   }
            if(p == "NUM_ITER"){          numiter = stoi(v);            continue;   }
            if(p == "RK_INDEX"){          rkindex = stoi(v);            continue;   }

            if(p == "GRID_DIMS"){
                sscanf(v.c_str(), "%d %d %d", &grid_dims[0], &grid_dims[1], &grid_dims[2]);
                continue;
            }
            if(p == "PERIODIC_DOMAIN"){
                sscanf(v.c_str(), "%d %d %d", &is_periodic[0], &is_periodic[1], &is_periodic[2]);
                continue;
            }

            if(p == "BADER_VOLUMES"){       bader = (Utils::toupper(v) == "TRUE");  continue;   }
            if(p == "THRESHOLD_VACUUM"){    vacuum_threshold = stof(v);             continue;   }
            if(p == "MOLECULAR_GRAPH"){     msc = (Utils::toupper(v) == "TRUE");    continue;   }
            if(p == "THRESHOLD_SIMPL"){     sval_threshold = stof(v);               continue;   }
            if(p == "THRESHOLD_FILT"){      fval_threshold = stof(v);               continue;   }

            if(p == "TRANSFER_FUNC"){       tfpath = v;                             continue;   }
        }

        // check if everytiong is in order
        if (infilename.length() == 0) {
            std::cerr << "\n\tConfig::parse() - Filename not specified. Aborting!\n";
            exit(1);
        }
        if (infiletype.length() == 0) {
            std::cerr << "\n\tConfig::parse() - Filetype not specified. Aborting!\n";
            exit(1);
        }

        if(msc){    bader = true;   }
        if(vacuum_threshold < 0){
            std::cout << "\n\tConfig::parse() - THRESHOLD_VACUUM not found. defaulting to 0.001";
            need_newline = true;
            vacuum_threshold = 0.1;
        }
        if(sval_threshold < 0){
            std::cout << "\n\tConfig::parse() - THRESHOLD_SIMPL not found. defaulting to 1.0";;
            need_newline = true;
            sval_threshold = 1.0;
        }
        if(fval_threshold < 0){
            std::cout << "\n\tConfig::parse() - THRESHOLD_FILT not found. defaulting to 1.0";
            need_newline = true;
            fval_threshold = 1.0;
        }
        if(need_newline)
            std::cout << std::endl;
        std::cout << " Done!\n";
    }
};

#endif
