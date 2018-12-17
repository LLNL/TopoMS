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
 *
 *  @brief This class handles the parsing of the input configuration file to TopoMS
 *
 *  @section DESCRIPTION
 *
 *  This class handles the parsing of the input configuration file to TopoMS.
 */

#ifndef _CONFIGPARSER_H_
#define _CONFIGPARSER_H_

#include <string>
#include <fstream>
#include <iostream>

#include "Utils.h"

class Config {

    const std::string configname;

public:
    std::string infiletype;     // VASP or CUBE
    std::string fieldtype;      // CHG or POT

    std::string infilename;
    std::string reffilename;

    bool is_periodic[3];

    bool do_bader;               // parameters for bader anlaysis
    double threshold_vacuum;

    bool do_msc;                 // parameters for msc
    double threshold_simp;
    double threshold_filt;

    double threshold_grad;      // other internal parameters
    double threshold_error;
    int numiter;
    int rkindex;

    std::string tfpath;         // optional

public:
    /**
      *   @brief  Initialize the constructor for a config file.
      *
      *   @param  fname is the name of the config file
      */
    Config(const std::string fname) : configname(fname),
        infilename(""), reffilename(""), infiletype(""), fieldtype(""), tfpath(""),
        do_bader(false), threshold_vacuum(-1),
        do_msc(false), threshold_simp(-1), threshold_filt(-1),
        threshold_grad(-1), threshold_error(-1), numiter(1), rkindex(1)
    {
        for(uint8_t i = 0; i < 3; i++){
            is_periodic[i] = 0;
            //grid_dims[i] = 0;
        }
    }

    /**
      *   @brief  Parses the config file specified to the constructor.
      *
      *   @return void
      */
    void parse() {

        bool need_newline = false;

        std::cout << " Parsing config file (" << configname << ")...";
        fflush(stdout);

        std::ifstream infile( configname.c_str() );
        std::string line;

        // read all known parameters
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

            if(p == "INFILE"){            infilename = v;                           continue;   }
            if(p == "REFCHG"){            reffilename = v;                          continue;   }
            if(p == "INFILETYPE"){        infiletype = Utils::toupper(v);           continue;   }
            if(p == "FIELDTYPE"){         fieldtype = Utils::toupper(v);            continue;   }

            if(p == "PERIODIC_DOMAIN"){
                sscanf(v.c_str(), "%d %d %d", &is_periodic[0], &is_periodic[1], &is_periodic[2]);
                continue;
            }

            if(p == "BADER_VOLUMES"){     do_bader = (Utils::toupper(v) == "TRUE"); continue;   }
            if(p == "THRESHOLD_VACUUM"){  threshold_vacuum = stof(v);               continue;   }

            if(p == "MOLECULAR_GRAPH"){   do_msc = (Utils::toupper(v) == "TRUE");   continue;   }
            if(p == "THRESHOLD_SIMPL"){   threshold_simp = stof(v);                 continue;   }
            if(p == "THRESHOLD_FILT"){    threshold_filt = stof(v);                 continue;   }

            if(p == "THRESHOLD_GRAD"){    threshold_grad = stof(v);                 continue;   }
            if(p == "THRESHOLD_ERROR"){   threshold_error = stof(v);                continue;   }
            if(p == "NUM_ITER"){          numiter = stoi(v);                        continue;   }
            if(p == "RK_INDEX"){          rkindex = stoi(v);                        continue;   }

            if(p == "TRANSFER_FUNC"){     tfpath = v;                             continue;   }

            std::cerr << "\n    Config::parse() - Ignoring unknown parameter: " << line;
            need_newline = true;
        }

        // check if everything is in order
        if (infilename.length() == 0) {
            std::cerr << "\n    Config::parse() - Filename not specified. Aborting!\n";
            exit(1);
        }
        if (infiletype.length() == 0) {
            std::cerr << "\n    Config::parse() - Filetype not specified. Aborting!\n";
            exit(1);
        }
        if (fieldtype.length() == 0) {
            std::cerr << "\n    Config::parse() - Fieldtype not specified. Aborting!\n";
            exit(1);
        }
        if (fieldtype != "CHG" && fieldtype != "POT" && fieldtype != "OTHER") {
            std::cerr << "\n    Config::parse() - Invalid Fieldtype ("<<fieldtype<<"). Expected CHG or POT or OTHER. Aborting!\n";
            exit(1);
        }

        if (fieldtype != "CHG" && reffilename.length() > 0) {
            reffilename = "";
            std::cerr << "\n    Config::parse() - Ignoring REFCHG because FIELDTYPE != CHG";
            need_newline = true;
        }
        if(threshold_vacuum < 0){
            threshold_vacuum = 0.1;
            std::cout << "\n    Config::parse() - THRESHOLD_VACUUM not found. defaulting to " << threshold_vacuum;
            need_newline = true;
        }
        if(threshold_simp < 0){
            threshold_simp = 1.0;
            std::cout << "\n    Config::parse() - THRESHOLD_SIMPL not found. defaulting to " << threshold_simp;
            need_newline = true;
        }
        if(threshold_filt < 0){
            threshold_filt = 1.0;
            std::cout << "\n    Config::parse() - THRESHOLD_FILT not found. defaulting to " << threshold_filt;
            need_newline = true;
        }
        if(need_newline)
            std::cout << std::endl;
        std::cout << " Done!\n";
    }
};

#endif
