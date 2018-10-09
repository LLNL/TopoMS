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
 *  @file    Utils.cpp
 *  @author  Harsh Bhatia (hbhatia@llnl.gov)
 *  @date    10/01/2017
 *
 *  @brief This file provides some basic utility functions
 *
 *  @section DESCRIPTION
 *
 *  This file provides some basic utility functions
 *
 */

#include <fstream>
#include <sstream>
#include <iostream>
#include <iterator>

#ifndef _WIN32
#include <unistd.h>
#include <sys/ioctl.h>
#endif

#include "Utils.h"

void Utils::print_separator(unsigned int n) {

#ifndef _WIN32
    struct winsize w;
    ioctl(STDOUT_FILENO, TIOCGWINSZ, &w);
    n = w.ws_col;
#else
    n = 40;
#endif
    for(int i=0; i<n; i++) printf("-");
    printf("\n");
}

// -----------------------------------------------------------
// -----------------------------------------------------------

std::string Utils::get_extn(std::string filename) {

    std::string::size_type idx = filename.rfind('.');

    if(idx != std::string::npos){
        return filename.substr(idx+1);
    }
    else
    {
        return "";
    }
}
std::string Utils::get_directory(std::string filename) {

    std::string::size_type idx = filename.rfind('/');

    if(idx != std::string::npos){
        return filename.substr(0, idx+1);
    }
    else
    {
        return "./";
    }
}
#if 0
std::string Utils::toupper(std::string& str) {
    for(uint32_t i=0; str[i]!=0; i++) {
        if(97 <= str[i] && str[i] <= 122){
            str[i]-=32;
        }
    }
    return str;
}

std::string Utils::trim(std::string& str) {
    str.erase(0, str.find_first_not_of(' '));       //prefixing spaces
    str.erase(str.find_last_not_of(' ')+1);         //surfixing spaces
    return str;
}

std::string Utils::remove_carriagereturn(std::string &str) {
    if (str[str.length()-1] == '\r')  str = str.erase(str.length()-1, 1);
}

std::string Utils::rtrim(std::string &str) {
    //remove_carriagereturn(str);
    str.erase(std::find_if(str.rbegin(), str.rend(), [](int ch) {
        return !std::isspace(ch);
    }).base(), str.end());
}
#endif
std::vector<std::string> Utils::tokenize(const std::string &line, char delim) {

    std::vector<std::string> tokens;

    std::stringstream ss;
    ss.str(line);

    std::string item;
    while (std::getline(ss, item, delim)) {
        tokens.push_back( item );
    }
    return tokens;
}

std::vector<std::string> Utils::tokenize(std::string line){

    // construct a stream from the string
    std::stringstream linestream(line);

    // use stream iterators to copy the stream to the vector as whitespace separated strings
    std::istream_iterator<std::string> it_line(linestream);
    std::istream_iterator<std::string> end_line;
    std::vector<std::string> tokens(it_line, end_line);
    return tokens;
}
