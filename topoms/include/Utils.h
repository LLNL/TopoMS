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
 *  @file    Utils.h
 *  @author  Harsh Bhatia (hbhatia@llnl.gov)
 *  @date    10/01/2017
 *
 *  @brief This header provides some basic utility functions
 *
 *  @section DESCRIPTION
 *
 *  This header provides some basic utility functions
 *
 */

#ifndef _UTILS_H_
#define _UTILS_H_

#include <cmath>
#include <numeric>
#include <algorithm>

#include <map>
#include <vector>
#include <string>

#include <cctype>

namespace Utils {

    // -------------------------------------------------------
    // basic mathematical utilities

    inline bool equals(double a, double b){     return (fabs(a-b) < powf(10,-5));   }

    inline float dist_sq(const float p[3], const float q[3]){
        return ((p[0]-q[0])*(p[0]-q[0]) + (p[1]-q[1])*(p[1]-q[1]) + (p[2]-q[2])*(p[2]-q[2]));
    }

    inline float dist(const float p[3], const float q[3]){  return std::sqrt(dist_sq(p, q));            }

    inline float magn_sq(const float a[3]){                 return (a[0]*a[0] + a[1]*a[1] + a[2]*a[2]); }

    inline float magn(const float a[3]){                    return std::sqrt(magn_sq(a));               }

    inline float dot(const float a[3], const float b[3]){   return (a[0]*b[0] + a[1]*b[1] + a[2]*b[2]); }

    inline float* cross(const float a[3], const float b[3]){
        float *cp = new float[3];
        cp[0] = a[1]*b[2] - a[2]*b[1];
        cp[1] = a[2]*b[0] - a[0]*b[2];
        cp[2] = a[0]*b[1] - a[1]*b[0];
        return cp;
    }

    template <typename T>
    inline T lerp(const T &a, const T &b, float l) {
        return (a + (b-a)*l);
    }

    template <typename T1, typename T2>
    inline T2 lerp (const T1 &x, const T1 &x0, const T1 &x1, const T2 &y0, const T2 &y1) {

        float l = (x0 == x1) ? 0.0f : (float) (x-x0) / (float) (x1-x0);
        return lerp (y0, y1, l);
    }

    template <typename T>
    inline T lerp2(const T &a, const T &b, const T &c, const T &d,
                   float l1, float l2) {
        return lerp( lerp(a,b,l1), lerp(c,d,l1), l2 );
    }

    template <typename T>
    inline T lerp3(const T &a, const T &b, const T &c, const T &d,
                   const T &e, const T &f, const T &g, const T &h,
                   float l1, float l2, float l3) {
        return lerp( lerp2(a, b, c, d, l1, l2), lerp2(e, f, g, h, l1, l2), l3 );
    }

    template <typename T>
    bool lies_value_in_range(T v, T a, T b){
        return !( v < a || v > b );
    }

    // -------------------------------------------------------------------------
    // -------------------------------------------------------------------------
    // Kahan summation

    namespace Kahan {

        template <typename T>
        struct KahanObject {
            T sum;
            T correction;
        };

        template <typename T>
        KahanObject<T> KahanSum(KahanObject<T> accumulation, T value) {
            KahanObject<T> result;
            T y = value - accumulation.correction;
            T t = accumulation.sum + y;
            result.correction = (t - accumulation.sum) - y;
            result.sum = t;
            return result;
        }

        template <typename T>
        T sum(const std::vector<T> &values) {
            KahanObject<T> k0 = {0};
            return std::accumulate(values.begin(), values.end(), k0, KahanSum<T>).sum;
        }
        template <typename T>
        T sum(const T *values, size_t sz) {
            KahanObject<T> k0 = {0};
            return std::accumulate(values, values+sz, k0, KahanSum<T>).sum;
        }
    }

    template <typename T>
    T sum(const std::vector<T> &values) {
        return std::accumulate(values.begin(), values.end(), T(0));
    }
    template <typename T>
    T sum(const T *values, size_t sz) {
        return std::accumulate(values, values+sz, T(0));
    }

    // -------------------------------------------------------------------------
    // -------------------------------------------------------------------------

    template <typename T>
    typename std::vector<T>::const_iterator get_lowerIterator(const std::vector<T> &vals, const T& value){

        // upper_bound returns an iterator pointing to the first element that is
        // greater than the value.
        typename std::vector<T>::const_iterator geq = std::upper_bound(vals.begin(), vals.end(), value);
        if(geq == vals.begin())     return geq;
        return --geq;
    }

    template <typename T>
    std::size_t get_lowerIndex(const std::vector<T> &vals, const T& value){
        return std::distance(vals.begin(), get_lowerIterator(vals, value));
    }

    template <typename T>
    inline std::size_t snap_to_axis(T val, const std::vector<T> &values){
        return get_closestIndex(values, val);
    }

    // -------------------------------------------------------------------------
    // handling periodic boundary


    // wrap around in [0,1]
    inline float wrap_around(float val, float min, float max){
        return (val < min) ? (val + (max-min)) :
               (val > max) ? (val - (max-min)) : val;
    }

    inline float get_periodicDisplacement(float p, float q, float span) {
        float d = p-q;
        return (d > 0.5f*span)  ? (d - span) :     // q is on the left end, and p is on the right end
               (d < -0.5f*span) ? (d + span) :     // p is on the left end, and q is on the right end
                                 d;                // else!
    }

    inline float get_periodicDist_sq(const float p[3], const float q[3], const float span[3]) {
        float r[3];
        for(unsigned int i = 0; i < 3; i++){
            r[i] = get_periodicDisplacement(p[i], q[i], span[i]);
        }
        return (r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);
    }

    inline float get_periodicDist(const float p[3], const float q[3], const float span[3]) {
        return std::sqrt(get_periodicDist_sq(p,q,span));
    }

    // -------------------------------------------------------------------------
    // IO related
    std::string get_extn(std::string filename);
    std::string get_directory(std::string filename);

    // read the files in a given directory
    std::vector<std::string> get_directory_listing(std::string searchstring);

    // -------------------------------------------------------------------------
    // strung utilities
    inline std::string remove_carriagereturn(std::string& str) {

        if (str.length() == 0)              return str;
        if (str[str.length()-1] == '\r')    str = str.erase(str.length()-1, 1);
        return str;
    }

    inline std::string trim(std::string& str) {
        str.erase(0, str.find_first_not_of(' '));       //prefixing spaces
        str.erase(str.find_last_not_of(' ')+1);         //surfixing spaces
        return str;
    }

    inline std::string rtrim(std::string& str) {
        str.erase(std::find_if(str.rbegin(), str.rend(), [](int ch) {
            return !std::isspace(static_cast<char>(ch)); 
        }).base(), str.end());
        return str;
    }

    inline std::string toupper(std::string& str) {
        for(uint32_t i=0; str[i]!=0; i++) {
            if(97 <= str[i] && str[i] <= 122){
                str[i]-=32;
            }
        }
        return str;
    }

    // tokenize a string
    std::vector<std::string> tokenize(std::string line);
    std::vector<std::string> tokenize(const std::string &line, char delim);


    // for formatting output
    void print_separator(unsigned int n = 120);
}








#endif
