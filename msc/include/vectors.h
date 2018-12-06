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

#ifndef VECTORS_H
#define VECTORS_H


#include <math.h>
#include <stdio.h>
#include <string.h>

#include "basic_types.h"



namespace MSC {


    template <typename DTYPE>
    class GenericVec2 {
    protected:
        DTYPE m_v[2]; //x,y,z;

    public:

        GenericVec2() {}
        GenericVec2(DTYPE x, DTYPE y) {
            m_v[0] = x; m_v[1] = y;
        }

        void operator=(const GenericVec2& other){
            this->m_v[0] = other.m_v[0];
            this->m_v[1] = other.m_v[1];
        }
        void operator+=(const GenericVec2& other) {
            this->m_v[0] += other.m_v[0];
            this->m_v[1] += other.m_v[1];
        }
        void operator-=(const GenericVec2& other) {
            this->m_v[0] -= other.m_v[0];
            this->m_v[1] -= other.m_v[1];
        }

        void operator-=(const DTYPE other) {
            this->m_v[0] -= other;
            this->m_v[1] -= other;
        }
        void operator*=(const DTYPE other) {
            this->m_v[0] *= other;
            this->m_v[1] *= other;
        }

        DTYPE MagSq() const {
            return (m_v[0] * m_v[0]) + (m_v[1] * m_v[1]) ;
        }

        DTYPE Mag() const {
            return sqrt(MagSq());
        }
        GenericVec2 operator+(const GenericVec2  &a) const {
            GenericVec2 res(*this);
            res += a;
            return res;
        }
        GenericVec2 operator-(const GenericVec2  &a) const {
            GenericVec2 res(*this);
            res -= a;
            return res;
        }

        GenericVec2 operator+(const DTYPE  val) const {
            GenericVec2 res(*this);
            res.m_v[0] += val; res.m_v[1] += val;
            return res;
        }
        GenericVec2 operator-(const DTYPE  val) const {
            GenericVec2 res(*this);
            res.m_v[0] -= val; res.m_v[1] -= val;
            return res;
        }

        GenericVec2 operator*(const DTYPE  val) const {
            GenericVec2 res(this->m_v[0], this->m_v[1]);// = (*this);
            res.m_v[0] *= val; res.m_v[1] *= val;
            return res;
        }
        bool operator==(const GenericVec2& other) const {
            return memcmp(this, &other, sizeof(GenericVec2)) == 0;
        }

        operator GenericVec2<int>() const { return GenericVec2<int>(this->m_v[0], this->m_v[1]); }
        operator GenericVec2<long long>() const { return GenericVec2<long long>(this->m_v[0], this->m_v[1]); }
        operator GenericVec2<double>() const { return GenericVec2<double>(this->m_v[0], this->m_v[1]); }

        // THIS GIVES REMAINDER, NOT MODULO
        GenericVec2 operator%(const GenericVec2& other) const{
            GenericVec2 res;
            res.m_v[0] = this->m_v[0] % other.m_v[0];
            res.m_v[1] = this->m_v[1] % other.m_v[1];
            return res;
        }

        DTYPE& operator[](int i) { return m_v[i]; }
        const DTYPE& operator[](int i) const { return m_v[i]; }

        void PrintInt() const {
            printf("(%d, %d)\n", m_v[0], m_v[1]);
        }
        void PrintFloat() const {
            printf("(%f, %f)\n", m_v[0], m_v[1]);
        }

        // this is a hack that rounds towards -infinity for "small" negative values
        GenericVec2<INDEX_TYPE> IntFloor() const {
            //(*this).print_vf();
            GenericVec2 dres = *this;
            //dres = dres + 1000.0;
            //dres.print_vf();
            GenericVec2<INDEX_TYPE> ires(floor(dres[0]), floor(dres[1]));
            //ires.print_vi();
            //ires = ires + -1000;
            //ires.print_vi();
            return ires;
        }
        GenericVec2<double> DoubleFloor() const {
            GenericVec2<double> ddres(floor(this->m_v[0]), floor(this->m_v[1]));
            return ddres;
        }
        // linear interpolation betwee 2 vectors: returns the vector a * (1-t) + b * t
        static GenericVec2<double> Lerp(const GenericVec2<double>& a, const GenericVec2<double>& b, const double t)  {
            return (a * (1.0 - t)) + (b * (t));
        }

    };


    typedef GenericVec2<int>			Vec2i;
    typedef GenericVec2<INDEX_TYPE>		Vec2l;
    typedef GenericVec2<double>			Vec2d;
    typedef GenericVec2<bool>			Vec2b;

    template <typename DTYPE>
    class GenericVec3 {
    protected:
        DTYPE m_v[3]; //x,y,z;

    public:

        GenericVec3() {}
        GenericVec3(DTYPE x, DTYPE y, DTYPE z) {
            m_v[0] = x; m_v[1] = y; m_v[2] = z;
        }

        void operator=(const GenericVec3& other){
            this->m_v[0] = other.m_v[0];
            this->m_v[1] = other.m_v[1];
            this->m_v[2] = other.m_v[2];
        }
        void operator+=(const GenericVec3& other) {
            this->m_v[0] += other.m_v[0];
            this->m_v[1] += other.m_v[1];
            this->m_v[2] += other.m_v[2];
        }
        void operator-=(const GenericVec3& other) {
            this->m_v[0] -= other.m_v[0];
            this->m_v[1] -= other.m_v[1];
            this->m_v[2] -= other.m_v[2];
        }

        void operator-=(const DTYPE other) {
            this->m_v[0] -= other;
            this->m_v[1] -= other;
            this->m_v[2] -= other;
        }
        void operator*=(const DTYPE other) {
            this->m_v[0] *= other;
            this->m_v[1] *= other;
            this->m_v[2] *= other;
        }

        DTYPE Dot(const GenericVec3& a) const {
            return (m_v[0] * a[0]) + (m_v[1] * a[1]) + (m_v[2] * a[2]);
        }

        const GenericVec3 Cross(const GenericVec3& b){
            GenericVec3 c;
            c[0] = m_v[1]*b[2] - m_v[2]*b[1];
            c[1] = m_v[2]*b[0] - m_v[0]*b[2];
            c[2] = m_v[0]*b[1] - m_v[1]*b[0];
            return c;
        }

        DTYPE MagSq() const {
            return (m_v[0] * m_v[0]) + (m_v[1] * m_v[1]) + (m_v[2] * m_v[2]);
        }

        DTYPE Mag() const {
            return sqrt(MagSq());
        }
        GenericVec3 operator+(const GenericVec3  &a) const {
            GenericVec3 res(*this);
            res += a;
            return res;
        }
        GenericVec3 operator-(const GenericVec3  &a) const {
            GenericVec3 res(*this);
            res -= a;
            return res;
        }

        GenericVec3 operator+(const DTYPE  val) const {
            GenericVec3 res(*this);
            res.m_v[0] += val; res.m_v[1] += val; res.m_v[2] += val;
            return res;
        }
        GenericVec3 operator-(const DTYPE  val) const {
            GenericVec3 res(*this);
            res.m_v[0] -= val; res.m_v[1] -= val; res.m_v[2] -= val;
            return res;
        }

        GenericVec3 operator*(const DTYPE  val) const {
            GenericVec3 res(this->m_v[0], this->m_v[1], this->m_v[2]);// = (*this);
            res.m_v[0] *= val; res.m_v[1] *= val; res.m_v[2] *= val;
            return res;
        }
        bool operator==(const GenericVec3& other) const {
            return memcmp(this, &other, sizeof(GenericVec3)) == 0;
        }

        GenericVec2<DTYPE> subset(int omit) {
            GenericVec2<DTYPE> v;
            int pos = 0;
            for (int i = 0; i < 3; i++)
                if (i != omit) v[pos++] = m_v[i];
            return v;
        }

        operator GenericVec3<int>() const { return GenericVec3<int>(this->m_v[0], this->m_v[1], this->m_v[2]); }
        operator GenericVec3<long long>() const { return GenericVec3<long long>(this->m_v[0], this->m_v[1], this->m_v[2]); }
        operator GenericVec3<double>() const { return GenericVec3<double>(this->m_v[0], this->m_v[1], this->m_v[2]); }

        // THIS GIVES REMAINDER, NOT MODULO
        GenericVec3 operator%(const GenericVec3& other) const{
            GenericVec3 res;
            res.m_v[0] = this->m_v[0] % other.m_v[0];
            res.m_v[1] = this->m_v[1] % other.m_v[1];
            res.m_v[2] = this->m_v[2] % other.m_v[2];
            return res;
        }

        DTYPE& operator[](int i) { return m_v[i]; }
        const DTYPE& operator[](int i) const { return m_v[i]; }

        void PrintInt() const {
            printf("(%d, %d, %d)\n", m_v[0], m_v[1], m_v[2]);
        }
        void PrintFloat() const {
            printf("(%f, %f, %f)\n", m_v[0], m_v[1], m_v[2]);
        }

        // this is a hack that rounds towards -infinity for "small" negative values
        GenericVec3<INDEX_TYPE> IntFloor() const {
            //(*this).print_vf();
            GenericVec3 dres = *this;
            //dres = dres + 1000.0;
            //dres.print_vf();
            GenericVec3<INDEX_TYPE> ires(floor(dres[0]), floor(dres[1]), floor(dres[2]));
            //ires.print_vi();
            //ires = ires + -1000;
            //ires.print_vi();
            return ires;
        }
        GenericVec3<double> DoubleFloor() const {
            GenericVec3<double> ddres(floor(this->m_v[0]), floor(this->m_v[1]), floor(this->m_v[2]));
            return ddres;
        }
        // linear interpolation betwee 2 vectors: returns the vector a * (1-t) + b * t
        static GenericVec3<double> Lerp(const GenericVec3<double>& a, const GenericVec3<double>& b, const double t)  {
            return (a * (1.0 - t)) + (b * (t));
        }



    };




    typedef GenericVec3<int> Vec3i;
    typedef GenericVec3<INDEX_TYPE> Vec3l;
    typedef GenericVec3<float> Vec3f;
    typedef GenericVec3<double> Vec3d;
    typedef GenericVec3<bool> Vec3b;



}

#endif
