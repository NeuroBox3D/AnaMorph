/* --------------------------------------------------------------------------------
 * 
 *                              THIS FILE IS PART OF                               
 * 
 * AnaMorph: A Framework for Geometric Modelling, Consistency Analysis and Surface
 * Mesh Generation of Anatomically Reconstructed Neuron Morphologies
 * 
 * Web site: http://www.anamorph.net
 * 
 * Copyright (c) 2013-2014, Konstantin Mörschel.
 * All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 * 3. All advertising materials mentioning features or use of this software
 *    must display the following acknowledgement:
 * 
 *    This product includes software developed by Konstantin Mörschel for
 *    AnaMorph (http://www.anamorph.net).
 * 
 * 4. Neither the name "AnaMorph" nor the names of its contributors may be
 *    used to endorse or promote products derived from this software without
 *    specific prior written permission.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE CONTRIBUTORS OF ANAMORPH ''AS IS'' AND ANY
 * EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE CONTRIBUTORS OF ANAMORPH BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 * 
 * -------------------------------------------------------------------------------- */

#ifndef VEC3_H
#define VEC3_H

#include "Vector.hh"
#include "debug.hh"

template <typename R>
class Vec3 : public Vector<R> {
    public:
                                            Vec3();
                                            Vec3(
                                                R const &x,
                                                R const &y,
                                                R const &z);
        /* compatibility size-based constructor and resize() that accepts size 3 only. useful if Vec3 is used as a
         * template argument that requires this type of constructor */
                                            Vec3(uint32_t n);

                                            Vec3(Vec3 const &v);
        Vec3                               &operator=(Vec3 const &v);

                                           ~Vec3();

        void                                resize(uint32_t size);

        /* component access */
        R                                  &operator()(uint32_t i);
        R const                            &operator()(uint32_t i) const;
        R                                  &operator[](uint32_t i);
        R const                            &operator[](uint32_t i) const;

        /* arithmetic */
        Vec3                                operator+(Vec3 const &v) const;
        Vec3                               &operator+=(Vec3 const &v);

        Vec3                                operator-(Vec3 const &v) const;
        Vec3                               &operator-=(Vec3 const &v);

        Vec3                                operator*(R const &x) const;
        Vec3                               &operator*=(R const &x);

        Vec3                                operator/(R const &x) const;
        Vec3                               &operator/=(R const &x);

        /* scalar product, inner product */
        R                                   operator*(Vec3 const &v) const;

        /* cross product */
        Vec3<R>                             cross(Vec3 const &v) const;

        /* relational operators */
        bool                                operator==(Vec3 const &v) const;
        bool                                operator!=(Vec3 const &v);

        /* tuple-like comparison operators: a < b iff a[0] < b[0] && .. && a[2] < b[2].
         *
         * a > b iff !(a == b) && !(a < b) */
        bool                                operator<(Vec3 const &v) const;
        bool                                operator>(Vec3 const &v) const;
        bool                                operator>=(Vec3 const &v) const;
        bool                                operator<=(Vec3 const &v) const;

        /* length for p2-norm (eucl) */
        R                                   len2(void) const;
        R                                   len2squared(void) const;

        Vec3<R>                            &normalize();

        void                                print() const;
        void                                print_debugl(uint32_t level) const;
};

#include "Vec3.tcc"

#endif
