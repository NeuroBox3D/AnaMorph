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
#ifndef VEC2_H
#define VEC2_H


#include "StaticVector.hh"

class Vec2
: public StaticVector<2, double>
{
    public:
        typedef StaticVector<2, double> base_type;

        // constructors
        Vec2();
        explicit Vec2(double v);
        Vec2(double x, double y);
        Vec2(const Vec2& v);

        // destructor
        ~Vec2();

        // assignment operator
        Vec2& operator=(const Vec2& v);

        // resizing (only for compatibility)
        void resize(uint32_t size);

        // arithmetic
        Vec2 operator+(const Vec2& v) const;
        Vec2& operator+=(const Vec2& v);
        Vec2 operator-(const Vec2& v) const;
        Vec2& operator-=(const Vec2& v);

        Vec2 operator*(double x) const;
        Vec2& operator*=(double x);
        Vec2 operator/(double x) const;
        Vec2& operator/=(double x);

        // scalar product, cross product
        double operator*(const Vec2& v) const;
        double cross(const Vec2& v) const;

        // relational operators
        bool operator==(const Vec2& v) const;
        bool operator!=(const Vec2& v);

        // tuple-like comparison operators
        bool operator<(const Vec2& v) const;
        bool operator>(const Vec2& v) const;
        bool operator>=(const Vec2& v) const;
        bool operator<=(const Vec2& v) const;

        // norm
        double len2(void) const;
        double len2squared(void) const;

        Vec2 &normalize();

        // output
        void print() const;
        void print_debugl(uint32_t level) const;

    private:
        using StaticVector<2, double>::v;
};


#endif
