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

#include "../include/Vec2.hh"
#include "debug.hh"


Vec2::Vec2()
{}

Vec2::Vec2(double _v)
{
    v[0] = _v;
    v[1] = _v;
}

Vec2::Vec2(double x, double y)
{
    v[0] = x;
    v[1] = y;
}


Vec2::Vec2(const Vec2& _v)
{
    v[0] = _v[0];
    v[1] = _v[1];
}


Vec2::~Vec2()
{}


Vec2& Vec2::operator=(const Vec2& _v)
{
    v[0] = _v[0];
    v[1] = _v[1];

    return *this;
}


void Vec2::resize(uint32_t size)
{
    if (size != 2)
        throw("Vec2::resize(size_t): given size != 2.");
}


Vec2 Vec2::operator+(const Vec2& _v) const
{
    Vec2 r;
    r[0] = v[0] + _v[0];
    r[1] = v[1] + _v[1];
    return r;
}

Vec2& Vec2::operator+=(const Vec2& _v)
{
    v[0] += _v[0];
    v[1] += _v[1];
    return *this;
}

Vec2 Vec2::operator-(const Vec2& _v) const
{
    Vec2 r;
    r[0] = v[0] - _v[0];
    r[1] = v[1] - _v[1];
    return r;
}

Vec2& Vec2::operator-=(const Vec2& _v)
{
    v[0] -= _v[0];
    v[1] -= _v[1];
    return *this;
}


Vec2 Vec2::operator*(double x) const
{
    Vec2 r;
    r[0] = v[0] * x;
    r[1] = v[1] * x;
    return r;
}

Vec2& Vec2::operator*=(double x)
{
    v[0] *= x;
    v[1] *= x;
    return *this;
}

Vec2 Vec2::operator/(double x) const
{
    Vec2 r;
    r[0] = v[0] / x;
    r[1] = v[1] / x;
    return r;
}

Vec2& Vec2::operator/=(double x)
{
    v[0] /= x;
    v[1] /= x;
    return *this;
}


double Vec2::operator*(const Vec2& _v) const
{
    return v[0]*_v[0] + v[1]*_v[1];
}

double Vec2::cross(const Vec2& _v) const
{
    return v[0]*_v[1] - v[1]*_v[0];
}


bool Vec2::operator==(const Vec2& _v) const
{
    return v[0] == _v[0] && v[1] == _v[1];
}

bool Vec2::operator!=(const Vec2& _v)
{
    return v[0] != _v[0] || v[1] != _v[1];
}


bool Vec2::operator<(const Vec2& _v) const
{
    return v[0] < _v[0] && v[1] < _v[1];
}

bool Vec2::operator>(const Vec2& _v) const
{
    return v[0] > _v[0] && v[1] > _v[1];
}

bool Vec2::operator>=(const Vec2& _v) const
{
    return v[0] >= _v[0] && v[1] >= _v[1];
}

bool Vec2::operator<=(const Vec2& _v) const
{
    return v[0] <= _v[0] && v[1] <= _v[1];
}


double Vec2::len2(void) const
{
    return sqrt(v[0]*v[0] + v[1]*v[1]);
}

double Vec2::len2squared(void) const
{
    return v[0]*v[0] + v[1]*v[1];
}


Vec2& Vec2::normalize()
{
    double n = sqrt(v[0]*v[0] + v[1]*v[1]);
    v[0] /= n;
    v[1] /= n;
    return *this;
}


void Vec2::print() const
{
    printf("(%f, %f)\n", v[0], v[1]);
}

void Vec2::print_debugl(uint32_t level) const
{
    debugTabInc();
    debugl(level, "(%f, %f)\n", v[0], v[1]);
    debugTabDec();
}


