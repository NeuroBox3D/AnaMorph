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

#include "../include/Vec3.hh"
#include "debug.hh"

template <typename R>
Vec3<R>::Vec3()
{}

template <typename R>
Vec3<R>::Vec3(R _v)
{
    v[0] = _v;
    v[1] = _v;
    v[2] = _v;
}

template <typename R>
Vec3<R>::Vec3(R x, R y, R z)
{
    v[0] = x;
    v[1] = y;
    v[2] = z;
}


template <typename R>
Vec3<R>::Vec3(const Vec3<R>& _v)
{
    v[0] = _v[0];
    v[1] = _v[1];
    v[2] = _v[2];
}


template <typename R>
Vec3<R>::~Vec3()
{}


template <typename R>
Vec3<R>& Vec3<R>::operator=(const Vec3<R>& _v)
{
    v[0] = _v[0];
    v[1] = _v[1];
    v[2] = _v[2];

    return *this;
}


template <typename R>
void Vec3<R>::resize(uint32_t size)
{
    if (size != 3)
        throw("Vec3<R>::resize(size_t): given size != 3.");
}


template <typename R>
Vec3<R> Vec3<R>::operator+(const Vec3<R>& _v) const
{
    Vec3<R> r;
    r[0] = v[0] + _v[0];
    r[1] = v[1] + _v[1];
    r[2] = v[2] + _v[2];
    return r;
}

template <typename R>
Vec3<R>& Vec3<R>::operator+=(const Vec3<R>& _v)
{
    v[0] += _v[0];
    v[1] += _v[1];
    v[2] += _v[2];
    return *this;
}

template <typename R>
Vec3<R> Vec3<R>::operator-(const Vec3<R>& _v) const
{
    Vec3<R> r;
    r[0] = v[0] - _v[0];
    r[1] = v[1] - _v[1];
    r[2] = v[2] - _v[2];
    return r;
}

template <typename R>
Vec3<R>& Vec3<R>::operator-=(const Vec3<R>& _v)
{
    v[0] -= _v[0];
    v[1] -= _v[1];
    v[2] -= _v[2];
    return *this;
}


template <typename R>
Vec3<R> Vec3<R>::operator*(R x) const
{
    Vec3<R> r;
    r[0] = v[0] * x;
    r[1] = v[1] * x;
    r[2] = v[2] * x;
    return r;
}

template <typename R>
Vec3<R>& Vec3<R>::operator*=(R x)
{
    v[0] *= x;
    v[1] *= x;
    v[2] *= x;
    return *this;
}

template <typename R>
Vec3<R> Vec3<R>::operator/(R x) const
{
    Vec3<R> r;
    r[0] = v[0] / x;
    r[1] = v[1] / x;
    r[2] = v[2] / x;
    return r;
}

template <typename R>
Vec3<R>& Vec3<R>::operator/=(R x)
{
    v[0] /= x;
    v[1] /= x;
    v[2] /= x;
    return *this;
}


template <typename R>
R Vec3<R>::operator*(const Vec3<R>& _v) const
{
    return v[0]*_v[0] + v[1]*_v[1] + v[2]*_v[2];
}

template <typename R>
Vec3<R> Vec3<R>::cross(const Vec3<R>& _v) const
{
    Vec3<R> r;
    r[0] = v[1]*_v[2] - v[2]*_v[1];
    r[1] = v[2]*_v[0] - v[0]*_v[2];
    r[2] = v[0]*_v[1] - v[1]*_v[0];
    return r;
}


template <typename R>
bool Vec3<R>::operator==(const Vec3<R>& _v) const
{
    return v[0] == _v[0] && v[1] == _v[1] && v[2] == _v[2];
}

template <typename R>
bool Vec3<R>::operator!=(const Vec3<R>& _v)
{
    return v[0] != _v[0] || v[1] != _v[1] || v[2] != _v[2];
}


template <typename R>
bool Vec3<R>::operator<(const Vec3<R>& _v) const
{
    return v[0] < _v[0] && v[1] < _v[1] && v[2] < _v[2];
}

template <typename R>
bool Vec3<R>::operator>(const Vec3<R>& _v) const
{
    return v[0] > _v[0] && v[1] > _v[1] && v[2] > _v[2];
}

template <typename R>
bool Vec3<R>::operator>=(const Vec3<R>& _v) const
{
    return v[0] >= _v[0] && v[1] >= _v[1] && v[2] >= _v[2];
}

template <typename R>
bool Vec3<R>::operator<=(const Vec3<R>& _v) const
{
    return v[0] <= _v[0] && v[1] <= _v[1] && v[2] <= _v[2];
}


template <typename R>
R Vec3<R>::len2(void) const
{
    return sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
}

template <typename R>
R Vec3<R>::len2squared(void) const
{
    return v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
}


template <typename R>
Vec3<R>& Vec3<R>::normalize()
{
    R n = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
    v[0] /= n;
    v[1] /= n;
    v[2] /= n;
    return *this;
}


template <typename R>
void Vec3<R>::print() const
{
    printf("(%f, %f, %f)\n", v[0], v[1], v[2]);
}

template <typename R>
void Vec3<R>::print_debugl(uint32_t level) const
{
    debugTabInc();
    debugl(level, "(%f, %f, %f)\n", v[0], v[1], v[2]);
    debugTabDec();
}



// make Vec3<double> usable
template class Vec3<double>;
