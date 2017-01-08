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
 *
 * author of this file: mbreit
 * creation: 04.01.2017
 * -------------------------------------------------------------------------------- */

#include <assert.h>
#include "StaticVector.hh"


template<uint32_t N, typename T>
StaticVector<N, T>::StaticVector()
{}

template<uint32_t N, typename T>
StaticVector<N, T>::StaticVector(const T& x)
{
    assign(x);
}

template<uint32_t N, typename T>
StaticVector<N, T>::StaticVector(const this_type& _v)
{
    for (size_t i = 0; i < N; ++i)
        v[i] = _v.v[i];
}

template<uint32_t N, typename T>
StaticVector<N, T>::~StaticVector()
{}

template<uint32_t N, typename T>
StaticVector<N, T>&
StaticVector<N, T>::operator=(const this_type& _v)
{
    for (size_t i = 0; i < N; ++i)
        v[i] = _v.v[i];
    return *this;
}

template<uint32_t N, typename T>
uint32_t
StaticVector<N, T>::size() const
{
    return N;
}

template<uint32_t N, typename T>
void
StaticVector<N, T>::assign(const T& x)
{
    for (size_t i = 0; i < N; ++i)
        v[i] = x;
}

template<uint32_t N, typename T>
T&
StaticVector<N, T>::operator()(uint32_t i)
{
    assert(i < N && "Invalid component requested.");
    return v[i];
}

template<uint32_t N, typename T>
T
StaticVector<N, T>::operator()(uint32_t i) const
{
    assert(i < N && "Invalid component requested.");
    return v[i];
}

template<uint32_t N, typename T>
T &
StaticVector<N, T>::operator[](uint32_t i)
{
    assert(i < N && "Invalid component requested.");
    return v[i];
}

template<uint32_t N, typename T>
T
StaticVector<N, T>::operator[](uint32_t i) const
{
    assert(i < N && "Invalid component requested.");
    return v[i];
}

template<uint32_t N, typename T>
StaticVector<N, T>
StaticVector<N, T>::operator+(const this_type& _v) const
{
    StaticVector r;
    for (size_t i = 0; i < N; ++i)
        r[i] = v[i] + _v[i];
    return r;
}

template<uint32_t N, typename T>
StaticVector<N, T> &
StaticVector<N, T>::operator+=(const this_type& _v)
{
    for (size_t i = 0; i < N; ++i)
        v[i] += _v[i];
    return *this;
}

template<uint32_t N, typename T>
StaticVector<N, T>
StaticVector<N, T>::operator-(const this_type& _v) const
{
    StaticVector r;
    for (size_t i = 0; i < N; ++i)
        r[i] = v[i] - _v[i];
    return r;
}

template<uint32_t N, typename T>
StaticVector<N, T> &
StaticVector<N, T>::operator-=(const this_type& _v)
{
    for (size_t i = 0; i < N; ++i)
        v[i] -= _v[i];
    return *this;
}

template<uint32_t N, typename T>
StaticVector<N, T>
StaticVector<N, T>::operator*(const T& alpha) const
{
    StaticVector r;
    for (size_t i = 0; i < N; ++i)
        r[i] = alpha * v[i];
    return r;
}

template<uint32_t N, typename T>
StaticVector<N, T> &
StaticVector<N, T>::operator*=(const T& alpha)
{
    for (size_t i = 0; i < N; ++i)
        v[i] *= alpha;
    return *this;
}
