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
#include "StaticMatrix.hh"


template<uint32_t M, uint32_t N, typename T>
StaticMatrix<M, N, T>::StaticMatrix()
{}

template<uint32_t M, uint32_t N, typename T>
StaticMatrix<M, N, T>::~StaticMatrix()
{}

template<uint32_t M, uint32_t N, typename T>
uint32_t
StaticMatrix<M, N, T>::numRows() const
{
    return M;
}

template<uint32_t M, uint32_t N, typename T>
uint32_t
StaticMatrix<M, N, T>::numCols() const
{
    return N;
}

template<uint32_t M, uint32_t N, typename T>
void StaticMatrix<M, N, T>::fill(const T& x)
{
    for (size_t i = 0; i < M; ++i)
        m[i].assign(x);
}

template<uint32_t M, uint32_t N, typename T>
T&
StaticMatrix<M, N, T>::operator()(uint32_t i, uint32_t j)
{
    assert(i < M && j < N && "Invalid component requested.");
    return m[i][j];
}

template<uint32_t M, uint32_t N, typename T>
T
StaticMatrix<M, N, T>::operator()(uint32_t i, uint32_t j) const
{
    assert(i < M && j < N && "Invalid component requested.");
    return m[i][j];
}

template<uint32_t M, uint32_t N, typename T>
StaticMatrix<M, N, T>
StaticMatrix<M, N, T>::operator+(const this_type& _m) const
{
    this_type r;
    for (size_t i = 0; i < M; ++i)
        r.m[i] = m[i] + _m.m[i];
    return r;
}

template<uint32_t M, uint32_t N, typename T>
StaticMatrix<M, N, T> &
StaticMatrix<M, N, T>::operator+=(const this_type& _m)
{
    for (size_t i = 0; i < M; ++i)
        m[i] += _m.m[i];
    return *this;
}

template<uint32_t M, uint32_t N, typename T>
StaticMatrix<M, N, T>
StaticMatrix<M, N, T>::operator-(const this_type& _m) const
{
    this_type r;
    for (size_t i = 0; i < M; ++i)
        r.m[i] = m[i] - _m.m[i];
    return r;
}

template<uint32_t M, uint32_t N, typename T>
StaticMatrix<M, N, T> &
StaticMatrix<M, N, T>::operator-=(const this_type& _m)
{
    for (size_t i = 0; i < M; ++i)
        m[i] -= _m.m[i];
    return *this;
}

template<uint32_t M, uint32_t N, typename T>
template <uint32_t L>
StaticMatrix<M, L, T>
StaticMatrix<M, N, T>::operator*(const StaticMatrix<N, L, T>& _m) const
{
    StaticMatrix<M, L, T> r;
    for (size_t i = 0; i < M; ++i)
    {
        for (size_t j = 0; j < L; ++j)
        {
            T& t = r(i,j);
            for (size_t k = 0; k < N; ++j)
                t += m[i][k] * _m(k,j);
        }
    }
    return r;
}

template<uint32_t M, uint32_t N, typename T>
StaticMatrix<M, N, T>
StaticMatrix<M, N, T>::operator*(const T& alpha) const
{
    this_type r;
    for (size_t i = 0; i < M; ++i)
        r[i] = alpha * m[i];
    return r;
}

template<uint32_t M, uint32_t N, typename T>
StaticMatrix<M, N, T> &
StaticMatrix<M, N, T>::operator*=(const T& alpha)
{
    for (size_t i = 0; i < M; ++i)
        m[i] *= alpha;
    return *this;
}

template<uint32_t M, uint32_t N, typename T>
const StaticVector<N, T>&
StaticMatrix<M, N, T>::getRow(uint32_t i) const
{
    return m[i];
}

template<uint32_t M, uint32_t N, typename T>
StaticVector<N, T>&
StaticMatrix<M, N, T>::getRow(uint32_t i)
{
    return m[i];
}

template<uint32_t M, uint32_t N, typename T>
StaticVector<M, T>
StaticMatrix<M, N, T>::getCol(uint32_t j) const
{
    StaticVector<M, T> r;
    for (size_t i = 0; i < M; ++i)
        r[i] = m[i][j];
    return r;
}

