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

#ifndef STATIC_MATRIX_HH
#define STATIC_MATRIX_HH

#include "StaticVector.hh"

template<uint32_t M, uint32_t N, typename T = double>
class StaticMatrix
{
    public:
        typedef StaticMatrix<M, N, T> this_type;
        typedef StaticVector<N, T> row_type;
        typedef StaticVector<M, T> col_type;

        StaticMatrix();
        ~StaticMatrix();

        uint32_t numRows() const;
        uint32_t numCols() const;

        void fill(const T& x);

        T& operator()(uint32_t i, uint32_t j);
        T operator()(uint32_t i, uint32_t j) const;

        this_type operator+(const this_type& m) const;
        this_type& operator+=(const this_type& m);

        this_type operator-(const this_type& m) const;
        this_type& operator-=(const this_type& m);

        template <uint32_t L>
        StaticMatrix<M, L, T> operator*(const StaticMatrix<N, L, T>& m) const;

        this_type operator*(const T& alpha) const;
        this_type& operator*=(const T& alpha);

        const StaticVector<N, T>& getRow(uint32_t i) const;
        StaticVector<N, T>& getRow(uint32_t i);
        StaticVector<M, T> getCol(uint32_t j) const;

    protected:
        StaticVector<N,T> m[M];
};

#include "../tsrc/StaticMatrix_impl.hh"

#endif // STATIC_MATRIX_HH
