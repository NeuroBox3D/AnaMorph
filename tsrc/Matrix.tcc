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

#include "common.hh"
#include "Matrix.hh"

template<typename T>
uint32_t
Matrix<T>::m() const
{
    return (this->dims[0]);
}

template<typename T>
uint32_t
Matrix<T>::n() const
{
    return (this->dims[1]);
}

template<typename T>
Matrix<T>::Matrix() : Tensor<2, T>()
{
}

template<typename T>
Matrix<T>::Matrix(uint32_t m, uint32_t n) : Tensor<2, T>({m, n})
{
}

template<typename T>
uint32_t
Matrix<T>::numRows() const
{
    return (this->m());
}

template<typename T>
uint32_t
Matrix<T>::numCols() const
{
    return (this->n());
}

template<typename T>
void
Matrix<T>::resize(uint32_t m, uint32_t n)
{
    Tensor<2, T>::resize({m, n});
}

template<typename T>
void
Matrix<T>::fill(T const &x)
{
    Tensor<2, T>::fill(x);
}

template<typename T>
T &
Matrix<T>::operator()(uint32_t i, uint32_t j)
{
    return Tensor<2, T>::operator()({i, j});
}

template<typename T>
T
Matrix<T>::operator()(uint32_t i, uint32_t j) const
{
    return Tensor<2, T>::operator()({i, j});
}

template<typename T>
Matrix<T>
Matrix<T>::operator+(Matrix<T> const &B) const
{
    return (Matrix<T>(*this) += B);
}

template<typename T>
Matrix<T> &
Matrix<T>::operator+=(Matrix<T> const &B)
{
    Tensor<2, T>::operator+=(B);
    return (*this);
}

template<typename T>
Matrix<T>
Matrix<T>::operator-(Matrix<T> const &B) const
{
    return (Matrix<T>(*this) -= B);
}

template<typename T>
Matrix<T> &
Matrix<T>::operator-=(Matrix<T> const &B)
{
    Tensor<2, T>::operator-=(B);
    return (*this);
}

template<typename T>
Matrix<T>
Matrix<T>::operator*(T const &alpha) const
{
    return (Matrix<T>(*this) *= alpha);
}

template<typename T>
Matrix<T> &
Matrix<T>::operator*=(T const &alpha)
{
    Tensor<2, T>::operator*=(alpha);
    return (*this);
}

template<typename T>
Matrix<T>
Matrix<T>::operator*(Matrix<T> const &B) const
{
    if (n() != B.m()) {
        throw("Matrix()::operator*: dimension mismatch.");
    }

    uint32_t i, j, k;
    Matrix<T> C(this->m(), B.n());

    for (i = 0; i < this->m(); i++) {
        for (j = 0; j < B.n(); j++) {
            C(i, j) = 0.0;
            for (k = 0; k < this->n(); k++) {
                C(i, j) += this->operator()(i, k) * B(k, j);
            }
        }
    }
    return C;
}

template<typename T>
Vector<T>
Matrix<T>::getRow(uint32_t i) const
{
    Vector<T> row(this->n());
    for (uint32_t j = 0; j < this->n(); j++) {
        row[j] = this->operator()(i, j);
    }
    return row;
}

template<typename T>
Vector<T>
Matrix<T>::getCol(uint32_t j) const
{
    Vector<T> col(this->m());
    for (uint32_t i = 0; i < this->m(); i++) {
        col[i] = this->operator()(i, j);
    }
    return col;
}
