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
#include "Vector.hh"

template<typename T>
uint32_t
Vector<T>::n() const
{
    return (this->dims[0]);
}

template<typename T>
Vector<T>::Vector() : Tensor<1, T>()
{
}

template<typename T>
Vector<T>::Vector(
    uint32_t    n,
    T const    &x) : Tensor<1, T>({n})
{
    this->fill(x);
}

template<typename T>
Vector<T>::Vector(Vector<T> const &v) : Tensor<1, T>(v)
{
}

template<typename T>
Vector<T> &
Vector<T>::operator=(Vector<T> const &v)
{
    Tensor<1, T>::operator=(v);
    return (*this);
}

template<typename T>
Vector<T>::~Vector()
{
}

template<typename T>
uint32_t
Vector<T>::size() const
{
    return (this->n());
}

template<typename T>
void
Vector<T>::resize(uint32_t n)
{
    Tensor<1, T>::resize({n});
}

template<typename T>
void
Vector<T>::fill(T const &x)
{
    Tensor<1, T>::fill(x);
}

template<typename T>
void
Vector<T>::assign(
    uint32_t n,
    T const &x)
{
    Tensor<1, T>::resize({n});
    Tensor<1, T>::fill(x);
}

template<typename T>
T &
Vector<T>::operator()(uint32_t i)
{
    return Tensor<1, T>::operator()({i});
}

template<typename T>
T
Vector<T>::operator()(uint32_t i) const
{
    return Tensor<1, T>::operator()({i});
}

template<typename T>
T &
Vector<T>::operator[](uint32_t i)
{
    return Tensor<1, T>::operator()({i});
}

template<typename T>
T
Vector<T>::operator[](uint32_t i) const
{
    return Tensor<1, T>::operator()({i});
}

template<typename T>
Vector<T>
Vector<T>::operator+(Vector<T> const &B) const
{
    return (Vector<T>(*this) += B);
}

template<typename T>
Vector<T> &
Vector<T>::operator+=(Vector<T> const &B)
{
    Tensor<1, T>::operator+=(B);
    return (*this);
}

template<typename T>
Vector<T>
Vector<T>::operator-(Vector<T> const &B) const
{
    return (Vector<T>(*this) -= B);
}

template<typename T>
Vector<T> &
Vector<T>::operator-=(Vector<T> const &B)
{
    Tensor<1, T>::operator-=(B);
    return (*this);
}

template<typename T>
Vector<T>
Vector<T>::operator*(T const &alpha) const
{
    return (Vector<T>(*this) *= alpha);
}

template<typename T>
Vector<T> &
Vector<T>::operator*=(T const &alpha)
{
    Tensor<1, T>::operator*=(alpha);
    return (*this);
}
