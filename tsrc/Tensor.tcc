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

template <uint32_t K, typename T>
void
Tensor<K, T>::checkDims(
    std::string         fn,
    Tensor<K, T> const &x) const
{
    if (this->dims != x.dims) {
        throw(std::string(fn + ": tensor dimension mismatch."));
    }
}


template <uint32_t K, typename T>
Tensor<K, T>::Tensor()
{
    /* set all dimensions to one, i.e. tensor has exactly one element */
    for (uint32_t i = 0; i < K; i++) {
        this->dims[i] = 1;
    }

    /* resize components to 1 with default constructor for type T */
    this->comps.resize(1);
}


template <uint32_t K, typename T>
Tensor<K, T>::Tensor(std::array<uint32_t, K> const  &dims)
{

    /* check if all dims are > 0, calculate number of elements, which is the product of all dims (assuming all are >= 1)
     * and the offsets for element position calulation:
     * 
     * if (n_0, .., n_{K-1}) are the dimensions, the offset o_i is defined as the product
     *
     *          o_i = (n_{0} * .. * n_{i-1})
     *
     * and o_0 = 1.
     *
     * if indices (i_0, .., i_{K-1}) are given, the position index of the corresponding tensor element in the
     * one-dimensional component std::vector this->comps is calculated as
     *
     * i_0*o_0 + .. i_{K-1} * o_{K-1}
     * */

    uint32_t nelements = 1;
    for (uint32_t i = 0; i < K; i++) {
        if (dims[i] == 0) {
            throw("Tensor::Tensor(dims, comps): zero-dimension discovered. all dimensions must be >= 1.");
        }

        nelements *= dims[i];

        if (i == 0) {
            this->offsets[i] = 1;
        }
        /* i >= 1 */
        else {
            this->offsets[i] = this->offsets[i-1] * dims[i-1];
        }
    }

    /* assign dimension array, resize (one-dimensional) components std::vector to the correct size using default
     * constructor for type T. */
    this->dims  = dims;
    this->comps.resize(nelements);
}


template <uint32_t K, typename T>
Tensor<K, T>::Tensor(Tensor<K, T> const &x)
{
    this->comps     = x.comps;
    this->dims      = x.dims;
    this->offsets   = x.offsets;
}

template <uint32_t K, typename T>
Tensor<K, T> &
Tensor<K, T>::operator=(Tensor<K, T> const &x)
{
    this->comps     = x.comps;
    this->dims      = x.dims;
    this->offsets   = x.offsets;

    return (*this);
}

template <uint32_t K, typename T>
Tensor<K, T>::~Tensor()
{
}

/* coefficient access: compute index in one-dimensional component vector from given index array */
template <uint32_t K, typename T>
T const &
Tensor<K, T>::operator()(std::array<uint32_t, K> const &indices) const
{
    uint32_t pos = 0;
    for (uint32_t i = 0; i < K; i++) {
        if (indices[i] <= this->dims[i]) {
            pos += indices[i] * this->offsets[i];
        }
        else {
            throw("Tensor::operator(): index out of range in at least one dimension.");
        }
    }

    return (this->comps[pos]);
}

template <uint32_t K, typename T>
T &
Tensor<K, T>::operator()(std::array<uint32_t, K> const &indices)
{
    uint32_t pos = 0;
    for (uint32_t i = 0; i < K; i++) {
        if (indices[i] <= this->dims[i]) {
            pos += indices[i] * this->offsets[i];
        }
        else {
            throw("Tensor::operator(): index out of range in at least one dimension.");
        }
    }

    return (this->comps[pos]);
}

template <uint32_t K, typename T>
void
Tensor<K, T>::resize(std::array<uint32_t, K> const &dims)
{
    (*this) = Tensor<K, T>(dims);
}

template <uint32_t K, typename T>
void
Tensor<K, T>::fill(T const &x)
{
    std::fill(this->comps.begin(), this->comps.end(), x);
}

template <uint32_t K, typename T>
Tensor<K, T>
Tensor<K, T>::operator+(Tensor<K, T> const &x) const
{
    return (Tensor<K, T>(*this) += x);
}

template <uint32_t K, typename T>
Tensor<K, T> &
Tensor<K, T>::operator+=(Tensor<K, T> const &x)
{
    this->checkDims("Tensor::operator+=()", x);
    for (uint32_t i = 0; i < this->comps.size(); i++) {
        this->comps[i] += x.comps[i];
    }

    return (*this);
}

template <uint32_t K, typename T>
Tensor<K, T>
Tensor<K, T>::operator-(Tensor<K, T> const &x) const
{
    return (Tensor<K, T>(*this) -= x);
}

template <uint32_t K, typename T>
Tensor<K, T> &
Tensor<K, T>::operator-=(Tensor<K, T> const &x)
{
    this->checkDims("Tensor::operator-=()", x);

    for (uint32_t i = 0; i < this->comps.size(); i++) {
        this->comps[i] -= x.comps[i];
    }

    return (*this);
}

template <uint32_t K, typename T>
Tensor<K, T>
Tensor<K, T>::operator*(T const &x) const
{
    return (Tensor<K, T>(*this) *= x);
}

template <uint32_t K, typename T>
Tensor<K, T> &
Tensor<K, T>::operator*=(T const &x)
{

    for (uint32_t i = 0; i < this->comps.size(); i++) {
        this->comps[i] *= x;
    }

    return (*this);
}

template <uint32_t K, typename T>
Tensor<K, T>
Tensor<K, T>::operator/(T const &x) const
{
    return (Tensor<K, T>(*this) /= x);
}

template <uint32_t K, typename T>
Tensor<K, T> &
Tensor<K, T>::operator/=(T const &x)
{

    for (uint32_t i = 0; i < this->comps.size(); i++) {
        this->comps[i] /= x;
    }

    return (*this);
}

template <uint32_t K, typename T>
bool
Tensor<K, T>::operator==(Tensor<K, T> const &x) const
{
    this->checkDims("Tensor::operator==()", x);
    return (this->comps == x.comps);
}

template <uint32_t K, typename T>
bool
Tensor<K, T>::operator!=(Tensor<K, T> const &x) const
{
    this->checkDims("Tensor::operator!=()", x);
    return (this->comps != x.comps);
}

template <uint32_t K, typename T>
bool
Tensor<K, T>::operator<(Tensor<K, T> const &x) const
{
    this->checkDims("Tensor::operator<()", x);
    return (this->comps < x.comps);
}

template <uint32_t K, typename T>
bool
Tensor<K, T>::operator<=(Tensor<K, T> const &x) const
{
    this->checkDims("Tensor::operator<=()", x);
    return (this->comps <= x.comps);
}

template <uint32_t K, typename T>
bool
Tensor<K, T>::operator>(Tensor<K, T> const &x) const
{
    this->checkDims("Tensor::operator>()", x);
    return (this->comps > x.comps);
}

template <uint32_t K, typename T>
bool
Tensor<K, T>::operator>=(Tensor<K, T> const &x) const
{
    this->checkDims("Tensor::operator>=()", x);
    return (this->comps >= x.comps);
}
