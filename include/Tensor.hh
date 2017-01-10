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

#ifndef TENSOR_HH
#define TENSOR_HH

#include "common.hh"

/*! @brief template class modelling tensors of arbitrary order.
 *
 * @tparam K non-negative order of tensor. must be known at compile-time.
 * @tparam T type of objects contained in tensor. */
template <
    uint32_t    K,
    typename    T
>
class Tensor {
    private:
        void                            checkDims(
                                            const std::string&         fn,
                                            Tensor<K, T> const &x) const; 
    protected:
        std::vector<T>                  comps;
        std::array<uint32_t, K>         dims;
        std::array<uint32_t, K>         offsets;

    public:
                                        Tensor();
        explicit                        Tensor(std::array<uint32_t, K> const  &dims);

        explicit                        Tensor(Tensor<K, T> const &x);
        Tensor<K, T>                   &operator=(Tensor<K, T> const &x);

        /* virtual dtor for polymorphism */
        virtual                        ~Tensor();


        /* component access, both as reference (non-const access) and value (const access).
           since operator[] cannot be overloaded to take more than one argument, operator() is used
           for the general case. */
        T const                        &operator()(std::array<uint32_t, K> const &indices) const;
        T                              &operator()(std::array<uint32_t, K> const &indices);

        void                            resize(std::array<uint32_t, K> const &dims);
        void                            fill(T const &x);

        /* arithmetic: derived classes can reuse this code, however the versions returning a new
         * object (operator{+, -, *, /} must be used with care, since they create objects of type
         * Tensor<K, T>, i.e. of base class type from the point of view of a derived class, say
         * Matrix<T>.  arithmetic operators usually need to be wrapped in derived classes to obtain
         * the desired semantics.  e.g.: Matrix<T> arithmetic should return Matrix<T> objects, which
         * might contain more specific information than Tensor<2, T>. however, an assignment
         * Matrix<T>::operator=(Tensor<2, T> const &x) will fix the issue as well. */
        Tensor<K, T>                    operator+(Tensor<K, T> const &x) const;
        Tensor<K, T>                   &operator+=(Tensor<K, T> const &x);

        Tensor<K, T>                    operator-(Tensor<K, T> const &x) const;
        Tensor<K, T>                   &operator-=(Tensor<K, T> const &x);

        Tensor<K, T>                    operator*(T const &x) const;
        Tensor<K, T>                   &operator*=(T const &x);

        Tensor<K, T>                    operator/(T const &x) const;
        Tensor<K, T>                   &operator/=(T const &x);

        /* relational operators */
        bool                            operator==(Tensor<K, T> const &x) const;
        bool                            operator!=(Tensor<K, T> const &x) const;
        bool                            operator<(Tensor<K, T> const &x) const;
        bool                            operator<=(Tensor<K, T> const &x) const;
        bool                            operator>(Tensor<K, T> const &x) const;
        bool                            operator>=(Tensor<K, T> const &x) const;
};

#include "../tsrc/Tensor.impl.hh"

#endif
