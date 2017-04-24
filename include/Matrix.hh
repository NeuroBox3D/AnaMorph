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

#ifndef MATRIX_H
#define MATRIX_H

#include "Tensor.hh"

template<typename T = double>
class Matrix : public Tensor<2, T> {
    private:
        inline uint32_t                     m() const;
        inline uint32_t                     n() const;

    public:
                                            Matrix();
                                            Matrix(uint32_t m, uint32_t n);

        uint32_t                            numRows() const;
        uint32_t                            numCols() const;

        void                                resize(uint32_t m, uint32_t n);
        void                                fill(T const &x);                            

        T                                  &operator()(uint32_t i, uint32_t j);
        T                                   operator()(uint32_t i, uint32_t j) const;

        Matrix<T>                           operator+(Matrix<T> const &B) const;
        Matrix<T>                          &operator+=(Matrix<T> const &B);

        Matrix<T>                           operator-(Matrix<T> const &B) const;
        Matrix<T>                          &operator-=(Matrix<T> const &B);

        Matrix<T>                           operator*(Matrix<T> const &B) const;

        Matrix<T>                           operator*(T const &alpha) const;
        Matrix<T>                          &operator*=(T const &alpha);

        Vector<T>                           getRow(uint32_t i) const;
        Vector<T>                           getCol(uint32_t j) const;
};

#include "../tsrc/Matrix_impl.hh"

#endif
