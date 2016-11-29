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

#ifndef BIVARIATE_POLYNOMIAL_H
#define BIVARIATE_POLYNOMIAL_H

/*! @brief abstract bivariate polynomial base class providing minimal interface */
template<
    typename F,
    typename R
>
class BivariatePolynomial {
    protected:
        uint32_t                            m, n;
        Matrix<F>                           coeff;

    public:
                                            BivariatePolynomial();
                                            BivariatePolynomial(uint32_t m, uint32_t n, F const &x);
                                            BivariatePolynomial(BivariatePolynomial<F, R> const &q);
        BivariatePolynomial<F, R>          &operator=(BivariatePolynomial const &q);

        virtual                            ~BivariatePolynomial();

        void                                getDegree(uint32_t &m, uint32_t &n) const;

        /* substitute for operator[], that unfortunately cannot be overloaded to have two arguments
         * */
        virtual F                          &operator()(uint32_t i, uint32_t j);
        virtual F                           operator()(uint32_t i, uint32_t j) const;

        Matrix<F>                           getCoeffs() const;
        void                                setCoeffs(Matrix<F> const &coeff_arg);
        F                                   getMaxAbsCoeff() const;

        void                                initConstant(F const &x);
        void                                zero();

        virtual F                           eval(R const &x, R const &y) const  = 0;

        /* arithmetic: make this protected just as for Polynomial to prevent misuse with two upcast pointers /
         * references from different derived classes */
    protected:
        /* addition in vector space \Pi^{m, n} */
        BivariatePolynomial<F, R>          &operator+=(BivariatePolynomial<F, R> const &q);
        BivariatePolynomial<F, R>          &operator-=(BivariatePolynomial<F, R> const &q);

        /* scalar multiplication in vector space \Pi^{m, n} */
        BivariatePolynomial<F, R>          &operator*=(F const &x);
        BivariatePolynomial<F, R>          &operator/=(F const &x);

        /* inner product with respect to L2 norm \int_0^1\int_0^1{p(x,y)q(x,y) dxdy} */
        //F                                   operator*(BivariatePolynomial<F, R> const &q) const;

    public:
        /* print coefficient matrix */
        void                                printCoeff() const;

        /* write plotfile suitable for gnuplot */
        void                                writePlotFile(uint32_t ticks, std::string filename) const;
};

/*! @brief class for bivariate polynomials represented in Bernstein-Bezier basis */
template<
    typename F,
    typename R
>
class BiBernsteinPolynomial : public BivariatePolynomial<F, R> {
    public:
                                            BiBernsteinPolynomial(); 
                                            BiBernsteinPolynomial(uint32_t m, uint32_t n, F const &x);
                                            BiBernsteinPolynomial(Matrix<F> const &coeff_arg);

                                            BiBernsteinPolynomial(BiBernsteinPolynomial<F, R> const &q);
        BiBernsteinPolynomial<F, R>        &operator=(BiBernsteinPolynomial<F, R> const &q);

        virtual                            ~BiBernsteinPolynomial();

        /* implementation of BivariatePolynomial interface */
        F                                   eval(
                                                R const &x,
                                                R const &y) const;

        /* addition in vector space \Pi^{m, n} */
        BiBernsteinPolynomial<F, R>         operator+(BiBernsteinPolynomial<F, R> const &q) const;
        BiBernsteinPolynomial<F, R>        &operator+=(BiBernsteinPolynomial<F, R> const &q);

        BiBernsteinPolynomial<F, R>         operator-(BiBernsteinPolynomial<F, R> const &q) const;
        BiBernsteinPolynomial<F, R>        &operator-=(BiBernsteinPolynomial<F, R> const &q);

        /* scalar multiplication in vector space \Pi^{m, n} */
        BiBernsteinPolynomial<F, R>         operator*(F const &x) const;
        BiBernsteinPolynomial<F, R>        &operator*=(F const &x);

        BiBernsteinPolynomial<F, R>         operator/(F const &x) const;
        BiBernsteinPolynomial<F, R>        &operator/=(F const &x);

        /* inner product with respect to L2 norm \int_0^1\int_0^1{p(x,y)q(x,y) dxdy} */
        F                                   operator*(BiBernsteinPolynomial<F, R> const &q) const;

        /* multiplication as polynomials over |R^2 or a suitably defined ring */
        BiBernsteinPolynomial                multiply(BiBernsteinPolynomial const &q) const;
        BiBernsteinPolynomial                square() const;

        /* other methods */
        void                                split_x(
                                                R const                        &x,
                                                BiBernsteinPolynomial<F, R>    *pleft,
                                                BiBernsteinPolynomial          *pright) const;

        void                                split_y(
                                                R const                        &y,
                                                BiBernsteinPolynomial<F, R>    *pdown,
                                                BiBernsteinPolynomial          *pup) const;

        void                                split_xy(
                                                R const                        &x,
                                                R const                        &y,
                                                BiBernsteinPolynomial<F, R>    *pleft_down,
                                                BiBernsteinPolynomial<F, R>    *pright_down,
                                                BiBernsteinPolynomial<F, R>    *pright_up,
                                                BiBernsteinPolynomial<F, R>    *pleft_up) const;

        void                                clipToInterval(
                                                R const                        &x0, 
                                                R const                        &x1,
                                                R const                        &y0,
                                                R const                        &y1,
                                                BiBernsteinPolynomial<F, R>    *pclip) const;
        /* elevate degree by (+r, +s) */
        void                                elevateDegree(
                                                uint32_t r_arg,
                                                uint32_t s_arg);

        /* static methods to precompute / free data needed for Bernstein-Bezier basis arithmetic
         * (binomial coefficients / inner products of basis polynomials). */
        static void                         initBernsteinPolyData();
        static void                         freeBernsteinPolyData();

        /* FIXME: add to PolyAlg::convertBasis set of polymorph conversion functions */
        void                                convertFromPowerBasis(Matrix<F> const &coeff_pow);
};

#include "../tsrc/BivariatePolynomial.impl.hh"

#endif
