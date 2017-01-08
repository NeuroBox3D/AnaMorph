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

#include "StaticMatrix.hh"

/*! @brief abstract bivariate polynomial base class providing minimal interface */
template <uint32_t deg1, uint32_t deg2, typename F, typename R>
class BivariatePolynomial
{
    protected:
        StaticMatrix<deg1+1, deg2+1, F> coeff;

    public:
        typedef BivariatePolynomial<deg1, deg2, F, R> this_type;
        typedef StaticMatrix<deg1+1, deg2+1, F> coeff_type;

        BivariatePolynomial();
        BivariatePolynomial(const F& x);
        BivariatePolynomial(const this_type& q);

        virtual ~BivariatePolynomial();

        this_type& operator=(const this_type& q);

        // substitute for operator[], that unfortunately cannot be overloaded to have two arguments
        virtual F& operator()(uint32_t i, uint32_t j);
        virtual F operator()(uint32_t i, uint32_t j) const;

        const coeff_type& getCoeffs() const;
        coeff_type& getCoeffs();
        void setCoeffs(const coeff_type& coeff);
        F getMaxAbsCoeff() const;

        void initConstant(const F& x);
        void zero();

        virtual F eval(const R& x, const R& y) const = 0;

        // arithmetic: make this protected just as for Polynomial to prevent misuse with two upcast pointers /
        // references from different derived classes
    protected:
        // addition in vector space \Pi^{m, n}
        this_type& operator+=(const this_type& q);
        this_type& operator-=(const this_type& q);

        // scalar multiplication in vector space \Pi^{m, n}
        this_type& operator*=(const F& x);
        this_type& operator/=(const F& x);

        // inner product with respect to L2 norm \int_0^1\int_0^1{p(x,y)q(x,y) dxdy}
        //F                                   operator*(BivariatePolynomial<F, R> const &q) const;

    public:
        // print coefficient matrix
        void printCoeff() const;

        // write plotfile suitable for gnuplot
        void writePlotFile(uint32_t ticks, const std::string& filename) const;

    protected:
        template <typename TF, typename TR, typename dummy = void>
        struct PrintCoeffImpl
        {
            PrintCoeffImpl(const BivariatePolynomial<deg1, deg2, F, R>& p);
        };
        template <typename dummy>
        struct PrintCoeffImpl<float, float, dummy>
        {
            PrintCoeffImpl(const BivariatePolynomial<deg1, deg2, F, R>& p);
        };
        template <typename dummy>
        struct PrintCoeffImpl<double, double, dummy>
        {
            PrintCoeffImpl(const BivariatePolynomial<deg1, deg2, F, R>& p);
        };

        template <typename TF, typename TR, typename dummy = void>
        struct WritePlotFileImpl
        {
            WritePlotFileImpl(uint32_t ticks, const std::string& filename,
                              const BivariatePolynomial<deg1, deg2, F, R>& p);
        };
        template <typename dummy>
        struct WritePlotFileImpl<float, float, dummy>
        {
            WritePlotFileImpl(uint32_t ticks, const std::string& filename,
                              const BivariatePolynomial<deg1, deg2, F, R>& p);
        };
        template <typename dummy>
        struct WritePlotFileImpl<double, double, dummy>
        {
            WritePlotFileImpl(uint32_t ticks, const std::string& filename,
                              const BivariatePolynomial<deg1, deg2, F, R>& p);
        };
};

/*! @brief class for bivariate polynomials represented in Bernstein-Bezier basis */
template <uint32_t deg1, uint32_t deg2, typename F, typename R>
class BiBernsteinPolynomial : public BivariatePolynomial<deg1, deg2, F, R>
{
    private:
        using BivariatePolynomial<deg1, deg2, F, R>::coeff;

    public:
        typedef BiBernsteinPolynomial<deg1, deg2, F, R> this_type;
        typedef StaticMatrix<deg1+1, deg2+1, F> coeff_type;

        BiBernsteinPolynomial();
        BiBernsteinPolynomial(const F& x);
        BiBernsteinPolynomial(const coeff_type& coeff);
        BiBernsteinPolynomial(const this_type& q);

        virtual ~BiBernsteinPolynomial();

        this_type& operator=(const this_type& q);

        // implementation of BivariatePolynomial interface
        F eval(const R& x, const R& y) const;

        // addition in vector space \Pi^{m, n} */
        this_type operator+(const this_type& q) const;
        this_type& operator+=(const this_type& q);

        this_type operator-(const this_type& q) const;
        this_type& operator-=(const this_type& q);

        // scalar multiplication in vector space \Pi^{m, n} */
        this_type operator*(const F& x) const;
        this_type& operator*=(const F& x);

        this_type operator/(const F& x) const;
        this_type& operator/=(const F& x);

        // inner product with respect to L2 norm \int_0^1\int_0^1{p(x,y)q(x,y) dxdy}
        template <uint32_t d1, uint32_t d2>
        F operator*(const BiBernsteinPolynomial<d1, d2, F, R>& q) const;

        // multiplication as polynomials over |R^2 or a suitably defined ring */
        template <uint32_t d1, uint32_t d2>
        BiBernsteinPolynomial<deg1+d1, deg2+d2, F, R> multiply(const BiBernsteinPolynomial<d1, d2, F, R>& q) const;

        BiBernsteinPolynomial<2*deg1, 2*deg2, F, R> square() const;

        // other methods
        void split_x(const R& x,this_type* pleft, this_type* pright) const;
        void split_y(const R& y, this_type* pdown, this_type* pup) const;
        void split_xy
        (
            const R& x,
            const R& y,
            this_type* pleft_down,
            this_type* pright_down,
            this_type* pright_up,
            this_type* pleft_up
        ) const;

        void clipToInterval(const R& x0, const R& x1, const R& y0, const R& y1, this_type* pclip) const;

        // elevate degree by (+d1, +d2)
        template <uint32_t d1, uint32_t d2>
        BiBernsteinPolynomial<deg1+d1, deg2+d2, F, R> elevateDegree();

        // static methods to precompute / free data needed for Bernstein-Bezier basis arithmetic
        // (binomial coefficients / inner products of basis polynomials).
        //static void initBernsteinPolyData();
        //static void freeBernsteinPolyData();

        // FIXME: add to PolyAlg::convertBasis set of polymorph conversion functions
        void convertFromPowerBasis(const coeff_type& coeff_pow);
};

#include "../tsrc/BivariatePolynomial.impl.hh"

#endif
