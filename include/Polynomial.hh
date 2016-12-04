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

/*! @file Polynomial.hh
 *  @brief common header file for all classes implementing polynomials (including arithmetic,
 *  conversion, ..) in different representations */
#ifndef POLYNOMIAL_H
#define POLYNOMIAL_H

#include "Tensor.hh"
#include "Vector.hh"
#include "Vec3.hh"
#include "aux.hh"

/*! @brief abstract univariate polynomial base class providing minimal interface */
template <
    typename F,
    typename R = double
>
class Polynomial {
    protected:
        uint32_t                            degree;
        Vector<F>                           coeff;

        void                                checkDegree(
                                                Polynomial<F, R> const &q,
                                                std::string const      &fn) const;

    public:
                                            Polynomial();
                                            Polynomial(uint32_t degree, F const &x);
                                            Polynomial(Vector<F> const &coeff);
                                            Polynomial(Polynomial<F, R> const &q);
        Polynomial<F, R>                   &operator=(Polynomial<F, R> const &q);
        virtual                            ~Polynomial();

        uint32_t                            getDegree() const;

        /* coefficient access, both as reference (non-const access) and value (const access).
           since operator[] cannot be overloaded to take more than one argument, operator() is used
           for bivariate and higher dimensional cases.  for completeness, it has been defined here
           as well, so operator[] and operator() behave identical. */
        F                                   operator()(uint32_t cidx) const;
        F                                  &operator()(uint32_t cidx);

        F                                   operator[](uint32_t cidx) const;
        F                                  &operator[](uint32_t cidx);

        Vector<F>                           getCoeffs() const;
        virtual void                        setCoeffs(Vector<F> const &coeff);
        F                                   getMaxAbsCoeff() const;

        virtual void                        initConstant(F const &x) = 0;
        virtual void                        zero();

        /* virtual evaluation for function and its first two derivatives are mandatory */
        virtual F                           eval(   R const &x) const               = 0;
        virtual F                           eval_d( R const &x) const               = 0;
        virtual F                           eval_d2(R const &x) const               = 0;
        //virtual F                           eval_dn(F const &x, uint32_t n) const   = 0;

        /* arithmetic: make this protected to prevent direct usage of the following operators
         * through upcasting of different derived polynomial pointers or references to abstract type
         * Polynomial. however, the implementations can and should be used in derived classes. */
    protected:
        Polynomial<F, R>                   &operator+=(Polynomial<F, R> const &q);
        Polynomial<F, R>                   &operator-=(Polynomial<F, R> const &q);

        Polynomial<F, R>                   &operator*=(F const &x);
        Polynomial<F, R>                   &operator/=(F const &x);

        /* scalar product. since this heavily depends on the scalar product (e.g. L2), and representation used (e.g.
         * concerning expressions for integrals), derived class may use this as a template for their more specific
         * implmentations.  */
        F                                   operator*(Polynomial<F, R> const &q) const;

    public:

        /* print coefficient vector */
        void                                printCoeff() const;

        /* write plotfile suitable for gnuplot */
        void                                writePlotFile(
                                                R               t0,
                                                R               t1,
                                                uint32_t        ticks,
                                                std::string     filename) const;
};


#define POWERPOLY_BASIS_IP_DEFAULT_SIZE 128u
/*! @brief class for univariate polynomials represented in power (aka monomial) basis */
template<
    typename F,
    typename R = double
>
class PowerPolynomial : public Polynomial<F, R> {
    private:
        static bool                         power_basis_inner_products_mutable;
        static Tensor<2, F>                 power_basis_inner_products;
        static uint32_t                     power_basis_inner_products_max_dim;

        static F                            computePowerBasisInnerProduct(
                                                uint32_t i,  
                                                uint32_t j);
    public:
                                            PowerPolynomial();
                                            PowerPolynomial(uint32_t degree, F const &x);
                                            PowerPolynomial(Vector<F> const &coeff);
                                            PowerPolynomial(PowerPolynomial<F, R> const &q);
        PowerPolynomial<F, R>              &operator=(PowerPolynomial const &q);
        virtual                            ~PowerPolynomial();

        static void                         initPowerBasisInnerProducts(uint32_t max_dim = POWERPOLY_BASIS_IP_DEFAULT_SIZE);
        static F                            getPowerBasisInnerProduct(
                                                uint32_t i,  
                                                uint32_t j);

        void                                setCoeffs(Vector<F> const &coeff);
        void                                initConstant(F const &x);

        F                                   eval(   R const &x) const;
        F                                   eval_d( R const &x) const;
        F                                   eval_d2(R const &x) const;

        /* arithmetic */

        /* addition in vector space \Pi^n */
        PowerPolynomial<F, R>               operator+(PowerPolynomial<F, R> const &q) const;
        PowerPolynomial<F, R>              &operator+=(PowerPolynomial<F, R> const &q);

        PowerPolynomial<F, R>               operator-(PowerPolynomial<F, R> const &q) const;
        PowerPolynomial<F, R>              &operator-=(PowerPolynomial<F, R> const &q);

        /* scalar multiplication in vector space \Pi^n */
        PowerPolynomial<F, R>               operator*(F const &x) const;
        PowerPolynomial<F, R>              &operator*=(F const &x);

        PowerPolynomial<F, R>               operator/(F const &x) const;
        PowerPolynomial<F, R>              &operator/=(F const &x);

        /* inner product in \Pi^n with respect to L2 norm \int_0^1{p(t)q(t) dt} */
        F                                   operator*(PowerPolynomial<F, R> const &q) const;

        /* multiplication as polynomials over |R or a suitably defined ring. */
        PowerPolynomial<F, R>               multiply(PowerPolynomial<F, R> const &q) const;
        PowerPolynomial<F, R>               square() const;

        PowerPolynomial<F, R>               getDerivative() const;

        void                                elevateDegree(uint32_t r); 

        /* static function that performs minimum elevation on either p or q such that their
         * representation degree matches, i.e. both p and q are represented in \Pi^n with
         * coefficient vectors of equal size (n + 1). */
        static void                         matchDegree(
                                                PowerPolynomial<F, R>  &p,
                                                PowerPolynomial<F, R>  &q);

        //void                                convertFromBernsteinBasis(Vector<F> const &coeff_arg);
};

/* define static members of template class PowerPolynomial<F, R> */
template<typename F, typename R> bool           PowerPolynomial<F, R>::power_basis_inner_products_mutable = true;
template<typename F, typename R> Tensor<2, F>   PowerPolynomial<F, R>::power_basis_inner_products( {1, 1} );
template<typename F, typename R> uint32_t       PowerPolynomial<F, R>::power_basis_inner_products_max_dim = 0;

#define BERNSTEINPOLY_BASIS_IP_DEFAULT_SIZE 32u
/*! @brief class for univariate polynomials represented in Bernstein-Bezier basis */
template<
    typename F,
    typename R = double
>
class BernsteinPolynomial : public Polynomial<F, R> {
    private:
        static bool                         bernstein_basis_inner_products_mutable;
        static Tensor<4, F>                 bernstein_basis_inner_products;
        static uint32_t                     bernstein_basis_inner_products_max_dim;

        static F                            computeBernsteinBasisInnerProduct(
                                                uint32_t m,  
                                                uint32_t n,  
                                                uint32_t i,  
                                                uint32_t j);

    public:
                                            BernsteinPolynomial();
                                            BernsteinPolynomial(uint32_t degree, F const &x);
                                            BernsteinPolynomial(Vector<F> const &coeff);
                                            BernsteinPolynomial(BernsteinPolynomial<F, R> const &q);
        BernsteinPolynomial<F, R>          &operator=(BernsteinPolynomial<F, R> const &q);
        virtual                            ~BernsteinPolynomial();

        static void                         initBernsteinBasisInnerProducts(uint32_t max_dim = POWERPOLY_BASIS_IP_DEFAULT_SIZE);
        static F                            getBernsteinBasisInnerProduct(
                                                uint32_t m,  
                                                uint32_t n,  
                                                uint32_t i,  
                                                uint32_t j);

        static void                         setInnerProductDataMutable();
        static void                         setInnerProductDataImmutable();

        void                                setCoeffs(Vector<F> const &coeff);
        void                                initConstant(F const &x);

        F                                   eval(   R const &x) const;
        F                                   eval_d( R const &x) const;
        F                                   eval_d2(R const &x) const;

        /* arithmetic */

        /* addition in vector space \Pi^n */
        BernsteinPolynomial<F, R>           operator+(BernsteinPolynomial<F, R> const &q) const;
        BernsteinPolynomial<F, R>          &operator+=(BernsteinPolynomial<F, R> const &q);

        BernsteinPolynomial<F, R>           operator-(BernsteinPolynomial<F, R> const &q) const;
        BernsteinPolynomial<F, R>          &operator-=(BernsteinPolynomial<F, R> const &q);

        /* scalar multiplication in vector space \Pi^n */
        BernsteinPolynomial<F, R>           operator*(F const &x) const;
        BernsteinPolynomial<F, R>          &operator*=(F const &x);

        BernsteinPolynomial<F, R>           operator/(F const &x) const;
        BernsteinPolynomial<F, R>          &operator/=(F const &x);

        /* inner product in \Pi^n with respect to L2 norm \int_0^1{p(t)q(t) dt} */
        F                                   operator*(BernsteinPolynomial<F, R> const &q) const;

        /* multiplication as polynomials over |R or a suitably defined ring */
        BernsteinPolynomial<F, R>           multiply(BernsteinPolynomial<F, R> const &q) const;
        BernsteinPolynomial<F, R>           square() const;

        /* derivative */
        BernsteinPolynomial<F, R>           getDerivative() const;

        /* degree elevation by r, effectively like multiplication with 1(t) in BB(r). this is done
         * in-place, though */
        void                                elevateDegree(uint32_t r_arg); 

        /* static function that performs minimum elevation on either p or q such that their
         * representation degree matches, i.e. both p and q are represented in BB(n) for some n
         * after the call. */
        static void                         matchDegree(
                                                BernsteinPolynomial<F, R>  &p,
                                                BernsteinPolynomial<F, R>  &q);

        void                                split(
                                                R const                    &t,
                                                BernsteinPolynomial<F, R>  *pleft,
                                                BernsteinPolynomial<F, R>  *pright) const;

        void                                clipToInterval(
                                                R const                    &t0,
                                                R const                    &t1,
                                                BernsteinPolynomial<F, R>  *pclip) const;
};

/* define static members of template class BernsteinPolynomial<F, R> */
template<typename F, typename R> bool           BernsteinPolynomial<F, R>::bernstein_basis_inner_products_mutable = true;
template<typename F, typename R> Tensor<4, F>   BernsteinPolynomial<F, R>::bernstein_basis_inner_products({1, 1, 1, 1});
template<typename F, typename R> uint32_t       BernsteinPolynomial<F, R>::bernstein_basis_inner_products_max_dim = 0;

#if 0
template <
    typename F,
    typename R = double
>
class LegendrePolynomial {
    private:
        /* custom init method, since legendre polynomials integrals do not depend on the dimension n
         * of \Pi^n and furthermore, legendre polynomials are orthogonal, i.e. the 4-array reduces
         * to a 1-array of length n for this case => reimplement initBasisInnerProducts efficiently
         * for this case. */
        F                                  *legendre_basis_inner_products;

        /* override generic virtual Polynomial<F>::initBasisInnerProducts() */
        void                                initBasisInnerProducts();
        
        /* implement pure virtual Polynomial<F>::computeBasisInnerProduct */
        F                                   computeBasisInnerProduct(
                                                uint32_t m,  
                                                uint32_t n,  
                                                uint32_t i,  
                                                uint32_t j);
    public:
                                            LegendrePolynomial();
                                            LegendrePolynomial(uint32_t degree, F const &x);
                                            LegendrePolynomial(Vector<F>      const &coeff);
                                            LegendrePolynomial(LegendrePolynomial<F, R> const &q);
        LegendrePolynomial<F, R>           &operator=(LegendrePolynomial<F, R> const &q);
        virtual                            ~LegendrePolynomial();


        F                                   eval(   R const &x) const;
        F                                   eval_d( R const &x) const;
        F                                   eval_d2(R const &x) const;

        /* arithmetic */

        /* addition in vector space \Pi^n */
        LegendrePolynomial<F, R>            operator+(LegendrePolynomial<F, R> const &q) const;
        LegendrePolynomial<F, R>           &operator+=(LegendrePolynomial<F, R> const &q);

        LegendrePolynomial<F, R>            operator-(LegendrePolynomial<F, R> const &q) const;
        LegendrePolynomial<F, R>           &operator-=(LegendrePolynomial<F, R> const &q);

        /* scalar multiplication in vector space \Pi^n */
        LegendrePolynomial<F, R>            operator*(F const &x) const;
        LegendrePolynomial<F, R>           &operator*=(F const &x);

        LegendrePolynomial<F, R>            operator/(F const &x) const;
        LegendrePolynomial<F, R>           &operator/=(F const &x);

        /* inner product in \Pi^n with respect to L2 norm \int_0^1{p(t)q(t) dt} */
        F                                   operator*(LegendrePolynomial<F, R> const &q) const;

        /* multiplication as polynomials over |R or a suitably defined ring */
        LegendrePolynomial<F, R>            multiply(LegendrePolynomial<F, R> const &q) const;
        LegendrePolynomial<F, R>            square() const;
};
#endif

#include "../tsrc/Polynomial.impl.hh"

#endif