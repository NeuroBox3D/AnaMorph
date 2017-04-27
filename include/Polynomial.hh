/*
 * This file is part of
 *
 * AnaMorph: a framework for geometric modelling, consistency analysis and surface
 * mesh generation of anatomically reconstructed neuron morphologies.
 * 
 * Copyright (c) 2013-2017: G-CSC, Goethe University Frankfurt - Queisser group
 * Author: Konstantin Mörschel
 * 
 * AnaMorph is free software: Redistribution and use in source and binary forms,
 * with or without modification, are permitted under the terms of the
 * GNU Lesser General Public License version 3 (as published by the
 * Free Software Foundation) with the following additional attribution
 * requirements (according to LGPL/GPL v3 §7):
 * 
 * (1) The following notice must be displayed in the Appropriate Legal Notices
 * of covered and combined works:
 * "Based on AnaMorph (https://github.com/NeuroBox3D/AnaMorph)."
 * 
 * (2) The following notice must be displayed at a prominent place in the
 * terminal output of covered works:
 * "Based on AnaMorph (https://github.com/NeuroBox3D/AnaMorph)."
 *
 * (3) Neither the name "AnaMorph" nor the names of its contributors may be
 * used to endorse or promote products derived from this software without
 * specific prior written permission.
 *
 * (4) The following bibliography is recommended for citation and must be
 * preserved in all covered files:
 * "Mörschel K, Breit M, Queisser G. Generating neuron geometries for detailed
 *   three-dimensional simulations using AnaMorph. Neuroinformatics (2017)"
 * "Grein S, Stepniewski M, Reiter S, Knodel MM, Queisser G.
 *   1D-3D hybrid modelling – from multi-compartment models to full resolution
 *   models in space and time. Frontiers in Neuroinformatics 8, 68 (2014)"
 * "Breit M, Stepniewski M, Grein S, Gottmann P, Reinhardt L, Queisser G.
 *   Anatomically detailed and large-scale simulations studying synapse loss
 *   and synchrony using NeuroBox. Frontiers in Neuroanatomy 10 (2016)"
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU Lesser General Public License for more details.
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

/*! @file Polynomial.hh
 *  @brief common header file for all classes implementing polynomials (including arithmetic,
 *  conversion, ..) in different representations */
#ifndef POLYNOMIAL_H
#define POLYNOMIAL_H

#include "StaticVector.hh"
#include "Vec3.hh"
#include "aux.hh"


// This implementation has deplorable runtime speed due to the implementation
// of Vector which makes too many dynamic allocations.
#if 0
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
        static bool                         bernstein_basis_inner_products_initialized;
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
template<typename F, typename R> bool           BernsteinPolynomial<F, R>::bernstein_basis_inner_products_initialized = true;
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

#endif



/*! @brief abstract univariate polynomial base class providing minimal interface */
template <uint32_t degree, typename F, typename R = double>
class Polynomial
{
    protected:
        StaticVector<degree+1, F> coeff;

    public:
        typedef Polynomial<degree, F, R> this_type;
        typedef StaticVector<degree+1, F> coeff_type;

        Polynomial();
        Polynomial(const F& x);
        Polynomial(const coeff_type& coeff);
        Polynomial(const this_type& q);

        virtual ~Polynomial();

        this_type& operator=(const this_type& q);

        uint32_t getDegree() const;

        // coefficient access, both as reference (non-const access) and value (const access).
        // since operator[] cannot be overloaded to take more than one argument, operator() is used
        // for bivariate and higher dimensional cases.  for completeness, it has been defined here
        // as well, so operator[] and operator() behave identical.
        F operator()(uint32_t cidx) const;
        F& operator()(uint32_t cidx);

        F operator[](uint32_t cidx) const;
        F& operator[](uint32_t cidx);

        const coeff_type& getCoeffs() const;
        coeff_type& getCoeffs();
        virtual void setCoeffs(const coeff_type& coeff);
        F getMaxAbsCoeff() const;

        virtual void initConstant(const F& x) = 0;
        virtual void zero();

        // virtual evaluation for function and its first two derivatives are mandatory
        virtual F eval(const R& x) const = 0;
        virtual F eval_d(const R& x) const = 0;
        virtual F eval_d2(const R& x) const = 0;

        // arithmetic: make this protected to prevent direct usage of the following operators
        // through upcasting of different derived polynomial pointers or references to abstract type
        // Polynomial. however, the implementations can and should be used in derived classes.
    protected:
        this_type& operator+=(const this_type& q);
        this_type& operator-=(const this_type& q);

        this_type& operator*=(const F& x);
        this_type& operator/=(const F& x);

        // scalar product. since this heavily depends on the scalar product (e.g. L2), and representation used (e.g.
        // concerning expressions for integrals), derived class may use this as a template for their more specific
        // implmentations.
        F operator*(const this_type& q) const;

    public:
        // print coefficient vector
        void printCoeff() const;

        // write plotfile suitable for gnuplot
        void writePlotFile(R t0, R t1, uint32_t ticks, const std::string& filename) const;

    protected:
        template <typename TF, typename TR, typename dummy = void>
        struct PrintCoeffImpl
        {
            PrintCoeffImpl(const Polynomial<degree, F, R>& p);
        };
        template <typename dummy>
        struct PrintCoeffImpl<float, float, dummy>
        {
            PrintCoeffImpl(const Polynomial<degree, F, R>& p);
        };
        template <typename dummy>
        struct PrintCoeffImpl<double, double, dummy>
        {
            PrintCoeffImpl(const Polynomial<degree, F, R>& p);
        };

        template <typename TF, typename TR, typename dummy = void>
        struct WritePlotFileImpl
        {
            WritePlotFileImpl(R t0, R t1, uint32_t ticks, const std::string& filename,
                              const Polynomial<degree, F, R>& p);
        };
        template <typename dummy>
        struct WritePlotFileImpl<float, float, dummy>
        {
            WritePlotFileImpl(R t0, R t1, uint32_t ticks, const std::string& filename,
                              const Polynomial<degree, F, R>& p);
        };
        template <typename dummy>
        struct WritePlotFileImpl<double, double, dummy>
        {
            WritePlotFileImpl(R t0, R t1, uint32_t ticks, const std::string& filename,
                              const Polynomial<degree, F, R>& p);
        };

};


#define POWERPOLY_BASIS_IP_DEFAULT_SIZE 128u
/*! @brief class for univariate polynomials represented in power (aka monomial) basis */
template <uint32_t degree, typename F, typename R = double>
class PowerPolynomial
: public Polynomial<degree, F, R>
{
    public:
        static void initPowerBasisInnerProducts();
        static F getPowerBasisInnerProduct(uint32_t i, uint32_t j);

    private:
        static StaticMatrix<degree, degree, F> power_basis_inner_products;
        static F computePowerBasisInnerProduct(uint32_t i, uint32_t j);

    public:
        typedef PowerPolynomial<degree, F, R> this_type;
        typedef Polynomial<degree, F, R> base_type;
        typedef PowerPolynomial<2*degree, F, R> squared_type;
        typedef StaticVector<degree+1, F> coeff_type;

        PowerPolynomial();
        PowerPolynomial(const F& x);
        PowerPolynomial(const coeff_type& coeff);
        PowerPolynomial(const this_type& q);

        virtual ~PowerPolynomial();

        this_type& operator=(const this_type& q);
        void initConstant(const F& x);

        F eval(const R& x) const;
        F eval_d(const R& x) const;
        F eval_d2(const R& x) const;

        // addition in vector space \Pi^n
        this_type operator+(const this_type& q) const;
        this_type& operator+=(const this_type& q);

        this_type operator-(const this_type& q) const;
        this_type& operator-=(const this_type& q);

        // scalar multiplication in vector space \Pi^n
        this_type operator*(const F& x) const;
        this_type& operator*=(const F& x);

        this_type operator/(const F& x) const;
        this_type& operator/=(const F& x);

        // inner product in \Pi^n with respect to L2 norm \int_0^1{p(t)q(t) dt}
        F operator*(const this_type& q) const;

        // multiplication as polynomials over |R or a suitably defined ring.
        template <uint32_t deg>
        PowerPolynomial<degree+deg, F, R> multiply(const PowerPolynomial<deg, F, R>& q) const;
        squared_type square() const;

        PowerPolynomial<(degree>0) ? degree-1 : 0, F, R> getDerivative() const;

        template <uint32_t deg>
        PowerPolynomial<degree+deg, F, R> elevateDegree();

    private:
        using base_type::coeff;
};

// define static members of template class PowerPolynomial<degree, F, R>
template <uint32_t degree, typename F, typename R>
StaticMatrix<degree, degree, F> PowerPolynomial<degree, F, R>::power_basis_inner_products;



#define BERNSTEINPOLY_BASIS_IP_DEFAULT_SIZE 32u
/*! @brief class for univariate polynomials represented in Bernstein-Bezier basis */
template<uint32_t degree, typename F, typename R = double>
class BernsteinPolynomial : public Polynomial<degree, F, R>
{
    public:
        static void initBernsteinBasisInnerProducts();
        static F getBernsteinBasisInnerProduct(uint32_t i, uint32_t j);

    private:
        // inner products of i-th and j-th basis functions of degree
        static StaticMatrix<degree+1, degree+1, F> bernstein_basis_inner_products;
        static F computeBernsteinBasisInnerProduct(uint32_t i, uint32_t j);

    public:
        typedef BernsteinPolynomial<degree, F, R> this_type;
        typedef Polynomial<degree, F, R> base_type;
        typedef BernsteinPolynomial<2*degree, F, R> squared_type;
        typedef StaticVector<degree+1, F> coeff_type;

        BernsteinPolynomial();
        BernsteinPolynomial(const F& x);
        BernsteinPolynomial(const coeff_type& coeff);
        BernsteinPolynomial(const this_type& q);

        virtual ~BernsteinPolynomial();

        this_type& operator=(const this_type& q);

        void setCoeffs(const coeff_type& coeff);
        void initConstant(const F& x);

        F eval(const R& x) const;
        F eval_d(const R& x) const;
        F eval_d2(const R& x) const;

        // addition in vector space \Pi^n
        this_type operator+(const this_type& q) const;
        this_type& operator+=(const this_type& q);

        this_type operator-(const this_type& q) const;
        this_type& operator-=(const this_type& q);

        // scalar multiplication in vector space \Pi^n
        this_type operator*(const F& x) const;
        this_type& operator*=(const F& x);

        this_type operator/(const F& x) const;
        this_type& operator/=(const F& x);

        // inner product in \Pi^n with respect to L2 norm \int_0^1{p(t)q(t) dt}
        F operator*(const this_type& q) const;

        // multiplication as polynomials over |R or a suitably defined ring
        template <uint32_t deg>
        BernsteinPolynomial<degree+deg, F, R> multiply(const BernsteinPolynomial<deg, F, R>& q) const;
        squared_type square() const;

        // derivative
        BernsteinPolynomial<(degree>0) ? degree-1 : 0, F, R> getDerivative() const;

        // degree elevation by r, effectively like multiplication with 1(t) in BB(r).
        template <uint32_t deg>
        BernsteinPolynomial<deg, F, R> elevateDegree();

        void split(const R& t, this_type* pleft, this_type* pright) const;

        void clipToInterval(const R& t0, const R& t1, this_type* pclip) const;

    private:
        using base_type::coeff;
};

// define static members of template class BernsteinPolynomial<F, R>
template<uint32_t degree, typename F, typename R>
StaticMatrix<degree+1, degree+1, F> BernsteinPolynomial<degree, F, R>::bernstein_basis_inner_products;


#include "../tsrc/Polynomial_impl.hh"

#endif
