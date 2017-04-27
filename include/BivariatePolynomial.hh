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

#include "../tsrc/BivariatePolynomial_impl.hh"

#endif
