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

#ifndef POLY_ALGORITHMS
#define POLY_ALGORITHMS

#include "Polynomial.hh"
#include "BivariatePolynomial.hh"

namespace PolyAlg {

    /* polynomial base conversion algorithms */
    template <typename R>
    void
    convertBasis(
        BernsteinPolynomial<R, R>      &b,
        PowerPolynomial<R, R> const    &p);

    template <typename R>
    void
    convertBasis(
        PowerPolynomial<R, R>              &p,
        BernsteinPolynomial<R, R> const    &b);

    /* tensor multiplication of two univariate BernsteinPolynomials to produce a bivariate BiBernsteinPolynomial */
    template<
        typename F = double,
        typename R = double
    >
    BiBernsteinPolynomial<F, R>
    BernsteinTensorMultiply(
        BernsteinPolynomial<F, R> const &p,
        BernsteinPolynomial<F, R> const &q);

    /* convert univariate BernsteinPolynomial in BB(k) to BiBernsteinPolynomial where one variable is just
     * unity: if as_x == true, use the univariate polynomial variable as x, otherwise as y => result basis is BB(k, d)
     * and BB(d, k), respectively, where k = deg(p) */
    template <
        typename F = double,
        typename R = double
    >
    BiBernsteinPolynomial<F, R>
    BernsteinConvertToBiPoly(
        BernsteinPolynomial<F, R> const    &p,
        uint32_t                            d,
        bool                                as_x);

    /* create bernstein basis polynomials in Bernstein-Bezier basis, which have a matching "kronecker delta" vector for
     * coefficients. */
    template <
        typename F,
        typename R
    >
    BernsteinPolynomial<F, R>
    computeBernsteinBasisPoly(uint32_t n, uint32_t i);

    /* compute the convex hull the control polygon of a BernsteinPolynomial over the reals,
     * where the polynomial may be interpreted as as a two-dimensional bezier-curve as described in
     * the thesis. */
    template <typename R>
    void
    computeBezierControlPolyConvexHull(
        BernsteinPolynomial<R, R> const    &p,
        std::vector<Vec2>                  &cvhull,
        R const                            &eps_slope);

    /* classes to store roots */
    template <typename R>
    class RealInterval {
        public:
            R       t0, t1;

            RealInterval(
                R   t0,
                R   t1)
            {
                this->t0    = t0;
                this->t1    = t1;
            }

            R 
            midpoint() const
            {
                return ( (this->t0 + this->t1) / (R)2 );
            }
    };

    template <typename R>
    class RealRectangle {
        public:
            R       x0, x1;
            R       y0, y1;
            int     d;

            RealRectangle() {
                using Aux::Numbers::inf;

                this->x0    = -inf<R>();
                this->x1    = -inf<R>();
                this->y0    = inf<R>();
                this->y1    = inf<R>();
                d           = 0;
            }

            RealRectangle(
                R x0,
                R x1,
                R y0,
                R y1,
                int d = 1)
            {
                this->x0    = x0;
                this->x1    = x1;
                this->y0    = y0;
                this->y1    = y1;
                this->d     = d;
            }

            bool
            subset(const RealRectangle &b) const
            {
                return (this->x0 >= b.x0 && this->x1 <= b.x1 && this->y0 >= b.y0 && this->y1 <= b.y1);
            }

            std::array<R, 2>
            midpoint() const
            {
                return std::array<R, 2>( { (this->x0 + this->x1) / (R)2, (this->y0 + this->y1) / (R)2 } );
            }
    };


    template <typename R>
    struct ComplexInterval {
        R       x0, x1;
        R       y0, y1;
        int     d;

        ComplexInterval()
        {
        }

        ComplexInterval(
            R x0,
            R x1,
            R y0,
            R y1,
            int d = 1)
        {
            this->x0    = x0;
            this->x1    = x1;
            this->y0    = y0;
            this->y1    = y1;
            this->d     = d;
        }
    };

    void
    initPolyAlgorithmData();
    
    void
    freePolyAlgorithmData();

    template <typename R = double>
    void
    BezClip_roots(
            BernsteinPolynomial<R, R> const    &pinput,
            R const                            &alpha,
            R const                            &beta,
            R const                            &tol,
            std::vector<RealInterval<R>>       &roots,
            R const                            &eps             = 1E-11,
            R const                            &eps_slope       = 1E-8);

    template <typename R>
    void
    BiLinClip_getApproximationData(
        uint32_t                                    m,
        uint32_t                                    n,
        bool                                        dynamic_recomputation_arg   = true,
        BiBernsteinPolynomial<R, R>                *BiLinClip_L00               = NULL,
        BiBernsteinPolynomial<R, R>                *BiLinClip_L10               = NULL,
        BiBernsteinPolynomial<R, R>                *BiLinClip_L01               = NULL,
        Matrix<R>                                  *BiLinClip_A00               = NULL,
        Matrix<R>                                  *BiLinClip_A10               = NULL,
        Matrix<R>                                  *BiLinClip_A01               = NULL);

    template <typename R = double>
    void
    BiLinClip_roots(
        BiBernsteinPolynomial<R, R> const  &pinput,
        BiBernsteinPolynomial<R, R> const  &qinput,
        R const                            &alpha0,
        R const                            &alpha1,
        R const                            &beta0,
        R const                            &beta1,
        R const                            &tol,
        std::vector<RealRectangle<R>>      &roots,
        bool                                approx_data_dynamic_recomputation   = false,
        bool                                use_blacklist                       = false,
        std::vector<RealRectangle<R>>      *blacklist                           = NULL,
        R const                            &eps                                 = 1E-11,
        R const                            &linsolve_eps                        = 1E-11,
        R const                            &rec_armax                           = 1E5);

    /*
    template <typename R = double>
    void
    TDB_roots(
            Polynomial<R, R>                   &p,
            R const                            &x0,
            R const                            &x1,
            R const                            &y0,
            R const                            &y1,
            uint32_t                            nsamp_per_edge,
            R const                            &tol,
            std::vector<ComplexInterval<R>>    &roots);
    */

    #include "../tsrc/PolyAlgorithms.impl.hh"
}

#endif