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

#ifndef PARAMETRIC_CURVE_H
#define PARAMETRIC_CURVE_H

#include "common.hh"
#include "Tensor.hh"
#include "Vec3.hh"
#include "BoundingBox.hh"
#include "Polynomial.hh"

template <
    typename    C2F,
    typename    T  = double,
    typename    Tr = StaticVector<3u, double>
>
class ParametricCurve {
    protected:
        /* dimension, component functors (functions) and domain [t0, t1] */
        std::vector<C2F>        component_functors;
        T                       t0, t1;

        /* private check method for matching dimension / domain of other ParametricCurve of matching
         * template parameters (which are implicitly used inside the declaration) */
        void                    checkDimDomain(
                                    ParametricCurve const  &x,
                                    std::string const      &fn) const;

        void                    checkEvalParameter(
                                    T const                &x,
                                    std::string const      &fn) const;

    public:
        /* constructors, assignment operator */
                                ParametricCurve();
                                ParametricCurve(
                                    std::vector<C2F> const         &component_functors,
                                    T                               t0,
                                    T                               t1);
                                ParametricCurve(ParametricCurve<C2F, T, Tr> const &x);
        ParametricCurve        &operator=(ParametricCurve const &x);

        /* virtual destructor for polymorphism */
        virtual                ~ParametricCurve();

        /* evaluation: virtual, may be overridden by derived classes */
        virtual Tr              eval(const T &x) const;
        virtual Tr              eval_d(const T &x) const;
        virtual Tr              eval_d2(const T &x) const;

        /* get dimension */    
        uint32_t                getDim() const;

        /* get /set domain */
        std::pair<T, T>         getDomain() const;
        /*
        void                    setDomain(
                                    T const &t0,
                                    T const &t1);*/

        /* access component functors */
        C2F                    &operator[](uint32_t i);
        C2F const              &operator[](uint32_t i) const;

        /* arithmetic */
        ParametricCurve         operator+(ParametricCurve const &x) const;
        ParametricCurve        &operator+=(ParametricCurve const &x);

        ParametricCurve         operator-(ParametricCurve const &x) const;
        ParametricCurve        &operator-=(ParametricCurve const &x);

        ParametricCurve         operator*(T const &x) const;
        ParametricCurve        &operator*=(T const &x);

        ParametricCurve         operator/(T const &x) const;
        ParametricCurve        &operator/=(T const &x);
};

template<
    typename C2F,
    typename R = double
>
class SpaceCurveReal : public ParametricCurve<C2F, R, Vec3<R> >
{
    protected:
        bool                    arclen_set;
        R                       arclen;
        R                       arclen_dt;

    public:
        /* ctors, dtor */
                                SpaceCurveReal();
                                SpaceCurveReal(
                                    std::array<C2F, 3> const   &component_functors,
                                    R const                    &t0,
                                    R const                    &t1);
                                SpaceCurveReal(SpaceCurveReal<C2F, R> const &x);
        SpaceCurveReal<C2F, R> &operator=(SpaceCurveReal<C2F, R> const &x);

        virtual                ~SpaceCurveReal();

        /* arithmetic */
        SpaceCurveReal          operator+(SpaceCurveReal const &x) const;
        SpaceCurveReal         &operator+=(SpaceCurveReal const &x);

        SpaceCurveReal          operator-(SpaceCurveReal const &x) const;
        SpaceCurveReal         &operator-=(SpaceCurveReal const &x);

        SpaceCurveReal          operator*(R const &x) const;
        SpaceCurveReal         &operator*=(R const &x);

        SpaceCurveReal          operator/(R const &x) const;
        SpaceCurveReal         &operator/=(R const &x);

        /* methods specific to real space curves */

        /* getDerivative: return space curve obtained from first derivatives of component functors */
        SpaceCurveReal<C2F, R>  getDerivative() const;

        /* get frenet-serret and render frame */
        void                    getFrenetFrame(
                                    R const    &t,
                                    Vec3<R>    &x,
                                    Vec3<R>    &y,
                                    Vec3<R>    &z) const;

        void                    getRenderFrame(
                                    R const    &t,
                                    Vec3<R>     rvec,
                                    Vec3<R>     &x,
                                    Vec3<R>     &y,
                                    Vec3<R>     &z) const;
        
        /* approximate arc length */
        R                       approxArcLength(
                                   R const     &tstart, 
                                   R const     &tend,
                                   R const     &dt = 1E-4) const;
        
        /* update internally stored arc length of entire curve in domain [t0, t1] */
        void                    updateArcLength(R const &dt = 1E-4);
};

template<uint32_t degree, typename R = double>
class BezierCurve : public SpaceCurveReal< BernsteinPolynomial<degree, R, R>, R >
{
    private:
        typedef BezierCurve<degree, R> this_type;
        typedef BezierCurve<(degree>0) ? degree-1 : 0, R> this_deriv_type;
        static const uint32_t derivDeg = degree>0 ? (degree-1) : 0;
        static const uint32_t deriv2Deg = degree>1 ? (degree-2) : 0;
        typedef BernsteinPolynomial<degree, R, R> pol_type;
        typedef BernsteinPolynomial<derivDeg, R, R>  pol_deriv_type;
        typedef BernsteinPolynomial<deriv2Deg, R, R>  pol_deriv2_type;

        /* also STORE first and second derivative of polynomial component functions.
         * component_functors has already been declared in template base class ParametricCurve */
        std::vector<pol_deriv_type>  d_component_functors;
        std::vector<pol_deriv2_type>  d2_component_functors;

    public:
        /* ctors, dtor */
        BezierCurve();
        BezierCurve(const std::array<pol_type, 3>& component_functors);
        BezierCurve(const std::vector<Vec3<R> >& control_points);
        BezierCurve(const this_type& x);

        ~BezierCurve();

        this_type& operator=(const this_type& x);

        /* arithmetic */
        this_type operator+(const this_type& x) const;
        this_type& operator+=(const this_type& x);

        this_type operator-(const this_type& x) const;
        this_type& operator-=(const this_type& x);

        this_type operator*(const R& x) const;
        this_type& operator*=(const R& x);

        this_type operator/(const R& x) const;
        this_type& operator/=(const R& x);

        /* methods specific to BezierCurve */
        this_deriv_type getDerivative() const;

        std::list<Vec3<R> > getControlPoints() const;

        void split(const R& xt, this_type* cleft, this_type* cright) const;

        void clipToInterval(const R& t0, const R& t1, this_type* cclip) const;

        /* get bounding box using recursive subdivision up to a certain level and then the bounding box of all control
         * points */
        BoundingBox<R> getBoundingBox(uint32_t subdivision_detph = 0) const;

        void computeRegularityPolynomial(BernsteinPolynomial<2*derivDeg, R, R>& p_reg) const;
        void computeStationaryPointDistPoly(const Vec3<R>& x, BernsteinPolynomial<degree+derivDeg, R, R>& p) const;
};

#include "../tsrc/ParametricCurve.impl.hh"

#endif
