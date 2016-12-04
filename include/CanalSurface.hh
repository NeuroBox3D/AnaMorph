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

/*! @file CanalSurface.h
 *  @brief header file of class CanalSurface
 *
 *  detailed CanalSurface.h description.
 * */
#ifndef CANAL_SURFACE_H
#define CANAL_SURFACE_H

#include "Vec3.hh"
#include "BoundingBox.hh"
#include "Mesh.hh"
#include "ParametricCurve.hh"

template<
    typename C2F,
    typename R
>
class ConstantRadiusFunctor {
    private:
        R                               radius;

    public:
                ConstantRadiusFunctor()
                {
                    this->radius = 0;
                }

                /* NOTE: implicitly generated copy ctor and assignment operator suffice here */

        void    setRadius(R const &r)
                {
                    this->radius = r;
                }

        void    clipToInterval(
                    R const                        &t0,
                    R const                        &t1,
                    SpaceCurveReal<C2F, R> const   &gamma) const
                {
                }

        R       getMaxRadius() const
                {
                    return (this->radius);
                }
         
        R       operator()(
                    R const                        &t,
                    SpaceCurveReal<C2F, R> const   &gamma) const
                {
                    return (this->radius);
                }
};

template<
    typename C2F,
    typename R
>
class LinearRadiusInterpolatorArcLen {
    private:
        R       r0, r1;
        R       arclen_dt;

    public:
                LinearRadiusInterpolatorArcLen()
                {
                    this->r0        = 0;
                    this->r1        = 0;
                    this->arclen_dt = 1E-2;
                }

                LinearRadiusInterpolatorArcLen(
                    R const    &r0,
                    R const    &r1)
                {
                    this->r0        = r0;
                    this->r1        = r1;
                    this->arclen_dt = 1E-2;
                }

                /* NOTE: implicitly generated copy ctor and assignment operator suffice here */

        std::pair<R, R>
                getRadii() const
                {
                    return (std::pair<R, R>(this->r0, this->r1));
                }

        R       getMaxRadius() const
                {
                    return (std::max(this->r0, this->r1));
                }

        R       getMinRadius() const
                {
                    return (std::min(this->r0, this->r1));
                }

        void    setRadii(
                    R const                        &r0,
                    R const                        &r1)
                {
                    this->r0    = r0;
                    this->r1    = r1;
                }

        void    clipToInterval(
                    R const                        &t0,
                    R const                        &t1,
                    SpaceCurveReal<C2F, R> const   &gamma) const
                {
                    if (t0 <= t1) {
                        this->r0 = gamma.radiusEval(t0);
                        this->r1 = gamma.radiusEval(t1);
                    }
                    else {
                        throw("LinearRadiusInterpolatorDomain::clipToInterval(): malformed interval [t0, t1]: t0 > t1.");
                    }
                }

        R       operator()(
                    R const                        &t,
                    SpaceCurveReal<C2F, R> const   &gamma) const
                {
                    /* linearly interpolate with respect to arc length */
                    R t0, t1, arclen_t0t, arclen_t0t1, ratio;
                    auto domain = gamma.getDomain();
                    t0          = domain.first;
                    t1          = domain.second;

                    if (t < t0 || t > t1) {
                        throw("LinearRadiusInterpolatorDomain::operator(): given parameter value t not in domain [t0, t1].");
                    }

                    arclen_t0t  = gamma.approxArcLength(t0, t, this->arclen_dt);
                    arclen_t0t1 = gamma.approxArcLength(t0, t1, this->arclen_dt);
                    ratio       = arclen_t0t / arclen_t0t1; 

                    debugl(2, "LinearRadiusInterpolatorArcLen::operator(): t: %5.4f, arclen in [t0, t] = %10.5f, total arclen in domain [t0, t1]: %10.5f, ratio: %10.5f\n", t, arclen_t0t, arclen_t0t1, ratio);
                    if (ratio < 0) ratio = 0;
                    if (ratio > 1) ratio = 1;

                    return (this->r0 + (this->r1 - this->r0)*ratio);
                }
};

template<
    typename C2F,
    typename R
>
class LinearRadiusInterpolatorDomain {
    public:
                LinearRadiusInterpolatorDomain()
                {
                    this->r0    = 0;
                    this->r1    = 0;
                }

                /* NOTE: implicitly generated copy ctor and assignment operator suffice here */

        std::pair<R, R>
                getRadii() const
                {
                    return (std::pair<R, R>(this->r0, this->r1));
                }

        R       getMaxRadius() const
                {
                    return (std::max(this->r0, this->r1));
                }

        R       getMinRadius() const
                {
                    return (std::min(this->r0, this->r1));
                }

        void    setRadii(
                    R const                        &r0,
                    R const                        &r1)
                {
                    this->r0    = r0;
                    this->r1    = r1;
                }

        void    clipToInterval(
                    R const                        &t0,
                    R const                        &t1,
                    SpaceCurveReal<C2F, R> const   &gamma) const
                {
                    if (t0 <= t1) {
                        this->r0 = gamma.radiusEval(t0);
                        this->r1 = gamma.radiusEval(t1);
                    }
                    else {
                        throw("LinearRadiusInterpolatorDomain::clipToInterval(): malformed interval [t0, t1]: t0 > t1.");
                    }
                }

        R       operator()(
                    R const                        &t,
                    SpaceCurveReal<C2F, R> const   &gamma) const
                {

                    /* linearly interpolate with respect to domain, this amounts to non-linear radius interpolation of
                     * the parametric curve has non-constant parametric speed, which holds in the general case. */
                    R t0, t1, ratio;

                    auto domain = gamma.getDomain();
                    t0          = domain.first;
                    t1          = domain.second;

                    if (t < t0 || t > t1) {
                        throw("LinearRadiusInterpolatorDomain::operator(): given parameter value t not in domain [t0, t1].");
                    }

                    ratio = (t - t0) / (t1 - t0);
                    if (ratio < 0) ratio = 0;
                    if (ratio > 1) ratio = 1;

                    return (this->r0 + (this->r1 - this->r0)*ratio);
                }
};

/*! @brief Canal Surface Class
 *
 * detailed Canal Surface class description */
template <
    typename C2F,
    typename RadF,
    typename R
>
class CanalSurface {
    protected:
        /* spine curve and its first two derivatives */
        SpaceCurveReal<C2F, R>          spine_curve;
        /*
        SpaceCurveReal<C2F, R>          spine_curve_d;
        SpaceCurveReal<C2F, R>          spine_curve_d2;
        */

        /* functor for evaluating radii */
        RadF                            radius_functor;

        /* domain */
        R                               t0, t1;

        void                            updateSpineCurveDerivatives();
        
    public:
                                        CanalSurface();
                                        CanalSurface(
                                                SpaceCurveReal<C2F, R> const   &spine_curve,
                                                RadF const                     &radius_functor);

                                        CanalSurface(
                                                std::array<C2F, 3> const       &component_functors,
                                                RadF const                     &radius_functor,
                                                R const                        &t0,
                                                R const                        &t1);

        /* copy ctor and assignment operator */
                                        CanalSurface(const CanalSurface<C2F, RadF, R> &delta);
        CanalSurface<C2F, RadF, R>     &operator=(const CanalSurface<C2F, RadF, R> &delta);

        /* virtual dtor for polymorphism */
        virtual                        ~CanalSurface();

        /* virtual spine curve / radius evaluation */
        virtual Vec3<R>                 spineCurveEval(   R const &t) const;
        virtual Vec3<R>                 spineCurveEval_d( R const &t) const;
        virtual Vec3<R>                 spineCurveEval_d2(R const &t) const;
        virtual R                       radiusEval(R const &t) const;

        /* orthonormal bases ("frames") for the spine curve */
        void                            spineCurveGetRenderFrame(
                                            R const    &t,
                                            Vec3<R>     rvec,
                                            Vec3<R>    &x,
                                            Vec3<R>    &y,
                                            Vec3<R>    &z) const;

        void                            spineCurveGetFrenetFrame(
                                            R const    &t,
                                            Vec3<R>    &x,
                                            Vec3<R>    &y,
                                            Vec3<R>    &z) const;

        /* get approximate arc length of spine curve */
        R                               spineCurveApproxArcLength(
                                            R const &tstart, 
                                            R const &tend,
                                            R const &dt) const;

        /* domain */
        std::pair<R, R>                 getDomain() const;
        /*
        void                            setDomain(
                                        R const &t0, 
                                        R const &t1); 
        */

        /* wrapper for spine curve / radius evaluation */
        SpaceCurveReal<C2F, R> const   &getSpineCurve() const;
        void                            setSpineCurve(SpaceCurveReal<C2F, R> const &gamma);

        RadF                            getRadiusFunctor() const;
        void                            setRadiusFunctor(RadF const &radius_functor);

        /* mesh generation methods, templated with mesh template parameters  */
        template <typename Tm, typename Tv, typename Tf>
        void                            generateMesh(
                                            Mesh<Tm, Tv, Tf, R>                                    &M,
                                            uint32_t                                                nphisegments,
                                            Vec3<R>                                                 rvec,
                                            R const                                                &phi_0,
                                            R const                                                &arclen_dt               = 1E-3,
                                            bool                                                    start_circle_offset     = false,
                                            std::vector<
                                                    typename Mesh<Tm, Tv, Tf, R>::vertex_iterator
                                                >                                                  *start_circle_its        = NULL,
                                            bool                                                   *end_circle_offset       = NULL,
                                            std::vector<
                                                    typename Mesh<Tm, Tv, Tf, R>::vertex_iterator
                                                >                                                  *end_circle_its          = NULL,
                                            typename Mesh<Tm, Tv, Tf, R>::vertex_iterator          *closing_vertex_it       = NULL,
                                            bool                                                    preserve_crease_edges   = false,
											R														triangle_height_factor	= 1.0) const ;
};

#if 0
template <
    typename C2F,
    typename R,
>
class PipeSurface : public CanalSurface<C2F, ConstantRadiusFunctor, R>
{
    public:
                                        PipeSurface();
                                        PipeSurface();
                                        PipeSurface();

                                        PipeSurface(PipeSurface<C2F, R> const &delta);
        PipeSurface<C2F, R>            &operator=(PipeSurface<C2F, R> const &delta);

        virtual                        ~PipeSurface(); 
};
#endif

#include "Polynomial.hh"
#include "BivariatePolynomial.hh"

template <
    typename RadF,
    typename R = double
>
class BezierCanalSurface : public CanalSurface<BernsteinPolynomial<R, R>, RadF, R> {
    protected:
        /* override spine curve to be BezierCurve, which is more specific than 
         * SpaceCurveReal<BernsteinPolynomial<R, R>, R> */
        BezierCurve<R>              spine_curve;

        uint32_t                    degree;

        bool                        bb_set;
        BoundingBox<R>              bb;

    public:
                                    BezierCanalSurface();

                                    BezierCanalSurface(
                                            BezierCurve<R> const   &spine_curve,
                                            RadF const             &radius_functor);

                                    BezierCanalSurface(
                                            std::array<
                                                    BernsteinPolynomial<R, R>,
                                                    3
                                                > const                        &component_functors,
                                            RadF const                         &radius_functor);

                                    BezierCanalSurface(
                                            std::vector<Vec3<R>> const &control_points,
                                            RadF const                 &radius_functor);

                                    BezierCanalSurface(const BezierCanalSurface<RadF, R> &delta);
        BezierCanalSurface         &operator=(const BezierCanalSurface<RadF, R> &delta);

        virtual                    ~BezierCanalSurface();

        BezierCurve<R>              getSpineCurve() const;
        //void                        setSpineCurve(BezierCurve<R> const &gamma);

        /* info methods */
        uint32_t                    getDegree() const;

        /* clip to interval */
        void                        clipToInterval(
                                        R const    &t0,
                                        R const    &t1);

        /* bounding box */
        void                        updateBoundingBox(uint32_t spine_curve_subdivision_depth = 8);

        BoundingBox<R>              getBoundingBox() const;

        /* check render vector */
        R                           checkRenderVector(Vec3<R> const &r) const;
};

template<typename R>
class BLRCanalSurface : public BezierCanalSurface<LinearRadiusInterpolatorArcLen<BernsteinPolynomial<R, R>, R>, R>
{
    private:
        /* precomputed data for construction of global self-intersection polynomial */
        static bool                                     gsi_data_mutable;
        static uint32_t                                 gsi_max_n;
        static Tensor<2, BiBernsteinPolynomial<R, R>>   F;
        static Tensor<3, BiBernsteinPolynomial<R, R>>   G;

    public:
                                    BLRCanalSurface();
                                    BLRCanalSurface(
                                            BezierCurve<R> const   &spine_curve,
                                            R const                &r0,
                                            R const                &r1);  

                                    BLRCanalSurface(
                                            std::array<
                                                    BernsteinPolynomial<R, R>,
                                                    3
                                                > const                        &component_functors,
                                            R const                            &r0,
                                            R const                            &r1);


                                    BLRCanalSurface(
                                            std::vector<Vec3<R>> const &control_points,
                                            R const                    &r0,
                                            R const                    &r1);

                                    BLRCanalSurface(const BLRCanalSurface &delta);
        BLRCanalSurface<R>         &operator=(const BLRCanalSurface &delta);


        /*! \brief static method initializing precomputable global self-intersection data for
         * BLRCanalSurface objects. */
        static void                 initGlobalSelfIntersectionData(uint32_t max_n);
        static void                 setGlobalSelfIntersectionDataMutable();
        static void                 setGlobalSelfIntersectionDataImmutable();

                                   ~BLRCanalSurface();                                
        
        /* get radii, min and max radius */
        std::pair<R, R>             getRadii() const;
        R                           getMinRadius() const;
        R                           getMaxRadius() const;

        /* spine curve regularity polynomial */
        void                        spineCurveComputeRegularityPolynomial(BernsteinPolynomial<R, R> &p_reg) const;

        /* spine curve stationary point distance polynomial */
        void                        spineCurveComputeStationaryPointDistPoly(
                                        Vec3<R> const              &x,
                                        BernsteinPolynomial<R, R>  &p) const;

        /* spine curve local self-intersection polynomial */
        void                        computeLocalSelfIntersectionPolynomial(BernsteinPolynomial<R, R> &p_lsi) const;

        /* global self-intersection system */
        void                        computeGlobalSelfIntersectionSystem(
                                        BiBernsteinPolynomial<R, R>    &p,
                                        BiBernsteinPolynomial<R, R>    &q,
                                        BernsteinPolynomial<R, R>      &p_edge_t0,
                                        BernsteinPolynomial<R, R>      &p_edge_t1) const;

        /* intersection system with other BLRCanalSurface */
        void                        computeIntersectionSystem(
                                        const BLRCanalSurface<R>      &delta,
                                        BiBernsteinPolynomial<R, R>    &p,
                                        BiBernsteinPolynomial<R, R>    &q,
                                        BernsteinPolynomial<R, R>      &p_edge_x0,
                                        BernsteinPolynomial<R, R>      &p_edge_x1,
                                        BernsteinPolynomial<R, R>      &p_edge_y0,
                                        BernsteinPolynomial<R, R>      &p_edge_y1) const;

};

/* initialize static members */
template <typename R> bool                                      BLRCanalSurface<R>::gsi_data_mutable    = true;
template <typename R> uint32_t                                  BLRCanalSurface<R>::gsi_max_n           = 0;
template <typename R> Tensor<2, BiBernsteinPolynomial<R, R>>    BLRCanalSurface<R>::F({ 1, 1 });
template <typename R> Tensor<3, BiBernsteinPolynomial<R, R>>    BLRCanalSurface<R>::G({ 1, 1, 1 });

#include "../tsrc/CanalSurface.impl.hh"

#endif 
