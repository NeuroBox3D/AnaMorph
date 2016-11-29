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

/* ------------------------------------------------------------------------------------------------------------------ *
 *                                                                                                                    *
 *                        CanalSurface implementation..                                                               *
 *                                                                                                                    *
 * ------------------------------------------------------------------------------------------------------------------ */
#include "PolyAlgorithms.hh"

/* default constructor: implicitly default constructs spine curve and radius functor.  these are
 * required to behave in a consistent way and default-init to trivial zero-valued parametric curve /
 * radius functor */
template <typename C2F, typename RadF, typename R>
CanalSurface<C2F, RadF, R>::CanalSurface()
{
}

/* ctor with spine curve and radius functor. domain [t0, t1] is extracted from given spine curve */
template <typename C2F, typename RadF, typename R>
CanalSurface<C2F, RadF, R>::CanalSurface(
    SpaceCurveReal<C2F, R> const   &spine_curve,
    RadF const                     &radius_functor)
{
    this->spine_curve       = spine_curve;
    /*
    this->spine_curve_d     = spine_curve.getDerivative();
    this->spine_curve_d2    = this->spine_curve_d.getDerivative();
    */

    this->radius_functor    = radius_functor;

    /* obtain domain from given spine curve and set as domain of canal surface */
    auto domain             = spine_curve.getDomain();
    this->t0                = domain.first;
    this->t1                = domain.second;
}

/* ctor taking array of three component functors, the radius functor and the domain */
template <typename C2F, typename RadF, typename R>
CanalSurface<C2F, RadF, R>::CanalSurface(
    std::array<C2F, 3> const       &component_functors,
    RadF const                     &radius_functor,
    R const                        &t0,
    R const                        &t1)
{
    this->spine_curve       = SpaceCurveReal<C2F, R>(component_functors, t0, t1);
    /*
    this->spine_curve_d     = spine_curve.getDerivative();
    this->spine_curve_d2    = this->spine_curve_d.getDerivative();
    */

    this->radius_functor    = radius_functor;

    this->t0                = t0;
    this->t1                = t1;
}

/* copy ctor, assignment operator */
template <typename C2F, typename RadF, typename R>
CanalSurface<C2F, RadF, R>::CanalSurface(const CanalSurface<C2F, RadF, R> &delta) 
{
    this->spine_curve       = delta.spine_curve;
    /*
    this->spine_curve_d     = delta.spine_curve_d;
    this->spine_curve_d2    = delta.spine_curve_d2;
    */

    this->radius_functor    = delta.radius_functor;

    this->t0                = delta.t0;
    this->t1                = delta.t1;
}

/* copy ctor, assignment operator */
template <typename C2F, typename RadF, typename R>
CanalSurface<C2F, RadF, R> &
CanalSurface<C2F, RadF, R>::operator=(const CanalSurface<C2F, RadF, R> &delta) 
{
    this->spine_curve       = delta.spine_curve;
    /*
    this->spine_curve_d     = delta.spine_curve_d;
    this->spine_curve_d2    = delta.spine_curve_d2;
    */

    this->radius_functor    = delta.radius_functor;

    this->t0                = delta.t0;
    this->t1                = delta.t1;

    return (*this);
}

/* virtual dtor */
template <typename C2F, typename RadF, typename R>
CanalSurface<C2F, RadF, R>::~CanalSurface()
{
}

/* virtual spine curve / radius evaluation. delegators in effect.. */
template <typename C2F, typename RadF, typename R>
Vec3<R>
CanalSurface<C2F, RadF, R>::spineCurveEval(R const &t) const
{
    return (this->spine_curve.eval(t));
}

template <typename C2F, typename RadF, typename R>
Vec3<R>
CanalSurface<C2F, RadF, R>::spineCurveEval_d(R const &t) const
{
    return (this->spine_curve.eval_d(t));
}

template <typename C2F, typename RadF, typename R>
Vec3<R>
CanalSurface<C2F, RadF, R>::spineCurveEval_d2(R const &t) const
{
    return (this->spine_curve.eval_d2(t));
}

template <typename C2F, typename RadF, typename R>
R
CanalSurface<C2F, RadF, R>::radiusEval(R const &t) const
{
    return (this->radius_functor(t, this->spine_curve));
}

template <typename C2F, typename RadF, typename R>
void
CanalSurface<C2F, RadF, R>::spineCurveGetRenderFrame(
    R const    &t,
    Vec3<R>     rvec,
    Vec3<R>    &x,
    Vec3<R>    &y,
    Vec3<R>    &z) const
{
    this->spine_curve.getRenderFrame(t, rvec, x, y, z);
}

template <typename C2F, typename RadF, typename R>
void
CanalSurface<C2F, RadF, R>::spineCurveGetFrenetFrame(
    R const    &t,
    Vec3<R>    &x,
    Vec3<R>    &y,
    Vec3<R>    &z) const
{
    this->spine_curve.getFrenetFrame(t, x, y, z);
}

template <typename C2F, typename RadF, typename R>
R
CanalSurface<C2F, RadF, R>::spineCurveApproxArcLength(
    R const &tstart, 
    R const &tend,
    R const &dt) const
{
    return (this->spine_curve.approxArcLength(tstart, tend, dt));
}

template <typename C2F, typename RadF, typename R>
std::pair<R, R>
CanalSurface<C2F, RadF, R>::getDomain() const
{
    return (std::pair<R, R>(this->t0, this->t1));
}

/*
template <typename C2F, typename RadF, typename R>
void
CanalSurface<C2F, RadF, R>::setDomain(
    R const &t0, 
    R const &t1)
{
    this->t0    = t0;
    this->t1    = t1;
}
*/

template <typename C2F, typename RadF, typename R>
SpaceCurveReal<C2F, R> const &
CanalSurface<C2F, RadF, R>::getSpineCurve() const
{
    return (this->spine_curve);
}

template <typename C2F, typename RadF, typename R>
void
CanalSurface<C2F, RadF, R>::setSpineCurve(SpaceCurveReal<C2F, R> const &gamma)
{
    this->spine_curve = gamma;
}

template <typename C2F, typename RadF, typename R>
RadF
CanalSurface<C2F, RadF, R>::getRadiusFunctor() const
{
    return (this->radius_functor);
}

template <typename C2F, typename RadF, typename R>
void
CanalSurface<C2F, RadF, R>::setRadiusFunctor(RadF const &radius_functor)
{
    this->radius_functor = radius_functor;
}

template <typename C2F, typename RadF, typename R>
template <typename Tm, typename Tv, typename Tf>
void
CanalSurface<C2F, RadF, R>::generateMesh(
    Mesh<Tm, Tv, Tf, R>                                    &M,
    uint32_t                                                n_phi_segments,
    Vec3<R>                                                 rvec,
    R const                                                &phi_0,
    R const                                                &arclen_dt,
    bool                                                    start_circle_offset,
    std::vector<
            typename Mesh<Tm, Tv, Tf, R>::vertex_iterator
        >                                                  *start_circle_its,
    bool                                                   *end_circle_offset,
    std::vector<
            typename Mesh<Tm, Tv, Tf, R>::vertex_iterator
        >                                                  *end_circle_its,
    typename Mesh<Tm, Tv, Tf, R>::vertex_iterator          *closing_vertex_it,
	bool                                                    preserve_crease_edges) const
{
    debugl(0, "CanalSurface::generateMesh().\n");
    debugTabInc();

    uint32_t    i, j, ntsegments;
    R           t, r, phi, dphi, phi_offset;
    Vec3<R>     c, p, n, last_p, a, px, py, pz, vpos;

    dphi = Common::twopi / (R)n_phi_segments;

    /* compute t values: start at t = 0.0 and compute an increment based on the radius. if the
     * radius is larger, we want to take a larger step than if it is smaller. the step is computed
     * in the following way: with n_phi_segments given, the length of a line segment on the
     * approximated circle is 
     *
     * c = r*sqrt(2*(1-cos(dphi)))
     *
     * and since we want equilateral triangles and the next circle is "offset" with dphi / 2.0, we
     * want to choose the step in to such that the height of the triangles is approx
     *
     * h = c/2 *sqrt(3)
     *
     * for a unit-length cylinder, this will result in only equilateral triangles. in practive, this
     * strategy results in close to equilateral triangles also for canal surfaces */
    R               l, h, t_new, tlast_slack;

    /* get arc length of spine curve */
    R               arclength  = this->spineCurveApproxArcLength(this->t0, this->t1, arclen_dt);

    std::vector<R>  t_values = { this->t0 };
    t = this->t0;
    while (t < this->t1) {
        /* compute and push next t-value */
        r       = this->radiusEval(t);
        l       = r*sqrt(2.0*(1.0 - cos(dphi))); 
        if (!preserve_crease_edges)
        	h   = sqrt(3.0)*l / 2.0;
        else
        	h   = l;

        /* delta_t depends on the arc length: we want to advance a step of h in arc length, which
         * approximately amounts to an advance of h / l in t (for close to constant parametric
         * speed). */
        t_new   = t + h / arclength;
        t_values.push_back(t_new);

        /* and shift */
        t       = t_new;
    }

    /* the last value will be larger than 1.0, pop it, but first compute tlast_slack, which is the
     * last computed value - this->t1. the last two circles will have a distance that is too small =>
     * distribute tlast_slack uniformly among the other t values by setting t(i) = t(i) -
     * i*tlast_slack / ntsegments */
    tlast_slack = t_values.back() - this->t1;
    t_values.pop_back();

    /* the final value t1 is not yet present in the vector of t values right now, therefore
     * ntsegments is the size of the vector, instead of size() - 1 */
    ntsegments   = t_values.size();

    /* offset to make room for last to circles */
    for (i = 1; i < t_values.size(); i++) {
        t_values[i] -= (R)i*tlast_slack/ntsegments;
    }

    /* append t1 to finish the list */
    t_values.push_back(this->t1);

    /* vectors storing vertices of current and last circle */
    std::vector<typename Mesh<Tm, Tv, Tf, R>::vertex_iterator>  last_circle;
    std::vector<typename Mesh<Tm, Tv, Tf, R>::vertex_iterator>  current_circle;

    /* resize vector */
    current_circle.resize(n_phi_segments);
    last_circle.resize(n_phi_segments);

    /* if start_circle_its == NULL, generate initial circle and close it with triangles */
    if (start_circle_its == NULL) {
        /* get starting point and radius */
        p   = this->spineCurveEval(this->t0);
        r   = this->radiusEval(this->t0);

        /* get curve base for t0 */
        this->spineCurveGetRenderFrame(this->t0, rvec, px, py, pz);

        /* frenet base is available to demonstrate its disadvantageous nature for meshing, which is due
         * to strong spinning of the frenet trihedron, which is quite bad for discretization, since
         * quads get twisted then.. */
        //this->getFrenetFrame(this->t0, px, py, pz);

        debugl(0, "start_circle_its == NULL => generating initial circle.\n"); 

        Vec3<R>                                         start_closing_vertex_pos;
        typename Mesh<Tm, Tv, Tf, R>::vertex_iterator   start_closing_vertex_it;

        /* start circle is never offset, overwrite wrong argument if necessary */
        phi_offset          = phi_0;
        start_circle_offset = false;

        /* generate starting circle vertices and save iterators */
        for (j = 0; j < n_phi_segments; j++) {
            phi                         = ( (R)j * Common::twopi) / (R)n_phi_segments;
            vpos                        = p + py*(r*cos(phi + phi_offset)) + pz*(r*sin(phi + phi_offset));
            current_circle[j]           = M.vertices.insert(vpos);
            start_closing_vertex_pos    += vpos;
        }

        /* get centroid of initial circle */
        start_closing_vertex_pos   *= ( 1.0 / (R)n_phi_segments);
        start_closing_vertex_it     = M.vertices.insert(start_closing_vertex_pos);

        /* closing triangles */
        M.faces.insert(start_closing_vertex_it, current_circle[0], current_circle[n_phi_segments - 1]);
        for (j = 0; j < n_phi_segments - 1; j++) {
            M.faces.insert(start_closing_vertex_it, current_circle[j+1], current_circle[j]);
        }
    }
    /* otherwise, this canal surface is not the start canal segment and the previous canal surface
     * has already generated the vertices of the start circle, which is the end circle of the
     * previous segment. the caller needs to take care to reopen the mesh around the end circle of the previous segment
     * by deleting the closing vertex. */
    else {
        /* NOTE: n_phi_segments has to match along the whole path for this to work properly.. */
        if (start_circle_its->size() == n_phi_segments) {
            current_circle = *start_circle_its;
            debugl(0, "start_circle_its != NULL => deleting closing vertex and taking given initial circle vertex ids..\n"); 
        }
        else {
            throw("CanalSurface::generateMesh(): supplied start_circle_its does not consist of exactly n_phi_segments vertices => invalid arguments supplied by the caller.");
        }
    }

    /* now the last two circles will, in general, have a distance that is too small compared to the
     * radii => compute the difference to the approximate good value and distribute this difference
     * equally among all the other segments */
    debugl(0, "rendering inner circles..\n");
    for (i = 1; i < ntsegments; i++) {
        /* if start is offset, then even values for i are not offset, odd values are offset.
         * otherwise, vice versa */
        if (!preserve_crease_edges)
        {
			if (start_circle_offset) {
				if (i % 2 == 0) {
					phi_offset = phi_0 + dphi / 2.0;
				}
				else {
					phi_offset = phi_0;
				}
			}
			else {
				if (i % 2 == 0) {
					phi_offset = phi_0;
				}
				else {
					phi_offset = phi_0 + dphi / 2.0;
				}
			}
        }

        //t = (R)i / (R)ntsegments;
        t = t_values[i];

        /* shift current_circle to last_circle. same for p and last_p */
        last_circle = current_circle;
        last_p      = p;

        /* calculate new current point and radius */
        p           = this->spineCurveEval(t);
        r           = this->radiusEval(t);

        this->spineCurveGetRenderFrame(t, rvec, px, py, pz);
        //this->spineCurveGetFrenetFrame(t, px, py, pz);

        /* generate current circle vertices and save iterators */
        for (j = 0; j < n_phi_segments; j++) {
            phi                 = ( (R)j * Common::twopi) / (R)n_phi_segments;
            current_circle[j]   = M.vertices.insert( p + py*(r*cos(phi + phi_offset)) + pz*(r*sin(phi + phi_offset)) );
        }

        /* generate quad faces between last_circle and current_circle */
        for (j = 0; j < n_phi_segments - 1; j++) {
            M.faces.insert(last_circle[j], last_circle[j+1], current_circle[j+1], current_circle[j]);
        }

        /* closing quad at index warp around */
        M.faces.insert(last_circle[n_phi_segments - 1], last_circle[0], current_circle[0], current_circle[n_phi_segments - 1]);
    }

    debugl(0, "rendering last circle..\n");
    /* and the same procedure for the endpoint, i.e. t = t1 */
    last_circle = current_circle;
    last_p      = p;
    p           = this->spineCurveEval(this->t1);
    r           = this->radiusEval(this->t1);

    this->spineCurveGetRenderFrame(this->t1, rvec, px, py, pz);
    //this->spineCurveGetFrenetFrame(this->t1, px, py, pz);

    Vec3<R>                                         end_closing_vertex_pos = Aux::VecMat::nullvec<R>();
    typename Mesh<Tm, Tv, Tf, R>::vertex_iterator   end_closing_vertex_it;
    
    /* last circle has index ntsegments, since we got (ntsegments + 1) circles */
    if (!preserve_crease_edges)
    {
    	if (start_circle_offset) {
			if (ntsegments % 2 == 0) {
				phi_offset = phi_0 + dphi / 2.0;
			}
			else {
				phi_offset = phi_0;
			}
		}
		else {
			if (ntsegments % 2 == 0) {
				phi_offset = phi_0;
			}
			else {
				phi_offset = phi_0 + dphi / 2.0;
			}
		}
    }

    px.print_debugl(0);
    py.print_debugl(0);
    pz.print_debugl(0);
    debugTabInc();
    for (j = 0; j < n_phi_segments; j++) {
        phi                     = ( (R)j * Common::twopi) / (R)n_phi_segments;
        vpos                    = p + py*(r*cos(phi + phi_offset)) + pz*(r*sin(phi + phi_offset));
        current_circle[j]       = M.vertices.insert(vpos);
        end_closing_vertex_pos += vpos;
    }
    debugTabDec();

    /* get centroid of last circle vertices and place closing vertex in the middle */
    end_closing_vertex_pos *= (1.0 / (R)n_phi_segments);
    end_closing_vertex_it   = M.vertices.insert(end_closing_vertex_pos);

    /* generate quad faces between last_circle and current_circle */
    for (j = 0; j < n_phi_segments - 1; j++) {
        M.faces.insert(last_circle[j], last_circle[j+1], current_circle[j+1], current_circle[j]);
    }

    /* closing quad at index warp around */
    M.faces.insert(last_circle[n_phi_segments - 1], last_circle[0], current_circle[0], current_circle[n_phi_segments - 1]);

    /* add closing triangles at the end */
    M.faces.insert(end_closing_vertex_it, current_circle[n_phi_segments - 1], current_circle[0]);
    for (j = 0; j < n_phi_segments - 1; j++) {
        M.faces.insert(end_closing_vertex_it, current_circle[j], current_circle[j + 1] );
    }

    /* if desired by the caller, write out information to referenced data */
    if (end_circle_its) {
        *end_circle_its = current_circle;
    }

    if (closing_vertex_it) {
        *closing_vertex_it = end_closing_vertex_it;
    }

    /* we got ntsegments + 1 circles with indices 0..ntsegments in total. so if ncircles is even,
     * then the offset bit will swap, otherwise it remains the same => XOR with ncircles_even */
    if (end_circle_offset) {
        *end_circle_offset   = Aux::Logic::lxor(start_circle_offset, (ntsegments + 1) % 2 == 0);
        //printf("start_circle_offset = %d, ntsegments: %5d, ncircles_even: %5d, neven: %d => end_circle_offset = %d\n", start_circle_offset, ntsegments, ncircles, ncircles_even, end_circle_offset);
    }

    debugTabDec();
    debugl(0, "CanalSurface::generateMesh(): done.\n");
}

/* ------------------------------------------------------------------------------------------------------------------ *
 *                                                                                                                    *
 *                        BezierCanalSurface implementation..                                                         *
 *                                                                                                                    *
 * ------------------------------------------------------------------------------------------------------------------ */
template <typename RadF, typename R>
BezierCanalSurface<RadF, R>::BezierCanalSurface()
    : CanalSurface<BernsteinPolynomial<R, R>, RadF, R>()
{
    /* domain is always [0,1] for BezierCurves and hence also for BezierCanalSurface */
    this->t0 = 0;
    this->t1 = 1;
}

template <typename RadF, typename R>
BezierCanalSurface<RadF, R>::BezierCanalSurface(
    BezierCurve<R> const   &spine_curve,
    RadF const             &radius_functor)
        : CanalSurface<BernsteinPolynomial<R, R>, RadF, R>(spine_curve, radius_functor)
{
    /* CanalSurface constructor sets [t0, t1] = [0,1] through domain of BezierCurve parameter spine_curve */
    this->spine_curve = spine_curve;
}

template <typename RadF, typename R>
BezierCanalSurface<RadF, R>::BezierCanalSurface(
    std::array<
            BernsteinPolynomial<R, R>,
            3
        > const                        &component_functors,
    RadF const                         &radius_functor)
    /* in-line constructed BezierCurve from component_functors arrays is passed to base class CanalSurface ctor. */
        : CanalSurface<BernsteinPolynomial<R, R>, RadF, R>(
                BezierCurve<R>(component_functors),
                radius_functor
            )
{
    this->spine_curve = BezierCurve<R>(component_functors);
}

/* same as previous constructor above for control-point based construction of BezierCurves */
template <typename RadF, typename R>
BezierCanalSurface<RadF, R>::BezierCanalSurface(
    std::vector<Vec3<R>> const &control_points,
    RadF const                 &radius_functor)
    /* in-line constructed BezierCurve from component_functors arrays is passed to base class CanalSurface ctor. */
        : CanalSurface<BernsteinPolynomial<R, R>, RadF, R>(
                BezierCurve<R>(control_points),
                radius_functor
            )
{
    this->spine_curve = BezierCurve<R>(control_points);
}

template <typename RadF, typename R>
BezierCanalSurface<RadF, R>::BezierCanalSurface(BezierCanalSurface const &delta)
    : CanalSurface<BernsteinPolynomial<R, R>, RadF, R>(delta)
{
    this->spine_curve   = delta.spine_curve;
    this->degree        = delta.degree;
    this->bb_set        = delta.bb_set;
    this->bb            = delta.bb;
}

template <typename RadF, typename R>
BezierCanalSurface<RadF, R> &
BezierCanalSurface<RadF, R>::operator=(BezierCanalSurface const &delta)
{
    CanalSurface<BernsteinPolynomial<R, R>, RadF, R>::operator=(delta);
    this->spine_curve   = delta.spine_curve;
    this->degree        = delta.degree;
    this->bb_set        = delta.bb_set;
    this->bb            = delta.bb;

    return (*this);
}

template <typename RadF, typename R>
BezierCanalSurface<RadF, R>::~BezierCanalSurface()
{
}

/* since spine_curve has been overridden with more specific type BezierCurve in BezierCanalSurface, also override
 * getSpineCurve() */
template <typename RadF, typename R>
BezierCurve<R>
BezierCanalSurface<RadF, R>::getSpineCurve() const
{
    return (this->spine_curve);
}

template <typename RadF, typename R>
uint32_t
BezierCanalSurface<RadF, R>::getDegree() const
{
    return (this->spine_curve.getDegree());
}

template <typename RadF, typename R>
void
BezierCanalSurface<RadF, R>::clipToInterval(
    R const    &t0,
    R const    &t1)
{
    if (t0 <= t1) {
        /* clip radius functor BEFORE clipping spine curve, since the old domain is needed to get radii */
        this->radius_functor.clipToInterval(t0, t1, this->spine_curve);

        /* clip this->spine_curve, which is a BezierCurve<R>. the spine_curve members from the base class
         * CanalSurface< BernsteinPolynomial<R, R>, RadF, R> has been overridden in BezierCanalSurface, which
         * necessitates the update below. */
        this->spine_curve.clipToInterval(t0, t1, &(this->spine_curve));

        /* update spine curve object in base class. note that BezierCanalSurface::spine_curve is a BezierCurve, which is
         * a derived class of SpaceCurveReal< BernsteinPolynomial<R, R>, R>, which is in turn the spine curve type of
         * the base class CanalSurface< BernsteinPolynomial<R, R>, RadF, R>.
         *
         * after clipping the BezierCurve in (this) BezierCanalSurface, use the resulting clipped BezierCurve to update
         * the base class member CanalSurface< BernsteinPolynomial<R, R>, RadF, R>::spine_curve. the following update
         * assignment "slices" the assigned BezierCurve, since the specifics of BezierCurve<R> are lost when assigning
         * to an object of the base class type SpaceCurveReal< BernsteinPolynomial<R, R>, R>.  this slicing is
         * unproblematic here however, since the CanalSurface<..> base class of BezierCanalSurface only require the
         * functionality at the "sliced" base level of abstraction. */
        CanalSurface< BernsteinPolynomial<R, R>, RadF, R>::spine_curve = this->spine_curve;
    }
    else {
        throw("BezierCanalSurface::clipToInterval(): malformed interval [t0, t1]: t0 > t1.");
    }
}

template <typename RadF, typename R>
void
BezierCanalSurface<RadF, R>::updateBoundingBox(uint32_t spine_curve_subdivision_depth)
{
    using Aux::VecMat::onesVec3;
    using Aux::VecMat::fabsVec3;
    using Aux::VecMat::maxVec3;

    this->bb        = this->spine_curve.getBoundingBox(spine_curve_subdivision_depth);
    
    /* extend bounding box of spine curve with maximum radius */
    R rmax          = this->radius_functor.getMaxRadius();
    Vec3<R> offset  = onesVec3<R>() * rmax;
    this->bb.extend(0, offset);

    /* extend resulting bounding box by 2,5%, but no less than 1E-3, in every component. */
    this->bb.extend(0.025, Vec3<R>(1E-3, 1E-3, 1E-3));
    this->bb_set    = true;
}

template <typename RadF, typename R>
BoundingBox<R>
BezierCanalSurface<RadF, R>::getBoundingBox() const
{
    if (this->bb_set) {
        return (this->bb);
    }
    else {
        throw("BezierCanalSurface::getBoundingBox(): bounding box has not been updated. use updateBoundingBox() first. this is a const method.");
    }
}

template <typename RadF, typename R>
R
BezierCanalSurface<RadF, R>::checkRenderVector(Vec3<R> const &r) const
{
    debugl(2, "BezierCanalSurface::checkRenderVector():\n");
    debugTabInc();

    BezierCurve<R> const       &gamma   = this->spine_curve;
    BezierCurve<R>              dgamma  = gamma.getDerivative();

    BernsteinPolynomial<R, R>   p, r_cross_dgamma[3], q, z;
   
    r_cross_dgamma[0]   = dgamma[2]*r[1] - dgamma[1]*r[2];
    r_cross_dgamma[1]   = dgamma[0]*r[2] - dgamma[2]*r[0];
    r_cross_dgamma[2]   = dgamma[1]*r[0] - dgamma[0]*r[1];
    p                   = r_cross_dgamma[0].square() + r_cross_dgamma[1].square() + r_cross_dgamma[2].square();
    q                   = dgamma[0].square() + dgamma[1].square() + dgamma[2].square();

    /* target polynomial, numerator of rational function describing || r \times t(t) || */
    z                   = (p.getDerivative()).multiply(q) - p.multiply(q.getDerivative());

    /* locate roots with bezier clipping */
    std::vector<PolyAlg::RealInterval<R>> z_roots;
    PolyAlg::BezClip_roots<R>(z, 0.0, 1.0, 1E-6, z_roots);

    R f_min = Aux::Numbers::inf<R>();

    /* boundary values */
    f_min = std::min(f_min, p.eval(0.0) / q.eval(0.0));
    f_min = std::min(f_min, p.eval(1.0) / q.eval(1.0));


    /* check all all roots in zroots */
    R root;
    debugTabInc();
    for (auto root_interval : z_roots) {
        root    = (root_interval.t0 + root_interval.t1) / 2.0;
        f_min   = std::min(f_min, p.eval(root) / q.eval(root));
        debugl(2, "function value of rational target function for root %5.4f of numerator polynomial: %5.4f\n",
            root, p.eval(root) / q.eval(root));
    }
    debugTabDec();
    debugl(2, "minimum target function value fmin: %5.4f\n", fmin);

    debugTabDec();
    debugl(2, "BezierCanalSurface::checkRenderVector(): done. f_min: %5.4f\n", f_min);

    return f_min;
}

/* ------------------------------------------------------------------------------------------------------------------ *
 *                                                                                                                    *
 *                        BLRCanalSurface implementation.                                                             *
 *                                                                                                                    *
 * ------------------------------------------------------------------------------------------------------------------ */
template <typename R>
BLRCanalSurface<R>::BLRCanalSurface()
    : BezierCanalSurface<LinearRadiusInterpolatorArcLen<BernsteinPolynomial<R, R>, R>, R>()
{
    /* domain is always [0,1], initialized by BezierCanalSurface ctor */
}

template <typename R>
BLRCanalSurface<R>::BLRCanalSurface(
    BezierCurve<R> const   &spine_curve,
    R const                &r0,
    R const                &r1)
        : BezierCanalSurface<LinearRadiusInterpolatorArcLen<BernsteinPolynomial<R, R>, R>, R>(
            spine_curve,
            LinearRadiusInterpolatorArcLen<BernsteinPolynomial<R,R>, R>(r0, r1)
        )
{
}

template <typename R>
BLRCanalSurface<R>::BLRCanalSurface(
    std::array<
            BernsteinPolynomial<R, R>,
            3
        > const                        &component_functors,
    R const                            &r0,
    R const                            &r1)
        : BezierCanalSurface<LinearRadiusInterpolatorArcLen<BernsteinPolynomial<R, R>, R>, R>(
            component_functors,
            LinearRadiusInterpolatorArcLen<BernsteinPolynomial<R,R>, R>(r0, r1)
        )
{
}

/* same as previous constructor above for control-point based construction of BezierCurves */
template <typename R>
BLRCanalSurface<R>::BLRCanalSurface(
    std::vector<Vec3<R>> const &control_points,
    R const                    &r0,
    R const                    &r1)
        : BezierCanalSurface<LinearRadiusInterpolatorArcLen<BernsteinPolynomial<R, R>, R>, R>(
            control_points,
            LinearRadiusInterpolatorArcLen<BernsteinPolynomial<R,R>, R>(r0, r1)
        )
{
}

template <typename R>
BLRCanalSurface<R>::BLRCanalSurface(BLRCanalSurface const &delta)
    : BezierCanalSurface<LinearRadiusInterpolatorArcLen<BernsteinPolynomial<R, R>, R>, R>(delta)
{
}

template <typename R>
BLRCanalSurface<R> &
BLRCanalSurface<R>::operator=(BLRCanalSurface const &delta)
{
    BezierCanalSurface<
            LinearRadiusInterpolatorArcLen<BernsteinPolynomial<R, R>, R>,
            R
        >::operator=(delta);

    return (*this);
}

template <typename R>
void
BLRCanalSurface<R>::initGlobalSelfIntersectionData(uint32_t gsi_max_n)
{
    if (!BLRCanalSurface<R>::gsi_data_mutable) {
        throw("(static) BLRCanalSurface<R>::initGlobalSelfIntersectionData(): can't process, since static global self-intersection data is set to immutable.");
    }

    uint32_t i, j, k, m, n;
    std::vector<PowerPolynomial<R, R>>  B_n_pow;
    Matrix<R>                           F_i_n_powercoeff;

    /* resize tensors F and G to the required sizes. */

    /*
    F   = alloc2dArray<BernsteinBasisBiPoly>(CanalSurface::gsi_max_n + 1, CanalSurface::gsi_max_n + 1);
    G   = alloc3dArray<BernsteinBasisBiPoly>(CanalSurface::gsi_max_n + 1, CanalSurface::gsi_max_n + 1, CanalSurface::gsi_max_n + 1);
    */

    BLRCanalSurface<R>::F.resize({ gsi_max_n + 1, gsi_max_n + 1});
    BLRCanalSurface<R>::G.resize({ gsi_max_n + 1, gsi_max_n + 1, gsi_max_n + 1 });

    /* up to the defined max value of n */
    for (n = 1; n < gsi_max_n + 1; n++) {

        /* compute the M(n) representation of B_i^n(t) for i = 0..n */
        B_n_pow.resize(n+1);

        for (i = 0; i < n + 1; i++) {
            //B_n_pow[i].convertFromBernsteinBasis( Aux::VecMat::kronecker_vec(n+1, i) );

            PolyAlg::convertBasis<R>(
                    B_n_pow[i],
                    PolyAlg::computeBernsteinBasisPoly<R, R>(n, i)
                );
        }

        debugl(2, "\t B_i^n in M(n) computed, 0 <= i <= n\n");

        /* compute all bivariate polynomials F_i^n(x,y), i = 0..n */
        for (i = 0; i < n + 1; i++) {
            /* first, compute coefficients of F_i_n in M(n-1, n-1) => nxn coefficient matrix */
            debugl(3, "\n\t\t i = %d, computing M(n-1, n-1) coefficients of F_i_n ..\n", i);

            /* init F_i_n_powercoeff to zerod (n, n) matrix */
            F_i_n_powercoeff.resize(n, n);
            F_i_n_powercoeff.fill(0.0);

            for (k = 1; k < n + 1; k++) {
                for (m = 0; m < k; m++) {
                    debugl(3, "\t\t\t F_i_n_powercoeff(%d, %d) += %f\n", m, k-1-m, B_n_pow[i][k]);

                    F_i_n_powercoeff(m, k - 1 - m) += B_n_pow[i][k];
                }
            }
            debugl(3, "\t\t done.\n");
            debugl(2, "\t\t i = %d, computing B(n-1, n-1) form of F_i_n through conversion .. ", i);

            /* now generate F_i^n from the power coefficient matrix F_i_n_powercoeff. */
            F({ n, i }).convertFromPowerBasis(F_i_n_powercoeff);

            /* ------ */
            uint32_t fm, fn;
            F({ n, i }).getDegree(fm, fn);
            debugl(2, "\t\tF[n][i = %d] computed in BB(%d, %d) ", i, fm, fn);
            /* ------ */
        }

        uint32_t tmp_m, tmp_n;
        
        /* compute all bivariate polynomials G_{ij}^n, i, j = 0..n */
        BiBernsteinPolynomial<R, R> B_in_y, B_jn_y;

        debugl(2, "\t------- Computing Gij..\n");
        for (i = 0; i < n + 1; i++) {
            for (j = 0; j < n + 1; j++) {
                debugl(3, "\t\t i = %d, j = %d\n", i, j);
                /*
                B_in_y = ( BB(n, i).convertToBiPoly(n, false) );
                B_jn_y = ( BB(n, j).convertToBiPoly(n, false) );
                */

                B_in_y = PolyAlg::BernsteinConvertToBiPoly(
                       PolyAlg::computeBernsteinBasisPoly<R, R>(n, i),
                       n,
                       false
                    );

                B_jn_y = PolyAlg::BernsteinConvertToBiPoly(
                        PolyAlg::computeBernsteinBasisPoly<R, R>(n, j),
                        n,
                        false
                    );

                /* ----- */
                B_in_y.getDegree(tmp_m, tmp_n);
                debugl(3, "\t\t B_i^n(y) converted to bipoly: BB(%d, %d)\n", tmp_m, tmp_n);
                B_jn_y.getDegree(tmp_m, tmp_n);
                debugl(3, "\t\t B_j^n(y) converted to bipoly: BB(%d, %d)\n", tmp_m, tmp_n);
                /* ----- */

                /* compute coefficient G(n, i, j) */
                //G[n][i][j] = B_jn_y.multiply( F[n][i] ) - B_in_y.multiply( F[n][j] );
                G({ n, i, j }) = B_jn_y.multiply(F({ n, i })) - B_in_y.multiply(F({ n, j }));

                /* ----- */
                G({ n, i, j }).getDegree(tmp_m, tmp_n);
                debugl(2, "\t G[n][i = %d][j = %d] computed in BB(%d, %d)\n", i, j, tmp_m, tmp_n);
                /* ----- */
            }
        }
    }

    /* update static members */
    BLRCanalSurface<R>::gsi_max_n = gsi_max_n;
}

template <typename R>
void
BLRCanalSurface<R>::setGlobalSelfIntersectionDataMutable()
{
    BLRCanalSurface<R>::gsi_data_mutable = true;
}

template <typename R>
void
BLRCanalSurface<R>::setGlobalSelfIntersectionDataImmutable()
{
    BLRCanalSurface<R>::gsi_data_mutable = false;
}

template <typename R>
BLRCanalSurface<R>::~BLRCanalSurface()
{
}

template <typename R>
std::pair<R, R>
BLRCanalSurface<R>::getRadii() const
{
    return (this->radius_functor.getRadii());
}

template <typename R>
R
BLRCanalSurface<R>::getMinRadius() const
{
    auto rpair = this->radius_functor.getRadii();
    return (std::min(rpair.first, rpair.second));
}

template <typename R>
R
BLRCanalSurface<R>::getMaxRadius() const
{
    auto rpair = this->radius_functor.getRadii();
    return (std::max(rpair.first, rpair.second));
}

template <typename R>
void
BLRCanalSurface<R>::spineCurveComputeRegularityPolynomial(BernsteinPolynomial<R, R> &p_reg) const
{
    this->spine_curve.computeRegularityPolynomial(p_reg);
}

template <typename R>
void
BLRCanalSurface<R>::spineCurveComputeStationaryPointDistPoly(
    Vec3<R> const              &x,
    BernsteinPolynomial<R, R>  &p) const
{
    this->spine_curve.computeStationaryPointDistPoly(x, p);
}

template <typename R>
void
BLRCanalSurface<R>::computeLocalSelfIntersectionPolynomial(BernsteinPolynomial<R, R> &p_lsi) const
{
    debugl(2, "BLRCanalSurface::computeLocalSelfIntersectionPolynomial()\n");
    debugTabInc();

    /* get degree, maximum radius, reference gamma to (this->spine_curve), first two derivatives dgamma and d2gamma of
     * spine_curve. */
    uint32_t const              n       = this->getDegree();
    R const                     rmax    = this->getMaxRadius();
    BezierCurve<R> const       &gamma   = this->spine_curve;
    BezierCurve<R> const        dgamma  = gamma.getDerivative();
    BezierCurve<R> const        d2gamma = dgamma.getDerivative();

    /* three components of the cross product (dgamma \cdot d2gamma). every summand is distinct
     * and has to be computed exactly once. */
    BernsteinPolynomial<R, R>   d_d2_yz_minus_zy, d_d2_zx_minus_xz,d_d2_xy_minus_yx;
    BernsteinPolynomial<R, R>   crossprod_square;
    BernsteinPolynomial<R, R>   dgamma_square;
    BernsteinPolynomial<R, R>   dgamma_sqcube;

    /* the derivative dgamma is represented in BB(n-1), d2gamma in BB(n-2), no degree elevantion was performed!
     * multiplying works nonetheless and yields a the summands of the cross product in in BB(2n-3) */
    d_d2_yz_minus_zy    = dgamma[1].multiply(d2gamma[2]) - dgamma[2].multiply(d2gamma[1]);
    d_d2_zx_minus_xz    = dgamma[2].multiply(d2gamma[0]) - dgamma[0].multiply(d2gamma[2]);
    d_d2_xy_minus_yx    = dgamma[0].multiply(d2gamma[1]) - dgamma[1].multiply(d2gamma[0]);

    /* the square of the cross product (that is, the square of the norm, inner product with itself) can be computed by
     * squaring the above terms and adding up. this yields a polynomial in BB(4n-6), which is then directly multiplied
     * by rmax^2 */
    crossprod_square    = (
            (d_d2_yz_minus_zy.multiply(d_d2_yz_minus_zy)) +
            (d_d2_zx_minus_xz.multiply(d_d2_zx_minus_xz)) + 
            (d_d2_xy_minus_yx.multiply(d_d2_xy_minus_yx))
        );

    crossprod_square    = crossprod_square * (rmax * rmax);

    /* the other term we need is (dgamma_x^2 + dgamma_y^2 + dgamma_z^2)^3, i.e.  the third power of the inner product of
     * the tangent vector with itself.  the tangent innner product yields coefficients in BB(2n-2), third power yields
     * BB(6n-6) => the other term, crossprod_square, is in BB(4n-6) and must therefore be degree-elevated by 2n. */
    dgamma_square       = 
            (dgamma[0].multiply(dgamma[0])) +
            (dgamma[1].multiply(dgamma[1])) +
            (dgamma[2].multiply(dgamma[2]));

    /* third power done via three multiplications, it should be possible in one step, but
     * we need multinomial coefficient then I think...*/
    dgamma_sqcube       = (dgamma_square.multiply(dgamma_square)).multiply(dgamma_square);

    /* elevante crossprod_square by 2n */
    crossprod_square.elevateDegree(2*n);

    /* self-intersection polynomial is now simply the difference. to check if its negative over its
     * entire domain [0,1], check corners values and compute roots */
    p_lsi               = crossprod_square - dgamma_sqcube;

    debugTabDec();
    debugl(2, "BLRCanalSurface::computeLocalSelfIntersectionPolynomial(): done.\n");
}

template <typename R>
void
BLRCanalSurface<R>::computeGlobalSelfIntersectionSystem(
    BiBernsteinPolynomial<R, R>    &p,
    BiBernsteinPolynomial<R, R>    &q,
    BernsteinPolynomial<R, R>      &p_edge_t0,
    BernsteinPolynomial<R, R>      &p_edge_t1) const
{
    debugl(1, "BLRCanalSurface::computeIntersectionSystem().\n");
    debugTabInc();

    /* for consistency with the above and the thesis, use const reference Gamma for "this" again */
    uint32_t const              n       = this->getDegree();
    BLRCanalSurface<R> const   &Gamma   = (*this);
    BezierCurve<R> const       &gamma   = Gamma.spine_curve;
    BezierCurve<R> const        dgamma  = gamma.getDerivative();

    /* distance vector with trivial solution factored out */
    std::vector<BiBernsteinPolynomial<R, R>> dist_nt(3, BiBernsteinPolynomial<R, R>(2*n-1, 2*n-1, 0));
    std::vector<BiBernsteinPolynomial<R, R>> dist_trivial(3);

    /*
    dist_nt[0].initConstant(2*n-1, 2*n-1, 0.0);
    dist_nt[1].initConstant(2*n-1, 2*n-1, 0.0);
    dist_nt[2].initConstant(2*n-1, 2*n-1, 0.0);
    */

    dist_trivial[0] = PolyAlg::BernsteinConvertToBiPoly(gamma[0], n, true) - 
                      PolyAlg::BernsteinConvertToBiPoly(gamma[0], n, false);

    dist_trivial[1] = PolyAlg::BernsteinConvertToBiPoly(gamma[1], n, true) - 
                      PolyAlg::BernsteinConvertToBiPoly(gamma[1], n, false);

    dist_trivial[2] = PolyAlg::BernsteinConvertToBiPoly(gamma[2], n, true) - 
                      PolyAlg::BernsteinConvertToBiPoly(gamma[2], n, false);

    /* dist_nt can be computed with the precomputed bivariate polynomials G_{ij}^n */
    uint32_t i, j;
    for (i = 0; i < n+1; i++) {
        for (j = 0; j < n+1; j++) {
            dist_nt[0] += (G({n, i, j}) * gamma[0][i]);
            dist_nt[1] += (G({n, i, j}) * gamma[1][i]);
            dist_nt[2] += (G({n, i, j}) * gamma[2][i]);
        }
    }

    /*
    double  x, y;
    Vec3    d_xy;
    Vec3    test;
    FILE    *df0, *df1, *df2;
    df0 = fopen("dist_div_x_y_0.plot", "w");
    df1 = fopen("dist_div_x_y_1.plot", "w");
    df2 = fopen("dist_div_x_y_2.plot", "w");

    uint32_t ticks = 100;
    for (i = 0; i <= ticks; i++) {
        x = (double)i / (double)ticks;
        for (j = 0; j <= ticks; j++) {
            y = (double)j / (double)ticks;
               
            d_xy = this->eval(x) - this->eval(y);
            
            if (x == y) {
                fprintf(df0, "%12.5E %12.5E %12.5E\n", x, y, 0.0);
                fprintf(df1, "%12.5E %12.5E %12.5E\n", x, y, 0.0);
                fprintf(df2, "%12.5E %12.5E %12.5E\n", x, y, 0.0);
            }
            else {
                fprintf(df0, "%12.5E %12.5E %12.5E\n", x, y, d_xy[0] / (x-y) - dist_nt[0].eval(x,y));
                fprintf(df1, "%12.5E %12.5E %12.5E\n", x, y, d_xy[1] / (x-y) - dist_nt[1].eval(x,y));
                fprintf(df2, "%12.5E %12.5E %12.5E\n", x, y, d_xy[2] / (x-y) - dist_nt[2].eval(x,y));
            }
        }
        fprintf(df0, "\n");
        fprintf(df1, "\n");
        fprintf(df2, "\n");
    }

    fclose(df0);
    fclose(df1);
    fclose(df2);

    dist_nt[0].writePlotFile(ticks, "dist_nt_0.plot");
    dist_nt[1].writePlotFile(ticks, "dist_nt_1.plot");
    dist_nt[2].writePlotFile(ticks, "dist_nt_2.plot");
    */

    /* p and q are computed as for the intersection of two pipe surfaces, only the specially
     * prepared "non trivial" distance vector dist_nt is used. */
    p = dist_nt[0].multiply( PolyAlg::BernsteinConvertToBiPoly(dgamma[0], 0, true) ) +
        dist_nt[1].multiply( PolyAlg::BernsteinConvertToBiPoly(dgamma[1], 0, true) ) +
        dist_nt[2].multiply( PolyAlg::BernsteinConvertToBiPoly(dgamma[2], 0, true) );

    q = dist_nt[0].multiply( PolyAlg::BernsteinConvertToBiPoly(dgamma[0], 0, false) ) +
        dist_nt[1].multiply( PolyAlg::BernsteinConvertToBiPoly(dgamma[1], 0, false) ) +
        dist_nt[2].multiply( PolyAlg::BernsteinConvertToBiPoly(dgamma[2], 0, false) );

    /*
    BiBernsteinPolynomial<R, R> p_trivial, q_trivial;

    p_trivial = 
        dist_trivial[0].multiply( PolyAlg::BernsteinConvertToBiPoly(dgamma[0], 0, true) ) +
        dist_trivial[1].multiply( PolyAlg::BernsteinConvertToBiPoly(dgamma[1], 0, true) ) +
        dist_trivial[2].multiply( PolyAlg::BernsteinConvertToBiPoly(dgamma[2], 0, true) );

    q_trivial =
        dist_trivial[0].multiply( PolyAlg::BernsteinConvertToBiPoly(dgamma[0], 0, false) ) +
        dist_trivial[1].multiply( PolyAlg::BernsteinConvertToBiPoly(dgamma[1], 0, false) ) +
        dist_trivial[2].multiply( PolyAlg::BernsteinConvertToBiPoly(dgamma[2], 0, false) );

    BiBernsteinPolynomial<R, R> ptmp = p, qtmp = q, r;
    ptmp.elevateDegree(0, n-1);
    qtmp.elevateDegree(n-1, 0);
    r = ptmp - qtmp;

    p.writePlotFile(200, "gsi_p.plot");
    q.writePlotFile(200, "gsi_q.plot");
    r.writePlotFile(200, "gsi_p_minus_q.plot");

    p_trivial.writePlotFile(200, "gsi_p_trivial.plot");
    q_trivial.writePlotFile(200, "gsi_q_trivial.plot");
    */

    uint32_t pm, pn, qm, qn;
    p.getDegree(pm, pn);
    q.getDegree(qm, qn);
    debugl(2, "CanalSurface::computeGlobalSelfIntersectionSystem: done. n = %d, p is in BB(%d, %d), q in BB(%d, %d).\n", n, pm, pn, qm, qn );

    /* there are only two edge systems here, which are called p_edge_t0 and p_edge_t1. since Gamma
     * == Delta in this case, the edges E_x0 and E_y0 logically describe the same search problem
     * and hence only two edge systems need to be considered. since dist(0,0) = 0 and dist(1,1) = 0
     * is trivial, only the two corners (0,1) and (1,0) need to be checked, which again coincide
     * logically, since || gamma(1) - gamma(0) || = || gamma(0) - gamma(1) ||. so only one 
     * candidate non-trvial candidate point from the four corners.. either one */
    Vec3<R>                     vec_gamma_t0, vec_gamma_t1;
    BernsteinPolynomial<R, R>   gamma_t0[3], gamma_t1[3];

    vec_gamma_t0 = gamma.eval(0.0);
    vec_gamma_t1 = gamma.eval(1.0);

    /* convert vectors to constant polynomials in BB(n) / BB(m) */
    for (j = 0; j < 3; j++) {
        gamma_t0[j] = BernsteinPolynomial<R, R>(n, vec_gamma_t0[j]);
        gamma_t1[j] = BernsteinPolynomial<R, R>(n, vec_gamma_t1[j]);
    }

    /* compute two edge polynomials */
    p_edge_t0 = 
        dgamma[0].multiply(gamma[0] - gamma_t0[0]) + 
        dgamma[1].multiply(gamma[1] - gamma_t0[1]) + 
        dgamma[2].multiply(gamma[2] - gamma_t0[2]);

    p_edge_t1 = 
        dgamma[0].multiply(gamma[0] - gamma_t1[0]) + 
        dgamma[1].multiply(gamma[1] - gamma_t1[1]) + 
        dgamma[2].multiply(gamma[2] - gamma_t1[2]);

    debugl(2, "BLRCanalSurface::computeGlobalSelfIntersectionSystem: p_edge_t0: BB(%d), p_edge_t1: BB(%d). coeffs vectors follow..\n", p_edge_t0.getDegree(), p_edge_t1.getDegree());

    debugTabDec();
    debugl(1, "BLRCanalSurface::computeIntersectionSystem(): done.\n");
}

template <typename R>
void
BLRCanalSurface<R>::computeIntersectionSystem(
    BLRCanalSurface<R> const       &Delta,
    BiBernsteinPolynomial<R, R>    &p,
    BiBernsteinPolynomial<R, R>    &q,
    BernsteinPolynomial<R, R>      &p_edge_x0,
    BernsteinPolynomial<R, R>      &p_edge_x1,
    BernsteinPolynomial<R, R>      &p_edge_y0,
    BernsteinPolynomial<R, R>      &p_edge_y1) const
{
    debugl(1, "BLRCanalSurface::computeIntersectionSystem().\n");
    debugTabInc();

   /* NOTE: (this) is Gamma within this context: a const reference is used for better readability. its spine curve
    * is gamma (+ derivative dgamma). same for Delta, delta, ddelta. */
    BLRCanalSurface const  &Gamma   = (*this);
    BezierCurve<R> const   &gamma   = Gamma.spine_curve;
    BezierCurve<R> const   dgamma   = gamma.getDerivative();

    BezierCurve<R> const   &delta   = Delta.spine_curve;
    BezierCurve<R> const   ddelta   = delta.getDerivative();

    uint32_t const          m       = Gamma.getDegree(), n = Delta.getDegree();
    uint32_t                i, j;

    debugl(1, "CanalSurface::computeIntersectionSystem(): m = %d, n = %d\n", m, n);
    debugTabInc();

    /* gamma[k] - delta[k] elevated to bidegree (m, n) */
    std::vector<BiBernsteinPolynomial<R, R>>    dist_gamma_delta(3, BiBernsteinPolynomial<R, R>(m, n, 0));

    /* compute common factor, the "distance vector" in BB(m, n) */
    for (i = 0; i < m + 1; i++) {
        for (j = 0; j < n + 1; j++) {
            dist_gamma_delta[0](i, j) = Gamma.spine_curve[0][i] - Delta.spine_curve[0][j];
            dist_gamma_delta[1](i, j) = Gamma.spine_curve[1][i] - Delta.spine_curve[1][j];
            dist_gamma_delta[2](i, j) = Gamma.spine_curve[2][i] - Delta.spine_curve[2][j];
        }
    } 

    /* p is the inner product of dgamma and the distance vector dist_gamma_delta. this is done
     * component-wise. for this purpose, the components of the derivatives dgamma and ddelta are
     * converted to bivariate polynomials in BB(m-1, 0) and BB(0, n-1), respectively */
    p = dist_gamma_delta[0].multiply( PolyAlg::BernsteinConvertToBiPoly(dgamma[0], 0, true) ) +
        dist_gamma_delta[1].multiply( PolyAlg::BernsteinConvertToBiPoly(dgamma[1], 0, true) ) +
        dist_gamma_delta[2].multiply( PolyAlg::BernsteinConvertToBiPoly(dgamma[2], 0, true) );

    q = dist_gamma_delta[0].multiply( PolyAlg::BernsteinConvertToBiPoly(ddelta[0], 0, false) ) +
        dist_gamma_delta[1].multiply( PolyAlg::BernsteinConvertToBiPoly(ddelta[1], 0, false) ) +
        dist_gamma_delta[2].multiply( PolyAlg::BernsteinConvertToBiPoly(ddelta[2], 0, false) );

    uint32_t pm, pn, qm, qn;
    p.getDegree(pm, pn);
    q.getDegree(qm, qn);
    debugl(2, "CanalSurface::computeIntersectionSystem: done. m = %d, n = %d, p is in BB(%d, %d), q in BB(%d, %d).\n", m, n, pm, pn, qm, qn );

    /* compute edge polynomials */
    /* the edge polynomial p_edge_x0 are the stationary points of the distance function restricted to the edge E_x0 =
     * {x0} x [y0, y1]. these are the roots of the polynomial
     *
     * p_edge_x0(y) = ddelta(y) \cdot (delta(y) - gamma(x0))
     *
     * the other three edge polynomials are defined accordingly.
     *
     * NOTE: we're working with the magnified curves -> x_0, y_0 = 0.0, x_1, y_1 = 1.0. */
    Vec3<R>                                 vec_gamma_x0, vec_gamma_x1, vec_delta_y0, vec_delta_y1; 
    std::vector<BernsteinPolynomial<R, R>>  gamma_x0(3), gamma_x1(3), delta_y0(3), delta_y1(3);

    vec_gamma_x0    = gamma.eval(0.0);
    vec_gamma_x1    = gamma.eval(1.0);
    vec_delta_y0    = delta.eval(0.0);
    vec_delta_y1    = delta.eval(1.0);
    
    /* convert vectors to constant polynomials in BB(n) / BB(m) */
    for (j = 0; j < 3; j++) {
        gamma_x0[j] = BernsteinPolynomial<R, R>(n, vec_gamma_x0[j]);
        gamma_x1[j] = BernsteinPolynomial<R, R>(n, vec_gamma_x1[j]);

        delta_y0[j] = BernsteinPolynomial<R, R>(m, vec_delta_y0[j]);
        delta_y1[j] = BernsteinPolynomial<R, R>(m, vec_delta_y1[j]);
    }

    /* compute four edge polynomials by component-wise evaluation of dot product */
    p_edge_x0 = 
        ddelta[0].multiply(delta[0] - gamma_x0[0]) + 
        ddelta[1].multiply(delta[1] - gamma_x0[1]) + 
        ddelta[2].multiply(delta[2] - gamma_x0[2]);

    p_edge_x1 = 
        ddelta[0].multiply(delta[0] - gamma_x1[0]) + 
        ddelta[1].multiply(delta[1] - gamma_x1[1]) + 
        ddelta[2].multiply(delta[2] - gamma_x1[2]);

    p_edge_y0 = 
        dgamma[0].multiply(gamma[0] - delta_y0[0]) + 
        dgamma[1].multiply(gamma[1] - delta_y0[1]) + 
        dgamma[2].multiply(gamma[2] - delta_y0[2]);

    p_edge_y1 = 
        dgamma[0].multiply(gamma[0] - delta_y1[0]) + 
        dgamma[1].multiply(gamma[1] - delta_y1[1]) + 
        dgamma[2].multiply(gamma[2] - delta_y1[2]);

    /* p_edge_{x0,x1} are in BB(2n-1), p_edge_{y0,y1} are in BB(2m-1) */
    debugl(2, "CanalSurface::computeIntersectionSystem: edge polynomials done. m = %d, n = %d, p_edge_x0: BB(%d), p_edge_x1: BB(%d), p_edge_y0: BB(%d), p_edge_y1: BB(%d)\n", 
            m, n, p_edge_x0.getDegree(), p_edge_x1.getDegree(), p_edge_y0.getDegree(), p_edge_y1.getDegree());

    debugTabDec();
    debugl(1, "CanalSurface::computeIntersectionSystem(): done.\n");
}
