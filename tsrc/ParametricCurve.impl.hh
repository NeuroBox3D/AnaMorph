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
 *                        ParametricCurve implementation..                                                            *
 *                                                                                                                    *
 * ------------------------------------------------------------------------------------------------------------------ */

template < typename C2F, typename T, typename Tr>
void
ParametricCurve<C2F, T, Tr>::checkDimDomain(
    ParametricCurve const  &x,
    std::string const      &fn) const
{
    if (this->getDim() != x.getDim() || this->t0 != x.t0 || this->t1 != x.t1) {
        throw(fn + " implement parametric curve domain mismatch exception please");
    }
}

template < typename C2F, typename T, typename Tr>
void
ParametricCurve<C2F, T, Tr>::checkEvalParameter(
    T const                &x,
    std::string const      &fn) const
{
    if (x < this->t0 || x > this->t1) {
        debugl(0, "ParametricCurve::checkEvalParameter(): x (%s) out of range[%s, %s]\n", std::to_string(x).c_str(), std::to_string(this->t0).c_str(), std::to_string(this->t1).c_str());
        throw(fn + " implement parameter out of domain exception please");
    }
}

/* default construction */
template < typename C2F, typename T, typename Tr>
ParametricCurve<C2F, T, Tr>::ParametricCurve()
    /* explicitly initialize domain variables with mem-initializer, since these will in
     * general be POD types (usually T = double or T = std::complex). POD types are are not
     * implicitly initialized by the compiler. for the scalar types to be used,
     * initialize to zero. */
    : t0(0), t1(0)
{
}

template < typename C2F, typename T, typename Tr>
ParametricCurve<C2F, T, Tr>::ParametricCurve(
    std::vector<C2F> const         &component_functors,
    T                               t0,
    T                               t1)
{
    printf("ParametricCurve::ParametricCurve(): component_functors.size(): %zu.\n", component_functors.size());
    this->component_functors    = component_functors;
    printf("ParametricCurve::ParametricCurve(): component_functors.size(): %zu.\n", this->component_functors.size());
    this->t0                    = t0;
    this->t1                    = t1;
}

/* copy construction and assigment works also on parametric curves with different dimensions and
 * domain, whereas arithmetic is constrained to parametric curves of matching dimension and domain
 * */
template < typename C2F, typename T, typename Tr>
ParametricCurve<C2F, T, Tr>::ParametricCurve(ParametricCurve<C2F, T, Tr> const &x)
{
    this->component_functors    = x.component_functors;
    this->t0                    = x.t0;
    this->t1                    = x.t1;
}

template < typename C2F, typename T, typename Tr>
ParametricCurve<C2F, T, Tr> &
ParametricCurve<C2F, T, Tr>::operator=(ParametricCurve<C2F, T, Tr> const &x)
{
    this->component_functors    = x.component_functors;
    this->t0                    = x.t0;
    this->t1                    = x.t1;

    return (*this);
}

/* dtor */
template < typename C2F, typename T, typename Tr>
ParametricCurve<C2F, T, Tr>::~ParametricCurve()
{
}

template < typename C2F, typename T, typename Tr>
Tr
ParametricCurve<C2F, T, Tr>::eval(T const &x) const
{
    this->checkEvalParameter(x, "ParametricCurve::eval()");
    Tr res;
    res.resize(this->getDim());
    for (uint32_t i = 0; i < this->getDim(); i++) {
        res[i] = this->component_functors[i].eval(x);
    }
    return res;
}

template < typename C2F, typename T, typename Tr>
Tr
ParametricCurve<C2F, T, Tr>::eval_d(T const &x) const
{
    this->checkEvalParameter(x, "ParametricCurve::eval_d()");
    Tr  res;
    res.resize(this->getDim());
    for (uint32_t i = 0; i < this->getDim(); i++) {
        res[i] = this->component_functors[i].eval_d(x);
    }
    return res;
}

template < typename C2F, typename T, typename Tr>
Tr
ParametricCurve<C2F, T, Tr>::eval_d2(T const &x) const
{
    this->checkEvalParameter(x, "ParametricCurve::eval_d2()");
    Tr  res;
    res.resize(this->getDim());
    for (uint32_t i = 0; i < this->getDim(); i++) {
        res[i] = this->component_functors[i].eval_d2(x);
    }
    return res;
}

template < typename C2F, typename T, typename Tr>
uint32_t
ParametricCurve<C2F, T, Tr>::getDim() const
{
    return (this->component_functors.size());
}

template < typename C2F, typename T, typename Tr>
std::pair<T, T>
ParametricCurve<C2F, T, Tr>::getDomain() const
{
    return (std::pair<T, T>(this->t0, this->t1));
}

/*
template < typename C2F, typename T, typename Tr>
void
ParametricCurve<C2F, T, Tr>::setDomain(
    T const &t0,
    T const &t1)
{
    if (!(t1 < t0)) {
        this->t0 = t0;
        this->t1 = t1;
    }
    else {
        throw("ParametricCurve::setDomain(): invalid domain [t0, t1] with t0 > t1 given.");
    }
}
*/

template < typename C2F, typename T, typename Tr>
C2F &
ParametricCurve<C2F, T, Tr>::operator[](uint32_t i)
{
    if (i < this->getDim()) {
        return (this->component_functors[i]);
    }
    else {
        throw("ParametricCurve::setComponent(): required component index larger than dimension.");
    }
}

template < typename C2F, typename T, typename Tr>
C2F const &
ParametricCurve<C2F, T, Tr>::operator[](uint32_t i) const
{
    if (i < this->getDim()) {
        return (this->component_functors[i]);
    }
    else {
        throw("ParametricCurve::setComponent(): required component index larger than dimension.");
    }
}

/* arithmetic */
template < typename C2F, typename T, typename Tr>
ParametricCurve<C2F, T, Tr>
ParametricCurve<C2F, T, Tr>::operator+(ParametricCurve<C2F, T, Tr> const &x) const
{
    this->checkDimDomain(x, "ParametricCurve::operator+(ParametricCurve const &)");

    ParametricCurve<C2F, T, Tr> res(*this);
    for (uint32_t i = 0; i < this->getDim(); i++) {
        res.component_functors[i] = this->component_functors[i] + x.component_functors[i];
    }

    return res;
}
template < typename C2F, typename T, typename Tr>
ParametricCurve<C2F, T, Tr> &
ParametricCurve<C2F, T, Tr>::operator+=(ParametricCurve<C2F, T, Tr> const &x)
{
    this->checkDimDomain(x, "ParametricCurve::operator+(ParametricCurve const &)");

    for (uint32_t i = 0; i < this->getDim(); i++) {
        this->component_functors[i] = this->component_functors[i] + x.component_functors[i];
    }

    return (*this);
}

template < typename C2F, typename T, typename Tr>
ParametricCurve<C2F, T, Tr>
ParametricCurve<C2F, T, Tr>::operator-(ParametricCurve<C2F, T, Tr> const &x) const
{
    this->checkDimDomain(x, "ParametricCurve::operator-(ParametricCurve const &)");

    ParametricCurve<C2F, T, Tr> res(*this);
    for (uint32_t i = 0; i < this->getDim(); i++) {
        res.component_functors[i] = this->component_functors[i] - x.component_functors[i];
    }

    return res;
}

template < typename C2F, typename T, typename Tr>
ParametricCurve<C2F, T, Tr> &
ParametricCurve<C2F, T, Tr>::operator-=(ParametricCurve<C2F, T, Tr> const &x)
{
    this->checkDimDomain(x, "ParametricCurve::operator-(ParametricCurve const &)");

    for (uint32_t i = 0; i < this->getDim(); i++) {
        this->component_functors[i] = this->component_functors[i] - x.component_functors[i];
    }

    return (*this);
}

template < typename C2F, typename T, typename Tr>
ParametricCurve<C2F, T, Tr>
ParametricCurve<C2F, T, Tr>::operator*(T const &x) const
{
    ParametricCurve<C2F, T, Tr> res(*this);
    for (uint32_t i = 0; i < this->getDim(); i++) {
        res.component_functors[i] = this->component_functors[i] * x;
    }

    return res;
}

template < typename C2F, typename T, typename Tr>
ParametricCurve<C2F, T, Tr> &
ParametricCurve<C2F, T, Tr>::operator*=(T const &x)
{
    for (uint32_t i = 0; i < this->getDim(); i++) {
        this->component_functors[i] = this->component_functors[i] * x;
    }

    return (*this);
}

template < typename C2F, typename T, typename Tr>
ParametricCurve<C2F, T, Tr>
ParametricCurve<C2F, T, Tr>::operator/(T const &x) const
{
    ParametricCurve<C2F, T, Tr> res(*this);
    for (uint32_t i = 0; i < this->getDim(); i++) {
        res.component_functors[i] = this->component_functors[i] / x;
    }

    return res;
}

template < typename C2F, typename T, typename Tr>
ParametricCurve<C2F, T, Tr> &
ParametricCurve<C2F, T, Tr>::operator/=(T const &x)
{
    for (uint32_t i = 0; i < this->getDim(); i++) {
        this->component_functors[i] = this->component_functors[i] / x;
    }

    return (*this);
}


/* ------------------------------------------------------------------------------------------------------------------ *
 *                                                                                                                    *
 *                        SpaceCurveReal implementation..                                                             *
 *                                                                                                                    *
 * ------------------------------------------------------------------------------------------------------------------ */

template <typename C2F, typename R>
SpaceCurveReal<C2F, R>::SpaceCurveReal()
    : ParametricCurve<C2F, R, Vec3<R>>()
{
    this->arclen_set    = false;
    this->arclen        = 0;
    this->arclen_dt     = 0;
    
    /* resize component functors to create three default-constructed objects of component function
     * template type C2F. */
    this->component_functors.resize(3);
}


template <typename C2F, typename R>
SpaceCurveReal<C2F, R>::SpaceCurveReal(
    std::array<C2F, 3> const   &component_functors,
    R const                    &t0,
    R const                    &t1)
        :   ParametricCurve<C2F, R, Vec3<R>>(
                { component_functors[0], component_functors[1], component_functors[2] },
                t0,
                t1
            ) 
{
    this->arclen_set            = false;
    this->arclen                = 0;
    this->arclen_dt             = 0;
}


template <typename C2F, typename R>
SpaceCurveReal<C2F, R>::SpaceCurveReal(SpaceCurveReal<C2F, R> const &x)
{
    ParametricCurve<C2F, R, Vec3<R>>::operator=(x);
    this->arclen_set    = x.arclen_set;
    this->arclen        = x.arclen;
    this->arclen_dt     = x.arclen_dt;
}

template <typename C2F, typename R>
SpaceCurveReal<C2F, R> &
SpaceCurveReal<C2F, R>::operator=(SpaceCurveReal<C2F, R> const &x)
{
    ParametricCurve<C2F, R, Vec3<R>>::operator=(x);
    this->arclen_set    = x.arclen_set;
    this->arclen        = x.arclen;
    this->arclen_dt     = x.arclen_dt;

    return (*this);
}

template <typename C2F, typename R>
SpaceCurveReal<C2F, R>::~SpaceCurveReal()
{
}

/* arithmetic */
template <typename C2F, typename R>
SpaceCurveReal<C2F, R>
SpaceCurveReal<C2F, R>::operator+(SpaceCurveReal<C2F, R> const &x) const
{
    return (SpaceCurveReal<C2F, R>(*this) += x);
}

template <typename C2F, typename R>
SpaceCurveReal<C2F, R> &
SpaceCurveReal<C2F, R>::operator+=(SpaceCurveReal<C2F, R> const &x)
{
    ParametricCurve<C2F, R, Vec3<R>>::operator+=(x);
    return (*this);
}

template <typename C2F, typename R>
SpaceCurveReal<C2F, R>
SpaceCurveReal<C2F, R>::operator-(SpaceCurveReal<C2F, R> const &x) const
{
    return (SpaceCurveReal<C2F, R>(*this) -= x);
}

template <typename C2F, typename R>
SpaceCurveReal<C2F, R> &
SpaceCurveReal<C2F, R>::operator-=(SpaceCurveReal<C2F, R> const &x)
{
    ParametricCurve<C2F, R, Vec3<R>>::operator-=(x);
    return (*this);
}

template <typename C2F, typename R>
SpaceCurveReal<C2F, R>
SpaceCurveReal<C2F, R>::operator*(R const &x) const
{
    return (SpaceCurveReal<C2F, R>(*this) * x);
}

template <typename C2F, typename R>
SpaceCurveReal<C2F, R> &
SpaceCurveReal<C2F, R>::operator*=(R const &x)
{
    ParametricCurve<C2F, R, Vec3<R>>::operator*=(x);
    return (*this);
}

template <typename C2F, typename R>
SpaceCurveReal<C2F, R>
SpaceCurveReal<C2F, R>::operator/(R const &x) const
{
    return (SpaceCurveReal<C2F, R>(*this) / x);
}

template <typename C2F, typename R>
SpaceCurveReal<C2F, R> &
SpaceCurveReal<C2F, R>::operator/=(R const &x)
{
    ParametricCurve<C2F, R, Vec3<R>>::operator/=(x);
    return (*this);
}

/* frenet-serret frame */
template <typename C2F, typename R>
void
SpaceCurveReal<C2F, R>::getFrenetFrame(
    R const    &t,
    Vec3<R>    &x,
    Vec3<R>    &y,
    Vec3<R>    &z) const
{
    x = this->eval_d(t);
    x.normalize();

    z = this->eval_d(t).cross(this->eval_d2(t));
    z.normalize();

    y = z.cross(x);
    y.normalize();
}

/* render frame */
template <typename C2F, typename R>
void
SpaceCurveReal<C2F, R>::getRenderFrame(
    R const    &t,
    Vec3<R>     r,
    Vec3<R>    &x,
    Vec3<R>    &y,
    Vec3<R>    &z) const
{
    x = this->eval_d(t);
    x.normalize();

    y = r.cross(x);
    y.normalize();

    z = x.cross(y);
    z.normalize();
}

template <typename C2F, typename R>
R
SpaceCurveReal<C2F, R>::approxArcLength(
    R const    &tstart,
    R const    &tend,
    R const    &dt) const
{
    R       t, arclen_t = 0;
    Vec3<R> last, current, disp;

    /* if tstart == tend or tstart = this->t0 && tend = this->t0 (and dt matched this->arclen_dt),
     * return 0.0 and this->arclength, respectively */
    if (tstart == tend) {
        debugl(4, "SpaceCurveReal::approxArcLength(): result trivial (tstart == tend). returning zero..\n");
        return 0.0;
    }
    else if (this->arclen_set && this->arclen_dt == dt && tstart == this->t0 && tend == this->t1) {
        debugl(4, "SpaceCurveReal::approxArcLength(): result cached: %10.5f. returning..\n", this->arclen);
        return this->arclen;
    }
    /* otherwise re-approximate for given settings */
    else {
        debugl(4, "SpaceCurveReal::approxArcLength(): result not trivial / cached: computing..\n");
        current = this->eval(tstart);
        debugTabInc();
        for (t = tstart + dt; t < tend; t += dt) {
            last        = current;
            current     = this->eval(t);
            disp        = current - last;
            arclen_t   += disp.len2();
            debugl(5, "t0: %f, t1: %f, t = %f, dt = %f, increment: %f\n", this->t0, this->t1, t, dt, disp.len2() );
        }
        debugTabDec();
        arclen_t += (this->eval(tend) - last).len2();
        debugl(4, "SpaceCurveReal::approxArcLength(): done. result: %10.5f.\n", arclen_t);
        return arclen_t;
    }
}

template <typename C2F, typename R>
void
SpaceCurveReal<C2F, R>::updateArcLength(R const &dt)
{
    if (!this->arclen_set || this->arclen_dt != dt) {
        this->arclen        = this->approxArcLength(dt, this->t0, this->t1);
        this->arclen_set    = true;
        this->arclen_dt     = dt;
    }
}


/* ------------------------------------------------------------------------------------------------------------------ *
 *                                                                                                                    *
 *                        implementation of real bezier curve class....                                               *
 *                                                                                                                    *
 * ------------------------------------------------------------------------------------------------------------------ */

template<typename R>
BezierCurve<R>::BezierCurve()
    : SpaceCurveReal<BernsteinPolynomial<R, R>, R>()
{
    /* domain is _always_ [0,1] for bezier curves. */
    this->t0    = 0;
    this->t1    = 1;

    this->d_component_functors.resize(3);
    this->d2_component_functors.resize(3);
}

template<typename R>
BezierCurve<R>::BezierCurve(
    std::array<BernsteinPolynomial<R, R>, 3> const &component_functors)
{
    /* check if all component polynomials have same degree */
    if (component_functors[0].getDegree() == component_functors[1].getDegree() &&
            component_functors[0].getDegree() == component_functors[2].getDegree()) {

        /* resize all component function vectors, assign component functions, compute and store
         * derivatives for faster evaluation */
        this->component_functors.resize(3);
        this->d_component_functors.resize(3);
        this->d2_component_functors.resize(3);

        for (int i = 0; i < 3; i++) {
            this->component_functors[i]     = component_functors[i];
            this->d_component_functors[i]   = this->component_functors[i].getDerivative();
            this->d2_component_functors[i]  = this->d_component_functors[i].getDerivative();
        }
    }
    else {
        throw("BezierCurve::BezierCurve(): given component polynomials not all of the same degree.");
    }

    /* set domain [0,1]*/
    this->t0    = 0;
    this->t1    = 1;
}

template<typename R>
BezierCurve<R>::BezierCurve(std::vector<Vec3<R>> const &control_points)
{
    /* resize vectors for component functions and its derivatives */
    this->component_functors.resize(3);
    this->d_component_functors.resize(3);
    this->d2_component_functors.resize(3);

    /* set bb component polynomial coefficients and compute derivatives from control points */
    Vector<R> gamma_i_coeff(4);
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 4; j++) {
            gamma_i_coeff[j] = control_points[j][i];
        }

        this->component_functors[i].setCoeffs(gamma_i_coeff);
        this->d_component_functors[i]   = this->component_functors[i].getDerivative();;
        this->d2_component_functors[i]  = this->d_component_functors[i].getDerivative();;
    }

    /* set domain */
    this->t0    = 0;
    this->t1    = 1;
}

template<typename R>
BezierCurve<R>::BezierCurve(BezierCurve const &x)
    : SpaceCurveReal<BernsteinPolynomial<R, R>, R>(x)
{
    this->d_component_functors  = x.d_component_functors;
    this->d2_component_functors = x.d2_component_functors;
}

template<typename R>
BezierCurve<R> &
BezierCurve<R>::operator=(BezierCurve const &x)
{
    SpaceCurveReal<BernsteinPolynomial<R, R>, R>::operator=(x);
    this->d_component_functors  = x.d_component_functors;
    this->d2_component_functors = x.d2_component_functors;

    return (*this);
}

template<typename R>
BezierCurve<R>::~BezierCurve()
{
}

/* arithmetic */
template <typename R>
BezierCurve<R>
BezierCurve<R>::operator+(BezierCurve<R> const &x) const
{
    return (BezierCurve<R>(*this) += x);
}

template <typename R>
BezierCurve<R> &
BezierCurve<R>::operator+=(BezierCurve<R> const &x)
{
    SpaceCurveReal< BernsteinPolynomial<R, R>, R>::operator+=(x);
    return (*this);
}

template <typename R>
BezierCurve<R>
BezierCurve<R>::operator-(BezierCurve<R> const &x) const
{
    return (BezierCurve<R>(*this) -= x);
}

template <typename R>
BezierCurve<R> &
BezierCurve<R>::operator-=(BezierCurve<R> const &x)
{
    SpaceCurveReal< BernsteinPolynomial<R, R>, R>::operator-=(x);
    return (*this);
}

template <typename R>
BezierCurve<R>
BezierCurve<R>::operator*(R const &x) const
{
    return (BezierCurve<R>(*this) * x);
}

template <typename R>
BezierCurve<R> &
BezierCurve<R>::operator*=(R const &x)
{
    SpaceCurveReal< BernsteinPolynomial<R, R>, R>::operator*=(x);
    return (*this);
}

template <typename R>
BezierCurve<R>
BezierCurve<R>::operator/(R const &x) const
{
    return (BezierCurve<R>(*this) / x);
}

template <typename R>
BezierCurve<R> &
BezierCurve<R>::operator/=(R const &x)
{
    SpaceCurveReal< BernsteinPolynomial<R, R>, R>::operator/=(x);
    return (*this);
}

/* get derivative curve */
template<typename R>
BezierCurve<R>
BezierCurve<R>::getDerivative() const
{
    return  BezierCurve<R>( {
                this->operator[](0).getDerivative(),
                this->operator[](1).getDerivative(),
                this->operator[](2).getDerivative()
            } );
}

/* get control points */
template<typename R>
std::list<Vec3<R>>
BezierCurve<R>::getControlPoints() const
{
    std::list<Vec3<R>> control_points;
    for (uint32_t i = 0; i < this->getDegree() + 1; i++) {
        control_points.push_back({
                this->component_functors[0][i],
                this->component_functors[1][i],
                this->component_functors[2][i]
            });
    }
    return control_points;
}

template<typename R>
void
BezierCurve<R>::split(
    R const        &t,
    BezierCurve<R> *cleft,
    BezierCurve<R> *cright) const
{
    if (t < 0 || t > 1) {
        throw("BezierCurve::split(): given value of t not in [0,1].");
    }
    std::array<BernsteinPolynomial<R, R>, 3> comps_left, comps_right;

    for (uint32_t i = 0; i < 3; i++) {
        this->operator[](i).split(t, &comps_left[i], &comps_right[i]);
    }

    if (cleft)  *cleft   = BezierCurve<R>(comps_left);
    if (cright) *cright  = BezierCurve<R>(comps_right);
}

template<typename R>
void
BezierCurve<R>::clipToInterval(
    R const        &t0,
    R const        &t1,
    BezierCurve<R> *cclip) const
{
    if (t0 > t1 || t0 < 0 || t0 > 1 || t1 < 0 || t1 > 1) {
        throw("BezierCurve::clipToInterval(): given interval [t0, t1] not subset of [0,1] or t1 < t0.");
    }

    if (cclip) {
        this->split(   t0,                      NULL,   cclip);
        cclip->split( (t1 - t0) / (1.0 - t0),   cclip,  NULL);
    }
}

template<typename R>
BoundingBox<R>
BezierCurve<R>::getBoundingBox(uint32_t subdivision_detph) const
{
    using Aux::Numbers::inf;
    using Aux::VecMat::onesVec3;
    using Aux::VecMat::fabsVec3;
    using Aux::VecMat::minVec3;
    using Aux::VecMat::maxVec3;

    std::list<BezierCurve<R>>   l1 = { *this }, l2;
    std::list<BezierCurve<R>>  *l = &l1, *lswap = &l2;
    BezierCurve<R>             *c, c1, c2; 
    for (uint32_t i = 0; i < subdivision_detph; i++) {
        /* clear swap list */
        lswap->clear();

        /* subdivide every bezier curve in l at 0.5 and push two resulting curves into swap list */
        for (auto cit = l->begin(); cit != l->end(); ++cit) {
            c = &(*cit);
            c->split(0.5, &c1, &c2);
            lswap->push_back(c1);
            lswap->push_back(c2);
        }

        /* swap pointers */
        std::swap(l, lswap);
    }

    /* l now contains (potentially a lot) of BezierCurve's from the iterated subdivision. get minimum / maximum
     * component vectgor from all control points */
    std::list<BezierCurve<R>> &list = *l;
    Vec3<R> bb_min                  =  onesVec3<R>() * inf<R>();
    Vec3<R> bb_max                  =  onesVec3<R>() * (-inf<R>());
    for (auto &bc : list) {
        auto cplist = bc.getControlPoints();
        for (auto &cp : cplist) {
            bb_min = minVec3<R>(bb_min, cp);
            bb_max = maxVec3<R>(bb_max, cp);
        }
    }

    /* extend bounding box by 2,5%, but no less than 1E-3, in every component. */
    return (BoundingBox<R>({ bb_min, bb_max}).extend(0.025, Vec3<R>(1E-3, 1E-3, 1E-3)));
}

template<typename R>
uint32_t
BezierCurve<R>::getDegree() const
{
    return (this->component_functors[0].getDegree());
}

template<typename R>
void
BezierCurve<R>::computeRegularityPolynomial(BernsteinPolynomial<R, R> &p_reg) const
{
    p_reg = this->d_component_functors[0].square() + 
            this->d_component_functors[1].square() +
            this->d_component_functors[2].square();
}

template<typename R>
void
BezierCurve<R>::computeStationaryPointDistPoly(
    Vec3<R> const              &x,
    BernsteinPolynomial<R, R>  &p) const
{
    uint32_t n = this->getDegree();

    /* convert components of x to three (constant) polynomials in BB(n), which can be subtracted
     * from gamma[j], j = 0..2 */
    BernsteinPolynomial<R, R> x_bb[3];

    x_bb[0] = BernsteinPolynomial<R, R>(n, x[0]);
    x_bb[1] = BernsteinPolynomial<R, R>(n, x[1]);
    x_bb[2] = BernsteinPolynomial<R, R>(n, x[2]);

    /* evaluate polynomial whose roots are stationary points as described in the thesis */
    p = this->d_component_functors[0].multiply(this->component_functors[0] - x_bb[0]) + 
        this->d_component_functors[1].multiply(this->component_functors[1] - x_bb[1]) + 
        this->d_component_functors[2].multiply(this->component_functors[2] - x_bb[2]);

    debugl(2, "CanalSurface::computeStationaryPointDistPoly(): n = %d, stationary point polynomial computed in BB(%d)\n", n, p.getDegree() );
}
