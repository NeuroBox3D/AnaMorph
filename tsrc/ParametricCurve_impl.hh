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
ParametricCurve<C2F, T, Tr>::checkDimDomain(
    ParametricCurve const  &x,
    const char* fn) const
{
    if (this->getDim() != x.getDim() || this->t0 != x.t0 || this->t1 != x.t1) {
        throw(std::string(fn) + " implement parametric curve domain mismatch exception please");
    }
}

template < typename C2F, typename T, typename Tr>
void
ParametricCurve<C2F, T, Tr>::checkEvalParameter(
    T const      &x,
    const char*   fn) const
{
    if (x < this->t0 || x > this->t1) {
        debugl(1, "ParametricCurve::checkEvalParameter(): x (%s) out of range[%s, %s]\n", std::to_string(x).c_str(), std::to_string(this->t0).c_str(), std::to_string(this->t1).c_str());
        throw(std::string(fn) + " implement parameter out of domain exception please");
    }
}

template < typename C2F, typename T, typename Tr>
void
ParametricCurve<C2F, T, Tr>::checkEvalParameter(
    T const                &x,
    std::string const      &fn) const
{
    if (x < this->t0 || x > this->t1) {
        debugl(1, "ParametricCurve::checkEvalParameter(): x (%s) out of range[%s, %s]\n", std::to_string(x).c_str(), std::to_string(this->t0).c_str(), std::to_string(this->t1).c_str());
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
: ParametricCurve<C2F, R, Vec3<R> >(x)
{
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
        debugl(5, "SpaceCurveReal::approxArcLength(): result trivial (tstart == tend). returning zero..\n");
        return 0.0;
    }
    else if (this->arclen_set && this->arclen_dt == dt && tstart == this->t0 && tend == this->t1) {
        debugl(5, "SpaceCurveReal::approxArcLength(): result cached: %10.5f. returning..\n", this->arclen);
        return this->arclen;
    }
    /* otherwise re-approximate for given settings */
    else {
        debugl(5, "SpaceCurveReal::approxArcLength(): result not trivial / cached: computing..\n");
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
        debugl(5, "SpaceCurveReal::approxArcLength(): done. result: %10.5f.\n", arclen_t);
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

template <uint32_t degree, typename R>
BezierCurve<degree, R>::BezierCurve()
    : SpaceCurveReal<pol_type, R>()
{
    /* domain is _always_ [0,1] for bezier curves. */
    this->t0    = 0;
    this->t1    = 1;

    this->d_component_functors.resize(3);
    this->d2_component_functors.resize(3);
}

template <uint32_t degree, typename R>
BezierCurve<degree, R>::BezierCurve(const std::array<pol_type, 3>& component_functors)
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

template <uint32_t degree, typename R>
BezierCurve<degree, R>::BezierCurve(const std::vector<Vec3<R> >& control_points)
{
    /* resize vectors for component functions and its derivatives */
    this->component_functors.resize(3);
    this->d_component_functors.resize(3);
    this->d2_component_functors.resize(3);

    /* set bb component polynomial coefficients and compute derivatives from control points */
    typename pol_type::coeff_type gamma_i_coeff(0);
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

template <uint32_t degree, typename R>
BezierCurve<degree, R>::BezierCurve(const this_type& x)
    : SpaceCurveReal<pol_type, R>(x)
{
    this->d_component_functors  = x.d_component_functors;
    this->d2_component_functors = x.d2_component_functors;
}

template <uint32_t degree, typename R>
BezierCurve<degree, R>&
BezierCurve<degree, R>::operator=(const this_type& x)
{
    SpaceCurveReal<pol_type, R>::operator=(x);
    this->d_component_functors  = x.d_component_functors;
    this->d2_component_functors = x.d2_component_functors;

    return *this;
}

template <uint32_t degree, typename R>
BezierCurve<degree, R>::~BezierCurve()
{
}

/* arithmetic */
template <uint32_t degree, typename R>
BezierCurve<degree, R>
BezierCurve<degree, R>::operator+(const this_type& x) const
{
    return this_type(*this) += x;
}

template <uint32_t degree, typename R>
BezierCurve<degree, R>&
BezierCurve<degree, R>::operator+=(const this_type& x)
{
    SpaceCurveReal<pol_type, R>::operator+=(x);
    return *this;
}

template <uint32_t degree, typename R>
BezierCurve<degree, R>
BezierCurve<degree, R>::operator-(const this_type& x) const
{
    return this_type(*this) -= x;
}

template <uint32_t degree, typename R>
BezierCurve<degree, R>&
BezierCurve<degree, R>::operator-=(const this_type& x)
{
    SpaceCurveReal<pol_type, R>::operator-=(x);
    return *this;
}

template <uint32_t degree, typename R>
BezierCurve<degree, R>
BezierCurve<degree, R>::operator*(const R& x) const
{
    return this_type(*this) * x;
}

template <uint32_t degree, typename R>
BezierCurve<degree, R>&
BezierCurve<degree, R>::operator*=(const R& x)
{
    SpaceCurveReal<pol_type, R>::operator*=(x);
    return *this;
}

template <uint32_t degree, typename R>
BezierCurve<degree, R>
BezierCurve<degree, R>::operator/(const R& x) const
{
    return this_type(*this) / x;
}

template <uint32_t degree, typename R>
BezierCurve<degree, R>&
BezierCurve<degree, R>::operator/=(const R& x)
{
    SpaceCurveReal<pol_type, R>::operator/=(x);
    return *this;
}

/* get derivative curve */
template <uint32_t degree, typename R>
BezierCurve<(degree>0) ? degree-1 : 0, R>
BezierCurve<degree, R>::getDerivative() const
{
    return this_deriv_type( {
                this->operator[](0).getDerivative(),
                this->operator[](1).getDerivative(),
                this->operator[](2).getDerivative()
            } );
}

/* get control points */
template <uint32_t degree, typename R>
std::list<Vec3<R> >
BezierCurve<degree, R>::getControlPoints() const
{
    std::list<Vec3<R>> control_points;
    for (uint32_t i = 0; i < degree+1; i++) {
        control_points.push_back({
                this->component_functors[0][i],
                this->component_functors[1][i],
                this->component_functors[2][i]
            });
    }
    return control_points;
}

template <uint32_t degree, typename R>
void
BezierCurve<degree, R>::split(
    const R& t,
    this_type* cleft,
    this_type* cright) const
{
    if (t < 0 || t > 1) {
        throw("BezierCurve::split(): given value of t not in [0,1].");
    }
    std::array<pol_type, 3> comps_left, comps_right;

    for (uint32_t i = 0; i < 3; i++) {
        this->operator[](i).split(t, &comps_left[i], &comps_right[i]);
    }

    if (cleft) *cleft = this_type(comps_left);
    if (cright) *cright = this_type(comps_right);
}

template <uint32_t degree, typename R>
void
BezierCurve<degree, R>::clipToInterval(const R& t0, const R& t1, this_type* cclip) const
{
    if (t0 > t1 || t0 < 0 || t0 > 1 || t1 < 0 || t1 > 1) {
        throw("BezierCurve::clipToInterval(): given interval [t0, t1] not subset of [0,1] or t1 < t0.");
    }

    if (cclip) {
        this->split(   t0,                      NULL,   cclip);
        cclip->split( (t1 - t0) / (1.0 - t0),   cclip,  NULL);
    }
}

template <uint32_t degree, typename R>
BoundingBox<R>
BezierCurve<degree, R>::getBoundingBox(uint32_t subdivision_detph) const
{
    using Aux::Numbers::inf;
    using Aux::VecMat::onesVec3;
    using Aux::VecMat::fabsVec3;
    using Aux::VecMat::minVec3;
    using Aux::VecMat::maxVec3;

    std::list<this_type>   l1 = { *this }, l2;
    std::list<this_type>  *l = &l1, *lswap = &l2;
    this_type             *c, c1, c2;
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
    std::list<this_type> &list = *l;
    Vec3<R> bb_min(inf<R>());
    Vec3<R> bb_max(-inf<R>());
    for (auto &bc : list) {
        auto cplist = bc.getControlPoints();
        for (auto &cp : cplist) {
            minVec3<R>(bb_min, bb_min, cp);
            maxVec3<R>(bb_max, bb_max, cp);
        }
    }

    /* extend bounding box by 2,5%, but no less than 1E-3, in every component. */
    return (BoundingBox<R>({ bb_min, bb_max}).extend(0.025, Vec3<R>(1E-3, 1E-3, 1E-3)));
}

template <uint32_t degree, typename R>
void
BezierCurve<degree, R>::computeRegularityPolynomial(BernsteinPolynomial<2*derivDeg, R, R>& p_reg) const
{
    p_reg = this->d_component_functors[0].square() + 
            this->d_component_functors[1].square() +
            this->d_component_functors[2].square();
}

template <uint32_t degree, typename R>
void
BezierCurve<degree, R>::computeStationaryPointDistPoly(
    Vec3<R> const              &x,
    BernsteinPolynomial<degree+derivDeg, R, R>& p) const
{
    /* convert components of x to three (constant) polynomials in BB(n), which can be subtracted
     * from gamma[j], j = 0..2 */
    pol_type x_bb[3];

    x_bb[0] = pol_type(x[0]);
    x_bb[1] = pol_type(x[1]);
    x_bb[2] = pol_type(x[2]);

    /* evaluate polynomial whose roots are stationary points as described in the thesis */
    p = this->d_component_functors[0].multiply(this->component_functors[0] - x_bb[0]) + 
        this->d_component_functors[1].multiply(this->component_functors[1] - x_bb[1]) + 
        this->d_component_functors[2].multiply(this->component_functors[2] - x_bb[2]);

    debugl(3, "CanalSurface::computeStationaryPointDistPoly(): n = %d, stationary point polynomial computed in BB(%d)\n", degree, degree);
}
