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

#include "common.hh"

#if 0
/* generate B_k^n(t) in BB(n) => delta_{kn} coefficient vector.. */
template <typename F, typename R>
BernsteinPolynomial<F, R>
BB(uint32_t n, uint32_t k);

/* deCasteljau's algorithm for bezier curve evaluation and splitting adapted for polynomials, i.e. only the "y"
 * coordinates of the control points are given (the "x'-coordinate equals the given splitpoint t).
 * with Vec2, this would of course work for general bezier curves, but that's unnecessary for its use here. */
template <typename F, typename R>
F
deCasteljau(
    uint32_t                degree,
    Vector<F> const        &coeff,
    R const                &t)
{
    uint32_t i, k;

    /* initialize Pi(0) to bezier coefficients */
    Vector<F>   Pi = coeff;
    F           plast_i;

    for (k = 1; k <= degree; k++) {
        /* save old value of Pi[0] as plast_i for first iteration */
        plast_i = Pi[0];

        /* compute new values of Pi, the old "lookahead" value plast_i is always saved */
        for (i = 0; i <= degree - k; i++) {
            Pi[i]   = (1.0 - t) * plast_i + t*Pi[i + 1];

            /* save old value of p[i+1] as pilast_i for next iteration */
            plast_i = Pi[i + 1];
        }
    }
    return Pi[0];
}

template <typename F, typename R>
void
deCasteljauSplit(
    uint32_t            degree,
    Vector<F> const    &coeff,
    R const            &t,
    Vector<F>          &coeff_left,
    Vector<F>          &coeff_right)
{
    uint32_t i, k;

    /* initialize Pi to b as above */
    Vector<F>       Pi = coeff;
    F               plast_i;

    /* resize and zero vectors for "coefficients", i.e. control point y-values, of left and right
     * split Bernstein polynomial. this works as expected also for F = std::complex<T>, since
     * std::complex<> has a matching assignment operator that takes a numeric value (in this case of
     * type T) as the REAL part, setting the complex part to 0.0 (T must be an scalar type). thus
     * this sets all coefficients to complex zero (0.0, 0.0), as desired. */
    coeff_left.assign(degree + 1, 0.0);
    coeff_right.assign(degree + 1, 0.0);

    coeff_left[0]       = Pi[0];
    coeff_right[degree] = Pi[degree];

    for (k = 1; k <= degree; k++) {
        /* save old value of Pi[0] as plast_i for first iteration */
        plast_i = Pi[0];

        /* compute new values of Pi, the old "lookahead" value plast_i is always saved */
        for (i = 0; i <= degree - k; i++) {
            Pi[i]   = (1.0 - t) * plast_i + t*Pi[i + 1];

            /* save old value of p[i+1] as pilast_i for next iteration */
            plast_i = Pi[i + 1];
        }

        /* Pi contains Pi(k). store relevant values for left and right polynomial */
        coeff_left[k]           = Pi[0];
        coeff_right[degree - k] = Pi[degree - k];
    }
}

/* ----------------------------------------------------------------------------------------------------------------- *
 *
 *                        implementation of abstract univariate polynomial class....                                  
 *
 * ----------------------------------------------------------------------------------------------------------------- */

template <typename F, typename R>
void
Polynomial<F, R>::checkDegree(
    Polynomial<F, R> const &q,
    std::string const      &fn) const
{
    if (this->getDegree() != q.getDegree()) {
        throw(fn + ": degree mismatch. operation can only performed on polynomials of matching degree.");
    }
}

template <typename F, typename R>
Polynomial<F, R>::Polynomial()
    /* default-initialize POD member degree, Vector<F> coeff is default-initalized before ctor
     * body is executed. */
    : degree(0)
{
    this->coeff.assign(1, 0);
}

template <typename F, typename R>
Polynomial<F, R>::Polynomial(
    uint32_t    degree,
    F const    &x) : degree(degree)
{
    this->coeff.assign(degree + 1, x);
}

template <typename F, typename R>
Polynomial<F, R>::Polynomial(Vector<F> const &coeff)
{
    this->degree    = coeff.size() - 1;
    this->coeff     = coeff;
}

template <typename F, typename R>
Polynomial<F, R>::Polynomial(Polynomial<F, R> const &q)
{
    this->degree    = q.degree;
    this->coeff     = q.coeff;
}

template <typename F, typename R>
Polynomial<F, R> &
Polynomial<F, R>::operator=(Polynomial const &q)
{
    this->degree = q.degree;
    this->coeff  = q.coeff;
    return (*this);
}

template <typename F, typename R>
Polynomial<F, R>::~Polynomial()
{
}

template <typename F, typename R>
uint32_t
Polynomial<F, R>::getDegree() const
{
    return this->degree;
}

template <typename F, typename R>
Vector<F>
Polynomial<F, R>::getCoeffs() const
{
    return (this->coeff);
}

template <typename F, typename R>
void
Polynomial<F, R>::setCoeffs(
    Vector<F> const &coeff_arg)
{
    this->degree    = coeff_arg.size() - 1;
    this->coeff     = coeff_arg;
}

template <typename F, typename R>
F
Polynomial<F, R>::getMaxAbsCoeff() const
{
    F cmax = std::abs(this->coeff[0]);
    F abs_i;
    for (uint32_t i = 1; i < this->degree + 1; i++) {
        abs_i = std::abs(this->coeff[i]);
        if (abs_i > cmax) {
            cmax = abs_i;
        }
    }
    return cmax;
}

template <typename F, typename R>
void
Polynomial<F, R>::zero()
{
    this->initConstant(0);
}

template <typename F, typename R>
F &
Polynomial<F, R>::operator()(uint32_t cidx)
{
    return (this->coeff[cidx]);
}

template <typename F, typename R>
F
Polynomial<F, R>::operator()(uint32_t cidx) const
{
    return (this->coeff[cidx]);
}

template <typename F, typename R>
F &
Polynomial<F, R>::operator[](uint32_t cidx)
{
    return (this->operator()(cidx));
}

template <typename F, typename R>
F
Polynomial<F, R>::operator[](uint32_t cidx) const
{
    return (this->operator()(cidx));
}

template <typename F, typename R>
Polynomial<F, R> &
Polynomial<F, R>::operator+=(Polynomial<F, R> const &q)
{
    this->checkDegree(q, "Polynomial<F, R>::operator+=(): ");
    for (uint32_t i = 0; i < this->degree + 1; i++) {
        this->coeff[i] = this->coeff[i] + q[i];
    }
    return (*this);
}

template <typename F, typename R>
Polynomial<F, R> &
Polynomial<F, R>::operator-=(Polynomial<F, R> const &q)
{
    this->checkDegree(q, "Polynomial<F, R>::operator-=(): ");
    for (uint32_t i = 0; i < this->degree + 1; i++) {
        this->coeff[i] = this->coeff[i] - q[i];
    }
    return (*this);
}

template <typename F, typename R>
Polynomial<F, R> &
Polynomial<F, R>::operator*=(F const &x) 
{
    for (uint32_t i = 0; i < this->degree + 1; i++) {
        this->coeff[i] = this->coeff[i] * x;
    }
    return (*this);
}

template <typename F, typename R>
Polynomial<F, R> &
Polynomial<F, R>::operator/=(F const &x) 
{
    for (uint32_t i = 0; i < this->degree + 1; i++) {
        this->coeff[i] = this->coeff[i] / x;
    }
    return (*this);
}

template <typename F, typename R>
F
Polynomial<F, R>::operator*(Polynomial<F, R> const &q) const
{
    /* check if degree n matches */
    this->checkDegree(q, std::string("Polynomial<F, R>::operator*(Polynomial<F, R> const &q) (inner / scalar product)"));

    /* p = const reference to (this) polynomial */
    const Polynomial<F, R> &p = (*this);
    int                     i, j, n = p.getDegree(); // == q.getDegree()
    F                       res = 0;

    for (i = 0; i < n + 1; i++) {
        for (j = 0; j < n + 1; j++) {
            /* use virtual method getBasisInnerProduct, which is pure-virtual in Polynomial<F, R>
             * and must be implemented by all derived classes. this enables the reuse of this method
             * from derived classes through late-binding. due to type safety however, every derived
             * class needs to reimplement operator* as a wrapper, which calls this method as
             * Polynomial<F, R>::operator*(..) (or implement a completely custom one obviously). */
            res += p(i)*q(j)*this->getBasisInnerProduct(n, n, i, j);
        }
    }
    return res;
}


/* ----------------------------------------------------------------------------------------------------------------- *
 *
 * implementation of univariate power basis (aka monomial basis) polynomial class, horner scheme for evaluation, etc.
 *
 * ----------------------------------------------------------------------------------------------------------------- */

/* the power basis inner products do not depend on the degree of the representation, but only on the degrees of the
 * monomials: x^i * x^j = \int_0^1{x^{i+j}dx} = 1 / (i + j + 1) 
 *
 * => use only 2d tensor (array) instead of 4d array with custom init method.
 *
 * note that for BernsteinPolynomial, the inner products do depend on the degree of the representation and hence the 4d
 * array is really necessary. for legendre polynomials on the other hand, even a 1d array will suffice due to
 * orthogonality. */
template<typename F, typename R>
void
PowerPolynomial<F, R>::initPowerBasisInnerProducts(uint32_t max_dim)
{
    static bool recompute = true;

    if (max_dim > power_basis_inner_products_max_dim) {
        recompute = true;
    }

    if (recompute) {

        uint32_t alloc_dim = std::max(max_dim, POWERPOLY_BASIS_IP_DEFAULT_SIZE);
        power_basis_inner_products.resize( { alloc_dim + 1, alloc_dim + 1 } );

        debugl(0, "(static) PowerPolynomial::initPowerBasisInnerProducts(): alloc_dim: %u\n", alloc_dim);

        for (uint32_t i = 0; i < alloc_dim + 1; i++) {
            for (uint32_t j = 0; j < alloc_dim + 1; j++) {
                power_basis_inner_products( {i, j} ) = computePowerBasisInnerProduct(i, j);
            }
        }

        power_basis_inner_products_max_dim  = alloc_dim;
        recompute                           = false;

        debugl(0, "(static) PowerPolynomial::initPowerBasisInnerProducts(): done..\n");
    }
}

/* indices m and n for representation degree are irrelevant: inner products depends only on i and j in
 * integral over x^i*x^j = x^{i + j}. */
template<typename F, typename R>
F
PowerPolynomial<F, R>::computePowerBasisInnerProduct(
    uint32_t i,  
    uint32_t j)
{
    debugl(4, "PowerPolynomial()::computePowerBasisInnerProduct(%d, %d)\n", i, j);

    /* int_0^1{x^{i+j}dx}  = 1 / (i+j+1)*/
    return ( (F)1 / ((F)(i + j + 1)) );
}

/* return pre-computed inner product value from static PowerPolynomial<F, R>::power_basis_inner_products. if value has
 * not been precomputed yet, do it now if dynamic recomputation is enabled */
template <typename F, typename R>
F
PowerPolynomial<F, R>::getPowerBasisInnerProduct(
    uint32_t i,  
    uint32_t j)
{
    uint32_t const n = std::max(i, j);
    if (PowerPolynomial<F, R>::power_basis_inner_products_max_dim < n) {
        if (PowerPolynomial<F, R>::power_basis_inner_products_mutable) {
            PowerPolynomial<F, R>::initPowerBasisInnerProducts(n);
        }
        else {
            throw("PowerPolynomial::getPowerBasisInnerProduct(): basis polynomial indices out of range and dynamic recomputation disabled.");
        }
    }
    return PowerPolynomial<F, R>::power_basis_inner_products( {i, j} );
}

template <typename F, typename R>
PowerPolynomial<F, R>::PowerPolynomial() : Polynomial<F, R>()
{
    if (PowerPolynomial<F, R>::power_basis_inner_products_max_dim == 0) {
        if (PowerPolynomial<F, R>::power_basis_inner_products_mutable) {
            PowerPolynomial<F, R>::initPowerBasisInnerProducts();
        }
        else {
            throw("PowerPolynomial: not all necessary inner products not initialized and dynamic inner product recomputation disabled.");
        }
    }
}

template <typename F, typename R>
PowerPolynomial<F, R>::PowerPolynomial(uint32_t degree, F const &x) : Polynomial<F, R>(degree, 0)
{
    if (PowerPolynomial<F, R>::power_basis_inner_products_max_dim < degree) {
        if (PowerPolynomial<F, R>::power_basis_inner_products_mutable) {
            PowerPolynomial<F, R>::initPowerBasisInnerProducts(degree);
        }
        else {
            throw("PowerPolynomial: not all necessary inner products not initialized and dynamic inner product recomputation disabled.");
        }
    }

    this->coeff[0] = x;
}

template <typename F, typename R>
PowerPolynomial<F, R>::PowerPolynomial(Vector<F> const &coeff) : Polynomial<F, R>(coeff)
{
    uint32_t const n = coeff.size() - 1;
    if (PowerPolynomial<F, R>::power_basis_inner_products_max_dim < n) {
        if (PowerPolynomial<F, R>::power_basis_inner_products_mutable) {
            PowerPolynomial<F, R>::initPowerBasisInnerProducts(n);
        }
        else {
            throw("PowerPolynomial: not all necessary inner products not initialized and dynamic inner product recomputation disabled.");
        }
    }
}

template <typename F, typename R>
PowerPolynomial<F, R>::PowerPolynomial(PowerPolynomial<F, R> const &q) : Polynomial<F, R>(q)
{
}

template <typename F, typename R>
PowerPolynomial<F, R> &
PowerPolynomial<F, R>::operator=(PowerPolynomial const &q)
{
    Polynomial<F, R>::operator=(q);
    return (*this);
}

template <typename F, typename R>
PowerPolynomial<F, R>::~PowerPolynomial()
{
}

/* set coefficients to vector. this changes degree, so call initPowerBasisInnerProducts() on changed
 * degree, which will allocate more inner products if need be */
template <typename F, typename R>
void
PowerPolynomial<F, R>::setCoeffs(Vector<F> const &coeff)
{
    Polynomial<F, R>::setCoeffs(coeff);

    /* update monomial basis inner products if necessary, degree has been changed manually outside
     * constructor. */
    uint32_t const n = this->getDegree();
    if (PowerPolynomial<F, R>::power_basis_inner_products_max_dim < n) {
        if (PowerPolynomial<F, R>::power_basis_inner_products_mutable) {
            PowerPolynomial<F, R>::initPowerBasisInnerProducts(n);
        }
        else {
            throw("PowerPolynomial: not all necessary inner products not initialized and dynamic inner product recomputation disabled.");
        }
    }
}

template <typename F, typename R>
void
PowerPolynomial<F, R>::initConstant(F const &alpha)
{
    this->coeff.assign(this->getDegree() + 1, 0);
    this->coeff[0] = alpha;
}

template <typename F, typename R>
F
PowerPolynomial<F, R>::eval(R const &x) const
{
    const uint32_t n = this->getDegree();

    int     i;
    F       p_x = this->coeff[n];

    for (i = n - 1; i >= 0; i--) {
        p_x = p_x*x + this->coeff[i];
    }
    return p_x;
}

template <typename F, typename R>
F
PowerPolynomial<F, R>::eval_d(R const &x) const
{
    const uint32_t n = this->getDegree();
    if (this->degree > 0) {
        int i;
        F   dp_x = this->coeff[n] * (R)n;

        for (i = n - 1; i >= 1; i--) {
            dp_x = dp_x*x + this->coeff[i] * (R)i;
        }
        return dp_x;
    }
    /* if degree is smaller than 1, return 0 */
    else return 0;
}

template <typename F, typename R>
F
PowerPolynomial<F, R>::eval_d2(R const &x) const
{
    const uint32_t n = this->getDegree();
    if (n > 1) {
        int i;
        F   dp2_x = (R)(n*(n - 1)) * this->coeff[n];

        for (i = n - 1; i >= 2; i--) {
            dp2_x = dp2_x*x + (R)(i*(i-1)) * this->coeff[i];
        }
        return dp2_x;
    }
    else {
        return 0;
    }
}

/* arithmetic */
template <typename F, typename R>
PowerPolynomial<F, R>
PowerPolynomial<F, R>::operator+(PowerPolynomial<F, R> const &q) const
{
    return ( PowerPolynomial(this->coeff) += q );
}

template <typename F, typename R>
PowerPolynomial<F, R> &
PowerPolynomial<F, R>::operator+=(PowerPolynomial<F, R> const &q)
{
    Polynomial<F, R>::operator+=(q);
    return (*this);
}

template <typename F, typename R>
PowerPolynomial<F, R>
PowerPolynomial<F, R>::operator-(PowerPolynomial<F, R> const &q) const
{
    return ( PowerPolynomial(this->coeff) -= q );
}

template <typename F, typename R>
PowerPolynomial<F, R> &
PowerPolynomial<F, R>::operator-=(PowerPolynomial<F, R> const &q)
{
    Polynomial<F, R>::operator-=(q);
    return (*this);
}

template <typename F, typename R>
PowerPolynomial<F, R>
PowerPolynomial<F, R>::operator*(F const &x) const
{
    return ( PowerPolynomial(this->coeff) *= x );
}

template <typename F, typename R>
PowerPolynomial<F, R> &
PowerPolynomial<F, R>::operator*=(F const &x)
{
    Polynomial<F, R>::operator*=(x);
    return (*this);
}

template <typename F, typename R>
PowerPolynomial<F, R>
PowerPolynomial<F, R>::operator/(F const &x) const
{
    return ( PowerPolynomial(this->coeff) /= x );
}

template <typename F, typename R>
PowerPolynomial<F, R> &
PowerPolynomial<F, R>::operator/=(F const &x)
{
    Polynomial<F, R>::operator/=(x);
    return (*this);
}

/* inner product with respect to L2 norm \int_0^1{p(t)q(t) dt}:
 * use individually computed 2d monomial basis inner products array. */
template <typename F, typename R>
F
PowerPolynomial<F, R>::operator*(PowerPolynomial<F, R> const &q) const
{
    /* check if degree n matches */
    this->checkDegree(q, std::string("Polynomial<F, R>::operator*(Polynomial<F, R> const &q) (inner / scalar product)"));

    /* p = const reference to (this) polynomial */
    const Polynomial<F, R> &p = (*this);
    uint32_t                i, j, n = p.getDegree(); // == q.getDegree()

    F                       res = 0;
    for (i = 0; i < n + 1; i++) {
        for (j = 0; j < n + 1; j++) {
            res += p(i)*q(j)*PowerPolynomial::power_basis_inner_products( {i, j} );
        }
    }
    return res;
}

template <typename F, typename R>
PowerPolynomial<F, R>
PowerPolynomial<F, R>::multiply(PowerPolynomial<F, R> const &q) const
{
    int k, i, m, n;

    PowerPolynomial<F, R> const &p = (*this);

    /* p is in BB(m), q is in BB(n) */
    n   = p.getDegree();
    m   = q.getDegree();

    /* result polynomial is in BB(m + n) */
    PowerPolynomial res(m + n, 0);

    for (k = 0; k < m + n + 1; k++) {
        res(k) = 0;
        for (i = std::max(0, k-m); i <= std::min(k, n); i++) {
            res(k) += p(i)*q(k-i);
        }
    }

    return res;
}

template <typename F, typename R>
PowerPolynomial<F, R>
PowerPolynomial<F, R>::square() const
{
    return (this->multiply( (*this) ));
}

template <typename F, typename R>
PowerPolynomial<F, R>
PowerPolynomial<F, R>::getDerivative() const
{
    uint32_t const      n = this->getDegree();
    uint32_t            k;

    /* only compute if degree of BB representation is > 0 */
    if (n > 0) {
        PowerPolynomial     d(n - 1, 0);

        for (k = 0; k < n; k++) {
            d(k) = (F)(k+1) * this->coeff[k + 1];
        }
        return d;
    }
    /* BB representation degree is zero, return zero polynomial */
    else {
        return (PowerPolynomial<F, R>(0, 0));
    }
}


/* elevate degree: pad coefficients with trailing zeros */
template <typename F, typename R>
void
PowerPolynomial<F, R>::elevateDegree(uint32_t r)
{
    this->degree += r;
    this->coeff.insert(this->coeff.end(), r, 0);

    /* update monomial basis inner products if necessary, degree has been changed manually outside
     * constructor. */
    uint32_t const n = this->getDegree();
    if (PowerPolynomial<F, R>::power_basis_inner_products_max_dim < n) {
        if (PowerPolynomial<F, R>::power_basis_inner_products_mutable) {
            PowerPolynomial<F, R>::initPowerBasisInnerProducts(n);
        }
        else {
            throw("PowerPolynomial: not all necessary inner products not initialized and dynamic inner product recomputation disabled.");
        }
    }
}

template<typename F, typename R>
void
PowerPolynomial<F, R>::matchDegree(
    PowerPolynomial<F, R>  &p,
    PowerPolynomial<F, R>  &q)
{
    uint32_t m, n;
    m = p.getDegree();
    n = q.getDegree();

    if (m > n) {
        q.elevateDegree(m - n);
    }
    else if (m < n) {
        p.elevateDegree(n - m);
    }
}

/*
template <typename F, typename R>
F
PowerPolynomial<F, R>::eval_dn(
    F const    &x,
    uint32_t    n) const
{
    if (this->degree > (n - 1) ) {
        int     i;
        double  dp2_x = this->degree * (this->degree - 1) * this->coeff[this->degree];

        for (i = this->degree - 1; i >= 2; i--) {
            dp2_x = dp2_x*x + i*(i-1)*this->coeff[i];
        }
        return dp2_x;
    }
    else {
        return 0;
    }
}
*/

/*
double
PowerPolynomial<F, R>::divide_monomial(double a)
{
    int     i;
    double  tmp, rem    = coeff[degree];
    coeff[degree]           = 0.0;

    for (i = degree - 1; i >= 0; i--) {
        tmp         = coeff[i];
        coeff[i]    = rem;
        rem         = tmp + rem*a;
    }

    degree--;

    return rem;
}

void
PowerPolynomial<F, R>::multiply_monomial(double a)
{
    int i;

    degree++;
    coeff.resize(degree);

    coeff[degree] = coeff[degree - 1];

    for (i = degree - 1; i >= 1; i--) {
        coeff[i] = coeff[i-1] - coeff[i]*a;
    }
    coeff[0] *= (-a);
}
*/

#if 0
double
PowerPolynomial<F, R>::root_ub()
{
    double  cmax        = fabs(coeff[0]);
    double  csum        = cmax;
    double  tmp, ub1, ub2;

    for (uint32_t i = 1; i <= degree; i++) {
        tmp     = fabs(coeff[i]);
        /* update max coefficient and index if necessary */
        if (tmp > cmax) {
            cmax        = tmp;
        }

        /* add fabs(coeff[i]) to csum */
        csum   += tmp;
    }

    /* use rouche's theorem for two upper bounds for all roots of p */
    ub1     = 1.0 + cmax / coeff[degree];

    tmp     = csum / coeff[degree];
    if (tmp > 1.0) {
        ub2 = tmp;
    }
    else {
        ub2 = 1.0;
    }

    if (ub1 < ub2) {
        return ub1 + 1E-5;
    }
    else {
        return ub2 + 1E-5;
    }
}
#endif

/* ----------------------------------------------------------------------------------------------------------------- *
 *
 * implementation of univariate bernstein-bezier basis polynomial class, deCasteljau algorithm for evaluation, etc..
 *
 * ----------------------------------------------------------------------------------------------------------------- */
template <typename F, typename R>
void
BernsteinPolynomial<F, R>::initBernsteinBasisInnerProducts(uint32_t max_dim)
{
    /* check if bernstein basis polynomial inner product have already been precomputed for a smaller
     * value than max_dim. if so, free it */
    bool recompute = false;
    if (BernsteinPolynomial<F, R>::bernstein_basis_inner_products_max_dim < max_dim) {
        recompute = true;
    }

    /* if recomputaton is necessary */
    if (recompute) {
        /* calculate allocation dimension, never less than default dimension */
        uint32_t alloc_dim = std::max(max_dim, BERNSTEINPOLY_BASIS_IP_DEFAULT_SIZE);

        debugl(0, "(static) BernsteinPolynomial::initBernsteinBasisInnerProducts(): alloc_dim: %u\n", alloc_dim);

        /* resize 4d inner product tensor */
        BernsteinPolynomial<F, R>::bernstein_basis_inner_products.resize(
                { alloc_dim + 1, alloc_dim + 1, alloc_dim + 1, alloc_dim + 1 }
            );

        /* and fill it */
        debugTabInc();
        for (uint32_t m = 0; m < alloc_dim + 1; m++) {
            for (uint32_t n = 0; n < alloc_dim + 1; n++) {
                for (uint32_t i = 0; i < m + 1; i++) {
                    for (uint32_t j = 0; j < n + 1; j++) {
                        debugl(4, "computing basis polynomial inner product element (%3d, %3d, %3d, %3d).\n", m, n, i, j);
                        BernsteinPolynomial<F, R>::bernstein_basis_inner_products( {m, n, i, j} )
                            = BernsteinPolynomial<F, R>::computeBernsteinBasisInnerProduct(m, n, i, j);
                    }
                }
            }
        }
        debugTabDec();

        /* write new allocated dimension back to caller */
        BernsteinPolynomial<F, R>::bernstein_basis_inner_products_max_dim = alloc_dim;

        debugl(0, "(static) BernsteinPolynomial::initBernsteinBasisInnerProducts(): done..\n");
    }
}

/* implementation of the virtual computation function for the inner product of the basis polynomials
 * with respect to the L2-norm on [0,1]. for the Bernstein-Bezier representation, this depends both
 * on the degree of the representation and the index of the basis function:
 * */
template <typename F, typename R>
F
BernsteinPolynomial<F, R>::computeBernsteinBasisInnerProduct(
        uint32_t m,
        uint32_t n,
        uint32_t i,
        uint32_t j)
{
    using Aux::Numbers::bicof;
    debugl(4, "BernsteinPolynomial()::computeBasisInnerProduct(%d, %d, %d, %d)\n", m, n, i, j);

    /* instead of SFINAE, use this (in my humble opinion much more readable) way of distinguishing
     * between different types. this has the advantage that the method type signature or return type
     * are not polluted by std::enable_if<..> and additional template parameters to make SFINAE
     * applicable.
     *
     * if the field type F is a floating point type (float, double, long double), use the following
     * implementation.  otherwise, throw an exception. further implementations (e.g. for custom
     * rational types) can be added with #elif */
    bool F_is_floating_point = std::is_floating_point<F>::value;
    if (F_is_floating_point) {
        /* compute new inner products up to given maximum degree if necessary */
        return (bicof<F>(m, i) * bicof<F>(n, j)) / ( (F)(m + n + 1) * bicof<F>(m + n, i + j) );
    }
    // else if (F_is_rational) {.. }
    else {
        throw("BernsteinPolynomial<F, R>::computeBasisInnerProduct(): not implemented for given type F.");
    }
}

/* return precomputed inner product value. this method relies on the fact that every BernsteinPolynomial constructor
 * (and degree-modifying method) precomputes bernstien basis polynomial inner products for higher degrees if necessary.
 * */
template <typename F, typename R>
F
BernsteinPolynomial<F, R>::getBernsteinBasisInnerProduct(
    uint32_t m,  
    uint32_t n,  
    uint32_t i,  
    uint32_t j)
{
    uint32_t max_dim = std::max(m, n);
    if (bernstein_basis_inner_products_max_dim < max_dim) {
        if (BernsteinPolynomial<F, R>::bernstein_basis_inner_products_initialized) {
            BernsteinPolynomial<F, R>::initBernsteinBasisInnerProducts(max_dim);
        }
        else {
            throw("BernsteinPolynomial()::getPowerBasisInnerProduct(): basis polynomial indices out of range with dynamic recomputation disabled.");
        }
    }
    return BernsteinPolynomial<F, R>::bernstein_basis_inner_products( {m, n, i, j} );
}

template <typename F, typename R>
void
BernsteinPolynomial<F, R>::setInnerProductDataMutable()
{
    BernsteinPolynomial<F, R>::bernstein_basis_inner_products_initialized = true;
}

template <typename F, typename R>
void
BernsteinPolynomial<F, R>::setInnerProductDataImmutable()
{
    BernsteinPolynomial<F, R>::bernstein_basis_inner_products_initialized = false;
}

template <typename F, typename R>
BernsteinPolynomial<F, R>::BernsteinPolynomial() : Polynomial<F, R>()
{
    if (BernsteinPolynomial<F, R>::bernstein_basis_inner_products_max_dim == 0) {
        if (BernsteinPolynomial<F, R>::bernstein_basis_inner_products_initialized) {
            BernsteinPolynomial<F, R>::initBernsteinBasisInnerProducts();
        }
        else {
            throw("BernsteinPolynomial: not all necessary inner products not initialized and dynamic inner product recomputation disabled.");
        }
    }
}

/* for bernstein-bezier polynomials, filling coeff with x amounts to setting the polynomial to the
 * BB(degree) representation of the constant x */
template <typename F, typename R>
BernsteinPolynomial<F, R>::BernsteinPolynomial(uint32_t degree, F const &x) : Polynomial<F, R>(degree, x)
{
    if (BernsteinPolynomial<F, R>::bernstein_basis_inner_products_max_dim < degree) {
        if (BernsteinPolynomial<F, R>::bernstein_basis_inner_products_initialized) {
            BernsteinPolynomial<F, R>::initBernsteinBasisInnerProducts(degree);
        }
        else {
            throw("BernsteinPolynomial: not all necessary inner products initialized and dynamic inner product recomputation disabled.");
        }
    }
}

template <typename F, typename R>
BernsteinPolynomial<F, R>::BernsteinPolynomial(Vector<F> const &coeff)
    : Polynomial<F, R>(coeff)
{
    uint32_t const n = coeff.size() - 1;
    if (BernsteinPolynomial<F, R>::bernstein_basis_inner_products_max_dim < n) {
        if (BernsteinPolynomial<F, R>::bernstein_basis_inner_products_initialized) {
            BernsteinPolynomial<F, R>::initBernsteinBasisInnerProducts(n);
        }
        else {
            throw("BernsteinPolynomial: not all necessary inner products initialized and dynamic inner product recomputation disabled.");
        }
    }
}

template <typename F, typename R>
BernsteinPolynomial<F, R>::BernsteinPolynomial(BernsteinPolynomial<F, R> const &q)
    : Polynomial<F, R>(q)
{
}

template <typename F, typename R>
BernsteinPolynomial<F, R> &
BernsteinPolynomial<F, R>::operator=(BernsteinPolynomial<F, R> const &q)
{
    /* call template base class Polynomial<F, R> (matching F) assignment operator on q, which can be
     * implicitly cast to suitable base class reference. */
    Polynomial<F, R>::operator=(q);

    /* return (*this) */
    return (*this);
}

template <typename F, typename R>
BernsteinPolynomial<F, R>::~BernsteinPolynomial()
{
}

template <typename F, typename R>
void
BernsteinPolynomial<F, R>::setCoeffs(Vector<F> const &coeff)
{
    Polynomial<F, R>::setCoeffs(coeff);

    /* update bernstein inner products if necessary, degree has been changed manually outside
     * constructor. */
    uint32_t const n = this->getDegree();
    if (BernsteinPolynomial<F, R>::bernstein_basis_inner_products_max_dim < n) {
        if (BernsteinPolynomial<F, R>::bernstein_basis_inner_products_initialized) {
            BernsteinPolynomial<F, R>::initBernsteinBasisInnerProducts(n);
        }
        else {
            throw("BernsteinPolynomial: not all necessary inner products initialized and dynamic inner product recomputation disabled.");
        }
    }
}

template <typename F, typename R>
void
BernsteinPolynomial<F, R>::initConstant(F const &alpha)
{
    this->coeff.assign(this->getDegree() + 1, alpha);
}

/* use de-Casteljau algorithm for bezier curves to evaluate polynomial. x-value is x, y-value is
 * evaluated using the algorithm */
template <typename F, typename R>
F
BernsteinPolynomial<F, R>::eval(R const &x) const
{
    return deCasteljau<F>(this->degree, this->coeff, x);
}

template <typename F, typename R>
F
BernsteinPolynomial<F, R>::eval_d(R const &x) const
{
    if (this->degree > 1) {
        return (this->getDerivative()).eval(x);
    }
    else return 0;
}

template <typename F, typename R>
F
BernsteinPolynomial<F, R>::eval_d2(R const &x) const
{
    if (this->degree > 2) {
        return (this->getDerivative().getDerivative()).eval(x);
    }
    else return 0;
}

/* arithmetic */
template <typename F, typename R>
BernsteinPolynomial<F, R>
BernsteinPolynomial<F, R>::operator+(BernsteinPolynomial<F, R> const &q) const
{
    return ( BernsteinPolynomial(this->coeff) += q );
}

template <typename F, typename R>
BernsteinPolynomial<F, R> &
BernsteinPolynomial<F, R>::operator+=(BernsteinPolynomial<F, R> const &q)
{
    Polynomial<F, R>::operator+=(q);
    return (*this);
}

template <typename F, typename R>
BernsteinPolynomial<F, R>
BernsteinPolynomial<F, R>::operator-(BernsteinPolynomial<F, R> const &q) const
{
    return ( BernsteinPolynomial(this->coeff) -= q );
}

template <typename F, typename R>
BernsteinPolynomial<F, R> &
BernsteinPolynomial<F, R>::operator-=(BernsteinPolynomial<F, R> const &q)
{
    Polynomial<F, R>::operator-=(q);
    return (*this);
}

template <typename F, typename R>
BernsteinPolynomial<F, R>
BernsteinPolynomial<F, R>::operator*(F const &x) const
{
    return ( BernsteinPolynomial(this->coeff) *= x );
}

template <typename F, typename R>
BernsteinPolynomial<F, R> &
BernsteinPolynomial<F, R>::operator*=(F const &x)
{
    Polynomial<F, R>::operator*=(x);
    return (*this);
}

template <typename F, typename R>
BernsteinPolynomial<F, R>
BernsteinPolynomial<F, R>::operator/(F const &x) const
{
    return ( BernsteinPolynomial(this->coeff) /= x );
}

template <typename F, typename R>
BernsteinPolynomial<F, R> &
BernsteinPolynomial<F, R>::operator/=(F const &x)
{
    Polynomial<F, R>::operator/=(x);
    return (*this);
}

/* inner product with respect to L2 norm \int_0^1{p(t)q(t) dt} */
template <typename F, typename R>
F
BernsteinPolynomial<F, R>::operator*(BernsteinPolynomial<F, R> const &q) const
{
    return (Polynomial<F, R>::operator*(q));
}

/* multiplication as polynomials over |R or a suitably defined ring */
template <typename F, typename R>
BernsteinPolynomial<F, R>
BernsteinPolynomial<F, R>::multiply(BernsteinPolynomial<F, R> const &q) const
{
    using Aux::Numbers::bicof;

    int i, k, m, n;

    BernsteinPolynomial<F, R> const &p = (*this);

    /* p is in BB(m), q is in BB(n) */
    n   = p.getDegree();
    m   = q.getDegree();

    /* result polynomial is in BB(m + n) */
    BernsteinPolynomial<F, R> res(m + n, 0);

    for (k = 0; k < m + n + 1; k++) {
        res(k) = 0;
        for (i = std::max(0, k - m); i < std::min(n, k) + 1; i++) {
            res(k) += 
                ( p(i) * q(k - i) * bicof<F>(n, i)*bicof<F>(m, k - i) ) / bicof<F>(m + n, k);
        }
    }

    return res;
}

template <typename F, typename R>
BernsteinPolynomial<F, R>
BernsteinPolynomial<F, R>::square() const
{
    return (this->multiply( (*this) ));
}

/* derivative computation */
template <typename F, typename R>
BernsteinPolynomial<F, R>
BernsteinPolynomial<F, R>::getDerivative() const
{
    uint32_t const      n = this->getDegree();
    uint32_t            k;

    /* only compute if degree of BB representation is > 0 */
    if (n > 0) {
        BernsteinPolynomial d(n - 1, 0);
        F                   n_F = (F)(n);

        for (k = 0; k < n; k++) {
            d(k) = n_F * (this->coeff[k + 1] - this->coeff[k]);
        }
        return d;
    }
    /* BB representation degree is zero, return zero polynomial */
    else {
        return (BernsteinPolynomial<F, R>(0, 0));
    }
}

/* degree elevation by r, effectively like multiplication with 1(t) in BB(r). this is done
 * in-place, though. update static inner products data since new degree (n+r) might be new maximum
 * degree ever seen. */
template <typename F, typename R>
void
BernsteinPolynomial<F, R>::elevateDegree(uint32_t r_arg)
{
    using Aux::Numbers::bicof;

    int i, k, n = this->getDegree(), r = r_arg;
    Vector<F> elev_coeff(n + r + 1, 0);

    /* equivalent: perform "multiplication" with 1(t) with coefficients (1, .., 1) in BB(r) */
    for (k = 0; k < n + r + 1; k++) {
        elev_coeff[k] = 0;
        for (i = std::max(0, k - r); i < std::min(n, k) + 1; i++) {
            elev_coeff[k] += ( this->coeff[i] * bicof<F>(n, i)*bicof<F>(r, k - i) ) / bicof<F>(n + r, k);
        }
    }

    /* commit new values */
    this->setCoeffs(elev_coeff);

    /* update inner products */
    uint32_t const d = this->getDegree();
    if (BernsteinPolynomial<F, R>::bernstein_basis_inner_products_max_dim < d) {
        if (BernsteinPolynomial<F, R>::bernstein_basis_inner_products_initialized) {
            BernsteinPolynomial<F, R>::initBernsteinBasisInnerProducts(d);
        }
        else {
            throw("BernsteinPolynomial: not all necessary inner products initialized and dynamic inner product recomputation disabled.");
        }
    }
}

template<typename F, typename R>
void
BernsteinPolynomial<F, R>::matchDegree(
    BernsteinPolynomial<F, R>  &p,
    BernsteinPolynomial<F, R>  &q)
{
    uint32_t m, n;
    m = p.getDegree();
    n = q.getDegree();

    if (m > n) {
        q.elevateDegree(m - n);
    }
    else if (m < n) {
        p.elevateDegree(n - m);
    }
}

template <typename F, typename R>
void
BernsteinPolynomial<F, R>::split(
        R const                        &t,
        BernsteinPolynomial<F, R>      *pleft,
        BernsteinPolynomial<F, R>      *pright) const
{
    Vector<F> coeff_left, coeff_right;
    deCasteljauSplit<F, R>(this->degree, this->coeff, t, coeff_left, coeff_right);
    
    /* update values for polynomials pleft and pright if present */
    if (pleft)  pleft->setCoeffs(coeff_left);
    if (pright) pright->setCoeffs(coeff_right);
}

template <typename F, typename R>
void
BernsteinPolynomial<F, R>::clipToInterval(
    R const                        &t0,
    R const                        &t1,
    BernsteinPolynomial<F, R>      *pclip) const
{
    /* first, clip at t0 and only save right part, then clip at position (t1 - t0) / (1 - t0), 
     * which is the relative position of t1 with respect to the intermediate interval [t0, 1].
     * save only left part then, and we're done.. */
    if (pclip) {
        this->split(   t0,                      NULL,   pclip);
        pclip->split( (t1 - t0) / (1.0 - t0),   pclip,  NULL);
    }
}

/* ----------------------------------------------------------------------------------------------------------------- *
 *
 *                                       auxiliary functions..                                                        
 *
 * ----------------------------------------------------------------------------------------------------------------- */

/* generate B_k^n(t) in BB(n) => delta_{kn} coefficient vector.. */
template <typename F, typename R>
BernsteinPolynomial<F, R>
BB(uint32_t n, uint32_t k)
{
    return ( BernsteinPolynomial<F, R>(Aux::VecMat::kronecker_vec<F>(n+1, k)) );
}

#endif


/* generate B_k^n(t) in BB(n) => delta_{kn} coefficient vector.. */
template <uint32_t degree, typename F, typename R>
BernsteinPolynomial<degree, F, R>
BB(uint32_t k);

/* deCasteljau's algorithm for bezier curve evaluation and splitting adapted for polynomials, i.e. only the "y"
 * coordinates of the control points are given (the "x'-coordinate equals the given splitpoint t).
 * with Vec2, this would of course work for general bezier curves, but that's unnecessary for its use here. */
template <uint32_t degree, typename F, typename R>
F
deCasteljau(const StaticVector<degree+1, F>& coeff, const R& t)
{
    uint32_t i, k;

    /* initialize Pi(0) to bezier coefficients */
    StaticVector<degree+1, F> Pi = coeff;
    F plast_i;

    for (k = 1; k <= degree; k++)
    {
        /* save old value of Pi[0] as plast_i for first iteration */
        plast_i = Pi[0];

        /* compute new values of Pi, the old "lookahead" value plast_i is always saved */
        for (i = 0; i <= degree - k; i++)
        {
            Pi[i] = (1.0 - t) * plast_i + t*Pi[i+1];

            /* save old value of p[i+1] as pilast_i for next iteration */
            plast_i = Pi[i+1];
        }
    }
    return Pi[0];
}

template <uint32_t degree, typename F, typename R>
void
deCasteljauSplit
(
    const StaticVector<degree+1, F>& coeff,
    const R& t,
    StaticVector<degree+1, F>& coeff_left,
    StaticVector<degree+1, F>& coeff_right)
{
    uint32_t i, k;

    // We need this copy as coeff might be the same reference as coeff_left or coeff_right.
    StaticVector<degree+1, F> Pi = coeff;
    F plast_i;

    coeff_left.assign(0);
    coeff_right.assign(0);

    coeff_left[0] = Pi[0];
    coeff_right[degree] = Pi[degree];

    for (k = 1; k <= degree; k++)
    {
        /* save old value of Pi[0] as plast_i for first iteration */
        plast_i = Pi[0];

        /* compute new values of Pi, the old "lookahead" value plast_i is always saved */
        for (i = 0; i <= degree - k; i++)
        {
            Pi[i] = (1.0 - t) * plast_i + t*Pi[i+1];

            /* save old value of p[i+1] as pilast_i for next iteration */
            plast_i = Pi[i+1];
        }

        /* Pi contains Pi(k). store relevant values for left and right polynomial */
        coeff_left[k] = Pi[0];
        coeff_right[degree - k] = Pi[degree - k];
    }
}



/* ----------------------------------------------------------------------------------------------------------------- *
 *
 *                        implementation of abstract univariate polynomial class....
 *
 * ----------------------------------------------------------------------------------------------------------------- */

template <uint32_t degree, typename F, typename R>
Polynomial<degree, F, R>::Polynomial()
: coeff(0.0)
{}

template <uint32_t degree, typename F, typename R>
Polynomial<degree, F, R>::Polynomial(const F& x)
: coeff(x)
{}

template <uint32_t degree, typename F, typename R>
Polynomial<degree, F, R>::Polynomial(const coeff_type& _coeff)
: coeff(_coeff)
{}

template <uint32_t degree, typename F, typename R>
Polynomial<degree, F, R>::Polynomial(const this_type& q)
: coeff(q.coeff)
{}

template <uint32_t degree, typename F, typename R>
Polynomial<degree, F, R>::~Polynomial()
{}

template <uint32_t degree, typename F, typename R>
Polynomial<degree, F, R>&
Polynomial<degree, F, R>::operator=(const this_type& q)
{
    coeff = q.coeff;
    return *this;
}

template <uint32_t degree, typename F, typename R>
uint32_t
Polynomial<degree, F, R>::getDegree() const
{
    return degree;
}

template <uint32_t degree, typename F, typename R>
const StaticVector<degree+1, F>&
Polynomial<degree, F, R>::getCoeffs() const
{
    return coeff;
}

template <uint32_t degree, typename F, typename R>
StaticVector<degree+1, F>&
Polynomial<degree, F, R>::getCoeffs()
{
    return coeff;
}

template <uint32_t degree, typename F, typename R>
void
Polynomial<degree, F, R>::setCoeffs(const coeff_type& _coeff)
{
    coeff = _coeff;
}

template <uint32_t degree, typename F, typename R>
F
Polynomial<degree, F, R>::getMaxAbsCoeff() const
{
    F cmax = std::abs(this->coeff[0]);
    F abs_i;
    for (uint32_t i = 1; i < degree + 1; i++)
    {
        abs_i = fabs(this->coeff[i]);
        if (abs_i > cmax)
            cmax = abs_i;
    }
    return cmax;
}

template <uint32_t degree, typename F, typename R>
void
Polynomial<degree, F, R>::zero()
{
    initConstant(0);
}

template <uint32_t degree, typename F, typename R>
F&
Polynomial<degree, F, R>::operator()(uint32_t cidx)
{
    return coeff[cidx];
}

template <uint32_t degree, typename F, typename R>
F
Polynomial<degree, F, R>::operator()(uint32_t cidx) const
{
    return coeff[cidx];
}

template <uint32_t degree, typename F, typename R>
F &
Polynomial<degree, F, R>::operator[](uint32_t cidx)
{
    return coeff[cidx];
}

template <uint32_t degree, typename F, typename R>
F
Polynomial<degree, F, R>::operator[](uint32_t cidx) const
{
    return coeff[cidx];
}

template <uint32_t degree, typename F, typename R>
Polynomial<degree, F, R> &
Polynomial<degree, F, R>::operator+=(const this_type& q)
{
    for (uint32_t i = 0; i < degree+1; i++)
        coeff[i] += q[i];

    return *this;
}

template <uint32_t degree, typename F, typename R>
Polynomial<degree, F, R> &
Polynomial<degree, F, R>::operator-=(const this_type& q)
{
    for (uint32_t i = 0; i < degree+1; i++)
        coeff[i] -= q[i];

    return *this;
}

template <uint32_t degree, typename F, typename R>
Polynomial<degree, F, R> &
Polynomial<degree, F, R>::operator*=(F const &x)
{
    for (uint32_t i = 0; i < degree+1; i++)
        coeff[i] *= x;

    return *this;
}

template <uint32_t degree, typename F, typename R>
Polynomial<degree, F, R> &
Polynomial<degree, F, R>::operator/=(F const &x)
{
    for (uint32_t i = 0; i < degree+1; i++)
        coeff[i] /= x;

    return *this;
}


// TODO: This method clearly cannot work, since the pure virtual getBasisInnerProduct
//       it calls has never been declared, let alone defined.
#if 0
template <uint32_t degree, typename F, typename R>
F
Polynomial<degree, F, R>::operator*(const this_type& q) const
{
    const Polynomial<degree, F, R>& p = *this;
    uint32_t i, j;
    F res = 0;

    for (i = 0; i < degree+1; i++) {
        for (j = 0; j < degree+1; j++) {
            /* use virtual method getBasisInnerProduct, which is pure-virtual in Polynomial<F, R>
             * and must be implemented by all derived classes. this enables the reuse of this method
             * from derived classes through late-binding. due to type safety however, every derived
             * class needs to reimplement operator* as a wrapper, which calls this method as
             * Polynomial<degree, F, R>::operator*(..) (or implement a completely custom one obviously). */
            res += p(i)*q(j)*this->getBasisInnerProduct(i, j);
        }
    }
    return res;
}
#endif

template <uint32_t degree, typename F, typename R>
void
Polynomial<degree, F, R>::printCoeff() const
{
    PrintCoeffImpl<F, R>(*this);
}

/* print coefficients: specializations for R = F = {float, double} */
template <uint32_t degree, typename F, typename R>
template <typename dummy>
Polynomial<degree, F, R>::PrintCoeffImpl<float, float, dummy>::PrintCoeffImpl
(const Polynomial<degree, F, R>& p)
{
    uint32_t        i;
    static char     s[1024];
    static char     tmp[256];

    sprintf(s, "[ ");
    for (i = 0; i < degree + 1; i++) {
        sprintf(tmp, "%+20.13E ", p.coeff[i]);
        strncat(s, tmp, 256);
    }
    strncat(s, "]\n", 2);
    debugl(0, "%s", s);
}

template <uint32_t degree, typename F, typename R>
template <typename dummy>
Polynomial<degree, F, R>::PrintCoeffImpl<double, double, dummy>::PrintCoeffImpl
(const Polynomial<degree, F, R>& p)
{
    uint32_t        i;
    static char     s[1024];
    static char     tmp[256];

    sprintf(s, "[ ");
    for (i = 0; i < degree + 1; i++) {
        sprintf(tmp, "%+20.13E ", p.coeff[i]);
        strncat(s, tmp, 256);
    }
    strncat(s, "]\n", 2);
    debugl(0, "%s", s);
    printf("%s", s);
}


template <uint32_t degree, typename F, typename R>
void
Polynomial<degree, F, R>::writePlotFile(R t0, R t1, uint32_t ticks, const std::string& filename) const
{
    WritePlotFileImpl<F, R>(t0, t1, ticks, filename, *this);
}

template <uint32_t degree, typename F, typename R>
template <typename dummy>
Polynomial<degree, F, R>::WritePlotFileImpl<float, float, dummy>::WritePlotFileImpl
(R t0, R t1, uint32_t ticks, const std::string& filename, const Polynomial<degree, F, R>& p)
{
    FILE *fout = fopen(filename.c_str(), "w");
    if (!fout) {
        fprintf(stderr, "Polynomial::writePlotFile(): can't open file \"%s\" for writing.\n", filename.c_str() );
        return;;
    }

    for (uint32_t i = 0; i <= ticks; i++) {
        float t = t0 + ((float)i / (float)ticks)*(t1 - t0);
        fprintf(fout, "%12.5E %12.5E\n", t, p.eval(t));
    }
    fclose(fout);
}

template <uint32_t degree, typename F, typename R>
template <typename dummy>
Polynomial<degree, F, R>::WritePlotFileImpl<double, double, dummy>::WritePlotFileImpl
(R t0, R t1, uint32_t ticks, const std::string& filename, const Polynomial<degree, F, R>& p)
{
    FILE *fout = fopen(filename.c_str(), "w");
    if (!fout) {
        fprintf(stderr, "Polynomial::writePlotFile(): can't open file \"%s\" for writing.\n", filename.c_str() );
        return;;
    }

    uint32_t    i;
    double      t;

    for (i = 0; i <= ticks; i++) {
        t = t0 + ((double)i / (double)ticks)*(t1 - t0);
        fprintf(fout, "%12.5E %12.5E\n", t, p.eval(t));
    }
    fclose(fout);
}



/* ----------------------------------------------------------------------------------------------------------------- *
 *
 * implementation of univariate power basis (aka monomial basis) polynomial class, horner scheme for evaluation, etc.
 *
 * ----------------------------------------------------------------------------------------------------------------- */

/* the power basis inner products do not depend on the degree of the representation, but only on the degrees of the
 * monomials: x^i * x^j = \int_0^1{x^{i+j}dx} = 1 / (i + j + 1)
 *
 * => use only 2d tensor (array) instead of 4d array with custom init method.
 *
 * note that for BernsteinPolynomial, the inner products do depend on the degree of the representation and hence the 4d
 * array is really necessary. for legendre polynomials on the other hand, even a 1d array will suffice due to
 * orthogonality. */
template <uint32_t degree, typename F, typename R>
void
PowerPolynomial<degree, F, R>::initPowerBasisInnerProducts()
{
    static bool recompute = true;

    if (recompute)
    {
        debugl(1, "(static) PowerPolynomial::initPowerBasisInnerProducts(): degree %u\n", degree);
        debugTabInc();

        for (uint32_t i = 0; i < degree + 1; ++i)
            for (uint32_t j = 0; j < degree + 1; ++j)
                power_basis_inner_products(i, j) = computePowerBasisInnerProduct(i, j);

        recompute = false;

        debugTabDec();
        debugl(1, "(static) PowerPolynomial::initPowerBasisInnerProducts(): done..\n");
    }
}

/* indices m and n for representation degree are irrelevant: inner products depends only on i and j in
 * integral over x^i*x^j = x^{i + j}. */
template <uint32_t degree, typename F, typename R>
F
PowerPolynomial<degree, F, R>::computePowerBasisInnerProduct(uint32_t i, uint32_t j)
{
    debugl(4, "PowerPolynomial()::computePowerBasisInnerProduct(%d, %d)\n", i, j);
    return (F)1 / (F)(i + j + 1);
}

/* return pre-computed inner product value from static PowerPolynomial<degree, F, R>::power_basis_inner_products. if value has
 * not been precomputed yet, do it now if dynamic recomputation is enabled */
template <uint32_t degree, typename F, typename R>
F
PowerPolynomial<degree, F, R>::getPowerBasisInnerProduct(uint32_t i, uint32_t j)
{
    if (i > degree || j > degree)
    {
        std::ostringstream oss;
        oss << "PowerPolynomial<" << degree << ">::getPowerBasisInnerProduct(): "
               "Requested illegal product (" << i << ", " << j << ").";
        throw(oss.str().c_str());
    }

    // initialize if necessary
    initPowerBasisInnerProducts();

    return power_basis_inner_products(i, j);
}



template <uint32_t degree, typename F, typename R>
PowerPolynomial<degree, F, R>::PowerPolynomial()
{}

template <uint32_t degree, typename F, typename R>
PowerPolynomial<degree, F, R>::PowerPolynomial(const F& x)
: Polynomial<degree, F, R>(0)
{}

template <uint32_t degree, typename F, typename R>
PowerPolynomial<degree, F, R>::PowerPolynomial(const coeff_type& coeff)
: Polynomial<degree, F, R>(coeff)
{}

template <uint32_t degree, typename F, typename R>
PowerPolynomial<degree, F, R>::PowerPolynomial(const this_type& q)
: Polynomial<degree, F, R>(q)
{}

template <uint32_t degree, typename F, typename R>
PowerPolynomial<degree, F, R>::~PowerPolynomial()
{}


template <uint32_t degree, typename F, typename R>
PowerPolynomial<degree, F, R>&
PowerPolynomial<degree, F, R>::operator=(const this_type& q)
{
    Polynomial<degree, F, R>::operator=(q);
    return (*this);
}

template <uint32_t degree, typename F, typename R>
void
PowerPolynomial<degree, F, R>::initConstant(const F& x)
{
    coeff.assign(0);
    coeff[0] = x;
}

template <uint32_t degree, typename F, typename R>
F
PowerPolynomial<degree, F, R>::eval(const R& x) const
{
    int i;
    F p_x = coeff[degree];

    for (i = degree-1; i >= 0; i--)
        p_x = p_x*x + coeff[i];

    return p_x;
}

template <uint32_t degree, typename F, typename R>
F
PowerPolynomial<degree, F, R>::eval_d(const R& x) const
{
    if (degree > 0)
    {
        int i;
        F dp_x = coeff[degree] * (R)degree;

        for (i = degree-1; i >= 1; i--)
            dp_x = dp_x*x + coeff[i] * (R)i;

        return dp_x;
    }
    else return 0;
}

template <uint32_t degree, typename F, typename R>
F
PowerPolynomial<degree, F, R>::eval_d2(const R& x) const
{
    if (degree > 1)
    {
        int i;
        F   dp2_x = (R)(degree*(degree-1)) * coeff[degree];

        for (i = degree - 1; i >= 2; i--)
            dp2_x = dp2_x*x + (R)(i*(i-1)) * coeff[i];

        return dp2_x;
    }
    else return 0;
}


template <uint32_t degree, typename F, typename R>
PowerPolynomial<degree, F, R>
PowerPolynomial<degree, F, R>::operator+(const this_type& q) const
{
    return PowerPolynomial(coeff) += q;
}

template <uint32_t degree, typename F, typename R>
PowerPolynomial<degree, F, R>&
PowerPolynomial<degree, F, R>::operator+=(const this_type& q)
{
    Polynomial<degree, F, R>::operator+=(q);
    return *this;
}

template <uint32_t degree, typename F, typename R>
PowerPolynomial<degree, F, R>
PowerPolynomial<degree, F, R>::operator-(const this_type& q) const
{
    return PowerPolynomial(coeff) -= q;
}

template <uint32_t degree, typename F, typename R>
PowerPolynomial<degree, F, R>&
PowerPolynomial<degree, F, R>::operator-=(const this_type& q)
{
    Polynomial<degree, F, R>::operator-=(q);
    return *this;
}

template <uint32_t degree, typename F, typename R>
PowerPolynomial<degree, F, R>
PowerPolynomial<degree, F, R>::operator*(const F& x) const
{
    return PowerPolynomial(coeff) *= x;
}

template <uint32_t degree, typename F, typename R>
PowerPolynomial<degree, F, R>&
PowerPolynomial<degree, F, R>::operator*=(const F& x)
{
    Polynomial<degree, F, R>::operator*=(x);
    return *this;
}

template <uint32_t degree, typename F, typename R>
PowerPolynomial<degree, F, R>
PowerPolynomial<degree, F, R>::operator/(const F& x) const
{
    return PowerPolynomial(coeff) /= x;
}

template <uint32_t degree, typename F, typename R>
PowerPolynomial<degree, F, R>&
PowerPolynomial<degree, F, R>::operator/=(const F& x)
{
    Polynomial<degree, F, R>::operator/=(x);
    return *this;
}

/* inner product with respect to L2 norm \int_0^1{p(t)q(t) dt}:
 * use individually computed 2d monomial basis inner products array. */
template <uint32_t degree, typename F, typename R>
F
PowerPolynomial<degree, F, R>::operator*(const this_type& q) const
{
    // p = const reference to (this) polynomial
    const Polynomial<degree, F, R>& p = *this;
    uint32_t i, j;

    F res = 0;
    for (i = 0; i < degree+1; i++)
        for (j = 0; j < degree+1; j++)
            res += p(i)*q(j)*PowerPolynomial::power_basis_inner_products( {i, j} );

    return res;
}

template <uint32_t degree, typename F, typename R>
template <uint32_t deg>
PowerPolynomial<degree+deg, F, R>
PowerPolynomial<degree, F, R>::multiply(const PowerPolynomial<deg, F, R>& q) const
{
    int k, i;
    const PowerPolynomial<degree, F, R>& p = *this;

    PowerPolynomial<degree+deg, F, R> res(0);
    for (k = 0; k < deg + degree + 1; ++k)
    {
        res(k) = 0;
        for (i = std::max(0, k-(int)deg); i <= std::min(k, (int)degree); ++i)
            res(k) += p(i)*q(k-i);
    }

    return res;
}

template <uint32_t degree, typename F, typename R>
PowerPolynomial<2*degree, F, R>
PowerPolynomial<degree, F, R>::square() const
{
    return multiply(*this);
}

template <uint32_t degree, typename F, typename R>
PowerPolynomial<(degree>0) ? degree-1 : 0, F, R>
PowerPolynomial<degree, F, R>::getDerivative() const
{
    uint32_t k;

    PowerPolynomial<(degree>0) ? degree-1 : 0, F, R> d(0);
    for (k = 0; k < degree; k++)
        d(k) = (F)(k+1) * coeff[k+1];

    return d;
}


/* elevate degree: pad coefficients with trailing zeros */
template <uint32_t degree, typename F, typename R>
template <uint32_t deg>
PowerPolynomial<degree+deg, F, R>
PowerPolynomial<degree, F, R>::elevateDegree()
{
    PowerPolynomial<degree+deg, F, R> p(0);
    for (uint32_t i = 0; i < deg+1; ++i)
        p(i) = coeff[i];
    return p;
}



/* ----------------------------------------------------------------------------------------------------------------- *
 *
 * implementation of univariate bernstein-bezier basis polynomial class, deCasteljau algorithm for evaluation, etc..
 *
 * ----------------------------------------------------------------------------------------------------------------- */
template <uint32_t degree, typename F, typename R>
void
BernsteinPolynomial<degree, F, R>::initBernsteinBasisInnerProducts()
{
    static bool recompute = true;

    // initialize if necessary
    if (recompute)
    {
        debugl(1, "(static) BernsteinPolynomial::initBernsteinBasisInnerProducts(): degree %u\n", degree);
        debugTabInc();

        for (uint32_t i = 0; i < degree + 1; ++i)
            for (uint32_t j = 0; j < degree + 1; ++j)
                bernstein_basis_inner_products(i, j) = computeBernsteinBasisInnerProduct(i, j);

        recompute = false;

        debugTabDec();
        debugl(0, "(static) BernsteinPolynomial::initBernsteinBasisInnerProducts(): done..\n");
    }
}

/* implementation of the virtual computation function for the inner product of the basis polynomials
 * with respect to the L2-norm on [0,1]. for the Bernstein-Bezier representation, this depends both
 * on the degree of the representation and the index of the basis function:
 * */
template <uint32_t degree, typename F, typename R>
F
BernsteinPolynomial<degree, F, R>::computeBernsteinBasisInnerProduct(uint32_t i, uint32_t j)
{
    using Aux::Numbers::bicof;
    debugl(4, "BernsteinPolynomial()::computeBasisInnerProduct(%d, %d, %d, %d)\n", degree, degree, i, j);

    /* instead of SFINAE, use this (in my humble opinion much more readable) way of distinguishing
     * between different types. this has the advantage that the method type signature or return type
     * are not polluted by std::enable_if<..> and additional template parameters to make SFINAE
     * applicable.
     *
     * if the field type F is a floating point type (float, double, long double), use the following
     * implementation.  otherwise, throw an exception. further implementations (e.g. for custom
     * rational types) can be added with #elif */
    bool F_is_floating_point = std::is_floating_point<F>::value;
    if (F_is_floating_point) {
        /* compute new inner products up to given maximum degree if necessary */
        return bicof<F>(degree, i) * bicof<F>(degree, j)
                / ((F)(degree + degree + 1) * bicof<F>(degree + degree, i + j));
    }
    // else if (F_is_rational) {.. }
    else {
        throw("BernsteinPolynomial<degree, F, R>::computeBasisInnerProduct(): not implemented for given type F.");
    }
}

/* return precomputed inner product value. this method relies on the fact that every BernsteinPolynomial constructor
 * (and degree-modifying method) precomputes bernstein basis polynomial inner products for higher degrees if necessary.
 * */
template <uint32_t degree, typename F, typename R>
F
BernsteinPolynomial<degree, F, R>::getBernsteinBasisInnerProduct(uint32_t i, uint32_t j)
{
    if (i > degree || j > degree)
    {
        std::ostringstream oss;
        oss << "PowerPolynomial<" << degree << ">::getPowerBasisInnerProduct(): "
               "Requested illegal product (" << i << ", " << j << ").";
        throw(oss.str().c_str());
    }

    // initialize if necessary
    initBernsteinBasisInnerProducts();

    return bernstein_basis_inner_products(i, j);
}

template <uint32_t degree, typename F, typename R>
BernsteinPolynomial<degree, F, R>::BernsteinPolynomial()
{}

/* for bernstein-bezier polynomials, filling coeff with x amounts to setting the polynomial to the
 * BB(degree) representation of the constant x */
template <uint32_t degree, typename F, typename R>
BernsteinPolynomial<degree, F, R>::BernsteinPolynomial(F const &x)
: base_type(x)
{}

template <uint32_t degree, typename F, typename R>
BernsteinPolynomial<degree, F, R>::BernsteinPolynomial(const coeff_type& coeff)
: base_type(coeff)
{}

template <uint32_t degree, typename F, typename R>
BernsteinPolynomial<degree, F, R>::BernsteinPolynomial(const this_type& q)
: base_type(q)
{}

template <uint32_t degree, typename F, typename R>
BernsteinPolynomial<degree, F, R>::~BernsteinPolynomial()
{}

template <uint32_t degree, typename F, typename R>
BernsteinPolynomial<degree, F, R>&
BernsteinPolynomial<degree, F, R>::operator=(const this_type& q)
{
    /* call template base class Polynomial<F, R> (matching F) assignment operator on q, which can be
     * implicitly cast to suitable base class reference. */
    base_type::operator=(q);

    /* return (*this) */
    return *this;
}

template <uint32_t degree, typename F, typename R>
void
BernsteinPolynomial<degree, F, R>::setCoeffs(const coeff_type& coeff)
{
    base_type::setCoeffs(coeff);
}

template <uint32_t degree, typename F, typename R>
void
BernsteinPolynomial<degree, F, R>::initConstant(const F& x)
{
    coeff.assign(x);
}

/* use de-Casteljau algorithm for bezier curves to evaluate polynomial. x-value is x, y-value is
 * evaluated using the algorithm */
template <uint32_t degree, typename F, typename R>
F
BernsteinPolynomial<degree, F, R>::eval(const R& x) const
{
    return deCasteljau<degree, F>(coeff, x);
}

template <uint32_t degree, typename F, typename R>
F
BernsteinPolynomial<degree, F, R>::eval_d(const R& x) const
{
    if (degree > 1)
        return getDerivative().eval(x);
    else return 0;
}

template <uint32_t degree, typename F, typename R>
F
BernsteinPolynomial<degree, F, R>::eval_d2(const R& x) const
{
    if (degree > 2) {
        return getDerivative().getDerivative().eval(x);
    }
    else return 0;
}


template <uint32_t degree, typename F, typename R>
BernsteinPolynomial<degree, F, R>
BernsteinPolynomial<degree, F, R>::operator+(const this_type& q) const
{
    return this_type(coeff) += q;
}

template <uint32_t degree, typename F, typename R>
BernsteinPolynomial<degree, F, R> &
BernsteinPolynomial<degree, F, R>::operator+=(const this_type& q)
{
    base_type::operator+=(q);
    return (*this);
}

template <uint32_t degree, typename F, typename R>
BernsteinPolynomial<degree, F, R>
BernsteinPolynomial<degree, F, R>::operator-(const this_type& q) const
{
    return this_type(coeff) -= q;
}

template <uint32_t degree, typename F, typename R>
BernsteinPolynomial<degree, F, R> &
BernsteinPolynomial<degree, F, R>::operator-=(const this_type& q)
{
    base_type::operator-=(q);
    return *this;
}

template <uint32_t degree, typename F, typename R>
BernsteinPolynomial<degree, F, R>
BernsteinPolynomial<degree, F, R>::operator*(const F& x) const
{
    return this_type(coeff) *= x;
}

template <uint32_t degree, typename F, typename R>
BernsteinPolynomial<degree, F, R> &
BernsteinPolynomial<degree, F, R>::operator*=(const F& x)
{
    base_type::operator*=(x);
    return *this;
}

template <uint32_t degree, typename F, typename R>
BernsteinPolynomial<degree, F, R>
BernsteinPolynomial<degree, F, R>::operator/(const F& x) const
{
    return this_type(coeff) /= x;
}

template <uint32_t degree, typename F, typename R>
BernsteinPolynomial<degree, F, R> &
BernsteinPolynomial<degree, F, R>::operator/=(const F& x)
{
    base_type::operator/=(x);
    return *this;
}

/* inner product with respect to L2 norm \int_0^1{p(t)q(t) dt} */
template <uint32_t degree, typename F, typename R>
F
BernsteinPolynomial<degree, F, R>::operator*(const this_type& q) const
{
    return (base_type::operator*(q));
}

/* multiplication as polynomials over |R or a suitably defined ring */
template <uint32_t degree, typename F, typename R>
template <uint32_t deg>
BernsteinPolynomial<degree+deg, F, R>
BernsteinPolynomial<degree, F, R>::multiply(const BernsteinPolynomial<deg, F, R>& q) const
{
    using Aux::Numbers::bicof;

    int i, k;
    const this_type& p = *this;

    BernsteinPolynomial<degree+deg, F, R> res(0);

    for (k = 0; k < (int)deg + (int)degree + 1; k++)
    {
        res(k) = 0;
        for (i = std::max(0, k - (int)deg); i < std::min((int)degree, k) + 1; i++)
            res(k) += p(i) * q(k-i) * bicof<F>(degree, i) * bicof<F>(deg, k-i) / bicof<F>(deg+degree, k);
    }

    return res;
}

template <uint32_t degree, typename F, typename R>
BernsteinPolynomial<2*degree, F, R>
BernsteinPolynomial<degree, F, R>::square() const
{
    return multiply(*this);
}

/* derivative computation */
template <uint32_t degree, typename F, typename R>
BernsteinPolynomial<(degree>0) ? degree-1 : 0, F, R>
BernsteinPolynomial<degree, F, R>::getDerivative() const
{
    uint32_t k;

    BernsteinPolynomial<(degree>0) ? degree-1 : 0, F, R> d(0);
    F n_F(degree);
    for (k = 0; k < degree; k++)
        d(k) = n_F * (coeff[k+1] - coeff[k]);
    return d;
}

/* degree elevation by r, effectively like multiplication with 1(t) in BB(r). this is done
 * in-place, though. update static inner products data since new degree (n+r) might be new maximum
 * degree ever seen. */
template <uint32_t degree, typename F, typename R>
template <uint32_t deg>
BernsteinPolynomial<deg, F, R>
BernsteinPolynomial<degree, F, R>::elevateDegree()
{
    using Aux::Numbers::bicof;

    int i, k;
    BernsteinPolynomial<deg, F, R> res(0);
    StaticVector<deg+1, F> elev_coeff;

    /* equivalent: perform "multiplication" with 1(t) with coefficients (1, .., 1) in BB(r) */
    for (k = 0; k < (int)deg + 1; k++)
    {
        elev_coeff[k] = 0;
        for (i = std::max(0, k - (int)deg + (int)degree); i < std::min((int)degree, k) + 1; i++)
            elev_coeff[k] += coeff[i] * bicof<F>(degree, i) * bicof<F>(deg-degree, k-i) / bicof<F>(deg, k);
    }

    res.setCoeffs(elev_coeff);
    return res;
}


template <uint32_t degree, typename F, typename R>
void
BernsteinPolynomial<degree, F, R>::split(const R& t, this_type* pleft, this_type* pright) const
{
    coeff_type coeff_left, coeff_right;
    deCasteljauSplit<degree, F, R>(coeff, t, coeff_left, coeff_right);

    /* update values for polynomials pleft and pright if present */
    if (pleft)  pleft->setCoeffs(coeff_left);
    if (pright) pright->setCoeffs(coeff_right);
}

template <uint32_t degree, typename F, typename R>
void
BernsteinPolynomial<degree, F, R>::clipToInterval(const R& t0, const R& t1, this_type* pclip) const
{
    /* first, clip at t0 and only save right part, then clip at position (t1 - t0) / (1 - t0),
     * which is the relative position of t1 with respect to the intermediate interval [t0, 1].
     * save only left part then, and we're done.. */
    if (pclip)
    {
        split(t0, NULL, pclip);
        pclip->split((t1 - t0) / (1.0 - t0), pclip, NULL);
    }
}

/* ----------------------------------------------------------------------------------------------------------------- *
 *
 *                                       auxiliary functions..
 *
 * ----------------------------------------------------------------------------------------------------------------- */

/* generate B_k^n(t) in BB(n) => delta_{kn} coefficient vector.. */
template <uint32_t degree, typename F, typename R>
BernsteinPolynomial<degree, F, R>
BB(uint32_t k)
{
    return BernsteinPolynomial<degree, F, R>(Aux::VecMat::kronecker_static_vec<degree+1, F>(k));
}
