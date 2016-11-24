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

/* ----------------------------------------------------------------------------------------------------------------- *
 *
 *                        implementation of abstract bivariate polynomial class.....                                  
 *
 * ----------------------------------------------------------------------------------------------------------------- */

template <typename F, typename R>
BivariatePolynomial<F, R>::BivariatePolynomial()
    : m(0), n(0)
{
    this->coeff.resize(1, 1);
    this->coeff.fill(0);
}

template <typename F, typename R>
BivariatePolynomial<F, R>::BivariatePolynomial(
    uint32_t    m,
    uint32_t    n,
    F const    &x)
{
    this->m = m;
    this->n = n;
    this->coeff.resize(m + 1, n + 1);
    this->coeff.fill(x);
}

template <typename F, typename R>
BivariatePolynomial<F, R>::BivariatePolynomial(BivariatePolynomial<F, R> const &q)
{
    this->m         = q.m;
    this->n         = q.n;
    this->coeff     = q.coeff;
}


template <typename F, typename R>
BivariatePolynomial<F, R> &
BivariatePolynomial<F, R>::operator=(BivariatePolynomial const &q)
{
    this->m     = q.m;
    this->n     = q.n;
    this->coeff = q.coeff;

    return (*this);
}

template <typename F, typename R>
BivariatePolynomial<F, R>::~BivariatePolynomial()
{
}

template <typename F, typename R>
void
BivariatePolynomial<F, R>::getDegree(uint32_t &m, uint32_t &n) const
{
    m = this->m;
    n = this->n;
}


template <typename F, typename R>
F &
BivariatePolynomial<F, R>::operator()(uint32_t i, uint32_t j)
{
    return this->coeff(i, j);
}

template <typename F, typename R>
F
BivariatePolynomial<F, R>::operator()(uint32_t i, uint32_t j) const
{
    return this->coeff(i, j);
}

template <typename F, typename R>
Matrix<F>
BivariatePolynomial<F, R>::getCoeffs() const
{
    return this->coeff;
}


template <typename F, typename R>
void
BivariatePolynomial<F, R>::setCoeffs(Matrix<F> const &coeff)
{
    uint32_t coeff_m, coeff_n;
    coeff_m = coeff.numRows();
    coeff_n = coeff.numCols();

    if (coeff_m > 0 && coeff_n > 0) {
        this->coeff   = coeff;
        this->m       = coeff_m - 1;
        this->n       = coeff_n - 1;
    }
    else {
        throw("BivariatePolynomial::setCoeffs(): (m x n) coefficient matrix with m = 0 || n = 0 given.");
    }
}

template <typename F, typename R>
F
BivariatePolynomial<F, R>::getMaxAbsCoeff() const
{
    uint32_t    i, j;
    F           cmax = std::abs(this->coeff(0, 0));
    F           abs_ij;
    for (i = 0; i < m + 1; i++) {
        for (j = 0; j < n + 1; j++) {
            abs_ij = std::abs(this->coeff(i, j));
            if (abs_ij > cmax) {
                cmax = abs_ij;
            }
        }
    }
    return cmax;
}

template <typename F, typename R>
void
BivariatePolynomial<F, R>::initConstant(F const &x)
{
    this->coeff.fill(x);
}

template <typename F, typename R>
void
BivariatePolynomial<F, R>::zero()
{
    this->coeff.fill(0);
}

/* arithmetic */
template <typename F, typename R>
BivariatePolynomial<F, R> &
BivariatePolynomial<F, R>::operator+=(BivariatePolynomial<F, R> const &q)
{
    this->coeff += q.coeff;
    return (*this);
}

template <typename F, typename R>
BivariatePolynomial<F, R> &
BivariatePolynomial<F, R>::operator-=(BivariatePolynomial<F, R> const &q)
{
    this->coeff -= q.coeff;
    return (*this);
}

template <typename F, typename R>
BivariatePolynomial<F, R> &
BivariatePolynomial<F, R>::operator*=(F const &x)
{
    this->coeff *= x;
    return (*this);
}

template <typename F, typename R>
BivariatePolynomial<F, R> &
BivariatePolynomial<F, R>::operator/=(F const &x)
{
    this->coeff /= x;
    return (*this);
}


/* ----------------------------------------------------------------------------------------------------------------- *
 *
 *                        implementation of bivariate bernstein polynomial class....                                  
 *
 * ----------------------------------------------------------------------------------------------------------------- */

template <typename F, typename R>
BiBernsteinPolynomial<F, R>::BiBernsteinPolynomial() : BivariatePolynomial<F, R>()
{
}


template <typename F, typename R>
BiBernsteinPolynomial<F, R>::BiBernsteinPolynomial(
    uint32_t    m,
    uint32_t    n,
    F const    &x) : BivariatePolynomial<F, R>(m, n, x)
{
}


template <typename F, typename R>
BiBernsteinPolynomial<F, R>::BiBernsteinPolynomial(Matrix<F> const &coeff_arg)
{
    uint32_t coeff_m, coeff_n;
    coeff_m = coeff_arg.numRows();
    coeff_n = coeff_arg.numCols();

    if (coeff_m > 0 && coeff_n > 0) {
        this->coeff = coeff_arg;
        this->m     = coeff_m - 1;
        this->n     = coeff_n - 1;
    }
    else {
        throw("BivariatePolynomial::setCoeffs(): (m x n) coefficient matrix with m = 0 || n = 0 given.");
    }
}

template <typename F, typename R>
BiBernsteinPolynomial<F, R>::BiBernsteinPolynomial(BiBernsteinPolynomial const &q) : 
    BivariatePolynomial<F, R>(q)
{
}

template <typename F, typename R>
BiBernsteinPolynomial<F, R> &
BiBernsteinPolynomial<F, R>::operator=(BiBernsteinPolynomial const &q)
{
    BivariatePolynomial<F, R>::operator=(q);
    return (*this);
}


template <typename F, typename R>
BiBernsteinPolynomial<F, R>::~BiBernsteinPolynomial()
{
}

template <typename F, typename R>
F
BiBernsteinPolynomial<F, R>::eval(
    R const &x,
    R const &y) const
{
    /* first, we insert y and obtain a univariate polynomial in x, whose coefficients
     * can be computed by applying the de-Casteljau algorithm for each of the (m + 1) rows of 
     * the coefficient matrix. we again apply the de-Casteljau algorithm for the resulting
     * polynomial for x to get the value */
    uint32_t    i;
    Vector<F>   d(this->m + 1);

    for (i = 0; i < this->m + 1; i++) {
        d[i] = deCasteljau<F, R>(this->n, this->coeff.getRow(i), y);
    }
    return deCasteljau(this->m, d, x);
}

/* arithmetic */
template <typename F, typename R>
BiBernsteinPolynomial<F, R>
BiBernsteinPolynomial<F, R>::operator+(BiBernsteinPolynomial<F, R> const &q) const
{
    return ( BiBernsteinPolynomial(this->coeff) += q );
}

template <typename F, typename R>
BiBernsteinPolynomial<F, R> &
BiBernsteinPolynomial<F, R>::operator+=(BiBernsteinPolynomial<F, R> const &q)
{
    BivariatePolynomial<F, R>::operator+=(q);
    return (*this);
}

template <typename F, typename R>
BiBernsteinPolynomial<F, R>
BiBernsteinPolynomial<F, R>::operator-(BiBernsteinPolynomial<F, R> const &q) const
{
    return ( BiBernsteinPolynomial(this->coeff) -= q );
}

template <typename F, typename R>
BiBernsteinPolynomial<F, R> &
BiBernsteinPolynomial<F, R>::operator-=(BiBernsteinPolynomial<F, R> const &q)
{
    BivariatePolynomial<F, R>::operator-=(q);
    return (*this);
}

template <typename F, typename R>
BiBernsteinPolynomial<F, R>
BiBernsteinPolynomial<F, R>::operator*(F const &x) const
{
    return ( BiBernsteinPolynomial(this->coeff) *= x );
}

template <typename F, typename R>
BiBernsteinPolynomial<F, R> &
BiBernsteinPolynomial<F, R>::operator*=(F const &x)
{
    BivariatePolynomial<F, R>::operator*=(x);
    return (*this);
}

template <typename F, typename R>
BiBernsteinPolynomial<F, R>
BiBernsteinPolynomial<F, R>::operator/(F const &x) const
{
    return ( BiBernsteinPolynomial(this->coeff) /= x );
}

template <typename F, typename R>
BiBernsteinPolynomial<F, R> &
BiBernsteinPolynomial<F, R>::operator/=(F const &x)
{
    BivariatePolynomial<F, R>::operator/=(x);
    return (*this);
}

/* FIXME: broken right now. */
/* inner product with respect to L2 norm \int_0^1\int_0^1{p(x,y)q(x,y) dxdy} */
template <typename F, typename R>
F
BiBernsteinPolynomial<F, R>::operator*(BiBernsteinPolynomial<F, R> const &q) const
{
    uint32_t                        i, j, k, l;
    F                               product = 0;
    BiBernsteinPolynomial<F, R> const &p = (*this);

    /*  quadruple loop. well... those polynomial coefficient matrices are effectively treated as
     *  vectors, and since BB(m, n) is no orthogonal base, we need to take all n^4 terms into
     *  account */
    for (i = 0; i < p.m + 1; i++) {
        for (j = 0; j < p.n + 1; j++) {
            for (k = 0; k < q.m + 1; k++) {
                for (l = 0; l < q.n + 1; l++) {
                    //product += p(i, j) * q(k, l) * BBInnerProducts[this->m][this->m][i][k] * BBInnerProducts[this->n][this->n][j][l];
                }
            }
        }
    }
    return product;
}


/* multiplication as polynomials over |R^2 or a suitably defined ring */
template <typename F, typename R>
BiBernsteinPolynomial<F, R>
BiBernsteinPolynomial<F, R>::multiply(BiBernsteinPolynomial const &q) const
{
    using Aux::Numbers::bicof;

    int                         i, j, k, l, m1 = this->m, m2 = q.m, n1 = this->n, n2 = q.n;
    Matrix<F>                   r_coeff(m1 + m2 + 1, n1 + n2 + 1);

    /* reference for better readability */
    BiBernsteinPolynomial const &p = (*this);

    for (k = 0; k < m1 + m2 + 1; k++) {
        for (l = 0; l < n1 + n2 + 1; l++) {
            //printf("calculating coeff: (k, l) = (%d, %d)... ", k, l);
            /* init coefficient to 0.0 */
            r_coeff(k, l) = 0.0;

            /* sum up all terms for (k, l), a direct generalization of the univariate case */
            for (i = std::max(0, k - m2); i < std::min(k , m1) + 1; i++) {
                for (j = std::max(0, l - n2); j < std::min(l, n1) + 1; j++) {
                    //printf("considering (i, j) = (%d, %d)\n", i, j);
                    r_coeff(k, l) +=
                        ( p(i, j) * q(k - i, l - j) * bicof<F>(m1, i) * bicof<F>(m2, k - i) * bicof<F>(n1, j) * bicof<F>(n2, l - j)) /
                        (bicof<F>(m1 + m2, k) * bicof<F>(n1 + n2, l));
                }
            }

            //printf("%f\n", r_coeff(k, l));
        }
    }

    return ( BiBernsteinPolynomial(r_coeff) );;
}


template <typename F, typename R>
BiBernsteinPolynomial<F, R>
BiBernsteinPolynomial<F, R>::square() const
{
    return (this->multiply( (*this) ));
}


/* to split the polynomial over [0,1]^2 at some value x into two polynomials representing it
 * in [0,x]x[0,1] and [x,1]x[0,1], we apply deCasteljauSplit() on the columns of the coefficient
 * matrix and get the two new coefficient matrices */
template <typename F, typename R>
void
BiBernsteinPolynomial<F, R>::split_x(
    R const                        &x,
    BiBernsteinPolynomial<F, R>    *pleft,
    BiBernsteinPolynomial<F, R>    *pright) const
{
    uint32_t    i, j;
    Vector<F>   col_j_left, col_j_right;
    Matrix<F>   coeff_left(this->m + 1, this->n + 1),
                coeff_right(this->m + 1, this->n + 1);

    for (j = 0; j < this->n + 1; j++) {
        deCasteljauSplit<F, R>(this->m, this->coeff.getCol(j), x, col_j_left, col_j_right);
        for (i = 0; i < this->m + 1; i++) {
            coeff_left(i, j)    = col_j_left[i];
            coeff_right(i, j)   = col_j_right[i];
        }
    }

    if (pleft)  pleft->setCoeffs(coeff_left);
    if (pright) pright->setCoeffs(coeff_right);
}

/* to split the polynomial over [0,1]^2 at some value y into two polynomials representing it
 * in [0,1]x[0,y] and [0,1]x[y,1], we apply deCasteljauSplit() on the rows of the coefficient
 * matrix and get the two new coefficient matrices */
template <typename F, typename R>
void
BiBernsteinPolynomial<F, R>::split_y(
    R const                        &y,
    BiBernsteinPolynomial<F, R>    *pdown,
    BiBernsteinPolynomial<F, R>    *pup) const
{
    uint32_t    i, j;
    Vector<F>   row_i_up, row_i_down;
    Matrix<F>           coeff_up(this->m + 1, this->n + 1),
                        coeff_down(this->m + 1, this->n + 1);

    for (i = 0; i < this->m + 1; i++) {
        deCasteljauSplit<F, R>(this->n, this->coeff.getRow(i), y, row_i_down, row_i_up);
        for (j = 0; j < this->n + 1; j++) {
            coeff_up(i, j)      = row_i_up[j];
            coeff_down(i, j)    = row_i_down[j];
        }
    }

    if (pdown)  pdown->setCoeffs(coeff_down);
    if (pup)    pup->setCoeffs(coeff_up);
}

/* combined method to split at (x,y) in [0,1]^2 */
template <typename F, typename R>
void
BiBernsteinPolynomial<F, R>::split_xy(
    R const                        &x,
    R const                        &y,
    BiBernsteinPolynomial<F, R>    *pleft_down,
    BiBernsteinPolynomial<F, R>    *pright_down,
    BiBernsteinPolynomial<F, R>    *pright_up,
    BiBernsteinPolynomial<F, R>    *pleft_up) const
{
    BiBernsteinPolynomial<F, R> pl_d, pr_d, pr_u, pl_u;

    this->split_x(x, &pl_d, &pr_d);
    pl_d.split_y(y, &pl_d, &pl_u); 
    pr_d.split_y(y, &pr_d, &pr_u); 

    if (pleft_down)     pleft_down->setCoeffs(pl_d.getCoeffs());
    if (pright_down)    pright_down->setCoeffs(pr_d.getCoeffs());
    if (pright_up)      pright_up->setCoeffs(pr_u.getCoeffs());
    if (pleft_up)       pleft_up->setCoeffs(pl_u.getCoeffs());
}


/* clip to interval [x0, x1]x[y0, y1] */
template <typename F, typename R>
void
BiBernsteinPolynomial<F, R>::clipToInterval(
    R const                        &x0, 
    R const                        &x1,
    R const                        &y0,
    R const                        &y1,
    BiBernsteinPolynomial<F, R>    *pclip) const
{
    this->split_x (     x0,                     NULL,   pclip);
    pclip->split_x(    (x1 - x0) / (1.0 - x0),  pclip,  NULL );
    pclip->split_y(     y0,                     NULL,   pclip);
    pclip->split_y(    (y1 - y0) / (1.0 - y0),  pclip,  NULL);
}


/* degree elevation, generalization of univariate case, effectively multiplying with 1(x, y) in
 * BB(r, s) */
template <typename F, typename R>
void
BiBernsteinPolynomial<F, R>::elevateDegree(
    uint32_t r_arg,
    uint32_t s_arg)
{
    using Aux::Numbers::bicof;

    int         i, j, k, l, r = r_arg, s = s_arg;

    int         m = this->m;
    int         n = this->n;

    Matrix<F>   elev_coeff(m + r + 1, n + s + 1);

    for (k = 0; k < m + r + 1; k++) {
        for (l = 0; l < n + s + 1; l++) {
            /* init coefficient to 0 */
            elev_coeff(k, l) = 0;

            /* sum up all terms for (k, l), a direct generalization of the univariate case */
            for (i = std::max(0, k - r); i < std::min(k , m) + 1; i++) {
                for (j = std::max(0, l - s); j < std::min(l, n) + 1; j++) {
                    elev_coeff(k, l) +=
                        ( this->coeff(i, j) * bicof<F>(m, i) * bicof<F>(r, k - i) * bicof<F>(n, j) * bicof<F>(s, l - j))
                        /
                        (bicof<R>(m + r, k) * bicof<F>(n + s, l));
                }
            }
        }
    }

    /* commit new values, qualified access is used as usual (so it doesnt matter that m and n have
     * been shadowed above) */
    this->m         = m + r;
    this->n         = n + s;
    this->coeff     = elev_coeff;
}

template <typename F, typename R>
void
BiBernsteinPolynomial<F, R>::convertFromPowerBasis(Matrix<F> const &coeff_pow)
{
    using Aux::Numbers::bicof;

    uint32_t i, j, k, l, coeff_m, coeff_n;
    coeff_m = coeff_pow.numRows();
    coeff_n = coeff_pow.numCols();

    /* resize matrix and init to 0.0 */
    this->coeff.resize(coeff_m, coeff_n);
    this->coeff.fill(0.0);

    this->m = coeff_m - 1;
    this->n = coeff_n - 1;

    for (k = 0; k <= this->m; k++) {
        for (l = 0; l <= this->n; l++) {
            this->coeff(k, l) = 0.0;
            for (i = 0; i <= k; i++) {
                for (j = 0; j <= l; j++) {
                    this->coeff(k, l) +=
                        coeff_pow(i, j) * bicof<F>(k, i) * bicof<F>(l, j) /
                        (bicof<F>(this->m, i) * bicof<F>(this->n, j));
                }
            }
        }
    }
}


/* ----------------------------------------------------------------------------------------------------------------- *
 *
 *                              auxiliary algorithms.                                                                 
 *
 * ----------------------------------------------------------------------------------------------------------------- */

/* B_{ij}^{m, n} in BB(m, n) => coefficient matrix C with only non-zero element C(i,j) = 1.0 */
template <typename F, typename R>
BiBernsteinPolynomial<F, R>
BiBB(uint32_t m, uint32_t n, uint32_t i, uint32_t j)
{
    return ( BiBernsteinPolynomial<F, R>(Aux::VecMat::kronecker_mat<F>(m+1, n+1, i, j)) );
}
