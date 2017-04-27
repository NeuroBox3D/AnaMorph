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

/* ----------------------------------------------------------------------------------------------------------------- *
 *
 *                        implementation of abstract bivariate polynomial class.....                                  
 *
 * ----------------------------------------------------------------------------------------------------------------- */

template <uint32_t deg1, uint32_t deg2, typename F, typename R>
BivariatePolynomial<deg1, deg2, F, R>::BivariatePolynomial()
{
    coeff.fill(0);
}

template <uint32_t deg1, uint32_t deg2, typename F, typename R>
BivariatePolynomial<deg1, deg2, F, R>::BivariatePolynomial(const F& x)
{
    coeff.fill(x);
}

template <uint32_t deg1, uint32_t deg2, typename F, typename R>
BivariatePolynomial<deg1, deg2, F, R>::BivariatePolynomial(const this_type& q)
{
    coeff = q.coeff;
}

template <uint32_t deg1, uint32_t deg2, typename F, typename R>
BivariatePolynomial<deg1, deg2, F, R>::~BivariatePolynomial()
{}


template <uint32_t deg1, uint32_t deg2, typename F, typename R>
BivariatePolynomial<deg1, deg2, F, R>&
BivariatePolynomial<deg1, deg2, F, R>::operator=(const this_type& q)
{
    coeff = q.coeff;
    return *this;
}

template <uint32_t deg1, uint32_t deg2, typename F, typename R>
F&
BivariatePolynomial<deg1, deg2, F, R>::operator()(uint32_t i, uint32_t j)
{
    return coeff(i, j);
}

template <uint32_t deg1, uint32_t deg2, typename F, typename R>
F
BivariatePolynomial<deg1, deg2, F, R>::operator()(uint32_t i, uint32_t j) const
{
    return coeff(i, j);
}

template <uint32_t deg1, uint32_t deg2, typename F, typename R>
const StaticMatrix<deg1+1, deg2+1, F>&
BivariatePolynomial<deg1, deg2, F, R>::getCoeffs() const
{
    return coeff;
}

template <uint32_t deg1, uint32_t deg2, typename F, typename R>
StaticMatrix<deg1+1, deg2+1, F>&
BivariatePolynomial<deg1, deg2, F, R>::getCoeffs()
{
    return coeff;
}


template <uint32_t deg1, uint32_t deg2, typename F, typename R>
void
BivariatePolynomial<deg1, deg2, F, R>::setCoeffs(const coeff_type& _coeff)
{
    coeff = _coeff;
}

template <uint32_t deg1, uint32_t deg2, typename F, typename R>
F
BivariatePolynomial<deg1, deg2, F, R>::getMaxAbsCoeff() const
{
    uint32_t    i, j;
    F           cmax = std::abs(coeff(0, 0));
    F           abs_ij;
    for (i = 0; i < deg1 + 1; i++) {
        for (j = 0; j < deg2 + 1; j++) {
            abs_ij = std::abs(this->coeff(i, j));
            if (abs_ij > cmax) {
                cmax = abs_ij;
            }
        }
    }
    return cmax;
}

template <uint32_t deg1, uint32_t deg2, typename F, typename R>
void
BivariatePolynomial<deg1, deg2, F, R>::initConstant(const F& x)
{
    coeff.fill(x);
}

template <uint32_t deg1, uint32_t deg2, typename F, typename R>
void
BivariatePolynomial<deg1, deg2, F, R>::zero()
{
    this->coeff.fill(0);
}

template <uint32_t deg1, uint32_t deg2, typename F, typename R>
BivariatePolynomial<deg1, deg2, F, R>&
BivariatePolynomial<deg1, deg2, F, R>::operator+=(const this_type& q)
{
    coeff += q.coeff;
    return *this;
}

template <uint32_t deg1, uint32_t deg2, typename F, typename R>
BivariatePolynomial<deg1, deg2, F, R>&
BivariatePolynomial<deg1, deg2, F, R>::operator-=(const this_type& q)
{
    coeff -= q.coeff;
    return *this;
}

template <uint32_t deg1, uint32_t deg2, typename F, typename R>
BivariatePolynomial<deg1, deg2, F, R>&
BivariatePolynomial<deg1, deg2, F, R>::operator*=(const F& x)
{
    coeff *= x;
    return *this;
}

template <uint32_t deg1, uint32_t deg2, typename F, typename R>
BivariatePolynomial<deg1, deg2, F, R>&
BivariatePolynomial<deg1, deg2, F, R>::operator/=(const F& x)
{
    coeff /= x;
    return *this;
}



template <uint32_t deg1, uint32_t deg2, typename F, typename R>
void
BivariatePolynomial<deg1, deg2, F, R>::printCoeff() const
{
    PrintCoeffImpl<F, R>(*this);
}


template <uint32_t deg1, uint32_t deg2, typename F, typename R>
template <typename dummy>
BivariatePolynomial<deg1, deg2, F, R>::PrintCoeffImpl<float, float, dummy>::PrintCoeffImpl
(const BivariatePolynomial<deg1, deg2, F, R>& p)
{
    uint32_t i, j;
    for (i = 0; i < deg1 + 1; ++i)
    {
        printf("( ");
        for (j = 0; j < deg2 + 1; ++j)
            printf("%+20.13E ", p.coeff(i, j));

        printf(")\n");
    }
}

template <uint32_t deg1, uint32_t deg2, typename F, typename R>
template <typename dummy>
BivariatePolynomial<deg1, deg2, F, R>::PrintCoeffImpl<double, double, dummy>::PrintCoeffImpl
(const BivariatePolynomial<deg1, deg2, F, R>& p)
{
    uint32_t i, j;
    for (i = 0; i < deg1 + 1; ++i)
    {
        printf("( ");
        for (j = 0; j < deg2 + 1; ++j)
            printf("%+20.13E ", p.coeff(i, j));

        printf(")\n");
    }
}


template <uint32_t deg1, uint32_t deg2, typename F, typename R>
void
BivariatePolynomial<deg1, deg2, F, R>::writePlotFile(uint32_t ticks, const std::string& filename) const
{
    WritePlotFileImpl<F, R>(ticks, filename, *this);
}


template <uint32_t deg1, uint32_t deg2, typename F, typename R>
template <typename dummy>
BivariatePolynomial<deg1, deg2, F, R>::WritePlotFileImpl<float, float, dummy>::WritePlotFileImpl
(uint32_t ticks, const std::string& filename, const BivariatePolynomial<deg1, deg2, F, R>& p)
{
    FILE *fout = fopen(filename.c_str(), "w");
    if (!fout)
    {
        fprintf(stderr, "BivariatePolynomial::writePlotFile(): can't open file \"%s\" for writing.\n", filename.c_str());
        return;
    }

    uint32_t i, j;

    for (i = 0; i <= ticks; ++i)
    {
        float x = (float)i / (float)ticks;
        for (j = 0; j <= ticks; ++j)
        {
            float y = (float)j / (float)ticks;
            fprintf(fout, "%12.5E %12.5E %12.5E\n", x, y, p.eval(x, y));
        }
        fprintf(fout, "\n");
    }
    fclose(fout);
}

template <uint32_t deg1, uint32_t deg2, typename F, typename R>
template <typename dummy>
BivariatePolynomial<deg1, deg2, F, R>::WritePlotFileImpl<double, double, dummy>::WritePlotFileImpl
(uint32_t ticks, const std::string& filename, const BivariatePolynomial<deg1, deg2, F, R>& p)
{
    FILE *fout = fopen(filename.c_str(), "w");
    if (!fout)
    {
        fprintf(stderr, "BivariatePolynomial::writePlotFile(): can't open file \"%s\" for writing.\n", filename.c_str());
        return;
    }

    uint32_t i, j;
    float x, y;

    for (i = 0; i <= ticks; ++i)
    {
        x = (float)i / (float)ticks;
        for (j = 0; j <= ticks; ++j)
        {
            y = (float)j / (float)ticks;
            fprintf(fout, "%12.5E %12.5E %12.5E\n", x, y, p.eval(x, y));
        }
        fprintf(fout, "\n");
    }
    fclose(fout);
}




/* ----------------------------------------------------------------------------------------------------------------- *
 *
 *                        implementation of bivariate bernstein polynomial class....                                  
 *
 * ----------------------------------------------------------------------------------------------------------------- */

template <uint32_t deg1, uint32_t deg2, typename F, typename R>
BiBernsteinPolynomial<deg1, deg2, F, R>::BiBernsteinPolynomial()
: BivariatePolynomial<deg1, deg2, F, R>()
{}


template <uint32_t deg1, uint32_t deg2, typename F, typename R>
BiBernsteinPolynomial<deg1, deg2, F, R>::BiBernsteinPolynomial(const F& x)
: BivariatePolynomial<deg1, deg2, F, R>(x)
{}


template <uint32_t deg1, uint32_t deg2, typename F, typename R>
BiBernsteinPolynomial<deg1, deg2, F, R>::BiBernsteinPolynomial(const coeff_type& _coeff)
{
    coeff = _coeff;
}

template <uint32_t deg1, uint32_t deg2, typename F, typename R>
BiBernsteinPolynomial<deg1, deg2, F, R>::BiBernsteinPolynomial(const this_type& q)
: BivariatePolynomial<deg1, deg2, F, R>(q)
{}

template <uint32_t deg1, uint32_t deg2, typename F, typename R>
BiBernsteinPolynomial<deg1, deg2, F, R>::~BiBernsteinPolynomial()
{}

template <uint32_t deg1, uint32_t deg2, typename F, typename R>
BiBernsteinPolynomial<deg1, deg2, F, R>&
BiBernsteinPolynomial<deg1, deg2, F, R>::operator=(const this_type& q)
{
    BivariatePolynomial<deg1, deg2, F, R>::operator=(q);
    return *this;
}


template <uint32_t deg1, uint32_t deg2, typename F, typename R>
F
BiBernsteinPolynomial<deg1, deg2, F, R>::eval(const R& x, const R& y) const
{
    /* first, we insert y and obtain a univariate polynomial in x, whose coefficients
     * can be computed by applying the de-Casteljau algorithm for each of the (m + 1) rows of 
     * the coefficient matrix. we again apply the de-Casteljau algorithm for the resulting
     * polynomial for x to get the value */
    uint32_t i;
    StaticVector<deg1+1,F> d;

    for (i = 0; i < deg1+1; ++i)
        d[i] = deCasteljau<deg2, F, R>(coeff.getRow(i), y);

    return deCasteljau<deg1>(d, x);
}

/* arithmetic */
template <uint32_t deg1, uint32_t deg2, typename F, typename R>
BiBernsteinPolynomial<deg1, deg2, F, R>
BiBernsteinPolynomial<deg1, deg2, F, R>::operator+(const this_type& q) const
{
    return BiBernsteinPolynomial(coeff) += q;
}

template <uint32_t deg1, uint32_t deg2, typename F, typename R>
BiBernsteinPolynomial<deg1, deg2, F, R>&
BiBernsteinPolynomial<deg1, deg2, F, R>::operator+=(const this_type& q)
{
    BivariatePolynomial<deg1, deg2, F, R>::operator+=(q);
    return *this;
}

template <uint32_t deg1, uint32_t deg2, typename F, typename R>
BiBernsteinPolynomial<deg1, deg2, F, R>
BiBernsteinPolynomial<deg1, deg2, F, R>::operator-(const this_type& q) const
{
    return BiBernsteinPolynomial(coeff) -= q;
}

template <uint32_t deg1, uint32_t deg2, typename F, typename R>
BiBernsteinPolynomial<deg1, deg2, F, R>&
BiBernsteinPolynomial<deg1, deg2, F, R>::operator-=(const this_type& q)
{
    BivariatePolynomial<deg1, deg2, F, R>::operator-=(q);
    return *this;
}

template <uint32_t deg1, uint32_t deg2, typename F, typename R>
BiBernsteinPolynomial<deg1, deg2, F, R>
BiBernsteinPolynomial<deg1, deg2, F, R>::operator*(const F& x) const
{
    return BiBernsteinPolynomial(coeff) *= x;
}

template <uint32_t deg1, uint32_t deg2, typename F, typename R>
BiBernsteinPolynomial<deg1, deg2, F, R>&
BiBernsteinPolynomial<deg1, deg2, F, R>::operator*=(const F& x)
{
    BivariatePolynomial<deg1, deg2, F, R>::operator*=(x);
    return *this;
}

template <uint32_t deg1, uint32_t deg2, typename F, typename R>
BiBernsteinPolynomial<deg1, deg2, F, R>
BiBernsteinPolynomial<deg1, deg2, F, R>::operator/(const F& x) const
{
    return BiBernsteinPolynomial(coeff) /= x;
}

template <uint32_t deg1, uint32_t deg2, typename F, typename R>
BiBernsteinPolynomial<deg1, deg2, F, R>&
BiBernsteinPolynomial<deg1, deg2, F, R>::operator/=(const F& x)
{
    BivariatePolynomial<deg1, deg2, F, R>::operator/=(x);
    return *this;
}

// FIXME: broken right now.
#if 0
template <uint32_t deg1, uint32_t deg2, typename F, typename R>
template <uint32_t d1, uint32_t d2>
F
BiBernsteinPolynomial<deg1, deg2, F, R>::operator*(const BiBernsteinPolynomial<d1, d2, F, R>& q) const
{
    uint32_t i, j, k, l;
    F product = 0;
    const this_type& p = *this;

    /*  quadruple loop. well... those polynomial coefficient matrices are effectively treated as
     *  vectors, and since BB(m, n) is no orthogonal base, we need to take all n^4 terms into
     *  account */
    for (i = 0; i < deg1+1; i++) {
        for (j = 0; j < deg2+1; j++) {
            for (k = 0; k < d1+1; k++) {
                for (l = 0; l < d2+1; l++) {
                    //product += p(i, j) * q(k, l) * BBInnerProducts[this->m][this->m][i][k] * BBInnerProducts[this->n][this->n][j][l];
                }
            }
        }
    }
    return product;
}
#endif

template <uint32_t deg1, uint32_t deg2, typename F, typename R>
template <uint32_t d1, uint32_t d2>
BiBernsteinPolynomial<deg1+d1, deg2+d2, F, R>
BiBernsteinPolynomial<deg1, deg2, F, R>::multiply(const BiBernsteinPolynomial<d1, d2, F, R>& q) const
{
    using Aux::Numbers::bicof;

    int i, j, k, l;
    StaticMatrix<deg1+d1+1, deg2+d2+1, F> r_coeff;

    // reference for better readability
    const this_type& p = *this;

    for (k = 0; k < (int)deg1 + (int)d1 + 1; ++k)
    {
        for (l = 0; l < (int)deg2 + (int)d2 + 1; ++l)
        {
            r_coeff(k, l) = 0;

            // sum up all terms for (k, l), a direct generalization of the univariate case
            for (i = std::max(0, k - (int)d1); i < std::min(k, (int)deg1) + 1; ++i)
            {
                for (j = std::max(0, l - (int)d2); j < std::min(l, (int)deg2) + 1; ++j)
                {
                    r_coeff(k, l) +=
                        p(i, j) * q(k - i, l - j) * bicof<F, deg1>(i) * bicof<F, d1>(k - i) * bicof<F, deg2>(j) * bicof<F, d2>(l - j)
                        / (bicof<F, deg1+d1>(k) * bicof<F, deg2+d2>(l));
                }
            }
        }
    }

    return BiBernsteinPolynomial<deg1+d1, deg2+d2, F, R>(r_coeff);
}


template <uint32_t deg1, uint32_t deg2, typename F, typename R>
BiBernsteinPolynomial<2*deg1, 2*deg2, F, R>
BiBernsteinPolynomial<deg1, deg2, F, R>::square() const
{
    return multiply(*this);
}


/* to split the polynomial over [0,1]^2 at some value x into two polynomials representing it
 * in [0,x]x[0,1] and [x,1]x[0,1], we apply deCasteljauSplit() on the columns of the coefficient
 * matrix and get the two new coefficient matrices */
template <uint32_t deg1, uint32_t deg2, typename F, typename R>
void
BiBernsteinPolynomial<deg1, deg2, F, R>::split_x(const R& x, this_type* pleft, this_type* pright) const
{
    // In contrast to split_y, we have to move some values around here
    // due to the fact that we cannot get a reference to a column of a matrix
    // in our row-storage implementation.
    uint32_t i, j;
    typename coeff_type::col_type col_j_left, col_j_right;

    if (pleft)
    {
        coeff_type& coeff_left = pleft->getCoeffs();

        if (pright)
        {
            coeff_type& coeff_right = pright->getCoeffs();
            for (j = 0; j < deg2 + 1; ++j)
            {
                deCasteljauSplit<deg1, F, R>(coeff.getCol(j), x, col_j_left, col_j_right);
                for (i = 0; i < deg1 + 1; ++i)
                {
                    coeff_left(i, j) = col_j_left[i];
                    coeff_right(i, j) = col_j_right[i];
                }
            }
        }
        else
        {
            coeff_type coeff_right;
            for (j = 0; j < deg2 + 1; ++j)
            {
                deCasteljauSplit<deg1, F, R>(coeff.getCol(j), x, col_j_left, col_j_right);
                for (i = 0; i < deg1 + 1; ++i)
                {
                    coeff_left(i, j) = col_j_left[i];
                    coeff_right(i, j) = col_j_right[i];
                }
            }
        }
    }
    else
    {
        if (pright)
        {
            coeff_type coeff_left;
            coeff_type& coeff_right = pright->getCoeffs();
            for (j = 0; j < deg2 + 1; ++j)
            {
                deCasteljauSplit<deg1, F, R>(coeff.getCol(j), x, col_j_left, col_j_right);
                for (i = 0; i < deg1 + 1; ++i)
                {
                    coeff_left(i, j) = col_j_left[i];
                    coeff_right(i, j) = col_j_right[i];
                }
            }
        }
        // else do nothing (meaningless)
    }
}

/* to split the polynomial over [0,1]^2 at some value y into two polynomials representing it
 * in [0,1]x[0,y] and [0,1]x[y,1], we apply deCasteljauSplit() on the rows of the coefficient
 * matrix and get the two new coefficient matrices */
template <uint32_t deg1, uint32_t deg2, typename F, typename R>
void
BiBernsteinPolynomial<deg1, deg2, F, R>::split_y(const R& y, this_type* pdown, this_type* pup) const
{
    if (pdown)
    {
        if (pup)
        {
            for (uint32_t i = 0; i < deg1 + 1; ++i)
            {
                const typename coeff_type::row_type& row = coeff.getRow(i);
                typename coeff_type::row_type& rowD = pdown->getCoeffs().getRow(i);
                typename coeff_type::row_type& rowU = pup->getCoeffs().getRow(i);

                deCasteljauSplit<deg2, F, R>(row , y, rowD, rowU);
            }
        }
        else
        {
            typename coeff_type::row_type rowU;
            for (uint32_t i = 0; i < deg1 + 1; ++i)
            {
                const typename coeff_type::row_type& row = coeff.getRow(i);
                typename coeff_type::row_type& rowD = pdown->getCoeffs().getRow(i);

                deCasteljauSplit<deg2, F, R>(row , y, rowD, rowU);
            }
        }
    }
    else
    {
        if (pup)
        {
            typename coeff_type::row_type rowD;
            for (uint32_t i = 0; i < deg1 + 1; ++i)
            {
               const typename coeff_type::row_type& row = coeff.getRow(i);
               typename coeff_type::row_type& rowU = pup->getCoeffs().getRow(i);

               deCasteljauSplit<deg2, F, R>(row , y, rowD, rowU);
            }
        }
        // else do nothing (meaningless)
    }

}

/* combined method to split at (x,y) in [0,1]^2 */
template <uint32_t deg1, uint32_t deg2, typename F, typename R>
void
BiBernsteinPolynomial<deg1, deg2, F, R>::split_xy
(
    const R& x,
    const R& y,
    this_type* pleft_down,
    this_type* pright_down,
    this_type* pright_up,
    this_type* pleft_up
) const
{
    if (pleft_down)
    {
        split_x(x, pleft_down, pright_down);
        pleft_down->split_y(y, pleft_down, pleft_up);
        pleft_down->split_y(y, pright_down, pright_up);
    }
    else
    {
        this_type tmp;
        split_x(x, &tmp, pright_down);
        tmp.split_y(y, &tmp, pleft_up);
        tmp.split_y(y, pright_down, pright_up);
    }
}


/* clip to interval [x0, x1]x[y0, y1] */
template <uint32_t deg1, uint32_t deg2, typename F, typename R>
void
BiBernsteinPolynomial<deg1, deg2, F, R>::clipToInterval
(
    const R& x0,
    const R& x1,
    const R& y0,
    const R& y1,
    this_type* pclip
) const
{
    split_x (x0, NULL, pclip);
    pclip->split_x((x1 - x0) / (1.0 - x0), pclip, NULL);
    pclip->split_y(y0, NULL, pclip);
    pclip->split_y((y1 - y0) / (1.0 - y0), pclip, NULL);
}


/* degree elevation, generalization of univariate case, effectively multiplying with 1(x, y) in
 * BB(r, s) */
template <uint32_t deg1, uint32_t deg2, typename F, typename R>
template <uint32_t d1, uint32_t d2>
BiBernsteinPolynomial<deg1+d1, deg2+d2, F, R>
BiBernsteinPolynomial<deg1, deg2, F, R>::elevateDegree()
{
    using Aux::Numbers::bicof;

    int i, j, k, l;
    StaticMatrix<deg1+d1+1, deg2+d2+1, F> elev_coeff;

    for (k = 0; k < (int)deg1 + (int)d1 + 1; ++k)
    {
        for (l = 0; l < (int)deg2 + (int)d2 + 1; ++l)
        {
            elev_coeff(k, l) = 0;

            for (i = std::max(0, k - (int)d1); i < std::min(k, (int)deg1) + 1; ++i)
            {
                for (j = std::max(0, l - (int)d2); j < std::min(l, (int)deg2) + 1; ++j)
                {
                    elev_coeff(k, l) +=
                        this->coeff(i, j) * bicof<F, deg1>(i) * bicof<F, d1>(k - i) * bicof<F, deg2>(j) * bicof<F, d2>(l - j)
                        / (bicof<R, deg1+d1>(k) * bicof<F, deg2+d2>(l));
                }
            }
        }
    }

    BiBernsteinPolynomial<deg1+d1, deg2+d2, F, R> r;
    r.setCoeffs(elev_coeff);
    return r;
}

template <uint32_t deg1, uint32_t deg2, typename F, typename R>
void
BiBernsteinPolynomial<deg1, deg2, F, R>::convertFromPowerBasis(const coeff_type& coeff_pow)
{
    uint32_t i, j, k, l;

    for (k = 0; k <= deg1; ++k)
    {
        for (l = 0; l <= deg2; ++l)
        {
            R ifac = (R)1;
            R jfac = (R)1;

            // first step (i==0)
            coeff(k, l) = coeff_pow(0, 0);
            for (j = 1; j <= l; ++j)
            {
                jfac *= (R)(l-j+1) / (R)(deg1-j+1);  // jfac = bicof(l,j) / bicof(deg2,j)
                coeff(k, l) += jfac * coeff_pow(0, j);
            }

            // rest
            for (i = 1; i <= k; ++i)
            {
                ifac *= (R)(k-i+1) / (R)(deg1-i+1); // ifac = bicof(k,i) / bicof(deg1,i)
                jfac = (R)1;
                coeff(k, l) += ifac * coeff_pow(i, 0);
                for (j = 1; j <= l; ++j)
                {
                    jfac *= (R)(l-j+1) / (R)(deg1-j+1);  // jfac = bicof(l,j) / bicof(deg2,j)
                    coeff(k, l) += ifac * jfac * coeff_pow(i, j);
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
template <uint32_t deg1, uint32_t deg2, typename F, typename R>
BiBernsteinPolynomial<deg1, deg2, F, R>
BiBB(uint32_t i, uint32_t j)
{
    return BiBernsteinPolynomial<deg1, deg2, F, R>(Aux::VecMat::kronecker_static_mat<deg1+1, deg2+1, F>(i, j));
}
