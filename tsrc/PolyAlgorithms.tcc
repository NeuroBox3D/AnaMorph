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
#include "Polynomial.hh"

/* ----------------------------------------------------------------------------------------------------------------- *
 *
 *          polynomial basis conversion                                                                                                    
 *
 * ----------------------------------------------------------------------------------------------------------------- */

template <typename R>
void
convertBasis(
    BernsteinPolynomial<R, R>      &b,
    PowerPolynomial<R, R> const    &p)
{
    using Aux::Numbers::bicof;

    uint32_t        i, j;
    const uint32_t  n  = p.getDegree();
    Vector<R>       b_coeff(n + 1, 0);

    //debugl("\n");
    for (i = 0; i < n + 1; i++) {
        b_coeff[i] = 0;
        for (j = 0; j < i + 1; j++) {
            //debugl("bicof[i = %d][j = %d]: %f, bicof[degree = %d][j = %d]: %f\n", i, j, bicof[i][j], degree, j, bicof[degree][j]);
            b_coeff[i] += (bicof<R>(i, j) * p[j]) / bicof<R>(n, j);
        }
        //debugl("b[i] = b[%d] = %f\n\n", i, b[i]);
    }

    b = BernsteinPolynomial<R, R>(b_coeff);
}

template <typename R>
void
convertBasis(
    PowerPolynomial<R, R>              &p,
    BernsteinPolynomial<R, R> const    &b)
{
    using Aux::Numbers::bicof;

    uint32_t        i, j;
    uint32_t const  n = b.getDegree();
    R               sign;
    Vector<R>       p_coeff(n + 1, 0);

    for (i = 0; i < n + 1; i++) {
        p_coeff[i] = 0;
        for (j = 0; j < i + 1; j++) {
            sign        = (R)( ( (i - j) % 2 == 0) ? 1 : -1);
            p_coeff[i] += sign * bicof<R>(n, i) * bicof<R>(i, j) * b[j];
        }
    }

    p = PowerPolynomial<R, R>(p_coeff);
}

/* ----------------------------------------------------------------------------------------------------------------- *
 *
 *          generate bivariate BiBernsteinPolynomials from univariate BernsteinPolynomials                                                   
 *
 * ----------------------------------------------------------------------------------------------------------------- */

template<
    typename F,
    typename R
>
BiBernsteinPolynomial<F, R>
BernsteinTensorMultiply(
    BernsteinPolynomial<F, R> const    &p,
    BernsteinPolynomial<F, R> const    &q)
{
    uint32_t i, j, m = p.getDegree(), n = q.getDegree();

    //debugl("tensor multiply: p.degree = %d, q.degree = %d => (%d, %d)\n", m, n, m, n);
    BiBernsteinPolynomial<F, R> r(m, n, 0);
    for (i = 0; i < m + 1; i++) {
        for (j = 0; j < n + 1; j++) {
            r(i, j) = p(i) * q(j);
        }
    }
    return r;
}

template <
    typename F,
    typename R
>
BiBernsteinPolynomial<F, R>
BernsteinConvertToBiPoly(
    BernsteinPolynomial<F, R> const    &p,
    uint32_t                            d,
    bool                                as_x)
{
    uint32_t i, j, k = p.getDegree();

    BiBernsteinPolynomial<F, R> r;
    if (as_x) {
        r = BiBernsteinPolynomial<F, R>(k, d, 0);
        for (i = 0; i < k + 1; i++) {
            for (j = 0; j < d + 1; j++) {
                r(i, j) = p(i);
            }
        }
    }
    else {
        r = BiBernsteinPolynomial<F, R>(d, k, 0);
        for (i = 0; i < d + 1; i++) {
            for (j = 0; j < k + 1; j++) {
                r(i, j) = p(j);
            }
        }
    }
    return r;
}

template <typename F, typename R>
BernsteinPolynomial<F, R>
computeBernsteinBasisPoly(
    uint32_t n,
    uint32_t i)
{
    return BernsteinPolynomial<F, R>(Aux::VecMat::kronecker_vec<F>(n + 1, i));
}

/* ----------------------------------------------------------------------------------------------------------------- *
 *
 *             root finding for univariate polynomials: bezier clipping and required auxiliary algorithms       
 *
 * ----------------------------------------------------------------------------------------------------------------- */
template <typename R>
void
computeBezierControlPolyConvexHull(
    BernsteinPolynomial<R, R> const    &p,
    std::vector<Vec2>                  &cvhull,
    R const                            &eps_slope)
{
    using Aux::Numbers::inf;

    int         i, j, n = p.getDegree();

    uint32_t    max_slope_idx;
    R           tmp, max_slope;

    /* first, compute "upper convex hull", starting with points p[0], which is always in the
     * convex hull. similarly, the point p[n] is also always in the convex hull, since it has
     * maximum x coordinate x = 1.0 */
    cvhull.clear();

    debugl(1, "\n\n Computing convex hull of control polygon.. control points are:\n");
    for (i = 1; i < n + 1; i++) {
        debugl(0, "%2d: %+20.13E %+20.13E\n", i, (R)i / (R)n, p[i]);
    }

    debugl(1, "\n\n---------------------- Computing upper convex hull by scanning to the right from i = 0 to n\n");
    cvhull.push_back( Vec2(0.0, p[0]) );
    i = 0;
    while (i < n) {
        debugl(1, "scanning from i = %d to the right..\n", i);
        max_slope       = -inf<R>();
        max_slope_idx   = i;

        /* find point with biggest "slope" with respect to current point (i/n, p[i]) by
        * computing (delta(y)) / (delta x) */
        for (j = i+1; j < n + 1; j++) {
            tmp = ( (R)n * (p[j] - p[i])) / (R)(j - i);
            debugl(1, "slope between (%d, %d) = %f\n", i, j, tmp);
            if (std::abs(tmp - max_slope) < eps_slope) {
                debugl(1, "two candidates pairs with same slope: (%d, %d) and (%d, %d). choosing current index %d, which is farther away..\n", i, max_slope_idx, i, j, j);
                max_slope       = tmp;
                max_slope_idx   = j;
            }
            else if (tmp > max_slope) {
                max_slope       = tmp;
                max_slope_idx   = j;
                debugl(1, "proper new max slope: max_slope = %f, max_slope_idx = %d\n", max_slope, max_slope_idx);
            }
        }

        /* next point on convex hull is ( (max_slope_idx / n), p[max_slope_idx]) */
        debugl(1, "\nscan finished. next i: %d\n", max_slope_idx);
        i = max_slope_idx;
        cvhull.push_back( Vec2( (R)i / (R)n, p[i]) );
    }

    /* i == degree == n here, last point was inserted. perform backwards scan for "lower" convex hull */
    if (i != n) {
        throw("PolyAlgorithms::computeBezierControlPolyConvexHull(): internal logic error.");
    }

    debugl(1, "\n\n---------------------- Upper convex hull done. Computing lower convex hull by scanning back from i = n to 0\n");

    /* construct "lower" convex hull, symmetric to above construction of the "upper" convex hull */
    i = n;
    while (i > 0) {
        debugl(1, "scanning from i = %d to the left..\n", i);
        max_slope           = -inf<R>();
        max_slope_idx       = i;

        /* find point with biggest slope left from i, symmetric to the case above */
        for (j = i - 1; j >= 0; j--) {
            tmp = ( (R)n * (p[i] - p[j])) / (R)(i - j);
            if (std::abs(tmp - max_slope) < eps_slope) {
                debugl(1, "two candidate pairs with numerically equal slope (vectors from last point to both colinear): (%d, %d) and (%d, %d). choosing current index %d, which is farther away..\n", i, max_slope_idx, i, j, j);
                max_slope       = tmp;
                max_slope_idx   = j;
            }
            else if (tmp > max_slope) {
                max_slope       = tmp;
                max_slope_idx   = j;
                debugl(1, "new proper max slope: max_slope = %f, max_slope_idx = %d\n", max_slope, max_slope_idx);
            }
        }

        /* next point has index max_slope_idx => ( (max_slope_idx / n), p[max_slope_idx] ) */
        i = max_slope_idx;
        cvhull.push_back( Vec2( (R)i / (R)n, p[i]) );
    }

    debugl(1, "\n");

    /* i == 0 here, and the first point (0.0, p[0]) has been inserted twice. pop it */
    cvhull.pop_back();
}


/* bezier clipping. scale input polynomial to [0, 1], represent with Bernstein basis
 * and iteratively construct convex hull, intersect with t-axis, find new interval,
 * bisect or not (splitting is done with de-Casteljau algorithm) until desired accuracy tol is
 * reached. order 2 for single roots.
 * might give false positive if graph of polynomial almost "touches" the t-axis (numerically or
 * tolerance too high) */

/* class to store triple (pointer to polynomial, interval bounds). a queue of such objects will be
 * worked off in BezClip_roots. the polynomial is in Bernstein basis with t in [0,1] (NOT [left, right],
 * linear transformation is applied by de-Casteljau splitting) and represents the input polynomial
 * for x in [left, right] with t in [0, 1] */
template <typename R>
struct BezClip_Triple {
    BernsteinPolynomial<R, R>  *p;
    R                           left, right;

    BezClip_Triple() {
        this->p     = NULL;
        this->left  = 0.0;
        this->right = 0.0;
    }

    BezClip_Triple(
        BernsteinPolynomial<R, R>  *p,
        R                          left,
        R                          right)
    {
        this->p     = p;
        this->left  = left;
        this->right = right;
    }
};

/* function that generates the new interval from given convex hull. makes algorithm more readable
 * indeed */

/* given the convex hull of p, scale it to [left, right], intersect with the x-axis to get
 * [new_left, new_right] if there are any intersection points. */
template <typename R = double>
void
BezClip_getNewInterval(
        std::vector<Vec2>   pcvhull,
        R const            &left,
        R const            &right,
        bool               &interval_relevant,
        R                  &new_left,
        R                  &new_right,
        bool               &pcvhull_subset_eps_strip,
        R const            &eps)
{
    using Aux::Numbers::inf;

    uint32_t    i, nisect = 0;
    Vec2        cp, cpnext;
    R           R_inf = inf<R>();
    R           tmp, x0, x1, y0, y1;
    R           xzero_min = R_inf, xzero_max = -R_inf;


    /* check for intersections, this is done using the y-coordinate (bi) only, scaling is
     * applied if intersections are found, since scaling does not change the y coordinate
     * anyway.. */
    pcvhull_subset_eps_strip = true;

    cpnext   = pcvhull[0];
    for (i = 0; i < pcvhull.size(); i++) {
        /* if i == pcvhull.size() - 1, init values to the "wrap around" pair (last, first). this is
         * done to avoid copying the body of the loop for the warp around pair.. might be unrolled
         * later */
        if(i == pcvhull.size() - 1) {
            cp      = pcvhull.back();
            cpnext  = pcvhull.front();
        }
        else {
            cp      = cpnext;
            cpnext  = pcvhull[i + 1];
        }
        //debugl(1, "\n\n i: %d, cp = (%+20.13E, %+20.13E), cpnext = (%+20.13E, %+20.13E)\n", i, cp[0], cp[1], cpnext[0], cpnext[1]);

        /* if current point of convex hull has y coordinate that is smaller than eps_strip */
        if ( std::abs(cp[1]) < eps) {
            /* update xzero_min and xzero_max if necessary */
            if (cp[0] < xzero_min) xzero_min = cp[0];
            if (cp[0] > xzero_max) xzero_max = cp[0];
        }
        /* otherwise, if both points have y coordinates > EPS in magnitude and there's an
         * intersection with the x axis */
        else if ( (cp[1] < -eps && cpnext[1] > -eps) || (cp[1] > eps && cpnext[1] < -eps) ) {
            //debugl(1, "intersection between current and next point... ");
            /* not all points with EPSilon strip of x axis */
            pcvhull_subset_eps_strip = false;

            /* perform intersection calculation and update xmin / xmax if necessary */

            /* check which point got the smaller x coordinate. this is well-defined, since x
             * coordinates are offset by 1 / n, which is ALWAYS >> eps for computable degrees */
            if ( cp[0] < cpnext[0] ) {
                x0     = cp[0];
                y0     = cp[1];

                x1     = cpnext[0];
                y1     = cpnext[1];
            }
            else {
                x0     = cpnext[0];
                y0     = cpnext[1];

                x1     = cp[0];
                y1     = cp[1];

            }

            //debugl(1, "(x0, y0) = (%+20.13E, %+20.13E), (x1, y1) = (%+20.13E, %+20.13E)\n", x0, y0, x1, y1);

            /* compute intersection x value  and update min/max zero value if necessary */
            tmp = x0 - (y0*(x1 - x0)) / (y1 - y0);

            //debugl(1, "point of intersection: %+20.13E\n", tmp);

            if (tmp < xzero_min) {
                xzero_min = tmp;
                /* cap to 0.0, might happen due to rounding error.. */
                if (xzero_min < 0.0) xzero_min = 0.0;
            }

            if (tmp > xzero_max) {
                /* cap to 1.0, might happen due to rounding error.. */
                xzero_max = tmp;
                if (xzero_max > 1.0) xzero_max = 1.0;
            }

            /* increment number of intersections */
            nisect++;
        }
        /* otherwise, one of the following cases have occured: let y = cp[1], ynext = cpnext[1]:
         *
         * 1. both y values are < EPS in magnitude, i.e. std::abs(y) < eps && std::abs(ynext) < EPS. the fact that
         * std::abs(ynext) < eps will be discovered in the next iteration and hence will not be dealt with here.. 
         *
         * 2. both y values are larger than eps or both smaller than -eps, i.e.
         *
         * (y > eps && ynext > eps) || (y < -eps && ynext < -eps) 
         *
         * there is no intersection here, moving along.. */
        else {
            //debugl("y: %+20.13E, ynext: %+20.13E. either (std::abs(y) < eps && std::abs(ynext) < eps) or ( (y > eps && ynext > eps) || (y < -eps && ynext < eps) )\n", cp[1], cpnext[1]);
            //debugl(1, "no intersection or both y values < eps\n");
        }
    }

    //debugl(1, "nisect = %d, xzero_min: %+20.13E, xzero_max: %+20.13E\n", nisect, xzero_min, xzero_max);

    if (nisect > 2) {
        throw("BezClip_getNewInterval(): nisect > 2, this must never happen.\n");
    }

    /* interval is relevant is xzero_min < INF and xzero_max > -INF. */
    /* if convex hull is contained in eps strip [0,1]x[-eps, eps], then so is the function itself.  this is indicated by
     * the bool pcvhull_subset_eps_strip */
    if (xzero_min < R_inf) {
        //debugl(1, "xzero_min: %+20.13E, xzero_max: %+20.13E\n", xzero_min, xzero_max);
        if (xzero_max == -R_inf) {
            throw("xzero_min < INF and xzero_max == -INF. must never happen\n");
        }

        new_left            = left + (xzero_min * (right - left));
        new_right           = left + (xzero_max * (right - left));
        interval_relevant   = true;
    }
    else {
        interval_relevant   = false;
    }
}

template <typename R>
void
BezClip_roots(
        BernsteinPolynomial<R, R> const    &pinput,
        R const                            &alpha,
        R const                            &beta,
        R const                            &tol,
        std::vector<RealInterval<R>>       &roots,
        R const                            &eps,
        R const                            &eps_slope)
{
    setDebugComponent(DBG_POLYSOLVERS);

    /* first, convert the input polynomial in power / monomial basis to Bernstein basis and
     * clip interval [alpha, beta] to [0, 1] with two de-Casteljau stepts, just as in the algorithm
     * below if the interval is not bisected. this is done only if alpha != 0.0 || beta != 1.0 (yes,
     * bitwise comparion on doubles). this enables us to exactly specify [0.0, 1.0] without any
     * rescaling taking place. */
    if (alpha < 0.0 || beta > 1.0 || alpha > beta) {
        throw("BezClip_roots(): given alpha / beta not sensible or not within [0, 1]");
    }

    debugl(1, "BezClip_roots(): welcome..\n");

    BernsteinPolynomial<R, R> *proot = new BernsteinPolynomial<R, R>(pinput);

    /* if not precisely [0.0, 1.0] has been specified, clip the interval to [0, 1] using
     * BernsteinPolynomial<R, R>::split(). notice that p itself is given as an argument and is changed by
     * its own split method */
    if (alpha != 0.0 || beta != 1.0) {
        debugl(1, "BezClip_roots(): alpha != 0.0 || beta != 1.0.. clipping interval to [0, 1]\n");

        /* left part is irrelevant, pass NULL, p is modified in-place */
        if (alpha != 0.0) {
            debugl(1, "BezClip_roots(): alpha = %+20.13E != 0.0\n", alpha);
            proot->split(alpha, NULL, proot);
        }
        /* right part is irrelevant, pass NULL, p is modified in-place */
        if (beta != 1.0) {
            debugl(1, "BezClip_roots(): beta  = %+20.13E != 1.0\n", beta);
            proot->split( (beta - alpha) / (1.0 - alpha), proot, NULL);
        }
    }

    /* queue to store triples(polynomial, intervall limits) */ 
    std::list<BezClip_Triple<R>>    S;
    BezClip_Triple<R>               T;
    BernsteinPolynomial<R, R>      *p, *pleft, *pright;
    R                               pmid, tol4 = tol / 4.0, tol4rel;
    R                               left, right, middle, isize, new_isize;
    R                               new_left = 0.0, new_right = 0.0;
    bool                            interval_relevant, pcvhull_subset_eps_strip, bisect;
    std::vector<Vec2>               pcvhull;

    /* insert root triple onto stack S */
    S.push_back( BezClip_Triple<R>(proot, alpha, beta) );

    /* main loop, work off stack */
    while ( !S.empty() ) {
        /* get top element of S, set variables and pop() */
        T       = S.front();
        p       = T.p;
        left    = T.left;
        right   = T.right;

        S.pop_front();

        debugl(1, "------------- new triple for interval: [%+20.13E, %+20.13E], size: %+20.13E\n", left, right, std::abs(right - left));

        /* while interval is relevant and keeps shrinking exponentially with a rate > 2 when new convex hull is 
         * intersected with the t-axis, keep going.. should that seize without convergence, break the loop and bisect
         * the interval */
        bisect  = false;
        while(1) {
            debugl(1, "\n\n------------- interval: [%+20.13E, %+20.13E], size: %+20.13E\n", left, right, std::abs(right - left));
            /* get convex hull */
            PolyAlg::computeBezierControlPolyConvexHull<R>(*p, pcvhull, 1E-10);

            /* compute new interval */
            PolyAlg::BezClip_getNewInterval(
                    pcvhull,
                    left, right,
                    interval_relevant,
                    new_left, new_right,
                    pcvhull_subset_eps_strip,
                    eps);

            /* if interval is relevant, process further */
            if (interval_relevant) {
                isize       = std::abs(right - left);
                new_isize   = std::abs(new_right - new_left);

                debugl(1, "\n\n interval [%+20.13E, %+20.13E] relevant. new interval: [%+20.13E, %+20.13E], new_isize: %+20.13E\n\n", left, right, new_left, new_right, new_isize);

                /* check if the convex hull, and hence also the function, is bounded by strip of length 2EPS around the x-axis */
                if (pcvhull_subset_eps_strip) {
                    debugl(1, "convex hull of function subset of strip [0,1]x[-EPS,EPS] => function < EPS in magnitude everywhere in [0,1]. convergence..\n");
                    debugl(1, "setting root interval:   [%+20.13E, %+20.13E]. f(left) = %+20.13E, f(right) = %+20.13E\n", new_left, new_right, pinput.eval(new_left), pinput.eval(new_right) );
                    roots.push_back( RealInterval<R>(new_left, new_right) );
                    break;
                }
                /* otherwise, check if shrink factor is > 2 and subdivide if not */
                else {
                    /* interval already small enough (i.e. converged)? if so, add new root and break */
                    if (new_isize < tol) {
                        debugl(1, "found root interval: [%+20.13E, %+20.13E]. f(left) = %+20.13E, f(right) = %+20.13E\n", new_left, new_right, pinput.eval(new_left), pinput.eval(new_right) );
                        roots.push_back( RealInterval<R>(new_left, new_right) );
                        break;
                    }
                    /* no convergence yet, restrict polynomial to new interval and check shrink factor */
                    else { 
                        /* clip the input polynomial to the interval [new_left, new_right] to avoid
                         * accumulation of roundoff errors */
                        pinput.clipToInterval(new_left, new_right, p);

                        /* commit new interval boundaries */
                        left    = new_left;
                        right   = new_right;

                        /* has interval shrunk at least a factor of two? if so, keep going, otherwise,
                         * break loop and bisect below */
                        if (new_isize > 0.5*isize) {
                            debugl(1, "interval hasn't shrunken by a factor of at least 2. subdividing..\n");
                            bisect = true;
                            break;
                        }
                        else {
                            debugl(1, "interval has shrunk by a factor of at least 2. keep shrinking..\n");
                        }
                    }
                }
            }
            /* interval irrelevant, break the loop. */
            else {
                debugl(1, "interval [%+20.13E, %+20.13E] irrelevant. moving along..\n", left, right);
                break;
            }
        }

        /* while loop has been broken. bisect if interval is relevant, still larger than tol and has
         * seized shrinking exponentially with a factor > 2 */
        if (bisect) {
            /* alloc new polys */
            pleft   = new BernsteinPolynomial<R, R>();
            pright  = new BernsteinPolynomial<R, R>();

            /* check if middle is a root */
            middle  = (left + right) / 2.0;
            pmid    = p->eval(0.5);
            if ( std::abs(pmid) < eps) {
                /* mid point is root. since we don't want to converge twice against the same root,
                 * the following approach is taken:
                 * we know the interval has size > tol, so we can generate the new intervals
                 * [left, middle - tol/4.0] and [middle + tol / 4.0, right] to prevent double
                 * convergence. notice that this might "swallow" roots inside [middle - tol/4,
                 * middle + tol/4], yet this interval is of size tol/2 and is guaranteed to contain
                 * a root, so this is more precise than the usual case, where different roots could
                 * be contained in one interval of size ~tol as well.. */
                debugl(1, "found root           [%+20.13E], function value: %+20.13E at middle of bisected interval.. \n", middle, pmid);
                roots.push_back( RealInterval<R>(middle - tol4, middle + tol4) );

                /* we need to know where to cut in [0, 1] */
                tol4rel = tol / (4.0 * (right - left));
                debugl(1, "tol4rel: %+20.13E\n", tol4rel);

                /* split p twice, once at 0.5 - tol4rel, once at 0.5 + tol4rel */
                p->split(0.5 - tol4rel, pleft , NULL);
                p->split(0.5 + tol4rel, NULL, pright);

                /* push two new intervals to consider onto the queue */
                S.push_front( BezClip_Triple<R>(pleft , left         , middle - tol4) );
                S.push_front( BezClip_Triple<R>(pright, middle + tol4, right        ) );
            }
            else {
                debugl(1, "bisecting interval: [%+20.13E, %+20.13E] and [%+20.13E, %+20.13E]\n", left, middle, middle, right);

                /* bisect interval and initialize left and right bernstein polys */
                p->split(0.5, pleft, pright);

                /* push two new intervals to consider onto stack */
                S.push_front( BezClip_Triple<R>(pright, middle, right ) );
                S.push_front( BezClip_Triple<R>(pleft , left  , middle) );
            }
        }

        /* delete old poly */
        delete p;

        debugl(1, "\n\n");
    }

    setDebugComponent(DBG_GLOBAL);
}


/* ----------------------------------------------------------------------------------------------------------------- *
 *
 *             root finding for bivariate polynomials: bivariate linear clipping and required auxiliary algorithms    
 *
 * ----------------------------------------------------------------------------------------------------------------- */

/* struct to store polynomials and rectanglular domain in [0,1]^2. this is the data strcture the mail loop
 * operates on, a queue of such structures is worked off (and in generalnew queue elements are generated along
 * the way) */
template <typename R = double>
struct BiLinClip_Tuple {
    BiBernsteinPolynomial<R, R>    *p, *q;
    R                               alpha0, alpha1, beta0, beta1;
    bool                            alpha_converged, beta_converged;
    uint32_t                        depth;

    BiLinClip_Tuple()
    {
        p       = q = NULL;
        depth   = 0;
        alpha0  = alpha1 = beta0 = beta1 = 0.0;
    } 

    BiLinClip_Tuple(
            BiBernsteinPolynomial<R, R>    *p,
            BiBernsteinPolynomial<R, R>    *q,
            R                               alpha0, 
            R                               alpha1, 
            bool                            alpha_converged,
            R                               beta0, 
            R                               beta1,
            bool                            beta_converged,
            uint32_t                        depth)
    {
        this->p                 = p;
        this->q                 = q;
        this->alpha0            = alpha0;
        this->alpha1            = alpha1;
        this->alpha_converged   = alpha_converged;
        this->beta0             = beta0;
        this->beta1             = beta1;
        this->beta_converged    = beta_converged;
        this->depth             = depth;
    } 
};

enum BiClinClip_errcodes {
    BLC_LINSOLVE_SINGULAR,
    BLC_LINSOLVE_SUCCESS
};

/* linear solver for 2x2 intersection matrices: gauss elimination with full pivoting, all right hand
 * sides given and manipulated in parallel. inverse of matrix, condition number and error bound on
 * solution is computed, which is used to increase robustness: we'd rather have a few false
 * positives than skipping a real root, since evaluation of resulsts is very cheap, a missed root
 * could be fatal on the other hand.. */
template <typename R>
int
BiLinClip_linsolve(
        R const         A_input[2][2],
        R const         b_input[4][2],
        R             (&s)[4][2],
        R const         eps = 1E-11)
{
    //int     i;
    bool    swap_cols, swap_rows;
    R       tmp;
    R       A[2][2];
    R       b[4][2];

    /* copy values by hand, faster for this small inputs than memcpy */
    A[0][0] = A_input[0][0];
    A[0][1] = A_input[0][1];
    A[1][0] = A_input[1][0];
    A[1][1] = A_input[1][1];

    b[0][0] = b_input[0][0];
    b[0][1] = b_input[0][1];

    b[1][0] = b_input[1][0];
    b[1][1] = b_input[1][1];

    b[2][0] = b_input[2][0];
    b[2][1] = b_input[2][1];

    b[3][0] = b_input[3][0];
    b[3][1] = b_input[3][1];

    /* determine largest element of matrix and set swap flags as necessary */
    tmp         = std::abs(A[0][0]);
    swap_rows   = false;
    swap_cols   = false;

    if (std::abs(A[0][1]) > tmp) {
        tmp         = std::abs(A[0][1]);
        swap_rows   = false;
        swap_cols   = true;
    }

    if (std::abs(A[1][0]) > tmp) {
        tmp         = std::abs(A[1][0]); 
        swap_rows   = true;
        swap_cols   = false;
    }

    if (std::abs(A[1][1]) > tmp) {
        tmp         = std::abs(A[1][1]);
        swap_rows   = true;
        swap_cols   = true;
    }

    if (swap_rows) {
        std::swap(A[0][0], A[1][0]);
        std::swap(A[0][1], A[1][1]);

        std::swap(b[0][0], b[0][1]);
        std::swap(b[1][0], b[1][1]);
        std::swap(b[2][0], b[2][1]);
        std::swap(b[3][0], b[3][1]);
    }

    /* if swap_cols is true, we need to swap variables back after solving, since column interchange
     * amounts to reordering variables */
    if (swap_cols) {
        std::swap(A[0][0], A[0][1]);
        std::swap(A[1][0], A[1][1]);
    }

    /* check if A[0][0], which now contains the element biggest in magnitude, is < eps */
    if ( std::abs(A[0][0]) < eps) {
        debugl(0, "BiLinClip_linsolve(): std::abs(Amax) < eps => system singular.\n");
        return BLC_LINSOLVE_SINGULAR;
    }

    /* upper triangular form */
    tmp     = -A[1][0] / A[0][0];
    A[1][0] += tmp*A[0][0];
    A[1][1] += tmp*A[0][1];

    b[0][1] += tmp*b[0][0];
    b[1][1] += tmp*b[1][0];
    b[2][1] += tmp*b[2][0];
    b[3][1] += tmp*b[3][0];

    if (std::abs(A[1][1]) < eps) {
        debugl(1, "BiLinClip_linsolve(): std::abs(A[1][1]) < EPS after elimination. matrix singular\n");
        return BLC_LINSOLVE_SINGULAR;
    }

    /* compute solutions to intersection system */
    s[0][1]    = b[0][1] / A[1][1];
    s[0][0]    = (b[0][0] - A[0][1] * (s[0][1])) / A[0][0];

    s[1][1]    = b[1][1] / A[1][1];
    s[1][0]    = (b[1][0] - A[0][1] * (s[1][1])) / A[0][0];

    s[2][1]    = b[2][1] / A[1][1];
    s[2][0]    = (b[2][0] - A[0][1] * (s[2][1])) / A[0][0];

    s[3][1]    = b[3][1] / A[1][1];
    s[3][0]    = (b[3][0] - A[0][1] * (s[3][1])) / A[0][0];

    /* swap solutions if necessary */
    if (swap_cols) {
        std::swap(s[0][0], s[0][1]);
        std::swap(s[1][0], s[1][1]);
        std::swap(s[2][0], s[2][1]);
        std::swap(s[3][0], s[3][1]);
    }

    return BLC_LINSOLVE_SUCCESS;
}


/* generate new rectangle for bivariate linear clipping. input are the polyomials p and q
 * and the current rectangle. 
 * get best linear approximants for p and q (or rather, their coefficients in the
 * legendre basis L(m, n), as well as the bounds delta_p and delta_q for the fat lines
 * in the xy-plane. compute the four points of intersection, found the smallest axis
 * aligned bounding rectangle for the resulting parallelogram and intersect with current
 * rectangle [alpha0, alpha1]x[beta0, beta1], that is to say with [0,1]x[0,1], since p and
 * q are always in Bernstein-Bezier form with respect to that interval, but rescaled to
 * the unit square, as usual. If the intersection is empty, the current rectangle is irrelevant,
 * because no root can lie inside the current domain. */
template <typename R>
void
BiLinClip_getNewRectangle(
        const BiBernsteinPolynomial<R, R>  &p,
        const BiBernsteinPolynomial<R, R>  &q,
        uint32_t                            m,
        uint32_t                            n,
        const Matrix<R>                    &BiLinClip_A00,
        const Matrix<R>                    &BiLinClip_A01,
        const Matrix<R>                    &BiLinClip_A10,
        const BiBernsteinPolynomial<R, R>  &BiLinClip_L00,
        const BiBernsteinPolynomial<R, R>  &BiLinClip_L01,
        const BiBernsteinPolynomial<R, R>  &BiLinClip_L10,
        R const                            &alpha0,
        R const                            &alpha1,
        bool                                alpha_frozen,
        R const                             beta0,
        R const                             beta1,
        bool                                beta_frozen,
        bool                               &rectangle_relevant,
        R                                  &new_alpha0,
        R                                  &new_alpha1,
        R                                  &new_beta0,
        R                                  &new_beta1,
        R const                            &eps,
        R const                            &linsolve_eps)
{
    using Aux::Numbers::inf;

    int         err;
    uint32_t    i, j;

    R           R_inf = inf<R>();
    R           tmp, offset;
    R           p_coeff, q_coeff;
    R           delta_p, delta_q, dalpha, dbeta;
    R           p_l00, p_l01, p_l10;
    R           q_l00, q_l01, q_l10;
    R           A[2][2] = { { 0 } };
    R           s[4][2] = { { 0 } }, b[4][2] = { { 0 } };
    R           s_xmin, s_xmax, s_ymin, s_ymax;

    int         free_axis;
    char        free_axis_name[16];
    R           p_newmin, p_newmax, p_newsize;
    R           q_newmin, q_newmax, q_newsize;
    bool        p_singular, q_singular;

    p_newmin = 0.0;
    p_newmax = 1.0;
    q_newmin = 0.0;
    q_newmax = 1.0;

    /* first, get the best linear fit legendre basis coefficients {p|q}_l{00|01|10} and bounds delta_{p|q} */
    /*
    BiLinClip_computeBestLinearApprox(
            p, q, m, n, 
            p_l00, p_l10, p_l01, delta_p,
            q_l00, q_l10, q_l01, delta_q
        );
    */
    /* ----------------- start of inlined function BiLinClip_computeBestLinearApprox ---------------------- */

    /* get coefficients for L(0, 0), L(1,0) and L(0,1) from the precomputed approximant matrices */
    p_l00 = 0.0;
    p_l10 = 0.0;
    p_l01 = 0.0;

    q_l00 = 0.0;
    q_l10 = 0.0;
    q_l01 = 0.0;

    for (i = 0; i < m + 1; i++) {
        for (j = 0; j < n + 1; j++) {
            p_coeff         = p(i, j);
            p_l00          += p_coeff * BiLinClip_A00(i, j);
            p_l10          += p_coeff * BiLinClip_A10(i, j);
            p_l01          += p_coeff * BiLinClip_A01(i, j);

            q_coeff         = q(i, j);
            q_l00          += q_coeff * BiLinClip_A00(i, j);
            q_l10          += q_coeff * BiLinClip_A10(i, j);
            q_l01          += q_coeff * BiLinClip_A01(i, j);
        }
    }

    /* scale with 1 / <L(i, j), L(i, j)> = (2i+1)(2j+1) */
    /* l00_coeff's facotor is 1, since (2*0 + 1)(2*0 + 1) = 1 */
    p_l10  *= 3.0;
    p_l01  *= 3.0;

    q_l10  *= 3.0;
    q_l01  *= 3.0;

    /* construct best linear approximant by orthogonal projection, but since legendre polynomials
     * are not orthonormal, their square norm <lij, lij> is taken into account above, which is just
     * 1/(2i+1)(2j+1), i.e. multiply with (2i+1)(2j+1) */
    BiBernsteinPolynomial<R, R> p_linfit = 
        BiLinClip_L00 * p_l00 +
        BiLinClip_L10 * p_l10 +
        BiLinClip_L01 * p_l01;

    BiBernsteinPolynomial<R, R> q_linfit = 
        BiLinClip_L00 * q_l00 +
        BiLinClip_L10 * q_l10 +
        BiLinClip_L01 * q_l01;

    /* compute the bound delta = max(i, j){|p(i, j) - linfit(i, j)|} */
    delta_p = 0.0;
    delta_q = 0.0;

    for (i = 0; i < m + 1; i++) {
        for (j = 0; j < n + 1; j++) {
            tmp = std::abs( p(i, j) - p_linfit(i, j) );
            if (tmp > delta_p) {
                delta_p = tmp;
            }

            tmp = std::abs( q(i, j) - q_linfit(i, j) );
            if (tmp > delta_q) {
                delta_q = tmp;
            }
        }
    } 

    /* ----------------- end of inlined function BiLinClip_computeBestLinearApprox ---------------------- */

    debugl(1, "BiLinClip_getNewRectangle: p_l00 = %20.13E, p_l10 = %20.13E, p_l01 = %20.13E, delta_p = %20.13E\n", p_l00, p_l10, p_l01, delta_p);
    debugl(1, "BiLinClip_getNewRectangle: q_l00 = %20.13E, q_l10 = %20.13E, q_l01 = %20.13E, delta_q = %20.13E\n", q_l00, q_l10, q_l01, delta_q);

    /* if neither alpha nor beta are frozen, try to solve intersection system for two fat lines.
     * that system can be singular because
     *
     * 1. at least one function (maybe both) have a best linear approximating plane that is parallel
     * to the xy-plane -> there is no intersection with the xy plane and hence there is no fat line.
     *
     * 2. the two intersection lines are parallel. this should not happen in general if the input
     * functions are different (but will of course happen if the input functions are identical),
     * and we have to deal with this in some non-obvious way. */
    if (!alpha_frozen && !beta_frozen) {
        /* fill matrix A */
        A[0][0] = 2.0 * p_l10; 
        A[0][1] = 2.0 * p_l01; 

        A[1][0] = 2.0 * q_l10; 
        A[1][1] = 2.0 * q_l01; 

        /* fill inhomogeneities b */
        b[0][0] = p_l10 + p_l01 - p_l00 + delta_p;
        b[0][1] = q_l10 + q_l01 - q_l00 + delta_q;

        b[1][0] = p_l10 + p_l01 - p_l00 + delta_p;
        b[1][1] = q_l10 + q_l01 - q_l00 - delta_q;

        b[2][0] = p_l10 + p_l01 - p_l00 - delta_p;
        b[2][1] = q_l10 + q_l01 - q_l00 + delta_q;

        b[3][0] = p_l10 + p_l01 - p_l00 - delta_p;
        b[3][1] = q_l10 + q_l01 - q_l00 - delta_q;

        /* call linear solver */
        err = BiLinClip_linsolve(A, b, s, linsolve_eps);

        /* print */
        debugl(1, "solutions back in BiLinClip_getNewRectangle(): \n");
        for (i = 0; i < 4; i++) {
            debugl(1, "s(%d): (%f %f)\n", i, s[i][0], s[i][1]);
        }

        /* solver returned success */
        if (err == BLC_LINSOLVE_SUCCESS) {
            /* determine minimum bounding box in [0,1] by finding minimum / maximum x / y coordinate */
            s_xmin = s_xmax = s[0][0];
            s_ymin = s_ymax = s[0][1];

            for (i = 1; i < 4; i++) {
                if (s[i][0] < s_xmin) {
                    s_xmin = s[i][0];
                }

                if (s[i][0] > s_xmax) {
                    s_xmax = s[i][0];
                }

                if (s[i][1] < s_ymin) {
                    s_ymin = s[i][1];
                }

                if (s[i][1] > s_ymax) {
                    s_ymax = s[i][1];
                }
            }

            debugl(1, "bounding box of parallelogram:                          [%10.5f (%20.13E), %10.5f (%20.13E)]x[%10.5f (%20.13E), %10.5f (%20.13E)]\n", s_xmin, s_xmin, s_xmax, s_xmax, s_ymin, s_ymin, s_ymax, s_ymax);

            /* intersection is empty iff x interval doesnt intersect [0,1] or y intervals doesnt intersect
             * [0,1]. */
            if (s_xmax < 0.0 || s_xmin > 1.0 || s_ymax < 0.0 || s_ymin > 1.0) {
                debugl(1, "setting rectangle to irrelevant.\n");
                rectangle_relevant = false;
                return;
            }
            /* rectangle is relevant <=> there is an intersection, which can be written as
             * [max(0, s_xmin), min(1, xsmax)]x[max(0,s_ymin),min(1,s_ymax)] */
            else {
                rectangle_relevant = true;
                
                /* calculate absolute values for new_alpha0. etc directly here */
                if (s_xmin > 0.0) {
                    new_alpha0  = alpha0 + (s_xmin * (alpha1 - alpha0));
                }
                else {
                    s_xmin      = 0.0;
                    new_alpha0  = alpha0;
                }

                if (s_xmax < 1.0) {
                    new_alpha1  = alpha0 + (s_xmax * (alpha1 - alpha0));
                }
                else {
                    s_xmax      = 1.0;
                    new_alpha1  = alpha1;
                }

                if (s_ymin > 0.0) {
                    new_beta0   = beta0 + (s_ymin * (beta1 - beta0) );
                }
                else {
                    s_ymin      = 0.0;
                    new_beta0   = beta0;
                }

                if (s_ymax < 1.0) {
                    new_beta1   = beta0 + (s_ymax * (beta1 - beta0) );
                }
                else {
                    s_ymax      = 1.0;
                    new_beta1   = beta1;
                }

                debugl(1, "relative: intersection of bounding box with [0,1]:      [%10.5f (%20.13E), %10.5f (%20.13E)]x[%10.5f (%20.13E), %10.5f (%20.13E)]\n", s_xmin, s_xmin, s_xmax, s_xmax, s_ymin, s_ymin, s_ymax, s_ymax);
                debugl(1, "absolute: new rectangle in absolute coordinates [0,1]:  [%10.5f (%20.13E), %10.5f (%20.13E)]x[%10.5f (%20.13E), %10.5f (%20.13E)]\n", new_alpha0, new_alpha0, new_alpha1, new_alpha1, new_beta0, new_beta0, new_beta1, new_beta1);
                debugl(1, "\n");
                return;
            }
        }
        else if (err == BLC_LINSOLVE_SINGULAR) {
            debugl(0, "BiLinClip_getNewRectangle(): intersection system for fat lines singular. checking if p or q best linear fit is parallel to xy plane...\n");

            tmp = std::max( std::abs(A[0][0]), std::abs(A[0][1]) );
            if (tmp < linsolve_eps) {
                debugl(0, "p plane nearly parallel to xy plane => checking if bounding planes are both above or below xy-plane..\n");
                offset = p_l00 - p_l10 - p_l01;
                if ( offset + delta_p < -eps) {
                    debugl(0, "p + delta_p has offset %f < -1E-14 => interval irrelevant\n", offset + delta_p);
                    rectangle_relevant = false;
                    return;
                }
                else if ( offset - delta_p > eps) {
                    debugl(0, "p - delta_p has offset %f > 1E-14 => interval irrelevant\n", offset - delta_p);
                    rectangle_relevant = false;
                    return;
                }
                else {
                    debugl(0, "p bounding planes probably sandwhich xy plane, offset: (%f, %f) => assuming relevant. leaving alpha and beta interval as is => subdivision..\n", offset - delta_q, offset + delta_p);
                    rectangle_relevant  = true;
                    new_alpha0          = alpha0;
                    new_alpha1          = alpha1;
                    new_beta0           = beta0;
                    new_beta1           = beta1;
                    return;
                }
            }

            tmp = std::max( std::abs(A[1][0]), std::abs(A[1][1]) );
            if (tmp < linsolve_eps) {
                debugl(0, "q plane nearly parallel to xy plane => checking if bounding planes are both above or below xy-plane..\n");
                offset = q_l00 - q_l10 - q_l01;
                if ( offset + delta_q < -eps) {
                    debugl(0, "q + delta_q has offset %f < -1E-14 => interval irrelevant\n", offset + delta_q);
                    rectangle_relevant = false;
                    return;
                }
                else if ( offset - delta_q > eps) {
                    debugl(0, "q - delta_q has offset %f > 1E-14 => interval irrelevant\n", offset - delta_q);
                    rectangle_relevant = false;
                    return;
                }
                else {
                    rectangle_relevant  = true;
                    new_alpha0          = alpha0;
                    new_alpha1          = alpha1;
                    new_beta0           = beta0;
                    new_beta1           = beta1;
                    debugl(0, "q bounding planes probably sandwhich xy plane, offset: (%f, %f) => assuming relevant. leaving alpha and beta interval as is => subdivision..\n", offset - delta_q, offset + delta_q);
                    return;
                }
            }

            /* if we ever reach this line, both rows of the matrix have elements > EPS, yet still
             * the system became singular => panic for now */
            throw("BiLinClip_getNewRectangle(): neither row small, yet intersection system singular... parallel fat lines?");
        }
        else {
            throw("BiLinClip_getNewRectangle(): invalid error code returned by linear solver.");
        }
    }
    /* either alpha or beta, but not both, are frozen: shrink rectangle for non-frozen variable only */
    else if ( (!alpha_frozen && beta_frozen) || (alpha_frozen && !beta_frozen) ) {
        /* intersect the two fat lines with the boundary lines y = 0 and y = 1 if alpha is frozen and with boundary
         * lines x = 0 and x = 1 if beta is frozen.
         *
         * after intersecting one fat line, we get two points from each boundary line. the intersection parallelogram of
         * both fat lines is definitely bounded by the bounding rectangles arising from intersecting each of the fat
         * lines with either of the two boundary lines of [0,1]^2.
         *
         * this is done for BOTH fat lines, and the smaller rectangle with respect to the unfrozen variable is taken. 
         *
         * it can happen that one of the systems solved is singular, because either the fat lines are numerically
         * parallel to the lines x-axis (y = 0, y = 1), or because there is no fat line, because the best linear
         * approximant does not intersect the xy-plane (parallel to it) */

        p_singular  = false;
        q_singular  = false;

        /* set free_axis to 0 (for x) if beta is frozen and to 1 (for y) if alpha is frozen */
        if (alpha_frozen) {
            free_axis = 1;
            sprintf(free_axis_name, "beta ");
        }
        else {
            free_axis = 0;
            sprintf(free_axis_name, "alpha");
        }

        /* intersect fat line of p with boundaries */

        /* fill matrix, first row for p, second row for parallel lines to x-axis */
        A[0][0] = 2*p_l10; 
        A[0][1] = 2*p_l01; 

        /* if alpha is frozen, we wish to intersect fat lines with x = 0 and x = 1, if beta is frozen,
         * intersect with y = 0 and y = 1 */
        if (alpha_frozen) {
            A[1][0] = 1.0;
            A[1][1] = 0.0;
        }
        else {
            A[1][0] = 0.0;
            A[1][1] = 1.0;
        }

        /* inhomogeneities define which part of the fat line we intersect (+-delta) with which
         * boundary (x/y = 0 or x/y = 1) */
        b[0][0] = p_l10 + p_l01 - p_l00 - delta_p;
        b[0][1] = 0.0;

        b[1][0] = p_l10 + p_l01 - p_l00 + delta_p;
        b[1][1] = 0.0;

        b[2][0] = p_l10 + p_l01 - p_l00 - delta_p;
        b[2][1] = 1.0;

        b[3][0] = p_l10 + p_l01 - p_l00 + delta_p;
        b[3][1] = 1.0;

        /* call linear solver */
        err = BiLinClip_linsolve(A, b, s);

        /* if system is singular, set new size to infinity */
        if (err == BLC_LINSOLVE_SINGULAR) {
            debugl(0, "system for p singular..\n");
            p_singular  = true;
            p_newsize   = R_inf;
        }
        else {
            /* system successfully solved. get min / max values for alpha/beta */
            p_newmin    = s[0][free_axis];

            if (s[1][free_axis] < p_newmin) p_newmin = s[1][free_axis];
            if (s[2][free_axis] < p_newmin) p_newmin = s[2][free_axis];
            if (s[3][free_axis] < p_newmin) p_newmin = s[3][free_axis];

            /* subtract 1% from p_newmin */
            p_newmin   -= 0.01 * p_newmin;

            /* compute maximum alpha/beta value */
            p_newmax    = s[0][free_axis];

            if (s[1][free_axis] > p_newmax) p_newmax = s[1][free_axis];
            if (s[2][free_axis] > p_newmax) p_newmax = s[2][free_axis];
            if (s[3][free_axis] > p_newmax) p_newmax = s[3][free_axis];

            /* add 1% from p_newmin */
            p_newmax   += 0.01 * p_newmax;

            /* compute new interval size */
            p_newsize   = p_newmax - p_newmin;

            debugl(1, "p fat line: new %s interval              [%f (%20.13E), %f (%20.13E)], size: %f (%20.13E)\n",
                    free_axis_name, p_newmin, p_newmin, p_newmax, p_newmax, p_newsize, p_newsize);

            /* determine if there is an intersectino with [0,1] */
            if (p_newmin > 1.0 || p_newmax < 0.0) {
                debugl(1, "new %s interval does not intersect [0,1] -> discard rectangle.\n", free_axis_name);
                rectangle_relevant = false;
                return;
            }
            /* there is an intersection, perform intersection and set newsize variable */
            else {
                rectangle_relevant = true;
                debugl(1, "new %s interval intersects [0,1] -> rectangle potentially relevant.\n", free_axis_name);
                if (p_newmin < 0.0) p_newmin = 0.0;
                if (p_newmax > 1.0) p_newmax = 1.0;

                p_newsize = p_newmax - p_newmin;

                debugl(1, "p fat line: intersection with [0,1]:        [%f (%20.13E), %f (%20.13E)], size: %f (%20.13E)\n",
                        p_newmin, p_newmin, p_newmax, p_newmax, p_newsize, p_newsize);
            }

            debugl(1, "p fat line: new %s interval              [%f (%20.13E), %f (%20.13E)], size: %f (%20.13E)\n",
                    free_axis_name, p_newmin, p_newmin, p_newmax, p_newmax, p_newsize, p_newsize);

        }

        /* same for q */

        /* fill matrix, first row for p, second row for parallel lines to x-axis */
        A[0][0] = 2*q_l10; 
        A[0][1] = 2*q_l01; 

        if (alpha_frozen) {
            A[1][0] = 1.0;
            A[1][1] = 0.0;
        }
        else {
            A[1][0] = 0.0;
            A[1][1] = 1.0;
        }

        b[0][0] = q_l10 + q_l01 - q_l00 - delta_q;
        b[0][1] = 0.0;

        b[1][0] = q_l10 + q_l01 - q_l00 + delta_q;
        b[1][1] = 0.0;

        b[2][0] = q_l10 + q_l01 - q_l00 - delta_q;
        b[2][1] = 1.0;

        b[3][0] = q_l10 + q_l01 - q_l00 + delta_q;
        b[3][1] = 1.0;

        /* call linear solver */
        err = BiLinClip_linsolve(A, b, s);

        /* if system is singular, set new size to infinity */
        if (err == BLC_LINSOLVE_SINGULAR) {
            debugl(0, "system for q singular..\n");
            q_singular  = true;
            q_newsize   = R_inf;
        }
        else {
            /* system successfully solved. get min / max values for alpha/beta */
            q_newmin    = s[0][free_axis];

            if (s[1][free_axis] < q_newmin) q_newmin = s[1][free_axis];
            if (s[2][free_axis] < q_newmin) q_newmin = s[2][free_axis];
            if (s[3][free_axis] < q_newmin) q_newmin = s[3][free_axis];

            q_newmin   -= 0.01 * q_newmin;

            /* compute maximum alpha/beta value */
            q_newmax    = s[0][free_axis];

            if (s[1][free_axis] > q_newmax) q_newmax = s[1][free_axis];
            if (s[2][free_axis] > q_newmax) q_newmax = s[2][free_axis];
            if (s[3][free_axis] > q_newmax) q_newmax = s[3][free_axis];

            /* add max_err to maximum value. this widens the interval and includes error analysis */
            q_newmax   += 0.01 * q_newmax;

            /* compute new interval size */
            q_newsize   = q_newmax - q_newmin;

            debugl(1, "q fat line: new %s interval              [%f (%20.13E), %f (%20.13E)], size: %f (%20.13E)\n",
                    free_axis_name, q_newmin, q_newmin, q_newmax, q_newmax, q_newsize, q_newsize);

            /* determine if there is an intersectino with [0,1] */
            if (q_newmin > 1.0 || q_newmax < 0.0) {
                debugl(1, "new %s interval does not intersect [0,1] -> discard rectangle.\n", free_axis_name);
                rectangle_relevant = false;
                return;
            }
            /* there is an intersection, perform intersection and set newsize variable */
            else {
                rectangle_relevant = true;
                debugl(1, "new %s interval intersectds [0,1] -> rectangle potentially relevant.\n", free_axis_name);
                if (q_newmin < 0.0) q_newmin = 0.0;
                if (q_newmax > 1.0) q_newmax = 1.0;

                q_newsize = q_newmax - q_newmin;

                debugl(1, "q fat line: intersection with [0,1]:        [%f (%20.13E), %f (%20.13E)], size: %f (%20.13E)\n",
                        q_newmin, q_newmin, q_newmax, q_newmax, q_newsize, q_newsize);
            }

            debugl(1, "q fat line: new %s interval              [%f (%20.13E), %f (%20.13E)], size: %f (%20.13E)\n",
                    free_axis_name, q_newmin, q_newmin, q_newmax, q_newmax, q_newsize, q_newsize);

        }

        debugl(1, "\n");

        /* compute new interval if possible */
        if (!p_singular || !q_singular) {
            if (alpha_frozen) {
                /* new values for beta */
                dbeta = beta1 - beta0;

                /* take the better of the two */
                if (p_newsize < q_newsize) {
                    debugl(1, "comitting new beta interval from p fat line..\n");
                    new_beta0  = beta0 + (p_newmin * dbeta);
                    new_beta1  = beta0 + (p_newmax * dbeta);
                }
                else {
                    debugl(1, "comitting new beta interval from q fat line..\n");
                    new_beta0  = beta0 + (q_newmin * dbeta);
                    new_beta1  = beta0 + (q_newmax * dbeta);
                }
            }
            else {
                /* new values for alpha */
                dalpha = alpha1 - alpha0;

                /* take the better of the two */
                if (p_newsize < q_newsize) {
                    debugl(1, "comitting new alpha interval from p fat line..\n");
                    new_alpha0  = alpha0 + (p_newmin * dalpha);
                    new_alpha1  = alpha0 + (p_newmax * dalpha);
                }
                else {
                    debugl(1, "comitting new alpha interval from q fat line..\n");
                    new_alpha0  = alpha0 + (q_newmin * dalpha);
                    new_alpha1  = alpha0 + (q_newmax * dalpha);
                }
            }
        }
        /* if both systems were singular, throw() exception. */
        else {
            throw("BiLinClip_getNewRectangle(): one axis frozen and both intersection systems for new free axis interval singular..");
        }
    }
    /* alpha and beta may never both be frozen */
    else {
        throw("BiLinClip_getNewRectangle(): both alpha and beta are frozen. internal logic error.");
    }
}

#define BILCLIP_APPROXDATA_DEFAULT_SIZE 8

template <typename R>
void
BiLinClip_getApproximationData(
    uint32_t                                    m,
    uint32_t                                    n,
    bool                                        dynamic_recomputation_arg,
    BiBernsteinPolynomial<R, R>                *BiLinClip_L00,
    BiBernsteinPolynomial<R, R>                *BiLinClip_L10,
    BiBernsteinPolynomial<R, R>                *BiLinClip_L01,
    Matrix<R>                                  *BiLinClip_A00,
    Matrix<R>                                  *BiLinClip_A10,
    Matrix<R>                                  *BiLinClip_A01)
{
    using Aux::Numbers::bicof;

    /* static variables to store precomputable legendre polynomials and approximation matrices. */
    static uint32_t                                 m_max = 0, n_max = 0;
    static bool                                     dynamic_recomputation = true;

    static Tensor<2, BernsteinPolynomial<R, R>>     LegendrePolBB;
    static Tensor<4, BiBernsteinPolynomial<R, R>>   LegendreBiPolBB;
    static Tensor<6, R>                             LegendreBiApproximantMatrices;

    /* recompute if necessary and dynamic_recomputation == true */
    if (m > m_max || n > n_max) {
        /* recomputation necessary */
        if (!dynamic_recomputation) {
            throw("BiLinClip_getApproximationData(): recomputation necessary but dynamic_recomputation flag set to false.");
        }

        debugl(0, "BiLinClip_getApproximationData(): recomputing approximation data for pair (m, n) = (%d, %d).\n", m, n);

        /* every degree is lower bounded by macro BILCLIP_APPROXDATA_DEFAULT_SIZE, get maximum of this lower bound and
         * the values m, n and the current m_max, n_max, respectively. */
        m_max = std::max(m_max, std::max(m, (uint32_t)BILCLIP_APPROXDATA_DEFAULT_SIZE));
        n_max = std::max(n_max, std::max(n, (uint32_t)BILCLIP_APPROXDATA_DEFAULT_SIZE));

        debugl(1, "initLegendrePolynomialsBB(): initializing univariate and bivariate shifted Legendre polynomials in Bernstein-Bezier basis... \n");
        fflush(stdout);

        /* resize 2d tensor for univariate legendre basis polynomials, size: (n_max + 1, n_max + 1) */
        LegendrePolBB.resize({ std::max(m_max, n_max) + 1, std::max(m_max, n_max) + 1});

        Vector<R>   lnk_coeff;
        R           sign;

        /* init univariate legendre polynomials in BB(n) for all degrees n = 0..max(n_max, m_max) */
        for (int n = 0; (uint32_t)n < std::max(m_max, n_max) + 1; n++) {
            //debugl("----------\n");
            lnk_coeff.assign(n + 1, 0.0);

            /* for all polynomials L_k, k = 0..n */
            for (int k = 0; k < n + 1; k++) {

                /* for all components c_j of L_k in BB(n), j = 0..n */
                for (int j = 0; j < n + 1; j++) {
                    lnk_coeff[j] = 0.0;

                    for (int i = std::max(0, j - n + k); i < std::min(j, k) + 1; i++) {
                        sign            = ( (k + i) % 2 == 0) ? 1.0 : -1.0;
                        lnk_coeff[j]   += sign * bicof<R>(k, i) * bicof<R>(k, i) * bicof<R>(n - k, j - i);
                    }

                    /* divide sum by 1 / (n choose j) */
                    lnk_coeff[j] /= bicof<R>(n, j);
                }

                /* assign coefficients for polynomial L_k in BB(n) */
                LegendrePolBB( {(uint32_t)n, (uint32_t)k} ).setCoeffs(lnk_coeff);

                //debugl("BB(%d): coefficients of legendre polynomial L_(%d): ", n, k), LegendrePolBB[n][k].printCoeff();
            }
        }

        /* now init bivariate legendre polynomials in BB(m, n), m = 0..m_max, n = 0..n_max using tensor product of the
         * univariate ones. */

        /* resize 4d tensor for bivariate legendre basis polynomials: size (m_max + 1, n_max + 1, m_max + 1, n_max + 1) */
        LegendreBiPolBB.resize({m_max + 1, n_max + 1, m_max + 1, n_max + 1});

        for (uint32_t m = 0; m < m_max + 1; m++) {
            for (uint32_t n = 0; n < n_max + 1; n++) {
                for (uint32_t i = 0; i < m + 1; i++) {
                    for (uint32_t j = 0; j < n + 1; j++) {
                        LegendreBiPolBB( {m, n, i, j} ) =
                            PolyAlg::BernsteinTensorMultiply(
                                LegendrePolBB( {m, i} ),
                                LegendrePolBB( {n, j} )
                            );
                    }
                }
            }
        }

        /* now for the legendre approximation matrices */

        /* resize to correct dimensions and zero */
        LegendreBiApproximantMatrices.resize(
            {
                m_max + 1,
                n_max + 1,
                3,
                3,
                m_max + 1,
                n_max + 1
            } );
        LegendreBiApproximantMatrices.fill(0);

        /* construct a univariate BernsteinPolynomial of degree max(max_m, max_n) to ensure that the Bernstein basis
         * polynomial inner products are precomputed at least up to the maximum degree required here. */
        //BernsteinPolynomial<R, R>       bb(std::max(m_max, n_max), 0);

        /* bivairate legendre polynomial L(i, j) represented in BB(m, n). */
        BiBernsteinPolynomial<R, R>     L_m_n_li_lj;

        /* for m, n in ranges given by m_max, n_max .. */
        for (uint32_t m = 0; m < m_max + 1; m++) {
            for (uint32_t n = 0; n < n_max + 1; n++) {

                /* compute the matrix that is needed to get the coefficient for legenrepolynomial L(l_i,
                 * l_j) in BB(m, n). */
                for (uint32_t l_i = 0; l_i < std::min(m + 1, 3u); l_i++) {
                    for (uint32_t l_j = 0; l_j < std::min(n + 1, 3u); l_j++) {

                        /* compute L(i, j) in BB(m, n) with tensor product */
                        //L_m_n_li_lj = LegendrePolBB[m][l_i].tensorMultiply(LegendrePolBB[n][l_j]);

                        L_m_n_li_lj = PolyAlg::BernsteinTensorMultiply(
                                LegendrePolBB( {m, l_i} ),
                                LegendrePolBB( {n, l_j} )
                            );

                        /* compute element (i, j) of approximation matix for L(l_i, l_j) in BB(m, n) */
                        for (uint32_t i = 0; i < m + 1; i++) {
                            for (uint32_t j = 0; j < n + 1; j++) {
                                //debugl("computing element LegendreBiApproximantMatrices[%d][%d][%d][%d][%d][%d]\n", m, n, l_i, l_j, i, j);

                                LegendreBiApproximantMatrices( {m, n, l_i, l_j, i, j} ) = 0.0;

                                for (uint32_t k = 0; k < m + 1; k++) {
                                    for (uint32_t l = 0; l < n + 1; l++) {
                                        /*
                                        LegendreBiApproximantMatrices[m][n][l_i][l_j][i][j] +=
                                            L_m_n_li_lj(k, l) * BBInnerProducts[m][m][i][k] * BBInnerProducts[n][n][j][l];
                                        */
                                        LegendreBiApproximantMatrices( {m, n, l_i, l_j, i, j} ) +=
                                            L_m_n_li_lj(k, l) *
                                            BernsteinPolynomial<R, R>::getBernsteinBasisInnerProduct(m, m, i, k) *
                                            BernsteinPolynomial<R, R>::getBernsteinBasisInnerProduct(n, n, j, l);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        debugl(0, "BiLinClip_getApproximationData(): approximation data computed for pair (m, n) = (%d, %d).\n", m, n);
    }

    /* if no exception has been thrown due to disabled dynamic recomputation, all approximation data has already been
     * available or recomputed if necessary.  write values desired by the caller. */

    //BiLinClip_L00 = LegendreBiPolBB[m][n][0][0];
    //BiLinClip_L10 = LegendreBiPolBB[m][n][1][0];
    //BiLinClip_L01 = LegendreBiPolBB[m][n][0][1];

    if (BiLinClip_L00) *BiLinClip_L00 = LegendreBiPolBB( {m, n, 0, 0} );
    if (BiLinClip_L10) *BiLinClip_L10 = LegendreBiPolBB( {m, n, 1, 0} );
    if (BiLinClip_L01) *BiLinClip_L01 = LegendreBiPolBB( {m, n, 0, 1} );

    if (BiLinClip_A00) BiLinClip_A00->resize(m + 1, n + 1);
    if (BiLinClip_A10) BiLinClip_A10->resize(m + 1, n + 1);
    if (BiLinClip_A01) BiLinClip_A01->resize(m + 1, n + 1);

    for (uint32_t i = 0; i < m + 1; i++) {
        for (uint32_t j = 0; j < n + 1; j++) {
            //BiLinClip_A00(i, j) = LegendreBiApproximantMatrices[m][n][0][0][i][j];
            //BiLinClip_A10(i, j) = LegendreBiApproximantMatrices[m][n][1][0][i][j];
            //BiLinClip_A01(i, j) = LegendreBiApproximantMatrices[m][n][0][1][i][j];

            if (BiLinClip_A00) BiLinClip_A00->operator()(i, j) = LegendreBiApproximantMatrices( {m, n, 0, 0, i, j} );
            if (BiLinClip_A10) BiLinClip_A10->operator()(i, j) = LegendreBiApproximantMatrices( {m, n, 1, 0, i, j} );
            if (BiLinClip_A01) BiLinClip_A01->operator()(i, j) = LegendreBiApproximantMatrices( {m, n, 0, 1, i, j} );
        }
    }

    /* set dynamic recomputation flag from parameters after update */
    dynamic_recomputation = dynamic_recomputation_arg;

    debugl(1, "done.\n");
}

/* NOTE: original algorithm again: diameter is numerically flawed as a metric,
 * since the the product of very small numbers numbers quickly approaches EPS, and the
 * sqrt of this is just meaningless noise. it will be almost certainly wiser to consider
 *
 * max(width, height) < tol
 *
 * as a termination criterion.
 * also, we have to guarantee that the rectangles don't degenerate too much, i.e. that they
 * don't become thin "stripes" rather than being approximately square. therefore, the
 * original algorithm was altered in such a way that both x (alpha) and y (beta) can be
 * "frozen", i.e. the width / height of the interval can be fixed, and only clipping
 * with respect to the other axis is performed. a good measure for this is the "aspect ratio"
 *
 * AR(rectangle) = max(width / height, height / width)
 *
 * where width = (alpha1 - alpha) and height = (beta1 - beta0), i.e. we consider the
 * alpha interval as width, the beta interval as height.
 *
 * AR is 1 iff the rectangle is a square. a large AR implies a rectangle with sides of
 * very different lengths, which is undesireable due to the following problems:
 *
 * 1. if width or height has already reached a value smaller than the tolerance, say 
 * width < tol, it is unwise to split with respect to alpha any more. the width could
 * then quickly become very small (< EPS) due to repeated bisection of the alpha interval,
 * while the height might still be huge in comparison. therefore, the width will be 
 * fixed and only clipping with respect to the beta axis, i.e. height will be performed.
 * same holds for the symmetrical case (for beta, height < tol) obviously.
 *
 * 2. if both width and height are still larger than tol, but the rectangle has a high
 * AR value, we might fix one of axes and shrink / split only in the other axis to
 * get rectangles with AR value closer to 1 again, since the numerical condition of
 * the intersection system is then much better.
 *
 * note that the intersection 2x2 system can quickly become singular if width and/or height are
 * very small, since the best approximant will be ~0, i.e. the xy plane for double roots
 * (or roots of even multiplicity). to prevent this, as said in 1., once width / height
 * drop below the tolerance, they are frozen and never unfrozen, and only clipping in
 * the other variable is performed until either the rectangles become irrelevant or
 * converge with respect to the other axis as well.. */
template <typename R>
void
BiLinClip_roots(
    BiBernsteinPolynomial<R, R> const  &pinput,
    BiBernsteinPolynomial<R, R> const  &qinput,
    R const                            &alpha0_input,
    R const                            &alpha1_input,
    R const                            &beta0_input,
    R const                            &beta1_input,
    R const                            &tol,
    std::vector<RealRectangle<R>>      &roots,
    bool                                approx_data_dynamic_recomputation,
    bool                                use_blacklist,
    std::vector<RealRectangle<R>>      *blacklist,
    R const                            &eps,
    R const                            &linsolve_eps,
    R const                            &rec_armax)
{
    using Aux::Numbers::inf;

    setDebugComponent(DBG_POLYSOLVERS);

    /* variables */
    uint32_t            m, n, p_m, p_n, q_m, q_n;

    pinput.getDegree(p_m, p_n);
    qinput.getDegree(q_m, q_n);

    if (p_m != q_m || p_n != q_n) {
        throw("BiLiClip_roots(): input polynomials not represented in same basis BB(m, n).");
    }
    else {
        m = p_m;
        n = p_n;
    }

    /* check for gargabe input rectangle */
    if (alpha0_input < 0.0 || alpha1_input > 1.0 ||
            alpha0_input > alpha1_input ||
            beta0_input < 0.0 || beta1_input > 1.0 ||
            beta0_input > beta1_input)
    {
        throw("BiLinClip_roots(): given rectangular domain [alpha0, alpha1] x [beta0, beta1] not sensible or not within [0,1]^2");
    }
    else if ( (alpha1_input - alpha0_input) < tol && (beta1_input - beta0_input) < tol) {
        throw("BiLinClip_roots(): input rectangle width and height below tolerance.\n");
    }

    /* extract precomputed data for relevant bidegree (m, n) if it does not match last seen bidegree */
    debugl(1, "BiLiClip_roots(): (m, n) = (%d, %d). fetching precomputed approximation data..\n", m, n);

    BiBernsteinPolynomial<R, R> BiLinClip_L00;
    BiBernsteinPolynomial<R, R> BiLinClip_L10;
    BiBernsteinPolynomial<R, R> BiLinClip_L01;

    Matrix<R>                   BiLinClip_A00;
    Matrix<R>                   BiLinClip_A10;
    Matrix<R>                   BiLinClip_A01;

    BiLinClip_getApproximationData(
            m, n,
            approx_data_dynamic_recomputation,
            &BiLinClip_L00,
            &BiLinClip_L10,
            &BiLinClip_L01,
            &BiLinClip_A00,
            &BiLinClip_A10,
            &BiLinClip_A01);

    /* necessary approximation data initialized. proceed with the algorithm */

    /* copy input polynomials for root input tuple */
    BiBernsteinPolynomial<R, R>    *proot = new BiBernsteinPolynomial<R, R>(pinput);
    BiBernsteinPolynomial<R, R>    *qroot = new BiBernsteinPolynomial<R, R>(qinput);

    /* check if input domain is PRECISELY [0,1]^2 (yes, bitwise), for generally this will be the
     * case and we don't need to cut around with deCasteljau to rescale to [0,1]^2 */
    if (alpha0_input != 0.0 || alpha1_input != 1.0 || beta0_input != 0.0 || beta1_input != 1.0) {
        debugl(1, "BiLiClip_roots(): input domain not [0,1]^2. cutting with deCasteljau\n");
        proot->clipToInterval(alpha0_input, alpha1_input, beta0_input, beta1_input, proot);
        qroot->clipToInterval(alpha0_input, alpha1_input, beta0_input, beta1_input, qroot);
    }
    else {
        debugl(1, "BiLiClip_roots(): input domain [0,1]^2: fine..\n");
    }     

    std::list<BiLinClip_Tuple<R>>   S;
    BiLinClip_Tuple<R>              T;
    BiBernsteinPolynomial<R, R>    *p, *q;
    BiBernsteinPolynomial<R, R>    *pleft_down, *pright_down, *pright_up, *pleft_up;
    BiBernsteinPolynomial<R, R>    *qleft_down, *qright_down, *qright_up, *qleft_up;
    BiBernsteinPolynomial<R, R>    *pup, *pdown, *pleft, *pright;
    BiBernsteinPolynomial<R, R>    *qup, *qdown, *qleft, *qright;

    RealRectangle<R>                current_rectangle;
    R                               alpha0, alpha1, beta0, beta1, dalpha, dbeta, alpha_middle, beta_middle;
    R                               new_alpha0, new_alpha1, new_dalpha, new_beta0, new_beta1, new_dbeta;
    R                               rec_width, rec_height, rec_ar;
    //R                               pval, qval;
    bool                            alpha_frozen, alpha_converged, beta_frozen, beta_converged;
    bool                            subdivide, rec_relevant, rec_blacklisted;
    uint32_t                        depth;

    /* push initial tuple with initial array and copies of input polynomials onto stack */
    S.push_front(
            BiLinClip_Tuple<R>(
                proot, qroot,
                alpha0_input, alpha1_input, (alpha1_input - alpha0_input) < tol,
                beta0_input,  beta1_input,  (beta1_input  - beta0_input)  < tol,
                0
            )
    );

    /* pre-loop variable init */
    rec_relevant    = false;
    new_alpha0      = 0.0;
    new_alpha1      = 1.0;
    new_beta0       = 0.0;
    new_beta1       = 1.0;

    while( !S.empty() ) {
        /* get front element of Q, set variables and pop() */
        T               = S.back();
        p               = T.p;
        q               = T.q;
        alpha0          = T.alpha0;
        alpha1          = T.alpha1;
        alpha_converged = T.alpha_converged;
        alpha_frozen    = T.alpha_converged;
        beta0           = T.beta0;
        beta1           = T.beta1;
        beta_converged  = T.beta_converged;
        beta_frozen     = T.beta_converged;
        depth           = T.depth;

        S.pop_back();

        /* reset flags */
        subdivide       = false;
        rec_blacklisted = false;

        /* while rectangle stays relevant, hasn't converged and keeps shrinking in area by a factor > 4 */
        while (1) {
            /* set current_rectangle */
            current_rectangle = RealRectangle<R>(alpha0, alpha1, beta0, beta1);

            /* if use_blacklist is true, check if rectangle is contained in a blacklisted
             * rectangle.. */
            if (use_blacklist) {
                for (auto &br : *blacklist) {
                    if (current_rectangle.subset(br)) {
                        debugl(1, "---------- DEPTH: %d -- rectangle: [%20.13E:%20.13E]x[%20.13E:%23.10E] --- dalpha: %20.13E, dbeta: %20.13E -- AR: %20.13E -- BLACKLISTED => setting to irrelevant.\n",
                            depth, alpha0, alpha1, beta0, beta1, (alpha1 - alpha0), (beta1 - beta0),
                            std::max( (alpha1 - alpha0) / (beta1 - beta0), (beta1 - beta0) / (alpha1 - alpha0) ) );

                        /* set blacklist flag, break loop */
                        rec_blacklisted = true;
                        break;
                    }
                }

                /* if rectangle has been blacklisted, break the while loop. subdivide == false holds here => interval is
                 * immediately discarded */
                if (rec_blacklisted) {
                    break;
                }
            }

            /* if neither alpha nor beta have converged, compute the aspect ratio of the rectangle and freeze if
             * necessary. notice that intervals considered here have never converged yet, i.e.  width, height or both
             * are still larger than tol <=> !alpha_converged or !beta_converged */
            if (!alpha_converged && !beta_converged) {
                rec_width   = (alpha1 - alpha0);
                rec_height  = (beta1  - beta0);
                rec_ar      = std::max( rec_width / rec_height, rec_height / rec_width);

                /* if rectangle ar is too large, freeze larger dimension variable */
                if (rec_ar > rec_armax) {
                    debugl(2, "AR: %f: maximum aspect %f ratio exceeded.... ", rec_ar, rec_armax);
                    /* freeze alpha / x */
                    if (rec_width > rec_height) {
                        debugl(2, "freezing beta.\n\n");
                        alpha_frozen    = false;
                        beta_frozen     = true;
                    }
                    /* freeze beta / y */
                    else {
                        debugl(2, "freezing alpha.\n\n");
                        alpha_frozen    = true;
                        beta_frozen     = false;
                    }
                }
                /* rectangle ar not too large, nothing frozen */
                else {
                    alpha_frozen    = false;
                    beta_frozen     = false;
                }

            }
            /* otherwise, either width or height (xor) has converged and that axis must be frozen. */
            else {
                alpha_frozen    = alpha_converged;
                beta_frozen     = beta_converged;
            }

            /* compute new interval */
            BiLinClip_getNewRectangle(
                   *p, *q,
                    m, n,
                    BiLinClip_A00,
                    BiLinClip_A01,
                    BiLinClip_A10,
                    BiLinClip_L00,
                    BiLinClip_L01,
                    BiLinClip_L10,
                    alpha0, alpha1, alpha_frozen,
                    beta0, beta1, beta_frozen,
                    rec_relevant,
                    new_alpha0, new_alpha1,
                    new_beta0, new_beta1,
                    eps, linsolve_eps
                );

            /* if rectangle is relevant, consider it further */
            if (rec_relevant) {
                /* if neither alpha nor beta are frozen */
                if (!alpha_frozen && !beta_frozen) {
                    /* compute old and new interval sizes */
                    dalpha      = (alpha1 - alpha0);
                    dbeta       = (beta1 - beta0);

                    new_dalpha  = (new_alpha1 - new_alpha0);
                    new_dbeta   = (new_beta1 - new_beta0);

                    debugl(1, "interval relevant: old interval sizes: (%20.13E)x(%20.13E), new interval sizes: (%20.13E)x(%20.13E)\n", dalpha, dbeta, new_dalpha, new_dbeta);

                    /* check if alpha or beta intervals have converged and freeze if so */
                    if (new_dalpha < tol) {
                        alpha_converged = true;
                        alpha_frozen    = true;
                    }
                    if (new_dbeta < tol) {
                        beta_converged  = true;
                        beta_frozen     = true;
                    }

                    /* check for convergence */
                    if (alpha_converged && beta_converged) {
                        debugl(1, "found root rectangle in depth %2d: [%20.13E, %20.13E]x[%20.13E, %20.13E]. p(midpoint) = %20.13E, q(midpoint) = %20.13E\n",
                                depth + 1, new_alpha0, new_alpha1, new_beta0, new_beta1,
                                pinput.eval( (new_alpha0 + new_alpha1) / 2.0, (new_beta0 + new_beta1) / 2.0),
                                qinput.eval( (new_alpha0 + new_alpha1) / 2.0, (new_beta0 + new_beta1) / 2.0)
                            );
                        roots.push_back( RealRectangle<R>(new_alpha0, new_alpha1, new_beta0, new_beta1) );
                        /*
                        printf("DEPTH: %d --ROOTZ: \\node[minimum width=0.1cm, circle, inner sep=0pt, draw=black, fill=orange!100] at (axis cs: %20.13E, %20.13E) {};\n",
                                depth + 1, new_alpha0, new_beta0);
                        */
                        break;
                    }
                    /* no convergence. restrict polynomials to new rectangle and check shrinkage */
                    else {
                        /* split with deCasteljau using the original polynomials pinput and qinput,
                         * so as not to accumulate round-off errors.. */
                        pinput.clipToInterval(new_alpha0, new_alpha1, new_beta0, new_beta1, p);
                        qinput.clipToInterval(new_alpha0, new_alpha1, new_beta0, new_beta1, q);

                        /* commit new rectangle boundaries */
                        alpha0  = new_alpha0;
                        alpha1  = new_alpha1;
                        beta0   = new_beta0;
                        beta1   = new_beta1;

                        /* check shrink factor. not via diameter, but via reduction of area. for a
                         * square, halfing the diameter amounts to quartering the area, so, as a
                         * heuristic setting, take this condition */
                        debugl(1, "area shrink factor: %f\n", (new_dalpha / dalpha) * (new_dbeta / dbeta));
                        if ( (new_dalpha / dalpha) * (new_dbeta / dbeta) > 0.9) {
                            debugl(1, "rectangle area hasn't shrunk by a factor of at least four... subdividing.\n");
                            subdivide = true;
                            break;
                        }
                        else {
                            debugl(1, "rectangle area has shrunk by a factor of at least four... keep shrinking..\n");
                            depth++;
                        }
                    }
                }
                /* if alpha is frozen but not beta */
                else if (alpha_frozen && !beta_frozen) {
                    /* we got a new beta interval. alpha is frozen and CAN have already converged
                     * (but not necessarily, might have been frozen due to AR). */
                    dbeta       = (beta1 - beta0);
                    new_dbeta   = (new_beta1 - new_beta0);

                    /* if beta has converged, it MUST be frozen. check if alpha interval has
                     * converged as well. if so, we got a root. if not, alpha must be unfrozen */
                    if (new_dbeta < tol) {
                        beta_converged  = true;
                        beta_frozen     = true;

                        /* both converged => got a root */
                        if (alpha_converged) {
                            debugl(1, "found root rectangle in depth %2d: [%20.13E, %20.13E]x[%20.13E, %20.13E]. p(midpoint) = %20.13E, q(midpoint) = %20.13E\n",
                                    depth + 1, alpha0, alpha1, new_beta0, new_beta1,
                                    pinput.eval( (alpha0 + alpha1) / 2.0, (new_beta0 + new_beta1) / 2.0),
                                    qinput.eval( (alpha0 + alpha1) / 2.0, (new_beta0 + new_beta1) / 2.0)
                                );
                            roots.push_back( RealRectangle<R>(alpha0, alpha1, new_beta0, new_beta1) );
                            /*
                            printf("DEPTH: %d --ROOTZ: \\node[minimum width=0.1cm, circle, inner sep=0pt, draw=black, fill=orange!100] at (axis cs: %20.13E, %20.13E) {};\n",
                                    depth + 1, alpha0, beta0);
                            */
                            break;
                        }
                        /* only beta interval converged => unfreeze alpha and commit new beta
                         * interval values, including the restriction of p and q to the new, final
                         * values for beta */
                        else {
                            debugl(1, "beta converged, alpha not converged but frozen => unfreezing alpha.\n");
                            alpha_frozen    = false;

                            pinput.clipToInterval(alpha0, alpha1, new_beta0, new_beta1, p);
                            qinput.clipToInterval(alpha0, alpha1, new_beta0, new_beta1, q);

                            beta0           = new_beta0;
                            beta1           = new_beta1;
                            depth++;
                        }
                    }
                    /* no convergence yet, beta is not frozen, alpha remains frozen. */
                    else {
                        /* restrict p and q to new beta interval with deCasteljau. again, use the input
                         * polynomials pinput and qinput to prevent accumulation of rounding errors. */
                        pinput.clipToInterval(alpha0, alpha1, new_beta0, new_beta1, p);
                        qinput.clipToInterval(alpha0, alpha1, new_beta0, new_beta1, q);

                        /* commit new beta interval boundaries */
                        beta0   = new_beta0;
                        beta1   = new_beta1;

                        /* check shrink factor. if beta interval has shrunk by a factor at least
                         * two, keep going, otherwise subdivide */
                        if ( new_dbeta / dbeta > 0.5) {
                            debugl(1, "beta shrink factor: %f. alpha frozen and beta interval hasn't shrunk by a factor of at least two... subdividing.\n", new_dbeta / dbeta);
                            subdivide = true;
                            break;
                        }
                        else {
                            debugl(1, "beta shrink factor: %f. alpha frozen and beta interval has shrunk by a factor at least two... keep shrinking.\n", new_dbeta / dbeta);
                            depth++;
                        }
                    }
                }
                /* if beta is frozen but not alpha */
                else if (!alpha_frozen && beta_frozen) {
                    /* we got a new alpha interval. beta is frozen and CAN have already converged
                     * (but not necessarily, might be frozen due to AR). */
                    dalpha      = (alpha1 - alpha0);
                    new_dalpha  = (new_alpha1 - new_alpha0);

                    /* if alpha has converged, it MUST be frozen. check if beta interval has
                     * converged as well. if so, we got a root. if not, beta must be unfrozen */
                    if (new_dalpha < tol) {
                        debugl(1, "alpha converged, new_dalpha = %20.13E\n", new_dalpha);
                        alpha_converged  = true;
                        alpha_frozen     = true;

                        /* both converged => got a root */
                        if (beta_converged) {
                            debugl(1, "found root rectangle in depth %2d: [%20.13E, %20.13E]x[%20.13E, %20.13E]. p(midpoint) = %20.13E, q(midpoint) = %20.13E\n",
                                    depth + 1, alpha0, alpha1, new_beta0, new_beta1,
                                    p->eval( (alpha0 + alpha1) / 2.0, (new_beta0 + new_beta1) / 2.0),
                                    q->eval( (alpha0 + alpha1) / 2.0, (new_beta0 + new_beta1) / 2.0)
                                );
                            roots.push_back( RealRectangle<R>(alpha0, alpha1, new_beta0, new_beta1) );
                            /*
                            printf("DEPTH: %d --ROOTZ: \\node[minimum width=0.1cm, circle, inner sep=0pt, draw=black, fill=orange!100] at (axis cs: %20.13E, %20.13E) {};\n",
                                    depth + 1, alpha0, beta0);
                                    */
                            break;
                        }
                        /* only alpha interval converged => unfreeze beta and commit new alpha
                         * interval values, including the restriction of p and q to the new, final
                         * values for alpha */
                        else {
                            debugl(1, "alpha converged, yet beta not converged but frozen => unfreezing beta.\n");

                            beta_frozen     = false;

                            pinput.clipToInterval(new_alpha0, new_alpha1, beta0, beta1, p);
                            qinput.clipToInterval(new_alpha0, new_alpha1, beta0, beta1, q);

                            alpha0          = new_alpha0;
                            alpha1          = new_alpha1;
                            depth++;
                        }
                    }
                    /* no convergence yet, alpha is not frozen, beta remains frozen. */
                    else {
                        /* restrict p and q to new alpha interval with deCasteljau. again, use the input
                         * polynomials pinput and qinput to prevent accumulation of rounding errors. */
                        pinput.clipToInterval(new_alpha0, new_alpha1, beta0, beta1, p);
                        qinput.clipToInterval(new_alpha0, new_alpha1, beta0, beta1, q);

                        /* commit new alpha interval boundaries */
                        alpha0   = new_alpha0;
                        alpha1   = new_alpha1;

                        /* check shrink factor. if alpha interval has shrunk by a factor at least
                         * two, keep going, otherwise subdivide */
                        if ( new_dalpha / dalpha > 0.5) {
                            debugl(1, "alpha shrink factor: %f. beta frozen and alpha interval hasn't shrunk by a factor of at least two... subdividing.\n", new_dalpha / dalpha);
                            subdivide = true;
                            break;
                        }
                        else {
                            debugl(1, "alpha shrink factor: %f. beta frozen and alpha interval has shrunk by a factor at least two... keep shrinking.\n", new_dalpha / dalpha);
                            depth++;
                        }
                    }
                }
                /* this should never happen, since only one dimension is frozen when aspect ratio gets
                 * too large and both are frozen in case of convergence (whence we should never reach
                 * this part) */
                else {
                    throw("BiLinClip_roots(): alpha and beta frozen in subdivision step. Now that must not happen..");
                }
            }
            /* current rectangle irrelevant, break the loop. here, subdivide is still false, hence polynomials
             * p and q are deleted and this branch of the search tree is left */
            else {
                debugl(1, "current rectangle [%f, %f]x[%f, %f] in depth %d irrelevant. moving along..\n", alpha0, alpha1, beta0, beta1, depth);
                break;
            }
        }

        /* while loop has been broken. if subdivide == true, the interval is relevant, still not
         * converged in alpha or beta (or both) and has seized shrinking exponentially in area..
         * subdivide rectangle with respect to unfrozen axes. 
         * if subdivide == false, it is irrelevant or has been blacklisted */
        if (subdivide) {
            /* if neither alpha nor beta are frozen, split p and q into four polynomials representing
             * the four subrectangles arising from the split at the middle (0.5, 0.5) with respect
             * to the unit square [0,1]^2, to which p and q are always affine linearly transformed by
             * deCasteljau. */
            if (!alpha_frozen && !beta_frozen) {
                //printf("neither alpha nor beta frozen: subdividing rectangle into four subrectangles.. \n");

                /* alloc new polys */
                pleft_down  = new BiBernsteinPolynomial<R, R>();
                pright_down = new BiBernsteinPolynomial<R, R>();
                pright_up   = new BiBernsteinPolynomial<R, R>();
                pleft_up    = new BiBernsteinPolynomial<R, R>();

                qleft_down  = new BiBernsteinPolynomial<R, R>();
                qright_down = new BiBernsteinPolynomial<R, R>();
                qright_up   = new BiBernsteinPolynomial<R, R>();
                qleft_up    = new BiBernsteinPolynomial<R, R>();

                /* check if midpoint (0.5, 0.5) is common root of p and q. if so..  it seems very unwise to cut into
                 * small slices of width / height 0.25*tol like in the univariate case, since those strips will have
                 * zero width / height long before the height / width is small. it seems better to keep the rectangles
                 * as mush squares as possible. maybe just push the root onto the root vector and leave it at that.  in
                 * the worst case, it will be found multiple times and we'll have to check that in the end, or just
                 * evaluate the same root multiple times outside this routine, which is very cheap compared to the
                 * routine itself.. */
                R pval = p->eval(0.5, 0.5);
                R qval = q->eval(0.5, 0.5);
                debugl(1, "subdividing: p(midpoint) = %f, q(midpoint) = %f\n", pval, qval);
                if ( std::abs(pval) < eps && std::abs(qval) < eps) {
                    debugl(0, "found root at midpoint of rectangle [%f, %f]x[%f, %f] to be subdividied.\n", alpha0, alpha1, beta0, beta1);
                }

                /* split with deCasteljau at (0.5, 0.5) */
                p->split_xy(0.5, 0.5, pleft_down, pright_down, pright_up, pleft_up);
                q->split_xy(0.5, 0.5, qleft_down, qright_down, qright_up, qleft_up);

                /* push new rectangles and polys onto stack */
                alpha_middle    = (alpha0 + alpha1) / 2.0;
                beta_middle     = (beta0  + beta1)  / 2.0;

                S.push_back( BiLinClip_Tuple<R>(pleft_down,     qleft_down,     alpha0,         alpha_middle,   alpha_converged,  beta0,          beta_middle,    beta_converged, depth + 1) );
                S.push_back( BiLinClip_Tuple<R>(pright_down,    qright_down,    alpha_middle,   alpha1,         alpha_converged,  beta0,          beta_middle,    beta_converged, depth + 1) );
                S.push_back( BiLinClip_Tuple<R>(pright_up,      qright_up,      alpha_middle,   alpha1,         alpha_converged,  beta_middle,    beta1,          beta_converged, depth + 1) );
                S.push_back( BiLinClip_Tuple<R>(pleft_up,       qleft_up,       alpha0,         alpha_middle,   alpha_converged,  beta_middle,    beta1,          beta_converged, depth + 1) );
            }
            /* if alpha is frozen but not beta */
            else if (alpha_frozen && !beta_frozen) {
                /* alloc new polys */
                pup     = new BiBernsteinPolynomial<R, R>();
                pdown   = new BiBernsteinPolynomial<R, R>();
                qup     = new BiBernsteinPolynomial<R, R>();
                qdown   = new BiBernsteinPolynomial<R, R>();

                /* split p and q with respect to beta at 0.5 */
                p->split_y(0.5, pdown, pup);
                q->split_y(0.5, qdown, qup);

                /* push new rectangles and polys onto stack */
                beta_middle = (beta0 + beta1) / 2.0;

                S.push_back( BiLinClip_Tuple<R>(pdown,    qdown,  alpha0,     alpha1,   alpha_converged,  beta0,          beta_middle,    beta_converged, depth + 1) );
                S.push_back( BiLinClip_Tuple<R>(pup,      qup,    alpha0,     alpha1,   alpha_converged,  beta_middle,    beta1,          beta_converged, depth + 1) );
            }
            /* if beta is frozen but not alpha */
            else if (!alpha_frozen && beta_frozen) {
                /* alloc new polys */
                pleft   = new BiBernsteinPolynomial<R, R>();
                pright  = new BiBernsteinPolynomial<R, R>();
                qleft   = new BiBernsteinPolynomial<R, R>();
                qright  = new BiBernsteinPolynomial<R, R>();

                /* split p and q with respect to alpha at 0.5 */
                p->split_x(0.5, pleft, pright);
                q->split_x(0.5, qleft, qright);

                /* push new rectangles and polys onto stack */
                alpha_middle = (alpha0 + alpha1) / 2.0;

                S.push_back( BiLinClip_Tuple<R>(pleft,    qleft,     alpha0,         alpha_middle,   alpha_converged,  beta0,    beta1,  beta_converged, depth + 1) );
                S.push_back( BiLinClip_Tuple<R>(pright,   qright,    alpha_middle,   alpha1,         alpha_converged,  beta0,    beta1,  beta_converged, depth + 1) );
            }
            /* this should never happen, since only one dimension is frozen when aspect ratio gets
             * too large and both are frozen in case of convergence (whence we should never reach
             * this part) */
            else {
                throw("BiLinClip_roots(): alpha and beta frozen in subdivision step. Now that must not happen..");
            }
        }

        /* delete polys, this branch of the search tree is left, either because a root has been found or
         * because the rectangle has been discarded as irrelevant or blacklisted. */
        delete p;
        delete q;
    }

    setDebugComponent(DBG_GLOBAL);
}
