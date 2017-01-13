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

#ifndef AUX_H
#define AUX_H

#ifndef TIME
#define TIME

#include <sys/time.h>
#include <time.h>

#define TIMER_REGISTERS 64

#endif

#include "common.hh"
#include "Vec2.hh"
#include "Vec3.hh"
#include "StaticVector.hh"
#include "StaticMatrix.hh"

#ifdef WITH_BOOST
	#include <boost/math/special_functions/binomial.hpp>
#endif

/* forward declaration of BoundingBox to break cyclical dependency. */
template <typename R> class BoundingBox;

namespace Aux {

    namespace Timing {
        void	tick(int i);
        double	tack(int i);
        double  doubletime();
    }

    namespace Numbers {
        /* FIXME: convert to templates, use std::{sin,cos,sqrt, ..} template specialization wrappers for arithmetic */
        double  frand(double min, double max);
        double  fmin3(double a, double b, double c);
        double  fmax3(double a, double b, double c);
        int     sign(int64_t d);
        double  deg2rad(double angle);
        double  rad2deg(double angle);

        /* template function inf, specializations for float, long and long double are in aux.cc */
        template<typename R>
        R
        inf();

        /* -------------- */

        /* numerically save init of binomial coefficients */
        template <typename R, uint32_t n>
        R
        bicof(uint32_t k)
        {
#ifdef WITH_BOOST
            return boost::math::binomial_coefficient<R>(n, k);
#else
            static StaticVector<n+1, R> bicof;
            static bool bicof_recompute = true;

            /* indicate necessity to (re)compute more binomial coefficients by toggling flag if required n is
             * larger than maximum precomputed one */
            if (bicof_recompute)
            {
                debugl(2, "initBinomialCoefficients(): ... \n");

                /* fill in binomial coefficients */
                bicof[0] = R(1);
                for (k = 1; k < n; ++k)
                    bicof[k] = floor(0.5 + exp(lgamma(n+1) - lgamma(k+1) - lgamma(n-k+1)));

                bicof[n] = (R)1;

                debugl(2, "done.\n");

                bicof_recompute = false;
            }
            return k < n+1 ? bicof[k] : (R)0;
#endif
        }

    }

    namespace Numerics {
        /* Thomson variant: extremely efficient, linear time, but of at least questionable stability. maybe better cholseky? */
        template <typename R>
        void
        solveSplineCoefficientSystem(
            std::vector<R> const   &t,
            std::vector<R> const   &x,
            R const                &ds0,
            R const                &dsn,
            std::vector<R>          &result_vector)
        {
            /*
            static bool narf_init = false;
            if (!narf_init) {
                NA_init();
                narf_init = true;
            }
            */

            int i, j, n = t.size() - 1; 

            /*
            printf("\n\n ------------------------------ \n\n solveSplineCoefficientSystem(): n = %d, ds0 = %f, dsn = %f\n", n, ds0, dsn);
            printf("t = [ ");
            for (i = 0; i <= n; i++) {
                printf("%+5.4e ", t[i]);
            }
            printf("]\n");
            printf("x = [ ");
            for (i = 0; i <= n; i++) {
                printf("%+5.4e ", x[i]);
            }
            printf("]\n");
            */

            /* parameter (knot) differences */
            std::vector<R>     h(n + 1);

            h[0] = 0.0;
            for (i = 0; i <= n - 1; i++) {
                h[i + 1] = t[i+1] - t[i];
            }

            /* moments of the spline, i.e. second derivatives. solution vector of the constructed system */
            std::vector<R>     M(n + 1, 0.0); 

            /* tridiagonal matrix encoding: every row is encoded as three values */
            std::vector<std::vector<R> >   A(n+1);
            for (i = 0; i < n+1; i++) {
                A[i].resize(3, 0.0);
            }

            /* inhomogeneity vector */
            std::vector<R> d(n+1);

            /* let A_real denote the real (n+1)x(n+1) matrix, not the tridiagonal matrix encoding.
             * first row of the matrix, A[0][0] is actually "A_real[0][-1]" and not part of the matrix, because
             * first row has only two entries. upper left element is 2.0, A_real[0][1] = lambda_0  = 1.0 for
             * hermite boundary conditions */
            A[0][0] = 0.0;
            A[0][1] = 2.0;
            A[0][2] = 1.0;

            /* value for d[0] is specially derived for hermite boundary conditions */
            d[0]    = 6.0 / h[1] * ( (x[1] - x[0]) / h[1] - ds0);
            
            /* last row. A_real[n][n-1] = mu_n = 1.0 for hermite boundary conditions. A_real[n][n] = 2.0.
             * A[n][2] = A_real[n][n+1] and is not part of the matrix */
            A[n][0] = 1.0;
            A[n][1] = 2.0;
            A[n][2] = 0.0;

            /* special value for hermite boundary conditions */
            d[n]    = 6.0 / h[n] * (dsn - (x[n] - x[n-1]) / h[n]);

            for (i = 1; i < n; i++) {
                /* mu_i, 2.0, lambda_i */
                A[i][0] = h[i] / (h[i] + h[i+1]);
                A[i][1] = 2.0;
                A[i][2] = h[i+1] / (h[i] + h[i+1]);

                /* inhomogeneity vector entry i */
                d[i]    = 6.0 / (h[i] + h[i+1]) * ( (x[i+1] - x[i]) / h[i+1] - (x[i] - x[i-1]) / h[i]);
            }

            /* ---------- debug print ----------- */
            /*
            printf("n= %d\n", n);
            NA_matrix A_real = NA_minit(n+1, n+1);
            A_real[0][0]    = A[0][1];
            A_real[0][1]    = A[0][2];
            for (i = 1; i < n; i++) {
                A_real[i][i-1] = A[i][0];
                A_real[i][i] = A[i][1];
                A_real[i][i + 1] = A[i][2];
            }
            A_real[n][n-1]  = A[n][0];
            A_real[n][n]    = A[n][1];

            NA_mprint(A_real, n + 1, n+1);
            printf("d = [ ");
            for (i = 0; i < n+1; i++) {
                printf("+%5.4e ", d[i]);
            }
            printf("]\n");
            */
            /* ---------- debug print ----------- */

            /* eliminate subdiagonal entries in tridiagonal system.
             * add a multiple of row i to row (i+1) to eliminate element A_real[i+1][i] */
            R lambda;
            for (i = 0; i < n; i++) {
                lambda = -A[i+1][0] / A[i][1]; 
                A[i+1][0] += lambda * A[i][1];
                A[i+1][1] += lambda * A[i][2];
                d[i+1]    += lambda * d[i];
            }

            /* ---------- debug print ----------- */
            /*
            A_real[0][0]    = A[0][1];
            A_real[0][1]    = A[0][2];
            for (i = 1; i < n; i++) {
                A_real[i][i-1] = A[i][0];
                A_real[i][i] = A[i][1];
                A_real[i][i + 1] = A[i][2];
            }
            A_real[n][n-1]  = A[n][0];
            A_real[n][n]    = A[n][1];

            NA_mprint(A_real, n + 1, n+1);
            */
            /* ---------- debug print ----------- */

            /* solve using backsubstitution */
            M[n]    = d[n] / A[n][1];
            for (i = n-1; i >= 0; i--) {
                M[i] = (d[i] - A[i][2] * M[i+1]) / A[i][1];
            } 

            /*
            printf("M = [ ");
            for (i = 0; i < n+1; i++) {
                printf("%+5.4e ", M[i]);
            }
            printf("]\n");
            */

            /* quantities used in the moment representation of the spline segment functions */
            std::vector<R> a(n), alpha(n), beta(n), gamma(n), delta(n);
            for (j = 0; j < n; j++) {
                //b[j]        = x[j] - ( M[j] * h[j+1] * h[j+1]) / 6.0;
                a[j]        = (x[j+1] - x[j]) / h[j+1] - h[j+1]/6.0*(M[j+1] - M[j]);
                alpha[j]    = x[j];
                gamma[j]    = M[j] / 2.0;
                beta[j]     = ( -M[j] * h[j+1])/2.0 + a[j];
                delta[j]    = (M[j+1] - M[j]) / (6.0 * h[j+1]);
            }

            /* the j-th spline is now given as:
             * delta[j]*(t - t[j])^3 + gamma[j]*(t - t[j])^2 + beta[j]*(t - t[j]) + alpha[j]
             * to get the monomial coefficients, this has to be expanded.. */
             
            /* output coefficients */
            result_vector.resize(4*n);

            R c_0, c_1, c_2, c_3;
            debugTabInc();
            for (j = 0; j < n; j++) {
                /* recover coefficients from solution */
                //c_0 = alpha[j] - beta[j]*t[j] + gamma[j]*t[j]*t[j] - delta[j]*t[j]*t[j]*t[j];
                c_0 = alpha[j] + t[j]*(- beta[j] + t[j]*(gamma[j] - delta[j]*t[j]) );
                //c_1 = beta[j] - 2.0*gamma[j]*t[j] + 3.0*delta[j]*t[j]*t[j];
                c_1 = beta[j] + t[j]*(- 2.0*gamma[j] + t[j]*(+ 3.0*delta[j]) );
                c_2 = gamma[j] + -3.0*delta[j]*t[j];
                c_3 = delta[j];

                debugl(2, "segment function = %d, c_0: %f, c_1: %f, c_2: %f, c_3: %f\n", j, c_0, c_1, c_2, c_3);
                debugl(2, "a[%d] = %f, alpha[%d] = %f, beta[%d] = %f, gamma[%d] = %f, delta[%d] = %f\n", j, a[j], j, alpha[j], j, beta[j], j, gamma[j], j, delta[j]);
                debugl(2, "t[%d] = %f, x[%d] = %f\n", j, t[j], j, x[j]);

                result_vector[4*j    ] = c_3;
                result_vector[4*j + 1] = c_2;
                result_vector[4*j + 2] = c_1;
                result_vector[4*j + 3] = c_0;
            }
            debugTabDec();
        }
    }

    namespace VecMat {
        template <typename R>
        R
        dist2(
            Vec3<R> const &u,
            Vec3<R> const &v)
        {
            return ((u-v).len2());
        }

        template<typename R>
        R
        dist2squared(
            Vec3<R> const &u,
            Vec3<R> const &v)
        {
            return ((u-v).len2squared());
        }

        template <typename R>
        Vec3<R>
        xbase()
        {
            return Vec3<R>((R)1.0, (R)0.0, (R)0.0);
        }

        template <typename R>
        Vec3<R>
        ybase()
        {
            return Vec3<R>((R)0.0, (R)1.0, (R)0.0);
        }

        template <typename R>
        Vec3<R>
        zbase()
        {
            return Vec3<R>((R)0.0, (R)0.0, (R)1.0);
        }

        template <typename R>
        Vec3<R>
        nullvec()
        {
            return Vec3<R>((R)0.0, (R)0.0, (R)0.0);
        }

        template <typename R>
        inline void
        minVec3(Vec3<R>& out, const Vec3<R>& a, const Vec3<R>& b)
        {
            out[0] = std::min(a[0], b[0]);
            out[1] = std::min(a[1], b[1]);
            out[2] = std::min(a[2], b[2]);
        }

        template <typename R>
        inline void
        maxVec3(Vec3<R>& out, const Vec3<R>& a, const Vec3<R>& b)
        {
            out[0] = std::max(a[0], b[0]);
            out[1] = std::max(a[1], b[1]);
            out[2] = std::max(a[2], b[2]);
        }

        template <typename R>
        inline Vec3<R>
        onesVec3()
        {
            return ( Vec3<R>(1.0, 1.0, 1.0) );
        }

        template <typename R>
        inline Vec3<R>
        randUnitVec3()
        {
            Vec3<R> r;

            r[0] = Aux::Numbers::frand(-1.1, 1.1);
            r[1] = Aux::Numbers::frand(-1.1, 1.1);
            r[2] = Aux::Numbers::frand(-1.1, 1.1);
            r.normalize();

            return r;
        }

        template <typename R>
        inline Vec3<R>
        fabsVec3(Vec3<R> const &x)
        {
            Vec3<R> r;

            r[0] = std::abs(x[0]);
            r[1] = std::abs(x[1]);
            r[2] = std::abs(x[2]);

            return r;
        }

        template <typename R>
        inline bool
        isfiniteVec3(Vec3<R> const &x)
        {
            return (std::isfinite(x[0]) &&
                    std::isfinite(x[1]) &&
                    std::isfinite(x[2])
                );
        }

        inline Vec2
        randUnitVec2()
        {
            Vec2    r;
            r[0] = Aux::Numbers::frand(-1.1, 1.1);
            r[1] = Aux::Numbers::frand(-1.1, 1.1);
            r.normalize();

            return r;

        }


        template <uint32_t N, typename T>
        inline StaticVector<N, T>
        kronecker_static_vec(uint32_t i)
        {
            StaticVector<N, T> r(0);
            r[i] = 1;
            return r;
        }

        template <uint32_t M, uint32_t N, typename T>
        inline StaticMatrix<M, N, T>
        kronecker_static_mat(uint32_t i, uint32_t j)
        {
            StaticMatrix<M, N, T> mat;
            mat.fill(0);
            mat(i,j) = 1;
            return mat;
        }

        template <typename R>
        void
        completeOrthonormalBase(
            Vec3<R> &x,
            Vec3<R> &y,
            Vec3<R> &z)
        {
            Vec3<R> nullvec = Aux::VecMat::nullvec<R>();

            /* if x is given */
            if (x != nullvec) {
                if (y == nullvec && z == nullvec) {
                    /* we need to construct y and z from x */
                    y.set(x[1], -x[0], 0.0);
                    y.normalize();
                    z = x.cross(y);
                }
                else if (y == nullvec) {
                    /* z == 0 -> we need to construct y from x and z */
                    y = z.cross(x);
                }
                else {
                    /* y == 0 -> we need to construct z from x and y */
                    z = x.cross(y);
                }
            }
            /* and so on .. */
            else if (y != nullvec) {
                if (x == nullvec && z == nullvec) {
                    /* we got y, construct x and z */
                    x.set(y[1], -y[0], 0.0);
                    x.normalize();
                    z = x.cross(y);
                }
                else if (x == nullvec) {
                    /* we got y and z, need x */
                    x = y.cross(z);
                }
                else {
                    /* we got y and x, need z */ 
                    z = x.cross(y);
                }
            }
            else if (z != nullvec) {
                if (y == nullvec && y == nullvec) {
                    x.set(z[1], -z[0], 0.0);
                    x.normalize();
                    y = z.cross(x);
                }
                else if (y == nullvec) {
                    y = z.cross(x);
                }
                else {
                    x = y.cross(z);
                }
            }
        }

        template <typename R>
        R
        det3x3(R A[3][3])
        {
            return (
                    A[0][0]*A[1][1]*A[2][2] +
                    A[0][1]*A[1][2]*A[2][0] +
                    A[0][2]*A[1][0]*A[2][1] -
                    A[0][2]*A[1][1]*A[2][0] -
                    A[0][1]*A[1][0]*A[2][2] -
                    A[0][0]*A[1][2]*A[2][1]
                   );
        }

        template <typename R>
        void
        print3x3(R A[3][3])
        {
            printf("(-------------------)\n");
            printf("( %9.4f %9.4f %9.4f )\n", A[0][0], A[0][1], A[0][2]);
            printf("( %9.4f %9.4f %9.4f )\n", A[1][0], A[1][1], A[1][2]);
            printf("( %9.4f %9.4f %9.4f )\n", A[2][0], A[2][1], A[2][2]);
            printf("(-------------------)\n");
        }
    }

    namespace Allocation {
        template<typename T>
        T **
        alloc2dArray(uint32_t n0, uint32_t n1)
        {
            uint32_t i0;

            T **A = new T * [n0];

            for (i0 = 0; i0 < n0; i0++) {
                A[i0] = new T [n1];
            }

            return A;
        }

        template<typename T>
        void
        free2dArray(T **A, uint32_t n0, uint32_t n1)
        {
            uint32_t i0;

            for (i0 = 0; i0 < n0; i0++) {
                delete[] A[i0];
            }

            delete[] A;
        }

        template<typename T>
        T ***
        alloc3dArray(uint32_t n0, uint32_t n1, uint32_t n2)
        {
            uint32_t i0, i1;

            T ***A = new T ** [n0];

            for (i0 = 0; i0 < n0; i0++) {
                A[i0] = new T* [n1];
                for (i1 = 0; i1 < n1; i1++) {
                    A[i0][i1] = new T [n2];
                }
            }

            return A;
        }

        template<typename T>
        void
        free3dArray(T ***A, uint32_t n0, uint32_t n1, uint32_t n2)
        {
            uint32_t i0, i1;

            for (i0 = 0; i0 < n0; i0++) {
                for (i1 = 0; i1 < n1; i1++) {
                    delete[] A[i0][i1];
                }

                delete[] A[i0];
            }

            delete [] A;
        }


        template<typename T>
        T ****
        alloc4dArray(uint32_t n0, uint32_t n1, uint32_t n2, uint32_t n3)
        {
            uint32_t i0, i1, i2;

            T ****A = new T *** [n0];

            for (i0 = 0; i0 < n0; i0++) {
                A[i0] = new T** [n1];
                
                for (i1 = 0; i1 < n1; i1++) {
                    A[i0][i1] = new T* [n2];

                    for (i2 = 0; i2 < n2; i2++) {

                        A[i0][i1][i2] = new T [n3];
                    }
                }
            }

            return A;
        }

        template<typename T>
        void
        free4dArray(T ****A, uint32_t n0, uint32_t n1, uint32_t n2, uint32_t n3)
        {
            uint32_t i0, i1, i2;

            for (i0 = 0; i0 < n0; i0++) {
                for (i1 = 0; i1 < n1; i1++) {
                    for (i2 = 0; i2 < n2; i2++) {
                        delete[] A[i0][i1][i2];
                    }

                    delete[] A[i0][i1];
                }

                delete[] A[i0];
            }

            delete[] A;
        }


        template<typename T>
        T ******
        alloc6dArray(uint32_t n0, uint32_t n1, uint32_t n2, uint32_t n3, uint32_t n4, uint32_t n5)
        {
            uint32_t i0, i1, i2, i3, i4;

            T ******A = new T ***** [n0];

            for (i0 = 0; i0 < n0; i0++) {
                A[i0] = new T **** [n1];

                for (i1 = 0; i1 < n1; i1++) {
                    A[i0][i1] = new T *** [n2];

                    for (i2 = 0; i2 < n2; i2++) {
                        A[i0][i1][i2] = new T ** [n3];

                        for (i3 = 0; i3 < n3; i3++) {
                            A[i0][i1][i2][i3] = new T * [n4];

                            for (i4 = 0; i4 < n4; i4++) {
                                A[i0][i1][i2][i3][i4] = new T [n5];
                            }
                        }
                    }
                }
            }

            return A;
        }

        template<typename T>
        void
        free6dArray(T ******A, uint32_t n0, uint32_t n1, uint32_t n2, uint32_t n3, uint32_t n4, uint32_t n5)
        {
            uint32_t i0, i1, i2, i3, i4;

            for (i0 = 0; i0 < n0; i0++) {
                for (i1 = 0; i1 < n1; i1++) {
                    for (i2 = 0; i2 < n2; i2++) {
                        for (i3 = 0; i3 < n3; i3++) {
                            for (i4 = 0; i4 < n4; i4++) {
                                delete[] A[i0][i1][i2][i3][i4];
                            }

                            delete[] A[i0][i1][i2][i3];
                        }

                        delete[] A[i0][i1][i2];
                    }

                    delete[] A[i0][i1];
                }

                delete[] A[i0];
            }

            delete[] A;
        }

    }

    namespace File {
        bool isEmpty(FILE *f);
    }

    namespace Geometry {
        namespace IntersectionTestResults {
            enum IntersectionTestResults_enum {
                A_SUBSET_B,
                B_SUBSET_A,
                DISJOINT,
                NOSUBSET_INTERSECT,
                INTERSECTION,
                POINT_IN,
                POINT_OUT,
                TEST_INCONCLUSIVE,
                EDGE_CASE,
                SUCCESS
            };
        }

        struct Tri2d {
            uint32_t    v0_id, v1_id, v2_id;

            Tri2d(uint32_t v0_id, uint32_t v1_id, uint32_t v2_id);
        };

        struct Vertex2d {
            uint32_t    id;
            Vec2        pos;

            Vertex2d() ;
            Vertex2d(uint32_t id, const Vec2& pos);
        };


        template <typename R>
        uint32_t
        simpleIntersect2AABB(
                Vec3<R> a_min,
                Vec3<R> a_max,
                Vec3<R> b_min,
                Vec3<R> b_max)
        {
            using namespace Aux::Geometry::IntersectionTestResults;

            if (
                a_min[0] > b_max[0] || a_max[0] < b_min[0] ||
                a_min[1] > b_max[1] || a_max[1] < b_min[1] ||
                a_min[2] > b_max[2] || a_max[2] < b_min[2] 
            )
            {
                return DISJOINT;
            }
            else {
                return INTERSECTION;
            }
        }


        template <typename R>
        uint32_t
        rayTriangle(
            const Vec3<R>&     p0,
            const Vec3<R>&     p1,
            const Vec3<R>&     v0,
            const Vec3<R>&     v1,
            const Vec3<R>&     v2,
            Vec3<R>    &x,
            R          &x_s,
            R          &x_t,
            R          &x_lambda,
            R           ieps = 1E-10)
        {
            using namespace Aux::Geometry::IntersectionTestResults;

            debugl(4, "Aux::Geometry::rayTriangle(): .. \n");
            debugTabInc();

            p0.print_debugl(3);
            p1.print_debugl(3);

            v0.print_debugl(3);
            v1.print_debugl(3);
            v2.print_debugl(3);

            /* calculate normal vector of the triangle plane */
            Vec3<R> u = v1 - v0;
            Vec3<R> v = v2 - v0;
            Vec3<R> p = p1 - p0;

            Vec3<R> n = u.cross(v);
            R denom = n * p;

            /* denom very small => infinite ray parallel to plane.. */
            if (std::abs(denom / n.len2() / p.len2()) < ieps) {
                debugl(1, "n.len2(): %5.4e. denom = %5.4e < ieps => segment parallel to plane.\n", n.len2(), denom);
                debugTabDec();
                return DISJOINT;
            }
            else {
                /* calculate parametric value lambda_plane of intersection point x: x = p_0 + lambda_plane(p_1 - p_0) */
                x_lambda = n * (v0 - p0) / denom;

                /* infinite ray intersects the plane definitely (not parallel if we reach this line), but the
                 * parametric value has to be in [0,1] for the segment (p0, p1) to intersect the plane */
                if (x_lambda < 0.0 || x_lambda > 1.0) {
                    debugl(2, "x_lambda = %10.5e not in [0,1] => finite segment does not intersect plane.\n", x_lambda);
                    debugTabDec();
                    return DISJOINT;
                }
                else {
                    debugl(2, "x_lambda = %10.5e in [0,1] => finite segment intersects plane.\n", x_lambda);
                    x = p0 + p*x_lambda;
                    Vec3<R> w = x - v0;

                    R uv = u*v;
                    R uu = u*u;
                    R vv = v*v;
                    R wu = w*u;
                    R wv = w*v;

                    denom = uv*uv - uu*vv;

                    x_s = (uv*wv - vv*wu) / denom;
                    x_t = (uv*wu - uu*wv) / denom;

                    /* barycentric parameters of triangle outside [0,1]^2 => point not in triangle */
                    if (x_s < 0.0 || x_t < 0.0 || (x_s + x_t > 1.0) ) {
                        debugl(4, "barycentric parametric coordinates x_s = %10.5e, x_t = %10.5e => point not on triangle\n", x_s, x_t);
                        debugl(4, "Aux::Geometry::rayTriangle(): DISJOINT.\n");
                        debugTabDec();
                        return DISJOINT;
                    }
                    /* got one */
                    else {
                        debugl(4, "rayTriangle(): hit. x_s: %5.4e, x_t: %5.4e, x_lambda: %5.4e\n", x_s, x_t, x_lambda);
                        debugl(4, "Aux::Geometry::rayTriangle(): INTERSECTION.\n");
                        debugTabDec();
                        return INTERSECTION;
                    }
                }
            }
        }

        template <typename R>
        bool
        raySphere(
            Vec3<R> const      &c,
            R const            &r,
            Vec3<R> const      &o,
            Vec3<R> const      &v,
            std::vector<R>     &lambda,
            R const            &eps = 1E-11)
        {
            Vec3<R> tmp;
            R       A, B, C, D;

            /* careful: '*' is overloaded and has additional scalar product meaning
             * here, '+' and '-' also have vector implementations */
            A       = v * v;
            tmp     = (o - c);
            B       = 2 * (tmp * v);
            C       = tmp * tmp - r * r;
            D       = B*B - 4*A*C;

            if (std::abs(D) < eps) {
                lambda.resize(1);
                lambda[0]  = -B / (2*A);
                return true;
            }
            else if (D < 0) {
                lambda.resize(0); 
                return false;
            }
            else {
                debug("two solutions. D = %f\n", D);
                lambda.resize(2);
                lambda[0]   = (-B - sqrt(D)) / (2*A);
                lambda[1]   = (-B + sqrt(D)) / (2*A);
                return true;
            }
        }

        template <typename R>
        void
        computeBaryCoordsOfProjectedPoint(
            const Vec3<R>&     p,
            const Vec3<R>&     v0,
            const Vec3<R>&     v1,
            const Vec3<R>&     v2,
            R          &s,
            R          &t)
        {
            debugl(4, "Aux::Geometry::computeBaryCoordsOfProjectedPoint().\n");
            debugTabInc();
          
            Vec3<R> u = v1 - v0;
            Vec3<R> v = v2 - v0;
            Vec3<R> w = p  - v0;
            Vec3<R> n = u.cross(v);

            R inv_denom, /* b0, */ b1, b2;

            inv_denom   = 1.0 / (n*n);
            b2          = ( (u.cross(w)) * n) * inv_denom;
            b1          = ( (w.cross(v)) * n) * inv_denom;
            /* b0          = 1.0 - b2 - b1; */

            debugl(5, "b0: %5.4f, b1: %5.4f, b2: %5.4f\n", 1.0 - b2 - b1, b1, b2);

            s = b1;
            t = b2;

            debugTabDec();
            debugl(4, "Aux::Geometry::computeBaryCoordsOfProjectedPoint(): done.\n");
        }


        uint32_t    lineSegmentLineSegment2d(
                        const Vec2&        p0,
                        const Vec2&        p1,
                        const Vec2&        q0,
                        const Vec2&        q1,
                        Vec2       &x,
                        double     &x_lambda,
                        double     &x_mu,
                        double      eps = 1E-15);

        uint32_t    rayLineSegment2d(
                        const Vec2&        p,
                        const Vec2&        dir,
                        const Vec2&        q0,
                        const Vec2&        q1,
                        Vec2       &x,
                        double     &x_lambda,
                        double     &x_mu,
                        double      eps = 1E-15);

        /* triangle - triangle intersection methods */
        template <typename R>
        bool
        triTriSharingEdge(
            const Vec3<R>& x,
            const Vec3<R>& u,
            const Vec3<R>& v,
            const Vec3<R>& y,
            R       eps = 1E-15)
        {
            using namespace Aux::Geometry::IntersectionTestResults;

            /* check if y lies in the plane of triangle (x, u, v) */
            Vec3<R> n, dn_uv, nz;
            R           dist_y;
            R           hp_x, hp_y;

            /* normal vector of triangle (x, u, v) */
            n       = (u - x).cross(v - x);
            n.normalize();

            /* normal vector pointing from u to v */
            dn_uv   = (v - u);
            dn_uv.normalize();

            /* normal vector on the plane defined by n and dn_uv => lies in the plane with normal vector
             * n and is orthogonal to dn_uv => points towards one half-plane */
            nz      = n.cross(dn_uv);

            dist_y  = (y - u) * n;

            if (std::abs(dist_y) < eps) {
                printf("y in plane of triangle (x, u, v)\n");
                /* the infinite line through u and v divides the plane pi with normal n in two half-planes.
                 * self-intersection occurs iff x and y lie in the same half-plane. */

                /* check dot product with nv: if either is numerically zero => point on edge uf => exception.
                 * otherwise, check signs */
                hp_x = (x - u) * nz;
                hp_y = (y - u) * nz;

                if ( std::abs(hp_x) < eps || std::abs(hp_y) < eps) {
                    printf("triTriSharingEdge(): x or y on ege (u, v): hp_x: %12.5e, hp_y: %12.5e\n", hp_x, hp_y);
                    throw("triTriSharingEdge(): x or y on ege (u, v)");
                }
                else {
                    /* check signs */
                    return ( (hp_x < 0.0 && hp_y < 0.0) || (hp_x > 0.0 && hp_y > 0.0) );
                }
            }
            else return false;
        }


        template <typename R>
        bool
        triTriSharingVertex(
            Vec3<R> a0,
            Vec3<R> a1,
            Vec3<R> s,
            Vec3<R> b0,
            Vec3<R> b1,
            R       eps = 1E-15)
        {
            using namespace Aux::Geometry::IntersectionTestResults;

            /* get normal of triangle (a0, a1, s) */
            Vec3<R> n;
            R       b0_dist, b1_dist;

            n = (a1 - a0).cross(s - a0);
            n.normalize();

            /* determine distances of b0 and b1 */
            b0_dist = (b0 - s) * n;
            b1_dist = (b1 - s) * n;
            //printf("\tb0_dist: %3.2e, b1_dist: %3.2e\n", b0_dist, b1_dist);

            /* four cases:
             *
             * 1. both on one side of the plane => no self-intersection possible
             * 2. one in the plane, one not => planar check / point containment
             * 3. both in the plane => entirely planar problem
             * 4. b0 on one side, b1 on the other => self-intersection iff 
             * (b0, b1) intersects (a0, a1, s) or (a0, a1) intersects (s, b0, b1) */
            if (std::abs(b0_dist) < eps || std::abs(b1_dist) < eps) {
                printf("triTriSharingEdge(): b0 and b1 both in plane pi of (a0, a1, b).\n");
                throw("triTriSharingEdge(): b0 and b1 both in plane pi of (a0, a1, b).");
            }
            else if ( (b0_dist > 0.0 && b1_dist > 0.0) || (b0_dist < 0.0 && b1_dist < 0.0) ) {
                //printf("\t both b0 and b1 on the same side of pi => no si.\n");
                return false;
            }
            else if ( (b0_dist > 0.0 && b1_dist < 0.0) || (b0_dist < 0.0 && b1_dist > 0.0) ) {
                //printf("\t b0 and b1 on different sides\n");
                Vec3<R> x;
                R       x_s, x_t, x_lambda;
                /* perform two intersection tests */

                /* edge (a0, a1) with triangle (s, b0, b1) */
                bool a0a1_sb0b1 = rayTriangle<R>(
                                a0, a1,
                                s,
                                b0,
                                b1,
                                x,
                                x_s,
                                x_t,
                                x_lambda) == INTERSECTION;

                /* edge (b0, b1) with triangle (a0, a1, s) */
                bool b0b1_a0a1s = rayTriangle<R>(
                                b0, b1,
                                a0,
                                a1,
                                s,
                                x,
                                x_s,
                                x_t,
                                x_lambda) == INTERSECTION;

                if (a0a1_sb0b1 || b0b1_a0a1s) {
                    printf("interset::triTriSharingEdge(): b0 on one side, b1 on the other. a0a1_sb0b1 || b0b1_a0a1s == true => si\n");
                    return true;
                }
                else return false;
            }
            /* impossible */
            else {
                throw("triTriSharingEdge(): impossible error, b0 and b1 must be somewhere..");
            }
            throw("triTriSharingEdge(): impossible case..");
        }

        /* declaration of moeller function */
        extern "C" {
            int
            tri_tri_intersect(
                double V0[3],
                double V1[3],
                double V2[3],
                double U0[3],
                double U1[3],
                double U2[3]);
        }

        template <typename R>
        bool
        triTri3d(
            Vec3<R> v0,
            Vec3<R> v1,
            Vec3<R> v2,
            Vec3<R> u0,
            Vec3<R> u1,
            Vec3<R> u2)
        {
            double  v0a[3] = { v0[0], v0[1], v0[2] };
            double  v1a[3] = { v1[0], v1[1], v1[2] };
            double  v2a[3] = { v2[0], v2[1], v2[2] };

            double  u0a[3] = { u0[0], u0[1], u0[2] };
            double  u1a[3] = { u1[0], u1[1], u1[2] };
            double  u2a[3] = { u2[0], u2[1], u2[2] };

            int ret = tri_tri_intersect(v0a, v1a, v2a, u0a, u1a, u2a);

            return (ret == 1);
        }

        /* check if a point is inside a simple polygon, given as a vector of vertices in counter
         * clockwise orientation */
        uint32_t    pointInSimplePolygon(
                        std::vector<Vertex2d>   vertices,
                        const Vec2&                    p,
                        double                  eps = 1E-15);

        /* triangulate a given simple planar polygon without holes, given as a vector (by value) of 2d vertices in
         * counter-clockwise order (when projected from 3d space, counter-clockwise when viewed from the
         * outside to uniquely determine orientation. the triangulation is stored in the vector
         * "triangulation" passed by reference. */
        uint32_t    triangulateSimplePlanarPolygon(
                        std::vector<Vertex2d>   vertices,
                        std::vector<Tri2d>     &triangulation);

        /* delaunay triangulation */
        bool        vertexInCircumCircle(Vec2 A, Vec2 B, Vec2 C, Vec2 D);
        bool        delaunayTryFlip(
                        Tri2d                          &A,
                        Tri2d                          &B,
                        std::map<uint32_t, Vertex2d>   &vertices);

        void        delaunay2d(
                        std::map<uint32_t, Vertex2d>   &vertices,
                        std::vector<Tri2d>             &triangles);
                        
        /* NOTE: since std::function / function pointers seems slow, compute bounding boxes once and pass them as
         * arguments */
        template <typename TA, typename TB, typename R = double>
        void
        computeSpatialIntersectionCandidatePairs(
            BoundingBox<R> const                       &bbox,
            std::vector<std::pair<TA, BoundingBox<R>>> &A_list,
            std::vector<std::pair<TB, BoundingBox<R>>> &B_list,
            uint32_t                                    rec_depth,
            uint32_t                                    max_elements,
            uint32_t                                    max_rec_depth,
            std::vector<std::pair<TA, TB>>             &candidate_pairs)
        {
            typedef std::vector<std::pair<TA, BoundingBox<R>>> APairListType;
            typedef std::vector<std::pair<TB, BoundingBox<R>>> BPairListType;

            using namespace Aux::Geometry::IntersectionTestResults;

            debugl(3, "partition_intersect(): depth %4d. A_list.size(): %6ld, B_list.size(): %6ld\n", rec_depth, A_list.size(), B_list.size());
            debugTabInc();

            /* if any of the facelists is empty => there cannot be any intersection */
            if (A_list.empty() || B_list.empty()) {
                debugl(1, "A_list or B_list empty => return\n");
            }
            /* leaf reached or number of components small enough => create leaf and compute intersection */
            else if (rec_depth >= max_rec_depth || A_list.size() + B_list.size() < max_elements) {
                /* check all pairs from A_list and B_list for intersection and insert all candidate
                 * pairs. */
                for (auto &A_tuple : A_list) {
                    TA &A_elem                      = A_tuple.first;
                    BoundingBox<R> const &A_elem_bb = A_tuple.second;

                    for (auto &B_tuple : B_list) {
                        TB &B_elem                      = B_tuple.first;
                        BoundingBox<R> const &B_elem_bb = B_tuple.second;
                        if (A_elem_bb && B_elem_bb) {
                            candidate_pairs.push_back(std::pair<TA, TB>(A_elem, B_elem));
                        }
                    }
                }
                debugl(3, "criterion reached => creating leaf.\n");
            }
            /* partition face list with current cube, recursive call */
            else {
                debugl(3, "max depth not yet reached => partitioning %8ld (A) and %8ld (B) faces among children..\n", A_list.size(), B_list.size());
                uint32_t                                    i;
                Vec3<R>                                     nc_min, nc_max, nc_m;
                Vec3<R>                                     face_bb_min, face_bb_max;
                Vec3<R>                                     xdisp, ydisp, zdisp;

                std::array<BoundingBox<R>, 8> sub_boxes;
                std::array<std::vector<std::pair<TA, BoundingBox<R> > >, 8> A_sub_lists, B_sub_lists;
                std::array<bool, 8> A_sub_box_relevant, B_sub_box_relevant;

                /* allocate all sub-box lists and default sub-boxes to irrelevant for both A and B */
                for (i = 0; i < 8; i++) {
                    A_sub_box_relevant[i]   = false;
                    B_sub_box_relevant[i]   = false;
                }

                /* compute sub-boxes */
                nc_min                  = bbox.min();
                nc_max                  = bbox.max();
                nc_m                    = (nc_min + nc_max) * 0.5;

                xdisp                   = Vec3<R>( (nc_max[0] - nc_min[0]) / 2.0, 0.0, 0.0);
                ydisp                   = Vec3<R>( 0.0, (nc_max[1] - nc_min[1]) / 2.0, 0.0);
                zdisp                   = Vec3<R>( 0.0, 0.0, (nc_max[2] - nc_min[2]) / 2.0);

                /* child 0: n_min and midpoint. offsetting from there with displacement vectors */
                /*
                sub_boxes[0].first      = nc_min;
                sub_boxes[0].second     = nc_m;
                */
                sub_boxes[0]            = BoundingBox<R>(nc_min, nc_m);

                /*
                sub_boxes[1].first      = nc_min + xdisp;
                sub_boxes[1].second     = nc_m + xdisp;
                */
                sub_boxes[1]            = BoundingBox<R>(nc_min + xdisp, nc_m + xdisp);

                /*
                sub_boxes[2].first      = nc_min + xdisp + ydisp;
                sub_boxes[2].second     = nc_m + xdisp + ydisp;
                */
                sub_boxes[2]            = BoundingBox<R>(nc_min + xdisp + ydisp, nc_m + xdisp + ydisp);

                /*
                sub_boxes[3].first      = nc_min + ydisp;
                sub_boxes[3].second     = nc_m + ydisp;
                */
                sub_boxes[3]            = BoundingBox<R>(nc_min + ydisp, nc_m + ydisp);

                /*
                sub_boxes[4].first      = nc_min + zdisp;
                sub_boxes[4].second     = nc_m + zdisp;
                */
                sub_boxes[4]            = BoundingBox<R>(nc_min + zdisp, nc_m + zdisp);

                /*
                sub_boxes[5].first      = nc_min + xdisp + zdisp;
                sub_boxes[5].second     = nc_m + xdisp + zdisp;
                */
                sub_boxes[5]            = BoundingBox<R>(nc_min + xdisp + zdisp, nc_m + xdisp + zdisp);

                /* specesial case: child 6 has min corner m and max corner nc_max */
                /*
                sub_boxes[6].first      = nc_m;
                sub_boxes[6].second     = nc_max;
                */
                sub_boxes[6]            = BoundingBox<R>(nc_m, nc_max);

                /*
                sub_boxes[7].first      = nc_min + ydisp + zdisp;
                sub_boxes[7].second     = nc_m + ydisp + zdisp;
                */
                sub_boxes[7]            = BoundingBox<R>(nc_min + ydisp + zdisp, nc_m + ydisp + zdisp);


                /* iterate over all faces in A_list and partition them into the sub-box lists.  if child i gets a face of
                 * A, set sub_box_relevant_A = true. same is done for B.  child i is relevant for the intersection iff
                 *
                 *  A_sub_box_relevant[i] && B_sub_box_relevant[i]
                 *  
                 *  */
                for (auto A_tuple : A_list) {
                    /* check for intersections with all eight sub-boxes. */
                    for (i = 0; i < 8; i++) {
                        /* get bounding box of element, check for intersection. */
                        BoundingBox<R> const &A_elem_bb = A_tuple.second;
                        if (A_elem_bb && sub_boxes[i]) {
                            A_sub_lists[i].push_back(A_tuple);
                            A_sub_box_relevant[i] = true;
                        }
                    }
                }

                /* same for B_list */
                for (auto B_tuple : B_list) {
                    /* check for intersections with all eight sub-boxes. */
                    for (i = 0; i < 8; i++) {
                        BoundingBox<R> const &B_elem_bb = B_tuple.second;
                        if (B_elem_bb && sub_boxes[i]) {
                            B_sub_lists[i].push_back(B_tuple);
                            B_sub_box_relevant[i] = true;
                        }
                    }
                }

                /* clear facelists from current call, and free the memory */
                APairListType().swap(A_list);
                BPairListType().swap(B_list);

                /* recursive calls for sub-boxes relevant to the intersection. */
                for (i = 0; i < 8; i++) {
                    if (A_sub_box_relevant[i] && B_sub_box_relevant[i]) {
                        /* recursive call */
                        computeSpatialIntersectionCandidatePairs(
                                sub_boxes[i],
                                A_sub_lists[i], B_sub_lists[i],
                                rec_depth + 1, max_elements, max_rec_depth,
                                candidate_pairs);
                    }
                }
            }
            debugTabDec();
            debugl(3, "all recursive calls finished => returning..\n");
        }

        /* version with std::function / function pointer. std::function version was extremely slow during testing.. */
        template <typename TA, typename TB, typename R = double>
        void
        computeSpatialIntersectionCandidatePairsWithFunctor(
            std::pair<Vec3<R>, Vec3<R>>                                     bbox,
            std::list<TA>                                                  &A_list,
            std::list<TB>                                                  &B_list,
            /* NOTE: std::function overhead seems extreme in this case, use function pointers for comparison.. */
            /*
            std::function<std::pair<Vec3<R>, Vec3<R>>(const TA &x)> const  &getBoundingBoxTA,
            std::function<std::pair<Vec3<R>, Vec3<R>>(const TB &x)> const  &getBoundingBoxTB,
            */
            BoundingBox<R>                                               (*getBoundingBoxTA)(const TA &x),
            BoundingBox<R>                                               (*getBoundingBoxTB)(const TB &x),
            uint32_t                                                        rec_depth,
            uint32_t                                                        max_elements,
            uint32_t                                                        max_rec_depth,
            std::list< std::pair<TA, TB> >                                 &candidate_pairs)
        {
            using namespace Aux::Geometry::IntersectionTestResults;

            debugl(3, "partition_intersect(): depth %4d. A_list.size(): %6ld, B_list.size(): %6ld\n", rec_depth, A_list.size(), B_list.size());
            debugTabInc();

            Vec3<R>                     A_elem_bb_min, A_elem_bb_max, B_elem_bb_min, B_elem_bb_max;
            std::pair<Vec3<R>, Vec3<R>> A_bb_pair, B_bb_pair;

            /* if any of the facelists is empty => there cannot be any intersection */
            if (A_list.empty() || B_list.empty()) {
                debugl(1, "A_list or B_list empty => return\n");
            }
            /* leaf reached or number of components small enough => create leaf and compute intersection */
            else if (rec_depth >= max_rec_depth || A_list.size() + B_list.size() < max_elements) {
                /* check all pairs from A_list and B_list for intersection and insert all candidate
                 * pairs. */
                for (auto &A_elem : A_list) {
                    auto A_bb = getBoundingBoxTA(A_elem);

                    for (auto &B_elem : B_list) {
                        auto B_bb = getBoundingBoxTB(B_elem);
                        if (A_bb && B_bb) {
                            candidate_pairs.push_back(std::pair<TA, TB>(A_elem, B_elem));
                        }
                    }
                }
                debugl(3, "criterion reached => creating leaf.\n");
            }
            /* partition face list with current cube, recursive call */
            else {
                debugl(3, "max depth not yet reached => partitioning %8ld (A) and %8ld (B) faces among children..\n", A_list.size(), B_list.size());
                uint32_t                        i;
                Vec3<R>                         nc_min, nc_max, nc_m;
                Vec3<R>                         face_bb_min, face_bb_max;
                Vec3<R>                         xdisp, ydisp, zdisp;

                std::array<BoundingBox<R>, 8>               sub_boxes;
                std::array<
                    std::list<
                            std::pair<
                                TA,
                                BoundingBox<R>
                            >
                        >,
                    8>                                      A_sub_lists, B_sub_lists;
                
                std::array<bool, 8>                         A_sub_box_relevant, B_sub_box_relevant;

                /* allocate all sub-box lists and default sub-boxes to irrelevant for both A and B */
                for (i = 0; i < 8; i++) {
                    A_sub_box_relevant[i]   = false;
                    B_sub_box_relevant[i]   = false;
                    A_sub_lists[i]          = {};
                    B_sub_lists[i]          = {};
                }

                /* compute sub-boxes */
                nc_min                  = bbox.min();
                nc_max                  = bbox.max();
                nc_m                    = (nc_min + nc_max) * 0.5;

                xdisp                   = Vec3<R>( (nc_max[0] - nc_min[0]) / 2.0, 0.0, 0.0);
                ydisp                   = Vec3<R>( 0.0, (nc_max[1] - nc_min[1]) / 2.0, 0.0);
                zdisp                   = Vec3<R>( 0.0, 0.0, (nc_max[2] - nc_min[2]) / 2.0);

                /* child 0: n_min and midpoint. offsetting from there with displacement vectors */
                /*
                sub_boxes[0].first      = nc_min;
                sub_boxes[0].second     = nc_m;
                */
                sub_boxes[0]            = BoundingBox<R>(nc_min, nc_m);

                /*
                sub_boxes[1].first      = nc_min + xdisp;
                sub_boxes[1].second     = nc_m + xdisp;
                */
                sub_boxes[1]            = BoundingBox<R>(nc_min + xdisp, nc_m + xdisp);

                /*
                sub_boxes[2].first      = nc_min + xdisp + ydisp;
                sub_boxes[2].second     = nc_m + xdisp + ydisp;
                */
                sub_boxes[2]            = BoundingBox<R>(nc_min + xdisp + ydisp, nc_m + xdisp + ydisp);

                /*
                sub_boxes[3].first      = nc_min + ydisp;
                sub_boxes[3].second     = nc_m + ydisp;
                */
                sub_boxes[3]            = BoundingBox<R>(nc_min + ydisp, nc_m + ydisp);

                /*
                sub_boxes[4].first      = nc_min + zdisp;
                sub_boxes[4].second     = nc_m + zdisp;
                */
                sub_boxes[4]            = BoundingBox<R>(nc_min + zdisp, nc_m + zdisp);

                /*
                sub_boxes[5].first      = nc_min + xdisp + zdisp;
                sub_boxes[5].second     = nc_m + xdisp + zdisp;
                */
                sub_boxes[5]            = BoundingBox<R>(nc_min + xdisp + zdisp, nc_m + xdisp + zdisp);

                /* specesial case: child 6 has min corner m and max corner nc_max */
                /*
                sub_boxes[6].first      = nc_m;
                sub_boxes[6].second     = nc_max;
                */
                sub_boxes[6]            = BoundingBox<R>(nc_m, nc_max);

                /*
                sub_boxes[7].first      = nc_min + ydisp + zdisp;
                sub_boxes[7].second     = nc_m + ydisp + zdisp;
                */
                sub_boxes[7]            = BoundingBox<R>(nc_min + ydisp + zdisp, nc_m + ydisp + zdisp);


                /* iterate over all faces in A_list and partition them into the sub-box lists.  if child i gets a face of
                 * A, set sub_box_relevant_A = true. same is done for B.  child i is relevant for the intersection iff
                 *
                 *  A_sub_box_relevant[i] && B_sub_box_relevant[i]
                 *  
                 *  */
                for (auto A_tuple : A_list) {
                    /* check for intersections with all eight sub-boxes. */
                    for (i = 0; i < 8; i++) {
                        /* get bounding box of element, check for intersection. */
                        BoundingBox<R> const &A_elem_bb = A_tuple.second;
                        if (A_elem_bb && sub_boxes[i]) {
                            A_sub_lists[i]->push_back(A_tuple);
                            A_sub_box_relevant[i] = true;
                        }
                    }
                }

                /* same for B_list */
                for (auto B_tuple : B_list) {
                    /* check for intersections with all eight sub-boxes. */
                    for (i = 0; i < 8; i++) {
                        BoundingBox<R> const &B_elem_bb = B_tuple.second;
                        if (B_elem_bb && sub_boxes[i]) {
                            B_sub_lists[i]->push_back(B_tuple);
                            B_sub_box_relevant[i] = true;
                        }
                    }
                }

                /* clear facelists from current call, they are not needed anymore */
                A_list.clear();
                B_list.clear();

                /* recursive calls for sub-boxes relevant to the intersection. */
                for (i = 0; i < 8; i++) {
                    if (A_sub_box_relevant[i] && B_sub_box_relevant[i]) {
                        /* recursive call */
                        computeSpatialIntersectionCandidatePairs(
                                sub_boxes[i],
                                *(A_sub_lists[i]), *(B_sub_lists[i]),
                                rec_depth + 1, max_elements, max_rec_depth,
                                candidate_pairs);
                    }

                    /* delete allocated face sub-lists */
                    delete A_sub_lists[i];
                    delete B_sub_lists[i];
                }
            }

            debugTabDec();
            debugl(3, "all recursive calls finished => returning..\n");
        }
    }

    namespace Alg {
        uint32_t
        stou(
            std::string const   &str,
            size_t             *idx     = 0,
            int                 base    = 10);

        template<typename T>
        bool
        vectorContains(const std::vector<T> &vec, const T &x)
        {
            return ( std::find(vec.begin(), vec.end(), x) != vec.end() );
        }

        template<typename T>
        typename std::list<T>::iterator
        listFind(
            std::list<T>   &l,
            const T        &x)
        {
            for (auto lit = l.begin(); lit != l.end(); ++lit) {
                if (x == *lit) {
                    return lit;
                }
            }
            return l.end();
        }

        template<typename T>
        bool
        listContains(
            const std::list<T> &l,
            const T            &x)
        {
            for (auto lit = l.begin(); lit != l.end(); ++lit) {
                if (x == *lit) {
                    return true;
                }
            }
            return false;
        }

        // EDIT: added pointer-* to parameter type of x and compare IDs rather than memory addresses
        // as this function is used exclusively in the case where x is a pointer to a geometric object
        template<typename T>
        bool
        listSortedInsert(
            std::list<T*>    &l,
            T*               x,
            bool             duplicates = false)
        {
            auto    lit = l.begin();
            while (lit != l.end() && (*lit)->id() < x->id()) {
                ++lit;
            }

            if (lit != l.end() && !duplicates && (*lit)->id() == x->id()) {
                return false;
            }
            else {
                l.insert(lit, x);
                return true;
            }
        }

        template<typename T>
        uint32_t
        listRemove(std::list<T> &l, const T &x)
        {
            typename std::list<T>::iterator lit;
            uint32_t nerased = 0;

            lit = l.begin();
            while (lit != l.end()) {
                if (x == *lit) {
                    lit = l.erase(lit);
                    ++nerased;
                }
                else ++lit;
            }
            return nerased;
        }

        template<template <typename, typename> class ListType, typename T, typename TAlloc>
        bool
        removeFirstOccurrenceFromList(ListType<T, TAlloc> &l, const T &x)
        {
            typename ListType<T, TAlloc>::iterator lit;

            for (lit = l.begin(); lit != l.end(); ++lit) {
                if (x == *lit) {
                    l.erase(lit);
                    return true;
                }
            }
            return false;
        }

        template<typename T>
        uint32_t
        listIntersection(std::list<T> l1, std::list<T> l2, std::list<T> &isec)
        {
            l1.sort();
            l1.unique();
            l2.sort();
            l2.unique();

            isec.clear();
            std::set_intersection(l1.begin(), l1.end(), l2.begin(), l2.end(), std::back_inserter(isec));

            return isec.size();
        }

        inline double
        sumReduce(const std::vector<double>  &values)
        {
            double  sum = 0.0;

            for (double v : values) {
                sum += v;
            }

            return sum;
        }

        /* generic binary searching for index on vectors, since C++ is so goddamn ugly it doesn't give you a efficient
         * way to do that. it WILL give you the iterator, but it won't give you its index without packing it
         * inside the key structure, which is just ugly */
        // And this is so goddamn important that it is never once used here!
        template <typename T>
        bool
        binarySearchVector(std::vector<T> v, T x, uint32_t &x_idx)
        {
            uint32_t left = 0, right = v.size() - 1;

            while (left <= right) {
                uint32_t middle = (left + right) / 2;
                if (v[middle] == x) {
                    x_idx = middle;
                    return true;
                }
                else if (x < v[middle]) {
                    right = middle - 1;     
                }
                else {
                    left = middle + 1;
                }
            }
            return false;
        }

        namespace Functional {
            template<
                typename    alpha,
                typename    beta
            >
            void
            map(
                std::list<alpha>       &list,
                beta                  (*op)(alpha const &x),
                std::list<beta>        &result)
            {
                result.resize(list.size());
                std::transform(list.begin(), list.end(), result.begin(), op);
            }

            template<
                typename    alpha,
                typename    beta
            >
            std::list<beta>
            map(
                std::list<alpha>       &list,
                beta                  (*op)(alpha const &x))
            {
                std::list<beta> result;
                Aux::Alg::Functional::map(list, op, result);
                return result;
            }

            template<
                typename alpha
            >
            void
            filter(
                std::list<alpha>       &list,
                bool                  (*pred)(const alpha &x))
            {
                /* remove_if predicate is FALSE => the container contains ever element for which
                 * predicate is TRUE after return. this is equivalent to usual filter() semantics. */
                auto it = std::remove_if(list.begin(), list.end(), [pred] (const alpha &x) -> bool { return !pred(x); } );
                list.erase(it, list.end());
            }

            template<
                typename alpha,
                typename beta
            >
            alpha
            foldl( 
                alpha             (*op)(const alpha &x, const beta &y),
                const alpha        &e,
                std::list<beta>    &list)
            {
                return std::accumulate(list.begin(), list.end(), e, op);
            }
                    
            template<
                typename alpha,
                typename beta,
                typename Ain, 
                typename Aout,
                template<typename, typename> class Cin,
                template<typename, typename> class Cout
            >
            void
            gMap(
                const Cin<alpha, Ain>  &container,
                beta                  (*op)(alpha const &x),
                Cout<beta, Aout>       &result)
            {
                std::transform(container.begin(), container.end(), result.begin(), op);
            }

            template<
                typename alpha,
                typename beta,
                typename Ain, 
                typename Aout,
                template<typename, typename> class Cin,
                template<typename, typename> class Cout
            >
            Cout<beta, Aout>
            gMap(
                const Cin<alpha, Ain>  &container,
                beta                  (*op)(alpha const &x))
            {
                Cout<beta, Aout> result(container.size());
                Aux::Alg::Functional::map<alpha, beta, Ain, Aout, Cin, Cout>(container, op, result);
                return result;
            }

            template<
                typename alpha,
                typename A,
                template<typename, typename> class C
            >
            void
            gFilter(
                C<alpha, A>            &container,
                bool                  (*pred)(const alpha &x))
            {
                /* remove_if predicate is FALSE => the container contains ever element for which
                 * predicate is TRUE after return. this is equivalent to usual filter() semantics. */
                auto it = std::remove_if(container.begin(), container.end(), [pred] (const alpha &x) -> bool { return !pred(x); } );
                container.erase(it, container.end());
            }

            template<
                typename beta,
                typename alpha,
                typename A = std::allocator<alpha>,
                template<typename, typename> class C
            >
            beta
            gFoldl(
                C<alpha, A>            &container,
                beta                  (*op)(const beta &x, const alpha &y),
                beta                    e)
            {
                return std::accumulate(container.begin(), container.end(), e, op);
            }
        }
    }

    namespace Stat {
        template <typename R>
        void
        computeMinMaxAvgSigma(
                std::vector<R> const       &values,
                R                          &min,
                R                          &max,
                R                          &avg,
                R                          &sigma)
        {
            R   n   = (R)values.size();
            R   sum = 0.0;

            /* set min = INF, max = -INF. note that Aux::Numbers::inf as well as std::{min, max} must have been
             * specialized for template type R. */
            min = Aux::Numbers::inf<R>();
            max = -Aux::Numbers::inf<R>();

            for (R x : values) {
                min     = std::min(min, x);
                max     = std::max(min, x);
                sum    += x;
            }
            avg = sum / n;

            /* sigma with bessel correction */
            R var = 0.0;

            for (R x : values) {
                var += (x - avg)*(x - avg);
            }

            sigma = std::sqrt( var / (n-1.0) );
        }
    }


    namespace Logic {
        inline bool
        lnot(bool a)
        {
            return (!a);
        }

        inline bool
        land(bool a, bool b)
        {
            return ( a && b);
        }

        inline bool
        lor(bool a, bool b)
        {
            return ( a || b);
        }

        inline bool
        lnor(bool a, bool b)
        {
            return !(a || b);
        }

        inline bool
        lnand(bool a, bool b)
        {
            return !(a && b);
        }

        inline bool
        lxor(bool a, bool b)
        {
            return ( (!a && b) || (a && !b) );
        }

        /* xnor == equivalence */
        inline bool
        lequiv(bool a, bool b)
        {
            return ( (a && b) || (!a && !b) );
        }

        inline bool
        getBit(uint32_t &flags, uint8_t bitpos)
        {
            return (flags & (1 << bitpos));
        }

        inline void
        setBit(uint32_t &flags, uint8_t bitpos, bool value)
        {
            if (value) {
                flags |= (1 << bitpos);
            }
            else {
                flags &= ~(1 << bitpos);
            }
        }

        inline uint32_t
        getBits(uint32_t &flags, uint8_t plow, uint8_t phigh)
        {
            uint32_t tmp = flags >> plow;
            tmp &= ( (1 << (phigh - plow + 1)) - 1);
            //printf("returning tmp = %d\n", tmp);
            return tmp;
        }

        inline void
        setBits(uint32_t &flags, uint8_t plow, uint8_t phigh, uint32_t value)
        {
            /* first, zero the range [plow..phigh]. to do that, generate (2^(phigh - plow + 1) - 1), which has the
             * first phigh bits set to 1. left shift that by plow, take one's complement and AND it over
             * flags */
            if (plow <= phigh) {
                //printf("value = %d\n", value);
                //uint32_t tmp = flags;
                flags &= ~( ( (1 << (phigh - plow + 1)) - 1) << plow);
                //printf("phigh: %d, plow: %d, 2^(phigh - plow + 1) - 1 : %d\n", phigh, plow, ( 1 << (phigh - plow + 1) ) - 1);

                /* then overwrite with value bits, very similar to resetting as done above. */
                flags |= (value & ( (1 << (phigh - plow + 1)) - 1)) << plow;
                //printf("flags after - flags before = %d\n", flags - tmp);
            }
            else throw("setBit(): range plow > phigh given.\n");
        }
    }
}

#endif
