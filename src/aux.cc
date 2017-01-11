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

#include <stdarg.h>
#include "debug.hh"
#include "Vec3.hh"
#include "aux.hh"

namespace Aux {
    namespace Timing {
        struct timeval 
            starttimes[TIMER_REGISTERS], 
            endtimes[TIMER_REGISTERS];

        /* timing functions. */

        /* tick: start time measurement for register i */
        void tick(int i) {
            gettimeofday(&starttimes[i], NULL);
        }

        /* tack: stop time measurement and return value for register i */
        double tack(int i) {
            double starttime, endtime;

            gettimeofday(&endtimes[i], NULL);

            starttime 		= starttimes[i].tv_sec 	+(starttimes[i].tv_usec/1000000.0);
            endtime 		= endtimes[i].tv_sec 	+(endtimes[i].tv_usec/1000000.0);

            return (endtime - starttime);
        }

        double doubletime()
        {
            struct timeval  time;
            gettimeofday(&time, NULL);
            return ( (double)time.tv_sec + (double)time.tv_usec / 1000000.0 );
        }
    }

    namespace Numbers {
        /* random double */
        double
        frand(double min, double max)
        {
            double f = (double)std::rand() / RAND_MAX;
            std::cout << " (ok)" << std::endl;
            return (min + f*(max - min));
        }

        double
        fmin3(double a, double b, double c)
        {
            return ( fmin(a, fmin(b, c)) );
        }

        double
        fmax3(double a, double b, double c)
        {
            return ( fmax(a, fmax(b, c)) );
        }

        int
        sign(int64_t d)
        {
            return ( (d > 0) - (d < 0) );
        }

        double 
        deg2rad(double angle)
        {
            return (M_PI * (angle / 180.0));
        }

        double 
        rad2deg(double angle)
        {
            return (180 * (angle / M_PI));
        }

        /* specializations of template inf<T>() for float, double and long double */
        template<>
        float
        inf() {
            return std::numeric_limits<float>::infinity();
        }

        template<>
        double
        inf() {
            return std::numeric_limits<double>::infinity();
        }

        template<>
        long double
        inf() {
            return std::numeric_limits<long double>::infinity();
        }


    }

    namespace Alg {
        uint32_t
        stou(
            std::string const   &str,
            size_t             *idx,
            int                 base)
        {
            unsigned long result = std::stoul(str, idx, base);
            if (result > std::numeric_limits<uint32_t>::max()) {
                throw std::out_of_range("stou");
            }
            return result;
        }
    }

    namespace File {
        bool isEmpty(FILE *f)
        {
            long offset = ftell(f);
            fseek(f, 0, SEEK_END);

            if (ftell(f) == 0)
            {
                return true;
            }
            else {
                fseek(f, offset, SEEK_SET);
                return false;
            }
        }
    }

    namespace Geometry {
        using namespace Aux::Geometry::IntersectionTestResults;

        Tri2d::Tri2d(
            uint32_t v0_id,
            uint32_t v1_id,
            uint32_t v2_id)
        {
            this->v0_id = v0_id;
            this->v1_id = v1_id;
            this->v2_id = v2_id;
        }


        Vertex2d::Vertex2d() 
        {
            this->id        = UINT32_MAX;
            this->pos[0]    = Aux::Numbers::inf<double>();
            this->pos[1]    = Aux::Numbers::inf<double>();
        }

        Vertex2d::Vertex2d(uint32_t _id, const Vec2& _pos)
        : id(_id), pos(_pos)
        {}


        uint32_t
        lineSegmentLineSegment2d(
            const Vec2&        p0,
            const Vec2&        p1,
            const Vec2&        q0,
            const Vec2&        q1,
            Vec2       &x,
            double     &x_lambda,
            double     &x_mu,
            double      eps)
        {
            /* get direction vectors: u = p1 - p0, v = q1 - q0 */
            Vec2    u, v, w;
            double  denom;

            u       = p1 - p0;
            v       = q1 - q0;
            w       = p0 - q0;
            denom   = u[0]*v[1] - u[1]*v[0];

            /* if lines are parallel, return EDGE_CASE for now, we need finite, non-coincident intersections */
            if ( fabs(denom) < eps) {
                return EDGE_CASE;
            }
            else {
                debugTabInc();

                double lambda_isec = (v[0]*w[1] - v[1]*w[0]) / denom;
                double mu_isec     = (u[0]*w[1] - u[1]*w[0]) / denom;

                x           = p0 + u*lambda_isec;
                x_lambda    = lambda_isec;
                x_mu        = mu_isec;

                if (lambda_isec < 0 || lambda_isec > 1.0 || mu_isec < 0 || mu_isec > 1.0) {
                    debugl(0, "lambda_isec: %5.4f, mu_isec: %5.4f for inifinite lines out of range.\n", lambda_isec, mu_isec);
                    debugTabDec();
                    return DISJOINT;
                }
                else {
                    debugl(0, "lambda_isec: %5.4f, mu_isec: %5.4f within range of finite segments => interseciton.\n", lambda_isec, mu_isec);
                    debugTabDec();
                    return INTERSECTION;
                }
            }
        }

        uint32_t
        rayLineSegment2d(
            const Vec2&        p,
            const Vec2&        dir,
            const Vec2&        q0,
            const Vec2&        q1,
            Vec2       &x,
            double     &x_lambda,
            double     &x_mu,
            double      eps)
        {
            /* get direction vectors: u = dir, v = q1 - q0 */
            Vec2    u, v, w;
            double  denom;

            u       = dir;
            v       = q1 - q0;
            w       = p - q0;
            denom   = u[0]*v[1] - u[1]*v[0];

            if ( fabs(denom) < eps) {
                return EDGE_CASE;
            }
            else {
                debugTabInc();

                x_lambda = (v[0]*w[1] - v[1]*w[0]) / denom;
                x_mu     = (u[0]*w[1] - u[1]*w[0]) / denom;

                if (x_mu < -eps || x_mu > 1.0 + eps) {
                    debugl(0, "rayLineSegment2d(): x_lambda: %5.4f (ray), x_mu: %5.4f (line segment) out of range.\n", x_lambda, x_mu);

                    debugTabDec();
                    return DISJOINT;
                }
                else {
                    debugl(0, "x_lambda: %5.4f (ray), x_mu: %5.4f (line segment) within range => interseciton.\n", x_lambda, x_mu);
                    x           = p + u*x_lambda;
                    
                    debugTabDec();
                    return INTERSECTION;
                }
            }
        }

        /* ray casting with random unit 2d vector. using a simple parity argument, point containment can be
         * decided efficiently */
        uint32_t
        pointInSimplePolygon(
            std::vector<Vertex2d>   vertices,
            const Vec2&             p,
            double                  eps)
        {
            Vec2                    d;
            uint32_t                i, attempts;
            std::vector<uint32_t>   isec_edges;
            Vec2                    x;
            Vec2                    v_i, v_ipo;
            double                  x_lambda, x_mu;

            /* append first vertex to vertex list to avoid wrap around case */
            vertices.push_back(vertices.front());

            debugl(0, "pointInSimplePolygon().\n");

            attempts    = 0;

            debugTabInc();
            while (attempts <= 1024) {
                debugl(0, "another attempt..");
                bool retry = false;

                /* fresh attempt */
                attempts++;
                isec_edges.clear();

                /* pick a random 2d vector d and cast the ray p + lambda*d, intersect it with
                 * all edges of the polygon. when close to an edge case, repeat with another random vector */
                std::cout << "    Calling rand for intersect direction (" << std::rand() << ")";
                d = Aux::VecMat::randUnitVec2();

                /* intersect all edges, interpreted as line segments in 2d. wrap around case (last, first)
                 * is handled by duplicating first vertex and appending it to the back of the list */
                for (i = 0; i < vertices.size() - 1; i++) {
                    v_i     = vertices[i].pos;
                    v_ipo   = vertices[i+1].pos;

                    uint32_t ray_line_result = rayLineSegment2d(p, d, v_i, v_ipo, x, x_lambda, x_mu);
                    if (ray_line_result == INTERSECTION) {
                        /* if lambda in [0,1] is smaller than eps in magnitude, the point lies on the edge =>
                         * edge case */
                        if (fabs(x_lambda) < eps) {
                            debugTabDec();
                            return EDGE_CASE;
                        }

                        /* if x is too close to the either one of the vertices, again.. */
                        if ( (x - v_i).len2squared() < 1E-8 || (x - v_ipo).len2squared() < 1E-8) {
                            debugTabDec();
                            return EDGE_CASE;
                        }
                        /* edge i = (v_i, v_{i+1}) is intersected with positive lambda (edge case caught
                         * above) */
                        else if (x_lambda > 0.0) {
                            debugl(0, "got one..\n");
                            isec_edges.push_back(i);
                        }
                    }
                    else if (ray_line_result == EDGE_CASE) {
                        /* parallel lines. retry */
                        retry = true;
                        break;
                    }
                }

                if (retry) {
                    continue;
                }

                /* if there are intersections, check if isec_edges contains consecutive edge indices */
                if (isec_edges.size() > 0) {
                    for (i = 0; i < isec_edges.size() - 1; i++) {
                        if (isec_edges[i] == isec_edges[i+1]) {
                            retry = true;
                        }
                    }
                }

                if (retry) {
                    continue;
                }

                /* point is inside iff number of intersections is odd */
                if (isec_edges.size() % 2 == 1) {
                    debugTabDec();
                    return POINT_IN;
                }
                else {
                    debugTabDec();
                    return POINT_OUT;
                }
            }
            debugTabDec();

            /* max number of attempts reached: indicate inconslusive test result */
            //printf("WARNING: pointInSimplePolygon(): maximum number of inconclusive attempts reached => indicating inconclusive result.\n");
            return TEST_INCONCLUSIVE;
        }

        uint32_t
        triangulateSimplePlanarPolygon(
            std::vector<Vertex2d>   vertices, 
            std::vector<Tri2d>     &triangulation)
        {
            uint32_t                        maxiter = 1024;
            uint32_t                        i, k, n = vertices.size();
            bool                            k_principal, got_ear = false;
            uint32_t                        u_id, v_id, w_id, x_id, y_id;
            Vec2                            u, v, w, x, y, p_isec, u_v_edgepoint, pisec;
            uint32_t                        pointcheck_result;
            double                          pisec_lambda, pisec_mu;
            std::vector<Vertex2d>           vertices_copy;
            std::vector<Vertex2d>::iterator vit;

            /* clear return list */
            triangulation.clear();

            if (n < 3) {
                throw("Aux::Geometry::triangulateSimplePlanarPolygon(): number of vertices < 3 => no polygon.");
            }
            else if (n == 3) {
                triangulation.push_back( Tri2d(vertices[0].id, vertices[1].id, vertices[2].id) );
                return SUCCESS;
            }
            else {
                /* find an ear, clip it, add triangle. repeat until vertices list is trivial, i.e. only one
                 * triangle remains */
                debugTabInc();
                while ( (n = vertices.size()) > 3) {
                    debugl(0, "triangulateSimplePlanarPolygon(): %ld vertices remaining.\n", n);

                    /* check all vertices 0..n-1 for principal property. to make the wrap-around issue
                     * easier, insert vertex 0 and 1 at the end, which results in vector of size n + 1.
                     * and loop with k from 1..n and check edge (u, v) = ( id(k-1), id(k+1) ) for intersections
                     * with all other edges. use the ids stored in the Vertex2d structs to filter out trivial intersections
                     * with any edge containing the id(k-1) and id(k+1), since they will always intersect in
                     * the common vertex. backup vertices first for later use. */
                    vertices_copy = vertices;
                    vertices.push_back(vertices[0]);
                    vertices.push_back(vertices[1]);

                    for (k = 1; k < n + 1; k++) {
                        /* retrieve ids */
                        u_id    = vertices[k-1].id;
                        u       = vertices[k-1].pos;

                        w_id    = vertices[k].id;
                        w       = vertices[k].pos;

                        v_id    = vertices[k+1].id;
                        v       = vertices[k+1].pos;

                        debugl(0, "checking vertex %5d for principal property by examining edge (%5d, %5d) of triangle (%5d, %5d, %5d) for intersections.\n",
                                w_id, u_id, v_id, u_id, w_id, v_id);

                        /* default to principal, if there's an intersection, this will be set to false.
                         * default to no ear, this will be teste below, if w is indeed a principal vertex */
                        k_principal = true;
                        got_ear     = false;

                        debugTabInc();

                        /* check for intersection with all edges (x_id, y_id) = ( vertices[i].id, vertices[i+1].id ) */
                        for (i = 0; i < n - 1; i++) {
                            x_id    = vertices[i].id;
                            x       = vertices[i].pos;

                            y_id    = vertices[i+1].id;
                            y       = vertices[i+1].pos;

                            /* filter out trivial intersections */
                            if (! (x_id == u_id || x_id == v_id || y_id == u_id || y_id == v_id) ) {
                                /* check */
                                uint32_t line_isec_result = lineSegmentLineSegment2d(u, v, x, y, pisec, pisec_lambda, pisec_mu);
                                if (line_isec_result == INTERSECTION) {
                                    k_principal = false;
                                    break;
                                }
                                else if (line_isec_result == EDGE_CASE) {
                                    return EDGE_CASE;
                                }
                            }
                        }

                        /* if vertex at position is a principal vertex, check if its a mouth or ear: it's an ear iff any
                         * point on the edge (k-1, k+1) is INSIDE the polygon (otherwise a mouth). since the
                         * edge (k-1, k+1) doesn't intersect any other edge, only one point must be checked,
                         * which is chosen randomly to avoid edge cases. */
                        if (k_principal) {
                            debugl(0, "got principal vertex %5d..\n", w_id);
                            /* if k is an ear, pop_back twice to get rid of copies, erase vertex k, insert
                             * triangle (k-1, k, k+1) into triangle return vector. note that vertices[0] is
                             * checked when k == n => catch that special case */

                            uint32_t iter = 0;
                            debugTabInc();
                            while (++iter < maxiter) {
                                debugl(0, "checking whether principal vertex is mouth or ear..\n");
                                std::cout << "    Calling rand in polygon triangulation (" << std::rand() << ")";
                                u_v_edgepoint       = u + (v - u) * Aux::Numbers::frand(0.1, 0.9);
                                pointcheck_result   = pointInSimplePolygon(vertices_copy, u_v_edgepoint);

                                /* point is in => w is an ear => add triangle (u, w, v) to triangulation and
                                 * delete vertex for w from vertices array */
                                if (pointcheck_result == POINT_IN) {
                                    debugl(0, "got an ear!\n");
                                    got_ear     = true;

                                    /* copy vertices from backup */
                                    vertices    = vertices_copy;

                                    /* add triangle (u, w, v) triangulation */
                                    triangulation.push_back( Tri2d(u_id, w_id, v_id) );

                                    /* special case if k == n -> set k = 0 */
                                    if (k == n) {
                                        k = 0;
                                    }

                                    /* delete vertex k from vertices. */
                                    vertices.erase(vertices.begin() + k);
                                    break;
                                }
                                /* mouth => keep on going */
                                else if (pointcheck_result == POINT_OUT) {
                                    debugl(0, "got a mouth.. :()\n");
                                    break;
                                }
                                /* edge case => this should not happen, but retry.. */
                                else if (pointcheck_result == EDGE_CASE) {
                                    debugl(0, "WARNING: triangulateSimplePlanarPolygon(): edge case during point containment check.\n");
                                }
                                /* inconclusive test, retry.. */
                                else if (pointcheck_result == TEST_INCONCLUSIVE) {
                                    debugl(0, "point containment check inconclusive => retrying..\n");
                                }
                                else {
                                    throw("Aux::Geometry::triangulateSimplePlanarPolygon(): point containment check returned invalid result. internal logic error.");
                                }
                            }
                            debugTabDec();

                            /* maxiter iterations have been performed without leaving the loop via return or exception
                             * => indicate inconclusive result after maximum number of iterations. */
                            if (iter == maxiter) {
                                debugTabDec(); debugTabDec();
                                return TEST_INCONCLUSIVE;
                            }
                        }
                        else {
                            debugl(0, "not principal..\n");
                        }

                        debugTabDec();

                        /* break if search was successful */
                        if (got_ear) break;
                    }

                    if (!got_ear) {
                        /* check if we have reached this point without finding an ear */
                        //printf("WARNING: triangulateSimplePlanarPolygon(): more than three vertices remaining, yet no ear left => simple / without holes?\n");
                        debugTabDec();
                        return EDGE_CASE;
                    }
                }
                debugTabDec();

                /* vertex.size() == 3 now. add the remaining triangle */
                debugl(0, "3 vertices remaining. adding final triangle..\n");
                triangulation.push_back( Tri2d(vertices[0].id, vertices[1].id, vertices[2].id) );
                return SUCCESS;
            }
        }

        /* 2d delaunay triangulation.
         * given a vector of 2d points, a vector of indices and a precomputed (usually bad) triangulation,
         * compute the delaunay triangulation, i.e. flip edges until every triangle is locally delaunay,
         * which means that the maximum minimum angle is maximixez.. o_O. the number of triangles will
         * generally be quite small (<< 100), so that a more sophisticated implementation seems overkill. */

        /* check if 2d vertex v is in circumcicle of triangle T */
        bool
        vertexInCircumCircle(Vec2 A, Vec2 B, Vec2 C, Vec2 D)
        {
            double  M[3][3];

            double  Dx_sq, Dy_sq;

            Dx_sq   = D[0]*D[0];
            Dy_sq   = D[1]*D[1];

            M[0][0] = A[0] - D[0];
            M[0][1] = A[1] - D[1];
            M[0][2] = (A[0]*A[0] - Dx_sq) + (A[1]*A[1] - Dy_sq);

            M[1][0] = B[0] - D[0];
            M[1][1] = B[1] - D[1];
            M[1][2] = (B[0]*B[0] - Dx_sq) + (B[1]*B[1] - Dy_sq);

            M[2][0] = C[0] - D[0];
            M[2][1] = C[1] - D[1];
            M[2][2] = (C[0]*C[0] - Dx_sq) + (C[1]*C[1] - Dy_sq);

            return (Aux::VecMat::det3x3<double>(M) > 0.0);
        }

        /* static aux function required below */
        static bool
        getTwoTrianglesSharedEdgeAndRemainingVertices(
            uint32_t    A_v0,
            uint32_t    A_v1,
            uint32_t    A_v2,
            uint32_t    B_v0,
            uint32_t    B_v1,
            uint32_t    B_v2,
            uint32_t   &fst_shared_vertex,
            uint32_t   &snd_shared_vertex,
            uint32_t   &A_remaining_vertex,
            uint32_t   &B_remaining_vertex)
        {
            /* first, check if triangles share exactly one common edge */
            std::vector<uint32_t>   A_vertices = { A_v0, A_v1, A_v2 };
            std::vector<uint32_t>   B_vertices = { B_v0, B_v1, B_v2 };
            std::vector<uint32_t>   AB_shared_vertices;
            std::vector<uint32_t>   A_remaining_vertices, B_remaining_vertices;

            std::sort(A_vertices.begin(), A_vertices.end());
            std::sort(B_vertices.begin(), B_vertices.end());
            std::set_intersection(
                    A_vertices.begin(),
                    A_vertices.end(),
                    B_vertices.begin(),
                    B_vertices.end(),
                    std::inserter(AB_shared_vertices, AB_shared_vertices.begin() )
                );

            /* if both triangles share exactly two vertices, i.e. a common edge */
            if (AB_shared_vertices.size() == 2) {
                /* get remaining vertices by set difference. the created lists should have size == 1 always
                 * */
                std::sort(AB_shared_vertices.begin(), AB_shared_vertices.end());
                std::set_difference(
                        A_vertices.begin(),
                        A_vertices.end(),
                        AB_shared_vertices.begin(),
                        AB_shared_vertices.end(),
                        std::inserter(A_remaining_vertices, A_remaining_vertices.begin())
                    );

                if (A_remaining_vertices.size() != 1) {
                    throw("(static) getTwoTrianglesSharedEdgeAndRemainingVertices(): two triangles share one edge, but A_remaining_vertices vector has size != 1. should be impossible..");
                }

                std::set_difference(
                        B_vertices.begin(),
                        B_vertices.end(),
                        AB_shared_vertices.begin(),
                        AB_shared_vertices.end(),
                        std::inserter(B_remaining_vertices, B_remaining_vertices.begin())
                    );

                if (B_remaining_vertices.size() != 1) {
                    throw("(static) getTwoTrianglesSharedEdgeAndRemainingVertices(): two triangles share one edge, but B_remaining_vertices vector has size != 1. should be impossible..");
                }

                /* write return variables */
                fst_shared_vertex   = AB_shared_vertices.front();
                snd_shared_vertex   = AB_shared_vertices.back();
                A_remaining_vertex  = A_remaining_vertices.front();
                B_remaining_vertex  = B_remaining_vertices.front();

                return true;
            }
            /* not edge-neighbours, return false */
            else return false;
        }

        /* static function required below */
        static bool
        getTriEdgeOrientationIndices(
            uint32_t v0,
            uint32_t v1,
            uint32_t v2,
            uint32_t u,
            uint32_t v)
        {
            /* check various cases */
            if (u == v0 && v == v1) {
                return true;
            }
            else if (u == v1 && v == v0) {
                return false;
            }
            else if (u == v1 && v == v2) {
                return true;
            }
            else if (u == v2 && v == v1) {
                return false;
            }
            else if (u == v2 && v == v0) {
                return true;
            }
            else if (u == v0 && v == v2) {
                return false;
            }
            else {
                throw ("getEdgeOrientationInTri(): given edge is not contained in supplied triangle.");
            }
        }

        bool
        delaunayTryFlip(
                Tri2d                          &A,
                Tri2d                          &B,
                std::map<uint32_t, Vertex2d>   &vertices)
        {
            debugl(4, "delaunayTryFlip(): triangles (%5d, %5d, %5d) - (%5d, %5d, %5d).\n", 
                    A.v0_id, A.v1_id, A.v2_id,
                    B.v0_id, B.v1_id, B.v2_id);
            debugTabInc();

            /* determine whether two triangles are neighbours. if so, check if the vertex of b that is not
             * in a (the remaining when the shared edge is taken out) lies in the circumcircle of a.
             * if so, flip edge and return true. in all other cases, return false */
            bool        neighbours;
            uint32_t    fst_shared_vertex, snd_shared_vertex;
            uint32_t    A_rem_vertex, B_rem_vertex;

            neighbours = getTwoTrianglesSharedEdgeAndRemainingVertices(
                    A.v0_id, A.v1_id, A.v2_id,
                    B.v0_id, B.v1_id, B.v2_id,
                    fst_shared_vertex, snd_shared_vertex,
                    A_rem_vertex, B_rem_vertex);

            /* if both triangles are "edge" neighbours, try to flip */
            if (neighbours) {
                debugl(4, "neighbours => checking local delaunay condition..\n");
                /* check if local delaunay condition is violated by numerically evaluating the determinant
                 * expression given in the paper cited in the thesis */
                if (vertexInCircumCircle(
                        vertices[A.v0_id].pos,
                        vertices[A.v1_id].pos,
                        vertices[A.v2_id].pos,
                        vertices[B_rem_vertex].pos))
                {
                    debugl(4, "local delaunay condition not satisifed => flipping\n");
                    /* get orientation of shared edge in A, this, together with the indices of the remaining
                     * vertices, uniquely defined the two "new" triangles and their orientation in case of a
                     * successful delaunay flip */
                    bool A_shared_edge_orientation = getTriEdgeOrientationIndices(
                            A.v0_id, A.v1_id, A.v2_id,
                            fst_shared_vertex, snd_shared_vertex);

                    /* if edge (fst_shared_vertex, snd_shared_vertex) is positively oriented in A, then the new
                     * triangles are: 
                     *
                     * (A_rem_vertex, fst_shared_vertex, B_remaining_vertex) and
                     * (A_rem_vertex, B_rem_vertex, snd_shared_vertex)
                     *
                     * and the roles of fst_shared_vertex and snd_shared_vertex get swapped if the orientation
                     * is negative in A. modify the two triangles in-place */
                    if (A_shared_edge_orientation) {
                        A.v0_id = A_rem_vertex;
                        A.v1_id = fst_shared_vertex;
                        A.v2_id = B_rem_vertex;

                        B.v0_id = A_rem_vertex;
                        B.v1_id = B_rem_vertex;
                        B.v2_id = snd_shared_vertex;
                    }
                    else {
                        A.v0_id = A_rem_vertex;
                        A.v1_id = snd_shared_vertex;
                        A.v2_id = B_rem_vertex;

                        B.v0_id = A_rem_vertex;
                        B.v1_id = B_rem_vertex;
                        B.v2_id = fst_shared_vertex;
                    }

                    /* flip has been applied, return true */
                    debugTabDec();
                    return true;
                } 
                /* two triangles are neighbours, but locally delaunay, don't flip */
                else {
                    debugl(4, "local delaunay condition satisifed. not flipping.\n");
                    debugTabDec();
                    return false;
                }
            }
            /* no neighbours, no flip is possible */
            else {
                debugl(4, "triangles are not neighbours => no flip possible.\n");
                debugTabDec();
                return false;
            }
        }

        /* naive algorithm: with an initial triangulation, flip edges of triangles violating the
         * local dalaunay condition until it converges or maximum number of iterations has been
         * reached (note that since the boundary is generally non-convex, convergence is NOT
         * guaranteed, but occurs quickly in most cases). */
        void
        delaunay2d(
            std::map<uint32_t, Vertex2d>   &vertices,
            std::vector<Tri2d>             &triangles)
        {
            bool                            again;
            const uint32_t                  max_iter = 2*vertices.size()*vertices.size();
            uint32_t                        iter;
            std::vector<Tri2d>::iterator    ta_it, tb_it;

            debugl(3, "delaunay2d(): flipping given triangulation until convergence..\n");
            debugTabInc();
            again   = true;
            iter    = 0;
            while (again && iter < max_iter)  {
                again = false;
                iter++;

                debugl(4, "no convergence yet => another flipping iteration..\n");
                debugTabInc();
                for (ta_it = triangles.begin(); ta_it != triangles.end(); ++ta_it) {
                    tb_it = ta_it;
                    ++tb_it;
                    for ( ; tb_it != triangles.end(); ++tb_it) {
                        if ( delaunayTryFlip(*ta_it, *tb_it, vertices) ) {
                            again = true;
                        }
                    }
                }
                debugTabDec();
            }
            if (iter == max_iter) {
                debugl(0, "delaunay2d() WARNING: returning due to maximum iteration limit, not due to convergence..\n");
                printf("delaunay2d() WARNING: returning due to maximum iteration limit, not due to convergence..\n");
            }
            debugTabDec();
            debugl(3, "delaunay2d(): done.\n");
        }


/* ---------------------------------------------------------------------------------------- */
/*                                                                                          */
/*                                                                                          */
/* FIXME : -------------- replace this with own clean implementaiton ---------- */
/*                                                                                          */
/*                                                                                          */
/* ---------------------------------------------------------------------------------------- */
#define IEPS 1E-10
        extern "C" {
            /* Triangle/triangle intersection test routine,
             * by Tomas Moller, 1997.
             * See article "A Fast Triangle-Triangle Intersection Test",
             * Journal of Graphics Tools, 2(2), 1997
             *
             * int tri_tri_intersect(double V0[3],double V1[3],double V2[3],
             *                         double U0[3],double U1[3],double U2[3])
             *
             * parameters: vertices of triangle 1: V0,V1,V2
             *             vertices of triangle 2: U0,U1,U2
             * result    : returns 1 if the triangles intersect, otherwise 0
             *
             */

#include <math.h>


            /* if USE_EPSILON_TEST is true then we do a check: 
                     if |dv|<EPSILON then dv=0.0;
               else no check is done (which is less robust)
            */
#define USE_EPSILON_TEST TRUE  
#define EPSILON IEPS


            /* some macros */
#define CROSS(dest,v1,v2)                      \
                          dest[0]=v1[1]*v2[2]-v1[2]*v2[1]; \
                          dest[1]=v1[2]*v2[0]-v1[0]*v2[2]; \
                          dest[2]=v1[0]*v2[1]-v1[1]*v2[0];

#define DOT(v1,v2) (v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2])

#define SUB(dest,v1,v2)          \
                        dest[0]=v1[0]-v2[0]; \
                        dest[1]=v1[1]-v2[1]; \
                        dest[2]=v1[2]-v2[2]; 

            /* sort so that a<=b */
#define SORT(a,b)       \
                         if(a>b)    \
                         {          \
                           double c; \
                           c=a;     \
                           a=b;     \
                           b=c;     \
                         }

#define ISECT(VV0,VV1,VV2,D0,D1,D2,isect0,isect1) \
                          isect0=VV0+(VV1-VV0)*D0/(D0-D1);    \
                          isect1=VV0+(VV2-VV0)*D0/(D0-D2);


#define COMPUTE_INTERVALS(VV0,VV1,VV2,D0,D1,D2,D0D1,D0D2,isect0,isect1) \
              if(D0D1>0.0f)                                         \
              {                                                     \
                /* here we know that D0D2<=0.0 */                   \
                /* that is D0, D1 are on the same side, D2 on the other or on the plane */ \
                ISECT(VV2,VV0,VV1,D2,D0,D1,isect0,isect1);          \
              }                                                     \
              else if(D0D2>0.0f)                                    \
              {                                                     \
                /* here we know that d0d1<=0.0 */                   \
                ISECT(VV1,VV0,VV2,D1,D0,D2,isect0,isect1);          \
              }                                                     \
              else if(D1*D2>0.0f || D0!=0.0f)                       \
              {                                                     \
                /* here we know that d0d1<=0.0 or that D0!=0.0 */   \
                ISECT(VV0,VV1,VV2,D0,D1,D2,isect0,isect1);          \
              }                                                     \
              else if(D1!=0.0f)                                     \
              {                                                     \
                ISECT(VV1,VV0,VV2,D1,D0,D2,isect0,isect1);          \
              }                                                     \
              else if(D2!=0.0f)                                     \
              {                                                     \
                ISECT(VV2,VV0,VV1,D2,D0,D1,isect0,isect1);          \
              }                                                     \
              else                                                  \
              {                                                     \
                /* triangles are coplanar */                        \
                return coplanar_tri_tri(N1,V0,V1,V2,U0,U1,U2);      \
              }



            /* this edge to edge test is based on Franlin Antonio's gem:
               "Faster Line Segment Intersection", in Graphics Gems III,
               pp. 199-202 */ 
#define EDGE_EDGE_TEST(V0,U0,U1)                      \
              Bx=U0[i0]-U1[i0];                                   \
              By=U0[i1]-U1[i1];                                   \
              Cx=V0[i0]-U0[i0];                                   \
              Cy=V0[i1]-U0[i1];                                   \
              f=Ay*Bx-Ax*By;                                      \
              d=By*Cx-Bx*Cy;                                      \
              if((f>0 && d>=0 && d<=f) || (f<0 && d<=0 && d>=f))  \
              {                                                   \
                e=Ax*Cy-Ay*Cx;                                    \
                if(f>0)                                           \
                {                                                 \
                  if(e>=0 && e<=f) return 1;                      \
                }                                                 \
                else                                              \
                {                                                 \
                  if(e<=0 && e>=f) return 1;                      \
                }                                                 \
              }                                

#define EDGE_AGAINST_TRI_EDGES(V0,V1,U0,U1,U2) \
            {                                              \
              double Ax,Ay,Bx,By,Cx,Cy,e,d,f;               \
              Ax=V1[i0]-V0[i0];                            \
              Ay=V1[i1]-V0[i1];                            \
              /* test edge U0,U1 against V0,V1 */          \
              EDGE_EDGE_TEST(V0,U0,U1);                    \
              /* test edge U1,U2 against V0,V1 */          \
              EDGE_EDGE_TEST(V0,U1,U2);                    \
              /* test edge U2,U1 against V0,V1 */          \
              EDGE_EDGE_TEST(V0,U2,U0);                    \
            }

#define POINT_IN_TRI(V0,U0,U1,U2)           \
            {                                           \
              double a,b,c,d0,d1,d2;                     \
              /* is T1 completly inside T2? */          \
              /* check if V0 is inside tri(U0,U1,U2) */ \
              a=U1[i1]-U0[i1];                          \
              b=-(U1[i0]-U0[i0]);                       \
              c=-a*U0[i0]-b*U0[i1];                     \
              d0=a*V0[i0]+b*V0[i1]+c;                   \
                                                        \
              a=U2[i1]-U1[i1];                          \
              b=-(U2[i0]-U1[i0]);                       \
              c=-a*U1[i0]-b*U1[i1];                     \
              d1=a*V0[i0]+b*V0[i1]+c;                   \
                                                        \
              a=U0[i1]-U2[i1];                          \
              b=-(U0[i0]-U2[i0]);                       \
              c=-a*U2[i0]-b*U2[i1];                     \
              d2=a*V0[i0]+b*V0[i1]+c;                   \
              if(d0*d1>0.0)                             \
              {                                         \
                if(d0*d2>0.0) return 1;                 \
              }                                         \
            }

            int coplanar_tri_tri(double N[3],double V0[3],double V1[3],double V2[3],
                                 double U0[3],double U1[3],double U2[3])
            {
               double A[3];
               short i0,i1;
               /* first project onto an axis-aligned plane, that maximizes the area */
               /* of the triangles, compute indices: i0,i1. */
               A[0]=fabs(N[0]);
               A[1]=fabs(N[1]);
               A[2]=fabs(N[2]);
               if(A[0]>A[1])
               {
                  if(A[0]>A[2])  
                  {
                      i0=1;      /* A[0] is greatest */
                      i1=2;
                  }
                  else
                  {
                      i0=0;      /* A[2] is greatest */
                      i1=1;
                  }
               }
               else   /* A[0]<=A[1] */
               {
                  if(A[2]>A[1])
                  {
                      i0=0;      /* A[2] is greatest */
                      i1=1;                                           
                  }
                  else
                  {
                      i0=0;      /* A[1] is greatest */
                      i1=2;
                  }
                }               
                            
                /* test all edges of triangle 1 against the edges of triangle 2 */
                EDGE_AGAINST_TRI_EDGES(V0,V1,U0,U1,U2);
                EDGE_AGAINST_TRI_EDGES(V1,V2,U0,U1,U2);
                EDGE_AGAINST_TRI_EDGES(V2,V0,U0,U1,U2);
                            
                /* finally, test if tri1 is totally contained in tri2 or vice versa */
                POINT_IN_TRI(V0,U0,U1,U2);
                POINT_IN_TRI(U0,V0,V1,V2);

                return 0;
            }


            int tri_tri_intersect(double V0[3],double V1[3],double V2[3],
                                  double U0[3],double U1[3],double U2[3])
            {
              double E1[3],E2[3];
              double N1[3],N2[3],d1,d2;
              double du0,du1,du2,dv0,dv1,dv2;
              double D[3];
              double isect1[2], isect2[2];
              double du0du1,du0du2,dv0dv1,dv0dv2;
              short index;
              double vp0,vp1,vp2;
              double up0,up1,up2;
              double b,c,max;

              /* compute plane equation of triangle(V0,V1,V2) */
              SUB(E1,V1,V0);
              SUB(E2,V2,V0);
              CROSS(N1,E1,E2);
              d1=-DOT(N1,V0);
              /* plane equation 1: N1.X+d1=0 */

              /* put U0,U1,U2 into plane equation 1 to compute signed distances to the plane*/
              du0=DOT(N1,U0)+d1;
              du1=DOT(N1,U1)+d1;
              du2=DOT(N1,U2)+d1;

              /* coplanarity robustness check */
#if USE_EPSILON_TEST==TRUE
              if(fabs(du0)<EPSILON) du0=0.0;
              if(fabs(du1)<EPSILON) du1=0.0;
              if(fabs(du2)<EPSILON) du2=0.0;
#endif
              du0du1=du0*du1;
              du0du2=du0*du2;

              if(du0du1>0.0f && du0du2>0.0f) /* same sign on all of them + not equal 0 ? */
                return 0;                    /* no intersection occurs */

              /* compute plane of triangle (U0,U1,U2) */
              SUB(E1,U1,U0);
              SUB(E2,U2,U0);
              CROSS(N2,E1,E2);
              d2=-DOT(N2,U0);
              /* plane equation 2: N2.X+d2=0 */

              /* put V0,V1,V2 into plane equation 2 */
              dv0=DOT(N2,V0)+d2;
              dv1=DOT(N2,V1)+d2;
              dv2=DOT(N2,V2)+d2;

#if USE_EPSILON_TEST==TRUE
              if(fabs(dv0)<EPSILON) dv0=0.0;
              if(fabs(dv1)<EPSILON) dv1=0.0;
              if(fabs(dv2)<EPSILON) dv2=0.0;
#endif

              dv0dv1=dv0*dv1;
              dv0dv2=dv0*dv2;
                    
              if(dv0dv1>0.0f && dv0dv2>0.0f) /* same sign on all of them + not equal 0 ? */
                return 0;                    /* no intersection occurs */

              /* compute direction of intersection line */
              CROSS(D,N1,N2);

              /* compute and index to the largest component of D */
              max=fabs(D[0]);
              index=0;
              b=fabs(D[1]);
              c=fabs(D[2]);
              if(b>max) max=b,index=1;
              if(c>max) index=2;

                    /* this is the simplified projection onto L*/
                    vp0=V0[index];
                    vp1=V1[index];
                    vp2=V2[index];

                    up0=U0[index];
                    up1=U1[index];
                    up2=U2[index];

              /* compute interval for triangle 1 */
              COMPUTE_INTERVALS(vp0,vp1,vp2,dv0,dv1,dv2,dv0dv1,dv0dv2,isect1[0],isect1[1]);

              /* compute interval for triangle 2 */
              COMPUTE_INTERVALS(up0,up1,up2,du0,du1,du2,du0du1,du0du2,isect2[0],isect2[1]);

              SORT(isect1[0],isect1[1]);
              SORT(isect2[0],isect2[1]);

              if(isect1[1]<isect2[0] || isect2[1]<isect1[0]) return 0;
              return 1;
            }
        }
/* ---------------------------------------------------------------------------------------- */

    }
}

/* --------- debug functionality ---------- */
static uint32_t debug_component                 = DBG_GLOBAL;
static bool     debug_component_enabled[128]    = { 0 };
static uint32_t debug_level                     = 0;
static uint32_t debug_max_level                 = 0;
static uint32_t debug_tab                       = 0;

void 
initDebug()
{
    int i;
    for (i = 0; i < 128; i++) {
        debug_component_enabled[i] = false;
    }
    debug_component                     = DBG_GLOBAL;
    debug_component_enabled[DBG_GLOBAL] = true;
    debug_level                         = 0;
    debug_max_level                     = 0;
}

void
setDebugComponent(uint32_t comp)
{
    debug_component = comp;
}

void
setDebugTab(uint32_t tab)
{
    debug_tab = tab;
}

uint32_t
getDebugTab()
{
    return debug_tab;
}



void
setDebugLevel(uint32_t level)
{
    debug_level = level;
}

void
setMaxDebugLevel(uint32_t max_level)
{
    debug_max_level = max_level;
}

void
enableComponentDebug(uint32_t comp)
{
    if (comp < 128) {
        debug_component_enabled[comp] = true;
    }
}

void
disableComponentDebug(uint32_t comp)
{
    if (comp < 128) {
        debug_component_enabled[comp] = false;
    }
}

void
debugTabIncrement()
{
    debug_tab++;
}   

void
debugTabDecrement()
{
    if (debug_tab == 0) {
        printf("debugTabDecrement(): already 0. mismatched tab calls.");
    }
    else {
        debug_tab--;
    }
}

int
debugprintf(std::string filename, int line, std::string fmt, ...)
{
    if (debug_component_enabled[debug_component]) {
        uint32_t    i;
        int         ret = 0;
        va_list     ap;

        char tabs[512]  = {0};
        const char* singletab  = "  ";
        /* construct tabs */
        for (i = 0; i < debug_tab; i++) {
            strncat(tabs, singletab, 2);
        }

        char foo[1024] = {0};
        strncpy(foo, filename.c_str(), 1024);

        char *bar = strrchr(foo, '/');
        if (bar) {
            /* get rid of last / */
            ++bar;

            printf("[ %20s | %5d ]: %s", bar, line, tabs);
        } 
        else {
            printf("[ %20s | %5d ]: %s", foo, line, tabs);
        }

        va_start(ap, fmt);
        ret = vprintf(fmt.c_str(), ap);
        va_end(ap);

        return ret;
    }
    else return 0;
}

int
debugprintf_wlevel(uint32_t level, std::string filename, int line, std::string fmt, ...)
{
    /* regular debug code. we can't call debugprintf() due to va_args macro hacking .. */
    if (debug_component_enabled[debug_component] && level <= debug_max_level) {
        uint32_t    i;
        va_list     ap;

        char foo[1024] = {0};
        strncpy(foo, filename.c_str(), 1024);
        char *bar = strrchr(foo, '/');

        char buffer[1024];
        int pos = 0;

        if (bar) {
            /* get rid of last / */
            ++bar;

            pos += sprintf(buffer, "[ %20s | %5d ]: ", bar, line);
            for (i = 0; i < debug_tab; i++) {
                pos += sprintf(&buffer[pos], "  ");
            }
        } 
        else {
            pos += sprintf(buffer, "[ %20s | %5d ]: ", foo, line);
            for (i = 0; i < debug_tab; i++) {
                pos += sprintf(&buffer[pos], "  ");
            }
        }

        va_start(ap, fmt);
        pos += vsprintf(&buffer[pos], fmt.c_str(), ap);
        printf("%s", buffer);
        va_end(ap);

        return pos;
    }
    else {
        return 0;
    }
}
