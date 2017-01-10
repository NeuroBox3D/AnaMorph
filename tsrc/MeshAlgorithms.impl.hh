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

#ifndef MESHALG_TCC
#define MESHALG_TCC

#include "common.hh"
#include "Octree.hh"
#include "Mesh.hh"
#include "MeshAlgorithms.hh"

#include "uTuple.hh"
#include "PriorityQueue.hh"


template <typename Tm, typename Tv, typename Tf, typename R>
void
refineTriangularUnitSphereApproximation(Mesh<Tm, Tv, Tf, R> &S)
{
    std::map<
            uPair<uint32_t>,
            typename Mesh<Tm, Tv, Tf, R>::vertex_iterator
        >                                                   id_cache;

    typename std::map<
            uPair<uint32_t>,
            typename Mesh<Tm, Tv, Tf, R>::vertex_iterator
        >::iterator                                         cachit;

    typedef std::tuple<
                typename Mesh<Tm, Tv, Tf, R>::vertex_iterator,
                typename Mesh<Tm, Tv, Tf, R>::vertex_iterator,
                typename Mesh<Tm, Tv, Tf, R>::vertex_iterator
            >                                               tri_tuple;
    std::vector<tri_tuple>                                  new_faces;

    /* iterate over all triangles of s, add new vertex for all three edges while avoiding double
     * vertices with id_cache. append four triples identifying the four new faces replacing the
     * currently considered old triangle in the produced refined mesh. */
    Vec3<R>                                         vnew_pos;
    uint32_t                                        v0_id, v1_id, v2_id;
    typename Mesh<Tm, Tv, Tf, R>::vertex_iterator   v0, v1, v2, v01, v12, v20;
    for (auto &f : S.faces) {
        /* get vertex indices and positions */
        f.getTriIndices(v0_id, v1_id, v2_id);
        f.getTriIterators(v0, v1, v2);

        /* add new vertices for (v0, v1), (v1, v2) and (v2, v0) */

        /* cache hit: get id */
        if ( (cachit = id_cache.find(uPair<uint32_t>(v0_id, v1_id))) != id_cache.end()) {
            v01 = cachit->second;
        }
        /* cache miss: add new vertex, insert into cache */
        else {
            vnew_pos            = (v0->pos() + v1->pos()) / 2.0;
            vnew_pos.normalize();

            v01                 = S.vertices.insert(vnew_pos);
            auto insres         = id_cache.insert( {uPair<uint32_t>(v0_id, v1_id), v01 } );
            if (!insres.second) {
                throw("refineTriangularUnitSphereApproximation(): after cache miss, id of newly created vertex on an edge has already been present in cache.");
            }
        }

        /* same for other two edges */
        if ( (cachit = id_cache.find(uPair<uint32_t>(v1_id, v2_id))) != id_cache.end()) {
            v12 = cachit->second;
        }
        else {
            vnew_pos            = (v1->pos() + v2->pos()) / 2.0;
            vnew_pos.normalize();

            v12                 = S.vertices.insert(vnew_pos);
            auto insres         = id_cache.insert( { uPair<uint32_t>(v1_id, v2_id), v12 } );
            if (!insres.second) {
                throw("refineTriangularUnitSphereApproximation(): after cache miss, id of newly created vertex on an edge has already been present in cache.");
            }
        }

        if ( (cachit = id_cache.find(uPair<uint32_t>(v2_id, v0_id))) != id_cache.end()) {
            v20 = cachit->second;
        }
        else {
            vnew_pos            = (v2->pos() + v0->pos()) / 2.0;
            vnew_pos.normalize();

            v20                 = S.vertices.insert(vnew_pos);
            auto insres         = id_cache.insert( { uPair<uint32_t>(v2_id, v0_id), v20 } );
            if (!insres.second) {
                throw("refineTriangularUnitSphereApproximation(): after cache miss, id of newly created vertex on an edge has already been present in cache.");
            }
        }

        /* add four triangles (v0, v01, v20), (v20, v01, v12), (v12, v01, v1), (v12, v2, v20) to
         * new_faces */
        new_faces.push_back( tri_tuple(v0,      v01,    v20 ) );
        new_faces.push_back( tri_tuple(v20,     v01,    v12 ) );
        new_faces.push_back( tri_tuple(v12,     v01,    v1  ) );
        new_faces.push_back( tri_tuple(v12,     v2,     v20 ) );
    }

    /* clear all faces, add all new triangles of the refined mesh */
    S.clearFaces();

    for (auto &nf : new_faces) {
        S.faces.insert( std::get<0>(nf), std::get<1>(nf), std::get<2>(nf) );
    }
}

template <typename Tm, typename Tv, typename Tf, typename R>
void
MeshAlg::generateOctSphere(
    Vec3<R>                 c,
    R const                &r,
    uint32_t                tessellation_depth,
    Mesh<Tm, Tv, Tf, R>    &S)
{
    /* clear input mesh */
    S.clear();

    R const                                         one_over_sqrt2 = (R)1.0 / std::sqrt((R)2.0);
    typename Mesh<Tm, Tv, Tf, R>::vertex_iterator   sqrt_00, sqrt_01, sqrt_10, sqrt_11, top, bottom;

    /* add 6 octahedron vertices and save their ids */
    sqrt_00 = S.vertices.insert( Vec3<R>(-one_over_sqrt2,  -one_over_sqrt2, 0.0) ); 
    sqrt_01 = S.vertices.insert( Vec3<R>(-one_over_sqrt2,   one_over_sqrt2, 0.0) ); 
    sqrt_10 = S.vertices.insert( Vec3<R>( one_over_sqrt2,  -one_over_sqrt2, 0.0) ); 
    sqrt_11 = S.vertices.insert( Vec3<R>( one_over_sqrt2,   one_over_sqrt2, 0.0) ); 
    top     = S.vertices.insert( Vec3<R>( 0.0,              0.0,            1.0) );
    bottom  = S.vertices.insert( Vec3<R>( 0.0,              0.0,           -1.0) );

    /* add 8 octahedron faces with correct orientation */
    S.faces.insert(sqrt_11,   top,        sqrt_10);
    S.faces.insert(sqrt_11,   sqrt_10,    bottom);
    S.faces.insert(sqrt_01,   sqrt_00,    top);
    S.faces.insert(sqrt_01,   bottom,     sqrt_00);
    S.faces.insert(sqrt_10,   top,        sqrt_00);
    S.faces.insert(sqrt_10,   sqrt_00,    bottom);
    S.faces.insert(sqrt_11,   sqrt_01,    top);
    S.faces.insert(sqrt_11,   bottom,     sqrt_01);

    /* refine tessellation_depth times */
    for (uint32_t i = 0; i < tessellation_depth; i++) {
        refineTriangularUnitSphereApproximation(S);
    }

    /* scale all vertices with radius r and translate origin to centre c */
    S.scale(r);
    S.translate(c);
}

/* ico sphere */
template <typename Tm, typename Tv, typename Tf, typename R>
void
MeshAlg::generateIcoSphere(
    Vec3<R>                 c,
    R const                &r,
    uint32_t                tessellation_depth,
    Mesh<Tm, Tv, Tf, R>    &S)
{
    /* clear input mesh */
    S.clear();

    R const                                         phi     = (1.0 + std::sqrt((R)5.0)) / 2.0;
    R const                                         norm    = std::sqrt((R)2.0 + phi);
    typename Mesh<Tm, Tv, Tf, R>::vertex_iterator   vits[12];

    /* add 12 icosahedron vertices and save their ids */
    vits[ 0]    = S.vertices.insert( Vec3<R>(  0.0, -1.0, -phi ) / norm );
    vits[ 1]    = S.vertices.insert( Vec3<R>(  0.0, -1.0, +phi ) / norm );
    vits[ 2]    = S.vertices.insert( Vec3<R>(  0.0, +1.0, -phi ) / norm );
    vits[ 3]    = S.vertices.insert( Vec3<R>(  0.0, +1.0, +phi ) / norm );

    vits[ 4]    = S.vertices.insert( Vec3<R>( -phi,  0.0, -1.0 ) / norm );
    vits[ 5]    = S.vertices.insert( Vec3<R>( -phi,  0.0, +1.0 ) / norm );
    vits[ 6]    = S.vertices.insert( Vec3<R>( +phi,  0.0, -1.0 ) / norm );
    vits[ 7]    = S.vertices.insert( Vec3<R>( +phi,  0.0, +1.0 ) / norm );

    vits[ 8]    = S.vertices.insert( Vec3<R>( -1.0, -phi,  0.0 ) / norm );
    vits[ 9]    = S.vertices.insert( Vec3<R>( -1.0, +phi,  0.0 ) / norm );
    vits[10]    = S.vertices.insert( Vec3<R>( +1.0, -phi,  0.0 ) / norm );
    vits[11]    = S.vertices.insert( Vec3<R>( +1.0, +phi,  0.0 ) / norm );

    /* add 20 icosahedron faces with correct orientation */
    S.faces.insert( vits[ 0], vits[ 8], vits[ 4]);
    S.faces.insert( vits[ 0], vits[ 2], vits[ 6]);
    S.faces.insert( vits[ 1], vits[ 7], vits[ 3]);
    S.faces.insert( vits[ 1], vits[ 5], vits[ 8]);
    S.faces.insert( vits[ 2], vits[ 0], vits[ 4]);
    S.faces.insert( vits[ 2], vits[11], vits[ 6]);
    S.faces.insert( vits[ 3], vits[ 7], vits[11]);
    S.faces.insert( vits[ 3], vits[ 5], vits[ 1]);
    S.faces.insert( vits[ 4], vits[ 8], vits[ 5]);
    S.faces.insert( vits[ 6], vits[11], vits[ 7]);
    S.faces.insert( vits[ 9], vits[ 5], vits[ 3]);
    S.faces.insert( vits[ 9], vits[ 3], vits[11]);
    S.faces.insert( vits[ 9], vits[11], vits[ 2]);
    S.faces.insert( vits[ 9], vits[ 2], vits[ 4]);
    S.faces.insert( vits[ 9], vits[ 4], vits[ 5]);
    S.faces.insert( vits[10], vits[ 7], vits[ 1]);
    S.faces.insert( vits[10], vits[ 1], vits[ 8]);
    S.faces.insert( vits[10], vits[ 8], vits[ 0]);
    S.faces.insert( vits[10], vits[ 0], vits[ 6]);
    S.faces.insert( vits[10], vits[ 6], vits[ 7]);

    /* refine tessellation_depth times */
    for (uint32_t i = 0; i < tessellation_depth; i++) {
        refineTriangularUnitSphereApproximation(S);
    }

    /* scale all vertices with radius r and translate origin to centre c */
    S.scale(r);
    S.translate(c);
}

template <typename Tm, typename Tv, typename Tf, typename R>
void
MeshAlg::appendHalfSphereToCanalSurfaceMesh(
    Mesh<Tm, Tv, Tf, R>                            &M,
    Vec3<R>                                         render_vector,
    Vec3<R>                                         start,
    R const                                        &radius,
    Vec3<R>                                         direction,
    uint32_t                                        nphisegments,
    R const                                        &phi_offset,
    std::vector<
            typename Mesh<Tm, Tv, Tf, R>
                ::vertex_iterator
        >                                           start_circle_its,
    typename Mesh<Tm, Tv, Tf, R>::vertex_iterator   closing_vertex_it)
{
    using Common::twopi;

    uint32_t    i, j;
    uint32_t    ntsegments;
    Vec3<R>     p, px, py, pz;
    R           t, phi, dt, dphi, r, alpha, dalpha;

    /* vectors storing vertex iterators of current and last circle */
    std::vector<typename Mesh<Tm, Tv, Tf, R>::vertex_iterator>  last_circle;
    std::vector<typename Mesh<Tm, Tv, Tf, R>::vertex_iterator>  current_circle;

    typename Mesh<Tm, Tv, Tf, R>::vertex_iterator               end_closing_vertex_it;

    /* experimentally conceived value for ntsegments */
    ntsegments  = nphisegments / 2;
    dt          = 1.0 / (R)ntsegments;
    dphi        = twopi / (R)nphisegments;
    dalpha      = (M_PI / 2.0) / (R)ntsegments;

    /* resize vectors, assign current_circle to start_circle_ids */
    current_circle.resize(nphisegments);
    last_circle.resize(nphisegments);
    current_circle = start_circle_its;

    /* delete the closing vertex, leaving the end "circle" open. render circles of the sphere and
     * connect them as for the canal surfaces. radius function: sqrt(1 - t^2) for t in [0,1].
     * stop before t = 1, and close the thing half-sphere up with a triangle fan. */
    M.vertices.erase(closing_vertex_it);

    /* normalize direction vector. we get the current centre of the circle as 
     *      p = start + t*r*direction
     * */
    direction.normalize();

    /* the "start" circle already exists and needs not be created again. the "end" circle of the
     * half-sphere coincides with a single point, which is specially connected to the second to last
     * circle, i.e. the last "real" non-degenerate circle, wiht a triangle fan */
    for (i = 1; i < ntsegments - 1; i++) {
        /* shift circle ids */
        last_circle = current_circle;

        /* get orthonormal base. direction is the x vector */
        px = direction;
        px.normalize();

        py = render_vector.cross(px);
        py.normalize();

        pz = px.cross(py);

        /* t = cos(alpha), where alpha = pi/2 - i*dalpha */
        alpha   = (M_PI / 2.0) - i*dalpha;
        t       = std::cos(alpha);

        /* calculate radius: r = radius * sqrt(1 - t^2) */
        r       = radius * std::sqrt(1.0 - t*t);

        /* the point p: start + (t*r)*direction */
        p   = start + direction*(t*radius);

        for (j = 0; j < nphisegments; j++) {
            phi                 = ( (R)j * twopi) / (R)nphisegments;
            current_circle[j]   = M.vertices.insert( p + py*(r*std::cos(phi + phi_offset)) + pz*(r*std::sin(phi + phi_offset)) );
        }

        /* generate quad faces between last_circle and current_circle */
        for (j = 0; j < nphisegments - 1; j++) {
            M.faces.insert(last_circle[j], last_circle[j+1], current_circle[j+1], current_circle[j]);
        }

        /* closing quad at index warp around */
        M.faces.insert(last_circle[nphisegments - 1], last_circle[0], current_circle[0], current_circle[nphisegments - 1]);
    }

    /* final point and closing triangle fan */
    p                       = start + direction*radius;
    end_closing_vertex_it   = M.vertices.insert(p);

    M.faces.insert(end_closing_vertex_it, current_circle[nphisegments - 1], current_circle[0]);
    for (j = 0; j < nphisegments - 1; j++) {
        M.faces.insert(end_closing_vertex_it, current_circle[j], current_circle[j + 1] );
    }
}


/* ---------------------------------------------------------------------------------------------- */
/*                                                                                                */
/*    Computation of potentially intersecting edge / face pairs with implicit Octree traversal    */
/*                                                                                                */
/* ---------------------------------------------------------------------------------------------- */
template <typename Tm, typename Tv, typename Tf, typename R>
void 
MeshAlg::getPotentiallyIntersectingEdgeFacePairs(
    Mesh<Tm, Tv, Tf, R>                        &X,
    Mesh<Tm, Tv, Tf, R>                        &Y,
    std::vector<EdgeFacePair<Mesh<Tm, Tv, Tf, R> > >&  X_edges_Y_faces_candidates,
    std::vector<EdgeFacePair<Mesh<Tm, Tv, Tf, R> > >&  Y_edges_X_faces_candidates,
    uint32_t                                    max_components,
    uint32_t                                    max_recursion_depth)
{
    typedef EdgeFacePair<Mesh<Tm, Tv, Tf, R> > EFPtype;
    typedef typename Mesh<Tm, Tv, Tf, R>::Face FaceType;
    using namespace Aux::Timing;
    tick(12);

    debugl(1, "MeshAlg::getPotentialEdgeFacePairs(): max_components: %5d, max_recursion_depth: %5d.\n", max_components, max_recursion_depth);
    debugTabInc();

    /* get bounding box surrounding both X and Y */
    auto XY_bb  = (X.getBoundingBox()).update(Y.getBoundingBox());

    /* store lists of all faces of X and Y in root lists for recursive top-down Octree-like algorithm 
     * Common::Geometry::computeSpatialIntersectionCandidatePairs() */
    std::vector<std::pair<FaceType*, BoundingBox<R> > > X_face_info_list, Y_face_info_list;
    X_face_info_list.reserve(X.numFaces());
    Y_face_info_list.reserve(Y.numFaces());
    for (auto &f : X.faces)
        X_face_info_list.push_back({&f, f.getBoundingBox()});
    for (auto &f : Y.faces)
        Y_face_info_list.push_back({&f, f.getBoundingBox()});

    std::vector<std::pair<FaceType*, FaceType*> > candidate_pairs;

    /* call recursive top-down spatial intersection algorithm */
    debugl(0, "MeshAlg::getPotentialEdgeFacePairs(): calling Aux::Geometry::computeSpatialIntersectionCandidatePairs().\n");
    Aux::Timing::tick(16);

    Aux::Geometry::computeSpatialIntersectionCandidatePairs<FaceType*, FaceType*, R>
    (
        XY_bb,
        X_face_info_list, Y_face_info_list,
        //face_bb_getter, face_bb_getter,
        //NULL, NULL,
        0, max_components, max_recursion_depth,
        candidate_pairs
    );

    debugl(0, "MeshAlg::getPotentialEdgeFacePairs(): done. time: %5.4f\n", Aux::Timing::tack(16));

    /* sort() and unique() candidate pairs. NOTE: comparison operators "<" and ">" are only well-defined for pointers if
     * certain criteria are met (see § 5.9 of the C++11 standard for further details), the cleanest and safest way is to
     * sort based on ids(). however: this requires the pointers to be "valid" in the sense that they pointer to properly
     * allocated data, which definitely holds for this particular situation. in general, this is unsatisfactory and
     * std::less<T *> for pointers must be used, which is guaranteed to provide a total ordering among pointers of the
     * same type (including NULL) by the standard. */
    auto face_ptr_pair_cmp = 
        [] (
            const std::pair<typename Mesh<Tm, Tv, Tf, R>::Face *, typename Mesh<Tm, Tv, Tf, R>::Face *> &x, 
            const std::pair<typename Mesh<Tm, Tv, Tf, R>::Face *, typename Mesh<Tm, Tv, Tf, R>::Face *> &y)
        -> bool 
        {
            if (x.first->id() < y.first->id()) return true;
            if (y.first->id() < x.first->id()) return false;
            return x.second->id() < y.second->id();
        };

    std::sort(candidate_pairs.begin(), candidate_pairs.end(), face_ptr_pair_cmp);
    auto newPairEnd = std::unique(candidate_pairs.begin(), candidate_pairs.end());
    candidate_pairs.erase(newPairEnd, candidate_pairs.end());

    debugl(1, "returned candidate pairs: %d\n", candidate_pairs.size());

    /* compute result from all pairs of potentially intersecting faces from X and Y */
    typename Mesh<Tm, Tv, Tf, R>::Face     *X_face, *Y_face;
    typename Mesh<Tm, Tv, Tf, R>::Vertex   *e_u, *e_v;
    std::list<
            std::pair<
                typename Mesh<Tm, Tv, Tf, R>::vertex_iterator,
                typename Mesh<Tm, Tv, Tf, R>::vertex_iterator
            >
        >                                   face_edges;

    for (auto &pair : candidate_pairs) {
        X_face  = pair.first;
        Y_face  = pair.second;

        /* get all edges e = {u, v} of X_face, insert pair (e, Y_face) into X_edges_Y_faces_candidates,
         * where care is taken to choose u as the "smaller" vertex pointer, using
         * Mesh<Tm, Tv, Tf>::Vertex::ptr_less for comparison. */
        X_face->getEdges(face_edges);
        for (auto &e : face_edges) {
            e_u = &(*e.first);
            e_v = &(*e.second);

            if (Mesh<Tm, Tv, Tf, R>::Vertex::ptr_less( e_v, e_u )) {
                std::swap(e_u, e_v);
            }

            X_edges_Y_faces_candidates.push_back(EFPtype(e_u, e_v, Y_face));
        }

        /* and vice versa for all edges e of Y_face and the corresponding pairs (e, X_face) */
        Y_face->getEdges(face_edges);
        for (auto &e : face_edges) {
            e_u = &(*e.first);
            e_v = &(*e.second);

            /* swap if necessary */
            if (Mesh<Tm, Tv, Tf, R>::Vertex::ptr_less( e_v, e_u )) {
                std::swap(e_u, e_v);
            }

            //Y_edges_X_faces_candidates.push_back( { { &(*e.first), &(*e.second) }, X_face } );
            Y_edges_X_faces_candidates.push_back(EFPtype(e_u, e_v, X_face));
        }
    }

    std::sort(X_edges_Y_faces_candidates.begin(), X_edges_Y_faces_candidates.end());
    auto newEnd = std::unique(X_edges_Y_faces_candidates.begin(), X_edges_Y_faces_candidates.end());
    X_edges_Y_faces_candidates.erase(newEnd, X_edges_Y_faces_candidates.end());

    std::sort(Y_edges_X_faces_candidates.begin(), Y_edges_X_faces_candidates.end());
    newEnd = std::unique(Y_edges_X_faces_candidates.begin(), Y_edges_X_faces_candidates.end());
    Y_edges_X_faces_candidates.erase(newEnd, Y_edges_X_faces_candidates.end());

    debugTabDec();

    debugl(1, "Mesh::getPotentialEdgeFacePairs(): time: %5.4f\n", tack(12));
}


/* ---------------------------------------------------------------------------------------------- */
/*                                                                                                */
/*                                Red-Blue-Union algorithm                                        */
/*                                                                                                */
/* ---------------------------------------------------------------------------------------------- */

/* red:                         flag. true: red point: point lies on an edge of M and properly inside a face of
 *                              A. false: blue point: lies on an edge of A and properly inside a face of
 *                              M. these conditions are enforced by the algorithm. edge cases, such
 *                              as as points that lie on edges of both M and A, are avoided by
 *                              unnoticeably tiny random perturbation of the vertex positions.
 *
 * x:                           position vector of the edge-face intersection point.
 *
 * edge_x_lambda:               
 * face_x_s, face_x_t:          parametric coordinates of x: lambda values for edge and barycentric
 *                              coordinates in face (again, everything depending on type).
 * u_it, v_it:              
 * face_it:                     edge (u_it, v_it) and id of face the point lies on. if type == true -> face_id
 *                              referring to A, edge_id referring to M. vice versa if type == false.
 *
 * backward_tri_it,
 * forward_tri_it:              ids of the two incident triangles in the mesh where point p lies on
 *                              an edge. again, if type == true => point lies on an edge of M => 
 *                              two ids of triangles incident ot that edge in M. symmetric for
 *                              false.
 *                              orientation is backward_tri -> common_edge -> forward_tr/
 *
 * R_new_it, B_new_it:          the corners of the intersection polygon will be added as new vertices in both meshes.
 *                              these two variables store the corresponding vertex iterators.
 * */
template <typename Tm, typename Tv, typename Tf, typename TR>
struct RB_Tuple {
    bool                                            red;
    Vec3<TR>                                        x;
    TR                                              edge_x_lambda, face_x_s, face_x_t;

    typename Mesh<Tm, Tv, Tf, TR>::vertex_iterator  edge_out_it, edge_in_it;
    typename Mesh<Tm, Tv, Tf, TR>::face_iterator    face_it;
    typename Mesh<Tm, Tv, Tf, TR>::face_iterator    backward_tri_it, forward_tri_it;
    typename Mesh<Tm, Tv, Tf, TR>::vertex_iterator  R_new_it, B_new_it;

    RB_Tuple(
        bool                                            red,
        Vec3<TR> const                                 &x,
        TR const                                       &edge_x_lambda,
        TR const                                       &face_x_s,
        TR const                                       &face_x_t,
        typename Mesh<Tm, Tv, Tf, TR>::vertex_iterator  edge_out_it,
        typename Mesh<Tm, Tv, Tf, TR>::vertex_iterator  edge_in_it,
        typename Mesh<Tm, Tv, Tf, TR>::face_iterator    face_it,
        typename Mesh<Tm, Tv, Tf, TR>::face_iterator    backward_tri_it,
        typename Mesh<Tm, Tv, Tf, TR>::face_iterator    forward_tri_it)
    {
        this->red                   = red;
        this->x                     = x;
        this->edge_x_lambda         = edge_x_lambda;
        this->face_x_s              = face_x_s;
        this->face_x_t              = face_x_t;
        this->edge_out_it           = edge_out_it,
        this->edge_in_it            = edge_in_it;
        this->face_it               = face_it;
        this->backward_tri_it       = backward_tri_it;
        this->forward_tri_it        = forward_tri_it;
    }


    /* use lexicographic ordering provided by std::tuple by creating temporary tuples from the
     * members and iterator->id()s and comparing those */
    bool
    operator<(const RB_Tuple &x) const
    {
        return (
                std::tuple<bool, uint32_t, uint32_t, uint32_t>(this->red, this->edge_out_it->id(), this->edge_in_it->id(), this->face_it->id()) < 
                std::tuple<bool, uint32_t, uint32_t, uint32_t>(x.red, x.edge_out_it->id(), x.edge_in_it->id(), x.face_it->id())  
           );
    }

    /* same reduction of equality comparison to std::tuple equality */
    bool
    operator==(const RB_Tuple &x) const
    {
        return (
                std::tuple<bool, uint32_t, uint32_t, uint32_t>(this->red, this->edge_out_it->id(), this->edge_in_it->id(), this->face_it->id()) == 
                std::tuple<bool, uint32_t, uint32_t, uint32_t>(x.red, x.edge_out_it->id(), x.edge_in_it->id(), x.face_it->id())  
           );
    }

    bool
    isNeighbour(const RB_Tuple &b) const {
        return (
            (backward_tri_it    == b.backward_tri_it)   ||
            (backward_tri_it    == b.forward_tri_it)    ||
            (forward_tri_it     == b.backward_tri_it)   ||
            (forward_tri_it     == b.forward_tri_it)
        );
    }
};

/* forward declarations of inline auxiliary functions used in Red-Blue-Union algorithm */
template <typename Tm, typename Tv, typename Tf, typename TR>
inline void
RedBlue_generateRBTupleList(
    Mesh<Tm, Tv, Tf, TR>                                   &X,
    Mesh<Tm, Tv, Tf, TR>                                   &Y,
    bool                                                    X_red,
    std::vector<MeshAlg::EdgeFacePair<Mesh<Tm, Tv, Tf, TR> > > &X_edges_Y_faces_candidates,
    std::list<RB_Tuple<Tm, Tv, Tf, TR>>                    &X_colour_tuples,
    std::list<RedBlue_EdgeIsecInfo<TR>>                    &complex_edge_info_list);

template <typename Tm, typename Tv, typename Tf, typename TR>
inline void
RedBlue_cyclicallyOrderTupleList(
    const Mesh<Tm, Tv, Tf, TR>             &M,
    bool                                    M_red,
    std::list<RB_Tuple<Tm, Tv, Tf, TR>>    &M_tuples);

template <typename Tm, typename Tv, typename Tf, typename TR>
inline void
RedBlue_cutHole(
    Mesh<Tm, Tv, Tf, TR>                                       &M,
    bool                                                        M_red,
    std::list<RB_Tuple<Tm, Tv, Tf, TR>>                        &isecpoly_tuples,
    bool                                                        keep_outside_part   = true,
    std::vector<
            typename Mesh<Tm, Tv, Tf, TR>::vertex_iterator
        >                                                      *update_its          = NULL,
    std::list<
            typename Mesh<Tm, Tv, Tf, TR>::face_iterator
        >                                                      *new_tri_its         = NULL);


/* Red-Blue algorithm implementation */
template <typename Tm, typename Tv, typename Tf, typename TR>
void 
MeshAlg::RedBlueAlgorithm(
    Mesh<Tm, Tv, Tf, TR>                                   &R,
    Mesh<Tm, Tv, Tf, TR>                                   &B,
    const bool                                             &keep_red_outside_part,
    const bool                                             &keep_blue_outside_part,
    std::vector<
            typename Mesh<Tm, Tv, Tf, TR>::vertex_iterator
        >                                                  *blue_update_its)
{
    debugl(1, "MeshAlg::RedBlueAlgorithm(): keep_red_outside_part: %d, keep_blue_outside_part: %d.\n",
            keep_red_outside_part, keep_blue_outside_part);
    debugTabInc();

    /* tuple lists used in the algorithm computation */
    std::list<RB_Tuple<Tm, Tv, Tf, TR>> red_tuples, blue_tuples, isecpoly_tuples;

    /* retrieve pairs of potentially intersecting eges and faces,  both
     *  
     *      (red edge / blue face)  and  (blue edge / red face)
     */
    std::vector<EdgeFacePair<Mesh<Tm, Tv, Tf, TR> > > R_edges_B_faces_candidates;
    std::vector<EdgeFacePair<Mesh<Tm, Tv, Tf, TR> > > B_edges_R_faces_candidates;

    /* retrieve candidate lists, where std::pair comparison operators are defined in the
     * straight-forward way form the comparison operators for the components. face pointers are
     * compared like every other integral type */
    Aux::Timing::tick(15);
    debugl(0, "RedBlueAlgorithm(): getting pairs of potentially intersecting edges / faces.\n");

    MeshAlg::getPotentiallyIntersectingEdgeFacePairs(R, B, R_edges_B_faces_candidates, B_edges_R_faces_candidates, 32, 8);

    debugl(0, "RedBlueAlgorithm(): done getting pairs of potentially intersecting edges / faces. time: %5.4f\n\n", Aux::Timing::tack(15));

#ifdef __DEBUG__
    debugl(2, "returned candidate pairs RED edges / BLUE faces: \n");
    debugTabInc();
    for (auto &R_B_cand : R_edges_B_faces_candidates) {
        debugl(2, "red edge:  (%5d, %5d), blue face: %5d\n",
                R_B_cand.vrt1->id(),
                R_B_cand.vrt2->id(),
                R_B_cand.f->id());
    }
    debugTabDec();

    debugl(2, "returned candidate pairs BLUE edges / RED faces: \n");
    debugTabInc();
    for (auto &B_R_cand : B_edges_R_faces_candidates) {
        debugl(2, "blue edge: (%5d, %5d), red face:  %5d\n",
                B_R_cand.vrt1->id(),
                B_R_cand.vrt2->id(),
                B_R_cand.f->id());
    }
    debugTabDec();
#endif

    /* generate tuple lists with extra function to avoid code repetition. pass list of edge intersection information
     * structs by reference to store information about potential complexly intersecting edges (both red and blue). */
    std::list<RedBlue_EdgeIsecInfo<TR>> complex_edge_info_list;
    RedBlue_generateRBTupleList<Tm, Tv, Tf, TR>(R, B, true,  R_edges_B_faces_candidates, red_tuples, complex_edge_info_list);
    RedBlue_generateRBTupleList<Tm, Tv, Tf, TR>(B, R, false, B_edges_R_faces_candidates, blue_tuples, complex_edge_info_list);
    debugl(2, "tuple lists generated..\n");

    /* if complex_edge_info_list is non-empty, there are complexly intersecting edges, maybe both red and blue
     * ones. throw exception to indicate this to the caller, which must split these edges or react with other
     * appropriate measures. */
    if (!complex_edge_info_list.empty()) {
        /* generate exception */
        throw RedBlue_Ex_ComplexEdges<TR>(
                "RedBlueAlgorithm(): two RedBlue_generateRBTupleList() calls discovered complexly intersecting edges.",
                complex_edge_info_list
            );
    }

#ifdef __DEBUG__
    debugl(2, "UNordered RED tuple list:\n");
    debugTabInc();
    for (auto &tuple : red_tuples) {
        debugl(2, " edge: (out, in) = (%5d, %5d) | face: %5d, backward_tri: %5d, forward_tri: %5d\n",
                tuple.edge_out_it->id(), tuple.edge_in_it->id(),
                tuple.face_it->id(),
                tuple.backward_tri_it->id(), tuple.forward_tri_it->id());
    }
    debugTabDec();
#endif

#ifdef __DEBUG__
    debugl(2, "UNordered BLUE tuple list:\n");
    debugTabInc();
    for (auto &tuple : blue_tuples) {
        debugl(2, " edge: (out, in) = (%5d, %5d) | face: %5d, backward_tri: %5d, forward_tri: %5d\n",
                tuple.edge_out_it->id(), tuple.edge_in_it->id(),
                tuple.face_it->id(),
                tuple.backward_tri_it->id(), tuple.forward_tri_it->id());
    }
    debugTabDec();
#endif

    /* if both lists are empty, there can't be any intersection. throw exception. */
    if ( red_tuples.empty() && blue_tuples.empty() ) {
        debugl(0, "redpoints and bluepoints both empty. nothing to do => returning.\n");
        throw RedBlue_Ex_Disjoint("RedBlue_Algorithm(): no red edge intersect a blue face and vice versa => red and blue meshes disjoint.");
    }
    /* if there is no red tuple but at least one blue tuple, then all blue edges intersect the SAME red face, resulting
     * the ''circle'' of affected red faces in the dual of the red mesh to be trivial, which violates an assumption of
     * the RedBlue (union) algorithm. throw exception to indicate this to the caller, who will have to react, e.g. by
     * splitting the red triangle properly. */
    else if ( red_tuples.empty() ) {
        debugl(0, "redpoints empty..\n");
        throw RedBlue_Ex_AffectedCircleTrivial(
                "RedBlue_Algorithm(): red tuple list empty => affected circle of red triangles trivial, i.e. consisting of exactly one red face => split it.", 
                true,
                blue_tuples.front().face_it->id()
            );
    }
    /* symmetric case: no blue tuples, yet at least one red tuple => trivial ''circle'' of blue
     * affected faces. throw exception */
    else if ( blue_tuples.empty() ) {
        debugl(0, "bluepoints empty..\n");
        throw RedBlue_Ex_AffectedCircleTrivial(
                "RedBlue_Algorithm(): blue tuple list empty => affected circle of blue triangles trivial, i.e. consisting of exactly one blue face => split it.", 
                false,
                red_tuples.front().face_it->id()
            );
    }

    /* now red_tuples contains all unordered red tuples and blue_tuples contains all unordered blue tuples. firstly, the
     * lists are each ordered / sorted individually according to the defined neighbourhood relation.  afterwards, the
     * two orderd lists are zipped together to form the completed tuple list uniquely describing the intersection
     * polygon.*/

    /*  cyclically order the computed tuple lists individually. */
    RedBlue_cyclicallyOrderTupleList(R, true,  red_tuples);
    RedBlue_cyclicallyOrderTupleList(B, false, blue_tuples);
    debugl(2, "tuple lists cyclically ordered..\n");

#ifdef __DEBUG__
    debugl(2, "cyclically ordered RED tuple list:\n");
    debugTabInc();
    for (auto &tuple : red_tuples) {
        debugl(2, " edge: (out, in) = (%5d, %5d) | face: %5d, backward_tri: %5d, forward_tri: %5d\n",
                tuple.edge_out_it->id(), tuple.edge_in_it->id(),
                tuple.face_it->id(),
                tuple.backward_tri_it->id(), tuple.forward_tri_it->id());
    }
    debugTabDec();
#endif

#ifdef __DEBUG__
    debugl(2, "cyclically ordered BLUE tuple list:\n");
    debugTabInc();
    for (auto &tuple : blue_tuples) {
        debugl(2, " edge: (out, in) = (%5d, %5d) | face: %5d, backward_tri: %5d, forward_tri: %5d\n",
                tuple.edge_out_it->id(), tuple.edge_in_it->id(),
                tuple.face_it->id(),
                tuple.backward_tri_it->id(), tuple.forward_tri_it->id());
    }
    debugTabDec();
#endif

    /* zip them together to produce the tuple list isecpoly_tuples, which is the desired representation of the
     * intersection polygon */
    debugl(2, "\"zipping\" together cyclically ordered red and blue tuple lists to produce cycically ordered tuple list representing the intersection polygon.\n");

    /* cycle through ordered red tuple list until front() is a red tuple whose intersected blue triangle is NOT the same
     * as intersected triangle of the next red tuple => for this front() tuple, there must be at least one blue tuple
     * in-between in the final list isecpoly_tuples, corresponding to at least blue corners of the intersection polygon
     * between the two corresponding two red points.  cycle the ordered blue tuple list to the correct place => zip them
     * together in a scan. */
    while ( red_tuples.front().face_it == (++red_tuples.begin())->face_it ) {
        red_tuples.push_back( red_tuples.front() );
        red_tuples.pop_front();
    }

    /* now cycle the ordered blue tuple list until the backward triangle of the current front() tuple is the intersected
     * triangle of the first red tuple (which has been cycled there above). */
    while ( blue_tuples.front().backward_tri_it != red_tuples.front().face_it) {
        blue_tuples.push_back(blue_tuples.front());
        blue_tuples.pop_front();
    }

    /* both lists have been prepared to enable linear-scan zipping: the front() red tuple R_f will
     * be followed by the front() blue tuple B_f in the final isecpoly_tuples list, since B_f has
     * R_f.face_it as its backward triangle. 
     *
     * the general invariant after every iteration:
     *
     * for the current front() tuple R_f of the red tuple list, the intersected face R_f.face_it is
     * the backward triangle of the current front() tuple B_f of the blue tuple list.
     *
     * Let R_n denote the first red tuple following R_f in the red tuple list that has a DIFFERENT
     * intersected face than R_f.
     *
     * the iteration creates moves the following sub-list
     *
     *  R_f -> contiguous sublist of all red tuples following R_f that have R_f.face_it as their 
     *          intersected face (can be empty)
     *      -> B_f (guaranteed to exist) 
     *      -> contiguous sublist of all "intermediate" blue tuples following B_f, i.e.  blue
     *          tuples following B_f whose FORWARD triangle is NOT the blue face R_n.face_it
     *          of the red tuple R_n (can be empty) 
     *      -> the last blue tuple B_n, whose foward triangle IS R_n.face_it
     *
     * to isecpoly_tuples in the given order. this process reinstates the invariant after the loop:
     * R_n is the front red() tuple, the front() blue tuple is the one following B_n and  has the
     * R_n.face_it as backward triangle as desired.
     *
     * note that, among other cases, the set of intermediate blue tuples can be empty.
     *
     * the while(1) loop below terminates once the lists run empty. */
    typename Mesh<Tm, Tv, Tf, TR>::face_iterator current_blue_tri, next_blue_tri;
    debugTabInc();
    while (1) {
        current_blue_tri = red_tuples.front().face_it;

        debugl(0, "current_blue_tri: %8d\n", current_blue_tri->id() );
        debugTabInc();
        /* insert all red tuples that have that current_blue_tri as face id */
        while ( red_tuples.front().face_it == current_blue_tri) {
            debugl(0, "red tuple has face_id %d == %d. pushing to final list.\n", red_tuples.front().face_it->id(), current_blue_tri->id() );
            isecpoly_tuples.push_back( red_tuples.front() );
            red_tuples.pop_front();
        }
        debugTabDec();

        /* if red tuple list is empty, append rest of blue tuples and break the loop */
        if ( red_tuples.empty() ) {
            debugl(0, "red tuple list empty => append rest of blue tuple list\n");
            debugTabInc();
            while (!blue_tuples.empty()) {
                debugl(0, "blue tuple: forward_tri_it: %8d | backward_tri_it: %8d\n", blue_tuples.front().forward_tri_it->id(), blue_tuples.front().backward_tri_it->id());
                isecpoly_tuples.push_back( blue_tuples.front() );
                blue_tuples.pop_front();
            }
            debugTabDec();

            /* done. lists have been zipped together */
            break;
        }

        /* now we got a red tuple with different face_id in front(). this is next_blue_tri */
        next_blue_tri = red_tuples.front().face_it;
        debugl(0, "next_blue_tri: %8d\n", next_blue_tri->id() );

        debugTabInc();

        /* now the first blue tuple B_f and all intermediate blue tuples. the first blue tuple has
         * backward_tri_it == current_blue_tri by invariant. the intermediate blue tuples are those
         * whose forward_tri_it is NOT next_blue_tri. while the first blue tuple is guaranteed to
         * exist, the list of intermediate blue tuples might well be empty. */
        while ( blue_tuples.front().forward_tri_it != next_blue_tri) {
            debugl(0, "blue tuple: forward_tri_it: %8d | backward_tri_it: %8d\n", blue_tuples.front().forward_tri_it->id(), blue_tuples.front().backward_tri_it->id());
            isecpoly_tuples.push_back( blue_tuples.front() );
            blue_tuples.pop_front();
        }

        /* and finally the last last blue tuple, whose forward_tri_it IS next_blue_tri */
        isecpoly_tuples.push_back( blue_tuples.front() );
        debugl(0, "blue tuple: forward_tri_it: %8d | backward_tri_it: %8d\n", blue_tuples.front().forward_tri_it->id(), blue_tuples.front().backward_tri_it->id());
        blue_tuples.pop_front();

        debugTabDec();
    }
    debugTabDec();
    debugl(0, "zipping done: ordered intersection polygon list fully constructed.\n");

    /* ---------------- cut holes in B and R -------- */

    /* first, the the red mesh R: cycle the isec poly tuple list until a RED vertex is in front () */
    debugl(2, "cutting hole in RED mesh R.\n");
    while (!isecpoly_tuples.front().red) {
        isecpoly_tuples.push_back( isecpoly_tuples.front() );
        isecpoly_tuples.pop_front();
    }

    /* cut hole in R */
    RedBlue_cutHole(R, true,  isecpoly_tuples, keep_red_outside_part);

    debugl(2, "RED mesh cut.\n");

    /* now the blue mesh B: cycle the isec poly tuple list until a BLUE vertex is in front () */
    debugl(2, "cutting hole in BLUE mesh B.\n");
    while (isecpoly_tuples.front().red) {
        isecpoly_tuples.push_back( isecpoly_tuples.front() );
        isecpoly_tuples.pop_front();
    }

    /* cut hole in B, pass blue_update_its parameters to have the cutting function update the given iterators. if a
     * vertex contained in blue_update_its is deleted, the corresponding iterator is explicitly invalidated, otherwise
     * the iterator is updated to refer to the vertex the cut mesh part of B. NOTE: after the merging step following
     * below, the iterators contained in blue_update_its have to be updated again to refer to the listed vertices as new
     * vertices of the union mesh. since all non-deleted blue vertices are moved from the (cut) rest of B to the cut
     * rest of R during merging, the container changes, and so all iterators would be invalid. to get around this, all
     * iterators are converted to pointers, all vertices are moved, and finally fresh red iterators are generated
     * afterwards. */
    RedBlue_cutHole(B, false, isecpoly_tuples, keep_blue_outside_part, blue_update_its);

    debugl(1, "BLUE mesh cut.\n");

    /* convert list of to-be-updated iterators to a list of vertex pointers if necessary (that is, if blue iterator
     * update is desired by the caller) */
    std::list<typename Mesh<Tm, Tv, Tf, TR>::Vertex *> blue_update_pointers;
    if (blue_update_its) {
        for (auto it = blue_update_its->begin(); it != blue_update_its->end(); ++it) {
            /* NOTE: "it" is a list iterator for a list of vertex iterators => dereferencing "it" yields a reference to
             * a vertex iterator. */
            if (it->explicitlyInvalid()) {
                blue_update_pointers.push_back(NULL);
            }
            else {
                blue_update_pointers.push_back(&(**it));
            }
        }
    }

    debugl(2, "BLUE mesh cut.\n");

    /* new way of doing things: MOVE append B to R while storing a std::list<std::pair<typename Mesh<Tm, Tv, Tf>::Vertex
     * *, typename Mesh<Tm, Tv, Tf>::Vertex *>> of associated corner vertices of the intersection polygon. every pair
     * identifies two vertices (v_r, v_b) with identical positions that represent a corner of the intersection polygon in
     * R and B.
     *
     * although Mesh::moveAppend() invalidates all iterators in the appendED mesh, the pointers remain intact. after the
     * moveAppend() call, valid iterators can be obtained simply by calling Vertex::iterator() on the old vertex
     * pointers from B, which are now vertices of R. 
     *
     * after the moveAppend() call, B is empty. from all pairs, the blue vertices have been moved to R: merge the
     * associated two corner vertices with identical positions topologically (union of adjacent vertices, union of
     * incident faces, replace all ids in neighbourhoods / faces) to produce the final mesh. this prevents space
     * overhead (no copying) and time overhead (no id map association). */

    /* RB_Tuples store iterators, which (at least for the blue mesh) are invalidated by the
     * following Mesh::moveAppend() call => compile the list of pairs of pointers to the red and
     * blue intersection polygon corner vertices. the isecpoly_tuples list is successively erase()d
     * during this process. */
    std::list<
            std::pair<
                typename Mesh<Tm, Tv, Tf, TR>::Vertex *,
                typename Mesh<Tm, Tv, Tf, TR>::Vertex *
            >
        > RB_border_vertex_pairs;

    for (auto it = isecpoly_tuples.begin(); it != isecpoly_tuples.end(); ) {
        /* push pair of pointers into border vertex pair list */
        RB_border_vertex_pairs.push_back( {&(*it->R_new_it), &(*it->B_new_it) } );

        /* delete tuple as we go so as not to waste space. */
        it = isecpoly_tuples.erase(it);
    }

    debugl(2, "moveAppend()ing BLUE mesh B to RED mesh R. B is empty afterwards.\n");
    R.moveAppend(B);

    /* merge all pairs of now "duplicate" intersection polygon corner vertices in R. */
    for (auto &p : RB_border_vertex_pairs) {
        R.mergeUnrelatedVertices(p.first->iterator(), p.second->iterator());
    }

    /* if desired by the caller, fill blue_update_its with (red) iterators of the corresponding vertices in the union
     * mesh (stored in R) or explicitly invalid iterators if the vertex has been deleted during cutting. */
    if (blue_update_its) {
        auto pit = blue_update_pointers.begin();
        auto lit = blue_update_its->begin();

        while (pit != blue_update_pointers.end()) {
            /* if pointer is not NULL, draw fresh red iterator from pointer */
            if (*pit) {
                *lit = (*pit)->iterator();
            }
            /* otherwise explicitly invalidate the iterator */
            else {
                lit->explicitlyInvalidate();
            }
            ++pit, ++lit;
        }
        if (lit != blue_update_its->end()) {
            throw RedBlue_Ex_InternalLogic("RedBlue_Algorithm(): blue_update_its list not yet finished after scan-update from blue_update_pointers list. internal logic error.");
        }
    }

    debugTabDec();
    debugl(1, "MeshAlg::RedBlueAlgorithm(): keep_red_outside_part: %d, keep_blue_outside_part: %d. done.\n",
            keep_red_outside_part, keep_blue_outside_part);

}


template <typename Tm, typename Tv, typename Tf, typename TR>
void
MeshAlg::RedBlueUnion(
    Mesh<Tm, Tv, Tf, TR>                               &R,
    Mesh<Tm, Tv, Tf, TR>                               &B,
    std::vector<
            typename Mesh<Tm, Tv, Tf, TR>::vertex_iterator
        >                                              *blue_update_its)
{
    debugl(0, "MeshAlg::RedBlueUnion()\n");
    debugTabInc();

    /* simple forward to RedBlueAlgorithm: keeping both OUTSIDE parts creates the union mesh */
    MeshAlg::RedBlueAlgorithm(R, B, true, true, blue_update_its);
    
    debugTabDec();
    debugl(0, "MeshAlg::RedBlueUnion(): done.\n");
}


template <typename Tm, typename Tv, typename Tf, typename TR>
void
MeshAlg::RedBlueRedMinusBlue(
    Mesh<Tm, Tv, Tf, TR>                               &R,
    Mesh<Tm, Tv, Tf, TR>                               &B,
    std::vector<
            typename Mesh<Tm, Tv, Tf, TR>::vertex_iterator
        >                                              *blue_update_its)
{
    debugl(0, "MeshAlg::RedBlueRedMinusBlue()\n");
    debugTabInc();

    /* simple forward to RedBlueAlgorithm: for set diffrence, keep outside part of R, keep inside
     * part of B */
    MeshAlg::RedBlueAlgorithm(R, B, true, false, blue_update_its);

    debugTabDec();
    debugl(0, "MeshAlg::RedBlueRedMinusBlue(): done.\n");
}


template <typename Tm, typename Tv, typename Tf, typename TR>
void
MeshAlg::RedBlueIntersection(
    Mesh<Tm, Tv, Tf, TR>                               &R,
    Mesh<Tm, Tv, Tf, TR>                               &B,
    std::vector<
            typename Mesh<Tm, Tv, Tf, TR>::vertex_iterator
        >                                              *blue_update_its)
{
    debugl(0, "MeshAlg::RedBlueIntersection()\n");
    debugTabInc();

    /* simple forward to RedBlueAlgorithm: keeping both INSIDE parts creates the intersection mesh.
     * during the cutting with RedBlue_cutHole(), the inside parts are reoriented consistently.
     * however, only in the case of "set" intersection, the orientation of the result mesh has to
     * be inverted again to produce the usually desired orientation. */
    MeshAlg::RedBlueAlgorithm(R, B, false, false, blue_update_its);
    R.invertOrientation();

    debugTabDec();
    debugl(0, "MeshAlg::RedBlueIntersection(): done.\n");
}


template <typename Tm, typename Tv, typename Tf, typename TR>
inline void
RedBlue_generateRBTupleList(
    Mesh<Tm, Tv, Tf, TR>                                   &X,
    Mesh<Tm, Tv, Tf, TR>                                   &Y,
    bool                                                    X_red,
    std::vector<MeshAlg::EdgeFacePair<Mesh<Tm, Tv, Tf, TR> > > &X_edges_Y_faces_candidates,
    std::list<RB_Tuple<Tm, Tv, Tf, TR>>                    &X_colour_tuples,
    /* RB edge info for exception handling */
    std::list<RedBlue_EdgeIsecInfo<TR>>                    &complex_edge_info_list)
{
    using namespace Aux::Geometry::IntersectionTestResults;

    debugl(2, "RedBlue_generateRBTupleList()\n");
    debugTabInc();

    /* NOTE: the entire function is commented with the assumption that X_red == true, i.e. X is the
     * RED mesh and Y is the BLUE mesh. if X_red == false, the roles of "red" and "blue" must be
     * interchanged. */
    typename Mesh<Tm, Tv, Tf, TR>::vertex_iterator          v_it, u_it, uv_in_it, uv_out_it;
    typename Mesh<Tm, Tv, Tf, TR>::face_iterator            uv_fst_tri_it, uv_snd_tri_it;

    Vec3<TR>                                                u, v, f_0, f_1, f_2;

    /* point of intersection / parametric value on segment lambda */
    Vec3<TR>                                                x;
    TR                                                      x_s, x_t, x_lambda;

    /* red edge e = (u, v) identified by iterators u_it, v_it */
    std::list<typename Mesh<Tm, Tv, Tf, TR>::Face *>        e_Y_candidate_tris;
    uint32_t                                                uv_nisec_faces;

    /* for all red edges e_r = (u, v): extract all potentially intersected blue faces from the
     * R_edges_B_faces_candidates list */
    debugTabInc();
    size_t sz = X_edges_Y_faces_candidates.size();
    for (size_t i = 0; i < sz;)
    {
        /* clear candidate faces and reset intersected face counter */
        e_Y_candidate_tris.clear();
        uv_nisec_faces = 0;

        /* extract the edge e = (u, v) from the current front() element of X_edges_Y_faces_candidates */ 
        debugl(0, "-----------------\n");
        u_it        = X_edges_Y_faces_candidates[i].vrt1->iterator();
        v_it        = X_edges_Y_faces_candidates[i].vrt2->iterator();

        u           = u_it->pos();
        v           = v_it->pos();

        /* extract blue face from front() element of X_edges_Y_faces_candidates and store in candidate list for e */
        e_Y_candidate_tris.push_back(X_edges_Y_faces_candidates[i].f);
        ++i;

        /* as long as the red edge e = (u, v) is the current edge in the SORTED candidate pair list, keep moving the
         * corresponding blue faces from Y into e_candidate_tris */
        while (i < sz && ((X_edges_Y_faces_candidates[i].vrt1->id() == u_it->id()
               && X_edges_Y_faces_candidates[i].vrt2->id() == v_it->id())
            || (X_edges_Y_faces_candidates[i].vrt1->id() == v_it->id()
                && X_edges_Y_faces_candidates[i].vrt2->id() == u_it->id())))
        {
            e_Y_candidate_tris.push_back(X_edges_Y_faces_candidates[i].f);
            ++i;
        }

#ifdef __DEBUG__
        /* pure debug: output candidate tris: */
        debugl(2, "candidate faces for edge (%5d, %5d)\n", u_it->id(), v_it->id());
        debugTabInc();
        for (auto &uv_cand_tri : e_Y_candidate_tris) {
            debugl(2, "%5d.\n", uv_cand_tri->id());
        }
        debugTabDec();
#endif

        /* perform ray-triangle intersections for every candidate and store intersected faces. if
         * the edge e = (u, v) intersects more than one candidate face, a corresponding exception is
         * thrown to indicate the violation of the invariant that every edge is simply intersecting
         * */
        std::vector<uint32_t>   e_isec_tri_ids;
        std::vector<TR>         e_isec_lambdas;

        /* check all blue candidate faces for e_r */
        debugTabInc();
        for (auto &candidate_tri : e_Y_candidate_tris) {
            debugl(2, "checking edge (%5d, %5d) vs candidate face %5d for intersection..\n", u_it->id(), v_it->id(), candidate_tri->id());
            if (candidate_tri->isQuad()) {
                throw RedBlue_Ex_InternalLogic("RedBlue_generateRBTupleList(): quad discovered. operatio on affected quads unsupported (and unsupportable) due to underlying mathematical structure.\n");
            }

            /* get vertex positions of blue candidate triangle */
            candidate_tri->getTriPositions(f_0, f_1, f_2);

            if (Aux::Geometry::rayTriangle(u_it->pos(), v_it->pos(), f_0, f_1, f_2, x, x_s, x_t, x_lambda) == INTERSECTION) {
                debugl(2, "intersection! => creating tuple for candidate face %5d.\n", candidate_tri->id());
                uv_nisec_faces++;

                /* check which endpoint of the red edge {u, v} lies on the outside of the blue mesh.
                 * also check if lambda is too close to zero or one, and throw edge case exception
                 * if it is. */
                if (std::abs(x_lambda) < 1E-10 || std::abs(x_lambda - 1.0) < 1E-10) {
                    debugTabDec();
                    debugTabDec();
                    throw RedBlue_Ex_NumericalEdgeCase(
                            "RedBlue_generateRBTupleList(): edge starting / ending on intersected face.",
                            /* both red and blue meshes are still intact at this point, since no cutting has been
                             * performed as yet. */
                            true,
                            true
                        );
                }

                e_isec_tri_ids.push_back(candidate_tri->id());
                e_isec_lambdas.push_back(x_lambda);
                
                /* get the normal n of the intersected blue triangle. since x is a point on the
                 * triangle, (u-x)*n > 0 implies that u is outside (and v inside), otherwise v is
                 * outside (and u inside). */
                if ( (u - x) * candidate_tri->getNormal() > 0.0) {
                    uv_out_it   = u_it;
                    uv_in_it    = v_it;
                }
                else {
                    uv_out_it   = v_it;
                    uv_in_it    = u_it;
                }

                /* get the two red triangles incident to edge e_r = (u, v). if e_r is non-manifold, this throws an
                 * exception */
                X.getFacesIncidentToManifoldEdge(u_it, v_it, uv_fst_tri_it, uv_snd_tri_it);

                /* compile and insert a new red tuple from the gathered information */
                X_colour_tuples.push_back( 
                        RB_Tuple<Tm, Tv, Tf, TR>(
                            X_red,
                            x,
                            x_lambda, x_s, x_t,
                            uv_out_it, uv_in_it,
                            candidate_tri->iterator(),
                            uv_fst_tri_it, uv_snd_tri_it
                        )
                    );
            }
            /* else: candidate face blue_candidate_tri has not been intersected by edge e_r = (u, v) and is not
             * considered any further.*/

        }
        debugTabDec();

        /* all candidate faces for red edge e = (u, v) have been checked. if more than one blue triangle has been
         * intersected, the edge e_r is not simply intersecting, but complexly intersecting, which violates an
         * assumption for the Red-Blue (Union) algorithm. in this case, the corresponding information is appended to
         * exception_cplx_info_list. an exception containing information about all complex edges from both the red and
         * the blue mesh is thrown by RedBlue_Algorithm if complex edges ever occur.  the caller must then split all
         * complex edge or react with other appropriate measures. */
        if (uv_nisec_faces > 1) {
            debugl(0, "Mesh::RedBlueUnion(): red edge {%5d %5d} is complexely intersecting: %5d > 1 blue faces intersected. throwing exception..\n", v_it->id(), u_it->id(), uv_nisec_faces);
            complex_edge_info_list.push_back(
                RedBlue_EdgeIsecInfo<TR>(
                    /* red edge iff X_red == true */
                    X_red,
                    /* vertex ids of the edge in the red mesh, in the ORIGINAL order, not sorted by orientation, since
                     * lambda depends on the order. original order is {v, u}, which is confusing, since in the struct
                     * its {u, v}, but well.. */
                    v_it->id(), u_it->id(),
                    e_isec_tri_ids,
                    e_isec_lambdas)
                );
        }
        /* red edge e = (u, v) intersects the blue mesh Y exactly once: okay */
        else if (uv_nisec_faces == 1) {
            debugl(3, "edge (%5d %5d) from X intersects exactly one triangle from Y.. okay\n", v_it->id(), u_it->id() );
        }
        /* red edge e = (u, v) does not intersect the blue mesh Y at all. also ok */
        else {
            debugl(3, "edge {%5d %5d} from X does not intersect Y at all.. okay\n", v_it->id(), u_it->id() );
        }
    }
    debugTabDec();

    /* sort() and unique() tuplelist. this uses RB_Tuple::operator<() and RB_Tuple::operator==() */
    X_colour_tuples.sort();
    X_colour_tuples.unique();

    debugTabDec();
    debugl(2, "RedBlue_generateRBTupleList(): done.\n");
}


template <typename Tm, typename Tv, typename Tf, typename TR>
inline void
RedBlue_cyclicallyOrderTupleList(
    const Mesh<Tm, Tv, Tf, TR>             &M,
    bool                                    M_red,
    std::list<RB_Tuple<Tm, Tv, Tf, TR>>    &M_tuples)
{
    /* NOTE: the entire function is commented with the assumption that M is the RED mesh.  If M is
     * the blue mesh, the roles of "red" and "blue" must be interchanged. */
    debugl(2, "RedBlue_cyclicallyOrderTupleList()\n");
    debugTabInc();

    /* create temporary "unordered" list and swap tuples into it */
    std::list<RB_Tuple<Tm, Tv, Tf, TR>> M_tuples_unordered;
    std::swap(M_tuples, M_tuples_unordered);

    /* insert first red tuple from unordered list into (ordered list) list and pop_back() */
    M_tuples.push_front( M_tuples_unordered.back() );
    M_tuples_unordered.pop_back();

    bool                                                    found_nb;
    typename Mesh<Tm, Tv, Tf, TR>::face_iterator            shared_tri;
    RB_Tuple<Tm, Tv, Tf, TR>                               *d;
    typename std::list<RB_Tuple<Tm, Tv, Tf, TR>>::iterator  pit;

    debugl(2, "sorting all tuples from unordere list into ordered list..\n");
    /* while there are still red tuples left to be sorted .. */
    debugTabInc();
    while (!M_tuples_unordered.empty()) {
        /* set d to be the last element of the partially completed ordered list */
        d           = &(M_tuples.back());
        found_nb    = false;

        debugl(2, "searching for correct neighbour for currently last tuple d of (partially completed) ordered list: (out, in) = (%5d, %5d)\n", d->edge_out_it->id(), d->edge_in_it->id() );
        /* find the neighbour nb of currently last added element d of ordered list */
        debugTabInc();
        for (pit = M_tuples_unordered.begin(); pit != M_tuples_unordered.end(); ++pit) {
            debugl(2, "checking (out, in) = (%5d, %5d)\n", pit->edge_out_it->id(), pit->edge_in_it->id() );
            /* if the edges of (*pit) and d share a common vertex, then they are part of the same
             * triangle and hence they are neighbours in the cyclically ordered list. however, there
             * are two neighbours and we need the right orientation in the polygonal line loop
             * describing the polygon of intersection. if we got a neighbour, retrieve their common
             * triangle and check if the edge (u,v) is in the direction of orientation in that
             * triangle. if so, accept the neighbour. otherwise, skip it.
             *
             * NOTE: this way, when the in-between blue points are added below, the intersection
             * polyogon is oriented in such a way that for all RED triangles, the polygon
             *
             * outside point -> red point -> in-between blue points -> next red point (-> other
             * outside point on triangle)
             *
             * is oriented CCW from the outside. */
            if (d->isNeighbour(*pit)) {
                /* get shared triangle */
                shared_tri = M.getSharedTri(d->edge_out_it, d->edge_in_it, pit->edge_out_it, pit->edge_in_it);

                /* get orientation of edge (out_it, in_it) = (u, v) in shared tri.
                 *
                 * RED case:    the forward triangle has an orientation that MATCHES,
                 *              i.e. contains (u, v), NOT (v, u).
                 *
                 * BLUE case:   the forward triangle has an orientation that OPPOSES,
                 *              i.e. contains (v, u), NOT (u, v).
                 *
                 * => M_red must be used appropriately in the predicate below. */
                if  (Aux::Logic::lequiv(
                            M_red,
                            shared_tri->getTriEdgeOrientation( d->edge_out_it->id(), d->edge_in_it->id() )
                        )
                    )
                {
                    debugl(4, "got neighbour, correct orientation..\n");

                    /* found the neighbour.. */
                    found_nb = true;

                    /* the "forward" triangle of d must be the shared triangle, which must also be
                     * the "backward" triangle of pit. swap if necessary. */
                    if (d->forward_tri_it != shared_tri) {
                        std::swap(d->forward_tri_it, d->backward_tri_it);
                    }
                    if (pit->backward_tri_it != shared_tri) {
                        std::swap(pit->forward_tri_it, pit->backward_tri_it);
                    }

                    /* append *pit to (ordered) tuple list and erase it from unordered tuple list*/
                    M_tuples.push_back(*pit);
                    M_tuples_unordered.erase(pit);
                    break;
                }
                else {
                    debugl(4, "got neighbour, but wrong orientation..\n");
                }
            }
        }
        debugTabDec();

        /* no neighbour found. check if redpoints already describes a closed ring of red faces => then the number of
         * intersection polygons is >= 2, throw matching exception to indicate this to the caller. */
        if (!found_nb) {
            if (M_tuples.back().isNeighbour(M_tuples.front())) {
                debugTabDec();
                debugTabDec();
                throw RedBlue_Ex_NumIsecPoly("RedBlue_cyclicallyOrderTupleList(): > 1 distinct intersection polygons discovered.");
            }
            else {
                debugTabDec();
                debugTabDec();
                throw RedBlue_Ex_InternalLogic("RedBlue_cyclicallyOrderTupleList(): discovered tuple without neighbour, but no intersection polygon has yet been completed => internal logic error.");
            }
        }
    }
    debugTabDec();

    /* check if first and last of cyclically ordered list are neighbours.. */
    if (!M_tuples.front().isNeighbour( M_tuples.back() ) ) {
        debugTabDec();
        throw RedBlue_Ex_InternalLogic("RedBlue_cyclicallyOrderTupleList(): first and last tuple of \"cyclically ordered\" red tuple list are not neighbours. this must not happen. internal logic error.");
    }

    debugTabDec();
    debugl(2, "RedBlue_cyclicallyOrderTupleList(): done.\n");
}


/* cut a hole in a mesh using the intersection curve "isecpoly_tuples" on the surface.  bool red indicated whether the
 * mesh is red or blue, i.e. whether the red points or the blue points are points on the edges of M (and on the faces of
 * the other mesh).
 *
 * isecpoly_tuples must be cyclically ordered and cycled in such a way that it starts with a red vertex iff red == true.
 * furthermore, the start_out_id identifies the vertex on the edge belonging to the first point that is supposed to be
 * used as the "outside" vertex, determining the orientation. by switching this, different boolean set operations can be
 * performed: union, intersection ,etc. the mesh M is cut in such a way that v_out_id remains on the produced mesh, and
 * everything beyond the border curve defined by isecpoly_tuples is deleted. */
template <typename Tm, typename Tv, typename Tf, typename TR>
inline void
RedBlue_cutHole(
    Mesh<Tm, Tv, Tf, TR>                                       &M,
    bool                                                        M_red,
    std::list<RB_Tuple<Tm, Tv, Tf, TR>>                        &isecpoly_tuples,
    bool                                                        keep_outside_part,
    std::vector<
            typename Mesh<Tm, Tv, Tf, TR>::vertex_iterator
        >                                                      *update_vertex_its,
    std::list<
            typename Mesh<Tm, Tv, Tf, TR>::face_iterator
        >                                                      *new_tri_its)
{
    using namespace Aux::Geometry;
    using namespace Aux::Geometry::IntersectionTestResults;

    debugl(3, "RedBlue_cutHole(): red mesh: %d.\n", M_red);
    debugTabInc();

    uint32_t                                        i, err;
    typename Mesh<Tm, Tv, Tf, TR>::vertex_iterator  v_out_it, v_in_it;


    /* check if front vertex is of the right colour */
    if ( !Aux::Logic::lequiv(M_red, isecpoly_tuples.front().red) ) {
        throw RedBlue_Ex_InternalLogic("RedBlue_cutHole(): first point in tuple list has different colour than mesh to be cut.");
    }

    /* if iterator update is desired by the caller, convert list of update iterators to list of pointers, which remain
     * intact even if the iterators are invalidated. */
    std::list<typename Mesh<Tm, Tv, Tf, TR>::Vertex *> update_vertex_pointers;
    if (update_vertex_its) {
        for (auto it = update_vertex_its->begin(); it != update_vertex_its->end(); ++it) {
            update_vertex_pointers.push_back(&(**it));
        }
    }

    /* get id of point on the edge belonging to the first point in the isecpoly_tuples list that is "inside" of the
     * created hole => his connected component will be deleted after the cutting and retriangulations are done */
    RB_Tuple<Tm, Tv, Tf, TR> &first_point = isecpoly_tuples.front();
    if (keep_outside_part) {
        v_out_it    = first_point.edge_out_it; 
        v_in_it     = first_point.edge_in_it;
    }
    else {
        v_out_it    = first_point.edge_in_it;
        v_in_it     = first_point.edge_out_it;
    }

    /* traversal id used during the function to set traversal states for the final "cutting" */
    uint32_t traversal_id = M.getFreshTraversalId();

    /* NOTE: for the following text, it is assumed that red == false. if red == true, substitute
     * "red" for "blue" :*/

    /* the blue points are on edges of A. walk from blue point to blue point and extract the generally non-convex planar
     * polygon that is created by the two consecutive blue points, the "outside" corners of the triangle both blue
     * points lie on and the intermittant red points.  the orientation is uniquely defined for the entire walk if for
     * only ONE edge of A that intersects M, the corresponding corner that is "inside" of M is known. from that point
     * onwards, everything is implied by the orientation of the triangles and the incidence relations.  the resulting
     * planar polyons are subsequently triangulated in 2d coordinates, maintaining the indices. the obtained triangles
     * must then only be added to M, while the original affected triangles and the "inside" connected component
     * encircled by these are deleted from the mesh. */

    /* add all corners of the intersection polygon as vetices of M and set iterators inside RB_Tuple struct, of course
     * depending on M_red. also compute number of such red and blue corner vertices, respectively. */
    uint32_t                                        npoints = 0, n_total_inbetween_points = 0;
    typename Mesh<Tm, Tv, Tf, TR>::vertex_iterator  vnew_it;
    for (auto &p : isecpoly_tuples) {
        /* add corner as vertex of M */
        vnew_it = M.vertices.insert(p.x);

        /* M is red */
        if (!M_red) {
            p.B_new_it = vnew_it;
        }
        else {
            p.R_new_it = vnew_it;
        }

        if ( Aux::Logic::lequiv(M_red, p.red) ) {
            npoints++;
        }
        else {
            n_total_inbetween_points++;
        }
    }

    debugl(0, "npoints: %d, n_total_inbetween_points: %d, size(): %d\n", npoints, n_total_inbetween_points, isecpoly_tuples.size());

    /* list to store red points between currently processed pair of points. if (M_red == false), this needs to be
     * REVERSED due to orientation after it is computed: the orientation of the intesection polygon was determined to
     * that consecutive red points are traversed in the direction of increasing index when orienting the polygons in M
     * in CCW fashion => for bue points, the csg points have to be traversed in reverse, i.e. in order of decreasing
     * index, to get ccw orientation */
    std::list<RB_Tuple<Tm, Tv, Tf, TR>>                     in_between_points;
    typename std::list<RB_Tuple<Tm, Tv, Tf, TR>>::iterator  pit, pit_current, pit_next;
    RB_Tuple<Tm, Tv, Tf, TR>                               *p, *pnext;

    typename Mesh<Tm, Tv, Tf, TR>::vertex_iterator          p_vin_it, p_vout_it, pnext_vout_it;
    typename Mesh<Tm, Tv, Tf, TR>::face_iterator            shared_tri_it;
    typename Mesh<Tm, Tv, Tf, TR>::vertex_iterator          shared_vertex_it, remaining_vertex_it;

    std::list< std::vector<Vertex2d> >                      planar_polygons;
    Vec3<TR>                                                shtri_x, shtri_y;

    Vec2                                                    v2d;
    Vec3<TR>                                                v3d;
    std::vector<Vertex2d>                                   planar_poly;
    std::vector<Tri2d>                                      planar_poly_triangulation;

    /* append current front() to back to avoid ugly wrap around case. both front() and back() are copies of the same
     * point now, so the wrap around pair (last, first) can be handled within the normal list iteration. */
    isecpoly_tuples.push_back( isecpoly_tuples.front() );

    /* set current point iterator */
    pit_current = isecpoly_tuples.begin();
    p           = &(isecpoly_tuples.front());

    /* the curve on the mesh surface of M defines two regions. v_out_id defines the vertex of the edge belonging to the
     * first point that is supposed to "survive" => cut the hole and delete the other region NOT containing v_out_id */
    p_vout_it = v_out_it;

    debugl(0, "splitting shared triangles for all consecutive pairs of points..\n");

    /* check pairs (i, i+1) of blue points for i = 0..(nblue - 1). note that we've appended a copy of the first point to
     * the end of the list to be able to scan linearly without wrap around. */
    debugTabInc();
    for (i = 0; i < npoints; i++) {
        debugl(0, "pair (%5d, %5d) of points. npoints: %5d.\n", i, i+1, npoints);
        debugTabInc();

        /* find next point, store in-between points */
        pit = pit_current;
        ++pit;
        in_between_points.clear();

        /* if we're cutting the blue mesh, skip red points, reverse them due to orientation */
        if (!M_red) {
            debugTabInc();

            debugl(0, "cutting blue mesh => collecting in-between red points and reversing..\n");
            while (pit->red) {
                debugl(0, "skipping red point, adding to in_between_points.\n");
                /* directly "reverse" by pushing at the front */
                in_between_points.push_front(*pit);
                ++pit;
            }
            debugl(0, "got next blue point. red points inbetween: %d\n", in_between_points.size() );

            debugTabDec();
        }
        /* if we're cutting the red mesh, skip the blue points. leave the orientation as is, see above */
        else {
            debugl(0, "cutting red mesh => collecting in-between blue points WITHOUT reversing..\n");

            debugTabInc();
            while (!pit->red) {
                debugl(0, "skipping blue point, adding to in_between_points.\n");
                /* push at back, keep orientation */
                in_between_points.push_back(*pit);
                ++pit;
            }
            debugl(0, "got next red point. blue points inbetween: %d\n", in_between_points.size() );

            debugTabDec();
        }


        /* got the next point */
        pit_next    = pit;
        pnext       = &(*pit_next);

        /* we already got p_vout_it from the last iteration or from the invariant init. only set p_vin_it */
        if (p->edge_out_it != p_vout_it) {
            p_vin_it = p->edge_out_it;
        }
        else {
            p_vin_it = p->edge_in_it;
        }

        /* get shared triangle id and its orthonormal basis */
        shared_tri_it = M.getSharedTri(
                            p->edge_out_it,
                            p->edge_in_it,
                            pnext->edge_out_it,
                            pnext->edge_in_it);

        shared_tri_it->getTriOrthonormalBase2d(shtri_x, shtri_y);

        debugl(0, "shared triangle id: %5d\n", shared_tri_it->id() );
        debugl(0, "triangle orthonormal basis..\n");
        shtri_x.print_debugl(0);
        shtri_y.print_debugl(0);

        /* determine the shared endpoint of the edges for b and bnext. two cases:
         *
         * i) M_red == false => blue, intersectED mesh. mind the orientation gap
         *
         *      1. the shared endpoint is b_vout_id => polygon is 
         *
         *          ( p_vout_id ++ pnext ++ red_points ++ p )
         *
         *      => the "outside point" on the edge of b_next remains the same, since it is shared.
         *
         *      2. the shared endpoint is b_vinid => polygon is 
         *
         *          (p_vout_id ++ remaining_vertex ++ p_next ++ red_points ++ p )
         *
         *      => the "outside point" on the edge of p_next is the other remaining vertex, since shared vertex is
         *      inside.
         *
         * i) M_red == true => red, intersectING mesh. mind the orientation gap.
         *
         *      1. the shared endpoint is b_vout_id => polygon is 
         *
         *          ( p_vout_id ++ p ++ blue_points ++ pnext )
         *
         *      => the "outside point" on the edge of p_next remains the same, since it is shared.
         *
         *      2. the shared endpoint is b_vin_id => polygon is 
         *
         *          (p_vout_id ++ p ++ blue_points ++ p_next ++ remaining_vertex)
         *
         *      => the "outside point" on the edge of b_next is the other remaining vertex, since shared vertex is
         *      inside.  */
        Mesh<Tm, Tv, Tf, TR>::Face::getTriSharedAndRemainingVertex(
                p->edge_out_it, p->edge_in_it,
                pnext->edge_out_it, pnext->edge_in_it,
                shared_vertex_it, remaining_vertex_it
            );

        debugl(0, "generating planar polygon..\n");

        /* generate planar polygon for current shared triangle in A */
        planar_poly.clear();
        planar_poly_triangulation.clear();

        /* if we're cutting the blue mesh */
        if (!M_red) {
            debugl(0, "cutting BLUE mesh.\n");
            debugTabInc();

            /* compute first 2d vertex by projecting the position of the outside vertex onto the triangle using the
             * computed 2s orthonormal base made up out of the two 3d vectors shtri_x and shtri_y */
            debugl(0, "BLUE MESH: first, v_out_id..\n");
            v3d     = p_vout_it->pos();
            v2d[0]  = v3d * shtri_x;
            v2d[1]  = v3d * shtri_y;
            planar_poly.push_back( Vertex2d(p_vout_it->id(), v2d) );

            /* insert remaining blue point if necessary, see comment above */
            if (shared_vertex_it == p_vin_it) {
                debugl(0, "BLUE MESH: shared vertex is INSIDE vertex on edge of p => add remaining vertex\n");
                /* remaining vertex */
                v3d     = remaining_vertex_it->pos();
                v2d[0]  = v3d * shtri_x;
                v2d[1]  = v3d * shtri_y;
                planar_poly.push_back( Vertex2d(remaining_vertex_it->id(), v2d) );

                /* set bnext_vout_id for shift */
                pnext_vout_it = remaining_vertex_it;
            }
            else {
                /* otherwise, the vertex p_vout_it is also the outside vertex on the edge of pnext */
                debugl(0, "BLUE MESH: shared vertex is OUTSIDE vertex on edge of p => fine.\n");
                pnext_vout_it = p_vout_it;
            }

            debugl(0, "BLUE MESH: now, pnext..\n");
           
            /* pnext. blue mesh! */
            v3d     = pnext->x;
            v2d[0]  = v3d * shtri_x;
            v2d[1]  = v3d * shtri_y;
            planar_poly.push_back( Vertex2d(pnext->B_new_it->id(), v2d) );

            debugl(0, "BLUE MESH: now, %5d in-between red vertices..\n", in_between_points.size() );

            /* now all red points, which are already reversed to be consistent with the orientation */
            for (auto &ibp : in_between_points) {
                v3d     = ibp.x;
                v2d[0]  = v3d * shtri_x;
                v2d[1]  = v3d * shtri_y;

                /* blue mesh, red points! */
                planar_poly.push_back( Vertex2d(ibp.B_new_it->id(), v2d) );
            }

            /* and finally p itself (blue mesh!) */
            debugl(0, "BLUE MESH: finally p..\n");
            v3d     = p->x;
            v2d[0]  = v3d * shtri_x;
            v2d[1]  = v3d * shtri_y;
            planar_poly.push_back( Vertex2d(p->B_new_it->id(), v2d) );
            
            debugTabDec();
        }
        /* if we're cutting the red mesh */
        else {
            /* first v_out_id. red mesh! */
            debugl(0, "cutting RED mesh..\n");
            debugTabInc();

            debugl(0, "RED MESH: first, v_out_id..\n");
            v3d     = p_vout_it->pos();
            v2d[0]  = v3d * shtri_x;
            v2d[1]  = v3d * shtri_y;
            planar_poly.push_back( Vertex2d(p_vout_it->id(), v2d) );

            /* then p. red mesh! */
            debugl(0, "RED MESH: p..\n");
            v3d     = p->x;
            v2d[0]  = v3d * shtri_x;
            v2d[1]  = v3d * shtri_y;
            planar_poly.push_back( Vertex2d(p->R_new_it->id(), v2d) );

            /* then all blue points, which haven't been reversed to be consistent with the
             * orientation */
            debugl(0, "RED MESH: now, %5d, in-between blue vertices..\n", in_between_points.size() );
            /* in-between points, which are already sorted in the right order. if in-between points are
             * red, they have been reversed, otherwise not. mind the orientation .. :D */
            for (auto &ibp : in_between_points) {
                v3d     = ibp.x;
                v2d[0]  = v3d * shtri_x;
                v2d[1]  = v3d * shtri_y;

                /* red mesh, blue points! */
                planar_poly.push_back( Vertex2d(ibp.R_new_it->id(), v2d) );
            }

            /* pnext. red mesh! */
            debugl(0, "RED MESH: pnext..\n");
            v3d     = pnext->x;
            v2d[0]  = v3d * shtri_x;
            v2d[1]  = v3d * shtri_y;
            planar_poly.push_back( Vertex2d(pnext->R_new_it->id(), v2d) );

            /* insert remaining blue point if necessary, see comment above */
            if (shared_vertex_it == p_vin_it) {
                debugl(0, "RED MESH: shared vertex is INSIDE vertex on edge of p => add remaining vertex\n");
                /* remaining vertex */
                v3d     = remaining_vertex_it->pos();
                v2d[0]  = v3d * shtri_x;
                v2d[1]  = v3d * shtri_y;
                planar_poly.push_back( Vertex2d(remaining_vertex_it->id(), v2d) );

                /* set bnext_vout_id for shift */
                pnext_vout_it = remaining_vertex_it;
            }
            else {
                /* otherwise, the vertex p_vout_it is also the outside vertex on the edge of bnext */
                debugl(0, "RED MESH: shared vertex is OUTSIDE vertex on edge of p => fine.\n");
                pnext_vout_it = p_vout_it;
            }
            debugl(0, "done..\n");

            debugTabDec();
        }

        /* call 2d triangulation method */
        debugl(2, "calling 2d triangulation algorithm on projected polygon. planar_poly.size(): %ld\n", planar_poly.size());
        err = Aux::Geometry::triangulateSimplePlanarPolygon(planar_poly, planar_poly_triangulation);
        switch (err) {
            case EDGE_CASE:
                /* the red mesh is never intact if this exception is thrown, since it has either been successfully cut
                 * or the exception happens during cutting the red mesh.  however, if M is the red mesh, i.e. M_red ==
                 * true, then the blue mesh is still intact.  otherwise both meshes are no longer intact when the
                 * exception is thrown. => B_intact is true iff M_red is true. */
                throw RedBlue_Ex_NumericalEdgeCase(
                        "edge case during (planar) outside polygon triangulation..",
                        false,
                        M_red
                        );
                break;

            case TEST_INCONCLUSIVE:
                throw RedBlue_Ex_NumericalEdgeCase(
                        "inconclusive test result during (planar) outside polygon triangulation..",
                        false,
                        M_red
                    );
                break;

            case SUCCESS:
                break;

            default:
                throw RedBlue_Ex_InternalLogic("RedBlue_cutHole(): Aux::Geometry::triangulateSimplePlanarPolygon retunred invalid result code.");
        }
        debugl(2, "2d triangulation algorithm done.\n");

        debugl(2, "calling 2d delaunay on initial triangulation..\n");

        /* convert the returned triangulation to a map */
        std::map<uint32_t, Vertex2d>    poly_map;
        for (auto v : planar_poly) {
            poly_map[v.id] = v;
        }
        delaunay2d(poly_map, planar_poly_triangulation);
        debugl(2, "2d delaunay done.\n");

        /* add all new triangles to the mesh and save them their ids return list new_tri_ids, if
         * present .. */
        debugl(2, "adding faces of of triangulation for affected triangle to mesh..\n");
        typename Mesh<Tm, Tv, Tf, TR>::face_iterator tri_it;
        debugTabInc();
        for (auto &tri : planar_poly_triangulation) {
            debugl(2, "triangle (%5d, %5d, %5d)\n", tri.v0_id, tri.v1_id, tri.v2_id);
            tri_it = M.faces.insert(tri.v0_id, tri.v1_id, tri.v2_id);

            /* if desired by the caller (new_tri_ids != NULL), return new triangle ids */
            if (new_tri_its) {
                new_tri_its->push_back(tri_it);
            }

            /* the new triangle has correct orientation, independent of whether the inside or
             * outside component will be kept. if the inside component should be kept however, all
             * other "inside" faces but the new ones have inverted orientation. mark the new faces
             * by setting their traversal state to DONE, so that these will not be returned by the
             * inside "connected component" traversal at the end of the function */
            tri_it->setTraversalState(traversal_id, Mesh<Tm, Tv, Tf, TR>::TRAV_DONE);
        }
        debugTabDec();

        /* shift values */
        pit_current = pit_next;
        p           = pnext;
        p_vout_it   = pnext_vout_it;

        debugTabDec();
    }
    debugTabDec();

    /* pop the front point copy from the back */
    isecpoly_tuples.pop_back();

    /* delete all old affected triangles and the entire connected component "encircled" by these triangles.  more
     * precisely: if parameter keep_outside_part == true, this means all vertices of the "connected component"
     * containing all INSIDE vertices of affected triangles and whose border is the ring of OUTSIDE vertices of all
     * affected triangles.  if keep_outside_part == false, roles of inside and outside are reversed.  the corresponding
     * vertices, and all faces containing any such vertex, are then deleted.
     *
     * as a matter of fact, that's also the way this "connected component" is computed: get a fresh traversal id, set
     * the traversal state of all inside / outside vertices (depending on keep_outside_part) of all affected triangles
     * to DONE, compute the "connected component" of an arbitrary outside / inside vertex with the outlined boundary
     * condition, delete all reachable vertices (this excludes the outside vertices, since these are initialized to
     * BFS_DONE prior to the actual traversal and are hence not counted as being reachable). note also that none of the
     * "new" corner vertices of the intersection polygon that have been above is reachable, and therefore no "new"
     * triangle created during outside polygon triangulation is ever touched here.  */


    debugl(2, "deleting \"%s\" connected component\" encircled by affected triangles (including the latter) with adapted traversal.\n", (keep_outside_part == true) ? "inside" : "outside");

    /* start traversal from inside vertex of an arbitrary tuple of matching colour, retrieve list of "reachable"
     * vertices with the above boundary condition and erase them from M, which also erases all faces containing any of
     * these vertices in the process. */
    std::list<typename Mesh<Tm, Tv, Tf, TR>::Vertex *> dead_vertices_walking;

    /* if we keep the outside part, delete all inside vertices */
    if (keep_outside_part) {
        debugTabInc();
        for (auto &p : isecpoly_tuples) {
            if (Aux::Logic::lequiv(M_red, p.red)) {
                debugl(3, "setting traversal state of OUTSIDE vertex %5d to DONE.\n", p.edge_out_it->id());
                p.edge_out_it->setTraversalState(traversal_id, Mesh<Tm, Tv, Tf, TR>::TRAV_DONE);
            }
        }
        debugTabDec();
        M.getConnectedComponent(isecpoly_tuples.front().edge_in_it, traversal_id, &dead_vertices_walking, NULL);
    }
    /* if we keep the inside part, delete all outside vertices AND invert the orientation of all faces on the inside
     * part. */
    else {
        debugTabInc();
        for (auto &p : isecpoly_tuples) {
            if (Aux::Logic::lequiv(M_red, p.red)) {
                debugl(3, "setting traversal state of INSIDE vertex %5d to DONE.\n", p.edge_out_it->id());
                p.edge_in_it->setTraversalState(traversal_id, Mesh<Tm, Tv, Tf, TR>::TRAV_DONE);
            }
        }
        debugTabDec();

        M.getConnectedComponent(isecpoly_tuples.front().edge_out_it, traversal_id, &dead_vertices_walking, NULL);

        /* get all ald "inside" faces, i.e. all faces of the "inside connected component" but the ones that have been
         * added during polygon triangulation. after their insertion above, all newly created inside triangles have been
         * marked as DONE for traversal_id.  after resetting the traversal states of all inside border vertices to
         * UNSEEN, a traversal from the inside vertex of an arbitrary tuple of matching colour will return only the
         * desired "old" inside faces. */ 
        debugTabInc();
        for (auto &p : isecpoly_tuples) {
            if (Aux::Logic::lequiv(M_red, p.red)) {
                debugl(3, "resetting traversal state of INSIDE vertex %5d to UNSEEN.\n", p.edge_out_it->id());
                p.edge_in_it->setTraversalState(traversal_id, Mesh<Tm, Tv, Tf, TR>::TRAV_UNSEEN);
            }
        }
        debugTabDec();

        std::list<typename Mesh<Tm, Tv, Tf, TR>::Face *> old_inside_faces;
        M.getConnectedComponent(isecpoly_tuples.front().edge_in_it, traversal_id, NULL, &old_inside_faces);

        /* invert face orientation on old faces */
        debugl(2, "keep_outside_part == false: inverting face orientation of all \"old\" inside faces.\n");
        debugTabInc();
        for (auto &f : old_inside_faces) {
            debugl(3, "inverting orientation of \"old\" inside face %5d.\n", f->id());
            f->invertOrientation();
        }
        debugTabDec();
        debugl(2, "orientation of \"old\" inside faces inverted.\n");
    }

    /* erase all "dead vertices walking". if desired by the caller, NULL all pointers of "dead vertices walking" in
     * update_vertex_pointers and compute update_vertex_its list: if a pointer is not null, generate fresh iterator,
     * otherwise explicitly invalidate the iterator. the order in update_vertex_its given by the caller is completely
     * preserved. */
    if (update_vertex_its) {
        for (auto v : dead_vertices_walking) {
            auto it = Aux::Alg::listFind(update_vertex_pointers, v);
            /* to-be erased vertex v has been found in update_vertex_pointers => NULL the entry of v */
            if (it != update_vertex_pointers.end()) {
                *it = NULL;
            }
        }

        /* write update_vertex_its by converting entries in update_vertex_pointers */
        auto pit = update_vertex_pointers.begin();
        auto lit = update_vertex_its->begin();
        while (pit != update_vertex_pointers.end()) {
            if (*pit) {
                *lit = (*pit)->iterator();
            }
            else {
                lit->explicitlyInvalidate();
            }
            ++pit, ++lit;
        }

        if (lit != update_vertex_its->end()) {
            throw RedBlue_Ex_InternalLogic("RedBlue_cutHole(): update_vertex_its list not yet finished after scan-update from update_vertex_pointers list. internal logic error.");
        }
    }

    /* delete all "dead vertices walking". */
    for (auto &v : dead_vertices_walking) {
        M.vertices.erase(v->iterator());
    }

    debugl(2, "\"%s\" connected component\" deleted.\n", (keep_outside_part == true) ? "inside" : "outside");

    debugTabDec();
    debugl(3, "RedBlue_cutHole(): done, red mesh: %d.\n", M_red);
}


/* ---------------------------------------------------------------------------------------------- */
/*                                                                                                */
/*      post-processing: greedy edge collapse mesh optimisation / simplification                  */
/*                                                                                                */
/* ---------------------------------------------------------------------------------------------- */
/* forward declaration of static function used in MeshAlg::greedyEdgeCollapsePostProcessing */
template <typename Tm, typename Tv, typename Tf, typename R>
inline R
GEC_getAvgAreaOfPermissibleSurroundingTriangles(
    typename Mesh<Tm, Tv, Tf, R>::face_iterator     tri_it,
    R const                                        &max_ar,
    uint32_t                                        depth);

/* greedy edge collapsing of shortest edge of triangles sorted by aspect ratio.  additionally,
 * before reinserting affected triangles, check if their size is within the average of the
 * surrounding trinagles that have a valid ar, i.e. an ar < max_ar. if the triangle is too small,
 * reinsert it. to prevent the collapsing from becoming uncontrolled, stop if the some factor times
 * the average area is exceeded. starting testing with factor 1.0, i.e. don't collpase any further
 * if the triangles have reached a typical area.. */
template <typename Tm, typename Tv, typename Tf, typename R>
void
MeshAlg::greedyEdgeCollapsePostProcessing(
    Mesh<Tm, Tv, Tf, R>    &M,
    R const                &alpha,
    R const                &lambda,
    R const                &mu,
    uint32_t                d)
{
    using namespace Aux::Timing;
    using Aux::Numbers::inf;

    debugl(0, "MeshAlg::greedyEdgeCollapsePostProcessing(): alpha: %f, lambda: %f, mu: %f, d: %d\n", alpha, lambda, mu, d);
    debugTabInc();

    /* triangulate quads, if there are any.. */
    M.triangulateQuads();

    tick(15);

    PriorityQueue<R, uint32_t>                  Q;
    std::pair<R, uint32_t>                      q_min;
    uint32_t                                    tri_id;
    typename Mesh<Tm, Tv, Tf, R>::face_iterator tri_it, uv_other_face_it;
    R                                           tri_ar, tri_area, tri_avg_surrounding_area;

    std::map<uint32_t, R>                       avg_surrounding_area;
    std::list<uint32_t>                         unsafe_tris;
    std::list<uint32_t>                         unsafe_tris_prev_iter;
    std::list<uint32_t>::iterator               lit;


    /* NOTE: collapsing edges only deletes faces => no face id can be freed and retaken by another face. when collapsing
     * an edge, two triangles get deleted => locate and delete them from the map. furthermore, the aspect ratios of all
     * faces incident to the new vertex will change => update them. */

    debugl(1, "searching for \"poor\" triangles / computing average neighbourhood triangle areas. this may take some time..\n");

    /* define processing predicate for convenience */
    #define proc(ar, area, avg_nbhd_area, alpha, lambda, mu) (area < mu * avg_nbhd_area && (ar >= alpha || area < lambda * avg_nbhd_area) )

    /* insert all "poor" triangles into Q, i.e. fill Q with the triangles that need processing */
    debugTabInc();
    for (auto &tri : M.faces) {
        /* get triangle's aspect ratio, area and average over permissible triangles in the
         * d-neighbourhood of currently processed triangle */
        tri_ar                      = tri.getTriAspectRatio();
        tri_area                    = tri.getTriArea();
        tri_avg_surrounding_area    = GEC_getAvgAreaOfPermissibleSurroundingTriangles<Tm, Tv, Tf>(tri.iterator(), alpha, d);
        avg_surrounding_area.insert( {tri.id(), tri_avg_surrounding_area} );

        debugl(3, "face %6d, ar: %10.5f, area: %10.5f, avg area in d-neighbourhood: %10.5f, d = %3d..\n", tri.id(), tri_ar, tri_area, tri_avg_surrounding_area, d);

        if ( proc(tri_ar, tri_area, tri_avg_surrounding_area, alpha, lambda, mu) ) {
            debugl(1, "adding face %6d, ar: %10.5f to min-heap Q..\n", tri.id(), tri_ar);
            Q.insert( { -tri_ar, tri.id() } );
        }

    }
    debugTabDec();
    debugl(1, "done. processing min-heap..\n");

    /* fix-point iteration: the queue might run out of elements but unsafe_tris is non-empty => reinsert them all and
     * rerun the algorithm until we reach a fixpoint, i.e. the unsafe-tri's before the run are equal to the unsafe tris
     * after the run. they are initialized to the empty list */
    typename Mesh<Tm, Tv, Tf, R>::vertex_iterator       u_it, v_it, w_it;
    typename Mesh<Tm, Tv, Tf, R>::face_iterator         us_tri_it;
    uint32_t                                            us_tri_id, w_inc_tri_id;
    std::list<
            typename Mesh<Tm, Tv, Tf, R>::face_iterator
        >                                               w_fstar;

    R                                                   w_inc_tri_ar, w_inc_tri_area, w_inc_tri_avg_surrounding_area;
    R                                                   us_tri_ar, us_tri_area, us_tri_avg_surrounding_area;

    bool                                                collapse_safe;
    bool                                                fixed_point;
    uint32_t                                            iter;

    /* enter fixed point loop */
    fixed_point = false;
    iter        = 0;
    while (!fixed_point) {
        debugl(2, "no fixed point reached yet => performing iteration %2d\n", iter);
        unsafe_tris_prev_iter = unsafe_tris;
        unsafe_tris.clear();

        /* while Q is not empty */
        debugTabInc();
        while (!Q.empty()) {
            /* get current "top" triangle T with maximum aspect ratio (note the - sign used when inserting AR values
             * into the MIN heap) from Q and deleteMin() */
            q_min       = Q.top();
            Q.deleteMin();

            /* get aspet ratio and id */
            tri_ar      = -q_min.first;
            tri_id      = q_min.second;

            /* if triangle T referred to by tri_id still exists in M, process it */
            tri_it = M.faces.find(tri_id);
            if (tri_it != M.faces.end()) {
                if ( tri_ar != tri_it->getTriAspectRatio()) {
                    debugl(0, "MeshAlg::greedyEdgeCollapsePostProcessing(): ar stored in top() element of Q (%10.5f) does not match return value of getTriAspectRatio() (%10.5f)\n", tri_ar, tri_it->getTriAspectRatio());
                    throw("tri_ar from Q doesn't match tri_ar from triangle..\n");
                }
                tri_area                    = tri_it->getTriArea();
                tri_avg_surrounding_area    = avg_surrounding_area[tri_it->id()];

                /* skip T if it does not require processing */
                if ( proc(tri_ar, tri_area, tri_avg_surrounding_area, alpha, lambda, mu) ) {
                //if (tri_area < mu * tri_avg_surrounding_area && ( tri_ar >= alpha || tri_area < lambda * tri_avg_surrounding_area) )
                    debugl(2, "triangle %d needs processing.. collapsing shortest edge..\n", tri_id);

                    /* get shortest edge e = {u, v} and the other face T' incident to the shortest edge {u, v} */
                    tri_it->getTriShortestEdge(u_it, v_it);
                    uv_other_face_it = M.getOtherFaceIncidentToManifoldEdge(u_it, v_it, tri_it);

                    /* collapse if topologically safe. if not, push triangle onto unsafe_tris list and refetch from
                     * queue. if queue runs empty, we're done. no safe collapse left. if not, update all unsafe_tris and
                     * reinsert them in the queue... */
                    collapse_safe = M.collapseTriEdge(u_it, v_it, &w_it, NULL);

                    /* collapse of edge {u, v} was topologically safe and has been performed */
                    if (collapse_safe) {

                        /* update all faces incident to the new vertex w, since their AR has changed in general */
                        debugl(2, "updating faces incident to new vertex %5d\n", w_it->id() );

                        /* get face star of w and iterate over all incident faces */
                        w_it->getFaceStarIterators(w_fstar);
                        debugTabInc();
                        for (auto w_inc_tri : w_fstar) {
                            if (!w_inc_tri->isTri()) {
                                throw ("MeshAlg::greedyEdgeCollapsePostProcessing(): discovered non-triangle face. only triangular meshes are supported.");
                            }
                            debugl(2, "updating neighbour face %5d..\n", w_inc_tri->id() );

                            /* erase w_inc_tri->id() from Q if it is currently present */
                            debugTabInc();
                            if ( Q.changeKey(w_inc_tri->id(), -inf<R>()) == true) {
                                /* changeKey() has returned true => R (ar) key of value w_inc_tri->id() has been set to
                                 * -INF and is therefore guaranteed to be the minimum element. */
                                debugl(2, "neighbour triangle found in Q, key decreased to -INF. getting Q.top() and comparing ids..\n");
                                q_min = Q.top();
                                if (q_min.second != w_inc_tri->id()) {
                                    throw MeshEx(MESH_LOGIC_ERROR, "MeshAlg::greedyEdgeCollapsePostProcessing(): after performing changeKey(value = id, new_key = -INF) on priority queue Q, element with value id is not top() element. internal logic error.");
                                }
                                else Q.deleteMin();
                            }
                            else {
                                debugl(2, "neighbour triangle NOT found in Q. checking if it needs reinsertion..\n");
                            }
                            debugTabDec();

                            /* reinsert w_inc_tri if necessary. compute new values and evaluate processing predicate */
                            w_inc_tri_id                    = w_inc_tri->id();
                            w_inc_tri_ar                    = w_inc_tri->getTriAspectRatio();
                            w_inc_tri_area                  = w_inc_tri->getTriArea();
                            w_inc_tri_avg_surrounding_area  = avg_surrounding_area[w_inc_tri_id];

                            if ( proc(w_inc_tri_ar, w_inc_tri_area, w_inc_tri_avg_surrounding_area, alpha, lambda, mu) ) {
                                debugl(2, "ar still too large, area ok => reinserting..\n");
                                Q.insert( { -w_inc_tri_ar, w_inc_tri->id() } );
                            }
                            else {
                                debugl(2, "ar / area ok => not reinserting\n");
                            }
                        }
                        debugTabDec();

                        debugl(2, "reinserting all \"unsafe\" triangles if still existent and necessary..\n");

                        /* reinsert all "unsafe" triangles, should they still exist and require processing */
                        debugTabInc();
                        while ( !unsafe_tris.empty() ) {
                            /* get and remove id of first remaining erstwhile unsafe triangle */
                            us_tri_id = unsafe_tris.front();
                            unsafe_tris.pop_front();

                            /* erase unsafe triangle if its part of the queue */
                            if ( Q.changeKey(us_tri_id, -inf<R>()) == true) {
                                /* changeKey() has returned true => R (ar) key of value w_inc_tri->id() has been set to
                                 * -INF and is therefore guaranteed to be the minimum element. */
                                debugl(3, "unsafe triangle found in Q, key decreased to -INF. getting Q.top() and comparing id to us_tri_id\n");
                                q_min = Q.top();
                                if (q_min.second != us_tri_id) {
                                    throw MeshEx(MESH_LOGIC_ERROR, "MeshAlg::greedyEdgeCollapsePostProcessing(): after performing changeKey(value = id, new_key = -INF) on priority queue Q, element with value id is not top() element. internal logic error.");
                                }
                                else Q.deleteMin();
                            }
                            else {
                                debugl(3, "unsafe triangle NOT found in Q. checking if it needs reinsertion..\n");
                            }

                            /* reinsert if necessary */
                            us_tri_it = M.faces.find(us_tri_id);
                            if ( us_tri_it != M.faces.end() ) {
                                /* unsafe triangle still exists.. compute necessary values to decide whether it needs to
                                 * be reinserted */
                                us_tri_ar                   = us_tri_it->getTriAspectRatio();
                                us_tri_area                 = us_tri_it->getTriArea();
                                us_tri_avg_surrounding_area = avg_surrounding_area[us_tri_it->id()];

                                /* evaluate processing predicate */
                                if ( proc(us_tri_ar, us_tri_area, us_tri_avg_surrounding_area, alpha, lambda, mu) ) {
                                //if ( us_tri_area < mu * us_tri_avg_surrounding_area && (us_tri_ar >= alpha || us_tri_area < lambda * us_tri_avg_surrounding_area) )
                                    debugl(3, "us_tri_ar: %5.4f, us_tri_area: %10.5f, us_tri_avg_surrounding_area: %10.5f\n", us_tri_ar, us_tri_area, us_tri_avg_surrounding_area);
                                    Q.insert( { -us_tri_ar, us_tri_id } );
                                }
                            }
                        }
                        debugTabDec();
                    }
                    /* collapse of edge {u, v} was topologically unsafe and has not been performed.  append tri_id to
                     * list of unsafe triangle ids */
                    else {
                        debugl(2, "collapse unsafe, pushing triangle %5d to unsafe tri list..\n", tri_id);
                        unsafe_tris.push_back(tri_id);
                    }
                }
                /* else: tri_id refers to still existent, but now permissible triangle that does not require further
                 * processing, i.e. proc(..) == false */

            }
            /* else: tri_id refers to no longer existent face in mesh => Q element is "lazily" discarded */

        }
        debugTabDec();

        unsafe_tris.sort();
        unsafe_tris_prev_iter.sort();

        /* check for fixed point, which is the case iff unsafe_tris == unsafe_tris_prev_iter. */
        if (unsafe_tris == unsafe_tris_prev_iter) {
            fixed_point = true;
        }
        /* no fixed point yet, reinsert existent unsafe_tris again and perform another loop
         * iteration. */
        else {
            fixed_point = false;
            for (uint32_t tri_id : unsafe_tris) {
                if ( (tri_it = M.faces.find(tri_id)) != M.faces.end()) {
                    tri_ar              = tri_it->getTriAspectRatio();
                    Q.insert( { -tri_ar, tri_it->id() } );
                }
            }
        }
        iter++;
    }
    debugl(2, "fixed point reached. returning...\n");

    debugTabDec();
    debugl(0, "MeshAlg::greedyEdgeCollapsePostProcessing(): done.\n");
}


template <typename Tm, typename Tv, typename Tf, typename R>
inline R
GEC_getAvgAreaOfPermissibleSurroundingTriangles(
    typename Mesh<Tm, Tv, Tf, R>::face_iterator     tri_it,
    R const                                        &max_ar,
    uint32_t                                        depth)
{
    debugl(3, "GEC_getAvgAreaOfPermissibleSurroundingTriangles()\n");
    debugTabInc();

    R           nbtri_ar, nbtri_area, avg_area;
    uint32_t    npermissible_triangles;

    /* get face-neighbourhood of depth "depth" for triangle tri_it */
    std::vector<typename Mesh<Tm, Tv, Tf, R>::Face *> surrounding_tris;
    tri_it->getFaceNeighbourhood(depth, surrounding_tris);

    /* get average area of all permissible triangles in the depth-neighbourhood of the face
     * tri_it */
    avg_area                = 0.0;
    npermissible_triangles  = 0;
    for (auto &nbtri : surrounding_tris) {
        nbtri_ar    = nbtri->getTriAspectRatio();
        nbtri_area  = nbtri->getTriArea();

        if (nbtri_ar < max_ar) {
            debugl(4, "got permissible triangle %5d with ar: %10.5f and area: %10.5f\n", nbtri->id(), nbtri_ar, nbtri_area);
            avg_area += nbtri_area;
            npermissible_triangles++;
        }
    }

    if (npermissible_triangles == 0) {
        debugl(0, "MeshAlg::GEC_getAvgAreaOfPermissibleSurroundingTriangles(): WARNING: triangle %d: can't compute average, since no %5.4f-permissible triangle found in the %d-neighbour of %d. returning area %5.4f as \"average\"\n.",
            tri_it->id(), max_ar, depth, tri_it->id(), tri_it->getTriArea());
        return (tri_it->getTriArea());
    }
    else {
        debugTabDec();
        debugl(3, "GEC_getAvgAreaOfPermissibleSurroundingTriangles(): done.\n");
        return ( avg_area / (R)npermissible_triangles );
    }
}

/* ---------------------------------------------------------------------------------------------- */
/*                                                                                                */
/*                                post-processing: mesh smoothing algorithms                      */
/*                                                                                                */
/* ---------------------------------------------------------------------------------------------- */
template <typename Tm, typename Tv, typename Tf, typename R>
void
MeshAlg::simpleLaplacianSmoothing(
    Mesh<Tm, Tv, Tf, R>    &M,
    R const                &lambda,
    uint32_t                maxiter)
{
    uint32_t                                            iter, i, m;

    /* build temporary vector of vertices to calculate offsets through discrete laplacian and
     * neighbouring information */
    std::map<uint32_t, Vec3<R>>                         offsets;
    typename std::map<uint32_t, Vec3<R>>::iterator      oit;

    /* calculate offsets for each vertex. we use a simple discrete version of the laplacian, the
     * umbrella operator. dxi = 1/num_neighbours SUM_{all neighbours x_}{x_i - x_j}, so, the
     * centroid the x_i's neighbours. */
    std::list<typename Mesh<Tm, Tv, Tf, R>::Vertex *>   vi_nbs;
    typename Mesh<Tm, Tv, Tf, R>::vertex_iterator       v_it;
    Vec3<R>                                             xi, dxi, offset;

    for (iter = 0; iter < maxiter; iter++) {
        for (auto &vi : M.vertices) {
            /* get xi, reset dxi to nullvec */
            i   = vi.id();
            xi  = vi.pos();
            dxi = Aux::VecMat::nullvec<R>;

            /* get vi's vertex star and their size */
            vi.getVertexStar(vi_nbs);
            m = vi_nbs.size();
             
            /* isolated vertex .. maybe throw exception? */
            if (m > 0) {
                /* use umbrella operator to calculate offset dxi */
                for (auto &w : vi_nbs) {
                    dxi += (w->pos() - xi);
                }
                dxi *= lambda / ( (R) m);

                /* store offset into map */
                offsets[i] = dxi;
            }
            else offsets[i] = Aux::VecMat::nullvec<R>;
        }

        /* displace all vertices with corresponding offsets */
        for (v_it = M.vertices.begin(), oit = offsets.begin(); v_it != M.vertices.end(); ++v_it, ++oit) {
            v_it->pos() += oit->second;
        }
    }
}

template <typename Tm, typename Tv, typename Tf, typename R>
void
MeshAlg::HCLaplacianSmoothing(
    Mesh<Tm, Tv, Tf, R>    &M,
    R const                &alpha,
    R const                &beta,
    uint32_t                maxiter)
{
    uint32_t                                iter, i, m;

    /* for vertex displacement calculation, a simple discrete version of the laplacian is used, the umbrella operator.
     * dxi = 1/num_neighbours SUM_{all neighbours x_}{x_i - x_j}, so, the centroid the x_i's neighbours. */
    std::list<typename Mesh<Tm, Tv, Tf, R>::Vertex *>   vi_nbs;
    std::list<uint32_t>                                 vi_nbs_ids;
    Vec3<R>                                             xi_new;

    /* original vertex coordinates */
    std::map<uint32_t, Vec3<R>> o;

    /* vertex coordinates before the step */
    std::map<uint32_t, Vec3<R>> q;

    /* vertex coordinates after the step */
    std::map<uint32_t, Vec3<R>> p;

    /* correction offsets, pushing back the vertices to a weighted sum of original and previous
     * position to avoid volume shrinkage */
    std::map<uint32_t, Vec3<R>> b;

    /* initialize original coordinates o and current coordinates q: prior to the first step,
     * q = o */
    for (auto &v : M.vertices) {
        i       = v.id();
        o[i]    = v.pos();
        q[i]    = v.pos();
    }

    for (iter = 0; iter < maxiter; iter++) {
        for (auto &vi : M.vertices) {
            /* get vertex id, vertex star and vertex star size */
            i = vi.id();
            vi.getVertexStar(vi_nbs);
            m = vi_nbs.size();
             
            /* throw exception here? isolated vertex .. */
            if (m > 0) {
                /* compute new position with "umbrella" discrete laplacian operator, store in
                 * temporary xi_new. */
                xi_new  = Aux::VecMat::nullvec<R>();
                for (auto &w : vi_nbs) {
                    xi_new += w->pos();
                }
                xi_new *= 1.0 / ( (R) m);

                /* copy computed value from temporary xi_new into map p: only one access */
                p[i]    = xi_new;

                /* compute offset value b[i] for vertex i */
                b[i]    = xi_new - ( o[i]*alpha + q[i]*(1.0 - alpha) );
            }
        }

        /* iterate over all vertices again, correct new position p[i] with offsets b[j] of the
         * neighbours AND the offset value for the currently processed "centre" vertex vi
         * itself.. */
        Vec3<R> p_i_correction;
        R       nbfactor;
        for (auto &vi : M.vertices) {
            i               = vi.id();
            vi.getVertexStarIndices(vi_nbs_ids);
            m               = vi_nbs_ids.size();

            /* throw exception here? isolated vertex .. */
            if (m > 0) {
                p_i_correction  = b[i] * beta;
                nbfactor        = (1.0 - beta) / (R)m;

                for (uint32_t w_id : vi_nbs_ids) {
                    p_i_correction += (b[w_id] * nbfactor);
                }

                /* set new vertex position (note: this invokes operator=() of typename Mesh<Tm, Tv, Tf>::Vertex, after dereferencing the
                 * typename Mesh<Tm, Tv, Tf>::vertex_iterator has produces a reference to typename Mesh<Tm, Tv, Tf>::Vertex as lvalue). */
                vi.pos() = p[i] - p_i_correction;
                if (!Aux::VecMat::isfiniteVec3(vi.pos())) {
                    printf("vertex: %5d. new coordinate not finite: p[i] = (%10.5f, %10.5f, %10.5f), p_i_correction: (%10.5f, %10.5f, %10.5f) \n",
                            i,
                            p[i][0], p[i][1], p[i][2],
                            p_i_correction[0], p_i_correction[1], p_i_correction[2]);
                    throw("MeshAlg::HCLaplacianSmoothing(): new vertex position vector computed during iteration has infinite component.");
                }
            }
        }

        /* shift: q = p. clear p. */
        q.clear();
        q.swap(p);
        /* q now contains contents of p before the swap(), p is empty. done for this iteration. */
    }
}

#endif

template <typename Tm, typename Tv, typename Tf, typename R>
void
MeshAlg::partialFlushToObjFile(
    Mesh<Tm, Tv, Tf, R>                                        &M,
    std::pair<FILE **, std::string> const                      &obj_file_info,
    std::list<typename Mesh<Tm, Tv, Tf, R>::Face *>            &face_list,
    std::list<
            std::pair<
                typename Mesh<Tm, Tv, Tf, R>::Vertex *,
                uint32_t
            >
        >                                                      &in_boundary_vertices,
    uint32_t const                                             &in_last_flush_vertex_id,
    std::list<
            std::pair<
                typename Mesh<Tm, Tv, Tf, R>::Vertex *,
                uint32_t          
            >
        >                                                      &out_boundary_vertices,
    uint32_t                                                   &out_last_flush_vertex_id)
{
    debugl(0, "MeshAlg::partialFlush().\n");
    debugTabInc();
    /* copy out information about all faces in face_list, the list of faces to be flushed, and subsequently delete them
     * in M. */

    /* struct to hold information about face */
    struct IdFace {
        bool                    quad;
        std::vector<uint32_t>   v_ids;

        IdFace(
            bool                            quad,
            std::vector<uint32_t> const    &v_ids)
        {
            this->quad  = quad;
            this->v_ids = v_ids;

            if (    (this->quad && v_ids.size() != 4) ||
                    (!this->quad && v_ids.size() != 3) )
            {
                debugl(0, "ERROR: quad: %d. v_ids.size(): %d\n", quad, v_ids.size());
                throw("MeshAlg::partialFlush(..)::IdFace::IdFace(): given vertex index vector has wrong size (neither tri / 3 nor quad / 4). internal logic error.");
            }
        }
    };

    debugl(0, "in_boundary_vertices.size(): %zu. in_last_flush_vertex_id: %d, face_list.size(): %zu\n",
        in_boundary_vertices.size(), in_last_flush_vertex_id, face_list.size());

#ifdef __DEBUG__
    debugTabInc();
    for (auto &vp : in_boundary_vertices) {
        debugl(1, "in_boundary_vertices element: mesh id: %d, flush id: %d\n", vp.first->id(), vp.second);
    }
    debugTabDec();
#endif

    std::list<IdFace> flush_face_list;
    for (auto &f : face_list) {
        if (f->isQuad()) {
            flush_face_list.push_back({ true, f->getIndices() });
        }
        else if (f->isTri()) {
            flush_face_list.push_back({ false, f->getIndices() });
        }
        else {
            throw("MeshAlg::partialFlush(): discovered face that is neither quad nor triangle. flushing not (yet) supported.");
        }

        M.faces.erase(f->iterator());
    }

    /* the deletion of all faces in face_list will generally create new isolated vertices and new non-isolated, but
     * non-manifold vertices (a little informally called boundary vertices). note however that not all boundary vertices
     * must be new: for example, a vertex in in_boundary_vertices might be incident to a face f from face_list and still
     * not become isolated. care must be taken to distinguish old from new boundary vertices and similarly, old and new
     * isolated vertices, where "old" isolated vertices are vertices that have been part of in_boundary_vertices at the
     * beginning of the call and have become isolated during face deletion. */

    /* firstly, find all isolated vertices and boundary vertices (i.e. non-isolated and non-manifold). set ids in pairs
     * to zero, they will be computed below. */
    std::list<
            std::pair<
                typename Mesh<Tm, Tv, Tf, R>::Vertex *,
                uint32_t
            >
        >                                                   isolated_vertices, new_isolated_vertices, boundary_vertices, new_boundary_vertices;

    debugTabInc();
    for (auto &v : M.vertices) {
        std::list<typename Mesh<Tm, Tv, Tf, R>::Face *> v_fstar;
        v.getFaceStar(v_fstar);
        debugl(1, "vertex %d. face star size: %d .. ", v.id(), v_fstar.size());

        if (v.isIsolated()) {
            isolated_vertices.push_back({&v, 0});
            debugl(1, "is isolated.\n");
        }
        else if (!v.isManifold()) {
            boundary_vertices.push_back({ &v, 0});
            debugl(1, "is not isolated but non-manifold => boundary vertex.\n");
        }
        else {
            debugl(1, "is regular.\n");
        }
    }
    debugTabDec();

    debugl(0, "isolated_vertices.size(): %zu, boundary_vertices.size(): %zu.\n", isolated_vertices.size(), boundary_vertices.size());

    /* sort all lists, unique input boundary vertex list for safety */
    auto cmp =
        [] (
            std::pair<
                typename Mesh<Tm, Tv, Tf, R>::Vertex *,
                uint32_t
            > const &x,
            std::pair<
                typename Mesh<Tm, Tv, Tf, R>::Vertex *,
                uint32_t
            > const &y)
        -> bool
        {
            return (Mesh<Tm, Tv, Tf, R>::Vertex::ptr_less(x.first, y.first)); // || (x.first == y.first && x.second < y.second) );
        };

    in_boundary_vertices.sort(cmp);
    in_boundary_vertices.unique();
    isolated_vertices.sort(cmp);
    boundary_vertices.sort(cmp);

    debugl(0, "after sort: isolated_vertices.size(): %zu, boundary_vertices.size(): %zu.\n", isolated_vertices.size(), boundary_vertices.size());

    /* compute sef-differences: isolated_vertices \ in_boundary_vertices and boundary_vertices \ in_boundary_vertices */ 
    std::set_difference(
        isolated_vertices.begin(), isolated_vertices.end(),
        in_boundary_vertices.begin(), in_boundary_vertices.end(),
        std::back_inserter(new_isolated_vertices), cmp);

    std::set_difference(
        boundary_vertices.begin(), boundary_vertices.end(),
        in_boundary_vertices.begin(), in_boundary_vertices.end(),
        std::back_inserter(new_boundary_vertices), cmp);


    debugl(0, "after sort: new_isolated_vertices.size(): %zu, new_boundary_vertices.size(): %zu.\n",
            new_isolated_vertices.size(), new_boundary_vertices.size());

    /* consecutively number all new isolated vertices and new boundary vertices. compile id replacement map and update
     * flush_face_list ids */
    uint32_t                        last_flush_vertex_id = in_last_flush_vertex_id;

    /* initialize id replacement map. associate mesh ids of all in_boundary_vertices with flush ids. number all new
     * vertices consecutively and associate new ids in the process */
    debugl(0, "initializing id replacement map.\n");

    std::map<uint32_t, uint32_t> id_replace_map;
    debugTabInc();
    for (auto vp : in_boundary_vertices) {
        debugl(1, "old in boundary vertex: %5d -> %d.\n", vp.first->id(), vp.second);
        id_replace_map[vp.first->id()] = vp.second;
    }
    for (auto it = new_isolated_vertices.begin(); it != new_isolated_vertices.end(); ++it) {
        debugl(1, "new isolated vertex: %d -> assigning id %d.\n", it->first->id(), last_flush_vertex_id);
        it->second = last_flush_vertex_id++;
        id_replace_map[it->first->id()] = it->second;
    }
    for (auto it = new_boundary_vertices.begin(); it != new_boundary_vertices.end(); ++it) {
        debugl(1, "new boundary vertex: %d -> assigning id %d.\n", it->first->id(), last_flush_vertex_id);
        it->second = last_flush_vertex_id++;
        id_replace_map[it->first->id()] = it->second;
    }
    debugTabDec();

    /* write out_last_vertex_id back to the caller. */
    out_last_flush_vertex_id = last_flush_vertex_id;

    /* replace all indices in flush_face_list using the generated id map above. */
    debugl(0, "replacing vertex indices in all flush faces with flush indices..\n"); 
    debugTabInc();
    for (auto &f : flush_face_list) {
        debugl(1, "replacing ids in face (%d, %d, %d)\n", f.v_ids[0], f.v_ids[1], f.v_ids[2]);
        for (auto &id : f.v_ids) {
            auto it = id_replace_map.find(id);
            if (it != id_replace_map.end()) {
                id = it->second;
            }
            else {
                throw("MeshAlgorithms::partialFlush(): failed to locate vertex index in id_replace_map. internal logic error.");
            }
        }
        debugl(1, "replaced ids: (%d, %d, %d)\n", f.v_ids[0], f.v_ids[1], f.v_ids[2]);
    }
    debugTabDec();

    /* write new isolated vertices (with correct id) and all NEW boundary vertices to obj file. as part of the
     * invariant, all old boundary vertices had already been written to the obj file when the call started. */

    /* since vertices should be the first block and the face definition block (which use vertex indices) should be
     * below, it is necessary to insert new vertex definition lines in the middle of obj_file, an operation that is
     * generally unsupported by most file systems. instead, "merge" the file and the new information into a temporary
     * file, rename to correct filename and adjust the file descriptor. 
     *
     * first, scan obj_file for the vertex block delimiter "# ____~V____". as long as it is not found, copy lines to
     * swap file. after delimiter has been found, insert new vertex lines into swap file, followed by the rest of
     * obj_file and finally the new faces. */
    debugl(0, "writing new partial mesh to swap file..\n");

    /* get "original" obj file */
    FILE               *obj_file    = *(obj_file_info.first);

    /* open swap file */
    std::string const &filename     = obj_file_info.second;
    std::string swap_filename       = filename + "_swap";
    FILE *swap_file                 = fopen( (swap_filename + ".obj").c_str(), "w");
    const char v_delim[]            = "# ____~V____";
    char line[1024];

    if (!swap_file) {
        throw("MeshAlg::partialFlush(): can't open swap file for writing.");
    }

    /* if orig file is empty, write new vertices / faces directly */
    if (Aux::File::isEmpty(obj_file)) {
        debugl(0, "given obj file empty..\n");

        fprintf(swap_file, "# %5zu flushed vertices\n", new_isolated_vertices.size() + new_boundary_vertices.size());
        Vec3<R> vpos;
        for (auto &vp : new_isolated_vertices) {
            vpos = vp.first->pos();
            fprintf(swap_file, "v %+.10e %+.10e %+.10e\n", vpos[0], vpos[1], vpos[2]);
        }
        for (auto &vp : new_boundary_vertices) {
            vpos = vp.first->pos();
            fprintf(swap_file, "v %+.10e %+.10e %+.10e\n", vpos[0], vpos[1], vpos[2]);
        }

        /* and write delimiter again */
        fprintf(swap_file, "%s\n", v_delim);

        debugl(0, "done writing new vertices and delimiter.\n");
    }
    /* otherwise assemble swap_file from new information and obj_file */
    else {
        debugl(0, "given obj file non-empty.. \"merging\" together old and new information into swap_file..\n");
        rewind(obj_file);
        debugTabInc();
        while ( fgets(line, sizeof(line), obj_file) == line) {
            if (line[strlen(line) - 1] == '\n') {
                /* overwrite newline with zero-termination */
                line[strlen(line) - 1] = '\0';
                debugl(1, "line from original file: \"%s\".\n", line);
                if (strncmp(line, v_delim, sizeof(v_delim)) == 0) {
                    debugl(2, "writing new vertices and delimiter.\n");

                    /* delimiter found. write new vertices */
                    fprintf(swap_file, "# %5zu flushed vertices\n", new_isolated_vertices.size() + new_boundary_vertices.size());
                    Vec3<R> vpos;
                    for (auto &vp : new_isolated_vertices) {
                        vpos = vp.first->pos();
                        fprintf(swap_file, "v %+.10e %+.10e %+.10e\n", vpos[0], vpos[1], vpos[2]);
                    }
                    for (auto &vp : new_boundary_vertices) {
                        vpos = vp.first->pos();
                        fprintf(swap_file, "v %+.10e %+.10e %+.10e\n", vpos[0], vpos[1], vpos[2]);
                    }

                    /* and write delimiter again */
                    fprintf(swap_file, "%s\n", v_delim);

                    debugl(2, "done writing new vertices and delimiter.\n");
                }
                /* copy line to swap_file */
                else {
                    fprintf(swap_file, "%s\n", line);
                }
            }
            else {
                throw("MeshAlg::partialFlush(): obj line longer than 1024 characters. please limit comments to reasonable line width.");
            }
        }
        debugTabDec();
    }

    /* append all faces to swap file */
    fprintf(swap_file, "# %5zu flushed faces.\n", flush_face_list.size());
    for (auto &f : flush_face_list) {
        if (f.quad) {
            fprintf(swap_file, "f %d//%d %d//%d %d//%d %d//%d\n",
                    f.v_ids[0] + 1, 
                    f.v_ids[0] + 1, 
                    f.v_ids[1] + 1,
                    f.v_ids[1] + 1,
                    f.v_ids[2] + 1,
                    f.v_ids[2] + 1,
                    f.v_ids[3] + 1,
                    f.v_ids[3] + 1);
        }
        else {
            fprintf(swap_file, "f %d//%d %d//%d %d//%d\n",
                    f.v_ids[0] + 1, 
                    f.v_ids[0] + 1, 
                    f.v_ids[1] + 1,
                    f.v_ids[1] + 1,
                    f.v_ids[2] + 1,
                    f.v_ids[2] + 1);
        }
    }

    debugl(1, "flushing / synching / closing obj_file\n");

    /* flush internal buffers, kernel buffers, close original file */
    fflush(obj_file);
    fsync(fileno(obj_file));
    fclose(obj_file);

    fflush(swap_file);
    fsync(fileno(swap_file));
    fclose(swap_file);

    debugl(1, "removing original file, rename swap file to original file.\n");

    /* remove original file, rename swap_file to original file, adjust FILE * reference in file_info */
    if ( remove( (filename + ".obj").c_str() ) != 0) {
        throw("MeshAlg::partialFlush(): can't remove old obj file before overwriting with swap file.");
    }
    if (rename( (swap_filename + ".obj").c_str(), (filename + ".obj").c_str() ) != 0) {
        throw("MeshAlg::partialFlush(): can't rename swap file to filename of obj file.");
    }

    debugl(1, "assigning new FILE * to obj_file_info pair by reopening filename, which now contains moved swap file...\n");
    /* update file pointer in obj_file_info: reopen new obj file (moved swap file) in append mode and rewind() */
    FILE *tmp = fopen( (filename + ".obj").c_str(), "r+");
    if (!tmp) {
        throw("MeshAlg::partialFlush(): can't re-open obj file after having removed and overwritten old one with swap file.");
    }
    else {
        rewind(tmp);
        *(obj_file_info).first = tmp;
    }
    
    debugl(1, "finishing invariants ..\n");

    /* write out_boundary_vertices for the caller: out_boundary_vertices is the union of new_boundary_vertices and all
     * old boundary vertices that have not become isolated. since a boundary vertex can either stay a boundary vertex or
     * become isolated (losing the boundary status) and no vertices have yet been deleted (and thus all pointers are
     * still intact), it's possible to simply iterate over in_boundary_vertices and extract all vertices that have not
     * become isolated.  note that the correct ids are copied as well: these are contained in in_boundary_vertices from
     * the beginning of the call. */

    /* copy in_boundary_vertices so as to enable the caller to use the same list for both references
     * in_boundary_vertices and out_boundary_vertices */
    auto in_boundary_vertices_copy = in_boundary_vertices;

    /* initialize (and overwrite) out_boundary_vertices with new_boundary_vertices. this might also overwrite
     * in_boundary_vertices if both references are identical. */
    out_boundary_vertices = new_boundary_vertices;

    /* handle input boundary vertices via copy. */
    for (auto &ibv : in_boundary_vertices_copy) {
        if (!ibv.first->isIsolated()) {
            out_boundary_vertices.push_back(ibv);
        }
    }

    /* sort */
    out_boundary_vertices.sort(cmp);
        
    /* delete all isolated vertices, old and new, from M. note that only old boundary vertices, which have become
     * isolated during the call, are thereby deleted. no other boundary vertex is deleted, but they have been written to
     * the obj file already to guarantee the invariant for the next flushing or the finalizing call. */
    debugl(1, "deleting all new isolated vertices.\n");
    debugTabInc();
    for (auto &vp : isolated_vertices) {
        debugl(2, "deleting isolated vertex %d.\n", vp.first->id());
        M.vertices.erase(vp.first->iterator());
    }
    debugTabDec();

    debugTabDec();
    debugl(0, "MeshAlg::partialFlush(): done.\n");
}

template <typename Tm, typename Tv, typename Tf, typename R>
void
MeshAlg::partialFlushToObjFile(
    Mesh<Tm, Tv, Tf, R>                                        &M,
    MeshObjFlushInfo<Tm, Tv, Tf, R>                            &M_flush_info,
    std::list<typename Mesh<Tm, Tv, Tf, R>::Face *>            &face_list)
{
    /* check if info has been prepared */
    if (!M_flush_info.obj_file) {
        throw("MeshAlg::partialFlushToObjFile(): given obj flush info struct not properly initialized. file handle is NULL.");
    }

    /* call above overload version with full parameter list */
    MeshAlg::partialFlushToObjFile<Tm, Tv, Tf, R>(
        M,
        { &M_flush_info.obj_file, M_flush_info.filename },
        face_list,
        M_flush_info.last_boundary_vertices,
        M_flush_info.last_flush_vertex_id,
        M_flush_info.last_boundary_vertices,
        M_flush_info.last_flush_vertex_id);
}
