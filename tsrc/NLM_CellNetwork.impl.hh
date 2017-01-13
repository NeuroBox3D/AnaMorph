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
 *                        NeuritePath implementation.
 *
 * ----------------------------------------------------------------------------------------------------------------- */
#include "PolyAlgorithms.hh"

/*
template<typename R>
NeuritePath<R>::NeuritePath()
{
    this->C                 = NULL;
    this->root_path         = false;;
    this->geometry_updated  = false;
    this->bb_set            = false;
    this->bb_min            = Vec3(0, 0, 0);
    this->bb_max            = Vec3(0, 0, 0);
}
*/

namespace NLM {
    template<typename R>
    NeuritePath<R>::NeuritePath(
        typename NLM_CellNetwork<R>::neurite_segment_iterator   start_ns_it,
        bool                                                    root_path)
            : neurite_segments { start_ns_it }, neurite_vertex_parameters { 0 }
    {
        NLM_CellNetwork<R> *C_downcast  = dynamic_cast<NLM_CellNetwork<R> *>(start_ns_it.getContainer());
        if (C_downcast) {
            this->C                     = C_downcast;
        }
        else {
            throw("NeuritePath::NeuritePath(): given start neurite segment does not belong to an NLM_CellNetwork as required.");
        }
        /*
        this->neurite_segments          = { start_ns_it };
        this->neurite_vertex_parameters = { 0 };
        */
        this->root_path                 = root_path;
        this->geometry_updated          = false;
        this->bb_set                    = false;
        this->bb                        = BoundingBox<R>({
                                                Aux::VecMat::onesVec3<R>() * -Aux::Numbers::inf<R>(),
                                                Aux::VecMat::onesVec3<R>() * +Aux::Numbers::inf<R>()
                                            });

        /* if root_path == true, check if the source vertex of ns_it is indeed
         * adjacent to exactly one soma */
        if (root_path) {
            auto tmp = start_ns_it->getSourceVertex()->
                template getFilteredInEdges<typename NLM_CellNetwork<R>::NeuriteRootEdge>();

            if (tmp.size() != 1) {
                throw("NeuritePath::NeuritePath(): path is supposed to be root path, yet supplied initial neurite segment is not incident to exactly one neurite root edge.");
            }
        }
    }

    template<typename R>
    NeuritePath<R>::NeuritePath(NeuritePath const &p)
    {
        this->C                         = p.C;
        this->neurite_vertex_parameters = p.neurite_vertex_parameters;
        this->neurite_segments          = p.neurite_segments;
        this->canal_segments_magnified  = p.canal_segments_magnified;
        this->root_path                 = p.root_path;
        this->geometry_updated          = p.geometry_updated;
        this->bb_set                    = p.bb_set;
        this->bb                        = p.bb;
    }

    template<typename R>
    NeuritePath<R> &
    NeuritePath<R>::operator=(NeuritePath const &p)
    {
        this->C                         = p.C;
        this->neurite_vertex_parameters = p.neurite_vertex_parameters;
        this->neurite_segments          = p.neurite_segments;
        this->canal_segments_magnified  = p.canal_segments_magnified;
        this->root_path                 = p.root_path;
        this->geometry_updated          = p.geometry_updated;
        this->bb_set                    = p.bb_set;
        this->bb                        = p.bb;

        return (*this);
    }

    template<typename R>
    NeuritePath<R>::~NeuritePath()
    {
    }

    template<typename R>
    uint32_t
    NeuritePath<R>::numVertices() const
    {
        if (!this->neurite_segments.empty()) {
            return (this->neurite_segments.size() + 1);
        }
        else {
            return 0;
        }
    }

    template<typename R>
    uint32_t
    NeuritePath<R>::numEdges() const
    {
        return (this->neurite_segments.size());
    }

    template<typename R>
    uint32_t
    NeuritePath<R>::extend(typename NLM_CellNetwork<R>::neurite_segment_iterator ns_it)
    {
        /* insert dummy parameter zero */
        if (ns_it.getContainer() == this->C) {
            /* check if source vertex of given neurite segment is destination vertex of the previous one (if any) */
            if (!this->neurite_segments.empty() && this->neurite_segments.back()->getDestinationVertex() != ns_it->getSourceVertex()) {
                throw("NeuritePath::extend(): can't extent path with neurite segment that is not incident to the previously inserted one.");
            }
            else {
                /* insert neurite segment */
                this->neurite_segments.push_back(ns_it);
                this->neurite_vertex_parameters.push_back(0);

                /* geometry needs update */
                this->geometry_updated = false;

                /* return id of segment */
                return (this->neurite_segments.size() - 1);
            }
        }
        else {
            throw("NeuritePath::extend(): cannot extend neurite path with iterator from different container.");
        }
    }

    template<typename R>
    typename NLM_CellNetwork<R>::neurite_iterator
    NeuritePath<R>::operator[](uint32_t idx) const
    {
        if (idx > this->neurite_segments.size()) {
            throw("NeuritePath::operator[]: given neurite vertex index out of range.");
        }
        /* get last vertex, which is the destination vertex of the back() neurite segment */
        else if (idx == this->neurite_segments.size()) {
            return (this->neurite_segments.back()->getDestinationVertex());
        }
        /* otherwise return source vertex of neurite segment with index idx */
        else {
            return (this->neurite_segments[idx]->getSourceVertex());
        }
    }

    template<typename R>
    void
    NeuritePath<R>::computeChordLengthParametrization()
    {
        debugl(3, "NeuritePath::computeChordLengthParametrization()\n");
        debugTabInc();

        uint32_t const  n = this->neurite_segments.size() + 1;
        uint32_t        i;

        /* init neurite_vertex_parameters to right size and zero */
        this->neurite_vertex_parameters.assign(n, 0);

        /* sum up neurite segment lentsh to get entire chord length and store partial chord lengths in
         * neurite_vertex_parameters in the process. */
        R cordlen = 0;
        for (i = 1; i < n; i++) {
            cordlen                            += this->neurite_segments[i-1]->getLength();
            this->neurite_vertex_parameters[i]  = cordlen;
        }

        /* cordlen now contains the entire chordlength. divide everything by cordlen to get relative chordlength in
         * [0,1] */
        this->neurite_vertex_parameters[0]      = 0;
        for (i = 1; i < n - 1; i++) {
            this->neurite_vertex_parameters[i] /= cordlen;
        }
        this->neurite_vertex_parameters[n-1]    = 1;

        debugTabDec();
        debugl(3, "NeuritePath::computeChordLengthParametrization(): done.\n");
    }

    template<typename R>
    void
    NeuritePath<R>::computeUniformParametrization()
    {
        debugl(3, "NeuritePath::computeUniformParametrization(): done.\n");
        debugTabInc();

        uint32_t i, n = this->neurite_segments.size() + 1;
        this->neurite_vertex_parameters.assign(n, 0.0);

        this->neurite_vertex_parameters.front() = 0.0;
        for (i = 1; i < n - 1; i++) {
            this->neurite_vertex_parameters[i] = (R)i/(R)n;
        }
        this->neurite_vertex_parameters.back() = 1.0;

        debugTabDec();
        debugl(3, "NeuritePath::computeUniformParametrization(): done.\n");
    }

    template<typename R>
    void
    NeuritePath<R>::computeCentripetalParametrization()
    {
        uint32_t const n = this->neurite_segments.size() + 1;
        this->neurite_vertex_parameters.assign(n, 0.0);
    }

    template<typename R>
    R
    NeuritePath<R>::getChordLength() const
    {
        R   clen = 0;
        for (auto &ns : this->neurite_segments) {
            clen += ns->getLength();
        }
        return clen;
    }

    template<typename R>
    void
    NeuritePath<R>::updateGeometry(
        std::function<
                std::tuple<
                    R,
                    std::vector<R>,
                    R
                >(NLM::NeuritePath<R> const &P)
            >  const                       &parametrization_algorithm)
    {
        using Aux::Numbers::inf;
        using Aux::VecMat::onesVec3;

        debugl(1, "NeuritePath::updateGeometry().\n");
        debugTabInc();

        debugl(2, "extracting vertex position components into vectors for linear solver..\n");
        /* set m to the number of segments on the path */
        uint32_t const m = this->neurite_segments.size();

        /* to compute all neurite canal segments, solve the spline system for all three space components {x,y,z}, the data
         * for which is extracted from all neurite vertices into the vectors {x,y,z}_components. */
        std::vector<R> vertex_position_components[3], result_powercoeff[3];
        for (uint32_t j = 0; j < 3; j++) {
            vertex_position_components[j].assign(m + 1, 0);
        }

        /* store coordinates of source vertex of first neurite segments */
        Vec3<R> npos = this->neurite_segments.front()->getSourceVertex()->getPosition();
        vertex_position_components[0][0] = npos[0];
        vertex_position_components[1][0] = npos[1];
        vertex_position_components[2][0] = npos[2];

        /* iterator over all neurite segments and handle destination vertices */
        for (uint32_t i = 1; i < m + 1; i++) {
            npos = this->neurite_segments[i - 1]->getDestinationVertex()->getPosition();
            vertex_position_components[0][i] = npos[0];
            vertex_position_components[1][i] = npos[1];
            vertex_position_components[2][i] = npos[2];
        }

        debugl(2, "computing hermite boundary condition vectors..\n");
        /* tangent boundary condition vectors for first and last vertex are obtained from the scaled unit direction vectors
         * between the two first and two last vertices: exception: if this is a neurite root path, set ds0 to be the normal
         * vector of the soma sphere */
        Vec3<R>                                                         ds0, dsn;
        typename NLM_CellNetwork<R>::neurite_iterator                   v0, v1;
        if (this->root_path) {
            v0          = this->neurite_segments[0]->getSourceVertex();
            auto slist  = v0->template getFilteredInNeighbours<typename NLM_CellNetwork<R>::SomaVertex>();
            if (slist.size() == 1) {
                ds0 = ( v0->getPosition() - slist.front()->getSinglePointPosition() );
            }
            else {
                debugl(1,"NeuritePath::updateGeometry(): root_path = %d, v0 = %d, number of soma in-neighbours: %d, number of in-neighbours: %d  => internal logic error.\n", this->root_path, v0->id(), slist.size(), v0->getInNeighbours().size());
                throw("NeuritePath::updateGeometry(): internal state indicates that (this) is a neurite root path, yet first vertex does not have exactly one soma as parent.");
            }
        }
        else {
            v0  = this->neurite_segments[0]->getSourceVertex();
            v1  = this->neurite_segments[0]->getDestinationVertex();
            ds0 = ( v1->getPosition() - v0->getPosition() );
        }

        dsn = ( this->neurite_segments.back()->getDestinationVertex()->getPosition() -
                this->neurite_segments.back()->getSourceVertex()->getPosition() );

        /* normalize tangent boundary condition vectors, they are scaled depending on the chosen parametrization below */
        ds0.normalize();
        dsn.normalize();

        debugl(2, "computing parametrization..\n");

        /* apply given parametrization algorithm, which returns a tuple containing
         *
         *  1. a scaling factor for ds0 as the first element
         *  2. a vector of parameters as the second element
         *  3. a scaling factor for dsn as the last element. */
        auto para_algo_tuple            = parametrization_algorithm((*this));
        ds0                            *= std::get<0>(para_algo_tuple);
        this->neurite_vertex_parameters = std::get<1>(para_algo_tuple);
        dsn                            *= std::get<2>(para_algo_tuple);

        debugl(1, "calling Thomson linear solver for spline system for three component splines ({x, y, z})..\n");
        /* get three result vectors for {x, y, z}-coefficient interpolation polynomial, each consisting of the 4n power
         * coefficients of the interpolating cubic spline with hermite boundary condition (ds0, dsn) for the respective
         * component {x,y,z}. for all (n-1) segments of the implicitly defined space curve, convert the component functions
         * to BernsteinPolynomials and init BezierCanalSurface of the respective neurite segment */
        for (uint32_t j = 0; j < 3; j++) {
            Aux::Numerics::solveSplineCoefficientSystem<R>(
                this->neurite_vertex_parameters,
                vertex_position_components[j],
                ds0[j],
                dsn[j],
                result_powercoeff[j]);
        }

        debugl(1, "systems solved. generating neurite segment curves and neurite canal segments..\n");
        /* for all n segments, assemble spine curve segment in power basis, convert to BezierCurve by converting
         * component functions to BernsteinPolynomials and finally create neurite canal segment BLRCanalSurface using the
         * radii of the two neurite vertices forming the neurite segment. */
        typename PowerPolynomial<3u, double, double>::coeff_type powercoeff_tmp(0);
        PowerPolynomial<3u, R, R>                                powerpoly_tmp;
        BernsteinPolynomial<3u, R, R>                            bbpoly_tmp;
        std::array<BernsteinPolynomial<3u, R, R>, 3>             segment_i_components;
        BezierCurve<3u, R>                                           segment_i_curve;
        Vec3<R>                                     segment_i_bb_min, segment_i_bb_max;

        /* reset bounding box of neurite path */
        this->bb = BoundingBox<R>({ onesVec3<R>() * inf<R>(), onesVec3<R>() * (-inf<R>()) });

        /* resize canal segments array to correct size m */
        this->canal_segments_magnified.resize(m);

        uint32_t coeff_index;
        debugTabInc();
        for (uint32_t i = 0; i < m; i++) {
            debugl(2, "segment %d.\n", i);
            debugTabInc();

            debugl(2, "extracting power coefficients for neurite segment curve and converting to BernsteinPolynomial.\n");
            coeff_index = 4*i;
            for (uint32_t j = 0; j < 3; j++) {
                /* read in result power coefficients for j-th component polynomial of segment i "backwards.." for historical
                 * reasons o_O */
                powercoeff_tmp[3] = result_powercoeff[j][coeff_index    ];
                powercoeff_tmp[2] = result_powercoeff[j][coeff_index + 1];
                powercoeff_tmp[1] = result_powercoeff[j][coeff_index + 2];
                powercoeff_tmp[0] = result_powercoeff[j][coeff_index + 3];

                /* create temporary power polynomial for component j of segment i */
                powerpoly_tmp = PowerPolynomial<3u, R, R>(powercoeff_tmp);

                /* convert to BernsteinPolynomial */
                PolyAlg::convertBasis(bbpoly_tmp, powerpoly_tmp);

                /* magnify bernstein polynomial to relevant domain [t_i, t_{i+1}] and store in segment_i_components
                 * array. */
                bbpoly_tmp.clipToInterval(
                        this->neurite_vertex_parameters[i],
                        this->neurite_vertex_parameters[i+1],
                        &(segment_i_components[j]));
            }

            debugl(2, "creating bezier curve from three component BernsteinPolynomials\n");
            /* create bezier spine curve for segment i */
            segment_i_curve = BezierCurve<3u, R>(segment_i_components);

            debugl(2, "creating BLRCanalSurface neurite canal segment..\n");
            /* finally, create BLRCanalSurface and store in NLM::NeuriteSegmentInfo of neurite segment i of (this) path. */
            NLM::NeuriteSegmentInfo<R> &segment_i_info  = this->neurite_segments[i]->neurite_segment_data;

            segment_i_info.canal_segment_magnified      =
                BLRCanalSurface<3u, R>(
                    segment_i_curve, 
                    this->neurite_segments[i]->getSourceVertex()->getRadius(),
                    this->neurite_segments[i]->getDestinationVertex()->getRadius()
                );

            /* store pointer to created magnified canal segment in NeuritePath::canal_segment_magnified. */
            this->canal_segments_magnified[i]           = &(segment_i_info.canal_segment_magnified);

            /* update bounding box of neurite canal segment and use it to calculate the bounding box of the entire neurite
             * path. */
            debugl(2, "updating canal segment bounding box and potentially bounding box of entire neurite path.\n");
            segment_i_info.canal_segment_magnified.updateBoundingBox();
            this->bb.update(segment_i_info.canal_segment_magnified.getBoundingBox());

            debugTabDec();
            debugl(2, "done for segment %d\n", i);
        }
        debugTabDec();

        /* if this is a root path, also compute the initial cylinder segment between first neurite vertex and the soma,
         * which is adjacent to the first neurite vertex. */
        if (this->root_path) {
            debugl(1, "(this) is a neurite ROOT path => computing initial cylinder segment between soma centre and neurite root vertex.\n");
            /* get NLM::NeuriteRootEdgeInfo */
            std::list<typename NLM_CellNetwork<R>::NeuriteRootEdge *> re_list;
            this->neurite_segments[0]->getSourceVertex()->
                template getFilteredInEdges<typename NLM_CellNetwork<R>::NeuriteRootEdge>(
                    re_list
                );

            if (re_list.size() == 1) {
                /* get NeuriteRootEdge between first neurite vertex of this path and soma */
                typename NLM_CellNetwork<R>::neurite_rootedge_iterator re_it = re_list.front()->iterator();

                /* get NLM::NeuriteRootEdgeInfo for re and iterators to source (soma) and destination (first neurite) vertex. */
                NLM::NeuriteRootEdgeInfo<R>                    &re_info = re_it->neurite_root_edge_data;
                typename NLM_CellNetwork<R>::neurite_iterator   v0_it   = re_it->getDestinationVertex();
                typename NLM_CellNetwork<R>::soma_iterator      s_it    = re_it->getSourceVertex();

                /* get position of soma and position of first vertex on the path. */
                Vec3<R> s_pos   = s_it->getSinglePointPosition();
                Vec3<R> v0_pos  = v0_it->getPosition();
                Vec3<R> d       = v0_pos - s_pos;
                R       r_v0    = v0_it->getRadius();

                re_info.initial_cylinder = BLRCanalSurface<3u, R>(
                    BezierCurve<3u, R>({ s_pos, s_pos + d / 3.0, s_pos + d * (2.0 / 3.0), v0_pos }),
                    r_v0,
                    r_v0);

                /* update bounding box of initial cylinder segment and use it to calculate the bounding box of the entire
                 * neurite path. */
                Vec3<R> cyl_bb_min, cyl_bb_max;
                re_info.initial_cylinder.updateBoundingBox();
                this->bb.update(re_info.initial_cylinder.getBoundingBox());
            }
            else {
                throw("NeuritePath::updateGeometry(): (this) path is supposed to be a root path, yet start vertex is not "\
                    "incident to exactly one neurite root edge.");
            }
        }

        /* extend bounding box of neurite path by 2,5%, but no less than 1E-3, in every component. */
        this->bb.extend(0.025, Vec3<R>(1E-3, 1E-3, 1E-3));
        this->bb_set            = true;

        /* geometry has been updated */
        this->geometry_updated  = true;

        debugTabDec();
        debugl(1, "NeuritePath::updateGeometry(): done.\n");
    }

    template<typename R>
    BoundingBox<R>
    NeuritePath<R>::getBoundingBox() const
    {
        if (this->bb_set) {
            return this->bb;
        }
        else {
            throw("NeuritePath::getBoundingBox(): bounding box has not been updated. use updateGeometry() first. this is a const method.");
        }
    }

    template<typename R>
    Vec3<R>
    NeuritePath<R>::findPermissibleRenderVector(
        R lambda,
        R mu) const
    {
        using Aux::Numbers::inf;

        /* check if geometry has been updated */
        if (!this->geometry_updated) {
            throw("NeuritePath::findPermissibleRenderVector(): geometry has not been updated. can't find render vector.\n");
        }

        Vec3<R> candidate = Aux::VecMat::zbase<R>();
        R                       r_min;
        Vec3<R>                 result;

        // We allow ourselves one educated guess before creating random candidates:
        // Check if global z vector is mu-permissible for all neurite segment curves of this path.
        r_min = inf<R>();
        for (auto &C : this->canal_segments_magnified)
            r_min = std::min(r_min, C->checkRenderVector(candidate));

        if (r_min > mu)
            return candidate;

        uint32_t const max_attempts = 1024;
        for (uint32_t attempt = 0; attempt < max_attempts; attempt++)
        {
            debugl(1, "NeuritePath::findPermissibleRenderVector(): attempt %d with %d candidates.\n", attempt + 1, 1024);

            candidate = Aux::VecMat::randUnitVec3<R>();

            /* check all candidates, pick alpha-permissible if possible, otherwise pick the best one if max_min is not
             * -inf */
			debugl(3, "NeuritePath()::findPermissibleRenderVector(): current candidate: ");
			#ifdef __DEBUG__
			candidate.print();
			#endif

			r_min = inf<R>();
			for (auto &C : this->canal_segments_magnified)
				r_min = std::min(r_min, C->checkRenderVector(candidate));

			debugl(3, "r_min: %5.4f", r_min);

            /* check if max_min > alpha, return result if that is the case, otherwise try the next candidate */
            if (r_min > lambda) {
                debugl(1, "attempt %d of %d, r_min = %5.4f => %5.4f-permissible for P. returning..\n", attempt + 1, max_attempts, r_min, lambda);
                return candidate;
            }
        }

        throw("NeuritePath::findPermissibleRenderVector(). no permissible render vector found.");
    }

    template <typename R>
    template <typename Tm, typename Tv, typename Tf>
    void
    NeuritePath<R>::generateInitialSegmentMesh(
        Mesh<Tm, Tv, Tf, R>                                    &M,
        uint32_t                                                n_phi_segments,
        Vec3<R>                                                 rvec,
        R const                                                &phi_0,
        R const                                                &arclen_dt,
        bool                                                   &end_circle_offset,
        std::vector<
                typename Mesh<Tm, Tv, Tf, R>::vertex_iterator
            >                                                  &end_circle_its,
        typename Mesh<Tm, Tv, Tf, R>::vertex_iterator          &closing_vertex_it,
        R const                                                &radius_reduction_factor,
		bool                                                    preserve_crease_edges) const
    {
        debugl(1, "NeuritePath::generateInitialSegmentMesh().\n");
        debugTabInc();

        /* check if geometry has been updated */
        if (!this->geometry_updated) {
            throw("NeuritePath::generateInitialSegmentMesh(): geometry has not been updated. can't find render vector.\n");
        }

        /* create new mesh for initial segment */
        /* if (this) is not a root path, generate mesh for initial segment mesh, but scale the radius r0 of the first
         * neurite canal segment by the factor radius_reduction_factor <= 1, which is the radius at the start vertex
         * of (this) path. details about this can be found in the thesis .. */
        if (!this->root_path) {
            debugl(2, "non-root path: getting copy of first neurite canal segment and applying radius reduction for start vertex. radius_reduction_factor = %f\n",
                radius_reduction_factor);

            /* copy first canal segment and reduce the start vertex radius with radius_reduction_factor */
            BLRCanalSurface<3u, R> const &C0    = *(this->canal_segments_magnified[0]);
            std::pair<R, R>     C0_radii    = C0.getRadii();
            BLRCanalSurface<3u, R>  C0_copy(
                    C0.getSpineCurve(),
                    C0_radii.first * radius_reduction_factor,
                    C0_radii.second
                );

            debugl(2, "generating mesh of reduced-radius initial segment..\n");
            C0_copy.template generateMesh<Tm, Tv, Tf>(
                /* append to given mesh M */
                M,
                /* pass n_phi_segments, render vector rvec, initial angle offset phi0 and arclen approximation precision
                 * from parameters */
                n_phi_segments,
                rvec,
                phi_0,
                arclen_dt,
                /* start circle not offset, supply no start circle iterators */
                false,
                NULL,
                /* retrieve end circle parity and iterators as well as closing vertex iterator */
                &end_circle_offset,
                &end_circle_its,
                &closing_vertex_it,
				preserve_crease_edges);

            debugl(2, "mesh generated.\n");
        }
        /* for root paths, also generate mesh for initial cylinder segment from neurite root edge incident to start
         * vertex of (this) path. */
        else {
            /* get neurite root segment incident to first vertex of this path */
            typename NLM_CellNetwork<R>::neurite_const_iterator         v0 = this->neurite_segments[0]->getSourceVertex();
            std::list<typename NLM_CellNetwork<R>::NeuriteRootEdge *>   v0_nres;

            v0->template getFilteredInEdges<typename NLM_CellNetwork<R>::NeuriteRootEdge>(v0_nres);

            if (v0_nres.size() == 1) {
                NLM::NeuriteRootEdgeInfo<R> const &re_info  = v0_nres.front()->neurite_root_edge_data;

                /* first, generate mesh for initial cylinder from the info about the neurite root edge, re_info. */
                re_info.initial_cylinder.template generateMesh<Tm, Tv, Tf>(
                    /* append to given mesh M */
                    M,
                    /* pass n_phi_segments, render vector rvec, initial angle offset phi0 and arclen approximation precision
                     * from parameters */
                    n_phi_segments,
                    rvec,
                    phi_0,
                    arclen_dt,
                    /* start circle not offset, supply no start circle iterators */
                    false,
                    NULL,
                    /* retrieve end circle parity and iterators as well as closing vertex iterator */
                    &end_circle_offset,
                    &end_circle_its,
                    &closing_vertex_it,
					preserve_crease_edges);

                /* erase the closing vertex to reopen the mesh */
                M.vertices.erase(closing_vertex_it);

                /* now append the mesh of canal segment 0 */
                this->canal_segments_magnified[0]->template generateMesh<Tm, Tv, Tf>(
                    /* append to given mesh M */
                    M,
                    /* pass n_phi_segments, render vector rvec, initial angle offset phi0 and arclen approximation precision
                     * from parameters */
                    n_phi_segments,
                    rvec,
                    phi_0,
                    arclen_dt,
                    /* hand over start circle parity / iterators */
                    end_circle_offset,
                    &end_circle_its,
                    /* pass pointers to the start circle variables above to seamlessly update these.  the semantics of
                     * CanalSurface::generateMesh() has been designed to support this scheme and performs sound updates.
                     * */
                    &end_circle_offset,
                    &end_circle_its,
                    /* save iterator to closing vertex */
                    &closing_vertex_it,
					preserve_crease_edges);
            }
            else {
                throw("NeuritePath::generateMesh(): (this) is a neurite root path whose first vertex not incident to "\
                    "exactly one neurite root edge. internal logic error.");
            }
            
            /*
            Mesh<Tm, Tv, Tf, R> CS;
            std::vector<
                    typename Mesh<Tm, Tv, Tf, R>::vertex_iterator
                >                                                   CS_end_circle_its;
            typename Mesh<Tm, Tv, Tf, R>::vertex_iterator           CS_closing_vertex_it;
            bool                                                    CS_end_circle_offset;
            */
        }

        debugTabDec();
        debugl(1, "NeuritePath::generateInitialSegmentMesh(): done.\n");
    }

    template <typename R>
    template<typename Tm, typename Tv, typename Tf>
    void
    NeuritePath<R>::appendTailMesh(
        Mesh<Tm, Tv, Tf, R>                                    &M,
        uint32_t                                                n_phi_segments,
        Vec3<R>                                                 rvec,
        R const                                                &phi_0,
        R const                                                &arclen_dt,
        bool const                                             &start_circle_offset,
        std::vector<
                typename Mesh<Tm, Tv, Tf, R>::vertex_iterator
            > const                                            &start_circle_its,
        typename Mesh<Tm, Tv, Tf, R>::vertex_iterator const    &start_circle_closing_vertex_it,
        bool                                                   *end_circle_offset,
        std::vector<
                typename Mesh<Tm, Tv, Tf, R>::vertex_iterator
            >                                                  *end_circle_its,
        typename Mesh<Tm, Tv, Tf, R>::vertex_iterator          *closing_vertex_it,
		bool                                                    preserve_crease_edges) const
    {
        debugl(1, "NeuritePath::appendTailMesh().\n");
        debugTabInc();

        /* check if geometry has been updated */
        if (!this->geometry_updated) {
            throw("NeuritePath::appendTailMesh(): geometry has not been updated. can't find render vector.\n");
        }

        uint32_t const m = this->neurite_segments.size();

        /* iterate through remaining segments 1, .. and append meshes on after the other. first, copy given start circle
         * parity and iterators into mutable local variables. */
        bool                                                        segment_i_start_circle_offset;
        std::vector<typename Mesh<Tm, Tv, Tf, R>::vertex_iterator>  segment_i_start_circle_its;
        typename Mesh<Tm, Tv, Tf, R>::vertex_iterator               segment_i_start_circle_closing_vertex_it;

        segment_i_start_circle_offset               = start_circle_offset;
        segment_i_start_circle_its                  = start_circle_its;
        segment_i_start_circle_closing_vertex_it    = start_circle_closing_vertex_it;


        debugl(2, "appending tail neurite canal segment meshes..\n");
        debugTabInc();
        for (uint32_t i = 1; i < m; i++) {
            debugl(2, "generating mesh for neurite canal segment %d\n", i);
            /* delete closing vertex of previous segment */
            M.vertices.erase(segment_i_start_circle_closing_vertex_it);

            this->canal_segments_magnified[i]->template generateMesh<Tm, Tv, Tf>(
                /* append to given mesh M */
                M,
                /* pass n_phi_segments, render vector rvec, initial angle offset phi0 and arclen approximation precision
                 * from parameters */
                n_phi_segments,
                rvec,
                phi_0,
                arclen_dt,
                /* start circle parity / iterators */
                segment_i_start_circle_offset,
                &segment_i_start_circle_its,
                /* if there is another segment, the end circle of the current segment is the start circle of the next.
                 * pass pointers to the start circle variables above to seamlessly update these for the next iteration.
                 * the semantics of CanalSurface::generateMesh() has been designed to support this scheme and performs
                 * sound updates. */
                &segment_i_start_circle_offset,
                &segment_i_start_circle_its,
                /* same scheme to update iterator to closing vertex */
                &segment_i_start_circle_closing_vertex_it,
				preserve_crease_edges);

            debugl(2, "done with mesh for neurite canal segment %d\n", i);
        }
        debugTabDec();

        /* all neurite canal segments have been meshed. start_circle_offset contains parity, start_circle_its contains
         * iterators of the END circle of the completed neurite path mesh. similarly,
         * segment_i_start_circle_closing_vertex_it identifies the closing vertex of the completed neurite path mesh. as
         * the closing vertex of the last neurite canal segment, it has not been deleted in the above loop, so the
         * neurite path mesh is closed.  write return values as desired by the caller. */
        if (end_circle_offset) {
            *end_circle_offset  = segment_i_start_circle_offset;
        }
        if (end_circle_its) {
            *end_circle_its     = segment_i_start_circle_its;
        }
        if (closing_vertex_it) {
            *closing_vertex_it  = segment_i_start_circle_closing_vertex_it;
        }

        debugTabDec();
        debugl(1, "NeuritePath::appendTailMesh(): done.\n");
    }
}

/* ----------------------------------------------------------------------------------------------------------------- *
 *
 *                        NLM_CellNetwork implementation.
 *
 * ----------------------------------------------------------------------------------------------------------------- */
template<typename R>
NLM_CellNetwork<R>::NLM_CellNetwork(std::string network_name) 
    : CellNetwork<
        NLM::CellNetworkInfo<R>,
        NLM::NeuronVertexInfo<R>,
        NLM::NeuronEdgeInfo<R>,
        NLM::SomaInfo<R>,
        NLM::NeuriteInfo<R>,
        NLM::AxonInfo<R>,
        NLM::DendriteInfo<R>,
        NLM::NeuriteSegmentInfo<R>,
        NLM::AxonSegmentInfo<R>,
        NLM::DendriteSegmentInfo<R>,
        NLM::NeuriteRootEdgeInfo<R>,
        NLM::AxonRootEdgeInfo<R>,
        NLM::DendriteRootEdgeInfo<R>,
        R
    >(network_name)
{
    this->nlm_neurite_root_vertices_adjusted        = false;
    this->nlm_network_info_updated                  = false;

    /* default settings for analysis and mesh generation */
    this->analysis_nthreads                         = 4;
    this->analysis_univar_solver_eps                = 1E-6;
    this->analysis_bivar_solver_eps                 = 1E-4;

    this->partition_filter_angle                    = M_PI / 2.0;
    this->partition_filter_max_ratio_ratio          = Aux::Numbers::inf<R>();
    this->partition_algo                            = this->partition_select_max_chordal_depth(
                                                            this->partition_filter_angle,
                                                            this->partition_filter_max_ratio_ratio);

    this->parametrization_algo                      = this->parametrization_chord_length();

    this->meshing_flush                             = true;
    this->meshing_flush_face_limit                  = 100000;

    this->meshing_canal_segment_n_phi_segments      = 12;
    this->meshing_outer_loop_maxiter                = 16;
    this->meshing_inner_loop_maxiter                = 8;

    this->meshing_preserve_crease_edges             = false;

    this->meshing_radius_factor_initial_value       = 0.975;
    this->meshing_radius_factor_decrement           = 0.01;
    this->meshing_complex_edge_max_growth_factor    = 2.0;
}

template <typename R>
NLM_CellNetwork<R>::~NLM_CellNetwork()
{
}

template <typename R>
typename NLM_CellNetwork<R>::Settings
NLM_CellNetwork<R>::getSettings() const
{
    NLM_CellNetwork<R>::Settings s;

    s.analysis_nthreads                         = this->analysis_nthreads;
    s.analysis_univar_solver_eps                = this->analysis_univar_solver_eps;
    s.analysis_bivar_solver_eps                 = this->analysis_bivar_solver_eps;
    /*
    s.partition_filter_angle                    = this->partition_filter_angle;
    s.partition_filter_max_ratio_ratio          = this->partition_filter_max_ratio_ratio;
    */
    s.partition_algo                            = this->partition_algo;
    s.parametrization_algo                      = this->parametrization_algo;;

    s.meshing_flush                             = this->meshing_flush;
    s.meshing_flush_face_limit                  = this->meshing_flush_face_limit;

    s.meshing_canal_segment_n_phi_segments      = this->meshing_canal_segment_n_phi_segments;
    s.meshing_outer_loop_maxiter                = this->meshing_outer_loop_maxiter;
    s.meshing_inner_loop_maxiter                = this->meshing_inner_loop_maxiter;

    s.meshing_preserve_crease_edges             = this->meshing_preserve_crease_edges;

    s.meshing_radius_factor_initial_value       = this->meshing_radius_factor_initial_value;
    s.meshing_radius_factor_decrement           = this->meshing_radius_factor_decrement;
    s.meshing_complex_edge_max_growth_factor    = this->meshing_complex_edge_max_growth_factor;

    return s;
}

template <typename R>
void
NLM_CellNetwork<R>::updateSettings(Settings const &s)
{

    this->analysis_nthreads                         = s.analysis_nthreads;
    this->analysis_univar_solver_eps                = s.analysis_univar_solver_eps;
    this->analysis_bivar_solver_eps                 = s.analysis_bivar_solver_eps;

    this->partition_algo                            = s.partition_algo;
    this->parametrization_algo                      = s.parametrization_algo;

    this->meshing_flush                             = s.meshing_flush;
    this->meshing_flush_face_limit                  = s.meshing_flush_face_limit;

    this->meshing_canal_segment_n_phi_segments      = s.meshing_canal_segment_n_phi_segments;
    this->meshing_outer_loop_maxiter                = s.meshing_outer_loop_maxiter;
    this->meshing_inner_loop_maxiter                = s.meshing_inner_loop_maxiter;

    this->meshing_preserve_crease_edges             = s.meshing_preserve_crease_edges;

    this->meshing_radius_factor_initial_value       = s.meshing_radius_factor_initial_value;
    this->meshing_radius_factor_decrement           = s.meshing_radius_factor_decrement;
    this->meshing_complex_edge_max_growth_factor    = s.meshing_complex_edge_max_growth_factor;

    printf("NLM_CellNetwork settings:\n"\
        "\t analysis_nthreads:                      %5d\n"\
        "\t analysis_univar_solver_eps:             %5.4e\n"\
        "\t analysis_bivar_solver_eps:              %5.4e\n"\
        "\t meshing_flush:                          %5d\n"\
        "\t meshing_flush_face_limit:               %5d\n"\
        "\t meshing_canal_segment_n_phi_segments:   %5d\n"\
        "\t meshing_outer_loop_maxiter:             %5d\n"\
        "\t meshing_inner_loop_maxiter:             %5d\n"\
		"\t meshing_preserve_crease_edges:          %s\n"\
        "\t meshing_radius_factor_initial_value:    %5.4f\n"\
        "\t meshing_radius_factor_decrement:        %5.4f\n"\
        "\t meshing_complex_edge_max_growth_factor: %5.4f\n\n",
        this->analysis_nthreads,
        this->analysis_univar_solver_eps,
        this->analysis_bivar_solver_eps,
        this->meshing_flush,
        this->meshing_flush_face_limit,
        this->meshing_canal_segment_n_phi_segments,
        this->meshing_outer_loop_maxiter,
        this->meshing_inner_loop_maxiter,
		this->meshing_preserve_crease_edges ? "true" : "false",
        this->meshing_radius_factor_initial_value,
        this->meshing_radius_factor_decrement,
        this->meshing_complex_edge_max_growth_factor);
}

template <typename R>
void
NLM_CellNetwork<R>::computeInitialNeuriteRootVertices()
{
    /* apply process to all somas of the network */
    for (auto &s : this->soma_vertices) {
        /* get reference to soma sphere S for s */
        NLM::SomaSphere<R> &S   = s.soma_data.soma_sphere;
        
        /* update info of soma sphere */
        S.centre()              = s.getSinglePointPosition();
        S.radius()              = s.getSinglePointRadius();

        /* consider every neurite of cell C_s rooted in s */
        std::list<NeuriteRootEdge *> s_nres;
        s.template getFilteredOutEdges<NeuriteRootEdge>(s_nres);

        for (auto &e : s_nres) {
            /* get current neurite root vertex u */
            auto u_it = e->getDestinationVertex();

            /* if u is not contained in S, insert new neurite root vertex x at the point of intersection between the
             * straight line connecting the centre of c and u and the sphere S. the radius for x is chosen as the radius
             * of u */
            if (!S.contains( u_it->getPosition() )) {
                Vec3<R> x_p = S.centre() + ( (u_it->getPosition() - S.centre()).normalize() ) * S.radius();
                R       x_r = u_it->getRadius();

                /* insert x, split edge (s, u) into (s, x) and (x, v). */
                this->neurite_root_edges.split( e->iterator(), { { x_p, x_r } } );
            }
            /* if u is contained in S, delete it and replace with its successor. this is done until the current neurite
             * root segment intersects S. the new neurite root vertex x is then inserted at the point of intersection
             * between the neurite root segment e and S. the radius chosen by linear interpolation between the radii of
             * e. */
            else while (true) {
                /* if u is NOT a branching vertex, proceed */
                if (!u_it->isNeuriteBranchingVertex()) {
                    /* get the child v of u and neurite segment (u, v) */
                    neurite_segment_iterator    uv_it   = u_it->getNeuriteSegmentFirstChild(); 
                    neurite_iterator            v_it    = uv_it->getDestinationVertex();

                    /* if v is contained in S as well, collapse e into v and set u to be newly created collapsed vertex,
                     * which is then processed in the next iteration of the while loop. */
                    if (S.contains(v_it->getPosition())) {
                        Vec3<R> const &c_p  = v_it->getPosition();
                        R       const &c_r  = v_it->getRadius();    
                        u_it                = this->neurite_segments.collapse(uv_it, { c_p, c_r });
                    }
                    /* otherwise v is outside of S. compute the point x of intersection between the straight line
                     * connecting (u -- v) and S and move u to x, where the radius is chosen by linear interpolation.
                     * break the loop afterwards. */
                    else {
                        Vec3<R> const p_u   = u_it->getPosition();
                        Vec3<R> const p_v   = v_it->getPosition();
                        Vec3<R> const d     = (p_v - p_u);

                        std::vector<R> lambda;
                        bool isec =
                            Aux::Geometry::raySphere(
                                S.centre(),
                                S.radius(),
                                p_u, 
                                d,
                                lambda);

                        R lambda_x;
                        /* NOTE: moderately subtle lazy evaluation here.. */
                        if (isec && lambda.size() == 2 && (lambda_x = std::max(lambda[0], lambda[1])) >= 0.0) {
                            Vec3<R> const p_x   = p_u + d * lambda_x;
                            R const r_x         = u_it->getRadius() + (v_it->getRadius() - u_it->getRadius())*lambda_x;

                            /* move / change u into x. */
                            u_it->setPosition(p_x);
                            u_it->setRadius(r_x);
                        }
                        else {
                            throw("NLM_CellNetwork::computeInitialNeuriteRootVertices(): in current neurite root "
                                "segment (u, v), u is contained in S, v is not, yet Aux::Geometry::raySphere() does "
                                "not report exactly two intersection points, one of which having a positive lambda value.");
                        }
                        break;
                    }
                }
                else {
                    throw("NLM_CellNetwork::computeInitialNeuriteRootVertices(): neurite branched prior to first vertex outside of soma sphere. won't split into multiple neurites.");
                }
            }
        }
    }

    /* update network info */
    this->network_info_initialized              = false;
    this->nlm_network_info_updated              = false;

    this->initializeNetworkInfo();
    this->updateNLMNetworkInfo();

    this->nlm_neurite_root_vertices_adjusted    = true;
}

template<typename R>
void
NLM_CellNetwork<R>::updateNLMNetworkInfo()
{
    if (!this->nlm_network_info_updated || !this->network_info_initialized) {
        debugl(1, "NLM_CellNetwork::updateNLMNetworkInfo().\n");
        debugTabInc();

        this->initializeNetworkInfo();

        /* for all cells (i.e. somas), init NLM_SomaInfo, traverse each neurite individually.  initialize
         * NLM::NeuriteInfo / NLM::NeuriteSegmentInfo containing e.g. information about the neurite path trees relevant
         * to the respective neurite vertex / neurite segment. this makes later access much more efficient and
         * comfortable.  */
        debugl(2, "iterating over all cells (somas)..\n");
        debugTabInc();
        for (auto &s : this->soma_vertices) {
            debugl(2, "processing soma %d\n", s.id());
            /* get reference to soma info and clear old data */
            NLM::SomaInfo<R>               &s_info  = s.soma_data;
            s_info.soma_sphere                      = NLM::SomaSphere<R>();
            s_info.neurite_path_trees.clear();

            std::list<NeuriteRootEdge *>    s_neurite_root_edges;

            s.template getFilteredOutEdges<NeuriteRootEdge>(s_neurite_root_edges);

            debugl(2, "processing all neurites connected to soma %d\n", s.id());
            debugTabInc();
            for (auto &nre : s_neurite_root_edges) {
                debugl(3, "processing neurite root edge %d with neurite root vertex r %d. retrieving connected component..\n", nre->id(), nre->getDestinationVertex()->id());

                debugl(3, "adding neurite path graph to NLM::SomaInfo for soma %d.\n", s.id());
                /* add neurite path graph for neurite root segment nre into soma info and get pointer to it. */
                s_info.neurite_path_trees.push_back(NeuritePathTree(nre->iterator()));

                NeuritePathTree *npt    = &(s_info.neurite_path_trees.back());

                /* extract neurite root vertex r of neurite starting with edge nre */
                neurite_iterator r      = nre->getDestinationVertex();

                debugl(3, "retrieving all neurite vertices / edges reachable from neurite roor vertex %d.\n", r->id());
                /* get all neurite vertices / segments reachable from r */
                std::list<NeuriteVertex *>  r_cc_nv;
                std::list<NeuriteSegment *> r_cc_ns;

                this->getNeuriteConnectedComponent(r, r_cc_nv, r_cc_ns);

                /* init NLM::NeuriteInfo inside all reachable vertices */
                debugl(3, "assigning basic NLM::NeuriteInfo to all %d neurite vertices reachable from r(%d).\n", r_cc_nv.size(), r->id());
                debugTabInc();
                for (auto &nv : r_cc_nv) {
                    debugl(4, "processing neurite vertex %d.\n", nv->id());
                    NLM::NeuriteInfo<R> &nv_info    = nv->neurite_data;
                    nv_info.npt                     = npt;
                    nv_info.npt_info.clear();
                }
                debugTabDec();
                debugl(3, "done assigning neurite vertex info.\n");

                /* init NLM::NeuriteSegmentInfo inside all reachable neurite segments */
                debugl(3, "assigning basic NLM::NeuriteSegmentInfo to all %d neurite segments reachable from r(%d).\n", r_cc_ns.size(), r->id());
                debugTabInc();
                for (auto &ns : r_cc_ns) {
                    debugl(4, "processing neurite segment %d.\n", ns->id());
                    NLM::NeuriteSegmentInfo<R> &ns_info = ns->neurite_segment_data;
                    ns_info.npt                         = npt;
                    ns_info.npt_it.explicitlyInvalidate();
                    ns_info.npt_ns_idx                  = 0;
                }
                debugTabDec();
            }
            debugTabDec();
            debugl(2, "done with soma %d.\n", s.id());
        }
        debugTabDec();
        debugl(1, "done processing somas.\n");

        this->nlm_network_info_updated = true;

        debugTabDec();
        debugl(1, "NLM_CellNetwork::updateNLMNetworkInfo(): done.\n");
    }
}

template <typename R>
bool
NLM_CellNetwork<R>::updateNsMDVInformation(neurite_segment_iterator ns_it)
{
    debugl(1, "NLM_CellNetwork::updateNsMDVInformation(): ns %d = (%d, %d)\n",
        ns_it->id(),
        ns_it->getSourceVertex()->id(),
        ns_it->getDestinationVertex()->id());

    debugTabInc();

    /* for all neurite segments, initialize minimum distance violation data / flags */
    neurite_const_iterator  const       u       = ns_it->getSourceVertex(),
                                        v       = ns_it->getDestinationVertex();
    R                                   ns_len  = ns_it->getLength();
    R                                   ns_rmax = ns_it->getMaxRadius();

    NLM::NeuriteSegmentInfo<R>         &ns_info = ns_it->neurite_segment_data;
    R                                   rmax_u_nb;
    R                                   rmax_v_nb;
    bool                                result  = false;

    debugl(1, "checking for PMDV..\n");
    /* check for PMDV: primary minimum distance violation */
    if (ns_len <= 2.0 * ns_rmax) {
        debugl(1, "PMDV ns: (%d, %d)\n", u->id(), v->id());

        ns_info.pmdv = true;
        debugTabInc();
        /* check for severe primary minimum distance violation */
        if (ns_len <= ns_rmax) {
            debugl(1, "SPMDV ns: (%d, %d)\n", u->id(), v->id());
            ns_info.spmdv = true;
        }
        else {
            ns_info.spmdv = false;
        }
        debugTabDec();

        result = true;
    }
    else {
        ns_info.pmdv    = false;
        ns_info.spmdv   = false;
    }

    debugl(1, "computing rmax_u_nb\n");
    if (!u->isNeuriteRootVertex()) {
        neurite_const_iterator const p = u->getNeuriteParent();
        rmax_u_nb   = std::max(u->getRadius(), p->getRadius());

        for (auto &nb : u->template getFilteredOutNeighbours<NeuriteVertex>()) {
            /* skip v */
            if (nb->iterator() != neurite_const_iterator(v)) {
                rmax_u_nb = std::max( rmax_u_nb, nb->getRadius() );
            }
        }
    }
    else {
        rmax_u_nb = 0.0;
    }

    debugl(1, "computing rmax_v_nb\n");
    /* get rmax_v_nb: v is the child of u, skip parent u */
    if (!v->isNeuriteTerminalVertex()) {
        rmax_v_nb = v->getRadius();
        for (auto &nb : u->template getFilteredOutNeighbours<NeuriteVertex>()) {
            rmax_v_nb = std::max( rmax_v_nb, nb->getRadius() );
        }
    }
    else {
        rmax_v_nb = 0.0;
    }

    debugl(1, "got neighbour radii quantities. checking for SMDV..\n");
    if (ns_len <= 1.0*(rmax_u_nb + rmax_v_nb) ) {
        debugl(1, "SMDV ns: (%d, %d)\n", u->id(), v->id());

        /* u corresponds to v_start, v to v_end in nsinfo */
        ns_info.smdv            = true;
        ns_info.v_src_rmax_nb   = rmax_u_nb;
        ns_info.v_dst_rmax_nb   = rmax_v_nb;

        result = true;
    }
    else {
        ns_info.smdv = false;
    }

    debugTabDec();
    debugl(1, "NLM_CellNetwork::updateNsMDVInformation(). done.\n");

    return result;
}

template <typename R>
void
NLM_CellNetwork<R>::updateAllMDVInformation()
{
    debugl(1, "NLM_CellNetwork::updateAllMDVInformation().\n");
    debugTabInc();

    for (auto &ns : this->neurite_segments) {
        this->updateNsMDVInformation(ns.iterator());
    }

    debugTabDec();
    debugl(1, "NLM_CellNetwork::updateAllMDVInformation(): done.\n");
}

template<typename R>
std::function<
    typename NLM_CellNetwork<R>::neurite_segment_iterator(
        NLM_CellNetwork<R> &,
        typename NLM_CellNetwork<R>::neurite_iterator
    )
>
NLM_CellNetwork<R>::partition_select_max_chordal_depth(
    R   filter_angle,
    R   max_radius_ratio)
{
    using Aux::Numbers::inf;

    return std::function<neurite_segment_iterator(NLM_CellNetwork<R> &, neurite_iterator)>(
        [filter_angle, max_radius_ratio] (
            NLM_CellNetwork    &C,
            neurite_iterator    it) -> neurite_segment_iterator 
        {
            neurite_segment_iterator    ret;
            R                           max_chordal_depth = -inf<R>(), tmp;

            std::list<NeuriteSegment *> nslist;
            it->template getFilteredOutEdges<NeuriteSegment>(nslist);

            for (auto &ns : nslist) {
                /* get chordal depth of neurite sub-tree rooted in destination vertex of ns */
                tmp     = (C.getNeuriteSubTreeDepths(ns->getDestinationVertex())).second;

                if (tmp > max_chordal_depth && ns->getAngle() < filter_angle && ns->getRadiusRatio() < max_radius_ratio) {
                    max_chordal_depth   = tmp;
                    ret                 = ns->iterator();
                }
            }

            if (max_chordal_depth == -inf<R>()) {
                ret.explicitlyInvalidate();
            }
            return ret;
        }
    );
}

template<typename R>
std::function<
    typename NLM_CellNetwork<R>::neurite_segment_iterator(
        NLM_CellNetwork<R> &,
        typename NLM_CellNetwork<R>::neurite_iterator
    )
>
NLM_CellNetwork<R>::partition_select_min_angle(
    R   filter_angle,
    R   max_radius_ratio)
{
    using Aux::Numbers::inf;

    return std::function<neurite_segment_iterator(NLM_CellNetwork<R> &, neurite_iterator)>(
        [filter_angle, max_radius_ratio] (
            NLM_CellNetwork    &C,
            neurite_iterator    it) -> neurite_segment_iterator 
        {
            neurite_segment_iterator    ret;
            R                           min_angle = inf<R>(), tmp;

            std::list<NeuriteSegment *> nslist;
            it->template getFilteredOutEdges<NeuriteSegment>(nslist);

            for (auto &ns : nslist) {
                /* get chordal depth of neurite sub-tree rooted in destination vertex of ns */
                tmp = ns->getAngle();

                /* if neurite segment is not filtered out and exhibits new minimum angle, use it */
                if (tmp < min_angle && tmp < filter_angle && ns->getRadiusRatio() < max_radius_ratio) {
                    min_angle   = tmp;
                    ret         = ns->iterator();
                }
            }
            
            /* if minimum angle is still infinity, no angle has been smaller than the filter angle. invalidate return iterator
             * so that comparison to all other iterators returns false */
            if (min_angle == inf<R>()) {
                ret.explicitlyInvalidate();
            }
            return ret;
        }
    );
}

template<typename R>
std::function<
    typename NLM_CellNetwork<R>::neurite_segment_iterator(
        NLM_CellNetwork<R> &,
        typename NLM_CellNetwork<R>::neurite_iterator
    )
>
NLM_CellNetwork<R>::partition_select_simple_neurite_paths()
{
    using Aux::Numbers::inf;

    return std::function<neurite_segment_iterator(NLM_CellNetwork<R> &, neurite_iterator)>(
        [] (
            NLM_CellNetwork    &C,
            neurite_iterator    it) -> neurite_segment_iterator 
        {
            neurite_segment_iterator    ret;
            ret.explicitlyInvalidate();
            return ret;
        }
    );
}

/* cell partitioning method */
template <typename R>
void
NLM_CellNetwork<R>::partitionCell(
    soma_iterator                       s_it,
    std::function<
            neurite_segment_iterator(
                NLM_CellNetwork<R> &,
                neurite_iterator
            )
        >  const                       &selection_algorithm)
{
    debugl(1, "NLM_CellNetwork::partitionCell(): partitioning cell tree rooted in soma %d\n", s_it->id());
    debugTabInc();

    /* initialize cell network information if needed (checked internally): NLM_SomaInfo, NLM::NeuriteInfo, .. annotated
     * to the various components of the network. in particular, this creates (empty) neurite path trees for every
     * neurite of the cell and hence also for the neurites connected to soma s. */
    this->initializeNetworkInfo();
    this->updateNLMNetworkInfo();

    /* iterator over all neurites, represented by neurite root edges, of soma s. */
    std::list<NeuriteRootEdge *> s_neurite_root_edges;
    s_it->template getFilteredOutEdges<NeuriteRootEdge>(s_neurite_root_edges);
    for (auto &nre : s_neurite_root_edges) {
        /* get neurite root vertex r */ 
        neurite_iterator r_it       = nre->getDestinationVertex();

        /* get reference to NLM::NeuriteInfo for neurite root vertex r. clear partitioning information for r (potentially
         * computed in a previous partitioning run). */
        NLM::NeuriteInfo<R> &r_info = r_it->neurite_data;
        r_info.npt_info.clear();

        /* get reference to NeuritePathTree for the current neurite by examining NLM::NeuriteInfo stored in r.  clear the
         * neurite path graph for the neurite root edge nre, which has potentially been computed previously during
         * another partitioning run. */
        NeuritePathTree &npt        = *(r_info.npt);
        npt.clear();

        /* stack of pairs (neurite_iterator, NeuritePathTree::vertex_iterator): the neurite_iterator points to a vertex v which is to
         * be added as the next vertex of the neurite path. the path might end at v, or it might be continued,
         * depending on the type of v and the selection algorithm */
        std::list<
                std::pair<
                    neurite_iterator,
                    typename NeuritePathTree::vertex_iterator
                >
            >                                                   S;

        /* for every neurite child c of r (usually 1, but v_r could be a neurite branching vertex, create a neurite
         * path P_c consisting only of v_r and push the pair (c, iterator to P_c in NeuritePathTree) onto the stack.
         * if computeNeuriteRootVertices has been called, every neurite root point lies on the soma sphere and the
         * corresponding neurite root vertex is never a branching vertex, so the following loop collapses to only one
         * out-going neurite root segment. */
        std::list<NeuriteSegment *> r_neurite_root_segments;
        r_it->template getFilteredOutEdges<NeuriteSegment>(r_neurite_root_segments);
        for (auto &r_ns : r_neurite_root_segments) {
            /* insert a vertex containing the one-edge neurite ROOT path { ns } into path NeuritePathTree npt and store
             * iterator. note the "true" argument of the NeuritePath<R> constructor, which identifies the path as a
             * neurite root path. */
            neurite_segment_iterator                    r_ns_it = r_ns->iterator();
            neurite_iterator                            c_it    = r_ns_it->getDestinationVertex();
            NLM::NeuritePath<R>                         P_c(r_ns_it, true);
            typename NeuritePathTree::vertex_iterator   npt_it = npt.vertices.insert(P_c);

            /* update NLM::NeuriteSegmentInfo for neurite root segment r_ns. */
            NLM::NeuriteSegmentInfo<R> &r_ns_info   = r_ns_it->neurite_segment_data;
            r_ns_info.npt                           = &npt;
            r_ns_info.npt_it                        = npt_it;
            r_ns_info.npt_ns_idx                    = 0;

            /* append information about the current path to NLM::NeuriteInfo for neurite root vertex r: r is the first
             * vertex of the path P_c contained in neurite path tree node referred to by npt_it. */
            r_info.npt_info.push_back({ npt_it, 0 });

            /* push pair (c_it, npt_it) onto stack */
            S.push_back( { c_it, npt_it } );
        }

        /* depth-first-style traversal of the neurite rooted in r */
        while (!S.empty()) {
            /* pop() currently visited neurite vertex v and corresponding neurite path tree iterator / neurite path */
            neurite_iterator                            v_it        = S.back().first;
            typename NeuritePathTree::vertex_iterator   npt_Pv_it   = S.back().second;
            NLM::NeuritePath<R>                        &P_v         = **npt_Pv_it;
            S.pop_back();

            /* get reference v_info to NLM::NeuriteInfo for v. v is currently the last vertex on the in general only
             * partially completed path P_v. assign that information to v_info, potentially overwriting data from a
             * previous partitioning run. */
            NLM::NeuriteInfo<R> &v_info = v_it->neurite_data;
            v_info.npt_info = { { npt_Pv_it, P_v.numVertices() - 1} };
             
            /* proceed depending on the properties of v: if v is a simple vertex, then the current path P_v will simply be
             * continued with the one and only child of v, i.e. P_v is extended with the neurite segment (v, c). */
            if (v_it->isNeuriteSimpleVertex()) {
                neurite_segment_iterator    vc_ns_it    = v_it->getNeuriteSegmentFirstChild();
                neurite_iterator            c_it        = vc_ns_it->getDestinationVertex();

                /* extend path with neurite segment (v, c) */
                uint32_t                    vc_path_idx = P_v.extend(vc_ns_it);

                /* retrieve and update neurite segment info of (v, c) */
                NLM::NeuriteSegmentInfo<R>  &vc_ns_info = vc_ns_it->neurite_segment_data;
                vc_ns_info.npt                          = &npt;
                vc_ns_info.npt_it                       = npt_Pv_it;
                vc_ns_info.npt_ns_idx                   = vc_path_idx;

                /* push onto stack */
                S.push_back( { c_it, npt_Pv_it } );
            }
            /* v is a branching vertex: decide which, if any, child vertex should continue the current path P_v using the
             * selection algorithm.
             *
             * create new paths, add the current vertex to all new paths and push all children along with their respective
             * paths onto the stack S. if the current path ends here, all paths will be new, otherwise exactly one child
             * will continue the current path P_v. */
            else if (v_it->isNeuriteBranchingVertex()) { 
                neurite_segment_iterator succ_ns_it = selection_algorithm(*this, v_it);
                /* otherwise the selection_algorithm has returned an explicitly invalidated neurite segment iterator,
                 * which indicates that no neurite segment has been selected => the path P_v ends here. start new paths
                 * for all neurite segments (v, c) connecting v to its children. */
                std::list<NeuriteSegment *> v_out_neurite_segments;
                v_it->template getFilteredOutEdges<NeuriteSegment>(v_out_neurite_segments);
                for (auto &vc_ns_ptr : v_out_neurite_segments) {
                    neurite_segment_iterator    vc_ns_it    = vc_ns_ptr->iterator();
                    neurite_iterator            c_it        = vc_ns_it->getDestinationVertex();

                    /* selection algorithm has returned neurite segment succ_ns == vc_ns: continue the path with the
                     * neurite segment succ_ns = (v, c) to the child c of v */
                    if (succ_ns_it == vc_ns_it) {
                        /* extend path with neurite segment succ_ns = (v, c) */
                        uint32_t                   succ_path_idx    = P_v.extend(succ_ns_it);

                        /* retrieve and update neurite segment info of succ_ns = (v, c) */
                        NLM::NeuriteSegmentInfo<R> &succ_ns_info    = succ_ns_it->neurite_segment_data;
                        succ_ns_info.npt                            = &npt;
                        succ_ns_info.npt_it                         = npt_Pv_it;
                        succ_ns_info.npt_ns_idx                     = succ_path_idx;

                        /* push onto stack */
                        S.push_back( { c_it, npt_Pv_it } );
                    }
                    /* create new path containing neurite segment (v, c), add edge to neurite path tree. */
                    else {
                        /* create non-root path { (n,s) } */
                        NLM::NeuritePath<R> P_vc(vc_ns_it, false);

                        /* create information struct for edge in neurite path tree. shared vertex is v, the index of v
                         * in P_v is the number of edges in P_v. */
                        NPTEdgeInfo e_info;
                        e_info.shared_vertex                = v_it;
                        e_info.shared_vertex_src_path_idx   = P_v.numEdges();

                        /* add new vertex for path P_vc and edge (P_v, P_vc) (containing information struct e_info) to neurite path tree */
                        typename NeuritePathTree::vertex_iterator   npt_Pvc_it  = npt.vertices.insert(P_vc);
                        auto npt_eins_rpair                                     = npt.edges.insert(npt_Pv_it, npt_Pvc_it, e_info);
                        if (!npt_eins_rpair.second) {
                            throw("NLM_CellNetwork::partitionCell(): failed to insert edge into neurite path tree. internal logic error.");
                        }
                        typename NeuritePathTree::edge_iterator     npt_e_it    = npt_eins_rpair.first;

                        /* retrieve and update neurite segment info of (v, c) */
                        NLM::NeuriteSegmentInfo<R> &vc_ns_info  = vc_ns_it->neurite_segment_data;
                        vc_ns_info.npt                          = &npt;
                        vc_ns_info.npt_it                       = npt_Pvc_it;
                        vc_ns_info.npt_ns_idx                   = 0;

                        /* currently visited vertex v it not only contained in P_v, but it is also the first vertex of
                         * the newly created path P_vc => append that information to v_info. */
                        v_info.npt_info.push_back( { npt_Pvc_it, 0 } );

                        /* push onto stack */
                        S.push_back( { c_it, npt_Pvc_it } );
                    }
                }
            }
            /* else v is a terminal vertex and the path ends here. depth-first style traversal will back-track. */
        }

        /* debug print all neurite paths for neurite nre */
#ifdef __DEBUG__
        debugl(1, "finished partitioning neurite with neurite root edge %d.\n", nre->id());
        debugl(1, "listing neurite paths:\n");
        debugTabInc();
        for (auto &npt_v : npt.vertices) {
            debugl(1, "neurite path tree vertex %d: contained path:..\n", npt_v.id());
            debugTabInc();
            for (auto &path_segment : npt_v->neurite_segments) {
                debugl(1, "neurite segment: %d\n", path_segment->id());
            }
            debugTabDec();
        }
        debugTabDec();
#endif
    }

    debugl(1, "checking integrity of partitioning..\n");
    /* all neurite segments belonging to soma s must be contained in exactly one neurite path. */
    debugTabInc();
    for (auto &ns : this->neurite_segments) {
        if (ns.getSoma() == s_it) {
            NLM::NeuriteSegmentInfo<R> &ns_info = ns.neurite_segment_data;
            if (ns_info.npt_it.explicitlyInvalid()) {
                debugl(1, "FATAL: NLM_CellNetwork::partitionCell(): neurite segment %d belonging to soma %d "\
                    "contains explicitly invalidated NeuritePathTree::vertex_iterator inside attached NLM_"\
                    "NeuriteInfo struct AFTER cell partitioning. internal logic error.", ns.id(), s_it->id());

                throw("NLM_CellNetwork::partitionCell(): discovered neurite segment with explicitly invalidated "\
                    "neurite path tree iterator contained in attached NLM::NeuriteInfo.");
            }
            debugl(2, "neurite segment %d: neurite: %d, path %d, segment index %d\n",
                ns.id(), ns.getNeurite()->id(), ns_info.npt_it->id(), ns_info.npt_ns_idx);
        }
    }

    for (auto &nv : this->neurite_vertices) {
        if (nv.getSoma() == s_it) {
            NLM::NeuriteInfo<R> &nv_info = nv.neurite_data;
            debugl(2, "neurite vertex %d (branching: %d): ..\n", nv.id(), nv.isNeuriteBranchingVertex());
            if (!nv.isNeuriteBranchingVertex() && nv_info.npt_info.size() != 1) {
                throw("CellNetwork::partitionCell(): found non-neurite branching vertex which is contained in more "\
                    "than one neurite path. internal logic error.");
            }

#ifdef __DEBUG__
            debugTabInc();
            for (auto &nv_npt_info : nv_info.npt_info) {
                debugl(2, "neurite: %d, path %d, vertex %d.\n",
                    nv.getNeurite()->id(), nv_npt_info.first->id(), nv_npt_info.second);
            }
            debugTabDec();
#endif
        }
    }
    debugTabDec();
    debugl(1, "partitioning seems ok.\n");

    debugTabDec();
    debugl(1, "NLM_CellNetwork::partitionCell(): done with cell tree rooted in soma %d\n", s_it->id());
}

template <typename R>
std::function<
        std::tuple<
            R,
            std::vector<R>,
            R
        >(NLM::NeuritePath<R> const &P)
    >
NLM_CellNetwork<R>::parametrization_chord_length()
{
    auto algo = std::function<std::tuple<R, std::vector<R>, R>(NLM::NeuritePath<R> const &P)>(
        [] (NLM::NeuritePath<R> const &P) -> std::tuple<R, std::vector<R>, R>
        {
            uint32_t const  n = P.numVertices();
            uint32_t        i;
            std::vector<R>  parameters(n, 0);

            /* sum up neurite segment lentsh to get entire chord length and store partial chord lengths in
             * neurite_vertex_parameters in the process. */
            R chordlen = 0;
            for (i = 1; i < n; i++) {
                chordlen       += P.neurite_segments[i-1]->getLength();
                parameters[i]   = chordlen;
            }

            /* cordlen now contains the entire chordlength. divide everything by cordlen to get relative chordlength in
             * [0,1] */
            parameters[0]   = 0;
            for (i = 1; i < n - 1; i++) {
                parameters[i] /= chordlen;
            }
            parameters[n-1] = 1;

            return (std::tuple<R, std::vector<R>, R>(chordlen, parameters, chordlen));
        }
    );

    return algo;
}

template <typename R>
std::function<
        std::tuple<
            R,
            std::vector<R>,
            R
        >(NLM::NeuritePath<R> const &P)
    >
NLM_CellNetwork<R>::parametrization_uniform()
{
    auto algo = std::function<std::tuple<R, std::vector<R>, R>(NLM::NeuritePath<R> const &P)>(
        [] (NLM::NeuritePath<R> const &P) -> std::tuple<R, std::vector<R>, R>
        {
            uint32_t const  n = P.numVertices();
            uint32_t        i;
            std::vector<R>  parameters(n, 0);

            parameters[0]   = 0;
            for (i = 1; i < n - 1; i++) {
                parameters[i]   = (R)i / R(n);
            }
            parameters[n-1] = 1;

            R ds0_scale = P.neurite_segments[0]->getLength() / parameters[1]; 
            R dsn_scale = P.neurite_segments[n-2]->getLength() / (1 - parameters[n-2]);

            return (std::tuple<R, std::vector<R>, R>(ds0_scale, parameters, dsn_scale));
        }
    );

    return algo;
}

template <typename R>
void
NLM_CellNetwork<R>::updateCellGeometry(
    soma_iterator   s_it,
    std::function<
            std::tuple<
                R,
                std::vector<R>,
                R
            >(NLM::NeuritePath<R> const &P)
        >  const                       &parametrization_algorithm)
{

    /* assumes that cell C_s identified by s has been partitioned. set soma sphere info, iterator over neurite path
     * trees and update canal surfaces for all of them. for neurite root paths, this also generates the initial cylinder
     * segments and stores them in the NLM::NeuriteRootEdgeInfo attached to all neurite root edges. */
    NLM::SomaInfo<R> &s_info    = s_it->soma_data;

    /* update information about soma sphere */
    s_info.soma_sphere.centre() = s_it->getSinglePointPosition();
    s_info.soma_sphere.radius() = s_it->getSinglePointRadius();

    /* update geometry of all neurites */
    for (auto &npt : s_info.neurite_path_trees) {
        for (auto &npt_v : npt.vertices) {
            npt_v->updateGeometry(parametrization_algorithm);
        }
    }
}

/* ----------------------------------------------------------------------------------------------------------------- *
 *
 *                                  static intersection analysis methods   
 *
 * ----------------------------------------------------------------------------------------------------------------- */
template <typename R>
bool
NLM_CellNetwork<R>::checkCanalSegmentRegularity(
    BLRCanalSurface<3u, R> const   &Gamma,
    R const                    &univar_solver_eps,
    std::vector<NLM::p2<R>>    &checkpoly_roots)
{
    bool                                    result;
    uint32_t                                i;
    R                                       val, feps, t_i;
    BernsteinPolynomial<4u, R, R>               gamma_reg;
    std::vector<PolyAlg::RealInterval<R>>   roots;

    /* compute "regularity" check polynomial, which is just the square of the parametric speed of
     * Gamma's spine curve */
    Gamma.spineCurveComputeRegularityPolynomial(gamma_reg);    

    /* get order of magnitude of regularity polynomial */
    feps   = gamma_reg.getMaxAbsCoeff();

    /* scale down to generous absolute error */
    feps   *= 1E-10;

    /* find roots of regularity polynomial with bezier clipping algorithm. */
    PolyAlg::BezClip_roots<4u, R>(gamma_reg, 0.0, 1.0, univar_solver_eps, roots);

    /* evaluate all candidate points and check against threshold */
    checkpoly_roots.clear();
    result = false;
    for (i = 0; i < roots.size(); i++) {
        t_i     = roots[i].midpoint();
        val     = gamma_reg.eval(t_i);

        if (val <= feps) {
            result = true;
            checkpoly_roots.push_back( NLM::p2<R>(roots[i], val) );
        }
    }
    return result;
}

template <typename R>
bool
NLM_CellNetwork<R>::checkSomaNeuriteIntersection(
    NLM::SomaSphere<R> const   &S,
    BLRCanalSurface<3u, R> const   &Gamma,
    bool                        neurite_root_segment,
    R const                    &univar_solver_eps,
    std::vector<NLM::p2<R>>    &isec_stat_points)
{
    //debugl(1, "(static) NLM_CellNetwork::checkSomaNeuriteIntersection()\n");
    bool                                    result;
    uint32_t                                i;
    R                                       r_S, r_Gamma_max, thres, dist, feps, t_i;
    Vec3<R>                                 S_c;
    BernsteinPolynomial<5u, R, R>               p;
    R const                                 offset = 2.0 * univar_solver_eps;
    std::vector<PolyAlg::RealInterval<R> >   roots;
    std::vector<PolyAlg::RealInterval<R> >   candidate_points;

    debugl(1, "getting data from soma sphere..\n");
    /* get soma centre, radius and max radius of pipe surface. compute threshold, which is (rmax + S_r)^2 */
    S_c         = S.centre();
    r_Gamma_max = Gamma.getMaxRadius();
    r_S         = S.radius();
    
    debugl(1, "computing threshold..\n");
    thres       = (r_Gamma_max + r_S);
    thres      *= thres;

    /* add boundary value candidate points t = 0.0 and t = 1.0 */
    candidate_points.push_back( PolyAlg::RealInterval<R>(0.0, 0.0));
    candidate_points.push_back( PolyAlg::RealInterval<R>(1.0, 1.0));

    debugl(1, "SONS: computing check polynomial..\n");

    /* get the check polynomial */
    Gamma.spineCurveComputeStationaryPointDistPoly(S_c, p);

    /* get order of magnitude of check polynomial */
    feps    = p.getMaxAbsCoeff();

    /* scale down to generous absolute error */
    feps   *= 1E-10;

    /* find roots of regularity polynomial and append roots to candidate_points. */
    PolyAlg::BezClip_roots<5u, R>(p, 0.0, 1.0, univar_solver_eps, roots);
    candidate_points.insert(candidate_points.end(), roots.begin(), roots.end());

    /* evaluate all candidate points and check against threshold */
    isec_stat_points.clear();
    result = false;
    for (i = 0; i < candidate_points.size(); i++) {
        /* if Gamma is a neurite root canal segment, skip all candidate points that are close to zero, since the canal
         * surface is supposed to intersect the soma aroung the starting point of Gamma */
        t_i     = candidate_points[i].midpoint();
        dist    = (S_c - Gamma.spineCurveEval(t_i)).len2squared();

        if (neurite_root_segment && t_i <= offset) {
            debugl(1, "skipping candidate point %5.4f for neurite root canal segment.\n", t_i); 
        } 
        else if (dist <= thres + feps) {
            result = true;
            isec_stat_points.push_back( NLM::p2<R>(candidate_points[i], dist) );
        }
    }
    return result;
}

/* local self intersection of neurite canal segment */
template <typename R>
bool
NLM_CellNetwork<R>::checkNeuriteLocalSelfIntersection(
    BLRCanalSurface<3u, R> const   &Gamma,
    R const                    &univar_solver_eps,
    std::vector<NLM::p2<R>>    &lsi_neg_points)
{
    bool                                    result;
    uint32_t                                i;
    R                                       feps, pval, t_i;
    BernsteinPolynomial<12u, R, R>               p_si;
    std::vector<PolyAlg::RealInterval<R>>   roots;
    std::vector<PolyAlg::RealInterval<R>>   candidate_points;

    /* add boundary value candidate points t = 0.0 and t = 1.0 */
    candidate_points.push_back( PolyAlg::RealInterval<R>(0.0, 0.0));
    candidate_points.push_back( PolyAlg::RealInterval<R>(1.0, 1.0));

    /* compute self-intersection polynomial of Gamma */
    Gamma.computeLocalSelfIntersectionPolynomial(p_si);

    /* get order of magnitude of self-intersection polynomial */
    feps    = p_si.getMaxAbsCoeff();

    /* scale down to generous absolute error */
    feps   *= 1E-10;

    /* find roots of regularity polynomial and append them to candidate_points */
    PolyAlg::BezClip_roots<12u, R>(p_si, 0.0, 1.0, univar_solver_eps, roots);
    candidate_points.insert(candidate_points.end(), roots.begin(), roots.end());

    /* check all candidate points. we got a focal point or a point between focal points if self-intersection polynomial
     * assumes a POSITIVE value.. */
    lsi_neg_points.clear();
    result = false;
    for (i = 0; i < candidate_points.size(); i++) {
        t_i     = candidate_points[i].midpoint();
        pval    = p_si.eval(t_i);

        if (pval > -feps) {
            result = true;
            lsi_neg_points.push_back( NLM::p2<R>(candidate_points[i], pval) );
        }
    }
    return result;
}

/* global self intersection of neurite canal segment */
template <typename R>
bool
NLM_CellNetwork<R>::checkNeuriteGlobalSelfIntersection(
    BLRCanalSurface<3u, R> const   &Gamma,
    R const                    &univar_solver_eps,
    R const                    &bivar_solver_eps,
    std::vector<NLM::p3<R>>    &gsi_stat_points)
{
    debugl(2, "NLM_CellNetwork::checkNeuriteGlobalSelfIntersection():\n");

    bool                                    result;
    uint32_t                                i;
    R                                       rmax_sum, thres, dist; //, feps;
    R const                                 offset = 2.0 * std::max(univar_solver_eps, bivar_solver_eps);

    BiBernsteinPolynomial<7u, 5u, R, R>             p;
    BiBernsteinPolynomial<5u, 7u, R, R>             q;
    std::vector<PolyAlg::RealRectangle<R>>  pq_roots;

    BernsteinPolynomial<5u, R, R>               pe_t0, pe_t1;
    std::vector<PolyAlg::RealInterval<R>>   edge_roots;

    std::vector<PolyAlg::RealRectangle<R>>  candidate_points;

    /* get sum of maximum radii */
    rmax_sum    = 2.0 * Gamma.getMaxRadius();
    thres       = rmax_sum * rmax_sum; 

    /* insert two corners of the unit square, represented by one corner candidate point due to symmetry in the GSI case. */
    candidate_points.push_back( { 0.0, 0.0, 1.0, 1.0} );

    /* construct intersection system */
    Gamma.computeGlobalSelfIntersectionSystem(p, q, pe_t0, pe_t1);

    /* obtain the order of magnitude of the system polynomials. get coefficient with max absolute
     * value over all coefficients of p and q */
    //feps = std::max( p.getMaxAbsCoeff(), q.getMaxAbsCoeff() );

    /* scale down to generous absolute error bound */
    //feps *= 1E-10;

    debugl(2, "solving bivariate system with bivariate linear clipping..\n");
    /* solve the system of bivariate polynomials with bivariate linear clipping. */

    /* unnecessarily, the bivariate linaer clipping implementation requires both polynomials to be in the same
     * basis, because the same legendre approximation matrices are used. this is not necessary, but inconvenient to
     * change right now, so: p is in BB(3n-2, 2n-1), q in BB(2n-1, 3n-2) => elevate p by (0, n-1) and q by (n-1, 0),
     * then both are in (3n-2, 3n-2) */
    BiBernsteinPolynomial<7u, 7u, R, R> p_elev, q_elev;
    p_elev = p.template elevateDegree<0,2u>();
    q_elev = q.template elevateDegree<2u,0>();
    try {
        PolyAlg::BiLinClip_roots<7u, 7u, R>(p_elev, q_elev, 0.0, 1.0, 0.0, 1.0, bivar_solver_eps, pq_roots);
    }
    catch (const char *err) {
        debugl(1, "checkNeuriteGlobalSelfIntersection(): caught exception from BiLinClip_roots: \'%s\'. outputting plot files of polynomial system and rethrowing.\n", err);
        p_elev.writePlotFile(200, "gsi_exception_p.plot");
        q_elev.writePlotFile(200, "gsi_exception_q.plot");

        /* exception caught, default to intersection */
        result = true;
    }

    /* append root rectangles returned by bivariate linear clipping to candidate points */
    candidate_points.insert(candidate_points.end(), pq_roots.begin(), pq_roots.end());

    debugl(2, "BilClip returned %ld roots.\n", pq_roots.size());

    /* solve _two_ edge polynomial systems and append respective roots, converted to rectangles, to candidate_points. */
    edge_roots.clear();
    PolyAlg::BezClip_roots<5u, R>(pe_t0, 0.0, 1.0, univar_solver_eps, edge_roots);
    for (i = 0; i < edge_roots.size(); i++) {
        //debugl(1, "t0 edge poly root interval (%20.13E, %20.13E)\n", edge_roots[i].x0, edge_roots[i].x1);
        candidate_points.push_back( { 0.0, 0.0, edge_roots[i].t0, edge_roots[i].t1 } );
                /*
                Vec2(
                    0.0,
                    (edge_roots[i].t0 + edge_roots[i].t1) / 2.0 
                ) 
                */
    }

    edge_roots.clear();
    PolyAlg::BezClip_roots<5u, R>(pe_t1, 0.0, 1.0, univar_solver_eps, edge_roots);
    for (i = 0; i < edge_roots.size(); i++) {
        //debugl(1, "t1 edge poly root interval (%20.13E, %20.13E)\n", edge_roots[i].x0, edge_roots[i].x1);
        candidate_points.push_back( { 1.0, 1.0, edge_roots[i].t0, edge_roots[i].t1 } );
                /*
                Vec2(
                    1.0,
                    (edge_roots[i].t0 + edge_roots[i].t1) / 2.0 
                ) 
                */
    }


    /* evaluate all candiate points, excepting the points that are equivalent to [1, 1] and [0,0]
     * within the numerical tolerance */
    gsi_stat_points.clear();
    result = false;
    for (i = 0; i < candidate_points.size(); i++) {
        auto z_i = candidate_points[i].midpoint();
        if ( (z_i[0] <= offset && z_i[1] <= offset) || (z_i[0] >= (1.0 - offset) && z_i[1] >= (1.0 - offset) ) ) {
            debugl(2, "\tgsi candidate point i = %d: (%20.13e, %20.13e) is being ignored..\n", i, z_i[0], z_i[1]);
        }
        else {
            //dist    = canalSurfacesSqDist(Gamma, Gamma, z_i[0], z_i[1]);
            dist    = (Gamma.spineCurveEval(z_i[0]) - Gamma.spineCurveEval(z_i[1])).len2squared();
            if (dist < thres) {
                debugl(1, "\tgsi candidate point i = %d: (%5.4e, %5.4e) has distance dist = %5.4f < thres = %f => gsi\n",
                        i, z_i[0], z_i[1], dist, thres);

                /* got a true candidate point. set result = true and append candidate points to gsi_stat_points. */
                result = true;
                gsi_stat_points.push_back( NLM::p3<R>(candidate_points[i], dist) );
            }
            else {
                debugl(2, "\tgsi candidate point i = %d: (%5.4e, %5.4e) has distance dist = %5.4f > thres = %f => no gsi\n",
                        i, z_i[0], z_i[1], dist, thres);
            }
        }
    }
    return result;
}

/* neurite / neurite intersection for non-adjacent neurite segment canal surfaces */
template <typename R>
bool
NLM_CellNetwork<R>::checkNeuriteNeuriteIntersection(
    BLRCanalSurface<3u, R> const   &Gamma,
    BLRCanalSurface<3u, R> const   &Delta,
    R const                    &univar_solver_eps,
    R const                    &bivar_solver_eps,
    std::vector<NLM::p3<R>>    &isec_stat_points)
{
    debugl(2, "NLM_CellNetwork::checkNeuriteNeuriteIntersection():\n");

    bool                                    result;
    uint32_t                                i;
    R                                       rmax_sum, thres, dist, feps;

    BiBernsteinPolynomial<5u, 3u, R, R>             p;
    BiBernsteinPolynomial<3u, 5u, R, R>             q;
    std::vector<PolyAlg::RealRectangle<R>>  pq_roots;

    BernsteinPolynomial<5u, R, R>               pe_x0, pe_x1, pe_y0, pe_y1;
    std::vector<PolyAlg::RealInterval<R>>   edge_roots;

    std::vector<PolyAlg::RealRectangle<R>>  candidate_points;

    /* get sum of maximum radii */
    rmax_sum    = Gamma.getMaxRadius() + Delta.getMaxRadius();
    thres       = rmax_sum * rmax_sum; 

    /* insert four corners of the unit square encoded as rectangles of diameter zero. */
    candidate_points.push_back( { 0.0, 0.0, 0.0, 0.0 } );
    candidate_points.push_back( { 0.0, 0.0, 1.0, 1.0 } );
    candidate_points.push_back( { 1.0, 1.0, 0.0, 0.0 } );
    candidate_points.push_back( { 1.0, 1.0, 1.0, 1.0 } );

    /* construct intersection system */
    Gamma.computeIntersectionSystem(Delta, p, q, pe_x0, pe_x1, pe_y0, pe_y1);

    /* obtain the order of magnitude of the system polynomials. get coefficient with max absolute
     * value over all coefficients of p and q */
    feps = std::max( p.getMaxAbsCoeff(), q.getMaxAbsCoeff() );

    /* scale down to generous absolute error bound */
    feps *= 1E-10;

    debugl(2, "solving bivariate system with bivariate linear clipping..\n");
    debugl(1, "NLM_CellNetwork::checkNeuriteNeuriteIntersection(): solving bivariate system with bivariate linear clipping..\n");
    /* solve the system of bivariate polynomials with bivariate linear clipping.
     * NOTE: unnecessarily, the bivariate linaer clipping implementation requires both polynomials
     * to be in the same basis, because the same legendre approximation matrices are used. this is
     * not necessary, but inconvenient to change right now, so:
     * p is in BB(2m-1, n), q in BB(m, 2n-1) => elevate p by (0, n-1) and q by (m-1, 0) */
    BiBernsteinPolynomial<5u, 5u, R, R> p_elev, q_elev;
    p_elev = p.template elevateDegree<0,2u>();
    q_elev = q.template elevateDegree<2u,0>();
    std::vector<PolyAlg::RealRectangle<R>> roots;
    try {
        PolyAlg::BiLinClip_roots<5u, 5u, R>(p_elev, q_elev, 0.0, 1.0, 0.0, 1.0, bivar_solver_eps, pq_roots);
    }
    catch (const char *err) {
        debugl(1, "checkNeuriteNeuriteIntersection(): caught exception from BiLinClip_roots: \'%s\'. outputting plot files of polynomial system and defaulting to intersection.\n", err);
        //p.writePlotFile(200, "gsi_exception_p.plot");
        //q.writePlotFile(200, "gsi_exception_q.plot");

        /* default to intersection */
        result = true;
        return true;
    }
    /* append root rectangles returned by bivariate linear clipping to candidate points */
    candidate_points.insert(candidate_points.end(), pq_roots.begin(), pq_roots.end());

    debugl(2, "BilClip returned %ld roots.\n", roots.size());
    debugl(1, "NLM_CellNetwork::checkNeuriteNeuriteIntersection(): bivariate linear clipping returned %ld roots.\n", roots.size());

    debugl(1, "NLM_CellNetwork::checkNeuriteNeuriteIntersection(): solving four univariate edge polynomial systems.\n");
    /* solve four edge polynomial systems and append respective roots, converted to rectangles, to candidate_points. */
    edge_roots.clear();
    PolyAlg::BezClip_roots<5u, R>(pe_x0, 0.0, 1.0, univar_solver_eps, edge_roots);
    for (i = 0; i < edge_roots.size(); i++) {
        candidate_points.push_back( { 0.0, 0.0, edge_roots[i].midpoint(), edge_roots[i].midpoint() } );
            /*
                Vec2(
                    0.0,
                    (edge_roots[i].x0 + edge_roots[i].x1) / 2.0 
                ) 
            */
        //debugl(1, "%20.15e %20.15e\n", candidate_points.back()[0], candidate_points.back()[1]);
    }

    edge_roots.clear();
    PolyAlg::BezClip_roots<5u, R>(pe_x1, 0.0, 1.0, univar_solver_eps, edge_roots);
    for (i = 0; i < edge_roots.size(); i++) {
        candidate_points.push_back( { 1.0, 1.0, edge_roots[i].midpoint(), edge_roots[i].midpoint() } );
            /*
                Vec2(
                    1.0,
                    (edge_roots[i].x0 + edge_roots[i].x1) / 2.0 
                ) 
            */
    }

    edge_roots.clear();
    PolyAlg::BezClip_roots<5u, R>(pe_y0, 0.0, 1.0, univar_solver_eps, edge_roots);
    for (i = 0; i < edge_roots.size(); i++) {
        candidate_points.push_back( { edge_roots[i].midpoint(), edge_roots[i].midpoint(), 0.0, 0.0 } );
            /*
                Vec2(
                    (edge_roots[i].x0 + edge_roots[i].x1) / 2.0,
                    0.0
                ) 
            */
    }

    edge_roots.clear();
    PolyAlg::BezClip_roots<5u, R>(pe_y1, 0.0, 1.0, univar_solver_eps, edge_roots);
    for (i = 0; i < edge_roots.size(); i++) {
        candidate_points.push_back( { edge_roots[i].midpoint(), edge_roots[i].midpoint(), 1.0, 1.0 } );
            /*
                Vec2(
                    (edge_roots[i].x0 + edge_roots[i].x1) / 2.0,
                    1.0
                ) 
            */
    }

    debugl(1, "NLM_CellNetwork::checkNeuriteNeuriteIntersection(): done. evaluating results..\n");

    /* evaluate all candiate points */
    isec_stat_points.clear();
    result = false;
    debugTabInc();
    for (i = 0; i < candidate_points.size(); i++) {
        auto z_i    = candidate_points[i].midpoint();
        //dist        = canalSurfacesSqDist(Gamma, Delta, z_i[0], z_i[1]);
        dist        = (Gamma.spineCurveEval(z_i[0]) - Delta.spineCurveEval(z_i[1])).len2squared();
        if (dist <= thres) {
            debugl(2, "candidate point i = %d: (%f, %f) has distance dist = %f < thres = %f => intersection\n",
                    i, z_i[0], z_i[1], dist, thres);

            debugl(1, "candidate point i = %d: (%f, %f) has distance dist = %f < thres = %f => intersection\n",
                    i, z_i[0], z_i[1], dist, thres);

            /* got a true candidate point. set result = true and append candidate points to isec_stat_points. */
            result = true;
            isec_stat_points.push_back( NLM::p3<R>(candidate_points[i], dist) );
        }
        else {
            debugl(2, "candidate point i = %d: (%f, %f) has distance dist = %f > thres = %f => no intersection\n",
                    i, z_i[0], z_i[1], dist, thres);
        }
    }
    debugTabDec();
    debugl(1, "NLM_CellNetwork::checkNeuriteNeuriteIntersection(): done...\n");
    return result;
}

/* same for adjacent neurite canal segments. if Gamma and Delta do not share their starting
 * point, then Gamma and Delta MUST be given in the order that  satisfies gamma(1.0) = delta(0.0),
 * i.e. the endpoint of Gamma's spine curve must be the starting point of Delta's spine curve. */
template <typename R>
bool
NLM_CellNetwork<R>::checkAdjacentNeuriteNeuriteIntersection(
    BLRCanalSurface<3u, R> const   &Gamma,
    BLRCanalSurface<3u, R> const   &Delta,
    R const                    &univar_solver_eps,
    R const                    &bivar_solver_eps,
    bool                        fst_end_snd_start,
    std::vector<NLM::p3<R>>    &isec_stat_points)
{
    debugl(2, "NLM_CellNetwork::checkAdjacentNeuriteNeuriteIntersection():\n");

    bool                                    result;
    uint32_t                                i;
    R                                       rmax_sum, thres, dist, feps;
    R const                                 offset = 1E-2; //2.0 * std::max(univar_solver_eps, bivar_solver_eps);

    BiBernsteinPolynomial<5u, 3u, R, R>     p;
    BiBernsteinPolynomial<3u, 5u, R, R>     q;
    std::vector<PolyAlg::RealRectangle<R>>  pq_roots;
    std::vector<PolyAlg::RealRectangle<R>>  pq_blacklist;

    BernsteinPolynomial<5u, R, R>               pe_x0, pe_x1, pe_y0, pe_y1;
    std::vector<PolyAlg::RealInterval<R>>   edge_roots;

    std::vector<PolyAlg::RealRectangle<R>>  candidate_points;

    /* get sum of maximum radii */
    rmax_sum    = Gamma.getMaxRadius() + Delta.getMaxRadius();
    thres       = rmax_sum * rmax_sum; 

    /* insert four corners of the unit square encoded as rectangles of diameter zero. */
    candidate_points.push_back( { 0.0, 0.0, 0.0, 0.0 } );
    candidate_points.push_back( { 0.0, 0.0, 1.0, 1.0 } );
    candidate_points.push_back( { 1.0, 1.0, 0.0, 0.0 } );
    candidate_points.push_back( { 1.0, 1.0, 1.0, 1.0 } );

    /* construct intersection system */
    Gamma.computeIntersectionSystem(Delta, p, q, pe_x0, pe_x1, pe_y0, pe_y1);

    /*
    p.writePlotFile(200, "p_adjacent.plot");
    q.writePlotFile(200, "q_adjacent.plot");
    */

    /* obtain the order of magnitude of the system polynomials. get coefficient with max absolute
     * value over all coefficients of p and q */
    feps = std::max( p.getMaxAbsCoeff(), q.getMaxAbsCoeff() );

    /* scale down to generous absolute error bound */
    feps *= 1E-10;

    debugl(2, "solving bivariate system with bivariate linear clipping..\n");
    /* solve the system of bivariate polynomials with bivariate linear clipping.  NOTE: unnecessarily, the bivariate
     * linaer clipping implementation requires both polynomials to be in the same basis, because the same legendre
     * approximation matrices are used. this is not necessary, but inconvenient to change right now, so: p is in
     * BB(2m-1, n), q in BB(m, 2n-1) => elevate p by (0, n-1) and q by (m-1, 0) */
    BiBernsteinPolynomial<5u, 5u, R, R> p_elev, q_elev;
    p_elev = p.template elevateDegree<0,2u>();
    q_elev = q.template elevateDegree<2u,0>();

    /* Gamma and Delta adjacent: blacklist a small rectangle around the corner point (1,0) or (0,0), depending on the
     * value of fst_end_snd_start:
     *
     * fst_end_snd_start == true:
     *
     *  (i) Gamma and Delta are consecutive segments on a path: both spine curves merge in a C2 fashion at the endpoint
     *  of gamma, which is the starting point of delta => tangent becomes identical, distance approaches zero => fat
     *  lines become parallel => bilinear clipping parallelogram system becomes singular => exception is thrown.  this
     *  can be avoided by blacklisting a small interval around (1, 0)
     *
     *  (ii) Gamma is the end segment of a path, Delta is the starting segment of another path, they share a branching
     *  vertex.  although there is in general no differentiable merging of the spine curves gamma and delta, blacklist
     *  interval anyway, because a single point on the corner of [0,1]^2 might cause numerical problems in any way,
     *  since there are no continously differentiable algebraic curves.
     *
     * fst_end_snd_start == false:
     *
     *  (i) both Gamma and Delta have the same starting point. blacklist small interval around (0,0)
     *
     *  */
    if (fst_end_snd_start) {
        return false;
        //pq_blacklist.push_back(PolyAlg::RealRectangle<R>( (1.0 - offset), 1.0, 0.0, offset) );
    }
    else {
        pq_blacklist.push_back(PolyAlg::RealRectangle<R>( 0.0, offset, 0.0, offset) );
    }

    /* call solver */
    try {
        PolyAlg::BiLinClip_roots<5u, 5u, R>(
                p_elev, q_elev,
                0.0, 1.0, 0.0, 1.0,
                bivar_solver_eps,
                pq_roots,
                /* disable dynamic recomputation of data to be thread-safe */
                false,
                /* use blacklist pq_blacklist */
                true, &pq_blacklist);
    }
    catch (const char *err) {
        debugl(1, "checkAdjacentNeuriteNeuriteIntersection(): caught exception from BiLinClip_roots: \'%s\'. outputting plot files of polynomial system..\n", err);
        debugl(1, "checkAdjacentNeuriteNeuriteIntersection(): caught exception from BiLinClip_roots: \'%s\'. outputting plot files of polynomial system..\n", err);

        //p_elev.writePlotFile(200, "consecutive_exception_p.plot");
        //q.writePlotFile(200, "consecutive_exception_q.plot");

        /* default to intersection */
        return !fst_end_snd_start;
    }

    /* append root rectangles returned by bivariate linear clipping to candidate points */
    candidate_points.insert(candidate_points.end(), pq_roots.begin(), pq_roots.end());

    debugl(2, "BilClip returned %ld roots.\n", pq_roots.size());

    /* solve four edge polynomial systems and append respective roots, converted to rectangles, to candidate_points. */
    edge_roots.clear();
    PolyAlg::BezClip_roots<5u, R>(pe_x0, 0.0, 1.0, univar_solver_eps, edge_roots);
    for (i = 0; i < edge_roots.size(); i++) {
        candidate_points.push_back( { 0.0, 0.0, edge_roots[i].midpoint(), edge_roots[i].midpoint() } );
            /*
                Vec2(
                    0.0,
                    (edge_roots[i].x0 + edge_roots[i].x1) / 2.0 
                ) 
            */
    }

    edge_roots.clear();
    PolyAlg::BezClip_roots<5u, R>(pe_x1, 0.0, 1.0, univar_solver_eps, edge_roots);
    for (i = 0; i < edge_roots.size(); i++) {
        candidate_points.push_back( { 1.0, 1.0, edge_roots[i].midpoint(), edge_roots[i].midpoint() } );
            /*
                Vec2(
                    1.0,
                    (edge_roots[i].x0 + edge_roots[i].x1) / 2.0 
                ) 
            */
    }

    edge_roots.clear();
    PolyAlg::BezClip_roots<5u, R>(pe_y0, 0.0, 1.0, univar_solver_eps, edge_roots);
    for (i = 0; i < edge_roots.size(); i++) {
        candidate_points.push_back( { edge_roots[i].midpoint(), edge_roots[i].midpoint(), 0.0, 0.0 } );
            /*
                Vec2(
                    (edge_roots[i].x0 + edge_roots[i].x1) / 2.0,
                    0.0
                ) 
            */
    }

    edge_roots.clear();
    PolyAlg::BezClip_roots<5u, R>(pe_y1, 0.0, 1.0, univar_solver_eps, edge_roots);
    for (i = 0; i < edge_roots.size(); i++) {
        candidate_points.push_back( { edge_roots[i].midpoint(), edge_roots[i].midpoint(), 1.0, 1.0 } );
            /*
                Vec2(
                    (edge_roots[i].x0 + edge_roots[i].x1) / 2.0,
                    1.0
                ) 
            */
    }

    /* evaluate all candiate points that are not too close to (1, 0) or (0,0), depending on
     * fst_end_snd_start flag, and compare them to the threshold */
    isec_stat_points.clear();
    result = false;
    debugTabInc();
    for (i = 0; i < candidate_points.size(); i++) {
        auto z_i = candidate_points[i].midpoint();

        if (( fst_end_snd_start && (z_i[0] >= 1.0 - offset   && z_i[1] <= offset) ) ||
            (!fst_end_snd_start && (z_i[0] <= offset         && z_i[1] <= offset) ) )
        {
            debugl(2, "candidate point i = %d: (%20.13e, %20.13e) is being ignored (fst_end_snd_start = %d)...\n", i, z_i[0], z_i[1], fst_end_snd_start);
        }
        else {
            debugl(2, "candidate point i = %d: (%20.13e, %20.13e) is considered (fst_end_snd_start = %d)...\n", i, z_i[0], z_i[1], fst_end_snd_start);
            //dist    = canalSurfacesSqDist(Gamma, Delta, z_i[0], z_i[1]);
            dist    = (Gamma.spineCurveEval(z_i[0]) - Delta.spineCurveEval(z_i[1])).len2squared();
            if (dist <= thres) {
                debugl(2, "candidate point i = %d: (%f, %f) has distance dist = %f < thres = %f => intersection\n",
                        i, z_i[0], z_i[1], dist, thres);

                /* got a true candidate point. set result = true and append candidate points to isec_stat_points. */
                result = true;
                isec_stat_points.push_back( NLM::p3<R>(candidate_points[i], dist) );
            }
            else {
                debugl(2, "candidate point i = %d: (%f, %f) has distance dist = %f > thres = %f => no intersection\n",
                        i, z_i[0], z_i[1], dist, thres);
            }
        }
    }
    debugTabDec();
    debugl(1, "NLM_CellNetwork::checkAdjacentNeuriteNeuriteIntersection(): done.\n");

    return result;
}

/* ----------------------------------------------------------------------------------------------------------------- *
 *
 *                                      thread related methods
 *
 * ----------------------------------------------------------------------------------------------------------------- */

template <typename R>
void
NLM_CellNetwork<R>::startWorkerThread(ThreadInfo *tinfo)
{
    using namespace Aux::Timing;

    IsecJob            *generic_job;
    REG_Job            *reg_job;
    LSI_Job            *lsi_job;
    GSI_Job            *gsi_job;
    SONS_Job           *sons_job;
    NSNS_Adj_Job       *nsns_adj_job;
    NSNS_NonAdj_Job    *nsns_nonadj_job;
   
    /* catch exceptions here */
    try {
        /* blocking lock: thread waits until it acquires the lock on the mutex in order to work safely. */
        tinfo->mutex.lock();

        /* measure time */
        tick(tinfo->thread_id + 1);

        /* loop over all jobs */
        #ifdef __DEBUG__
        uint32_t job_index = 0;
        #endif
        for (auto generic_job_shared_ptr : tinfo->job_list) {
            generic_job         = generic_job_shared_ptr.get();

            debugl(1, "Thread %2d: processing job %5d: type: %2d\n", tinfo->thread_id, job_index, generic_job->type());

            /* reset all specialized pointers to NULL */
            reg_job             = NULL;
            lsi_job             = NULL;
            gsi_job             = NULL;
            sons_job            = NULL;
            nsns_adj_job        = NULL;
            nsns_nonadj_job     = NULL;

            /* depending on the job type, down-cast to specialized job class and call solver with the
             * stored arguments */
            switch (generic_job->type()) {
                case JOB_REG:
                    reg_job                         = dynamic_cast<REG_Job *>(generic_job);
                    if (reg_job) {
                        BLRCanalSurface<3u, R> const &Gamma = reg_job->ns_it->neurite_segment_data.canal_segment_magnified;

                        reg_job->result                 = checkCanalSegmentRegularity(
                                Gamma,
                                reg_job->univar_solver_eps,
                                reg_job->checkpoly_roots);
                    }
                    else {
                        throw("(static) NLM_CellNetwork::startWorkerThread(): failed to down-cast generic job to specialized job of indicated type.");
                    }
                    break;
                    
                case JOB_LSI:
                    lsi_job = dynamic_cast<LSI_Job *>(generic_job);
                    if (lsi_job) {
                        BLRCanalSurface<3u, R> const &Gamma = lsi_job->ns_it->neurite_segment_data.canal_segment_magnified;

                        lsi_job->result                 = checkNeuriteLocalSelfIntersection(
                                Gamma,
                                lsi_job->univar_solver_eps,
                                lsi_job->lsi_neg_points);
                    }
                    else {
                        throw("(static) NLM_CellNetwork::startWorkerThread(): failed to down-cast generic job to specialized job of indicated type.");
                    }
                    break;
                    
                case JOB_GSI:
                    gsi_job = dynamic_cast<GSI_Job *>(generic_job);

                    if (gsi_job) {
                        BLRCanalSurface<3u, R> const &Gamma = gsi_job->ns_it->neurite_segment_data.canal_segment_magnified;

                        gsi_job->result                 = checkNeuriteGlobalSelfIntersection(
                                Gamma,
                                gsi_job->univar_solver_eps,
                                gsi_job->bivar_solver_eps,
                                gsi_job->gsi_stat_points);
                    }
                    else {
                        throw("(static) NLM_CellNetwork::startWorkerThread(): failed to down-cast generic job to specialized job of indicated type.");
                    }
                    break;
                    
                case JOB_SONS:
                    sons_job = dynamic_cast<SONS_Job *>(generic_job);

                    if (sons_job) {
                        BLRCanalSurface<3u, R> const &Gamma = sons_job->ns_it->neurite_segment_data.canal_segment_magnified;

                        sons_job->result                = checkSomaNeuriteIntersection(
                                sons_job->s_it->soma_data.soma_sphere,
                                Gamma,
                                sons_job->neurite_root_segment,
                                sons_job->univar_solver_eps,
                                sons_job->isec_stat_points);
                    }
                    else {
                        throw("(static) NLM_CellNetwork::startWorkerThread(): failed to down-cast generic job to specialized job of indicated type.");
                    }
                    break;
                    
                case JOB_NS_NS_ADJ:
                    nsns_adj_job            = dynamic_cast<NSNS_Adj_Job *>(generic_job);

                    if (nsns_adj_job) {
                        BLRCanalSurface<3u, R> const &Gamma = nsns_adj_job->ns_first_it->neurite_segment_data.canal_segment_magnified;
                        BLRCanalSurface<3u, R> const &Delta = nsns_adj_job->ns_second_it->neurite_segment_data.canal_segment_magnified;

                        nsns_adj_job->result            = checkAdjacentNeuriteNeuriteIntersection(
                                Gamma,
                                Delta,
                                nsns_adj_job->univar_solver_eps,
                                nsns_adj_job->bivar_solver_eps,
                                nsns_adj_job->fst_end_snd_start,
                                nsns_adj_job->isec_stat_points);
                    }
                    else {
                        throw("(static) NLM_CellNetwork::startWorkerThread(): failed to down-cast generic job to specialized job of indicated type.");
                    }
                    break;

                case JOB_NS_NS_NONADJ:
                    nsns_nonadj_job = dynamic_cast<NSNS_NonAdj_Job *>(generic_job);

                    if (nsns_nonadj_job) {
                        BLRCanalSurface<3u, R> const &Gamma = nsns_nonadj_job->ns_first_it->neurite_segment_data.canal_segment_magnified;
                        BLRCanalSurface<3u, R> const &Delta = nsns_nonadj_job->ns_second_it->neurite_segment_data.canal_segment_magnified;

                        nsns_nonadj_job->result         = checkNeuriteNeuriteIntersection(
                                Gamma,
                                Delta,
                                nsns_nonadj_job->univar_solver_eps,
                                nsns_nonadj_job->bivar_solver_eps,
                                nsns_nonadj_job->isec_stat_points);
                    }
                    else {
                        throw("(static) NLM_CellNetwork::startWorkerThread(): failed to down-cast generic job to specialized job of indicated type.");
                    }
                    break;

                default:
                    throw("NLM_CellNetwork::processIntersectionJobs(): unknown job type encountered.\n");
            }
            /* set job state */
            generic_job->job_state = JOB_DONE;
        }

        debugl(1, "Thread slot %02d: thread finished after %5.4f seconds.\n", tinfo->thread_id, tack(tinfo->thread_id + 1));

        /* set thread state to DONE */
        tinfo->thread_state = THREAD_DONE;

        /* unlock mutex */
        tinfo->mutex.unlock();
    }
    catch (const char *errmsg) {
        /* panic for now */
        printf("NLM_CellNetwork::startWorkerThread(): caught exception \"%s\". shutting down...\n", errmsg);
        exit(1);
    }
}


template <typename R>
void
NLM_CellNetwork<R>::getAllNeuritePaths(std::list<NLM::NeuritePath<R> *> &neurite_paths)
{
    neurite_paths.clear();

    /* iterate over all cells, then over all neurites within the cell and extract paths */
    for (auto &s : this->soma_vertices) {
        NLM::SomaInfo<R> &s_info = s.soma_data;
        for (auto &npt : s_info.neurite_path_trees) {
            for (auto &npt_v : npt.vertices) {
                neurite_paths.push_back( &(*npt_v) );   
            }
        }
    }
}

template <typename R>
void
NLM_CellNetwork<R>::getAllNeuritePaths(std::list<NLM::NeuritePath<R> const *> &neurite_paths) const
{
    neurite_paths.clear();

    /* iterate over all cells, then over all neurites within the cell and extract paths */
    for (auto &s : this->soma_vertices) {
        NLM::SomaInfo<R> const &s_info = s.soma_data;
        for (auto &npt : s_info.neurite_path_trees) {
            for (auto &npt_v : npt.vertices) {
                neurite_paths.push_back( &(*npt_v) );   
            }
        }
    }
}

/* compute all intersection jobs for one full analysis cycle */
template <typename R>
void
NLM_CellNetwork<R>::computeFullAnalysisIntersectionJobs(
    std::list<std::shared_ptr<IsecJob>> &job_queue) const
{
    using namespace Aux::Geometry::IntersectionTestResults;

    /* clear job queue passed by reference */
    job_queue.clear();

    /* independent of the soma or neurite it belongs to, every neurite segment of the network is checked for
     *
     *  1. regularity of its spine curve (REG)
     *
     *  2. local self-intersection (LSI)
     *
     *  3. global self-intersection (GSI)
     *
     * */
    for (auto &ns : this->neurite_segments) {
        job_queue.push_back( std::shared_ptr<IsecJob>(new REG_Job( ns.iterator(), this->analysis_univar_solver_eps)) );
        job_queue.push_back( std::shared_ptr<IsecJob>(new LSI_Job( ns.iterator(), this->analysis_univar_solver_eps)) );
        job_queue.push_back( std::shared_ptr<IsecJob>(new GSI_Job( ns.iterator(), this->analysis_univar_solver_eps, this->analysis_bivar_solver_eps)) );
    }

    /* compute list of all neurite paths. check every neurite path P 
     *
     * 1. check for bb intersection with all somas. is positive, check all neurite segments against that soma.
     *
     * 1. internally: special "adjacent" jobs for consecutive neurite segments. "normal non-adjacent" jobs for non-consecutive neurite segments.
     *
     * 2. check against all other paths Q for bb intersection. if positive, check all pairs of canal segments from P and
     * Q, taking care to catch special cases of incident neurite segments (intra-cell-intra-neurite but not intra-path)
     * */

    std::list<NLM::NeuritePath<R> const *>    np_list;
    this->getAllNeuritePaths(np_list);

    for (auto pit = np_list.begin(); pit != np_list.end(); ++pit) {
        NLM::NeuritePath<R> const &P    = **pit;
        auto P_bb                       = P.getBoundingBox();

        /* ----- analyze intersection between pairs of canal segments that are both from the SAME path P ---- */
        /* check all pairs (Gamma_i, Gamma_j) of neurite canal segments from path P for bounding box
         * intersection.  in the special case of consecutive segments (Gamma_i, Gamma_{i+1}), generate an
         * NSNS_Adj_Job, which is processed  using a special check that tolerates intersections around common
         * points. for j > (i+1), the canal segment (Gamma_i, Gamma_j) must not intersect at all =>
         * NSNS_NonAdj_Jobs are generated. */
        for (uint32_t i = 0; i < P.canal_segments_magnified.size() - 1; i++) {
            /* get reference to canal segment C_i as well as its bounding box */
            BLRCanalSurface<3u, R> const &P_Gamma_i = *(P.canal_segments_magnified[i]);
            auto P_Gamma_i_bb                   = P_Gamma_i.getBoundingBox();

            /* check for intersection of P_Gamma_i with all somas and generate SONS_Jobs if necessary. */
            for (auto &s_m : this->soma_vertices) {
                NLM::SomaInfo<R> const &s_m_info    = s_m.soma_data;
                NLM::SomaSphere<R> const &S_m       = s_m_info.soma_sphere;
                auto S_m_bb                         = S_m.getBoundingBox();

                if (S_m_bb && P_Gamma_i_bb) {
                    job_queue.push_back(std::shared_ptr<IsecJob>(
                            new SONS_Job(
                                s_m.iterator(),
                                P.neurite_segments[i],
                                this->analysis_univar_solver_eps
                            )
                        ));
                }
            }

            /* generate special-case job for adjacent neurite canal segments (Gamma_i, Gamma_{i+1}), which have to be
             * checked anyway, since their bounding boxes always intersect. */
            job_queue.push_back(std::shared_ptr<IsecJob>(
                    new NSNS_Adj_Job(
                        P.neurite_segments[i],
                        P.neurite_segments[i+1],
                        /* fst_end_snd_start == true, since endpoint of Gamma_i is starting point of Gamma_{i+1} */
                        true,
                        this->analysis_univar_solver_eps,
                        this->analysis_bivar_solver_eps
                    )
                ));
            
            /* check the rest: j = i + 2, .. */
            for (uint32_t j = i + 2; j < P.canal_segments_magnified.size(); j++) {
                BLRCanalSurface<3u, R> const &P_Gamma_j = *(P.canal_segments_magnified[j]);
                auto P_Gamma_j_bb                   = P_Gamma_j.getBoundingBox();

                /* check for bounding box intersection */
                if (P_Gamma_i_bb && P_Gamma_j_bb) {
                    debugl(1, "creating non-adj nsns job (from within one path P): (%d, %d) - (%d, %d)\n",
                            P.neurite_segments[i]->getSourceVertex()->id(),
                            P.neurite_segments[i]->getDestinationVertex()->id(),
                            P.neurite_segments[j]->getSourceVertex()->id(),
                            P.neurite_segments[j]->getDestinationVertex()->id());

                    job_queue.push_back(std::shared_ptr<IsecJob>(
                            new NSNS_NonAdj_Job(
                                P.neurite_segments[i],
                                P.neurite_segments[j],
                                this->analysis_univar_solver_eps,
                                this->analysis_bivar_solver_eps
                            )
                        ));
                }
            }
        }

        /* ----------------------------------------- iterate over all other paths Q --------------------------------- */
        auto qit = pit;
        ++qit;
        for ( ; qit != np_list.end(); ++qit) {
            NLM::NeuritePath<R> const &Q    = **qit;
            auto Q_bb                       = Q.getBoundingBox();

            /* only if bounding boxes of paths P and Q (P != Q) intersect, perform check on the canal segment level. */
            if (P_bb && Q_bb) {
                /* iterate over all pairs (P_e, Q_d) of neurite segments from P and Q. if (P_e, Q_d) are adjacent, i.e.
                 * are incident to a common neurite vertex, generate an NSNS_Adj_Job.  the corresponding special cases
                 * are handled in detail below. 
                 *
                 * otherwise, check the neurite canal segments P_Gamma(for P_c) and Q_Delta(of Q_d) for bounding box
                 * intersection and generate an NSNS_NonAdj_Job if necessary. */
                for (auto &P_c : P.neurite_segments) {
                    NLM::NeuriteSegmentInfo<R> const &P_c_info  = P_c->neurite_segment_data;
                    BLRCanalSurface<3u, R> const &P_Gamma           = P_c_info.canal_segment_magnified;
                    auto P_Gamma_bb                             = P_Gamma.getBoundingBox();

                    for (auto &Q_d : Q.neurite_segments) {
                        NLM::NeuriteSegmentInfo<R> const &Q_d_info  = Q_d->neurite_segment_data;

                        /* P_c and Q_d share the same starting (source) vertex */
                        if (P_c->getSourceVertex() == Q_d->getSourceVertex()) {
                            job_queue.push_back(std::shared_ptr<IsecJob>(
                                    new NSNS_Adj_Job(
                                        P_c,
                                        Q_d,
                                        /* endpoint of P_c is not start point of Q_d */
                                        false,
                                        this->analysis_univar_solver_eps,
                                        this->analysis_bivar_solver_eps
                                    )
                                ));
                        }
                        /* end (destination) vertex of P_c is the start (source) vertex of Q_d */
                        else if (P_c->getDestinationVertex() == Q_d->getSourceVertex()) {
                            job_queue.push_back(std::shared_ptr<IsecJob>(
                                    new NSNS_Adj_Job(
                                        P_c,
                                        Q_d,
                                        /* endpoint of P_c is start point of Q_d */
                                        true,
                                        this->analysis_univar_solver_eps,
                                        this->analysis_bivar_solver_eps
                                    )
                                ));
                        }
                        /* other way around: start (source) vertex of P_c is the end (destination) vertex of Q_d */
                        else if (P_c->getSourceVertex() == Q_d->getDestinationVertex()) {
                            job_queue.push_back(std::shared_ptr<IsecJob>(
                                    new NSNS_Adj_Job(
                                        /* reversed order! */
                                        Q_d,
                                        P_c,
                                        /* endpoint of Q_d is start point of P_c */
                                        true,
                                        this->analysis_univar_solver_eps,
                                        this->analysis_bivar_solver_eps
                                    )
                                ));
                        }
                        /* this must never happen in a cell-tree that exhibits the proper tree topology */
                        else if (P_c->getDestinationVertex() == Q_d->getDestinationVertex()) {
                            throw("NLM_CellNetwork::computeFullAnalysisIntersectionJobs(): discovered two "\
                                "neurite segments P_c and Q_d from same cell C_n and same neurite N_n_i that"
                                "have the same destination vertex => invalid topology of cell tree.");
                        }
                        /* P_c and Q_d are non-adjacent => check for bounding box intersection */
                        else {
                            BLRCanalSurface<3u, R> const &Q_Delta   = Q_d_info.canal_segment_magnified;
                            auto Q_Delta_bb                     = Q_Delta.getBoundingBox();

                            if (P_Gamma_bb && Q_Delta_bb) {
                                debugl(1, "creating non-adj nsns job (from two paths P !- Q): (%d, %d) - (%d, %d)\n",
                                        P_c->getSourceVertex()->id(),
                                        P_c->getDestinationVertex()->id(),
                                        Q_d->getSourceVertex()->id(),
                                        Q_d->getDestinationVertex()->id());

                                job_queue.push_back(std::shared_ptr<IsecJob>(
                                        new NSNS_NonAdj_Job(
                                            P_c,
                                            Q_d,
                                            this->analysis_univar_solver_eps,
                                            this->analysis_bivar_solver_eps
                                        )
                                    ));
                            }
                        }
                    }
                }
            }
        }
    }
}

/* thread-related methods */
template <typename R>
void
NLM_CellNetwork<R>::processIntersectionJobsMultiThreaded(
    uint32_t const                             &nthreads,
    std::list<std::shared_ptr<IsecJob>> const  &job_queue,
    std::list<std::shared_ptr<IsecJob>>        &results)
{
    /* resize thread slot vector */
    //this->thread_slots.assign(nthreads, { std::thread(), ThreadInfo() });
    this->thread_slots.clear();
    for (uint32_t i = 0; i < nthreads; i++) {
        this->thread_slots.push_back(ThreadInfo(i));
    }

    debugl(1, "NLM_CellNetwork::processIntersectionJobs(). number of jobs: %ld\n", job_queue.size());
    //debugl(1, "NLM_CellNetwork::processIntersectionJobs(). number of jobs: %ld\n", job_queue.size());

    /* if not single-threaded, set njobs_per_thread to 1/(5*nthreads) of the number of total jobs, but not less than 500. just a
     * heuristic setting to avoid too low or too high work load */
    uint32_t njobs_per_thread;
    if (nthreads <= 1) {
        njobs_per_thread = job_queue.size();
    }
    else {
        njobs_per_thread = std::max( 500u, (uint32_t)(job_queue.size() / (5 * nthreads) ) );
    }

    debugTabInc();
    uint32_t    njobs_remaining = job_queue.size();
    auto        jit             = job_queue.begin();
    while (jit != job_queue.end()) {
        /* check for free thread capacity */
        for (uint32_t i = 0; i < nthreads; i++) {
            /* only proceed if end of queue hasn't been reached yet, otherwise break the loop. this can happen because
             * for loop iterates over nthreads slots. if, e.g., only nthreads - 1 jobs are remaining, the queue will be
             * empty on the last slot check.. */
            if (jit != job_queue.end()) {
                /* try to lock the mutex of thread slot i. if this succeeds and the slot is not processing, get another
                 * job from the job queue and launch a thread for it */
                if (thread_slots[i].mutex.try_lock()) {
                    /* check if thread is IDLE or DONE. otherwise, if thread is in state PROCESSING, the main loop here was
                     * too quick to relock the mutex .. */
                    if (thread_slots[i].thread_state == THREAD_IDLE || thread_slots[i].thread_state == THREAD_DONE) {
                        /* if thread has been running, join with it to ensure that it is fully done. */
                        if (thread_slots[i].thread_state == THREAD_DONE) {
                            thread_slots[i].thread.join();
                        }

                        debugl(1, "Thread slot %2d: free. fetching results and preparing another thread..\n", i);

                        /* push results */
                        for (auto job : thread_slots[i].job_list) {
                            if (job->result) {
                                results.push_back(job);
                            }
                        }

                        /* reset thread_slots for slot i. this clear()s the contained job list and causes the reference
                         * count of all shared pointers contained in it to be decreased by one. */
                        thread_slots[i].reset();

                        /* get next max(job_queue.size(), njobs_per_thread) jobs from queue and load them into
                         * thread_slots struct */
                        for (uint32_t k = 0; k < njobs_per_thread && jit != job_queue.end(); k++, ++jit) {
                            thread_slots[i].job_list.push_back(*jit);
                        }

                        /* set state to processing. */
                        thread_slots[i].thread_state = THREAD_PROCESSING;

                        /* decrease remaining jobs counter */
                        njobs_remaining -= thread_slots[i].job_list.size();

                        /* create new thread */
                        debugl(1, "Thread slot %2d: launching new thread for %5zu jobs. job queue contains another %u waiting jobs.\n",
                                i, thread_slots[i].job_list.size(), njobs_remaining );

                        try {
                            thread_slots[i].thread =
                                std::move(
                                    std::thread(NLM_CellNetwork<R>::startWorkerThread, &(thread_slots[i]))
                                );
                        }
                        catch (std::system_error& err) {
                            throw("(static) NLM_CellNetwork::processIntersectionJobsMultiThreaded(): caught std::system-error from thread() constructor => system could not spawn thread.");
                        }
                    }
                    else if (thread_slots[i].thread_state == THREAD_PROCESSING) {
                        debugl(1, "main loop too quick to relock thread mutex %d. unlocking..\n", i);
                    }
                    else {
                        throw("NLM_CellNetwork::processIntersectionJobsMultiThreaded(): discovered thread with unknown value for ThreadInfo::thread_state. internal logic error.");
                    }

                    /* unlock the mutex */
                    thread_slots[i].mutex.unlock();
                }
                /* thread is still active */
                else if (thread_slots[i].thread_state == THREAD_PROCESSING) {
                    debugl(5, "Thread slot %2d: still active..\n", i);
                }
                else {
                    throw("NLM_CellNetwork::processIntersectionJobsMultiThreaded(): discovered thread with unknown value for ThreadInfo::thread_state. internal logic error.");
                }
            }
            /* end of queue reached in slot-check for loop. break the loop. */
            else {
                break; 
            }
        }

        /* sleep a few microseconds so as not to burn cpu time on the main thread with polling */
        usleep(50000);
    }
    debugTabDec();

    debugl(1, "NLM_CellNetwork::processIntersectionJobs(): job queue empty. joining..\n");

    /* lock all slot mutices. if slot i has ever been in use, join with it and retrieve results.
     * */
    for (uint32_t i = 0; i < nthreads; i++) {
        /* blocking lock of mutex for thread slot i */
        thread_slots[i].mutex.lock();

        /* if thread slot i has been used at all, then its state is now DONE, because the thread
         * releases the mutex only after being done and setting the state to done. the mutex is
         * locked in blocking fashion above. if thread slot i has never been used, it's state is
         * idle. the state processing must not occur here. */
        if (thread_slots[i].thread_state == THREAD_DONE) {
            /* the thread returns directly after releasing the mutex, but it's all asynchronous here. perform blocking
             * join() to prevent undefined behaviour.. */
            thread_slots[i].thread.join();

            /* push results */
            for (auto job : thread_slots[i].job_list) {
                if (job->result) {
                    results.push_back(job);
                }
            }

            /* reset to IDLE and clear job list */
            thread_slots[i].reset();
        }
        else if (thread_slots[i].thread_state == THREAD_IDLE) {
            debugl(1, "Thread slot %2d: never been used..\n", i);
        }
        else if (thread_slots[i].thread_state == THREAD_PROCESSING) {
            throw("NLM_CellNetwork::processIntersectionJobs(): discovered thread with state PROCESSING after locking mutex. internal logic error.");
        }
        else {
            throw("NLM_CellNetwork::processIntersectionJobsMultiThreaded(): discovered thread with unknown value for ThreadInfo::thread_state. internal logic error.");
        }
    }
    fflush(stdout);
    debugl(1, "NLM_CellNetwork::processIntersectionJobs(): joined. done.\n");
}

template <typename R>
void
NLM_CellNetwork<R>::readFromNeuroMorphoSWCFile(
    std::string     filename,
    bool const     &check_coincident_positions)
{
    NLM_CellNetworkBaseType::readFromNeuroMorphoSWCFile(filename, check_coincident_positions);
    this->computeInitialNeuriteRootVertices();
}

template <typename R>
void
NLM_CellNetwork<R>::partitionNetwork()
{
    /* update network info */
    this->initializeNetworkInfo();
    this->updateNLMNetworkInfo();

    /* partition all cells and update geometry */
    for (auto &s : this->soma_vertices) {
        this->partitionCell(s.iterator(), this->partition_algo );
    }
}

template <typename R>
void
NLM_CellNetwork<R>::updateNetworkGeometry()
{
    for (auto &s : this->soma_vertices) {
        this->updateCellGeometry(s.iterator(), this->parametrization_algo );
    }
}

/* one full analysis iteration on the entire cell network */
template <typename R>
bool
NLM_CellNetwork<R>::performFullAnalysis()
{
    debugl(1, "NLM_CellNetwork::performFullAnalysis().\n");
    debugTabInc();

    /* precompute data for polynomials / numerical solvers.  disable dynamic recomputation of inner products,
     * approximation data, etc.. i.e. set them to immutable. otherwise, multiple concurrent write-access during e.g.
     * dynamic recomputation of required values would result in undefined behaviour in multi-tnhreaded environments. all
     * quired values are precomputed beforehand and never touched during analysis. maximum BB degree 24 and BB bi-degree
     * (8, 8) is sufficient for the currently considered modelling settings (mainly cubic bezier curves as spine
     * curves). should the settings be insufficient, the solvers will throw exceptions. in single-threaded environments,
     * it is safe to dynamically recompute the data as needed. */
    // EDIT (mbreit, 06-01-2017): This is no longer necessary as the polynomials are now templated by their degree
    // and sizes of the precomputed arrays do no longer change (once initialized).
    // Same goes for BLRCanalSurfaces.
    /*
    BernsteinPolynomial<3u, R, R>::initBernsteinBasisInnerProducts(24);
    BernsteinPolynomial<3u, R, R>::setInnerProductDataImmutable();

    debugl(1, "various approximation data for univariate and bivariate numerical solvers..\n");
    PolyAlg::BiLinClip_getApproximationData<8u, 8u, R>(
            // disable dynamic recomputation after computation to be thread-safe.
            false);

    debugl(1, "global self-intersection data for maximum radius pipe surface approximation..\n");
    BLRCanalSurface<3u, R>::initGlobalSelfIntersectionData();
    */

    /* update mdv information */
    this->updateAllMDVInformation();

    debugl(1, "computing all intersection jobs for one full analysis interation..\n");
    /* compute all intersection jobs */
    std::list<std::shared_ptr<IsecJob>> job_list, intersection_list;
    this->computeFullAnalysisIntersectionJobs(job_list);

    debugl(1, "got %zu jobs => shuffling..\n", job_list.size());

    /* move job list into vector, shuffle, move back to job list */
    std::vector<std::shared_ptr<IsecJob>> job_vec(job_list.size());

    std::copy(job_list.begin(), job_list.end(), job_vec.begin());
    job_list.clear();

    //std::random_shuffle(job_vec.begin(), job_vec.end());

    std::copy(job_vec.begin(), job_vec.end(), std::back_inserter(job_list));
    job_vec.clear();

    printf("processing %zu intersection jobs using %d worker threads.\n", job_list.size(), this->analysis_nthreads); 

    /* shuffle jobs to generate more even load with high probability */
    //std::random_shuffle(job_list.begin(), job_list.end());

    /* process all jobs multi-threaded */
    this->processIntersectionJobsMultiThreaded(this->analysis_nthreads, job_list, intersection_list);
    bool clean = intersection_list.empty();

    if (!clean) {
    printf("intersection jobs processed: number of positive intersection results returned by solvers: %5zu. results in detail:\n",
        intersection_list.size());
    }
    else {
        printf("\t intersection jobs processed: network CLEAN.\n");
    }

    /* process intersections. get non-const iterators from down-cast job objects and attach intersection info in neurite
     * segments. */
    for (auto &job_sptr : intersection_list) {
        IsecJob *generic_job    = job_sptr.get();
        uint32_t type           = generic_job->type();
        if (!generic_job->result) {
            debugTabDec();
            throw("NLM_CellNetwork<R>::performFullAnalysis(): list of positive intersection jobs contains job whose "\
                "result value indicates no intersection. internal logic error.");
        }

        if (type == JOB_REG) {
            REG_Job        *reg_job             = dynamic_cast<REG_Job *>(generic_job);
            REG_IsecInfo   *reg_isec_info       = new REG_IsecInfo(*this, *reg_job); 

            /* invoke update method of neurite segment info, which stores the created REG_IsecInfo in smart pointer and
             * sets the lsi flag. same for other per-neurite-segment intersections below.. */
            NLM::NeuriteSegmentInfo<R> &ns_info = reg_isec_info->ns_it->neurite_segment_data;

            ns_info.updateREGStatus(reg_isec_info);

            printf("\t REG:       (%5d).\n", reg_job->ns_it->id());
        }
        else if (type == JOB_LSI) {
            LSI_Job        *lsi_job                     = dynamic_cast<LSI_Job *>(generic_job);
            LSI_IsecInfo   *lsi_isec_info               = new LSI_IsecInfo(*this, *lsi_job); 
            NLM::NeuriteSegmentInfo<R> &ns_info         = lsi_isec_info->ns_it->neurite_segment_data;

            ns_info.updateLSIStatus(lsi_isec_info);

            printf("\t LSI:       (%5d).\n", lsi_isec_info->ns_it->id());
        }
        else if (type == JOB_GSI) {
            GSI_Job        *gsi_job                     = dynamic_cast<GSI_Job *>(generic_job);
            GSI_IsecInfo   *gsi_isec_info               = new GSI_IsecInfo(*this, *gsi_job); 
            NLM::NeuriteSegmentInfo<R> &ns_info         = gsi_isec_info->ns_it->neurite_segment_data;
            ns_info.updateGSIStatus(gsi_isec_info);

            printf("\t GSI:       (%5d)\n", gsi_isec_info->ns_it->id());
        }
        else if (type == JOB_SONS) {
            SONS_Job            *sons_job               = dynamic_cast<SONS_Job *>(generic_job);
            IC_SONS_IsecInfo    *sons_isec_info         = new IC_SONS_IsecInfo(*this, *sons_job); 

            NLM::NeuriteSegmentInfo<R> &ns_info         = sons_isec_info->ns_it->neurite_segment_data;
            ns_info.ic_sons                             = true;

            printf("\t SONS:      (%5d, %5d).\n", sons_isec_info->s_it->id(), sons_isec_info->ns_it->id());
        }
        else if (type == JOB_NS_NS_ADJ || type == JOB_NS_NS_NONADJ) {
            NSNS_Job       *nsns_job                    = dynamic_cast<NSNS_Job *>(generic_job);

            std::shared_ptr<NSNS_IsecInfo>  sp          = this->generateNSNSIsecInfo(*nsns_job);  
            NSNS_IsecInfo  *nsns_isec_info              = sp.get();

            /* get references to neurite segments and attached info */
            neurite_segment_iterator &e_it      = nsns_isec_info->ns_first_it;
            neurite_segment_iterator &f_it      = nsns_isec_info->ns_second_it;
            NLM::NeuriteSegmentInfo<R> &e_info  = e_it->neurite_segment_data;
            NLM::NeuriteSegmentInfo<R> &f_info  = f_it->neurite_segment_data;

            /* set flags depending on the type of NSNS intersection. */
            uint32_t isec_type = nsns_isec_info->type();
            if (isec_type == RC_NSNS) {
                printf("\t RC-NSNS:   (%5d, %5d).\n", nsns_job->ns_first_it->id(), nsns_job->ns_second_it->id());
                e_info.rc_nsns      = true;
                f_info.rc_nsns      = true;
            }
            else if (isec_type == ICRN_NSNS) {
                printf("\t ICRN-NSNS: (%5d, %5d).\n", nsns_job->ns_first_it->id(), nsns_job->ns_second_it->id());
                e_info.icrn_nsns    = true;
                f_info.icrn_nsns    = true;
            }
            else if (isec_type == ICIN_NSNS || isec_type == ICIN_NSNSA) {
                printf("\t ICIN-NSNS: (%5d, %5d).\n", nsns_job->ns_first_it->id(), nsns_job->ns_second_it->id());
                e_info.icin_nsns    = true;
                f_info.icin_nsns    = true;
            }
            else {
                debugTabDec();
                throw("NLM_CellNetwork::performFullAnalysis(): intersection analysis returned NSNS_Job with unknown type. internal logic error.");
            }
        }
        else {
            debugTabDec();
            throw("NLM_CellNetwork::performFullAnalysis(): unknown job type discovered. internal logic error.\n");
        }
    }

    debugTabDec();
    debugl(1, "NLM_CellNetwork::performFullAnalysis(): done.\n");

    return clean;
}


/* ----------------------------------------------------------------------------------------------------------------- *
 *
 *                                      mesh generation methods 
 *
 * ----------------------------------------------------------------------------------------------------------------- */
#include "MeshAlgorithms.hh"

template <typename R>
template <typename Tm, typename Tv, typename Tf>
void
NLM_CellNetwork<R>::renderCellNetwork(std::string filename)
{
    debugl(1, "NLM_CellNetwork<R>::renderCellNetwork(): \"%s\".\n", filename.c_str());
    debugTabInc();

    using namespace RedBlue_ExCodes;

    std::list<typename NeuritePathTree::vertex_iterator>    npt_vertices_bfs_ordered;

    Mesh<Tm, Tv, Tf, R>                                     M_cell, M_cell_backup, M_S, M_P;

    bool                                                    end_circle_offset;
    std::vector<
            typename Mesh<Tm, Tv, Tf, R>::vertex_iterator
        >                                                   end_circle_its, circle_its_update,
                                                            circle_its_update_original;
    typename Mesh<Tm, Tv, Tf, R>::vertex_iterator           closing_vertex_it;

    /* initialize the cell mesh to consist of all soma spheres. while iterating over all somas, get breadth-first
     * ordering of neurite paths for all neurites and append to global list */
    debugl(1, "initializing soma sphere meshes and computing BFS ordering among neurite paths.\n");
    for (auto &s : this->soma_vertices) {
        NLM::SomaInfo<R> &s_info  = s.soma_data;

        s_info.soma_sphere.template generateMesh<Tm, Tv, Tf>(M_S);

        M_cell.moveAppend(M_S);

        /* iterate over all neurite path trees, i.e. all neurites, of the soma. */
        for (NeuritePathTree &npt : s_info.neurite_path_trees) {
            std::list<typename NeuritePathTree::Vertex const *>  source_vertices;
            for (auto &v : npt.vertices) {
                if (v.indeg() == 0) {
                    source_vertices.push_back(&v);
                }
            }

            for (auto &sv : source_vertices) {
                /* get connected component of source vertex in breadth-first order and append corresponding
                 * NeuritePath pointers to neurite path list. */
                std::list<typename NeuritePathTree::Vertex *> sv_cc;  

                uint32_t tid = npt.getFreshTraversalId();

                npt.getConnectedComponentBreadthFirst(
                    sv->iterator(),
                    tid,
                    &sv_cc,
                    NULL);

                for (auto &u : sv_cc) {
                    npt_vertices_bfs_ordered.push_back(u->iterator());
                }
            }
        }
    }

    /* initialize flush info */
    MeshAlg::MeshObjFlushInfo<Tm, Tv, Tf, R>    M_cell_flushinfo(filename);
    std::list<uint32_t>                         M_cell_flush_last_boundary_vertices_ids_backup;                          

    /* in the computed bread-first ordering, inductively append neurite path meshes */
    debugl(1, "processing neurite paths in BFS order.\n");
    debugTabInc();
    uint32_t np_idx = 1;
    for (auto npt_vit = npt_vertices_bfs_ordered.begin(); npt_vit != npt_vertices_bfs_ordered.end(); ++npt_vit) {
        debugl(2, "processing neurite path %d\n", (*npt_vit)->id());
        printf("\t meshing neurite path %5u of %5zu.\n", np_idx, npt_vertices_bfs_ordered.size() );

        /* get reference to neurite path */
        NLM::NeuritePath<R> const &P    = (*npt_vit)->vertex_data;

        /* find permissible render vector for neurite path */
        Vec3<R> render_vector;
        render_vector = P.findPermissibleRenderVector();

        /* variables for angular offset */
        R           phi_0;

        /* merge neurite path initial mesh segment with RedBlueUnion, catch exceptions, re-randomize / decrease radius
         * factor / handle errors as required */
        uint32_t    outer_loop_iter             = 0;
        uint32_t    inner_loop_iter             = 0;
        bool        done                        = false;
        bool        restore_M_cell              = false;
        bool        new_outer_iteration         = false;
        bool        break_inner_meshing_loop    = true;

        R           complex_edge_growth_factor  = 0.0;
        uint32_t    complex_edge_initial_count  = 0;

        R           radius_factor           = this->meshing_radius_factor_initial_value;

        /* compute a radius factor the corresponds to a generous lower bound on a safe radius, which can be obtained
         * by analysing the local neighbourhood of P's start vertex. */
        NeuriteVertex const &P_vstart       = *(P.neurite_segments[0]->getSourceVertex());
        R           P_vstart_radius         = P_vstart.getRadius();
        auto        P_vstart_vertex_nbs     = P_vstart.template getFilteredNeighbours<NeuriteVertex>();
        R           radius_factor_safe_lb   = 0.5;

        for (auto &x : P_vstart_vertex_nbs) {
            radius_factor_safe_lb = std::min(radius_factor_safe_lb, x->getRadius() / (2.0 * P_vstart_radius));
        }
        debugl(1, "radius_factor_safe_lb = %5.4f\n", radius_factor_safe_lb);

        /* if mesh has more than meshing_flush_face_limit faces, flush all vertices and faces to disk which
         * definitely don't participate in any merging operation that remains to be done. the set of respective
         * faces is computed as follows: for every remaining neurite path and all neurite canal segments of this
         * path, get the bounding box, extend it generously and locate all faces of M_cell that properly intersect
         * that bounding box. compute the union of all such sets of faces (over all neurite canal segments yet to be
         * meshed) and invert the set of resulting faces.
         *
         * all faces in the obtained set (list) will never be affected during a merging operation and can safely be
         * flushed to disk. */
        if (this->meshing_flush && M_cell.numFaces() > this->meshing_flush_face_limit) {
            printf("\t partial cell mesh has %5d > %5d (flush face limit)) faces. flushing definitely no longer needed parts of to disk.. ",
                    M_cell.numFaces(), this->meshing_flush_face_limit);
            fflush(stdout);

            std::list<typename Mesh<Tm, Tv, Tf, R>::Face *>  flush_faces, tmp;
            for (auto npt_wit = npt_vit; npt_wit != npt_vertices_bfs_ordered.end(); ++npt_wit) {
                NLM::NeuritePath<R> const &Q = (*npt_wit)->vertex_data;
                for (auto &Gamma : Q.canal_segments_magnified) {
                    auto Gamma_search_bb = (Gamma->getBoundingBox()).extend(0.1, Vec3<R>(1E-2, 1E-2, 1E-2));
                    M_cell.findFaces(Gamma_search_bb, tmp);
                    flush_faces.insert(flush_faces.end(), tmp.begin(), tmp.end());
                }
            }

            auto sortFct = [] (const typename Mesh<Tm, Tv, Tf, R>::Face* x, const typename Mesh<Tm, Tv, Tf, R>::Face* y) -> bool {return (x->id() < y->id());};
            auto uniqueFct = [] (const typename Mesh<Tm, Tv, Tf, R>::Face* x, const typename Mesh<Tm, Tv, Tf, R>::Face* y) -> bool {return (x->id() == y->id());};
            flush_faces.sort(sortFct);
            flush_faces.unique(uniqueFct);

            /* invert face selection to get set of definitely not affected faces for all merging operations yet to
             * be performed. */
            M_cell.invertFaceSelection(flush_faces);

            /* .. and perform the flush */
            try {MeshAlg::partialFlushToObjFile(M_cell, M_cell_flushinfo, flush_faces);}
            catch (...) {debugTabDec(); debugTabDec(); throw;}

            printf("done.\n");
        }

        /* backup (potentially just partially flushed) cell mesh */
        M_cell_backup                   = M_cell;

        /* backup ids of all boundary_vertices (referring to M_cell) from M_cell_flushinfo */
        M_cell_flush_last_boundary_vertices_ids_backup.clear();
        for (auto &vp : M_cell_flushinfo.last_boundary_vertices) {
            M_cell_flush_last_boundary_vertices_ids_backup.push_back(vp.first->id());
        }

        /*
        tmp = M_cell;
        tmp.writeObjFile("M_cell_before_merge");
        */

        debugl(1, "entering outer meshing loop for path %d.\n", (*npt_vit)->id());
        debugTabInc();
        while (!done) {
            /* if necessary, use back M_cell_backup to restore original M_cell before the RedBlueUnion call if
             * necessary. this is the case iff the exception that lead to the necessity of another run indicated
             * R_intact == true (M_cell has been used as the red mesh). */
            if (restore_M_cell) {
                /* NOTE: the mesh flush info struct contains _pointers_ to boundary vertices of the old M_cell. once the
                 * backup is used to restore, the current M_cell is effectively clear()ed and overwritten with the
                 * backup copy by Mesh::operator=(). this process frees all vertices and faces of the mesh to which the
                 * backup copy is assigned, which in turn causes dangling invalid pointers inside M_cell_flushinfo. it
                 * therefore becomes necessary to update all last_boundary_vertices pointer to point to the correct
                 * vertices of M_cell after the restore. however, mesh assignment with Mesh::operator=() preserves
                 * integer ids, and so the update can be achieved as follows:
                 *
                 * 1. restore M_cell from M_cell_backup with assignment.
                 *
                 * 2. locate all boundary vertex ids in the restored M_cell and update pointers in M_cell_flushinfo. the
                 * reference to M_cell stays intact, since the object itself is never deleted but only overwritten by
                 * assignment. the integer ids of all boundary vertices are stored in the list
                 * M_cell_boundary_vertex_ids_backup.
                 *
                 * */
                debugl(1, "restoring M_cell from backup via assignment.\n");
                M_cell = M_cell_backup;

                debugl(1, "updating pointers in M_flush_info.last_boundary_vertices via id lookup..\n");
                auto fp_it = M_cell_flushinfo.last_boundary_vertices.begin();
                for (auto &v_id : M_cell_flush_last_boundary_vertices_ids_backup) {
                    auto v_it = M_cell.vertices.find(v_id);
                    if (v_it != M_cell.vertices.end()) {
                        fp_it->first = &(*v_it);
                        ++fp_it; 
                    }
                    else {
                        debugTabDec(); debugTabDec(); debugTabDec();
                        throw("NLM_CellNetwork::renderCellNetwork(): failed to locate flush boundary vertex via id during cell mesh backup.");
                    }
                }

                if (fp_it != M_cell_flushinfo.last_boundary_vertices.end()) {
                    debugTabDec(); debugTabDec(); debugTabDec();
                    throw("NLM_CellNetwork::renderCellNetwork(): iterator to MeshObjFlushInfo::last_boundary_vertices has not reached end() after cell mesh backup. internal logic error.");
                }

                debugl(1, "M_cell and flush info fully restored.\n");
            }

            new_outer_iteration         = false;
            restore_M_cell              = false;
            outer_loop_iter++;

            if (outer_loop_iter % this->meshing_outer_loop_maxiter == 0) {
                radius_factor   = std::max(radius_factor_safe_lb, radius_factor - this->meshing_radius_factor_decrement);
            }

            debugl(1, "outer meshing loop: iteration %d. radius factor: %5.4f\n", outer_loop_iter, radius_factor);
            debugl(1, "choosing random phi_0 and generating initial mesh segment..\n");

            /* compute random angular offset phi_0 */
            phi_0 = Aux::Numbers::frand(0.0, (2*(R)M_PI) / (R)this->meshing_canal_segment_n_phi_segments);

            /* clear path mesh */
            M_P.clear();

            /* generate P's initial segment mesh and append to M */
            try
            {
                P.template generateInitialSegmentMesh<Tm, Tv, Tf>(
                    /* append to mesh M_P for path P*/
                    M_P,
                    /* n_phi_segments default to 16 for testing */
                    this->meshing_canal_segment_n_phi_segments,
                    /* render vector */
                    render_vector,
                    /* phi_0, arclen_dt = 1E-3 */
                    phi_0,
                    1E-3,
                    /* end circle info */
                    end_circle_offset,
                    end_circle_its,
                    closing_vertex_it,
                    /* radius factor, which is being ignored for neurite root paths. */
                    radius_factor,
                    this->meshing_preserve_crease_edges);
            }
            catch (...) {debugTabDec(); debugTabDec(); debugTabDec(); throw;}

            /* triangulate M_P for RedBlueAlgorithm */
            M_P.triangulateQuads();

            /*
            tmp = M_P;
            tmp.writeObjFile("M_P_before_merge");
            */

            /* merge initial path mesh segment into cell mesh with RedBlueUnion. to append the tail of the path mesh, it
             * is required to know the end circle and closing vertex iterators AFTER merging, so these are assembled
             * into the update iterator vector circle_its_update for the RedBlueUnion call. */
            circle_its_update_original = end_circle_its;
            circle_its_update_original.push_back(closing_vertex_it);

            /* inner meshing loop: while the initial mesh segment generated above might still be usable (e.g. by
             * splitting complex edges), try to use it. as soon as exception handling sets new_outer_iteration or
             * break_inner_meshing_loop is set to false at the end of the try {..} block, the inner loop breaks.
             * directly after the body of the inner meshing loop, it is checked whether new_outer_iteration == true. */
            break_inner_meshing_loop    = false;
            inner_loop_iter             = 0;
            complex_edge_growth_factor  = 0.0;
            complex_edge_initial_count  = 0;

            debugl(1, "entering inner meshing loop..\n");
            debugTabInc();
            while (!new_outer_iteration && !break_inner_meshing_loop && inner_loop_iter < this->meshing_inner_loop_maxiter) {
                /* copy circle_its_update_original into circle_its_update for current meshing run */
                circle_its_update = circle_its_update_original;
                inner_loop_iter++;

                debugl(2, "calling RedBlueUnion algorithm to merge initial mesh segment into partially completed cell mesh.\n");

                Aux::Timing::tick(14);

                try {
                    MeshAlg::RedBlueUnion<Tm, Tv, Tf, R>(
                        /* R = M_cell, which is to be union mesh afterwards */
                        M_cell,
                        /* B = M_P, the mesh for the initial segment of P */
                        M_P,
                        /* list of end circle iterators from M_P which are updated to reflect the corresponding vertices in
                         * the union mesh */
                       &circle_its_update);

                    /* RedBlueUnion call has been succcessful. break inner meshing loop */
                    break_inner_meshing_loop = true;
                }
                catch (RedBlue_Ex_InternalLogic& logic_ex) {
                    debugTabDec(); debugTabDec(); debugTabDec(); debugTabDec();
                    throw;
                }
                catch (RedBlue_Ex_Disjoint& disjoint_ex) {
                    debugTabDec(); debugTabDec(); debugTabDec(); debugTabDec();
                    throw;
                }
                catch (RedBlue_Ex_ComplexEdges<R>& complex_ex) {
                    debugl(0, "NLM_CellNetwork::renderCellNetwork(): RedBlueAlgorithm returned exception: %d complexly intersecting edges.. splitting.\n", complex_ex.edge_isec_info.size());
                    debugTabInc();

                    uint32_t const nce = complex_ex.edge_isec_info.size();

                    /* if this is the first complex edge exception, set initial complex edge count */
                    if (complex_edge_growth_factor == 0.0) {
                        complex_edge_initial_count  = nce;
                        complex_edge_growth_factor  = 1.0;
                        debugl(1, "setting complex edge initial count to %d, factor to %5.4f\n", nce, complex_edge_growth_factor);
                    }
                    /* otherwise calculate "growth factor". in certain situations, splitting all complex edges creates
                     * even more complex edges. this process can amplify itself exponentially. to prevent this, check
                     * if the growth factor exceeds a certain limit */
                    else {
                        complex_edge_growth_factor  = (R)nce / (R)complex_edge_initial_count;
                        debugl(1, "setting complex edge growth factor to %5.4f\n", complex_edge_growth_factor);
                    }

                    /* if complex edge growth factor is too large, start new outer meshing iteration */
                    if (complex_edge_growth_factor > this->meshing_complex_edge_max_growth_factor) {
                        debugl(0, "Complex edge growth factor too large.\n");
                        radius_factor   = std::max(radius_factor_safe_lb, radius_factor - this->meshing_radius_factor_decrement);
                        new_outer_iteration = true;
                        restore_M_cell      = true;
                    }
                    else for (auto e_info : complex_ex.edge_isec_info) {
                        debugl(2, "complex edge: e = (%d, %d)\n", e_info.u_id, e_info.v_id);

                        /* sort and check lambdas */
                        debugl(2, "sorting complex exception edge lambdas. size(): %zu..\n", e_info.edge_lambdas.size() );
                        std::sort(e_info.edge_lambdas.begin(), e_info.edge_lambdas.end(), std::less<R>() );
                        debugTabInc();
                        for (auto &lambda : e_info.edge_lambdas) {
                            debugl(3, "lambda: %5.4f\n", lambda);
                            if (lambda <= 0 || lambda >= 1) {
                                debugTabDec(); debugTabDec(); debugTabDec(); debugTabDec(); debugTabDec(); debugTabDec();
                                throw("NLM_CellNetwork::renderCellNetwork(): got exceptino about complex edge indicating fractional edge intersection value lambda outside ]0, 1[."\
                                    " internal logic error.");
                            }
                        }
                        debugTabDec();

                        /* compute split points */
                        uint32_t const n = e_info.edge_lambdas.size();
                        std::vector<R> e_split_points(n + 1);

                        e_split_points[0] = e_info.edge_lambdas[0] / 2.0;
                        debugl(3, "split_points[0] = %5.4f\n", e_split_points[0]);
                        for (uint32_t i = 1; i < n; i++) {
                            e_split_points[i] = (e_info.edge_lambdas[i] + e_info.edge_lambdas[i - 1]) / 2.0;
                            debugl(3, "split_points[%d] = %5.4f\n", i, e_split_points[i]);
                        }
                        e_split_points[n] = (1.0 + e_info.edge_lambdas[n-1]) / 2.0;
                        debugl(3, "split_points[%d] = %5.4f\n", n, e_split_points[n]);

                        /* split in red mesh, i.e. M_cell */
                        if (e_info.red) {
                            auto u_it = M_cell.vertices.find(e_info.u_id);
                            auto v_it = M_cell.vertices.find(e_info.v_id);

                            if (u_it != M_cell.vertices.end() && v_it != M_cell.vertices.end()) {
                                M_cell.splitEdge(u_it, v_it, e_split_points);
                            }
                            else {
                                debugTabDec(); debugTabDec(); debugTabDec(); debugTabDec(); debugTabDec();
                                throw("NLM_CellNetwork::renderCellNetwork(): discovered invalid vertex id contained "\
                                    " in complex edge information of RedBlue_Ex_ComplexEdges (id not found). internal"\
                                    " logic error.");
                            }
                        }
                        /* split in blue mesh, i.e. M_P */
                        else {
                            auto u_it = M_P.vertices.find(e_info.u_id);
                            auto v_it = M_P.vertices.find(e_info.v_id);

                            if (u_it != M_P.vertices.end() && v_it != M_P.vertices.end()) {
                                M_P.splitEdge(u_it, v_it, e_split_points);
                            }
                            else {
                                debugTabDec(); debugTabDec(); debugTabDec(); debugTabDec(); debugTabDec();
                                throw("NLM_CellNetwork::renderCellNetwork(): discovered invalid vertex id contained "\
                                    " in complex edge information of RedBlue_Ex_ComplexEdges (id not found). internal"\
                                    " logic error.");
                            }
                        }
                    }
                    debugTabDec();

                    /* all splits performed. at this point, new_outer_iteration == false, break_inner_meshing_loop == false =>
                     * another iteration of the inner meshing loop is performed, which issues another RedBlueUnion call
                     * on the same initial mesh segment after splitting complex edges. */
                    /*
                    tmp = M_cell;
                    tmp.writeObjFile("M_cell_split");

                    tmp = M_P;
                    tmp.writeObjFile("M_P_split");
                    */
                }
                catch (RedBlue_Ex_NumericalEdgeCase& numerical_ex) {
                    debugl(0, "NLM_CellNetwork::renderCellNetwork(): RedBlueAlgorithm returned exception: numerical edge case => retry..\n");
                    new_outer_iteration = true;
                    restore_M_cell      = !numerical_ex.R_intact;
                }
                catch (RedBlue_Ex_Triangulation<R>& tri_ex) {
                    debugl(0, "NLM_CellNetwork::renderCellNetwork(): RedBlueAlgorithm returned exception: error during triangulation of outside / inside polygons. => retry..\n");

                    /* decrease radius factor, but lower bound by radius_factor_safe_lb. */
                    radius_factor       = std::max(radius_factor_safe_lb, radius_factor - this->meshing_radius_factor_decrement);
                    new_outer_iteration = true;
                    restore_M_cell      = !tri_ex.R_intact;
                }
                catch (RedBlue_Ex_NumIsecPoly& isecpoly_ex) {
                    debugl(0, "NLM_CellNetwork::renderCellNetwork(): RedBlueAlgorithm returned exception: number of intersection polygons != 1.\n");
                    radius_factor       = std::max(radius_factor_safe_lb, radius_factor - this->meshing_radius_factor_decrement);
                    new_outer_iteration = true;
                    restore_M_cell      = !isecpoly_ex.R_intact;
                }
                catch (RedBlue_Ex_AffectedCircleTrivial& trivcircle_ex) {
                    debugTabDec(); debugTabDec(); debugTabDec(); debugTabDec();
                    throw;
                }
                debugl(1, "inner meshing loop time: %5.4f\n\n", Aux::Timing::tack(14));
            }
            debugTabDec();
            debugl(1, "inner meshing loop left..\n");

            /* if maximum number of inner meshing loop iterations has been reached, restart outer meshing loop. */
            if (inner_loop_iter == this->meshing_inner_loop_maxiter) {
                debugl(0, "inner meshing loop broken because maximum number of iterations has been reached => restart outer meshing loop.\n");
                radius_factor   = std::max(radius_factor_safe_lb, radius_factor - this->meshing_radius_factor_decrement);
                new_outer_iteration = true;
                restore_M_cell      = true;
            }

            /* start fresh iteration of outer meshing loop if required */
            if (new_outer_iteration) {
                debugl(1, "re-iteration of outer meshing loop necessary..\n");
                /* if radius_factor has reached radius_factor_safe_lb, throw exception, since a definitely safe radius
                 * should already have been reached. */
                if (radius_factor == radius_factor_safe_lb) {
                    debugTabDec(); debugTabDec(); debugTabDec();
                    throw("NLM_CellNetwork::renderCellNetwork(): reached safe lower bound radius factor for current neurite path. this must not happen for clean cell networks.");
                }
                /* otherwise start another meshing run */
                else {
                    continue;
                }
            }

            /* if this line is reached, the initial segment has been successfully merged into M_cell, where
             * circle_its_update contain iterators to the end circle and closing vertex as new vertices of the partial
             * union mesh M_cell.  check if any iterator has been explicitly invalidated by the RedBlueUnion algorithm.
             * if so, throw an exception. if the network is clean in the sense of the definition in the thesis, this
             * must not happen.  otherwise, unpack the updated iterators back into end_circle_its and closing_vertex_it
             * and append rest of path to M_cell. */
            for (auto &it : circle_its_update) {
                if (it.explicitlyInvalid()) {
                    debugTabDec(); debugTabDec(); debugTabDec();
                    throw("NLM_CellNetwork::renderCellNetwork(): RedBlueUnion algorithm has explicitly invalidated an "\
                        "end circle iterator or the closing vertex iterator of the current path P's initial mesh "
                        "segment. in a clean network, this should be impossible. numerical edge case due to tight "
                        "PMDV / SMDV constants?");

                }
                else if (!it.checkContainer(M_cell)) {
                    debugl(1, "it.container: %p, M_cell (ptr): %p, M_P (ptr): %p.\n", it.getContainer(), &M_cell, &M_P);
                    debugTabDec(); debugTabDec(); debugTabDec();
                    throw("NLM_CellNetwork::renderCellNetwork(): RedBlueUnion algorithm has returned an "\
                        "updated end circle iterator that does not refer to the partially completed cell mesh. "\
                        "internal logic error.");
                }
            }

            debugl(2, "initial mesh segment successfully merged. appending path tail mesh..\n");

            /* unpack updated iterators referring to M_cell. */
            closing_vertex_it   = circle_its_update.back();
            circle_its_update.pop_back();
            end_circle_its      = circle_its_update;

            /* append P's tail mesh (for neurite canal segments 1, .., m) to M */
            try
            {
                P.template appendTailMesh<Tm, Tv, Tf>(
                    /* append to the partial cell mesh M_cell */
                    M_cell,
                    /* n_phi_segments default to 16 for testing */
                    this->meshing_canal_segment_n_phi_segments,
                    /* render vector */
                    render_vector,
                    /* phi_0, arclen_dt = 1E-3 */
                    phi_0,
                    1E-3,
                    /* end circle information from RedBlue merged initial segment as start circle information for tail */
                    end_circle_offset,
                    end_circle_its,
                    closing_vertex_it,
                    /* store return iterators for subsequent generation of terminal half-sphere */
                   &end_circle_offset,
                   &end_circle_its,
                   &closing_vertex_it,
                    this->meshing_preserve_crease_edges);
            }
            catch (...) {debugTabDec(); debugTabDec(); debugTabDec(); throw;}

            debugl(2, "tail path mesh appended. appending terminal half-sphere.\n");

            /* append a terminal "half-sphere" at the end neurite point of P */
            BLRCanalSurface<3u, R> &C_end   = *(P.canal_segments_magnified.back());
            Vec3<R> start               = C_end.spineCurveEval(1.0);
            Vec3<R> direction           = C_end.spineCurveEval_d(1.0);
            R       radius              = C_end.radiusEval(1.0);

            MeshAlg::appendHalfSphereToCanalSurfaceMesh<Tm, Tv, Tf, R>(
                    M_cell,
                    render_vector,
                    start,
                    radius,
                    direction,
                    this->meshing_canal_segment_n_phi_segments,
                    phi_0,
                    end_circle_its,
                    closing_vertex_it);

            debugl(2, "half-sphere appended. triangulating quads..\n");

            /* triangulate quads in M_cell */
            M_cell.triangulateQuads();

            debugl(1, "path %d completely processed. M_cell.numVertices(): %d\n", (*npt_vit)->id(), M_cell.numVertices());

            /* done for path P */
            done = true;
        }
        debugTabDec();
        np_idx++;
    }
    debugTabDec();
    debugl(1, "all neurite paths processed. finalizing obj file..\n");

    /* select all faces from cell mesh and flush them .. if no flush has been performed before, this is semantically
     * equivalent to writeObjFile(), otherwise it completes partially flushed cell meshes that are yet incomplete in the
     * flush obj file. */
    std::list<typename Mesh<Tm, Tv, Tf, R>::Face *> remaining_faces = {};
    M_cell.invertFaceSelection(remaining_faces);
    try {MeshAlg::partialFlushToObjFile(M_cell, M_cell_flushinfo, remaining_faces);}
    catch (...) {debugTabDec(); debugTabDec(); throw;}

    debugTabDec();
    debugl(1, "NLM_CellNetwork<R>::renderCellNetwork(): done.\n");
}

template <typename R>
template <typename Tm, typename Tv, typename Tf>
void
NLM_CellNetwork<R>::renderModellingMeshesIndividually(std::string filename) const
{
    Mesh<Tm, Tv, Tf, R>                                     M_cell, M_S, M_P;
    bool                                                    initial_segment_end_circle_offset;
    std::vector<
            typename Mesh<Tm, Tv, Tf, R>::vertex_iterator
        >                                                   initial_segment_end_circle_its;
    typename Mesh<Tm, Tv, Tf, R>::vertex_iterator           initial_segment_closing_vertex_it;

    /* iterate over all somas */
    for (auto &s : this->soma_vertices) {
        /* generate mesh for soma sphere and append to cell mesh. */
        NLM::SomaInfo<R> const &s_info  = s.soma_data;
        s_info.soma_sphere.template generateMesh<Tm, Tv, Tf>(M_S);
        M_cell.moveAppend(M_S);
        
        /* iterate over all neurite path trees and append meshes from all neurite paths */
        for (auto &npt : s_info.neurite_path_trees) {
            for (auto &npt_v : npt.vertices) {
                NLM::NeuritePath<R> const  &P = *npt_v;
                Vec3<R>                     render_vector;

                /* find permissible render vector for neurite path */
                render_vector = P.findPermissibleRenderVector();

                /* clear path mesh */
                M_P.clear();

                /* generate P's initial segment mesh and append to M */
                P.template generateInitialSegmentMesh<Tm, Tv, Tf>(
                    /* append to mesh M_P for path P*/
                    M_P,
                    /* n_phi_segments default to 16 for testing */
                    16, 
                    /* render vector */
                    render_vector,
                    /* phi_0 = 0, arclen_dt = 1E-3 */
                    0,
                    1E-3,
                    /* end circle info */
                    initial_segment_end_circle_offset,
                    initial_segment_end_circle_its,
                    initial_segment_closing_vertex_it,
                    /* no radius reduction => factor 1 */
                    1);

                /*
                Mesh<Tm, Tv, Tf, R> M_initial = M_P;
                M_initial.writeObjFile("M_initial", false);
                */

                /* append P's tail mesh (for neurite canal segments 1, .., m) to M */
                P.template appendTailMesh<Tm, Tv, Tf>(
                    /* append to mesh M_P for path P */
                    M_P,
                    /* n_phi_segments default to 16 for testing */
                    16, 
                    /* render vector */
                    render_vector,
                    /* phi_0 = 0, arclen_dt = 1E-3 */
                    0,
                    1E-3,
                    /* end circle information from initial segment as start circle information for tail */
                    initial_segment_end_circle_offset,
                    initial_segment_end_circle_its,
                    initial_segment_closing_vertex_it,
                    /* no return information required */
                    NULL,
                    NULL,
                    NULL);

                M_cell.moveAppend(M_P);
            }
        }
    }

    /* output M */
    M_cell.triangulateQuads();
    M_cell.writeObjFile(filename.c_str());
}

template <typename R>
void
NLM_CellNetwork<R>::writeMorphViewFile(std::string filename) const
{
    debugl(1, "NLM_CellNetwork::writeMorphViewFile(): writing to filename \"%s\".\n", filename.c_str() );
    debugTabInc();

    FILE *f = fopen(filename.c_str(), "w");
    if (!f) {
        throw("NLM_CellNetwork()::writeBinaryDump(): can't open output filename.");
    }

    std::list<NLM::NeuritePath<R> const *>  neurite_paths;
    this->getAllNeuritePaths(neurite_paths);

    uint8_t flag;

    /* write total number of neurite vertices, coordinates and radii */
    uint32_t num_nv = this->neurite_vertices.size();
    fwrite(&num_nv, sizeof(uint32_t), 1, f);

    debugl(1, "writing position / radius / branching flag for %d neurite vertices.\n", num_nv);

    /* iterate over all neurite vertices */
    Vec3<R> nv_pos;
    R       nv_r, tmp;
    uint8_t nv_branch;
    for (auto &nv : this->neurite_vertices) {
        nv_pos      = nv.getPosition();
        nv_r        = nv.getRadius();
        nv_branch   = nv.isNeuriteBranchingVertex();

        fwrite(&nv_pos[0],  sizeof(R), 1, f);
        fwrite(&nv_pos[1],  sizeof(R), 1, f);
        fwrite(&nv_pos[2],  sizeof(R), 1, f);
        fwrite(&nv_r,       sizeof(R), 1, f);
        fwrite(&nv_branch,  sizeof(uint8_t), 1, f);
    }

    /* write number of neurite paths */
    uint32_t num_neurite_paths = neurite_paths.size();
    fwrite(&num_neurite_paths, sizeof(uint32_t), 1, f);

    debugl(1, "writing information about %d neurite paths.\n", num_neurite_paths);

    /* iterate over all neurite paths */
    uint32_t    nP, i, j, k, degree;
    Vec3<R>     bb_min, bb_max, bb_len;
    debugTabInc();
    for (auto &P : neurite_paths) {
        debugl(2, " ----- Dumping Path ----- \n");

        /* write number of neurite segment curves of P in file */
        nP                  = P->numVertices();
        uint32_t tmp_int    = nP - 1;
        fwrite(&tmp_int, sizeof(uint32_t), 1, f);

        /* write the nP parameter values P's vertices to file */
        for (i = 0; i < nP; i++) {
            fwrite(&P->neurite_vertex_parameters[i], sizeof(R), 1, f);
        }

        /* now the radii of the nP vertices */
        for (i = 0; i < nP; i++) {
            tmp = P->operator[](i)->getRadius();
            fwrite(&tmp, sizeof(double), 1, f);
        }

        /* the bounding box of the whole path */
        auto P_bb = P->getBoundingBox();
        bb_min  = P_bb.min();
        bb_max  = P_bb.max();
        bb_len  = bb_max - bb_min;

        /* minimum vector of the bounding box */
        fwrite(&bb_min[0], sizeof(R), 1, f);
        fwrite(&bb_min[1], sizeof(R), 1, f);
        fwrite(&bb_min[2], sizeof(R), 1, f);

        /* bounding box length in x, y and z */
        fwrite(&bb_len[0], sizeof(R), 1, f);
        fwrite(&bb_len[1], sizeof(R), 1, f);
        fwrite(&bb_len[2], sizeof(R), 1, f);

        /* write all spine curve component polynomial coefficients, bounding boxes and flags of neurite segment curves
         * */
        debugTabInc();
        for (i = 0; i < nP - 1; i++) {
            /* get const reference to neurite segment curve and neurite segment info */
            neurite_segment_const_iterator ns_it        = P->neurite_segments[i];
            NLM::NeuriteSegmentInfo<R> const &ns_info   = ns_it->neurite_segment_data;
            BLRCanalSurface<3u, R> const &Gamma             = ns_info.canal_segment_magnified;
            BezierCurve<3u, R> const &gamma                 = Gamma.getSpineCurve();

            /* dump coefficients of component functions to file */
            /* for all three component functions */
            debugTabInc();
            for (j = 0; j < 3; j++) {
                /* write degree of component function j, which will all be equal in general, but
                 * this allows more flexibility */
                degree = gamma[j].getDegree();
                fwrite(&degree, sizeof(uint32_t), 1, f);

                /* now write the (degree + 1) coefficients of the component function j */
                for (k = 0; k < degree + 1; k++) {
                    tmp = gamma[j](k);
                    fwrite(&tmp, sizeof(double), 1, f);
                }
            }
            debugTabDec();

            /* write bounding box of neurite segment: min vector and length */
            auto Gamma_bb = Gamma.getBoundingBox();
            bb_min  = Gamma_bb.min();
            bb_max  = Gamma_bb.max();
            bb_len  = bb_max - bb_min;

            fwrite(&bb_min[0], sizeof(R), 1, f);
            fwrite(&bb_min[1], sizeof(R), 1, f);
            fwrite(&bb_min[2], sizeof(R), 1, f);

            /* bounding box length in x, y and z */
            fwrite(&bb_len[0], sizeof(R), 1, f);
            fwrite(&bb_len[1], sizeof(R), 1, f);
            fwrite(&bb_len[2], sizeof(R), 1, f);

            /* write all flags */
            flag    = ns_info.pmdv;
            fwrite(&flag, sizeof(uint8_t), 1, f);

            flag    = ns_info.smdv;
            fwrite(&flag, sizeof(uint8_t), 1, f);

            flag    = ns_info.lsi;;
            fwrite(&flag, sizeof(uint8_t), 1, f);

            flag    = ns_info.gsi;;
            fwrite(&flag, sizeof(uint8_t), 1, f);

            flag    = ns_info.ic_sons;
            fwrite(&flag, sizeof(uint8_t), 1, f);

            flag    = ns_info.icrn_nsns;
            fwrite(&flag, sizeof(uint8_t), 1, f);

            flag   = ns_info.icin_nsns;
            fwrite(&flag, sizeof(uint8_t), 1, f);
        }
        debugTabDec();
    }
    debugTabDec();


    /* write information about first soma, only applicable for cells from NeuroMorpho where there's only one soma.. */
    soma_const_iterator s       = this->soma_vertices.begin();
    NLM::SomaInfo<R>    s_info  = s->soma_data;
    Vec3<R> s_sphere_centre     = s_info.soma_sphere.centre();
    R s_sphere_radius           = s_info.soma_sphere.radius();

    auto s_bb                   = s_info.soma_sphere.getBoundingBox();
    bb_min                      = s_bb.min();
    bb_max                      = s_bb.max();
    bb_len                      = bb_max - bb_min;

    debugl(2, "writing soma sphere position vector / radius / bounding box for first soma %d\n", s->id());

    /* write soma centre, radius and bounding box */
    fwrite(&s_sphere_centre[0], sizeof(R), 1, f);
    fwrite(&s_sphere_centre[1], sizeof(R), 1, f);
    fwrite(&s_sphere_centre[2], sizeof(R), 1, f);
    fwrite(&s_sphere_radius,    sizeof(R), 1, f);

    /* minimum vector of the bounding box */
    fwrite(&bb_min[0], sizeof(R), 1, f);
    fwrite(&bb_min[1], sizeof(R), 1, f);
    fwrite(&bb_min[2], sizeof(R), 1, f);

    /* bounding box length in x, y and z */
    fwrite(&bb_len[0], sizeof(R), 1, f);
    fwrite(&bb_len[1], sizeof(R), 1, f);
    fwrite(&bb_len[2], sizeof(R), 1, f);

    fclose(f);

    debugTabDec();
    debugl(1, "NLM_CellNetwork::writeMorphViewFile(): done writing file \"%s\".\n", filename.c_str() );
}
