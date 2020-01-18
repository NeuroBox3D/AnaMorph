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

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax,typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, // typename Tsyn = Common::UnitType
    typename R
>
void
preliminaryPreconditioning(
    CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr>  &C,
    R const                                                                    &alpha,
    R const                                                                    &beta,
    R const                                                                    &gamma)
{
    debugl(1, "CellNetwork::preliminaryPreconditioning(): alpha: %2.3f, beta: %2.3f, gamma: %2.3f.\n", alpha, beta, gamma);
    debugTabInc();

    typedef CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr> CellNetworkType;

    /* first, get all neurite segments ns whose length is longer than gamma * maxRadius(ns) */
    debugl(1, "extracting all \"gamma\"-long neurite segments.\n");
    std::list<typename CellNetworkType::neurite_segment_iterator>   long_neurite_segments;
    debugTabInc();
    for (auto &ns : C.neurite_segments) {
        if (ns.getLength() > gamma * ns.getMaxRadius()) {
            debugl(2, "neurite segment: %d: len: %5.4f, gamma * ns.getMaxRadius(): %5.4f\n",
                    ns.id(), ns.getLength(), gamma*ns.getMaxRadius());

            long_neurite_segments.push_back(ns.iterator());
        }
    }
    debugTabDec();

    debugl(1, "splitting all \"gamma\"-long neurite segments ..\n");
    /* now split every long neurite segment ns individually with the following algorithm:
     *
     *  1. choose an n >= 2 (see 2.) and split ns into n neurite segments of equal length with linearly interpolated
     *  radii
     *
     *  2. n is chosen as follows: for given value k of segments to be created, calculate a quadratic penalty function
     *  that is minimized whenever all k neurite segments (m_1, .., m_k) all exactly satisfy:
     *
     *      len(m_i) = gamma*maxRadius(m_i)
     *
     *  the quadratic error function is calculated as follows:
     *
     *      sum_{i = 0}^k{ (len(m_i) - gamma*maxRadius(m_i))^2 }.
     *
     *  with an initial upper bound nmax, all values n = 2,.., nmax are tested and the value n with minimum quadratic
     *  error function is chosen.  */
    typename CellNetworkType::neurite_iterator  u, v;
    Vec3<R>                                     p_u, p_v;
    R                                           uv_len, r_u, r_v, rmax;
    uint32_t                                    k, nmax, m;
    R                                           k_penalty, n_penalty, ratio;
    std::vector<std::pair<Vec3<R>, R>>          k_vertex_info, n_vertex_info;

    auto
    penalty_function = [gamma] (std::vector< std::pair<Vec3<R>, R>> vertices) -> R
    {
        /* sum individual penalties over all neurite segments.. */
        uint32_t const  n = vertices.size();
        R               p, penalty_total = 0, len, rmax;

        for (uint32_t i = 0; n > 0 && i < n - 1; i++) {
            len             = (vertices[i + 1].first - vertices[i].first).len2();
            rmax            = std::max(vertices[i+1].second, vertices[i].second);
            p               = (len - gamma * rmax);
            penalty_total  += p*p;
        }

        return penalty_total;
    };

    debugTabInc();
    for (auto &lns : long_neurite_segments) {
        debugl(2, "handling \"gamma\"-long neurite segment %d\n", lns->id() );
        debugTabInc();

        uv_len  = lns->getLength();
        u       = lns->getSourceVertex();
        v       = lns->getDestinationVertex();
        p_u     = u->getPosition();
        p_v     = v->getPosition();
        r_u     = u->getRadius();
        r_v     = v->getRadius();
        rmax    = std::max(r_u, r_v);

        debugl(2, "len: %5.4f, rmax: %5.4f\n", uv_len, rmax);

        /* calculate upper bound: nmax is the minimum number of segments to be created such that the length of one
         * segment falls beneath rmax => 
         *
         *  uv_len / nmax < rmax where uv_len / (nmax - 1) > rmax =>
         *
         *  nmax = ceil(uv_len / rmax)
         *
         *  */
        nmax        = std::ceil(uv_len / rmax);
        n_penalty   = Aux::Numbers::inf<R>();

        debugl(2, "nmax: %d\n", nmax);

        debugl(2, "testing possible values for n..\n");
        /* test all possible values 2, .., nmax */
        debugTabInc();
        for (k = 2; k <= nmax; k++) {
            debugl(3, "testing n = %d.\n", k);
            /* construct vertex info vector: neurite segment is partitioned into k segments => k + 1 vertices, the first
             * and last of which are u and v, respectively. */
            k_vertex_info.resize(k + 1);  
            k_vertex_info[0] = { p_u, r_u };
            k_vertex_info[k] = { p_v, r_v };

            for (m = 1; m < k; m++) {
                ratio = (R)m / (R)k;

                k_vertex_info[m] = {
                        p_u + (p_v - p_u) * ratio,
                        r_u + (r_v - r_u) * ratio
                    };
            }

            /* calculate penalty */
            k_penalty = penalty_function(k_vertex_info);
            if (k_penalty < n_penalty) {
                n_penalty       = k_penalty;
                n_vertex_info   = k_vertex_info;
            }
        }
        debugTabDec();

        /* value for n has been chosen. split the segment using n_vertex_info, but remove u (front) and v (back) to get
         * information about intermediate vertices. */
        n_vertex_info.erase(n_vertex_info.begin());
        n_vertex_info.pop_back();

        C.neurite_segments.split(lns, n_vertex_info);

        debugTabDec();
    }
    debugTabDec();

    /* greedy neurite segment collapsing with priority queue indexing on reduced neurite segment length */
    PriorityQueue<double, uint32_t>             Q;
    auto
    calculate_ns_weight = [alpha, beta] (typename CellNetworkType::neurite_segment_const_iterator const &ns_it) -> R
    {
        R ns_len        = ns_it->getLength();            
        auto smdv_radii = ns_it->getSMDVRadii();

        return (
                std::min(
                    ns_len - alpha*ns_it->getMaxRadius(),
                    ns_len - beta*(smdv_radii.first + smdv_radii.second)
                )
            );
    };


    debugl(1, "performing greedy collapse fixed-point iteration.\n");
    debugTabInc();
    R           e_len, e_weight;
    uint32_t    e_id;
    bool        fixed_point = false;

    while (!fixed_point) {
        debugl(2, "performing another greedy collapse iteration. inserting all necessary neurite segments into Q..\n");
        fixed_point = true;
        Q.clear();

        /* using the length as key, insert all neurite segments into the initial Q for this iteration of the fixed-point
         * scheme. */
        for (auto &ns : C.neurite_segments) {
            Q.insert(ns.getLength(), ns.id());
        }

        debugl(2, "processing Q.\n");
        debugTabInc();
        while (!Q.empty()) {
            /* get neurite segment e = (u, v) with minimum key off priority queue Q */
            auto min_pair   = Q.top();
            e_len           = min_pair.first;
            e_id            = min_pair.second;
            Q.deleteMin();

            debugl(2, "processing neurite segment %d, key from Q, i.e. length: %5.4f\n", e_id, e_len);

            /* only process if e still exists and its key is < 0, i.e. if e is an alpha-PMDV or a beta-SMDV. although
             * this holds for all neurite segments in Q, this can change due to updates after collapses. */
            auto e_it = C.neurite_segments.find(e_id);
            debugTabInc();
            if (e_it != C.neurite_segments.end()) {
                /* check if key-value just dequeued indeed matches the length of e, since it must have been
                 * enqueued with this length and updated on every change. */
                if ( std::abs(e_it->getLength() - e_len) > 1E-5 ) {
                    debugl(1, "CellNetworkAlg::preliminaryPreconditioning(): current Q.top() ns: %d. length: %5.4f, "\
                        "key from Q: %5.4f => mismatch. internal logic error.\n",
                        e_id, e_it->getLength(), e_len);

                    throw("CellNetworkAlg::preliminaryPreconditioning(): currently considered top() neurite segment "\
                        "from Q has a length that differs from the key value it has been enqueued with. these should "\
                        "always be identical. internal logic error.");
                }

                /* act depending upon the types of e's vertices u and v */
                u           = e_it->getSourceVertex();
                v           = e_it->getDestinationVertex();

                /* get weight */
                e_weight    = calculate_ns_weight(e_it);

                debugl(2, "e = (u, v) | %d = (%d, %d). weight: %5.4f.\n", e_id, u->id(), v->id(), e_weight);

                /* only process e if weight is <= 0, i.e. e exhibits an alpha-PMDV or beta-SMDV */
                if (e_weight <= 0.0) {
                    /* variables for the vertex that e will (potentially) be collapsed into. */
                    Vec3<R>                                     v_collapsed_pos;
                    R                                           v_collapsed_radius;
                    typename CellNetworkType::neurite_iterator  v_collapsed_it;

                    /* u is neurite root vertex */
                    if (u->isNeuriteRootVertex()) {
                        /* neurite root vertices are constrained to the surface of the soma sphere, or at least its vicinity.
                         * u's position is not to be altered. can't collapse if neurite ends or branches directly in v
                         * (would lead to an immediately branching neurite directly "on the soma", which is unrealistic). */
                        if (v->isNeuriteTerminalVertex() || v->isNeuriteBranchingVertex()) {
                            debugl(3, "u is neurite root vertex, v is {terminal, branching} vertex. skipping e..\n");
                            debugTabDec();
                            continue;
                        }
                        else if (v->isNeuriteSimpleVertex()) {
                            /* delete v, i.e. collapse e into u. */
                            v_collapsed_pos     = u->getPosition();
                            v_collapsed_radius  = u->getRadius();
                        }
                        /* this should never happen */
                        else {
                            throw("CellNetworkAlg::preliminaryPreconditioning(): logically impossible case for "\
                                "destination vertex v of currently considered neurite segment e = (u, v). internal logic "\
                                "error.");
                        }
                    }
                    /* u is a neurite branching vertex */
                    else if (u->isNeuriteBranchingVertex()) {
                        /* can't collapse or delete if v is a terminal vertex */
                        if (v->isNeuriteTerminalVertex()) {
                            debugl(3, "u is neurite branching vertex, v is terminal vertex. skipping e..\n");
                            debugTabDec();
                            continue;
                        }
                        /* branching vertex */
                        else if (v->isNeuriteBranchingVertex()) {
                            /* collapse e */
                            v_collapsed_pos     = (u->getPosition() + v->getPosition()) * 0.5;
                            v_collapsed_radius  = (u->getRadius() + v->getRadius()) * 0.5;
                        }
                        else if (v->isNeuriteSimpleVertex()) {
                            /* delete v, i.e. collapse e into u. */
                            v_collapsed_pos     = u->getPosition();
                            v_collapsed_radius  = u->getRadius();
                        }
                        /* this should never happen */
                        else {
                            throw("CellNetworkAlg::preliminaryPreconditioning(): logically impossible case for "\
                                "destination vertex v of currently considered neurite segment e = (u, v). internal logic "\
                                "error.");
                        }
                    }
                    /* otherwise u is a simple neurite vertex on interior of some neurite path. */
                    else if (u->isNeuriteSimpleVertex()) {
                        if (v->isNeuriteTerminalVertex() || v->isNeuriteBranchingVertex()) {
                            /* delete u, i.e. collapse e into v */
                            v_collapsed_pos     = v->getPosition();
                            v_collapsed_radius  = v->getRadius();
                        }
                        else if (v->isNeuriteSimpleVertex()) {
                            /* collapse e */
                            v_collapsed_pos     = (u->getPosition() + v->getPosition()) * 0.5;
                            v_collapsed_radius  = (u->getRadius() + v->getRadius()) * 0.5;
                        }
                        else {
                            throw("CellNetworkAlg::preliminaryPreconditioning(): logically impossible case for "\
                                "destination vertex v of currently considered neurite segment e = (u, v). internal logic "\
                                "error.");
                        }
                    }
                    else {
                        throw("CellNetworkAlg::preliminaryPreconditioning(): logically impossible case for "\
                            "source vertex u of currently considered neurite segment e = (u, v). internal logic "\
                            "error.");
                    }

                    debugl(2, "collapsing e..\n");
                    /* perform the collapse, store return iterator. */
                    v_collapsed_it = C.neurite_segments.collapse(e_it, { v_collapsed_pos, v_collapsed_radius } );
                    debugl(2, "collapse performed. collapsed vertex: %d\n", v_collapsed_it->id());

                    /* set fixed_point to false, since at least one additional iteration is required. */
                    fixed_point = false;

                    /* the neurite segment has been collapsed in one way or the other and v_collapsed_it points to the
                     * vertex into which the neurite segment has been collapsed. inspect all edges incident to v_collapsed,
                     * change their keys or insert them into Q if necessary. */
                    debugl(2, "updating edge-neighbours of collapsed e / collapsed vertex.\n");
                    std::list<typename CellNetworkType::NeuriteSegment *> v_collapsed_incident_edges, tmp;

                    v_collapsed_it->template getFilteredInEdges<typename CellNetworkType::NeuriteSegment>(v_collapsed_incident_edges);
                    v_collapsed_it->template getFilteredOutEdges<typename CellNetworkType::NeuriteSegment>(tmp);
                    v_collapsed_incident_edges.insert(v_collapsed_incident_edges.end(), tmp.begin(), tmp.end());

                    debugTabInc();
                    for (auto f : v_collapsed_incident_edges) {
                        debugl(3, "updating edge-neighbour %d\n", f->id());
                        R f_weight = calculate_ns_weight(f->iterator());
                        
                        /* NOTE: moderately subtle lazy evaluation here.. if changeKey() returns false (<=> value not
                         * present) AND calculated weight of f is <= 0, (re)insert f into Q. */
                        if (!Q.changeKey(f->id(), f->getLength()) && f_weight <= 0) {
                            debugl(3, "edge-neighbour %d: not in Q and weight %5.4f <= 0 => inserting with key (new length) %5.4f.\n",
                                    f->id(), f_weight, f->getLength());

                            Q.insert(f->getLength(), f->id());
                        }
                        else {
                            debugl(3, "edge-neighbour %d: found in Q, key updated to %5.4f.\n", f->id(), f->getLength());
                        }
                    }
                    debugTabDec();
                }
                else {
                    debugl(2, "neurite segment %d not relevant for processing.\n", e_id);
                }
            }
            else {
                debugl(2, "neurite segment %d no longer existent.\n");
            }
            debugTabDec();
        }
        debugTabDec();
    }
    debugTabDec();
    debugl(1, "fixed point reached. returning..\n");

    debugTabDec();
    debugl(1, "CellNetwork::preliminaryPreconditioning(): done.\n");
}



template <typename network_type, typename R>
void scale_radii(network_type& C, const R& scale)
{
    // scaling
    typename network_type::neuron_iterator nvit = C.neuron_vertices.begin();
    typename network_type::neuron_iterator nvit_end = C.neuron_vertices.end();
    for (; nvit != nvit_end; ++nvit)
        nvit->scaleSinglePointRadius(scale);

    // adjust root vertex positions to preserve intersection with soma
    C.computeInitialNeuriteRootVertices();
}

